// Copyright 2019 Google LLC
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//     https://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.
// -----------------------------------------------------------------------------
//
// WP2 residual encoding.
//

#include <limits>

#include "src/enc/wp2_enc_i.h"

namespace WP2 {

WP2Status ResidualWriter::Init(bool use_aom_coeffs, bool has_alpha) {
  num_channels_ = (has_alpha ? 4 : 3);
  aom_writer_.Init();
  ResidualIO::Init(use_aom_coeffs);
  return WP2_STATUS_OK;
}

WP2Status ResidualWriter::WriteHeaderForResidualSymbols(
    Channel channel, uint32_t num_coeffs_max, const SymbolRecorder& recorder,
    ANSEncBase* const enc, SymbolWriter* const sw,
    ANSDictionaries* const dicts_in) {
  bool is_maybe_used[kResidual1Max + 1];
  for (EncodingMethod method :
       {EncodingMethod::kMethod0, EncodingMethod::kMethod1}) {
    std::string label;

    // Write the dictionaries for the number of consecutive zeros.
    const uint32_t method_index = GetMethodIndex(method);
    (void)method_index;
    WP2SPrint(&label, "%s_num_zeros_%d", kChannelStr[channel], method_index);
    const uint32_t cluster = GetCluster(channel, num_channels_, method);
    WP2_CHECK_STATUS(sw->WriteHeader(cluster, num_coeffs_max,
                                     kSymbolResNumZeros, recorder,
                                     label.c_str(), enc, dicts_in));

    // Write the dictionaries for small residuals.
    WP2SPrint(&label, "%s_bits0_%d", kChannelStr[channel], method_index);
    WP2_CHECK_STATUS(sw->WriteHeader(cluster, num_coeffs_max, kSymbolBits0,
                                     recorder, label.c_str(), enc, dicts_in));
    sw->GetPotentialUsage(cluster, kSymbolBits0, is_maybe_used,
                          kResidual1Max + 1);

    if (is_maybe_used[kResidual1Max]) {
      // Write the dictionaries for prefixes of big residuals.
      WP2SPrint(&label, "%s_bits1_%d", kChannelStr[channel], method_index);
      WP2_CHECK_STATUS(sw->WriteHeader(cluster, num_coeffs_max, kSymbolBits1,
                                       recorder, label.c_str(), enc, dicts_in));
    } else {
      const WP2::ANSDictionary& d =
          recorder.GetRecordedDict(cluster, kSymbolBits1);
      (void)d;
      assert(*std::max_element(d.Counts().begin(), d.Counts().end()) == 0);
    }
  }
  return WP2_STATUS_OK;
}

WP2Status ResidualWriter::WriteHeader(
    uint32_t num_coeffs_max_y, uint32_t num_coeffs_max_uv,
    uint32_t num_transforms, bool has_lossy_alpha,
    const SymbolRecorder& recorder, ANSDictionaries* const dicts,
    ANSEncBase* const enc, SymbolWriter* const sw) {
  if (use_aom_coeffs_) {
    for (Symbol sym : kSymbolsForAOMCoeffs) {
      WP2_CHECK_STATUS(sw->WriteHeader(num_transforms, sym, recorder,
                                       "aom_symbols", enc, dicts));
    }
  } else {
    for (Channel channel : {kYChannel, kUChannel, kVChannel, kAChannel}) {
      if (channel == kAChannel && !has_lossy_alpha) continue;

      WP2_CHECK_STATUS(WriteHeaderForResidualSymbols(
          channel,
          (channel == kYChannel || channel == kAChannel) ? num_coeffs_max_y
                                                         : num_coeffs_max_uv,
          recorder, enc, sw, dicts));
    }

    for (Symbol sym : kSymbolsForCoeffs) {
      // For U/V, num_coeffs_max_y should be num_coeffs_max_uv instead, but in
      // practice it makes no difference and this is simpler.
      WP2_CHECK_STATUS(sw->WriteHeader(num_coeffs_max_y, sym, recorder,
                                       "residual_symbols", enc, dicts));
    }

    for (Symbol sym : kSymbolsForCoeffsPerTf) {
      WP2_CHECK_STATUS(sw->WriteHeader(num_transforms, sym, recorder,
                                       "block_symbols", enc, dicts));
    }
  }

  if (has_lossy_alpha) {
    WP2_CHECK_STATUS(sw->WriteHeader(num_transforms, kSymbolEncodingMethodA,
                                     recorder, "coeff_method_alpha", enc,
                                     dicts));
  }
  WP2_CHECK_STATUS(sw->WriteHeader(num_transforms, kSymbolHasCoeffs, recorder,
                                   "has_coeffs", enc, dicts));
  WP2_CHECK_STATUS(sw->WriteHeader(num_transforms, kSymbolEncodingMethod,
                                   recorder, "coeff_method", enc, dicts));

  return WP2_STATUS_OK;
}

void ResidualWriter::StoreDC(Channel channel, uint32_t num_channels,
                             ANSEncBase* const enc, SymbolManager* const sm,
                             int16_t v, bool can_be_zero) {
  if (!can_be_zero) assert(v != 0);
  // Converting the signed value to an unsigned one.
  uint32_t n = (v >= 0) ? 2 * v : -2 * v - 1;
  if (!can_be_zero) --n;
  sm->Process(GetCluster(channel, num_channels), kSymbolDC, n, "DC", enc);
}

inline void WriteBounds(Channel channel, uint32_t num_channels,
                        EncodingMethod method, TrfSize tdim, bool is_x_first,
                        uint32_t val1, uint32_t min1, uint32_t max1,
                        bool use_bounds2, uint32_t val2, ANSEncBase* const enc,
                        SymbolManager* const sm) {
  enc->PutRange(val1, min1, max1, "bound1");
  uint8_t min2, max2;
  if (is_x_first) {
    ResidualBoxAnalyzer::GetRangePerX(tdim, val1, &min2, &max2);
  } else {
    ResidualBoxAnalyzer::GetRangePerY(tdim, val1, &min2, &max2);
  }
  if (min2 != 255 && min2 <= max2 &&
      sm->Process(ResidualWriter::GetClusterMergedUV(channel, num_channels,
                                                     method, tdim),
                  kSymbolResidualUseBound2, use_bounds2, "use_bound2", enc)) {
    enc->PutRange(val2, min2, max2, "bound2");
  }
}

// 'dicts' is a pointer to an array of pointers to dictionary-like structs.
// The residuals are stored in 'coeffs' of size 'num_coeffs'.
// TODO(vrabaud) 'is_pseudo_rate' should be removed but results are slightly
// worse (on an RD-curve) without it. We should have better default values
// for ANSBinSymbol for the bound checking as those are usually very biased
// (in comparison, the ANSBinSymbol for is_zero and such are very close to 1).
void ResidualWriter::StoreCoeffs(const int16_t* const coeffs,
                                 uint32_t num_coeffs, bool first_is_dc,
                                 TrfSize tdim, Channel channel,
                                 uint32_t num_channels, EncodingMethod method,
                                 ANSEncBase* const enc, SymbolManager* const sm,
                                 bool is_pseudo_rate) {
  // Figure out the rectangle in which the residuals fit.
  uint32_t max_x, max_y;
  ResidualBoxAnalyzer::FindBoundingBox(coeffs, num_coeffs, tdim, &max_x,
                                       &max_y);
  const uint32_t bw = TrfWidth[tdim];
  const uint32_t bh = TrfHeight[tdim];

  // Figure out the number of zeros we would save by storing the boundaries of
  // the residuals.
  uint32_t num_zeros_x = 0, num_zeros_y = 0;
  uint32_t last_zigzag_ind = 0;
  {
    // Find the last index in zigzag order.
    const uint16_t* const zigzag = ResidualIterator::GetZigzag(tdim);
    for (uint32_t i = kNumCoeffs[tdim] - 1; i > 0; --i) {
      if (coeffs[zigzag[i]] != 0) {
        last_zigzag_ind = i;
        break;
      }
    }
    // Count the saved zeros.
    ResidualIterator iter(tdim);
    uint32_t zigzag_ind = 0;
    if (first_is_dc) {
      ++iter;
      ++zigzag_ind;
    }
    for (; zigzag_ind <= last_zigzag_ind; ++iter, ++zigzag_ind) {
      if (coeffs[iter.Index()] == 0) {
        if (iter.x() > max_x) ++num_zeros_x;
        if (iter.y() > max_y) ++num_zeros_y;
      } else if (zigzag_ind == last_zigzag_ind) {
        break;
      }
    }
  }
  bool use_bounds_x, use_bounds_y;
  uint32_t ind_min = 0u;
  bool can_use_bounds_x, can_use_bounds_y;
  ResidualBoxAnalyzer::CanUseBounds(tdim, &can_use_bounds_x, &can_use_bounds_y);
  bool use_bounds;
  if ((can_use_bounds_x || can_use_bounds_y) && !is_pseudo_rate) {
    ResidualBoxAnalyzer::ShouldUseBounds(tdim, last_zigzag_ind, max_x, max_y,
                                         &use_bounds_x, &use_bounds_y);
    use_bounds = (use_bounds_x || use_bounds_y);
  } else {
    use_bounds = false;
  }

  if (can_use_bounds_x || can_use_bounds_y) {
    sm->Process(GetClusterMergedUV(channel, num_channels, method, tdim),
                kSymbolResidualUseBounds, use_bounds, "use_bounds", enc);
  }
  if (use_bounds) {
    uint8_t range_min_x, range_max_x, range_min_y, range_max_y;
    ResidualBoxAnalyzer::GetRangeX(tdim, &range_min_x, &range_max_x);
    ResidualBoxAnalyzer::GetRangeY(tdim, &range_min_y, &range_max_y);
    bool is_x_first;
    if (can_use_bounds_x && !can_use_bounds_y) {
      is_x_first = true;
    } else if (!can_use_bounds_x && can_use_bounds_y) {
      is_x_first = false;
    } else {
      if (use_bounds_x && !use_bounds_y) {
        is_x_first = true;
      } else if (!use_bounds_x && use_bounds_y) {
        is_x_first = false;
      } else {
        is_x_first = (range_max_x - range_min_x <= range_max_y - range_min_y);
      }

      sm->Process(GetClusterMergedUV(channel, num_channels, method, tdim),
                  kSymbolResidualBound1IsX, is_x_first, "is_x_first", enc);
    }
    if (is_x_first) {
      WriteBounds(channel, num_channels, method, tdim, /*is_x_first=*/true,
                  max_x, range_min_x, range_max_x, use_bounds_y, max_y, enc,
                  sm);
    } else {
      WriteBounds(channel, num_channels, method, tdim, /*is_x_first=*/false,
                  max_y, range_min_y, range_max_y, use_bounds_x, max_x, enc,
                  sm);
    }
    if (!use_bounds_x) max_x = bw - 1;
    if (!use_bounds_y) max_y = bh - 1;

    // Figure out the minimal index we will reach (and therefore the one from
    // which we need to store EOB).
    uint32_t min_zig_zag_ind_x, min_zig_zag_ind_y;
    ResidualBoxAnalyzer::FindBounds(tdim, max_x, max_y, &min_zig_zag_ind_x,
                                    &min_zig_zag_ind_y);
    if (use_bounds_x) ind_min = min_zig_zag_ind_x;
    if (use_bounds_y) ind_min = std::max(ind_min, min_zig_zag_ind_y);
  } else {
    use_bounds_x = use_bounds_y = false;
    max_x = bw - 1;
    max_y = bh - 1;
  }

  // Figure out the index of the first one from which all coefficients after are
  // just 1 or 0.
  uint32_t first_ending_one = num_coeffs;
  uint32_t last_non_zero_index = 0;
  BoundedResidualIterator iter(tdim, use_bounds_x, use_bounds_y, max_x, max_y);
  if (first_is_dc) ++iter;
  for (; !iter.IsDone(); ++iter) {
    const uint32_t ind = iter.Index();
    if (coeffs[ind] == 0) continue;
    last_non_zero_index = ind;
    if (std::abs(coeffs[ind]) == 1) {
      if (first_ending_one == num_coeffs) first_ending_one = ind;
    } else  {
      first_ending_one = num_coeffs;
    }
  }
  if (first_is_dc) assert(last_non_zero_index != 0);

  // Store the residuals.
  bool previous_is_a_one = false;
  bool has_written_zeros = false;
  bool has_only_ones_left = false;
  uint32_t sector_cluster;
  iter.Reset();
  if (first_is_dc) ++iter;  // Skip the DC.
  while (true) {
    const uint32_t x = iter.x();
    const uint32_t y = iter.y();
    const uint32_t i = iter.Index();
    // Set the probas to the ones of the current sector.
    const uint32_t sector = ResidualIO::GetSector(x, y, tdim);
    sector_cluster =
        GetClusterMergedUV(channel, num_channels, method, tdim, sector);

    const int16_t v = coeffs[i];
    // Count the number of consecutive zeros.
    if (v == 0) {
      // Write down the number of consecutive zeros by batches of
      // kResidualCons0Max.
      sm->Process(sector_cluster, kSymbolResidualIsZero, 1, "is_zero", enc);
      ++iter;
      do {
        uint32_t num_zeros_tmp = 0;
        const uint32_t max_num_coeffs_left = iter.MaxNumCoeffsLeft();
        while (num_zeros_tmp < kResidualCons0Max && coeffs[iter.Index()] == 0) {
          ++iter;
          ++num_zeros_tmp;
        }

        // Check if by adding the max number of 0s and one non-zero element, we
        // get further than allowed.
        if (max_num_coeffs_left < kResidualCons0Max + 1) {
          enc->PutRange(num_zeros_tmp, 0, max_num_coeffs_left - 1, "num_zeros");
        } else {
          // TODO(vrabaud) ANSAdaptiveSymbol performs 4% worse, investigate why.
          sm->Process(GetCluster(channel, num_channels, method),
                      kSymbolResNumZeros, num_zeros_tmp, "num_zeros", enc);
        }
        if (num_zeros_tmp != kResidualCons0Max) break;
      } while (true);

      has_written_zeros = true;
      continue;
    } else {
      // If we are the last element, we know we are not 0.
      if (!has_written_zeros && iter.MaxNumCoeffsLeft() > 1) {
        sm->Process(sector_cluster, kSymbolResidualIsZero, 0, "is_zero", enc);
      }
    }
    has_written_zeros = false;
    iter.SetAsNonZero();

    // If we don't know yet we only have 1s left.
    const uint16_t abs_v = std::abs(v);
    if (!has_only_ones_left &&
        !sm->Process(sector_cluster, kSymbolResidualIsOne, abs_v <= 1, "is_one",
                     enc)) {
      if (!sm->Process(sector_cluster, kSymbolResidualIsTwo, abs_v <= 2,
                       "is_two", enc)) {
        const uint16_t residual1 =
            std::min((uint32_t)(abs_v - 3), kResidual1Max);
        sm->Process(GetCluster(channel, num_channels, method), kSymbolBits0,
                    residual1, "residual1", enc);
        if (residual1 == kResidual1Max) {
          const uint16_t residual2 = abs_v - 3 - kResidual1Max;
          const Golomb golomb(residual2, /*prefix_size=*/1);
          sm->Process(GetCluster(channel, num_channels, method), kSymbolBits1,
                      golomb.prefix, "residual2_prefix", enc);
          enc->PutUValue(golomb.extra_bits_value, golomb.extra_bits_num,
                         "residual2_extra");
        }
      }
    }
    enc->PutBool(v < 0, "is_negative");
    // Exit if we are at the last element.
    const uint32_t zigzag_ind = iter.ZigZagIndex();
    ++iter;
    if (iter.IsDone()) break;
    if (zigzag_ind >= ind_min && iter.CanEOB() &&
        sm->Process(sector_cluster, kSymbolResidualEndOfBlock,
                    (i == last_non_zero_index), "eob", enc)) {
      break;
    }
    if (abs_v == 1 && !has_only_ones_left && !previous_is_a_one) {
      has_only_ones_left = (i == first_ending_one);
      sm->Process(sector_cluster, kSymbolResidualHasOnlyOnesLeft,
                  has_only_ones_left, "has_only_ones_left", enc);
    }
    previous_is_a_one = (abs_v == 1);
  }
}

WP2Status ResidualWriter::FindBestEncodingMethod(
    TrfSize tdim, const int16_t* const coeffs, uint32_t num_coeffs,
    bool first_is_dc, Channel channel, uint32_t num_channels,
    SymbolCounter* const counter, EncodingMethod* const encoding_method,
    float* const cost) {
  // We ignore the first value as it is DC coded.
  if (num_coeffs == 0 || (num_coeffs == 1 && first_is_dc)) {
    *encoding_method =
        (num_coeffs == 0) ? EncodingMethod::kAllZero : EncodingMethod::kDCOnly;
    if (cost != nullptr) *cost = 0.f;
    return WP2_STATUS_OK;
  }

  // Compute the best of the methods, and decide.
  float cost_best = std::numeric_limits<float>::max();
  *encoding_method = EncodingMethod::kMethod0;
  for (EncodingMethod method :
       {EncodingMethod::kMethod0, EncodingMethod::kMethod1}) {
    ANSEncCounter enc;
    counter->Clear();
    StoreCoeffs(coeffs, num_coeffs, first_is_dc, tdim, channel, num_channels,
                method, &enc, counter, /*is_pseudo_rate=*/false);
    const float cost_diff = enc.GetCost();
    if (cost_diff <= cost_best) {
      cost_best = cost_diff;
      *encoding_method = method;
    }
  }
  if (cost != nullptr) {
    *cost = cost_best;
  }
  return WP2_STATUS_OK;
}

WP2Status ResidualWriter::RecordCoeffs(const CodedBlock& cb, Channel channel,
                                       SymbolRecorder* const recorder,
                                       libgav1::AOMContext* const aom_context) {
  // Using a counter because it updates adaptive bits, which we do want.
  // TODO(maryla): be more consistant with using counter vs noop.
  ANSEncCounter enc;
  return WriteCoeffs(cb, channel, &enc, recorder, aom_context);
}

WP2Status ResidualWriter::WriteCoeffs(const CodedBlock& cb, Channel channel,
                                      ANSEncBase* enc, SymbolManager* const sm,
                                      libgav1::AOMContext* const aom_context) {
  const ANSDebugPrefix prefix(enc, (channel == kYChannel)   ? "Y"
                                   : (channel == kAChannel) ? "A"
                                                            : "UV");

  const CodedBlock::CodingParams& coding_params = cb.GetCodingParams(channel);
  const BlockSize split_size = GetSplitSize(cb.dim(), coding_params.split_tf);
  const uint32_t split_w = BlockWidthPix(split_size);
  const uint32_t split_h = BlockHeightPix(split_size);
  uint32_t tf_i = 0;
  for (uint32_t split_y = 0; split_y < cb.h_pix(); split_y += split_h) {
    for (uint32_t split_x = 0; split_x < cb.w_pix(); split_x += split_w) {
      const int16_t* const coeffs = cb.coeffs_[channel][tf_i];
      const uint32_t num_coeffs = cb.num_coeffs_[channel][tf_i];
      if (use_aom_coeffs_) {
        WP2_CHECK_STATUS(aom_writer_.WriteCoeffs(
            cb.x_pix() + split_x, cb.y_pix() + split_y, split_size,
            coding_params.tf, cb.IsFirstCoeffDC(channel), cb.is420_, coeffs,
            channel, /*do_update=*/true, sm, enc, aom_context));
      } else {
        const EncodingMethod method = cb.method_[channel][tf_i];
        if (method == EncodingMethod::kAllZero) {
          assert(cb.num_coeffs_[channel][tf_i] == 0);
        } else {
          const bool first_is_dc = cb.IsFirstCoeffDC(channel);
          if (first_is_dc) {
            const bool can_be_zero = (method != EncodingMethod::kDCOnly);
            StoreDC(channel, num_channels_, enc, sm, coeffs[0], can_be_zero);
          } else {
            assert(method != EncodingMethod::kDCOnly);
          }

          if (method != EncodingMethod::kDCOnly) {
            // Write the coefficients after DC.
            std::string label;
            WP2SPrint(&label, "C%d", GetMethodIndex(method));
            const ANSDebugPrefix coeff_prefix(enc, label.c_str());

            StoreCoeffs(coeffs, num_coeffs, first_is_dc, cb.tdim(channel),
                        channel, num_channels_, method, enc, sm,
                        /*is_pseudo_rate=*/false);
          }
        }
      }
      ++tf_i;
    }
  }
  return WP2_STATUS_OK;
}

WP2Status ResidualWriter::GetRate(Channel channel, uint32_t num_channels,
                                  TrfSize tdim, const int16_t* const coeffs,
                                  uint32_t num_coeffs, bool first_is_dc,
                                  SymbolCounter* const counter,
                                  float* const rate,
                                  EncodingMethod* const encoding_method) {
  EncodingMethod method;
  float residual_cost;
  counter->Clear();
  WP2_CHECK_STATUS(FindBestEncodingMethod(tdim, coeffs, num_coeffs, first_is_dc,
                                          channel, num_channels, counter,
                                          &method, &residual_cost));
  if (!first_is_dc) assert(method != EncodingMethod::kDCOnly);
  if (encoding_method != nullptr) *encoding_method = method;

  float dc = 0.f;
  if (first_is_dc && method != EncodingMethod::kAllZero) {
    ANSEncCounter enc;
    counter->Clear();
    StoreDC(channel, num_channels, &enc, counter, coeffs[0],
            /*can_be_zero=*/(method != EncodingMethod::kDCOnly));
    dc = enc.GetCost();
  }

  *rate = residual_cost + dc;
  return WP2_STATUS_OK;
}

WP2Status ResidualWriter::GetRateAOM(const CodedBlock& cb, Channel channel,
                                     const libgav1::AOMContext& aom_context,
                                     SymbolCounter* const counter,
                                     float* const rate) {
  counter->Clear();
  ANSEncCounter enc;

  const CodedBlock::CodingParams& coding_params = cb.GetCodingParams(channel);
  const BlockSize split_size = GetSplitSize(cb.dim(), coding_params.split_tf);
  const uint32_t split_w = BlockWidthPix(split_size);
  const uint32_t split_h = BlockHeightPix(split_size);
  uint32_t tf_i = 0;
  for (uint32_t split_y = 0; split_y < cb.h_pix(); split_y += split_h) {
    for (uint32_t split_x = 0; split_x < cb.w_pix(); split_x += split_w) {
      WP2_CHECK_STATUS(libgav1::AOMResidualWriter::WriteCoeffs(
          cb.x_pix() + split_x, cb.y_pix() + split_y, split_size,
          coding_params.tf, cb.IsFirstCoeffDC(channel), cb.is420_,
          cb.coeffs_[channel][tf_i], channel, /*do_update=*/false, counter,
          &enc, const_cast<libgav1::AOMContext*>(&aom_context)));
      ++tf_i;
    }
  }

  *rate = enc.GetCost();
  return WP2_STATUS_OK;
}

WP2Status ResidualWriter::GetPseudoRate(Channel channel, uint32_t num_channels,
                                        TrfSize tdim,
                                        const int16_t* const coeffs,
                                        uint32_t num_coeffs, bool first_is_dc,
                                        SymbolCounter* const counter,
                                        float* const pseudo_rate) {
  *pseudo_rate = 0;

  // Return right away if there are no elements (kAllZero).
  if (num_coeffs == 0) return WP2_STATUS_OK;

  if (first_is_dc) {
    // Add cost of storing DC.
    // 11 is the number of yuv bits of YCoCg. It differs for other transforms
    // but it's in the same ballpark.
    static constexpr uint32_t kChannelBits = 11;
    const uint32_t extra_bits_num =
        Golomb((coeffs[0] >= 0) ? 2 * coeffs[0] : -2 * coeffs[0] - 1,
               /*prefix_size=*/1)
            .extra_bits_num;
    // We assume a Golomb storage for DC and store DC as in StoreDC.
    *pseudo_rate +=
        WP2Log2(kChannelBits) /* prefix */ + extra_bits_num /* abs_val */;

    if (num_coeffs == 1) return WP2_STATUS_OK;
  }

  counter->Clear();
  ANSEncCounter enc;
  StoreCoeffs(coeffs, num_coeffs, first_is_dc, tdim, channel, num_channels,
              EncodingMethod::kMethod0, &enc, counter, /*is_pseudo_rate=*/true);
  const float residual_cost = enc.GetCost();

  // Also consider the cost for a proba of ANSBinSymbol(8, 2).
  const uint32_t start_i = first_is_dc ? 1 : 0;
  uint32_t num_zeros_streaks = 0;
  uint32_t num_non_zeros = 0;
  for (uint32_t i = start_i; i < num_coeffs; ++i) {
    if (coeffs[i] != 0) {
      if (i > start_i && coeffs[i - 1] == 0) ++num_zeros_streaks;
      ++num_non_zeros;
    }
  }
  constexpr float p = 1.f / 5.f;
  const float is_not_zeros_costs_2_8 =
      num_zeros_streaks * std::log2(1.f - p) + num_non_zeros * std::log2(p);
  const float is_not_zeros_costs_8_2 =
      num_zeros_streaks * std::log2(p) + num_non_zeros * std::log2(1.f - p);

  *pseudo_rate +=
      std::min(residual_cost,
               residual_cost - is_not_zeros_costs_2_8 + is_not_zeros_costs_8_2);
  return WP2_STATUS_OK;
}

}  // namespace WP2
