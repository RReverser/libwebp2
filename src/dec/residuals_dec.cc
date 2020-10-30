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
// WP2 lossy decoding of residuals.
//
// Author: Skal (pascal.massimino@gmail.com)

#include "src/common/constants.h"
#include "src/dec/wp2_dec_i.h"
#include "src/utils/ans_utils.h"
#include "src/utils/utils.h"
#include "src/wp2/decode.h"

namespace WP2 {

//------------------------------------------------------------------------------

WP2Status ResidualReader::ReadHeaderForResidualSymbols(uint32_t num_coeffs_max,
                                                       Channel channel,
                                                       SymbolReader* const sr) {
  bool is_maybe_used[kResidual1Max + 1];
  for (EncodingMethod method :
       {EncodingMethod::kMethod0, EncodingMethod::kMethod1}) {
    std::string label;
    const uint32_t cluster = GetCluster(channel, num_channels_, method);

    // Read the dictionaries for the number of consecutive zeros.
    const uint32_t method_index = GetMethodIndex(method);
    (void)method_index;
    WP2SPrint(&label, "%s_num_zeros_%d", kChannelStr[channel], method_index);
    WP2_CHECK_STATUS(sr->ReadHeader(cluster, num_coeffs_max, kSymbolResNumZeros,
                                    label.c_str()));

    // Read the dictionaries for small residuals.
    WP2SPrint(&label, "%s_bits0_%d", kChannelStr[channel], method_index);
    WP2_CHECK_STATUS(
        sr->ReadHeader(cluster, num_coeffs_max, kSymbolBits0, label.c_str()));
    sr->GetPotentialUsage(cluster, kSymbolBits0, is_maybe_used,
                          kResidual1Max + 1);

    if (is_maybe_used[kResidual1Max]) {
      // Read the dictionaries for prefixes of big residuals.
      WP2SPrint(&label, "%s_bits1_%d", kChannelStr[channel], method_index);
      WP2_CHECK_STATUS(
          sr->ReadHeader(cluster, num_coeffs_max, kSymbolBits1, label.c_str()));
    }
  }
  return WP2_STATUS_OK;
}

WP2Status ResidualReader::ReadHeader(SymbolReader* const sr,
                                     uint32_t num_coeffs_max_y,
                                     uint32_t num_coeffs_max_uv,
                                     uint32_t num_transforms, bool has_alpha,
                                     bool has_lossy_alpha) {
  num_channels_ = (has_alpha ? 4 : 3);
  if (use_aom_coeffs_) {
    for (Symbol sym : kSymbolsForAOMCoeffs) {
      WP2_CHECK_STATUS(sr->ReadHeader(num_transforms, sym, "aom_symbols"));
    }
  } else {
    for (Channel channel : {kYChannel, kUChannel, kVChannel, kAChannel}) {
      if (channel == kAChannel && !has_lossy_alpha) continue;

      WP2_CHECK_STATUS(ReadHeaderForResidualSymbols(
          (channel == kYChannel || channel == kAChannel) ? num_coeffs_max_y
                                                         : num_coeffs_max_uv,
          channel, sr));
    }

    for (Symbol sym : kSymbolsForCoeffs) {
      WP2_CHECK_STATUS(
          sr->ReadHeader(num_coeffs_max_y, sym, "residual_symbols"));
    }

    for (Symbol sym : kSymbolsForCoeffsPerTf) {
      WP2_CHECK_STATUS(sr->ReadHeader(num_transforms, sym, "block_symbols"));
    }
  }

  if (has_lossy_alpha) {
    WP2_CHECK_STATUS(sr->ReadHeader(num_transforms, kSymbolEncodingMethodA,
                                    "coeff_method_alpha"));
  }
  WP2_CHECK_STATUS(
      sr->ReadHeader(num_transforms, kSymbolHasCoeffs, "has_coeffs"));
  WP2_CHECK_STATUS(
      sr->ReadHeader(num_transforms, kSymbolEncodingMethod, "coeff_method"));
  return WP2_STATUS_OK;
}

//------------------------------------------------------------------------------

static int16_t ReadDC(Channel channel, uint32_t num_channels, bool can_be_zero,
                      SymbolReader* const sr) {
  int n =
      sr->Read(ResidualIO::GetCluster(channel, num_channels), kSymbolDC, "DC");
  if (!can_be_zero) ++n;
  // Converting the unsigned value to a signed one.
  return (n % 2 == 0) ? n / 2 : -((n + 1) / 2);
}

inline void ReadBounds(Channel channel, uint32_t num_channels,
                       EncodingMethod method, TrfSize tdim, bool is_x_first,
                       ANSDec* const dec, SymbolReader* const sr,
                       bool* const use_bounds2, uint32_t* const val1,
                       uint32_t* const val2) {
  uint8_t min1, max1;
  if (is_x_first) {
    ResidualBoxAnalyzer::GetRangeX(tdim, &min1, &max1);
  } else {
    ResidualBoxAnalyzer::GetRangeY(tdim, &min1, &max1);
  }
  *val1 = dec->ReadRange(min1, max1, "bound1");
  uint8_t min2, max2;
  if (is_x_first) {
    ResidualBoxAnalyzer::GetRangePerX(tdim, *val1, &min2, &max2);
  } else {
    ResidualBoxAnalyzer::GetRangePerY(tdim, *val1, &min2, &max2);
  }
  if (min2 != 255 && min2 <= max2) {
    *use_bounds2 = sr->Read(
        ResidualReader::GetClusterMergedUV(channel, num_channels, method, tdim),
        kSymbolResidualUseBound2, "use_bound2");
  } else {
    *use_bounds2 = false;
  }
  if (*use_bounds2) *val2 = dec->ReadRange(min2, max2, "bound2");
}

WP2Status ResidualReader::ReadCoeffs(Channel channel, ANSDec* const dec,
                                     SymbolReader* const sr,
                                     CodedBlock* const cb,
                                     libgav1::AOMContext* const aom_context,
                                     BlockInfo* const info) {
  const bool is_uv = (channel == kUChannel || channel == kVChannel);
  ANSDebugPrefix prefix(dec, is_uv ? "UV" : (channel == kAChannel) ? "A" : "Y");

  CodedBlock::CodingParams* const params = cb->GetCodingParams(channel);
  const BlockSize split_size = GetSplitSize(cb->dim(), params->split_tf);
  const uint32_t split_w = BlockWidthPix(split_size);
  const uint32_t split_h = BlockHeightPix(split_size);
  uint32_t tf_i = 0;
  for (uint32_t split_y = 0; split_y < cb->h_pix(); split_y += split_h) {
    for (uint32_t split_x = 0; split_x < cb->w_pix(); split_x += split_w) {
      int16_t* const coeffs = cb->coeffs_[channel][tf_i];
      std::fill(coeffs, coeffs + cb->NumCoeffsPerTransform(channel), 0);

      if (use_aom_coeffs_) {
        const int max_num_coeffs = (cb->num_coeffs_[channel][tf_i] > 0)
                                       ? WP2::kNumCoeffs[cb->tdim(channel)]
                                       : 0;
        WP2_CHECK_STATUS(aom_reader_.ReadCoeffs(
            channel, cb->x_pix() + split_x, cb->y_pix() + split_y, split_size,
            cb->is420_, params->tf, cb->IsFirstCoeffDC(channel), max_num_coeffs,
            dec, sr, aom_context, coeffs, &cb->num_coeffs_[channel][tf_i]));
        SetGeometry(cb->num_coeffs_[channel][tf_i],
                    &cb->method_[channel][tf_i]);
      } else if (cb->method_[channel][tf_i] == EncodingMethod::kAllZero) {
        cb->num_coeffs_[channel][tf_i] = 0;
      } else {
        WP2_CHECK_STATUS(ReadCoeffsMethod01(channel, tf_i, dec, sr, cb, info));
      }
      ++tf_i;
    }
  }
  return WP2_STATUS_OK;
}

WP2Status ResidualReader::ReadCoeffsMethod01(Channel channel, uint32_t tf_i,
                                             ANSDec* const dec,
                                             SymbolReader* const sr,
                                             CodedBlock* const cb,
                                             BlockInfo* const info) {
  int16_t* const coeffs = cb->coeffs_[channel][tf_i];
  const EncodingMethod method = cb->method_[channel][tf_i];
  const TrfSize tdim = cb->tdim(channel);
  const bool first_is_dc = cb->IsFirstCoeffDC(channel);

  if (first_is_dc) {
    // Read the DC.
    const bool can_be_zero = (method != EncodingMethod::kDCOnly);
    coeffs[0] = ReadDC(channel, num_channels_, can_be_zero, sr);
    // [-2047:2047] if zero is included, -2048 is possible if it is not.
    // TODO(yguyon): Or is it GlobalParams::GetMaxAbsDC()?
    assert(coeffs[0] >= -kMaxDcValue - (can_be_zero ? 0 : 1) &&
           coeffs[0] <= kMaxDcValue);

    if (method == EncodingMethod::kDCOnly) {
      cb->num_coeffs_[channel][tf_i] = 1;
      return WP2_STATUS_OK;
    }
  } else {
    assert(method != EncodingMethod::kDCOnly);
  }

  std::string label;
  WP2SPrint(&label, "C%d", GetMethodIndex(method));
  ANSDebugPrefix coeff_prefix(dec, label.c_str());

  const uint32_t bw = TrfWidth[tdim];
  const uint32_t bh = TrfHeight[tdim];
  bool use_bounds_x, use_bounds_y;
  uint32_t max_x, max_y;
  uint32_t ind_min = 0u;
  bool can_use_bounds_x, can_use_bounds_y;
  ResidualBoxAnalyzer::CanUseBounds(tdim, &can_use_bounds_x, &can_use_bounds_y);
  if ((can_use_bounds_x || can_use_bounds_y) &&
      sr->Read(GetClusterMergedUV(channel, num_channels_, method, tdim),
               kSymbolResidualUseBounds, "use_bounds")) {
    bool is_x_first;
    if (can_use_bounds_x && !can_use_bounds_y) {
      is_x_first = true;
    } else if (!can_use_bounds_x && can_use_bounds_y) {
      is_x_first = false;
    } else {
      is_x_first =
          sr->Read(GetClusterMergedUV(channel, num_channels_, method, tdim),
                   kSymbolResidualBound1IsX, "is_x_first");
    }
    if (is_x_first) {
      use_bounds_x = true;
      ReadBounds(channel, num_channels_, method, tdim, /*is_x_first=*/true, dec,
                 sr, &use_bounds_y, &max_x, &max_y);
      if (!use_bounds_y) max_y = bh - 1;
    } else {
      use_bounds_y = true;
      ReadBounds(channel, num_channels_, method, tdim, /*is_x_first=*/false,
                 dec, sr, &use_bounds_x, &max_y, &max_x);
      if (!use_bounds_x) max_x = bw - 1;
    }

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

  // Debug.
  if (info != nullptr) {
    if (use_bounds_x) {
      info->residual_info[channel][tf_i].push_back("use x bound " +
                                                   std::to_string(max_x));
    }
    if (use_bounds_y) {
      info->residual_info[channel][tf_i].push_back("use y bound " +
                                                   std::to_string(max_y));
    }
  }

  bool has_written_zeros = false;
  bool has_only_ones_left = false;
  bool previous_is_a_one = false;
  uint32_t sector_cluster;
  BoundedResidualIterator iter(tdim, use_bounds_x, use_bounds_y, max_x, max_y);
  if (first_is_dc) ++iter;  // Skip the DC.
  for (; !iter.IsDone();) {
    const uint32_t x = iter.x();
    const uint32_t y = iter.y();
    const uint32_t i = iter.Index();
    const uint32_t sector = ResidualIO::GetSector(x, y, tdim);
    sector_cluster =
        GetClusterMergedUV(channel, num_channels_, method, tdim, sector);

    // If we have more than one element left and not written 0s before, check if
    // there is any 0 coming.
    if (!has_written_zeros && iter.MaxNumCoeffsLeft() > 1 &&
        sr->Read(sector_cluster, kSymbolResidualIsZero, "is_zero")) {
      // Read the number of consecutive 0s we have, by batches of
      // kResidualCons0Max.
      ++iter;
      int32_t num_zeros_tmp;
      do {
        if (iter.MaxNumCoeffsLeft() < kResidualCons0Max + 1u) {
          // If the number of elements left is smaller than the max number of
          // possible 0s plus one non-zero element, go with a ramge to force
          // feasibility.
          num_zeros_tmp = dec->ReadRValue(iter.MaxNumCoeffsLeft(), "num_zeros");
        } else {
          num_zeros_tmp = sr->Read(GetCluster(channel, num_channels_, method),
                                   kSymbolResNumZeros, "num_zeros");
        }
        for (int32_t j = 0; j < num_zeros_tmp; ++j) ++iter;
      } while (num_zeros_tmp == kResidualCons0Max);
      has_written_zeros = true;
      continue;
    }

    has_written_zeros = false;
    iter.SetAsNonZero();

    uint32_t abs_v;
    if (has_only_ones_left ||
        sr->Read(sector_cluster, kSymbolResidualIsOne, "is_one")) {
      abs_v = 1;
    } else {
      if (sr->Read(sector_cluster, kSymbolResidualIsTwo, "is_two")) {
        abs_v = 2;
      } else {
        const uint32_t residual1 =
            sr->Read(GetCluster(channel, num_channels_, method), kSymbolBits0,
                     "residual1");
        if (residual1 == kResidual1Max) {
          abs_v = 3 + kResidual1Max;
          const uint32_t prefix =
              sr->Read(GetCluster(channel, num_channels_, method), kSymbolBits1,
                       "residual2_prefix");
          const uint32_t extra =
              dec->ReadUValue(Golomb::NumExtraBits(prefix, /*prefix_size=*/1),
                              "residual2_extra");
          abs_v += Golomb::Merge(prefix, /*prefix_size=*/1, extra);
        } else {
          abs_v = 3 + residual1;
        }
      }
    }
    coeffs[i] =
        (int16_t)(dec->ReadBool("is_negative") ? -(int16_t)abs_v : abs_v);

    const uint32_t zigzag_ind = iter.ZigZagIndex();
    ++iter;
    // Exit if we are at the last element.
    if (iter.IsDone()) break;
    // Read an End Of Block if we have reached both sides of the box.
    if (zigzag_ind >= ind_min && iter.CanEOB() &&
        sr->Read(sector_cluster, kSymbolResidualEndOfBlock, "eob")) {
      break;
    }
    if (abs_v == 1 && !has_only_ones_left && !previous_is_a_one) {
      has_only_ones_left = (bool)sr->Read(
          sector_cluster, kSymbolResidualHasOnlyOnesLeft, "has_only_ones_left");
    }
    previous_is_a_one = (abs_v == 1);
  }
  cb->num_coeffs_[channel][tf_i] = 0;
  const uint32_t start_i = (first_is_dc ? 1 : 0);
  for (uint32_t i = cb->NumCoeffsPerTransform(channel); i-- > start_i;) {
    if (coeffs[i] != 0) {
      cb->num_coeffs_[channel][tf_i] = i + 1;
      break;
    }
  }
  assert(cb->num_coeffs_[channel][tf_i] > start_i);
  for (uint32_t i = start_i; i < cb->num_coeffs_[channel][tf_i]; ++i) {
    assert(coeffs[i] >= -kMaxCoeffValue);
    assert(coeffs[i] <= kMaxCoeffValue);
  }

  // Debug.
  if (info != nullptr) {
    // Compute the real box bounds.
    uint32_t real_max_x = 0u, real_max_y = 0u;
    const uint32_t max_i = cb->NumCoeffsPerTransform(channel) - 1;
    for (uint32_t i = 0u; i <= max_i; ++i) {
      const uint32_t x = i % bw;
      const uint32_t y = i / bw;
      if (coeffs[i] != 0) {
        real_max_x = std::max(real_max_x, x);
        real_max_y = std::max(real_max_y, y);
      }
    }
    // last_i contains the index of the last element in the whole block.
    // last_j contains the index of the last element in the chosen box (there
    // can be none).
    // last_k contains the index of the last element in the bounding box.
    uint32_t last_i = 0u, last_j = 0u, last_k = 0u;
    for (uint32_t i = 0u, j = 0u, k = 0u; i <= max_i; ++i) {
      const uint32_t x = i % bw;
      const uint32_t y = i / bw;
      // No need to store anything if we are out of bounds: it is 0s.
      if (x > max_x || y > max_y) continue;
      if (coeffs[i] != 0) {
        last_i = i;
        last_j = j;
        last_k = k;
      }
      ++j;
      if (x <= real_max_x && y <= real_max_y) ++k;
    }
    assert(last_j <= last_i);
    assert(last_k <= last_j);
    info->residual_info[channel][tf_i].push_back("last index: " +
                                                 std::to_string(last_i));
    if (use_bounds_x && use_bounds_y) {
      info->residual_info[channel][tf_i].push_back(
          "last index in chosen full box: " + std::to_string(last_j));
    } else if (use_bounds_x || use_bounds_y) {
      info->residual_info[channel][tf_i].push_back(
          "last index in chosen box: " + std::to_string(last_j));
      info->residual_info[channel][tf_i].push_back("last index in full box: " +
                                                   std::to_string(last_k));
    } else {
      info->residual_info[channel][tf_i].push_back("last index in full box: " +
                                                   std::to_string(last_k));
    }
  }

  return WP2_STATUS_OK;
}

//------------------------------------------------------------------------------

}  // namespace WP2
