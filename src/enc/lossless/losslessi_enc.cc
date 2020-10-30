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
// main entry for the lossless encoder.
//
// Author: Vincent Rabaud (vrabaud@google.com)

#include "src/enc/lossless/losslessi_enc.h"

#include <algorithm>
#include <cassert>
#include <cstdlib>
#include <memory>
#include <string>

#include "src/common/color_precision.h"
#include "src/common/symbols.h"
#include "src/dsp/lossless/lossless_common.h"
#include "src/enc/lossless/backward_references_enc.h"
#include "src/enc/lossless/histogram_enc.h"
#include "src/enc/symbols_enc.h"
#include "src/utils/ans.h"
#include "src/utils/ans_utils.h"
#include "src/utils/plane.h"
#include "src/utils/vector.h"

namespace WP2L {

// -----------------------------------------------------------------------------

static void swap(EncodeInfo& a, EncodeInfo& b) {
  swap(a.line_tokens, b.line_tokens);
  swap(a.bits_per_pixel, b.bits_per_pixel);
}

// Make sure the array contains indices in the green channel.
// From one to the next, those should be the same or less, and if they increase,
// they increase only by one.
static bool IsArrayMadeOfLabels(const int16_t v[], uint32_t size) {
  if (size == 0) return true;
  assert(size % 4 == 0);
  int16_t max_value = v[2];
  if (max_value != 0) return false;
  for (uint32_t i = 6; i < size; i += 4) {
    if (v[i] <= max_value) continue;
    if (v[i] - max_value > 1) return false;
    max_value = v[i];
  }
  return true;
}

WP2Status StorePixels(uint32_t width, uint32_t height, uint32_t histo_bits,
                      const BackwardRefs& refs,
                      const uint16_t* const histogram_symbols,
                      WP2::ANSEncBase* const enc, WP2::SymbolManager* const sw,
                      EncodeInfo* const encode_info) {
  const uint32_t histo_xsize =
      histo_bits ? WP2LSubSampleSize(width, histo_bits) : 1;
  const uint32_t tile_mask = (histo_bits == 0) ? 0 : ~((1u << histo_bits)-1);
  // x and y trace the position in the image.
  uint32_t x = 0;
  uint32_t y = 0;
  uint32_t tile_x = x & tile_mask;
  uint32_t tile_y = y & tile_mask;
  uint32_t histogram_ix =
      (histogram_symbols == nullptr) ? 0 : histogram_symbols[0];

  if (encode_info != nullptr) {
    assert(encode_info->line_tokens.empty() ||
           encode_info->line_tokens.size() == height);
    assert(encode_info->bits_per_pixel.empty() ||
           encode_info->bits_per_pixel.size() == width * height);
  }
  const bool get_cost =
      encode_info != nullptr && !encode_info->bits_per_pixel.empty();

  // Define unused symbols.
  const bool has_alpha = LosslessSymbolsInfo::HasAlpha(sw->symbols_info());
  const bool is_label_image =
      LosslessSymbolsInfo::IsLabelImage(sw->symbols_info());

  double cost;
  double* const cost_ptr = get_cost ? &cost : nullptr;
  for (RefsCursor c(refs); c.Ok(); c.Next()) {
    cost = 0.;
    if (histo_bits != 0u) {
      if ((tile_x != (x & tile_mask)) || (tile_y != (y & tile_mask))) {
        tile_x = x & tile_mask;
        tile_y = y & tile_mask;
        const uint32_t ind =
            (y >> histo_bits) * histo_xsize + (x >> histo_bits);
        histogram_ix =
            (histogram_symbols == nullptr) ? ind : histogram_symbols[ind];
      }
    }
    const PixelMode& v = *c.cur_pos_;
    sw->ProcessWithCost(histogram_ix, kSymbolType, v.GetMode(), "lossless_mode",
                        enc, cost_ptr);
    switch (v.GetMode()) {
      case kSymbolTypeLiteral: {
        WP2::ANSDebugPrefix prefix(enc, "literal");
        if (has_alpha) {
          sw->ProcessWithCost(histogram_ix, kSymbolA, v.GetLiteral(0), "A", enc,
                              cost_ptr);
        }
        if (!is_label_image) {
          sw->ProcessWithCost(histogram_ix, kSymbolR, v.GetLiteral(1), "R", enc,
                              cost_ptr);
          sw->ProcessWithCost(histogram_ix, kSymbolB, v.GetLiteral(3), "B", enc,
                              cost_ptr);
        }
        sw->ProcessWithCost(histogram_ix, kSymbolG, v.GetLiteral(2), "G", enc,
                            cost_ptr);
        break;
      }
      case kSymbolTypeCacheIdx: {
        const uint32_t code = v.GetCacheIdx();
        sw->ProcessWithCost(histogram_ix, kSymbolCache, code,
                            v.GetCacheIdxRange() - 1, "cache_idx", enc,
                            cost_ptr);
        break;
      }
      case kSymbolTypeCopy: {
        WP2::ANSDebugPrefix prefix(enc, "copy");
        int bits, n_bits;
        int code;
        WP2LPrefixEncode(v.GetLength(), &code, &n_bits, &bits);
        sw->ProcessWithCost(histogram_ix, kSymbolLen, code, "len", enc,
                            cost_ptr);
        if (n_bits) {
          enc->PutUValue(bits, n_bits, "len_bits");
          cost += n_bits;
        }

        const uint32_t distance = DistanceToPlaneCode(width, v.GetDistance());
        WP2LPrefixEncode(distance, &code, &n_bits, &bits);
        sw->ProcessWithCost(histogram_ix, kSymbolDist, code, "dist", enc,
                            cost_ptr);
        if (n_bits) {
          if (n_bits <= IO_BITS) {
            enc->PutUValue(bits, n_bits, "dist_bits");
          } else {
            enc->PutUValue(bits & IO_MASK, IO_BITS, "dist_bits1");
            enc->PutUValue(bits >> IO_BITS, n_bits - IO_BITS, "dist_bits2");
          }
          cost += n_bits;
        }

        break;
      }
      default:
        assert(false);
    }
    if (get_cost) {
      const float bits_per_pixel = cost / v.GetLength();
      for (uint32_t i = 0; i < v.GetLength(); ++i) {
        encode_info->bits_per_pixel[y * width + x + i] = bits_per_pixel;
      }
    }
    if (histo_bits != 0u) {
      x += v.GetLength();
      while (x >= width) {
        x -= width;
        if (encode_info != nullptr && !encode_info->line_tokens.empty()) {
          encode_info->line_tokens[y] = enc->NumTokens();
        }
        ++y;
      }
    }
  }
  return enc->GetStatus();
}

WP2Status StorePixels(uint32_t width, const BackwardRefs& refs,
                      WP2::ANSEncBase* const enc,
                      WP2::SymbolManager* const sw) {
  WP2_CHECK_STATUS(StorePixels(width, /*height=*/0, /*histo_bits=*/0, refs,
                               /*histogram_symbols=*/nullptr, enc, sw,
                               /*encode_info=*/nullptr));
  return WP2_STATUS_OK;
}

// Store the ANS statistics in the bitstream, so that they can be re-used
// when decoding the symbols.
// n_pixels is the number of pixels in the image, and use_palette a flag
// indicating whether the image uses a palette and therefore has no information
// in the A,R,B channels.
WP2Status WriteHeaders(const WP2::SymbolRecorder& recorder,
                       const LosslessSymbolsInfo& symbols_info,
                       uint32_t num_pixels, WP2::ANSEncBase* const enc,
                       WP2::ANSDictionaries* const dicts,
                       WP2::SymbolWriter* const sw) {
  WP2::ANSDebugPrefix prefix(enc, "symbols");
  WP2_CHECK_STATUS(sw->Init(symbols_info));
  WP2_CHECK_STATUS(sw->Allocate());

  // Iterate over all histograms and get the aggregate number of codes used.
  for (uint32_t i = 0; i < symbols_info.NumClusters(kSymbolType); ++i) {
    // Deal with the first symbol type.
    WP2_CHECK_STATUS(sw->WriteHeader(i, num_pixels, kSymbolType, recorder,
                                     kSymbolNames[kSymbolType], enc, dicts));
    bool is_maybe_used[kSymbolTypeNum];
    // TODO(vrabaud) For kSymbolType, it's either a dictionary, or a trivial
    // symbol (not a range) so we could limit the storage to those two methods.
    // We would get exact counts when decoding and we would therefore know
    // exactly what is used or not. That would lead to a very small gain in
    // compressibility.
    sw->GetPotentialUsage(i, kSymbolType, is_maybe_used, kSymbolTypeNum);

    // Figure out the number of elements and the max values in bits.
    for (int s = 1; s < kSymbolNum; ++s) {
      const Symbol sym = static_cast<Symbol>(s);
      // Ignore symbols we don't even use.
      if (symbols_info.IsSymbolUseless(sym, is_maybe_used)) continue;
      WP2_CHECK_STATUS(sw->WriteHeader(i, num_pixels, sym, recorder,
                                       kSymbolNames[s], enc, dicts));
    }
  }

  return WP2_STATUS_OK;
}

WP2Status EncodeImageNoClusters(
    WP2::ANSEncBase* const enc, WP2::ANSDictionaries* const dicts,
    const int16_t* const argb, HashChain* const hash_chain,
    BackwardRefsPool* const ref_pool, uint32_t width, uint32_t height,
    const LosslessSymbolsInfo& symbols_info, int speed) {
  WP2::SymbolWriter sw;

  // Calculate backward references from ARGB image.
  WP2_CHECK_STATUS(hash_chain->Fill(speed, argb, width, height));
  LosslessSymbolsInfo symbols_info_new;
  WP2_CHECK_STATUS(symbols_info_new.CopyFrom(symbols_info));
  BackwardRefsPool::RefsPtr refs = BackwardRefsPool::GetEmptyBackwardRefs();
  CacheConfig cache_config;
  WP2_CHECK_STATUS(GetBackwardReferences(
      width, height, argb, speed, kLZ77Standard | kLZ77RLE,
      /*cache_bits_max=*/0, *hash_chain, symbols_info_new, &cache_config,
      ref_pool, &refs));

  // Record symbol statistics.
  WP2::SymbolRecorder recorder;
  WP2::ANSEncNoop noop;
  WP2_CHECK_STATUS(recorder.Allocate(symbols_info_new, /*num_records=*/0));
  WP2_CHECK_STATUS(StorePixels(width, *refs, &noop, &recorder));

  // Store headers.
  WP2_CHECK_STATUS(WriteHeaders(recorder, symbols_info_new, width * height, enc,
                                dicts, &sw));

  // Store actual literals.
  return StorePixels(width, height, /*histo_bits=*/0, *refs,
                     /*histogram_symbols=*/nullptr, enc, &sw);
}

static WP2Status EncodeImageInternal(
    const uint16_t* const argb_in, const CrunchConfig& config, uint32_t width,
    uint32_t height, uint32_t histogram_bits, uint32_t cache_bits_max,
    uint32_t speed, WP2::ANSEnc* const enc, WP2::ANSDictionaries* const dicts,
    HashChain* const hash_chain, BackwardRefsPool* ref_pool,
    LosslessSymbolsInfo* const symbols_info_best,
    EncodeInfo* const encode_info) {
  const uint32_t histogram_width = WP2LSubSampleSize(width, histogram_bits);
  const uint32_t histogram_height = WP2LSubSampleSize(height, histogram_bits);
  const uint32_t histogram_image_xysize = histogram_width * histogram_height;
  HistogramSet histogram_image;
  WP2::Vector_u16 histogram_symbols;
  WP2::ANSEnc enc_best;
  WP2::ANSDictionaries dicts_best;
  WP2::ANSEnc enc_init;
  WP2::ANSDictionaries dicts_init;
  float best_cost = std::numeric_limits<float>::max();
  WP2::SymbolWriter sw;

  assert(width * height >= 1 && histogram_image_xysize >= 1);
  assert(histogram_bits >= kHistogramBitsMin);
  assert(histogram_bits <= kHistogramBitsMax);

  // If the value is different from zero, it has been set during the
  // palette analysis.
  LosslessSymbolsInfo symbols_info_init;
  WP2_CHECK_STATUS(symbols_info_init.CopyFrom(*symbols_info_best));

  // Convert input to int16_t.
  WP2::Vector_s16 argb_vec;
  WP2_CHECK_ALLOC_OK(argb_vec.resize(4 * width * height));
  int16_t* const argb = argb_vec.data();
  std::copy(argb_in, argb_in + 4 * width * height, argb);

  // Calculate backward references from ARGB image.
  WP2_CHECK_STATUS(enc_init.Clone(*enc));
  WP2_CHECK_STATUS(dicts_init.CopyFrom(*dicts));
  WP2_CHECK_STATUS(enc_best.Clone(*enc));
  WP2_CHECK_STATUS(dicts_best.CopyFrom(*dicts));
  WP2_CHECK_ALLOC_OK(histogram_symbols.resize(histogram_image_xysize));
  WP2_CHECK_STATUS(hash_chain->Fill(speed, argb, width, height));
  for (uint32_t lz77s_idx = 0; lz77s_idx < config.lz77s_types_to_try_size;
       ++lz77s_idx) {
    LosslessSymbolsInfo symbols_info;
    WP2_CHECK_STATUS(symbols_info.CopyFrom(symbols_info_init));
    BackwardRefsPool::RefsPtr refs_best =
        BackwardRefsPool::GetEmptyBackwardRefs();
    CacheConfig cache_config;
    WP2_CHECK_STATUS(GetBackwardReferences(
        width, height, argb, speed, config.lz77s_types_to_try[lz77s_idx],
        cache_bits_max, *hash_chain, symbols_info, &cache_config, ref_pool,
        &refs_best));
    symbols_info.SetCacheRange(GetColorCacheRange(cache_config));
    WP2_CHECK_ALLOC_OK(
        histogram_image.Allocate(histogram_image_xysize, symbols_info));

    // Build histogram image and symbols from backward references.
    WP2_CHECK_STATUS(GetHistoImageSymbols(
        width, height, *refs_best, speed, histogram_bits, symbols_info,
        &histogram_image, histogram_symbols.data()));

    {
      WP2::ANSDebugPrefix prefix_gh(enc, "GlobalHeader");
      // Color Cache parameters.
      if (enc->PutBool(cache_config.type != CacheType::kNone, "color_cache")) {
        assert(cache_config.cache_bits > 0);
        enc->PutRange(cache_config.cache_bits, 1, kMaxCacheBits,
                      "color_cache_bits");
      }

      // Find the maximum subsampling that would result in >= 2 histograms.
      uint32_t max_sampling;
      bool more_than_one_hist = false;
      for (max_sampling = kHistogramBitsMax + 1;
           !more_than_one_hist && max_sampling-- > kHistogramBitsMin;) {
        more_than_one_hist = (WP2LSubSampleSize(width, max_sampling) > 1 ||
                              WP2LSubSampleSize(height, max_sampling) > 1);
      }

      // Histogram image.
      const uint32_t num_clusters = histogram_image.histograms_.size();
      if (num_clusters > 1) assert(more_than_one_hist);
      const bool write_histogram_image =
          (more_than_one_hist &&
           enc->PutBool(num_clusters > 1, "write_histogram_image"));
      if (write_histogram_image) {
        WP2::Vector_s16 labels;
        WP2_CHECK_ALLOC_OK(labels.resize(4 * histogram_image_xysize));
        for (size_t i = 0; i < histogram_image_xysize; ++i) {
          labels[4 * i + 0] = 0;
          labels[4 * i + 1] = 0;
          labels[4 * i + 2] = histogram_symbols[i];
          labels[4 * i + 3] = 0;
        }
        if (!IsArrayMadeOfLabels(labels.data(), labels.size())) {
          assert(false);
        }

        WP2::ANSDebugPrefix prefix_histogram_image(enc, "histogram_image");
        enc->PutRange(histogram_bits, kHistogramBitsMin, max_sampling,
                      "histogram_bits");
        enc->PutRange(num_clusters, 2,
                      std::min(kMaxHistogramImageSize, histogram_image_xysize),
                      "num_histograms_m2");
        LosslessSymbolsInfo symbols_info_tmp;
        WP2_CHECK_STATUS(symbols_info_tmp.CopyFrom(symbols_info));
        symbols_info_tmp.SetCacheRange(0);
        symbols_info_tmp.SetAsLabelImage();
        symbols_info_tmp.SetMinMax(kSymbolG, /*min=*/0,
                                   /*max=*/num_clusters - 1);
        WP2_CHECK_STATUS(EncodeImageNoClusters(
            enc, dicts, labels.data(), hash_chain, ref_pool, histogram_width,
            histogram_height, symbols_info_tmp, speed));
      }

      // Get the symbols stats.
      const uint32_t num_pixels = width * height;
      WP2::ANSEncNoop enc_noop;
      WP2::SymbolRecorder recorder;
      symbols_info.SetNumClusters(num_clusters);
      WP2_CHECK_STATUS(recorder.Allocate(symbols_info, /*num_records=*/0));
      WP2_CHECK_STATUS(StorePixels(width, height, histogram_bits, *refs_best,
                                   histogram_symbols.data(), &enc_noop,
                                   &recorder));
      // Store symbol headers.
      WP2_CHECK_STATUS(
          WriteHeaders(recorder, symbols_info, num_pixels, enc, dicts, &sw));
    }

    // Store actual literals.
    EncodeInfo current_encode_info;
    if (encode_info != nullptr) {
      WP2_CHECK_ALLOC_OK(current_encode_info.line_tokens.resize(
          encode_info->line_tokens.size()));
      WP2_CHECK_ALLOC_OK(current_encode_info.bits_per_pixel.resize(
          encode_info->bits_per_pixel.size()));
    }
    WP2_CHECK_STATUS(StorePixels(width, height, histogram_bits, *refs_best,
                                 histogram_symbols.data(), enc, &sw,
                                 &current_encode_info));

    // Keep track of the smallest image so far.
    if (enc->GetCost(*dicts) < best_cost) {
      WP2_CHECK_STATUS(symbols_info_best->CopyFrom(symbols_info));
      best_cost = enc->GetCost(*dicts);
      WP2::swap(*enc, enc_best);
      WP2::swap(*dicts, dicts_best);
      if (encode_info != nullptr) {
        swap(*encode_info, current_encode_info);
      }
    }
    // Reset the bit writer for the following iteration if any.
    if (config.lz77s_types_to_try_size > 1) {
      WP2_CHECK_STATUS(enc->Clone(enc_init));
      dicts->DeepClear();
      WP2_CHECK_STATUS(dicts->CopyFrom(dicts_init));
    }
  }
  WP2::swap(*enc, enc_best);
  WP2::swap(*dicts, dicts_best);

  return WP2_STATUS_OK;
}

// -----------------------------------------------------------------------------
// Transforms

static void ApplySubtractGreen(Encoder* const encoder, uint32_t width,
                               uint32_t height, WP2::ANSEnc* const enc) {
  WP2LSubtractGreenFromBlueAndRed(
      encoder->argb_, encoder->symbols_info_.ChannelBits(), width * height);
}

// Finds the best spatial predictor for each prediction block (of size
// 2^(enc->transform_bits_)), saves the residuals in enc->argb_ and writes the
// transform for each block in the form of an image.
static WP2Status ApplyPredictFilter(const Encoder* const encoder,
                                    uint32_t width, uint32_t height, int speed,
                                    WP2::ANSEnc* const enc,
                                    WP2::ANSDictionaries* const dicts) {
  const int pred_bits = encoder->transform_bits_;
  const int transform_width = WP2LSubSampleSize(width, pred_bits);
  const int transform_height = WP2LSubSampleSize(height, pred_bits);

  WP2_CHECK_STATUS(ResidualImage(
      width, height, pred_bits, encoder->symbols_info_.ChannelBits(), speed,
      encoder->argb_, encoder->has_alpha_, encoder->argb_scratch_,
      encoder->transform_data_));
  const WP2::ANSDebugPrefix prefix(enc, "predictor");
  enc->PutRange(pred_bits, kTransformBitsMin, kTransformBitsMax,
                "transform_bits");
  LosslessSymbolsInfo symbols_info;
  symbols_info.Init(/*has_alpha=*/false, encoder->symbols_info_.SampleFormat());
  symbols_info.SetMinMax(kSymbolG, /*min=*/0, /*max=*/kNumPredictors - 1);
  symbols_info.SetAsLabelImage();
  const WP2Status status = EncodeImageNoClusters(
      enc, dicts, encoder->transform_data_, (HashChain*)&encoder->hash_chain_,
      (BackwardRefsPool*)&encoder->ref_pool_, transform_width, transform_height,
      symbols_info, speed);
  return status;
}

// Finds the best cross channel transform for each prediction block (of size
// 2^(enc->transform_bits_)), saves the residuals in enc->argb_ and writes the
// transform for each block in the form of an image.
static WP2Status ApplyCrossColorFilter(const Encoder* const encoder,
                                       uint32_t width, uint32_t height,
                                       int speed, WP2::ANSEnc* const enc,
                                       WP2::ANSDictionaries* const dicts) {
  const int ccolor_transform_bits = encoder->transform_bits_;
  const int transform_width = WP2LSubSampleSize(width, ccolor_transform_bits);
  const int transform_height = WP2LSubSampleSize(height, ccolor_transform_bits);

  WP2_CHECK_STATUS(ColorSpaceTransform(width, height, ccolor_transform_bits,
                                       encoder->symbols_info_.SampleFormat(),
                                       speed, encoder->argb_,
                                       encoder->transform_data_));
  const WP2::ANSDebugPrefix prefix(enc, "cross_color");
  enc->PutRange(ccolor_transform_bits, kTransformBitsMin, kTransformBitsMax,
                "transform_bits");
  LosslessSymbolsInfo symbols_info;
  symbols_info.Init(/*has_alpha=*/false, encoder->symbols_info_.SampleFormat());
  symbols_info.SetAsCrossColorImage();
  WP2_CHECK_STATUS(EncodeImageNoClusters(
      enc, dicts, encoder->transform_data_, (HashChain*)&encoder->hash_chain_,
      (BackwardRefsPool*)&encoder->ref_pool_, transform_width, transform_height,
      symbols_info, speed));
  return WP2_STATUS_OK;
}

// -----------------------------------------------------------------------------

// Allocates the memory for argb (W x H) buffer, 2 rows of context for
// prediction and transform data.
// Flags influencing the memory allocated:
//  transform_bits_
//  use_predict_, use_cross_color_
WP2Status Encoder::AllocateTransformBuffer(const CrunchConfig& config,
                                           uint32_t width, uint32_t height) {
  const uint64_t image_size = width * height;
  // ResidualImage needs room for 2 scanlines of uint32 pixels with an extra
  // pixel in each, plus 2 regular scanlines of bytes.
  // TODO(skal): Clean up by using arithmetic in bytes instead of words.
  const uint64_t argb_scratch_size =
      config.use_predict
          ? (width + 1) * 2 +
                (width * 2 + sizeof(uint32_t) - 1) / sizeof(uint32_t)
          : 0;
  const uint64_t transform_data_size =
      (config.use_predict || config.use_cross_color)
          ? WP2LSubSampleSize(width, transform_bits_) *
                WP2LSubSampleSize(height, transform_bits_)
          : 0;
  const uint64_t mem_size =
      image_size + argb_scratch_size + transform_data_size;
  uint16_t* mem;
  if (4 * mem_size > transform_mem_.size()) {
    WP2_CHECK_ALLOC_OK(transform_mem_.resize(4 * mem_size));
    mem = transform_mem_.data();
    argb_content_ = kEncoderNone;
  } else {
    mem = transform_mem_.data();
  }
  argb_ = mem;
  mem += 4 * image_size;
  argb_scratch_ = mem;
  mem += 4 * argb_scratch_size;
  transform_data_ = reinterpret_cast<int16_t*>(mem);

  return WP2_STATUS_OK;
}

static WP2Status MakeInputImageCopy(const CrunchConfig& config,
                                    Encoder* const encoder) {
  const WP2::ArgbBuffer& pic = encoder->pic_;
  const uint32_t width = pic.width;
  const uint32_t height = pic.height;
  WP2_CHECK_STATUS(encoder->AllocateTransformBuffer(config, width, height));
  if (encoder->argb_content_ == kEncoderARGB) return WP2_STATUS_OK;
  assert(WP2FormatBpc(pic.format) == 1 || WP2FormatBpc(pic.format) == 2);
  for (uint32_t y = 0; y < height; ++y) {
    if (WP2FormatBpc(pic.format) == 1) {
      std::copy(pic.GetRow8(y), pic.GetRow8(y) + 4 * width,
                &encoder->argb_[4 * y * width]);
    } else {
      std::copy(pic.GetRow16(y), pic.GetRow16(y) + 4 * width,
                &encoder->argb_[4 * y * width]);
    }
  }
  encoder->argb_content_ = kEncoderARGB;
  return WP2_STATUS_OK;
}

// -----------------------------------------------------------------------------

// Note: Expects "enc->palette_" to be set properly.
static WP2Status MapImageFromPalette(const CrunchConfig& config,
                                     Encoder* const encoder) {
  const WP2::ArgbBuffer& pic = encoder->pic_;
  const uint32_t width = pic.width;
  const uint32_t height = pic.height;

  WP2_CHECK_STATUS(encoder->AllocateTransformBuffer(config, width, height));

  encoder->argb_content_ = kEncoderPalette;
  WP2Status status = encoder->palette_.Apply(encoder->pic_, encoder->argb_);
  return status;
}

// Save palette_[] to bitstream.
static WP2Status EncodePalette(WP2::ANSEnc* const enc,
                               WP2::ANSDictionaries* const dicts,
                               Encoder* const encoder) {
  return encoder->palette_.Write(enc, dicts, encoder);
}

// Applies the different transforms decided by the cruncher and writes them to
// the header.
static WP2Status ApplyAndWriteTransforms(
    uint32_t width, uint32_t height, uint32_t speed, const CrunchConfig& config,
    Encoder* const encoder, LosslessSymbolsInfo* const symbols_info,
    WP2::ANSEnc* const enc, WP2::ANSDictionaries* const dicts) {
  WP2_CHECK_STATUS(symbols_info->CopyFrom(encoder->symbols_info_));
  symbols_info->SetCacheRange(0);

  encoder->argb_content_ = kEncoderNone;

  const WP2::ANSDebugPrefix prefix1(enc, "GlobalHeader");
  const WP2::ANSDebugPrefix prefix2(enc, "transforms");
  if (enc->PutBool(config.use_palette || config.use_subtract_green ||
                       config.use_predict || config.use_cross_color ||
                       config.use_group4,
                   "transform_present")) {
    if (!enc->PutBool(config.use_palette, "use_color_indexing") &&
        !enc->PutBool(config.use_group4, "use_group4")) {
      // Color indexing cannot be used with other transforms.
      enc->PutBool(config.use_subtract_green, "use_subtract_green");
      enc->PutBool(config.use_predict, "use_predictor");
      enc->PutBool(config.use_cross_color, "use_cross_color");
    }
  }

  // Encode palette
  if (config.use_palette) {
    symbols_info->SetAsLabelImage();
    const uint32_t palette_size = encoder->palette_.Size();
    symbols_info->SetMinMax(kSymbolG, /*min=*/0, /*max=*/palette_size - 1);
    WP2_CHECK_STATUS(EncodePalette(enc, dicts, encoder));
    WP2_CHECK_STATUS(MapImageFromPalette(config, encoder));
  }

  if (config.use_group4) {
    assert(encoder->palette_.Size() >= 2);
    WP2_CHECK_STATUS(EncodePalette(enc, dicts, encoder));
    WP2_CHECK_STATUS(MapImageFromPalette(config, encoder));
    // Green channel of first pixel.
    const uint16_t first_color = encoder->argb_[2];
    enc->PutRValue(first_color, encoder->palette_.Size(), "first_color");
    if (encoder->palette_.Size() > 3) {
      // MoveToFront is only relevant when there are more than 3 colors. With
      // 2 colors, any color change is necessarily the other color (different
      // from the current one). With 3 colors, we use a boolean to encode
      // whether the new color is the one expected (from the row above), or the
      // other one, and that is enough.
      enc->PutBool(config.group4_use_move_to_front, "use_move_to_front");
    } else {
      assert(!config.group4_use_move_to_front);
    }
  }

  // In case image is not packed, copy it to the argb_ buffer.
  if (encoder->argb_content_ != kEncoderNearLossless &&
      encoder->argb_content_ != kEncoderPalette) {
    WP2_CHECK_STATUS(MakeInputImageCopy(config, encoder));
  }

  // Apply transforms and write transform data.
  if (config.use_subtract_green) {
    ApplySubtractGreen(encoder, width, height, enc);
  }

  if (config.use_predict) {
    WP2_CHECK_STATUS(
        ApplyPredictFilter(encoder, width, height, speed, enc, dicts));
  }

  if (config.use_cross_color) {
    WP2_CHECK_STATUS(
        ApplyCrossColorFilter(encoder, width, height, speed, enc, dicts));
  }

  return WP2_STATUS_OK;
}

// -----------------------------------------------------------------------------
// Encoder

Encoder::Encoder(const WP2::EncoderConfig& config,
                 const WP2::ArgbBuffer& picture, bool has_alpha)
    : config_(config), pic_(picture), has_alpha_(has_alpha) {
  argb_ = nullptr;
  argb_content_ = kEncoderNone;
  argb_scratch_ = nullptr;
  transform_data_ = nullptr;

  histo_bits_ = 0;
  transform_bits_ = 0;
  // TODO(vrabaud) Deal with channel_bits
  symbols_info_.Init(has_alpha, picture.format);

  argb_content_ = kEncoderNone;
  WP2LEncDspInit();
}

WP2Status Encoder::Allocate() {
  const uint32_t num_pixels = pic_.width * pic_.height;
  WP2_CHECK_ALLOC_OK(hash_chain_.Allocate(num_pixels));
  palette_.Init(symbols_info_.SampleFormat(), has_alpha_);
  ref_pool_.Init(num_pixels);

  return WP2_STATUS_OK;
}

// -----------------------------------------------------------------------------
// Main call

typedef struct {
  const WP2::EncoderConfig* config;
  WP2::ANSEnc* enc;
  Encoder* encoder;
  CrunchConfig crunch_configs[kCrunchNumMax];
  int num_crunch_configs;
} StreamEncodeContext;

static WP2Status EncodeStreamHook(void* input, void* data2,
                                  EncodeInfo* const encode_info) {
  auto* const params = (StreamEncodeContext*)input;
  WP2::ANSEnc* const enc_init = params->enc;
  Encoder* const encoder = params->encoder;
  const CrunchConfig* const crunch_configs = params->crunch_configs;
  const uint32_t num_crunch_configs = params->num_crunch_configs;
  const uint32_t speed = params->config->speed;
  const uint32_t width = encoder->pic_.width;
  const uint32_t height = encoder->pic_.height;
  uint32_t best_size = std::numeric_limits<uint32_t>::max();
  WP2::ANSEnc enc;
  WP2::ANSEnc enc_best;
  (void)data2;

  for (const CrunchConfig* config = &crunch_configs[0];
       config < &crunch_configs[0] + num_crunch_configs; ++config) {
    // Reset the bit writer.
    WP2_CHECK_STATUS(enc.Clone(*enc_init));
    WP2::ANSDictionaries dicts;

    // Reset any parameter in the encoder that is set in the previous iteration.
    LosslessSymbolsInfo symbols_info;
    encoder->ref_pool_.Reset();

    WP2_CHECK_STATUS(ApplyAndWriteTransforms(
        width, height, speed, *config, encoder, &symbols_info, &enc, &dicts));

    // -------------------------------------------------------------------------
    // Encode and write the transformed image.
    EncodeInfo current_encode_info;
    if (encode_info != nullptr) {
      WP2_CHECK_ALLOC_OK(current_encode_info.line_tokens.resize(height));
      WP2_CHECK_ALLOC_OK(
          current_encode_info.bits_per_pixel.resize(width * height));
    }

    if (config->use_group4) {
      WP2_CHECK_STATUS(Group4Encode(
          encoder->argb_, width, height, encoder->palette_.Size(),
          config->group4_use_move_to_front, &enc, &current_encode_info));
    } else {
      // If using a color cache, do not have it bigger than the number of
      // colors.
      const uint32_t cache_bits_max =
          config->use_palette
              ? (1u + (uint32_t)WP2Log2Floor(encoder->palette_.Size()))
              : kMaxCacheBits;
      WP2_CHECK_STATUS(EncodeImageInternal(
          encoder->argb_, *config, width, height, encoder->histo_bits_,
          cache_bits_max, speed, &enc, &dicts, &encoder->hash_chain_,
          &encoder->ref_pool_, &symbols_info, &current_encode_info));
    }

    // If we are better than what we already have.
    WP2_CHECK_STATUS(enc.Assemble());
    if (enc.BufferSize() < best_size) {
      best_size = enc.BufferSize();
      // Store the BitWriter.
      WP2::swap(enc, enc_best);
      if (encode_info != nullptr) {
        swap(*encode_info, current_encode_info);
      }
    }
  }
  WP2::swap(*enc_init, enc_best);

  return WP2_STATUS_OK;
}

WP2Status EncodeImage(const WP2::EncoderConfig& config,
                      const WP2::ArgbBuffer& picture, bool has_alpha,
                      WP2::ANSEnc* const enc, EncodeInfo* const encode_info) {
  assert(!picture.IsEmpty());
  assert(enc != nullptr);
  Encoder encoder(config, picture, has_alpha);
  CrunchConfig crunch_configs[kCrunchNumMax];
  uint32_t num_crunch_configs = 0;

  // Analyze image (entropy, num_palettes etc)
  WP2_CHECK_STATUS(encoder.Allocate());
  WP2_CHECK_STATUS(
      EncoderAnalyze(&encoder, crunch_configs, &num_crunch_configs));

  StreamEncodeContext params;
  for (uint32_t idx = 0; idx < num_crunch_configs; ++idx) {
    params.crunch_configs[idx] = crunch_configs[idx];
  }
  params.num_crunch_configs = num_crunch_configs;

  // Fill in the parameters for the thread workers.
  params.config = &config;
  params.enc = enc;
  params.encoder = &encoder;
  return EncodeStreamHook(&params, nullptr, encode_info);
}

//------------------------------------------------------------------------------

bool EncodeInfo::CopyFrom(const EncodeInfo& other) {
  return line_tokens.copy_from(other.line_tokens) &&
         bits_per_pixel.copy_from(other.bits_per_pixel);
}

//------------------------------------------------------------------------------

}  // namespace WP2L
