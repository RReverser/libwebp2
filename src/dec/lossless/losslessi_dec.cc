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
// main entry for the decoder
//
// Author: Vincent Rabaud (vrabaud@google.com)

#include "src/dec/lossless/losslessi_dec.h"

#include <algorithm>
#include <string>

#include "src/common/color_precision.h"
#include "src/common/constants.h"
#include "src/common/lossless/color_cache.h"
#include "src/common/progress_watcher.h"
#include "src/common/symbols.h"
#include "src/dec/wp2_dec_i.h"
#include "src/dsp/lossless/lossless.h"
#include "src/dsp/lossless/lossless_common.h"
#include "src/dsp/math.h"
#include "src/utils/ans_utils.h"
#include "src/utils/plane.h"
#include "src/utils/vector.h"
#include "src/wp2/format_constants.h"

namespace WP2L {

// -----------------------------------------------------------------------------

static const uint8_t kCodeToPlane[kCodeToPlaneCodes] = {
  0x18, 0x07, 0x17, 0x19, 0x28, 0x06, 0x27, 0x29, 0x16, 0x1a,
  0x26, 0x2a, 0x38, 0x05, 0x37, 0x39, 0x15, 0x1b, 0x36, 0x3a,
  0x25, 0x2b, 0x48, 0x04, 0x47, 0x49, 0x14, 0x1c, 0x35, 0x3b,
  0x46, 0x4a, 0x24, 0x2c, 0x58, 0x45, 0x4b, 0x34, 0x3c, 0x03,
  0x57, 0x59, 0x13, 0x1d, 0x56, 0x5a, 0x23, 0x2d, 0x44, 0x4c,
  0x55, 0x5b, 0x33, 0x3d, 0x68, 0x02, 0x67, 0x69, 0x12, 0x1e,
  0x66, 0x6a, 0x22, 0x2e, 0x54, 0x5c, 0x43, 0x4d, 0x65, 0x6b,
  0x32, 0x3e, 0x78, 0x01, 0x77, 0x79, 0x53, 0x5d, 0x11, 0x1f,
  0x64, 0x6c, 0x42, 0x4e, 0x76, 0x7a, 0x21, 0x2f, 0x75, 0x7b,
  0x31, 0x3f, 0x63, 0x6d, 0x52, 0x5e, 0x00, 0x74, 0x7c, 0x41,
  0x4f, 0x10, 0x20, 0x62, 0x6e, 0x30, 0x73, 0x7d, 0x51, 0x5f,
  0x40, 0x72, 0x7e, 0x61, 0x6f, 0x50, 0x71, 0x7f, 0x60, 0x70
};

//------------------------------------------------------------------------------

static inline int GetCopyDistance(uint32_t distance_symbol,
                                  WP2::ANSDec* const dec, double* const cost,
                                  const char* const label = "dist_bits") {
  if (distance_symbol < 4) {
    return distance_symbol + 1;
  }
  const uint32_t extra_bits = (distance_symbol - 2) >> 1;
  const int offset = (2 + (distance_symbol & 1)) << extra_bits;

  if (extra_bits > 0) {
    int extra_val;
    if (extra_bits <= IO_BITS) {
      extra_val = dec->ReadUValue(extra_bits, label);
    } else {
      extra_val = dec->ReadUValue(IO_BITS, "dist_bits1");
      extra_val |= dec->ReadUValue(extra_bits - IO_BITS, "dist_bits2")
                   << IO_BITS;
    }
    if (cost != nullptr) {
      *cost += extra_bits;
    }
    return offset + extra_val + 1;
  } else {
    return offset + 1;
  }
}

static inline int GetCopyLength(uint32_t length_symbol, WP2::ANSDec* const dec,
                                double* const cost) {
  // Length and distance prefixes are encoded the same way.
  return GetCopyDistance(length_symbol, dec, cost, "len_bits");
}

static inline int PlaneCodeToDistance(int xsize, int plane_code) {
  if (plane_code > kCodeToPlaneCodes) {
    return plane_code - kCodeToPlaneCodes;
  } else {
    const int dist_code = kCodeToPlane[plane_code - 1];
    const int yoffset = dist_code >> 4;
    const int xoffset = 8 - (dist_code & 0xf);
    const int dist = yoffset * xsize + xoffset;
    return (dist >= 1) ? dist : 1;  // dist<1 can happen if xsize is very small
  }
}

//------------------------------------------------------------------------------

// Reads packed symbol depending on GREEN channel
WP2Status Decoder::ReadANSStats(uint32_t width, uint32_t height,
                                const LosslessSymbolsInfo& symbols_info) {
  WP2::ANSDebugPrefix prefix(dec_, "symbols");

  WP2::SymbolReader* const sr = &hdr_.sr_;
  WP2_CHECK_STATUS(sr->Init(symbols_info, dec_));
  WP2_CHECK_STATUS(sr->Allocate());
  const uint32_t num_clusters = symbols_info.NumClusters(0);
  for (uint32_t cluster = 0; cluster < num_clusters; ++cluster) {
    // Deal with the first symbol.
    WP2_CHECK_STATUS(sr->ReadHeader(cluster, width * height, kSymbolType,
                                    kSymbolNames[kSymbolType]));
    bool is_maybe_used[kSymbolTypeNum];
    sr->GetPotentialUsage(cluster, kSymbolType, is_maybe_used, kSymbolTypeNum);
    for (int s = 1; s < kSymbolNum; ++s) {
      const Symbol sym = (Symbol)(s);
      if (symbols_info.IsSymbolUseless(sym, is_maybe_used)) continue;
      WP2_CHECK_STATUS(
          sr->ReadHeader(cluster, width * height, sym, kSymbolNames[s]));
    }
  }
  return WP2_STATUS_OK;
}

//------------------------------------------------------------------------------

// Emit rows without any scaling.
template <typename T>
uint32_t EmitRows(const uint16_t* row_in, uint32_t mb_w, uint32_t mb_h,
                  bool has_alpha, T* const out, uint32_t out_stride) {
  int lines = mb_h;
  T* row_out = out;
  while (lines-- > 0) {
    // TODO(vrabaud) support multiple output color spaces.
    for (size_t i = 0; i < 4 * mb_w; i += 4, row_in += 4) {
      const uint8_t alpha = (uint8_t)row_in[0];
      if (!has_alpha) assert(alpha == WP2::kAlphaMax);
      row_out[i] = alpha;
      if (alpha > 0) {
        for (int j = 1; j <= 3; ++j) row_out[i + j] = (T)row_in[j];
      } else {
        // To optimize encoding, the encoder can store arbitrary rgb values when
        // a = 0 (fully transparent pixel, see GetResidual in predictor_enc.cc).
        // To make sure we return proper pre-multiplied Argb, we sanitize these
        // values by setting rgb to 0.
        for (int j = 1; j <= 3; ++j) row_out[i + j] = 0;
      }
    }
    row_out = (T*)(((uint8_t*)row_out) + out_stride);
  }
  return mb_h;  // Num rows out == num rows in.
}

//------------------------------------------------------------------------------

static inline int GetMetaIndex(const WP2::Vector_s16& image, uint32_t width,
                               int bits, int x, int y) {
  if (bits == 0 || image.empty()) return 0;
  return image[4 * (width * (y >> bits) + (x >> bits)) + 2];
}

//------------------------------------------------------------------------------
// Main loop, with custom row-processing function

WP2Status Decoder::ApplyInverseTransforms(uint32_t num_rows, uint32_t num_bits,
                                          const int16_t* const rows) {
  uint32_t n = next_transform_;
  const uint32_t cache_pixs = tile_->rect.width * num_rows;

  if (n == 0) {
    // No transform called, hence just copy.
    std::copy(rows, rows + 4 * cache_pixs, argb_cache_);
    return WP2_STATUS_OK;
  }

  const uint32_t start_row = last_row_;
  const uint32_t end_row = start_row + num_rows;

  // TODO(vrabaud) remove when switching to int16_t.
  WP2::Vector_u16 data_vec;
  WP2_CHECK_ALLOC_OK(data_vec.resize(4 * cache_pixs));
  std::copy(rows, rows + 4 * cache_pixs, data_vec.data());

  const uint16_t* rows_in = data_vec.data();
  uint16_t* const rows_out = argb_cache_;
  // Inverse transforms.
  while (n-- > 0) {
    const Transform& transform = transforms_[n];
    if (transform.type_ == WP2L::COLOR_INDEXING_TRANSFORM ||
        transform.type_ == WP2L::GROUP4) {
      // Check that all color map indices are valid here rather than in dsp.
      // TODO(yguyon): Remove pessimization.
      const uint32_t color_map_size = (uint32_t)transform.data_.size() / 4;
      const uint16_t* const r_end = rows_in + 4 * num_rows * transform.width_;
      for (const uint16_t* r = rows_in; r < r_end; r += 4) {
        WP2_CHECK_OK(r[2] < color_map_size, WP2_STATUS_BITSTREAM_ERROR);
      }
    }
    WP2LInverseTransform(&transform, start_row, end_row, num_bits, rows_in,
                         gparams_->has_alpha_, rows_out);
    rows_in = rows_out;
  }
  return WP2_STATUS_OK;
}

WP2Status Decoder::ProcessRows(uint32_t row,
                               WP2::DecoderInfo* const decoder_info) {
  const int16_t* const rows =
      pixels_.data() + 4 * tile_->rect.width * last_row_;

  if (row > last_row_) {    // Emit output.
    const uint32_t num_rows = row - last_row_;
    RegisterUnprocessedRowForVDebug(num_rows, num_bits_, rows, decoder_info);
    // We can't process more than kNumARGBCacheRows at a time (that's the size
    // of argb_cache_), but we currently don't need more than that.
    assert(num_rows <= kNumARGBCacheRows);
    WP2::ArgbBuffer* const buf = &tile_->rgb_output;
    assert(!buf->IsEmpty());
    uint16_t* const rows_data = argb_cache_;

    WP2_CHECK_STATUS(ApplyInverseTransforms(num_rows, num_bits_, rows));

    uint32_t num_rows_out;
    if (WP2FormatBpc(buf->format) == 1) {
      num_rows_out =
          EmitRows(rows_data, buf->width, num_rows, gparams_->has_alpha_,
                   buf->GetRow8(last_out_row_), buf->stride);
    } else {
      num_rows_out =
          EmitRows(rows_data, buf->width, num_rows, gparams_->has_alpha_,
                   buf->GetRow16(last_out_row_), buf->stride);
    }
    // Update 'last_out_row_'.
    last_out_row_ += num_rows_out;
    tile_->num_decoded_rows += num_rows_out;
    if (progress_ != nullptr) {
      WP2_CHECK_STATUS(progress_->AdvanceBy(num_rows_out * buf->width));
    }
  }

  // Update 'last_row_'.
  last_row_ = row;
  assert(last_row_ <= tile_->rect.height);
  return WP2_STATUS_OK;
}

//------------------------------------------------------------------------------
// Helper functions for fast pattern copy.

static inline void CopyBlock(int16_t* dst, uint32_t dist, uint32_t length) {
  uint32_t length_written = 0;
  const int16_t* const dst_ini = dst - 4 * dist;
  int16_t* dst_end = dst;
  while (length_written < length) {
    // We cannot write more than what is left, and we cannot write further than
    // dst_end.
    const uint32_t length_to_write =
        std::min(length - length_written, (uint32_t)((dst_end - dst_ini) / 4));
    std::copy(dst_ini, dst_ini + 4 * length_to_write, dst_end);
    dst_end += 4 * length_to_write;
    length_written += length_to_write;
  }
}

//------------------------------------------------------------------------------

typedef void (*ProcessRowsFunc)(Decoder* const dec, uint32_t row);

WP2Status Decoder::DecodeAndProcess(uint32_t width, uint32_t last_row,
                                    bool do_process,
                                    WP2::Vector_s16* const data_vec,
                                    WP2::Planef* const bits_per_pixel,
                                    WP2::DecoderInfo* const decoder_info,
                                    int16_t** src_out) {
  int16_t* const data = data_vec->data();
  uint32_t row = last_pixel_ / width;
  uint32_t col = last_pixel_ % width;
  int16_t* src = data + 4 * last_pixel_;
  const int16_t* last_cached = src;
  const int16_t* const src_end = &data_vec->back() + 1;  // End of data
  int16_t* const src_last =
      data + 4 * width * last_row;  // Last pixel to decode
  ColorCache* const color_cache =
      do_process ? hdr_.color_cache_.get() : nullptr;
  const uint32_t mask = hdr_.histogram_mask_;
  WP2::SymbolReader* const sr = &hdr_.sr_;

  uint32_t cluster =
      (src < src_last)
          ? GetMetaIndex(hdr_.histogram_image_, hdr_.histogram_xsize_,
                         hdr_.histogram_subsample_bits_, col, row)
          : 0;
  assert(last_row_ < last_row);
  assert(src_last <= src_end);

  const bool has_alpha = LosslessSymbolsInfo::HasAlpha(sr->symbols_info());
  const bool is_label_image =
      LosslessSymbolsInfo::IsLabelImage(sr->symbols_info());

  const bool get_cost = (bits_per_pixel != nullptr || decoder_info != nullptr);
  uint32_t i_pixel = 0;
  while (src < src_last) {
    double cost = 0.;
    double* cost_ptr = get_cost ? &cost : nullptr;
    uint32_t length;
    const auto type =
        (SymbolType)sr->Read(cluster, kSymbolType, "lossless_mode", cost_ptr);
    WP2_CHECK_STATUS(dec_->GetStatus());
    if (type == kSymbolTypeLiteral) {
      WP2::ANSDebugPrefix prefix(dec_, "literal");
      length = 1;
      if (has_alpha) {
        src[0] = sr->Read(cluster, kSymbolA, "A", cost_ptr);
      } else {
        // For transform images/entropy image, alpha is not used: any value
        // works. The default value only matters for real ARGB samples.
        src[0] = WP2::kAlphaMax;
      }
      if (!gparams_->has_alpha_) assert(src[0] == WP2::kAlphaMax);
      if (is_label_image) {
        src[1] = src[3] = 0;
      } else {
        src[1] = sr->Read(cluster, kSymbolR, "R", cost_ptr);
        src[3] = sr->Read(cluster, kSymbolB, "B", cost_ptr);
      }
      src[2] = sr->Read(cluster, kSymbolG, "G", cost_ptr);
      WP2_CHECK_STATUS(dec_->GetStatus());
    } else if (type == kSymbolTypeCopy) {
      WP2::ANSDebugPrefix prefix(dec_, "copy");
      int dist_code, dist;
      const uint32_t length_sym =
          sr->Read(cluster, kSymbolLen, "len", cost_ptr);
      length = GetCopyLength(length_sym, dec_, cost_ptr);
      const uint32_t dist_symbol =
          sr->Read(cluster, kSymbolDist, "dist", cost_ptr);
      dist_code = GetCopyDistance(dist_symbol, dec_, cost_ptr);
      dist = PlaneCodeToDistance(width, dist_code);
      WP2_CHECK_STATUS(dec_->GetStatus());
      WP2_CHECK_OK(src - data >= (ptrdiff_t)(4 * dist) &&
                       src_end - src >= (ptrdiff_t)(4 * length),
                   WP2_STATUS_BITSTREAM_ERROR);
      CopyBlock(src, dist, length);
      if (!gparams_->has_alpha_) assert(src[0] == WP2::kAlphaMax);
    } else if (type == kSymbolTypeCacheIdx) {
      length = 1;
      assert(color_cache != nullptr);
      while (last_cached < src) {
        color_cache->Insert(last_cached, /*index_ptr=*/nullptr);
        last_cached += 4;
      }
      int32_t cache_index;
      WP2_CHECK_STATUS(sr->TryRead(cluster, kSymbolCache,
                                   color_cache->IndexRange() - 1, "cache_idx",
                                   &cache_index, cost_ptr));
      color_cache->Lookup(cache_index, src);
      if (!gparams_->has_alpha_) assert(src[0] == WP2::kAlphaMax);
    } else {                                            // Not reached
      WP2_CHECK_OK(false, WP2_STATUS_BITSTREAM_ERROR);  // For error reporting.
    }
    if (get_cost) {
      RegisterSymbolForVDebug(type, /*is_group4=*/false, last_pixel_ + i_pixel,
                              length, cost, decoder_info);
      if (bits_per_pixel != nullptr) {
        const double bits = cost / length;
        for (uint32_t j = 0; j < length; ++j) {
          const uint32_t pixel = last_pixel_ + i_pixel + j;
          bits_per_pixel->Row(pixel / width)[pixel % width] = bits;
        }
      }
    }
    i_pixel += length;
    src += 4 * length;
    col += length;
    while (col >= width) {
      col -= width;
      ++row;
      if (do_process) {
        if (row <= last_row && (row % kNumARGBCacheRows == 0)) {
          WP2_CHECK_STATUS(ProcessRows(row, config_.info));
        }
      }
    }

    // Only update when changing tile. Note we could use this test:
    // if "((((prev_col ^ col) | prev_row ^ row)) > mask)" -> tile changed
    // but that's actually slower and needs storing the previous col/row.
    if ((type == kSymbolTypeCopy || (col & mask) == 0) && src < src_end) {
      cluster = GetMetaIndex(hdr_.histogram_image_, hdr_.histogram_xsize_,
                             hdr_.histogram_subsample_bits_, col, row);
    }
  }

  // Finish filling the color cache.
  if (color_cache != nullptr) {
    while (last_cached < src) {
      color_cache->Insert(last_cached, /*index_ptr=*/nullptr);
      last_cached += 4;
    }
  }

  WP2_CHECK_STATUS(dec_->GetStatus());

  if (src_out != nullptr) *src_out = src;

  return WP2_STATUS_OK;
}

WP2Status Decoder::DecodeNoProcess(uint32_t width, uint32_t height,
                                   uint32_t last_row,
                                   WP2::Vector_s16* const data_vec) {
  WP2_CHECK_ALLOC_OK(data_vec->resize(4 * (uint64_t)width * height));
  WP2_CHECK_STATUS(
      DecodeAndProcess(width, last_row, /*do_process=*/false, data_vec));
  return WP2_STATUS_OK;
}

WP2Status Decoder::DecodeImageData(uint32_t width, uint32_t height,
                                   uint32_t last_row,
                                   WP2::Planef* const bits_per_pixel) {
  WP2_CHECK_ALLOC_OK(pixels_.resize(4 * (uint64_t)width * height));

  int16_t* src;
  if (use_group4_) {
    WP2_CHECK_STATUS(DecodeGroup4(width, last_row, bits_per_pixel, &src));
  } else {
    WP2_CHECK_STATUS(DecodeAndProcess(width, last_row, /*do_process=*/true,
                                      &pixels_, bits_per_pixel, config_.info,
                                      &src));
  }

  // Process the remaining rows corresponding to last row-block.
  const uint32_t row = (src - pixels_.data()) / 4 / width;
  WP2_CHECK_STATUS(ProcessRows(row > last_row ? last_row : row, config_.info));

  WP2_CHECK_STATUS(dec_->GetStatus());

  // When decoding the main image, update the last pixel, otherwise reset.
  int16_t* const data = pixels_.data();
  last_pixel_ = (uint32_t)(src - data) / 4;  // end-of-scan marker
  return WP2_STATUS_OK;
}

// -----------------------------------------------------------------------------
// Transform

static WP2Status ReadPaletteVerbatim(uint32_t num_colors,
                                     WP2SampleFormat format, bool has_alpha,
                                     bool is_grayscale, WP2::ANSDec* const dec,
                                     Transform* const transform) {
  WP2_CHECK_ALLOC_OK(transform->data_.resize(4 * num_colors));
  uint32_t channel_max[4];
  for (uint32_t c = 0; c < 4; ++c) {
    channel_max[c] = WP2::FormatMax(format, /*channel=*/c);
  }
  for (uint32_t i = 0; i < num_colors; ++i) {
    if (has_alpha) {
      transform->data_[4 * i + 0] =
          dec->ReadRValue(channel_max[0] + 1, "color_verbatim");
      for (uint32_t c = 1; c < 4; ++c) {
        const uint32_t max = WP2::DivRound(
            (uint32_t)(channel_max[c] * transform->data_[4 * i + 0]),
            channel_max[0]);
        transform->data_[4 * i + c] =
            (is_grayscale && c > 1)
                ? transform->data_[4 * i + 1]
                : dec->ReadRValue(max + 1, "color_verbatim");
      }
    } else {
      transform->data_[4 * i + 0] = channel_max[0];
      for (uint32_t c = 1; c < 4; ++c) {
        transform->data_[4 * i + c] =
            (is_grayscale && c > 1)
                ? transform->data_[4 * i + 1]
                : dec->ReadRValue(channel_max[c] + 1, "color_verbatim");
      }
    }
  }
  return WP2_STATUS_OK;
}

WP2Status Decoder::ReadPalette(Transform* const transform) {
  transform->bits_ = 0;

  dec_->AddDebugPrefix("palette");
  // There are at least 2 colors in the palette.
  num_colors_ = ReadGolomb(
      /*min=*/2, MaxPaletteSize(transform->width_ * transform->height_), 1,
      dec_, "num_colors");

  const bool is_grayscale = dec_->ReadBool("is_grayscale");
  const uint32_t num_channels = (is_grayscale ? 2 : 4);

  // Modify the main image ranges.
  hdr_.symbols_info_.SetAsLabelImage();
  hdr_.symbols_info_.SetMinMax(kSymbolG, /*min=*/0, /*max=*/num_colors_ - 1);

  if (num_colors_ <= kSmallPaletteLimit) {
    WP2_CHECK_STATUS(ReadPaletteVerbatim(
        num_colors_, hdr_.symbols_info_.SampleFormat(), gparams_->has_alpha_,
        is_grayscale, dec_, transform));
    dec_->PopDebugPrefix();
    return WP2_STATUS_OK;
  }

  // Fill the sign array.
  WP2::Vector_u8 is_positive[4];
  for (WP2::Vector_u8& vec : is_positive) {
    WP2_CHECK_ALLOC_OK(vec.resize(num_colors_ - 1));
  }
  for (size_t c = 0; c < num_channels; ++c) {
    if (c == 0 && !gparams_->has_alpha_) continue;
    const bool is_monotonous = dec_->ReadBool("is_monotonous");
    const bool use_opposite_sign = dec_->ReadBool("use_opposite_sign");
    if (is_monotonous) {
      std::fill(is_positive[c].begin(), is_positive[c].end(),
                !use_opposite_sign);
    } else {
      dec_->AddDebugPrefix("transform_index_signs");
      WP2::ReadVector(dec_, 1, is_positive[c]);
      dec_->PopDebugPrefix();
      if (use_opposite_sign) {
        for (auto& is_pos : is_positive[c]) is_pos = !is_pos;
      }
    }
  }

  if (dec_->ReadBool("use_image")) {
    // Read the header of the palette image.
    LosslessSymbolsInfo symbols_info;
    symbols_info.Init(gparams_->has_alpha_, hdr_.symbols_info_.SampleFormat());
    WP2_CHECK_STATUS(ReadANSStats(num_colors_, 1, symbols_info));

    // Read the palette image.
    WP2_CHECK_STATUS(
        DecodeNoProcess(num_colors_, 1, /*last_row=*/1, &transform->data_));
    WP2_CHECK_STATUS(dec_->GetStatus());

    // Convert the diffs to proper colors.
    int16_t* const data = transform->data_.data();
    // data[0] stays unchanged and serves as a basis to adding successive
    // differences.
    if (!gparams_->has_alpha_) assert(data[0] == WP2::kAlphaMax);
    for (size_t i = 1; i < num_colors_; ++i) {
      for (size_t c = 0; c < num_channels; ++c) {
        if (c == 0 && !gparams_->has_alpha_) {
          assert(data[4 * i + c] == WP2::kAlphaMax);
          continue;
        }
        const int16_t alpha = data[4 * i + 0];
        int16_t channel_prev = data[4 * (i - 1) + c];
        if (c > 0 && channel_prev > alpha) channel_prev = alpha;
        const auto channel_diff = data[4 * i + c];
        if (is_positive[c][i - 1]) {
          data[4 * i + c] = channel_prev + channel_diff;
        } else {
          data[4 * i + c] = channel_prev - channel_diff;
        }
      }
    }
  } else {
    WP2_CHECK_ALLOC_OK(transform->data_.resize(4 * num_colors_));
    std::fill(transform->data_.begin(), transform->data_.end(), 0);
    // Read the channels one by one.
    const uint32_t alpha_max =
        WP2::FormatMax(hdr_.symbols_info_.SampleFormat(), 0);
    for (size_t c = 0; c < num_channels; ++c) {
      const uint32_t channel_max =
          WP2::FormatMax(hdr_.symbols_info_.SampleFormat(), c);
      if (c == 0 && !gparams_->has_alpha_) {
        for (size_t i = 0; i < num_colors_; ++i) {
          transform->data_[4 * i + 0] = channel_max;
        }
        continue;
      }
      const bool use_vector = dec_->ReadBool("use_vector");
      // Read the first channel.
      const uint32_t max_first = (c == 0) ? channel_max : transform->data_[0];
      int16_t channel = dec_->ReadRValue(max_first + 1, "first_channel");
      transform->data_[c] = channel;

      if (use_vector) {
        WP2::Vector_u16 diffs;
        WP2_CHECK_ALLOC_OK(diffs.resize(num_colors_ - 1));
        dec_->AddDebugPrefix("transform_index_vals");
        WP2::ReadVector(dec_, channel_max, diffs);
        dec_->PopDebugPrefix();
        for (size_t i = 1; i < num_colors_; ++i) {
          if (gparams_->has_alpha_) {
            const int16_t alpha = transform->data_[4 * i + 0];
            if (c > 0 && channel > alpha) channel = alpha;
          }
          channel = WP2::Clamp<int32_t>(channel +
              (is_positive[c][i - 1] ? diffs[i - 1] : -diffs[i - 1]),
              0, channel_max);
          transform->data_[4 * i + c] = channel;
        }
      } else {
        for (size_t i = 1; i < num_colors_; ++i) {
          const int16_t alpha = transform->data_[4 * i + 0];
          assert(alpha >= 0 && alpha <= (int16_t)alpha_max);
          if (gparams_->has_alpha_ && c > 0 && channel > alpha) channel = alpha;
          if (is_positive[c][i - 1]) {
            const uint32_t max =
                (c == 0) ? channel_max
                         : WP2::DivRound(channel_max * alpha, alpha_max);
            channel += dec_->ReadRValue(max + 1 - channel, "store_diff");
          } else {
            channel -= dec_->ReadRValue(channel + 1, "store_diff");
          }
          transform->data_[4 * i + c] = channel;
        }
      }
    }
  }
  // Duplicate the G/B channels if needed.
  if (is_grayscale) {
    for (size_t i = 0; i < num_colors_; ++i) {
      const uint32_t channel = transform->data_[4 * i + 1];
      transform->data_[4 * i + 2] = channel;
      transform->data_[4 * i + 3] = channel;
    }
  }
  dec_->PopDebugPrefix();
  return WP2_STATUS_OK;
}

WP2Status Decoder::ReadTransform(uint32_t width, uint32_t height,
                                ImageTransformType type) {
  Transform* const transform = &transforms_[next_transform_];

  transform->type_ = type;
  transform->width_ = width;
  transform->height_ = height;
  ++next_transform_;
  assert(next_transform_ <= NUM_TRANSFORMS);

  switch (type) {
    case PREDICTOR_TRANSFORM:
    case CROSS_COLOR_TRANSFORM: {
      if (type == CROSS_COLOR_TRANSFORM) {
        dec_->AddDebugPrefix("cross_color");
      } else {
        assert(type == PREDICTOR_TRANSFORM);
        dec_->AddDebugPrefix("predictor");
      }
      transform->bits_ = dec_->ReadRange(kTransformBitsMin, kTransformBitsMax,
                                         "transform_bits");

      // Read the transform header.
      LosslessSymbolsInfo symbols_info;
      WP2_CHECK_STATUS(symbols_info.CopyFrom(hdr_.symbols_info_));
      symbols_info.SetInfo(kSymbolA, /*min=*/WP2::SymbolsInfo::kInvalidBound,
                           /*max=*/WP2::SymbolsInfo::kInvalidBound,
                           /*num_clusters=*/0,
                           WP2::SymbolsInfo::StorageMethod::kUnused);
      if (type == PREDICTOR_TRANSFORM) {
        symbols_info.SetAsLabelImage();
        symbols_info.SetMinMax(kSymbolG, /*min=*/0, /*max=*/kNumPredictors - 1);
      } else {
        symbols_info.SetAsCrossColorImage();
      }
      symbols_info.SetNumClusters(1);
      symbols_info.SetCacheRange(0);
      const uint32_t
          sub_width = WP2LSubSampleSize(transform->width_, transform->bits_);
      const uint32_t
          sub_height = WP2LSubSampleSize(transform->height_, transform->bits_);
      WP2_CHECK_STATUS(ReadANSStats(sub_width, sub_height, symbols_info));

      // Read the transform image.
      WP2_CHECK_STATUS(DecodeNoProcess(
          sub_width, sub_height, /*last_row=*/
          WP2LSubSampleSize(transform->height_, transform->bits_),
          &transform->data_));
      WP2_CHECK_STATUS(dec_->GetStatus());

      dec_->PopDebugPrefix();
      break;
    }
    case COLOR_INDEXING_TRANSFORM: {
      WP2_CHECK_STATUS(ReadPalette(transform));
      break;
    }
    case GROUP4: {
      WP2_CHECK_STATUS(ReadPalette(transform));
      const uint32_t color_map_size = (uint32_t)transform->data_.size() / 4;
      group4_first_color_ = dec_->ReadRValue(color_map_size, "first_color");
      const bool use_move_to_front =
          (num_colors_ > 3) ? dec_->ReadBool("use_move_to_front") : false;
      WP2_CHECK_STATUS(mtf_.Init(use_move_to_front, num_colors_));

      WP2::SymbolsInfo info;
      InitGroup4(width, color_map_size, &info);

      // Read the headers.
      WP2::SymbolReader* const sr = &hdr_.sr_;
      WP2_CHECK_STATUS(sr->Init(info, dec_));
      WP2_CHECK_STATUS(sr->Allocate());

      {
        WP2::ANSDebugPrefix prefix(dec_, "group4");
        for (uint32_t s = 0; s < (uint32_t)kSymbolG4Num; ++s) {
          if (color_map_size <= 2 &&
              (s == kSymbolG4ColorChange || s == kSymbolG4NewColor)) {
            continue;
          }
          const uint32_t max_nnz = width * height;
          WP2_CHECK_STATUS(sr->ReadHeader(max_nnz, s, kSymbolGroup4Names[s]));
        }
      }
      break;
    }
    case SUBTRACT_GREEN:
      break;
    default:
      assert(0);    // can't happen
      break;
  }

  return WP2_STATUS_OK;
}

// -----------------------------------------------------------------------------
// Metadata

void Metadata::Init(bool has_alpha, WP2SampleFormat format) {
  symbols_info_.Init(has_alpha, format);
  histogram_mask_ = 0;
  histogram_subsample_bits_ = 0;
  histogram_xsize_ = 0;
}

// -----------------------------------------------------------------------------
// Decoder

Decoder::Decoder() {
  WP2LDspInit();  // Init critical function pointers.
}

void Decoder::Init(const WP2::DecoderConfig& config,
                   const WP2::GlobalParams& gparams,
                   WP2::ProgressWatcher* const progress, WP2::ANSDec* const dec,
                   WP2::Tile* const tile) {
  config_ = config;
  gparams_ = &gparams;
  progress_ = progress;
  dec_ = dec;
  last_out_row_ = 0;
  tile_ = tile;
  tile->num_decoded_rows = 0;
  hdr_.Init(gparams.has_alpha_, tile->rgb_output.format);
}

WP2Status Decoder::AllocateInternalBuffers(uint32_t final_width) {
  const uint64_t num_pixels = (uint64_t)tile_->rect.GetArea();
  // Scratch buffer corresponding to top-prediction row for transforming the
  // first row in the row-blocks. Not needed for paletted alpha.
  const uint64_t cache_top_pixels = final_width;
  // Scratch buffer for temporary BGRA storage. Not needed for paletted alpha.
  const uint64_t cache_pixels = (uint64_t)final_width * kNumARGBCacheRows;
  const uint64_t total_num_pixels =
      num_pixels + cache_top_pixels + cache_pixels;

  assert(tile_->rect.width <= final_width);
  argb_cache_ = nullptr;  // for sanity check
  WP2_CHECK_ALLOC_OK(pixels_.resize(4 * total_num_pixels));
  // TODO(vrabaud) Remove the cast once switched to int16_t.
  argb_cache_ = reinterpret_cast<uint16_t*>(pixels_.data()) +
                4 * (num_pixels + cache_top_pixels);

  if (config_.info != nullptr && config_.info->bits_per_pixel != nullptr) {
    WP2_CHECK_STATUS(
        bits_per_pixel_.SetView(*config_.info->bits_per_pixel, tile_->rect));
  } else if (VDMatch(config_, "bits-per-pixel")) {
    WP2_CHECK_STATUS(
        bits_per_pixel_.Resize(tile_->rect.width, tile_->rect.height));
  }

  return WP2_STATUS_OK;
}

WP2Status Decoder::DecodeHeader() {
  state_ = READ_HDR;

  const WP2::ANSDebugPrefix prefix1(dec_, "GlobalHeader");

  // Read the transforms.
  use_group4_ = false;
  {
    const WP2::ANSDebugPrefix prefix2(dec_, "transforms");
    if (dec_->ReadBool("transform_present")) {
      if (dec_->ReadBool("use_color_indexing")) {
        WP2_CHECK_STATUS(ReadTransform(tile_->rect.width, tile_->rect.height,
                                       COLOR_INDEXING_TRANSFORM));
      } else if (dec_->ReadBool("use_group4")) {
        // TODO(vrabaud) Merge with color indexing as it is a subcase.
        use_group4_ = true;
        WP2_CHECK_STATUS(ReadTransform(tile_->rect.width, tile_->rect.height,
                                       GROUP4));
        return WP2_STATUS_OK;
      } else {
        // Color indexing cannot be used with other transforms.
        const bool use_subtract_green = dec_->ReadBool("use_subtract_green");
        const bool use_predictor = dec_->ReadBool("use_predictor");
        const bool use_cross_color = dec_->ReadBool("use_cross_color");
        if (use_subtract_green) {
          WP2_CHECK_STATUS(ReadTransform(tile_->rect.width, tile_->rect.height,
                                         SUBTRACT_GREEN));
        }
        if (use_predictor) {
          WP2_CHECK_STATUS(ReadTransform(tile_->rect.width, tile_->rect.height,
                                         PREDICTOR_TRANSFORM));
        }
        if (use_cross_color) {
          WP2_CHECK_STATUS(ReadTransform(tile_->rect.width, tile_->rect.height,
                                         CROSS_COLOR_TRANSFORM));
        }
      }
    }
  }

  // Color cache
  hdr_.cache_config_.type = CacheType::kNone;
  hdr_.cache_config_.cache_bits = 0;
  if (dec_->ReadBool("color_cache")) {
    hdr_.cache_config_.type = CacheType::kHash;
    hdr_.cache_config_.cache_bits =
        dec_->ReadRange(1, kMaxCacheBits, "color_cache_bits");
    WP2_CHECK_STATUS(hdr_.color_cache_.Init(hdr_.cache_config_));
  }

  // Find the maximum subsampling that would result in at least 2 histograms.
  uint32_t max_sampling;
  bool more_than_one_hist = false;
  for (max_sampling = kHistogramBitsMax + 1;
       !more_than_one_hist && max_sampling-- > kHistogramBitsMin;) {
    more_than_one_hist =
        (WP2LSubSampleSize(tile_->rect.width, max_sampling) > 1 ||
         WP2LSubSampleSize(tile_->rect.height, max_sampling) > 1);
  }

  // Read the entropy image if needed.
  uint32_t num_histo;
  if (more_than_one_hist && dec_->ReadBool("write_histogram_image")) {
    WP2::ANSDebugPrefix prefix_histogram_image(dec_, "histogram_image");
    const uint32_t histogram_precision =
        dec_->ReadRange(kHistogramBitsMin, max_sampling, "histogram_bits");
    const uint32_t histogram_width =
        WP2LSubSampleSize(tile_->rect.width, histogram_precision);
    const uint32_t histogram_height =
        WP2LSubSampleSize(tile_->rect.height, histogram_precision);
    assert(histogram_width > 1 || histogram_height > 1);  // more_than_one_hist

    num_histo = dec_->ReadRange(
        2, std::min(kMaxHistogramImageSize, histogram_width * histogram_height),
        "num_histograms_m2");

    // Read the entropy image probabilities.
    LosslessSymbolsInfo symbols_info;
    symbols_info.Init(gparams_->has_alpha_, hdr_.symbols_info_.SampleFormat());
    symbols_info.SetAsLabelImage();
    symbols_info.SetMinMax(kSymbolG, /*min=*/0, /*max=*/num_histo - 1);
    WP2_CHECK_STATUS(
        ReadANSStats(histogram_width, histogram_height, symbols_info));

    // Read the entropy image.
    WP2_CHECK_STATUS(DecodeNoProcess(histogram_width, histogram_height,
                                     /*last_row=*/histogram_height,
                                     &hdr_.histogram_image_));
    WP2_CHECK_STATUS(dec_->GetStatus());

    hdr_.histogram_subsample_bits_ = histogram_precision;
  } else {
    num_histo = 1;
    hdr_.histogram_image_.clear();
  }

  WP2_CHECK_STATUS(dec_->GetStatus());

  // Read the ANS probabilities.
  LosslessSymbolsInfo symbols_info;
  WP2_CHECK_STATUS(symbols_info.CopyFrom(hdr_.symbols_info_));
  symbols_info.SetNumClusters(num_histo);
  symbols_info.SetCacheRange(GetColorCacheRange(hdr_.cache_config_));
  WP2_CHECK_STATUS(
      ReadANSStats(tile_->rect.width, tile_->rect.height, symbols_info));

  // Finish setting up the color-cache
  if (hdr_.cache_config_.type != CacheType::kNone) {
    hdr_.symbols_info_.SetCacheRange(GetColorCacheRange(hdr_.cache_config_));
  }

  // Update information about the image.
  const int num_bits = hdr_.histogram_subsample_bits_;

  hdr_.histogram_xsize_ = WP2LSubSampleSize(tile_->rect.width, num_bits);
  hdr_.histogram_mask_ = (num_bits == 0) ? ~0 : (1 << num_bits) - 1;

  HeaderVDebug();

  return WP2_STATUS_OK;
}

//------------------------------------------------------------------------------

WP2Status Decoder::DecodeImage() { return DecodeLines(tile_->rect.height); }

WP2Status Decoder::DecodeLines(uint32_t n_lines) {
  assert(last_row_ + n_lines <= tile_->rect.height);
  // Initialization.
  if (state_ != READ_DATA) {
    WP2_CHECK_STATUS(AllocateInternalBuffers(tile_->rect.width));
    state_ = READ_DATA;
  }

  WP2::Planef* const bits_per_pixel =
      bits_per_pixel_.IsEmpty() ? nullptr : &bits_per_pixel_;

  // Decode.
  WP2_CHECK_STATUS(DecodeImageData(tile_->rect.width, tile_->rect.height,
                                   last_row_ + n_lines, bits_per_pixel));

  if (VDMatch(config_, "bits-per-pixel") && last_row_ == tile_->rect.height) {
    WP2::ArgbBuffer debug_output;
    WP2_CHECK_STATUS(
        debug_output.SetView(config_.info->debug_output, tile_->rect));
    WP2_CHECK_STATUS(
        bits_per_pixel_.ToGray(&debug_output, /*num_bits=*/5, /*shift=*/false));
  }

  return WP2_STATUS_OK;
}

//------------------------------------------------------------------------------

}  // namespace WP2L
