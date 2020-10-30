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
// Intratile directional filter.
//
// Author: Yannis Guyon (yguyon@google.com)

#include "src/dec/filters/directional_filter.h"

#include "src/dec/wp2_dec_i.h"
#include "src/dsp/dsp.h"
#include "src/dsp/math.h"
#include "src/utils/utils.h"
#include "src/wp2/format_constants.h"

namespace WP2 {

//------------------------------------------------------------------------------

constexpr uint32_t DirectionalFilter::kMaxPriStr;
constexpr uint32_t DirectionalFilter::kMaxSecStr;

static uint32_t GetStrengthFromQuality(uint32_t quality_hint,
                                       uint32_t max_strength) {
  if (quality_hint == kLosslessQualityHint) return 0;
  assert(quality_hint <= kMaxLossyQualityHint);
  // Full filtering at low/medium qualities, no filtering at high qualities.
  constexpr uint32_t kThresh = kMaxLossyQualityHint - 1;  // quality ~88
  return std::min(
      DivRound(SafeSub(kThresh, quality_hint) * max_strength * 3u, kThresh),
      max_strength);
}

bool IsDirectionalFilterEnabled(const DecoderConfig& config,
                                const BitstreamFeatures& features) {
  const uint32_t strength = GetStrengthFromQuality(
      features.quality_hint,
      std::max(DirectionalFilter::kMaxPriStr, DirectionalFilter::kMaxSecStr));

  return ((config.enable_directional_filter && strength > 0) ||
          VDMatch(config, "directional-filter"));
}

//------------------------------------------------------------------------------

DirectionalFilter::DirectionalFilter(const DecoderConfig& config,
                                     const BitstreamFeatures& features,
                                     const FilterBlockMap& blocks)
    : config_(config),
      features_(features),
      blocks_(blocks),
      enabled_(IsDirectionalFilterEnabled(config, features)) {
  assert(blocks_.num_precision_bits_ >= 8u);  // To match AV1 algorithms.
}

WP2Status DirectionalFilter::Allocate() {
  if (!enabled_) return WP2_STATUS_OK;

  WP2_CHECK_STATUS(top_taps_.Resize(blocks_.tile_rect_.width, kDrctFltTapDist));
  WP2_CHECK_STATUS(bottom_backup_.Resize(top_taps_.Y.w_, top_taps_.Y.h_));

  WP2_CHECK_STATUS(left_taps_.Resize(kDrctFltTapDist, kDrctFltSize));
  WP2_CHECK_STATUS(right_backup_.Resize(left_taps_.Y.w_, left_taps_.Y.h_));

  WP2MathInit();
  DrctFilterInit();
  return WP2_STATUS_OK;
}

//------------------------------------------------------------------------------

// Reduces the filter primary strength for low 'variance' values (when there is
// uncertainty about the detected direction), only for the luma channel.
static uint32_t AdjustPriStr(Channel channel, uint32_t pri_str,
                             uint32_t variance) {
  if (channel == kYChannel) {
    if (variance == 0u) return 0u;
    const uint32_t variance_strength =
        std::min((uint32_t)WP2Log2Floor(variance >> 6), 12u);
    return (pri_str * (4u + variance_strength) + 8u) >> 4;
  }
  return pri_str;
}

// The difference between the original value of the pixel to filter and the taps
// will be right-shifted by this damping. Higher values mean weaker filtering.
static uint32_t GetDamping(Channel channel, uint32_t precision_shift) {
  return 3u + precision_shift - ((channel == kYChannel) ? 0u : 1u);
}

// Converts available stats to an AV1 CDEF primary or secondary strength.
static uint32_t GetCdefStr(uint32_t quality_hint, uint32_t min_bpp,
                           uint32_t res_den, uint32_t max_cdef_str) {
  assert(min_bpp <= 255 && res_den <= 255);
  min_bpp = (min_bpp * min_bpp) / 255;
  res_den = (res_den * res_den) / 255;
  // More details (bpp, res) need less filtering.
  const uint32_t strength = GetStrengthFromQuality(quality_hint, max_cdef_str);
  return ((255 - min_bpp) * (255 - res_den) * strength + 127 * 127) /
         (255 * 255);
}

uint32_t DirectionalFilter::Smooth(uint32_t num_rows) {
  if (!enabled_) return num_rows;

  YUVPlane* const pixels = blocks_.pixels_;
  assert(num_rows <= blocks_.tile_rect_.height);
  assert(num_filtered_rows_ <= num_rows);

  uint32_t num_filterable_rows;
  if (num_rows == blocks_.tile_rect_.height) {
    // All pixels are available: filter everything left.
    num_filterable_rows = num_rows;
  } else {
    // Wait for bottom taps to be available.
    num_filterable_rows = SafeSub(num_rows, kDrctFltTapDist);
  }

  if (VDMatch(config_, "directional-filter")) SavePixelsForVDebug(num_rows);
  if (num_filterable_rows <= num_filtered_rows_) return num_filtered_rows_;

  const uint32_t block_step = blocks_.max_num_blocks_x_;
  const uint32_t precision_shift = (blocks_.num_precision_bits_ - 8u);

  // Temporary padded buffer for CDEF computation.
  // 'kDrctFltSize x kDrctFltSize' pixels + margins of 'kDrctFltTapDist' pixels.
  static constexpr int32_t kTmpStep =
      kDrctFltTapDist + kDrctFltSize + kDrctFltTapDist;
  int16_t tmp_buf[kTmpStep * kTmpStep];
  int16_t* const tmp = tmp_buf + kDrctFltTapDist * kTmpStep + kDrctFltTapDist;

  uint32_t y = num_filtered_rows_;
  while (y < blocks_.tile_rect_.height) {
    const uint32_t height =
        std::min(kDrctFltSize, blocks_.tile_rect_.height - y);
    if (y + height > num_filterable_rows) break;

    // Backup bottom pixels for the top taps of the next row of blocks, if any.
    const bool backup_bottom = (y + height < blocks_.tile_rect_.height);
    if (backup_bottom) {
      assert(height >= bottom_backup_.Y.h_);
      YUVPlane bottom;
      WP2_ASSERT_STATUS(
          bottom.SetView(*pixels, {0, y + height - bottom_backup_.Y.h_,
                                   bottom_backup_.Y.w_, bottom_backup_.Y.h_}));
      bottom.A.Clear();  // Alpha is not filtered here.
      WP2_ASSERT_STATUS(
          bottom_backup_.Copy(bottom, /*resize_if_needed=*/false));
    }

    for (uint32_t x = 0; x < blocks_.tile_rect_.width; x += kDrctFltSize) {
      const uint32_t width =
          std::min(kDrctFltSize, blocks_.tile_rect_.width - x);

      // TODO(yguyon): Filter partial blocks (bottom- or right-most in a tile)
      if (height != kDrctFltSize || width != kDrctFltSize) break;

      const uint32_t block_x = x / kMinBlockSizePix;
      const uint32_t block_y = y / kMinBlockSizePix;
      const uint32_t num_blocks_x = SizeBlocks(width);
      const uint32_t num_blocks_y = SizeBlocks(height);

      // Backup right pixels for the left taps of the next block.
      const bool backup_right = (x + width < blocks_.tile_rect_.width);
      if (backup_right) {
        assert(width >= right_backup_.Y.w_);
        YUVPlane right;
        WP2_ASSERT_STATUS(
            right.SetView(*pixels, {x + width - right_backup_.Y.w_, y,
                                    right_backup_.Y.w_, right_backup_.Y.h_}));
        right.A.Clear();  // Alpha is not filtered here.
        WP2_ASSERT_STATUS(
            right_backup_.Copy(right, /*resize_if_needed=*/false));
      }

      // Get average num-bits-per-pixel and residual density.
      uint32_t min_bpp = 0, res_den[3]{0};
      const FilterBlockMap::BlockFeatures* block_row =
          &blocks_.features_[block_y * block_step + block_x];
      for (uint32_t offset_y = 0; offset_y < num_blocks_y; ++offset_y) {
        for (uint32_t offset_x = 0; offset_x < num_blocks_x; ++offset_x) {
          min_bpp += block_row[offset_x].min_bpp;
          for (Channel c : {kYChannel, kUChannel, kVChannel}) {
            res_den[c] += block_row[offset_x].res_den[c];
          }
        }
        block_row += block_step;
      }
      const uint32_t num_blocks = num_blocks_x * num_blocks_y;
      min_bpp = (min_bpp + num_blocks / 2) / num_blocks;
      for (Channel c : {kYChannel, kUChannel, kVChannel}) {
        res_den[c] = (res_den[c] + num_blocks / 2) / num_blocks;
      }

      // Get overall direction and its variance.
      uint32_t direction, variance;
      CdefDirection8x8(&pixels->Y.At(x, y), (int32_t)pixels->Y.Step(),
                       blocks_.num_precision_bits_, &direction, &variance);

      // Margins for taps. Top-left is the origin.
      const uint32_t n_left = x;
      const uint32_t n_right = blocks_.tile_rect_.width - (x + width);
      const uint32_t n_top = y;
      const uint32_t n_bottom = blocks_.tile_rect_.height - (y + height);

      const bool margin_left = (n_left >= left_taps_.Y.w_);
      const bool margin_top = (n_top >= top_taps_.Y.h_);
      const uint32_t backup_x = margin_left ? (x - left_taps_.Y.w_) : x;

      for (Channel c : {kYChannel, kUChannel, kVChannel}) {
        Plane16* const pixels_plane = &pixels->GetChannel(c);
        const Plane16& left_plane = left_taps_.GetChannel(c);
        const Plane16& top_plane = top_taps_.GetChannel(c);

        const int16_t* const left =
            margin_left ? &left_plane.At(0, 0) : nullptr;
        const int16_t* const top =
            margin_top ? &top_plane.At(backup_x, 0) : nullptr;

        const uint32_t pri_str =
            GetCdefStr(features_.quality_hint, min_bpp, res_den[c], kMaxPriStr);
        const uint32_t sec_str =
            GetCdefStr(features_.quality_hint, min_bpp, res_den[c], kMaxSecStr);
        const uint32_t shifted_pri_str = pri_str << precision_shift;
        const uint32_t shifted_sec_str = sec_str << precision_shift;

        // Copy pixels from source plane and backups to 'tmp'.
        CdefPad(&pixels_plane->At(x, y), (int32_t)pixels_plane->Step(),
                left, (int32_t)left_plane.Step(),
                top, (int32_t)top_plane.Step(),
                (int32_t)width, (int32_t)height, (int32_t)n_left,
                (int32_t)n_right, (int32_t)n_top, (int32_t)n_bottom,
                tmp, kTmpStep);

        CdefFiltering(tmp, kTmpStep, blocks_.num_precision_bits_, width, height,
                      AdjustPriStr(c, shifted_pri_str, variance),
                      shifted_sec_str, GetDamping(c, precision_shift),
                      direction, &pixels_plane->At(x, y),
                      (int32_t)pixels_plane->Step());
      }

      if (VDMatch(config_, "directional-filter")) {
        RegisterPixelsForVDebug(
            x, x + width - 1u, y, y + height - 1u,
            GetCdefStr(features_.quality_hint, min_bpp, res_den[0], kMaxPriStr),
            GetCdefStr(features_.quality_hint, min_bpp, res_den[0], kMaxSecStr),
            direction, variance);
      }

      // One's right column is another one's left taps.
      if (backup_right) {
        using std::swap;
        swap(left_taps_, right_backup_);
      }
    }

    if (backup_bottom) {
      using std::swap;
      swap(top_taps_, bottom_backup_);
    }

    y += height;
  }
  num_filtered_rows_ = y;

  if (VDMatch(config_, "directional-filter")) ApplyVDebug(num_rows);
  return num_filtered_rows_;
}

//------------------------------------------------------------------------------

}  // namespace WP2
