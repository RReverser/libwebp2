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
// Intratile restoration filter.
//
// Author: Yannis Guyon (yguyon@google.com)

#include "src/dec/filters/restoration_filter.h"

#include "src/dec/wp2_dec_i.h"
#include "src/dsp/dsp.h"
#include "src/dsp/math.h"
#include "src/utils/utils.h"

namespace WP2 {

//------------------------------------------------------------------------------

bool IsRestorationFilterEnabled(const DecoderConfig& config) {
  return (config.enable_restoration_filter ||
          VDMatch(config, "restoration-filter"));
}

//------------------------------------------------------------------------------

RestorationFilter::RestorationFilter(const DecoderConfig& config,
                                     const FilterBlockMap& blocks,
                                     const RstrFltParams& params)
    : config_(config),
      blocks_(blocks),
      enabled_(IsRestorationFilterEnabled(config)),
      params_(params) {
  assert(blocks_.num_precision_bits_ >= 8u);  // To match AV1 algorithms.
}

WP2Status RestorationFilter::Allocate() {
  if (!enabled_) return WP2_STATUS_OK;

  WP2_CHECK_STATUS(top_taps_.Resize(blocks_.tile_rect_.width, kWieFltTapDist));
  WP2_CHECK_STATUS(bottom_backup_.Resize(top_taps_.Y.w_, top_taps_.Y.h_));
  WP2_CHECK_ALLOC_OK(filter_strength_map_.resize(
      std::min(kWieFltWidth, blocks_.tile_rect_.width) *
      std::min(kWieFltHeight, blocks_.tile_rect_.height)));

  WP2MathInit();
  WienerFilterInit();
  return WP2_STATUS_OK;
}

//------------------------------------------------------------------------------

uint32_t RestorationFilter::Enhances(uint32_t num_rows) {
  if (!enabled_) return num_rows;

  YUVPlane* const pixels = blocks_.pixels_;
  assert(num_rows <= blocks_.tile_rect_.height);
  assert(num_filtered_rows_ <= num_rows);

  const bool last_rows = (num_rows == blocks_.tile_rect_.height);
  uint32_t num_filterable_rows;
  if (last_rows) {
    // All pixels are available: filter everything left.
    num_filterable_rows = num_rows;
  } else {
    // Wait for bottom taps to be available.
    num_filterable_rows = SafeSub(num_rows, kWieFltTapDist);
  }

  if (VDMatch(config_, "restoration-filter")) SavePixelsForVDebug(num_rows);
  if (num_filterable_rows <= num_filtered_rows_) return num_filtered_rows_;

  const size_t block_step = blocks_.max_num_blocks_x_;

  uint32_t y = num_filtered_rows_;
  while (y < blocks_.tile_rect_.height) {
    const uint32_t height =
        std::min(kWieFltHeight, blocks_.tile_rect_.height - y);
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

    // Whole rows can be processed in batches.
    assert(blocks_.tile_rect_.width <= kWieFltWidth);
    const uint32_t width = blocks_.tile_rect_.width;

    // Copy filter strength to a layout usable by dsp functions.
    bool apply_filtering = false;
    for (uint32_t offset_y = 0; offset_y < height; ++offset_y) {
      const FilterBlockMap::BlockFeatures* block_row =
          &blocks_.features_[((y + offset_y) / kMinBlockSizePix) * block_step];
      for (uint32_t offset_x = 0; offset_x < blocks_.tile_rect_.width;
           ++offset_x) {
        // TODO(yguyon): Adjust, take res_den into account
        const uint8_t strength =
            63 * (255 - block_row[offset_x / kMinBlockSizePix].min_bpp) / 255;
        filter_strength_map_[offset_y * blocks_.tile_rect_.width + offset_x] =
            strength;
        apply_filtering |= (strength > 0);
      }
    }

    if (apply_filtering) {
      // Margins for taps. Top-left is the origin.
      const uint32_t n_top = std::min(kWieFltTapDist, y);
      const uint32_t n_bottom =
          std::min(kWieFltTapDist, blocks_.tile_rect_.height - (y + height));

      const uint32_t area = y / kWieFltHeight;
      for (Channel c : {kYChannel, kUChannel, kVChannel}) {
        // Construct the tap weights from those extracted from the bitstream.
        int32_t tap_weights_h[kWieFltNumTaps], tap_weights_v[kWieFltNumTaps];
        bool identity = true;
        for (uint32_t i = 0; i < kWieFltTapDist; ++i) {
          tap_weights_h[i] = params_.half_tap_weights[c][area][0][i];
          tap_weights_v[i] = params_.half_tap_weights[c][area][1][i];
          identity &= (tap_weights_h[i] == 0 && tap_weights_v[i] == 0);
        }

        // If the filter is made of 0s, no need to apply it.
        if (!identity) {
          Plane16* const pixels_plane = &pixels->GetChannel(c);
          const Plane16& top_plane = top_taps_.GetChannel(c);

          const int16_t* const unfiltered_top =
              (n_top > 0) ? top_plane.Row(kWieFltTapDist - n_top) : nullptr;
          const int16_t* const unfiltered_bottom =
              (n_bottom > 0) ? pixels_plane->Row(y + height) : nullptr;

          WienerHalfToFullWgts(tap_weights_h, tap_weights_h);
          WienerHalfToFullWgts(tap_weights_v, tap_weights_v);

          WP2_CHECK_STATUS(WienerFilter(
              width, height,
              /*left=*/nullptr, /*left_step=*/0,
              /*right=*/nullptr, /*right_step=*/0,
              /*top=*/unfiltered_top, top_plane.Step(),
              /*bottom=*/unfiltered_bottom, pixels_plane->Step(),
              /*n_left=*/0, /*n_right=*/0, n_top, n_bottom, tap_weights_h,
              tap_weights_v, filter_strength_map_.data(),
              blocks_.tile_rect_.width, blocks_.num_precision_bits_,
              pixels_plane->Step(), pixels_plane->Row(y)));
        }
      }
    }
    if (VDMatch(config_, "restoration-filter")) {
      RegisterPixelsForVDebug(0, width - 1u, y, y + height - 1u);
    }

    if (backup_bottom) {
      using std::swap;
      swap(top_taps_, bottom_backup_);
    }

    y += height;
  }
  num_filtered_rows_ = y;

  if (VDMatch(config_, "restoration-filter")) ApplyVDebug(num_rows);
  return num_filtered_rows_;
}

//------------------------------------------------------------------------------

WP2Status ReadRestorationFilterParams(ANSDec* const dec,
                                      RstrFltParams* const params) {
  for (Channel channel : {kYChannel, kUChannel, kVChannel}) {
    for (uint32_t area = 0; area < params->num_areas; ++area) {
      for (uint32_t h_or_v : {0, 1}) {
        int8_t* const half = params->half_tap_weights[channel][area][h_or_v];
        for (uint32_t i = 0; i < kWieFltTapDist; ++i) {
          // TODO(yguyon): Read it from bitstream when computed by encoder.
          half[i] = (i == 0) ? 3 : (i == 1) ? -7 : 15;  // Mock values.
          // const uint32_t num_bits = kWieNumBitsPerTapWgt[i];
          // const int32_t min = GetMinWienerTapWeight(i);
          // half[i] = dec->ReadUValue(num_bits, "wiener_weights") + min;
        }
        assert(WienerVerifyHalfTapWeights(half));
      }
    }
  }
  return dec->GetStatus();
}

//------------------------------------------------------------------------------

}  // namespace WP2
