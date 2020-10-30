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
// Intertile filter.
//
// Author: Yannis Guyon (yguyon@google.com)

#include "src/dec/filters/intertile_filter.h"

#include "src/dec/wp2_dec_i.h"
#include "src/dec/filters/deblocking_filter.h"
#include "src/dsp/dsp.h"
#include "src/dsp/math.h"
#include "src/utils/utils.h"

namespace WP2 {

typedef FilterBlockMap::BlockFeatures BlockFeatures;

//------------------------------------------------------------------------------

bool IsIntertileFilterEnabled(const DecoderConfig& config,
                              const BitstreamFeatures& features,
                              const GlobalParams& gparams,
                              const TilesLayout& tiles_layout) {
  return (config.enable_deblocking_filter &&
          gparams.type_ != GlobalParams::GP_LOSSLESS &&
          DeblockingFilter::GetStrengthFromQuality(features.quality_hint) > 0 &&
          (tiles_layout.num_tiles_x > 1 || tiles_layout.num_tiles_y > 1)) ||
         VDMatch(config, "intertile-filter");
  // TODO(yguyon): There might be tiles that are lossy but output RGB (neural)
  // TODO(yguyon): There might be no pair of adjacent lossy tiles
}

//------------------------------------------------------------------------------

void IntertileFilter::Init(const DecoderConfig& config,
                           const BitstreamFeatures& features,
                           const GlobalParams& gparams,
                           const TilesLayout& tiles_layout,
                           YUVPlane* const canvas) {
  assert(!canvas->IsEmpty());
  assert(tiles_layout.num_tiles_x >= 1 && tiles_layout.num_tiles_y >= 1);
  const Tile& last_tile = tiles_layout.tiles.back();
  const uint32_t width = last_tile.rect.x + last_tile.rect.width;
  const uint32_t height = last_tile.rect.y + last_tile.rect.height;
  // 'tiles_layout' gives access to the 'canvas' in each tile but it's
  // convenient to have the whole buffer. Make sure it's the same padded size:
  assert(Pad(width, kPredWidth) == canvas->Y.w_);
  assert(Pad(height, kPredWidth) == canvas->Y.h_);

  Clear();
  config_ = &config;
  gparams_ = &gparams;
  tiles_layout_ = &tiles_layout;
  // Non-padded view.
  WP2_ASSERT_STATUS(canvas_.SetView(*canvas, {0, 0, width, height}));
  quality_hint_ = features.quality_hint;
  enabled_ = IsIntertileFilterEnabled(config, features, gparams, tiles_layout);
  if (enabled_) DblkFilterInit();
}

void IntertileFilter::Clear() {
  enabled_ = false;
  tiles_layout_ = nullptr;
  canvas_.Clear();
  quality_hint_ = 0;
  num_horizontally_deblocked_rows_ = 0;
  num_vertically_deblocked_rows_ = 0;

  num_debug_rows_ = 0;
}

//------------------------------------------------------------------------------

void IntertileFilter::Deblock(uint32_t num_rows) {
  assert(num_horizontally_deblocked_rows_ <= num_rows);
  assert(num_vertically_deblocked_rows_ <= num_horizontally_deblocked_rows_);

  if (!enabled_) {
    num_horizontally_deblocked_rows_ = num_rows;
    num_vertically_deblocked_rows_ = num_horizontally_deblocked_rows_;
  } else {
    assert(num_rows <= canvas_.Y.h_);
    if (VDMatch(*config_, "intertile-filter")) SavePixelsForVDebug(num_rows);

    DeblockHorizontally(num_rows);
    DeblockVertically(num_horizontally_deblocked_rows_);

    if (VDMatch(*config_, "intertile-filter")) ApplyVDebug(num_rows);
  }
}

//------------------------------------------------------------------------------
// For the given tiles A, B, C and D:     AAA | BB
//                                        AAA | BB
//                                        ----+---
//                                        CCC | DD

// Deblocks vertical edges (A | B, C | D) horizontally.
void IntertileFilter::DeblockHorizontally(uint32_t num_rows) {
  const Rectangle frame_rect = {0, 0, canvas_.GetWidth(), canvas_.GetHeight()};
  const uint32_t num_bits = gparams_->transf_.GetYUVPrecisionBits() + 1;
  const uint32_t yuv_min = gparams_->transf_.GetYUVMin();
  const uint32_t yuv_max = gparams_->transf_.GetYUVMax();

  uint32_t y = num_horizontally_deblocked_rows_;
  while (y < num_rows) {
    const uint32_t tile_width = tiles_layout_->tile_width;
    const uint32_t tile_height = tiles_layout_->tile_height;
    const uint32_t y_in_tile = y % tile_height;
    const uint32_t height = std::min(num_rows - (y - y_in_tile), tile_height);

    for (uint32_t x = tile_width; x < canvas_.Y.w_; x += tile_width) {
      const Tile& p_tile = tiles_layout_->GetTileAt(x - tile_width, y);
      const Tile& q_tile = tiles_layout_->GetTileAt(x, y);
      if (!p_tile.output_is_yuv) continue;  // This filter only applies
      if (!q_tile.output_is_yuv) continue;  // to two lossy tiles.

      const uint32_t p_features_step = p_tile.block_map.max_num_blocks_x_;
      const uint32_t q_features_step = q_tile.block_map.max_num_blocks_x_;
      const BlockFeatures* p_features =
          p_tile.block_map.GetRow(y_in_tile / kMinBlockSizePix) +
          (p_features_step - 1);
      const BlockFeatures* q_features =
          q_tile.block_map.GetRow(y_in_tile / kMinBlockSizePix);

      // For each pixel along the edge = for each line accross the edge:
      uint32_t line = y_in_tile;
      while (line < height) {
        // Number of lines accross the edge sharing the same block features.
        const uint32_t num_same_lines =
            std::min(height - line, kMinBlockSizePix - line % kMinBlockSizePix);

        const Rectangle edge_rect = {x, q_tile.rect.y + line, 0,
                                     num_same_lines};
        DeblockingFilter::DeblockEdge(*config_, /*intertile=*/true, frame_rect,
                                      canvas_.HasAlpha(), quality_hint_,
                                      num_bits, yuv_min, yuv_max, *p_features,
                                      *q_features, edge_rect, &canvas_);

        line += num_same_lines;
        p_features += p_features_step;
        q_features += q_features_step;
      }
    }

    y += (height - y_in_tile);
  }
  num_horizontally_deblocked_rows_ = y;
}

//                             A    B
// Deblocks horizontal edges (---, ---) vertically.
//                             C    D
void IntertileFilter::DeblockVertically(uint32_t num_rows) {
  const uint32_t max_half_n = kMaxBlockSizePix / 2;
  const uint32_t tile_width = tiles_layout_->tile_width;
  const uint32_t tile_height = tiles_layout_->tile_height;
  const Rectangle frame_rect = {0, 0, canvas_.GetWidth(), canvas_.GetHeight()};
  const uint32_t num_bits = gparams_->transf_.GetYUVPrecisionBits() + 1;
  const uint32_t yuv_min = gparams_->transf_.GetYUVMin();
  const uint32_t yuv_max = gparams_->transf_.GetYUVMax();

  uint32_t y = num_vertically_deblocked_rows_;
  while (y < num_rows) {
    uint32_t height;
    const uint32_t tile_y = y / tile_height;
    const uint32_t y_in_tile = y % tile_height;
    const bool is_first_tile = (tile_y == 0);
    const bool is_last_tile = (tile_y + 1 == tiles_layout_->num_tiles_y);

    // Check whether y is near the edges between two consecutive rows of tiles.
    if ((is_first_tile || y_in_tile >= max_half_n) &&
        (is_last_tile || y_in_tile < tile_height - max_half_n)) {
      // Within a tile. These rows of pixels are not vertically filtered.
      if (is_last_tile) {
        height = num_rows - y;
      } else {
        height = std::min(num_rows - y, tile_height - max_half_n - y_in_tile);
      }
    } else {
      assert(tile_height >= kMaxBlockSizePix);
      // Near the junction of two rows of tiles (at 'y + max_half_n' exactly).
      height = std::min(canvas_.Y.h_ - y, 2 * max_half_n);
      // Wait for all pixels on each side of the horizontal edge.
      if (y + height > num_rows) break;
      const uint32_t edge_y = y + max_half_n;
      assert(edge_y % tile_height == 0);  // Edge expected here.

      for (uint32_t x = 0; x < canvas_.Y.w_; x += tile_width) {
        const Tile& p_tile = tiles_layout_->GetTileAt(x, y);
        const Tile& q_tile = tiles_layout_->GetTileAt(x, y + max_half_n);
        if (!p_tile.output_is_yuv) continue;  // This filter only applies
        if (!q_tile.output_is_yuv) continue;  // to two lossy tiles.

        const uint32_t features_step = 1;
        const BlockFeatures* p_features =
            p_tile.block_map.GetRow(p_tile.block_map.max_num_blocks_y_ - 1);
        const BlockFeatures* q_features = q_tile.block_map.GetRow(0);

        const uint32_t width = std::min(tile_width, canvas_.Y.w_ - x);

        // For each pixel along the edge = for each line accross the edge:
        uint32_t line = 0;
        while (line < width) {
          // Number of lines accross the edge sharing the same block features.
          const uint32_t num_same_lines = std::min(
              width - line, kMinBlockSizePix - line % kMinBlockSizePix);

          const Rectangle edge_rect = {x + line, edge_y, num_same_lines, 0};
          DeblockingFilter::DeblockEdge(
              *config_, /*intertile=*/true, frame_rect, canvas_.HasAlpha(),
              quality_hint_, num_bits, yuv_min, yuv_max, *p_features,
              *q_features, edge_rect, &canvas_);

          line += num_same_lines;
          p_features += features_step;
          q_features += features_step;
        }
      }
    }

    y += height;
  }
  num_vertically_deblocked_rows_ = y;
}

//------------------------------------------------------------------------------

}  // namespace WP2
