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
// Deblocking filter.
//
// Author: Yannis Guyon (yguyon@google.com)

#include "src/dec/filters/deblocking_filter.h"

#include "src/dec/wp2_dec_i.h"
#include "src/dsp/dsp.h"
#include "src/dsp/math.h"
#include "src/utils/utils.h"
#include "src/wp2/format_constants.h"

namespace WP2 {

//------------------------------------------------------------------------------

uint32_t DeblockingFilter::GetStrengthFromQuality(uint32_t quality_hint) {
  if (quality_hint == kLosslessQualityHint) return 0;
  assert(quality_hint <= kMaxLossyQualityHint);
  // Full filtering below mid quality, no filtering at top quality.
  return std::min(
      DivRound((kMaxLossyQualityHint - quality_hint) * kDblkMaxStrength * 2,
               kMaxLossyQualityHint),
      kDblkMaxStrength);
}

uint32_t DeblockingFilter::GetStrengthFromBpp(uint32_t num_bits_per_pixel) {
  assert(num_bits_per_pixel <= 255);
  // Full filtering at 0bpp, no filtering at 2bpp (=255).
  return std::min((255 - num_bits_per_pixel) * kDblkMaxStrength / 255,
                  kDblkMaxStrength);
}

uint32_t DeblockingFilter::GetSharpnessFromQuality(uint32_t quality_hint) {
  if (quality_hint == kLosslessQualityHint) return kDblkMaxSharpness;
  return DivRound(quality_hint * kDblkMaxSharpness, kMaxLossyQualityHint);
}

uint32_t DeblockingFilter::GetSharpnessFromResidualDensity(uint32_t res_den) {
  assert(res_den <= 255);
  // Heavy filtering at 10% (=25) of non-zero coeffs, slight filtering above 30%
  // (=75).
  return std::min(SafeSub(res_den, 25u) * kDblkMaxSharpness / (75 - 25),
                  kDblkMaxSharpness);
}

//------------------------------------------------------------------------------

bool IsDeblockingFilterEnabled(const DecoderConfig& config,
                               const BitstreamFeatures& features) {
  const uint32_t strength =
      DeblockingFilter::GetStrengthFromQuality(features.quality_hint);

  return (config.enable_deblocking_filter && strength > 0) ||
         VDMatch(config, "deblocking-filter");
}

//------------------------------------------------------------------------------

DeblockingFilter::DeblockingFilter(const DecoderConfig& config,
                                   const BitstreamFeatures& features,
                                   const FilterBlockMap& blocks)
    : config_(config),
      features_(features),
      blocks_(blocks),
      enabled_(IsDeblockingFilterEnabled(config, features)) {
  assert(blocks_.num_precision_bits_ >= 8u);  // To match AV1 algorithms.
}

WP2Status DeblockingFilter::Allocate() {
  if (!enabled_) return WP2_STATUS_OK;

  // 'x_to_next_up_block_y_' is filled with zeros (by ctor).
  WP2_CHECK_ALLOC_OK(x_to_next_up_block_y_.resize(blocks_.max_num_blocks_x_));
  // "Block position should fit in x_to_next_up_block_y_
  static_assert(kMaxTileSize / kMinBlockSizePix <
                (1u << (sizeof(x_to_next_up_block_y_[0]) * 8)),
                "Incorrect kMaxTileSize value");
  // No vertical deblocking will happen before at least 1 min block and a half.
  min_num_rows_to_vdblk_ = std::min(kMinBlockSizePix + kMinBlockSizePix / 2,
                                    blocks_.tile_rect_.height);

  DblkFilterInit();
  return WP2_STATUS_OK;
}

//------------------------------------------------------------------------------

uint32_t DeblockingFilter::Deblock(uint32_t num_rows) {
  if (!enabled_) return num_rows;

  assert(num_rows <= blocks_.tile_rect_.height);
  assert(num_horizontally_deblocked_rows_ <= num_rows);
  assert(num_vertically_deblocked_rows_ <= num_horizontally_deblocked_rows_);

  uint32_t num_deblockable_rows = num_rows;

  if (VDMatch(config_, "deblocking-filter")) SavePixelsForVDebug(num_rows);

  // All vertical edges up to 'num_deblockable_rows' can be horizontally
  // deblocked now.
  const uint32_t from_y = num_horizontally_deblocked_rows_;
  if (from_y >= num_deblockable_rows) return num_vertically_deblocked_rows_;
  num_horizontally_deblocked_rows_ =
      DeblockHorizontally(from_y, /*to_y=*/num_deblockable_rows - 1u);

  // Horizontal edges before row index 'num_horizontally_deblocked_rows_ -
  // half_num_filtered_pixels' can be vertically deblocked now; vertical
  // deblocking must be applied only on pixels that are already horizontally
  // deblocked, for gradient consistency.
  num_deblockable_rows = num_horizontally_deblocked_rows_;
  if (num_deblockable_rows < min_num_rows_to_vdblk_) {
    // It's unsure whether at least one edge per column is ready to be deblocked
    // so wait before maybe uselessly browsing 'max_num_blocks_x_'.
    return num_vertically_deblocked_rows_;
  }

  num_vertically_deblocked_rows_ =
      DeblockVertically(/*to_y=*/num_deblockable_rows - 1u);
  min_num_rows_to_vdblk_ =
      std::min(num_vertically_deblocked_rows_ + blocks_.max_registered_height_,
               blocks_.tile_rect_.height);

  if (VDMatch(config_, "deblocking-filter")) ApplyVDebug(num_rows);
  return num_vertically_deblocked_rows_;
}

//------------------------------------------------------------------------------

uint32_t DeblockingFilter::DeblockHorizontally(uint32_t from_y, uint32_t to_y) {
  const uint32_t block_step = blocks_.max_num_blocks_x_;

  uint32_t y = from_y;
  uint32_t block_y = y / kMinBlockSizePix;
  const FilterBlockMap::BlockFeatures* row_block_features =
      &blocks_.features_[block_y * block_step];
  while (y <= to_y) {
    // At most 'kMinBlockSizePix' rows can be deblocked with the same pattern.
    uint32_t num_pixel_rows =
        std::min(kMinBlockSizePix - y % kMinBlockSizePix, to_y - y + 1u);

    // In 'kMinBlockSizePix' units.
    uint32_t left_block_x = 0;
    uint32_t left_block_width =
        BlockWidth[row_block_features[left_block_x].size];
    uint32_t right_block_x = left_block_width;

    // All vertical edges between adjacent pairs of blocks belonging to these
    // 'num_pixel_rows' are parsed, from the left.
    while (right_block_x < blocks_.max_num_blocks_x_) {
      // The X coordinate of the left-most pixel from right block.
      const uint32_t q0_x = right_block_x * kMinBlockSizePix;
      assert(q0_x < blocks_.tile_rect_.width);

      DeblockEdge(config_, /*intertile=*/false, blocks_.tile_rect_,
                  blocks_.pixels_->HasAlpha(), features_.quality_hint,
                  blocks_.num_precision_bits_, blocks_.yuv_min_,
                  blocks_.yuv_max_, row_block_features[left_block_x],
                  row_block_features[right_block_x],
                  /*edge_rect=*/{q0_x, y, 0, num_pixel_rows}, blocks_.pixels_);

      // Advance to next pair of adjacent blocks.
      left_block_x = right_block_x;
      right_block_x += BlockWidth[row_block_features[right_block_x].size];
    }
    // Advance to next batch of at most 'kMinBlockSizePix' pixel rows.
    y += num_pixel_rows;
    ++block_y;
    row_block_features += block_step;
  }
  return y;
}

//------------------------------------------------------------------------------

uint32_t DeblockingFilter::DeblockVertically(uint32_t to_y) {
  // Includes partial blocks.
  const uint32_t to_block_y = to_y / kMinBlockSizePix;
  const uint32_t block_step = blocks_.max_num_blocks_x_;

  // For each column, wide of kMinBlockSizePix pixels.
  for (uint32_t block_x = 0; block_x < blocks_.max_num_blocks_x_; ++block_x) {
    const uint32_t x = block_x * kMinBlockSizePix;  // In pixel units.
    const uint32_t num_pixel_cols =
        std::min(kMinBlockSizePix, blocks_.tile_rect_.width - x);

    // In 'kMinBlockSizePix' units.
    const uint32_t first_block_y = x_to_next_up_block_y_[block_x];
    const FilterBlockMap::BlockFeatures* up_block_features =
        &blocks_.features_[first_block_y * block_step + block_x];
    uint32_t up_block_height = BlockHeight[up_block_features->size];
    uint32_t down_block_y = first_block_y + up_block_height;

    // Deblock as many edges as available starting from the top.
    while (down_block_y <= to_block_y) {
      const FilterBlockMap::BlockFeatures* down_block_features =
          &blocks_.features_[down_block_y * block_step + block_x];

      // In 'kMinBlockSizePix' units.
      const uint32_t down_block_height = BlockHeight[down_block_features->size];
      // The Y coordinate of the up-most pixel from the block below the edge.
      const uint32_t q0_y = down_block_y * kMinBlockSizePix;

      // Up block is never cropped.
      const uint32_t up_block_half_pixel_height =
          up_block_height * (kMinBlockSizePix / 2);
      // Down block might be cropped if it's the last one of the column.
      const uint32_t down_block_half_pixel_height =
          std::min(blocks_.tile_rect_.height - q0_y,
                   down_block_height * (kMinBlockSizePix / 2));

      // Only symmetrical filtering.
      const uint32_t min_half =
          std::min(up_block_half_pixel_height, down_block_half_pixel_height);
      const uint32_t min_num_deblockable_rows = q0_y + min_half;
      if (min_num_deblockable_rows > to_y + 1u) {
        break;  // Not enough horizontally deblocked pixel rows to continue.
      }

      DeblockEdge(config_, /*intertile=*/false, blocks_.tile_rect_,
                  blocks_.pixels_->HasAlpha(), features_.quality_hint,
                  blocks_.num_precision_bits_, blocks_.yuv_min_,
                  blocks_.yuv_max_, *up_block_features, *down_block_features,
                  /*edge_rect=*/{x, q0_y, num_pixel_cols, 0}, blocks_.pixels_);

      // Up block of next edge is down block of current edge.
      x_to_next_up_block_y_[block_x] = (uint16_t)down_block_y;
      // Advance to next pair of adjacent blocks.
      up_block_features = down_block_features;
      up_block_height = down_block_height;
      down_block_y += down_block_height;
    }
  }
  // Decoded all lines, all edges should have been deblocked.
  if (to_y + 1u == blocks_.tile_rect_.height) return blocks_.tile_rect_.height;
  // Return the position of the up-most not yet deblocked edge.
  const uint32_t min_edge_y = *std::min_element(x_to_next_up_block_y_.begin(),
                                                x_to_next_up_block_y_.end());
  return std::min(min_edge_y * kMinBlockSizePix, blocks_.tile_rect_.height);
}

//------------------------------------------------------------------------------

void DeblockingFilter::DeblockEdge(const DecoderConfig& config, bool intertile,
                                   const Rectangle& tile_rect, bool has_alpha,
                                   uint32_t quality_hint, uint32_t yuv_num_bits,
                                   int32_t yuv_min, int32_t yuv_max,
                                   const FilterBlockMap::BlockFeatures& p_block,
                                   const FilterBlockMap::BlockFeatures& q_block,
                                   const Rectangle& edge_rect,
                                   YUVPlane* const pixels) {
  assert((edge_rect.width == 0) ^ (edge_rect.height == 0));
  const bool vertical_edge = (edge_rect.width == 0);

  // The strength is determined by the block using the fewest bits per px.
  const uint32_t yuv_strength = std::min(
      GetStrengthFromQuality(quality_hint),
      GetStrengthFromBpp(std::min(p_block.min_bpp, q_block.min_bpp)));
  // TODO(maryla): this should be computed based on alpha_quality and
  //               alpha bpp instead!
  const uint32_t alpha_strength = yuv_strength;
  if ((!has_alpha || alpha_strength == 0) && yuv_strength == 0) {
    return;  // Early no-op exit.
  }

  // Left/up block is never cropped.
  const uint32_t p_block_half =
      (vertical_edge ? BlockWidth[p_block.size] : BlockHeight[p_block.size]) *
      (kMinBlockSizePix / 2);
  // Right/down block might be cropped if it's the last one of the row/column.
  const uint32_t q_block_half =
      vertical_edge
          ? std::min(tile_rect.width - edge_rect.x,
                     BlockWidth[q_block.size] * (kMinBlockSizePix / 2))
          : std::min(tile_rect.height - edge_rect.y,
                     BlockHeight[q_block.size] * (kMinBlockSizePix / 2));
  const uint32_t half = std::min({p_block_half, q_block_half, kDblkMaxHalf});
  const uint32_t edge_length =
      vertical_edge ? edge_rect.height : edge_rect.width;

  const bool has_lossless_alpha =
      (!p_block.has_lossy_alpha || !q_block.has_lossy_alpha);
  const uint32_t sharpness_from_quality = GetSharpnessFromQuality(quality_hint);

  assert(edge_length <= kMinBlockSizePix);
  // Store whether each alpha pixel line across the edge was deblocked.
  bool deblocked_A[kMinBlockSizePix];

  const bool vd_match =
      VDMatch(config, intertile ? "intertile-filter" : "deblocking-filter");

  for (Channel c : {kAChannel, kYChannel, kUChannel, kVChannel}) {
    if (c == kAChannel && !has_alpha) continue;

    const uint32_t strength = (c == kAChannel) ? alpha_strength : yuv_strength;
    if (strength == 0) {
      // If A is not deblocked, the remaining channels neither.
      if (c == kAChannel) break;
      continue;
    }

    const uint32_t num_bits =
        (c == kAChannel) ? (kAlphaBits + 1) : yuv_num_bits;
    const int32_t min = (c == kAChannel) ? 0 : yuv_min;
    const int32_t max = (c == kAChannel) ? (int32_t)kAlphaMax : yuv_max;
    const bool do_deblock = (c != kAChannel || !has_lossless_alpha);
    // TODO(maryla): also filter lossy blocks next to lossless ones?

    uint32_t sharpness;
    if (c == kAChannel && has_lossless_alpha) {
      // We can't compute sharpness from residual density since we
      // don't have it for lossless alpha blocks.
      sharpness = sharpness_from_quality;
    } else {
      // The sharpness is determined by the block with more details.
      sharpness = std::max(sharpness_from_quality,
                           GetSharpnessFromResidualDensity(std::max(
                               p_block.res_den[c], q_block.res_den[c])));
    }
    const int32_t deblock_threshold =
        DeblockThresholdFromSharpness(sharpness, num_bits);

    Plane16* const plane = &pixels->GetChannel(c);
    int16_t* q0 = &plane->At(edge_rect.x, edge_rect.y);
    const uint32_t step_along_edge = vertical_edge ? plane->Step() : 1u;
    const uint32_t step_across_edge = vertical_edge ? 1u : plane->Step();
    const bool is_chroma = (c == kUChannel || c == kVChannel);

    // Parse pixels of either a row or a column.
    for (uint32_t line = 0; line < edge_length; ++line, q0 += step_along_edge) {
      // If A is not deblocked, the remaining channels neither.
      if (c != kAChannel && has_alpha && !deblocked_A[line]) continue;

      const uint32_t flat_half =
          MeasureFlatLength(sharpness, half, q0, step_across_edge);

      if (vd_match) {
        RegisterPixelsForVDebug(
            config, tile_rect.x + edge_rect.x + (vertical_edge ? 0 : line),
            tile_rect.y + edge_rect.y + (vertical_edge ? line : 0), half,
            flat_half, vertical_edge, c, strength, sharpness,
            /*before_deblocking=*/true, do_deblock, num_bits, q0,
            step_across_edge);
      }

      const bool deblocked =
          do_deblock
              ? DeblockLine(strength, deblock_threshold, flat_half, is_chroma,
                            min, max, step_across_edge, q0)
              : WouldDeblockLine(strength, deblock_threshold, flat_half,
                                 is_chroma, min, max, step_across_edge, q0);
      if (c == kAChannel) deblocked_A[line] = deblocked;

      if (vd_match) {
        RegisterPixelsForVDebug(
            config, tile_rect.x + edge_rect.x + (vertical_edge ? 0 : line),
            tile_rect.y + edge_rect.y + (vertical_edge ? line : 0), half,
            flat_half, vertical_edge, c, strength, sharpness,
            /*before_deblocking=*/false, do_deblock && deblocked, num_bits, q0,
            step_across_edge);
      }
    }
  }
}

//------------------------------------------------------------------------------

}  // namespace WP2
