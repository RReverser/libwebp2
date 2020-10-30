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
// Grain 'filter'.
//
// Author: Skal (pascal.massimino@gmail.com)

#include "src/dec/filters/grain_filter.h"

#include "src/dec/wp2_dec_i.h"
#include "src/dsp/math.h"

namespace WP2 {

//------------------------------------------------------------------------------

bool IsGrainFilterEnabled(const DecoderConfig& config,
                          const Vector<Segment>& segments) {
  for (const Segment& segment : segments) {
    if (segment.grain_.IsUsed()) return true;
  }
  return (VDMatch(config, "grain-filter"));
}

//------------------------------------------------------------------------------

GrainFilter::GrainFilter(const DecoderConfig& config,
                         const Vector<Segment>& segments,
                         const FilterBlockMap& blocks)
    : config_(config),
      segments_(segments),
      blocks_(blocks),
      enabled_(IsGrainFilterEnabled(config, segments)),
      last_row_(0) {}

WP2Status GrainFilter::Allocate() {
  if (!enabled_) return WP2_STATUS_OK;

  GrainFilterInit();
  return WP2_STATUS_OK;
}

//------------------------------------------------------------------------------

uint32_t GrainFilter::Apply(uint32_t num_rows) {
  if (!enabled_) return num_rows;

  if (VDMatch(config_, "grain-filter")) SavePixelsForVDebug(num_rows);

  assert(num_rows <= blocks_.tile_rect_.height);
  YUVPlane* const pixels = blocks_.pixels_;
  assert(pixels->GetHeight() % kMinBlockSizePix == 0);
  const uint32_t block_step = blocks_.max_num_blocks_x_;
  uint32_t y = last_row_;
  assert(y % kMinBlockSizePix == 0);

  const uint32_t num_padded_rows =
      (num_rows == blocks_.tile_rect_.height) ? pixels->GetHeight() : num_rows;
  for (; y + kMinBlockSizePix <= num_padded_rows; y += kMinBlockSizePix) {
    const uint32_t block_y = y / kMinBlockSizePix;
    const FilterBlockMap::BlockFeatures* row_block_features =
        &blocks_.features_[block_y * block_step];
    for (uint32_t block_x = 0; block_x < block_step; ++block_x) {
      const FilterBlockMap::BlockFeatures& blk = row_block_features[block_x];
      const Segment& seg = segments_[blk.segment_id];
      const GrainParams& grain = seg.grain_;
      const uint32_t x_pix = block_x * kMinBlockSizePix;
      AddGrain4x4(&pixels->Y.At(x_pix, y), pixels->Y.Step(), &rng_,
                  grain.y_, grain.cut_y_);
      AddGrain4x4(&pixels->U.At(x_pix, y), pixels->U.Step(), &rng_,
                  grain.uv_, grain.cut_uv_);
      AddGrain4x4(&pixels->V.At(x_pix, y), pixels->V.Step(), &rng_,
                  grain.uv_, grain.cut_uv_);
    }
  }
  last_row_ = y;

  if (VDMatch(config_, "grain-filter")) ApplyVDebug(num_rows);
  return std::min(last_row_, num_rows);
}

//------------------------------------------------------------------------------

}  // namespace WP2
