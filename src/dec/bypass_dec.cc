// Copyright 2020 Google LLC
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
//  Bypass with raw pixel coding
//
// Author: Yannis Guyon (yguyon@google.com)

#include "src/dec/tile_dec.h"
#include "src/dec/wp2_dec_i.h"
#include "src/dsp/math.h"
#include "src/wp2/base.h"
#include "src/wp2/format_constants.h"
#include "src/common/global_params.h"
#include "src/common/lossy/block.h"
#include "src/common/lossy/block_size.h"

namespace WP2 {

//------------------------------------------------------------------------------

WP2Status TileDecoder::BypassTileDec() {
  const uint32_t max_num_bytes = GetTileMaxNumBytes(*gparams_, tile_->rect);
  if (gparams_->type_ == GlobalParams::GP_LOSSLESS) {
    tile_->output_is_yuv = false;
    assert(!tile_->rgb_output.IsEmpty());
    const uint32_t bytes_per_row =
        tile_->rect.width * (gparams_->has_alpha_ ? 4 : 3);
    assert(bytes_per_row * tile_->rect.height == max_num_bytes);

    uint8_t* row = (uint8_t*)tile_->rgb_output.GetRow(0);
    tile_->num_decoded_rows = 0;
    while (tile_->num_decoded_rows < tile_->rect.height) {
      DataSource::DataHandle handle;
      WP2_CHECK_OK(tile_->input->TryReadNext(bytes_per_row, &handle),
                   WP2_STATUS_NOT_ENOUGH_DATA);  // TODO(yguyon): Or OOM?
      if (gparams_->has_alpha_) {
        std::memcpy(row, handle.GetBytes(), handle.GetSize());
      } else {
        for (uint32_t x = 0; x < tile_->rect.width; ++x) {
          row[x * 4 + 0] = kAlphaMax;
          row[x * 4 + 1] = handle.GetBytes()[x * 3 + 0];
          row[x * 4 + 2] = handle.GetBytes()[x * 3 + 1];
          row[x * 4 + 3] = handle.GetBytes()[x * 3 + 2];
        }
      }
      row += tile_->rgb_output.stride;
      ++tile_->num_decoded_rows;
      if (progress_ != nullptr) {
        WP2_CHECK_STATUS(progress_->AdvanceBy(tile_->rect.width));
      }
    }
  } else {  // !GP_LOSSLESS
    tile_->output_is_yuv = true;
    assert(!tile_->yuv_output.IsEmpty());
    const uint32_t num_bits = gparams_->transf_.GetYUVPrecisionBits() + 1u;
    const uint32_t min_num_bytes_per_px =
        ((gparams_->has_alpha_ ? kAlphaBits : 0) + num_bits * 3) / 8;

    tile_->block_map.Init(tile_->rect, num_bits, gparams_->transf_.GetYUVMin(),
                          gparams_->transf_.GetYUVMax(), &tile_->yuv_output);
    if (IsFilterBlockMapNeeded(*config_, *features_, *gparams_,
                               *tiles_layout_)) {
      WP2_CHECK_STATUS(tile_->block_map.Allocate());
    }
    // Simulate a grid of blocks for the IntertileFilter settings.
    // Even though this is lossless it is likely to border a lossy tile which
    // might benefit from slight deblocking.
    constexpr BlockSize kSimulateBlockSize = BLK_4x4;
    const uint32_t sim_block_height = BlockHeightPix(kSimulateBlockSize);
    const uint32_t sim_block_width = BlockWidthPix(kSimulateBlockSize);
    CodedBlock sim_block;  // For FilterBlockMap::RegisterBlock().
    sim_block.alpha_mode_ = kBlockAlphaLossless;
    sim_block.id_ = 0;  // Set segment id as the least quantized one.
    sim_block.is420_ = false;
    for (Channel c : {kYChannel, kUChannel, kVChannel, kAChannel}) {
      std::fill(&sim_block.coeffs_[c][0][0],  // Maximum residual density.
                &sim_block.coeffs_[c][0][0] + kMaxBlockSizePix2, 1);
    }

    DataSource::DataHandle handle;
    WP2_CHECK_OK(tile_->input->TryReadNext(max_num_bytes, &handle),
                 WP2_STATUS_NOT_ENOUGH_DATA);
    HeaderDec bit_unpacker(handle.GetBytes(), handle.GetSize());
    bit_unpacker.SetDebugPrefix("raw_px");
    if (!gparams_->has_alpha_ && !tile_->yuv_output.A.IsEmpty()) {
      tile_->yuv_output.A.Fill(kAlphaMax);
    }

    tile_->num_decoded_rows = 0;
    while (tile_->num_decoded_rows < tile_->rect.height) {
      const uint32_t y = tile_->num_decoded_rows;
      // Read interleaved channels to output samples line by line.
      for (Channel channel : {kYChannel, kUChannel, kVChannel}) {
        int16_t* const row = tile_->yuv_output.GetChannel(channel).Row(y);
        for (uint32_t x = 0; x < tile_->rect.width; ++x) {
          row[x] = (int16_t)bit_unpacker.ReadSBits(num_bits, "yuv");
        }
      }
      if (gparams_->has_alpha_) {
        int16_t* const row = tile_->yuv_output.A.Row(y);
        for (uint32_t x = 0; x < tile_->rect.width; ++x) {
          row[x] = (int16_t)bit_unpacker.ReadBits(kAlphaBits, "alpha");
        }
      }
      ++tile_->num_decoded_rows;
      if (tile_->num_decoded_rows % sim_block_height == 0 ||
          tile_->num_decoded_rows == tile_->rect.height) {
        // Register simulated blocks for intertile filtering when complete.
        for (uint32_t x = 0; x < tile_->rect.width; x += sim_block_width) {
          sim_block.SetDimDefault(
              {x / kMinBlockSizePix, y / kMinBlockSizePix, kSimulateBlockSize});
          tile_->block_map.RegisterBlock(sim_block, min_num_bytes_per_px);
        }
      }
      if (progress_ != nullptr) {
        WP2_CHECK_STATUS(progress_->AdvanceBy(tile_->rect.width));
      }
    }

    WP2_CHECK_OK(bit_unpacker.Ok(), WP2_STATUS_BAD_READ);
    assert(bit_unpacker.Used() == max_num_bytes);
  }
  return WP2_STATUS_OK;
}

//------------------------------------------------------------------------------

}   // namespace WP2
