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
// Functions related to tile decoding.
//
// Author: Skal (pascal.massimino@gmail.com)

#include "src/dec/tile_dec.h"

#include <algorithm>
#include <cstddef>

#include "src/common/global_params.h"
#include "src/common/lossy/block.h"
#include "src/common/lossy/block_size.h"
#include "src/dec/wp2_dec_i.h"
#include "src/dsp/math.h"
#include "src/wp2/base.h"
#include "src/wp2/format_constants.h"

namespace WP2 {

//------------------------------------------------------------------------------

void Tile::Clear() {
  chunk_size_is_known = false;
  chunk_size = 0;
  input = nullptr;
  private_input.Reset();
  data_handle = DataSource::DataHandle();

  rect = {};
  output_is_yuv = false;
  // Output is a view, no memory dealloc/realloc happens.
  assert(rgb_output.IsEmpty() || rgb_output.IsView());
  rgb_output.Deallocate();
  assert(yuv_output.IsEmpty() || yuv_output.IsView());
  yuv_output.Clear();
  num_decoded_rows = 0;
  block_map.Clear();
}

//------------------------------------------------------------------------------

const Tile& TilesLayout::GetTileAt(uint32_t x, uint32_t y) const {
  return tiles[(y / tile_height) * num_tiles_x + x / tile_width];
}

//------------------------------------------------------------------------------
// TileDecoder definitions.

WP2Status TileDecoder::TryAssignNextTile() {
  tile_ = nullptr;

  if (self_assigning_) {
    // Get first unassigned tile, and increment the counter if valid.
    WP2_CHECK_STATUS(tiles_layout_->assignment_lock.Acquire());
    uint32_t tile_index = tiles_layout_->first_unassigned_tile_index;
    if (tile_index < tiles_layout_->num_assignable_tiles) {
      tile_ = &tiles_layout_->tiles[tile_index];
      ++tiles_layout_->first_unassigned_tile_index;
    }
    tiles_layout_->assignment_lock.Release();
  }
  return WP2_STATUS_OK;
}

WP2Status TileDecoder::Execute() {
  while (tile_ != nullptr) {
    const uint32_t max_num_bytes = GetTileMaxNumBytes(*gparams_, tile_->rect);

    if (tile_->chunk_size == max_num_bytes) {
      // Raw samples are expected instead of the regular ANS decoding.
      WP2_CHECK_STATUS(BypassTileDec());
    } else {
      ANSDec dec(tile_->input);
      bool is_lossless;
      {
        ANSDebugPrefix prefix(&dec, "GlobalHeader");
        is_lossless =
            (gparams_->type_ == GlobalParams::GP_LOSSY) ? false :
            (gparams_->type_ == GlobalParams::GP_LOSSLESS) ? true :
            dec.ReadBool("use_lossless");
      }
      WP2_CHECK_STATUS(dec.GetStatus());
      WP2_CHECK_STATUS(is_lossless
                           ? LosslessDecode(*features_, *config_, *gparams_,
                                            progress_, &dec, tile_)
                           : LossyDecode(*features_, *config_, tiles_layout_,
                                         progress_, &dec, tile_));

#if defined(WP2_BITTRACE)
      if (config_->info != nullptr) {
        WP2_CHECK_STATUS(tiles_layout_->assignment_lock.Acquire());
        for (const auto& it : dec.GetBitTraces()) {
          config_->info->bit_traces[it.first].bits += it.second.bits;
          config_->info->bit_traces[it.first].num_occurrences +=
              it.second.num_occurrences;
          config_->info->bit_traces[it.first].type = it.second.type;
          for (const auto& p : it.second.histo) {
            config_->info->bit_traces[it.first].histo[p.first] += p.second;
          }
        }
        tiles_layout_->assignment_lock.Release();
      }
#endif
    }

    if (VDMatch(*config_, "blocks/partition")) tile_->Draw(*config_);

    WP2_CHECK_STATUS(TryAssignNextTile());
  }
  return WP2_STATUS_OK;
}

//------------------------------------------------------------------------------

WP2Status SetupWorkers(const BitstreamFeatures& features,
                       const DecoderConfig& config,
                       ProgressWatcher* const progress,
                       TilesLayout* const tiles_layout,
                       Vector<TileDecoder>* const workers) {
  WP2_CHECK_OK(config.thread_level <= 2147483647u,
               WP2_STATUS_INVALID_CONFIGURATION);
  // The number of workers is limited by the number of threads and tiles.
  const uint32_t num_workers =
      std::min(config.thread_level + 1, (uint32_t)tiles_layout->tiles.size());
  WP2_CHECK_ALLOC_OK(workers->resize(num_workers));
  for (TileDecoder& worker : *workers) {
    worker.features_ = &features;
    worker.config_ = &config;
    worker.gparams_ = tiles_layout->gparams;
    worker.progress_ = progress;
    worker.tiles_layout_ = tiles_layout;
    worker.self_assigning_ = true;
    WP2_CHECK_STATUS(worker.TryAssignNextTile());
  }
  return WP2_STATUS_OK;
}

//------------------------------------------------------------------------------

uint32_t GetNumTiles(uint32_t frame_width, uint32_t frame_height,
                     uint32_t tile_width, uint32_t tile_height,
                     uint32_t* const num_tiles_x, uint32_t* const num_tiles_y) {
  const uint32_t num_x = DivCeil(frame_width, tile_width);
  const uint32_t num_y = DivCeil(frame_height, tile_height);
  if (num_tiles_x != nullptr) *num_tiles_x = num_x;
  if (num_tiles_y != nullptr) *num_tiles_y = num_y;
  return num_x * num_y;
}

Rectangle GetTileRect(uint32_t frame_width, uint32_t frame_height,
                      uint32_t tile_width, uint32_t tile_height,
                      uint32_t tile_index) {
  uint32_t num_tiles_x;
  if (tile_index >= GetNumTiles(frame_width, frame_height, tile_width,
                                tile_height, &num_tiles_x)) {
    assert(false);
  }
  const uint32_t tile_x = (tile_index % num_tiles_x) * tile_width;
  const uint32_t tile_y = (tile_index / num_tiles_x) * tile_height;
  return {tile_x, tile_y, std::min(frame_width - tile_x, tile_width),
          std::min(frame_height - tile_y, tile_height)};
}

WP2Status GetTilesLayout(uint32_t frame_width, uint32_t frame_height,
                         uint32_t tile_width, uint32_t tile_height,
                         ArgbBuffer* const rgb_output,
                         YUVPlane* const yuv_output,
                         TilesLayout* const tiles_layout) {
  WP2_CHECK_OK(frame_width > 0 && frame_height > 0, WP2_STATUS_BITSTREAM_ERROR);
  WP2_CHECK_OK(frame_width <= kImageDimMax && frame_height <= kImageDimMax,
               WP2_STATUS_BITSTREAM_ERROR);
  tiles_layout->tile_width = tile_width;
  tiles_layout->tile_height = tile_height;
  const uint32_t num_tiles =
      GetNumTiles(frame_width, frame_height, tile_width, tile_height,
                  &tiles_layout->num_tiles_x, &tiles_layout->num_tiles_y);
  WP2_CHECK_ALLOC_OK(tiles_layout->tiles.resize(num_tiles));

  for (uint32_t tile_index = 0; tile_index < num_tiles; ++tile_index) {
    Tile& tile = tiles_layout->tiles[tile_index];
    tile.Clear();
    tile.rect = GetTileRect(frame_width, frame_height, tile_width, tile_height,
                            tile_index);
    // Set views on RGB and YUV output buffers, both can be used.
    if (!rgb_output->IsEmpty()) {
      assert(rgb_output->width == frame_width &&
             rgb_output->height == frame_height);
      WP2_CHECK_STATUS(tile.rgb_output.SetFormat(rgb_output->format));
      WP2_CHECK_STATUS(tile.rgb_output.SetView(*rgb_output, tile.rect));
      tile.output_is_yuv = false;
    }
    if (!yuv_output->IsEmpty()) {
      assert(yuv_output->Y.w_ == Pad(frame_width, kPredWidth));
      assert(yuv_output->Y.h_ == Pad(frame_height, kPredWidth));
      WP2_CHECK_STATUS(tile.yuv_output.SetView(
          *yuv_output,
          {tile.rect.x, tile.rect.y, Pad(tile.rect.width, kPredWidth),
           Pad(tile.rect.height, kPredWidth)}));
      tile.output_is_yuv = true;  // Overridden in Lossy/LosslessDecode() anyway
    }
  }
  tiles_layout->num_assignable_tiles = 0;
  tiles_layout->first_unassigned_tile_index = 0;
  tiles_layout->num_decoded_tiles = 0;
  return WP2_STATUS_OK;
}

uint32_t GetTileMaxNumBytes(const GlobalParams& gparams,
                            const Rectangle& tile_rect) {
#if defined(WP2_ENC_DEC_MATCH)
  (void)gparams, (void)tile_rect;
  // Make sure it does not fall back to raw pixel coding in CodeTiles().
  constexpr uint32_t kPlentyEnoughNumBytes = 2147483647u;
  static_assert(kPlentyEnoughNumBytes < kUpperVarInt, "kPlentyEnoughNumBytes");
  return kPlentyEnoughNumBytes;
#else
  const bool lossless = (gparams.type_ == GlobalParams::GP_LOSSLESS);
  const uint32_t channel_bits =
      lossless ? 8 : (gparams.transf_.GetYUVPrecisionBits() + 1u);
  const uint32_t num_pixels = tile_rect.width * tile_rect.height;
  const uint32_t alpha_bits = gparams.has_alpha_ ? kAlphaBits : 0;
  const uint32_t max_num_bits = (alpha_bits + 3 * channel_bits) * num_pixels;
  return DivCeil(max_num_bits, 8);
#endif  // defined(WP2_ENC_DEC_MATCH)
}

//------------------------------------------------------------------------------

// Reads the size of each tile and skips them.
WP2Status SkipTiles(DataSource* const data_source,
                    const BitstreamFeatures& features,
                    const Rectangle& window) {
  GlobalParams gparams;
  WP2_CHECK_STATUS(
      DecodeGLBL(data_source, DecoderConfig::kDefault, features, &gparams));
  const uint32_t num_tiles = GetNumTiles(
      window.width, window.height, features.tile_width, features.tile_height);
  for (uint32_t tile_index = 0; tile_index < num_tiles; ++tile_index) {
    const Rectangle& tile_rect =
        GetTileRect(window.width, window.height, features.tile_width,
                    features.tile_height, tile_index);
    size_t tile_chunk_size;
    WP2_CHECK_STATUS(
        DecodeTileChunkSize(gparams, tile_rect, data_source, &tile_chunk_size));
    data_source->MarkNumBytesAsRead(tile_chunk_size);
  }
  return WP2_STATUS_OK;
}

//------------------------------------------------------------------------------

}   // namespace WP2
