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
//  incremental decoding
//
// Author: Skal (pascal.massimino@gmail.com)

#include <algorithm>

#include "src/dec/incr/decoder_state.h"
#include "src/dec/wp2_dec_i.h"
#include "src/utils/data_source.h"
#include "src/wp2/decode.h"
#include "src/common/global_params.h"

namespace WP2 {

//------------------------------------------------------------------------------

namespace {

bool TryReadTileChunkSize(const GlobalParams& gparams,
                          DataSource* const data_source,
                          const DecoderConfig& config, Tile* const tile) {
  assert(tile->rect.GetArea() > 0);
  if (!tile->chunk_size_is_known) {  // If not previously decoded.
    const size_t size_before = data_source->GetNumReadBytes();
    const WP2Status chunk_size_status = DecodeTileChunkSize(
        gparams, tile->rect, data_source, &tile->chunk_size);
    if (chunk_size_status != WP2_STATUS_OK) {
      // Tile chunk size is unknown, don't decode this tile for now.
      // No byte is marked as read by DecodeTileChunkSize() if not STATUS_OK,
      // so UnmarkAllReadBytes() is not needed.
      return false;
    }
    tile->chunk_size_is_known = true;
    RegisterChunkSize(config, tile->chunk_size,
                      data_source->GetNumReadBytes() - size_before);
  }
  return true;
}

// Returns true if the next tile has all its chunk data available.
bool TryReadTile(DataSource* const data_source,
                 const BitstreamFeatures& features, const DecoderConfig& config,
                 TilesLayout* const tiles_layout, Tile* const tile,
                 PartialTileDecoder* const partial_tile_decoder) {
  if (!TryReadTileChunkSize(*tiles_layout->gparams, data_source, config,
                            tile)) {
    return false;
  }

  if (data_source->TryReadNext(tile->chunk_size, &tile->data_handle)) {
    // All of the tile chunk data is available, decode it in one shot.
    ++tiles_layout->num_assignable_tiles;
  } else {
    if (partial_tile_decoder != nullptr) {
      // TryGetNext() should always be able to return GetNumNextBytes(), but
      // just in case, decode previous tiles if possible instead of returning
      // an error right now.
      if (data_source->TryGetNext(data_source->GetNumNextBytes(),
                                  &tile->data_handle)) {
        if (partial_tile_decoder->Init(features, config,
                                       tiles_layout, tile) != WP2_STATUS_OK) {
          return false;  // this is most probably a memory alloc failure
        }
      }
    }
    return false;  // Part of the data is missing.
  }
  return true;
}

bool TrySkipTile(DataSource* const data_source, const DecoderConfig& config,
                 TilesLayout* const tiles_layout, Tile* const tile) {
  if (!TryReadTileChunkSize(*tiles_layout->gparams, data_source, config,
                            tile)) {
    return false;
  }
  data_source->MarkNumBytesAsRead(tile->chunk_size);
  ++tiles_layout->num_assignable_tiles;
  return true;
}

}  // namespace

//------------------------------------------------------------------------------

WP2Status Decoder::DecodeAvailableTiles() {
  TilesLayout* const tiles_layout = &state_->tiles_layout;
  PartialTileDecoder* const partial_tile_decoder =
      state_->partial_tile_decoder.get();
  if ((partial_tile_decoder != nullptr) &&
      (partial_tile_decoder->GetPartialTile() != nullptr)) {
    Tile* const partial_tile = partial_tile_decoder->GetPartialTile();
    // DataSource::GetNumNextBytes() gets updated on TryGetNext() so try to get
    // the whole tile chunk first.
    if (!state_->data_source->TryGetNext(partial_tile->chunk_size,
                                         &partial_tile->data_handle)) {
      WP2_CHECK_OK(state_->data_source->TryGetNext(
                       state_->data_source->GetNumNextBytes(),
                       &partial_tile->data_handle),
                   WP2_STATUS_BITSTREAM_OUT_OF_MEMORY);
    }

    // If the partial tile is not fully decoded (WP2_STATUS_OK), don't look
    // further.
    WP2_CHECK_STATUS(partial_tile_decoder->Continue());

    state_->data_source->MarkNumBytesAsRead(partial_tile->chunk_size);
    state_->data_source->Discard(state_->data_source->GetNumReadBytes());
    partial_tile_decoder->Clear();
    ++tiles_layout->num_decoded_tiles;
  }

  while (tiles_layout->num_assignable_tiles < tiles_layout->tiles.size()) {
    const uint32_t tile_index = tiles_layout->num_assignable_tiles;
    if (!TryReadTile(state_->data_source, features_, config_,
                     &state_->tiles_layout,
                     &state_->tiles_layout.tiles[tile_index],
                     state_->partial_tile_decoder.get())) {
      break;
    }
  }

  // Now that the primary 'state_->data_source' was read as much as possible and
  // won't be used till next DecodeAvailableTiles(), use ExternalDataSources as
  // input for each tile, to prevent slowness and concurrency issues.
  for (uint32_t tile_index = tiles_layout->num_decoded_tiles;
       tile_index < tiles_layout->num_assignable_tiles; ++tile_index) {
    Tile* const tile = &tiles_layout->tiles[tile_index];
    const uint8_t* const tile_data = tile->data_handle.GetBytes();
    // Maybe a DataSource::TryReadNext() invalidated previously read data.
    WP2_CHECK_OK(tile_data != nullptr, WP2_STATUS_BITSTREAM_OUT_OF_MEMORY);
    tile->private_input = ExternalDataSource(tile_data, tile->chunk_size);
    tile->input = &tile->private_input;
  }

  // Assign first tiles to available workers.
  const uint32_t num_full_tiles =
      tiles_layout->num_assignable_tiles - tiles_layout->num_decoded_tiles;
  const uint32_t num_workers =
      std::min((uint32_t)state_->workers.size(), num_full_tiles);
  for (uint32_t i = 0; i < num_workers; ++i) {
    TileDecoder* const worker = &state_->workers[i];
    WP2_CHECK_STATUS(worker->TryAssignNextTile());
    assert(worker->tile_ != nullptr);  // Has a job.
  }

  // Decode full tiles in parallel or sequentially.
  WP2Status status = WP2_STATUS_OK;
  for (uint32_t i = 0; i < num_workers; ++i) {
    TileDecoder* const worker = &state_->workers[i];

    // The last worker runs in the main thread, after starting others.
    const bool threaded = (i + 1 != num_workers);
    status = worker->Start(threaded);
    // Only the 'partial_tile_decoder' is allowed to return missing data.
    assert(status != WP2_STATUS_NOT_ENOUGH_DATA);
    if (status != WP2_STATUS_OK) break;
  }

  if ((status == WP2_STATUS_OK) && (partial_tile_decoder != nullptr) &&
      (partial_tile_decoder->GetPartialTile() != nullptr)) {
    // Start partial tile decoding, worker will suspend itself.
    status = partial_tile_decoder->Start();
  }

  // Loop over all threads and close them properly. No early return.
  for (uint32_t i = 0; i < num_workers; ++i) {
    TileDecoder* const worker = &state_->workers[i];
    WP2Status worker_status = worker->End();
    worker->tile_ = nullptr;
    // Only the 'partial_tile_decoder' is allowed to return missing data.
    assert(worker_status != WP2_STATUS_NOT_ENOUGH_DATA);

    if ((worker_status != WP2_STATUS_OK) && (status == WP2_STATUS_OK)) {
      // The first encountered error is returned (by tile order, not time).
      status = worker_status;
    }
  }

  if ((partial_tile_decoder != nullptr) &&
      (partial_tile_decoder->GetPartialTile() != nullptr)) {
    if (status == WP2_STATUS_OK) {
      // Don't increase these earlier or regular workers might try to decode the
      // partial tile.
      ++tiles_layout->num_assignable_tiles;
      ++tiles_layout->first_unassigned_tile_index;
    } else {
      partial_tile_decoder->Clear();
    }
  }

  WP2_CHECK_STATUS(status);

  state_->data_source->Discard(state_->data_source->GetNumReadBytes());
  // Increase num_decoded_tiles only if there is no error. We don't know here
  // which tile failed.
  tiles_layout->num_decoded_tiles += num_full_tiles;

  state_->intertile_filter.Deblock(state_->GetFrameNumDecodedRows());

  const uint32_t num_converted_rows = state_->frame_num_converted_rows;
  const uint32_t num_convertible_rows =
      state_->intertile_filter.GetNumFilteredRows();

  if (num_converted_rows < num_convertible_rows) {
    // Per-tile YUV to RGB conversion, if any.
    for (Tile& tile : tiles_layout->tiles) {
      // This Decoder does not support YUV output, final buffer is always RGB.
      assert(!tile.rgb_output.IsEmpty());
      if (!tile.output_is_yuv) continue;

      if (tile.rect.y + tile.rect.height <= num_converted_rows) continue;
      if (tile.rect.y >= num_convertible_rows) break;
      const uint32_t num_converted_rows_in_tile =
          std::max(tile.rect.y, num_converted_rows) - tile.rect.y;
      const uint32_t num_convertible_rows_in_tile =
          std::min(tile.rect.height, num_convertible_rows - tile.rect.y);
      assert(num_converted_rows_in_tile < num_convertible_rows_in_tile);
      assert(num_convertible_rows_in_tile <= tile.num_decoded_rows);

      const Rectangle new_rows_rect = {
          0, num_converted_rows_in_tile, tile.rect.width,
          num_convertible_rows_in_tile - num_converted_rows_in_tile};

      ArgbBuffer rgb_output_rows(tile.rgb_output.format);
      WP2_CHECK_STATUS(rgb_output_rows.SetView(tile.rgb_output, new_rows_rect));
      YUVPlane yuv_output_rows;
      WP2_CHECK_STATUS(yuv_output_rows.SetView(tile.yuv_output, new_rows_rect));
      WP2_CHECK_STATUS(yuv_output_rows.Export(tiles_layout->gparams->transf_,
                                              /*resize_if_needed=*/false,
                                              &rgb_output_rows));
    }
  }

  state_->frame_num_converted_rows = num_convertible_rows;
  state_->frame_num_final_rows = state_->frame_num_converted_rows;

  if (tiles_layout->num_decoded_tiles < tiles_layout->tiles.size()) {
    return WP2_STATUS_NOT_ENOUGH_DATA;
  }
  return WP2_STATUS_OK;
}

WP2Status Decoder::DiscardAvailableTiles() {
  TilesLayout* const tiles_layout = &state_->tiles_layout;
  PartialTileDecoder* const partial_tile_decoder =
      state_->partial_tile_decoder.get();
  if ((partial_tile_decoder != nullptr) &&
      (partial_tile_decoder->GetPartialTile() != nullptr)) {
    state_->data_source->MarkNumBytesAsRead(
        partial_tile_decoder->GetPartialTile()->chunk_size);
    state_->data_source->Discard(state_->data_source->GetNumReadBytes());
    partial_tile_decoder->Clear();
  }

  while (tiles_layout->num_assignable_tiles < tiles_layout->tiles.size()) {
    const uint32_t tile_index = tiles_layout->num_assignable_tiles;
    const bool skipped_tile =
        TrySkipTile(state_->data_source, config_, &state_->tiles_layout,
                    &state_->tiles_layout.tiles[tile_index]);
    state_->data_source->Discard(state_->data_source->GetNumReadBytes());
    if (!skipped_tile) break;
  }

  if (tiles_layout->num_assignable_tiles < tiles_layout->tiles.size()) {
    return WP2_STATUS_NOT_ENOUGH_DATA;
  }
  return WP2_STATUS_OK;
}

//------------------------------------------------------------------------------

}  // namespace WP2
