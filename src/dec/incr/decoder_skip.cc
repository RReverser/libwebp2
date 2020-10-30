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
// Functions used to skip frames and glimpse at future frame features.
//
// Author: Yannis Guyon (yguyon@google.com)

#include "src/dec/incr/decoder_skip.h"

#include <algorithm>

#include "src/dec/incr/decoder_info.h"
#include "src/dec/incr/decoder_state.h"
#include "src/dec/wp2_dec_i.h"
#include "src/utils/orientation.h"
#include "src/utils/utils.h"
#include "src/wp2/decode.h"

namespace WP2 {

//------------------------------------------------------------------------------

void CopyGlimpseFromCurrent(Decoder::State* const state) {
  state->glimpsed_bitstream_position =
      state->data_source->GetNumDiscardedBytes();
  state->glimpsed_frame_index = state->current_frame_index;
  if (state->stage < Decoder::State::Stage::WAITING_FOR_TILE_DATA) {
    state->glimpsed_frame_anmf =
        (state->stage > Decoder::State::Stage::WAITING_FOR_ANMF_CHUNK);
    state->glimpsed_frame_glbl = false;
    state->glimpsed_frame_gparams.Reset();
    state->num_glimpsed_tiles = 0;
  } else {
    state->glimpsed_frame_anmf = true;
    state->glimpsed_frame_glbl = true;
    // DecodeTileChunkSize() needs only GlobalParamas::type_,transf_,has_alpha_.
    state->glimpsed_frame_gparams.type_ = state->gparams.type_;
    state->glimpsed_frame_gparams.transf_ = state->gparams.transf_;
    state->glimpsed_frame_gparams.has_alpha_ = state->gparams.has_alpha_;
    state->num_glimpsed_tiles = state->tiles_layout.num_decoded_tiles;
    if (state->tiles_layout.num_decoded_tiles <
        state->tiles_layout.tiles.size()) {
      const Tile& partial_tile =
          state->tiles_layout.tiles[state->tiles_layout.num_decoded_tiles];
      if (partial_tile.chunk_size_is_known) {
        ++state->num_glimpsed_tiles;
        state->glimpsed_bitstream_position += partial_tile.chunk_size;
      }
    }
  }
}

WP2Status GlimpseTillNextANMFChunk(const DecoderConfig& config,
                                   const BitstreamFeatures& features,
                                   Decoder::State* const state) {
  // Only the next ANMF chunk is wanted; skip any tile before it.
  if (state->glimpsed_frame_anmf && state->glimpsed_frame_glbl) {
    const Decoder::State::InternalFrameFeatures& glimpsed_frame_features =
        state->frames[state->glimpsed_frame_index];
    const Rectangle& window = glimpsed_frame_features.info.window;
    const uint32_t num_tiles = GetNumTiles(
        window.width, window.height, features.tile_width, features.tile_height);
    while (state->num_glimpsed_tiles < num_tiles) {
      const Rectangle& tile_rect =
          GetTileRect(window.width, window.height, features.tile_width,
                      features.tile_height, state->num_glimpsed_tiles);
      size_t tile_chunk_size = 0;
      WP2_CHECK_STATUS(DecodeTileChunkSize(state->glimpsed_frame_gparams,
                                           tile_rect, state->data_source,
                                           &tile_chunk_size));
      state->data_source->MarkNumBytesAsRead(tile_chunk_size);
      ++state->num_glimpsed_tiles;
    }
    WP2_CHECK_STATUS(
        SaveFrameNumBytes(state->data_source->GetNumDiscardedBytes() +
                              state->data_source->GetNumReadBytes(),
                          &state->frames[state->glimpsed_frame_index]));
    if (glimpsed_frame_features.duration_ms > 0) {  // Regular frame.
      state->num_available_final_frames =
          std::max(state->num_available_final_frames,
                   glimpsed_frame_features.final_frame_index + 1u);
    }
    if (!glimpsed_frame_features.info.is_last) {
      ++state->glimpsed_frame_index;
      state->glimpsed_frame_anmf = false;
      state->glimpsed_frame_glbl = false;
      state->glimpsed_frame_gparams.Reset();
      state->num_glimpsed_tiles = 0;
    }
  }
  if (!state->glimpsed_frame_anmf) {
    WP2_CHECK_STATUS(DecodeFrameInfo(state->data_source, features,
                                     /*glimpse=*/true, state));
    state->glimpsed_frame_anmf = true;
    assert(!state->glimpsed_frame_glbl);
  }
  if (!state->glimpsed_frame_glbl) {
    WP2_CHECK_STATUS(DecodeGLBL(state->data_source, config, features,
                                &state->glimpsed_frame_gparams));
    state->glimpsed_frame_glbl = true;
    assert(state->glimpsed_frame_anmf);
  }
  return WP2_STATUS_OK;
}

WP2Status GlimpseAtNextFrame(const DecoderConfig& config,
                             const BitstreamFeatures& features,
                             Decoder::State* const state) {
  assert(state->data_source->GetNumReadBytes() == 0);
  const size_t current_bitstream_position =
      state->data_source->GetNumDiscardedBytes();
  if (state->glimpsed_bitstream_position < current_bitstream_position) {
    CopyGlimpseFromCurrent(state);
  }

  state->data_source->MarkNumBytesAsRead(state->glimpsed_bitstream_position -
                                         current_bitstream_position);
  WP2Status status = GlimpseTillNextANMFChunk(config, features, state);
  state->glimpsed_bitstream_position =
      current_bitstream_position + state->data_source->GetNumReadBytes();
  state->data_source->UnmarkAllReadBytes();
  return status;
}

//------------------------------------------------------------------------------

bool FrameIsAlreadyInOutputBuffer(uint32_t tile_width, uint32_t tile_height,
                                  Decoder::State* const state) {
  uint32_t target_frame_index = state->current_frame_index;
  if (state->skip_final_frames_before_index <
      state->final_to_regular_index.size()) {
    target_frame_index = std::max(
        target_frame_index,
        state->final_to_regular_index[state->skip_final_frames_before_index]);
  }
  if (state->index_of_frame_in_output_buffer == target_frame_index) {
    assert(state->frames[target_frame_index].num_bytes !=
           Decoder::State::kUnknown);
    const size_t target_bitstream_position =
        state->frames[target_frame_index].bitstream_position +
        state->frames[target_frame_index].num_bytes;
    const size_t bitstream_position =
        state->data_source->GetNumDiscardedBytes();

    if (bitstream_position > target_bitstream_position) {
      state->status = WP2_STATUS_BAD_READ;  // Data source changed.
    } else {
      state->data_source->Discard(target_bitstream_position -
                                  bitstream_position);
      state->current_frame_index = target_frame_index;

      // Simulate that all tiles are decoded.
      const Rectangle& window = state->frames[target_frame_index].info.window;
      state->status =
          GetTilesLayout(window.width, window.height, tile_width, tile_height,
                         &state->frame_output_buffer, &state->yuv_output_buffer,
                         &state->tiles_layout);
      state->tiles_layout.num_assignable_tiles =
          state->tiles_layout.tiles.size();
      state->tiles_layout.first_unassigned_tile_index =
          state->tiles_layout.num_assignable_tiles;
      state->tiles_layout.num_decoded_tiles =
          state->tiles_layout.num_assignable_tiles;
      state->frame_num_converted_rows = window.height;
      state->frame_num_final_rows = window.height;
      state->intertile_filter.Clear();
      // TODO(yguyon): Store one or two pixels when setting
      //               'index_of_frame_in_output_buffer' and assert them here.
    }
    return true;  // Going directly to WAITING_FOR_TILE_DATA is expected.
  }

  if (state->index_of_frame_in_output_buffer >=
          state->discard_frames_before_index &&
      state->index_of_frame_in_output_buffer < target_frame_index) {
    // Even more frames can be discarded.
    state->discard_frames_before_index =
        state->index_of_frame_in_output_buffer + 1;
    // Rely on the following DiscardWholeFrames() to actually do the work.
  }
  return false;
}

//------------------------------------------------------------------------------

bool DiscardWholeFrames(Decoder::State* const state) {
  if (state->frames.empty()) return false;
  uint32_t target_frame_index = std::min(state->discard_frames_before_index,
                                         (uint32_t)state->frames.size());
  if (target_frame_index <= state->current_frame_index) return false;

  size_t target_bitstream_position =
      state->frames[target_frame_index - 1].bitstream_position;
  if (state->frames[target_frame_index - 1].num_bytes !=
      Decoder::State::kUnknown) {
    target_bitstream_position +=
        state->frames[target_frame_index - 1].num_bytes;
  } else {
    // Previous 'size' is not yet available so one less discarded frame.
    target_frame_index -= 1;
    if (target_frame_index <= state->current_frame_index) return false;
  }

  const size_t bitstream_position = state->data_source->GetNumDiscardedBytes();
  if (bitstream_position > target_bitstream_position) {
    state->status = WP2_STATUS_BAD_READ;  // Data source changed.
  } else {
    state->data_source->Discard(target_bitstream_position - bitstream_position);
    state->current_frame_index = target_frame_index;
    state->frame_num_converted_rows = 0;
    state->frame_num_final_rows = 0;
    state->intertile_filter.Clear();
  }
  return true;  // Going directly to WAITING_FOR_ANMF_CHUNK is expected.
}

bool DiscardTillMetadata(Decoder::State* const state) {
  if (state->ShouldDiscardAllFrames() && !state->frames.empty() &&
      state->frames.back().is_last &&
      state->frames.back().num_bytes != Decoder::State::kUnknown) {
    const size_t metadata_bitstream_position =
        state->frames.back().bitstream_position +
        state->frames.back().num_bytes;
    const size_t bitstream_position =
        state->data_source->GetNumDiscardedBytes();
    if (bitstream_position > metadata_bitstream_position) {
      // This function is called only during frame decoding, so being past the
      // first byte of metadata can only mean that the data source changed.
      state->status = WP2_STATUS_BAD_READ;
      return false;
    }
    state->data_source->Discard(metadata_bitstream_position -
                                bitstream_position);
    // If 'metadata_bitstream_position' is known, so are all frames.
    assert(!state->frames.empty() && state->frames.back().info.is_last);
    state->current_frame_index = (uint32_t)state->frames.size() - 1;
    state->frame_num_converted_rows = 0;
    state->frame_num_final_rows = 0;
    state->intertile_filter.Clear();
    return true;  // Going directly to WAITING_FOR_METADATA_CHUNKS is expected.
  }
  return false;
}

//------------------------------------------------------------------------------

}  // namespace WP2
