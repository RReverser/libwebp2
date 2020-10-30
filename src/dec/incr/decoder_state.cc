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
//  Incremental decoding state getters
//
// Author: Yannis Guyon (yguyon@google.com)

#include "src/dec/incr/decoder_state.h"

#include <algorithm>

namespace WP2 {

//------------------------------------------------------------------------------

void Decoder::State::Rewind(bool output_buffer_changed) {
  if (partial_tile_decoder != nullptr) partial_tile_decoder->Clear();
  stage = Decoder::State::Stage::WAITING_FOR_HEADER_CHUNK;
  status = WP2_STATUS_NOT_ENOUGH_DATA;
  current_frame_index = 0;
  current_frame_is_ready = false;
  tiles_layout.num_decoded_tiles = 0;
  frame_num_converted_rows = 0;
  frame_num_final_rows = 0;
  intertile_filter.Clear();
  data_source->Reset();
  if (output_buffer_changed) index_of_frame_in_output_buffer = kMaxNumFrames;
  skip_final_frames_before_index = 0;
  discard_frames_before_index = 0;
  // Keep the glimpsed frame position and the already decoded frame features.
}

void DecoderStateArray::Rewind(bool output_buffer_changed) {
  Decoder::State::Rewind(output_buffer_changed);
  array_data_source.Update(last_input_set.bytes, last_input_set.size);
}

//------------------------------------------------------------------------------

bool Decoder::State::InternalFrameFeatures::operator==(
    const InternalFrameFeatures& other) const {
  return duration_ms == other.duration_ms &&
         window.x == other.window.x &&
         window.y == other.window.y &&
         window.width == other.window.width &&
         window.height == other.window.height &&
         last_dispose_frame_index == other.last_dispose_frame_index &&
         bitstream_position == other.bitstream_position &&
         num_bytes == other.num_bytes &&
         info.dispose == other.info.dispose &&
         info.blend == other.info.blend &&
         info.duration_ms == other.info.duration_ms &&
         info.window.x == other.info.window.x &&
         info.window.y == other.info.window.y &&
         info.window.width == other.info.window.width &&
         info.window.height == other.info.window.height &&
         info.is_last == other.info.is_last &&
         is_last == other.is_last &&
         final_frame_index == other.final_frame_index;
}

bool Decoder::State::InternalFrameFeatures::operator!=(
    const InternalFrameFeatures& other) const {
  return !operator==(other);
}

//------------------------------------------------------------------------------

bool Decoder::State::NextStage() {
  if (status != WP2_STATUS_OK) {
    data_source->UnmarkAllReadBytes();
    return false;
  }
  data_source->Discard(data_source->GetNumReadBytes());
  assert(stage != DECODED);
  stage = (Stage)(stage + 1);
  return true;
}

//------------------------------------------------------------------------------

uint32_t Decoder::State::GetFinalFrameIndex(uint32_t frame_index) const {
  if (frame_index < frames.size()) {
    return frames[frame_index].final_frame_index;
  } else if (frame_index == frames.size()) {
    if (frame_index == 0) {
      return 0;
    } else if (frames.back().info.duration_ms == 0) {
      // The previous one is a preframe, so it's the same final frame index.
      return frames.back().final_frame_index;
    } else {
      return frames.back().final_frame_index + 1;
    }
  }
  assert(false);  // Should not happen, frame_index <= frames.size()
  return 0;
}

bool Decoder::State::TryGetFirstFrameBelongingToFinalFrame(
    uint32_t final_frame_index, uint32_t* const frame_index) const {
  if (final_frame_index == 0) {
    *frame_index = 0;
    return true;
  }
  const uint32_t previous_final_frame_index = final_frame_index - 1;
  if (previous_final_frame_index < final_to_regular_index.size()) {
    const uint32_t previous_regular_frame_index =
        final_to_regular_index[previous_final_frame_index];
    *frame_index = previous_regular_frame_index + 1;
    return true;
  }
  return false;
}

WP2Status Decoder::State::ComputeNumFramesToDiscard() {
  // Quick return for most cases (no SkipNumNextFrames() was called).
  if (skip_final_frames_before_index == 0) return WP2_STATUS_OK;

  if (GetFinalFrameIndex(current_frame_index) >=
      skip_final_frames_before_index) {
    // There's no point in computing 'discard_frames_before_index', it will be
    // before the current frame anyway.
    return WP2_STATUS_OK;
  }
  if (skip_final_frames_before_index >= kMaxNumFrames) {
    discard_frames_before_index = kMaxNumFrames;
    return WP2_STATUS_OK;
  }

  if (skip_final_frames_before_index < final_to_regular_index.size()) {
    const uint32_t target_frame_index =
        final_to_regular_index[skip_final_frames_before_index];

    const uint32_t index_of_first_final_frame_to_keep =
        frames[target_frame_index].last_dispose_frame_index;

    uint32_t index_of_first_frame_to_keep;
    if (TryGetFirstFrameBelongingToFinalFrame(
            index_of_first_final_frame_to_keep,
            &index_of_first_frame_to_keep)) {
      discard_frames_before_index = index_of_first_frame_to_keep;
      return WP2_STATUS_OK;
    }
  }
  if (!frames.empty() && frames.back().info.is_last) {
    discard_frames_before_index = (uint32_t)frames.size();
    return WP2_STATUS_OK;
  }
  // TODO(yguyon): Assumptions can be made with kMaxNumDependentFrames to
  //               discard some frames.
  return WP2_STATUS_NOT_ENOUGH_DATA;
}

// If this returns true once, it will do so until Decoder::Rewind().
bool Decoder::State::ShouldDiscardAllFrames() const {
  if (!frames.empty() && frames.back().info.is_last) {
    return (discard_frames_before_index >= frames.size());
  }
  return (skip_final_frames_before_index >= kMaxNumFrames);
}

//------------------------------------------------------------------------------

uint32_t Decoder::State::GetFrameNumDecodedRows() const {
  uint32_t num_decoded_rows = 0;
  if (stage >= State::Stage::WAITING_FOR_TILE_DATA &&
      current_frame_index < frames.size()) {
    const State::InternalFrameFeatures& frame = frames[current_frame_index];

    assert(tiles_layout.num_tiles_x > 0);
    const uint32_t num_decoded_tiles_y =
        tiles_layout.num_decoded_tiles / tiles_layout.num_tiles_x;

    if (num_decoded_tiles_y < tiles_layout.num_tiles_y) {
      num_decoded_rows = num_decoded_tiles_y * tiles_layout.tile_height;

      const uint32_t partial_tile_index = tiles_layout.num_decoded_tiles;
      assert(partial_tile_index < tiles_layout.tiles.size());
      if (((partial_tile_index + 1) % tiles_layout.num_tiles_x) == 0) {
        // If the partial tile is the right-most one of a row, it's sure that
        // all other tiles in that row are fully decoded, and thus the number
        // of decoded pixel rows is increased by as many as the partial tile.
        const Tile& partial_tile = tiles_layout.tiles[partial_tile_index];
        num_decoded_rows += partial_tile.num_decoded_rows;
      }
    } else {
      num_decoded_rows = frame.info.window.height;
    }
  }
  return num_decoded_rows;
}

//------------------------------------------------------------------------------

}  // namespace WP2
