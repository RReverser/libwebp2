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
// Frame header (ANMF chunk) decoding and verification.
//
// Author: Yannis Guyon (yguyon@google.com)

#include "src/dec/incr/decoder_info.h"

#include <algorithm>

#include "src/dec/incr/decoder_state.h"
#include "src/dec/wp2_dec_i.h"
#include "src/utils/orientation.h"
#include "src/utils/utils.h"
#include "src/wp2/decode.h"

namespace WP2 {

//------------------------------------------------------------------------------

static Rectangle BoundingRectangle(const Rectangle& a, const Rectangle& b) {
  const uint32_t left_x = std::min(a.x, b.x);
  const uint32_t top_y = std::min(a.y, b.y);
  const uint32_t right_x = std::max(a.x + a.width, b.x + b.width);
  const uint32_t bottom_y = std::max(a.y + a.height, b.y + b.height);
  return {left_x, top_y, right_x - left_x, bottom_y - top_y};
}

//------------------------------------------------------------------------------

WP2Status DecodeFrameInfo(DataSource* const data_source,
                          const BitstreamFeatures& features, bool glimpse,
                          Decoder::State* const state) {
  const uint32_t frame_index =
      glimpse ? state->glimpsed_frame_index : state->current_frame_index;
  assert(frame_index <= state->frames.size());
  const size_t bitstream_position =
      data_source->GetNumDiscardedBytes() + data_source->GetNumReadBytes();
  AnimationFrame info;
  // Fills 'info' with default values if it's not an animation.
  WP2_CHECK_STATUS(DecodeANMF(data_source, features, frame_index, &info));

  Decoder::State::InternalFrameFeatures frame;
  frame.bitstream_position = bitstream_position;
  frame.num_bytes = Decoder::State::kUnknown;
  frame.info = info;
  frame.duration_ms = frame.info.duration_ms;
  frame.window = RotateRectangle(features.orientation, features.raw_width,
                                 features.raw_height, frame.info.window);
  frame.is_last = frame.info.is_last;
  frame.final_frame_index = state->GetFinalFrameIndex(frame_index);

  const bool is_regular_frame = (frame.info.duration_ms > 0);
  // 'last_dispose_frame_index' is actually a final frame index.
  if (frame.info.dispose) {
    frame.last_dispose_frame_index = frame.final_frame_index;
  } else {
    assert(frame_index > 0);  // Otherwise it would be 'frame.info.dispose'.
    frame.last_dispose_frame_index =
        state->frames[frame_index - 1].last_dispose_frame_index;
    if (!is_regular_frame) {
      const uint32_t previous_regular_frame_index =
          (frame.final_frame_index == 0)
              ? 0
              : state->final_to_regular_index[frame.final_frame_index - 1];
      WP2_CHECK_OK(
          (frame_index - previous_regular_frame_index) <= kMaxNumPreframes,
          WP2_STATUS_BITSTREAM_ERROR);
    }
  }

  if (!frame.info.dispose) {
    uint32_t previous_dispose = frame_index;
    while (previous_dispose > 0 &&
           !state->frames[--previous_dispose].info.dispose) {
    }
    WP2_CHECK_OK((frame.final_frame_index -
                  state->frames[previous_dispose].final_frame_index) <=
                     kMaxNumDependentFrames + 1,
                 WP2_STATUS_BITSTREAM_ERROR);
  }

  if (is_regular_frame) {
    // Merge previous preframe windows with the current regular frame window.
    for (uint32_t i = frame_index; i > 0; --i) {
      if (state->frames[i - 1].info.duration_ms > 0) break;
      frame.window =
          BoundingRectangle(frame.window, state->frames[i - 1].window);
    }
  }

  if (frame_index > 0) {
    assert(frame_index - 1 < state->frames.size());
    const Decoder::State::InternalFrameFeatures& previous_frame =
        state->frames[frame_index - 1];
    assert(previous_frame.num_bytes != Decoder::State::kUnknown);
    WP2_CHECK_OK(frame.bitstream_position == previous_frame.bitstream_position +
                                                 previous_frame.num_bytes,
                 WP2_STATUS_BAD_READ);
  }

  if (frame_index < state->frames.size()) {
    // Features are already known, make sure the bitstream did not change.
    frame.num_bytes = state->frames[frame_index].num_bytes;  // Not known yet.
    WP2_CHECK_OK(frame == state->frames[frame_index], WP2_STATUS_BAD_READ);
  } else {
    WP2_CHECK_ALLOC_OK(state->frames.push_back(frame));
    if (is_regular_frame) {
      assert(frame.final_frame_index == state->final_to_regular_index.size());
      WP2_CHECK_ALLOC_OK(state->final_to_regular_index.push_back(frame_index));
    }
  }
  return WP2_STATUS_OK;
}

//------------------------------------------------------------------------------

WP2Status SaveFrameNumBytes(
    size_t bitstream_position_at_end_of_frame,
    Decoder::State::InternalFrameFeatures* const frame) {
  WP2_CHECK_OK(bitstream_position_at_end_of_frame > frame->bitstream_position,
               WP2_STATUS_BAD_READ);
  const size_t num_bytes =
      bitstream_position_at_end_of_frame - frame->bitstream_position;
  if (frame->num_bytes != Decoder::State::kUnknown) {
    WP2_CHECK_OK(frame->num_bytes == num_bytes, WP2_STATUS_BAD_READ);
  }
  frame->num_bytes = num_bytes;
  return WP2_STATUS_OK;
}

//------------------------------------------------------------------------------

}  // namespace WP2
