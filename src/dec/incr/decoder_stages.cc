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

#include "src/common/color_precision.h"
#include "src/dec/incr/decoder_info.h"
#include "src/dec/incr/decoder_skip.h"
#include "src/dec/incr/decoder_state.h"
#include "src/dec/wp2_dec_i.h"
#include "src/utils/orientation.h"
#include "src/utils/vector.h"
#include "src/wp2/decode.h"

#if __cplusplus >= 201703L
#define WP2_FALLTHROUGH [[fallthrough]]
#else
#define WP2_FALLTHROUGH do { } while (0)
#endif

namespace WP2 {

//------------------------------------------------------------------------------

static WP2Status UpdateNumFramesToDiscard(const DecoderConfig& config,
                                          const BitstreamFeatures& features,
                                          Decoder::State* const state) {
  WP2Status status = state->ComputeNumFramesToDiscard();
  while (status == WP2_STATUS_NOT_ENOUGH_DATA) {
    WP2_CHECK_STATUS(GlimpseAtNextFrame(config, features, state));
    status = state->ComputeNumFramesToDiscard();
  }
  return status;
}

static WP2Status PrepareTiles(const DecoderConfig& config,
                              const BitstreamFeatures& features,
                              Decoder::State* const state) {
  WP2_CHECK_STATUS(
      DecodeGLBL(state->data_source, config, features, &state->gparams));
  FillDecoderInfo(state->gparams, config);

  const bool discard_frame =
      (state->current_frame_index < state->discard_frames_before_index);
  const Rectangle& window =
      state->frames[state->current_frame_index].info.window;

  if (!discard_frame) {
    WP2_CHECK_STATUS(
        state->frame_output_buffer.SetView(state->raw_output_buffer, window));
    if (state->gparams.type_ != GlobalParams::GP_LOSSLESS) {
      WP2_CHECK_STATUS(state->yuv_output_buffer.Resize(
          window.width, window.height, /*pad=*/kPredWidth,
          state->gparams.has_alpha_));
    } else {
      state->yuv_output_buffer.Clear();
    }
  }

  WP2_CHECK_STATUS(
      GetTilesLayout(window.width, window.height, features.tile_width,
                     features.tile_height, &state->frame_output_buffer,
                     &state->yuv_output_buffer, &state->tiles_layout));
  state->tiles_layout.gparams = &state->gparams;

  if (!discard_frame) {
    if (state->gparams.type_ != GlobalParams::GP_LOSSLESS) {
      state->intertile_filter.Init(config, features, state->gparams,
                                   state->tiles_layout,
                                   &state->yuv_output_buffer);
    }
    // TODO(maryla): support multiple frames (add the bpps?). This will
    // only contain information for the last frame.
    if (config.info != nullptr && config.info->bits_per_pixel != nullptr) {
      WP2_CHECK_STATUS(
          config.info->bits_per_pixel->Resize(window.width, window.height));
    }

    WP2_CHECK_STATUS(SetupWorkers(features, config, /*progress=*/nullptr,
                                  &state->tiles_layout, &state->workers));
  }
  return WP2_STATUS_OK;
}

static WP2Status SetupOutputBuffer(const AnimationFrame& frame_info,
                                   ArgbBuffer* const output,
                                   Decoder::State* const state) {
  /*
   * BUFFER USE.
   * If 'must_use_raw_output_buffer' is set, we cannot render to 'output'
   * directly, we must use an intermediate buffer. This happens if
   * 'decoding_orientation' is set, or the 'output' buffer is slow, or the pixel
   * format differs. 'must_use_raw_output_buffer' is global to the animation.
   * Moreover, individual frames can require a temporary buffer if their 'blend'
   * bit is set.
   *
   * There are four cases. For simplicity, we'll consider that
   * 'must_use_raw_output_buffer' is set because of rotation, but it could be
   * because of is_slow or pixel format as well.
   *
   * - No rotation, no blending: we render to 'raw_output_buffer' which is set
   *   up as a view to 'output', so we basically render to output directly.
   *   using_raw_output_buffer == false
   *   using_blended_buffer == false
   *
   * - No rotation, with blending: we allocate 'raw_output_buffer' and render to
   *   it, then blend onto 'output'.
   *   using_raw_output_buffer == true
   *   using_blended_buffer == false
   *
   * - With rotation, no blending: we allocate 'raw_output_buffer' and render to
   *   it, then copy with rotation to 'output'
   *   using_raw_output_buffer == true
   *   using_blended_buffer == false
   *
   * - With rotation, with blending: we allocate 'raw_output_buffer' and
   *   'blended_buffer', we render on 'raw_output_buffer', then blend onto
   *   'blended_buffer', then copy with rotation to 'output'.
   *   using_raw_output_buffer == true
   *   using_blended_buffer == true
   */
  if (frame_info.blend) {
    if (state->using_blended_buffer || (state->using_raw_output_buffer &&
                                        !state->must_use_raw_output_buffer)) {
      // If we're already using the blended buffer, or if we're using the
      // raw_output_buffer even if we don't necessarily have to, it means the
      // previous frame was also blending and buffers are already correctly set
      // up.
    } else if (state->using_raw_output_buffer) {
      // raw_output_buffer is already being used (for example because
      // a rotation is set). We need an extra buffer for blending.
      assert(!state->raw_output_buffer.IsView());
      // Swap so blended_buffer contains the pixels that we will blend onto.
      WP2_CHECK_STATUS(state->blended_buffer.Swap(&state->raw_output_buffer));
      if (state->raw_output_buffer.IsEmpty()) {
        WP2_CHECK_STATUS(state->raw_output_buffer.Resize(
            state->blended_buffer.width, state->blended_buffer.height));
      }
      state->using_blended_buffer = true;
    } else {
      assert(!state->using_blended_buffer);
      // Use raw_output_buffer. It will later be blended onto 'output' directly.
      if (state->raw_output_buffer.IsView()) {
        state->raw_output_buffer.Deallocate();
      }
      WP2_CHECK_STATUS(
          state->raw_output_buffer.Resize(output->width, output->height));
      state->using_raw_output_buffer = true;
    }
  } else {
    if (state->must_use_raw_output_buffer) {
      assert(state->using_raw_output_buffer);
      if (state->using_blended_buffer) {
        // The previous frame was being blended. We don't need the extra
        // 'blended_buffer' anymore but we have to copy the result.
        WP2_CHECK_STATUS(state->raw_output_buffer.Swap(&state->blended_buffer));
        state->using_blended_buffer = false;
      }
    } else {
      assert(!state->using_blended_buffer);
      state->status = state->raw_output_buffer.SetView(*output);
      state->using_raw_output_buffer = false;
    }
  }
  return WP2_STATUS_OK;
}

WP2Status Decoder::BlendFrame(Rectangle window, bool dispose) {
  ArgbBuffer* const bg =
      state_->using_blended_buffer ? &state_->blended_buffer : output_;
  if (dispose) {
    bg->Fill(ToArgb32b(features_.background_color));
  }
  WP2_CHECK_STATUS(bg->CompositeUnder(state_->raw_output_buffer, window));
  return WP2_STATUS_OK;
}

//------------------------------------------------------------------------------

void Decoder::DecodeAvailableData() {
  if (Failed()) return;

  WP2DecDspInit();
  assert(state_->data_source != nullptr);

  Rectangle previous_decoded_area = GetDecodedArea();

  bool keep_decoding;
  do {
    keep_decoding = false;
    if (state_->current_frame_is_ready) {  // If a frame was just fully decoded:
      state_->current_frame_is_ready = false;
      if (state_->frames[state_->current_frame_index].info.is_last) {
        state_->stage = State::Stage::WAITING_FOR_METADATA_CHUNKS;
      } else {
        state_->stage = State::Stage::WAITING_FOR_ANMF_CHUNK;
        ++state_->current_frame_index;
        previous_decoded_area = {0, 0, 0, 0};
        state_->frame_num_converted_rows = 0;
        state_->frame_num_final_rows = 0;
        state_->intertile_filter.Clear();
      }
    }

    switch (state_->stage) {
      case State::Stage::WAITING_FOR_HEADER_CHUNK: {
        // TODO(yguyon): If Rewind(), verify it's the same bitstream/header.
        state_->status = DecodeHeader(state_->data_source, &features_);
        if (state_->status == WP2_STATUS_OK) {
          state_->status = SetupDecoderInfo(features_, config_);
        }
        if (!state_->NextStage()) break;
        WP2_FALLTHROUGH;
      }
      // TODO(maryla): merge with SetupOutputBuffer?
      case State::Stage::ALLOCATE_OUTPUT_BUFFER: {
        if (!state_->ShouldDiscardAllFrames()) {
          if (state_->index_of_frame_in_output_buffer < kMaxNumFrames &&
              (output_->IsEmpty() || output_->width != features_.width ||
               output_->height != features_.height)) {
            // Either Rewind() was given a false negative parameter
            // 'output_buffer_changed', or the bitstream changed.
            state_->status = WP2_STATUS_BAD_READ;
          } else {
            state_->status = output_->Resize(features_.width, features_.height);
          }
          if (state_->status == WP2_STATUS_OK) {
            state_->must_use_raw_output_buffer =
                (features_.orientation != Orientation::kOriginal ||
                 output_->format != state_->raw_output_buffer.format ||
                 output_->IsSlow());
            if (state_->must_use_raw_output_buffer) {
              // Temporary buffer to store the raw (non-oriented Argb32) output.
              state_->status = state_->raw_output_buffer.Resize(
                  features_.raw_width, features_.raw_height);
              state_->using_raw_output_buffer = true;
            }
          }
        }
        if (!state_->NextStage()) break;
        WP2_FALLTHROUGH;
      }
      case State::Stage::WAITING_FOR_PREVIEW_CHUNK: {
        state_->status = SkipPreview(state_->data_source, features_);
        if (!state_->NextStage()) break;
        WP2_FALLTHROUGH;
      }
      case State::Stage::WAITING_FOR_ICC_CHUNK: {
        state_->status =
            DecodeICC(state_->data_source, features_, &output_->metadata.iccp);
        if (!state_->NextStage()) break;
        WP2_FALLTHROUGH;
      }
      case State::Stage::WAITING_FOR_ANMF_CHUNK: {
        state_->status = UpdateNumFramesToDiscard(config_, features_, state_);
        if (state_->status == WP2_STATUS_OK) {
          if (DiscardTillMetadata(state_)) {
            state_->stage = State::Stage::WAITING_FOR_METADATA_CHUNKS;
            keep_decoding = true;
            break;  // Go to WAITING_FOR_METADATA_CHUNKS part.
          }
          if (FrameIsAlreadyInOutputBuffer(features_.tile_width,
                                           features_.tile_height, state_)) {
            state_->stage = State::Stage::WAITING_FOR_TILE_DATA;
            keep_decoding = true;
            break;  // Go to WAITING_FOR_TILE_DATA part.
          }
          if (DiscardWholeFrames(state_)) {
            // Maybe the output buffer contains 'current_frame_index' now.
            if (FrameIsAlreadyInOutputBuffer(features_.tile_width,
                                             features_.tile_height, state_)) {
              state_->stage = State::Stage::WAITING_FOR_TILE_DATA;
              keep_decoding = true;
              break;  // Go to WAITING_FOR_TILE_DATA part.
            }
          }
        }

        if (state_->status == WP2_STATUS_OK) {
          state_->status = DecodeFrameInfo(state_->data_source, features_,
                                           /*glimpse=*/false, state_);
          if (state_->status == WP2_STATUS_OK &&
              state_->current_frame_index >=
                  state_->discard_frames_before_index) {
            // Make sure the (supposed) content of the output buffer is correct.
            const AnimationFrame& frame_info =
                state_->frames[state_->current_frame_index].info;
            assert(frame_info.dispose ||
                   (state_->index_of_frame_in_output_buffer + 1 ==
                    state_->current_frame_index));
            // The content of the output buffer is no longer an entire frame.
            state_->index_of_frame_in_output_buffer = kMaxNumFrames;

            state_->status = SetupOutputBuffer(frame_info, output_, state_);
            if (state_->status == WP2_STATUS_OK) {
              state_->status = FillBorders(
                  features_, state_->frames[state_->current_frame_index].info,
                  &state_->raw_output_buffer);
            }
          }
        }
        if (!state_->NextStage()) break;
        WP2_FALLTHROUGH;
      }
      case State::Stage::PREPARE_TILES: {
        state_->tiles_layout.gparams = nullptr;   // for sanity
        state_->status = PrepareTiles(config_, features_, state_);
        if (!state_->NextStage()) break;
        WP2_FALLTHROUGH;
      }
      case State::Stage::WAITING_FOR_TILE_DATA: {
        state_->status = UpdateNumFramesToDiscard(config_, features_, state_);
        if (state_->status == WP2_STATUS_OK) {
          PartialTileDecoder* const partial_tile_decoder =
              state_->partial_tile_decoder.get();
          if (DiscardTillMetadata(state_)) {
            if (partial_tile_decoder != nullptr) partial_tile_decoder->Clear();
            state_->stage = State::Stage::WAITING_FOR_METADATA_CHUNKS;
            keep_decoding = true;
            break;
          }
          if (DiscardWholeFrames(state_)) {
            if (partial_tile_decoder != nullptr) partial_tile_decoder->Clear();
            state_->stage = State::Stage::WAITING_FOR_ANMF_CHUNK;
            keep_decoding = true;
            break;
          }
        }

        if (state_->status == WP2_STATUS_OK) {
          if (state_->current_frame_index <
              state_->discard_frames_before_index) {
            state_->status = DiscardAvailableTiles();
          } else {
            state_->status = DecodeAvailableTiles();
          }
        }

        if (state_->status == WP2_STATUS_OK) {
          state_->status =
              SaveFrameNumBytes(state_->data_source->GetNumDiscardedBytes() +
                                    state_->data_source->GetNumReadBytes(),
                                &state_->frames[state_->current_frame_index]);
        }

        if (state_->status == WP2_STATUS_OK) {
          const uint32_t current_frame_index = state_->current_frame_index;
          const State::InternalFrameFeatures& current_frame_features =
              state_->frames[current_frame_index];
          if (current_frame_features.info.duration_ms > 0) {  // Regular frame.
            state_->num_available_final_frames =
                std::max(state_->num_available_final_frames,
                         current_frame_features.final_frame_index + 1u);
          }

          if (current_frame_index < state_->discard_frames_before_index) {
            keep_decoding = true;
          } else {
            state_->index_of_frame_in_output_buffer = current_frame_index;
            const bool is_preframe =
                current_frame_features.info.duration_ms == 0;
            // Blend the preframe if necessary.
            if (is_preframe && current_frame_features.info.blend) {
              const Rectangle window = RotateRectangle(
                  GetInverseOrientation(features_.orientation), features_.width,
                  features_.height, current_frame_features.window);
              state_->status =
                  BlendFrame(window, current_frame_features.info.dispose);
            }
            // Merge all consecutive preframes with the next regular frame,
            // forming a final frame. Don't notify for preframe readiness.
            keep_decoding = is_preframe;
          }
          if (state_->status == WP2_STATUS_OK) {
            state_->current_frame_is_ready = true;
            state_->status = WP2_STATUS_NOT_ENOUGH_DATA;
            break;  // If !keep_decoding, halts the decoding so that
                    // the current frame can be displayed or copied.
          }
        }
        if (!state_->NextStage()) break;
        WP2_FALLTHROUGH;
      }
      case State::Stage::WAITING_FOR_METADATA_CHUNKS: {
        state_->status =
            DecodeMetadata(state_->data_source, features_,
                           &output_->metadata.exif, &output_->metadata.xmp);
        if (!state_->NextStage()) break;
        WP2_FALLTHROUGH;
      }
      case State::Stage::DECODED:
      default:
        break;
    }
  } while (keep_decoding && !Failed());

  // Blend/rotate/copy to output, if needed.
  if (TryGetDecodedFeatures() != nullptr && state_->using_raw_output_buffer &&
      !output_->IsEmpty()) {
    const Rectangle previous_raw_decoded_area = RotateRectangle(
        GetInverseOrientation(features_.orientation),
        features_.width, features_.height, previous_decoded_area);
    const Rectangle raw_decoded_area = RotateRectangle(
        GetInverseOrientation(features_.orientation),
        features_.width, features_.height, GetDecodedArea());
    assert(raw_decoded_area.x == 0 && raw_decoded_area.y == 0);
    assert(raw_decoded_area.height >= previous_raw_decoded_area.height);
    assert(previous_raw_decoded_area.width == 0 ||
           (raw_decoded_area.width == previous_raw_decoded_area.width));

    const Rectangle new_raw_decoded_area =
        Rectangle(raw_decoded_area.x,
                  raw_decoded_area.y + previous_raw_decoded_area.height,
                  raw_decoded_area.width,
                  raw_decoded_area.height - previous_raw_decoded_area.height);

    if (new_raw_decoded_area.GetArea() > 0) {
      WP2Status postprocess_status = WP2_STATUS_OK;
      assert(state_->current_frame_index < state_->frames.size());
      const AnimationFrame& frame =
          state_->frames[state_->current_frame_index].info;
      if (frame.blend) {
        const Rectangle clipped = new_raw_decoded_area.ClipWith(frame.window);
        if (clipped.GetArea() > 0) {
          postprocess_status = BlendFrame(clipped, frame.dispose);
        }
      }

      if (postprocess_status == WP2_STATUS_OK &&
          state_->must_use_raw_output_buffer) {
        // Also handles the conversion from Argb32 to something else, if needed.
        postprocess_status = RotateSubBuffer(features_.orientation,
                                             state_->using_blended_buffer
                                                 ? state_->blended_buffer
                                                 : state_->raw_output_buffer,
                                             new_raw_decoded_area, output_);
      }
      if (postprocess_status != WP2_STATUS_OK) {
        state_->status = postprocess_status;
      }
    }
  }
}

//------------------------------------------------------------------------------

}  // namespace WP2
