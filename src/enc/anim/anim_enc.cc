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
//  Animation encoding functions.
//
// Author: Skal (pascal.massimino@gmail.com)

#include "./anim_enc.h"

#include "src/enc/preview/preview_enc.h"
#include "src/enc/wp2_enc_i.h"
#include "src/utils/utils.h"
#include "src/wp2/encode.h"
#include "src/wp2/format_constants.h"

namespace WP2 {

//------------------------------------------------------------------------------
// AnimationEncoder

// Pointer to implementation for AnimationEncoder state.
struct AnimationEncoder::State : public WP2Allocable {
  Vector<Frame> frames;
};

// If PImpl allocation failed, AddFrame() and Encode() will return OUT_OF_MEM.
AnimationEncoder::AnimationEncoder()
    : state_(new (WP2Allocable::nothrow) State()) {}

AnimationEncoder::AnimationEncoder(AnimationEncoder&& other) noexcept
    : state_(other.state_) {
  other.state_ = nullptr;
}

AnimationEncoder::~AnimationEncoder() { delete state_; }

WP2Status AnimationEncoder::AddFrame(const ArgbBuffer& pixels,
                                     uint32_t duration_ms, bool force_dispose,
                                     PictureHint picture_hint) {
  (void)picture_hint;
  WP2_CHECK_ALLOC_OK(state_ != nullptr);
  Vector<Frame>& frames = state_->frames;

  WP2_CHECK_OK(!pixels.IsEmpty(), WP2_STATUS_NULL_PARAMETER);
  WP2_CHECK_OK(duration_ms > 0 && duration_ms <= kMaxFrameDurationMs,
               WP2_STATUS_INVALID_PARAMETER);

  // Checks dimension validity only for the first frame. For the following ones,
  // checks dimension equality to the first frame.
  if (frames.empty()) {
    WP2_CHECK_OK(pixels.format == WP2_Argb_32, WP2_STATUS_INVALID_COLORSPACE);
    WP2_CHECK_OK((pixels.width > 0) && (pixels.height > 0) &&
                     (pixels.width <= kImageDimMax) &&
                     (pixels.height <= kImageDimMax),
                 WP2_STATUS_BAD_DIMENSION);
  } else {
    WP2_CHECK_OK(frames.size() < kMaxNumFrames, WP2_STATUS_INVALID_PARAMETER);
    WP2_CHECK_OK(pixels.width == frames.front().GetWidth() &&
                     pixels.height == frames.front().GetHeight(),
                 WP2_STATUS_BAD_DIMENSION);

    // A mix of color spaces is currently unsupported.
    WP2_CHECK_OK(frames.front().ccsp_pixels.IsEmpty(),
                 WP2_STATUS_UNSUPPORTED_FEATURE);
  }

  WP2_CHECK_ALLOC_OK(frames.resize(frames.size() + 1));
  WP2Status status = frames.back().rgb_pixels.CopyFrom(pixels);
  if (status != WP2_STATUS_OK) frames.pop_back();  // No-op if problem.
  WP2_CHECK_STATUS(status);
  frames.back().duration_ms = duration_ms;
  frames.back().force_dispose = force_dispose;
  return WP2_STATUS_OK;
}

WP2Status AnimationEncoder::AddFrame(
    uint32_t width, uint32_t height,
    const int16_t* c0_buffer, uint32_t c0_step,
    const int16_t* c1_buffer, uint32_t c1_step,
    const int16_t* c2_buffer, uint32_t c2_step,
    const int16_t* a_buffer, uint32_t a_step,
    const int16_t ccsp_to_rgb_matrix[9],
    uint32_t ccsp_to_rgb_shift, uint32_t duration_ms, bool force_dispose) {
  WP2_CHECK_ALLOC_OK(state_ != nullptr);
  WP2_CHECK_OK(ccsp_to_rgb_matrix != nullptr, WP2_STATUS_NULL_PARAMETER);
  const CSPMtx ccsp_to_rgb(ccsp_to_rgb_matrix, ccsp_to_rgb_shift);
  Vector<Frame>& frames = state_->frames;

  const bool has_alpha = (a_buffer != nullptr);

  const int16_t* const buffers[] = {c0_buffer, c1_buffer, c2_buffer, a_buffer};
  uint32_t steps[] = {c0_step, c1_step, c2_step, a_step};
  for (Channel channel : {kYChannel, kUChannel, kVChannel, kAChannel}) {
    if (channel == kAChannel && !has_alpha) continue;
    WP2_CHECK_OK(buffers[channel] != nullptr, WP2_STATUS_NULL_PARAMETER);
    WP2_CHECK_OK(steps[channel] >= width, WP2_STATUS_BAD_DIMENSION);
  }

  WP2_CHECK_OK(duration_ms > 0 && duration_ms <= kMaxFrameDurationMs,
               WP2_STATUS_INVALID_PARAMETER);

  // Checks dimension validity only for the first frame. For the following ones,
  // checks dimension equality to the first frame.
  if (frames.empty()) {
    WP2_CHECK_OK((width > 0 && height > 0) &&
                     (width <= kImageDimMax && height <= kImageDimMax),
                 WP2_STATUS_BAD_DIMENSION);
    WP2_CHECK_OK(ccsp_to_rgb.shift <= 16, WP2_STATUS_INVALID_PARAMETER);
  } else {
    const Frame& first_frame = frames.front();
    WP2_CHECK_OK(frames.size() < kMaxNumFrames, WP2_STATUS_INVALID_PARAMETER);
    WP2_CHECK_OK(
        width == first_frame.GetWidth() && height == first_frame.GetHeight(),
        WP2_STATUS_BAD_DIMENSION);

    // A mix of color spaces is currently unsupported.
    WP2_CHECK_OK(first_frame.rgb_pixels.IsEmpty(),
                 WP2_STATUS_UNSUPPORTED_FEATURE);
    WP2_CHECK_OK(ccsp_to_rgb == first_frame.ccsp_to_rgb,
                 WP2_STATUS_UNSUPPORTED_FEATURE);
  }

  WP2_CHECK_ALLOC_OK(frames.resize(frames.size() + 1));
  WP2Status success = WP2_STATUS_OK;
  for (Channel channel : {kYChannel, kUChannel, kVChannel, kAChannel}) {
    if (channel == kAChannel && !has_alpha) continue;
    Plane16* const plane = &frames.back().ccsp_pixels.GetChannel(channel);
    success = plane->Resize(width, height);
    if (success != WP2_STATUS_OK) break;
    const int16_t* from_row = buffers[channel];
    for (uint32_t y = 0; y < plane->h_; ++y) {
      std::copy(from_row, from_row + width, plane->Row(y));
      from_row += steps[channel];
    }
  }
  if (success != WP2_STATUS_OK) frames.pop_back();  // No-op if problem.
  WP2_CHECK_STATUS(success);
  frames.back().ccsp_to_rgb = ccsp_to_rgb;
  frames.back().duration_ms = duration_ms;
  frames.back().force_dispose = force_dispose;
  return WP2_STATUS_OK;
}

WP2Status AnimationEncoder::Encode(Writer* output, const EncoderConfig& config,
                                   uint8_t loop_count,
                                   const Metadata& metadata) const {
  WP2_CHECK_ALLOC_OK(state_ != nullptr);
  const Vector<Frame>& frames = state_->frames;

  WP2_CHECK_OK(!frames.empty(), WP2_STATUS_INVALID_PARAMETER);
  WP2_CHECK_OK(output != nullptr, WP2_STATUS_NULL_PARAMETER);
  WP2_CHECK_OK(config.IsValid(), WP2_STATUS_INVALID_CONFIGURATION);
  WP2_CHECK_OK(loop_count <= kMaxLoopCount, WP2_STATUS_INVALID_PARAMETER);

  const bool has_alpha = !AllFramesAreOpaque(frames);
  const Frame& first_frame = frames.front();
  const RGB12b preview_color =
      first_frame.ccsp_pixels.IsEmpty()
          ? GetPreviewColor(first_frame.rgb_pixels)
          : GetPreviewColor(first_frame.ccsp_pixels, first_frame.ccsp_to_rgb);
  const bool has_icc = (metadata.iccp.size > 0);
  const bool has_trailing_data =
      (metadata.xmp.size > 0) || (metadata.exif.size > 0);

  // TODO(yguyon): Actually encode frames and compute background color etc.
  WP2_CHECK_STATUS(EncodeHeader(
      config, first_frame.GetWidth(), first_frame.GetHeight(), has_alpha,
      /*is_anim=*/true, loop_count, /*background_color=*/kTransparentArgb38b,
      preview_color, has_icc, has_trailing_data, output));

  if (config.create_preview && first_frame.rgb_pixels.IsEmpty()) {
    ArgbBuffer converted;
    // TODO(yguyon): Maybe also converted in FrameEncoder::EncodeSubframe()
    WP2_CHECK_STATUS(first_frame.ccsp_pixels.Export(
        first_frame.ccsp_to_rgb, /*resize_if_needed=*/true, &converted));
    WP2_CHECK_STATUS(EncodePreview(converted, config, output));
  } else {
    WP2_CHECK_STATUS(EncodePreview(first_frame.rgb_pixels, config, output));
  }
  WP2_CHECK_STATUS(
      EncodeICC({metadata.iccp.bytes, metadata.iccp.size}, output));

  {
    FrameEncoder frame_encoder(config, /*background_color=*/kTransparentArgb38b,
                               frames);
    WP2_CHECK_STATUS(frame_encoder.Encode(output));
  }

  WP2_CHECK_STATUS(EncodeMetadata(metadata, output));
  return WP2_STATUS_OK;
}

//------------------------------------------------------------------------------

}  // namespace WP2
