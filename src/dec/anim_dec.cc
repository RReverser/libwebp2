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
// Functions related to animation decoding.
//
// Author: Yannis Guyon (yguyon@google.com)

#include <cassert>

#include "src/common/color_precision.h"
#include "src/dec/wp2_dec_i.h"
#include "src/utils/ans_utils.h"
#include "src/utils/data_source.h"
#include "src/wp2/base.h"
#include "src/wp2/decode.h"
#include "src/wp2/format_constants.h"

namespace WP2 {

//------------------------------------------------------------------------------

WP2Status DecodeANMF(DataSource* const data_source,
                     const BitstreamFeatures& features, uint32_t frame_index,
                     AnimationFrame* const frame) {
  if (features.is_animation) {
    const uint8_t* data;
    WP2_CHECK_OK(data_source->TryReadNext(kANMFHeaderSize, &data),
                 WP2_STATUS_NOT_ENOUGH_DATA);

    HeaderDec dec(data, kANMFHeaderSize);
    const uint32_t tag = dec.ReadBits(kANMFTagNumBits, "anmf_tag");
    frame->dispose = (tag == kANMFTagDispose);
    frame->blend = (dec.ReadBits(1, "blend") == 1);
    frame->is_last = (dec.ReadBits(1, "is_last") == 1);
    frame->duration_ms = dec.ReadBits(kFrameDurationNumBits, "duration");
    frame->window.x = dec.ReadBits(kImageDimNumBits, "frame_x");
    frame->window.y = dec.ReadBits(kImageDimNumBits, "frame_y");
    frame->window.width = dec.ReadBits(kImageDimNumBits, "frame_w") + 1;
    frame->window.height = dec.ReadBits(kImageDimNumBits, "frame_h") + 1;

    WP2_CHECK_OK((tag == kANMFTagDispose) || (tag == kANMFTagFrame),
                 WP2_STATUS_BITSTREAM_ERROR);
    WP2_CHECK_OK(frame->window.x + frame->window.width <= features.raw_width,
                 WP2_STATUS_BITSTREAM_ERROR);
    WP2_CHECK_OK(frame->window.y + frame->window.height <= features.raw_height,
                 WP2_STATUS_BITSTREAM_ERROR);
    if (frame_index == 0) {  // The first frame MUST dispose.
      WP2_CHECK_OK(frame->dispose, WP2_STATUS_BITSTREAM_ERROR);
    } else {
      WP2_CHECK_OK(frame_index < kMaxNumFrames, WP2_STATUS_BITSTREAM_ERROR);
      if (frame->duration_ms == 0) {
        WP2_CHECK_OK(!frame->is_last, WP2_STATUS_BITSTREAM_ERROR);
      }
    }
    // TODO(yguyon): Enforce these checks by design
  } else {
    frame->dispose = true;
    frame->blend = false;
    frame->duration_ms = kMaxFrameDurationMs;
    frame->window.x = 0;
    frame->window.y = 0;
    frame->window.width = features.raw_width;
    frame->window.height = features.raw_height;
    frame->is_last = true;
  }
  return WP2_STATUS_OK;
}

//------------------------------------------------------------------------------

// Fills the area outside of 'frame.window' with 'features.background_color' if
// 'frame.dispose' is true. The destination is the non-null buffer among
// 'rgb_output' and 'yuv_output'.
static WP2Status FillBorders(const BitstreamFeatures& features,
                             const AnimationFrame& frame,
                             ArgbBuffer* const rgb_output,
                             const CSPTransform* const csp_transform,
                             YUVPlane* const yuv_output) {
  if (!frame.dispose) return WP2_STATUS_OK;

  WP2_CHECK_OK((rgb_output != nullptr) != (yuv_output != nullptr),
               WP2_STATUS_INVALID_PARAMETER);
  const uint32_t width =
      (rgb_output != nullptr) ? rgb_output->width : yuv_output->GetWidth();
  const uint32_t height =
      (rgb_output != nullptr) ? rgb_output->height : yuv_output->GetHeight();
  const Rectangle& window = frame.window;
  assert(window.x + window.width <= width);
  assert(window.y + window.height <= height);

  const Rectangle top(0, 0, width, window.y);
  const Rectangle bot(0, window.y + window.height, width,
                      height - (window.y + window.height));
  const Rectangle lft(0, window.y, window.x, window.height);
  const Rectangle rgt(window.x + window.width, window.y,
                      width - (window.x + window.width), window.height);

  if (rgb_output != nullptr) {
    if (WP2FormatBpc(rgb_output->format) == 1) {
      const Argb32b color = ToArgb32b(features.background_color);

      rgb_output->Fill(top, color);
      rgb_output->Fill(bot, color);
      rgb_output->Fill(lft, color);
      rgb_output->Fill(rgt, color);
    } else {
      rgb_output->Fill(top, features.background_color);
      rgb_output->Fill(bot, features.background_color);
      rgb_output->Fill(lft, features.background_color);
      rgb_output->Fill(rgt, features.background_color);
    }
  } else {
    assert(csp_transform != nullptr);
    const Ayuv38b color =
        csp_transform->ToYUV(ToArgb32b(features.background_color));

    yuv_output->Fill(top, color);
    yuv_output->Fill(bot, color);
    yuv_output->Fill(lft, color);
    yuv_output->Fill(rgt, color);
  }
  return WP2_STATUS_OK;
}

WP2Status FillBorders(const BitstreamFeatures& features,
                      const AnimationFrame& frame, ArgbBuffer* const output) {
  return FillBorders(features, frame, output,
                     /*csp_transform=*/nullptr, /*yuv_output=*/nullptr);
}

WP2Status FillBorders(const BitstreamFeatures& features,
                      const AnimationFrame& frame,
                      const CSPTransform& csp_transform,
                      YUVPlane* const output) {
  return FillBorders(features, frame, /*rgb_output=*/nullptr,
                     /*csp_transform=*/&csp_transform, /*yuv_output=*/output);
}

//------------------------------------------------------------------------------

}  // namespace WP2
