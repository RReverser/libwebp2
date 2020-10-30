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
// Misc. common utility functions
//
// Author: Skal (pascal.massimino@gmail.com)

#include "src/utils/plane.h"

#include <algorithm>
#include <numeric>

#include "src/dsp/dsp.h"
#include "src/utils/csp.h"

namespace WP2 {

//------------------------------------------------------------------------------

bool SamplingTaps::Taps::operator==(const SamplingTaps::Taps& rhs) const {
  return radius == rhs.radius &&
         std::equal(taps - radius, taps + radius + 1, rhs.taps - radius);
}

bool SamplingTaps::operator==(const SamplingTaps& rhs) const {
  return t1 == rhs.t1 && t2 == rhs.t2;
}

bool SamplingTaps::Taps::IsValid() const {
  return std::accumulate(taps - radius, taps + radius + 1, 0) == 16;
}

//------------------------------------------------------------------------------

// Lightly sharpen while downscaling.
// Centered on the top-left pixel for a 2x2 block.
static constexpr tap_t kDownSharpCoeffs[2 * 2 + 1] = {0, -2, 10, 10, -2};
const SamplingTaps SamplingTaps::kDownSharp = {
    SamplingTaps::Taps(2, kDownSharpCoeffs + 2),
    SamplingTaps::Taps(2, kDownSharpCoeffs + 2)};

// Just average a block of 2x2 pixels into a single pixel.
// Centered on the top-left pixel for a 2x2 block.
static constexpr tap_t kDownAvgCoeffs[2 * 1 + 1] = {0, 8, 8};
const SamplingTaps SamplingTaps::kDownAvg = {
    SamplingTaps::Taps(1, kDownAvgCoeffs + 1),
    SamplingTaps::Taps(1, kDownAvgCoeffs + 1)};

//------------------------------------------------------------------------------

// Slightly smooths while upscaling, should compensate for the sharp downscale.
static constexpr tap_t kUpSmoothFull[2 * 2 + 1] = {1, 3, 11, 1, 0};
static constexpr tap_t kUpSmoothHalf[2 * 2 + 1] = {0, 1, 11, 3, 1};
const SamplingTaps SamplingTaps::kUpSmooth = {
    SamplingTaps::Taps(2, kUpSmoothFull + 2),
    SamplingTaps::Taps(2, kUpSmoothHalf + 2)};

// Just scale by two with no interpolation.
static constexpr tap_t kUpNearestCoeffs[2 * 0 + 1] = {16};
const SamplingTaps SamplingTaps::kUpNearest = {
    SamplingTaps::Taps(0, kUpNearestCoeffs),
    SamplingTaps::Taps(0, kUpNearestCoeffs)};

//------------------------------------------------------------------------------

WP2Status PlaneBase::CheckRect(const Rectangle& rect) const {
  WP2_CHECK_OK(rect.GetArea() > 0, WP2_STATUS_BAD_DIMENSION);
  WP2_CHECK_OK(rect.x + rect.width <= w_, WP2_STATUS_BAD_DIMENSION);
  WP2_CHECK_OK(rect.y + rect.height <= h_, WP2_STATUS_BAD_DIMENSION);
  return WP2_STATUS_OK;
}

//------------------------------------------------------------------------------

template <>
WP2Status Plane16::ToGray(WP2::ArgbBuffer* dst, uint32_t num_bits,
                          bool shift) const {
  assert(num_bits >= 8 && num_bits <= 16);
  const uint32_t H = std::min(h_, dst->height);
  const uint32_t W = std::min(w_, dst->width);
  for (uint32_t y = 0; y < H; ++y) {
    uint8_t* const dst_y = (uint8_t*)dst->GetRow(y);
    const int16_t* const src_y = Row(y);
    for (uint32_t x = 0; x < W; ++x) {
      dst_y[4 * x + 0] = 0xff;
      dst_y[4 * x + 1] = dst_y[4 * x + 2] = dst_y[4 * x + 3] = (uint8_t)Clamp(
          (shift ? 128 : 0) + RightShiftRound(src_y[x], num_bits - 8), 0, 255);
    }
  }
  return WP2_STATUS_OK;
}

template <>
WP2Status Planef::ToGray(WP2::ArgbBuffer* dst, uint32_t num_bits,
                         bool shift) const {
  const uint32_t H = std::min(h_, dst->height);
  const uint32_t W = std::min(w_, dst->width);
  const float scale =
      (num_bits >= 8) ? 1.f / (1 << (num_bits - 8)) : (1 << (8 - num_bits));
  for (uint32_t y = 0; y < H; ++y) {
    uint8_t* const dst_y = (uint8_t*)dst->GetRow(y);
    const float* const src_y = Row(y);
    for (uint32_t x = 0; x < W; ++x) {
      dst_y[4 * x + 0] = 0xff;
      dst_y[4 * x + 1] = dst_y[4 * x + 2] = dst_y[4 * x + 3] = (uint8_t)Clamp(
          (shift ? 128.f : 0.f) + (src_y[x] * scale), 0.f, 255.f);
    }
  }
  return WP2_STATUS_OK;
}

template <>
WP2Status Plane8u::ToGray(WP2::ArgbBuffer* dst, uint32_t num_bits,
                          bool shift) const {
  assert(num_bits <= 8);
  (void)shift;
  const uint32_t h = std::min(h_, dst->height);
  const uint32_t w = std::min(w_, dst->width);
  for (uint32_t y = 0; y < h; ++y) {
    uint8_t* const dst_y = (uint8_t*)dst->GetRow(y);
    const uint8_t* const src_y = Row(y);
    for (uint32_t x = 0; x < w; ++x) {
      dst_y[4 * x + 0] = 0xff;
      dst_y[4 * x + 1] = dst_y[4 * x + 2] = dst_y[4 * x + 3] =
          (src_y[x] << (8 - num_bits));
    }
  }
  return WP2_STATUS_OK;
}

//------------------------------------------------------------------------------

template <>
WP2Status Plane16::GetSSE(const Plane16& src, uint64_t* const score) const {
  WP2_CHECK_OK(w_ == src.w_ && h_ == src.h_, WP2_STATUS_BAD_DIMENSION);
  WP2PSNRInit();
  *score = 0;
  for (uint32_t y = 0; y < h_; ++y) {
    *score += WP2SumSquaredError16s(Row(y), src.Row(y), w_);
  }
  return WP2_STATUS_OK;
}

template <>
WP2Status Plane8u::GetSSE(const Plane8u& src, uint64_t* const score) const {
  WP2_CHECK_OK(w_ == src.w_ && h_ == src.h_, WP2_STATUS_BAD_DIMENSION);
  WP2PSNRInit();
  *score = 0;
  for (uint32_t y = 0; y < h_; ++y) {
    *score += WP2SumSquaredError8u(Row(y), src.Row(y), w_);
  }
  return WP2_STATUS_OK;
}

//------------------------------------------------------------------------------

bool YUVPlane::IsEmpty() const {
  assert(U.IsEmpty() == Y.IsEmpty() && V.IsEmpty() == Y.IsEmpty());
  if (Y.IsEmpty()) assert(A.IsEmpty());
  return Y.IsEmpty();
}

bool YUVPlane::IsView() const {
  if (IsEmpty()) return false;
  assert(U.IsView() == Y.IsView() && V.IsView() == Y.IsView());
  if (!A.IsEmpty()) assert(A.IsView() == Y.IsView());
  return Y.IsView();
}

uint32_t YUVPlane::GetWidth() const {
  assert(U.w_ == Y.w_ || U.w_ == DivCeil(Y.w_, 2));  // 4:4:4 or 4:2:0
  assert(V.w_ == U.w_);
  if (!A.IsEmpty()) assert(A.w_ == Y.w_);
  return Y.w_;
}
uint32_t YUVPlane::GetHeight() const {
  assert(U.h_ == Y.h_ || U.h_ == DivCeil(Y.h_, 2));  // 4:4:4 or 4:2:0
  assert(V.h_ == U.h_);
  if (!A.IsEmpty()) assert(A.h_ == Y.h_);
  return Y.h_;
}

void YUVPlane::Clear() {
  A.Clear();
  Y.Clear();
  U.Clear();
  V.Clear();
}

const Plane16& YUVPlane::GetChannel(Channel channel) const {
  if (channel == kYChannel) return Y;
  if (channel == kUChannel) return U;
  if (channel == kVChannel) return V;
  assert(channel == kAChannel);
  return A;
}

Plane16& YUVPlane::GetChannel(Channel channel) {
  if (channel == kYChannel) return Y;
  if (channel == kUChannel) return U;
  if (channel == kVChannel) return V;
  assert(channel == kAChannel);
  return A;
}

//------------------------------------------------------------------------------

WP2Status YUVPlane::Resize(uint32_t width, uint32_t height, uint32_t pad,
                           bool has_alpha, bool as_yuv420) {
  WP2_CHECK_OK(pad > 0, WP2_STATUS_INVALID_PARAMETER);
  const uint32_t w = Pad(width, pad);
  const uint32_t h = Pad(height, pad);
  const uint32_t uv_w = as_yuv420 ? Pad((width  + 1u) >> 1, pad) : w;
  const uint32_t uv_h = as_yuv420 ? Pad((height + 1u) >> 1, pad) : h;
  WP2_CHECK_STATUS(GetChannel(kYChannel).Resize(w, h));
  WP2_CHECK_STATUS(GetChannel(kUChannel).Resize(uv_w, uv_h));
  WP2_CHECK_STATUS(GetChannel(kVChannel).Resize(uv_w, uv_h));
  if (has_alpha) {
    WP2_CHECK_STATUS(GetChannel(kAChannel).Resize(w, h));
  } else {
    GetChannel(kAChannel).Clear();
  }
  return WP2_STATUS_OK;
}

WP2Status YUVPlane::Copy(const YUVPlane& src, bool resize_if_needed) {
  for (Channel channel : {kYChannel, kUChannel, kVChannel, kAChannel}) {
    WP2_CHECK_STATUS(
        GetChannel(channel).Copy(src.GetChannel(channel), resize_if_needed));
  }
  return WP2_STATUS_OK;
}

void YUVPlane::Fill(const Rectangle& rect, const Ayuv38b& color) {
  if (!A.IsEmpty()) A.Fill(rect, color.a);
  Y.Fill(rect, color.y);
  U.Fill(rect, color.u);
  V.Fill(rect, color.v);
}

void YUVPlane::Fill(const Ayuv38b& color) {
  if (!A.IsEmpty()) A.Fill(color.a);
  Y.Fill(color.y);
  U.Fill(color.u);
  V.Fill(color.v);
}

WP2Status YUVPlane::FillPad(uint32_t non_padded_width,
                            uint32_t non_padded_height) {
  WP2_CHECK_OK(!IsDownsampled(), WP2_STATUS_INVALID_PARAMETER);
  for (Channel channel : {kYChannel, kUChannel, kVChannel, kAChannel}) {
    if (channel == kAChannel && A.IsEmpty()) continue;
    WP2_CHECK_STATUS(
        GetChannel(channel).FillPad(non_padded_width, non_padded_height));
  }
  return WP2_STATUS_OK;
}

WP2Status YUVPlane::SetView(const YUVPlane& src, const Rectangle& rect) {
  WP2_CHECK_OK(!IsDownsampled(), WP2_STATUS_INVALID_PARAMETER);
  for (Channel channel : {kYChannel, kUChannel, kVChannel, kAChannel}) {
    if (channel == kAChannel && src.A.IsEmpty()) {
      A.Clear();
    } else {
      WP2_CHECK_STATUS(
          GetChannel(channel).SetView(src.GetChannel(channel), rect));
    }
  }
  return WP2_STATUS_OK;
}

WP2Status YUVPlane::SetView(const YUVPlane& src) {
  for (Channel channel : {kYChannel, kUChannel, kVChannel, kAChannel}) {
    if (channel == kAChannel && src.A.IsEmpty()) {
      A.Clear();
    } else {
      WP2_CHECK_STATUS(GetChannel(channel).SetView(src.GetChannel(channel)));
    }
  }
  return WP2_STATUS_OK;
}

//------------------------------------------------------------------------------
// One-liners.

static WP2Status ResizeIfAllowed(uint32_t width, uint32_t height, uint32_t pad,
                                 bool has_alpha, bool allowed,
                                 YUVPlane* const buffer) {
  if (allowed) {
    WP2_CHECK_STATUS(buffer->Resize(width, height, pad, has_alpha));
  } else {
    WP2_CHECK_OK(pad > 0, WP2_STATUS_INVALID_PARAMETER);
    for (Channel channel : {kYChannel, kUChannel, kVChannel, kAChannel}) {
      Plane16* const plane = &buffer->GetChannel(channel);
      if (channel == kAChannel && !has_alpha) {
        WP2_CHECK_OK(plane->IsEmpty(), WP2_STATUS_BAD_DIMENSION);
      } else {
        WP2_CHECK_OK(plane->w_ == Pad(width, pad), WP2_STATUS_BAD_DIMENSION);
        WP2_CHECK_OK(plane->h_ == Pad(height, pad), WP2_STATUS_BAD_DIMENSION);
      }
    }
  }
  return WP2_STATUS_OK;
}

//------------------------------------------------------------------------------

WP2Status YUVPlane::Import(const ArgbBuffer& rgb, bool import_alpha,
                           const CSPTransform& csp_transform,
                           bool resize_if_needed, size_t pad) {
  WP2_CHECK_STATUS(ResizeIfAllowed(rgb.width, rgb.height, pad, import_alpha,
                                   resize_if_needed, this));
  const uint32_t x_stride = WP2FormatBpp(rgb.format);
  for (size_t y = 0; y < rgb.height; ++y) {
    int16_t* const y_row = Y.Row(y);
    int16_t* const u_row = U.Row(y);
    int16_t* const v_row = V.Row(y);
    int16_t* const a_row = import_alpha ? A.Row(y) : nullptr;
    const uint8_t* rgb_row = (const uint8_t*)rgb.GetRow(y);
    for (size_t x = 0; x < rgb.width; ++x, rgb_row += x_stride) {
      csp_transform.ToYUV(rgb_row[1], rgb_row[2], rgb_row[3],
                          &y_row[x], &u_row[x], &v_row[x]);
      if (import_alpha) a_row[x] = rgb_row[0];
    }
  }
  WP2_CHECK_STATUS(FillPad(rgb.width, rgb.height));
  return WP2_STATUS_OK;
}

WP2Status YUVPlane::Import(const ArgbBuffer& rgb, bool import_alpha,
                           const CSPMtx& rgb_to_ccsp, bool resize_if_needed,
                           size_t pad) {
  int16_t m[9];
  for (uint32_t i = 0; i < 9; ++i) {
    m[i] = ChangePrecision(rgb_to_ccsp.matrix[i], rgb_to_ccsp.shift,
                           CSPTransform::kMtxShift);
    WP2_CHECK_OK(std::abs(m[i]) < (1 << CSPTransform::kMtxBits),
                 WP2_STATUS_INVALID_PARAMETER);
  }
  WP2_CHECK_STATUS(ResizeIfAllowed(rgb.width, rgb.height, pad, import_alpha,
                                   resize_if_needed, this));

  for (size_t y = 0; y < rgb.height; ++y) {
    int16_t* const y_row = Y.Row(y);
    int16_t* const u_row = U.Row(y);
    int16_t* const v_row = V.Row(y);
    int16_t* const a_row = import_alpha ? A.Row(y) : nullptr;
    const uint8_t* const rgb_row = (const uint8_t*)rgb.GetRow(y);
    for (size_t x = 0; x < rgb.width; ++x) {
      const uint32_t rgb_x = WP2FormatBpp(rgb.format) * x;
      if (import_alpha) a_row[x] = rgb_row[rgb_x + 0];
      Multiply(rgb_row[rgb_x + 1], rgb_row[rgb_x + 2], rgb_row[rgb_x + 3], m,
               CSPTransform::kMtxShift, &y_row[x], &u_row[x], &v_row[x]);
    }
  }
  WP2_CHECK_STATUS(FillPad(rgb.width, rgb.height));
  return WP2_STATUS_OK;
}

WP2Status YUVPlane::Import(const YUVPlane& ccsp, const CSPMtx& ccsp_to_rgb,
                           const CSPTransform& csp_transform,
                           bool resize_if_needed, size_t pad) {
  WP2_CHECK_OK(!ccsp.IsDownsampled(), WP2_STATUS_INVALID_PARAMETER);
  const uint32_t width = ccsp.GetWidth(), height = ccsp.GetHeight();
  WP2_CHECK_STATUS(ResizeIfAllowed(width, height, pad, ccsp.HasAlpha(),
                                   resize_if_needed, this));

  WP2_CHECK_STATUS(csp_transform.CustomToYUV(
      width, height, ccsp.Y.Row(0), ccsp.Y.Step(), ccsp.U.Row(0), ccsp.U.Step(),
      ccsp.V.Row(0), ccsp.V.Step(), ccsp_to_rgb, Y.Row(0), Y.Step(), U.Row(0),
      U.Step(), V.Row(0), V.Step()));
  if (!ccsp.A.IsEmpty()) {
    Plane16 non_padded_alpha;
    WP2_CHECK_STATUS(non_padded_alpha.SetView(A, {0, 0, width, height}));
    WP2_CHECK_STATUS(non_padded_alpha.Copy(ccsp.A, /*resize_if_needed=*/false));
  }

  WP2_CHECK_STATUS(FillPad(width, height));
  return WP2_STATUS_OK;
}

//------------------------------------------------------------------------------

WP2Status YUVPlane::Export(const CSPTransform& transform, bool resize_if_needed,
                           ArgbBuffer* const dst,
                           const SamplingTaps* const upsample_if_needed) const {
  if (IsDownsampled()) {
    WP2_CHECK_OK(upsample_if_needed != nullptr, WP2_STATUS_INVALID_COLORSPACE);
    YUVPlane upsampled;
    WP2_CHECK_STATUS(upsampled.UpsampleFrom(*this, *upsample_if_needed,
                                            /*use_views_for_ya=*/true));
    WP2_CHECK_STATUS(upsampled.Export(transform, resize_if_needed, dst));
    return WP2_STATUS_OK;
  }
  if (dst->width != GetWidth() || dst->height != GetHeight()) {
    WP2_CHECK_OK(resize_if_needed, WP2_STATUS_BAD_DIMENSION);
    WP2_CHECK_STATUS(dst->Resize(GetWidth(), GetHeight()));
  }
  // we only handle XRGB32 or ARGB32 for now
  WP2_CHECK_OK(WP2FormatBpp(dst->format) == 4, WP2_STATUS_INVALID_COLORSPACE);
  const bool has_alpha = HasAlpha();
  static_assert(kRgbBits == 8, "YUVPlane::Export() handles only 8b RGB");
  WP2CSPConverterInit();
  WP2YUVToArgbFunc cvrt_func = has_alpha ? WP2YUVToArgb : WP2YUVToXRGB;
  for (size_t y = 0; y < dst->height; ++y) {
    const int16_t* const y_row = Y.Row(y);
    const int16_t* const u_row = U.Row(y);
    const int16_t* const v_row = V.Row(y);
    const int16_t* const a_row = has_alpha ? A.Row(y) : nullptr;
    uint8_t* const argb = (uint8_t*)dst->GetRow(y);
    cvrt_func(y_row, u_row, v_row, a_row,
              transform.GetRgbAverage(), transform.GetYuvToRgbMatrix(),
              argb, dst->width);
  }
  return WP2_STATUS_OK;
}

WP2Status YUVPlane::Export(const CSPMtx& ccsp_to_rgb, bool resize_if_needed,
                           ArgbBuffer* const dst,
                           const SamplingTaps* const upsample_if_needed) const {
  if (IsDownsampled()) {
    WP2_CHECK_OK(upsample_if_needed != nullptr, WP2_STATUS_INVALID_COLORSPACE);
    YUVPlane upsampled;
    WP2_CHECK_STATUS(upsampled.UpsampleFrom(*this, *upsample_if_needed,
                                            /*use_views_for_ya=*/true));
    WP2_CHECK_STATUS(upsampled.Export(ccsp_to_rgb, resize_if_needed, dst));
    return WP2_STATUS_OK;
  }
  if (dst->width != GetWidth() || dst->height != GetHeight()) {
    WP2_CHECK_OK(resize_if_needed, WP2_STATUS_BAD_DIMENSION);
    WP2_CHECK_STATUS(dst->Resize(GetWidth(), GetHeight()));
  }

  WP2_CHECK_STATUS(CSPTransform::CustomToArgb(
      dst->width, dst->height, Y.Row(0), Y.Step(), U.Row(0), U.Step(), V.Row(0),
      V.Step(), HasAlpha() ? A.Row(0) : nullptr, A.Step(), ccsp_to_rgb, dst));
  return WP2_STATUS_OK;
}

//------------------------------------------------------------------------------

WP2Status YUVPlane::DownsampleFrom(const YUVPlane& src, const SamplingTaps& f,
                                   bool use_views_for_ya) {
  WP2_CHECK_OK(!src.IsDownsampled(), WP2_STATUS_INVALID_PARAMETER);
  if (use_views_for_ya) {
    WP2_CHECK_STATUS(Y.SetView(src.Y));
    WP2_CHECK_STATUS(A.SetView(src.A));
  } else {
    WP2_CHECK_STATUS(Y.Copy(src.Y, /*resize_if_needed=*/true));
    WP2_CHECK_STATUS(A.Copy(src.A, /*resize_if_needed=*/true));
  }
  WP2_CHECK_STATUS(U.DownsampleFrom(src.U, f));
  WP2_CHECK_STATUS(V.DownsampleFrom(src.V, f));
  return WP2_STATUS_OK;
}

WP2Status YUVPlane::Downsample(const SamplingTaps& f) {
  WP2_CHECK_OK(!IsDownsampled(), WP2_STATUS_INVALID_PARAMETER);
  WP2_CHECK_STATUS(U.Downsample(f));
  WP2_CHECK_STATUS(V.Downsample(f));
  return WP2_STATUS_OK;
}

WP2Status YUVPlane::UpsampleFrom(const YUVPlane& src, const SamplingTaps& f,
                                 bool use_views_for_ya) {
  if (src.GetWidth() == 1 && src.GetHeight() == 1) {
    if (use_views_for_ya) return SetView(src);
    return Copy(src, /*resize_if_needed=*/true);
  }
  WP2_CHECK_OK(src.IsDownsampled(), WP2_STATUS_INVALID_PARAMETER);
  if (use_views_for_ya) {
    WP2_CHECK_STATUS(Y.SetView(src.Y));
    if (src.HasAlpha()) {
      WP2_CHECK_STATUS(A.SetView(src.A));
    } else {
      A.Clear();
    }
  } else {
    WP2_CHECK_STATUS(Y.Copy(src.Y, /*resize_if_needed=*/true));
    WP2_CHECK_STATUS(A.Copy(src.A, /*resize_if_needed=*/true));
  }
  WP2_CHECK_STATUS(U.UpsampleFrom(src.U, Y.w_, Y.h_, f));
  WP2_CHECK_STATUS(V.UpsampleFrom(src.V, Y.w_, Y.h_, f));
  return WP2_STATUS_OK;
}

WP2Status YUVPlane::Upsample(const SamplingTaps& f) {
  if (GetWidth() == 1 && GetHeight() == 1) return WP2_STATUS_OK;
  WP2_CHECK_OK(IsDownsampled(), WP2_STATUS_INVALID_PARAMETER);
  WP2_CHECK_STATUS(U.Upsample(Y.w_, Y.h_, f));
  WP2_CHECK_STATUS(V.Upsample(Y.w_, Y.h_, f));
  return WP2_STATUS_OK;
}

bool YUVPlane::IsDownsampled() const {
  // GetWidth() and GetHeight() already verify the dimensions of the planes.
  return (U.w_ != GetWidth() || U.h_ != GetHeight());
}

// Checks that Y, U, V and A planes (if present) are all the same size.
static WP2Status CheckPlanesSameSize(const YUVPlane& yuv) {
  WP2_CHECK_OK(yuv.Y.IsSameSize(yuv.U), WP2_STATUS_BAD_DIMENSION);
  WP2_CHECK_OK(yuv.Y.IsSameSize(yuv.V), WP2_STATUS_BAD_DIMENSION);
  if (yuv.HasAlpha()) {
    WP2_CHECK_OK(yuv.Y.IsSameSize(yuv.A), WP2_STATUS_BAD_DIMENSION);
  }
  return WP2_STATUS_OK;
}

static WP2Status Composite(const YUVPlane& background,
                           const YUVPlane& foreground,
                           const CSPTransform& transform, YUVPlane* const dst) {
  WP2_CHECK_STATUS(CheckPlanesSameSize(background));
  WP2_CHECK_STATUS(CheckPlanesSameSize(foreground));
  WP2_CHECK_OK(background.Y.IsSameSize(foreground.Y), WP2_STATUS_BAD_DIMENSION);

  const Ayuv38b shift = transform.ToYUV({255, 0, 0, 0});
  for (uint32_t y = 0; y < foreground.Y.h_; ++y) {
    int16_t* const dst_row[4] = {dst->Y.Row(y), dst->U.Row(y), dst->V.Row(y),
                                 dst->A.Row(y)};
    const int16_t* const fg_row[4] = {foreground.Y.Row(y), foreground.U.Row(y),
                                      foreground.V.Row(y), foreground.A.Row(y)};
    const int16_t* const bg_row[4] = {background.Y.Row(y), background.U.Row(y),
                                      background.V.Row(y), background.A.Row(y)};
    for (uint32_t x = 0; x < foreground.Y.w_; ++x) {
      const uint8_t b = kAlphaMax - fg_row[kAChannel][x];
      dst_row[kYChannel][x] =
          fg_row[kYChannel][x] + DivBy255(b * (bg_row[kYChannel][x] - shift.y));
      dst_row[kUChannel][x] =
          fg_row[kUChannel][x] + DivBy255(b * (bg_row[kUChannel][x] - shift.u));
      dst_row[kVChannel][x] =
          fg_row[kVChannel][x] + DivBy255(b * (bg_row[kVChannel][x] - shift.v));
      dst_row[kAChannel][x] =
          fg_row[kAChannel][x] + DivBy255(b * bg_row[kAChannel][x]);
    }
  }
  return WP2_STATUS_OK;
}

WP2Status YUVPlane::CompositeOver(const YUVPlane& background,
                                  const CSPTransform& transform) {
  if (!HasAlpha()) return WP2_STATUS_OK;
  return Composite(background, *this, transform, this);
}

WP2Status YUVPlane::CompositeUnder(const YUVPlane& foreground,
                                   const CSPTransform& transform) {
  if (!foreground.HasAlpha()) {
    return Copy(foreground, /*resize_if_needed=*/false);
  }
  return Composite(*this, foreground, transform, this);
}

//------------------------------------------------------------------------------

void swap(YUVPlane& a, YUVPlane& b) {
  using std::swap;
  swap(a.Y, b.Y);
  swap(a.U, b.U);
  swap(a.V, b.V);
  swap(a.A, b.A);
}

//------------------------------------------------------------------------------

}    // namespace WP2
