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
// Colorspace utilities
//
// Author: Skal (pascal.massimino@gmail.com)

#include "src/utils/csp.h"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdio>

#include "src/dsp/math.h"
#include "src/utils/ans.h"
#include "src/utils/ans_utils.h"
#include "src/utils/utils.h"

namespace WP2 {

//------------------------------------------------------------------------------

static int16_t ToInt16(int32_t v, int32_t min_v, int32_t max_v) {
  assert(std::abs(v) < max_v + 8);  // Rounding errors may happen.
  return (int16_t)Clamp(v, min_v, max_v);
}
static int16_t RoundToInt16(double v, int16_t max_v) {
  assert(v >= -(double)max_v && v <= (double)max_v);
  return (int16_t)lround(v);
}

static bool TryInvertMatrix(const double m[9], double out[9]) {
  const double det = m[0] * (m[4] * m[8] - m[5] * m[7]) -
                     m[3] * (m[1] * m[8] - m[7] * m[2]) +
                     m[6] * (m[1] * m[5] - m[4] * m[2]);
  if (det == 0.) return false;
  const double idet = 1.0 / det;
  out[0] = idet * (m[4] * m[8] - m[5] * m[7]);
  out[3] = idet * (m[6] * m[5] - m[3] * m[8]);
  out[6] = idet * (m[3] * m[7] - m[4] * m[6]);
  out[1] = idet * (m[7] * m[2] - m[1] * m[8]);
  out[4] = idet * (m[0] * m[8] - m[2] * m[6]);
  out[7] = idet * (m[1] * m[6] - m[0] * m[7]);
  out[2] = idet * (m[1] * m[5] - m[2] * m[4]);
  out[5] = idet * (m[2] * m[3] - m[0] * m[5]);
  out[8] = idet * (m[0] * m[4] - m[1] * m[3]);
  return true;
}

static double Norm2(const double a[3]) {
  return a[0] * a[0] + a[1] * a[1] + a[2] * a[2];
}

static void CrossProduct(const double a[3], const double b[3], double out[3]) {
  out[0] = a[1] * b[2] - a[2] * b[1];
  out[1] = a[2] * b[0] - a[0] * b[2];
  out[2] = a[0] * b[1] - a[1] * b[0];
  double norm = Norm2(out);
  if (norm > 0.) {
    norm = 1. / std::sqrt(norm);
    out[0] *= norm;
    out[1] *= norm;
    out[2] *= norm;
  }
}

static float Distance(const double a[3], const double b[3]) {
  float dot = std::fabs(a[0] * b[0] + a[1] * b[1] + a[2] * b[2]);
  if (dot != 0.) dot /= std::sqrt(Norm2(a) * Norm2(b));
  return 1.f - dot;
}

static void MultiplyMatrix(double m[9], double by) {
  for (int i = 0; i < 9; ++i) m[i] *= by;
}

//------------------------------------------------------------------------------

constexpr uint32_t CSPTransform::kMtxBits;
constexpr uint32_t CSPTransform::kMtxShift;
constexpr int16_t CSPTransform::kDefaultRgbAvg[3];

static constexpr int32_t kMtxFixedPointMul = 1 << CSPTransform::kMtxShift;
static constexpr int32_t kMtxRound = kMtxFixedPointMul >> 1;
static constexpr int16_t kMaxMtxValue = (1 << CSPTransform::kMtxBits) - 1;

//------------------------------------------------------------------------------

void CSPTransform::StoreCustomYUVBounds() {  // Compute min/max YUV values.
  int32_t min[3] = {0}, max[3] = {0};
  for (uint32_t i = 0; i < 3; ++i) {
    for (uint32_t j = 0; j < 3; ++j) {
      const int16_t rgb_min = 0 - avg_[j];
      const int16_t rgb_max = kRgbMax - avg_[j];

      const int32_t coeff = m_[i * 3 + j];
      min[i] += coeff * ((coeff > 0) ? rgb_min : rgb_max);
      max[i] += coeff * ((coeff > 0) ? rgb_max : rgb_min);
    }
    min[i] = (min[i] + kMtxRound) >> kMtxShift;
    max[i] = (max[i] + kMtxRound) >> kMtxShift;
  }
  const int32_t global_min = std::min({min[0], min[1], min[2]});
  const int32_t global_max = std::max({max[0], max[1], max[2]});
  SetYUVBounds(global_min, global_max);
}

bool CSPTransform::Init(const int16_t yuv_to_rgb_matrix[9],
                        const int16_t rgb_avg[3]) {
  for (size_t i = 0; i < 9; ++i) im_[i] = yuv_to_rgb_matrix[i];
  for (uint32_t i : {0, 1, 2}) avg_[i] = rgb_avg[i];
  if (!ComputeRGBToYUVMatrix()) return false;
  StoreCustomYUVBounds();
  type_ = Csp::kCustom;
  return true;
}

WP2Status CSPTransform::Init(Csp csp_type, const ArgbBuffer& argb_buffer) {
  if (csp_type == Csp::kCustom) {
    assert(!argb_buffer.IsEmpty());
    const WP2Status status = Optimize(argb_buffer);
    // Stop here if success or unrecoverable error. Otherwise fallback on YCoCg.
    if (status != WP2_STATUS_INVALID_COLORSPACE) return status;
    csp_type = Csp::kYCoCg;
  }
  if (csp_type == Csp::kYCoCg) {
    InitYCoCg();
  } else if (csp_type == Csp::kYCbCr) {
    InitYCbCr();
  } else {
    assert(csp_type == Csp::kYIQ);
    InitYIQ();
  }
  return WP2_STATUS_OK;
}

void CSPTransform::InitYCoCg() {
  // The RGB to YUV matrix has to be scaled by 4 to be lossless, so we divide
  // the inverse by 4.
  static constexpr int16_t I = (kMtxFixedPointMul >> 2);
  // clang-format off
  static constexpr int16_t kYCoCgMatrix[] = {
      I,  I, -I,
      I,  0,  I,
      I, -I, -I
  };
  // clang-format on
  if (!Init(kYCoCgMatrix)) assert(false);
  type_ = Csp::kYCoCg;
  // With 10 bits par channel ([-512:511] if 'rgb_avg' is 128), the YCoCg
  // transformation is lossless.
  assert(yuv_bits_ <= kMaxYuvBits + 1);
}

void CSPTransform::InitYCbCr() {
#if 0  // reference code used for generating kYCbCrMatrix[]
  double yuv_to_rgb_matrix[] = {
      1.16,  0.00,  1.60,
      1.16, -0.39, -0.81,
      1.16,  2.01,  0.00
  };
  // This is not lossless: use the maximum range (yuv_bits_ == kMaxYuvBits + 1).
  MultiplyMatrix(yuv_to_rgb_matrix, kMtxFixedPointMul / 4);
  int16_t kYCbCrMatrix[9];
  for (int i = 0; i < 9; ++i) {
    kYCbCrMatrix[i] = RoundToInt16(yuv_to_rgb_matrix[i], kMaxMtxValue);
  }
  for (auto v : kYCbCrMatrix) printf("%6d,\n", v);
#else
  static constexpr int16_t kYCbCrMatrix[] = {
      1188,    0, 1638,
      1188, -399, -829,
      1188, 2058,    0,
  };
#endif
  if (!Init(kYCbCrMatrix)) assert(false);
  type_ = Csp::kYCbCr;
  assert(yuv_bits_ == kMaxYuvBits + 1);
}

void CSPTransform::InitYIQ() {
#if 0  // reference code used for generating kIQMatrix[]
  static constexpr double rgb_to_yuv_matrix[] = {
    0.299,  0.587,  0.114,
    0.596, -0.274, -0.322,
    0.211, -0.523,  0.312
  };
  double yuv_to_rgb_matrix[9];
  TryInvertMatrix(rgb_to_yuv_matrix, yuv_to_rgb_matrix);
  // This is not lossless: use the maximum range (yuv_bits_ == kMaxYuvBits + 1).
  MultiplyMatrix(yuv_to_rgb_matrix, kMtxFixedPointMul / 2);
  int16_t kYIQMatrix[9];
  for (int i = 0; i < 9; ++i) {
    kYIQMatrix[i] = RoundToInt16(yuv_to_rgb_matrix[i], kMaxMtxValue);
  }
  for (auto v : kYIQMatrix) printf("%6d,\n", v);
#else
  static constexpr int16_t kYIQMatrix[] = {
      2048,  1958,  1273,
      2048,  -558, -1325,
      2048, -2260,  3483,
  };
#endif
  if (!Init(kYIQMatrix)) assert(false);
  type_ = Csp::kYIQ;
  assert(yuv_bits_ == kMaxYuvBits + 1);
}

void CSPTransform::Print() const {
  printf("CSP type %u\n", (uint32_t)type_);
  double rgb_to_yuv_matrix[9], yuv_to_rgb_matrix[9];
  for (int i = 0; i < 9; ++i) {
    rgb_to_yuv_matrix[i] = (double)m_[i] / kMtxFixedPointMul;
    yuv_to_rgb_matrix[i] = (double)im_[i] / kMtxFixedPointMul;
  }

  printf("=== RGB to YUV: ===\n");
  for (int j = 0; j < 3; ++j) {
    for (int i = 0; i < 3; ++i) {
      printf("%10.4f ", rgb_to_yuv_matrix[i + j * 3]);
    }
    printf("\n");
  }

  printf("=== YUV to RGB: ===\n");
  for (int j = 0; j < 3; ++j) {
    for (int i = 0; i < 3; ++i) {
      printf("%10.4f ", yuv_to_rgb_matrix[i + j * 3]);
    }
    printf("\n");
  }

  printf("=== avg: ");
  printf("%3d %3d %3d", avg_[0], avg_[1], avg_[2]);
  printf("\n");
}

//------------------------------------------------------------------------------

void CSPTransform::SetYUVBounds(int16_t yuv_min, int16_t yuv_max) {
  yuv_bits_ =
      std::ceil(std::log2(std::max((int16_t)std::abs(yuv_min), yuv_max))) + 1;
  min_yuv_value_ = yuv_min;
  max_yuv_value_ = yuv_max;
}

bool CSPTransform::ComputeRGBToYUVMatrix() {
  double rgb_to_yuv_matrix[9], yuv_to_rgb_matrix[9];
  constexpr double norm = 1. / kMtxFixedPointMul;
  for (int i = 0; i < 9; ++i) {
    if (std::abs(im_[i]) > kMaxMtxValue) return false;
    yuv_to_rgb_matrix[i] = norm * im_[i];
  }

  if (!TryInvertMatrix(yuv_to_rgb_matrix, rgb_to_yuv_matrix)) return false;
  for (int i = 0; i < 9; ++i) {
    const double v = rgb_to_yuv_matrix[i] * kMtxFixedPointMul;
    if (std::abs(v) > kMaxMtxValue) return false;
    m_[i] = RoundToInt16(v, kMaxMtxValue);
  }
  return true;
}

//------------------------------------------------------------------------------

template <typename Trgb>
static inline void RGBToYUV(const int16_t avg[3], const int16_t m[9],
                            int16_t min_yuv_value, int16_t max_yuv_value,
                            Trgb r, Trgb g, Trgb b, int16_t* const y,
                            int16_t* const u, int16_t* const v) {
  // kRgbBits to (kRgbBits+1)
  const int32_t tmp_r = (int32_t)r - avg[0];
  const int32_t tmp_g = (int32_t)g - avg[1];
  const int32_t tmp_b = (int32_t)b - avg[2];
  // (kRgbBits+1) to (yuv_bits_ + kMtxShift)
  int32_t tmp_y = m[0] * tmp_r + m[1] * tmp_g + m[2] * tmp_b;
  int32_t tmp_u = m[3] * tmp_r + m[4] * tmp_g + m[5] * tmp_b;
  int32_t tmp_v = m[6] * tmp_r + m[7] * tmp_g + m[8] * tmp_b;
  // (yuv_bits_ + kMtxShift) to yuv_bits_
  tmp_y = RightShift(tmp_y + kMtxRound, CSPTransform::kMtxShift);
  tmp_u = RightShift(tmp_u + kMtxRound, CSPTransform::kMtxShift);
  tmp_v = RightShift(tmp_v + kMtxRound, CSPTransform::kMtxShift);
  *y = ToInt16(tmp_y, min_yuv_value, max_yuv_value);
  *u = ToInt16(tmp_u, min_yuv_value, max_yuv_value);
  *v = ToInt16(tmp_v, min_yuv_value, max_yuv_value);
}

template <typename Trgb>
static inline void YUVToRGB(const int16_t avg[3], const int16_t im[9],
                            int16_t y, int16_t u, int16_t v,
                            Trgb* const r, Trgb* const g, Trgb* const b) {
  // yuv_bits_ to (kRgbBits+1 + kMtxShift)
  int32_t tmp_r = im[0] * y + im[1] * u + im[2] * v;
  int32_t tmp_g = im[3] * y + im[4] * u + im[5] * v;
  int32_t tmp_b = im[6] * y + im[7] * u + im[8] * v;
  // (kRgbBits+1 + kMtxShift) to kRgbBits+1
  tmp_r = RightShift(tmp_r + kMtxRound, CSPTransform::kMtxShift);
  tmp_g = RightShift(tmp_g + kMtxRound, CSPTransform::kMtxShift);
  tmp_b = RightShift(tmp_b + kMtxRound, CSPTransform::kMtxShift);
  // (kRgbBits+1) to kRgbBits
  *r = (Trgb)Clamp<int32_t>(tmp_r + avg[0], 0, kRgbMax);
  *g = (Trgb)Clamp<int32_t>(tmp_g + avg[1], 0, kRgbMax);
  *b = (Trgb)Clamp<int32_t>(tmp_b + avg[2], 0, kRgbMax);
  // Don't assert that rgb fits in kRgbBits, just clamp it (lossy conversion).
}

void CSPTransform::ToYUV(int16_t r, int16_t g, int16_t b, int16_t* const y,
                         int16_t* const u, int16_t* const v) const {
  RGBToYUV(avg_, m_, min_yuv_value_, max_yuv_value_, r, g, b, y, u, v);
}

void CSPTransform::ToRGB(int16_t y, int16_t u, int16_t v, int16_t* const r,
                         int16_t* const g, int16_t* const b) const {
  YUVToRGB(avg_, im_, y, u, v, r, g, b);
}

void CSPTransform::ToYUV(uint8_t r, uint8_t g, uint8_t b, int16_t* const y,
                         int16_t* const u, int16_t* const v) const {
  RGBToYUV(avg_, m_, min_yuv_value_, max_yuv_value_, r, g, b, y, u, v);
}

void CSPTransform::ToRGB(int16_t y, int16_t u, int16_t v, uint8_t* const r,
                         uint8_t* const g, uint8_t* const b) const {
  YUVToRGB(avg_, im_, y, u, v, r, g, b);
}

void CSPTransform::ToYUV(const int16_t rgb[3], int16_t yuv[3]) const {
  RGBToYUV(avg_, m_, min_yuv_value_, max_yuv_value_, rgb[0], rgb[1], rgb[2],
           &yuv[0], &yuv[1], &yuv[2]);
  // CheckRoundTrip(yuv, rgb);
}

void CSPTransform::ToRGB(const int16_t yuv[3], int16_t rgb[3]) const {
  YUVToRGB(avg_, im_, yuv[0], yuv[1], yuv[2], &rgb[0], &rgb[1], &rgb[2]);
}

void CSPTransform::ToYUV(const uint8_t rgb[3], int16_t yuv[3]) const {
  RGBToYUV(avg_, m_, min_yuv_value_, max_yuv_value_, rgb[0], rgb[1], rgb[2],
           &yuv[0], &yuv[1], &yuv[2]);
}

void CSPTransform::ToRGB(const int16_t yuv[3], uint8_t rgb[3]) const {
  YUVToRGB(avg_, im_, yuv[0], yuv[1], yuv[2], &rgb[0], &rgb[1], &rgb[2]);
}

Ayuv38b CSPTransform::ToYUV(const Argb32b& color) const {
  Ayuv38b result;
  result.a = color.a;
  RGBToYUV(avg_, m_, min_yuv_value_, max_yuv_value_, color.r, color.g, color.b,
           &result.y, &result.u, &result.v);
  return result;
}

Argb32b CSPTransform::ToRGB(const Ayuv38b& color) const {
  Argb32b result;
  result.a = color.a;
  YUVToRGB(avg_, im_, color.y, color.u, color.v, &result.r, &result.g,
           &result.b);
  return result;
}

//------------------------------------------------------------------------------

WP2Status CSPTransform::CustomToArgb(uint32_t width, uint32_t height,
                                     const int16_t* c0_buffer, uint32_t c0_step,
                                     const int16_t* c1_buffer, uint32_t c1_step,
                                     const int16_t* c2_buffer, uint32_t c2_step,
                                     const int16_t* a_buffer, uint32_t a_step,
                                     const CSPMtx& ccsp_to_rgb,
                                     ArgbBuffer* const argb) {
  WP2_CHECK_OK(argb->width == width && argb->height == height,
               WP2_STATUS_BAD_DIMENSION);
  int16_t m[9];
  for (uint32_t i = 0; i < 9; ++i) {
    const int32_t e =
        ChangePrecision(ccsp_to_rgb.matrix[i], ccsp_to_rgb.shift, kMtxShift);
    WP2_CHECK_OK(std::abs(e) < kMaxMtxValue, WP2_STATUS_INVALID_PARAMETER);
    m[i] = (int16_t)e;
  }

  for (uint32_t y = 0; y < height; ++y) {
    uint8_t* const dst_row = (uint8_t*)argb->GetRow(y);
    if (a_buffer != nullptr) {
      for (uint32_t x = 0; x < width; ++x) {
        uint8_t* const dst_pixel = &dst_row[x * 4];
        dst_pixel[0] =
            (uint8_t)Clamp<int32_t>(a_buffer[y * a_step + x], 0, kAlphaMax);
        int32_t rgb[3];
        Multiply(c0_buffer[x], c1_buffer[x], c2_buffer[x], m, &rgb[0], &rgb[1],
                 &rgb[2]);
        for (uint32_t i : {0, 1, 2}) {
          dst_pixel[1 + i] = (uint8_t)Clamp<int32_t>(
              RightShiftRound(rgb[i], kMtxShift), 0, dst_pixel[0]);
        }
      }
    } else {
      for (uint32_t x = 0; x < width; ++x) {
        uint8_t* const dst_pixel = &dst_row[x * 4];
        dst_pixel[0] = kAlphaMax;
        int32_t rgb[3];
        Multiply(c0_buffer[x], c1_buffer[x], c2_buffer[x], m, &rgb[0], &rgb[1],
                 &rgb[2]);
        for (uint32_t i : {0, 1, 2}) {
          dst_pixel[1 + i] = (uint8_t)Clamp<int32_t>(
              RightShiftRound(rgb[i], kMtxShift), 0, kAlphaMax);
        }
      }
    }
    c0_buffer += c0_step;
    c1_buffer += c1_step;
    c2_buffer += c2_step;
  }
  return WP2_STATUS_OK;
}

WP2Status CSPTransform::CustomToYUV(uint32_t width, uint32_t height,
                                    const int16_t* c0_buffer, uint32_t c0_step,
                                    const int16_t* c1_buffer, uint32_t c1_step,
                                    const int16_t* c2_buffer, uint32_t c2_step,
                                    const CSPMtx& ccsp_to_rgb,
                                    int16_t* y_buffer, uint32_t y_step,
                                    int16_t* u_buffer, uint32_t u_step,
                                    int16_t* v_buffer, uint32_t v_step) const {
  // 'ccsp_to_rgb_matrix' being A, 'm_' (from RGB to YUV) being M:
  // A.ccsp = rgb    and    M.(rgb - avg) = yuv    so    M.A.ccsp - M.avg = yuv

  const int16_t* const a = ccsp_to_rgb.mtx();  // matrix A
  const uint32_t a_shift = ccsp_to_rgb.shift;   // fixed point bits of A

  int16_t ma[9] = {};  // matrix M.A
  {
    int32_t ma_tmp[9] = {};  // 32 bits necessary for matrix multiplication.
    Multiply(m_, a, ma_tmp);
    for (uint32_t i = 0; i < 9; ++i) {
      ma_tmp[i] = ChangePrecision(ma_tmp[i], kMtxShift + a_shift, kMtxShift);
      WP2_CHECK_OK(std::abs(ma_tmp[i]) <= kMaxMtxValue,
                   WP2_STATUS_INVALID_PARAMETER);
      ma[i] = (int16_t)ma_tmp[i];
    }
  }

  int32_t offset[3];  // M.avg
  Multiply<int32_t>(avg_[0], avg_[1], avg_[2], m_,
                    &offset[0], &offset[1], &offset[2]);
  // 'offset' is now scaled by kMtxShift so it can be subtracted as is to 'yuv'
  // below which will also be scaled by kMtxShift.

  for (uint32_t y = 0; y < height; ++y) {
    for (uint32_t x = 0; x < width; ++x) {
      // TODO(skal): this loop should go in dsp/
      int32_t yuv[3];
      Multiply<int32_t>(c0_buffer[x], c1_buffer[x], c2_buffer[x], ma,
                        &yuv[0], &yuv[1], &yuv[2]);
      for (uint32_t i : {0, 1, 2}) {
        yuv[i] = RightShiftRound(yuv[i] - offset[i], kMtxShift);
        yuv[i] = Clamp<int32_t>(yuv[i], min_yuv_value_, max_yuv_value_);
      }
      y_buffer[x] = (int16_t)yuv[0];
      u_buffer[x] = (int16_t)yuv[1];
      v_buffer[x] = (int16_t)yuv[2];
    }
    c0_buffer += c0_step;
    c1_buffer += c1_step;
    c2_buffer += c2_step;
    y_buffer += y_step;
    u_buffer += u_step;
    v_buffer += v_step;
  }
  return WP2_STATUS_OK;
}

WP2Status CSPTransform::YUVToCustom(
    uint32_t width, uint32_t height, const int16_t* y_buffer, uint32_t y_step,
    const int16_t* u_buffer, uint32_t u_step, const int16_t* v_buffer,
    uint32_t v_step, const CSPMtx& rgb_to_ccsp, int16_t* c0_buffer,
    uint32_t c0_step, int16_t* c1_buffer, uint32_t c1_step, int16_t* c2_buffer,
    uint32_t c2_step) const {
  // 'rgb_to_ccsp_matrix' being B, 'm_' (from RGB to YUV) being M and
  // 'im_' (from YUV to RGB) being N:
  // B.rgb = ccsp    and    N.yuv + avg = rgb    so    B.N.(yuv + M.avg) = ccsp

  const int16_t* const b = rgb_to_ccsp.mtx();  // matrix B
  const uint32_t b_shift = rgb_to_ccsp.shift;   // fixed point bits of B

  int16_t bn[9] = {};  // matrix B.N
  {
    int32_t bn_tmp[9] = {};  // 32 bits necessary for matrix multiplication.
    Multiply(b, im_, bn_tmp);
    for (uint32_t i = 0; i < 9; ++i) {
      bn_tmp[i] = ChangePrecision(bn_tmp[i], b_shift + kMtxShift, kMtxShift);
      WP2_CHECK_OK(std::abs(bn_tmp[i]) < kMaxMtxValue,
                   WP2_STATUS_INVALID_PARAMETER);
      bn[i] = (int16_t)bn_tmp[i];
    }
  }

  int16_t offset[3];  // M.avg
  Multiply(avg_[0], avg_[1], avg_[2], m_, kMtxShift,
           &offset[0], &offset[1], &offset[2]);

  WP2CSPConverterInit();
  for (uint32_t y = 0; y < height; ++y) {
    WP2YUVToCustom(y_buffer, u_buffer, v_buffer, offset,
                   bn, c0_buffer, c1_buffer, c2_buffer, width);
    y_buffer += y_step;
    u_buffer += u_step;
    v_buffer += v_step;
    c0_buffer += c0_step;
    c1_buffer += c1_step;
    c2_buffer += c2_step;
  }
  return WP2_STATUS_OK;
}

//------------------------------------------------------------------------------

bool CSPTransform::CheckRoundTrip(const int16_t yuv[3],
                                  const int16_t rgb[3]) const {
  int16_t out[3];
  ToRGB(yuv, out);
  const int16_t delta =
      abs(rgb[0] - out[0]) + abs(rgb[1] - out[1]) + abs(rgb[2] - out[2]);
  const int16_t kErrorThreshold = 5;
  if (delta > kErrorThreshold) {
    static int cnt = 0;
    printf("#%d: %d %d %d / %d %d %d (yuv=%d %d %d) err=%d\n",
           cnt, rgb[0], rgb[1], rgb[2], out[0], out[1], out[2],
           yuv[0], yuv[1], yuv[2], delta);
    if (cnt == 0) Print();
    ++cnt;
    if (cnt == 10) exit(0);
    return false;
  }
  return true;
}

//------------------------------------------------------------------------------

// Returns |v|^2. Result is <= 3 * 2^(15+15).
static int64_t GetVecSqLen(int16_t v0, int16_t v1, int16_t v2) {
  return (int64_t)v0 * v0 + (int64_t)v1 * v1 + (int64_t)v2 * v2;
}

// Predicts a vector element from the others and the magnitude.
static int16_t CompleteVector(int16_t v0, int16_t v1, int64_t v_sq_len) {
  // No need to worry about fixed point.
  int64_t v2 =
      SqrtFloor(SafeSub(v_sq_len, (int64_t)v0 * v0 + (int64_t)v1 * v1));
  if (std::abs(GetVecSqLen(v0, v1, v2) - v_sq_len) >
      std::abs(GetVecSqLen(v0, v1, v2 + 1) - v_sq_len)) {
    v2 += 1;  // In case sqrt floor was "worse" than ceil.
  }
  return (int16_t)Clamp<int64_t>(v2, std::numeric_limits<int16_t>::min(),
                                     std::numeric_limits<int16_t>::max());
}

//------------------------------------------------------------------------------

static void PutSValue(int32_t v, uint32_t bits, ANSEnc* const enc) {
  if (enc->PutBool(v != 0, "non_zero")) {
    const bool sign = (v < 0);
    enc->PutUValue((uint32_t)abs(v) - 1, bits - 1, "abs_val");
    enc->PutBool(sign, "sign");
  }
}

static int32_t ReadSValue(uint32_t bits, ANSDec* const dec) {
  int32_t v = 0;
  if (dec->ReadBool("non_zero")) {
    v = 1 + (int32_t)dec->ReadUValue(bits - 1, "abs_val");
    if (dec->ReadBool("sign")) v = -v;
  }
  return v;
}

// The matrix needs to be precise enough but encodable through ANS.
static constexpr int32_t kEigenBits = 14;
static_assert(kEigenBits <= CSPTransform::kMtxBits &&
                  kEigenBits <= kANSMaxUniformBits &&
                  kEigenBits <= kANSMaxRangeBits,
              "kEigenBits");

static constexpr int32_t kMaxDelta = 5;
static constexpr int32_t kMaxDeltaBits = WP2Log2Ceil_k(kMaxDelta) + 1;

static inline int32_t PutDelta(int16_t original, int16_t predicted,
                        ANSEnc* const enc) {
  const int32_t delta = Clamp(original - predicted, -kMaxDelta, kMaxDelta);
  PutSValue(delta, kMaxDeltaBits, enc);
  return delta;
}

static inline int32_t ReadDelta(ANSDec* const dec) {
  return ReadSValue(kMaxDeltaBits, dec);
}

// For each of these vectors of indices, the last element will be predicted.
static uint32_t kPredictedVectors[3][3] = {{1, 4, 7}, {0, 1, 2}, {3, 4, 5}};

//------------------------------------------------------------------------------

// Encodes 'm' assuming that:
//  - all column- and row-vectors have the same magnitude,
//  - the first row and column contain only positive values,
//  - column- and row-vectors are mutually orthogonal.
static void WriteEigenMatrix(const int16_t m[9], ANSEnc* const enc) {
  ANSDebugPrefix prefix(enc, "eigen_mtx");
  //                                                  0 1 2
  // Final matrix that will be read. For reference:   3 4 5
  //                                                  6 7 8
  int16_t dm[9] = {0};

  for (uint32_t i : {0, 3, 6}) dm[i] = std::abs(m[i]);
  int64_t sq_len = GetVecSqLen(dm[0], dm[3], dm[6]);
  uint32_t num_sq_len = 1;

  // Find the maximum value of the middle top element based on top row.
  const int16_t max_dm1 = CompleteVector(dm[0], 0, sq_len) + kMaxDelta;
  dm[1] = std::min((int16_t)std::abs(m[1]), max_dm1);

  // Find the maximum value of the middle element based on middle row, column.
  const int16_t max_dm4 = RightShiftRound(CompleteVector(dm[1], 0, sq_len) +
                                              CompleteVector(dm[3], 0, sq_len),
                                          1) +
                          kMaxDelta;
  dm[4] = std::min((int16_t)std::abs(m[4]), max_dm4);

  // Find out and signal the maximum number of bits needed to encode 0,1,2,3,4.
  const int32_t max_dm = std::max({dm[0], dm[3], dm[6], dm[1], dm[4]});
  assert(max_dm < (1 << kEigenBits));
  const uint32_t num_bits = std::max(1u, (uint32_t)WP2Log2Ceil_k(max_dm + 1));
  enc->PutRange(num_bits, 1, kEigenBits, "num_bits");

  // Encode 0,3,6 as is and 1,4 as clamped.
  for (uint32_t i : {0, 3, 6}) enc->PutUValue(dm[i], num_bits, "value");
  (max_dm1 + 1 < (1 << num_bits)) ? enc->PutRValue(dm[1], max_dm1 + 1, "value")
                                  : enc->PutUValue(dm[1], num_bits, "value");
  (max_dm4 + 1 < (1 << num_bits)) ? enc->PutRValue(dm[4], max_dm4 + 1, "value")
                                  : enc->PutUValue(dm[4], num_bits, "value");

  // Other elements can be predicted; adjust by +-kMaxDelta to avoid any loss.
  for (const auto& vec : kPredictedVectors) {
    const uint32_t a = vec[0], b = vec[1], c = vec[2];  // Indices.
    dm[c] = CompleteVector(dm[a], dm[b], sq_len);
    dm[c] += PutDelta(std::abs(m[c]), dm[c], enc);
    // Improve 'sq_len' accuracy by taking into account the corrected values.
    sq_len = DivRound<int64_t>(
        sq_len * num_sq_len + GetVecSqLen(dm[a], dm[b], dm[c]), num_sq_len + 1);
    ++num_sq_len;
  }

  // The last one reuses the predictions and deltas so compute it after them.
  dm[8] = RightShiftRound((CompleteVector(dm[2], dm[5], sq_len) +
                           CompleteVector(dm[6], dm[7], sq_len)),
                          1);
  dm[8] += PutDelta(std::abs(m[8]), dm[8], enc);
}

static void ReadEigenMatrix(ANSDec* const dec, int16_t m[9]) {
  ANSDebugPrefix prefix(dec, "eigen_mtx");
  const int32_t num_bits = dec->ReadRange(1u, kEigenBits, "num_bits");

  for (uint32_t i : {0, 3, 6}) m[i] = dec->ReadUValue(num_bits, "value");
  int64_t sq_len = GetVecSqLen(m[0], m[3], m[6]);
  uint32_t num_sq_len = 1;

  const int16_t range_m1 = CompleteVector(m[0], 0, sq_len) + kMaxDelta + 1;
  m[1] = (range_m1 < (1 << num_bits)) ? dec->ReadRValue(range_m1, "value")
                                      : dec->ReadUValue(num_bits, "value");

  const int16_t range_m4 = RightShiftRound(CompleteVector(m[1], 0, sq_len) +
                                               CompleteVector(m[3], 0, sq_len),
                                           1) +
                           kMaxDelta + 1;
  m[4] = (range_m4 < (1 << num_bits)) ? dec->ReadRValue(range_m4, "value")
                                      : dec->ReadUValue(num_bits, "value");

  for (const auto& vec : kPredictedVectors) {
    const uint32_t a = vec[0], b = vec[1], c = vec[2];
    m[c] = CompleteVector(m[a], m[b], sq_len) + ReadDelta(dec);
    sq_len = DivRound<int64_t>(
        sq_len * num_sq_len + GetVecSqLen(m[a], m[b], m[c]), num_sq_len + 1);
    ++num_sq_len;
  }

  m[8] = RightShiftRound(
      CompleteVector(m[2], m[5], sq_len) + CompleteVector(m[6], m[7], sq_len),
      1);
  m[8] += ReadDelta(dec);

  // Find signs for the elements 4,5,7,8. As row-vectors are orthogonal, the
  // element generating the highest sub-dot-product has an opposite sign.
  for (uint32_t row : {3, 6}) {  // Middle and bottom vectors.
    int32_t subdot[3];           // Partial dot-product results with top vector.
    for (uint32_t col : {0, 1, 2}) subdot[col] = (int32_t)m[col] * m[col + row];
    const auto highest = std::max_element(subdot, subdot + 3) - subdot;
    for (int32_t col : {1, 2}) {
      if (highest == 0 || highest == col) m[col + row] = -m[col + row];
    }
  }
}

WP2Status CSPTransform::MakeEigenMatrixEncodable(const int16_t m[9],
                                                 int16_t dm[9],
                                                 int32_t* const error) {
  WP2MathInit();
  ANSEnc enc;
  WriteEigenMatrix(m, &enc);
  WP2_CHECK_STATUS(enc.Assemble());
  ExternalDataSource data(enc.Buffer(), enc.BufferSize());
  ANSDec dec(&data);
  ReadEigenMatrix(&dec, dm);
  WP2_CHECK_STATUS(dec.GetStatus());
  *error = 0;
  for (uint32_t i = 0; i < 9; ++i) {
    *error = std::max(*error, std::abs(m[i] - dm[i]));
  }
  return WP2_STATUS_OK;
}

//------------------------------------------------------------------------------

WP2Status CSPTransform::Write(ANSEnc* const enc)  const {
  enc->PutRValue((uint32_t)type_, kNumCspTypes, "csp_type");
  if (type_ == Csp::kCustom) {
    int32_t error;
    int16_t encodable_matrix[9];
    WP2_CHECK_STATUS(MakeEigenMatrixEncodable(im_, encodable_matrix, &error));
    if (enc->PutBool(error == 0, "eigen")) {
      WriteEigenMatrix(im_, enc);  // Matrix can be encoded with fewer bits.
    } else {
      for (auto& v : im_) PutSValue(v, kANSMaxUniformBits + 1, enc);
    }
    for (auto& v : avg_) enc->PutUValue(v, kRgbBits, "avg");
  }
  return WP2_STATUS_OK;
}

WP2Status CSPTransform::Read(ANSDec* const dec) {
  type_ = (Csp)dec->ReadRValue(kNumCspTypes, "csp_type");
  if (type_ == Csp::kCustom) {
    SetYUVBounds(-(1 << kMaxYuvBits), (1 << kMaxYuvBits) - 1);
    assert(yuv_bits_ <= kMaxYuvBits + 1);
    if (dec->ReadBool("eigen")) {
      ReadEigenMatrix(dec, im_);
    } else {
      // Completely custom matrices can be decoded.
      for (auto& v : im_) v = ReadSValue(kANSMaxUniformBits + 1, dec);
    }
    for (auto& v : avg_) v = (uint16_t)dec->ReadUValue(kRgbBits, "avg");
  } else {
    WP2_CHECK_STATUS(Init(type_));
  }
  WP2_CHECK_STATUS(dec->GetStatus());
  return WP2_STATUS_OK;
}

//------------------------------------------------------------------------------

namespace {

template<typename A> constexpr double Dot(const A a[3], const A b[3]) {
  return (double)a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

constexpr uint32_t kAvgBlk = 8;   // size of the averaging block

// returns the number of samples added to avg_sum[] average accumulator.
uint32_t ExtractCov(const uint8_t* argb, uint32_t w, uint32_t h, size_t stride,
                    double cov[3][3], double avg_sum[3]) {
  double avg[3] = { 0., 0., 0. };
  double tmp[3 * kAvgBlk * kAvgBlk];
  assert(w <= kAvgBlk && h <= kAvgBlk && w * h > 0);
  for (uint32_t idx = 0, j = 0; j < h; ++j, argb += stride) {
    for (uint32_t i = 0; i < w; ++i) {
      const double r = argb[4 * i + 1];
      const double g = argb[4 * i + 2];
      const double b = argb[4 * i + 3];
      avg[0] += r;
      avg[1] += g;
      avg[2] += b;
      tmp[idx++] = r;
      tmp[idx++] = g;
      tmp[idx++] = b;
    }
  }
  avg_sum[0] += avg[0];
  avg_sum[1] += avg[1];
  avg_sum[2] += avg[2];
  avg[0] /= w * h;
  avg[1] /= w * h;
  avg[2] /= w * h;
  for (uint32_t idx = 0, i = 0; i < w * h; ++i, idx += 3) {
    const double r = tmp[idx + 0] - avg[0];
    const double g = tmp[idx + 1] - avg[1];
    const double b = tmp[idx + 2] - avg[2];
    cov[0][0] += r * r;
    cov[0][1] += g * r;
    cov[0][2] += b * r;
    cov[1][0] += r * g;
    cov[1][1] += g * g;
    cov[1][2] += b * g;
    cov[2][0] += r * b;
    cov[2][1] += g * b;
    cov[2][2] += b * b;
  }
  return w * h;
}

bool ExtractCovarianceMtx(const uint8_t* argb,
                          uint32_t w, uint32_t h, size_t stride,
                          double cov[3][3], int16_t avg[3]) {
  if ((uint64_t)w * h <= 64 * 64) {
    return false;   // not enough sample to build reliable estimate!
  }
  double avg_sum[3];
  for (uint32_t j = 0; j < 3; ++j) {
    for (uint32_t i = 0; i < 3; ++i) cov[j][i] = 0.;
    avg_sum[j] = 0.;
  }
  uint32_t num_samples = 0;
  for (uint32_t j = 0; j < h; j += kAvgBlk, argb += kAvgBlk * stride) {
    for (uint32_t i = 0; i < w; i += kAvgBlk) {
      const uint32_t W = std::min(i + kAvgBlk, w) - i;
      const uint32_t H = std::min(j + kAvgBlk, h) - j;
      num_samples += ExtractCov(argb + i, W, H, stride, cov, avg_sum);
    }
  }
  // Normalize a bit the covariance matrix (not totally necessary, just comfort)
  for (uint32_t j = 0; j < 3; ++j) {
    for (uint32_t i = 0; i < 3; ++i) cov[j][i] /= num_samples;
  }
  // store the average
  for (uint32_t i = 0; i < 3; ++i) {
    if (i == 0) {
      // We only centralize the 'luma' channel. Empirically, the U/V channels
      // are most unstable to centralization and renormalization (during later
      // calls NormalizeFromInput()).
      avg_sum[i] /= num_samples;
      avg[i] = (int16_t)Clamp<double>(avg_sum[i], 0, kRgbMax);
    } else {
      avg[i] = (1u << (kRgbBits - 1));
    }
  }
  return true;
}

bool ExtractEigenMtx(double cov[3][3], double out[9]) {
  // Good seeding matrix for starting the iteration. We expect the
  // final eigenvectors to have a similar signature.
  double v[3][3] = {
    { .5,  0.0,  0.5 },
    { .5,  0.5, -0.5 },
    { .5, -0.5,  0.5 }
  };
  double lambdas[3] = { 1., 1., 1. };
  constexpr double kStopThreshold = 0.001;
  for (uint32_t n = 0; n < 3; ++n) {
    double lambda;
    const uint32_t num_iters = 20;
    uint32_t k = 0;
    for (; k < num_iters; ++k) {
      lambda = 0.;
      double vv[3];
      for (uint32_t i = 0; i < 3; ++i) {
        vv[i] = Dot(cov[i], v[n]);
        lambda += vv[i] * vv[i];
      }
      if (lambda < kStopThreshold) break;
      lambda = std::sqrt(lambda);
      double lnorm = (lambda > 0.) ? 1. / lambda : 1.;
      if (vv[0] < 0.) lnorm = -lnorm;
      double epsilon = 0.;
      for (uint32_t i = 0; i < 3; ++i) {
        vv[i] *= lnorm;
        epsilon += std::abs(v[n][i] - vv[i]);
        v[n][i] = vv[i];
      }
      if (epsilon < kStopThreshold) break;
    }
    if (n == 0) {
      // for the partitioning algo (based on variance) to work correctly, we
      // need a primary axis that is close to a 'luma' one.
      const double kLumaAxis[3] = {.299, .587, .114};
      if (k == num_iters || Distance(v[0], kLumaAxis) > 0.4f) {
        return false;  // didn't converge on luma -> probably problematic
      }
    }

    if (lambda < kStopThreshold) {  // can't converge -> force vector(s)
      if (n == 0) {
        return false;       // very bad case
      } else if (n == 1) {  // monochrome input?
        const double kAxis[3] = {0., 1., -1.};
        CrossProduct(kAxis, v[0], v[1]);
      } else {  // bichrome input?
        assert(n == 2);
        CrossProduct(v[0], v[1], v[2]);
        if (v[2][0] < 0.) {
          for (uint32_t i : {0, 1, 2}) v[2][i] = -v[2][i];
        }
      }
    }
    // store the column vector in the inverse matrix
    out[n + 0] = v[n][0];
    out[n + 3] = v[n][1];
    out[n + 6] = v[n][2];
    lambdas[n] = lambda;

    // remove eigen vector from matrix
    if (n < 2) {
      for (uint32_t i = 0; i < 3; ++i) {
        const double c = Dot(cov[i], v[n]);
        for (uint32_t j = 0; j < 3; ++j) cov[i][j] -= c * v[n][j];
      }
    }
  }
  // very simple and conservative classifier
  bool is_good =
      (lambdas[0] >= 2500.f && lambdas[1] >= 20. && lambdas[2] >= 1.) ||
      (lambdas[0] <= 5000. && lambdas[1] >= 300. * (1. - lambdas[0] / 5000.));
  is_good &= (lambdas[1] < 1000.);
  // TODO(skal): we need ML and NN here!
  return is_good;
}

bool NormalizeFromInput(const uint8_t* argb,
                        uint32_t w, uint32_t h, size_t stride,
                        double max_value,
                        const int16_t avg[3],
                        double custom_to_rgb_matrix[9],
                        int16_t out[9]) {
  double rgb_to_custom_matrix[9];
  if (!TryInvertMatrix(custom_to_rgb_matrix, rgb_to_custom_matrix)) {
    return false;
  }

  double max_custom = 1.;
  for (uint32_t y = 0; y < h; ++y) {
    for (uint32_t x = 0; x < w; ++x) {
      const double rgb[3] = {(double)argb[4 * x + 1] - avg[0],
                             (double)argb[4 * x + 2] - avg[1],
                             (double)argb[4 * x + 3] - avg[2]};
      const double X = std::abs(Dot(rgb_to_custom_matrix + 0, rgb));
      const double Y = std::abs(Dot(rgb_to_custom_matrix + 3, rgb));
      const double Z = std::abs(Dot(rgb_to_custom_matrix + 6, rgb));
      max_custom = std::max(max_custom, std::max({X, Y, Z}));
    }
    argb += stride;
  }

  // Scale rgb_to_custom_matrix so it goes from ~kRgbBits to yuv_bits_ (so it
  // uses the whole yuv_bits_ range).
  max_value = max_value / max_custom;
  for (uint32_t i = 0; i < 9; ++i) rgb_to_custom_matrix[i] *= max_value;

  // Compute custom_to_rgb_matrix, prepare it for fixed point.
  if (!TryInvertMatrix(rgb_to_custom_matrix, custom_to_rgb_matrix)) {
    return false;
  }
  MultiplyMatrix(custom_to_rgb_matrix, kMtxFixedPointMul);

  constexpr double kMin = 10., kMax = (1 << kEigenBits) - 1.;
  const double max =
      std::max(std::abs(*std::min_element(custom_to_rgb_matrix,
                                          custom_to_rgb_matrix + 9)),
               std::abs(*std::max_element(custom_to_rgb_matrix,
                                          custom_to_rgb_matrix + 9)));
  if (max < kMin || max > kMax) return false;

  for (uint32_t i = 0; i < 9; ++i) {
    out[i] = (int16_t)lround(custom_to_rgb_matrix[i]);
  }
  return true;
}

}    // anonymous namespace

//------------------------------------------------------------------------------

WP2Status CSPTransform::Optimize(const ArgbBuffer& image) {
  WP2_CHECK_OK(image.format == WP2_Argb_32, WP2_STATUS_INVALID_PARAMETER);
  WP2_CHECK_OK(!image.IsEmpty(), WP2_STATUS_INVALID_PARAMETER);
  const uint8_t* const argb = (const uint8_t*)image.GetRow(0);
  const uint32_t w = image.width;
  const uint32_t h = image.height;
  const size_t stride = image.stride;
  SetYUVBounds(-(1 << kMaxYuvBits), (1 << kMaxYuvBits) - 1);

  double cov[3][3];
  double custom_to_rgb_matrix[9];
  int16_t im[9];
  if (ExtractCovarianceMtx(argb, w, h, stride, cov, avg_) &&
      ExtractEigenMtx(cov, custom_to_rgb_matrix) &&
      NormalizeFromInput(
          argb, w, h, stride,
          std::max((int16_t)std::abs(min_yuv_value_), max_yuv_value_), avg_,
          custom_to_rgb_matrix, im)) {
    int32_t error;
    WP2_CHECK_STATUS(MakeEigenMatrixEncodable(im, im_, &error));
    if (error > 0) {
      // 'kMaxDelta' should be big enough to correct most errors, though it may
      // happen that some ExtractEigenMtx() outputs are too imprecise. In this
      // case, encode those as is.
      std::copy(im, im + 9, im_);
      // Reducing 'kMaxDelta' and keeping the "quantized" matrix (assuming
      // 'error' is reasonable) is possible but it leads to poor results.
    }
    if (ComputeRGBToYUVMatrix()) {
      type_ = Csp::kCustom;
      assert(yuv_bits_ <= kMaxYuvBits + 1);
      return WP2_STATUS_OK;
    }
  }
  return WP2_STATUS_INVALID_COLORSPACE;
}

//------------------------------------------------------------------------------

void CSPTransform::Apply(ArgbBuffer* const image, uint32_t num_precision_bits,
                         uint32_t x, uint32_t y, uint32_t w, uint32_t h) const {
  if (image == nullptr || image->IsEmpty()) return;
  for (uint32_t j = 0; j < h; ++j) {
    if (y + j >= image->height) break;
    uint8_t* const dst = (uint8_t*)image->GetRow(y + j);
    for (uint32_t i = x; i < std::min(x + w, image->width); ++i) {
      // Set alpha to opaque.
      dst[4 * i + 0] = 255;
      int16_t yuv[3];
      ToYUV(&dst[4 * i + 1], yuv);
      for (uint32_t c = 0; c < 3; ++c) {
        const int32_t max_yuv_value = (1 << (num_precision_bits - 1));
        // Convert to 8b-unsigned.
        dst[4 * i + 1 + c] =
            (uint8_t)Clamp(DivRound(yuv[c] * 128, max_yuv_value) + 128, 0, 255);
      }
    }
  }
}

//------------------------------------------------------------------------------

}  // namespace WP2
