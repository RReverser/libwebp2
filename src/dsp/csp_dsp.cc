// Copyright 2020 Google LLC
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
// Custom-CSP critical functions
//
// Author: Skal (pascal.massimino@gmail.com)

#include <cassert>

#include "src/dsp/dsp.h"
#include "src/dsp/math.h"
#include "src/wp2/format_constants.h"

namespace {

constexpr uint32_t kMtxShift = 12u;

//------------------------------------------------------------------------------

void YUVToCustom_C(const int16_t* y, const int16_t* u, const int16_t* v,
                   const int16_t offset[3], const int16_t mtx[9],
                   int16_t* dst0, int16_t* dst1, int16_t* dst2,
                   uint32_t width) {
  for (uint32_t x = 0; x < width; ++x) {
    // TODO(yguyon): No overflow guaranteed only if custom ccsp is kMaxYuvBits
    int32_t ccsp[3];
    WP2::Multiply(y[x] + offset[0], u[x] + offset[1], v[x] + offset[2], mtx,
                  &ccsp[0], &ccsp[1], &ccsp[2]);
    dst0[x] = (int16_t)WP2::RightShiftRound(ccsp[0], kMtxShift);
    dst1[x] = (int16_t)WP2::RightShiftRound(ccsp[1], kMtxShift);
    dst2[x] = (int16_t)WP2::RightShiftRound(ccsp[2], kMtxShift);
  }
}

//------------------------------------------------------------------------------


inline void YUVToXRGB(const int16_t avg[3], const int16_t im[9],
                      int16_t y, int16_t u, int16_t v,
                      uint8_t* const argb) {
  int32_t tmp_r, tmp_g, tmp_b;
  WP2::Multiply<int32_t>(y, u, v, im, &tmp_r, &tmp_g, &tmp_b);
  tmp_r = avg[0] + WP2::RightShiftRound(tmp_r, kMtxShift);
  tmp_g = avg[1] + WP2::RightShiftRound(tmp_g, kMtxShift);
  tmp_b = avg[2] + WP2::RightShiftRound(tmp_b, kMtxShift);
  argb[0] = WP2::kAlphaMax;
  argb[1] = (uint8_t)WP2::Clamp<int32_t>(tmp_r, 0, 255);
  argb[2] = (uint8_t)WP2::Clamp<int32_t>(tmp_g, 0, 255);
  argb[3] = (uint8_t)WP2::Clamp<int32_t>(tmp_b, 0, 255);
}

inline void YUVToARGB(const int16_t avg[3], const int16_t im[9],
                      int16_t y, int16_t u, int16_t v, int32_t a,
                      uint8_t* const argb) {
  int32_t tmp_r, tmp_g, tmp_b;
  WP2::Multiply<int32_t>(y, u, v, im, &tmp_r, &tmp_g, &tmp_b);
  tmp_r = avg[0] + WP2::RightShiftRound(tmp_r, kMtxShift);
  tmp_g = avg[1] + WP2::RightShiftRound(tmp_g, kMtxShift);
  tmp_b = avg[2] + WP2::RightShiftRound(tmp_b, kMtxShift);
  a = WP2::Clamp<int32_t>(a, 0, 255);
  argb[0] = (uint8_t)a;
  if (a == 0) {
    argb[1] = argb[2] = argb[3] = 0;
  } else {
    argb[1] = (uint8_t)WP2::Clamp<int32_t>(tmp_r, 0, a);
    argb[2] = (uint8_t)WP2::Clamp<int32_t>(tmp_g, 0, a);
    argb[3] = (uint8_t)WP2::Clamp<int32_t>(tmp_b, 0, a);
    const uint32_t M = WP2::kAlphaDiv[a];
    argb[1] = WP2::DivByAlphaDiv(argb[1], M);
    argb[2] = WP2::DivByAlphaDiv(argb[2], M);
    argb[3] = WP2::DivByAlphaDiv(argb[3], M);
  }
}

inline void YUVToArgb(const int16_t avg[3], const int16_t im[9],
                      int16_t y, int16_t u, int16_t v, int32_t a,
                      uint8_t* const argb) {
  int32_t tmp_r, tmp_g, tmp_b;
  WP2::Multiply<int32_t>(y, u, v, im, &tmp_r, &tmp_g, &tmp_b);
  tmp_r = avg[0] + WP2::RightShiftRound(tmp_r, kMtxShift);
  tmp_g = avg[1] + WP2::RightShiftRound(tmp_g, kMtxShift);
  tmp_b = avg[2] + WP2::RightShiftRound(tmp_b, kMtxShift);
  a = WP2::Clamp<int32_t>(a, 0, 255);
  argb[0] = (uint8_t)a;
  if (a > 0) {
    argb[1] = (uint8_t)WP2::Clamp<int32_t>(tmp_r, 0, a);
    argb[2] = (uint8_t)WP2::Clamp<int32_t>(tmp_g, 0, a);
    argb[3] = (uint8_t)WP2::Clamp<int32_t>(tmp_b, 0, a);
  } else {
    argb[1] = 0;
    argb[2] = 0;
    argb[3] = 0;
  }
}

void YUVToXRGB_C(const int16_t* y, const int16_t* u, const int16_t* v,
                 const int16_t* a, const int16_t avg[3], const int16_t mtx[9],
                 uint8_t* argb, uint32_t width) {
  (void)a;
  for (uint32_t x = 0; x < width; ++x) {
    YUVToXRGB(avg, mtx, y[x], u[x], v[x], argb + 4 * x);
  }
}

void YUVToArgb_C(const int16_t* y, const int16_t* u, const int16_t* v,
                 const int16_t* a, const int16_t avg[3], const int16_t mtx[9],
                 uint8_t* argb, uint32_t width) {
  for (uint32_t x = 0; x < width; ++x) {
    YUVToArgb(avg, mtx, y[x], u[x], v[x], a[x], argb + 4 * x);
  }
}

void YUVToARGB_C(const int16_t* y, const int16_t* u, const int16_t* v,
                 const int16_t* a, const int16_t avg[3], const int16_t mtx[9],
                 uint8_t* argb, uint32_t width) {
  for (uint32_t x = 0; x < width; ++x) {
    YUVToARGB(avg, mtx, y[x], u[x], v[x], a[x], argb + 4 * x);
  }
}

//------------------------------------------------------------------------------
// SSE4.1 implementation

#if defined(WP2_USE_SSE)

//------------------------------------------------------------------------------
// some useful macros used extensively below

// load the matrix elements mtx[] into variables mtx00, ... mtx21
#define MAKE32b(A, B)  ((uint32_t)((A) & 0xffff) | ((uint32_t)(B) << 16))

#define LOAD_MTX_ELEMENTS(MTX)                                    \
  constexpr uint32_t kRnd = 1u << kMtxShift >> 1;                 \
  const __m128i mtx00 = _mm_set1_epi32(MAKE32b(MTX[0], MTX[1]));  \
  const __m128i mtx01 = _mm_set1_epi32(MAKE32b(MTX[2], kRnd));    \
  const __m128i mtx10 = _mm_set1_epi32(MAKE32b(MTX[3], MTX[4]));  \
  const __m128i mtx11 = _mm_set1_epi32(MAKE32b(MTX[5], kRnd));    \
  const __m128i mtx20 = _mm_set1_epi32(MAKE32b(MTX[6], MTX[7]));  \
  const __m128i mtx21 = _mm_set1_epi32(MAKE32b(MTX[8], kRnd));    \
  const __m128i kCst = _mm_set1_epi16(1)

// load the y[x],u[x],v[x] samples in variables Y0, U0, V0
#define LOAD_YUV(Y0, U0, V0)                                     \
  const __m128i Y0 = _mm_loadu_si128((const __m128i*)(y + x));   \
  const __m128i U0 = _mm_loadu_si128((const __m128i*)(u + x));   \
  const __m128i V0 = _mm_loadu_si128((const __m128i*)(v + x))

// unpack Y/U/V to variables 16b-pairs Y|U and V|K
#define UNPACK_YUV(Y, U, V)                          \
  const __m128i YU_lo = _mm_unpacklo_epi16(Y, U);    \
  const __m128i YU_hi = _mm_unpackhi_epi16(Y, U);    \
  const __m128i VK_lo = _mm_unpacklo_epi16(V, kCst); \
  const __m128i VK_hi = _mm_unpackhi_epi16(V, kCst)

#define MULT_ROW(MTX0, MTX1, OUT)                      \
  __m128i OUT;                                         \
  do {                                                 \
    const __m128i a0 = _mm_madd_epi16(YU_lo, MTX0);    \
    const __m128i a1 = _mm_madd_epi16(YU_hi, MTX0);    \
    const __m128i b0 = _mm_madd_epi16(VK_lo, MTX1);    \
    const __m128i b1 = _mm_madd_epi16(VK_hi, MTX1);    \
    const __m128i c0 = _mm_add_epi32(a0, b0);          \
    const __m128i c1 = _mm_add_epi32(a1, b1);          \
    const __m128i d0 = _mm_srai_epi32(c0, kMtxShift);  \
    const __m128i d1 = _mm_srai_epi32(c1, kMtxShift);  \
    OUT = _mm_packs_epi32(d0, d1);                     \
  } while (false)

#define PACK_RGB(A, R, G, B, DST) do {                       \
  const __m128i out_rb = _mm_packus_epi16(R, B);             \
  const __m128i out_ag = _mm_packus_epi16(A, G);             \
  const __m128i out_ar = _mm_unpacklo_epi8(out_ag, out_rb);  \
  const __m128i out_gb = _mm_unpackhi_epi8(out_ag, out_rb);  \
  const __m128i argb0 = _mm_unpacklo_epi16(out_ar, out_gb);  \
  const __m128i argb1 = _mm_unpackhi_epi16(out_ar, out_gb);  \
  _mm_storeu_si128((__m128i*)(DST +  0), argb0);             \
  _mm_storeu_si128((__m128i*)(DST + 16), argb1);             \
} while (false)

//------------------------------------------------------------------------------

void YUVToCustom_SSE(const int16_t* y, const int16_t* u, const int16_t* v,
                     const int16_t offset[3], const int16_t mtx[9],
                     int16_t* dst0, int16_t* dst1, int16_t* dst2,
                     uint32_t width) {
  const __m128i Off0 = _mm_set1_epi16(offset[0]);
  const __m128i Off1 = _mm_set1_epi16(offset[1]);
  const __m128i Off2 = _mm_set1_epi16(offset[2]);
  LOAD_MTX_ELEMENTS(mtx);

  uint32_t x = 0;
  for (; x + 8 <= width; x += 8) {
    // load source
    LOAD_YUV(Y0, U0, V0);
    // add offset
    const __m128i Y1 = _mm_adds_epi16(Y0, Off0);
    const __m128i U1 = _mm_adds_epi16(U0, Off1);
    const __m128i V1 = _mm_adds_epi16(V0, Off2);
    // unpack to YU + VK
    UNPACK_YUV(Y1, U1, V1);
    // multiply with mtx
    MULT_ROW(mtx00, mtx01, out0);
    MULT_ROW(mtx10, mtx11, out1);
    MULT_ROW(mtx20, mtx21, out2);
    // and store
    _mm_storeu_si128((__m128i*)(dst0 + x), out0);
    _mm_storeu_si128((__m128i*)(dst1 + x), out1);
    _mm_storeu_si128((__m128i*)(dst2 + x), out2);
  }

  if (x < width) {  // left-over
    YUVToCustom_C(y + x, u + x, v + x, offset, mtx,
                  dst0 + x, dst1 + x, dst2 + x, width - x);
  }
}

//------------------------------------------------------------------------------

void YUVToXRGB_SSE(const int16_t* y, const int16_t* u, const int16_t* v,
                   const int16_t* a, const int16_t avg[3], const int16_t mtx[9],
                   uint8_t* argb, uint32_t width) {
  (void)a;
  const __m128i Avg0 = _mm_set1_epi16(avg[0]);
  const __m128i Avg1 = _mm_set1_epi16(avg[1]);
  const __m128i Avg2 = _mm_set1_epi16(avg[2]);
  LOAD_MTX_ELEMENTS(mtx);
  const __m128i alpha = _mm_set1_epi16(255u);
  uint32_t x = 0;
  for (; x + 8 <= width; x += 8) {
    // load source
    LOAD_YUV(Y0, U0, V0);
    UNPACK_YUV(Y0, U0, V0);

    // multiply with mtx and add average
    MULT_ROW(mtx00, mtx01, R0);
    MULT_ROW(mtx10, mtx11, G0);
    MULT_ROW(mtx20, mtx21, B0);

    const __m128i R1 = _mm_adds_epi16(Avg0, R0);
    const __m128i G1 = _mm_adds_epi16(Avg1, G0);
    const __m128i B1 = _mm_adds_epi16(Avg2, B0);

    // repack to 8b, with alpha
    PACK_RGB(alpha, R1, G1, B1, argb + 4 * x);
  }

  if (x < width) {  // left-over
    YUVToXRGB_C(y + x, u + x, v + x, nullptr, avg, mtx,
                argb + 4 * x, width - x);
  }
}

//------------------------------------------------------------------------------

void YUVToArgb_SSE(const int16_t* y, const int16_t* u, const int16_t* v,
                   const int16_t* a, const int16_t avg[3], const int16_t mtx[9],
                   uint8_t* argb, uint32_t width) {
  const __m128i Avg0 = _mm_set1_epi16(avg[0]);
  const __m128i Avg1 = _mm_set1_epi16(avg[1]);
  const __m128i Avg2 = _mm_set1_epi16(avg[2]);
  LOAD_MTX_ELEMENTS(mtx);
  uint32_t x = 0;
  for (; x + 8 <= width; x += 8) {
    LOAD_YUV(Y0, U0, V0);
    UNPACK_YUV(Y0, U0, V0);
    MULT_ROW(mtx00, mtx01, R0);
    MULT_ROW(mtx10, mtx11, G0);
    MULT_ROW(mtx20, mtx21, B0);
    const __m128i A1 = _mm_loadu_si128((const __m128i*)(a + x));
    const __m128i R1 = _mm_min_epi16(A1, _mm_adds_epi16(Avg0, R0));
    const __m128i G1 = _mm_min_epi16(A1, _mm_adds_epi16(Avg1, G0));
    const __m128i B1 = _mm_min_epi16(A1, _mm_adds_epi16(Avg2, B0));
    PACK_RGB(A1, R1, G1, B1, argb + 4 * x);
  }

  if (x < width) {  // left-over
    YUVToArgb_C(y + x, u + x, v + x, a + x, avg, mtx, argb + 4 * x, width - x);
  }
}

void YUVToARGB_SSE(const int16_t* y, const int16_t* u, const int16_t* v,
                   const int16_t* a, const int16_t avg[3], const int16_t mtx[9],
                   uint8_t* argb, uint32_t width) {
  const __m128i Avg0 = _mm_set1_epi16(avg[0]);
  const __m128i Avg1 = _mm_set1_epi16(avg[1]);
  const __m128i Avg2 = _mm_set1_epi16(avg[2]);
  LOAD_MTX_ELEMENTS(mtx);
  uint32_t x = 0;
  for (; x + 8 <= width; x += 8) {
    LOAD_YUV(Y0, U0, V0);
    UNPACK_YUV(Y0, U0, V0);
    MULT_ROW(mtx00, mtx01, R0);
    MULT_ROW(mtx10, mtx11, G0);
    MULT_ROW(mtx20, mtx21, B0);
    const __m128i A1 = _mm_loadu_si128((const __m128i*)(a + x));
    const __m128i R1 = _mm_min_epi16(A1, _mm_adds_epi16(Avg0, R0));
    const __m128i G1 = _mm_min_epi16(A1, _mm_adds_epi16(Avg1, G0));
    const __m128i B1 = _mm_min_epi16(A1, _mm_adds_epi16(Avg2, B0));
    uint8_t tmp[4 * 8];
    PACK_RGB(A1, R1, G1, B1, tmp);
    for (uint32_t i = 0; i < 8; ++i) {
      const uint8_t A = tmp[4 * i + 0];
      if (A == 255) {
        argb[4 * (x + i) + 0] = tmp[4 * i + 0];
        argb[4 * (x + i) + 1] = tmp[4 * i + 1];
        argb[4 * (x + i) + 2] = tmp[4 * i + 2];
        argb[4 * (x + i) + 3] = tmp[4 * i + 3];
      } else  {
        const uint32_t M = WP2::kAlphaDiv[A];
        argb[4 * (x + i) + 0] = A;
        argb[4 * (x + i) + 1] = WP2::DivByAlphaDiv(tmp[4 * i + 1], M);
        argb[4 * (x + i) + 2] = WP2::DivByAlphaDiv(tmp[4 * i + 2], M);
        argb[4 * (x + i) + 3] = WP2::DivByAlphaDiv(tmp[4 * i + 3], M);
      }
    }
  }

  if (x < width) {  // left-over
    YUVToARGB_C(y + x, u + x, v + x, a + x, avg, mtx, argb + 4 * x, width - x);
  }
}

//------------------------------------------------------------------------------

#undef LOAD_MTX_ELEMENTS
#undef LOAD_YUV
#undef UNPACK_YUV
#undef MULT_ROW
#undef PACK_RGB

WP2_TSAN_IGNORE_FUNCTION void WP2CSPConverterInitSSE() {
  WP2YUVToCustom = YUVToCustom_SSE;

  WP2YUVToXRGB = YUVToXRGB_SSE;
  WP2YUVToArgb = YUVToArgb_SSE;
  WP2YUVToARGB = YUVToARGB_SSE;
}

#endif  // WP2_USE_SSE

//------------------------------------------------------------------------------
// NEON implementation (TODO)

//------------------------------------------------------------------------------

}  // namespace

void (*WP2YUVToCustom)(const int16_t* y, const int16_t* u, const int16_t* v,
                       const int16_t offset[3], const int16_t mtx[9],
                       int16_t* dst0, int16_t* dst1, int16_t* dst2,
                       uint32_t width) = nullptr;

WP2YUVToArgbFunc WP2YUVToXRGB = nullptr;
WP2YUVToArgbFunc WP2YUVToArgb = nullptr;
WP2YUVToArgbFunc WP2YUVToARGB = nullptr;

static volatile WP2CPUInfo csp_converter_last_cpuinfo_used =
    (WP2CPUInfo)&csp_converter_last_cpuinfo_used;

WP2_TSAN_IGNORE_FUNCTION void WP2CSPConverterInit() {
  if (csp_converter_last_cpuinfo_used == WP2GetCPUInfo) return;

  WP2YUVToCustom = YUVToCustom_C;
  WP2YUVToXRGB = YUVToXRGB_C;
  WP2YUVToArgb = YUVToArgb_C;
  WP2YUVToARGB = YUVToARGB_C;

  if (WP2GetCPUInfo != nullptr) {
#if defined(WP2_USE_SSE)
    if (WP2GetCPUInfo(kSSE)) WP2CSPConverterInitSSE();
#endif
  }
  csp_converter_last_cpuinfo_used = WP2GetCPUInfo;
}

//------------------------------------------------------------------------------
