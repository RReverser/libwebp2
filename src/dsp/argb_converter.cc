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
// Speed-critical convert-to-Argb functions
//
// Author: Skal (pascal.massimino@gmail.com)

#include <cstring>

#include "src/dsp/dsp.h"
#include "src/dsp/math.h"

namespace {

#define MAKE_COLOR(A, R, G, B, DST) do {                           \
  (DST)[0] = (A); (DST)[1] = (R); (DST)[2] = (G); (DST)[3] = (B);  \
} while (0)

#define MAKE_ARGB(S, D) MAKE_COLOR((S)[0], (S)[1], (S)[2], (S)[3], (D))
#define MAKE_XRGB(S, D) MAKE_COLOR(0xff,   (S)[1], (S)[2], (S)[3], (D))
#define MAKE_RGBA(S, D) MAKE_COLOR((S)[3], (S)[0], (S)[1], (S)[2], (D))
#define MAKE_RGBX(S, D) MAKE_COLOR(0xff,   (S)[0], (S)[1], (S)[2], (D))
#define MAKE_BGRA(S, D) MAKE_COLOR((S)[3], (S)[2], (S)[1], (S)[0], (D))
#define MAKE_BGRX(S, D) MAKE_COLOR(0xff,   (S)[2], (S)[1], (S)[0], (D))

//------------------------------------------------------------------------------
// ConvertFrom...

void Premultiply(uint8_t* ARGB, uint32_t width) {
  uint32_t i;
  for (i = 0; i < width; ++i, ARGB += 4) {
    const uint32_t A = ARGB[0];
    if (A < 255) {
      ARGB[1] = WP2::DivBy255(ARGB[1] * A);
      ARGB[2] = WP2::DivBy255(ARGB[2] * A);
      ARGB[3] = WP2::DivBy255(ARGB[3] * A);
    }
  }
}
#define CONVERT_FROM_FUNC(NAME, CALL)                                        \
void NAME(const uint8_t* src, uint32_t width, uint8_t* dst) {                \
  for (uint32_t i = 0; i < width; ++i, src += 4, dst += 4) CALL(src, dst);   \
}

CONVERT_FROM_FUNC(ConvertFrom_Argb_C, MAKE_ARGB)
CONVERT_FROM_FUNC(ConvertFrom_XRGB_C, MAKE_XRGB)
void ConvertFrom_ARGB_C(const uint8_t* src, uint32_t width, uint8_t* dst) {
  ConvertFrom_Argb_C(src, width, dst);
  Premultiply(dst, width);
}

CONVERT_FROM_FUNC(ConvertFrom_rgbA_C, MAKE_RGBA)
CONVERT_FROM_FUNC(ConvertFrom_RGBX_C, MAKE_RGBX)
void ConvertFrom_RGBA_C(const uint8_t* src, uint32_t width, uint8_t* dst) {
  ConvertFrom_rgbA_C(src, width, dst);
  Premultiply(dst, width);
}

CONVERT_FROM_FUNC(ConvertFrom_bgrA_C, MAKE_BGRA)
CONVERT_FROM_FUNC(ConvertFrom_BGRX_C, MAKE_BGRX)
void ConvertFrom_BGRA_C(const uint8_t* src, uint32_t width, uint8_t* dst) {
  ConvertFrom_bgrA_C(src, width, dst);
  Premultiply(dst, width);
}

void ConvertFrom_RGB_C(const uint8_t* src, uint32_t width, uint8_t* dst) {
  uint32_t i;
  for (i = 0; i < width; ++i, src += 3, dst += 4) {
    dst[0] = 0xff;
    dst[1] = src[0];
    dst[2] = src[1];
    dst[3] = src[2];
  }
}

void ConvertFrom_BGR_C(const uint8_t* src, uint32_t width, uint8_t* dst) {
  uint32_t i;
  for (i = 0; i < width; ++i, src += 3, dst += 4) {
    dst[0] = 0xff;
    dst[1] = src[2];
    dst[2] = src[1];
    dst[3] = src[0];
  }
}

#undef CONVERT_FROM_FUNC

//------------------------------------------------------------------------------
// ConvertTo...

void Unmult(uint32_t a, uint8_t* const v) {
  if (a > 0) {
    const uint32_t M = WP2::kAlphaDiv[a];
    v[0] = std::min(WP2::DivByAlphaDiv(v[0], M), 255u);
    v[1] = std::min(WP2::DivByAlphaDiv(v[1], M), 255u);
    v[2] = std::min(WP2::DivByAlphaDiv(v[2], M), 255u);
  } else {
    v[0] = v[1] = v[2] = 0;
  }
}

#define CONVERT_TO_FUNC(NAME, A, R, G, B, UNMULT, FILL)                 \
void NAME(const uint8_t* src, uint32_t width, uint8_t* dst) {           \
  for (uint32_t i = 0; i < width; ++i, src += 4, dst += 4) {            \
    dst[(A)] = (FILL) ? 0xff : src[0];                                  \
    dst[(R)] = src[1];                                                  \
    dst[(G)] = src[2];                                                  \
    dst[(B)] = src[3];                                                  \
    if (!(FILL) && (UNMULT) >= 0) Unmult(dst[(A)], dst + (UNMULT));     \
  }                                                                     \
}

CONVERT_TO_FUNC(ConvertTo_Argb_C, 0, 1, 2, 3, -1, 0);
CONVERT_TO_FUNC(ConvertTo_ARGB_C, 0, 1, 2, 3, 1, 0);
CONVERT_TO_FUNC(ConvertTo_XRGB_C, 0, 1, 2, 3, 0, 1);
CONVERT_TO_FUNC(ConvertTo_rgbA_C, 3, 0, 1, 2, -1, 0);
CONVERT_TO_FUNC(ConvertTo_RGBA_C, 3, 0, 1, 2, 0, 0);
CONVERT_TO_FUNC(ConvertTo_RGBX_C, 3, 0, 1, 2, 0, 1);
CONVERT_TO_FUNC(ConvertTo_bgrA_C, 3, 2, 1, 0, -1, 0);
CONVERT_TO_FUNC(ConvertTo_BGRA_C, 3, 2, 1, 0, 0, 0);
CONVERT_TO_FUNC(ConvertTo_BGRX_C, 3, 2, 1, 0, 0, 1);

void ConvertTo_RGB_C(const uint8_t* src, uint32_t width, uint8_t* dst) {
  for (uint32_t i = 0; i < width; ++i, src += 4, dst += 3) {
    dst[0] = src[1];
    dst[1] = src[2];
    dst[2] = src[3];
    Unmult(src[0], dst);
  }
}

void ConvertTo_BGR_C(const uint8_t* src, uint32_t width, uint8_t* dst) {
  for (uint32_t i = 0; i < width; ++i, src += 4, dst += 3) {
    dst[2] = src[1];
    dst[1] = src[2];
    dst[0] = src[3];
    Unmult(src[0], dst);
  }
}

#undef CONVERT_TO_FUNC

//------------------------------------------------------------------------------
// SSE4.1 implementation

#if defined(WP2_USE_SSE)

#define PSHUFB_CST32(CST)  \
  _mm_set_epi32(CST | 0x0c0c0c0c, CST | 0x08080808, CST | 0x04040404, CST)

#define CONVERT_FUNC(NAME, CST)                                              \
void NAME##_SSE(const uint8_t* src, uint32_t width, uint8_t* dst) {          \
  const __m128i kShuffle = PSHUFB_CST32(CST);                                \
  uint32_t x;                                                                \
  for (x = 0; x + 8 <= width; x += 8) {                                      \
    const __m128i A0 = _mm_loadu_si128((const __m128i*)(src + 4 * x +  0));  \
    const __m128i A1 = _mm_loadu_si128((const __m128i*)(src + 4 * x + 16));  \
    const __m128i B0 = _mm_shuffle_epi8(A0, kShuffle);                       \
    const __m128i B1 = _mm_shuffle_epi8(A1, kShuffle);                       \
    _mm_storeu_si128((__m128i*)(dst + 4 * x +  0), B0);                      \
    _mm_storeu_si128((__m128i*)(dst + 4 * x + 16), B1);                      \
  }                                                                          \
  if (x < width) NAME##_C(src + 4 * x, width - x, dst + 4 * x);              \
}

#define CONVERT_FUNC_FILL(NAME, CST, CST_A)                                  \
void NAME##_SSE(const uint8_t* src, uint32_t width, uint8_t* dst) {          \
  const __m128i kShuffle = PSHUFB_CST32(CST);                                \
  const __m128i kFill = _mm_set1_epi32(CST_A);                               \
  uint32_t x;                                                                \
  for (x = 0; x + 8 <= width; x += 8) {                                      \
    const __m128i A0 = _mm_loadu_si128((const __m128i*)(src + 4 * x +  0));  \
    const __m128i A1 = _mm_loadu_si128((const __m128i*)(src + 4 * x + 16));  \
    const __m128i B0 = _mm_or_si128(_mm_shuffle_epi8(A0, kShuffle), kFill);  \
    const __m128i B1 = _mm_or_si128(_mm_shuffle_epi8(A1, kShuffle), kFill);  \
    _mm_storeu_si128((__m128i*)(dst + 4 * x +  0), B0);                      \
    _mm_storeu_si128((__m128i*)(dst + 4 * x + 16), B1);                      \
  }                                                                          \
  if (x < width) NAME##_C(src + 4 * x, width - x, dst + 4 * x);              \
}

void ConvertFrom_Argb_SSE(const uint8_t* src, uint32_t width, uint8_t* dst) {
  memcpy(dst, src, 4 * width);
}

CONVERT_FUNC_FILL(ConvertFrom_XRGB, 0x030201ff, 0x000000ffu)
void ConvertFrom_ARGB_SSE(const uint8_t* src, uint32_t width, uint8_t* dst) {
  ConvertFrom_Argb_SSE(src, width, dst);
  Premultiply(dst, width);
}

CONVERT_FUNC(ConvertFrom_rgbA, 0x02010003)
CONVERT_FUNC_FILL(ConvertFrom_RGBX, 0x020100ff, 0x000000ffu)
void ConvertFrom_RGBA_SSE(const uint8_t* src, uint32_t width, uint8_t* dst) {
  ConvertFrom_rgbA_SSE(src, width, dst);
  Premultiply(dst, width);
}

CONVERT_FUNC(ConvertFrom_bgrA, 0x00010203)
CONVERT_FUNC_FILL(ConvertFrom_BGRX, 0x000102ff, 0x000000ffu)
void ConvertFrom_BGRA_SSE(const uint8_t* src, uint32_t width, uint8_t* dst) {
  ConvertFrom_bgrA_SSE(src, width, dst);
  Premultiply(dst, width);
}

// ConvertTo

void Unmultiply(uint8_t* const dst, const uint8_t* dst_a, uint32_t width) {
  for (uint32_t x = 0; x < width; ++x) Unmult(dst_a[4 * x], dst + 4 * x);
}

CONVERT_FUNC_FILL(ConvertTo_XRGB, 0x030201ff, 0x000000ffu)
void ConvertTo_ARGB_SSE(const uint8_t* src, uint32_t width, uint8_t* dst) {
  ConvertFrom_Argb_SSE(src, width, dst);
  Unmultiply(dst + 1, dst + 0, width);
}

CONVERT_FUNC(ConvertTo_rgbA, 0x00030201)
CONVERT_FUNC_FILL(ConvertTo_RGBX, 0xff030201, 0xff000000u)
void ConvertTo_RGBA_SSE(const uint8_t* src, uint32_t width, uint8_t* dst) {
  ConvertTo_rgbA_SSE(src, width, dst);
  Unmultiply(dst + 0, dst + 3, width);
}

CONVERT_FUNC(ConvertTo_bgrA, 0x00010203)
CONVERT_FUNC_FILL(ConvertTo_BGRX, 0xff010203, 0xff000000u)
void ConvertTo_BGRA_SSE(const uint8_t* src, uint32_t width, uint8_t* dst) {
  ConvertTo_bgrA_SSE(src, width, dst);
  Unmultiply(dst + 0, dst + 3, width);
}

#undef CONVERT_FUNC
#undef PSHUFB_CST32

WP2_TSAN_IGNORE_FUNCTION void ArgbConverterDspInitSSE() {
  WP2ArgbConvertFrom[WP2_Argb_32] = ConvertFrom_Argb_SSE;
  WP2ArgbConvertFrom[WP2_ARGB_32] = ConvertFrom_ARGB_SSE;
  WP2ArgbConvertFrom[WP2_XRGB_32] = ConvertFrom_XRGB_SSE;

  WP2ArgbConvertFrom[WP2_rgbA_32] = ConvertFrom_rgbA_SSE;
  WP2ArgbConvertFrom[WP2_RGBA_32] = ConvertFrom_RGBA_SSE;
  WP2ArgbConvertFrom[WP2_RGBX_32] = ConvertFrom_RGBX_SSE;

  WP2ArgbConvertFrom[WP2_bgrA_32] = ConvertFrom_bgrA_SSE;
  WP2ArgbConvertFrom[WP2_BGRA_32] = ConvertFrom_BGRA_SSE;
  WP2ArgbConvertFrom[WP2_BGRX_32] = ConvertFrom_BGRX_SSE;

  WP2ArgbConvertTo[WP2_Argb_32] = ConvertFrom_Argb_SSE;
  WP2ArgbConvertTo[WP2_ARGB_32] = ConvertTo_ARGB_SSE;
  WP2ArgbConvertTo[WP2_XRGB_32] = ConvertTo_XRGB_SSE;

  WP2ArgbConvertTo[WP2_rgbA_32] = ConvertTo_rgbA_SSE;
  WP2ArgbConvertTo[WP2_RGBA_32] = ConvertTo_RGBA_SSE;
  WP2ArgbConvertTo[WP2_RGBX_32] = ConvertTo_RGBX_SSE;

  WP2ArgbConvertTo[WP2_bgrA_32] = ConvertTo_bgrA_SSE;
  WP2ArgbConvertTo[WP2_BGRA_32] = ConvertTo_BGRA_SSE;
  WP2ArgbConvertTo[WP2_BGRX_32] = ConvertTo_BGRX_SSE;
}

#endif  // WP2_USE_SSE

//------------------------------------------------------------------------------
// NEON implementation (TODO)

//------------------------------------------------------------------------------

}  // namespace

WP2ArgbConverterF WP2ArgbConvertFrom[WP2_FORMAT_NUM];
WP2ArgbConverterF WP2ArgbConvertTo[WP2_FORMAT_NUM];

static volatile WP2CPUInfo argb_converter_last_cpuinfo_used =
    (WP2CPUInfo)&argb_converter_last_cpuinfo_used;

WP2_TSAN_IGNORE_FUNCTION void WP2ArgbConverterInit() {
  if (argb_converter_last_cpuinfo_used == WP2GetCPUInfo) return;

  WP2ArgbConvertFrom[WP2_Argb_32] = ConvertFrom_Argb_C;
  WP2ArgbConvertFrom[WP2_ARGB_32] = ConvertFrom_ARGB_C;
  WP2ArgbConvertFrom[WP2_XRGB_32] = ConvertFrom_XRGB_C;

  WP2ArgbConvertFrom[WP2_rgbA_32] = ConvertFrom_rgbA_C;
  WP2ArgbConvertFrom[WP2_RGBA_32] = ConvertFrom_RGBA_C;
  WP2ArgbConvertFrom[WP2_RGBX_32] = ConvertFrom_RGBX_C;

  WP2ArgbConvertFrom[WP2_bgrA_32] = ConvertFrom_bgrA_C;
  WP2ArgbConvertFrom[WP2_BGRA_32] = ConvertFrom_BGRA_C;
  WP2ArgbConvertFrom[WP2_BGRX_32] = ConvertFrom_BGRX_C;

  WP2ArgbConvertFrom[WP2_RGB_24] = ConvertFrom_RGB_C;
  WP2ArgbConvertFrom[WP2_BGR_24] = ConvertFrom_BGR_C;

  WP2ArgbConvertTo[WP2_Argb_32] = ConvertTo_Argb_C;
  WP2ArgbConvertTo[WP2_ARGB_32] = ConvertTo_ARGB_C;
  WP2ArgbConvertTo[WP2_XRGB_32] = ConvertTo_XRGB_C;

  WP2ArgbConvertTo[WP2_rgbA_32] = ConvertTo_rgbA_C;
  WP2ArgbConvertTo[WP2_RGBA_32] = ConvertTo_RGBA_C;
  WP2ArgbConvertTo[WP2_RGBX_32] = ConvertTo_RGBX_C;

  WP2ArgbConvertTo[WP2_bgrA_32] = ConvertTo_bgrA_C;
  WP2ArgbConvertTo[WP2_BGRA_32] = ConvertTo_BGRA_C;
  WP2ArgbConvertTo[WP2_BGRX_32] = ConvertTo_BGRX_C;

  WP2ArgbConvertTo[WP2_RGB_24] = ConvertTo_RGB_C;
  WP2ArgbConvertTo[WP2_BGR_24] = ConvertTo_BGR_C;

  if (WP2GetCPUInfo != nullptr) {
#if defined(WP2_USE_SSE)
    if (WP2GetCPUInfo(kSSE)) ArgbConverterDspInitSSE();
#endif
  }
  argb_converter_last_cpuinfo_used = WP2GetCPUInfo;
}

//------------------------------------------------------------------------------
