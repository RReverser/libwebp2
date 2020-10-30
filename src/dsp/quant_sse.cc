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
// Quantization in SSE4.1
//
// Author: Skal (pascal.massimino@gmail.com)

#include <algorithm>
#include <cassert>

#include "src/dsp/dsp.h"

#if defined(WP2_USE_SSE)

//------------------------------------------------------------------------------

namespace {

void Quantize_SSE(const uint32_t iq[], const uint32_t bias[],
                  const int32_t res[], int16_t coeffs[], uint32_t len) {
  for (uint32_t i = 0; i < len; i += 8) {
    const __m128i IQ0 = _mm_loadu_si128((const __m128i*)&iq[i + 0]);
    const __m128i IQ1 = _mm_loadu_si128((const __m128i*)&iq[i + 4]);
    const __m128i BIAS0 = _mm_loadu_si128((const __m128i*)&bias[i + 0]);
    const __m128i BIAS1 = _mm_loadu_si128((const __m128i*)&bias[i + 4]);
    const __m128i r0 = _mm_loadu_si128((const __m128i*)&res[i + 0]);
    const __m128i r1 = _mm_loadu_si128((const __m128i*)&res[i + 4]);
    const __m128i abs_r0 = _mm_abs_epi32(r0);
    const __m128i abs_r1 = _mm_abs_epi32(r1);
    // |res| * iq + bias
    const __m128i A0 = _mm_mullo_epi32(abs_r0, IQ0);
    const __m128i A1 = _mm_mullo_epi32(abs_r1, IQ1);
    const __m128i B0 = _mm_add_epi32(A0, BIAS0);
    const __m128i B1 = _mm_add_epi32(A1, BIAS1);
    // >> WP2QBits
    const __m128i C0 = _mm_srli_epi32(B0, WP2QBits);
    const __m128i C1 = _mm_srli_epi32(B1, WP2QBits);
    // coeffs[] = sign ? -B : B;
    const __m128i D0 = _mm_sign_epi32(C0, r0);
    const __m128i D1 = _mm_sign_epi32(C1, r1);
    const __m128i E = _mm_packs_epi32(D0, D1);
    _mm_storeu_si128((__m128i*)&coeffs[i], E);
  }
}

void Dequantize_SSE(const int16_t in[], const int16_t dequants[],
                    int32_t out[], uint32_t len , uint32_t max_len) {
  uint32_t i = 0;
  for (; i < len; i += 8) {
    const __m128i A0 = _mm_loadu_si128((const __m128i*)&in[i]);
    const __m128i B0 = _mm_loadu_si128((const __m128i*)&dequants[i]);
    const __m128i C0 = _mm_mullo_epi16(A0, B0);
    const __m128i C1 = _mm_mulhi_epi16(A0, B0);
    const __m128i D0 = _mm_unpacklo_epi16(C0, C1);
    const __m128i D1 = _mm_unpackhi_epi16(C0, C1);
    _mm_storeu_si128((__m128i*)&out[i + 0], D0);
    _mm_storeu_si128((__m128i*)&out[i + 4], D1);
  }
  const __m128i zero = _mm_setzero_si128();
  for (; i < max_len; i += 4) _mm_storeu_si128((__m128i*)&out[i], zero);
}

}  // namespace

//------------------------------------------------------------------------------

extern void WP2QuantizeInitSSE();

WP2_TSAN_IGNORE_FUNCTION void WP2QuantizeInitSSE() {
  WP2Quantize = Quantize_SSE;
  WP2Dequantize = Dequantize_SSE;
}

#else  // !WP2_USE_SSE

WP2_DSP_INIT_STUB(WP2QuantizeInitSSE);

#endif  // WP2_USE_SSE
