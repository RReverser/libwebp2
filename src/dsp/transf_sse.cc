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
// SSE implementation of Transforms
//
// Author: Skal (pascal.massimino@gmail.com)
//

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdio>
#include <cstring>

#include "src/dsp/dsp.h"
#include "src/dsp/math.h"

#if defined(WP2_USE_SSE)

namespace {
//------------------------------------------------------------------------------

inline void Transpose_4x4(const __m128i* const in0, const __m128i* const in1,
                          const __m128i* const in2, const __m128i* const in3,
                          __m128i* const out0, __m128i* const out1,
                          __m128i* const out2, __m128i* const out3) {
  const __m128i t_0 = _mm_unpacklo_epi32(*in0, *in1);
  const __m128i t_1 = _mm_unpacklo_epi32(*in2, *in3);
  const __m128i t_2 = _mm_unpackhi_epi32(*in0, *in1);
  const __m128i t_3 = _mm_unpackhi_epi32(*in2, *in3);
  *out0 = _mm_unpacklo_epi64(t_0, t_1);
  *out1 = _mm_unpackhi_epi64(t_0, t_1);
  *out2 = _mm_unpacklo_epi64(t_2, t_3);
  *out3 = _mm_unpackhi_epi64(t_2, t_3);
}

void Transpose_SSE(const int32_t* in, uint32_t w, uint32_t h, int32_t* out) {
  assert(in != out);
  if (w == 2 || h == 2) {
    // fallback to C-impl for these cases
    // TODO(skal): remove need for 2x2 transforms
    for (size_t j = 0; j < h; ++j) {
      for (size_t i = 0; i < w; ++i) {
        out[j + h * i] = in[i + w * j];
      }
    }
    return;
  }
  for (size_t j = 0; j < h; j += 4) {
    const __m128i* ptr0 = (const __m128i*)&in[w * (j + 0)];
    const __m128i* ptr1 = (const __m128i*)&in[w * (j + 1)];
    const __m128i* ptr2 = (const __m128i*)&in[w * (j + 2)];
    const __m128i* ptr3 = (const __m128i*)&in[w * (j + 3)];
    for (size_t i = 0; i < w; i += 4) {
      __m128i out0, out1, out2, out3;
      const __m128i in0 = _mm_loadu_si128(ptr0++);
      const __m128i in1 = _mm_loadu_si128(ptr1++);
      const __m128i in2 = _mm_loadu_si128(ptr2++);
      const __m128i in3 = _mm_loadu_si128(ptr3++);
      Transpose_4x4(&in0, &in1, &in2, &in3, &out0, &out1, &out2, &out3);
      _mm_storeu_si128((__m128i*)&out[j + h * (i + 0)], out0);
      _mm_storeu_si128((__m128i*)&out[j + h * (i + 1)], out1);
      _mm_storeu_si128((__m128i*)&out[j + h * (i + 2)], out2);
      _mm_storeu_si128((__m128i*)&out[j + h * (i + 3)], out3);
    }
  }
}

//------------------------------------------------------------------------------

// double horizontal-add
int32_t hadd(const __m128 v0, const __m128 v1) {
  float tmp[4];
  _mm_storeu_ps(tmp, _mm_add_ps(v0, v1));
  return (int32_t)lrintf((tmp[0] + tmp[1] + tmp[2] + tmp[3]) * 256.f);
}

int32_t SlowDct8x8_SSE(const int32_t in[64],
                       const float cos_x[8], const float cos_y[8]) {
  __m128 sum0 = _mm_set1_ps(0.f);
  __m128 sum4 = _mm_set1_ps(0.f);
  const __m128 cos_x0  = _mm_loadu_ps(cos_x + 0);
  const __m128 cos_x4  = _mm_loadu_ps(cos_x + 4);
  const __m128i* src = (const __m128i*)in;
  for (uint32_t j = 0; j < 8; ++j, src += 2) {
    const __m128 C = _mm_set1_ps(cos_y[j]);
    const __m128 in0  = _mm_cvtepi32_ps(_mm_loadu_si128(src + 0));
    const __m128 in4  = _mm_cvtepi32_ps(_mm_loadu_si128(src + 1));
    const __m128 a0 = _mm_mul_ps(in0, cos_x0);
    const __m128 a4 = _mm_mul_ps(in4, cos_x4);
    sum0 = _mm_add_ps(sum0, _mm_mul_ps(a0, C));
    sum4 = _mm_add_ps(sum4, _mm_mul_ps(a4, C));
  }
  return hadd(sum0, sum4);
}

//------------------------------------------------------------------------------
}  // namespace

extern void WP2TransformInitSSE();

WP2_TSAN_IGNORE_FUNCTION void WP2TransformInitSSE() {
  WP2Transpose = Transpose_SSE;
  WP2SlowDct8x8 = SlowDct8x8_SSE;
}

#else  // !WP2_USE_SSE

WP2_DSP_INIT_STUB(WP2TransformInitSSE);

#endif  // WP2_USE_SSE
