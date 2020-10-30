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
// ANS DSP
//
// Author: Skal (pascal.massimino@gmail.com)

#include <cassert>

#include "src/dsp/dsp.h"
#include "src/dsp/math.h"

namespace WP2 {
namespace {

//------------------------------------------------------------------------------
// C version

void ANSUpdateCDF_C(uint32_t n, const uint16_t cdf_base[],
                    const uint16_t cdf_var[], uint32_t mult, uint16_t cumul[]) {
  if (mult == 65535) {  // special case 'copy'
    for (uint32_t i = 1; i < n; ++i) {
      const int32_t v = (int32_t)cdf_var[i] + cdf_base[i] - cumul[i];
      cumul[i] += v;
    }
  } else {
    for (uint32_t i = 1; i < n; ++i) {
      const int32_t v = (int32_t)cdf_var[i] + cdf_base[i] - cumul[i];
      cumul[i] += RightShift(v * mult, 16);
    }
  }
  for (uint32_t i = 1; i < n; ++i) assert(cumul[i] >= cumul[i - 1]);
}

//------------------------------------------------------------------------------
// SSE version

#if defined(WP2_USE_SSE)

void ANSUpdateCDF_SSE(uint32_t n, const uint16_t cdf_base[],
                      const uint16_t cdf_var[], uint32_t mult,
                      uint16_t cumul[]) {
  if (mult == 65535) {  // special case 'copy'
    for (uint32_t i = 0; i < n; i += 8) {
      const __m128i A = _mm_loadu_si128((const __m128i*)&cdf_var[i]);
      const __m128i B = _mm_loadu_si128((const __m128i*)&cdf_base[i]);
      const __m128i D = _mm_add_epi16(A, B);
      _mm_storeu_si128((__m128i*)&cumul[i], D);
    }
  } else {
    assert(mult < 32768);
    const __m128i M = _mm_set1_epi16(mult);
    for (uint32_t i = 0; i < n; i += 8) {
      const __m128i A = _mm_loadu_si128((const __m128i*)&cdf_var[i]);
      const __m128i B = _mm_loadu_si128((const __m128i*)&cdf_base[i]);
      const __m128i C = _mm_loadu_si128((const __m128i*)&cumul[i]);
      const __m128i D = _mm_add_epi16(A, B);
      const __m128i E = _mm_sub_epi16(D, C);
      const __m128i F = _mm_mulhi_epi16(E, M);
      const __m128i I = _mm_add_epi16(C, F);
      _mm_storeu_si128((__m128i*)&cumul[i], I);
    }
  }
}

#endif  // WP2_USE_SSE

}  // namespace

//------------------------------------------------------------------------------

ANSUpdateCDFFunc ANSUpdateCDF = nullptr;

static volatile WP2CPUInfo ans_last_cpuinfo_used =
    (WP2CPUInfo)&ans_last_cpuinfo_used;

WP2_TSAN_IGNORE_FUNCTION void ANSInit() {
  if (ans_last_cpuinfo_used == WP2GetCPUInfo) return;

  ANSUpdateCDF = ANSUpdateCDF_C;

  if (WP2GetCPUInfo != nullptr) {
#if defined(WP2_USE_SSE)
    if (WP2GetCPUInfo(kSSE)) ANSUpdateCDF = ANSUpdateCDF_SSE;
#endif
  }
  ans_last_cpuinfo_used = WP2GetCPUInfo;
}

//------------------------------------------------------------------------------

}  // namespace WP2
