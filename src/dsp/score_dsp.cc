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
// Scoring functions
//
// Author: Skal (pascal.massimino@gmail.com)

#include <algorithm>
#include <cmath>
#include <numeric>

#include "src/dsp/dsp.h"
#include "src/dsp/math.h"

namespace WP2 {

void GetBlockMinMaxGeneric(const int16_t* src, uint32_t step,
                           uint32_t w, uint32_t h,
                           int16_t* const min, int16_t* const max) {
  *min = src[0];
  *max = src[0];
  for (uint32_t j = 0; j < h; ++j) {
    const auto p = std::minmax_element(src, src + w);
    *min = std::min(*min, *p.first);
    *max = std::max(*max, *p.second);
    src += step;
  }
}

namespace {

//------------------------------------------------------------------------------

void GetBlockMinMax_C(const int16_t* src, uint32_t step,
                      uint32_t num_blocks,
                      int16_t min[], int16_t max[]) {
  for (uint32_t n = 0; n < num_blocks; ++n) {
    const int16_t* p = src + n * 4;
    min[n] = std::min({p[0], p[1], p[2], p[3]});
    max[n] = std::max({p[0], p[1], p[2], p[3]});
    for (uint32_t i = 1; i < 4; ++i) {
      p += step;
      min[n] = std::min({min[n], p[0], p[1], p[2], p[3]});
      max[n] = std::max({max[n], p[0], p[1], p[2], p[3]});
    }
  }
}

void GetBlockMinMax_5x5_C(const int16_t* src, uint32_t step,
                          int16_t* const min, int16_t* const max) {
  GetBlockMinMaxGeneric(src, step, 5, 5, min, max);
}

//------------------------------------------------------------------------------
// SSE4.1 implementation

#if defined(WP2_USE_SSE)

static void DoMinMax8x4(const __m128i& A, const __m128i& B, const __m128i& C,
                        const __m128i& D, __m128i* const m, __m128i* const M) {
  const __m128i m0 = _mm_min_epi16(A, B);
  const __m128i M0 = _mm_max_epi16(A, B);
  const __m128i m1 = _mm_min_epi16(C, D);
  const __m128i M1 = _mm_max_epi16(C, D);
  const __m128i m2 = _mm_min_epi16(m0, m1);
  const __m128i M2 = _mm_max_epi16(M0, M1);
  const __m128i m3 = _mm_min_epi16(m2, _mm_srli_epi64(m2, 32));
  const __m128i M3 = _mm_max_epi16(M2, _mm_srli_epi64(M2, 32));
  *m = _mm_min_epi16(m3, _mm_srli_epi64(m3, 16));
  *M = _mm_max_epi16(M3, _mm_srli_epi64(M3, 16));
}

void GetBlockMinMax_SSE(const int16_t* src, uint32_t step,
                        uint32_t num_blocks,
                        int16_t min[], int16_t max[]) {
  uint32_t n = 0;
  for (; n + 2 <= num_blocks; n += 2, src += 8) {
    const __m128i A = _mm_loadu_si128((const __m128i*)(src + 0 * step));
    const __m128i B = _mm_loadu_si128((const __m128i*)(src + 1 * step));
    const __m128i C = _mm_loadu_si128((const __m128i*)(src + 2 * step));
    const __m128i D = _mm_loadu_si128((const __m128i*)(src + 3 * step));
    __m128i m, M;
    DoMinMax8x4(A, B, C, D, &m, &M);
    min[n + 0] = _mm_extract_epi16(m, 0);
    min[n + 1] = _mm_extract_epi16(m, 4);
    max[n + 0] = _mm_extract_epi16(M, 0);
    max[n + 1] = _mm_extract_epi16(M, 4);
  }
  if (num_blocks & 1) {
    const __m128i A = _mm_loadl_epi64((const __m128i*)(src + 0 * step));
    const __m128i B = _mm_loadl_epi64((const __m128i*)(src + 1 * step));
    const __m128i C = _mm_loadl_epi64((const __m128i*)(src + 2 * step));
    const __m128i D = _mm_loadl_epi64((const __m128i*)(src + 3 * step));
    __m128i m, M;
    DoMinMax8x4(A, B, C, D, &m, &M);
    min[n + 0] = _mm_extract_epi16(m, 0);
    max[n + 0] = _mm_extract_epi16(M, 0);
  }
}

void GetBlockMinMax_5x5_SSE(const int16_t* src, uint32_t step,
                            int16_t* const min, int16_t* const max) {
  const __m128i A = _mm_loadu_si128((const __m128i*)(src + 0 * step));
  const __m128i B = _mm_loadu_si128((const __m128i*)(src + 1 * step));
  const __m128i C = _mm_loadu_si128((const __m128i*)(src + 2 * step));
  const __m128i D = _mm_loadu_si128((const __m128i*)(src + 3 * step));
  const __m128i E = _mm_loadu_si128((const __m128i*)(src + 4 * step));
  const __m128i m0 = _mm_min_epi16(A, B);
  const __m128i M0 = _mm_max_epi16(A, B);
  const __m128i m1 = _mm_min_epi16(C, D);
  const __m128i M1 = _mm_max_epi16(C, D);
  const __m128i m2 = _mm_min_epi16(m0, m1);
  const __m128i M2 = _mm_max_epi16(M0, M1);
  const __m128i m3 = _mm_min_epi16(m2, E);
  const __m128i M3 = _mm_max_epi16(M2, E);
  min[0] = _mm_extract_epi16(m3, 4);
  max[0] = _mm_extract_epi16(M3, 4);
  const __m128i m4 = _mm_min_epi16(m3, _mm_srli_epi64(m3, 32));
  const __m128i M4 = _mm_max_epi16(M3, _mm_srli_epi64(M3, 32));
  const __m128i m5 = _mm_min_epi16(m4, _mm_srli_epi64(m4, 16));
  const __m128i M5 = _mm_max_epi16(M4, _mm_srli_epi64(M4, 16));
  min[0] = std::min(min[0], (int16_t)_mm_extract_epi16(m5, 0));
  max[0] = std::max(max[0], (int16_t)_mm_extract_epi16(M5, 0));
}

WP2_TSAN_IGNORE_FUNCTION void ScoreDspInitSSE() {
  GetBlockMinMax = GetBlockMinMax_SSE;
  GetBlockMinMax_5x5 = GetBlockMinMax_5x5_SSE;
}

#endif  // WP2_USE_SSE

//------------------------------------------------------------------------------
}  // namespace

void (*GetBlockMinMax)(const int16_t* src, uint32_t step, uint32_t num_blocks,
                       int16_t min[], int16_t max[]) = nullptr;

void (*GetBlockMinMax_5x5)(const int16_t* src, uint32_t step,
                           int16_t* const min, int16_t* const max) = nullptr;

static volatile WP2CPUInfo score_last_cpuinfo_used =
    (WP2CPUInfo)&score_last_cpuinfo_used;

WP2_TSAN_IGNORE_FUNCTION void ScoreDspInit() {
  if (score_last_cpuinfo_used == WP2GetCPUInfo) return;

  GetBlockMinMax = GetBlockMinMax_C;
  GetBlockMinMax_5x5 = GetBlockMinMax_5x5_C;
  if (WP2GetCPUInfo != nullptr) {
#if defined(WP2_USE_SSE)
    if (WP2GetCPUInfo(kSSE)) ScoreDspInitSSE();
#endif
  }
  score_last_cpuinfo_used = WP2GetCPUInfo;
}

}  // namespace WP2

//------------------------------------------------------------------------------
