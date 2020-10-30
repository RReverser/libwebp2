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
// SSIM and PSNR calculations
//
// Author: Skal (pascal.massimino@gmail.com)

#include <cassert>

#include "src/dsp/dsp.h"

#if defined(WP2_USE_SSE)
#include <stdio.h>  // for Print()
#include <initializer_list>
#endif

//------------------------------------------------------------------------------

// hat-shaped filter. Sum of coefficients is equal to 16.
static constexpr uint32_t kWeight[2 * kWP2SSIMKernel + 1] = {
  1, 2, 3, 4, 3, 2, 1
};

double WP2SSIMCalculation(uint32_t bit_depth, const WP2DistoStats& stats) {
  const uint32_t N = stats.w;
  assert(bit_depth >= 8);
  // Scaling w2 helps scaling C1, C2, C3.
  const uint32_t w2 = (N * N) << ((bit_depth - 8) * 2);
  const uint32_t C1 = 20 * w2;
  const uint32_t C2 = 60 * w2;
  const uint32_t C3 = 8 * 8 * w2;   // 'dark' limit ~= 6
  const uint64_t xmxm = (uint64_t)stats.xm * stats.xm;
  const uint64_t ymym = (uint64_t)stats.ym * stats.ym;
  if (xmxm + ymym >= C3) {
    const int64_t xmym = (int64_t)stats.xm * stats.ym;
    const int64_t sxy = stats.xym * N - xmym;    // can be negative
    const uint64_t sxx = stats.xxm * N - xmxm;
    const uint64_t syy = stats.yym * N - ymym;
    // we descale by 8 to prevent overflow during the fnum/fden multiply.
    const uint64_t num_S = (2 * (uint64_t)(sxy < 0 ? 0 : sxy) + C2) >> 8;
    const uint64_t den_S = (sxx + syy + C2) >> 8;
    const double fnum = (double)(2 * xmym + C1) * (double)num_S;
    const double fden = (double)(xmxm + ymym + C1) * (double)den_S;
    const double r = (fnum < 0. ? 0. : fnum) / fden;
    assert(r >= 0. && r <= 1.0);
    return r;
  }
  return 1.;   // area is too dark to contribute meaningfully
}

double WP2CsSSIMCalculation(uint32_t bit_depth, const WP2DistoStats& stats) {
  const uint32_t N = stats.w;
  assert(bit_depth >= 8);
  // Scaling w2 helps scaling C2, C3.
  const uint32_t w2 = (N * N) << ((bit_depth - 8) * 2);
  const uint32_t C2 = 60 * w2;
  const uint32_t C3 = 8 * 8 * w2;   // 'dark' limit ~= 6
  const uint64_t xmxm = (uint64_t)stats.xm * stats.xm;
  const uint64_t ymym = (uint64_t)stats.ym * stats.ym;
  if (xmxm + ymym >= C3) {
    const int64_t xmym = (int64_t)stats.xm * stats.ym;
    const int64_t sxy = stats.xym * N - xmym;    // can be negative
    const uint64_t sxx = stats.xxm * N - xmxm;
    const uint64_t syy = stats.yym * N - ymym;
    // we descale by 8 to prevent overflow during the fnum/fden multiply.
    const uint64_t num_S = (2 * (uint64_t)(sxy < 0 ? 0 : sxy) + C2) >> 8;
    const uint64_t den_S = (sxx + syy + C2) >> 8;
    const double r = (double)num_S / den_S;
    assert(r >= 0. && r <= 1.0);
    return r;
  }
  return 1.;   // area is too dark to contribute meaningfully
}

//------------------------------------------------------------------------------

namespace {

template <class T, uint32_t channel_step>
void SSIMGetClipped_C(const T* src1, size_t step1,
                      const T* src2, size_t step2, uint32_t xo, uint32_t yo,
                      uint32_t W, uint32_t H, WP2DistoStats* const stats) {
  const uint32_t ymin = (yo < kWP2SSIMKernel) ? 0 : yo - kWP2SSIMKernel;
  const uint32_t ymax = (yo + kWP2SSIMKernel > H - 1) ? H - 1
                                                      : yo + kWP2SSIMKernel;
  const uint32_t xmin = (xo < kWP2SSIMKernel) ? 0 : xo - kWP2SSIMKernel;
  const uint32_t xmax = (xo + kWP2SSIMKernel > W - 1) ? W - 1
                                                      : xo + kWP2SSIMKernel;
  src1 += ymin * step1;
  src2 += ymin * step2;
  for (uint32_t y = ymin; y <= ymax; ++y) {
    for (uint32_t x = xmin; x <= xmax; ++x) {
      const int32_t w =
          kWeight[kWP2SSIMKernel + x - xo] * kWeight[kWP2SSIMKernel + y - yo];
      const int32_t s1 = src1[channel_step * x];
      const int32_t s2 = src2[channel_step * x];
      stats->w += w;
      stats->xm += w * s1;
      stats->ym += w * s2;
      stats->xxm += w * s1 * s1;
      stats->xym += w * s1 * s2;
      stats->yym += w * s2 * s2;
    }
    src1 += step1;
    src2 += step2;
  }
}

template <class T, uint32_t channel_step>
void SSIMGet_C(const T* src1, size_t step1, const T* src2, size_t step2,
               WP2DistoStats* const stats) {
  for (uint32_t y = 0; y <= 2 * kWP2SSIMKernel; ++y) {
    for (uint32_t x = 0; x <= 2 * kWP2SSIMKernel; ++x) {
      const int32_t w = kWeight[x] * kWeight[y];
      const int32_t s1 = src1[channel_step * x];
      const int32_t s2 = src2[channel_step * x];
      stats->xm  += w * s1;
      stats->ym  += w * s2;
      stats->xxm += w * s1 * s1;
      stats->xym += w * s1 * s2;
      stats->yym += w * s2 * s2;
    }
    src1 += step1;
    src2 += step2;
  }
  stats->w += kWP2SSIMWeightSum;
}

//------------------------------------------------------------------------------
// SSE implementation

#if defined(WP2_USE_SSE)

// debug helper (TODO(skal): move to dsp.h ?)
template<typename T> void Print(const __m128i A, const char name[]) {
  constexpr size_t size = 16 / sizeof(T);
  T tmp[size];
  _mm_storeu_si128((__m128i*)tmp, A);
  printf("%s :", name);
  for (auto c : tmp) {
    printf((size == 16) ? " %.2x" : (size == 8) ? " %.4x" : " %.8x", c);
  }
  printf("\n");
}
void Print32b(const __m128i A, const char name[]) { Print<uint32_t>(A, name); }
void Print16b(const __m128i A, const char name[]) { Print<uint16_t>(A, name); }
void Print8b(const __m128i A, const char name[]) { Print<uint8_t>(A, name); }

//------------------------------------------------------------------------------

__m128i MakeWeight16b(int w) {
  return _mm_set_epi16(0, 1 * w, 2 * w, 3 * w, 4 * w, 3 * w, 2 * w, 1 * w);
}
const __m128i kWeight1_16b = MakeWeight16b(1);
const __m128i kWeight2_16b = MakeWeight16b(2);
const __m128i kWeight3_16b = MakeWeight16b(3);
const __m128i kWeight4_16b = MakeWeight16b(4);

__m128i MakeWeight32b(int w) {
  return _mm_set_epi32(4 * w, 3 * w, 2 * w, 1 * w);
}
const __m128i kWeight1_32b = MakeWeight32b(1);
const __m128i kWeight2_32b = MakeWeight32b(2);
const __m128i kWeight3_32b = MakeWeight32b(3);
const __m128i kWeight4_32b = MakeWeight32b(4);

uint32_t HorizontalAdd16b_SSE(const __m128i* const m) {
  uint16_t tmp[8];
  _mm_storeu_si128((__m128i*)tmp, *m);
  uint16_t res = 0;  // use a 16b accumulator to help clang find optims!
  for (uint32_t i = 0; i < 8; ++i) res += tmp[i];
  return (uint32_t)res;
}

uint32_t HorizontalAdd32b_SSE(const __m128i* const m) {
  const __m128i a = _mm_srli_si128(*m, 8);
  const __m128i b = _mm_add_epi32(*m, a);
  const __m128i c = _mm_add_epi32(b, _mm_srli_si128(b, 4));
  return (uint32_t)_mm_cvtsi128_si32(c);
}

#define ACCUMULATE_ROW_4x8u(W) do {                                 \
  /* process 7 32b-pixels at a time */                              \
  const __m128i a0 = _mm_loadu_si128((const __m128i*)(src1 +  0));  \
  const __m128i b0 = _mm_loadu_si128((const __m128i*)(src2 +  0));  \
  const __m128i c0 = _mm_loadu_si128((const __m128i*)(src1 + 16));  \
  const __m128i d0 = _mm_loadu_si128((const __m128i*)(src2 + 16));  \
  /* select the lowest 8b channel */                                \
  const __m128i a1 = _mm_and_si128(a0, mask);                       \
  const __m128i b1 = _mm_and_si128(b0, mask);                       \
  const __m128i c1 = _mm_and_si128(c0, mask);                       \
  const __m128i d1 = _mm_and_si128(d0, mask);                       \
  /* convert to 16b and multiply by weight */                       \
  const __m128i a2 = _mm_packs_epi32(a1, c1);                       \
  const __m128i b2 = _mm_packs_epi32(b1, d1);                       \
  const __m128i wa2 = _mm_mullo_epi16(a2, W);                       \
  const __m128i wb2 = _mm_mullo_epi16(b2, W);                       \
  /* accumulate */                                                  \
  xm  = _mm_add_epi16(xm, wa2);                                     \
  ym  = _mm_add_epi16(ym, wb2);                                     \
  xxm = _mm_add_epi32(xxm, _mm_madd_epi16(a2, wa2));                \
  xym = _mm_add_epi32(xym, _mm_madd_epi16(a2, wb2));                \
  yym = _mm_add_epi32(yym, _mm_madd_epi16(b2, wb2));                \
  src1 += stride1;                                                  \
  src2 += stride2;                                                  \
} while (0)

void SSIMGet_SSE_4x(const uint8_t* src1, size_t stride1,
                    const uint8_t* src2, size_t stride2,
                    WP2DistoStats* const stats) {
  const __m128i zero = _mm_setzero_si128();
  const __m128i mask = _mm_set1_epi32(0xffu);
  __m128i xm = zero, ym = zero;                 // 16b accums
  __m128i xxm = zero, yym = zero, xym = zero;   // 32b accums
  ACCUMULATE_ROW_4x8u(kWeight1_16b);
  ACCUMULATE_ROW_4x8u(kWeight2_16b);
  ACCUMULATE_ROW_4x8u(kWeight3_16b);
  ACCUMULATE_ROW_4x8u(kWeight4_16b);
  ACCUMULATE_ROW_4x8u(kWeight3_16b);
  ACCUMULATE_ROW_4x8u(kWeight2_16b);
  ACCUMULATE_ROW_4x8u(kWeight1_16b);

  stats->w   += kWP2SSIMWeightSum;
  stats->xm  += HorizontalAdd16b_SSE(&xm);
  stats->ym  += HorizontalAdd16b_SSE(&ym);
  stats->xxm += HorizontalAdd32b_SSE(&xxm);
  stats->xym += HorizontalAdd32b_SSE(&xym);
  stats->yym += HorizontalAdd32b_SSE(&yym);
}
#undef ACCUMULATE_ROW_4x8u

__m128i Load8x8bTo16s(const uint8_t* const ptr) {
  return _mm_cvtepu8_epi16(_mm_loadl_epi64((const __m128i*)ptr));
}

#define ACCUMULATE_ROW(a0, b0, W1, W2) do {                         \
  /* convert to 32b */                                              \
  const __m128i wa = _mm_madd_epi16(a0, W1);                        \
  const __m128i wb = _mm_madd_epi16(b0, W1);                        \
  const __m128i a1 = _mm_shuffle_epi8(a0, kShuffleCst);             \
  const __m128i b1 = _mm_shuffle_epi8(b0, kShuffleCst);             \
  const __m128i aa = _mm_madd_epi16(a1, a1);                        \
  const __m128i bb = _mm_madd_epi16(b1, b1);                        \
  const __m128i ab = _mm_madd_epi16(a1, b1);                        \
  /* accumulate */                                                  \
  xm  = _mm_add_epi32(xm, wa);                                      \
  ym  = _mm_add_epi32(ym, wb);                                      \
  xxm = _mm_add_epi32(xxm, _mm_mullo_epi32(aa, W2));                \
  xym = _mm_add_epi32(xym, _mm_mullo_epi32(ab, W2));                \
  yym = _mm_add_epi32(yym, _mm_mullo_epi32(bb, W2));                \
  src1 += stride1;                                                  \
  src2 += stride2;                                                  \
} while (0)

#define ACCUMULATE_ROW_8b(W1, W2) do {                              \
  /* process 7 8b-pixels at a time */                               \
  const __m128i a0 = Load8x8bTo16s(src1);                           \
  const __m128i b0 = Load8x8bTo16s(src2);                           \
  ACCUMULATE_ROW(a0, b0, W1, W2);                                   \
} while (0)

#define ACCUMULATE_ROW_16b(W1, W2) do {                             \
  /* process 7 16b-pixels at a time */                              \
  const __m128i a0 = _mm_loadu_si128((const __m128i*)(src1));       \
  const __m128i b0 = _mm_loadu_si128((const __m128i*)(src2));       \
  ACCUMULATE_ROW(a0, b0, W1, W2);                                   \
} while (0)

void SSIMGet_SSE_8u(const uint8_t* src1, size_t stride1,
                    const uint8_t* src2, size_t stride2,
                    WP2DistoStats* const stats) {
  // Constant to shuffle 16b values [A B C D E F G H] into [A G|B F|C E|D 0].
  const __m128i kShuffleCst = _mm_set_epi16(
      0xffff, 0x0706, 0x0908, 0x0504, 0x0b0a, 0x0302, 0x0d0c, 0x0100);
  const __m128i zero = _mm_setzero_si128();
  __m128i xm = zero, ym = zero;                 // 32b accums
  __m128i xxm = zero, yym = zero, xym = zero;   // 32b accums
  ACCUMULATE_ROW_8b(kWeight1_16b, kWeight1_32b);
  ACCUMULATE_ROW_8b(kWeight2_16b, kWeight2_32b);
  ACCUMULATE_ROW_8b(kWeight3_16b, kWeight3_32b);
  ACCUMULATE_ROW_8b(kWeight4_16b, kWeight4_32b);
  ACCUMULATE_ROW_8b(kWeight3_16b, kWeight3_32b);
  ACCUMULATE_ROW_8b(kWeight2_16b, kWeight2_32b);
  ACCUMULATE_ROW_8b(kWeight1_16b, kWeight1_32b);
  stats->w   += kWP2SSIMWeightSum;
  stats->xm  += HorizontalAdd32b_SSE(&xm);
  stats->ym  += HorizontalAdd32b_SSE(&ym);
  stats->xxm += HorizontalAdd32b_SSE(&xxm);
  stats->xym += HorizontalAdd32b_SSE(&xym);
  stats->yym += HorizontalAdd32b_SSE(&yym);
}

void SSIMGet_SSE_16s(const int16_t* src1, size_t stride1,
                     const int16_t* src2, size_t stride2,
                     WP2DistoStats* const stats) {
  // Constant to shuffle 16b values [A B C D E F G H] into [A G|B F|C E|D 0].
  const __m128i kShuffleCst = _mm_set_epi16(
      0xffff, 0x0706, 0x0908, 0x0504, 0x0b0a, 0x0302, 0x0d0c, 0x0100);
  const __m128i zero = _mm_setzero_si128();
  __m128i xm = zero, ym = zero;                 // 32b accums
  __m128i xxm = zero, yym = zero, xym = zero;   // 32b accums
  ACCUMULATE_ROW_16b(kWeight1_16b, kWeight1_32b);
  ACCUMULATE_ROW_16b(kWeight2_16b, kWeight2_32b);
  ACCUMULATE_ROW_16b(kWeight3_16b, kWeight3_32b);
  ACCUMULATE_ROW_16b(kWeight4_16b, kWeight4_32b);
  ACCUMULATE_ROW_16b(kWeight3_16b, kWeight3_32b);
  ACCUMULATE_ROW_16b(kWeight2_16b, kWeight2_32b);
  ACCUMULATE_ROW_16b(kWeight1_16b, kWeight1_32b);
  stats->w   += kWP2SSIMWeightSum;
  stats->xm  += HorizontalAdd32b_SSE(&xm);
  stats->ym  += HorizontalAdd32b_SSE(&ym);
  stats->xxm += HorizontalAdd32b_SSE(&xxm);
  stats->xym += HorizontalAdd32b_SSE(&xym);
  stats->yym += HorizontalAdd32b_SSE(&yym);
}
#undef ACCUMULATE_ROW_8b
#undef ACCUMULATE_ROW_16b

WP2_TSAN_IGNORE_FUNCTION void WP2SSIMInitSSE() {
  WP2SSIMGet4x8u = SSIMGet_SSE_4x;
  WP2SSIMGet8u = SSIMGet_SSE_8u;
  WP2SSIMGet16s = SSIMGet_SSE_16s;

  // to avoid compile warnings:
  (void)Print8b;
  (void)Print16b;
  (void)Print32b;
}

#endif  // WP2_USE_SSE

}  // namespace

//------------------------------------------------------------------------------

WP2SSIMGet8uFunc WP2SSIMGet4x8u = nullptr;
WP2SSIMGet8uFunc WP2SSIMGet8u = nullptr;
WP2SSIMGet16sFunc WP2SSIMGet16s = nullptr;
WP2SSIMGetClipped8uFunc WP2SSIMGetClipped4x8u = nullptr;
WP2SSIMGetClipped8uFunc WP2SSIMGetClipped8u = nullptr;
WP2SSIMGetClipped16sFunc WP2SSIMGetClipped16s = nullptr;

static volatile WP2CPUInfo ssim_last_cpuinfo_used =
    (WP2CPUInfo)&ssim_last_cpuinfo_used;

WP2_TSAN_IGNORE_FUNCTION void WP2SSIMInit() {
  if (ssim_last_cpuinfo_used == WP2GetCPUInfo) return;

  WP2SSIMGet4x8u = SSIMGet_C<uint8_t, 4>;
  WP2SSIMGet8u = SSIMGet_C<uint8_t, 1>;
  WP2SSIMGet16s = SSIMGet_C<int16_t, 1>;
  WP2SSIMGetClipped4x8u = SSIMGetClipped_C<uint8_t, 4>;
  WP2SSIMGetClipped8u = SSIMGetClipped_C<uint8_t, 1>;
  WP2SSIMGetClipped16s = SSIMGetClipped_C<int16_t, 1>;

  if (WP2GetCPUInfo != nullptr) {
#if defined(WP2_USE_SSE)
    if (WP2GetCPUInfo(kSSE)) WP2SSIMInitSSE();
#endif
  }

  ssim_last_cpuinfo_used = WP2GetCPUInfo;
}

//------------------------------------------------------------------------------
// PSNR

namespace {

template<typename T> uint64_t Accumulate(const T* const src1,
                                         const T* const src2,
                                         uint32_t len) {
  assert(sizeof(T) <= 2);
  uint64_t sum = 0;
  for (uint32_t i = 0; i < len; ++i) {
    const int32_t diff = (int32_t)src1[i] - src2[i];
    sum += (uint64_t)(diff * diff);
  }
  return sum;
}

uint64_t WP2SumSquaredError8u_C(const uint8_t* src1, const uint8_t* src2,
                                uint32_t len) {
  return Accumulate(src1, src2, len);
}

uint64_t WP2SumSquaredError16s_C(const int16_t* src1, const int16_t* src2,
                                 uint32_t len) {
  return Accumulate(src1, src2, len);
}

void WP2SumSquaredError4x8u_C(const uint8_t* src1, const uint8_t* src2,
                              uint32_t len, uint64_t result[4]) {
  for (uint32_t i = 0; i < len; ++i) {
    for (uint32_t c = 0; c < 4; ++c) {
      const int32_t diff = (int32_t)src1[4 * i + c] - src2[4 * i + c];
      result[c] += diff * diff;
    }
  }
}

void WP2SumSquaredError2x16u_C(const uint16_t* src1, const uint16_t* src2,
                               uint32_t len, uint64_t result[4]) {
  for (uint32_t i = 0; i < len; ++i) {
    for (uint32_t c = 0; c < 4; ++c) {
      const int32_t diff = (int32_t)src1[4 * i + c] - src2[4 * i + c];
      result[c] += diff * diff;
    }
  }
}

uint64_t WP2SumSquaredErrorBlock_C(const int16_t* src1, uint32_t step1,
                                   const int16_t* src2, uint32_t step2,
                                   uint32_t w, uint32_t h) {
  uint64_t sum = 0;
  for (uint32_t j = 0; j < h; ++j) {
    sum += WP2SumSquaredError16s_C(src1, src2, w);
    src1 += step1;
    src2 += step2;
  }
  return sum;
}

//------------------------------------------------------------------------------
// SSE implementation

#if defined(WP2_USE_SSE)

inline void SubtractAndAccumulate_SSE(const __m128i a, const __m128i b,
                                      __m128i* const sum) {
  const __m128i a_b = _mm_subs_epu8(a, b);
  const __m128i b_a = _mm_subs_epu8(b, a);
  const __m128i abs_a_b = _mm_or_si128(a_b, b_a);   // |a - b|
  const __m128i zero = _mm_setzero_si128();
  const __m128i C0 = _mm_unpacklo_epi8(abs_a_b, zero);
  const __m128i C1 = _mm_unpackhi_epi8(abs_a_b, zero);
  const __m128i sum1 = _mm_madd_epi16(C0, C0);  // sum{|a - b|^2}_lo
  const __m128i sum2 = _mm_madd_epi16(C1, C1);  // sum{|a - b|^2}_hi
  *sum = _mm_add_epi32(sum1, sum2);
}

uint64_t SumSquaredError8u_SSE(const uint8_t* A, const uint8_t* B,
                               uint32_t len) {
  __m128i sum = _mm_setzero_si128();
  uint32_t i = 0;
  for (; i + 16 <= len; i += 16) {
    const __m128i a = _mm_loadu_si128((const __m128i*)&A[i]);
    const __m128i b = _mm_loadu_si128((const __m128i*)&B[i]);
    __m128i tmp_sum;
    SubtractAndAccumulate_SSE(a, b, &tmp_sum);
    sum = _mm_add_epi32(sum, tmp_sum);
  }

  uint64_t sum_i = 0;
  for (; i < len; ++i) sum_i += (A[i] - B[i]) * (A[i] - B[i]);

  uint32_t tmp[4];
  _mm_storeu_si128((__m128i*)tmp, sum);
  sum_i += tmp[3] + tmp[2] + tmp[1] + tmp[0];
  return sum_i;
}

inline void SubtractAndAccumulate16s_SSE(const __m128i* const a,
                                         const __m128i* const b,
                                         __m128i* const sum1,
                                         __m128i* const sum2) {
  const __m128i A0 = _mm_sub_epi16(_mm_loadu_si128(a + 0),
                                   _mm_loadu_si128(b + 0));
  const __m128i A1 = _mm_sub_epi16(_mm_loadu_si128(a + 1),
                                   _mm_loadu_si128(b + 1));
  const __m128i B0 = _mm_madd_epi16(A0, A0);
  const __m128i B1 = _mm_madd_epi16(A1, A1);
  const __m128i C0 = _mm_add_epi32(B0, B1);
  // horizontal add 32b->64b (unsigned, since it's a squared value)
  const __m128i zero = _mm_setzero_si128();
  const __m128i D0 = _mm_unpacklo_epi32(C0, zero);
  const __m128i D1 = _mm_unpackhi_epi32(C0, zero);
  // need to accumulate in 64b, to avoid 32b overflow for long runs
  *sum1 = _mm_add_epi64(*sum1, D0);
  *sum2 = _mm_add_epi64(*sum2, D1);
}

uint64_t SumSquaredError16s_SSE(const int16_t* A, const int16_t* B,
                                uint32_t len) {
  __m128i sum1 = _mm_setzero_si128(), sum2 = sum1;
  uint32_t i = 0;
  for (; i + 16 <= len; i += 16) {
    const __m128i* const a = (const __m128i*)&A[i];
    const __m128i* const b = (const __m128i*)&B[i];
    SubtractAndAccumulate16s_SSE(a, b, &sum1, &sum2);
  }
  sum1 = _mm_add_epi64(sum1, sum2);
  uint64_t tmp[2];
  _mm_storeu_si128((__m128i*)&tmp[0], sum1);
  uint64_t sum_i = tmp[1] + tmp[0];

  for (; i < len; ++i) sum_i += (A[i] - B[i]) * (A[i] - B[i]);
  return sum_i;
}

// We must sum four lanes of four samples in parallel. So we can't
// use _mm_madd_epi16 here.
inline void SubtractAndAccumulate4x_SSE(const __m128i a, const __m128i b,
                                        __m128i* const sum) {
  const __m128i a_b = _mm_subs_epu8(a, b);
  const __m128i b_a = _mm_subs_epu8(b, a);
  const __m128i abs_a_b = _mm_or_si128(a_b, b_a);   // |a - b|
  const __m128i zero = _mm_setzero_si128();
  const __m128i C0 = _mm_unpacklo_epi8(abs_a_b, zero);
  const __m128i C1 = _mm_unpackhi_epi8(abs_a_b, zero);
  const __m128i D0 = _mm_mullo_epi16(C0, C0);
  const __m128i D1 = _mm_mulhi_epi16(C0, C0);
  const __m128i D2 = _mm_mullo_epi16(C1, C1);
  const __m128i D3 = _mm_mulhi_epi16(C1, C1);
  const __m128i E0 = _mm_unpacklo_epi16(D0, D1);
  const __m128i E1 = _mm_unpacklo_epi16(D2, D3);
  const __m128i E2 = _mm_unpackhi_epi16(D0, D1);
  const __m128i E3 = _mm_unpackhi_epi16(D2, D3);
  const __m128i F0 = _mm_add_epi32(E0, E2);
  const __m128i F1 = _mm_add_epi32(E1, E3);
  *sum = _mm_add_epi32(F0, F1);
}

void SumSquaredError4x8u_SSE(const uint8_t* A, const uint8_t* B,
                             uint32_t len, uint64_t result[4]) {
  __m128i sum = _mm_setzero_si128();
  uint32_t i = 0;
  for (; i + 4 <= len; i += 4) {
    const __m128i a = _mm_loadu_si128((const __m128i*)&A[4 * i]);
    const __m128i b = _mm_loadu_si128((const __m128i*)&B[4 * i]);
    __m128i tmp_sum;
    SubtractAndAccumulate4x_SSE(a, b, &tmp_sum);
    sum = _mm_add_epi32(sum, tmp_sum);
  }
  uint32_t tmp[4];
  _mm_storeu_si128((__m128i*)tmp, sum);
  result[0] += tmp[0];
  result[1] += tmp[1];
  result[2] += tmp[2];
  result[3] += tmp[3];
  for (; i < len; ++i) {
    for (uint32_t c = 0; c < 4; ++c) {
      const uint32_t diff = A[4 * i + c] - B[4 * i + c];
      result[c] += diff * diff;
    }
  }
}

uint64_t SumSquaredErrorBlock_SSE(const int16_t* src1, uint32_t step1,
                                  const int16_t* src2, uint32_t step2,
                                  uint32_t bw, uint32_t bh) {
  uint64_t sum_i = 0;
  uint32_t tmp[4];
  for (uint32_t j = 0; j < bh; ++j) {
    uint32_t i = 0;
    for (; i + 8 <= bw; i += 8) {
      const __m128i a0 = _mm_loadu_si128((const __m128i*)&src1[i]);
      const __m128i b0 = _mm_loadu_si128((const __m128i*)&src2[i]);
      const __m128i c0 = _mm_sub_epi16(a0, b0);
      const __m128i sum_tmp = _mm_madd_epi16(c0, c0);
      _mm_storeu_si128((__m128i*)tmp, sum_tmp);
      sum_i += tmp[3] + tmp[2] + tmp[1] + tmp[0];
    }
    for (; i + 4 <= bw; i += 4) {
      const __m128i a0 = _mm_loadl_epi64((const __m128i*)&src1[i]);
      const __m128i b0 = _mm_loadl_epi64((const __m128i*)&src2[i]);
      const __m128i c0 = _mm_sub_epi16(a0, b0);
      const __m128i sum_tmp = _mm_madd_epi16(c0, c0);
      _mm_storeu_si128((__m128i*)tmp, sum_tmp);
      sum_i += tmp[1] + tmp[0];
    }
    for (; i < bw; ++i) {
      const int32_t diff = (int32_t)src1[i] - src2[i];
      sum_i += (uint64_t)(diff * diff);
    }
    src1 += step1;
    src2 += step2;
  }
  return sum_i;
}

static WP2_TSAN_IGNORE_FUNCTION void WP2PSNRInitSSE() {
  WP2SumSquaredErrorBlock = SumSquaredErrorBlock_SSE;
  WP2SumSquaredError8u = SumSquaredError8u_SSE;
  WP2SumSquaredError4x8u = SumSquaredError4x8u_SSE;
  WP2SumSquaredError16s = SumSquaredError16s_SSE;
}

#endif  // WP2_USE_SSE

//------------------------------------------------------------------------------

}  // namespace

uint64_t (*WP2SumSquaredError8u)(const uint8_t* src1, const uint8_t* src2,
                                uint32_t len) = nullptr;
uint64_t (*WP2SumSquaredError16s)(const int16_t* src1, const int16_t* src2,
                                  uint32_t len) = nullptr;
void (*WP2SumSquaredError4x8u)(const uint8_t* src1, const uint8_t* src2,
                               uint32_t len, uint64_t result[4]) = nullptr;
void (*WP2SumSquaredError2x16u)(const uint16_t* src1, const uint16_t* src2,
                                uint32_t len, uint64_t result[4]) = nullptr;

uint64_t (*WP2SumSquaredErrorBlock)(const int16_t* src1, uint32_t step1,
                                    const int16_t* src2, uint32_t step2,
                                    uint32_t w, uint32_t h);

static volatile WP2CPUInfo psnr_last_cpuinfo_used =
    (WP2CPUInfo)&psnr_last_cpuinfo_used;

WP2_TSAN_IGNORE_FUNCTION void WP2PSNRInit() {
  if (psnr_last_cpuinfo_used == WP2GetCPUInfo) return;

  WP2SumSquaredError8u = WP2SumSquaredError8u_C;
  WP2SumSquaredError16s = WP2SumSquaredError16s_C;
  WP2SumSquaredError4x8u = WP2SumSquaredError4x8u_C;
  WP2SumSquaredError2x16u = WP2SumSquaredError2x16u_C;
  WP2SumSquaredErrorBlock = WP2SumSquaredErrorBlock_C;

  if (WP2GetCPUInfo != nullptr) {
#if defined(WP2_USE_SSE)
    if (WP2GetCPUInfo(kSSE)) WP2PSNRInitSSE();
#endif
  }

  psnr_last_cpuinfo_used = WP2GetCPUInfo;
}
