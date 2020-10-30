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
// Speed-critical deblocking functions
//
// Author: Yannis Guyon (yguyon@google.com)

#include <algorithm>
#include <numeric>

#include "src/dsp/dsp.h"
#include "src/dsp/math.h"

namespace WP2 {
namespace {

//------------------------------------------------------------------------------

constexpr uint32_t kWeightLog2 = 8;
constexpr uint32_t kStrengthShift = WP2Log2Ceil_k(kDblkMaxStrength);

static_assert(kDblkMaxHalf <= 8u, "kDblkMaxHalf should be adjusted");

// Weights representing how much to blend a pixel with the average.
// Per line length, per distance to edge. Fixed point precision is kWeightLog2.
// Values are chosen experimentally.
constexpr uint16_t kBlurWeights[kDblkMaxHalf + 1][8] = {
    {},
    {160},
    {180, 70},
    {190, 89, 22},
    {200, 113, 50, 13},
    {200, 128, 72, 32, 8},
    {200, 139, 89, 50, 22, 6},
    {200, 147, 102, 65, 37, 16, 4},
    {200, 153, 113, 78, 50, 28, 13, 3},
};

// Per sharpness. Threshold to decide whether to filter a line or not-at-all.
// Values are chosen experimentally with a fixed-point precision of 10 bits.
constexpr int32_t kSharpnessThreshold[kDblkMaxSharpness + 1] = {
  220, 212, 204, 168, 108, 72, 40, 20
};

// Per sharpness, per distance to edge. Threshold between this pixel and the
// previous one (hence index 0 is irrelevant). Values are chosen experimentally.
// TODO(yguyon): Take 'num_precision_bits' into account.
// clang-format off
constexpr int8_t kFlatnessThreshold[kDblkMaxSharpness + 1][kDblkMaxHalf] = {
  // 0x7f is the max int8_t value for thresh, to guaranty that the comparison
  // "abs_all_8b > thresh" will always be false (and result in a leading 0 bit).
  {0x7f, 82, 15, 7, 4, 3, 2, 2},
  {0x7f, 80, 10, 6, 4, 3, 2, 2},
  {0x7f, 78,  7, 5, 4, 3, 2, 2},
  {0x7f, 66,  6, 4, 4, 3, 2, 1},
  {0x7f, 54,  6, 4, 3, 2, 1, 1},
  {0x7f, 40,  6, 4, 3, 2, 1, 1},
  {0x7f, 37,  4, 3, 2, 2, 1, 1},
  {0x7f, 25,  3, 2, 2, 1, 1, 1}
};
// clang-format on

//------------------------------------------------------------------------------
// plain-C implementation

// returns ...p0|q0... edge difference
static int32_t EdgeDelta(uint32_t half, const int16_t* const q0, int32_t step) {
  assert(half > 0 && half <= kDblkMaxHalf);
  const int16_t* const p0 = q0 - step;

  // Extrapolate the edge color from the two closest pixels (consider it to be
  // between p0 and q0, at a distance of 0.5 pixels from each).
  int32_t edge_according_to_p = p0[0];
  int32_t edge_according_to_q = q0[0];
  if (half > 1) {
    edge_according_to_p += RightShift(p0[0] - p0[-step], 1);
    edge_according_to_q += RightShift(q0[0] - q0[ step], 1);
  }
  return (edge_according_to_q - edge_according_to_p);
}

bool WouldDeblockLine_C(uint32_t filter_strength, int32_t threshold,
                        uint32_t half, bool is_chroma, int32_t min, int32_t max,
                        int32_t step, int16_t* q0) {
  (void)min, (void)max, (void)is_chroma, (void)filter_strength;
  assert(filter_strength != 0);
  return (std::abs(EdgeDelta(half, q0, step)) <= threshold);
}

bool DeblockLine_C(uint32_t filter_strength, int32_t threshold, uint32_t half,
                   bool is_chroma, int32_t min, int32_t max,
                   int32_t step, int16_t* q0) {
  assert(filter_strength > 0);
  int32_t delta =  EdgeDelta(half, q0, step);
  if (std::abs(delta) > threshold) return false;
  // Do not filter invisible block edges, consider them already deblocked.
  if (std::abs(delta) <= 1) return true;
  // Color offset to go from the edge extrapolated by p0 to the adjusted edge.
  delta *= (int32_t)filter_strength;

  int16_t* p = q0 - step;
  int16_t* q = q0;
  for (uint32_t i = 0; i < half; ++i, p -= step, q += step) {
    const int32_t weighted_delta = RightShiftRound(
        delta * kBlurWeights[half][i], kWeightLog2 + kStrengthShift + 1);
    // since kBlurWeights[] is decreasing, check if we can early-out
    if (weighted_delta == 0) break;
    *p = (int16_t)Clamp(*p + weighted_delta, min, max);
    *q = (int16_t)Clamp(*q - weighted_delta, min, max);
  }
  return true;
}

uint32_t MeasureFlatLength_2(uint32_t filter_sharpness,
                             const int16_t* const q0, int32_t step) {
  const int16_t* const p = q0 - step;
  const int16_t* const q = q0;
  const int32_t diff0 = q[0] - p[0];
  const int32_t diff_q = q[step] - q[0];
  const int32_t diff_p = p[0] - p[-step];
  const int32_t t = kFlatnessThreshold[filter_sharpness][1];
  if (std::abs(diff_q - diff0) > t || std::abs(diff_p - diff0) > t) return 1;
  return 2;
}

uint32_t MeasureFlatLength_C(uint32_t filter_sharpness, uint32_t half,
                             const int16_t* q0, int32_t step) {
  assert(half > 0 && half <= kDblkMaxHalf);
  if (half == 1) return 1;
  if (half == 2) return MeasureFlatLength_2(filter_sharpness, q0, step);
  const int16_t* p = q0 - step;
  const int16_t* q = q0;
  int32_t previous_diff_p = *p - *q;
  int32_t previous_diff_q = -previous_diff_p;
  for (uint32_t i = 1; i < half; ++i) {
    const int32_t previous_p = *p, previous_q = *q;
    p -= step;
    q += step;
    const int32_t diff_p = *p - previous_p, diff_q = *q - previous_q;
    const int32_t t = kFlatnessThreshold[filter_sharpness][i];
    if (std::abs(diff_p - previous_diff_p) > t ||
        std::abs(diff_q - previous_diff_q) > t) {
      return i;
    }
    previous_diff_p = diff_p;
    previous_diff_q = diff_q;
  }
  return half;
}

}  // namespace

//------------------------------------------------------------------------------
// SSE4.1 implementation

#if defined(WP2_USE_SSE)

namespace {

inline __m128i Load_s16x4(const int16_t p[], int step) {
  return _mm_set_epi16(0, 0, 0, 0, p[3 * step], p[2 * step], p[1 * step], p[0]);
}
inline __m128i Load_s16x8(const int16_t p[], int step) {
  return _mm_set_epi16(p[7 * step], p[6 * step], p[5 * step], p[4 * step],
                       p[3 * step], p[2 * step], p[1 * step], p[0]);
}

// Loads P0/Q0 with minimal memory access.
// Q0 = [q0 q1 q2 q3 q4 q5 q6 q7]
// P0 = [p0 p1 p2 p3 p4 p5 p6 p7]
template<bool half_larger_than_4>
void GetP0Q0(const __m128i thresh, const int16_t* const q0, int32_t step,
             __m128i* const P0, __m128i* const Q0) {
  const int16_t* const p0 = q0 - step;
  if (step == 1) {
    if (half_larger_than_4) {
      *Q0 = _mm_loadu_si128((const __m128i*)(q0 + 0));
      *P0 = _mm_loadu_si128((const __m128i*)(q0 - 8));
      const __m128i kRev =
        _mm_set_epi32(0x01000302, 0x05040706, 0x09080b0a, 0x0d0c0f0e);
      *P0 = _mm_shuffle_epi8(*P0, kRev);
    } else {
      *Q0 = _mm_loadl_epi64((const __m128i*)(q0 + 0));
      *P0 = _mm_loadl_epi64((const __m128i*)(q0 - 4));
      const __m128i kRev =
        _mm_set_epi32(0, 0, 0x01000302, 0x05040706);
      *P0 = _mm_shuffle_epi8(*P0, kRev);
    }
  } else {
    if (half_larger_than_4) {
      *Q0 = Load_s16x8(q0,  step);
      *P0 = Load_s16x8(p0, -step);
    } else {
      *Q0 = Load_s16x4(q0,  step);
      *P0 = Load_s16x4(p0, -step);
    }
  }
}

uint32_t MeasureFlatLength_SSE(uint32_t filter_sharpness, uint32_t half,
                               const int16_t* q0, int32_t step) {
  assert(half > 0 && half <= kDblkMaxHalf);
  if (half == 1) return 1;
  if (half <= 2) return MeasureFlatLength_2(filter_sharpness, q0, step);
  const __m128i thresh =
      _mm_loadl_epi64((const __m128i*)kFlatnessThreshold[filter_sharpness]);
  // Q0 = [q0 q1 q2 q3 q4 q5 q6 q7]
  // P0 = [p0 p1 p2 p3 p4 p5 p6 p7]
  __m128i Q0, P0;
  if (half <= 4) {
    GetP0Q0</*half_larger_than_4=*/false>(thresh, q0, step, &P0, &Q0);
  } else {
    GetP0Q0</*half_larger_than_4=*/true>(thresh, q0, step, &P0, &Q0);
  }
  // Q1 = [p0 q0 q1 q2 q3 q4 q5 q6]
  // P1 = [q0 p0 p1 p2 p3 p4 p5 p6]
  const __m128i Q1 = _mm_insert_epi16(_mm_slli_si128(Q0, 2), q0[-step], 0);
  const __m128i P1 = _mm_insert_epi16(_mm_slli_si128(P0, 2), q0[0], 0);
  // dQ10 = [p0-q0 q0-q1 q1-q2 q2-q3 q3-q4 q4-q5 q5-q6 q6-q7]
  // dP10 = [q0-p0 p0-p1 p1-p2 p2-p3 p3-p4 p4-p5 p5-p6 p6-p7]
  const __m128i dQ10 = _mm_subs_epi16(Q1, Q0);
  const __m128i dP10 = _mm_subs_epi16(P1, P0);
  // dQ21 = [0     p0-q0 q0-q1 q1-q2 q2-q3 q3-q4 q4-q5 q5-q6]
  // dP21 = [0     q0-p0 p0-p1 p1-p2 p2-p3 p3-p4 p4-p5 p5-p6]
  const __m128i dQ21 = _mm_slli_si128(dQ10, 2);
  const __m128i dP21 = _mm_slli_si128(dP10, 2);
  // ddQ = [|q0-p0| |2q0-q1-p0| |2q1-q0-q2| |2q2-q1-q3| ... |2q6-q5-q7|]
  // ddP = [|p0-q0| |2p0-p1-q0| |2p1-p0-p2| |2p2-p1-p3| ... |2p6-p5-p7|]
  const __m128i abs_ddQ = _mm_abs_epi16(_mm_subs_epi16(dQ21, dQ10));
  const __m128i abs_ddP = _mm_abs_epi16(_mm_subs_epi16(dP21, dP10));
  const __m128i abs_all = _mm_max_epu16(abs_ddQ, abs_ddP);
  const __m128i abs_all_8b = _mm_packs_epi16(abs_all, abs_all);
  const __m128i bits_PQ = _mm_cmpgt_epi8(abs_all_8b, thresh);
  const uint32_t bits = _mm_movemask_epi8(bits_PQ);
  return WP2Ctz(bits | /*sentinel=*/(1 << half));
}

}  // namespace

static WP2_TSAN_IGNORE_FUNCTION void DblkFilterInitSSE() {
  MeasureFlatLength = MeasureFlatLength_SSE;
}

#endif  // WP2_USE_SSE
//------------------------------------------------------------------------------

int32_t DeblockThresholdFromSharpness(uint32_t filter_sharpness,
                                      uint32_t num_precision_bits) {
  return ChangePrecision(kSharpnessThreshold[filter_sharpness], 10u,
                         num_precision_bits);
}

DeblockLineF DeblockLine = nullptr;
DeblockLineF WouldDeblockLine = nullptr;
uint32_t (*MeasureFlatLength)(uint32_t filter_sharpness, uint32_t half,
                              const int16_t* const q0, int32_t step) = nullptr;

static volatile WP2CPUInfo dblk_filter_last_cpuinfo_used =
    (WP2CPUInfo)&dblk_filter_last_cpuinfo_used;

WP2_TSAN_IGNORE_FUNCTION void DblkFilterInit() {
  if (dblk_filter_last_cpuinfo_used == WP2GetCPUInfo) return;

  DeblockLine = DeblockLine_C;
  WouldDeblockLine = WouldDeblockLine_C;

  MeasureFlatLength = MeasureFlatLength_C;

  if (WP2GetCPUInfo != nullptr) {
#if defined(WP2_USE_SSE)
    if (WP2GetCPUInfo(kSSE)) DblkFilterInitSSE();
#endif
  }

  dblk_filter_last_cpuinfo_used = WP2GetCPUInfo;
}

//------------------------------------------------------------------------------

}  // namespace WP2
