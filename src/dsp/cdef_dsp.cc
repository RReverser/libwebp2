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
// Speed-critical CDEF functions
// Designed to adaptively filter blocks based on identifying the direction.
//
// Implementation of AV1's Constrained Directional Enhancement Filter.
// Specification at https://aomediacodec.github.io/av1-spec/#cdef-process
// Copied from libgav1 commit bf190c43e5c7cc81751867c917a81bc2920be079 at
// https://chromium.googlesource.com/codecs/libgav1/+/refs/heads/master/src/dsp/
// The algorithm is also nicely explained on hacks.mozilla.org.
//
// Author: Yannis Guyon (yguyon@google.com)

#include <algorithm>
#include <limits>

#include "src/dsp/dsp.h"
#include "src/dsp/math.h"

namespace WP2 {
namespace {

//------------------------------------------------------------------------------

inline void Fill(int16_t src_value, int32_t width, int32_t height,
                 size_t dst_step, int16_t* dst) {
  for (int32_t y = 0; y < height; ++y, dst += dst_step) {
    std::fill(dst, dst + width, src_value);
  }
}

inline void Copy(const int16_t* src, size_t src_step,
                 int32_t width, int32_t height, size_t dst_step, int16_t* dst) {
  for (int32_t y = 0; y < height; ++y, dst += dst_step, src += src_step) {
    std::copy(src, src + width, dst);
  }
}

//------------------------------------------------------------------------------

// Given 4x4 'src' samples, returns the overall 'direction' and the 'variance'
// (certainty) of that orientation. The 'num_precision_bits' include the sign.
// Directions are expressed in [0:7], from pi/4, clockwise till -3*pi/8.
// Variance is in [0:63+], from "noise or flat area" to "most confident angle".
void CdefDirection4x4_C(const int16_t* src, int32_t step, uint32_t bitdepth,
                        uint32_t* const direction, uint32_t* const variance) {
  // TODO(yguyon): could be worth doing something more accurate at some point
  //               (but that requires way more computation, like an arctan):
  //               https://docs.opencv.org/master/d4/d70/
  //                 tutorial_anisotropic_image_segmentation_by_a_gst.html
  constexpr uint32_t kMaxNumPrecisionBits = 12;
  assert(bitdepth <= kMaxNumPrecisionBits);

  // Direction in [0:7] and y,x in [0:3] to line index in [0:6].
  // Based on https://hacks.mozilla.org/2018/06/
  //   av1-next-generation-video-the-constrained-directional-enhancement-filter/
  constexpr uint8_t kLineIndex[kDrctFltNumDirs][4][4] = {
      {{0, 1, 2, 3}, {1, 2, 3, 4}, {2, 3, 4, 5}, {3, 4, 5, 6}},
      {{0, 0, 1, 1}, {1, 1, 2, 2}, {2, 2, 3, 3}, {3, 3, 4, 4}},
      {{0, 0, 0, 0}, {1, 1, 1, 1}, {2, 2, 2, 2}, {3, 3, 3, 3}},
      {{1, 1, 0, 0}, {2, 2, 1, 1}, {3, 3, 2, 2}, {4, 4, 3, 3}},
      {{3, 2, 1, 0}, {4, 3, 2, 1}, {5, 4, 3, 2}, {6, 5, 4, 3}},
      {{3, 2, 1, 0}, {3, 2, 1, 0}, {4, 3, 2, 1}, {4, 3, 2, 1}},
      {{0, 1, 2, 3}, {0, 1, 2, 3}, {0, 1, 2, 3}, {0, 1, 2, 3}},
      {{0, 1, 2, 3}, {0, 1, 2, 3}, {1, 2, 3, 4}, {1, 2, 3, 4}}};
  constexpr uint32_t kNumLines = 4 + 4 - 1;
  constexpr uint8_t kNumPixelsPerLine[kDrctFltNumDirs][kNumLines] = {
      {1, 2, 3, 4, 3, 2, 1}, {2, 4, 4, 4, 2, 0, 0}, {4, 4, 4, 4, 0, 0, 0},
      {2, 4, 4, 4, 2, 0, 0}, {1, 2, 3, 4, 3, 2, 1}, {2, 4, 4, 4, 2, 0, 0},
      {4, 4, 4, 4, 0, 0, 0}, {2, 4, 4, 4, 2, 0, 0}};

  // For each direction, compute the average of all lines of pixels parallel to
  // that direction.
  const int16_t* row = src;
  int32_t line_avg[kDrctFltNumDirs][kNumLines] = {{0}};
  for (uint32_t y = 0; y < 4; ++y) {
    for (uint32_t x = 0; x < 4; ++x) {
      const int32_t value = ChangePrecision(row[x], bitdepth, 10);
      for (uint32_t d = 0; d < kDrctFltNumDirs; ++d) {
        const uint32_t line_index = kLineIndex[d][y][x];
        line_avg[d][line_index] += value;
      }
    }
    row += step;
  }

  for (uint32_t d = 0; d < kDrctFltNumDirs; ++d) {
    for (uint32_t l = 0; l < kNumLines; ++l) {
      const int32_t num = kNumPixelsPerLine[d][l];
      line_avg[d][l] = (num > 0) ? DivRound(line_avg[d][l], num) : 0;
    }
  }

  constexpr uint32_t kNumPx = 4 * 4;
  constexpr uint32_t kSumNumBits = kMaxNumPrecisionBits + WP2Log2Ceil_k(kNumPx);
  static_assert(kSumNumBits < sizeof(**line_avg) * 8 - 1, "Possible overflow");

  // For each pixel, compute its difference with the average of the line it is
  // in. Compute the overall variance of each direction.
  int64_t variances[kDrctFltNumDirs] = {0};
  for (uint32_t d = 0; d < kDrctFltNumDirs; ++d) {
    row = src;
    for (uint32_t y = 0; y < 4; ++y) {
      for (uint32_t x = 0; x < 4; ++x) {
        const int32_t value = ChangePrecision(row[x], bitdepth, 10);
        const uint32_t line_index = kLineIndex[d][y][x];
        const int64_t diff = value - line_avg[d][line_index];
        variances[d] += diff * diff;
      }
      row += step;
    }
  }

  constexpr uint32_t kDiffNumBits = kMaxNumPrecisionBits + kSumNumBits;
  constexpr uint32_t kVarianceNumBits =
      kDiffNumBits * 2 + WP2Log2Ceil_k(kNumPx);
  static_assert(kVarianceNumBits + 1 < sizeof(variances[0]) * 8 - 1,
                "Possible overflow");

  // The best 'direction' is the one with the lowest variance (pixels blending
  // very well along the lines).
  *direction =
      std::min_element(variances, variances + kDrctFltNumDirs) - variances;
  const uint32_t orthogonal_direction =
      (*direction + kDrctFltNumDirs / 2) % kDrctFltNumDirs;
  // The 'variance' is the gap between the variances of the best and orthogonal
  // directions.
  *variance = (uint32_t)(
      (variances[orthogonal_direction] - variances[*direction]) >> 10);
}

//------------------------------------------------------------------------------
// Constants taken from libgav1 files.

constexpr int kCdefSecondaryTap0 = 2;
constexpr int kCdefSecondaryTap1 = 1;

constexpr uint8_t kCdefPrimaryTaps[2][2] = {{4, 2}, {3, 3}};

// This is Cdef_Directions (section 7.15.3) with 2 padding entries at the
// beginning and end of the table. The cdef direction range is [0, 7] and the
// first index is offset +/-2. This removes the need to constrain the first
// index to the same range using e.g., & 7.
constexpr int8_t kCdefDirectionsPadded[12][2][2] = {
    {{1, 0}, {2, 0}},    // Padding: Cdef_Directions[6]
    {{1, 0}, {2, -1}},   // Padding: Cdef_Directions[7]
    {{-1, 1}, {-2, 2}},  // Begin Cdef_Directions
    {{0, 1}, {-1, 2}},   //
    {{0, 1}, {0, 2}},    //
    {{0, 1}, {1, 2}},    //
    {{1, 1}, {2, 2}},    //
    {{1, 0}, {2, 1}},    //
    {{1, 0}, {2, 0}},    //
    {{1, 0}, {2, -1}},   // End Cdef_Directions
    {{-1, 1}, {-2, 2}},  // Padding: Cdef_Directions[0]
    {{0, 1}, {-1, 2}},   // Padding: Cdef_Directions[1]
};

const int8_t (*const kCdefDirections)[2][2] = kCdefDirectionsPadded + 2;

//------------------------------------------------------------------------------
// CdefDirection_C() from libgav1/src/dsp/cdef.cc adapted for int16_t.

constexpr uint32_t kDivisionTable[] = {
  840, 420, 280, 210, 168, 140, 120, 105
};
constexpr uint32_t kDivisionTableOdd[] = {  // = kDivisionTable[2 * i + 1] + pad
  420, 210, 140, 105, 105, 105, 105, 105
};

constexpr int32_t Square(int32_t x) { return x * x; }

void CdefDirection8x8_C(const int16_t* src, int32_t step, uint32_t bitdepth,
                        uint32_t* const direction, uint32_t* const variance) {
  assert(direction != nullptr);
  assert(variance != nullptr);
  const uint32_t src_shift = bitdepth - 8;
  int32_t cost[8] = {};
  int32_t partial[8][15] = {};
  for (uint32_t i = 0; i < 8; ++i) {
    for (uint32_t j = 0; j < 8; ++j) {
      const int x = RightShift(src[j], src_shift);
      partial[0][i + j] += x;
      partial[1][i + j / 2] += x;
      partial[2][i] += x;
      partial[3][3 + i - j / 2] += x;
      partial[4][7 + i - j] += x;
      partial[5][3 - i / 2 + j] += x;
      partial[6][j] += x;
      partial[7][i / 2 + j] += x;
    }
    src += step;
  }
  for (uint32_t i = 0; i < 8; ++i) {
    cost[2] += Square(partial[2][i]);
    cost[6] += Square(partial[6][i]);
  }
  cost[2] *= kDivisionTable[7];
  cost[6] *= kDivisionTable[7];

  for (uint32_t i = 0; i < 7; ++i) {
    cost[0] += (Square(partial[0][i]) + Square(partial[0][14 - i])) *
               kDivisionTable[i];
    cost[4] += (Square(partial[4][i]) + Square(partial[4][14 - i])) *
               kDivisionTable[i];
  }
  cost[0] += Square(partial[0][7]) * kDivisionTable[7];
  cost[4] += Square(partial[4][7]) * kDivisionTable[7];

  for (uint32_t i = 1; i < 8; i += 2) {
    for (uint32_t j = 3; j < 8; ++j) {
      cost[i] += Square(partial[i][j]);
    }
    cost[i] *= kDivisionTable[7];
    for (uint32_t j = 0; j < 3; ++j) {
      cost[i] += (Square(partial[i][j]) + Square(partial[i][10 - j])) *
                 kDivisionTableOdd[j];
    }
  }
  int32_t best_cost = 0;
  *direction = 0;
  for (uint32_t i = 0; i < 8; ++i) {
    if (cost[i] > best_cost) {
      best_cost = cost[i];
      *direction = i;
    }
  }
  *variance = (uint32_t)(best_cost - cost[(*direction + 4) & 7]) >> 10;
}

//------------------------------------------------------------------------------

// This constant should:
//  . be larger than any sample values (10bit signed)
//  . result in a 0 value from calling Constain()
//  . and map to a large negative value when calling MapMax()
static constexpr int16_t kUnknown = 0x4000;

// Fills the extended input buffer with unfiltered pixels or unknown values.
// Area is sliced into corners, sides, and center areas as follows ->  ABBBC
// 'tmp[0]' and 'src[0]' represent the top-left E.                     DEEEF
// Margins are kDrctFltTapDist thick.                                  GHHHI
// 'top' contains [A]B[C] and 'left' contains D.
void CdefPad_C(const int16_t* const src, int32_t src_step,
               const int16_t* const left, int32_t left_step,
               const int16_t* const top, int32_t top_step,
               int32_t width, int32_t height,
               int32_t n_left, int32_t n_right,
               int32_t n_top, int32_t n_bottom,
               int16_t* const tmp, int32_t tmp_step) {
  static constexpr int32_t gap = (int32_t)kDrctFltTapDist;
  int32_t x_start = -gap, x_end = width + gap;
  int32_t y_start = -gap, y_end = height + gap;

  // Set flagged values if the margin is not available.
  if (n_top < gap) {  // ABC
    Fill(kUnknown, /*width=*/gap + width + gap, /*height=*/gap,
         tmp_step, /*dst=*/tmp - gap * tmp_step - gap);
    y_start = 0;
  }
  if (n_bottom < gap) {  // GHI
    Fill(kUnknown, /*width=*/gap + width + gap, /*height=*/gap,
         tmp_step, /*dst=*/tmp + height * tmp_step - gap);
    y_end -= gap;
  }
  if (n_left < gap) {  // D
    Fill(kUnknown, /*width=*/gap, /*height=*/y_end - y_start,
         tmp_step, /*dst=*/tmp + y_start * tmp_step - gap);
    x_start = 0;
  }
  if (n_right < gap) {  // F
    Fill(kUnknown, /*width=*/gap, /*height=*/y_end - y_start,
         tmp_step, /*dst=*/tmp + y_start * tmp_step + width);
    x_end -= gap;
  }

  // Copy pixels if they exist.
  if (n_top >= gap) {  // [A]B[C]
    Copy(top, top_step, /*width=*/x_end - x_start, /*height=*/gap,
         tmp_step, /*dst=*/tmp + y_start * tmp_step + x_start);
  }
  if (n_bottom >= gap) {  // [G]H[I]
    Copy(src + height * src_step + x_start, src_step,
         /*width=*/x_end - x_start, /*height=*/gap,
         tmp_step, /*dst=*/tmp + height * tmp_step + x_start);
  }
  if (n_left >= gap) {  // D
    Copy(left, left_step, /*width=*/gap, /*height=*/height,
         tmp_step, /*dst=*/tmp + x_start);
  }
  // E[F]
  Copy(src, src_step, /*width=*/x_end, /*height=*/height,
       tmp_step, /*dst=*/tmp);
}

//------------------------------------------------------------------------------
// CdefFilter_C() from libgav1/src/dsp/cdef.cc adapted for int16_t.

int Constrain(int diff, int threshold, int damping) {
  assert(threshold != 0);
  damping = std::max(0, damping - WP2Log2Floor(threshold));
  const int sign = (diff < 0) ? -1 : 1, abs_diff = std::abs(diff);
  return sign * Clamp(threshold - (abs_diff >> damping), 0, abs_diff);
}

// Filters the source block. It doesn't check whether the candidate pixel is
// inside the frame. However it requires the source input to be padded with a
// constant large value (kUnknown) if at the boundary.
template <bool enable_primary, bool enable_secondary>
void CdefFilter_C(const int16_t* src, int32_t src_step, uint32_t bitdepth,
                  int block_width, int block_height, int primary_strength,
                  int secondary_strength, int damping, int direction,
                  int16_t* dst, int32_t dst_step) {
  static_assert(enable_primary || enable_secondary, "Wrong template params");
  assert(block_width == 4 || block_width == 8);
  assert(block_height == 4 || block_height == 8);
  assert(direction >= 0 && direction <= 7);

  const int coeff_shift = bitdepth - 8u;
  // Section 5.9.19. CDEF params syntax.
  assert(primary_strength >= 0 && primary_strength <= 15 << coeff_shift);
  assert(secondary_strength >= 0 && secondary_strength <= 4 << coeff_shift);
  // assert(secondary_strength != 3 << coeff_shift);  // Weird so disabled.
  assert(primary_strength != 0 || secondary_strength != 0);
  // damping is decreased by 1 for chroma.
  assert((damping >= 3 && damping <= 6 + coeff_shift) ||
         (damping >= 2 && damping <= 5 + coeff_shift));
  // When only primary_strength or secondary_strength are non-zero the number of
  // pixels inspected (4 for primary_strength, 8 for secondary_strength) and the
  // taps used don't exceed the amount the sum is descaled by (16) so we can
  // skip tracking and clipping to the minimum and maximum value observed.
  constexpr bool clipping_required = enable_primary && enable_secondary;
  static constexpr int kCdefSecondaryTaps[2] = {kCdefSecondaryTap0,
                                                kCdefSecondaryTap1};
  const int pri_str_shifted_1 = (primary_strength >> coeff_shift) & 1;
  for (int y = block_height; y > 0; --y) {
    for (int x = 0; x < block_width; ++x) {
      int16_t sum = 0;
      const int16_t pixel_value = src[x];
      int16_t max_value = pixel_value;
      int16_t min_value = pixel_value;
      for (int k = 0; k <= 1; ++k) {
        for (int sign : {-1, 1}) {
          if (enable_primary) {
            const int dy = sign * kCdefDirections[direction][k][0];
            const int dx = sign * kCdefDirections[direction][k][1];
            const int16_t value = src[dy * src_step + dx + x];
            if (value != kUnknown) {
              sum += Constrain(value - pixel_value, primary_strength, damping)
                   * kCdefPrimaryTaps[pri_str_shifted_1][k];
              if (clipping_required) {
                max_value = std::max(value, max_value);
                min_value = std::min(value, min_value);
              }
            }
          }

          if (enable_secondary) {
            for (int offset : {-2, 2}) {
              const int dy = sign * kCdefDirections[direction + offset][k][0];
              const int dx = sign * kCdefDirections[direction + offset][k][1];
              const int16_t value = src[dy * src_step + dx + x];
              if (value != kUnknown) {
                sum += Constrain(value - pixel_value, secondary_strength, damping)
                     * kCdefSecondaryTaps[k];
                if (clipping_required) {
                  max_value = std::max(value, max_value);
                  min_value = std::min(value, min_value);
                }
              }
            }
          }
        }
      }

      dst[x] = pixel_value + RightShift(8 + sum - (sum < 0), 4);
      if (clipping_required) {
        dst[x] = Clamp(dst[x], min_value, max_value);
      }
    }
    src += src_step;
    dst += dst_step;
  }
}

// Dispatches to the right template instanciation based on pri/sec strengths.
void CdefFiltering_C(const int16_t* src, int32_t src_step, uint32_t bitdepth,
                     int block_width, int block_height, int primary_strength,
                     int secondary_strength, int damping, int direction,
                     int16_t* const dst, int32_t dst_step) {
  if (primary_strength > 0) {
    if (secondary_strength > 0) {
      CdefFilter_C</*enable_primary=*/true, /*enable_secondary=*/true>(
          src, src_step, bitdepth, block_width, block_height, primary_strength,
          secondary_strength, damping, direction, dst, dst_step);
    } else {
      CdefFilter_C</*enable_primary=*/true, /*enable_secondary=*/false>(
          src, src_step, bitdepth, block_width, block_height, primary_strength,
          secondary_strength, damping, direction, dst, dst_step);
    }
  } else {
    if (secondary_strength > 0) {
      CdefFilter_C</*enable_primary=*/false, /*enable_secondary=*/true>(
          src, src_step, bitdepth, block_width, block_height, primary_strength,
          secondary_strength, damping, direction, dst, dst_step);
    } else {
      // Do nothing.
    }
  }
}

//------------------------------------------------------------------------------
// SSE4.1 implementation

#if defined(WP2_USE_SSE)

template<typename T> __m128i LoadU(const T* a) {
  return _mm_loadu_si128((const __m128i*)a);
}
template<typename T> __m128i LoadUx2(const T* a, const T* b) {
  const __m128i lo = _mm_loadl_epi64((const __m128i*)a);
  const __m128i hi = _mm_loadl_epi64((const __m128i*)b);
  return _mm_unpacklo_epi64(lo, hi);
}
template<typename T> void StoreU(const __m128i v, T* dst) {
  return _mm_storeu_si128((__m128i*)dst, v);
}
template<typename T> void StoreUx2(const __m128i v, T* a, T* b) {
  _mm_storel_epi64((__m128i*)a, v);
  _mm_storel_epi64((__m128i*)b, _mm_srli_si128(v, 8));
}

// horizontal sum 32bit
inline uint32_t SumVector_S32(__m128i a) {
  a = _mm_hadd_epi32(a, a);
  a = _mm_add_epi32(a, _mm_srli_si128(a, 4));
  return _mm_cvtsi128_si32(a);
}
// Sum of squared elements.
inline uint32_t SquareSum_S16(const __m128i a) {
  const __m128i square = _mm_madd_epi16(a, a);
  return SumVector_S32(square);
}

//   A0 A1 A2 A3 E0 E1 E2 E3 |
//      B0 B1 B2 B3 F0 F1 F2 | F3
//         C0 C1 C2 C3 G0 G1 | G2 G3
//            D0 D1 D2 D3 H0 | H1 H2 H3
//  <--lo ----------------->  <-- hi --------...>
inline void StaggeredSum(const __m128i v[4], __m128i* lo, __m128i* hi) {
  const __m128i tmp0 = _mm_add_epi16(v[0],
                                     _mm_slli_si128(v[1], 2));
  const __m128i tmp1 = _mm_add_epi16(_mm_slli_si128(v[2], 4),
                                     _mm_slli_si128(v[3], 6));
  *lo = _mm_add_epi16(tmp0, tmp1);

  const __m128i tmp2 = _mm_add_epi16(_mm_srli_si128(v[1], 14),
                                     _mm_srli_si128(v[2], 12));
  *hi = _mm_add_epi16(tmp2, _mm_srli_si128(v[3], 10));
}

// ----------------------------------------------------------------------------
// 00 01 02 03 04 05 06 07                          |
//    10 11 12 13 14 15 16  17                      |-> StaggeredSum()
//       20 21 22 23 24 25  26 27                   |
//          30 31 32 33 34  35 36 37                |
//             40 41 42 43  44 45 46 47          xx
//                50 51 52  53 54 55 56 57       xx
//                   60 61  62 63 64 65 66 67    xx
//                      70  71 72 73 74 75 76 77 xx
inline void AddPartial_D0_D4(__m128i* v_src, __m128i* lo, __m128i* hi) {
  StaggeredSum(v_src, lo, hi);
  const __m128i tmp0 = _mm_add_epi16(_mm_slli_si128(v_src[4], 8),
                                     _mm_slli_si128(v_src[5], 10));
  const __m128i tmp1 = _mm_add_epi16(_mm_slli_si128(v_src[6], 12),
                                     _mm_slli_si128(v_src[7], 14));
  *lo = _mm_add_epi16(*lo, _mm_add_epi16(tmp0, tmp1));
  const __m128i tmp2 = _mm_add_epi16(_mm_srli_si128(v_src[4], 8),
                                     _mm_srli_si128(v_src[5], 6));
  const __m128i tmp3 = _mm_add_epi16(_mm_srli_si128(v_src[6], 4),
                                     _mm_srli_si128(v_src[7], 2));
  *hi = _mm_add_epi16(*hi, _mm_add_epi16(tmp2, tmp3));
}

// ----------------------------------------------------------------------------
// These two functions only differ by _mm_hadd_epi16 vs _mm_add_epi16 !

inline void AddPartial_D1_D2_D3(__m128i* v_src, __m128i* lo, __m128i* hi) {
  const __m128i kRev2 =
      _mm_set_epi32(0x09080b0a, 0x0d0c0f0e, 0x01000302, 0x05040706);
  __m128i v[4] = {
    _mm_hadd_epi16(v_src[0], v_src[4]), _mm_hadd_epi16(v_src[1], v_src[5]),
    _mm_hadd_epi16(v_src[2], v_src[6]), _mm_hadd_epi16(v_src[3], v_src[7])
  };
  StaggeredSum(v, lo + 1, hi + 1);
  for (uint32_t i = 0; i < 4; ++i) v[i] = _mm_shuffle_epi8(v[i], kRev2);
  StaggeredSum(v, lo + 3, hi + 3);
  const __m128i B0 = _mm_hadd_epi16(v[0], v[1]);
  const __m128i B1 = _mm_hadd_epi16(v[2], v[3]);
  lo[2] = _mm_hadd_epi16(B0, B1);  // hi[2] is zero
}

inline void AddPartial_D7_D6_D5(__m128i* v_src, __m128i* lo, __m128i* hi) {
  const __m128i kRev1 =
      _mm_set_epi32(0x01000302, 0x05040706, 0x09080b0a, 0x0d0c0f0e);
  __m128i v[4] = {
    _mm_add_epi16(v_src[0], v_src[1]), _mm_add_epi16(v_src[2], v_src[3]),
    _mm_add_epi16(v_src[4], v_src[5]), _mm_add_epi16(v_src[6], v_src[7])
  };
  StaggeredSum(v, lo + 7, hi + 7);
  for (uint32_t i = 0; i < 4; ++i) v[i] = _mm_shuffle_epi8(v[i], kRev1);
  StaggeredSum(v, lo + 5, hi + 5);
  const __m128i B0 = _mm_add_epi16(v[0], v[1]);
  const __m128i B1 = _mm_add_epi16(v[2], v[3]);
  lo[6] = _mm_add_epi16(B0, B1);  // hi[6] is zero
}

// ----------------------------------------------------------------------------

inline void AddPartial(const int16_t* src, uint32_t step,
                       __m128i* const partial_lo, __m128i* const partial_hi,
                       const __m128i shift) {
  __m128i v_src[8];
  for (auto& v : v_src) {
    v = _mm_sra_epi16(LoadU(src), shift);  // TODO(skal): remove the descaling?
    src += step;
  }

  AddPartial_D1_D2_D3(v_src, partial_lo, partial_hi);
  AddPartial_D7_D6_D5(v_src, partial_lo, partial_hi);
  AddPartial_D0_D4(v_src, partial_lo, partial_hi);
  const __m128i kRev1 =
      _mm_set_epi32(0x01000302, 0x05040706, 0x09080b0a, 0x0d0c0f0e);
  for (auto& v : v_src) v = _mm_shuffle_epi8(v, kRev1);
  AddPartial_D0_D4(v_src, &partial_lo[4], &partial_hi[4]);
}


inline uint32_t CostSSE(const __m128i a, const __m128i b,
                        const __m128i shuffle, const uint32_t mult[8]) {
  const __m128i b_reversed = _mm_shuffle_epi8(b, shuffle);
  const __m128i ab_lo = _mm_unpacklo_epi16(a, b_reversed);
  const __m128i ab_hi = _mm_unpackhi_epi16(a, b_reversed);
  const __m128i square_lo = _mm_madd_epi16(ab_lo, ab_lo);
  const __m128i square_hi = _mm_madd_epi16(ab_hi, ab_hi);
  const __m128i kMult0 = LoadU(mult + 0);
  const __m128i kMult4 = LoadU(mult + 4);

  const __m128i c = _mm_mullo_epi32(square_lo, kMult0);
  const __m128i d = _mm_mullo_epi32(square_hi, kMult4);
  return SumVector_S32(_mm_add_epi32(c, d));
}

void CdefDirection8x8_SSE(const int16_t* src, int32_t step, uint32_t bitdepth,
                          uint32_t* const direction, uint32_t* const variance) {
  assert(direction != nullptr);
  assert(variance != nullptr);

  const uint32_t src_shift = bitdepth - 8;
  __m128i partial_lo[8], partial_hi[8];
  const __m128i shift = _mm_cvtsi32_si128(src_shift);
  AddPartial(src, step, partial_lo, partial_hi, shift);

  uint32_t cost[8];
  cost[2] = kDivisionTable[7] * SquareSum_S16(partial_lo[2]);
  cost[6] = kDivisionTable[7] * SquareSum_S16(partial_lo[6]);

  // reverse 7 words and clear upper 2 bytes
  const __m128i kRev7 =
      _mm_set_epi32(0x80800100, 0x03020504, 0x07060908, 0x0b0a0d0c);
  cost[0] = CostSSE(partial_lo[0], partial_hi[0], kRev7, kDivisionTable);
  cost[4] = CostSSE(partial_lo[4], partial_hi[4], kRev7, kDivisionTable);

  // reverse 3 words and clear upper 10 bytes
  const __m128i kRev3 =
      _mm_set_epi32(0x80808080, 0x80808080, 0x80800100, 0x03020504);
  cost[1] = CostSSE(partial_lo[1], partial_hi[1], kRev3, kDivisionTableOdd);
  cost[3] = CostSSE(partial_lo[3], partial_hi[3], kRev3, kDivisionTableOdd);
  cost[5] = CostSSE(partial_lo[5], partial_hi[5], kRev3, kDivisionTableOdd);
  cost[7] = CostSSE(partial_lo[7], partial_hi[7], kRev3, kDivisionTableOdd);

  uint32_t best_cost = 0;
  *direction = 0;
  for (int i = 0; i < 8; ++i) {
    if (cost[i] > best_cost) {
      best_cost = cost[i];
      *direction = i;
    }
  }
  *variance = (best_cost - cost[(*direction + 4) & 7]) >> 10;
}

//------------------------------------------------------------------------------

inline void LoadDir(const int16_t src[], int stride,
                    __m128i* output, const int8_t directions[][2]) {
  const int y_0 = directions[0][0];
  const int x_0 = directions[0][1];
  const int y_1 = directions[1][0];
  const int x_1 = directions[1][1];
  output[0] = LoadU(src - y_0 * stride - x_0);
  output[1] = LoadU(src + y_0 * stride + x_0);
  output[2] = LoadU(src - y_1 * stride - x_1);
  output[3] = LoadU(src + y_1 * stride + x_1);
}

inline void LoadDirx2(const int16_t src1[], int stride,
                      __m128i* output, const int8_t directions[][2]) {
  const int16_t* const src2 = src1 + stride;
  const int y_0 = directions[0][0];
  const int x_0 = directions[0][1];
  const int y_1 = directions[1][0];
  const int x_1 = directions[1][1];
  output[0] = LoadUx2(src1 - y_0 * stride - x_0, src2 - y_0 * stride - x_0);
  output[1] = LoadUx2(src1 + y_0 * stride + x_0, src2 + y_0 * stride - x_0);
  output[2] = LoadUx2(src1 - y_1 * stride - x_1, src2 - y_1 * stride - x_1);
  output[3] = LoadUx2(src1 + y_1 * stride + x_1, src2 + y_1 * stride + x_1);
}

inline __m128i Constrain_SSE(const __m128i pixel, const __m128i ref,
                             const __m128i shift, const __m128i threshold) {
  const __m128i diff = _mm_sub_epi16(pixel, ref);
  const __m128i abs_diff = _mm_abs_epi16(diff);
  const __m128i shifted_diff = _mm_srl_epi16(abs_diff, shift);
  const __m128i max_abs_diff = _mm_subs_epu16(threshold, shifted_diff);
  const __m128i final_diff = _mm_min_epi16(max_abs_diff, abs_diff);
  return _mm_sign_epi16(final_diff, diff);
}

// map any kUnknown value to a large negative value, that will not interfere
// with the tracking of the max pixel value
inline __m128i MapMax(const __m128i max) {
  return _mm_srai_epi16(_mm_slli_epi16(max, 1), 1);
}

template <bool enable_primary, bool enable_secondary, int block_width>
void CdefFilter_SSE(const int16_t* src, int32_t src_step, uint32_t bitdepth,
                    int block_height, int primary_strength,
                    int secondary_strength, int damping, int direction,
                    int16_t* dst, int32_t dst_step) {
  assert(block_width == 4 || block_width == 8);
  assert(direction >= 0 && direction <= 7);

  __m128i primary_damping_shift, secondary_damping_shift;
  if (enable_primary) {
    primary_damping_shift =
      _mm_cvtsi32_si128(std::max(0, damping - WP2Log2Floor(primary_strength)));
  }
  if (enable_secondary) {
    secondary_damping_shift =
        _mm_cvtsi32_si128(damping - WP2Log2Floor(secondary_strength));
  }

  const int coeff_shift = bitdepth - 8u;
  const int pri_str_shifted_1 = (primary_strength >> coeff_shift) & 1;
  const __m128i primary_tap_0 =
      _mm_set1_epi16(kCdefPrimaryTaps[pri_str_shifted_1][0]);
  const __m128i primary_tap_1 =
      _mm_set1_epi16(kCdefPrimaryTaps[pri_str_shifted_1][1]);

  const __m128i primary_threshold = _mm_set1_epi16(primary_strength);
  const __m128i secondary_threshold = _mm_set1_epi16(secondary_strength);
  constexpr bool clipping_required = enable_primary && enable_secondary;

  int y = block_height;
  while (y > 0) {
    const __m128i pixel =
        (block_width == 8) ? LoadU(src) : LoadUx2(src, src + src_step);
    __m128i min = pixel;
    __m128i max = pixel;
    __m128i sum;
    __m128i vals[8];
    if (enable_primary) {
      if (block_width == 8) {
        LoadDir(src, src_step, vals, kCdefDirections[direction]);
      } else {
        LoadDirx2(src, src_step, vals, kCdefDirections[direction]);
      }
      for (uint32_t i = 0; i < 4; ++i) {
        if (clipping_required) {
          min = _mm_min_epi16(min, vals[i]);
          max = _mm_max_epi16(max, MapMax(vals[i]));
        }
        vals[i] = Constrain_SSE(vals[i], pixel,
                                primary_damping_shift, primary_threshold);
      }
      const __m128i tmp01 = _mm_add_epi16(vals[0], vals[1]);
      const __m128i tmp23 = _mm_add_epi16(vals[2], vals[3]);
      sum = _mm_add_epi16(_mm_mullo_epi16(tmp01, primary_tap_0),
                          _mm_mullo_epi16(tmp23, primary_tap_1));
    } else {
      sum = _mm_setzero_si128();
    }

    if (enable_secondary) {
      if (block_width == 8) {
        LoadDir(src, src_step, vals + 0, kCdefDirections[direction + 2]);
        LoadDir(src, src_step, vals + 4, kCdefDirections[direction - 2]);
      } else {
        LoadDirx2(src, src_step, vals + 0, kCdefDirections[direction + 2]);
        LoadDirx2(src, src_step, vals + 4, kCdefDirections[direction - 2]);
      }
      for (uint32_t i = 0; i < 8; ++i) {
        if (clipping_required) {
          min = _mm_min_epi16(min, vals[i]);
          max = _mm_max_epi16(max, MapMax(vals[i]));
        }
        vals[i] = Constrain_SSE(vals[i], pixel,
                                secondary_damping_shift, secondary_threshold);
      }

      // cascading adds
      const __m128i tmp0 = _mm_add_epi16(vals[0], vals[4]);
      const __m128i tmp1 = _mm_add_epi16(vals[1], vals[5]);
      const __m128i tmp2 = _mm_add_epi16(vals[2], vals[6]);
      const __m128i tmp3 = _mm_add_epi16(vals[3], vals[7]);
      const __m128i tmp4 = _mm_add_epi16(tmp0, tmp1);
      const __m128i tmp5 = _mm_add_epi16(tmp2, tmp3);
      const __m128i tmp6 = _mm_slli_epi16(tmp4, 1);  // <- secondary_tap0 = 2
      sum = _mm_add_epi16(sum, _mm_add_epi16(tmp5, tmp6));
    }
    const __m128i sum_lt_0 = _mm_srai_epi16(sum, 15);
    sum = _mm_add_epi16(sum, _mm_set1_epi16(8));
    sum = _mm_add_epi16(sum, sum_lt_0);
    sum = _mm_srai_epi16(sum, 4);
    sum = _mm_add_epi16(sum, pixel);
    if (clipping_required) {
      sum = _mm_min_epi16(sum, max);
      sum = _mm_max_epi16(sum, min);
    }
    if (block_width == 8) {
      StoreU(sum, dst);
      src += src_step;
      dst += dst_step;
      y -= 1;
    } else {
      StoreUx2(sum, dst, dst + dst_step);
      src += 2 * src_step;
      dst += 2 * dst_step;
      y -= 2;
    }
  }
}

// Dispatches to the right template instanciation based on pri/sec strengths.
// We use the fact that block_width is 4 or 8.
void CdefFiltering_SSE(const int16_t* src, int32_t src_step, uint32_t bitdepth,
                       int block_width, int block_height, int primary_strength,
                       int secondary_strength, int damping, int direction,
                       int16_t* const dst, int32_t dst_step) {
  if (primary_strength > 0) {
    if (secondary_strength > 0) {
      if (block_width == 4) {
        CdefFilter_SSE</*enable_primary=*/true, /*enable_secondary=*/true, 4>(
          src, src_step, bitdepth, block_height, primary_strength,
          secondary_strength, damping, direction, dst, dst_step);
      } else {
        CdefFilter_SSE</*enable_primary=*/true, /*enable_secondary=*/true, 8>(
          src, src_step, bitdepth, block_height, primary_strength,
          secondary_strength, damping, direction, dst, dst_step);
      }
    } else {
      if (block_width == 4) {
        CdefFilter_SSE</*enable_primary=*/true, /*enable_secondary=*/false, 4>(
            src, src_step, bitdepth, block_height, primary_strength,
            secondary_strength, damping, direction, dst, dst_step);
      } else {
        CdefFilter_SSE</*enable_primary=*/true, /*enable_secondary=*/false, 8>(
          src, src_step, bitdepth, block_height, primary_strength,
          secondary_strength, damping, direction, dst, dst_step);
      }
    }
  } else {
    if (secondary_strength > 0) {
      if (block_width == 4) {
        CdefFilter_SSE</*enable_primary=*/false, /*enable_secondary=*/true, 4>(
            src, src_step, bitdepth, block_height, primary_strength,
            secondary_strength, damping, direction, dst, dst_step);
      } else {
        CdefFilter_SSE</*enable_primary=*/false, /*enable_secondary=*/true, 8>(
            src, src_step, bitdepth, block_height, primary_strength,
            secondary_strength, damping, direction, dst, dst_step);
      }
    } else {
      // Do nothing.
    }
  }
}

//------------------------------------------------------------------------------

WP2_TSAN_IGNORE_FUNCTION void DrctFilterInitSSE() {
  CdefDirection8x8 = CdefDirection8x8_SSE;
  CdefFiltering = CdefFiltering_SSE;
}

#endif  // WP2_USE_SSE

}  // namespace

//------------------------------------------------------------------------------

CdefDirectionFunc CdefDirection4x4 = nullptr;
CdefDirectionFunc CdefDirection8x8 = nullptr;
CdefPadFunc CdefPad = nullptr;
CdefFilteringFunc CdefFiltering = nullptr;

static volatile WP2CPUInfo drct_filter_last_cpuinfo_used =
    (WP2CPUInfo)&drct_filter_last_cpuinfo_used;

WP2_TSAN_IGNORE_FUNCTION void DrctFilterInit() {
  if (drct_filter_last_cpuinfo_used == WP2GetCPUInfo) return;

  CdefDirection4x4 = CdefDirection4x4_C;
  CdefDirection8x8 = CdefDirection8x8_C;
  CdefPad = CdefPad_C;
  CdefFiltering = CdefFiltering_C;

  if (WP2GetCPUInfo != nullptr) {
#if defined(WP2_USE_SSE)
    if (WP2GetCPUInfo(kSSE)) DrctFilterInitSSE();
#endif
  }

  drct_filter_last_cpuinfo_used = WP2GetCPUInfo;
}

//------------------------------------------------------------------------------

}  // namespace WP2
