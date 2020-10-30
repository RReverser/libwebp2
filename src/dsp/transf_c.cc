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
// AV1 Transforms

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstring>

#include "src/dsp/dsp.h"
#include "src/dsp/math.h"

//------------------------------------------------------------------------------
// Constants for 10bit precision

// kCosPi[j] = (int)round(cos(M_PI*j/128) * (1<<kPrecision));
static constexpr uint32_t kPrecision = 10u;
static constexpr int32_t kCosPi[64] = {
  1024, 1024, 1023, 1021, 1019, 1016, 1013, 1009, 1004, 999, 993, 987, 980,
  972,  964,  955,  946,  936,  926,  915,  903,  891,  878, 865, 851, 837,
  822,  807,  792,  775,  759,  742,  724,  706,  688,  669, 650, 630, 610,
  590,  569,  548,  526,  505,  483,  460,  438,  415,  392, 369, 345, 321,
  297,  273,  249,  224,  200,  175,  150,  125,  100,  75,  50,  25
};
// Pre-calculated constants kCosPi[32] * kCosPi[16/48]
constexpr int32_t kCos32Cos48 = 277;   // cos(pi/4) * cos(3 * pi / 8)
constexpr int32_t kCos32Cos16 = 669;   // cos(pi/4) * cos(pi / 8)
// These are the same constants but with one extra precision bit (p=pre-scaled).
constexpr int32_t kCos32Cos48p = 554;    // 2 * cos(pi/4) * cos(3 * pi / 8)
constexpr int32_t kCos32Cos16p = 1337;   // 2 * cos(pi/4) * cos(pi / 8)

//  kSinPi[j] =
//    (int)round((sqrt(2) * sin(j*Pi/9) * 2 / 3) * (1 << kPrecision))
//  modified so that elements j=1,2 sum to element j=4.
static constexpr int32_t kSinPi[5] = { 0, 330, 621, 836, 951 };
static_assert(kSinPi[1] + kSinPi[2] == kSinPi[4],
              "invalid pre-condition on kSinPi[]");

//------------------------------------------------------------------------------

static inline int32_t round_shift(int64_t value, uint32_t bit = kPrecision) {
  return (int32_t)((value + ((int64_t)1 << bit >> 1)) >> bit);
}
static inline int32_t noround_shift(int64_t value, uint32_t bit = kPrecision) {
  return (value >> bit);
}

static inline int32_t half_btf(int32_t w0, int32_t in0, int32_t w1, int32_t in1,
                               uint32_t precision = kPrecision) {
  const int64_t tmp = (int64_t)w0 * in0 + (int64_t)w1 * in1;
  return round_shift(tmp, precision);
}

static inline int32_t div_sqrt2(int32_t v, uint32_t shift = 0u) {
  return round_shift((int64_t)v * kCosPi[32], kPrecision + shift);
}

//------------------------------------------------------------------------------
// Inverse transforms

static void idct2_C(const int32_t* input, int32_t* output) {
  const int32_t A = input[0], B = input[1];
  output[0] = div_sqrt2(A + B);
  output[1] = div_sqrt2(A - B);
}

static void idct4_C(const int32_t* input, int32_t* output) {
  int32_t *bf0, *bf1;
  int32_t step[4];

  // stage 1;
  assert(output != input);
  bf1 = output;
  bf1[0] = input[0];
  bf1[1] = input[2];
  bf1[2] = input[1];
  bf1[3] = input[3];

  // stage 2
  bf0 = output;
  bf1 = step;
  bf1[0] = bf0[0] + bf0[1];
  bf1[1] = bf0[0] - bf0[1];
  bf1[2] = half_btf(kCos32Cos48p, bf0[2], -kCos32Cos16p, bf0[3]);
  bf1[3] = half_btf(kCos32Cos16p, bf0[2],  kCos32Cos48p, bf0[3]);

  // stage 3
  bf0 = step;
  bf1 = output;
  bf1[0] = noround_shift(bf0[0] + bf0[3], 1);
  bf1[1] = noround_shift(bf0[1] + bf0[2], 1);
  bf1[2] = noround_shift(bf0[1] - bf0[2], 1);
  bf1[3] = noround_shift(bf0[0] - bf0[3], 1);
}

static void idct8_C(const int32_t* input, int32_t* output) {
  constexpr uint32_t shift = kPrecision + 1;
  int32_t *bf0, *bf1;
  int32_t step[8];

  // stage 1;
  assert(output != input);
  bf1 = output;
  bf1[0] = input[0];
  bf1[1] = input[4];
  bf1[2] = input[2];
  bf1[3] = input[6];
  bf1[4] = input[1];
  bf1[5] = input[5];
  bf1[6] = input[3];
  bf1[7] = input[7];

  // stage 2
  bf0 = output;
  bf1 = step;
  bf1[0] = bf0[0];
  bf1[1] = bf0[1];
  bf1[2] = bf0[2];
  bf1[3] = bf0[3];
  bf1[4] = half_btf(kCosPi[56], bf0[4], -kCosPi[8], bf0[7]);
  bf1[5] = half_btf(kCosPi[24], bf0[5], -kCosPi[40], bf0[6]);
  bf1[6] = half_btf(kCosPi[40], bf0[5], kCosPi[24], bf0[6]);
  bf1[7] = half_btf(kCosPi[8], bf0[4], kCosPi[56], bf0[7]);

  // stage 3
  bf0 = step;
  bf1 = output;
  bf1[0] = div_sqrt2(bf0[0] + bf0[1]) + 1;   // +1=rounder
  bf1[1] = div_sqrt2(bf0[0] - bf0[1]) + 1;
  bf1[2] = half_btf(kCosPi[48], bf0[2], -kCosPi[16], bf0[3]);
  bf1[3] = half_btf(kCosPi[16], bf0[2], kCosPi[48], bf0[3]);
  bf1[4] = bf0[4] + bf0[5];
  bf1[5] = bf0[4] - bf0[5];
  bf1[6] = -bf0[6] + bf0[7];
  bf1[7] = bf0[6] + bf0[7];

  // stage 4
  bf0 = output;
  bf1 = step;
  bf1[0] = bf0[0] + bf0[3];
  bf1[1] = bf0[1] + bf0[2];
  bf1[2] = bf0[1] - bf0[2];
  bf1[3] = bf0[0] - bf0[3];
  bf1[4] = bf0[4];
  bf1[5] = div_sqrt2(bf0[6] - bf0[5]);
  bf1[6] = div_sqrt2(bf0[6] + bf0[5]);
  bf1[7] = bf0[7];

  // stage 5
  bf0 = step;
  bf1 = output;
  bf1[0] = noround_shift(bf0[0] + bf0[7], shift - kPrecision);
  bf1[1] = noround_shift(bf0[1] + bf0[6], shift - kPrecision);
  bf1[2] = noround_shift(bf0[2] + bf0[5], shift - kPrecision);
  bf1[3] = noround_shift(bf0[3] + bf0[4], shift - kPrecision);
  bf1[4] = noround_shift(bf0[3] - bf0[4], shift - kPrecision);
  bf1[5] = noround_shift(bf0[2] - bf0[5], shift - kPrecision);
  bf1[6] = noround_shift(bf0[1] - bf0[6], shift - kPrecision);
  bf1[7] = noround_shift(bf0[0] - bf0[7], shift - kPrecision);
}

static void idct16_C(const int32_t* input, int32_t* output) {
  int32_t *bf0, *bf1;
  int32_t step[16];

  // stage 1;
  assert(output != input);
  bf1 = output;
  bf1[0] = input[0];
  bf1[1] = input[8];
  bf1[2] = input[4];
  bf1[3] = input[12];
  bf1[4] = input[2];
  bf1[5] = input[10];
  bf1[6] = input[6];
  bf1[7] = input[14];
  bf1[8] = input[1];
  bf1[9] = input[9];
  bf1[10] = input[5];
  bf1[11] = input[13];
  bf1[12] = input[3];
  bf1[13] = input[11];
  bf1[14] = input[7];
  bf1[15] = input[15];

  // stage 2
  bf0 = output;
  bf1 = step;
  bf1[0] = bf0[0];
  bf1[1] = bf0[1];
  bf1[2] = bf0[2];
  bf1[3] = bf0[3];
  bf1[4] = bf0[4];
  bf1[5] = bf0[5];
  bf1[6] = bf0[6];
  bf1[7] = bf0[7];
  bf1[8] = half_btf(kCosPi[60], bf0[8], -kCosPi[4], bf0[15]);
  bf1[9] = half_btf(kCosPi[28], bf0[9], -kCosPi[36], bf0[14]);
  bf1[10] = half_btf(kCosPi[44], bf0[10], -kCosPi[20], bf0[13]);
  bf1[11] = half_btf(kCosPi[12], bf0[11], -kCosPi[52], bf0[12]);
  bf1[12] = half_btf(kCosPi[52], bf0[11], kCosPi[12], bf0[12]);
  bf1[13] = half_btf(kCosPi[20], bf0[10], kCosPi[44], bf0[13]);
  bf1[14] = half_btf(kCosPi[36], bf0[9], kCosPi[28], bf0[14]);
  bf1[15] = half_btf(kCosPi[4], bf0[8], kCosPi[60], bf0[15]);

  // stage 3
  bf0 = step;
  bf1 = output;
  bf1[0] = bf0[0];
  bf1[1] = bf0[1];
  bf1[2] = bf0[2];
  bf1[3] = bf0[3];
  bf1[4] = half_btf(kCosPi[56], bf0[4], -kCosPi[8], bf0[7]);
  bf1[5] = half_btf(kCosPi[24], bf0[5], -kCosPi[40], bf0[6]);
  bf1[6] = half_btf(kCosPi[40], bf0[5], kCosPi[24], bf0[6]);
  bf1[7] = half_btf(kCosPi[8], bf0[4], kCosPi[56], bf0[7]);
  bf1[8] = bf0[8] + bf0[9];
  bf1[9] = bf0[8] - bf0[9];
  bf1[10] = -bf0[10] + bf0[11];
  bf1[11] = bf0[10] + bf0[11];
  bf1[12] = bf0[12] + bf0[13];
  bf1[13] = bf0[12] - bf0[13];
  bf1[14] = -bf0[14] + bf0[15];
  bf1[15] = bf0[14] + bf0[15];

  // stage 4
  bf0 = output;
  bf1 = step;
  bf1[0] = half_btf(kCosPi[32], bf0[0], kCosPi[32], bf0[1]);
  bf1[1] = half_btf(kCosPi[32], bf0[0], -kCosPi[32], bf0[1]);
  bf1[2] = half_btf(kCosPi[48], bf0[2], -kCosPi[16], bf0[3]);
  bf1[3] = half_btf(kCosPi[16], bf0[2], kCosPi[48], bf0[3]);
  bf1[4] = bf0[4] + bf0[5];
  bf1[5] = bf0[4] - bf0[5];
  bf1[6] = -bf0[6] + bf0[7];
  bf1[7] = bf0[6] + bf0[7];
  bf1[8] = bf0[8];
  bf1[9] = half_btf(-kCosPi[16], bf0[9], kCosPi[48], bf0[14]);
  bf1[10] = half_btf(-kCosPi[48], bf0[10], -kCosPi[16], bf0[13]);
  bf1[11] = bf0[11];
  bf1[12] = bf0[12];
  bf1[13] = half_btf(-kCosPi[16], bf0[10], kCosPi[48], bf0[13]);
  bf1[14] = half_btf(kCosPi[48], bf0[9], kCosPi[16], bf0[14]);
  bf1[15] = bf0[15];

  // stage 5
  bf0 = step;
  bf1 = output;
  bf1[0] = bf0[0] + bf0[3];
  bf1[1] = bf0[1] + bf0[2];
  bf1[2] = bf0[1] - bf0[2];
  bf1[3] = bf0[0] - bf0[3];
  bf1[4] = bf0[4];
  bf1[5] = half_btf(-kCosPi[32], bf0[5], kCosPi[32], bf0[6]);
  bf1[6] = half_btf(kCosPi[32], bf0[5], kCosPi[32], bf0[6]);
  bf1[7] = bf0[7];
  bf1[8] = bf0[8] + bf0[11];
  bf1[9] = bf0[9] + bf0[10];
  bf1[10] = bf0[9] - bf0[10];
  bf1[11] = bf0[8] - bf0[11];
  bf1[12] = -bf0[12] + bf0[15];
  bf1[13] = -bf0[13] + bf0[14];
  bf1[14] = bf0[13] + bf0[14];
  bf1[15] = bf0[12] + bf0[15];

  // stage 6
  bf0 = output;
  bf1 = step;
  bf1[0] = bf0[0] + bf0[7];
  bf1[1] = bf0[1] + bf0[6];
  bf1[2] = bf0[2] + bf0[5];
  bf1[3] = bf0[3] + bf0[4];
  bf1[4] = bf0[3] - bf0[4];
  bf1[5] = bf0[2] - bf0[5];
  bf1[6] = bf0[1] - bf0[6];
  bf1[7] = bf0[0] - bf0[7];
  bf1[8] = bf0[8];
  bf1[9] = bf0[9];
  bf1[10] = half_btf(-kCosPi[32], bf0[10], kCosPi[32], bf0[13]);
  bf1[11] = half_btf(-kCosPi[32], bf0[11], kCosPi[32], bf0[12]);
  bf1[12] = half_btf(kCosPi[32], bf0[11], kCosPi[32], bf0[12]);
  bf1[13] = half_btf(kCosPi[32], bf0[10], kCosPi[32], bf0[13]);
  bf1[14] = bf0[14];
  bf1[15] = bf0[15];

  // stage 7
  bf0 = step;
  bf1 = output;
  constexpr uint32_t shift = 1;
  bf1[0] = div_sqrt2(bf0[0] + bf0[15], shift);
  bf1[1] = div_sqrt2(bf0[1] + bf0[14], shift);
  bf1[2] = div_sqrt2(bf0[2] + bf0[13], shift);
  bf1[3] = div_sqrt2(bf0[3] + bf0[12], shift);
  bf1[4] = div_sqrt2(bf0[4] + bf0[11], shift);
  bf1[5] = div_sqrt2(bf0[5] + bf0[10], shift);
  bf1[6] = div_sqrt2(bf0[6] + bf0[9], shift);
  bf1[7] = div_sqrt2(bf0[7] + bf0[8], shift);
  bf1[8] = div_sqrt2(bf0[7] - bf0[8], shift);
  bf1[9] = div_sqrt2(bf0[6] - bf0[9], shift);
  bf1[10] = div_sqrt2(bf0[5] - bf0[10], shift);
  bf1[11] = div_sqrt2(bf0[4] - bf0[11], shift);
  bf1[12] = div_sqrt2(bf0[3] - bf0[12], shift);
  bf1[13] = div_sqrt2(bf0[2] - bf0[13], shift);
  bf1[14] = div_sqrt2(bf0[1] - bf0[14], shift);
  bf1[15] = div_sqrt2(bf0[0] - bf0[15], shift);
}

static void idct32_C(const int32_t* input, int32_t* output) {
  constexpr uint32_t shift = kPrecision + 2;
  int32_t *bf0, *bf1;
  int32_t step[32];

  // stage 1;
  assert(output != input);
  bf1 = output;
  bf1[0] = input[0];
  bf1[1] = input[16];
  bf1[2] = input[8];
  bf1[3] = input[24];
  bf1[4] = input[4];
  bf1[5] = input[20];
  bf1[6] = input[12];
  bf1[7] = input[28];
  bf1[8] = input[2];
  bf1[9] = input[18];
  bf1[10] = input[10];
  bf1[11] = input[26];
  bf1[12] = input[6];
  bf1[13] = input[22];
  bf1[14] = input[14];
  bf1[15] = input[30];
  bf1[16] = input[1];
  bf1[17] = input[17];
  bf1[18] = input[9];
  bf1[19] = input[25];
  bf1[20] = input[5];
  bf1[21] = input[21];
  bf1[22] = input[13];
  bf1[23] = input[29];
  bf1[24] = input[3];
  bf1[25] = input[19];
  bf1[26] = input[11];
  bf1[27] = input[27];
  bf1[28] = input[7];
  bf1[29] = input[23];
  bf1[30] = input[15];
  bf1[31] = input[31];

  // stage 2
  bf0 = output;
  bf1 = step;
  bf1[0] = bf0[0];
  bf1[1] = bf0[1];
  bf1[2] = bf0[2];
  bf1[3] = bf0[3];
  bf1[4] = bf0[4];
  bf1[5] = bf0[5];
  bf1[6] = bf0[6];
  bf1[7] = bf0[7];
  bf1[8] = bf0[8];
  bf1[9] = bf0[9];
  bf1[10] = bf0[10];
  bf1[11] = bf0[11];
  bf1[12] = bf0[12];
  bf1[13] = bf0[13];
  bf1[14] = bf0[14];
  bf1[15] = bf0[15];
  bf1[16] = half_btf(kCosPi[62], bf0[16], -kCosPi[2], bf0[31]);
  bf1[17] = half_btf(kCosPi[30], bf0[17], -kCosPi[34], bf0[30]);
  bf1[18] = half_btf(kCosPi[46], bf0[18], -kCosPi[18], bf0[29]);
  bf1[19] = half_btf(kCosPi[14], bf0[19], -kCosPi[50], bf0[28]);
  bf1[20] = half_btf(kCosPi[54], bf0[20], -kCosPi[10], bf0[27]);
  bf1[21] = half_btf(kCosPi[22], bf0[21], -kCosPi[42], bf0[26]);
  bf1[22] = half_btf(kCosPi[38], bf0[22], -kCosPi[26], bf0[25]);
  bf1[23] = half_btf(kCosPi[6], bf0[23], -kCosPi[58], bf0[24]);
  bf1[24] = half_btf(kCosPi[58], bf0[23], kCosPi[6], bf0[24]);
  bf1[25] = half_btf(kCosPi[26], bf0[22], kCosPi[38], bf0[25]);
  bf1[26] = half_btf(kCosPi[42], bf0[21], kCosPi[22], bf0[26]);
  bf1[27] = half_btf(kCosPi[10], bf0[20], kCosPi[54], bf0[27]);
  bf1[28] = half_btf(kCosPi[50], bf0[19], kCosPi[14], bf0[28]);
  bf1[29] = half_btf(kCosPi[18], bf0[18], kCosPi[46], bf0[29]);
  bf1[30] = half_btf(kCosPi[34], bf0[17], kCosPi[30], bf0[30]);
  bf1[31] = half_btf(kCosPi[2], bf0[16], kCosPi[62], bf0[31]);

  // stage 3
  bf0 = step;
  bf1 = output;
  bf1[0] = bf0[0];
  bf1[1] = bf0[1];
  bf1[2] = bf0[2];
  bf1[3] = bf0[3];
  bf1[4] = bf0[4];
  bf1[5] = bf0[5];
  bf1[6] = bf0[6];
  bf1[7] = bf0[7];
  bf1[8] = half_btf(kCosPi[60], bf0[8], -kCosPi[4], bf0[15]);
  bf1[9] = half_btf(kCosPi[28], bf0[9], -kCosPi[36], bf0[14]);
  bf1[10] = half_btf(kCosPi[44], bf0[10], -kCosPi[20], bf0[13]);
  bf1[11] = half_btf(kCosPi[12], bf0[11], -kCosPi[52], bf0[12]);
  bf1[12] = half_btf(kCosPi[52], bf0[11], kCosPi[12], bf0[12]);
  bf1[13] = half_btf(kCosPi[20], bf0[10], kCosPi[44], bf0[13]);
  bf1[14] = half_btf(kCosPi[36], bf0[9], kCosPi[28], bf0[14]);
  bf1[15] = half_btf(kCosPi[4], bf0[8], kCosPi[60], bf0[15]);
  bf1[16] = bf0[16] + bf0[17];
  bf1[17] = bf0[16] - bf0[17];
  bf1[18] = -bf0[18] + bf0[19];
  bf1[19] = bf0[18] + bf0[19];
  bf1[20] = bf0[20] + bf0[21];
  bf1[21] = bf0[20] - bf0[21];
  bf1[22] = -bf0[22] + bf0[23];
  bf1[23] = bf0[22] + bf0[23];
  bf1[24] = bf0[24] + bf0[25];
  bf1[25] = bf0[24] - bf0[25];
  bf1[26] = -bf0[26] + bf0[27];
  bf1[27] = bf0[26] + bf0[27];
  bf1[28] = bf0[28] + bf0[29];
  bf1[29] = bf0[28] - bf0[29];
  bf1[30] = -bf0[30] + bf0[31];
  bf1[31] = bf0[30] + bf0[31];

  // stage 4
  bf0 = output;
  bf1 = step;
  bf1[0] = bf0[0];
  bf1[1] = bf0[1];
  bf1[2] = bf0[2];
  bf1[3] = bf0[3];
  bf1[4] = half_btf(kCosPi[56], bf0[4], -kCosPi[8], bf0[7]);
  bf1[5] = half_btf(kCosPi[24], bf0[5], -kCosPi[40], bf0[6]);
  bf1[6] = half_btf(kCosPi[40], bf0[5], kCosPi[24], bf0[6]);
  bf1[7] = half_btf(kCosPi[8], bf0[4], kCosPi[56], bf0[7]);
  bf1[8] = bf0[8] + bf0[9];
  bf1[9] = bf0[8] - bf0[9];
  bf1[10] = -bf0[10] + bf0[11];
  bf1[11] = bf0[10] + bf0[11];
  bf1[12] = bf0[12] + bf0[13];
  bf1[13] = bf0[12] - bf0[13];
  bf1[14] = -bf0[14] + bf0[15];
  bf1[15] = bf0[14] + bf0[15];
  bf1[16] = bf0[16];
  bf1[17] = half_btf(-kCosPi[8], bf0[17], kCosPi[56], bf0[30]);
  bf1[18] = half_btf(-kCosPi[56], bf0[18], -kCosPi[8], bf0[29]);
  bf1[19] = bf0[19];
  bf1[20] = bf0[20];
  bf1[21] = half_btf(-kCosPi[40], bf0[21], kCosPi[24], bf0[26]);
  bf1[22] = half_btf(-kCosPi[24], bf0[22], -kCosPi[40], bf0[25]);
  bf1[23] = bf0[23];
  bf1[24] = bf0[24];
  bf1[25] = half_btf(-kCosPi[40], bf0[22], kCosPi[24], bf0[25]);
  bf1[26] = half_btf(kCosPi[24], bf0[21], kCosPi[40], bf0[26]);
  bf1[27] = bf0[27];
  bf1[28] = bf0[28];
  bf1[29] = half_btf(-kCosPi[8], bf0[18], kCosPi[56], bf0[29]);
  bf1[30] = half_btf(kCosPi[56], bf0[17], kCosPi[8], bf0[30]);
  bf1[31] = bf0[31];

  // stage 5
  bf0 = step;
  bf1 = output;
  bf1[0] = half_btf(kCosPi[32], bf0[0], kCosPi[32], bf0[1]);
  bf1[1] = half_btf(kCosPi[32], bf0[0], -kCosPi[32], bf0[1]);
  bf1[2] = half_btf(kCosPi[48], bf0[2], -kCosPi[16], bf0[3]);
  bf1[3] = half_btf(kCosPi[16], bf0[2], kCosPi[48], bf0[3]);
  bf1[4] = bf0[4] + bf0[5];
  bf1[5] = bf0[4] - bf0[5];
  bf1[6] = -bf0[6] + bf0[7];
  bf1[7] = bf0[6] + bf0[7];
  bf1[8] = bf0[8];
  bf1[9] = half_btf(-kCosPi[16], bf0[9], kCosPi[48], bf0[14]);
  bf1[10] = half_btf(-kCosPi[48], bf0[10], -kCosPi[16], bf0[13]);
  bf1[11] = bf0[11];
  bf1[12] = bf0[12];
  bf1[13] = half_btf(-kCosPi[16], bf0[10], kCosPi[48], bf0[13]);
  bf1[14] = half_btf(kCosPi[48], bf0[9], kCosPi[16], bf0[14]);
  bf1[15] = bf0[15];
  bf1[16] = bf0[16] + bf0[19];
  bf1[17] = bf0[17] + bf0[18];
  bf1[18] = bf0[17] - bf0[18];
  bf1[19] = bf0[16] - bf0[19];
  bf1[20] = -bf0[20] + bf0[23];
  bf1[21] = -bf0[21] + bf0[22];
  bf1[22] = bf0[21] + bf0[22];
  bf1[23] = bf0[20] + bf0[23];
  bf1[24] = bf0[24] + bf0[27];
  bf1[25] = bf0[25] + bf0[26];
  bf1[26] = bf0[25] - bf0[26];
  bf1[27] = bf0[24] - bf0[27];
  bf1[28] = -bf0[28] + bf0[31];
  bf1[29] = -bf0[29] + bf0[30];
  bf1[30] = bf0[29] + bf0[30];
  bf1[31] = bf0[28] + bf0[31];

  // stage 6
  bf0 = output;
  bf1 = step;
  bf1[0] = bf0[0] + bf0[3];
  bf1[1] = bf0[1] + bf0[2];
  bf1[2] = bf0[1] - bf0[2];
  bf1[3] = bf0[0] - bf0[3];
  bf1[4] = bf0[4];
  bf1[5] = half_btf(-kCosPi[32], bf0[5], kCosPi[32], bf0[6]);
  bf1[6] = half_btf(kCosPi[32], bf0[5], kCosPi[32], bf0[6]);
  bf1[7] = bf0[7];
  bf1[8] = bf0[8] + bf0[11];
  bf1[9] = bf0[9] + bf0[10];
  bf1[10] = bf0[9] - bf0[10];
  bf1[11] = bf0[8] - bf0[11];
  bf1[12] = -bf0[12] + bf0[15];
  bf1[13] = -bf0[13] + bf0[14];
  bf1[14] = bf0[13] + bf0[14];
  bf1[15] = bf0[12] + bf0[15];
  bf1[16] = bf0[16];
  bf1[17] = bf0[17];
  bf1[18] = half_btf(-kCosPi[16], bf0[18], kCosPi[48], bf0[29]);
  bf1[19] = half_btf(-kCosPi[16], bf0[19], kCosPi[48], bf0[28]);
  bf1[20] = half_btf(-kCosPi[48], bf0[20], -kCosPi[16], bf0[27]);
  bf1[21] = half_btf(-kCosPi[48], bf0[21], -kCosPi[16], bf0[26]);
  bf1[22] = bf0[22];
  bf1[23] = bf0[23];
  bf1[24] = bf0[24];
  bf1[25] = bf0[25];
  bf1[26] = half_btf(-kCosPi[16], bf0[21], kCosPi[48], bf0[26]);
  bf1[27] = half_btf(-kCosPi[16], bf0[20], kCosPi[48], bf0[27]);
  bf1[28] = half_btf(kCosPi[48], bf0[19], kCosPi[16], bf0[28]);
  bf1[29] = half_btf(kCosPi[48], bf0[18], kCosPi[16], bf0[29]);
  bf1[30] = bf0[30];
  bf1[31] = bf0[31];

  // stage 7
  bf0 = step;
  bf1 = output;
  bf1[0] = bf0[0] + bf0[7];
  bf1[1] = bf0[1] + bf0[6];
  bf1[2] = bf0[2] + bf0[5];
  bf1[3] = bf0[3] + bf0[4];
  bf1[4] = bf0[3] - bf0[4];
  bf1[5] = bf0[2] - bf0[5];
  bf1[6] = bf0[1] - bf0[6];
  bf1[7] = bf0[0] - bf0[7];
  bf1[8] = bf0[8];
  bf1[9] = bf0[9];
  bf1[10] = half_btf(-kCosPi[32], bf0[10], kCosPi[32], bf0[13]);
  bf1[11] = half_btf(-kCosPi[32], bf0[11], kCosPi[32], bf0[12]);
  bf1[12] = half_btf(kCosPi[32], bf0[11], kCosPi[32], bf0[12]);
  bf1[13] = half_btf(kCosPi[32], bf0[10], kCosPi[32], bf0[13]);
  bf1[14] = bf0[14];
  bf1[15] = bf0[15];
  bf1[16] = bf0[16] + bf0[23];
  bf1[17] = bf0[17] + bf0[22];
  bf1[18] = bf0[18] + bf0[21];
  bf1[19] = bf0[19] + bf0[20];
  bf1[20] = bf0[19] - bf0[20];
  bf1[21] = bf0[18] - bf0[21];
  bf1[22] = bf0[17] - bf0[22];
  bf1[23] = bf0[16] - bf0[23];
  bf1[24] = -bf0[24] + bf0[31];
  bf1[25] = -bf0[25] + bf0[30];
  bf1[26] = -bf0[26] + bf0[29];
  bf1[27] = -bf0[27] + bf0[28];
  bf1[28] = bf0[27] + bf0[28];
  bf1[29] = bf0[26] + bf0[29];
  bf1[30] = bf0[25] + bf0[30];
  bf1[31] = bf0[24] + bf0[31];

  // stage 8
  bf0 = output;
  bf1 = step;
  bf1[0] = bf0[0] + bf0[15];
  bf1[1] = bf0[1] + bf0[14];
  bf1[2] = bf0[2] + bf0[13];
  bf1[3] = bf0[3] + bf0[12];
  bf1[4] = bf0[4] + bf0[11];
  bf1[5] = bf0[5] + bf0[10];
  bf1[6] = bf0[6] + bf0[9];
  bf1[7] = bf0[7] + bf0[8];
  bf1[8] = bf0[7] - bf0[8];
  bf1[9] = bf0[6] - bf0[9];
  bf1[10] = bf0[5] - bf0[10];
  bf1[11] = bf0[4] - bf0[11];
  bf1[12] = bf0[3] - bf0[12];
  bf1[13] = bf0[2] - bf0[13];
  bf1[14] = bf0[1] - bf0[14];
  bf1[15] = bf0[0] - bf0[15];
  bf1[16] = bf0[16];
  bf1[17] = bf0[17];
  bf1[18] = bf0[18];
  bf1[19] = bf0[19];
  bf1[20] = half_btf(-kCosPi[32], bf0[20], kCosPi[32], bf0[27]);
  bf1[21] = half_btf(-kCosPi[32], bf0[21], kCosPi[32], bf0[26]);
  bf1[22] = half_btf(-kCosPi[32], bf0[22], kCosPi[32], bf0[25]);
  bf1[23] = half_btf(-kCosPi[32], bf0[23], kCosPi[32], bf0[24]);
  bf1[24] = half_btf(kCosPi[32], bf0[23], kCosPi[32], bf0[24]);
  bf1[25] = half_btf(kCosPi[32], bf0[22], kCosPi[32], bf0[25]);
  bf1[26] = half_btf(kCosPi[32], bf0[21], kCosPi[32], bf0[26]);
  bf1[27] = half_btf(kCosPi[32], bf0[20], kCosPi[32], bf0[27]);
  bf1[28] = bf0[28];
  bf1[29] = bf0[29];
  bf1[30] = bf0[30];
  bf1[31] = bf0[31];

  // stage 9
  bf0 = step;
  bf1 = output;
  bf1[0] = round_shift(bf0[0] + bf0[31], shift - kPrecision);
  bf1[1] = round_shift(bf0[1] + bf0[30], shift - kPrecision);
  bf1[2] = round_shift(bf0[2] + bf0[29], shift - kPrecision);
  bf1[3] = round_shift(bf0[3] + bf0[28], shift - kPrecision);
  bf1[4] = round_shift(bf0[4] + bf0[27], shift - kPrecision);
  bf1[5] = round_shift(bf0[5] + bf0[26], shift - kPrecision);
  bf1[6] = round_shift(bf0[6] + bf0[25], shift - kPrecision);
  bf1[7] = round_shift(bf0[7] + bf0[24], shift - kPrecision);
  bf1[8] = round_shift(bf0[8] + bf0[23], shift - kPrecision);
  bf1[9] = round_shift(bf0[9] + bf0[22], shift - kPrecision);
  bf1[10] = round_shift(bf0[10] + bf0[21], shift - kPrecision);
  bf1[11] = round_shift(bf0[11] + bf0[20], shift - kPrecision);
  bf1[12] = round_shift(bf0[12] + bf0[19], shift - kPrecision);
  bf1[13] = round_shift(bf0[13] + bf0[18], shift - kPrecision);
  bf1[14] = round_shift(bf0[14] + bf0[17], shift - kPrecision);
  bf1[15] = round_shift(bf0[15] + bf0[16], shift - kPrecision);
  bf1[16] = round_shift(bf0[15] - bf0[16], shift - kPrecision);
  bf1[17] = round_shift(bf0[14] - bf0[17], shift - kPrecision);
  bf1[18] = round_shift(bf0[13] - bf0[18], shift - kPrecision);
  bf1[19] = round_shift(bf0[12] - bf0[19], shift - kPrecision);
  bf1[20] = round_shift(bf0[11] - bf0[20], shift - kPrecision);
  bf1[21] = round_shift(bf0[10] - bf0[21], shift - kPrecision);
  bf1[22] = round_shift(bf0[9] - bf0[22], shift - kPrecision);
  bf1[23] = round_shift(bf0[8] - bf0[23], shift - kPrecision);
  bf1[24] = round_shift(bf0[7] - bf0[24], shift - kPrecision);
  bf1[25] = round_shift(bf0[6] - bf0[25], shift - kPrecision);
  bf1[26] = round_shift(bf0[5] - bf0[26], shift - kPrecision);
  bf1[27] = round_shift(bf0[4] - bf0[27], shift - kPrecision);
  bf1[28] = round_shift(bf0[3] - bf0[28], shift - kPrecision);
  bf1[29] = round_shift(bf0[2] - bf0[29], shift - kPrecision);
  bf1[30] = round_shift(bf0[1] - bf0[30], shift - kPrecision);
  bf1[31] = round_shift(bf0[0] - bf0[31], shift - kPrecision);
}

static void iadst4_C(const int32_t* input, int32_t* output) {
  int32_t s0, s1, s2, s3, s4, s5, s6, s7;
  // stage 0;

  int32_t x0 = input[0];
  int32_t x1 = input[1];
  int32_t x2 = input[2];
  int32_t x3 = input[3];

  if (!(x0 | x1 | x2 | x3)) {
    output[0] = output[1] = output[2] = output[3] = 0;
    return;
  }

  // stage 1
  s0 = round_shift((int64_t)kSinPi[1] * x0);
  s1 = round_shift((int64_t)kSinPi[2] * x0);
  s2 = round_shift((int64_t)kSinPi[3] * x1);
  s3 = round_shift((int64_t)kSinPi[4] * x2);
  s4 = round_shift((int64_t)kSinPi[1] * x2);
  s5 = round_shift((int64_t)kSinPi[2] * x3);
  s6 = round_shift((int64_t)kSinPi[4] * x3);
  s7 = x0 - x2;

  // stage 2
  s7 = s7 + x3;

  // stage 3
  s0 = s0 + s3;
  s1 = s1 - s4;
  s3 = s2;
  s2 = round_shift((int64_t)kSinPi[3] * s7);

  // stage 4
  s0 = s0 + s5;
  s1 = s1 - s6;

  // stage 5
  x0 = s0 + s3;
  x1 = s1 + s3;
  x2 = s2;
  x3 = s0 + s1;

  // stage 6
  x3 = x3 - s3;

  output[0] = div_sqrt2(x0);
  output[1] = div_sqrt2(x1);
  output[2] = div_sqrt2(x2);
  output[3] = div_sqrt2(x3);
}

static void iadst8_C(const int32_t* input, int32_t* output) {
  constexpr uint32_t shift = kPrecision + 1;
  int32_t *bf0, *bf1;
  int32_t step[8];

  // stage 1;
  bf1 = output;
  bf1[0] = input[7];
  bf1[1] = input[0];
  bf1[2] = input[5];
  bf1[3] = input[2];
  bf1[4] = input[3];
  bf1[5] = input[4];
  bf1[6] = input[1];
  bf1[7] = input[6];

  // stage 2
  bf0 = output;
  bf1 = step;
  bf1[0] = half_btf(kCosPi[4], bf0[0], kCosPi[60], bf0[1]);
  bf1[1] = half_btf(-kCosPi[4], bf0[1], kCosPi[60], bf0[0]);
  bf1[2] = half_btf(kCosPi[20], bf0[2], kCosPi[44], bf0[3]);
  bf1[3] = half_btf(-kCosPi[20], bf0[3], kCosPi[44], bf0[2]);
  bf1[4] = half_btf(kCosPi[36], bf0[4], kCosPi[28], bf0[5]);
  bf1[5] = half_btf(-kCosPi[36], bf0[5], kCosPi[28], bf0[4]);
  bf1[6] = half_btf(kCosPi[52], bf0[6], kCosPi[12], bf0[7]);
  bf1[7] = half_btf(-kCosPi[52], bf0[7], kCosPi[12], bf0[6]);

  // stage 3
  bf0 = step;
  bf1 = output;
  bf1[0] = bf0[0] + bf0[4];
  bf1[1] = bf0[1] + bf0[5];
  bf1[2] = bf0[2] + bf0[6];
  bf1[3] = bf0[3] + bf0[7];
  bf1[4] = -bf0[4] + bf0[0];
  bf1[5] = -bf0[5] + bf0[1];
  bf1[6] = -bf0[6] + bf0[2];
  bf1[7] = -bf0[7] + bf0[3];

  // stage 4
  bf0 = output;
  bf1 = step;
  bf1[0] = bf0[0];
  bf1[1] = bf0[1];
  bf1[2] = bf0[2];
  bf1[3] = bf0[3];
  bf1[4] = half_btf(kCosPi[16], bf0[4], kCosPi[48], bf0[5]);
  bf1[5] = half_btf(-kCosPi[16], bf0[5], kCosPi[48], bf0[4]);
  bf1[6] = half_btf(-kCosPi[48], bf0[6], kCosPi[16], bf0[7]);
  bf1[7] = half_btf(kCosPi[48], bf0[7], kCosPi[16], bf0[6]);

  // stage 5
  bf0 = step;
  bf1 = output;
  bf1[0] = bf0[0] + bf0[2];
  bf1[1] = bf0[1] + bf0[3];
  bf1[2] = -bf0[2] + bf0[0];
  bf1[3] = -bf0[3] + bf0[1];
  bf1[4] = bf0[4] + bf0[6];
  bf1[5] = bf0[5] + bf0[7];
  bf1[6] = -bf0[6] + bf0[4];
  bf1[7] = -bf0[7] + bf0[5];

  // stage 6
  bf0 = output;
  bf1 = step;
  bf1[0] = round_shift(bf0[0], shift - kPrecision);
  bf1[1] = round_shift(bf0[1], shift - kPrecision);
  bf1[2] = half_btf(kCosPi[32], bf0[2], kCosPi[32], bf0[3], shift);
  bf1[3] = half_btf(-kCosPi[32], bf0[3], kCosPi[32], bf0[2], shift);
  bf1[4] = round_shift(bf0[4], shift - kPrecision);
  bf1[5] = round_shift(bf0[5], shift - kPrecision);
  bf1[6] = half_btf(kCosPi[32], bf0[6], kCosPi[32], bf0[7], shift);
  bf1[7] = half_btf(-kCosPi[32], bf0[7], kCosPi[32], bf0[6], shift);

  // stage 7
  bf0 = step;
  bf1 = output;
  bf1[0] = bf0[0];
  bf1[1] = -bf0[4];
  bf1[2] = bf0[6];
  bf1[3] = -bf0[2];
  bf1[4] = bf0[3];
  bf1[5] = -bf0[7];
  bf1[6] = bf0[5];
  bf1[7] = -bf0[1];
}

static void iadst16_C(const int32_t* input, int32_t* output) {
  constexpr uint32_t shift = kPrecision + 1;
  int32_t *bf0, *bf1;
  int32_t step[16];

  // stage 1;
  bf1 = output;
  bf1[0] = input[15];
  bf1[1] = input[0];
  bf1[2] = input[13];
  bf1[3] = input[2];
  bf1[4] = input[11];
  bf1[5] = input[4];
  bf1[6] = input[9];
  bf1[7] = input[6];
  bf1[8] = input[7];
  bf1[9] = input[8];
  bf1[10] = input[5];
  bf1[11] = input[10];
  bf1[12] = input[3];
  bf1[13] = input[12];
  bf1[14] = input[1];
  bf1[15] = input[14];

  // stage 2
  bf0 = output;
  bf1 = step;
  bf1[0] = half_btf(kCosPi[2], bf0[0], kCosPi[62], bf0[1]);
  bf1[1] = half_btf(-kCosPi[2], bf0[1], kCosPi[62], bf0[0]);
  bf1[2] = half_btf(kCosPi[10], bf0[2], kCosPi[54], bf0[3]);
  bf1[3] = half_btf(-kCosPi[10], bf0[3], kCosPi[54], bf0[2]);
  bf1[4] = half_btf(kCosPi[18], bf0[4], kCosPi[46], bf0[5]);
  bf1[5] = half_btf(-kCosPi[18], bf0[5], kCosPi[46], bf0[4]);
  bf1[6] = half_btf(kCosPi[26], bf0[6], kCosPi[38], bf0[7]);
  bf1[7] = half_btf(-kCosPi[26], bf0[7], kCosPi[38], bf0[6]);
  bf1[8] = half_btf(kCosPi[34], bf0[8], kCosPi[30], bf0[9]);
  bf1[9] = half_btf(-kCosPi[34], bf0[9], kCosPi[30], bf0[8]);
  bf1[10] = half_btf(kCosPi[42], bf0[10], kCosPi[22], bf0[11]);
  bf1[11] = half_btf(-kCosPi[42], bf0[11], kCosPi[22], bf0[10]);
  bf1[12] = half_btf(kCosPi[50], bf0[12], kCosPi[14], bf0[13]);
  bf1[13] = half_btf(-kCosPi[50], bf0[13], kCosPi[14], bf0[12]);
  bf1[14] = half_btf(kCosPi[58], bf0[14], kCosPi[6], bf0[15]);
  bf1[15] = half_btf(-kCosPi[58], bf0[15], kCosPi[6], bf0[14]);

  // stage 3
  bf0 = step;
  bf1 = output;
  bf1[0] = bf0[0] + bf0[8];
  bf1[1] = bf0[1] + bf0[9];
  bf1[2] = bf0[2] + bf0[10];
  bf1[3] = bf0[3] + bf0[11];
  bf1[4] = bf0[4] + bf0[12];
  bf1[5] = bf0[5] + bf0[13];
  bf1[6] = bf0[6] + bf0[14];
  bf1[7] = bf0[7] + bf0[15];
  bf1[8] = -bf0[8] + bf0[0];
  bf1[9] = -bf0[9] + bf0[1];
  bf1[10] = -bf0[10] + bf0[2];
  bf1[11] = -bf0[11] + bf0[3];
  bf1[12] = -bf0[12] + bf0[4];
  bf1[13] = -bf0[13] + bf0[5];
  bf1[14] = -bf0[14] + bf0[6];
  bf1[15] = -bf0[15] + bf0[7];

  // stage 4
  bf0 = output;
  bf1 = step;
  bf1[0] = bf0[0];
  bf1[1] = bf0[1];
  bf1[2] = bf0[2];
  bf1[3] = bf0[3];
  bf1[4] = bf0[4];
  bf1[5] = bf0[5];
  bf1[6] = bf0[6];
  bf1[7] = bf0[7];
  bf1[8] = half_btf(kCosPi[8], bf0[8], kCosPi[56], bf0[9]);
  bf1[9] = half_btf(-kCosPi[8], bf0[9], kCosPi[56], bf0[8]);
  bf1[10] = half_btf(kCosPi[40], bf0[10], kCosPi[24], bf0[11]);
  bf1[11] = half_btf(-kCosPi[40], bf0[11], kCosPi[24], bf0[10]);
  bf1[12] = half_btf(-kCosPi[56], bf0[12], kCosPi[8], bf0[13]);
  bf1[13] = half_btf(kCosPi[56], bf0[13], kCosPi[8], bf0[12]);
  bf1[14] = half_btf(-kCosPi[24], bf0[14], kCosPi[40], bf0[15]);
  bf1[15] = half_btf(kCosPi[24], bf0[15], kCosPi[40], bf0[14]);

  // stage 5
  bf0 = step;
  bf1 = output;
  bf1[0] = bf0[0] + bf0[4];
  bf1[1] = bf0[1] + bf0[5];
  bf1[2] = bf0[2] + bf0[6];
  bf1[3] = bf0[3] + bf0[7];
  bf1[4] = -bf0[4] + bf0[0];
  bf1[5] = -bf0[5] + bf0[1];
  bf1[6] = -bf0[6] + bf0[2];
  bf1[7] = -bf0[7] + bf0[3];
  bf1[8] = bf0[8] + bf0[12];
  bf1[9] = bf0[9] + bf0[13];
  bf1[10] = bf0[10] + bf0[14];
  bf1[11] = bf0[11] + bf0[15];
  bf1[12] = -bf0[12] + bf0[8];
  bf1[13] = -bf0[13] + bf0[9];
  bf1[14] = -bf0[14] + bf0[10];
  bf1[15] = -bf0[15] + bf0[11];

  // stage 6
  bf0 = output;
  bf1 = step;
  bf1[0] = bf0[0];
  bf1[1] = bf0[1];
  bf1[2] = bf0[2];
  bf1[3] = bf0[3];
  bf1[4] = half_btf(kCosPi[16], bf0[4], kCosPi[48], bf0[5]);
  bf1[5] = half_btf(-kCosPi[16], bf0[5], kCosPi[48], bf0[4]);
  bf1[6] = half_btf(-kCosPi[48], bf0[6], kCosPi[16], bf0[7]);
  bf1[7] = half_btf(kCosPi[48], bf0[7], kCosPi[16], bf0[6]);
  bf1[8] = bf0[8];
  bf1[9] = bf0[9];
  bf1[10] = bf0[10];
  bf1[11] = bf0[11];
  bf1[12] = half_btf(kCosPi[16], bf0[12], kCosPi[48], bf0[13]);
  bf1[13] = half_btf(-kCosPi[16], bf0[13], kCosPi[48], bf0[12]);
  bf1[14] = half_btf(-kCosPi[48], bf0[14], kCosPi[16], bf0[15]);
  bf1[15] = half_btf(kCosPi[48], bf0[15], kCosPi[16], bf0[14]);

  // stage 7
  bf0 = step;
  bf1 = output;
  bf1[0] = bf0[0] + bf0[2];
  bf1[1] = bf0[1] + bf0[3];
  bf1[2] = -bf0[2] + bf0[0];
  bf1[3] = -bf0[3] + bf0[1];
  bf1[4] = bf0[4] + bf0[6];
  bf1[5] = bf0[5] + bf0[7];
  bf1[6] = -bf0[6] + bf0[4];
  bf1[7] = -bf0[7] + bf0[5];
  bf1[8] = bf0[8] + bf0[10];
  bf1[9] = bf0[9] + bf0[11];
  bf1[10] = -bf0[10] + bf0[8];
  bf1[11] = -bf0[11] + bf0[9];
  bf1[12] = bf0[12] + bf0[14];
  bf1[13] = bf0[13] + bf0[15];
  bf1[14] = -bf0[14] + bf0[12];
  bf1[15] = -bf0[15] + bf0[13];

  // stage 8
  bf0 = output;
  bf1 = step;
  bf1[0] = div_sqrt2(bf0[0], shift - kPrecision);
  bf1[1] = div_sqrt2(bf0[1], shift - kPrecision);
  bf1[2] = round_shift(bf0[2] + bf0[3], shift - kPrecision + 1);
  bf1[3] = round_shift(bf0[2] - bf0[3], shift - kPrecision + 1);
  bf1[4] = div_sqrt2(bf0[4], shift - kPrecision);
  bf1[5] = div_sqrt2(bf0[5], shift - kPrecision);
  bf1[6] = round_shift(bf0[6] + bf0[7], shift - kPrecision + 1);
  bf1[7] = round_shift(bf0[6] - bf0[7], shift - kPrecision + 1);
  bf1[8] = div_sqrt2(bf0[8], shift - kPrecision);
  bf1[9] = div_sqrt2(bf0[9], shift - kPrecision);
  bf1[10] = round_shift(bf0[10] + bf0[11], shift - kPrecision + 1);
  bf1[11] = round_shift(bf0[10] - bf0[11], shift - kPrecision + 1);
  bf1[12] = div_sqrt2(bf0[12], shift - kPrecision);
  bf1[13] =  div_sqrt2(bf0[13], shift - kPrecision);
  bf1[14] = round_shift(bf0[14] + bf0[15], shift - kPrecision + 1);
  bf1[15] = round_shift(bf0[14] - bf0[15], shift - kPrecision + 1);

  // stage 9
  bf0 = step;
  bf1 = output;
  bf1[0] = bf0[0];
  bf1[1] = -bf0[8];
  bf1[2] = bf0[12];
  bf1[3] = -bf0[4];
  bf1[4] = bf0[6];
  bf1[5] = -bf0[14];
  bf1[6] = bf0[10];
  bf1[7] = -bf0[2];
  bf1[8] = bf0[3];
  bf1[9] = -bf0[11];
  bf1[10] = bf0[15];
  bf1[11] = -bf0[7];
  bf1[12] = bf0[5];
  bf1[13] = -bf0[13];
  bf1[14] = bf0[9];
  bf1[15] = -bf0[1];
}

static void iadst32_C(const int32_t* input, int32_t* output) {
  constexpr uint32_t shift = kPrecision + 2;
  int32_t *bf0, *bf1;
  int32_t step[32];

  // stage 1;
  assert(output != input);
  bf1 = output;
  bf1[0] = input[0];
  bf1[1] = -input[31];
  bf1[2] = -input[15];
  bf1[3] = input[16];
  bf1[4] = -input[7];
  bf1[5] = input[24];
  bf1[6] = input[8];
  bf1[7] = -input[23];
  bf1[8] = -input[3];
  bf1[9] = input[28];
  bf1[10] = input[12];
  bf1[11] = -input[19];
  bf1[12] = input[4];
  bf1[13] = -input[27];
  bf1[14] = -input[11];
  bf1[15] = input[20];
  bf1[16] = -input[1];
  bf1[17] = input[30];
  bf1[18] = input[14];
  bf1[19] = -input[17];
  bf1[20] = input[6];
  bf1[21] = -input[25];
  bf1[22] = -input[9];
  bf1[23] = input[22];
  bf1[24] = input[2];
  bf1[25] = -input[29];
  bf1[26] = -input[13];
  bf1[27] = input[18];
  bf1[28] = -input[5];
  bf1[29] = input[26];
  bf1[30] = input[10];
  bf1[31] = -input[21];

  // stage 2
  bf0 = output;
  bf1 = step;
  bf1[0] = bf0[0];
  bf1[1] = bf0[1];
  bf1[2] = half_btf(kCosPi[32], bf0[2], kCosPi[32], bf0[3]);
  bf1[3] = half_btf(kCosPi[32], bf0[2], -kCosPi[32], bf0[3]);
  bf1[4] = bf0[4];
  bf1[5] = bf0[5];
  bf1[6] = half_btf(kCosPi[32], bf0[6], kCosPi[32], bf0[7]);
  bf1[7] = half_btf(kCosPi[32], bf0[6], -kCosPi[32], bf0[7]);
  bf1[8] = bf0[8];
  bf1[9] = bf0[9];
  bf1[10] = half_btf(kCosPi[32], bf0[10], kCosPi[32], bf0[11]);
  bf1[11] = half_btf(kCosPi[32], bf0[10], -kCosPi[32], bf0[11]);
  bf1[12] = bf0[12];
  bf1[13] = bf0[13];
  bf1[14] = half_btf(kCosPi[32], bf0[14], kCosPi[32], bf0[15]);
  bf1[15] = half_btf(kCosPi[32], bf0[14], -kCosPi[32], bf0[15]);
  bf1[16] = bf0[16];
  bf1[17] = bf0[17];
  bf1[18] = half_btf(kCosPi[32], bf0[18], kCosPi[32], bf0[19]);
  bf1[19] = half_btf(kCosPi[32], bf0[18], -kCosPi[32], bf0[19]);
  bf1[20] = bf0[20];
  bf1[21] = bf0[21];
  bf1[22] = half_btf(kCosPi[32], bf0[22], kCosPi[32], bf0[23]);
  bf1[23] = half_btf(kCosPi[32], bf0[22], -kCosPi[32], bf0[23]);
  bf1[24] = bf0[24];
  bf1[25] = bf0[25];
  bf1[26] = half_btf(kCosPi[32], bf0[26], kCosPi[32], bf0[27]);
  bf1[27] = half_btf(kCosPi[32], bf0[26], -kCosPi[32], bf0[27]);
  bf1[28] = bf0[28];
  bf1[29] = bf0[29];
  bf1[30] = half_btf(kCosPi[32], bf0[30], kCosPi[32], bf0[31]);
  bf1[31] = half_btf(kCosPi[32], bf0[30], -kCosPi[32], bf0[31]);

  // stage 3
  bf0 = step;
  bf1 = output;
  bf1[0] = bf0[0] + bf0[2];
  bf1[1] = bf0[1] + bf0[3];
  bf1[2] = bf0[0] - bf0[2];
  bf1[3] = bf0[1] - bf0[3];
  bf1[4] = bf0[4] + bf0[6];
  bf1[5] = bf0[5] + bf0[7];
  bf1[6] = bf0[4] - bf0[6];
  bf1[7] = bf0[5] - bf0[7];
  bf1[8] = bf0[8] + bf0[10];
  bf1[9] = bf0[9] + bf0[11];
  bf1[10] = bf0[8] - bf0[10];
  bf1[11] = bf0[9] - bf0[11];
  bf1[12] = bf0[12] + bf0[14];
  bf1[13] = bf0[13] + bf0[15];
  bf1[14] = bf0[12] - bf0[14];
  bf1[15] = bf0[13] - bf0[15];
  bf1[16] = bf0[16] + bf0[18];
  bf1[17] = bf0[17] + bf0[19];
  bf1[18] = bf0[16] - bf0[18];
  bf1[19] = bf0[17] - bf0[19];
  bf1[20] = bf0[20] + bf0[22];
  bf1[21] = bf0[21] + bf0[23];
  bf1[22] = bf0[20] - bf0[22];
  bf1[23] = bf0[21] - bf0[23];
  bf1[24] = bf0[24] + bf0[26];
  bf1[25] = bf0[25] + bf0[27];
  bf1[26] = bf0[24] - bf0[26];
  bf1[27] = bf0[25] - bf0[27];
  bf1[28] = bf0[28] + bf0[30];
  bf1[29] = bf0[29] + bf0[31];
  bf1[30] = bf0[28] - bf0[30];
  bf1[31] = bf0[29] - bf0[31];

  // stage 4
  bf0 = output;
  bf1 = step;
  bf1[0] = bf0[0];
  bf1[1] = bf0[1];
  bf1[2] = bf0[2];
  bf1[3] = bf0[3];
  bf1[4] = half_btf(kCosPi[16], bf0[4], kCosPi[48], bf0[5]);
  bf1[5] = half_btf(kCosPi[48], bf0[4], -kCosPi[16], bf0[5]);
  bf1[6] = half_btf(-kCosPi[48], bf0[6], kCosPi[16], bf0[7]);
  bf1[7] = half_btf(kCosPi[16], bf0[6], kCosPi[48], bf0[7]);
  bf1[8] = bf0[8];
  bf1[9] = bf0[9];
  bf1[10] = bf0[10];
  bf1[11] = bf0[11];
  bf1[12] = half_btf(kCosPi[16], bf0[12], kCosPi[48], bf0[13]);
  bf1[13] = half_btf(kCosPi[48], bf0[12], -kCosPi[16], bf0[13]);
  bf1[14] = half_btf(-kCosPi[48], bf0[14], kCosPi[16], bf0[15]);
  bf1[15] = half_btf(kCosPi[16], bf0[14], kCosPi[48], bf0[15]);
  bf1[16] = bf0[16];
  bf1[17] = bf0[17];
  bf1[18] = bf0[18];
  bf1[19] = bf0[19];
  bf1[20] = half_btf(kCosPi[16], bf0[20], kCosPi[48], bf0[21]);
  bf1[21] = half_btf(kCosPi[48], bf0[20], -kCosPi[16], bf0[21]);
  bf1[22] = half_btf(-kCosPi[48], bf0[22], kCosPi[16], bf0[23]);
  bf1[23] = half_btf(kCosPi[16], bf0[22], kCosPi[48], bf0[23]);
  bf1[24] = bf0[24];
  bf1[25] = bf0[25];
  bf1[26] = bf0[26];
  bf1[27] = bf0[27];
  bf1[28] = half_btf(kCosPi[16], bf0[28], kCosPi[48], bf0[29]);
  bf1[29] = half_btf(kCosPi[48], bf0[28], -kCosPi[16], bf0[29]);
  bf1[30] = half_btf(-kCosPi[48], bf0[30], kCosPi[16], bf0[31]);
  bf1[31] = half_btf(kCosPi[16], bf0[30], kCosPi[48], bf0[31]);

  // stage 5
  bf0 = step;
  bf1 = output;
  bf1[0] = bf0[0] + bf0[4];
  bf1[1] = bf0[1] + bf0[5];
  bf1[2] = bf0[2] + bf0[6];
  bf1[3] = bf0[3] + bf0[7];
  bf1[4] = bf0[0] - bf0[4];
  bf1[5] = bf0[1] - bf0[5];
  bf1[6] = bf0[2] - bf0[6];
  bf1[7] = bf0[3] - bf0[7];
  bf1[8] = bf0[8] + bf0[12];
  bf1[9] = bf0[9] + bf0[13];
  bf1[10] = bf0[10] + bf0[14];
  bf1[11] = bf0[11] + bf0[15];
  bf1[12] = bf0[8] - bf0[12];
  bf1[13] = bf0[9] - bf0[13];
  bf1[14] = bf0[10] - bf0[14];
  bf1[15] = bf0[11] - bf0[15];
  bf1[16] = bf0[16] + bf0[20];
  bf1[17] = bf0[17] + bf0[21];
  bf1[18] = bf0[18] + bf0[22];
  bf1[19] = bf0[19] + bf0[23];
  bf1[20] = bf0[16] - bf0[20];
  bf1[21] = bf0[17] - bf0[21];
  bf1[22] = bf0[18] - bf0[22];
  bf1[23] = bf0[19] - bf0[23];
  bf1[24] = bf0[24] + bf0[28];
  bf1[25] = bf0[25] + bf0[29];
  bf1[26] = bf0[26] + bf0[30];
  bf1[27] = bf0[27] + bf0[31];
  bf1[28] = bf0[24] - bf0[28];
  bf1[29] = bf0[25] - bf0[29];
  bf1[30] = bf0[26] - bf0[30];
  bf1[31] = bf0[27] - bf0[31];

  // stage 6
  bf0 = output;
  bf1 = step;
  bf1[0] = bf0[0];
  bf1[1] = bf0[1];
  bf1[2] = bf0[2];
  bf1[3] = bf0[3];
  bf1[4] = bf0[4];
  bf1[5] = bf0[5];
  bf1[6] = bf0[6];
  bf1[7] = bf0[7];
  bf1[8] = half_btf(kCosPi[8], bf0[8], kCosPi[56], bf0[9]);
  bf1[9] = half_btf(kCosPi[56], bf0[8], -kCosPi[8], bf0[9]);
  bf1[10] = half_btf(kCosPi[40], bf0[10], kCosPi[24], bf0[11]);
  bf1[11] = half_btf(kCosPi[24], bf0[10], -kCosPi[40], bf0[11]);
  bf1[12] = half_btf(-kCosPi[56], bf0[12], kCosPi[8], bf0[13]);
  bf1[13] = half_btf(kCosPi[8], bf0[12], kCosPi[56], bf0[13]);
  bf1[14] = half_btf(-kCosPi[24], bf0[14], kCosPi[40], bf0[15]);
  bf1[15] = half_btf(kCosPi[40], bf0[14], kCosPi[24], bf0[15]);
  bf1[16] = bf0[16];
  bf1[17] = bf0[17];
  bf1[18] = bf0[18];
  bf1[19] = bf0[19];
  bf1[20] = bf0[20];
  bf1[21] = bf0[21];
  bf1[22] = bf0[22];
  bf1[23] = bf0[23];
  bf1[24] = half_btf(kCosPi[8], bf0[24], kCosPi[56], bf0[25]);
  bf1[25] = half_btf(kCosPi[56], bf0[24], -kCosPi[8], bf0[25]);
  bf1[26] = half_btf(kCosPi[40], bf0[26], kCosPi[24], bf0[27]);
  bf1[27] = half_btf(kCosPi[24], bf0[26], -kCosPi[40], bf0[27]);
  bf1[28] = half_btf(-kCosPi[56], bf0[28], kCosPi[8], bf0[29]);
  bf1[29] = half_btf(kCosPi[8], bf0[28], kCosPi[56], bf0[29]);
  bf1[30] = half_btf(-kCosPi[24], bf0[30], kCosPi[40], bf0[31]);
  bf1[31] = half_btf(kCosPi[40], bf0[30], kCosPi[24], bf0[31]);

  // stage 7
  bf0 = step;
  bf1 = output;
  bf1[0] = bf0[0] + bf0[8];
  bf1[1] = bf0[1] + bf0[9];
  bf1[2] = bf0[2] + bf0[10];
  bf1[3] = bf0[3] + bf0[11];
  bf1[4] = bf0[4] + bf0[12];
  bf1[5] = bf0[5] + bf0[13];
  bf1[6] = bf0[6] + bf0[14];
  bf1[7] = bf0[7] + bf0[15];
  bf1[8] = bf0[0] - bf0[8];
  bf1[9] = bf0[1] - bf0[9];
  bf1[10] = bf0[2] - bf0[10];
  bf1[11] = bf0[3] - bf0[11];
  bf1[12] = bf0[4] - bf0[12];
  bf1[13] = bf0[5] - bf0[13];
  bf1[14] = bf0[6] - bf0[14];
  bf1[15] = bf0[7] - bf0[15];
  bf1[16] = bf0[16] + bf0[24];
  bf1[17] = bf0[17] + bf0[25];
  bf1[18] = bf0[18] + bf0[26];
  bf1[19] = bf0[19] + bf0[27];
  bf1[20] = bf0[20] + bf0[28];
  bf1[21] = bf0[21] + bf0[29];
  bf1[22] = bf0[22] + bf0[30];
  bf1[23] = bf0[23] + bf0[31];
  bf1[24] = bf0[16] - bf0[24];
  bf1[25] = bf0[17] - bf0[25];
  bf1[26] = bf0[18] - bf0[26];
  bf1[27] = bf0[19] - bf0[27];
  bf1[28] = bf0[20] - bf0[28];
  bf1[29] = bf0[21] - bf0[29];
  bf1[30] = bf0[22] - bf0[30];
  bf1[31] = bf0[23] - bf0[31];

  // stage 8
  bf0 = output;
  bf1 = step;
  bf1[0] = bf0[0];
  bf1[1] = bf0[1];
  bf1[2] = bf0[2];
  bf1[3] = bf0[3];
  bf1[4] = bf0[4];
  bf1[5] = bf0[5];
  bf1[6] = bf0[6];
  bf1[7] = bf0[7];
  bf1[8] = bf0[8];
  bf1[9] = bf0[9];
  bf1[10] = bf0[10];
  bf1[11] = bf0[11];
  bf1[12] = bf0[12];
  bf1[13] = bf0[13];
  bf1[14] = bf0[14];
  bf1[15] = bf0[15];
  bf1[16] = half_btf(kCosPi[4], bf0[16], kCosPi[60], bf0[17]);
  bf1[17] = half_btf(kCosPi[60], bf0[16], -kCosPi[4], bf0[17]);
  bf1[18] = half_btf(kCosPi[20], bf0[18], kCosPi[44], bf0[19]);
  bf1[19] = half_btf(kCosPi[44], bf0[18], -kCosPi[20], bf0[19]);
  bf1[20] = half_btf(kCosPi[36], bf0[20], kCosPi[28], bf0[21]);
  bf1[21] = half_btf(kCosPi[28], bf0[20], -kCosPi[36], bf0[21]);
  bf1[22] = half_btf(kCosPi[52], bf0[22], kCosPi[12], bf0[23]);
  bf1[23] = half_btf(kCosPi[12], bf0[22], -kCosPi[52], bf0[23]);
  bf1[24] = half_btf(-kCosPi[60], bf0[24], kCosPi[4], bf0[25]);
  bf1[25] = half_btf(kCosPi[4], bf0[24], kCosPi[60], bf0[25]);
  bf1[26] = half_btf(-kCosPi[44], bf0[26], kCosPi[20], bf0[27]);
  bf1[27] = half_btf(kCosPi[20], bf0[26], kCosPi[44], bf0[27]);
  bf1[28] = half_btf(-kCosPi[28], bf0[28], kCosPi[36], bf0[29]);
  bf1[29] = half_btf(kCosPi[36], bf0[28], kCosPi[28], bf0[29]);
  bf1[30] = half_btf(-kCosPi[12], bf0[30], kCosPi[52], bf0[31]);
  bf1[31] = half_btf(kCosPi[52], bf0[30], kCosPi[12], bf0[31]);

  // stage 9
  bf0 = step;
  bf1 = output;
  bf1[0] = bf0[0] + bf0[16];
  bf1[1] = bf0[1] + bf0[17];
  bf1[2] = bf0[2] + bf0[18];
  bf1[3] = bf0[3] + bf0[19];
  bf1[4] = bf0[4] + bf0[20];
  bf1[5] = bf0[5] + bf0[21];
  bf1[6] = bf0[6] + bf0[22];
  bf1[7] = bf0[7] + bf0[23];
  bf1[8] = bf0[8] + bf0[24];
  bf1[9] = bf0[9] + bf0[25];
  bf1[10] = bf0[10] + bf0[26];
  bf1[11] = bf0[11] + bf0[27];
  bf1[12] = bf0[12] + bf0[28];
  bf1[13] = bf0[13] + bf0[29];
  bf1[14] = bf0[14] + bf0[30];
  bf1[15] = bf0[15] + bf0[31];
  bf1[16] = bf0[0] - bf0[16];
  bf1[17] = bf0[1] - bf0[17];
  bf1[18] = bf0[2] - bf0[18];
  bf1[19] = bf0[3] - bf0[19];
  bf1[20] = bf0[4] - bf0[20];
  bf1[21] = bf0[5] - bf0[21];
  bf1[22] = bf0[6] - bf0[22];
  bf1[23] = bf0[7] - bf0[23];
  bf1[24] = bf0[8] - bf0[24];
  bf1[25] = bf0[9] - bf0[25];
  bf1[26] = bf0[10] - bf0[26];
  bf1[27] = bf0[11] - bf0[27];
  bf1[28] = bf0[12] - bf0[28];
  bf1[29] = bf0[13] - bf0[29];
  bf1[30] = bf0[14] - bf0[30];
  bf1[31] = bf0[15] - bf0[31];

  // stage 10
  bf0 = output;
  bf1 = step;
  bf1[0] = half_btf(kCosPi[1], bf0[0], kCosPi[63], bf0[1], shift);
  bf1[1] = half_btf(kCosPi[63], bf0[0], -kCosPi[1], bf0[1], shift);
  bf1[2] = half_btf(kCosPi[5], bf0[2], kCosPi[59], bf0[3], shift);
  bf1[3] = half_btf(kCosPi[59], bf0[2], -kCosPi[5], bf0[3], shift);
  bf1[4] = half_btf(kCosPi[9], bf0[4], kCosPi[55], bf0[5], shift);
  bf1[5] = half_btf(kCosPi[55], bf0[4], -kCosPi[9], bf0[5], shift);
  bf1[6] = half_btf(kCosPi[13], bf0[6], kCosPi[51], bf0[7], shift);
  bf1[7] = half_btf(kCosPi[51], bf0[6], -kCosPi[13], bf0[7], shift);
  bf1[8] = half_btf(kCosPi[17], bf0[8], kCosPi[47], bf0[9], shift);
  bf1[9] = half_btf(kCosPi[47], bf0[8], -kCosPi[17], bf0[9], shift);
  bf1[10] = half_btf(kCosPi[21], bf0[10], kCosPi[43], bf0[11], shift);
  bf1[11] = half_btf(kCosPi[43], bf0[10], -kCosPi[21], bf0[11], shift);
  bf1[12] = half_btf(kCosPi[25], bf0[12], kCosPi[39], bf0[13], shift);
  bf1[13] = half_btf(kCosPi[39], bf0[12], -kCosPi[25], bf0[13], shift);
  bf1[14] = half_btf(kCosPi[29], bf0[14], kCosPi[35], bf0[15], shift);
  bf1[15] = half_btf(kCosPi[35], bf0[14], -kCosPi[29], bf0[15], shift);
  bf1[16] = half_btf(kCosPi[33], bf0[16], kCosPi[31], bf0[17], shift);
  bf1[17] = half_btf(kCosPi[31], bf0[16], -kCosPi[33], bf0[17], shift);
  bf1[18] = half_btf(kCosPi[37], bf0[18], kCosPi[27], bf0[19], shift);
  bf1[19] = half_btf(kCosPi[27], bf0[18], -kCosPi[37], bf0[19], shift);
  bf1[20] = half_btf(kCosPi[41], bf0[20], kCosPi[23], bf0[21], shift);
  bf1[21] = half_btf(kCosPi[23], bf0[20], -kCosPi[41], bf0[21], shift);
  bf1[22] = half_btf(kCosPi[45], bf0[22], kCosPi[19], bf0[23], shift);
  bf1[23] = half_btf(kCosPi[19], bf0[22], -kCosPi[45], bf0[23], shift);
  bf1[24] = half_btf(kCosPi[49], bf0[24], kCosPi[15], bf0[25], shift);
  bf1[25] = half_btf(kCosPi[15], bf0[24], -kCosPi[49], bf0[25], shift);
  bf1[26] = half_btf(kCosPi[53], bf0[26], kCosPi[11], bf0[27], shift);
  bf1[27] = half_btf(kCosPi[11], bf0[26], -kCosPi[53], bf0[27], shift);
  bf1[28] = half_btf(kCosPi[57], bf0[28], kCosPi[7], bf0[29], shift);
  bf1[29] = half_btf(kCosPi[7], bf0[28], -kCosPi[57], bf0[29], shift);
  bf1[30] = half_btf(kCosPi[61], bf0[30], kCosPi[3], bf0[31], shift);
  bf1[31] = half_btf(kCosPi[3], bf0[30], -kCosPi[61], bf0[31], shift);

  // stage 11
  bf0 = step;
  bf1 = output;
  bf1[0] = bf0[1];
  bf1[1] = bf0[30];
  bf1[2] = bf0[3];
  bf1[3] = bf0[28];
  bf1[4] = bf0[5];
  bf1[5] = bf0[26];
  bf1[6] = bf0[7];
  bf1[7] = bf0[24];
  bf1[8] = bf0[9];
  bf1[9] = bf0[22];
  bf1[10] = bf0[11];
  bf1[11] = bf0[20];
  bf1[12] = bf0[13];
  bf1[13] = bf0[18];
  bf1[14] = bf0[15];
  bf1[15] = bf0[16];
  bf1[16] = bf0[17];
  bf1[17] = bf0[14];
  bf1[18] = bf0[19];
  bf1[19] = bf0[12];
  bf1[20] = bf0[21];
  bf1[21] = bf0[10];
  bf1[22] = bf0[23];
  bf1[23] = bf0[8];
  bf1[24] = bf0[25];
  bf1[25] = bf0[6];
  bf1[26] = bf0[27];
  bf1[27] = bf0[4];
  bf1[28] = bf0[29];
  bf1[29] = bf0[2];
  bf1[30] = bf0[31];
  bf1[31] = bf0[0];
}

//------------------------------------------------------------------------------
// Hadamard transform

// sequency reordering (=bit-reversal + Gray code)
static const uint32_t kHadamardR2[2] = { 0, 1 };
static const uint32_t kHadamardR4[4] = { 0, 3, 1, 2 };
static const uint32_t kHadamardR8[8] = { 0, 7, 3, 4, 1, 6, 2, 5 };
static const uint32_t kHadamardR16[16] = {  // https://oeis.org/A240909
  0, 15, 7, 8, 3, 12, 4, 11, 1, 14, 6, 9, 2, 13, 5, 10
};
static const uint32_t kHadamardR32[32] = {  // https://oeis.org/A240910
  0, 31, 15, 16, 7, 24, 8, 23, 3, 28, 12, 19, 4, 27, 11, 20,
  1, 30, 14, 17, 6, 25, 9, 22, 2, 29, 13, 18, 5, 26, 10, 21
};

static const uint32_t* const kHadamardReorder[5] = {
  kHadamardR2, kHadamardR4, kHadamardR8, kHadamardR16, kHadamardR32
};

static void CoreHadamard_C(const int32_t* src, int32_t* dst, uint32_t size,
                           int8_t divider = 1) {
  assert(size == 2 || size == 4 || size == 8 || size == 16 || size == 32);
  int32_t tmp[32];
  for (uint32_t half_size = 1; half_size < size; half_size *= 2) {
    for (uint32_t i = 0; i < size; i += 2 * half_size) {
      for (uint32_t j = 0; j < half_size; ++j) {
        const int32_t A = src[i + j +         0];
        const int32_t B = src[i + j + half_size];
        tmp[i + j +         0] = A + B;
        tmp[i + j + half_size] = A - B;
      }
    }
    src = tmp;
  }
  // Sequency-reordering
  const uint32_t* const reorder = kHadamardReorder[WP2Log2Floor(size) - 1];
  // note we could have a different norm for fwd and bwd transform, as long as
  // their product is equal to 1/2. Using ">>1" for descale gives an exact
  // transform. Using 1/sqrt(2) for both is a good compromise of dynamic range
  // compared to the other transforms (dct, adst, ...).
  const float norm = 0.70710678118f / divider;
  for (uint32_t i = 0; i < size; ++i) dst[i] = tmp[reorder[i]] * norm;
}

#define GENERATE_HADAMARD_FUNC(SIZE)                               \
static void iHadamard##SIZE##_C(const int32_t* src, int32_t* dst) {   \
  CoreHadamard_C(src, dst, SIZE, SIZE / 2);                        \
}                                                                  \
static void fHadamard##SIZE##_C(const int32_t* src, int32_t* dst) {   \
  CoreHadamard_C(src, dst, SIZE, 1);                               \
}

GENERATE_HADAMARD_FUNC(2)
GENERATE_HADAMARD_FUNC(4)
GENERATE_HADAMARD_FUNC(8)
GENERATE_HADAMARD_FUNC(16)
GENERATE_HADAMARD_FUNC(32)

#undef GENERATE_HADAMARD_FUNC

//------------------------------------------------------------------------------
// Forward transforms

static void fdct2_C(const int32_t* input, int32_t* output) {
  const int32_t A = input[0], B = input[1];
  output[0] = div_sqrt2(A + B);
  output[1] = div_sqrt2(A - B);
}

static void fdct4_C(const int32_t* input, int32_t* output) {
  int32_t *bf0, *bf1;
  int32_t step[4];

  // stage 1;
  bf1 = output;
  bf1[0] = input[0] + input[3];
  bf1[1] = input[1] + input[2];
  bf1[2] = input[1] - input[2];
  bf1[3] = input[0] - input[3];

  // stage 2
  bf0 = output;
  bf1 = step;
  bf1[0] = round_shift(bf0[0] + bf0[1], 1);
  bf1[1] = round_shift(bf0[0] - bf0[1], 1);
  bf1[2] = half_btf( kCos32Cos48, bf0[2], kCos32Cos16, bf0[3]);
  bf1[3] = half_btf(-kCos32Cos16, bf0[2], kCos32Cos48, bf0[3]);

  // stage 3
  bf0 = step;
  bf1 = output;
  bf1[0] = bf0[0];
  bf1[1] = bf0[2];
  bf1[2] = bf0[1];
  bf1[3] = bf0[3];
}

static void fdct8_C(const int32_t* input, int32_t* output) {
  constexpr uint32_t shift = kPrecision + 1;
  int32_t *bf0, *bf1;
  int32_t step[8];

  // stage 1;
  bf1 = output;
  bf1[0] = input[0] + input[7];
  bf1[1] = input[1] + input[6];
  bf1[2] = input[2] + input[5];
  bf1[3] = input[3] + input[4];
  bf1[4] = -input[4] + input[3];
  bf1[5] = -input[5] + input[2];
  bf1[6] = -input[6] + input[1];
  bf1[7] = -input[7] + input[0];

  // stage 2
  bf0 = output;
  bf1 = step;
  bf1[0] = bf0[0] + bf0[3];
  bf1[1] = bf0[1] + bf0[2];
  bf1[2] = -bf0[2] + bf0[1];
  bf1[3] = -bf0[3] + bf0[0];
  bf1[4] = bf0[4];
  bf1[5] = half_btf(-kCosPi[32], bf0[5], kCosPi[32], bf0[6]);
  bf1[6] = half_btf(kCosPi[32], bf0[6], kCosPi[32], bf0[5]);
  bf1[7] = bf0[7];

  // stage 3
  bf0 = step;
  bf1 = output;
  bf1[0] = half_btf(kCosPi[32], bf0[0], kCosPi[32], bf0[1], shift);
  bf1[1] = half_btf(-kCosPi[32], bf0[1], kCosPi[32], bf0[0], shift);
  bf1[2] = half_btf(kCosPi[48], bf0[2], kCosPi[16], bf0[3], shift);
  bf1[3] = half_btf(kCosPi[48], bf0[3], -kCosPi[16], bf0[2], shift);
  bf1[4] = bf0[4] + bf0[5];
  bf1[5] = -bf0[5] + bf0[4];
  bf1[6] = -bf0[6] + bf0[7];
  bf1[7] = bf0[7] + bf0[6];

  // stage 4
  bf0 = output;
  bf1 = step;
  bf1[0] = bf0[0];
  bf1[1] = bf0[1];
  bf1[2] = bf0[2];
  bf1[3] = bf0[3];
  bf1[4] = half_btf(kCosPi[56], bf0[4], kCosPi[8], bf0[7], shift);
  bf1[5] = half_btf(kCosPi[24], bf0[5], kCosPi[40], bf0[6], shift);
  bf1[6] = half_btf(kCosPi[24], bf0[6], -kCosPi[40], bf0[5], shift);
  bf1[7] = half_btf(kCosPi[56], bf0[7], -kCosPi[8], bf0[4], shift);

  // stage 5
  bf0 = step;
  bf1 = output;
  bf1[0] = bf0[0];
  bf1[1] = bf0[4];
  bf1[2] = bf0[2];
  bf1[3] = bf0[6];
  bf1[4] = bf0[1];
  bf1[5] = bf0[5];
  bf1[6] = bf0[3];
  bf1[7] = bf0[7];
}

static void fdct16_C(const int32_t* input, int32_t* output) {
  int32_t *bf0, *bf1;
  int32_t step[16];

  // stage 1;
  bf1 = output;
  bf1[0] = input[0] + input[15];
  bf1[1] = input[1] + input[14];
  bf1[2] = input[2] + input[13];
  bf1[3] = input[3] + input[12];
  bf1[4] = input[4] + input[11];
  bf1[5] = input[5] + input[10];
  bf1[6] = input[6] + input[9];
  bf1[7] = input[7] + input[8];
  bf1[8] = -input[8] + input[7];
  bf1[9] = -input[9] + input[6];
  bf1[10] = -input[10] + input[5];
  bf1[11] = -input[11] + input[4];
  bf1[12] = -input[12] + input[3];
  bf1[13] = -input[13] + input[2];
  bf1[14] = -input[14] + input[1];
  bf1[15] = -input[15] + input[0];

  // stage 2
  bf0 = output;
  bf1 = step;
  bf1[0] = bf0[0] + bf0[7];
  bf1[1] = bf0[1] + bf0[6];
  bf1[2] = bf0[2] + bf0[5];
  bf1[3] = bf0[3] + bf0[4];
  bf1[4] = -bf0[4] + bf0[3];
  bf1[5] = -bf0[5] + bf0[2];
  bf1[6] = -bf0[6] + bf0[1];
  bf1[7] = -bf0[7] + bf0[0];
  bf1[8] = bf0[8];
  bf1[9] = bf0[9];
  bf1[10] = half_btf(-kCosPi[32], bf0[10], kCosPi[32], bf0[13]);
  bf1[11] = half_btf(-kCosPi[32], bf0[11], kCosPi[32], bf0[12]);
  bf1[12] = half_btf(kCosPi[32], bf0[12], kCosPi[32], bf0[11]);
  bf1[13] = half_btf(kCosPi[32], bf0[13], kCosPi[32], bf0[10]);
  bf1[14] = bf0[14];
  bf1[15] = bf0[15];

  // stage 3
  bf0 = step;
  bf1 = output;
  bf1[0] = bf0[0] + bf0[3];
  bf1[1] = bf0[1] + bf0[2];
  bf1[2] = -bf0[2] + bf0[1];
  bf1[3] = -bf0[3] + bf0[0];
  bf1[4] = bf0[4];
  bf1[5] = half_btf(-kCosPi[32], bf0[5], kCosPi[32], bf0[6]);
  bf1[6] = half_btf(kCosPi[32], bf0[6], kCosPi[32], bf0[5]);
  bf1[7] = bf0[7];
  bf1[8] = bf0[8] + bf0[11];
  bf1[9] = bf0[9] + bf0[10];
  bf1[10] = -bf0[10] + bf0[9];
  bf1[11] = -bf0[11] + bf0[8];
  bf1[12] = -bf0[12] + bf0[15];
  bf1[13] = -bf0[13] + bf0[14];
  bf1[14] = bf0[14] + bf0[13];
  bf1[15] = bf0[15] + bf0[12];

  // stage 4
  bf0 = output;
  bf1 = step;
  bf1[0] = half_btf(kCosPi[32], bf0[0], kCosPi[32], bf0[1]);
  bf1[1] = half_btf(-kCosPi[32], bf0[1], kCosPi[32], bf0[0]);
  bf1[2] = half_btf(kCosPi[48], bf0[2], kCosPi[16], bf0[3]);
  bf1[3] = half_btf(kCosPi[48], bf0[3], -kCosPi[16], bf0[2]);
  bf1[4] = bf0[4] + bf0[5];
  bf1[5] = -bf0[5] + bf0[4];
  bf1[6] = -bf0[6] + bf0[7];
  bf1[7] = bf0[7] + bf0[6];
  bf1[8] = bf0[8];
  bf1[9] = half_btf(-kCosPi[16], bf0[9], kCosPi[48], bf0[14]);
  bf1[10] = half_btf(-kCosPi[48], bf0[10], -kCosPi[16], bf0[13]);
  bf1[11] = bf0[11];
  bf1[12] = bf0[12];
  bf1[13] = half_btf(kCosPi[48], bf0[13], -kCosPi[16], bf0[10]);
  bf1[14] = half_btf(kCosPi[16], bf0[14], kCosPi[48], bf0[9]);
  bf1[15] = bf0[15];

  // stage 5
  bf0 = step;
  bf1 = output;
  bf1[0] = bf0[0];
  bf1[1] = bf0[1];
  bf1[2] = bf0[2];
  bf1[3] = bf0[3];
  bf1[4] = half_btf(kCosPi[56], bf0[4], kCosPi[8], bf0[7]);
  bf1[5] = half_btf(kCosPi[24], bf0[5], kCosPi[40], bf0[6]);
  bf1[6] = half_btf(kCosPi[24], bf0[6], -kCosPi[40], bf0[5]);
  bf1[7] = half_btf(kCosPi[56], bf0[7], -kCosPi[8], bf0[4]);
  bf1[8] = bf0[8] + bf0[9];
  bf1[9] = -bf0[9] + bf0[8];
  bf1[10] = -bf0[10] + bf0[11];
  bf1[11] = bf0[11] + bf0[10];
  bf1[12] = bf0[12] + bf0[13];
  bf1[13] = -bf0[13] + bf0[12];
  bf1[14] = -bf0[14] + bf0[15];
  bf1[15] = bf0[15] + bf0[14];

  // stage 6
  bf0 = output;
  bf1 = step;
  bf1[0] = bf0[0];
  bf1[1] = bf0[1];
  bf1[2] = bf0[2];
  bf1[3] = bf0[3];
  bf1[4] = bf0[4];
  bf1[5] = bf0[5];
  bf1[6] = bf0[6];
  bf1[7] = bf0[7];
  bf1[8] = half_btf(kCosPi[60], bf0[8], kCosPi[4], bf0[15]);
  bf1[9] = half_btf(kCosPi[28], bf0[9], kCosPi[36], bf0[14]);
  bf1[10] = half_btf(kCosPi[44], bf0[10], kCosPi[20], bf0[13]);
  bf1[11] = half_btf(kCosPi[12], bf0[11], kCosPi[52], bf0[12]);
  bf1[12] = half_btf(kCosPi[12], bf0[12], -kCosPi[52], bf0[11]);
  bf1[13] = half_btf(kCosPi[44], bf0[13], -kCosPi[20], bf0[10]);
  bf1[14] = half_btf(kCosPi[28], bf0[14], -kCosPi[36], bf0[9]);
  bf1[15] = half_btf(kCosPi[60], bf0[15], -kCosPi[4], bf0[8]);

  // stage 7
  bf0 = step;
  bf1 = output;
  constexpr uint32_t shift = 1;
  bf1[0] = div_sqrt2(bf0[0], shift);
  bf1[1] = div_sqrt2(bf0[8], shift);
  bf1[2] = div_sqrt2(bf0[4], shift);
  bf1[3] = div_sqrt2(bf0[12], shift);
  bf1[4] = div_sqrt2(bf0[2], shift);
  bf1[5] = div_sqrt2(bf0[10], shift);
  bf1[6] = div_sqrt2(bf0[6], shift);
  bf1[7] = div_sqrt2(bf0[14], shift);
  bf1[8] = div_sqrt2(bf0[1], shift);
  bf1[9] = div_sqrt2(bf0[9], shift);
  bf1[10] = div_sqrt2(bf0[5], shift);
  bf1[11] = div_sqrt2(bf0[13], shift);
  bf1[12] = div_sqrt2(bf0[3], shift);
  bf1[13] = div_sqrt2(bf0[11], shift);
  bf1[14] = div_sqrt2(bf0[7], shift);
  bf1[15] = div_sqrt2(bf0[15], shift);
}

static void fdct32_C(const int32_t* input, int32_t* output) {
  constexpr uint32_t shift = kPrecision + 2;
  int32_t *bf0, *bf1;
  int32_t step[32];

  // stage 1;
  bf1 = output;
  bf1[0] = input[0] + input[31];
  bf1[1] = input[1] + input[30];
  bf1[2] = input[2] + input[29];
  bf1[3] = input[3] + input[28];
  bf1[4] = input[4] + input[27];
  bf1[5] = input[5] + input[26];
  bf1[6] = input[6] + input[25];
  bf1[7] = input[7] + input[24];
  bf1[8] = input[8] + input[23];
  bf1[9] = input[9] + input[22];
  bf1[10] = input[10] + input[21];
  bf1[11] = input[11] + input[20];
  bf1[12] = input[12] + input[19];
  bf1[13] = input[13] + input[18];
  bf1[14] = input[14] + input[17];
  bf1[15] = input[15] + input[16];
  bf1[16] = -input[16] + input[15];
  bf1[17] = -input[17] + input[14];
  bf1[18] = -input[18] + input[13];
  bf1[19] = -input[19] + input[12];
  bf1[20] = -input[20] + input[11];
  bf1[21] = -input[21] + input[10];
  bf1[22] = -input[22] + input[9];
  bf1[23] = -input[23] + input[8];
  bf1[24] = -input[24] + input[7];
  bf1[25] = -input[25] + input[6];
  bf1[26] = -input[26] + input[5];
  bf1[27] = -input[27] + input[4];
  bf1[28] = -input[28] + input[3];
  bf1[29] = -input[29] + input[2];
  bf1[30] = -input[30] + input[1];
  bf1[31] = -input[31] + input[0];

  // stage 2
  bf0 = output;
  bf1 = step;
  bf1[0] = bf0[0] + bf0[15];
  bf1[1] = bf0[1] + bf0[14];
  bf1[2] = bf0[2] + bf0[13];
  bf1[3] = bf0[3] + bf0[12];
  bf1[4] = bf0[4] + bf0[11];
  bf1[5] = bf0[5] + bf0[10];
  bf1[6] = bf0[6] + bf0[9];
  bf1[7] = bf0[7] + bf0[8];
  bf1[8] = -bf0[8] + bf0[7];
  bf1[9] = -bf0[9] + bf0[6];
  bf1[10] = -bf0[10] + bf0[5];
  bf1[11] = -bf0[11] + bf0[4];
  bf1[12] = -bf0[12] + bf0[3];
  bf1[13] = -bf0[13] + bf0[2];
  bf1[14] = -bf0[14] + bf0[1];
  bf1[15] = -bf0[15] + bf0[0];
  bf1[16] = bf0[16];
  bf1[17] = bf0[17];
  bf1[18] = bf0[18];
  bf1[19] = bf0[19];
  bf1[20] = half_btf(-kCosPi[32], bf0[20], kCosPi[32], bf0[27]);
  bf1[21] = half_btf(-kCosPi[32], bf0[21], kCosPi[32], bf0[26]);
  bf1[22] = half_btf(-kCosPi[32], bf0[22], kCosPi[32], bf0[25]);
  bf1[23] = half_btf(-kCosPi[32], bf0[23], kCosPi[32], bf0[24]);
  bf1[24] = half_btf(kCosPi[32], bf0[24], kCosPi[32], bf0[23]);
  bf1[25] = half_btf(kCosPi[32], bf0[25], kCosPi[32], bf0[22]);
  bf1[26] = half_btf(kCosPi[32], bf0[26], kCosPi[32], bf0[21]);
  bf1[27] = half_btf(kCosPi[32], bf0[27], kCosPi[32], bf0[20]);
  bf1[28] = bf0[28];
  bf1[29] = bf0[29];
  bf1[30] = bf0[30];
  bf1[31] = bf0[31];

  // stage 3
  bf0 = step;
  bf1 = output;
  bf1[0] = bf0[0] + bf0[7];
  bf1[1] = bf0[1] + bf0[6];
  bf1[2] = bf0[2] + bf0[5];
  bf1[3] = bf0[3] + bf0[4];
  bf1[4] = -bf0[4] + bf0[3];
  bf1[5] = -bf0[5] + bf0[2];
  bf1[6] = -bf0[6] + bf0[1];
  bf1[7] = -bf0[7] + bf0[0];
  bf1[8] = bf0[8];
  bf1[9] = bf0[9];
  bf1[10] = half_btf(-kCosPi[32], bf0[10], kCosPi[32], bf0[13]);
  bf1[11] = half_btf(-kCosPi[32], bf0[11], kCosPi[32], bf0[12]);
  bf1[12] = half_btf(kCosPi[32], bf0[12], kCosPi[32], bf0[11]);
  bf1[13] = half_btf(kCosPi[32], bf0[13], kCosPi[32], bf0[10]);
  bf1[14] = bf0[14];
  bf1[15] = bf0[15];
  bf1[16] = bf0[16] + bf0[23];
  bf1[17] = bf0[17] + bf0[22];
  bf1[18] = bf0[18] + bf0[21];
  bf1[19] = bf0[19] + bf0[20];
  bf1[20] = -bf0[20] + bf0[19];
  bf1[21] = -bf0[21] + bf0[18];
  bf1[22] = -bf0[22] + bf0[17];
  bf1[23] = -bf0[23] + bf0[16];
  bf1[24] = -bf0[24] + bf0[31];
  bf1[25] = -bf0[25] + bf0[30];
  bf1[26] = -bf0[26] + bf0[29];
  bf1[27] = -bf0[27] + bf0[28];
  bf1[28] = bf0[28] + bf0[27];
  bf1[29] = bf0[29] + bf0[26];
  bf1[30] = bf0[30] + bf0[25];
  bf1[31] = bf0[31] + bf0[24];

  // stage 4
  bf0 = output;
  bf1 = step;
  bf1[0] = bf0[0] + bf0[3];
  bf1[1] = bf0[1] + bf0[2];
  bf1[2] = -bf0[2] + bf0[1];
  bf1[3] = -bf0[3] + bf0[0];
  bf1[4] = bf0[4];
  bf1[5] = half_btf(-kCosPi[32], bf0[5], kCosPi[32], bf0[6]);
  bf1[6] = half_btf(kCosPi[32], bf0[6], kCosPi[32], bf0[5]);
  bf1[7] = bf0[7];
  bf1[8] = bf0[8] + bf0[11];
  bf1[9] = bf0[9] + bf0[10];
  bf1[10] = -bf0[10] + bf0[9];
  bf1[11] = -bf0[11] + bf0[8];
  bf1[12] = -bf0[12] + bf0[15];
  bf1[13] = -bf0[13] + bf0[14];
  bf1[14] = bf0[14] + bf0[13];
  bf1[15] = bf0[15] + bf0[12];
  bf1[16] = bf0[16];
  bf1[17] = bf0[17];
  bf1[18] = half_btf(-kCosPi[16], bf0[18], kCosPi[48], bf0[29]);
  bf1[19] = half_btf(-kCosPi[16], bf0[19], kCosPi[48], bf0[28]);
  bf1[20] = half_btf(-kCosPi[48], bf0[20], -kCosPi[16], bf0[27]);
  bf1[21] = half_btf(-kCosPi[48], bf0[21], -kCosPi[16], bf0[26]);
  bf1[22] = bf0[22];
  bf1[23] = bf0[23];
  bf1[24] = bf0[24];
  bf1[25] = bf0[25];
  bf1[26] = half_btf(kCosPi[48], bf0[26], -kCosPi[16], bf0[21]);
  bf1[27] = half_btf(kCosPi[48], bf0[27], -kCosPi[16], bf0[20]);
  bf1[28] = half_btf(kCosPi[16], bf0[28], kCosPi[48], bf0[19]);
  bf1[29] = half_btf(kCosPi[16], bf0[29], kCosPi[48], bf0[18]);
  bf1[30] = bf0[30];
  bf1[31] = bf0[31];

  // stage 5
  bf0 = step;
  bf1 = output;
  bf1[0] = half_btf(kCosPi[32], bf0[0], kCosPi[32], bf0[1], shift);
  bf1[1] = half_btf(-kCosPi[32], bf0[1], kCosPi[32], bf0[0], shift);
  bf1[2] = half_btf(kCosPi[48], bf0[2], kCosPi[16], bf0[3], shift);
  bf1[3] = half_btf(kCosPi[48], bf0[3], -kCosPi[16], bf0[2], shift);
  bf1[4] = bf0[4] + bf0[5];
  bf1[5] = -bf0[5] + bf0[4];
  bf1[6] = -bf0[6] + bf0[7];
  bf1[7] = bf0[7] + bf0[6];
  bf1[8] = bf0[8];
  bf1[9] = half_btf(-kCosPi[16], bf0[9], kCosPi[48], bf0[14]);
  bf1[10] = half_btf(-kCosPi[48], bf0[10], -kCosPi[16], bf0[13]);
  bf1[11] = bf0[11];
  bf1[12] = bf0[12];
  bf1[13] = half_btf(kCosPi[48], bf0[13], -kCosPi[16], bf0[10]);
  bf1[14] = half_btf(kCosPi[16], bf0[14], kCosPi[48], bf0[9]);
  bf1[15] = bf0[15];
  bf1[16] = bf0[16] + bf0[19];
  bf1[17] = bf0[17] + bf0[18];
  bf1[18] = -bf0[18] + bf0[17];
  bf1[19] = -bf0[19] + bf0[16];
  bf1[20] = -bf0[20] + bf0[23];
  bf1[21] = -bf0[21] + bf0[22];
  bf1[22] = bf0[22] + bf0[21];
  bf1[23] = bf0[23] + bf0[20];
  bf1[24] = bf0[24] + bf0[27];
  bf1[25] = bf0[25] + bf0[26];
  bf1[26] = -bf0[26] + bf0[25];
  bf1[27] = -bf0[27] + bf0[24];
  bf1[28] = -bf0[28] + bf0[31];
  bf1[29] = -bf0[29] + bf0[30];
  bf1[30] = bf0[30] + bf0[29];
  bf1[31] = bf0[31] + bf0[28];

  // stage 6
  bf0 = output;
  bf1 = step;
  bf1[0] = bf0[0];
  bf1[1] = bf0[1];
  bf1[2] = bf0[2];
  bf1[3] = bf0[3];
  bf1[4] = half_btf(kCosPi[56], bf0[4], kCosPi[8], bf0[7], shift);
  bf1[5] = half_btf(kCosPi[24], bf0[5], kCosPi[40], bf0[6], shift);
  bf1[6] = half_btf(kCosPi[24], bf0[6], -kCosPi[40], bf0[5], shift);
  bf1[7] = half_btf(kCosPi[56], bf0[7], -kCosPi[8], bf0[4], shift);
  bf1[8] = bf0[8] + bf0[9];
  bf1[9] = -bf0[9] + bf0[8];
  bf1[10] = -bf0[10] + bf0[11];
  bf1[11] = bf0[11] + bf0[10];
  bf1[12] = bf0[12] + bf0[13];
  bf1[13] = -bf0[13] + bf0[12];
  bf1[14] = -bf0[14] + bf0[15];
  bf1[15] = bf0[15] + bf0[14];
  bf1[16] = bf0[16];
  bf1[17] = half_btf(-kCosPi[8], bf0[17], kCosPi[56], bf0[30]);
  bf1[18] = half_btf(-kCosPi[56], bf0[18], -kCosPi[8], bf0[29]);
  bf1[19] = bf0[19];
  bf1[20] = bf0[20];
  bf1[21] = half_btf(-kCosPi[40], bf0[21], kCosPi[24], bf0[26]);
  bf1[22] = half_btf(-kCosPi[24], bf0[22], -kCosPi[40], bf0[25]);
  bf1[23] = bf0[23];
  bf1[24] = bf0[24];
  bf1[25] = half_btf(kCosPi[24], bf0[25], -kCosPi[40], bf0[22]);
  bf1[26] = half_btf(kCosPi[40], bf0[26], kCosPi[24], bf0[21]);
  bf1[27] = bf0[27];
  bf1[28] = bf0[28];
  bf1[29] = half_btf(kCosPi[56], bf0[29], -kCosPi[8], bf0[18]);
  bf1[30] = half_btf(kCosPi[8], bf0[30], kCosPi[56], bf0[17]);
  bf1[31] = bf0[31];

  // stage 7
  bf0 = step;
  bf1 = output;
  bf1[0] = bf0[0];
  bf1[1] = bf0[1];
  bf1[2] = bf0[2];
  bf1[3] = bf0[3];
  bf1[4] = bf0[4];
  bf1[5] = bf0[5];
  bf1[6] = bf0[6];
  bf1[7] = bf0[7];
  bf1[8] = half_btf(kCosPi[60], bf0[8], kCosPi[4], bf0[15], shift);
  bf1[9] = half_btf(kCosPi[28], bf0[9], kCosPi[36], bf0[14], shift);
  bf1[10] = half_btf(kCosPi[44], bf0[10], kCosPi[20], bf0[13], shift);
  bf1[11] = half_btf(kCosPi[12], bf0[11], kCosPi[52], bf0[12], shift);
  bf1[12] = half_btf(kCosPi[12], bf0[12], -kCosPi[52], bf0[11], shift);
  bf1[13] = half_btf(kCosPi[44], bf0[13], -kCosPi[20], bf0[10], shift);
  bf1[14] = half_btf(kCosPi[28], bf0[14], -kCosPi[36], bf0[9], shift);
  bf1[15] = half_btf(kCosPi[60], bf0[15], -kCosPi[4], bf0[8], shift);
  bf1[16] = bf0[16] + bf0[17];
  bf1[17] = -bf0[17] + bf0[16];
  bf1[18] = -bf0[18] + bf0[19];
  bf1[19] = bf0[19] + bf0[18];
  bf1[20] = bf0[20] + bf0[21];
  bf1[21] = -bf0[21] + bf0[20];
  bf1[22] = -bf0[22] + bf0[23];
  bf1[23] = bf0[23] + bf0[22];
  bf1[24] = bf0[24] + bf0[25];
  bf1[25] = -bf0[25] + bf0[24];
  bf1[26] = -bf0[26] + bf0[27];
  bf1[27] = bf0[27] + bf0[26];
  bf1[28] = bf0[28] + bf0[29];
  bf1[29] = -bf0[29] + bf0[28];
  bf1[30] = -bf0[30] + bf0[31];
  bf1[31] = bf0[31] + bf0[30];

  // stage 8
  bf0 = output;
  bf1 = step;
  bf1[0] = bf0[0];
  bf1[1] = bf0[1];
  bf1[2] = bf0[2];
  bf1[3] = bf0[3];
  bf1[4] = bf0[4];
  bf1[5] = bf0[5];
  bf1[6] = bf0[6];
  bf1[7] = bf0[7];
  bf1[8] = bf0[8];
  bf1[9] = bf0[9];
  bf1[10] = bf0[10];
  bf1[11] = bf0[11];
  bf1[12] = bf0[12];
  bf1[13] = bf0[13];
  bf1[14] = bf0[14];
  bf1[15] = bf0[15];
  bf1[16] = half_btf(kCosPi[62], bf0[16], kCosPi[2], bf0[31], shift);
  bf1[17] = half_btf(kCosPi[30], bf0[17], kCosPi[34], bf0[30], shift);
  bf1[18] = half_btf(kCosPi[46], bf0[18], kCosPi[18], bf0[29], shift);
  bf1[19] = half_btf(kCosPi[14], bf0[19], kCosPi[50], bf0[28], shift);
  bf1[20] = half_btf(kCosPi[54], bf0[20], kCosPi[10], bf0[27], shift);
  bf1[21] = half_btf(kCosPi[22], bf0[21], kCosPi[42], bf0[26], shift);
  bf1[22] = half_btf(kCosPi[38], bf0[22], kCosPi[26], bf0[25], shift);
  bf1[23] = half_btf(kCosPi[6], bf0[23], kCosPi[58], bf0[24], shift);
  bf1[24] = half_btf(kCosPi[6], bf0[24], -kCosPi[58], bf0[23], shift);
  bf1[25] = half_btf(kCosPi[38], bf0[25], -kCosPi[26], bf0[22], shift);
  bf1[26] = half_btf(kCosPi[22], bf0[26], -kCosPi[42], bf0[21], shift);
  bf1[27] = half_btf(kCosPi[54], bf0[27], -kCosPi[10], bf0[20], shift);
  bf1[28] = half_btf(kCosPi[14], bf0[28], -kCosPi[50], bf0[19], shift);
  bf1[29] = half_btf(kCosPi[46], bf0[29], -kCosPi[18], bf0[18], shift);
  bf1[30] = half_btf(kCosPi[30], bf0[30], -kCosPi[34], bf0[17], shift);
  bf1[31] = half_btf(kCosPi[62], bf0[31], -kCosPi[2], bf0[16], shift);

  // stage 9
  bf0 = step;
  bf1 = output;
  bf1[0] = bf0[0];
  bf1[1] = bf0[16];
  bf1[2] = bf0[8];
  bf1[3] = bf0[24];
  bf1[4] = bf0[4];
  bf1[5] = bf0[20];
  bf1[6] = bf0[12];
  bf1[7] = bf0[28];
  bf1[8] = bf0[2];
  bf1[9] = bf0[18];
  bf1[10] = bf0[10];
  bf1[11] = bf0[26];
  bf1[12] = bf0[6];
  bf1[13] = bf0[22];
  bf1[14] = bf0[14];
  bf1[15] = bf0[30];
  bf1[16] = bf0[1];
  bf1[17] = bf0[17];
  bf1[18] = bf0[9];
  bf1[19] = bf0[25];
  bf1[20] = bf0[5];
  bf1[21] = bf0[21];
  bf1[22] = bf0[13];
  bf1[23] = bf0[29];
  bf1[24] = bf0[3];
  bf1[25] = bf0[19];
  bf1[26] = bf0[11];
  bf1[27] = bf0[27];
  bf1[28] = bf0[7];
  bf1[29] = bf0[23];
  bf1[30] = bf0[15];
  bf1[31] = bf0[31];
}

static void fadst4_C(const int32_t* input, int32_t* output) {
  int64_t x0, x1, x2, x3;
  int64_t s0, s1, s2, s3, s4, s5, s6, s7;

  // stage 0
  x0 = input[0];
  x1 = input[1];
  x2 = input[2];
  x3 = input[3];

  if (!(x0 | x1 | x2 | x3)) {
    output[0] = output[1] = output[2] = output[3] = 0;
    return;
  }

  // stage 1
  s0 = (int64_t)kSinPi[1] * x0;
  s1 = (int64_t)kSinPi[4] * x0;
  s2 = (int64_t)kSinPi[2] * x1;
  s3 = (int64_t)kSinPi[1] * x1;
  s4 = (int64_t)kSinPi[3] * x2;
  s5 = (int64_t)kSinPi[4] * x3;
  s6 = (int64_t)kSinPi[2] * x3;
  s7 = x0 + x1;

  // stage 2
  s7 = s7 - x3;

  // stage 3
  x0 = s0 + s2;
  x1 = (int64_t)kSinPi[3] * s7;
  x2 = s1 - s3;
  x3 = s4;

  // stage 4
  x0 = x0 + s5;
  x2 = x2 + s6;

  // stage 5
  s0 = x0 + x3;
  s1 = x1;
  s2 = x2 - x3;
  s3 = x2 - x0;

  // stage 6
  s3 = s3 + x3;

  output[0] = div_sqrt2(round_shift(s0));
  output[1] = div_sqrt2(round_shift(s1));
  output[2] = div_sqrt2(round_shift(s2));
  output[3] = div_sqrt2(round_shift(s3));
}

static void fadst8_C(const int32_t* input, int32_t* output) {
  const uint32_t shift = kPrecision + 1;
  int32_t *bf0, *bf1;
  int32_t step[8];

  // stage 1;
  assert(output != input);
  bf1 = output;
  bf1[0] = input[0];
  bf1[1] = -input[7];
  bf1[2] = -input[3];
  bf1[3] = input[4];
  bf1[4] = -input[1];
  bf1[5] = input[6];
  bf1[6] = input[2];
  bf1[7] = -input[5];

  // stage 2
  bf0 = output;
  bf1 = step;
  bf1[0] = bf0[0];
  bf1[1] = bf0[1];
  bf1[2] = half_btf(kCosPi[32], bf0[2], kCosPi[32], bf0[3]);
  bf1[3] = half_btf(kCosPi[32], bf0[2], -kCosPi[32], bf0[3]);
  bf1[4] = bf0[4];
  bf1[5] = bf0[5];
  bf1[6] = half_btf(kCosPi[32], bf0[6], kCosPi[32], bf0[7]);
  bf1[7] = half_btf(kCosPi[32], bf0[6], -kCosPi[32], bf0[7]);

  // stage 3
  bf0 = step;
  bf1 = output;
  bf1[0] = bf0[0] + bf0[2];
  bf1[1] = bf0[1] + bf0[3];
  bf1[2] = bf0[0] - bf0[2];
  bf1[3] = bf0[1] - bf0[3];
  bf1[4] = bf0[4] + bf0[6];
  bf1[5] = bf0[5] + bf0[7];
  bf1[6] = bf0[4] - bf0[6];
  bf1[7] = bf0[5] - bf0[7];

  // stage 4
  bf0 = output;
  bf1 = step;
  bf1[0] = bf0[0];
  bf1[1] = bf0[1];
  bf1[2] = bf0[2];
  bf1[3] = bf0[3];
  bf1[4] = half_btf(kCosPi[16], bf0[4], kCosPi[48], bf0[5]);
  bf1[5] = half_btf(kCosPi[48], bf0[4], -kCosPi[16], bf0[5]);
  bf1[6] = half_btf(-kCosPi[48], bf0[6], kCosPi[16], bf0[7]);
  bf1[7] = half_btf(kCosPi[16], bf0[6], kCosPi[48], bf0[7]);

  // stage 5
  bf0 = step;
  bf1 = output;
  bf1[0] = bf0[0] + bf0[4];
  bf1[1] = bf0[1] + bf0[5];
  bf1[2] = bf0[2] + bf0[6];
  bf1[3] = bf0[3] + bf0[7];
  bf1[4] = bf0[0] - bf0[4];
  bf1[5] = bf0[1] - bf0[5];
  bf1[6] = bf0[2] - bf0[6];
  bf1[7] = bf0[3] - bf0[7];

  // stage 6
  bf0 = output;
  bf1 = step;
  bf1[0] = half_btf(kCosPi[4], bf0[0], kCosPi[60], bf0[1], shift);
  bf1[1] = half_btf(kCosPi[60], bf0[0], -kCosPi[4], bf0[1], shift);
  bf1[2] = half_btf(kCosPi[20], bf0[2], kCosPi[44], bf0[3], shift);
  bf1[3] = half_btf(kCosPi[44], bf0[2], -kCosPi[20], bf0[3], shift);
  bf1[4] = half_btf(kCosPi[36], bf0[4], kCosPi[28], bf0[5], shift);
  bf1[5] = half_btf(kCosPi[28], bf0[4], -kCosPi[36], bf0[5], shift);
  bf1[6] = half_btf(kCosPi[52], bf0[6], kCosPi[12], bf0[7], shift);
  bf1[7] = half_btf(kCosPi[12], bf0[6], -kCosPi[52], bf0[7], shift);

  // stage 7
  bf0 = step;
  bf1 = output;
  bf1[0] = bf0[1];
  bf1[1] = bf0[6];
  bf1[2] = bf0[3];
  bf1[3] = bf0[4];
  bf1[4] = bf0[5];
  bf1[5] = bf0[2];
  bf1[6] = bf0[7];
  bf1[7] = bf0[0];
}

static void fadst16_C(const int32_t* input, int32_t* output) {
  const uint32_t shift = kPrecision + 1;
  int32_t *bf0, *bf1;
  int32_t step[16];

  // stage 1;
  assert(output != input);
  bf1 = output;
  bf1[0] = input[0];
  bf1[1] = -input[15];
  bf1[2] = -input[7];
  bf1[3] = input[8];
  bf1[4] = -input[3];
  bf1[5] = input[12];
  bf1[6] = input[4];
  bf1[7] = -input[11];
  bf1[8] = -input[1];
  bf1[9] = input[14];
  bf1[10] = input[6];
  bf1[11] = -input[9];
  bf1[12] = input[2];
  bf1[13] = -input[13];
  bf1[14] = -input[5];
  bf1[15] = input[10];

  // stage 2
  bf0 = output;
  bf1 = step;
  bf1[0] = bf0[0];
  bf1[1] = bf0[1];
  bf1[2] = half_btf(kCosPi[32], bf0[2], kCosPi[32], bf0[3]);
  bf1[3] = half_btf(kCosPi[32], bf0[2], -kCosPi[32], bf0[3]);
  bf1[4] = bf0[4];
  bf1[5] = bf0[5];
  bf1[6] = half_btf(kCosPi[32], bf0[6], kCosPi[32], bf0[7]);
  bf1[7] = half_btf(kCosPi[32], bf0[6], -kCosPi[32], bf0[7]);
  bf1[8] = bf0[8];
  bf1[9] = bf0[9];
  bf1[10] = half_btf(kCosPi[32], bf0[10], kCosPi[32], bf0[11]);
  bf1[11] = half_btf(kCosPi[32], bf0[10], -kCosPi[32], bf0[11]);
  bf1[12] = bf0[12];
  bf1[13] = bf0[13];
  bf1[14] = half_btf(kCosPi[32], bf0[14], kCosPi[32], bf0[15]);
  bf1[15] = half_btf(kCosPi[32], bf0[14], -kCosPi[32], bf0[15]);

  // stage 3
  bf0 = step;
  bf1 = output;
  bf1[0] = bf0[0] + bf0[2];
  bf1[1] = bf0[1] + bf0[3];
  bf1[2] = bf0[0] - bf0[2];
  bf1[3] = bf0[1] - bf0[3];
  bf1[4] = bf0[4] + bf0[6];
  bf1[5] = bf0[5] + bf0[7];
  bf1[6] = bf0[4] - bf0[6];
  bf1[7] = bf0[5] - bf0[7];
  bf1[8] = bf0[8] + bf0[10];
  bf1[9] = bf0[9] + bf0[11];
  bf1[10] = bf0[8] - bf0[10];
  bf1[11] = bf0[9] - bf0[11];
  bf1[12] = bf0[12] + bf0[14];
  bf1[13] = bf0[13] + bf0[15];
  bf1[14] = bf0[12] - bf0[14];
  bf1[15] = bf0[13] - bf0[15];

  // stage 4
  bf0 = output;
  bf1 = step;
  bf1[0] = bf0[0];
  bf1[1] = bf0[1];
  bf1[2] = bf0[2];
  bf1[3] = bf0[3];
  bf1[4] = half_btf(kCosPi[16], bf0[4], kCosPi[48], bf0[5]);
  bf1[5] = half_btf(kCosPi[48], bf0[4], -kCosPi[16], bf0[5]);
  bf1[6] = half_btf(-kCosPi[48], bf0[6], kCosPi[16], bf0[7]);
  bf1[7] = half_btf(kCosPi[16], bf0[6], kCosPi[48], bf0[7]);
  bf1[8] = bf0[8];
  bf1[9] = bf0[9];
  bf1[10] = bf0[10];
  bf1[11] = bf0[11];
  bf1[12] = half_btf(kCosPi[16], bf0[12], kCosPi[48], bf0[13]);
  bf1[13] = half_btf(kCosPi[48], bf0[12], -kCosPi[16], bf0[13]);
  bf1[14] = half_btf(-kCosPi[48], bf0[14], kCosPi[16], bf0[15]);
  bf1[15] = half_btf(kCosPi[16], bf0[14], kCosPi[48], bf0[15]);

  // stage 5
  bf0 = step;
  bf1 = output;
  bf1[0] = bf0[0] + bf0[4];
  bf1[1] = bf0[1] + bf0[5];
  bf1[2] = bf0[2] + bf0[6];
  bf1[3] = bf0[3] + bf0[7];
  bf1[4] = bf0[0] - bf0[4];
  bf1[5] = bf0[1] - bf0[5];
  bf1[6] = bf0[2] - bf0[6];
  bf1[7] = bf0[3] - bf0[7];
  bf1[8] = bf0[8] + bf0[12];
  bf1[9] = bf0[9] + bf0[13];
  bf1[10] = bf0[10] + bf0[14];
  bf1[11] = bf0[11] + bf0[15];
  bf1[12] = bf0[8] - bf0[12];
  bf1[13] = bf0[9] - bf0[13];
  bf1[14] = bf0[10] - bf0[14];
  bf1[15] = bf0[11] - bf0[15];

  // stage 6
  bf0 = output;
  bf1 = step;
  bf1[0] = bf0[0];
  bf1[1] = bf0[1];
  bf1[2] = bf0[2];
  bf1[3] = bf0[3];
  bf1[4] = bf0[4];
  bf1[5] = bf0[5];
  bf1[6] = bf0[6];
  bf1[7] = bf0[7];
  bf1[8] = half_btf(kCosPi[8], bf0[8], kCosPi[56], bf0[9]);
  bf1[9] = half_btf(kCosPi[56], bf0[8], -kCosPi[8], bf0[9]);
  bf1[10] = half_btf(kCosPi[40], bf0[10], kCosPi[24], bf0[11]);
  bf1[11] = half_btf(kCosPi[24], bf0[10], -kCosPi[40], bf0[11]);
  bf1[12] = half_btf(-kCosPi[56], bf0[12], kCosPi[8], bf0[13]);
  bf1[13] = half_btf(kCosPi[8], bf0[12], kCosPi[56], bf0[13]);
  bf1[14] = half_btf(-kCosPi[24], bf0[14], kCosPi[40], bf0[15]);
  bf1[15] = half_btf(kCosPi[40], bf0[14], kCosPi[24], bf0[15]);

  // stage 7
  bf0 = step;
  bf1 = output;
  bf1[0] = bf0[0] + bf0[8];
  bf1[1] = bf0[1] + bf0[9];
  bf1[2] = bf0[2] + bf0[10];
  bf1[3] = bf0[3] + bf0[11];
  bf1[4] = bf0[4] + bf0[12];
  bf1[5] = bf0[5] + bf0[13];
  bf1[6] = bf0[6] + bf0[14];
  bf1[7] = bf0[7] + bf0[15];
  bf1[8] = bf0[0] - bf0[8];
  bf1[9] = bf0[1] - bf0[9];
  bf1[10] = bf0[2] - bf0[10];
  bf1[11] = bf0[3] - bf0[11];
  bf1[12] = bf0[4] - bf0[12];
  bf1[13] = bf0[5] - bf0[13];
  bf1[14] = bf0[6] - bf0[14];
  bf1[15] = bf0[7] - bf0[15];

  // stage 8
  bf0 = output;
  bf1 = step;
  bf1[0] = half_btf(kCosPi[2], bf0[0], kCosPi[62], bf0[1]);
  bf1[1] = half_btf(kCosPi[62], bf0[0], -kCosPi[2], bf0[1]);
  bf1[2] = half_btf(kCosPi[10], bf0[2], kCosPi[54], bf0[3]);
  bf1[3] = half_btf(kCosPi[54], bf0[2], -kCosPi[10], bf0[3]);
  bf1[4] = half_btf(kCosPi[18], bf0[4], kCosPi[46], bf0[5]);
  bf1[5] = half_btf(kCosPi[46], bf0[4], -kCosPi[18], bf0[5]);
  bf1[6] = half_btf(kCosPi[26], bf0[6], kCosPi[38], bf0[7]);
  bf1[7] = half_btf(kCosPi[38], bf0[6], -kCosPi[26], bf0[7]);
  bf1[8] = half_btf(kCosPi[34], bf0[8], kCosPi[30], bf0[9]);
  bf1[9] = half_btf(kCosPi[30], bf0[8], -kCosPi[34], bf0[9]);
  bf1[10] = half_btf(kCosPi[42], bf0[10], kCosPi[22], bf0[11]);
  bf1[11] = half_btf(kCosPi[22], bf0[10], -kCosPi[42], bf0[11]);
  bf1[12] = half_btf(kCosPi[50], bf0[12], kCosPi[14], bf0[13]);
  bf1[13] = half_btf(kCosPi[14], bf0[12], -kCosPi[50], bf0[13]);
  bf1[14] = half_btf(kCosPi[58], bf0[14], kCosPi[6], bf0[15]);
  bf1[15] = half_btf(kCosPi[6], bf0[14], -kCosPi[58], bf0[15]);

  // stage 9
  bf0 = step;
  bf1 = output;
  bf1[0] = div_sqrt2(bf0[1], shift - kPrecision);
  bf1[1] = div_sqrt2(bf0[14], shift - kPrecision);
  bf1[2] = div_sqrt2(bf0[3], shift - kPrecision);
  bf1[3] = div_sqrt2(bf0[12], shift - kPrecision);
  bf1[4] = div_sqrt2(bf0[5], shift - kPrecision);
  bf1[5] = div_sqrt2(bf0[10], shift - kPrecision);
  bf1[6] = div_sqrt2(bf0[7], shift - kPrecision);
  bf1[7] = div_sqrt2(bf0[8], shift - kPrecision);
  bf1[8] = div_sqrt2(bf0[9], shift - kPrecision);
  bf1[9] = div_sqrt2(bf0[6], shift - kPrecision);
  bf1[10] = div_sqrt2(bf0[11], shift - kPrecision);
  bf1[11] = div_sqrt2(bf0[4], shift - kPrecision);
  bf1[12] = div_sqrt2(bf0[13], shift - kPrecision);
  bf1[13] = div_sqrt2(bf0[2], shift - kPrecision);
  bf1[14] = div_sqrt2(bf0[15], shift - kPrecision);
  bf1[15] = div_sqrt2(bf0[0], shift - kPrecision);
}

static void fadst32_C(const int32_t* input, int32_t* output) {
  constexpr uint32_t shift = kPrecision + 2;
  int32_t *bf0, *bf1;
  int32_t step[32];

  // stage 1;
  bf1 = output;
  bf1[0] = input[31];
  bf1[1] = input[0];
  bf1[2] = input[29];
  bf1[3] = input[2];
  bf1[4] = input[27];
  bf1[5] = input[4];
  bf1[6] = input[25];
  bf1[7] = input[6];
  bf1[8] = input[23];
  bf1[9] = input[8];
  bf1[10] = input[21];
  bf1[11] = input[10];
  bf1[12] = input[19];
  bf1[13] = input[12];
  bf1[14] = input[17];
  bf1[15] = input[14];
  bf1[16] = input[15];
  bf1[17] = input[16];
  bf1[18] = input[13];
  bf1[19] = input[18];
  bf1[20] = input[11];
  bf1[21] = input[20];
  bf1[22] = input[9];
  bf1[23] = input[22];
  bf1[24] = input[7];
  bf1[25] = input[24];
  bf1[26] = input[5];
  bf1[27] = input[26];
  bf1[28] = input[3];
  bf1[29] = input[28];
  bf1[30] = input[1];
  bf1[31] = input[30];

  // stage 2
  bf0 = output;
  bf1 = step;
  bf1[0] = half_btf(kCosPi[1], bf0[0], kCosPi[63], bf0[1]);
  bf1[1] = half_btf(-kCosPi[1], bf0[1], kCosPi[63], bf0[0]);
  bf1[2] = half_btf(kCosPi[5], bf0[2], kCosPi[59], bf0[3]);
  bf1[3] = half_btf(-kCosPi[5], bf0[3], kCosPi[59], bf0[2]);
  bf1[4] = half_btf(kCosPi[9], bf0[4], kCosPi[55], bf0[5]);
  bf1[5] = half_btf(-kCosPi[9], bf0[5], kCosPi[55], bf0[4]);
  bf1[6] = half_btf(kCosPi[13], bf0[6], kCosPi[51], bf0[7]);
  bf1[7] = half_btf(-kCosPi[13], bf0[7], kCosPi[51], bf0[6]);
  bf1[8] = half_btf(kCosPi[17], bf0[8], kCosPi[47], bf0[9]);
  bf1[9] = half_btf(-kCosPi[17], bf0[9], kCosPi[47], bf0[8]);
  bf1[10] = half_btf(kCosPi[21], bf0[10], kCosPi[43], bf0[11]);
  bf1[11] = half_btf(-kCosPi[21], bf0[11], kCosPi[43], bf0[10]);
  bf1[12] = half_btf(kCosPi[25], bf0[12], kCosPi[39], bf0[13]);
  bf1[13] = half_btf(-kCosPi[25], bf0[13], kCosPi[39], bf0[12]);
  bf1[14] = half_btf(kCosPi[29], bf0[14], kCosPi[35], bf0[15]);
  bf1[15] = half_btf(-kCosPi[29], bf0[15], kCosPi[35], bf0[14]);
  bf1[16] = half_btf(kCosPi[33], bf0[16], kCosPi[31], bf0[17]);
  bf1[17] = half_btf(-kCosPi[33], bf0[17], kCosPi[31], bf0[16]);
  bf1[18] = half_btf(kCosPi[37], bf0[18], kCosPi[27], bf0[19]);
  bf1[19] = half_btf(-kCosPi[37], bf0[19], kCosPi[27], bf0[18]);
  bf1[20] = half_btf(kCosPi[41], bf0[20], kCosPi[23], bf0[21]);
  bf1[21] = half_btf(-kCosPi[41], bf0[21], kCosPi[23], bf0[20]);
  bf1[22] = half_btf(kCosPi[45], bf0[22], kCosPi[19], bf0[23]);
  bf1[23] = half_btf(-kCosPi[45], bf0[23], kCosPi[19], bf0[22]);
  bf1[24] = half_btf(kCosPi[49], bf0[24], kCosPi[15], bf0[25]);
  bf1[25] = half_btf(-kCosPi[49], bf0[25], kCosPi[15], bf0[24]);
  bf1[26] = half_btf(kCosPi[53], bf0[26], kCosPi[11], bf0[27]);
  bf1[27] = half_btf(-kCosPi[53], bf0[27], kCosPi[11], bf0[26]);
  bf1[28] = half_btf(kCosPi[57], bf0[28], kCosPi[7], bf0[29]);
  bf1[29] = half_btf(-kCosPi[57], bf0[29], kCosPi[7], bf0[28]);
  bf1[30] = half_btf(kCosPi[61], bf0[30], kCosPi[3], bf0[31]);
  bf1[31] = half_btf(-kCosPi[61], bf0[31], kCosPi[3], bf0[30]);

  // stage 3
  bf0 = step;
  bf1 = output;
  bf1[0] = bf0[0] + bf0[16];
  bf1[1] = bf0[1] + bf0[17];
  bf1[2] = bf0[2] + bf0[18];
  bf1[3] = bf0[3] + bf0[19];
  bf1[4] = bf0[4] + bf0[20];
  bf1[5] = bf0[5] + bf0[21];
  bf1[6] = bf0[6] + bf0[22];
  bf1[7] = bf0[7] + bf0[23];
  bf1[8] = bf0[8] + bf0[24];
  bf1[9] = bf0[9] + bf0[25];
  bf1[10] = bf0[10] + bf0[26];
  bf1[11] = bf0[11] + bf0[27];
  bf1[12] = bf0[12] + bf0[28];
  bf1[13] = bf0[13] + bf0[29];
  bf1[14] = bf0[14] + bf0[30];
  bf1[15] = bf0[15] + bf0[31];
  bf1[16] = -bf0[16] + bf0[0];
  bf1[17] = -bf0[17] + bf0[1];
  bf1[18] = -bf0[18] + bf0[2];
  bf1[19] = -bf0[19] + bf0[3];
  bf1[20] = -bf0[20] + bf0[4];
  bf1[21] = -bf0[21] + bf0[5];
  bf1[22] = -bf0[22] + bf0[6];
  bf1[23] = -bf0[23] + bf0[7];
  bf1[24] = -bf0[24] + bf0[8];
  bf1[25] = -bf0[25] + bf0[9];
  bf1[26] = -bf0[26] + bf0[10];
  bf1[27] = -bf0[27] + bf0[11];
  bf1[28] = -bf0[28] + bf0[12];
  bf1[29] = -bf0[29] + bf0[13];
  bf1[30] = -bf0[30] + bf0[14];
  bf1[31] = -bf0[31] + bf0[15];

  // stage 4
  bf0 = output;
  bf1 = step;
  bf1[0] = bf0[0];
  bf1[1] = bf0[1];
  bf1[2] = bf0[2];
  bf1[3] = bf0[3];
  bf1[4] = bf0[4];
  bf1[5] = bf0[5];
  bf1[6] = bf0[6];
  bf1[7] = bf0[7];
  bf1[8] = bf0[8];
  bf1[9] = bf0[9];
  bf1[10] = bf0[10];
  bf1[11] = bf0[11];
  bf1[12] = bf0[12];
  bf1[13] = bf0[13];
  bf1[14] = bf0[14];
  bf1[15] = bf0[15];
  bf1[16] = half_btf(kCosPi[4], bf0[16], kCosPi[60], bf0[17]);
  bf1[17] = half_btf(-kCosPi[4], bf0[17], kCosPi[60], bf0[16]);
  bf1[18] = half_btf(kCosPi[20], bf0[18], kCosPi[44], bf0[19]);
  bf1[19] = half_btf(-kCosPi[20], bf0[19], kCosPi[44], bf0[18]);
  bf1[20] = half_btf(kCosPi[36], bf0[20], kCosPi[28], bf0[21]);
  bf1[21] = half_btf(-kCosPi[36], bf0[21], kCosPi[28], bf0[20]);
  bf1[22] = half_btf(kCosPi[52], bf0[22], kCosPi[12], bf0[23]);
  bf1[23] = half_btf(-kCosPi[52], bf0[23], kCosPi[12], bf0[22]);
  bf1[24] = half_btf(-kCosPi[60], bf0[24], kCosPi[4], bf0[25]);
  bf1[25] = half_btf(kCosPi[60], bf0[25], kCosPi[4], bf0[24]);
  bf1[26] = half_btf(-kCosPi[44], bf0[26], kCosPi[20], bf0[27]);
  bf1[27] = half_btf(kCosPi[44], bf0[27], kCosPi[20], bf0[26]);
  bf1[28] = half_btf(-kCosPi[28], bf0[28], kCosPi[36], bf0[29]);
  bf1[29] = half_btf(kCosPi[28], bf0[29], kCosPi[36], bf0[28]);
  bf1[30] = half_btf(-kCosPi[12], bf0[30], kCosPi[52], bf0[31]);
  bf1[31] = half_btf(kCosPi[12], bf0[31], kCosPi[52], bf0[30]);

  // stage 5
  bf0 = step;
  bf1 = output;
  bf1[0] = bf0[0] + bf0[8];
  bf1[1] = bf0[1] + bf0[9];
  bf1[2] = bf0[2] + bf0[10];
  bf1[3] = bf0[3] + bf0[11];
  bf1[4] = bf0[4] + bf0[12];
  bf1[5] = bf0[5] + bf0[13];
  bf1[6] = bf0[6] + bf0[14];
  bf1[7] = bf0[7] + bf0[15];
  bf1[8] = -bf0[8] + bf0[0];
  bf1[9] = -bf0[9] + bf0[1];
  bf1[10] = -bf0[10] + bf0[2];
  bf1[11] = -bf0[11] + bf0[3];
  bf1[12] = -bf0[12] + bf0[4];
  bf1[13] = -bf0[13] + bf0[5];
  bf1[14] = -bf0[14] + bf0[6];
  bf1[15] = -bf0[15] + bf0[7];
  bf1[16] = bf0[16] + bf0[24];
  bf1[17] = bf0[17] + bf0[25];
  bf1[18] = bf0[18] + bf0[26];
  bf1[19] = bf0[19] + bf0[27];
  bf1[20] = bf0[20] + bf0[28];
  bf1[21] = bf0[21] + bf0[29];
  bf1[22] = bf0[22] + bf0[30];
  bf1[23] = bf0[23] + bf0[31];
  bf1[24] = -bf0[24] + bf0[16];
  bf1[25] = -bf0[25] + bf0[17];
  bf1[26] = -bf0[26] + bf0[18];
  bf1[27] = -bf0[27] + bf0[19];
  bf1[28] = -bf0[28] + bf0[20];
  bf1[29] = -bf0[29] + bf0[21];
  bf1[30] = -bf0[30] + bf0[22];
  bf1[31] = -bf0[31] + bf0[23];

  // stage 6
  bf0 = output;
  bf1 = step;
  bf1[0] = bf0[0];
  bf1[1] = bf0[1];
  bf1[2] = bf0[2];
  bf1[3] = bf0[3];
  bf1[4] = bf0[4];
  bf1[5] = bf0[5];
  bf1[6] = bf0[6];
  bf1[7] = bf0[7];
  bf1[8] = half_btf(kCosPi[8], bf0[8], kCosPi[56], bf0[9]);
  bf1[9] = half_btf(-kCosPi[8], bf0[9], kCosPi[56], bf0[8]);
  bf1[10] = half_btf(kCosPi[40], bf0[10], kCosPi[24], bf0[11]);
  bf1[11] = half_btf(-kCosPi[40], bf0[11], kCosPi[24], bf0[10]);
  bf1[12] = half_btf(-kCosPi[56], bf0[12], kCosPi[8], bf0[13]);
  bf1[13] = half_btf(kCosPi[56], bf0[13], kCosPi[8], bf0[12]);
  bf1[14] = half_btf(-kCosPi[24], bf0[14], kCosPi[40], bf0[15]);
  bf1[15] = half_btf(kCosPi[24], bf0[15], kCosPi[40], bf0[14]);
  bf1[16] = bf0[16];
  bf1[17] = bf0[17];
  bf1[18] = bf0[18];
  bf1[19] = bf0[19];
  bf1[20] = bf0[20];
  bf1[21] = bf0[21];
  bf1[22] = bf0[22];
  bf1[23] = bf0[23];
  bf1[24] = half_btf(kCosPi[8], bf0[24], kCosPi[56], bf0[25]);
  bf1[25] = half_btf(-kCosPi[8], bf0[25], kCosPi[56], bf0[24]);
  bf1[26] = half_btf(kCosPi[40], bf0[26], kCosPi[24], bf0[27]);
  bf1[27] = half_btf(-kCosPi[40], bf0[27], kCosPi[24], bf0[26]);
  bf1[28] = half_btf(-kCosPi[56], bf0[28], kCosPi[8], bf0[29]);
  bf1[29] = half_btf(kCosPi[56], bf0[29], kCosPi[8], bf0[28]);
  bf1[30] = half_btf(-kCosPi[24], bf0[30], kCosPi[40], bf0[31]);
  bf1[31] = half_btf(kCosPi[24], bf0[31], kCosPi[40], bf0[30]);

  // stage 7
  bf0 = step;
  bf1 = output;
  bf1[0] = bf0[0] + bf0[4];
  bf1[1] = bf0[1] + bf0[5];
  bf1[2] = bf0[2] + bf0[6];
  bf1[3] = bf0[3] + bf0[7];
  bf1[4] = -bf0[4] + bf0[0];
  bf1[5] = -bf0[5] + bf0[1];
  bf1[6] = -bf0[6] + bf0[2];
  bf1[7] = -bf0[7] + bf0[3];
  bf1[8] = bf0[8] + bf0[12];
  bf1[9] = bf0[9] + bf0[13];
  bf1[10] = bf0[10] + bf0[14];
  bf1[11] = bf0[11] + bf0[15];
  bf1[12] = -bf0[12] + bf0[8];
  bf1[13] = -bf0[13] + bf0[9];
  bf1[14] = -bf0[14] + bf0[10];
  bf1[15] = -bf0[15] + bf0[11];
  bf1[16] = bf0[16] + bf0[20];
  bf1[17] = bf0[17] + bf0[21];
  bf1[18] = bf0[18] + bf0[22];
  bf1[19] = bf0[19] + bf0[23];
  bf1[20] = -bf0[20] + bf0[16];
  bf1[21] = -bf0[21] + bf0[17];
  bf1[22] = -bf0[22] + bf0[18];
  bf1[23] = -bf0[23] + bf0[19];
  bf1[24] = bf0[24] + bf0[28];
  bf1[25] = bf0[25] + bf0[29];
  bf1[26] = bf0[26] + bf0[30];
  bf1[27] = bf0[27] + bf0[31];
  bf1[28] = -bf0[28] + bf0[24];
  bf1[29] = -bf0[29] + bf0[25];
  bf1[30] = -bf0[30] + bf0[26];
  bf1[31] = -bf0[31] + bf0[27];

  // stage 8
  bf0 = output;
  bf1 = step;
  bf1[0] = bf0[0];
  bf1[1] = bf0[1];
  bf1[2] = bf0[2];
  bf1[3] = bf0[3];
  bf1[4] = half_btf(kCosPi[16], bf0[4], kCosPi[48], bf0[5]);
  bf1[5] = half_btf(-kCosPi[16], bf0[5], kCosPi[48], bf0[4]);
  bf1[6] = half_btf(-kCosPi[48], bf0[6], kCosPi[16], bf0[7]);
  bf1[7] = half_btf(kCosPi[48], bf0[7], kCosPi[16], bf0[6]);
  bf1[8] = bf0[8];
  bf1[9] = bf0[9];
  bf1[10] = bf0[10];
  bf1[11] = bf0[11];
  bf1[12] = half_btf(kCosPi[16], bf0[12], kCosPi[48], bf0[13]);
  bf1[13] = half_btf(-kCosPi[16], bf0[13], kCosPi[48], bf0[12]);
  bf1[14] = half_btf(-kCosPi[48], bf0[14], kCosPi[16], bf0[15]);
  bf1[15] = half_btf(kCosPi[48], bf0[15], kCosPi[16], bf0[14]);
  bf1[16] = bf0[16];
  bf1[17] = bf0[17];
  bf1[18] = bf0[18];
  bf1[19] = bf0[19];
  bf1[20] = half_btf(kCosPi[16], bf0[20], kCosPi[48], bf0[21]);
  bf1[21] = half_btf(-kCosPi[16], bf0[21], kCosPi[48], bf0[20]);
  bf1[22] = half_btf(-kCosPi[48], bf0[22], kCosPi[16], bf0[23]);
  bf1[23] = half_btf(kCosPi[48], bf0[23], kCosPi[16], bf0[22]);
  bf1[24] = bf0[24];
  bf1[25] = bf0[25];
  bf1[26] = bf0[26];
  bf1[27] = bf0[27];
  bf1[28] = half_btf(kCosPi[16], bf0[28], kCosPi[48], bf0[29]);
  bf1[29] = half_btf(-kCosPi[16], bf0[29], kCosPi[48], bf0[28]);
  bf1[30] = half_btf(-kCosPi[48], bf0[30], kCosPi[16], bf0[31]);
  bf1[31] = half_btf(kCosPi[48], bf0[31], kCosPi[16], bf0[30]);

  // stage 9
  bf0 = step;
  bf1 = output;
  bf1[0] = bf0[0] + bf0[2];
  bf1[1] = bf0[1] + bf0[3];
  bf1[2] = -bf0[2] + bf0[0];
  bf1[3] = -bf0[3] + bf0[1];
  bf1[4] = bf0[4] + bf0[6];
  bf1[5] = bf0[5] + bf0[7];
  bf1[6] = -bf0[6] + bf0[4];
  bf1[7] = -bf0[7] + bf0[5];
  bf1[8] = bf0[8] + bf0[10];
  bf1[9] = bf0[9] + bf0[11];
  bf1[10] = -bf0[10] + bf0[8];
  bf1[11] = -bf0[11] + bf0[9];
  bf1[12] = bf0[12] + bf0[14];
  bf1[13] = bf0[13] + bf0[15];
  bf1[14] = -bf0[14] + bf0[12];
  bf1[15] = -bf0[15] + bf0[13];
  bf1[16] = bf0[16] + bf0[18];
  bf1[17] = bf0[17] + bf0[19];
  bf1[18] = -bf0[18] + bf0[16];
  bf1[19] = -bf0[19] + bf0[17];
  bf1[20] = bf0[20] + bf0[22];
  bf1[21] = bf0[21] + bf0[23];
  bf1[22] = -bf0[22] + bf0[20];
  bf1[23] = -bf0[23] + bf0[21];
  bf1[24] = bf0[24] + bf0[26];
  bf1[25] = bf0[25] + bf0[27];
  bf1[26] = -bf0[26] + bf0[24];
  bf1[27] = -bf0[27] + bf0[25];
  bf1[28] = bf0[28] + bf0[30];
  bf1[29] = bf0[29] + bf0[31];
  bf1[30] = -bf0[30] + bf0[28];
  bf1[31] = -bf0[31] + bf0[29];

  // stage 10
  bf0 = output;
  bf1 = step;
  bf1[0] = round_shift(bf0[0], shift - kPrecision);
  bf1[1] = round_shift(bf0[1], shift - kPrecision);
  bf1[2] = half_btf(kCosPi[32], bf0[2], kCosPi[32], bf0[3], shift);
  bf1[3] = half_btf(-kCosPi[32], bf0[3], kCosPi[32], bf0[2], shift);
  bf1[4] = round_shift(bf0[4], shift - kPrecision);
  bf1[5] = round_shift(bf0[5], shift - kPrecision);
  bf1[6] = half_btf(kCosPi[32], bf0[6], kCosPi[32], bf0[7], shift);
  bf1[7] = half_btf(-kCosPi[32], bf0[7], kCosPi[32], bf0[6], shift);
  bf1[8] = round_shift(bf0[8], shift - kPrecision);
  bf1[9] = round_shift(bf0[9], shift - kPrecision);
  bf1[10] = half_btf(kCosPi[32], bf0[10], kCosPi[32], bf0[11], shift);
  bf1[11] = half_btf(-kCosPi[32], bf0[11], kCosPi[32], bf0[10], shift);
  bf1[12] = round_shift(bf0[12], shift - kPrecision);
  bf1[13] = round_shift(bf0[13], shift - kPrecision);
  bf1[14] = half_btf(kCosPi[32], bf0[14], kCosPi[32], bf0[15], shift);
  bf1[15] = half_btf(-kCosPi[32], bf0[15], kCosPi[32], bf0[14], shift);
  bf1[16] = round_shift(bf0[16], shift - kPrecision);
  bf1[17] = round_shift(bf0[17], shift - kPrecision);
  bf1[18] = half_btf(kCosPi[32], bf0[18], kCosPi[32], bf0[19], shift);
  bf1[19] = half_btf(-kCosPi[32], bf0[19], kCosPi[32], bf0[18], shift);
  bf1[20] = round_shift(bf0[20], shift - kPrecision);
  bf1[21] = round_shift(bf0[21], shift - kPrecision);
  bf1[22] = half_btf(kCosPi[32], bf0[22], kCosPi[32], bf0[23], shift);
  bf1[23] = half_btf(-kCosPi[32], bf0[23], kCosPi[32], bf0[22], shift);
  bf1[24] = round_shift(bf0[24], shift - kPrecision);
  bf1[25] = round_shift(bf0[25], shift - kPrecision);
  bf1[26] = half_btf(kCosPi[32], bf0[26], kCosPi[32], bf0[27], shift);
  bf1[27] = half_btf(-kCosPi[32], bf0[27], kCosPi[32], bf0[26], shift);
  bf1[28] = round_shift(bf0[28], shift - kPrecision);
  bf1[29] = round_shift(bf0[29], shift - kPrecision);
  bf1[30] = half_btf(kCosPi[32], bf0[30], kCosPi[32], bf0[31], shift);
  bf1[31] = half_btf(-kCosPi[32], bf0[31], kCosPi[32], bf0[30], shift);

  // stage 11
  bf0 = step;
  bf1 = output;
  bf1[0] = bf0[0];
  bf1[1] = -bf0[16];
  bf1[2] = bf0[24];
  bf1[3] = -bf0[8];
  bf1[4] = bf0[12];
  bf1[5] = -bf0[28];
  bf1[6] = bf0[20];
  bf1[7] = -bf0[4];
  bf1[8] = bf0[6];
  bf1[9] = -bf0[22];
  bf1[10] = bf0[30];
  bf1[11] = -bf0[14];
  bf1[12] = bf0[10];
  bf1[13] = -bf0[26];
  bf1[14] = bf0[18];
  bf1[15] = -bf0[2];
  bf1[16] = bf0[3];
  bf1[17] = -bf0[19];
  bf1[18] = bf0[27];
  bf1[19] = -bf0[11];
  bf1[20] = bf0[15];
  bf1[21] = -bf0[31];
  bf1[22] = bf0[23];
  bf1[23] = -bf0[7];
  bf1[24] = bf0[5];
  bf1[25] = -bf0[21];
  bf1[26] = bf0[29];
  bf1[27] = -bf0[13];
  bf1[28] = bf0[9];
  bf1[29] = -bf0[25];
  bf1[30] = bf0[17];
  bf1[31] = -bf0[1];
}

//------------------------------------------------------------------------------

static void Identity(const int32_t* input, int32_t* output, int8_t size,
                     int8_t multiplier) {
  if (multiplier == 1) {
    memcpy(output, input, sizeof(int32_t) * size);
  } else {
    for (int8_t i = 0; i < size; ++i) {
      output[i] = input[i] * multiplier;
    }
  }
}

static void iIdentity(const int32_t* input, int32_t* output, int8_t size,
                     int8_t divider) {
  if (divider == 1) {
    memcpy(output, input, sizeof(int32_t) * size);
  } else {
    for (int8_t i = 0; i < size; ++i) {
      const int32_t v = input[i];
      output[i] = (v < 0) ? -((-v) / divider) : (v / divider);
    }
  }
}

//------------------------------------------------------------------------------

static void Transpose_C(const int32_t* in, uint32_t w, uint32_t h,
                        int32_t* out) {
  assert(in != out);
  for (size_t j = 0; j < h; ++j) {
    for (size_t i = 0; i < w; ++i) {
      out[j + h * i] = in[i + w * j];
    }
  }
}

//------------------------------------------------------------------------------

static int32_t SlowDct8x8_C(const int32_t in[64],
                            const float cos_x[8], const float cos_y[8]) {
  float sum = 0.;
  for (uint32_t j = 0; j < 8; ++j) {
    for (uint32_t i = 0; i < 8; ++i) {
      sum += in[i + j * 8] * cos_x[i] * cos_y[j];
    }
  }
  return (int32_t)lrintf(256.f * sum);
}

//------------------------------------------------------------------------------

WP2TransposeF WP2Transpose = nullptr;
int32_t (*WP2SlowDct8x8)(const int32_t in[64],
                         const float cos_x[8], const float cos_y[8]) = nullptr;
WP2TransformF WP2InvDct[5] = { nullptr }, WP2FwdDct[5] = { nullptr };
WP2TransformF WP2InvAdst[5] = { nullptr }, WP2FwdAdst[5] = { nullptr };
WP2TransformF WP2InvHadamard[5] = { nullptr }, WP2FwdHadamard[5] = { nullptr };

static volatile WP2CPUInfo transform_last_cpuinfo_used =
    (WP2CPUInfo)&transform_last_cpuinfo_used;

extern void WP2TransformInitSSE();

WP2_TSAN_IGNORE_FUNCTION void WP2TransformInit() {
  if (transform_last_cpuinfo_used == WP2GetCPUInfo) return;

  WP2Transpose = Transpose_C;

  WP2InvDct[0] = idct2_C;
  WP2InvDct[1] = idct4_C;
  WP2InvDct[2] = idct8_C;
  WP2InvDct[3] = idct16_C;
  WP2InvDct[4] = idct32_C;
  WP2FwdDct[0] = fdct2_C;
  WP2FwdDct[1] = fdct4_C;
  WP2FwdDct[2] = fdct8_C;
  WP2FwdDct[3] = fdct16_C;
  WP2FwdDct[4] = fdct32_C;
  WP2InvAdst[0] = idct2_C;
  WP2InvAdst[1] = iadst4_C;
  WP2InvAdst[2] = iadst8_C;
  WP2InvAdst[3] = iadst16_C;
  WP2InvAdst[4] = iadst32_C;
  WP2FwdAdst[0] = fdct2_C;
  WP2FwdAdst[1] = fadst4_C;
  WP2FwdAdst[2] = fadst8_C;
  WP2FwdAdst[3] = fadst16_C;
  WP2FwdAdst[4] = fadst32_C;
  WP2InvHadamard[0] = iHadamard2_C;
  WP2InvHadamard[1] = iHadamard4_C;
  WP2InvHadamard[2] = iHadamard8_C;
  WP2InvHadamard[3] = iHadamard16_C;
  WP2InvHadamard[4] = iHadamard32_C;
  WP2FwdHadamard[0] = fHadamard2_C;
  WP2FwdHadamard[1] = fHadamard4_C;
  WP2FwdHadamard[2] = fHadamard8_C;
  WP2FwdHadamard[3] = fHadamard16_C;
  WP2FwdHadamard[4] = fHadamard32_C;

  WP2SlowDct8x8 = SlowDct8x8_C;

  if (WP2GetCPUInfo != nullptr) {
#if defined(WP2_USE_SSE)
    if (WP2GetCPUInfo(kSSE)) WP2TransformInitSSE();
#endif
  }

  transform_last_cpuinfo_used = WP2GetCPUInfo;
}

//------------------------------------------------------------------------------
// Entry points

static inline uint32_t SizeToIdx(uint32_t size) {
  assert(size == 2 || size == 4 || size == 8 || size == 16 || size == 32);
  return (size == 2) ? 0 : (size == 4) ? 1 : (size == 8) ? 2 :
         (size == 16) ? 3 : 4;
}

//------------------------------------------------------------------------------

void WP2InvTransform(const int32_t src[], WP2TransformType type, uint32_t size,
                     int32_t dst[]) {
  const uint32_t idx = SizeToIdx(size);
  // clang-format off
  switch (type) {
    case WP2TransformType::kDct: WP2InvDct[idx](src, dst); break;
    case WP2TransformType::kAdst: WP2InvAdst[idx](src, dst); break;
    case WP2TransformType::kHadamard: WP2InvHadamard[idx](src, dst); break;
    case WP2TransformType::kIdentity:
      iIdentity(src, dst, size, /*divider=*/size / 2);
      break;
    default:
      assert(false);
  }
  // clang-format on
}

void WP2Transform(const int32_t src[], WP2TransformType type, uint32_t size,
                  int32_t dst[]) {
  const uint32_t idx = SizeToIdx(size);
  // clang-format off
  switch (type) {
    case WP2TransformType::kDct: WP2FwdDct[idx](src, dst); break;
    case WP2TransformType::kAdst: WP2FwdAdst[idx](src, dst); break;
    case WP2TransformType::kHadamard: WP2FwdHadamard[idx](src, dst); break;
    case WP2TransformType::kIdentity:
      Identity(src, dst, size, /*multiplier=*/(size / 2));
      break;
    default:
      assert(false);
  }
  // clang-format on
}

static void Transform2D(const int32_t src[], WP2TransformType tf_x,
                        WP2TransformType tf_y, uint32_t size_x,
                        uint32_t size_y, int32_t dst[]) {
  int32_t tmp[32 * 32];
  const uint32_t tw = size_x;
  const uint32_t th = size_y;
  for (size_t j = 0; j < th; ++j) {
    WP2Transform(src + j * tw, tf_x, size_x, tmp + j * tw);
  }
  WP2Transpose(tmp, tw, th, dst);  // tw x th -> th x tw
  for (size_t j = 0; j < tw; ++j) {
    WP2Transform(dst + j * th, tf_y, size_y, tmp + j * th);
  }
  WP2Transpose(tmp, th, tw, dst);  // th x tw -> tw x th
}

// 'x_end' is the first x for which all coeffs are 0.
// TODO(vrabaud) the DC only case can be sped up.
static void InvTransform2D(const int32_t src[], WP2TransformType tf_x,
                           WP2TransformType tf_y, uint32_t size_x,
                           uint32_t x_end, uint32_t size_y, int32_t dst[]) {
  int32_t tmp[32 * 32];
  const uint32_t tw = size_x;
  const uint32_t th = size_y;
  WP2Transpose(src, tw, th, tmp);  // tw x th -> th x tw
  for (uint32_t j = 0; j < x_end; ++j) {
    WP2InvTransform(tmp + j * th, tf_y, size_y, dst + j * th);
  }
  std::fill(dst + x_end * th, dst + tw * th, 0);
  WP2Transpose(dst, th, tw, tmp);  // th x tw -> tw x th
  for (size_t j = 0; j < th; ++j) {
    WP2InvTransform(tmp + j * tw, tf_x, size_x, dst + j * tw);
  }
}

//------------------------------------------------------------------------------

void WP2Transform2D(const int32_t src[], WP2TransformType tf_x,
                    WP2TransformType tf_y, uint32_t size_x, uint32_t size_y,
                    int32_t dst[], bool reduced) {
  if (!reduced) return Transform2D(src, tf_x, tf_y, size_x, size_y, dst);
  const uint32_t tw = size_x / 2;
  const uint32_t th = size_y / 2;

  int32_t tmp[32 * 32];
  const int32_t* coeffs;
  if (src == dst) {
    memcpy(tmp, src, 2 * 2 * tw * th * sizeof(tmp[0]));
    coeffs = tmp;
  } else {
    coeffs = src;
  }

  for (uint32_t y = 0; y < th; ++y) {
    for (uint32_t x = 0; x < tw; ++x) {
      const size_t offset = 2 * y * 2 * tw + 2 * x;
      dst[x + tw * y] = (coeffs[offset +          0] +
                         coeffs[offset +          1] +
                         coeffs[offset + 2 * tw + 0] +
                         coeffs[offset + 2 * tw + 1] + 2) >> 2;
    }
  }
  Transform2D(dst, tf_x, tf_y, tw, th, dst);
}

// Upsampling taps
static constexpr int32_t kTaps[4][9] = {
  {   64,     9,    -9,     9,    -9,     2,    -1,    -1,     0 },
  {   64,    -9,     9,     9,    -9,    -1,     2,     0,    -1 },
  {   64,     9,    -9,    -9,     9,    -1,     0,     2,    -1 },
  {   64,    -9,     9,    -9,     9,     0,    -1,    -1,     2 },
};

void WP2InvTransform2D(const int32_t src[], WP2TransformType tf_x,
                       WP2TransformType tf_y, uint32_t size_x, uint32_t size_y,
                       int32_t dst[], bool reduced) {
  const uint32_t tw = reduced ? size_x / 2 : size_x;
  const uint32_t th = reduced ? size_y / 2 : size_y;
  // Figure out the bounding column.
  uint32_t x_end = 0;
  for (uint32_t j = 0; j < th; ++j) {
    for (uint32_t i = x_end; i < tw; ++i) {
      if (src[j * tw + i] != 0) x_end = i + 1;
    }
  }
  // If the residuals are null.
  if (x_end == 0) {
    std::fill(dst, dst + size_x * size_y, 0);
    return;
  }

  if (!reduced) {
    return InvTransform2D(src, tf_x, tf_y, size_x, x_end, size_y, dst);
  }
  int32_t tmp[16 * 16];
  assert(tw <= 16 && th <= 16);
  InvTransform2D(src, tf_x, tf_y, tw, x_end, th, tmp);

  // Check if we are DC only.
  if (x_end == 1 && tf_x == kDct && tf_y == kDct) {
    uint32_t y_end = th;
    while (y_end > 0 && src[(y_end - 1) * tw] == 0) --y_end;
    assert(y_end > 0);  // the all zero case has been dealt with before
    if (y_end == 1) {
      std::fill(dst, dst + size_x * size_y, tmp[0]);
      return;
    }
  }

  // Upsample
  for (uint32_t y = 0; y < th; ++y) {
    for (uint32_t x = 0; x < tw; ++x) {
      const uint32_t xp = (x + 1 < tw) ? x + 1 : x;
      const uint32_t yp = (y + 1 < th) ? y + 1 : y;
      const uint32_t xm = (x > 0) ? x - 1 : 0;
      const uint32_t ym = (y > 0) ? y - 1 : 0;
      // Source sample numbering (samples[]):
      // 5...3...6
      // .   .   .
      // .   .   .
      // .  x.x  .
      // 1...0...2
      // .  x.x  .
      // .   .   .
      // .   .   .
      // 7...4...8
      // (interpolated samples v[] marked with 'x')
      const int32_t samples[9] = {
          tmp[x  + tw * y],                      // #0
          tmp[xm + tw * y],  tmp[xp + tw * y],   // #1 #2
          tmp[x  + tw * ym], tmp[x  + tw * yp],  // #3 #4
          tmp[xm + tw * ym], tmp[xp + tw * ym],  // #5 #6
          tmp[xm + tw * yp], tmp[xp + tw * yp]   // #7 #8
      };
      int32_t v[4] = {32, 32, 32, 32};
      for (uint32_t i = 0; i < 9; ++i) {
        v[0] += kTaps[0][i] * samples[i];
        v[1] += kTaps[1][i] * samples[i];
        v[2] += kTaps[2][i] * samples[i];
        v[3] += kTaps[3][i] * samples[i];
      }
      const size_t offset = 2 * y * 2 * tw + 2 * x;
      dst[offset +          0] = v[0] >> 6;
      dst[offset +          1] = v[1] >> 6;
      dst[offset + 2 * tw + 0] = v[2] >> 6;
      dst[offset + 2 * tw + 1] = v[3] >> 6;
    }
  }
}

//------------------------------------------------------------------------------

void WP2ReduceCoeffs(const int32_t src[], uint32_t size_x, uint32_t size_y,
                     int32_t dst[]) {
  const uint32_t half_size_x = size_x / 2;
  const uint32_t half_size_y = size_y / 2;
  for (uint32_t y = 0; y < half_size_y; ++y) {
    for (uint32_t x = 0; x < half_size_x; ++x) {
      const int32_t v0 = src[x + size_x * y];
      // rounds toward 0 with a +1 bias, using >>'s instead of /2
      dst[x + half_size_x * y] = (v0 + (v0 < 0)) >> 1;
    }
  }
}

//------------------------------------------------------------------------------
