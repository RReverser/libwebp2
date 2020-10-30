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
// Image transform methods for lossless encoder.
//
// Authors: Vincent Rabaud (vrabaud@google.com)

#include "src/dsp/lossless/lossless.h"

#include <cstddef>
#include <cstdint>

#include "src/dsp/dsp.h"
#include "src/dsp/lossless/lossless_common.h"

const WP2LPrefixCode kWP2PrefixEncodeCode[PREFIX_LOOKUP_IDX_MAX] = {
  { 0, 0}, { 0, 0}, { 1, 0}, { 2, 0}, { 3, 0}, { 4, 1}, { 4, 1}, { 5, 1},
  { 5, 1}, { 6, 2}, { 6, 2}, { 6, 2}, { 6, 2}, { 7, 2}, { 7, 2}, { 7, 2},
  { 7, 2}, { 8, 3}, { 8, 3}, { 8, 3}, { 8, 3}, { 8, 3}, { 8, 3}, { 8, 3},
  { 8, 3}, { 9, 3}, { 9, 3}, { 9, 3}, { 9, 3}, { 9, 3}, { 9, 3}, { 9, 3},
  { 9, 3}, {10, 4}, {10, 4}, {10, 4}, {10, 4}, {10, 4}, {10, 4}, {10, 4},
  {10, 4}, {10, 4}, {10, 4}, {10, 4}, {10, 4}, {10, 4}, {10, 4}, {10, 4},
  {10, 4}, {11, 4}, {11, 4}, {11, 4}, {11, 4}, {11, 4}, {11, 4}, {11, 4},
  {11, 4}, {11, 4}, {11, 4}, {11, 4}, {11, 4}, {11, 4}, {11, 4}, {11, 4},
  {11, 4}, {12, 5}, {12, 5}, {12, 5}, {12, 5}, {12, 5}, {12, 5}, {12, 5},
  {12, 5}, {12, 5}, {12, 5}, {12, 5}, {12, 5}, {12, 5}, {12, 5}, {12, 5},
  {12, 5}, {12, 5}, {12, 5}, {12, 5}, {12, 5}, {12, 5}, {12, 5}, {12, 5},
  {12, 5}, {12, 5}, {12, 5}, {12, 5}, {12, 5}, {12, 5}, {12, 5}, {12, 5},
  {12, 5}, {13, 5}, {13, 5}, {13, 5}, {13, 5}, {13, 5}, {13, 5}, {13, 5},
  {13, 5}, {13, 5}, {13, 5}, {13, 5}, {13, 5}, {13, 5}, {13, 5}, {13, 5},
  {13, 5}, {13, 5}, {13, 5}, {13, 5}, {13, 5}, {13, 5}, {13, 5}, {13, 5},
  {13, 5}, {13, 5}, {13, 5}, {13, 5}, {13, 5}, {13, 5}, {13, 5}, {13, 5},
  {13, 5}, {14, 6}, {14, 6}, {14, 6}, {14, 6}, {14, 6}, {14, 6}, {14, 6},
  {14, 6}, {14, 6}, {14, 6}, {14, 6}, {14, 6}, {14, 6}, {14, 6}, {14, 6},
  {14, 6}, {14, 6}, {14, 6}, {14, 6}, {14, 6}, {14, 6}, {14, 6}, {14, 6},
  {14, 6}, {14, 6}, {14, 6}, {14, 6}, {14, 6}, {14, 6}, {14, 6}, {14, 6},
  {14, 6}, {14, 6}, {14, 6}, {14, 6}, {14, 6}, {14, 6}, {14, 6}, {14, 6},
  {14, 6}, {14, 6}, {14, 6}, {14, 6}, {14, 6}, {14, 6}, {14, 6}, {14, 6},
  {14, 6}, {14, 6}, {14, 6}, {14, 6}, {14, 6}, {14, 6}, {14, 6}, {14, 6},
  {14, 6}, {14, 6}, {14, 6}, {14, 6}, {14, 6}, {14, 6}, {14, 6}, {14, 6},
  {14, 6}, {15, 6}, {15, 6}, {15, 6}, {15, 6}, {15, 6}, {15, 6}, {15, 6},
  {15, 6}, {15, 6}, {15, 6}, {15, 6}, {15, 6}, {15, 6}, {15, 6}, {15, 6},
  {15, 6}, {15, 6}, {15, 6}, {15, 6}, {15, 6}, {15, 6}, {15, 6}, {15, 6},
  {15, 6}, {15, 6}, {15, 6}, {15, 6}, {15, 6}, {15, 6}, {15, 6}, {15, 6},
  {15, 6}, {15, 6}, {15, 6}, {15, 6}, {15, 6}, {15, 6}, {15, 6}, {15, 6},
  {15, 6}, {15, 6}, {15, 6}, {15, 6}, {15, 6}, {15, 6}, {15, 6}, {15, 6},
  {15, 6}, {15, 6}, {15, 6}, {15, 6}, {15, 6}, {15, 6}, {15, 6}, {15, 6},
  {15, 6}, {15, 6}, {15, 6}, {15, 6}, {15, 6}, {15, 6}, {15, 6}, {15, 6},
  {15, 6}, {16, 7}, {16, 7}, {16, 7}, {16, 7}, {16, 7}, {16, 7}, {16, 7},
  {16, 7}, {16, 7}, {16, 7}, {16, 7}, {16, 7}, {16, 7}, {16, 7}, {16, 7},
  {16, 7}, {16, 7}, {16, 7}, {16, 7}, {16, 7}, {16, 7}, {16, 7}, {16, 7},
  {16, 7}, {16, 7}, {16, 7}, {16, 7}, {16, 7}, {16, 7}, {16, 7}, {16, 7},
  {16, 7}, {16, 7}, {16, 7}, {16, 7}, {16, 7}, {16, 7}, {16, 7}, {16, 7},
  {16, 7}, {16, 7}, {16, 7}, {16, 7}, {16, 7}, {16, 7}, {16, 7}, {16, 7},
  {16, 7}, {16, 7}, {16, 7}, {16, 7}, {16, 7}, {16, 7}, {16, 7}, {16, 7},
  {16, 7}, {16, 7}, {16, 7}, {16, 7}, {16, 7}, {16, 7}, {16, 7}, {16, 7},
  {16, 7}, {16, 7}, {16, 7}, {16, 7}, {16, 7}, {16, 7}, {16, 7}, {16, 7},
  {16, 7}, {16, 7}, {16, 7}, {16, 7}, {16, 7}, {16, 7}, {16, 7}, {16, 7},
  {16, 7}, {16, 7}, {16, 7}, {16, 7}, {16, 7}, {16, 7}, {16, 7}, {16, 7},
  {16, 7}, {16, 7}, {16, 7}, {16, 7}, {16, 7}, {16, 7}, {16, 7}, {16, 7},
  {16, 7}, {16, 7}, {16, 7}, {16, 7}, {16, 7}, {16, 7}, {16, 7}, {16, 7},
  {16, 7}, {16, 7}, {16, 7}, {16, 7}, {16, 7}, {16, 7}, {16, 7}, {16, 7},
  {16, 7}, {16, 7}, {16, 7}, {16, 7}, {16, 7}, {16, 7}, {16, 7}, {16, 7},
  {16, 7}, {16, 7}, {16, 7}, {16, 7}, {16, 7}, {16, 7}, {16, 7}, {16, 7},
  {16, 7}, {17, 7}, {17, 7}, {17, 7}, {17, 7}, {17, 7}, {17, 7}, {17, 7},
  {17, 7}, {17, 7}, {17, 7}, {17, 7}, {17, 7}, {17, 7}, {17, 7}, {17, 7},
  {17, 7}, {17, 7}, {17, 7}, {17, 7}, {17, 7}, {17, 7}, {17, 7}, {17, 7},
  {17, 7}, {17, 7}, {17, 7}, {17, 7}, {17, 7}, {17, 7}, {17, 7}, {17, 7},
  {17, 7}, {17, 7}, {17, 7}, {17, 7}, {17, 7}, {17, 7}, {17, 7}, {17, 7},
  {17, 7}, {17, 7}, {17, 7}, {17, 7}, {17, 7}, {17, 7}, {17, 7}, {17, 7},
  {17, 7}, {17, 7}, {17, 7}, {17, 7}, {17, 7}, {17, 7}, {17, 7}, {17, 7},
  {17, 7}, {17, 7}, {17, 7}, {17, 7}, {17, 7}, {17, 7}, {17, 7}, {17, 7},
  {17, 7}, {17, 7}, {17, 7}, {17, 7}, {17, 7}, {17, 7}, {17, 7}, {17, 7},
  {17, 7}, {17, 7}, {17, 7}, {17, 7}, {17, 7}, {17, 7}, {17, 7}, {17, 7},
  {17, 7}, {17, 7}, {17, 7}, {17, 7}, {17, 7}, {17, 7}, {17, 7}, {17, 7},
  {17, 7}, {17, 7}, {17, 7}, {17, 7}, {17, 7}, {17, 7}, {17, 7}, {17, 7},
  {17, 7}, {17, 7}, {17, 7}, {17, 7}, {17, 7}, {17, 7}, {17, 7}, {17, 7},
  {17, 7}, {17, 7}, {17, 7}, {17, 7}, {17, 7}, {17, 7}, {17, 7}, {17, 7},
  {17, 7}, {17, 7}, {17, 7}, {17, 7}, {17, 7}, {17, 7}, {17, 7}, {17, 7},
  {17, 7}, {17, 7}, {17, 7}, {17, 7}, {17, 7}, {17, 7}, {17, 7}, {17, 7},
};

const uint8_t kWP2PrefixEncodeExtraBitsValue[PREFIX_LOOKUP_IDX_MAX] = {
  0,  0,  0,  0,  0,  0,  1,  0,  1,  0,  1,  2,  3,  0,  1,  2,  3,
  0,  1,  2,  3,  4,  5,  6,  7,  0,  1,  2,  3,  4,  5,  6,  7,
  0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14, 15,
  0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14, 15,
  0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14, 15,
 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31,
  0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14, 15,
 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31,
  0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14, 15,
 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31,
 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47,
 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63,
  0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14, 15,
 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31,
 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47,
 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63,
  0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14, 15,
 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31,
 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47,
 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63,
 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79,
 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95,
 96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111,
 112, 113, 114, 115, 116, 117, 118, 119, 120, 121, 122, 123, 124, 125, 126,
 127,
  0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14, 15,
 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31,
 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47,
 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63,
 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79,
 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95,
 96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111,
 112, 113, 114, 115, 116, 117, 118, 119, 120, 121, 122, 123, 124, 125, 126
};

//------------------------------------------------------------------------------
// Methods to calculate Entropy (Shannon).

// Compute the combined Shanon's entropy for distribution {X} and {X+Y}
static float CombinedShannonEntropy_C(const uint32_t* X, const uint32_t* Y,
                                      uint32_t length) {
  double retval = 0.;
  uint32_t sumX = 0, sumXY = 0;
  for (uint32_t i = 0; i < length; ++i) {
    const uint32_t x = X[i];
    if (x != 0) {
      const uint32_t xy = x + Y[i];
      sumX += x;
      retval -= WP2SLog2(x);
      sumXY += xy;
      retval -= WP2SLog2(xy);
    } else if (Y[i] != 0) {
      sumXY += Y[i];
      retval -= WP2SLog2(Y[i]);
    }
  }
  retval += WP2SLog2(sumX) + WP2SLog2(sumXY);
  return (float)retval;
}

//------------------------------------------------------------------------------

void WP2LSubtractGreenFromBlueAndRed_C(uint16_t* argb_data,
                                       uint32_t channel_bits,
                                       uint32_t num_pixels) {
  uint16_t mask = (1u << channel_bits) - 1;
  for (size_t i = 0; i < num_pixels; ++i) {
    const int green = argb_data[4 * i + 2];
    argb_data[4 * i + 1] = (argb_data[4 * i + 1] - green) & mask;
    argb_data[4 * i + 3] = (argb_data[4 * i + 3] - green) & mask;
  }
}

static inline int32_t ColorTransformDelta(int16_t color_pred, int16_t color) {
  return (((int32_t)color_pred * (int8_t)color) >> 5);
}

void WP2LTransformColor_C(const WP2LMultipliers* const m, int num_pixels,
                          int16_t* const data) {
  int i;
  for (i = 0; i < num_pixels; ++i) {
    const int16_t green = data[4 * i + 2];
    const int16_t red = data[4 * i + 1];
    int new_red = red;
    int new_blue = data[4 * i + 3];
    new_red -= ColorTransformDelta(m->green_to_red, green);
    new_blue -= ColorTransformDelta(m->green_to_blue, green);
    new_blue -= ColorTransformDelta(m->red_to_blue, red);
    data[4 * i + 1] = new_red;
    data[4 * i + 3] = new_blue;
  }
}

static inline int16_t TransformColorRed(int16_t green_to_red,
                                        const int16_t* const argb) {
  const int16_t green = argb[2];
  int16_t new_red = argb[1];
  new_red -= ColorTransformDelta(green_to_red, green);
  return new_red;
}

static inline int16_t TransformColorBlue(int16_t green_to_blue,
                                         int16_t red_to_blue,
                                         const int16_t* const argb) {
  const int16_t green = argb[2];
  const int16_t red = argb[1];
  int16_t new_blue = argb[3];
  new_blue -= ColorTransformDelta(green_to_blue, green);
  new_blue -= ColorTransformDelta(red_to_blue, red);
  return new_blue;
}

void WP2LCollectColorRedTransforms_C(const int16_t* argb, uint32_t width,
                                     uint32_t tile_width, uint32_t tile_height,
                                     int16_t green_to_red,
                                     uint32_t* const histo) {
  while (tile_height-- > 0) {
    for (uint32_t x = 0; x < tile_width; ++x) {
      ++histo[MakeIndex(TransformColorRed(green_to_red, &argb[4 * x]))];
    }
    argb += 4 * width;
  }
}

void WP2LCollectColorBlueTransforms_C(const int16_t* argb, uint32_t width,
                                      uint32_t tile_width, uint32_t tile_height,
                                      int16_t green_to_blue,
                                      int16_t red_to_blue,
                                      uint32_t* const histo) {
  while (tile_height-- > 0) {
    for (uint32_t x = 0; x < tile_width; ++x) {
      ++histo[MakeIndex(
          TransformColorBlue(green_to_blue, red_to_blue, &argb[4 * x]))];
    }
    argb += 4 * width;
  }
}

//------------------------------------------------------------------------------

static double ExtraCost_C(const uint32_t* population, int length) {
  int i;
  double cost = 0.;
  for (i = 2; i < length - 2; ++i) cost += (i >> 1) * population[i + 2];
  return cost;
}

static double ExtraCostCombined_C(const uint32_t* X, const uint32_t* Y,
                                  int length) {
  int i;
  double cost = 0.;
  for (i = 2; i < length - 2; ++i) {
    const int xy = X[i + 2] + Y[i + 2];
    cost += (i >> 1) * xy;
  }
  return cost;
}

//------------------------------------------------------------------------------

static void BufferAdd_C(const uint32_t* const a, const uint32_t* const b,
                        uint32_t size, uint32_t* const out,
                        uint32_t* const nonzeros) {
  assert(out != nullptr);
  if (nonzeros != nullptr) *nonzeros = 0;
  if (b != out) {
    for (uint32_t i = 0; i < size; ++i) {
      out[i] = a[i] + b[i];
      if (nonzeros != nullptr && out[i]) ++*nonzeros;
    }
  } else {
    for (uint32_t i = 0; i < size; ++i) {
      out[i] += a[i];
      if (nonzeros != nullptr && out[i]) ++*nonzeros;
    }
  }
}

//------------------------------------------------------------------------------
// Predictors

static void Predictor2(const uint16_t* left, const uint16_t* top,
                       uint16_t max_value, uint16_t* out) {
  (void)left;
  (void)max_value;
  out[0] = top[0];
  out[1] = top[1];
  out[2] = top[2];
  out[3] = top[3];
}
static void Predictor3(const uint16_t* left, const uint16_t* top,
                       uint16_t max_value, uint16_t* out) {
  (void)left;
  (void)max_value;
  out[0] = top[4];
  out[1] = top[5];
  out[2] = top[6];
  out[3] = top[7];
}
static void Predictor4(const uint16_t* left, const uint16_t* top,
                       uint16_t max_value, uint16_t* out) {
  (void)left;
  (void)max_value;
  out[0] = top[-4];
  out[1] = top[-3];
  out[2] = top[-2];
  out[3] = top[-1];
}
static void Predictor5(const uint16_t* left, const uint16_t* top,
                       uint16_t max_value, uint16_t* out) {
  Average3(left, top, top + 4, max_value, out);
}
static void Predictor6(const uint16_t* left, const uint16_t* top,
                       uint16_t max_value, uint16_t* out) {
  Average2(left, top - 4, max_value, out);
}
static void Predictor7(const uint16_t* left, const uint16_t* top,
                       uint16_t max_value, uint16_t* out) {
  Average2(left, top, max_value, out);
}
static void Predictor8(const uint16_t* left, const uint16_t* top,
                       uint16_t max_value, uint16_t* out) {
  (void)left;
  Average2(top - 4, top, max_value, out);
}
static void Predictor9(const uint16_t* left, const uint16_t* top,
                       uint16_t max_value, uint16_t* out) {
  (void)left;
  Average2(top, top + 4, max_value, out);
}
static void Predictor10(const uint16_t* left, const uint16_t* top,
                        uint16_t max_value, uint16_t* out) {
  Average4(left, top - 4, top, top + 4, max_value, out);
}
static void Predictor11(const uint16_t* left, const uint16_t* top,
                        uint16_t max_value, uint16_t* out) {
  (void)max_value;
  Select(top, left, top - 4, out);
}
static void Predictor12(const uint16_t* left, const uint16_t* top,
                        uint16_t max_value, uint16_t* out) {
  ClampedAddSubtractFull(left, top, top - 4, max_value, out);
}
static void Predictor13(const uint16_t* left, const uint16_t* top,
                        uint16_t max_value, uint16_t* out) {
  ClampedAddSubtractHalf(left, top, top - 4, max_value, out);
}

//------------------------------------------------------------------------------

static void PredictorSub0_C(const uint16_t* in, bool has_alpha,
                            const uint16_t* const upper, uint32_t num_pixels,
                            uint32_t num_bits, uint16_t* out) {
  (void)upper;
  const uint16_t max_value = (1u << num_bits) - 1;
  for (const uint16_t* const in_end = in + 4 * num_pixels; in < in_end;
       in += 4, out += 4) {
    out[0] = has_alpha ? ((in[0] - max_value) & max_value) : WP2::kAlphaMax;
    out[1] = in[1];
    out[2] = in[2];
    out[3] = in[3];
  }
}

static void PredictorSub1_C(const uint16_t* in, bool has_alpha,
                            const uint16_t* const upper, uint32_t num_pixels,
                            uint32_t num_bits, uint16_t* out) {
  (void)upper;
  const uint16_t mask = (1u << num_bits) - 1;
  for (const uint16_t* const in_end = in + 4 * num_pixels; in < in_end;
       in += 4, out += 4) {
    WP2LSubPixels(in, has_alpha, &in[-4], mask, out);
  }
}

GENERATE_PREDICTOR_SUB(Predictor2, PredictorSub2_C)
GENERATE_PREDICTOR_SUB(Predictor3, PredictorSub3_C)
GENERATE_PREDICTOR_SUB(Predictor4, PredictorSub4_C)
GENERATE_PREDICTOR_SUB(Predictor5, PredictorSub5_C)
GENERATE_PREDICTOR_SUB(Predictor6, PredictorSub6_C)
GENERATE_PREDICTOR_SUB(Predictor7, PredictorSub7_C)
GENERATE_PREDICTOR_SUB(Predictor8, PredictorSub8_C)
GENERATE_PREDICTOR_SUB(Predictor9, PredictorSub9_C)
GENERATE_PREDICTOR_SUB(Predictor10, PredictorSub10_C)
GENERATE_PREDICTOR_SUB(Predictor11, PredictorSub11_C)
GENERATE_PREDICTOR_SUB(Predictor12, PredictorSub12_C)
GENERATE_PREDICTOR_SUB(Predictor13, PredictorSub13_C)

//------------------------------------------------------------------------------

WP2LProcessEncBlueAndRedFunc WP2LSubtractGreenFromBlueAndRed;

WP2LTransformColorFunc WP2LTransformColor;

WP2LCollectColorBlueTransformsFunc WP2LCollectColorBlueTransforms;
WP2LCollectColorRedTransformsFunc WP2LCollectColorRedTransforms;

WP2LCostFunc WP2LExtraCost;
WP2LCostCombinedFunc WP2LExtraCostCombined;
WP2LCombinedShannonEntropyFunc WP2LCombinedShannonEntropy;

WP2LBufferAddFunc WP2LBufferAdd;

WP2LPredictorAddSubFunc WP2LPredictorsSub[kNumPredictors];
WP2LPredictorAddSubFunc WP2LPredictorsSub_C[kNumPredictors];

static volatile WP2CPUInfo lossless_enc_last_cpuinfo_used =
    (WP2CPUInfo)&lossless_enc_last_cpuinfo_used;

WP2_TSAN_IGNORE_FUNCTION void WP2LEncDspInit() {
  if (lossless_enc_last_cpuinfo_used == WP2GetCPUInfo) return;

  WP2LDspInit();

  WP2LSubtractGreenFromBlueAndRed = WP2LSubtractGreenFromBlueAndRed_C;

  WP2LTransformColor = WP2LTransformColor_C;

  WP2LCollectColorBlueTransforms = WP2LCollectColorBlueTransforms_C;
  WP2LCollectColorRedTransforms = WP2LCollectColorRedTransforms_C;

  WP2LExtraCost = ExtraCost_C;
  WP2LExtraCostCombined = ExtraCostCombined_C;
  WP2LCombinedShannonEntropy = CombinedShannonEntropy_C;

  WP2LBufferAdd = BufferAdd_C;

  WP2LPredictorsSub[0] = PredictorSub0_C;
  WP2LPredictorsSub[1] = PredictorSub1_C;
  WP2LPredictorsSub[2] = PredictorSub2_C;
  WP2LPredictorsSub[3] = PredictorSub3_C;
  WP2LPredictorsSub[4] = PredictorSub4_C;
  WP2LPredictorsSub[5] = PredictorSub5_C;
  WP2LPredictorsSub[6] = PredictorSub6_C;
  WP2LPredictorsSub[7] = PredictorSub7_C;
  WP2LPredictorsSub[8] = PredictorSub8_C;
  WP2LPredictorsSub[9] = PredictorSub9_C;
  WP2LPredictorsSub[10] = PredictorSub10_C;
  WP2LPredictorsSub[11] = PredictorSub11_C;
  WP2LPredictorsSub[12] = PredictorSub12_C;
  WP2LPredictorsSub[13] = PredictorSub13_C;
  WP2LPredictorsSub[14] = PredictorSub0_C;  // <- padding security sentinels
  WP2LPredictorsSub[15] = PredictorSub0_C;

  WP2LPredictorsSub_C[0] = PredictorSub0_C;
  WP2LPredictorsSub_C[1] = PredictorSub1_C;
  WP2LPredictorsSub_C[2] = PredictorSub2_C;
  WP2LPredictorsSub_C[3] = PredictorSub3_C;
  WP2LPredictorsSub_C[4] = PredictorSub4_C;
  WP2LPredictorsSub_C[5] = PredictorSub5_C;
  WP2LPredictorsSub_C[6] = PredictorSub6_C;
  WP2LPredictorsSub_C[7] = PredictorSub7_C;
  WP2LPredictorsSub_C[8] = PredictorSub8_C;
  WP2LPredictorsSub_C[9] = PredictorSub9_C;
  WP2LPredictorsSub_C[10] = PredictorSub10_C;
  WP2LPredictorsSub_C[11] = PredictorSub11_C;
  WP2LPredictorsSub_C[12] = PredictorSub12_C;
  WP2LPredictorsSub_C[13] = PredictorSub13_C;
  WP2LPredictorsSub_C[14] = PredictorSub0_C;  // <- padding security sentinels
  WP2LPredictorsSub_C[15] = PredictorSub0_C;

  // If defined, use CPUInfo() to overwrite some pointers with faster versions.
  if (WP2GetCPUInfo != NULL) {
  // TODO(skal): SSE2, etc.
  }
  lossless_enc_last_cpuinfo_used = WP2GetCPUInfo;
}

//------------------------------------------------------------------------------
