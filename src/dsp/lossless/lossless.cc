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
// Image transforms and color space conversion methods for lossless decoder.
//
// Authors: Vikas Arora (vikaas.arora@gmail.com)
//          Jyrki Alakuijala (jyrki@google.com)
//          Urvang Joshi (urvang@google.com)

#include "src/dsp/lossless/lossless.h"

#include "src/common/constants.h"
#include "src/dec/lossless/losslessi_dec.h"
#include "src/dsp/dsp.h"
#include "src/dsp/lossless/lossless_common.h"

#define MAX_DIFF_COST (1e30f)

//------------------------------------------------------------------------------
// Predictors

static void Predictor0_C(const uint16_t* left, const uint16_t* top,
                         uint16_t max_value, uint16_t* out) {
  (void)top;
  (void)left;
  out[0] = max_value;
  out[1] = 0;
  out[2] = 0;
  out[3] = 0;
}
static void Predictor1_C(const uint16_t* left, const uint16_t* top,
                         uint16_t max_value, uint16_t* out) {
  (void)top;
  (void)max_value;
  out[0] = left[0];
  out[1] = left[1];
  out[2] = left[2];
  out[3] = left[3];
}
static void Predictor2_C(const uint16_t* left, const uint16_t* top,
                         uint16_t max_value, uint16_t* out) {
  (void)left;
  (void)max_value;
  out[0] = top[0];
  out[1] = top[1];
  out[2] = top[2];
  out[3] = top[3];
}
static void Predictor3_C(const uint16_t* left, const uint16_t* top,
                         uint16_t max_value, uint16_t* out) {
  (void)left;
  (void)max_value;
  out[0] = top[4];
  out[1] = top[5];
  out[2] = top[6];
  out[3] = top[7];
}
static void Predictor4_C(const uint16_t* left, const uint16_t* top,
                         uint16_t max_value, uint16_t* out) {
  (void)left;
  (void)max_value;
  out[0] = top[-4];
  out[1] = top[-3];
  out[2] = top[-2];
  out[3] = top[-1];
}
static void Predictor5_C(const uint16_t* left, const uint16_t* top,
                         uint16_t max_value, uint16_t* out) {
  Average3(left, top, top + 4, max_value, out);
}
static void Predictor6_C(const uint16_t* left, const uint16_t* top,
                         uint16_t max_value, uint16_t* out) {
  Average2(left, top - 4, max_value, out);
}
static void Predictor7_C(const uint16_t* left, const uint16_t* top,
                         uint16_t max_value, uint16_t* out) {
  Average2(left, top, max_value, out);
}
static void Predictor8_C(const uint16_t* left, const uint16_t* top,
                         uint16_t max_value, uint16_t* out) {
  (void)left;
  Average2(top - 4, top, max_value, out);
}
static void Predictor9_C(const uint16_t* left, const uint16_t* top,
                         uint16_t max_value, uint16_t* out) {
  (void)left;
  Average2(top, top + 4, max_value, out);
}
static void Predictor10_C(const uint16_t* left, const uint16_t* top,
                          uint16_t max_value, uint16_t* out) {
  Average4(left, top - 4, top, top + 4, max_value, out);
}
static void Predictor11_C(const uint16_t* left, const uint16_t* top,
                          uint16_t max_value, uint16_t* out) {
  (void)left;
  Select(top, left, top - 4, out);
}
static void Predictor12_C(const uint16_t* left, const uint16_t* top,
                          uint16_t max_value, uint16_t* out) {
  ClampedAddSubtractFull(left, top, top - 4, max_value, out);
}
static void Predictor13_C(const uint16_t* left, const uint16_t* top,
                          uint16_t max_value, uint16_t* out) {
  ClampedAddSubtractHalf(left, top, top - 4, max_value, out);
}

typedef void (*PredictoAddFunc)(const uint16_t* in, const uint16_t* upper,
                                uint32_t num_pixels, uint32_t num_bits,
                                uint16_t* out);

static void PredictorAdd0_C(const uint16_t* in, bool has_alpha, const uint16_t*,
                            uint32_t num_pixels, uint32_t num_bits,
                            uint16_t* out) {
  uint16_t pred[4];
  const uint16_t max_value = (1u << num_bits) - 1;
  Predictor0_C(/*left=*/nullptr, /*top=*/nullptr, max_value, pred);
  for (const uint16_t* in_end = in + 4 * num_pixels; in < in_end;
       in += 4, out += 4) {
    WP2LAddPixels(in, has_alpha, pred, max_value, out);
  }
}
static void PredictorAdd1_C(const uint16_t* in, bool has_alpha, const uint16_t*,
                            uint32_t num_pixels, uint32_t num_bits,
                            uint16_t* out) {
  const uint16_t max_value = (1u << num_bits) - 1;
  for (const uint16_t* in_end = in + 4 * num_pixels; in < in_end;
       in += 4, out += 4) {
    WP2LAddPixels(in, has_alpha, &out[-4], max_value, out);
  }
}
GENERATE_PREDICTOR_ADD(Predictor2_C, PredictorAdd2_C)
GENERATE_PREDICTOR_ADD(Predictor3_C, PredictorAdd3_C)
GENERATE_PREDICTOR_ADD(Predictor4_C, PredictorAdd4_C)
GENERATE_PREDICTOR_ADD(Predictor5_C, PredictorAdd5_C)
GENERATE_PREDICTOR_ADD(Predictor6_C, PredictorAdd6_C)
GENERATE_PREDICTOR_ADD(Predictor7_C, PredictorAdd7_C)
GENERATE_PREDICTOR_ADD(Predictor8_C, PredictorAdd8_C)
GENERATE_PREDICTOR_ADD(Predictor9_C, PredictorAdd9_C)
GENERATE_PREDICTOR_ADD(Predictor10_C, PredictorAdd10_C)
GENERATE_PREDICTOR_ADD(Predictor11_C, PredictorAdd11_C)
GENERATE_PREDICTOR_ADD(Predictor12_C, PredictorAdd12_C)
GENERATE_PREDICTOR_ADD(Predictor13_C, PredictorAdd13_C)

//------------------------------------------------------------------------------

// Inverse prediction.
static void PredictorInverseTransform_C(const WP2L::Transform* const transform,
                                        int y_start, int y_end,
                                        uint32_t num_bits, const uint16_t* in,
                                        bool has_alpha, uint16_t* out) {
  const int width = transform->width_;
  if (y_start == 0) {  // First Row follows the L (mode=1) mode.
    PredictorAdd0_C(in, has_alpha, nullptr, 1, num_bits, out);
    PredictorAdd1_C(in + 4, has_alpha, nullptr, width - 1, num_bits, out + 4);
    in += 4 * width;
    out += 4 * width;
    ++y_start;
  }

  {
    int y = y_start;
    const int tile_width = 1 << transform->bits_;
    const int mask = tile_width - 1;
    const int tiles_per_row = WP2LSubSampleSize(width, transform->bits_);
    const int16_t* pred_mode_base =
        transform->data_.data() + 4 * (y >> transform->bits_) * tiles_per_row;

    while (y < y_end) {
      const int16_t* pred_mode_src = pred_mode_base;
      int x = 1;
      // First pixel follows the T (mode=2) mode.
      PredictorAdd2_C(in, has_alpha, out - 4 * width, 1, num_bits, out);
      // .. the rest:
      while (x < width) {
        const WP2LPredictorAddSubFunc pred_func =
            WP2LPredictorsAdd[pred_mode_src[2]];
        int x_end = (x & ~mask) + tile_width;
        if (x_end > width) x_end = width;
        pred_func(in + 4 * x, has_alpha, out + 4 * x - 4 * width, x_end - x,
                  num_bits, out + 4 * x);
        x = x_end;
        pred_mode_src += 4;
      }
      in += 4 * width;
      out += 4 * width;
      ++y;
      if ((y & mask) == 0) {   // Use the same mask, since tiles are squares.
        pred_mode_base += 4 * tiles_per_row;
      }
    }
  }
}

// Add green to blue and red channels (i.e. perform the inverse transform of
// 'subtract green').
void WP2LAddGreenToBlueAndRed_C(const uint16_t* src, uint32_t num_pixels,
                                uint32_t channel_bits, uint16_t* dst) {
  const uint16_t mask = (uint16_t)((1u << channel_bits) - 1);
  for (const uint16_t* src_end = src + 4 * num_pixels; src < src_end;
       src += 4, dst += 4) {
    const uint16_t green = src[2];
    dst[0] = src[0];
    dst[1] = (src[1] + green) & mask;
    dst[2] = src[2];
    dst[3] = (src[3] + green) & mask;
  }
}

static inline int32_t ColorTransformDelta(int16_t color_pred,
                                          uint16_t max_value_half,
                                          uint16_t color) {
  const int color_int =
      (color >= max_value_half) ? (int)color - 2 * max_value_half : (int)color;
  return (((int32_t)color_pred * color_int) >> 5);
}

static inline void ColorCodeToMultipliers(const int16_t* const color_code,
                                          WP2LMultipliers* const m) {
  m->green_to_red = color_code[3];
  m->green_to_blue = color_code[2];
  m->red_to_blue = color_code[1];
}

static void WP2LTransformColorInverse_C(const WP2LMultipliers* const m,
                                        const uint16_t* src,
                                        uint32_t num_pixels, uint32_t num_bits,
                                        uint16_t* dst) {
  const uint16_t max_value = (1u << num_bits) - 1;
  const uint16_t max_value_half = (max_value + 1) / 2;
  for (const uint16_t* src_end = src + 4 * num_pixels; src < src_end;
       src += 4, dst += 4) {
    const uint16_t red = src[1];
    const uint16_t green = src[2];
    int new_red = red;
    int new_blue = src[3];
    new_red += ColorTransformDelta(m->green_to_red, max_value_half, green);
    new_red &= max_value;
    new_blue += ColorTransformDelta(m->green_to_blue, max_value_half, green);
    new_blue += ColorTransformDelta(m->red_to_blue, max_value_half, new_red);
    new_blue &= max_value;
    dst[0] = src[0];
    dst[1] = (uint16_t)new_red;
    dst[2] = src[2];
    dst[3] = (uint16_t)new_blue;
  }
}

// Color space inverse transform.
static void ColorSpaceInverseTransform_C(const WP2L::Transform* const transform,
                                         int y_start, int y_end,
                                         uint32_t num_bits, const uint16_t* src,
                                         uint16_t* dst) {
  const uint32_t width = transform->width_;
  const int tile_width = 1 << transform->bits_;
  const int mask = tile_width - 1;
  const int safe_width = width & ~mask;
  const int remaining_width = width - safe_width;
  const int tiles_per_row = WP2LSubSampleSize(width, transform->bits_);
  int y = y_start;
  const int16_t* pred_row =
      transform->data_.data() + 4 * (y >> transform->bits_) * tiles_per_row;

  while (y < y_end) {
    const int16_t* pred = pred_row;
    WP2LMultipliers m = { 0, 0, 0 };
    const uint16_t* const src_safe_end = src + 4 * safe_width;
    const uint16_t* const src_end = src + 4 * width;
    while (src < src_safe_end) {
      ColorCodeToMultipliers(pred, &m);
      WP2LTransformColorInverse(&m, src, tile_width, num_bits, dst);
      src += 4 * tile_width;
      dst += 4 * tile_width;
      pred += 4;
    }
    if (src < src_end) {  // Left-overs using C-version.
      ColorCodeToMultipliers(pred, &m);
      WP2LTransformColorInverse(&m, src, remaining_width, num_bits, dst);
      src += 4 * remaining_width;
      dst += 4 * remaining_width;
      pred += 4;
    }
    ++y;
    if ((y & mask) == 0) pred_row += 4 * tiles_per_row;
  }
}

static void MapARGB_C(const uint16_t* src, const int16_t* const color_map,
                      uint32_t color_map_size, uint32_t y_start, uint32_t y_end,
                      uint32_t width, uint16_t* dst) {
  for (const uint16_t* const src_end = src + 4 * (y_end - y_start) * width;
       src < src_end; src += 4, dst += 4) {
    if (src[2] >= color_map_size) assert(false);
    // TODO(vrabaud) use the following once switched to int16_t;
    // WP2L::ColorCopy(&color_map[4 * src[2]], dst);
    std::copy(&color_map[4 * src[2]], &color_map[4 * src[2]] + 4, dst);
  }
}
static void ColorIndexInverseTransform_C(const WP2L::Transform* const transform,
                                         uint32_t y_start, uint32_t y_end,
                                         const uint16_t* const src,
                                         uint16_t* const dst) {
  const uint32_t width = transform->width_;
  WP2LMapColor(src, /*color_map=*/transform->data_.data(),
               /*color_map_size=*/transform->data_.size() / 4, y_start, y_end,
               width, dst);
}

void WP2LInverseTransform(const WP2L::Transform* const transform,
                          uint32_t row_start, uint32_t row_end,
                          uint32_t channel_bits, const uint16_t* const in,
                          bool has_alpha, uint16_t* const out) {
  const uint32_t width = transform->width_;
  assert(row_start < row_end);
  assert(row_end <= transform->height_);

  switch (transform->type_) {
    case WP2L::SUBTRACT_GREEN:
      WP2LAddGreenToBlueAndRed(in, (row_end - row_start) * width, channel_bits,
                               out);
      break;
    case WP2L::PREDICTOR_TRANSFORM: {
      PredictorInverseTransform_C(transform, row_start, row_end, channel_bits,
                                  in, has_alpha, out);
      if (row_end != transform->height_) {
        // The last predicted row in this iteration will be the top-pred row
        // for the first row in next iteration.
        std::copy(out + 4 * (row_end - row_start - 1) * width,
                  out + 4 * (row_end - row_start) * width, out - 4 * width);
      }
      break;
    }
    case WP2L::CROSS_COLOR_TRANSFORM:
      ColorSpaceInverseTransform_C(transform, row_start, row_end, channel_bits,
                                   in, out);
      break;
    case WP2L::GROUP4:
    case WP2L::COLOR_INDEXING_TRANSFORM:
      ColorIndexInverseTransform_C(transform, row_start, row_end, in, out);
      break;
    default:
      assert(false);
  }
}

//------------------------------------------------------------------------------

WP2LProcessDecBlueAndRedFunc WP2LAddGreenToBlueAndRed;
WP2LPredictorAddSubFunc WP2LPredictorsAdd[kNumPredictors];
WP2LPredictorFunc WP2LPredictors[kNumPredictors];

// exposed plain-C implementations
WP2LPredictorAddSubFunc WP2LPredictorsAdd_C[kNumPredictors];
WP2LPredictorFunc WP2LPredictors_C[kNumPredictors];

WP2LTransformColorInverseFunc WP2LTransformColorInverse;

WP2LMapARGBFunc WP2LMapColor;

static volatile WP2CPUInfo lossless_last_cpuinfo_used =
    (WP2CPUInfo)&lossless_last_cpuinfo_used;

#define COPY_PREDICTOR_ARRAY(IN, OUT) do {                \
  (OUT)[0] = IN##0_C;                                     \
  (OUT)[1] = IN##1_C;                                     \
  (OUT)[2] = IN##2_C;                                     \
  (OUT)[3] = IN##3_C;                                     \
  (OUT)[4] = IN##4_C;                                     \
  (OUT)[5] = IN##5_C;                                     \
  (OUT)[6] = IN##6_C;                                     \
  (OUT)[7] = IN##7_C;                                     \
  (OUT)[8] = IN##8_C;                                     \
  (OUT)[9] = IN##9_C;                                     \
  (OUT)[10] = IN##10_C;                                   \
  (OUT)[11] = IN##11_C;                                   \
  (OUT)[12] = IN##12_C;                                   \
  (OUT)[13] = IN##13_C;                                   \
  (OUT)[14] = IN##0_C; /* <- padding security sentinels*/ \
  (OUT)[15] = IN##0_C;                                    \
} while (0);

WP2_TSAN_IGNORE_FUNCTION void WP2LDspInit() {
  if (lossless_last_cpuinfo_used == WP2GetCPUInfo) return;

  WP2MathInit();

  COPY_PREDICTOR_ARRAY(Predictor, WP2LPredictors)
  COPY_PREDICTOR_ARRAY(Predictor, WP2LPredictors_C)
  COPY_PREDICTOR_ARRAY(PredictorAdd, WP2LPredictorsAdd)
  COPY_PREDICTOR_ARRAY(PredictorAdd, WP2LPredictorsAdd_C)

  WP2LAddGreenToBlueAndRed = WP2LAddGreenToBlueAndRed_C;

  WP2LTransformColorInverse = WP2LTransformColorInverse_C;

  WP2LMapColor = MapARGB_C;

  // If defined, use CPUInfo() to overwrite some pointers with faster versions.
  if (WP2GetCPUInfo != NULL) {
  // TODO(skal): SSE2, etc.
  }
  lossless_last_cpuinfo_used = WP2GetCPUInfo;
}
#undef COPY_PREDICTOR_ARRAY

//------------------------------------------------------------------------------
