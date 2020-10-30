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
// Intra predictions
//
// Author: Skal (pascal.massimino@gmail.com)

#include <algorithm>
#include <cmath>
#include <numeric>

#include "src/dsp/dsp.h"
#include "src/dsp/math.h"
#include "src/utils/utils.h"

namespace {

//------------------------------------------------------------------------------

inline void FillBlock(int16_t value, uint32_t bw, uint32_t bh,
                      int16_t* dst, size_t step) {
  for (uint32_t j = 0; j < bh; ++j) {
    for (uint32_t i = 0; i < bw; ++i) dst[i] = value;
    dst += step;
  }
}

void DC_C(const int16_t* ctx, uint32_t bw, uint32_t bh, int16_t min_value,
          int16_t max_value, int16_t* dst, size_t step) {
  // Only smooth over the left, top-left and top contexts.
  const int32_t size = bh + 1 + bw;
  const int32_t sum = std::accumulate(ctx, ctx + size, 0);
  FillBlock(WP2::DivRound(sum, size), bw, bh, dst, step);
}

void DC_L_C(const int16_t* ctx, uint32_t bw, uint32_t bh, int16_t min_value,
            int16_t max_value, int16_t* dst, size_t step) {
  const int32_t sum = std::accumulate(ctx, ctx + bh, 0);
  FillBlock(WP2::DivRound(sum, (int32_t)bh), bw, bh, dst, step);
}

void DC_T_C(const int16_t* ctx, uint32_t bw, uint32_t bh, int16_t min_value,
            int16_t max_value, int16_t* dst, size_t step) {
  const int32_t sum = std::accumulate(ctx + bh + 1, ctx + bh + 1 + bw, 0);
  FillBlock(WP2::DivRound(sum, (int32_t)bw), bw, bh, dst, step);
}

//------------------------------------------------------------------------------

// The Sm_Weights_Tx_* in AV1.
static constexpr int16_t kSmoothWeights[] = {
  0, 0, 0, 0,   // <- dummy entries to make the offsets easy
  255, 149, 85, 64,  // 4
  255, 197, 146, 105, 73, 50, 37, 32,  // 8
  255, 225, 196, 170, 145, 123, 102, 84, 68, 54, 43, 33, 26, 20, 17, 16,  // 16
  255, 240, 225, 210, 196, 182, 169, 157, 145, 133, 122, 111, 101, 92, 83, 74,
  66, 59, 52, 45, 39, 34, 29, 25, 21, 17, 14, 12, 10, 9, 8, 8,  // 32
};
STATIC_ASSERT_ARRAY_SIZE(kSmoothWeights, 4 + 4 + 8 + 16 + 32);
static constexpr int32_t kSmoothShift = 8;

static const int16_t* GetWeights(uint32_t d) { return &kSmoothWeights[d]; }

void Smooth_C(const int16_t* ctx, uint32_t bw, uint32_t bh, int16_t min_value,
              int16_t max_value, int16_t* dst, size_t step) {
  const int16_t* const wLR = GetWeights(bw);  // left / right weights
  const int16_t* const wAB = GetWeights(bh);  // above / below weights
  const int16_t below = ctx[0];
  for (uint32_t j = 0; j < bh; ++j) {
    const int16_t left = ctx[bh - 1 - j];
    const int16_t right = ctx[bh + 1 + bw + 1 + j];
    for (uint32_t i = 0; i < bw; ++i) {
      const int16_t above = ctx[bh + 1 + i];
      const int32_t v = above * wAB[j] + below * (256 - wAB[j])
                      + left * wLR[i] + right * (256 - wLR[i]);
      dst[i] = WP2::RightShiftRound(v, kSmoothShift + 1);
    }
    dst += step;
  }
}

void Smooth_H_C(const int16_t* ctx, uint32_t bw, uint32_t bh, int16_t min_value,
                int16_t max_value, int16_t* dst, size_t step) {
  const int16_t* const wLR = GetWeights(bw);  // left / right weights
  for (uint32_t j = 0; j < bh; ++j) {
    const int16_t left = ctx[bh - 1 - j];
    const int16_t right = ctx[bh + 1 + bw + 1 + j];
    for (uint32_t i = 0; i < bw; ++i) {
      const int32_t v = left * wLR[i] + right * (256 - wLR[i]);
      dst[i] = WP2::RightShiftRound(v, kSmoothShift);
    }
    dst += step;
  }
}

void Smooth_V_C(const int16_t* ctx, uint32_t bw, uint32_t bh, int16_t min_value,
                int16_t max_value, int16_t* dst, size_t step) {
  const int16_t* const wAB = GetWeights(bh);  // above / below weights
  const int16_t below = ctx[0];
  for (uint32_t j = 0; j < bh; ++j) {
    for (uint32_t i = 0; i < bw; ++i) {
      const int16_t above = ctx[bh + 1 + i];
      const int32_t v = above * wAB[j] + below * (256 - wAB[j]);
      dst[i] = WP2::RightShiftRound(v, kSmoothShift);
    }
    dst += step;
  }
}

//------------------------------------------------------------------------------

void TrueMotion_C(const int16_t* ctx, uint32_t bw, uint32_t bh,
                  int16_t min_value, int16_t max_value, int16_t* dst,
                  size_t step) {
  const int16_t* const top = ctx + bh + 1;
  const int16_t top_left = ctx[bh];
  for (uint32_t y = 0; y < bh; ++y) {
    const int16_t base = ctx[bh - 1 - y] - top_left;
    for (uint32_t x = 0; x < bw; ++x) {
      dst[x] = WP2::Clamp<int16_t>(top[x] + base, min_value, max_value);
    }
    dst += step;
  }
}

//------------------------------------------------------------------------------

static float FilteredProjWeight(float x0, float y0, float c, float s, float x,
                                float y) {
  const float w = WP2::ProjWeight(x0, y0, c, s, x, y);
  constexpr float kMinWeight = 0.01f;  // "too-far" threshold
  return (w < kMinWeight) ? kMinWeight : w;
}

static int16_t WeightedSum(float weights[], const int16_t ctx[], uint32_t size,
                           int16_t min_value, int16_t max_value) {
  float sum = 0.f, v = 0.f;
  for (uint32_t m = 0; m < size; ++m) {
    v += ctx[m] * weights[m];
    sum += weights[m];
  }
  assert(sum > 0.f);
  return WP2::Clamp((int16_t)(v / sum), min_value, max_value);
}

}  // namespace

namespace WP2 {

//------------------------------------------------------------------------------
// Paeth

static inline int16_t PaethPredictOnePixel(int16_t left, int16_t top,
                                           int16_t top_left) {
  // abs((top + left - top_left) - left)
  const uint16_t p_left = std::abs(top - top_left);
  // abs((top + left - top_left) - top)
  const uint16_t p_top = std::abs(left - top_left);
  // abs((top + left - top_left) - top_left)
  const uint16_t p_top_left = std::abs((top - top_left) + (left - top_left));
  // Return nearest to base of left, top and top_left.
  return (p_left <= p_top && p_left <= p_top_left)
             ? left
             : (p_top <= p_top_left) ? top : top_left;
}

void BasePaethPredictor(const int16_t* ctx, uint32_t bw, uint32_t bh,
                        int16_t min_value, int16_t max_value, int16_t* dst,
                        size_t step) {
  const int16_t top_left = ctx[bh];
  for (uint32_t y = 0; y < bh; ++y) {
    const int16_t left = ctx[bh - 1 - y];
    for (uint32_t x = 0; x < bw; ++x) {
      const int16_t top = ctx[bh + 1 + x];
      dst[x] = PaethPredictOnePixel(left, top, top_left);
    }
    dst += step;
  }
}

//------------------------------------------------------------------------------

float ProjWeight(float x0, float y0, float c, float s, float x, float y) {
  // The axes we use here are: Y going down, X going right.
  // But the angle is defined using the VP9/AV1 convention: Y going up, X going
  // right.
  // Compute distance from block point to context point, projected onto
  // the angle axis of vector (-c, s).
  const float d = -(x - x0) * c + (y - y0) * s;
  if (d > 0.f) return 0.f;  // not in the correct direction
  // Compute vector from angle axis to context point.
  const float rx = (x - x0) + d * c;
  const float ry = (y - y0) - d * s;
  const float w =
      1.f / (1.f + 0.2f * std::pow(rx * rx + ry * ry, 2.f));  // NOLINT
  return w;
}

void BaseAnglePredictor(float angle_deg, const int16_t ctx[], uint32_t bw,
                        uint32_t bh, int16_t* dst, size_t step,
                        int16_t min_value, int16_t max_value) {
  const float alpha = M_PI / 180.f * angle_deg;
  const float c = std::cos(alpha);
  const float s = std::sin(alpha);
  for (uint32_t y = 0; y < bh; ++y) {
    for (uint32_t x = 0; x < bw; ++x) {
      float weights[kMaxContextSize];
      uint32_t m = 0;
      // left
      for (uint32_t j = 0; j < bh; ++j) {
        weights[m++] = FilteredProjWeight(x, y, c, s, -1, bh - 1 - j);
      }
      // top-left
      weights[m++] = FilteredProjWeight(x, y, c, s, -1, -1);
      // top
      for (uint32_t i = 0; i < bw; ++i) {
        weights[m++] = FilteredProjWeight(x, y, c, s, i, -1);
      }
      // top-right
      weights[m++] = FilteredProjWeight(x, y, c, s, bw, -1);
      // right
      for (uint32_t j = 0; j < bh; ++j) {
        weights[m++] = FilteredProjWeight(x, y, c, s, bw, j);
      }
      dst[x] = WeightedSum(weights, ctx, m, min_value, max_value);
    }
    dst += step;
  }
}

static inline void AssertInt32_TBounds(int64_t i) {
  assert(i <= std::numeric_limits<int32_t>::max());
  assert(i >= std::numeric_limits<int32_t>::min());
}

// Fixed-point precision in bits, maximized to have all asserts pass.
constexpr uint32_t kAnglePredPrecision = 15;
static_assert(kAnglePredPrecision + 1 + 10 + 1 <= 32,
              "Too many precision bits");
constexpr int32_t kOne = (1u << kAnglePredPrecision);
constexpr uint32_t kMask = kOne - 1;
// (int32_t)(1. / std::tan(alpha) * (1u << kAnglePredPrecision))
static constexpr int32_t kCoTanTable[] = {
  46182, 41089, 36667, 32768, 29283, 26131, 23250, 20589,
  18110, 15780, 13572, 11466, 9440, 7479, 5567, 3692, 1840,
  0,  // 90 degree
  -1840, -3692, -5567, -7479, -9440, -11466,  -13572, -15780, -18110, -20589,
  -23250, -26131, -29283, -32767, -36667, -41089, -46182, -52149, -59289,
  -68043, -79108, -93645, -113740, -143565, -192858, -290824, -583488,
  -2147483648,  // 180 degree
  583488, 290824, 192858, 143565, 113740, 93645, 79108, 68043, 59289, 52149
};
STATIC_ASSERT_ARRAY_SIZE(kCoTanTable, kNumDirectionalAngles);

// (int32_t)(std::tan(alpha) * (1u << kAnglePredPrecision))
static constexpr int32_t kTanTable[] = {
  // (these 'angle < 90' entries are not actually used)
  23250, 26131, 29283, 32767, 36667, 41089, 46182, 52149,
  59289, 68043, 79108, 93645, 113740, 143565, 192858, 290824, 583488,
  // 90 degrees and up...
  -2147483648, -583488, -290824, -192858, -143565, -113740, -93645,
  -79108, -68043, -59289, -52149, -46182, -41089, -36667, -32768,
  -29283, -26131, -23250, -20589, -18110, -15780, -13572, -11466,
  -9440, -7479, -5567, -3692, -1840, 0, 1840, 3692,
  5567, 7479, 9440, 11466, 13572, 15780, 18110, 20589,
};
STATIC_ASSERT_ARRAY_SIZE(kTanTable, kNumDirectionalAngles);

static void AnglePredInterpolate_C(const int16_t* src,
                                   int32_t frac, int16_t* dst, uint32_t len) {
  assert(frac >= 0 && frac <= kOne);
  uint32_t x = 0;
  for (; x < len; ++x) {
    // Computes src[x] * (kOne - frac) + src[x + 1] * frac
    dst[x] = src[x] + (((src[x + 1] - src[x]) * frac) >> kAnglePredPrecision);
  }
}

void SimpleAnglePredictor(uint8_t angle_idx, const int16_t ctx[],
                          uint32_t log2_bw, uint32_t log2_bh,
                          int16_t* dst, size_t step) {
  const uint32_t bw = 1u << log2_bw;
  const uint32_t bh = 1u << log2_bh;
  const uint32_t context_size =
      (angle_idx < kAngle_90) ? ContextWithTRSize(bw, bh) : ContextSize(bw, bh);
  if (angle_idx == kAngle_45) {
    for (uint32_t y = 0; y < bh; ++y, dst += step) {
      for (uint32_t x = 0; x < bw; ++x) dst[x] = ctx[bh + 1 + y + 1 + x];
    }
  } else if (angle_idx == kAngle_90) {
    for (uint32_t y = 0; y < bh; ++y, dst += step) {
      for (uint32_t x = 0; x < bw; ++x) dst[x] = ctx[bh + 1 + x];
    }
  } else if (angle_idx == kAngle_135) {
    for (uint32_t y = 0; y < bh; ++y, dst += step) {
      for (uint32_t x = 0; x < bw; ++x) dst[x] = ctx[bh - y + x];
    }
  } else if (angle_idx == kAngle_180) {
    for (uint32_t y = 0; y < bh; ++y, dst += step) {
      for (uint32_t x = 0; x < bw; ++x) dst[x] = ctx[bh - 1 - y];
    }
  } else {
    const int32_t cot = kCoTanTable[angle_idx];
    // Line equation in X and Y of slope 'angle' going through pixel (x, y):
    // X*sin(-angle)-Y*cos(-angle)=x*sin(-angle)-y*cos(-angle) x axis going
    // right, y axis going down (hence the angle opposition).
    // Let's simplify it to:
    // X*sin(angle)+Y*cos(angle)=x*sin(angle)+y*cos(angle)
    if (angle_idx < kAngle_90) {
      for (uint32_t y = 0; y < bh; ++y, dst += step) {
        // Hit the top context. The first index X at Y = -1, for x = 0 is:
        assert((uint64_t)(y + 1) * cot <= std::numeric_limits<uint32_t>::max());
        const uint32_t X = (y + 1) * (uint32_t)cot;
        // Round down to get the first index.
        uint32_t ind = bh + 1 + (X >> kAnglePredPrecision);
        if (ind >= context_size) {
          // If we reach outside the context, repeat the last element of the
          // line above.
          assert(y > 0);
          std::fill(dst, dst + bw, dst[-step + bw - 1]);
        } else if (ind + 1 >= context_size) {
          // If we cannot interpolate with the next element because it is out
          // of context, repeat the last context element.
          std::fill(dst, dst + bw, ctx[ind]);
        } else {
          const int32_t decimal = (X & kMask);
          assert(context_size - ind - 1 >= bw);
          AnglePredInterpolate(ctx + ind, decimal, dst, bw);
        }
      }
    } else {        // angle >= 90.f
      const int32_t t = kTanTable[angle_idx];
      const int32_t bh_minus_1 =
          ChangePrecision((int32_t)bh - 1, 0, kAnglePredPrecision);
      AssertInt32_TBounds(t);
      for (int32_t y = 0; y < (int32_t)bh; ++y, dst += step) {
        int32_t x_max;
        // Figure out the x that hits a boundary of the left context.
        if (angle_idx > kAngle_180) {
          // Hits the left context for X = -1 and Y = bh - 1;
          x_max = -kOne + (bh - 1) * cot - y * cot;
        } else {
          // Hits the left context for X = -1 and Y = -1;
          x_max = -kOne + (-1) * cot - y * cot;
        }
        assert(x_max >= -kOne);
        // Round down.
        x_max = (x_max < 0) ? -1 : (x_max >> kAnglePredPrecision);
        if (x_max >= (int32_t)bw) x_max = bw - 1;

        int32_t x = 0;
        // Hit the left context. X = -1.
        for (; x <= x_max; ++x) {
          AssertInt32_TBounds((int64_t)t + (x * t + y * kOne));
          const int32_t Y = t + (x * t + y * kOne);
          assert(Y >= -kOne && Y <= bh_minus_1);
          uint32_t ind = bh_minus_1 - Y;
          const int32_t decimal = (ind & kMask);
          ind >>= kAnglePredPrecision;
          assert(ind < context_size);
          // Compute ctx[ind] * (kOne - decimal) + ctx[ind + 1] * decimal.
          dst[x] =
              ctx[ind] + ChangePrecision((ctx[ind + 1] - ctx[ind]) * decimal,
                                         kAnglePredPrecision, 0u);
        }
        if (x == (int32_t)bw) continue;
        if (angle_idx > kAngle_180) {
          // Repeat the bottom of the left context.
          std::fill(dst + x, dst + bw, ctx[0]);
        } else {
          // Hit the top context. Y = -1.
          AssertInt32_TBounds(((int64_t)x * kOne + y * cot) - (-1) * cot);
          const int32_t X = (x * kOne + y * cot) - (-1) * cot;
          assert(X >= -kOne && X <= 0);
          // Add kOne to have a positive value.
          const uint32_t X_plus_1 = X + kOne;
          const int32_t decimal = (X_plus_1 & kMask);
          const uint32_t ind = bh + (X_plus_1 >> kAnglePredPrecision);
          AnglePredInterpolate(ctx + ind, decimal, dst + x, bw - x);
        }
      }
    }
  }
}

//------------------------------------------------------------------------------

// Weighted distance between predicted pixel (x0, y0) and context pixel (x, y).
// The larger 'strength' is, the more faster the weight decreases with distance.
static constexpr float kMinWeight = 0.01f;     // "too-far" threshold
static float WDistance(int32_t x0, int32_t y0, int32_t x, int32_t y,
                       float strength) {
  // Threshold below which strength is considered flat.
  const float kStrengthThreshold = 0.01f;
  if (strength <= kStrengthThreshold) return 1.f;
  const int32_t d2 = (x - x0) * (x - x0) + (y - y0) * (y - y0);
  if (d2 == 0) return 0.;
  const float weight = std::pow(d2, -6.f * strength);
  return (weight < kMinWeight) ? kMinWeight : weight;
}

// Table-based version
void PrecomputeLargeWeightTable(float strength, LargeWeightTable table) {
  for (uint32_t y = 0; y < LargeWeightTableDim; ++y) {
    for (uint32_t x = 0; x < LargeWeightTableDim; ++x) {
      table[x + y * LargeWeightTableDim] = WDistance(0, 0, x, y, strength);
    }
  }
}

static inline float WDistance(int32_t x0, int32_t y0, int32_t x, int32_t y,
                              int32_t w, const LargeWeightTable table) {
  static constexpr int32_t limit = LargeWeightTableDim - 1u;
  const int32_t dx = std::min(std::abs(x - x0), limit);
  const int32_t dy = std::min(std::abs(y - y0), limit);
  const float weight = table[dx + LargeWeightTableDim * dy];
  if (weight == kMinWeight) {
    const int32_t min_distance = std::min({x0 + 1, w - x0, y0 + 1});
    // If the weight is low but this is the closest we'll ever be anyway,
    // set a weight of 1, so the pixel is the average of its closest neighbors.
    // Otherwise, ignore this context pixel completely.
    return (dx == min_distance || dy == min_distance) ? 1.f : kMinWeight;
  }
  return weight;
}

uint32_t ComputeFuseWeights(uint32_t w, uint32_t h, int32_t x, int32_t y,
                            const LargeWeightTable table,
                            float weights[kMaxContextSize]) {
  uint32_t m = 0;
  // left
  for (uint32_t j = 0; j < h; ++j) {
    weights[m++] = WDistance(x, y, -1, h - 1 - j, w, table);
  }
  // top-left
  weights[m++] = WDistance(x, y, -1, -1, w, table);
  // top
  for (uint32_t i = 0; i < w; ++i) {
    weights[m++] = WDistance(x, y, i, -1, w, table);
  }
  // top-right
  weights[m++] = WDistance(x, y, w, -1, w, table);
  // right
  for (uint32_t j = 0; j < h; ++j) {
    weights[m++] = WDistance(x, y, w, j, w, table);
  }
  assert(m <= kMaxContextSize);
  return m;
}

void BaseFusePredictor(const LargeWeightTable table, const int16_t ctx[],
                       uint32_t bw, uint32_t bh, int16_t* dst, size_t step,
                       int16_t min_value, int16_t max_value) {
  float weights[kMaxContextSize];
  for (uint32_t idx = 0, y = 0; y < bh; ++y) {
    for (uint32_t x = 0; x < bw; ++x, ++idx) {
      const uint32_t num_weights =
          ComputeFuseWeights(bw, bh, x, y, table, weights);
      dst[x] = WeightedSum(weights, ctx, num_weights, min_value, max_value);
    }
    dst += step;
  }
}

//------------------------------------------------------------------------------
// Add / Sub predictions

static void AddRow_C(const int16_t src[], const int32_t res[],
                     int32_t min, int32_t max, int16_t dst[], uint32_t len) {
  assert((len & 3) == 0);
  for (uint32_t x = 0; x < len; ++x) dst[x] = Clamp(res[x] + src[x], min, max);
}

static void SubtractRow_C(const int16_t src[], const int16_t pred[],
                          int32_t dst[], uint32_t len) {
  assert((len & 3) == 0);
  for (uint32_t x = 0; x < len; ++x) dst[x] = src[x] - pred[x];
}

#define ADD_SUB_BLOCK_FUNCS_DCL(WIDTH, ADD_FUNC, SUB_FUNC, EXT)                \
static void AddBlock_## WIDTH ## EXT(const int16_t src[], uint32_t src_step,   \
                                     const int32_t res[], uint32_t res_step,   \
                                     int32_t min, int32_t max,                 \
                                     int16_t dst[], uint32_t dst_step,         \
                                     uint32_t height) {                        \
  for (uint32_t y = 0; y < height; ++y) {                                      \
    ADD_FUNC(src, res, min, max, dst, WIDTH);                                  \
    src += src_step;                                                           \
    res += res_step;                                                           \
    dst += dst_step;                                                           \
  }                                                                            \
}                                                                              \
static void SubBlock_## WIDTH ## EXT(const int16_t src[], uint32_t src_step,   \
                                     const int16_t pred[], uint32_t pred_step, \
                                     int32_t dst[], uint32_t dst_step,         \
                                     uint32_t height) {                        \
  for (uint32_t y = 0; y < height; ++y) {                                      \
    SUB_FUNC(src, pred, dst, WIDTH);                                           \
    src += src_step;                                                           \
    pred += pred_step;                                                         \
    dst += dst_step;                                                           \
  }                                                                            \
}

ADD_SUB_BLOCK_FUNCS_DCL( 4, AddRow_C, SubtractRow_C, _C)
ADD_SUB_BLOCK_FUNCS_DCL( 8, AddRow_C, SubtractRow_C, _C)
ADD_SUB_BLOCK_FUNCS_DCL(16, AddRow_C, SubtractRow_C, _C)
ADD_SUB_BLOCK_FUNCS_DCL(32, AddRow_C, SubtractRow_C, _C)

//------------------------------------------------------------------------------
// SSE4.1 implementation

#if defined(WP2_USE_SSE)

static void AnglePredInterpolate_SSE(const int16_t* src,
                                     int32_t frac, int16_t* dst, uint32_t len) {
  assert(frac >= 0 && frac <= kOne);
  static_assert(kAnglePredPrecision == 15, "loop is tuned for 15bit!");
  uint32_t x = 0;
  const __m128i M = _mm_set1_epi16(frac);
  for (; x + 8 <= len; x += 8) {
    const __m128i A = _mm_loadu_si128((const __m128i*)(src + x + 0));
    const __m128i B = _mm_loadu_si128((const __m128i*)(src + x + 1));
    const __m128i C = _mm_sub_epi16(B, A);
    // _mm_mulhi_epi16 will perform the >>16 automatically. Since we want
    // diff >> 15, we pre-multiply diff by 2: (2*diff)>>16 = diff>>15
    const __m128i D = _mm_add_epi16(C, C);   // diff *= 2
    const __m128i E = _mm_mulhi_epi16(D, M);
    const __m128i F = _mm_add_epi16(A, E);
    _mm_storeu_si128((__m128i*)(dst + x), F);
  }
  if (x < len) AnglePredInterpolate_C(src + x, frac, dst + x, len - x);
}

static void AddRow_SSE(const int16_t src[], const int32_t res[],
                       int32_t min, int32_t max, int16_t dst[], uint32_t len) {
  assert((len & 3) == 0);
  const __m128i m_max = _mm_set1_epi16(max);
  const __m128i m_min = _mm_set1_epi16(min);
  uint32_t x = 0;
  for (; x + 8 <= len; x += 8) {
    const __m128i A0 = _mm_loadu_si128((const __m128i*)(res + x + 0));
    const __m128i A1 = _mm_loadu_si128((const __m128i*)(res + x + 4));
    const __m128i B = _mm_loadu_si128((const __m128i*)(src + x));
    const __m128i C = _mm_packs_epi32(A0, A1);
    const __m128i D = _mm_add_epi16(C, B);
    const __m128i E = _mm_min_epi16(D, m_max);
    const __m128i F = _mm_max_epi16(E, m_min);
    _mm_storeu_si128((__m128i*)(dst + x), F);
  }
  if (x < len) {
    const __m128i A = _mm_loadu_si128((const __m128i*)(res + x));
    const __m128i B = _mm_loadl_epi64((const __m128i*)(src + x));
    const __m128i C = _mm_add_epi16(_mm_packs_epi32(A, A), B);
    const __m128i D = _mm_min_epi16(C, m_max);
    const __m128i E = _mm_max_epi16(D, m_min);
    _mm_storel_epi64((__m128i*)(dst + x), E);
  }
}

static void SubtractRow_SSE(const int16_t src[], const int16_t pred[],
                            int32_t dst[], uint32_t len) {
  assert((len & 3) == 0);
  uint32_t x = 0;
  for (; x + 8 <= len; x += 8) {
    const __m128i A = _mm_loadu_si128((const __m128i*)(src + x));
    const __m128i B = _mm_loadu_si128((const __m128i*)(pred + x));
    const __m128i C = _mm_sub_epi16(A, B);
    const __m128i D0 = _mm_cvtepi16_epi32(C);
    const __m128i D1 = _mm_cvtepi16_epi32(_mm_srli_si128(C, 8));
    _mm_storeu_si128((__m128i*)(dst + x + 0), D0);
    _mm_storeu_si128((__m128i*)(dst + x + 4), D1);
  }
  if (x < len) {   // left over 4 pixels
    const __m128i A = _mm_loadl_epi64((const __m128i*)(src + x));
    const __m128i B = _mm_loadl_epi64((const __m128i*)(pred + x));
    const __m128i C = _mm_sub_epi16(A, B);
    const __m128i D = _mm_cvtepi16_epi32(C);
    _mm_storeu_si128((__m128i*)(dst + x), D);
  }
}

ADD_SUB_BLOCK_FUNCS_DCL( 4, AddRow_SSE, SubtractRow_SSE, _SSE)
ADD_SUB_BLOCK_FUNCS_DCL( 8, AddRow_SSE, SubtractRow_SSE, _SSE)
ADD_SUB_BLOCK_FUNCS_DCL(16, AddRow_SSE, SubtractRow_SSE, _SSE)
ADD_SUB_BLOCK_FUNCS_DCL(32, AddRow_SSE, SubtractRow_SSE, _SSE)

static WP2_TSAN_IGNORE_FUNCTION void PredictionInitSSE() {
  AnglePredInterpolate = AnglePredInterpolate_SSE;
  SubtractRow = SubtractRow_SSE;
  AddRow = AddRow_SSE;

  SubtractBlock[0] = SubBlock_4_SSE;
  SubtractBlock[1] = SubBlock_8_SSE;
  SubtractBlock[2] = SubBlock_16_SSE;
  SubtractBlock[3] = SubBlock_32_SSE;
  AddBlock[0] = AddBlock_4_SSE;
  AddBlock[1] = AddBlock_8_SSE;
  AddBlock[2] = AddBlock_16_SSE;
  AddBlock[3] = AddBlock_32_SSE;
}

#endif  // WP2_USE_SSE

//------------------------------------------------------------------------------
#undef ADD_SUB_BLOCK_FUNCS_DCL

LPredictF BasePredictors[BPRED_LAST];

void (*AnglePredInterpolate)(const int16_t* src, int32_t frac,
                             int16_t* dst, uint32_t len) = nullptr;

void (*SubtractRow)(const int16_t src[], const int16_t pred[],
                    int32_t dst[], uint32_t len) = nullptr;

void (*AddRow)(const int16_t src[], const int32_t res[],
               int32_t min, int32_t max,
               int16_t dst[], uint32_t len) = nullptr;

SubtractBlockFunc SubtractBlock[4] = { nullptr };
AddBlockFunc AddBlock[4] = { nullptr };

static volatile WP2CPUInfo intra_last_cpuinfo_used =
    (WP2CPUInfo)&intra_last_cpuinfo_used;

WP2_TSAN_IGNORE_FUNCTION void PredictionInit() {
  if (intra_last_cpuinfo_used == WP2GetCPUInfo) return;

  BasePredictors[BPRED_DC] = DC_C;
  BasePredictors[BPRED_DC_L] = DC_L_C;
  BasePredictors[BPRED_DC_T] = DC_T_C;
  BasePredictors[BPRED_SMOOTH] = Smooth_C;
  BasePredictors[BPRED_SMOOTH_H] = Smooth_H_C;
  BasePredictors[BPRED_SMOOTH_V] = Smooth_V_C;
  BasePredictors[BPRED_TM] = TrueMotion_C;

  AnglePredInterpolate = AnglePredInterpolate_C;

  SubtractRow = SubtractRow_C;
  AddRow = AddRow_C;

  SubtractBlock[0] = SubBlock_4_C;
  SubtractBlock[1] = SubBlock_8_C;
  SubtractBlock[2] = SubBlock_16_C;
  SubtractBlock[3] = SubBlock_32_C;
  AddBlock[0] = AddBlock_4_C;
  AddBlock[1] = AddBlock_8_C;
  AddBlock[2] = AddBlock_16_C;
  AddBlock[3] = AddBlock_32_C;

  if (WP2GetCPUInfo != nullptr) {
#if defined(WP2_USE_SSE)
    if (WP2GetCPUInfo(kSSE)) PredictionInitSSE();
#endif
  }
  intra_last_cpuinfo_used = WP2GetCPUInfo;
}

}  // namespace WP2

//------------------------------------------------------------------------------
