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

#include <limits>

#include "src/common/color_precision.h"
#include "src/dsp/lossless/lossless.h"
#include "src/dsp/lossless/lossless_common.h"
#include "src/dsp/math.h"
#include "src/enc/lossless/losslessi_enc.h"

using WP2::Vector_u32;

static const float kSpatialPredictorBias = 15.f;
static const int kPredLowEffort = 11;

//------------------------------------------------------------------------------
// Methods to calculate Entropy (Shannon).

static float PredictionCostSpatial(const Vector_u32& counts, int weight_0,
                                   double exp_val) {
  // TODO(vrabaud) make this constant be dependent on counts.size().
  //               This will matter for 10-bit.
  const uint32_t counts_size = 256;
  const uint32_t significant_symbols = counts_size >> 4;
  const double exp_decay_factor = 0.6;
  double bits = weight_0 * counts[0];
  for (uint32_t i = 1; i < significant_symbols; ++i) {
    // TODO(vrabaud) once subsequent usage of argb is signed, use MakeIndex
    //               to compute the cost.
    bits += exp_val * (counts[i] + counts[counts_size - i]);
    exp_val *= exp_decay_factor;
  }
  return (float)(-0.1 * bits);
}

static float PredictionCostSpatialHistogram(const Vector_u32 accumulated[4],
                                            const Vector_u32 tile[4]) {
  int i;
  double retval = 0;
  for (i = 0; i < 4; ++i) {
    const double kExpValue = 0.94;
    retval += PredictionCostSpatial(tile[i], 1, kExpValue);
    retval += WP2LCombinedShannonEntropy(tile[i].data(), accumulated[i].data(),
                                         tile[i].size());
  }
  return (float)retval;
}

static inline void UpdateHisto(Vector_u32 histo_argb[4],
                               const uint16_t* const argb) {
  for (size_t i = 0; i < 4; ++i) ++histo_argb[i][argb[i]];
}

//------------------------------------------------------------------------------
// Spatial transform functions.

static inline void PredictBatch(int mode, int x_start, int y,
                                uint32_t num_pixels, uint32_t channel_bits,
                                const uint16_t* const current, bool has_alpha,
                                const uint16_t* const upper, uint16_t* out) {
  if (x_start == 0) {
    if (y == 0) {
      // ARGB_BLACK.
      WP2LPredictorsSub[0](current, has_alpha, nullptr, 1, channel_bits, out);
    } else {
      // Top one.
      WP2LPredictorsSub[2](current, has_alpha, upper, 1, channel_bits, out);
    }
    ++x_start;
    out += 4;
    --num_pixels;
  }
  if (y == 0) {
    // Left one.
    WP2LPredictorsSub[1](&current[4 * x_start], has_alpha, nullptr, num_pixels,
                         channel_bits, out);
  } else {
    WP2LPredictorsSub[mode](&current[4 * x_start], has_alpha,
                            &upper[4 * x_start], num_pixels, channel_bits, out);
  }
}

namespace WP2L {
// Stores the difference between the pixel and its prediction in "out".
// In case of a lossy encoding, updates the source image to avoid propagating
// the deviation further to pixels which depend on the current pixel for their
// predictions.
static inline void GetResidual(int width, uint32_t channel_bits,
                               uint16_t* const upper_row,
                               uint16_t* const current_row, bool has_alpha,
                               int mode, int x_start, int x_end, int y,
                               uint16_t* const out) {
  const uint16_t mask = (1 << channel_bits) - 1;
  const WP2LPredictorFunc pred_func = WP2LPredictors[mode];
  int x;
  for (x = x_start; x < x_end; ++x) {
    uint16_t predict[4], residual[4];
    if (y == 0) {
      if (x == 0) {
        predict[0] = mask;
        predict[1] = 0;
        predict[2] = 0;
        predict[3] = 0;
      } else {
        ColorCopy(&current_row[4 * (x - 1)], predict);  // Left.
      }
    } else if (x == 0) {
      ColorCopy(&upper_row[4 * x], predict);  // Top.
    } else {
      pred_func(&current_row[4 * (x - 1)], &upper_row[4 * x], mask, predict);
    }
    WP2LSubPixels(&current_row[4 * x], has_alpha, predict, mask, residual);
    if (current_row[4 * x] == 0) {
      // If alpha is 0, cleanup RGB. We can choose the RGB values of the
      // residual for best compression. The prediction of alpha itself can be
      // non-zero and must be kept though. We choose RGB of the residual to be
      // 0.
      residual[1] = 0;
      residual[2] = 0;
      residual[3] = 0;
      // Update the source image.
      current_row[4 * x + 0] = has_alpha ? 0 : WP2::kAlphaMax;
      current_row[4 * x + 1] = predict[1];
      current_row[4 * x + 2] = predict[2];
      current_row[4 * x + 3] = predict[3];
      // The prediction for the rightmost pixel in a row uses the leftmost
      // pixel in that row as its top-right context pixel. Hence if we change
      // the leftmost pixel of current_row, the corresponding change must be
      // applied to upper_row as well where top-right context is being read
      // from.
      if (x == 0 && y != 0) {
        ColorCopy(current_row, &upper_row[4 * width]);
      }
    }
    ColorCopy(residual, &out[4 * (x - x_start)]);
  }
}

// Returns best predictor and updates the accumulated histogram.
static int GetBestPredictorForTile(int width, int height, int tile_x,
                                   int tile_y, uint32_t bits,
                                   uint16_t* const argb_scratch,
                                   const uint16_t* const argb, bool has_alpha,
                                   const int16_t* const modes,
                                   Vector_u32 accumulated[4],
                                   Vector_u32 tmp_histo[][4]) {
  const int kNumPredModes = 14;
  const int start_x = tile_x << bits;
  const int start_y = tile_y << bits;
  const int tile_size = 1 << bits;
  const int max_y = std::min(tile_size, height - start_y);
  const int max_x = std::min(tile_size, width - start_x);
  // Whether there exist columns just outside the tile.
  const int have_left = (start_x > 0);
  // Position and size of the strip covering the tile and adjacent columns if
  // they exist.
  const int context_start_x = start_x - have_left;
  const int tiles_per_row = WP2LSubSampleSize(width, bits);
  // Prediction modes of the left and above neighbor tiles.
  const int left_mode =
      (tile_x > 0) ? modes[4 * (tile_y * tiles_per_row + tile_x - 1) + 2]
                   : 0xff;
  const int above_mode =
      (tile_y > 0) ? modes[4 * ((tile_y - 1) * tiles_per_row + tile_x) + 2]
                   : 0xff;
  // The width of upper_row and current_row is one pixel larger than image width
  // to allow the top right pixel to point to the leftmost pixel of the next row
  // when at the right edge.
  uint16_t* upper_row = argb_scratch;
  uint16_t* current_row = &upper_row[4 * (width + 1)];
  float best_diff = std::numeric_limits<float>::max();
  int best_mode = 0;
  int mode;
  // Need pointers to be able to swap arrays.
  Vector_u32* histo_argb = &tmp_histo[0][0];
  Vector_u32* best_histo = &tmp_histo[1][0];
  uint16_t residuals[4 * (1 << kTransformBitsMax)];
  assert(bits >= kTransformBitsMin && bits <= kTransformBitsMax);
  assert(max_x <= (1 << kTransformBitsMax));

  // TODO(vrabaud) fix channel_bits.
  const uint32_t channel_bits = 8;
  for (mode = 0; mode < kNumPredModes; ++mode) {
    float cur_diff;
    int relative_y;
    for (uint32_t i = 0; i < 4; ++i) {
      std::fill(histo_argb[i].begin(), histo_argb[i].end(), 0);
    }
    if (start_y > 0) {
      // Read the row above the tile which will become the first upper_row.
      // Include a pixel to the left if it exists; include a pixel to the right
      // in all cases (wrapping to the leftmost pixel of the next row if it does
      // not exist).
      std::copy(&argb[4 * ((start_y - 1) * width + context_start_x)],
                &argb[4 * ((start_y - 1) * width + context_start_x) +
                    4 * (max_x + have_left + 1)],
                &current_row[4 * context_start_x]);
    }
    for (relative_y = 0; relative_y < max_y; ++relative_y) {
      const int y = start_y + relative_y;
      int relative_x;
      std::swap(upper_row, current_row);
      // Read current_row. Include a pixel to the left if it exists; include a
      // pixel to the right in all cases except at the bottom right corner of
      // the image (wrapping to the leftmost pixel of the next row if it does
      // not exist in the current row).
      std::copy(&argb[4 * (y * width + context_start_x)],
                &argb[4 * (y * width + context_start_x) +
                    4 * (max_x + have_left + (y + 1 < height))],
                &current_row[4 * context_start_x]);

      GetResidual(width, channel_bits, upper_row, current_row, has_alpha, mode,
                  start_x, start_x + max_x, y, residuals);
      for (relative_x = 0; relative_x < max_x; ++relative_x) {
        UpdateHisto(histo_argb, &residuals[4 * relative_x]);
      }
    }
    cur_diff = PredictionCostSpatialHistogram(accumulated, histo_argb);
    // Favor keeping the areas locally similar.
    if (mode == left_mode) cur_diff -= kSpatialPredictorBias;
    if (mode == above_mode) cur_diff -= kSpatialPredictorBias;

    if (cur_diff < best_diff) {
      std::swap(histo_argb, best_histo);
      best_diff = cur_diff;
      best_mode = mode;
    }
  }

  for (uint32_t i = 0; i < 4; i++) {
    for (uint32_t j = 0; j < best_histo[i].size(); j++) {
      accumulated[i][j] += best_histo[i][j];
    }
  }

  return best_mode;
}

// Converts pixels of the image to residuals with respect to predictions.
// If max_quantization > 1, applies near lossless processing, quantizing
// residuals to multiples of quantization levels up to max_quantization
// (the actual quantization level depends on smoothness near the given pixel).
static void CopyImageWithPrediction(int width, int height, int bits,
                                    uint32_t channel_bits,
                                    const int16_t* const modes,
                                    uint16_t* const argb_scratch,
                                    uint16_t* const argb, bool has_alpha,
                                    int speed) {
  const int tiles_per_row = WP2LSubSampleSize(width, bits);
  // The width of upper_row and current_row is one pixel larger than image width
  // to allow the top right pixel to point to the leftmost pixel of the next row
  // when at the right edge.
  uint16_t* upper_row = argb_scratch;
  uint16_t* current_row = &upper_row[4 * (width + 1)];
  int y;

  for (y = 0; y < height; ++y) {
    int x;
    std::swap(upper_row, current_row);
    std::copy(&argb[4 * y * width],
              &argb[4 * y * width + 4 * (width + (y + 1 < height))],
              current_row);

    if (speed == 0) {
      PredictBatch(kPredLowEffort, 0, y, width, channel_bits, current_row,
                   has_alpha, upper_row, &argb[4 * y * width]);
    } else {
      for (x = 0; x < width;) {
        const int mode =
            modes[4 * ((y >> bits) * tiles_per_row + (x >> bits)) + 2];
        int x_end = x + (1 << bits);
        if (x_end > width) x_end = width;
        GetResidual(width, channel_bits, upper_row, current_row, has_alpha,
                    mode, x, x_end, y, &argb[4 * (y * width + x)]);
        x = x_end;
      }
    }
  }
}

// Finds the best predictor for each tile, and converts the image to residuals
// with respect to predictions. If near_lossless_quality < 100, applies
// near lossless processing, shaving off more bits of residuals for lower
// qualities.
WP2Status ResidualImage(uint32_t width, uint32_t height, uint32_t bits,
                        uint32_t channel_bits, int speed, uint16_t* const argb,
                        bool has_alpha, uint16_t* const argb_scratch,
                        int16_t* const image) {
  const uint32_t tiles_per_row = WP2LSubSampleSize(width, bits);
  const uint32_t tiles_per_col = WP2LSubSampleSize(height, bits);
  if (speed == 0) {
    for (size_t i = 0; i < tiles_per_row * tiles_per_col; ++i) {
      image[4 * i + 0] = 0;
      image[4 * i + 1] = 0;
      image[4 * i + 2] = kPredLowEffort;
      image[4 * i + 3] = 0;
    }
  } else {
    const uint32_t histo_size = (1 << channel_bits);
    Vector_u32 histo[4];
    Vector_u32 tmp_histo[2][4];
    for (uint32_t i = 0; i < 4; ++i) {
      WP2_CHECK_ALLOC_OK(histo[i].resize(histo_size));
      std::fill(histo[i].begin(), histo[i].end(), 0);
      for (uint32_t j = 0; j < 2; ++j) {
        WP2_CHECK_ALLOC_OK(tmp_histo[j][i].resize(histo_size));
      }
    }
    for (size_t tile_y = 0; tile_y < tiles_per_col; ++tile_y) {
      for (size_t tile_x = 0; tile_x < tiles_per_row; ++tile_x) {
        const int pred = GetBestPredictorForTile(
            width, height, tile_x, tile_y, bits, argb_scratch, argb, has_alpha,
            image, histo, tmp_histo);
        image[4 * (tile_y * tiles_per_row + tile_x) + 0] = 0;
        image[4 * (tile_y * tiles_per_row + tile_x) + 1] = 0;
        image[4 * (tile_y * tiles_per_row + tile_x) + 2] = pred;
        image[4 * (tile_y * tiles_per_row + tile_x) + 3] = 0;
      }
    }
  }

  CopyImageWithPrediction(width, height, bits, channel_bits, image,
                          argb_scratch, argb, has_alpha, speed);
  return WP2_STATUS_OK;
}

//------------------------------------------------------------------------------
// Color transform functions.

static inline void MultipliersClear(WP2LMultipliers* const m) {
  m->green_to_red = 0;
  m->green_to_blue = 0;
  m->red_to_blue = 0;
}

static inline void ColorCodeToMultipliers(const int16_t* const color,
                                          WP2LMultipliers* const m) {
  m->green_to_red = color[3];
  m->green_to_blue = color[2];
  m->red_to_blue = color[1];
}

static inline void MultipliersToColorCode(const WP2LMultipliers* const m,
                                          int16_t* const res) {
  res[0] = 0;
  res[1] = m->red_to_blue;
  res[2] = m->green_to_blue;
  res[3] = m->green_to_red;
}

static float PredictionCostCrossColor(const Vector_u32& accumulated,
                                      const Vector_u32& counts) {
  // Favor low entropy, locally and globally.
  // Favor small absolute values for PredictionCostSpatial
  static const double kExpValue = 2.4;
  assert(accumulated.size() == counts.size());
  return WP2LCombinedShannonEntropy(counts.data(), accumulated.data(),
                                    counts.size()) +
         PredictionCostSpatial(counts, 3, kExpValue);
}

static float GetPredictionCostCrossColorRed(
    const int16_t* argb, uint32_t width, int tile_width, int tile_height,
    WP2LMultipliers prev_x, WP2LMultipliers prev_y, int16_t green_to_red,
    const Vector_u32& accumulated_red_histo, Vector_u32* const tmp_histo) {
  std::fill(tmp_histo->begin(), tmp_histo->end(), 0);
  float cur_diff;

  WP2LCollectColorRedTransforms(argb, width, tile_width, tile_height,
                                green_to_red, tmp_histo->data());

  cur_diff =
      PredictionCostCrossColor(accumulated_red_histo, *tmp_histo);
  // cur_diff -= 3 favors keeping the areas locally similar
  if (green_to_red == prev_x.green_to_red) {
    cur_diff -= 3;
  }
  if (green_to_red == prev_y.green_to_red) {
    cur_diff -= 3;  // favor keeping the areas locally similar
  }
  if (green_to_red == 0) {
    cur_diff -= 3;
  }
  return cur_diff;
}

static void GetBestGreenToRed(const int16_t* argb, uint32_t width,
                              int tile_width, int tile_height,
                              WP2LMultipliers prev_x, WP2LMultipliers prev_y,
                              int speed,
                              const Vector_u32& accumulated_red_histo,
                              WP2LMultipliers* const best_tx,
                              Vector_u32* const tmp_histo) {
  const int kMaxIters = (speed < 5) ? 4 + speed * 2 / 5 : 6;  // in range [4..6]
  int green_to_red_best = 0;
  int iter, offset;
  float best_diff = GetPredictionCostCrossColorRed(
      argb, width, tile_width, tile_height, prev_x, prev_y, green_to_red_best,
      accumulated_red_histo, tmp_histo);
  for (iter = 0; iter < kMaxIters; ++iter) {
    // ColorTransformDelta is a 3.5 bit fixed point, so 32 is equal to
    // one in color computation. Having initial delta here as 1 is sufficient
    // to explore the range of (-2, 2).
    const int delta = 32 >> iter;
    // Try a negative and a positive delta from the best known value.
    for (offset = -delta; offset <= delta; offset += 2 * delta) {
      const int green_to_red_cur =
          WP2::Clamp(offset + green_to_red_best, -(int32_t)kGreenToRedMax,
                     (int32_t)kGreenToRedMax);
      const float cur_diff = GetPredictionCostCrossColorRed(
          argb, width, tile_width, tile_height, prev_x, prev_y,
          green_to_red_cur, accumulated_red_histo, tmp_histo);
      if (cur_diff < best_diff) {
        best_diff = cur_diff;
        green_to_red_best = green_to_red_cur;
      }
    }
  }
  best_tx->green_to_red = green_to_red_best;
}

static float GetPredictionCostCrossColorBlue(
    const int16_t* argb, int stride, int tile_width, int tile_height,
    WP2LMultipliers prev_x, WP2LMultipliers prev_y, int16_t green_to_blue,
    int16_t red_to_blue, const Vector_u32& accumulated_blue_histo,
    Vector_u32* tmp_histo) {
  std::fill(tmp_histo->begin(), tmp_histo->end(), 0);
  float cur_diff;

  WP2LCollectColorBlueTransforms(argb, stride, tile_width, tile_height,
                                 green_to_blue, red_to_blue, tmp_histo->data());

  cur_diff = PredictionCostCrossColor(accumulated_blue_histo, *tmp_histo);
  // cur_diff -= 3 favors keeping the areas locally similar
  if (green_to_blue == prev_x.green_to_blue) {
    cur_diff -= 3;
  }
  if (green_to_blue == prev_y.green_to_blue) {
    cur_diff -= 3;
  }
  if (red_to_blue == prev_x.red_to_blue) {
    cur_diff -= 3;
  }
  if (red_to_blue == prev_y.red_to_blue) {
    cur_diff -= 3;
  }
  if (green_to_blue == 0) {
    cur_diff -= 3;
  }
  if (red_to_blue == 0) {
    cur_diff -= 3;
  }
  return cur_diff;
}

static constexpr uint32_t kGreenRedToBlueNumAxis = 8;
static constexpr uint32_t kGreenRedToBlueMaxIters = 7;
static void GetBestGreenRedToBlue(const int16_t* argb, uint32_t width,
                                  int tile_width, int tile_height,
                                  WP2LMultipliers prev_x,
                                  WP2LMultipliers prev_y, int speed,
                                  const Vector_u32& accumulated_blue_histo,
                                  WP2LMultipliers* const best_tx,
                                  Vector_u32* tmp_histo) {
  const int8_t offset[kGreenRedToBlueNumAxis][2] =
      {{0, -1}, {0, 1}, {-1, 0}, {1, 0}, {-1, -1}, {-1, 1}, {1, -1}, {1, 1}};
  const int8_t delta_lut[kGreenRedToBlueMaxIters] = { 16, 16, 8, 4, 2, 2, 2 };
  const int iters =
      (speed < 3) ? 1 : (speed >= 5) ? kGreenRedToBlueMaxIters : 4;
  int green_to_blue_best = 0;
  int red_to_blue_best = 0;
  int iter;
  // Initial value at origin:
  float best_diff = GetPredictionCostCrossColorBlue(
      argb, width, tile_width, tile_height, prev_x, prev_y, green_to_blue_best,
      red_to_blue_best, accumulated_blue_histo, tmp_histo);
  for (iter = 0; iter < iters; ++iter) {
    const int delta = delta_lut[iter];
    for (uint32_t axis = 0; axis < kGreenRedToBlueNumAxis; ++axis) {
      const int green_to_blue_cur =
          WP2::Clamp(offset[axis][0] * delta + green_to_blue_best,
                     -(int32_t)kGreenToBlueMax, (int32_t)kGreenToBlueMax);
      const int red_to_blue_cur =
          WP2::Clamp(offset[axis][1] * delta + red_to_blue_best,
                     -(int32_t)kRedToBlueMax, (int32_t)kRedToBlueMax);
      const float cur_diff = GetPredictionCostCrossColorBlue(
          argb, width, tile_width, tile_height, prev_x, prev_y,
          green_to_blue_cur, red_to_blue_cur, accumulated_blue_histo,
          tmp_histo);
      if (cur_diff < best_diff) {
        best_diff = cur_diff;
        green_to_blue_best = green_to_blue_cur;
        red_to_blue_best = red_to_blue_cur;
      }
      if (speed < 3 && iter == 4) {
        // Only axis aligned diffs for lower speed.
        break;  // next iter.
      }
    }
    if (delta == 2 && green_to_blue_best == 0 && red_to_blue_best == 0) {
      // Further iterations would not help.
      break;  // out of iter-loop.
    }
  }
  best_tx->green_to_blue = green_to_blue_best;
  best_tx->red_to_blue = red_to_blue_best;
}

static WP2LMultipliers GetBestColorTransformForTile(
    int tile_x, int tile_y, int bits, WP2LMultipliers prev_x,
    WP2LMultipliers prev_y, int speed, int xsize, int ysize,
    const Vector_u32& accumulated_red_histo,
    const Vector_u32& accumulated_blue_histo, const int16_t* const argb,
    Vector_u32* const tmp_histo) {
  const int max_tile_size = 1 << bits;
  const int tile_y_offset = tile_y * max_tile_size;
  const int tile_x_offset = tile_x * max_tile_size;
  const int all_x_max = std::min(tile_x_offset + max_tile_size, xsize);
  const int all_y_max = std::min(tile_y_offset + max_tile_size, ysize);
  const int tile_width = all_x_max - tile_x_offset;
  const int tile_height = all_y_max - tile_y_offset;
  const int16_t* const tile_argb =
      &argb[4 * (tile_y_offset * xsize + tile_x_offset)];
  WP2LMultipliers best_tx;
  MultipliersClear(&best_tx);

  GetBestGreenToRed(tile_argb, xsize, tile_width, tile_height, prev_x, prev_y,
                    speed, accumulated_red_histo, &best_tx, tmp_histo);
  GetBestGreenRedToBlue(tile_argb, xsize, tile_width, tile_height, prev_x,
                        prev_y, speed, accumulated_blue_histo, &best_tx,
                        tmp_histo);
  return best_tx;
}

static void CopyTileWithColorTransform(uint32_t width, uint32_t height,
                                       int tile_x, int tile_y,
                                       uint32_t max_tile_size,
                                       WP2LMultipliers color_transform,
                                       int16_t* argb) {
  const int xscan = std::min(max_tile_size, width - tile_x);
  int yscan = std::min(max_tile_size, height - tile_y);
  argb += 4 * (tile_y * width + tile_x);
  while (yscan-- > 0) {
    WP2LTransformColor(&color_transform, xscan, argb);
    argb += 4 * width;
  }
}

WP2Status ColorSpaceTransform(uint32_t width, uint32_t height, uint32_t bits,
                              WP2SampleFormat format, uint32_t speed,
                              uint16_t* const argb_in, int16_t* image) {
  const int max_tile_size = 1 << bits;
  const int tile_xsize = WP2LSubSampleSize(width, bits);
  const int tile_ysize = WP2LSubSampleSize(height, bits);
  int tile_x, tile_y;
  WP2LMultipliers prev_x, prev_y;
  MultipliersClear(&prev_y);
  MultipliersClear(&prev_x);
  Vector_u32 tmp_histo, accumulated_red_histo, accumulated_blue_histo;
  // At most, we have a delta of twice a color.
  const uint32_t max_index_red = 3 * 2 * WP2::FormatMax(format, /*channel=*/1);
  const uint32_t max_index_blue = 3 * 2 * WP2::FormatMax(format, /*channel=*/3);
  WP2_CHECK_ALLOC_OK(
      tmp_histo.resize(std::max(max_index_red, max_index_blue) + 1));
  WP2_CHECK_ALLOC_OK(accumulated_red_histo.resize(max_index_red + 1));
  std::fill(accumulated_red_histo.begin(), accumulated_red_histo.end(), 0);
  WP2_CHECK_ALLOC_OK(accumulated_blue_histo.resize(max_index_blue + 1));
  std::fill(accumulated_blue_histo.begin(), accumulated_blue_histo.end(), 0);
  WP2::VectorNoCtor<int16_t> argb_signed;
  WP2_CHECK_ALLOC_OK(argb_signed.resize(4 * width * height));
  // Convert data to signed first.
  int16_t* const argb = argb_signed.data();
  std::copy(argb_in, argb_in + 4 * width * height, argb);
  for (tile_y = 0; tile_y < tile_ysize; ++tile_y) {
    for (tile_x = 0; tile_x < tile_xsize; ++tile_x) {
      int y;
      const uint32_t tile_x_offset = tile_x * max_tile_size;
      const uint32_t tile_y_offset = tile_y * max_tile_size;
      const int all_x_max = std::min(tile_x_offset + max_tile_size, width);
      const int all_y_max = std::min(tile_y_offset + max_tile_size, height);
      const int offset = tile_y * tile_xsize + tile_x;
      if (tile_y != 0) {
        ColorCodeToMultipliers(&image[4 * (offset - tile_xsize)], &prev_y);
      }
      prev_x = GetBestColorTransformForTile(
          tile_x, tile_y, bits, prev_x, prev_y, speed, width, height,
          accumulated_red_histo, accumulated_blue_histo, argb, &tmp_histo);
      MultipliersToColorCode(&prev_x, &image[4 * offset]);
      CopyTileWithColorTransform(width, height, tile_x_offset, tile_y_offset,
                                 max_tile_size, prev_x, argb);

      // Gather accumulated histogram data.
      for (y = tile_y_offset; y < all_y_max; ++y) {
        uint32_t ix = y * width + tile_x_offset;
        const uint32_t ix_end = ix + all_x_max - tile_x_offset;
        for (; ix < ix_end; ++ix) {
          if (ix >= 2 &&
              ColorEq(&argb[4 * ix], &argb[4 * (ix - 2)]) &&
              ColorEq(&argb[4 * ix], &argb[4 * (ix - 1)])) {
            continue;  // repeated pixels are handled by backward references
          }
          if (ix >= width + 2 &&
              ColorEq(&argb[4 * (ix - 2)], &argb[4 * (ix - width - 2)]) &&
              ColorEq(&argb[4 * (ix - 1)], &argb[4 * (ix - width - 1)]) &&
              ColorEq(&argb[4 * ix], &argb[4 * (ix - width)])) {
            continue;  // repeated pixels are handled by backward references
          }
          assert(MakeIndex(argb[4 * ix + 1]) < accumulated_red_histo.size());
          assert(MakeIndex(argb[4 * ix + 3]) < accumulated_blue_histo.size());
          ++accumulated_red_histo[MakeIndex(argb[4 * ix + 1])];
          ++accumulated_blue_histo[MakeIndex(argb[4 * ix + 3])];
        }
      }
    }
  }
  // Copy back data and wrap around.
  for (uint32_t i = 0; i < 4 * width * height; i += 4) {
    argb_in[i + 1] = argb[i + 1] & 0xff;
    argb_in[i + 3] = argb[i + 3] & 0xff;
  }

  return WP2_STATUS_OK;
}
}  // namespace WP2L
