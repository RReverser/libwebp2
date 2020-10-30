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
// main entry for the lossless encoder.
//
// Author: Vincent Rabaud (vrabaud@google.com)

#include "src/common/color_precision.h"
#include "src/dsp/lossless/lossless_common.h"
#include "src/enc/lossless/losslessi_enc.h"

namespace WP2L {

// -----------------------------------------------------------------------------
// Entropy analysis

// These five modes are evaluated and their respective entropy is computed.
enum class EntropyMode {
  // No transform, colors are stored as is.
  kDirect = 0,
  // Spatial prediction. We predict the next pixel (using one of many prediction
  // modes) and store the difference between the prediction and the actual value
  // (= "the residual").
  kSpatial = 1,
  // Instead of storing red, green, and blue, we store red minus green, green,
  // and blue minus green (takes advantage of correlation between the three
  // channels).
  kSubGreen = 2,
  // We use both kSpatial and kSubGreen a the same time.
  kSpatialSubGreen = 3,
  // We create a palette of colors and store one palette index per pixel.
  kPalette = 4,
  kNum = 5
};

// Histogram types.
typedef enum {
  kHistoAlpha = 0,
  kHistoAlphaPred,  // Alpha difference with previous pixel.
  kHistoGreen,
  kHistoGreenPred,  // Green difference with previous pixel.
  kHistoRed,
  kHistoRedPred,  // Red difference with previous pixel.
  kHistoBlue,
  kHistoBluePred,          // Blue difference with previous pixel.
  kHistoRedSubGreen,       // Red minus green.
  kHistoRedPredSubGreen,   // Red minus green difference with previous pixel.
  kHistoBlueSubGreen,      // Blue minus green.
  kHistoBluePredSubGreen,  // Blue minus green difference with previous pixel.
  kHistoPalette,
  kHistoTotal  // Must be last.
} HistoIx;

// Adds the red minus green, and blue minus green values from the given argb
// pixel 'p' to the given histograms 'r' and 'b'.
template <typename T>
static void AddSingleSubGreen(const T* const p, uint32_t mask,
                              uint32_t* const r, uint32_t* const b) {
  const int green = p[2];  // The upper bits are masked away later.
  ++r[((int)p[1] - green) & mask];
  ++b[((int)p[3] - green) & mask];
}

// Adds the values from the given argb pixel 'p' to the given histograms 'a',
// 'r', 'g', and 'b'.
template <typename T>
static void AddSingle(const T* const p, uint32_t* const a, uint32_t* const r,
                      uint32_t* const g, uint32_t* const b) {
  ++a[p[0]];
  ++r[p[1]];
  ++g[p[2]];
  ++b[p[3]];
}

// Analyzes the distribution of colors on the image to estimate the best
// transform, returned in 'min_entropy_ix'. Also sets 'red_and_blue_always_zero'
// to true if both red and blue channels are always zero when applying this
// transform.
template <typename T>
static WP2Status AnalyzeEntropyImpl(const T* curr_row, uint32_t stride,
                                    uint32_t width, uint32_t height,
                                    WP2SampleFormat format, bool use_palette,
                                    Encoder* const encoder,
                                    EntropyMode* const min_entropy_ix,
                                    bool* const red_and_blue_always_zero) {
  if (use_palette && encoder->palette_.Size() <= 16) {
    // In practice, small palettes are better than any other transform.
    *min_entropy_ix = EntropyMode::kPalette;
    *red_and_blue_always_zero = true;
    return WP2_STATUS_OK;
  }
  uint32_t max_channel = 0;
  uint32_t masks[4];
  uint32_t max_rgb = 0;
  for (uint32_t c = 0; c < 4; ++c) {
    masks[c] = WP2::FormatMax(format, c);
    assert(((masks[c] + 1) & masks[c]) == 0);
    if (c != 0) max_rgb = masks[c];
    max_channel = std::max(max_channel, masks[c]);
  }
  for (uint32_t c = 0; c < 4; ++c) {
    if (c != 0) {
      // We only support RGB of same depth in this code for now as it gets
      // tricky to subtract green from red if the depth is deifferent.
      assert(max_rgb == masks[c]);
    }
  }

  // Allocate histogram set with cache_bits = 0.
  WP2::Vector_u32 histo;
  const uint32_t histo_len = max_channel + 1;
  WP2_CHECK_ALLOC_OK(histo.resize(kHistoTotal * histo_len));
  std::fill(histo.begin(), histo.end(), 0);

  const T* prev_row = nullptr;
  // Skip the first pixel.
  T pix_prev[4] = {curr_row[0], curr_row[1], curr_row[2], curr_row[3]};
  for (uint32_t y = 0; y < height; ++y) {
    for (uint32_t x = 0; x < width; ++x) {
      // Compute the difference between the two pixels.
      const T* pix = &curr_row[4 * x];
      uint32_t pix_diff[4];
      for (uint32_t c = 0; c < 4; ++c) {
        pix_diff[c] = ((int32_t)pix[c] - (int32_t)pix_prev[c]) & masks[c];
      }
      std::copy(pix, pix + 4, pix_prev);
      if ((pix_diff[0] + pix_diff[1] + pix_diff[2] + pix_diff[3] == 0) ||
          (prev_row != nullptr
              && memcmp(pix, &prev_row[4 * x], 4 * sizeof(uint8_t)) == 0)) {
        continue;
      }
      AddSingle(pix,
                &histo[kHistoAlpha * histo_len],
                &histo[kHistoRed * histo_len],
                &histo[kHistoGreen * histo_len],
                &histo[kHistoBlue * histo_len]);
      AddSingle(pix_diff,
                &histo[kHistoAlphaPred * histo_len],
                &histo[kHistoRedPred * histo_len],
                &histo[kHistoGreenPred * histo_len],
                &histo[kHistoBluePred * histo_len]);
      AddSingleSubGreen(pix, max_rgb, &histo[kHistoRedSubGreen * histo_len],
                        &histo[kHistoBlueSubGreen * histo_len]);
      AddSingleSubGreen(pix_diff, max_rgb,
                        &histo[kHistoRedPredSubGreen * histo_len],
                        &histo[kHistoBluePredSubGreen * histo_len]);
      // Approximate the palette by the entropy of the multiplicative hash.
      const uint32_t hash = HashPix(pix, WP2Log2Floor(max_channel + 1));
      ++histo[kHistoPalette * histo_len + hash];
    }
    prev_row = curr_row;
    curr_row = (const T*)((const uint8_t*)curr_row + stride);
  }

  double entropy_comp[kHistoTotal];
  double entropy[(uint32_t)EntropyMode::kNum];
  const EntropyMode last_mode_to_analyze =
      use_palette ? EntropyMode::kPalette : EntropyMode::kSpatialSubGreen;
  // Let's add one zero to the predicted histograms. The zeros are removed
  // too efficiently by the pix_diff == 0 comparison, at least one of the
  // zeros is likely to exist.
  ++histo[kHistoRedPredSubGreen * histo_len];
  ++histo[kHistoBluePredSubGreen * histo_len];
  ++histo[kHistoRedPred * histo_len];
  ++histo[kHistoGreenPred * histo_len];
  ++histo[kHistoBluePred * histo_len];
  ++histo[kHistoAlphaPred * histo_len];

  for (uint32_t j = 0; j < kHistoTotal; ++j) {
    entropy_comp[j] = WP2::ANSCountsQuantizedCost(
        &histo[j * histo_len], &histo[j * histo_len], histo_len);
  }
  entropy[(uint32_t)EntropyMode::kDirect] =
      entropy_comp[kHistoAlpha] + entropy_comp[kHistoRed] +
      entropy_comp[kHistoGreen] + entropy_comp[kHistoBlue];
  entropy[(uint32_t)EntropyMode::kSpatial] =
      entropy_comp[kHistoAlphaPred] + entropy_comp[kHistoRedPred] +
      entropy_comp[kHistoGreenPred] + entropy_comp[kHistoBluePred];
  entropy[(uint32_t)EntropyMode::kSubGreen] =
      entropy_comp[kHistoAlpha] + entropy_comp[kHistoRedSubGreen] +
      entropy_comp[kHistoGreen] + entropy_comp[kHistoBlueSubGreen];
  entropy[(uint32_t)EntropyMode::kSpatialSubGreen] =
      entropy_comp[kHistoAlphaPred] + entropy_comp[kHistoRedPredSubGreen] +
      entropy_comp[kHistoGreenPred] + entropy_comp[kHistoBluePredSubGreen];
  entropy[(uint32_t)EntropyMode::kPalette] = entropy_comp[kHistoPalette];

  // When including transforms, there is an overhead in bits from
  // storing them. This overhead is small but matters for small images.
  // For spatial, there are 14 transformations.
  const int transform_bits = encoder->transform_bits_;
  entropy[(uint32_t)EntropyMode::kSpatial] +=
      WP2LSubSampleSize(width, transform_bits) *
      WP2LSubSampleSize(height, transform_bits) * WP2Log2(14);
  // For color transforms: 24 as only 3 channels are considered in a
  // ColorTransformElement.
  entropy[(uint32_t)EntropyMode::kSpatialSubGreen] +=
      WP2LSubSampleSize(width, transform_bits) *
      WP2LSubSampleSize(height, transform_bits) * WP2Log2(24);
  // For palettes, add the cost of storing the palette.
  if (use_palette) {
    const float palette_cost = encoder->palette_.GetCost();
    entropy[(uint32_t)EntropyMode::kPalette] += palette_cost;
  }

  *min_entropy_ix = EntropyMode::kDirect;
  for (uint32_t k = (uint32_t)EntropyMode::kDirect + 1;
       k <= (uint32_t)last_mode_to_analyze; ++k) {
    if (entropy[(uint32_t)*min_entropy_ix] > entropy[k]) {
      *min_entropy_ix = (EntropyMode)k;
    }
  }
  assert((uint32_t)*min_entropy_ix <= (uint32_t)last_mode_to_analyze);
  *red_and_blue_always_zero = true;
  // Let's check if the histogram of the chosen entropy mode has
  // non-zero red and blue values. If all are zero, we can later skip
  // the cross color optimization.

  static const uint8_t kHistoPairs[5][2] = {
      {kHistoRed, kHistoBlue},
      {kHistoRedPred, kHistoBluePred},
      {kHistoRedSubGreen, kHistoBlueSubGreen},
      {kHistoRedPredSubGreen, kHistoBluePredSubGreen},
      {kHistoRed, kHistoBlue}};
  const uint32_t* const red_histo =
      &histo[histo_len * kHistoPairs[(uint32_t)*min_entropy_ix][0]];
  const uint32_t* const blue_histo =
      &histo[histo_len * kHistoPairs[(uint32_t)*min_entropy_ix][1]];
  for (size_t i = 1; i < histo_len; ++i) {
    if ((red_histo[i] | blue_histo[i]) != 0) {
      *red_and_blue_always_zero = false;
      break;
    }
  }

  return WP2_STATUS_OK;
}

static WP2Status AnalyzeEntropy(const WP2::ArgbBuffer& pic, uint32_t width,
                                uint32_t height, bool use_palette,
                                Encoder* const encoder,
                                EntropyMode* const min_entropy_ix,
                                bool* const red_and_blue_always_zero) {
  if (WP2FormatBpc(pic.format) == 1) {
    WP2_CHECK_STATUS(AnalyzeEntropyImpl(
        pic.GetRow8(0), pic.stride, width, height, pic.format, use_palette,
        encoder, min_entropy_ix, red_and_blue_always_zero));
  } else {
    WP2_CHECK_STATUS(AnalyzeEntropyImpl(
        pic.GetRow16(0), pic.stride, width, height, pic.format, use_palette,
        encoder, min_entropy_ix, red_and_blue_always_zero));
  }

  return WP2_STATUS_OK;
}

static int GetHistoBits(int speed, bool use_palette, int width, int height) {
  // Make tile size a function of encoding method (Range: 0 to 6).
  // TODO(vrabaud) use rounding. We leave it as is so that for the default
  // value of 5, it matches the old default value of 3.
  uint32_t histo_bits = (use_palette ? 9. : 7.) - speed * 6. / 9.;
  while (true) {
    const uint32_t huff_image_size = WP2LSubSampleSize(width, histo_bits) *
                                     WP2LSubSampleSize(height, histo_bits);
    if (huff_image_size <= kMaxHistogramImageSize) break;
    ++histo_bits;
  }
  return (histo_bits < kHistogramBitsMin) ? kHistogramBitsMin :
         (histo_bits > kHistogramBitsMax) ? kHistogramBitsMax : histo_bits;
}

static uint32_t GetTransformBits(uint32_t speed, uint32_t histo_bits) {
  const uint32_t max_transform_bits = (speed < 5) ? 6 : (speed > 5) ? 4 : 5;
  const uint32_t res =
      (histo_bits > max_transform_bits) ? max_transform_bits : histo_bits;
  assert(res <= kTransformBitsMax);
  return res;
}

// Translates an entropy mode to a CrunchConfig.
static void TranslateEntropyMode(EntropyMode entropy_idx, uint32_t speed,
                                 bool red_and_blue_always_zero,
                                 CrunchConfig* const config) {
  config->use_group4 = false;
  config->use_palette = (entropy_idx == EntropyMode::kPalette);
  config->use_subtract_green = (entropy_idx == EntropyMode::kSubGreen) ||
                               (entropy_idx == EntropyMode::kSpatialSubGreen);
  config->use_predict = (entropy_idx == EntropyMode::kSpatial) ||
                        (entropy_idx == EntropyMode::kSpatialSubGreen);
  if (speed == 0) {
    config->use_cross_color = false;
  } else {
    config->use_cross_color =
        red_and_blue_always_zero ? false : config->use_predict;
  }
}

// Fills in 'crunch_configs' with a set of transforms and LZ77 variants to try
// based on heuristics. Also initializes the palette in enc->palette_ and
// sets red_and_blue_always_zero to true if the red and blue components are
// always zero after applying the best transform.
// TODO(maryla): the logic for red_and_blue_always_zero seems broken: it's only
// valid if there's only one transform (speed != 9)
WP2Status EncoderAnalyze(Encoder* const encoder,
                         CrunchConfig crunch_configs[kCrunchNumMax],
                         uint32_t* const crunch_configs_size) {
  const WP2::ArgbBuffer& pic = encoder->pic_;
  const uint32_t width = pic.width;
  const uint32_t height = pic.height;
  const int speed = encoder->config_.speed;
  uint32_t n_lz77s;
  assert(!pic.IsEmpty());

  WP2_CHECK_STATUS(encoder->palette_.AnalyzeAndCreate(pic));
  WP2_CHECK_STATUS(encoder->palette_.FindBestMethod(encoder, speed));

  // Do not use a palette when there is only one color as ANS will optimize it
  // anyway.
  const bool use_palette = (encoder->palette_.Size() > 1);

  // TODO(jyrki): replace the decision to be based on an actual estimate
  // of entropy, or even spatial variance of entropy.
  encoder->histo_bits_ =
      GetHistoBits(speed, use_palette, pic.width, pic.height);
  encoder->transform_bits_ = GetTransformBits(speed, encoder->histo_bits_);

  *crunch_configs_size = 0;
  if (false) {
    // Try group4 for images with a palette.
    CrunchConfig* config;
    config = &crunch_configs[*crunch_configs_size];
    config->use_group4 = true;
    config->group4_use_move_to_front = false;
    config->use_palette = false;
    config->use_subtract_green = false;
    config->use_predict = false;
    config->use_cross_color = false;
    ++(*crunch_configs_size);

    if (encoder->palette_.Size() > 3 && speed > 0) {
      config = &crunch_configs[*crunch_configs_size];
      // Copy the previous config.
      *config = crunch_configs[*crunch_configs_size - 1];
      // Try using a MoveToFrontCache for images with more than 3 colors.
      config->group4_use_move_to_front = true;
      ++(*crunch_configs_size);
    }
  }
  if (speed == 0) {
    if (*crunch_configs_size != 1) {
      // AnalyzeEntropy is somewhat slow. 'red_and_blue_always_zero' is unused
      // for speed == 0.
      TranslateEntropyMode(
          use_palette ? EntropyMode::kPalette : EntropyMode::kSpatialSubGreen,
          speed, /*red_and_blue_always_zero=*/false, &crunch_configs[0]);
      *crunch_configs_size = 1;
    }
    n_lz77s = 1;
  } else {
    EntropyMode min_entropy_ix = EntropyMode::kDirect;
    // Try out multiple LZ77 on images with few colors.
    n_lz77s = (use_palette && encoder->palette_.Size() <= 16) ? 2 : 1;
    // TODO(maryla): if speed == 9 this doesn't seem useful.
    bool red_and_blue_always_zero = false;
    WP2_CHECK_STATUS(AnalyzeEntropy(pic, width, height, use_palette, encoder,
                                    &min_entropy_ix,
                                    &red_and_blue_always_zero));
    if (speed == 9) {
      // Go brute force on all transforms.
      static_assert(1 + (uint32_t)EntropyMode::kNum <= kCrunchNumMax,
                    "Too many entropy types");
      for (uint32_t i = 0; i < (uint32_t)EntropyMode::kNum; ++i) {
        if (i != (uint32_t)EntropyMode::kPalette || use_palette) {
          assert(*crunch_configs_size < kCrunchNumMax);
          TranslateEntropyMode((EntropyMode)i, speed, red_and_blue_always_zero,
                               &crunch_configs[(*crunch_configs_size)++]);
        }
      }
    } else {
      // Only choose the guessed best transform.
      TranslateEntropyMode(min_entropy_ix, speed, red_and_blue_always_zero,
                           &crunch_configs[(*crunch_configs_size)++]);
    }
  }
  // Fill in the different LZ77s.
  assert(n_lz77s <= kCrunchLZ77Max);
  for (uint32_t i = 0; i < *crunch_configs_size; ++i) {
    for (uint32_t j = 0; j < n_lz77s; ++j) {
      crunch_configs[i].lz77s_types_to_try[j] =
          (j == 0) ? kLZ77Standard | kLZ77RLE : kLZ77Box;
    }
    crunch_configs[i].lz77s_types_to_try_size = n_lz77s;
  }
  return WP2_STATUS_OK;
}

//------------------------------------------------------------------------------

}  // namespace WP2L
