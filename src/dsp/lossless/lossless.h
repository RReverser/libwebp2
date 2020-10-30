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
// Authors: Vincent Rabaud (vrabaud@google.com)

#ifndef WP2_DSP_LOSSLESS_LOSSLESS_H_
#define WP2_DSP_LOSSLESS_LOSSLESS_H_

#include <algorithm>
#include <cstdint>

#ifdef WP2_EXPERIMENTAL_FEATURES
#include "src/enc/delta_palettization_enc.h"
#endif  // WP2_EXPERIMENTAL_FEATURES


//------------------------------------------------------------------------------
// Common

// Transforms a signed value into an unsigned index.
static inline uint16_t MakeIndex(int16_t val) {
  // TODO(vrabaud) remove this wrap around.
  return (val & 0xff);
}

//------------------------------------------------------------------------------
// Decoding

constexpr uint32_t kNumPredictors = 16;

typedef void (*WP2LPredictorFunc)(const uint16_t* const left,
                                  const uint16_t* const top, uint16_t max_value,
                                  uint16_t* const out);
extern WP2LPredictorFunc WP2LPredictors[kNumPredictors];
extern WP2LPredictorFunc WP2LPredictors_C[kNumPredictors];
// These Add/Sub function expects upper[-1] and out[-1] to be readable.
typedef void (*WP2LPredictorAddSubFunc)(const uint16_t* in, bool has_alpha,
                                        const uint16_t* upper,
                                        uint32_t num_pixels,
                                        uint32_t channel_bits, uint16_t* out);
extern WP2LPredictorAddSubFunc WP2LPredictorsAdd[kNumPredictors];
extern WP2LPredictorAddSubFunc WP2LPredictorsAdd_C[kNumPredictors];

typedef void (*WP2LProcessDecBlueAndRedFunc)(const uint16_t* src,
                                             uint32_t num_pixels,
                                             uint32_t n_bits, uint16_t* dst);
extern WP2LProcessDecBlueAndRedFunc WP2LAddGreenToBlueAndRed;

typedef struct {
  int16_t green_to_red;
  int16_t green_to_blue;
  int16_t red_to_blue;
} WP2LMultipliers;
typedef void (*WP2LTransformColorInverseFunc)(const WP2LMultipliers* const m,
                                              const uint16_t* src,
                                              uint32_t num_pixels,
                                              uint32_t channel_bits,
                                              uint16_t* dst);
extern WP2LTransformColorInverseFunc WP2LTransformColorInverse;

namespace WP2L {
class Transform;  // Defined in dec/lossless/losslessi_dec.h.
}

// Performs inverse transform of data given transform information, start and end
// rows. Transform will be applied to rows [row_start, row_end[.
// The *in and *out pointers refer to source and destination data respectively
// corresponding to the intermediate row (row_start).
void WP2LInverseTransform(const WP2L::Transform* const transform,
                          uint32_t row_start, uint32_t row_end,
                          uint32_t channel_bits, const uint16_t* const in,
                          bool has_alpha, uint16_t* const out);

// Fills 'dst' with the samples from the 'color_map' at the indices specified by
// the green channel of the pixels in 'src'.
typedef void (*WP2LMapARGBFunc)(const uint16_t* src,
                                const int16_t* const color_map,
                                uint32_t color_map_size, uint32_t y_start,
                                uint32_t y_end, uint32_t width, uint16_t* dst);

extern WP2LMapARGBFunc WP2LMapColor;

// Expose some C-only fallback functions
void WP2LTransformColorInverse_C(const WP2LMultipliers* const m,
                                 const uint32_t* src, int num_pixels,
                                 uint32_t* dst);


void WP2LAddGreenToBlueAndRed_C(const uint16_t* src, uint32_t num_pixels,
                                uint32_t channel_bits, uint16_t* dst);

// Must be called before calling any of the above methods.
void WP2LDspInit();

//------------------------------------------------------------------------------
// Encoding

typedef void (*WP2LProcessEncBlueAndRedFunc)(uint16_t* dst,
                                             uint32_t channel_bits,
                                             uint32_t num_pixels);
extern WP2LProcessEncBlueAndRedFunc WP2LSubtractGreenFromBlueAndRed;
typedef void (*WP2LTransformColorFunc)(const WP2LMultipliers* const m,
                                       int num_pixels, int16_t* const data);
extern WP2LTransformColorFunc WP2LTransformColor;
typedef void (*WP2LCollectColorBlueTransformsFunc)(
    const int16_t* argb, uint32_t width, uint32_t tile_width,
    uint32_t tile_height, int16_t green_to_blue, int16_t red_to_blue,
    uint32_t* const histo);
extern WP2LCollectColorBlueTransformsFunc WP2LCollectColorBlueTransforms;

typedef void (*WP2LCollectColorRedTransformsFunc)(
    const int16_t* argb, uint32_t width, uint32_t tile_width,
    uint32_t tile_height, int16_t green_to_red, uint32_t* const histo);
extern WP2LCollectColorRedTransformsFunc WP2LCollectColorRedTransforms;

// Expose some C-only fallback functions
void WP2LTransformColor_C(const WP2LMultipliers* const m, int num_pixels,
                          int16_t* const data);
void WP2LSubtractGreenFromBlueAndRed_C(uint16_t* argb_data,
                                       uint32_t channel_bits,
                                       uint32_t num_pixels);
void WP2LCollectColorRedTransforms_C(const int16_t* argb, uint32_t width,
                                     uint32_t tile_width, uint32_t tile_height,
                                     int16_t green_to_red,
                                     uint32_t* const histo);
void WP2LCollectColorBlueTransforms_C(const int16_t* argb, uint32_t width,
                                      uint32_t tile_width, uint32_t tile_height,
                                      int16_t green_to_blue,
                                      int16_t red_to_blue,
                                      uint32_t* const histo);

extern WP2LPredictorAddSubFunc WP2LPredictorsSub[kNumPredictors];
extern WP2LPredictorAddSubFunc WP2LPredictorsSub_C[kNumPredictors];

// -----------------------------------------------------------------------------
// Entropy-cost related functions.

typedef void (*WP2LBufferAddFunc)(const uint32_t* const a,
                                  const uint32_t* const b, uint32_t size,
                                  uint32_t* const out,
                                  uint32_t* const nonzeros);
extern WP2LBufferAddFunc WP2LBufferAdd;

typedef double (*WP2LCostFunc)(const uint32_t* population, int length);
typedef double (*WP2LCostCombinedFunc)(const uint32_t* X, const uint32_t* Y,
                                       int length);
typedef float (*WP2LCombinedShannonEntropyFunc)(const uint32_t* X,
                                                const uint32_t* Y,
                                                uint32_t length);

extern WP2LCostFunc WP2LExtraCost;
extern WP2LCostCombinedFunc WP2LExtraCostCombined;
extern WP2LCombinedShannonEntropyFunc WP2LCombinedShannonEntropy;

// -----------------------------------------------------------------------------

// Must be called before calling any of the above methods.
void WP2LEncDspInit();

//------------------------------------------------------------------------------

#endif  // WP2_DSP_LOSSLESS_LOSSLESS_H_
