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
// Internal constants.
//
// Authors: Vincent Rabaud (vrabaud@google.com)

#ifndef WP2_COMMON_CONSTANTS_H_
#define WP2_COMMON_CONSTANTS_H_

#include <cmath>

#include "src/dsp/math.h"
#include "src/utils/utils.h"
#include "src/wp2/base.h"
#include "src/wp2/format_constants.h"

//------------------------------------------------------------------------------
// Lossy and common constants.

namespace WP2 {

static constexpr uint32_t kLosslessQualityHint = (1 << kQualityHintNumBits) - 1;
static constexpr uint32_t kMaxLossyQualityHint = kLosslessQualityHint - 1;
static constexpr uint32_t GetQualityHint(float quality) {
  return (quality <= kMaxLossyQuality)
             ? DivRound<uint32_t>(quality * kMaxLossyQualityHint,
                                  kMaxLossyQuality)
             : kLosslessQualityHint;
}

static constexpr uint32_t kMaxVarIntLength = 5;
// Upper value expressible with a var int with kMaxVarIntLength values.
static constexpr uint64_t kUpperVarInt =
    (1ull << (7 * (kMaxVarIntLength - 1) + 8));  // == 64G
static_assert(kChunkSizeMax < kUpperVarInt,
              "The maximum chunk size cannot be represented with var ints.");

// Maximum (inclusive) bytes per pixel.
static constexpr uint64_t kMaxNumBytesPerPixel =
#if defined(WP2_ENC_DEC_MATCH)
    // Each channel is sizeof(uint16_t) and its ANS hash sizeof(uint16_t).
    4 * 2 * sizeof(uint16_t);
#else
    4 * sizeof(uint16_t);
#endif

static constexpr char kChannelStr[4][2] = {"Y", "U", "V", "A"};
const char* const kPartitionMethodString[] = {
    "Multipass", "Block-enc",    "Split-enc",   "Tile-enc",
    "Area-enc",  "Sub-area-enc", "Exhaustive",  "Fixed-4x4",
    "Fixed-8x8", "Fixed-16x16",  "Fixed-32x32", "Auto"};
STATIC_ASSERT_ARRAY_SIZE(kPartitionMethodString, WP2::NUM_PARTITION_METHODS);

// Number of sets of residual statistics. For now, there is one set per method.
// kMethod0 and kMethod1.
static constexpr uint32_t kNumResidualStats = 2;
// Maximum value used by StoreCoeffs to use the first dictionary.
// A value of kResidual1Max is an escape code to get to the second
// dictionary.
static constexpr uint32_t kResidual1Max = 12;
// Maximum value of consecutive 0s we store for residuals.
static constexpr uint32_t kResidualCons0Max = 6;
// Maximum quantized DC value that can be coded.
static constexpr int32_t kMaxDcValue = (1 << 12) - 1;
// Maximum quantized value that can be coded for AC coeffs.
static constexpr int32_t kMaxCoeffValue =
    ((1 << kMaxCoeffBits) - 1) + 3 + kResidual1Max;

static constexpr Argb38b kTransparentArgb38b{0x00u, 0x0000u, 0x0000u, 0x0000u};

// Decoding progress proportion in [0:1] range. The remaining is left to pixels.
static constexpr float kProgressDecHeader =   0.01f * 1;  // % for readability
static constexpr float kProgressDecPreview =  0.01f * 3;
static constexpr float kProgressDecICC =      0.01f * 2;
static constexpr float kProgressDecANMF =     0.01f * 1;
static constexpr float kProgressDecMetadata = 0.01f * 4;
static constexpr float kProgressDecEnd =      0.01f * 0.1f;

// Returns the tile width/height in pixels for a given TileShape and a given
// image width.
static inline uint32_t TileWidth(TileShape tile_shape, uint32_t image_width) {
  switch (tile_shape) {
    case TILE_SHAPE_SQUARE_128:
      return 128;
    case TILE_SHAPE_SQUARE_256:
      return 256;
    case TILE_SHAPE_SQUARE_512:
      return 512;
    case TILE_SHAPE_WIDE:
      return std::min(Pad(image_width, kMinBlockSizePix), kMaxTileSize);
    default:
      assert(false);
      return 0;
  }
}
static inline uint32_t TileHeight(TileShape tile_shape, uint32_t image_width) {
  const uint32_t tile_width = TileWidth(tile_shape, image_width);
  switch (tile_shape) {
    case TILE_SHAPE_SQUARE_128:
    case TILE_SHAPE_SQUARE_256:
    case TILE_SHAPE_SQUARE_512:
      return tile_width;
    case TILE_SHAPE_WIDE:
      // Should this be smaller than width?
      return std::min(Pad(kMaxTilePixels / tile_width, kMinBlockSizePix),
                      kMaxTileSize);
    default:
      assert(false);
      return 0;
  }
}

}  // namespace WP2

//------------------------------------------------------------------------------
// Lossless constants.

namespace WP2L {

typedef enum {
  PREDICTOR_TRANSFORM      = 0,
  CROSS_COLOR_TRANSFORM    = 1,
  SUBTRACT_GREEN           = 2,
  COLOR_INDEXING_TRANSFORM = 3,
  GROUP4                   = 4,
  NUM_TRANSFORMS           = 5
} ImageTransformType;

// Modes used in Group4 compression: pass is just about skipping the definition
// of a segment, horizontal mode is about defining the next two segments,
// vertical is about defining a segment with respect to the line above.
enum class Group4Mode { kVertical, kHorizontal, kPass };

// Gives the maximum possible palette size for a given image size.
static inline uint32_t MaxPaletteSize(uint32_t num_pixels) {
  // We only allow palettes of at least 2 elements.
  if (num_pixels < 2) return 2;
  // The palette cannot be bigger than the number of pixels.
  return std::min(kMaxPaletteSize, num_pixels);
}

}  // namespace WP2L

//------------------------------------------------------------------------------

#endif /* WP2_COMMON_CONSTANTS_H_ */
