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
//  Internal header for constants related to WP2 file format.
//
// Author: Skal (pascal.massimino@gmail.com)

#ifndef WP2_WP2_FORMAT_CONSTANTS_H_
#define WP2_WP2_FORMAT_CONSTANTS_H_

#include <cstddef>
#include <cstdint>

namespace WP2 {

// maximum width/height allowed (inclusive), in pixels
// width/height limit
static constexpr uint32_t kMaxBufferDimension = (1u << 24);
// width * height limit
static constexpr uint32_t kMaxBufferArea = (1u << 28);

static constexpr uint32_t kSignature = 0x6ffff4;
static constexpr uint32_t kMaxPreviewSize = 333;
static constexpr uint32_t kHeaderMinSize = 10;
static constexpr uint32_t kHeaderMaxSize = 17;
static constexpr uint32_t kTagMask = 0x60ffe9;

// Maximum (inclusive) size of a buffer chunk. -1 to fit in our var ints.
static constexpr uint64_t kChunkSizeMax = (1ull << 36) - 1;  // ~=64G

static constexpr const char kFileExtension[] = "wp2";

// Number of bits used to store image dimension or frame offset.
static constexpr uint32_t kImageDimNumBits = 14;
// Maximum width/height allowed (inclusive), in pixels. We do not need to remove
// 1 to fit in kImageDimNumBits bits as we always remove one when storing
// dimensions (as they cannot be 0).
static constexpr uint32_t kImageDimMax = (1 << kImageDimNumBits);

// Qualities up to 95 included are considered lossy; above is lossless.
static constexpr uint32_t kMaxLossyQuality = 95;
static constexpr uint32_t kMaxQuality = 100;
static constexpr uint32_t kQualityHintNumBits = 4;
static constexpr uint32_t kAlphaQualityRange = kMaxLossyQuality + 2;

// For animations:
// Header chunk
static constexpr uint8_t kLoopCountNumBits = 6;
static constexpr uint32_t kMaxLoopCount =
    (1u << kLoopCountNumBits) - 1u;  // = 63, inclusive
static constexpr uint8_t kInfiniteLoop = 0;
static constexpr uint8_t kBackgroundNumBits = 2;
static constexpr uint8_t kBackgroundTransparent = 0x0;
static constexpr uint8_t kBackgroundWhite = 0x1;
static constexpr uint8_t kBackgroundBlack = 0x2;
static constexpr uint8_t kBackgroundCustom = 0x3;
// ANMF chunk
static constexpr size_t kANMFHeaderSize = 10;
static constexpr uint8_t kANMFTagNumBits = 6;
static constexpr uint32_t kANMFTagDispose = 0x33;
static constexpr uint32_t kANMFTagFrame = 0x15;
static constexpr uint8_t kFrameDurationNumBits = 16;
static constexpr uint32_t kMaxFrameDurationMs =
    (1u << kFrameDurationNumBits) - 1u;  // ~= 65 seconds
// Maximum number of frames in an animation, including preframes.
// Compute kMaxNumFrames so that the animation's total duration fits in 31 bits.
static constexpr uint32_t kMaxNumFrames =
    (1u << (31 - kFrameDurationNumBits)) - 1u;  // = 32767, inclusive
// Number of allowed consecutive preframes (meaning with a duration of 0 ms).
static constexpr uint32_t kMaxNumPreframes = 4;
// Number of allowed consecutive non-disposed regular frames (ignore preframes).
// To take preframes into account:
//   (1 + kMaxNumDependentFrames) * (1 + kMaxNumPreframes) - 1
static constexpr uint32_t kMaxNumDependentFrames = 16;

// Maximum inclusive angle step magnitude. MAX_ANGLE_DELTA in AV1.
// An angle predictor at angle 'angle' will generate sub-predictors at angles
// 'angle' +/- kDirectionalMaxAngleDelta * 'angle_unit' where 'angle_unit'
// subdivides 22.5 degrees in (2 * kDirectionalMaxAngleDeltaYA + 1) angles.
constexpr uint32_t kDirectionalMaxAngleDeltaYA = 3u;
constexpr uint32_t kDirectionalMaxAngleDeltaUV = 1u;
constexpr uint32_t kYBasePredNum = 8;
constexpr uint32_t kABasePredNum = 8;
constexpr uint32_t kUVBasePredNum = 4;
constexpr uint32_t kAnglePredNum = 8;

enum BasePredictor {
  BPRED_DC = 0,
  BPRED_DC_L,  // dc-left
  BPRED_DC_T,  // dc-top
  BPRED_SMOOTH,
  BPRED_SMOOTH_H,
  BPRED_SMOOTH_V,
  BPRED_TM,
  BPRED_GRADIENT,
  BPRED_LAST = kYBasePredNum
};

// Algorithm used for image partitioning.
enum PartitionMethod {
  // Combine several heuristics in successive block size passes:
  MULTIPASS_PARTITIONING,        // (fast)
  // For each pos starting from top left, take the best block size depending on:
  BLOCK_ENCODE_PARTITIONING,     // block encoding score (slow)
  SPLIT_RECURSE_PARTITIONING,    // whole or split block encoding score (slow)
  TILE_ENCODE_PARTITIONING,      // tile encoding score (super slow)
  // Per 32x32 area, take the best one of a few block layouts:
  AREA_ENCODE_PARTITIONING,      // (slow)
  SUB_AREA_ENCODE_PARTITIONING,  // same but for each block in area (very slow)
  // Per tile, take the best one of all possible block layouts:
  EXHAUSTIVE_PARTITIONING,       // (extremely slow)
  // Fixed block size (except on edges). Also depends on the partition set.
  ALL_4X4_PARTITIONING,
  ALL_8X8_PARTITIONING,
  ALL_16X16_PARTITIONING,
  ALL_32X32_PARTITIONING,
  // Choose one above based on compression effort and image size.
  AUTO_PARTITIONING,
  NUM_PARTITION_METHODS
};

// The set of allowed block sizes for image partitioning.
enum PartitionSet {  // The smallest block size is 4x4.
  SMALL_SQUARES,     // Up to 8x8
  SMALL_RECTS,       // Up to 16x16
  ALL_RECTS,         // Up to 32x32, ratio at most 4:1
  THICK_RECTS,       // Up to 32x32, ratio at most 2:1
  MEDIUM_SQUARES,    // Up to 16x16
  ALL_SQUARES,       // Up to 32x32
  SOME_RECTS,        // Up to 32x32, subset of frequently used rects
  NUM_PARTITION_SETS
};

static constexpr uint32_t kMaxNumSegments = 8;

typedef enum {
  WP2_TF_ITU_R_BT2020_10BIT = 0,
  WP2_TF_ITU_R_BT709,
  WP2_TF_UNSPECIFIED,
  WP2_TF_GAMMA_22,
  WP2_TF_GAMMA_28,
  WP2_TF_SMPTE_170M,
  WP2_TF_SMPTE_240M,
  WP2_TF_LINEAR,
  WP2_TF_LOG,
  WP2_TF_SQRT,
  WP2_TF_IEC_61966_2_4,
  WP2_TF_ITU_R_BT1361_EXTENDED,
  WP2_TF_IEC_61966_2_1,
  WP2_TF_ITU_R_BT2020_12BIT,
  WP2_TF_SMPTE_ST_2084,
  WP2_TF_SMPTE_ST_428_1,
  WP2_TF_ARIB_STD_B67_HLG
} TransferFunction;

typedef enum { kYChannel = 0, kUChannel, kVChannel, kAChannel } Channel;

static constexpr uint32_t kRgbBits = 8u;
static constexpr uint32_t kRgbMax = (1 << kRgbBits) - 1;  // Inclusive.
static constexpr uint32_t kAlphaBits = 8u;
static constexpr uint32_t kAlphaMax = (1 << kAlphaBits) - 1;  // Inclusive.

// Max precision for YUV values, excluding the sign.
static constexpr uint32_t kMaxYuvBits = 9u;

// Prediction block dimensions:
static constexpr uint32_t kMinBlockSizePix = 4u;   // min block size in pixels
static constexpr uint32_t kMaxBlockSizePix = 32u;  // max block size in pixels
// maximum block size  (for stack reservation)
static constexpr uint32_t kMaxBlockSizePix2 =
    (kMaxBlockSizePix * kMaxBlockSizePix);
// blocks dimensions (w or h) in kMinBlockSizePix units are [1..kMaxBlockSize]
// That's 64 possible block types.
static constexpr uint32_t kMaxBlockSize = (kMaxBlockSizePix / kMinBlockSizePix);
static constexpr uint32_t kMaxBlockSize2 = (kMaxBlockSize * kMaxBlockSize);

// Tile sizes.
static constexpr uint32_t kTileShapeBits = 2;
static constexpr uint32_t kMaxTilePixels = 512 * 512;
static constexpr uint32_t kMaxTileSize =
    kMaxTilePixels / (2 * kMaxBlockSizePix);
// Extreme aspect ratios are undesirable as they create more borders with less
// context.
static constexpr uint32_t kMaxAspectRatio = 64;
static_assert(kMaxTileSize * kMaxTileSize / kMaxTilePixels <= kMaxAspectRatio,
              "Extreme aspect ratios are no good");

typedef enum {
  TILE_SHAPE_SQUARE_128,
  TILE_SHAPE_SQUARE_256,
  TILE_SHAPE_SQUARE_512,
  // As wide as possible (with a max of kMaxTileSize) and height set so that
  // it has at most kMaxTilePixels.
  TILE_SHAPE_WIDE,
  TILE_SHAPE_AUTO,
  NUM_TILE_SHAPES,
} TileShape;

// max coded level. number of prefix bits for the dictionaries.
static constexpr uint32_t kMaxCoeffBits = 10;

// convert pixels units to min-block-size units (rounding up)
static constexpr uint32_t SizeBlocks(uint32_t size_pix) {
  return ((size_pix + kMinBlockSizePix - 1u) / kMinBlockSizePix);
}

// Prediction block size
static constexpr uint32_t kPredWidth = 4;
static constexpr uint32_t kPredHeight = 4;
static constexpr uint32_t kPredSize = (kPredWidth * kPredHeight);
static constexpr uint32_t kLinearPredBits = 10;  // For calculation.
// Size of the usual context for prediction (it goes all around a block).
static constexpr uint32_t ContextSize(uint32_t bw, uint32_t bh) {
  return bw + 2 * (bh + 1);
}
// Size of the context going left of a block, top, and top right extending the
// top context. Extending by kMaxBlockSizePix is good enough for 45 degree
// predictors but we need more for lower angles, hence 2 * kMaxBlockSizePix.
static constexpr uint32_t kMaxContextTRExtent = 2 * kMaxBlockSizePix;
static constexpr uint32_t ContextWithTRSize(uint32_t bw, uint32_t bh) {
  return bh + 1 + bw + kMaxContextTRExtent;
}
static constexpr uint32_t kContextSize = ContextSize(kPredWidth, kPredHeight);
// Maximum context size (for stack reservation).
static constexpr uint32_t kMaxContextSize =
    ContextWithTRSize(kMaxBlockSizePix, kMaxBlockSizePix);
static_assert(kMaxContextSize >=
                  ContextSize(kMaxBlockSizePix, kMaxBlockSizePix),
              "Wrong kMaxContextSize");

// Maximum (inclusive) number of bits to use to represent a frequency in ANS
// quantization.
static constexpr uint32_t kMaxFreqBits = 14;

// Index to Adaptation-Speed mapping.
// TODO(skal): fine-tune, learn.
// Natural photo prefer slower adaptation (speed ~= 2..4)
// Web pictures prefer 3..6
static constexpr uint32_t kAdaptationSpeeds[] = {
    500, 800, 1000, 1300, 1800, 3000, 5000, 65535  /* <~unique symbol */
};
static constexpr uint32_t kNumAdaptationSpeeds =
    sizeof(kAdaptationSpeeds) / sizeof(kAdaptationSpeeds[0]);
static constexpr uint32_t kMinAOMAdaptSpeed = 1024u;
static constexpr uint32_t kMaxAOMAdaptSpeed = 8192u;
static constexpr uint32_t kAOMAdaptSpeedShift = 5;

// quality-factor scale (used for fixed-point precision)
static constexpr uint16_t kQFactorScale = 4;
static inline constexpr uint16_t kQualityToQFactor(float q) {
  return (uint16_t)(q * kQFactorScale);
}
// maximum quality factor (inclusive)
static constexpr uint16_t kQFactorMax = kQualityToQFactor(kMaxLossyQuality);

}  // namespace WP2

//------------------------------------------------------------------------------
// Lossless

namespace WP2L {

// Maximum number of colors in the palette.
// TODO(vrabaud) replace with something related to the image size (ratio?).
static constexpr uint32_t kMaxPaletteSize = 256;
static constexpr size_t kMaxCacheBits = 10;
// Small palettes are stored verbatim.
static constexpr size_t kSmallPaletteLimit = 3;

// Constant defining the special codes for LZ77 distances, for pixels in a
// neighboring window.
// TODO(vrabaud) get rid of it by defining the most common distances/length as
// part of a dictionary.
static constexpr int kCodeToPlaneCodes = 120;

// Maximum length and distance (included) used in LZ77.
// WebP had (1<<12)-1 for length and (1<<20)-120 for distance.
// Any value can be used, it does not have to be a power of two.
static constexpr uint32_t kMaxLZ77Length = (1u << 14);
// TODO(vrabaud) get the code to work for a higher value.
static constexpr uint32_t kMinLZ77Length = 1;
static constexpr uint32_t kMaxLZ77Distance = (1u << 20);

// Min and max number of sampling bits for histograms.
static constexpr uint32_t kHistogramBitsMin = 2;
static constexpr uint32_t kHistogramBitsMax = 9;
// Maximum number of histogram images (sub-blocks). This impacts the compression
// but most images require at the very most 200 histograms. It also impacts the
// amount of RAM used for decoding.
static constexpr uint32_t kMaxHistogramImageSize = 1500;

// Minimum and maximum value of transform bits (to downscale the image).
static constexpr uint32_t kTransformBitsMin = 2;
static constexpr uint32_t kTransformBitsMax = 6;

// TODO(vrabaud) Play with this parameter. This value is the official spec one.
// Window size to consider the vertical mode in Group4 compression.
static constexpr uint32_t kGroup4Window = 3;

// TODO(vrabaud) Play with those values. kRedToBlueMax could be 128, it is just
// that the code is computing as far as this value.
static constexpr uint32_t kRedToBlueMax = 63;
static constexpr uint32_t kGreenToBlueMax = 128;
static constexpr uint32_t kGreenToRedMax = 128;

}  // namespace WP2L

#endif  /* WP2_WP2_FORMAT_CONSTANTS_H_ */
