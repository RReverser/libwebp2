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
// Lossless encoder: internal header.
//
// Author: Vincent Rabaud (vrabaud@google.com)

#ifndef WP2_ENC_LOSSLESS_LOSSLESSI_ENC_H_
#define WP2_ENC_LOSSLESS_LOSSLESSI_ENC_H_

#include "src/common/symbols.h"
#include "src/dsp/lossless/lossless.h"
#include "src/enc/lossless/backward_references_enc.h"
#include "src/enc/lossless/histogram_enc.h"
#include "src/enc/lossless/palette.h"
#include "src/utils/ans.h"
#include "src/utils/ans_utils.h"
#include "src/wp2/encode.h"
#include "src/wp2/format_constants.h"

namespace WP2 {
class SymbolManager;
class SymbolRecorder;
class SymbolWriter;
}

namespace WP2L {

typedef enum {
  kEncoderNone = 0,
  kEncoderARGB,
  kEncoderNearLossless,  // Not implemented yet.
  kEncoderPalette
} EncoderARGBContent;

struct CrunchConfig;
class Encoder {
 public:
  Encoder(const WP2::EncoderConfig& config, const WP2::ArgbBuffer& picture,
          bool has_alpha);
  WP2Status Allocate();
  WP2Status AllocateTransformBuffer(const CrunchConfig& config, uint32_t width,
                                    uint32_t height);

  const WP2::EncoderConfig& config_;  // user configuration and parameters
  const WP2::ArgbBuffer& pic_;        // input picture (tile view)
  bool has_alpha_;  // false if the whole image is opaque (not only this tile)

  uint16_t* argb_;                       // Transformed argb image data.
  EncoderARGBContent argb_content_;      // Content type of the argb buffer.
  uint16_t* argb_scratch_;               // Scratch memory for argb rows
                                         // (used for prediction).
  int16_t* transform_data_;              // Scratch memory for transform data.
  WP2::Vector_u16 transform_mem_;        // Currently allocated memory.

  // Encoding parameters derived from quality parameter.
  uint32_t histo_bits_;
  uint32_t transform_bits_;
  LosslessSymbolsInfo symbols_info_;

  Palette palette_;

  // Some 'scratch' (potentially large) objects.
  BackwardRefsPool ref_pool_;   // Backward Refs array for temporaries.
  HashChain hash_chain_;        // HashChain data for constructing
                                // backward references.
};

//------------------------------------------------------------------------------
// internal functions. Not public.

WP2Status EncodeImageNoClusters(
    WP2::ANSEncBase* const enc, WP2::ANSDictionaries* const dicts,
    const int16_t* const argb, HashChain* const hash_chain,
    BackwardRefsPool* const ref_pool, uint32_t width, uint32_t height,
    const LosslessSymbolsInfo& symbols_info, int speed);

// Optional extra info that can be exported by EncodeImage.
// Callers should allocate themselves the vectors that they wish filled in.
// Empty vectors will not be filled in.
struct EncodeInfo {
  WP2_NO_DISCARD bool CopyFrom(const EncodeInfo& other);

  // The number of tokens in the ANSEnc at which point each line has finished
  // being written is returned in it. If non-empty, its size should match the
  // height of the image.
  WP2::Vector_u32 line_tokens;
  // Cost in bits of storing each pixel of the image (not including bits of
  // headers like palettes, ans stats...), in row major order. If non-empty,
  // its size should match the number of pixels of the image.
  WP2::Vector_f bits_per_pixel;
};

// Encodes the picture.
// 'has_alpha' might be true for another tile even if this one does not.
WP2Status EncodeImage(const WP2::EncoderConfig& config,
                      const WP2::ArgbBuffer& picture, bool has_alpha,
                      WP2::ANSEnc* const enc,
                      EncodeInfo* const encode_info = nullptr);

// Stores the pixels of the image, without caring about the header.
WP2Status StorePixels(uint32_t width, uint32_t height, uint32_t histo_bits,
                      const BackwardRefs& refs,
                      const uint16_t* const histogram_symbols,
                      WP2::ANSEncBase* const enc, WP2::SymbolManager* const sw,
                      EncodeInfo* const encode_info = nullptr);
// Same as above except for an image were there are no entropy clusters.
WP2Status StorePixels(uint32_t width, const BackwardRefs& refs,
                      WP2::ANSEncBase* const enc, WP2::SymbolManager* const sw);
// Stores the different symbol headers and sets up 'sw' with the right symbol
// info.
WP2Status WriteHeaders(const WP2::SymbolRecorder& recorder,
                       const LosslessSymbolsInfo& symbols_info,
                       uint32_t num_pixels, WP2::ANSEncBase* const enc,
                       WP2::ANSDictionaries* const dicts,
                       WP2::SymbolWriter* const sw);

// Converts a distance to a code taking into account the proximity to the
// current pixel.
uint32_t DistanceToPlaneCode(uint32_t width, uint32_t dist);

//------------------------------------------------------------------------------
// Image transforms in predictor_enc.cc

WP2Status ResidualImage(uint32_t width, uint32_t height, uint32_t bits,
                        uint32_t channel_bits, int speed, uint16_t* const argb,
                        bool has_alpha, uint16_t* const argb_scratch,
                        int16_t* const image);

WP2Status ColorSpaceTransform(uint32_t width, uint32_t height, uint32_t bits,
                              WP2SampleFormat format, uint32_t speed,
                              uint16_t* const argb, int16_t* image);

//------------------------------------------------------------------------------

// Set of parameters to be used in each iteration of the cruncher.
static constexpr uint32_t kCrunchLZ77Max = 2;
struct CrunchConfig {
  uint32_t lz77s_types_to_try[kCrunchLZ77Max];
  uint32_t lz77s_types_to_try_size;
  // Encoding parameters derived from image characteristics.
  bool use_cross_color;
  bool use_subtract_green;
  bool use_predict;
  bool use_palette;
  bool use_group4;
  // Use a MoveToFrontCache in group4 to signal new colors.
  bool group4_use_move_to_front;
};

static constexpr uint32_t kCrunchNumMax = 7;

WP2Status EncoderAnalyze(Encoder* const encoder,
                         CrunchConfig crunch_configs[kCrunchNumMax],
                         uint32_t* const crunch_configs_size);

//------------------------------------------------------------------------------
// Group 4.

WP2Status Group4Encode(const uint16_t* const argb, uint32_t width,
                       uint32_t height, uint32_t num_colors,
                       bool use_move_to_front, WP2::ANSEncBase* const enc,
                       EncodeInfo* const encode_info);

}  // namespace WP2L

//------------------------------------------------------------------------------

#endif  /* WP2_ENC_LOSSLESS_LOSSLESSI_ENC_H_ */
