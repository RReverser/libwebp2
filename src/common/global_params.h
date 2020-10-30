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
//  Global parameters used by lossy / lossless tiles
//
// Author: Skal (pascal.massimino@gmail.com)
//
#ifndef WP2_COMMON_GLOBAL_PARAMS_H_
#define WP2_COMMON_GLOBAL_PARAMS_H_

#include "src/wp2/base.h"
#include "src/wp2/format_constants.h"
#include "src/common/lossy/predictor.h"
#include "src/common/lossy/rnd_mtx.h"
#include "src/common/lossy/segment.h"
#include "src/utils/vector.h"
#include "src/utils/csp.h"

namespace WP2 {

class ANSEnc;
class ANSDec;
class FeatureMap;
struct DecoderConfig;

//------------------------------------------------------------------------------
// Global parameter

// Quant multiplier allows adjusting U/V quantization compared to Y.
// With 4 bits we allow multipliers between 1/8 and 2 in increments of 1/8.
static constexpr uint32_t kQuantMultiplierBits = 4;
static constexpr uint32_t kNumQuantMultiplierValues = 1 << kQuantMultiplierBits;
// This value means multiply by 1 basically.
static constexpr uint32_t kNeutralMultiplier =
    (1 << (kQuantMultiplierBits - 1));
static constexpr uint32_t kDefaultQuantMultiplier = 7;  // 0.875

class GlobalParams {
 public:
  WP2Status Write(uint32_t quality_hint, bool image_has_alpha,
                  ANSEnc* const enc) const;
  WP2Status Read(uint32_t quality_hint, bool image_has_alpha,
                 ANSDec* const dec);

  void Reset();

  // Permanently modify the params according to the decoder config. Call Init()
  // afterward. Only used by the decoder.
  WP2Status ApplyDecoderConfig(const DecoderConfig& config);

  // finish initializing the satellite data
  WP2Status Init();

 public:
  // TODO(skal): should we have a GP_NONE for super-small images with default
  //             global params?
  enum Type { GP_LOSSY = 0, GP_LOSSLESS, GP_BOTH, GP_LAST };
  Type type_ = GP_LOSSY;

  // True if the frame contains alpha (either lossless or lossy).
  // Can only be true if 'BitstreamFeatures::has_alpha' is true as well.
  bool has_alpha_ = false;

  // lossy part
  CSPTransform transf_;

  // True if segment ids are signaled per block. Otherwise they are
  // deduced from the block sizes.
  bool explicit_segment_ids_ = false;

  PartitionSet partition_set_ = PartitionSet::ALL_SQUARES;
  bool partition_snapping_ = false;

  bool use_rnd_mtx_ = false;
  RndMtxSet mtx_set_;

  // UV quantization multipliers. [1 .. 16]
  uint32_t u_quant_multiplier_ = kDefaultQuantMultiplier;
  uint32_t v_quant_multiplier_ = kDefaultQuantMultiplier;

  YPredictors y_preds_;
  UVPredictors uv_preds_;
  APredictors a_preds_;
  Vector<Segment> segments_;

  // True if the image might contain lossy alpha.
  bool maybe_use_lossy_alpha_ = false;
  bool enable_alpha_filter_ = false;  // For lossy alpha only.
  Vector<Segment> a_segments_;
  FeatureMap* features_ = nullptr;  // only used by the encoder

 public:   // lossless part
  // TODO(skal)

 public:
  const Predictors& predictors(Channel channel) const {
    if (channel == kAChannel) return a_preds_;
    if (channel != kYChannel) return uv_preds_;
    return y_preds_;
  }

  // initialize y_preds_ / uv_preds_ / a_preds_
  WP2Status InitFixedPredictors();

  WP2Status InitRndMtxSet();

  bool IsOk() const;

  uint16_t GetMaxAbsDC(Channel channel) const;
  uint32_t GetMaxDCRange(Channel channel) const;

  // These functions are responsible for mapping the segment's risk_ factor
  // to a set of quant factors, and applying them on segments.
  WP2Status AssignQuantizations(const EncoderConfig& config);
  WP2Status AssignAlphaQuantizations(const YUVPlane& yuv,
                                     const EncoderConfig& config);
};

//------------------------------------------------------------------------------

}  // namespace WP2

#endif  // WP2_COMMON_GLOBAL_PARAMS_H_
