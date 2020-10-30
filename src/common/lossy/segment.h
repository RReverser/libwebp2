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
//  Segment
//
// Author: Skal (pascal.massimino@gmail.com)
//
#ifndef WP2_COMMON_LOSSY_SEGMENT_H_
#define WP2_COMMON_LOSSY_SEGMENT_H_

#include <cmath>

#include "src/common/lossy/quant_mtx.h"
#include "src/utils/ans_utils.h"
#include "src/wp2/base.h"
#include "src/wp2/format_constants.h"

namespace WP2 {

struct EncoderConfig;
class CSPTransform;
class CodedBlock;
class SymbolManager;
class SymbolReader;
class GlobalParams;

//------------------------------------------------------------------------------

// Maximum number of segments.
uint32_t GetMaxNumSegments(bool explicit_segment_ids, uint32_t quality_hint,
                           PartitionSet set);

//------------------------------------------------------------------------------
// Grain

class GrainParams {
 public:
  void Write(ANSEnc* const enc) const;
  void Read(ANSDec* const dec);

  bool operator==(const GrainParams& other) const;
  bool IsUsed() const { return (y_ > 0 || uv_ > 0); }
  void Reset();

  void Print() const;  // debug

 public:
  uint8_t y_ = 0, uv_ = 0;          // grain level [0=off, 15=full]
  uint8_t cut_y_ = 0, cut_uv_ = 0;  // grain cut-off [0..7]
};

//------------------------------------------------------------------------------
// Segment

class Segment {
 public:
  WP2Status WriteHeader(ANSEnc* const enc, bool write_grain) const;
  WP2Status ReadHeader(PartitionSet partition_set, uint32_t quality_hint,
                       uint32_t u_quant_multiplier, uint32_t v_quant_multiplier,
                       bool read_grain, ANSDec* const dec);

 public:
  // quantization
  QuantMtx quant_y_;
  QuantMtx quant_u_;
  QuantMtx quant_v_;
  QuantMtx quant_a_;

  const QuantMtx& GetQuant(Channel channel) const {
    switch (channel) {
      case kYChannel:
        return quant_y_;
      case kUChannel:
        return quant_u_;
      case kVChannel:
        return quant_v_;
      default:
        assert(channel == kAChannel);
        return quant_a_;
    }
  }

  // Returns the maximum absolute value of the DC coefficient over all possible
  // block dimensions.
  uint16_t GetMaxAbsDC(Channel channel) const;

  uint16_t quant_steps_[4][kNumQuantZones] = {{0}, {0}, {0}, {0}};  // Y/U/V/A
  float lambda_f_ = 0.f;
  bool use_quality_factor_ = true;
  uint16_t quality_factor_ = kQualityToQFactor(75.f);

  // Sets the range of YUV values to [yuv_min, yuv_max] inclusive.
  // q_scale (quantization scale) is updated based on it.
  void SetYUVBounds(int16_t yuv_min, int16_t yuv_max) {
    // Because residuals are the subtraction of the actual and the predicted
    // values, the range is doubled.
    max_residual_ = std::max(std::abs(yuv_min), std::abs(yuv_max)) * 2;
    const uint32_t range = 1 + yuv_max - yuv_min;
    // Adjust quantization based on the actual yuv/alpha range, compared to
    // the max yuv range which is what the quantizer is tuned for by default.
    q_scale_ = (float)(1 << (kMaxYuvBits + 1)) / range;
    q_scale_a_ = (float)(1 << (kMaxYuvBits + 1)) / (kAlphaMax + 1);
  }

  // Fills quant_steps_[][] with flat values based on 'quality_hint' and
  // 'quality_factor' in [0,kQFactorMax]. Note 'quality_factor' is backwards:
  // the higher it is, the more the image is compressed.
  void SetQuality(uint32_t quality_hint, uint32_t quality_factor,
                  uint32_t u_quant_multiplier, uint32_t v_quant_multiplier);
  // Fills the corresponding quant_steps_[]
  void SetQuantSteps(uint16_t quants[kNumQuantZones], Channel channel);
  void GetQuantSteps(uint16_t quants[kNumQuantZones], Channel channel) const;

  // Initializes the quant_mtx based on quant_steps_[].
  // q_scale must be set before calling this.
  void FinalizeQuant(float lambda = 0.f);

  WP2Status AllocateForEncoder();   // must be called before FinalizeQuant()

  // Adjusts the quantization steps if some stats were recorded.
  // TODO(skal): right now it only adjusts the Y-Channel.
  void AdjustQuantSteps();

  // Sets the YUV bounds, the quality and the quant mtx based on 'qfactor'.
  // 'qfactor' should be in [0, kFactorMax] range.
  void SetQuantizationFactor(const CSPTransform& transform,
                             uint32_t quality_hint, uint32_t u_quant_multiplier,
                             uint32_t v_quant_multiplier, uint32_t qfactor,
                             float lambda_mult, PartitionSet partition_set);

  // Returns true if this segment is "close enough" to 'other' that they can be
  // merged (some fields might still be different).
  bool IsMergeableWith(const Segment& other) const;
  WP2Status CopyFrom(const Segment& other);

  // Filtering and noise/grain generation.
  bool use_grain_ = false;
  GrainParams grain_;

  // Encoding methods.
  WP2Status ReadEncodingMethods(SymbolReader* const sr,
                                CodedBlock* const cb) const;

  void StoreEncodingMethods(const CodedBlock& cb, SymbolManager* const sm,
                            ANSEncBase* const enc) const;

  // Basic characteristics of a segment, which can be use to determine
  // quantization and other compression parameters. Used by encoder only.
  // 'risk_' is heuristic a value between 0 (high risk of visible artifacts)
  // and 1 (low risk of artifacts). The value impacts the choice of quantizers.
  // 'risk_class_' is a heuristic index mapping to risk_ values.
  float risk_ = 0.f;
  int risk_class_ = 0;  // class 0 = most difficult

  static constexpr float kDefaultAvgLambda = 1.f;
  float avg_lambda_ = kDefaultAvgLambda;  // complexity multiplier for lambda

 private:
  float q_scale_ = 0.f;  // y/u/v quantization precision
  float q_scale_a_ = 0.f;  // a quantization precision
  uint32_t max_residual_;
  PartitionSet partition_set_ = ALL_RECTS;
};

//------------------------------------------------------------------------------

// Colors of segments for visual debugging.
extern const Argb32b kSegmentColors[kMaxNumSegments];

//------------------------------------------------------------------------------

}  // namespace WP2

#endif  // WP2_COMMON_LOSSY_SEGMENT_H_
