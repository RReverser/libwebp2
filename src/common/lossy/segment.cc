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
//   Segment
//
// Author: Skal (pascal.massimino@gmail.com)

#include "src/common/lossy/segment.h"

#include "src/common/constants.h"
#include "src/common/lossy/block.h"
#include "src/dec/symbols_dec.h"
#include "src/enc/symbols_enc.h"
#include "src/utils/ans_utils.h"
#include "src/utils/utils.h"

namespace WP2 {

//------------------------------------------------------------------------------

uint32_t GetMaxNumSegments(bool explicit_segment_ids, uint32_t quality_hint,
                           PartitionSet set) {
  assert(quality_hint <= kMaxLossyQualityHint);
  const uint32_t max_num_segments =
      (explicit_segment_ids ? kMaxNumSegments : kMaxNumSegments / 2);
  // Use fewer segments as quality decreases.
  const uint32_t res = 1u + DivRound((max_num_segments - 1) * quality_hint,
                                     kMaxLossyQualityHint);
  // TODO(vrabaud) Try for other partition sets, only tested on the default one.
  if (set == SOME_RECTS) {
    return std::max(2u, res);
  } else {
    return res;
  }
}

//------------------------------------------------------------------------------
// Grain read/write

bool GrainParams::operator==(const GrainParams& other) const {
  return (y_ == other.y_) &&
         (uv_ == other.uv_) &&
         (cut_y_ == other.cut_y_) &&
         (cut_uv_ == other.cut_uv_);
}

void GrainParams::Write(ANSEnc* const enc) const {
  enc->PutUValue(y_,  4, "grain");
  enc->PutUValue(uv_, 4, "grain");
  enc->PutUValue(cut_y_,  3, "grain");
  enc->PutUValue(cut_uv_, 3, "grain");
}

void GrainParams::Read(ANSDec* const dec) {
  y_ = dec->ReadUValue(4,  "grain");
  uv_ = dec->ReadUValue(4, "grain");
  cut_y_ = dec->ReadUValue(3,  "grain");
  cut_uv_ = dec->ReadUValue(3, "grain");
}

void GrainParams::Reset() {
  y_ = uv_ = 0;
  cut_y_ = cut_uv_ = 0;
}

void GrainParams::Print() const {
  printf("[y=%d uv=%d, cut y=%d uv=%d]\n", y_, uv_, cut_y_, cut_uv_);
}

//------------------------------------------------------------------------------
// Segment

bool Segment::IsMergeableWith(const Segment& other) const {
  const bool same =
      (q_scale_ == other.q_scale_) &&
      (use_quality_factor_ == other.use_quality_factor_) &&
      (use_grain_ == other.use_grain_) &&
      (!use_grain_ || grain_ == other.grain_);
  if (!same) return false;
  for (uint32_t c = 0; c < 4; ++c) {
    for (uint32_t i = 0; i < kNumQuantZones; ++i) {
      if (quant_steps_[c][i] != other.quant_steps_[c][i]) return false;
    }
  }
  return true;
}

WP2Status Segment::WriteHeader(ANSEnc* const enc, bool write_grain) const {
  enc->PutBool(use_quality_factor_, "quant-factor");

  if (use_quality_factor_) {
    enc->PutRValue(quality_factor_, kQFactorMax + 1, "quant-factor");
  } else {
    // TODO(skal): use better delta-coding.
    for (uint32_t i = 0; i < kNumQuantZones; ++i) {
      for (auto c : {kYChannel, kUChannel, kVChannel, kAChannel}) {
        enc->PutRange(quant_steps_[c][i],
                      kMinQuantStep, kMaxQuantStep, "quant-factor");
      }
    }
  }

  if (write_grain && enc->PutUValue(use_grain_, 1, "grain")) grain_.Write(enc);
  return enc->GetStatus();
}

void Segment::StoreEncodingMethods(const CodedBlock& cb,
                                   SymbolManager* const sm,
                                   ANSEncBase* const enc) const {
  for (Channel channel : {kYChannel, kUChannel, kVChannel, kAChannel}) {
    if (channel == kAChannel && !cb.HasLossyAlpha()) continue;
    for (uint32_t tf_i = 0; tf_i < cb.GetNumTransforms(channel); ++tf_i) {
      if (cb.num_coeffs_[channel][tf_i] == 0) {
        // TODO(yguyon): At least one method != kAllZero in block if split
        assert(cb.method_[channel][tf_i] == EncodingMethod::kAllZero);
      } else {
        const uint32_t cluster =
            (channel == kAChannel) ? 0 : (3 * cb.id_ + (uint32_t)channel);
        const Symbol symbol = (channel == kAChannel) ? kSymbolEncodingMethodA
                                                     : kSymbolEncodingMethod;
        if (cb.IsFirstCoeffDC(channel)) {
          sm->Process(cluster, symbol, (uint32_t)cb.method_[channel][tf_i],
                      kChannelStr[channel], enc);
        } else {
          // Use 'max_value' to remove the EncodingMethod::kDCOnly possibility.
          sm->Process(cluster, symbol, (uint32_t)cb.method_[channel][tf_i],
                      /*max_value=*/(uint32_t)EncodingMethod::kMethod1,
                      kChannelStr[channel], enc);
        }
      }
    }
  }
}

WP2Status Segment::ReadHeader(PartitionSet partition_set, uint32_t quality_hint,
                              uint32_t u_quant_multiplier,
                              uint32_t v_quant_multiplier, bool read_grain,
                              ANSDec* const dec) {
  partition_set_ = partition_set;
  const bool use_quality_factor = dec->ReadBool("quant-factor");
  if (use_quality_factor) {
    const uint32_t quality_factor =
        dec->ReadRValue(kQFactorMax + 1, "quant-factor");
    SetQuality(quality_hint, quality_factor,
               u_quant_multiplier, v_quant_multiplier);
  } else {
    for (uint32_t i = 0; i < kNumQuantZones; ++i) {
      for (auto c : {kYChannel, kUChannel, kVChannel, kAChannel}) {
        quant_steps_[c][i] =
            dec->ReadRange(kMinQuantStep, kMaxQuantStep, "quant-factor");
        SetQuantSteps(quant_steps_[c], c);
      }
    }
  }
  FinalizeQuant();

  if (read_grain) {
    use_grain_ = dec->ReadUValue(1, "grain");
    if (use_grain_) grain_.Read(dec);
  } else {
    use_grain_ = false;
  }
  return dec->GetStatus();
}

WP2Status Segment::ReadEncodingMethods(SymbolReader* const sr,
                                       CodedBlock* const cb) const {
  for (Channel channel : {kYChannel, kUChannel, kVChannel, kAChannel}) {
    if (channel == kAChannel && !cb->HasLossyAlpha()) continue;
    for (uint32_t tf_i = 0; tf_i < cb->GetNumTransforms(channel); ++tf_i) {
      if (cb->num_coeffs_[channel][tf_i] == 0) {
        cb->method_[channel][tf_i] = EncodingMethod::kAllZero;
      } else {
        const uint32_t cluster =
            (channel == kAChannel) ? 0 : (3 * cb->id_ + (uint32_t)channel);
        const Symbol symbol = (channel == kAChannel) ? kSymbolEncodingMethodA
                                                     : kSymbolEncodingMethod;
        if (cb->IsFirstCoeffDC(channel)) {
          cb->method_[channel][tf_i] =
              (EncodingMethod)sr->Read(cluster, symbol, kChannelStr[channel]);
        } else {
          // Use 'max_value' to remove the EncodingMethod::kDCOnly possibility.
          int32_t method;
          WP2_CHECK_STATUS(
              sr->TryRead(cluster, symbol,
                          /*max_value=*/(uint32_t)EncodingMethod::kMethod1,
                          kChannelStr[channel], &method));
          cb->method_[channel][tf_i] = (EncodingMethod)method;
        }
      }
    }
  }
  return WP2_STATUS_OK;
}

//------------------------------------------------------------------------------

uint16_t Segment::GetMaxAbsDC(Channel channel) const {
  const QuantMtx& mtx = GetQuant(channel);
  uint16_t max_abs_dc = 0;
  // TODO(vrabaud) only use the TrfSizes defined by -ps.
  for (const TrfSize tdim : kAllTrfSizes) {
    max_abs_dc = std::max(max_abs_dc, mtx.GetMaxAbsResDC(tdim));
  }
  return max_abs_dc;
}

//------------------------------------------------------------------------------


static constexpr uint32_t kQuantizerFrac = 16u;   // 4bit fractional precision
// Table mapping quality factor to quantizer steps, discretized using a reduced
// set of values. Linear interpolation is performed for in-between values.
// The first value (highest quality) is always 128. Subsequent ones are
// exponentially increasing.
static constexpr uint32_t kNumQuantizers = 32u;
static constexpr uint16_t kQuantizers[kNumQuantizers + 1] = {
  kMinQuantStep,
  592, 637, 685, 738, 794, 855, 920, 991, 1066, 1148, 1236, 1330, 1432,
  1541, 1659, 1786, 1922, 2069, 2227, 2397, 2580, 2777, 2990, 3218, 3464,
  3729, 4013, 4320, 4650, 5005, 5388,
  kMaxQuantStep
};
STATIC_ASSERT_ARRAY_SIZE(kQuantizers, kNumQuantizers+1);

static uint32_t InterpolateQuant(int32_t qf, float mult) {
  qf = Clamp(qf, 0, (int32_t)kQFactorMax);
  uint32_t idx = qf * kNumQuantizers * kQuantizerFrac / kQFactorMax;
  const uint32_t eps = idx % kQuantizerFrac;
  idx /= kQuantizerFrac;
  assert(idx <= kNumQuantizers);
  uint32_t q0 = kQuantizers[idx];
  if (eps > 0) {   // interpolate
    q0 = q0 * (kQuantizerFrac - eps) + kQuantizers[idx + 1] * eps;
    q0 = (q0 + kQuantizerFrac / 2) / kQuantizerFrac;
  }
  const int32_t q = std::lround(q0 * mult);
  return (uint32_t)std::max(std::min(q, (int32_t)kMaxQuantStep),
                                        (int32_t)kMinQuantStep);
}

void Segment::SetQuality(uint32_t quality_hint, uint32_t quality_factor,
                         uint32_t u_quant_multiplier,
                         uint32_t v_quant_multiplier) {
  // The following ad-hoc formulae maps the [0-kQFactorMax] range of quality
  // factor to approximately 0.05 bpp - 3 bpp usable range. Note that
  // quality_factor 0 is the highest quality (= smaller quant-step).
  quality_factor = std::min(quality_factor, (uint32_t)kQFactorMax);
  const int32_t q_y = quality_factor;
  const int32_t q_u = quality_factor + u_quant_multiplier - kNeutralMultiplier;
  const int32_t q_v = quality_factor + v_quant_multiplier - kNeutralMultiplier;
  const int32_t q_a = quality_factor;
  // Quantize DC less than AC for UV, especially at lower qualities.
  constexpr float kMinScale = 0.15f;
  constexpr float kMaxScale = 0.70f;
  const float qf =
      1.f * std::min(quality_hint, kMaxLossyQualityHint) / kMaxLossyQualityHint;
  const float uv_dc_scale = kMinScale + (kMaxScale - kMinScale) * qf;
  const float kUVScales[kNumQuantZones] = {uv_dc_scale, 1.0, 1.0, 1.0};

  for (uint32_t i = 0; i < kNumQuantZones; ++i) {
    quant_steps_[kYChannel][i] = InterpolateQuant(q_y, 1.0f);
    quant_steps_[kUChannel][i] = InterpolateQuant(q_u, kUVScales[i]);
    quant_steps_[kVChannel][i] = InterpolateQuant(q_v, kUVScales[i]);
    quant_steps_[kAChannel][i] = InterpolateQuant(q_a, 1.0f);
  }
  quality_factor_ = quality_factor;
  use_quality_factor_ = true;
}

void Segment::SetQuantSteps(uint16_t quants[kNumQuantZones], Channel channel) {
  for (uint32_t i = 0; i < kNumQuantZones; ++i) {
    if (quants[i]) {
      quant_steps_[channel][i] = std::min((uint32_t)quants[i], kMaxQuantStep);
    }
  }
  use_quality_factor_ = false;
}

void Segment::GetQuantSteps(uint16_t quants[kNumQuantZones],
                            Channel channel) const {
  for (uint32_t i = 0; i < kNumQuantZones; ++i) {
    quants[i] = quant_steps_[channel][i];
  }
}

WP2Status Segment::AllocateForEncoder() {
  WP2_CHECK_STATUS(quant_y_.Allocate());
  WP2_CHECK_STATUS(quant_u_.Allocate());
  WP2_CHECK_STATUS(quant_v_.Allocate());
  WP2_CHECK_STATUS(quant_a_.Allocate());
  return WP2_STATUS_OK;
}

void Segment::FinalizeQuant(float lambda) {
  assert(q_scale_ > 0.f);   // must have been set prior to calling this
  quant_y_.Init(max_residual_, q_scale_, partition_set_,
                quant_steps_[kYChannel]);
  quant_u_.Init(max_residual_, q_scale_, partition_set_,
                quant_steps_[kUChannel]);
  quant_v_.Init(max_residual_, q_scale_, partition_set_,
                quant_steps_[kVChannel]);
  // Max residual value is different for alpha since it doesn't go through the
  // CSP transform. Alpha is between 0 and kAlphaMax. After subtracting
  // predicted values from actual values, the residuals are between -kAlphaMax
  // and kAlphaMax. So the max absolute value is still kAlphaMax.
  quant_a_.Init(/*max_residual=*/kAlphaMax, q_scale_a_, partition_set_,
                quant_steps_[kAChannel]);
  if (lambda > 0) {
    lambda_f_ = lambda;
    quant_y_.SetLambda(lambda / q_scale_);
    quant_u_.SetLambda(lambda / q_scale_);
    quant_v_.SetLambda(lambda / q_scale_);
    quant_a_.SetLambda(lambda / q_scale_a_);
  }
}

void Segment::AdjustQuantSteps() {
  if (quant_y_.qstats != nullptr) {
    uint16_t quant_steps[kNumQuantZones];
    GetQuantSteps(quant_steps, kYChannel);
    quant_y_.qstats->AdjustQuantSteps(quant_steps);
    SetQuantSteps(quant_steps, kYChannel);
  }
  // TODO(skal): quant_uv_

  FinalizeQuant();
}

WP2Status GlobalParams::AssignQuantizations(const EncoderConfig& config) {
  const uint32_t num_segments = segments_.size();
  assert(num_segments <= (uint32_t)config.segments);
  assert(num_segments <= kMaxNumSegments);
  uint32_t quants[kMaxNumSegments];
  const float quant_factor0 =
      kQFactorMax * (kMaxLossyQuality - config.quality) / kMaxLossyQuality;
  for (uint32_t idx = 0; idx < num_segments; ++idx) {
    const Segment& s = segments_[idx];
    float quant = config.segment_factors[idx];
    if (quant == 0.f) {
      const float adjust = std::min(s.risk_, 1.f);

      // Now we map from 'adjust' to 'quant'.
      // General slope of the curve.
      const float amp = 1.0f * Clamp(config.sns / 100.f, 0.f, 1.f);
      // Exponent for 'adjust' (1 = linear).
      const float exponent = 1.0;
      // Compute slope and intercept values so that an 'adjust' value
      // of 'fixedX' will have a value of 'fixedY'.
      const float fixedX = 0.5;
      const float fixedY = 1.;
      const float intercept = fixedY + amp;
      const float slope = -amp / std::pow(fixedX, exponent);
      quant = quant_factor0 * (intercept + slope * std::pow(adjust, exponent));
    } else {
      quant = kQFactorMax * (kMaxLossyQuality - quant) / kMaxLossyQuality;
    }
    quants[idx] = (uint32_t)std::lround(quant);
  }
  for (uint32_t idx = 0; idx < num_segments; ++idx) {
    WP2_CHECK_STATUS(segments_[idx].AllocateForEncoder());
    segments_[idx].SetQuantizationFactor(
        transf_, GetQualityHint(config.quality), u_quant_multiplier_,
        v_quant_multiplier_, quants[idx], 1.f, config.partition_set);
  }
  return WP2_STATUS_OK;
}

WP2Status GlobalParams::AssignAlphaQuantizations(const YUVPlane& yuv,
                                                 const EncoderConfig& config) {
  WP2AlphaInit();
  has_alpha_ = false;
  if (!yuv.A.IsEmpty()) {
    for (uint32_t y = 0; y < yuv.A.h_; ++y) {
      const int16_t* const row = yuv.A.Row(y);
      has_alpha_ = has_alpha_ || WP2HasOtherValue16b(row, yuv.A.w_, kAlphaMax);
      if (has_alpha_) break;
    }
  }
  maybe_use_lossy_alpha_ =
      has_alpha_ && (config.alpha_quality <= kMaxLossyQuality);
  if (maybe_use_lossy_alpha_) {
    // TODO(maryla): we could detect automatically cases where the alpha
    // filter should not be applied.
    enable_alpha_filter_ = config.enable_alpha_filter;
    WP2_CHECK_ALLOC_OK(a_segments_.resize(1));
    const uint32_t max_quality = kMaxLossyQuality;
    const uint32_t quant_factor0 =
        (uint32_t)std::lround(
            kQFactorMax * (max_quality - config.alpha_quality) / max_quality);
    const float lambda = std::max(85.f, -42.f * config.alpha_quality + 3400.f);
    for (Segment& s : a_segments_) {
      WP2_CHECK_STATUS(s.AllocateForEncoder());
      s.SetQuantizationFactor(
          transf_, GetQualityHint(config.alpha_quality),
          /*u_quant_multiplier=*/0, /*v_quant_multiplier=*/0, quant_factor0,
          1.f + lambda / kMaxLossyQuality, config.partition_set);
    }
  }
  return WP2_STATUS_OK;
}

void Segment::SetQuantizationFactor(const CSPTransform& transform,
                                    uint32_t quality_hint,
                                    uint32_t u_quant_multiplier,
                                    uint32_t v_quant_multiplier,
                                    uint32_t qfactor, float lambda_mult,
                                    PartitionSet partition_set) {
  partition_set_ = partition_set;
  SetYUVBounds(transform.GetYUVMin(), transform.GetYUVMax());
  SetQuality(quality_hint, qfactor, u_quant_multiplier, v_quant_multiplier);
  FinalizeQuant(lambda_mult);
}

//------------------------------------------------------------------------------

}   // namespace WP2
