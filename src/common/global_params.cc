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

#include "src/common/constants.h"
#include "src/common/global_params.h"
#include "src/common/lossy/segment.h"
#include "src/wp2/decode.h"

namespace WP2 {

WP2Status GlobalParams::Write(uint32_t quality_hint, bool image_has_alpha,
                              ANSEnc* const enc) const {
  ANSDebugPrefix prefix(enc, "GlobalParams");
  if (quality_hint == kLosslessQualityHint) {
    assert(type_ == GP_LOSSLESS);
  } else {
    enc->PutRValue(type_, GP_LAST, "gp_type");
  }
  if (image_has_alpha) {
    enc->PutBool(has_alpha_, "has_alpha");
  } else {
    assert(!has_alpha_);
  }
  if (type_ != GP_LOSSLESS) {   // write lossy
    if (has_alpha_) {
      enc->PutBool(maybe_use_lossy_alpha_, "maybe_use_lossy_alpha");
    } else {
      assert(!maybe_use_lossy_alpha_);
    }
    WP2_CHECK_STATUS(transf_.Write(enc));  // write inverse matrix
    enc->PutBool(explicit_segment_ids_, "explicit_segment_ids");
    enc->PutRValue(partition_set_, NUM_PARTITION_SETS, "block_size");
    enc->PutRange(segments_.size(), 1,
        GetMaxNumSegments(explicit_segment_ids_, quality_hint, partition_set_),
        "num_segments_id");
    enc->PutBool(partition_snapping_, "block_size");
    bool use_grain = false;
    for (const Segment& s : segments_) {
      use_grain |= s.grain_.IsUsed();
      if (use_grain) break;
    }
    enc->PutBool(use_grain, "grain");

    const bool has_uv_quant_multiplier =
        (u_quant_multiplier_ != kDefaultQuantMultiplier) ||
        (v_quant_multiplier_ != kDefaultQuantMultiplier);
    enc->PutBool(has_uv_quant_multiplier, "uv_quant_multiplier");
    if (has_uv_quant_multiplier) {
      enc->PutRange(u_quant_multiplier_, 1, kNumQuantMultiplierValues,
                    "uv_quant_multiplier");
      if (u_quant_multiplier_ == kDefaultQuantMultiplier) {
        enc->PutRange((v_quant_multiplier_ > kDefaultQuantMultiplier)
                          ? (v_quant_multiplier_ - 1)
                          : v_quant_multiplier_,
                      1, kNumQuantMultiplierValues - 1, "uv_quant_multiplier");
      } else {
        enc->PutRange(v_quant_multiplier_, 1, kNumQuantMultiplierValues,
                      "uv_quant_multiplier");
      }
    }
    for (const Segment& s : segments_) {
      WP2_CHECK_STATUS(s.WriteHeader(enc, use_grain));
    }
    // TODO(skal): use some fixed-value probas
    enc->PutBool(use_rnd_mtx_, "rnd_mtx");  // TODO(skal): write num_mtx?

    if (maybe_use_lossy_alpha_) {
      enc->PutBool(enable_alpha_filter_, "enable_alpha_filter");
      assert(a_segments_.size() == 1);
      for (const Segment& s : a_segments_) {
        WP2_CHECK_STATUS(s.WriteHeader(enc, /*write_grain=*/false));
      }
    } else {
      assert(!enable_alpha_filter_);
    }
  }
  if (type_ != GP_LOSSY) {
    // write lossless
  }

  return WP2_STATUS_OK;
}

// TODO(skal): use delta-coding instead, for saving more bits
WP2Status GlobalParams::Read(uint32_t quality_hint, bool image_has_alpha,
                             ANSDec* const dec) {
  ANSDebugPrefix prefix(dec, "GlobalParams");
  if (quality_hint == kLosslessQualityHint) {
    type_ = GP_LOSSLESS;
  } else {
    type_ = (Type)dec->ReadRValue(GP_LAST, "gp_type");
  }
  has_alpha_ = image_has_alpha && dec->ReadBool("has_alpha");
  if (type_ != GP_LOSSLESS) {   // read lossy
    maybe_use_lossy_alpha_ =
        has_alpha_ && dec->ReadBool("maybe_use_lossy_alpha");

    WP2_CHECK_STATUS(transf_.Read(dec));
    explicit_segment_ids_ = dec->ReadBool("explicit_segment_ids");
    // segments
    partition_set_ =
        (PartitionSet)dec->ReadRValue(NUM_PARTITION_SETS, "block_size");
    const uint32_t num_segments = dec->ReadRange(
        1,
        GetMaxNumSegments(explicit_segment_ids_, quality_hint, partition_set_),
        "num_segments_id");
    WP2_CHECK_ALLOC_OK(segments_.resize(num_segments));
    partition_snapping_ = dec->ReadBool("block_size");

    const bool use_grain = dec->ReadBool("grain");

    const bool has_uv_quant_multiplier = dec->ReadBool("uv_quant_multiplier");
    if (has_uv_quant_multiplier) {
      u_quant_multiplier_ =
          dec->ReadRange(1, kNumQuantMultiplierValues, "uv_quant_multiplier");
      if (u_quant_multiplier_ == kDefaultQuantMultiplier) {
        v_quant_multiplier_ = dec->ReadRange(1, kNumQuantMultiplierValues - 1,
                                             "uv_quant_multiplier");
        if (v_quant_multiplier_ >= kDefaultQuantMultiplier) {
          ++v_quant_multiplier_;
        }
      } else {
        v_quant_multiplier_ =
            dec->ReadRange(1, kNumQuantMultiplierValues, "uv_quant_multiplier");
      }
    }

    for (Segment& segment : segments_) {
      segment.SetYUVBounds(transf_.GetYUVMin(), transf_.GetYUVMax());
      WP2_CHECK_STATUS(segment.ReadHeader(partition_set_, quality_hint,
                                          u_quant_multiplier_,
                                          v_quant_multiplier_, use_grain, dec));
    }
    use_rnd_mtx_ = dec->ReadBool("rnd_mtx");
    WP2_CHECK_STATUS(InitRndMtxSet());
    WP2_CHECK_STATUS(InitFixedPredictors());
    enable_alpha_filter_ =
        maybe_use_lossy_alpha_ && dec->ReadBool("enable_alpha_filter");
    if (maybe_use_lossy_alpha_) {
      WP2_CHECK_ALLOC_OK(a_segments_.resize(1));
      for (Segment& segment : a_segments_) {
        segment.SetYUVBounds(transf_.GetYUVMin(), transf_.GetYUVMax());
        WP2_CHECK_STATUS(segment.ReadHeader(
            partition_set_, quality_hint, /*u_quant_multiplier=*/0,
            /*v_quant_multiplier=*/0, /*read_grain=*/false, dec));
      }
    }
  }
  if (type_ != GP_LOSSY) {  // read lossless
  }
  return WP2_STATUS_OK;
}

WP2Status GlobalParams::InitFixedPredictors() {
  y_preds_.reset();
  WP2_CHECK_STATUS(y_preds_.Fill(transf_.GetYUVMin(), transf_.GetYUVMax()));
  uv_preds_.reset();
  WP2_CHECK_STATUS(uv_preds_.Fill(transf_.GetYUVMin(), transf_.GetYUVMax()));

  if (maybe_use_lossy_alpha_) {
    WP2_CHECK_STATUS(InitAlphaPredictors(transf_, &a_preds_));
  }

  return WP2_STATUS_OK;
}

WP2Status GlobalParams::InitRndMtxSet() {
  if (use_rnd_mtx_) {
    const uint32_t num_mtx = kMaxNumRndMtx;
    assert(!segments_.empty());
    WP2_CHECK_STATUS(mtx_set_.Init(num_mtx, segments_[0].quant_y_));
  } else {
    mtx_set_.Reset();
  }
  return WP2_STATUS_OK;
}

bool GlobalParams::IsOk() const {
  if (y_preds_.size() > kYPredNum) return false;
  if (uv_preds_.size() > kUVPredNum) return false;
  if (u_quant_multiplier_ < 1 ||
      u_quant_multiplier_ > kNumQuantMultiplierValues)
    return false;
  if (v_quant_multiplier_ < 1 ||
      v_quant_multiplier_ > kNumQuantMultiplierValues)
    return false;
  if (a_preds_.size() > kAPredNum) return false;
  return true;
}

uint16_t GlobalParams::GetMaxAbsDC(Channel channel) const {
  uint16_t max_abs_dc = 0;
  const Vector<Segment>& segments =
      (channel == kAChannel) ? a_segments_ : segments_;
  for (const Segment& s : segments) {
    // Figure out the maximum range.
    max_abs_dc = std::max(max_abs_dc, s.GetMaxAbsDC(channel));
  }
  max_abs_dc = std::min(max_abs_dc, (uint16_t)kMaxDcValue);
  return max_abs_dc;
}

uint32_t GlobalParams::GetMaxDCRange(Channel channel) const {
  // 2 * max_abs_dc + 1 as DC is in [-max_abs_dc, max_abs_dc], zero included.
  return 2u * GetMaxAbsDC(channel) + 1u;
}

void GlobalParams::Reset() {
  type_ = GP_LOSSY;
  has_alpha_ = false;
  transf_.InitYCoCg();
  segments_.reset();
  a_segments_.reset();
  y_preds_.reset();
  uv_preds_.reset();
  maybe_use_lossy_alpha_ = false;
  enable_alpha_filter_ = false;
  a_preds_.reset();
  explicit_segment_ids_ = false;
  partition_set_ = PartitionSet::ALL_SQUARES;
  partition_snapping_ = false;
  use_rnd_mtx_ = false;
  mtx_set_.Reset();
  features_ = nullptr;
  u_quant_multiplier_ = kDefaultQuantMultiplier;
  v_quant_multiplier_ = kDefaultQuantMultiplier;
}

WP2Status GlobalParams::Init() {
  if (type_ != GP_LOSSLESS) {
    WP2_CHECK_STATUS(InitFixedPredictors());
  }
  if (type_ != GP_LOSSY) {
    // TODO(skal): lossless
  }
  return WP2_STATUS_OK;
}

WP2Status GlobalParams::ApplyDecoderConfig(const DecoderConfig& config) {
  const uint32_t grain_amp = 2000u * config.grain_amplitude / 100u;
  for (Segment& s : segments_) {
    s.grain_.y_ = (s.grain_.y_ * grain_amp) >> 8;
    s.grain_.uv_ = (s.grain_.uv_ * grain_amp) >> 8;
  }
  return Init();  // Finish initialization
}

}    // namespace WP2
