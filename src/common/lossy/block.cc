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
//   Coded blocks
//
// Author: Skal (pascal.massimino@gmail.com)

#include "src/common/lossy/block.h"

#include <algorithm>
#include <cassert>

#include "src/common/lossy/block_size.h"
#include "src/common/lossy/rnd_mtx.h"
#include "src/common/lossy/transforms.h"
#include "src/common/symbols.h"
#include "src/dsp/dsp.h"
#include "src/dsp/math.h"
#include "src/enc/wp2_enc_i.h"
#include "src/utils/ans.h"
#include "src/utils/front_mgr.h"
#include "src/utils/utils.h"
#include "src/wp2/format_constants.h"

namespace WP2 {

//------------------------------------------------------------------------------

WP2Status Counters::Init(const SymbolRecorder& recorder) {
  recorder_ = &recorder;

  transform_.reset(new (WP2Allocable::nothrow) SymbolCounter(recorder_));
  WP2_CHECK_ALLOC_OK(transform_ != nullptr);
  WP2_CHECK_STATUS(transform_->Allocate(
      {kSymbolHasCoeffs, kSymbolTransform, kSymbolSplitTransform}));

  residuals_.reset(new (WP2Allocable::nothrow) SymbolCounter(recorder_));
  WP2_CHECK_ALLOC_OK(residuals_ != nullptr);
  WP2_CHECK_STATUS(residuals_->Allocate(
      {kSymbolDC, kSymbolResidualUseBounds, kSymbolResidualBound1IsX,
       kSymbolResidualUseBound2, kSymbolResidualIsZero, kSymbolResidualIsOne,
       kSymbolResidualIsTwo, kSymbolResidualEndOfBlock,
       kSymbolResidualHasOnlyOnesLeft}));

  residuals_aom_.reset(new (WP2Allocable::nothrow) SymbolCounter(recorder_));
  WP2_CHECK_ALLOC_OK(residuals_aom_ != nullptr);
  WP2_CHECK_STATUS(residuals_aom_->Allocate(
      {kSymbolTransform, kAOMAllZero, kAOMEOBPT16, kAOMEOBPT32, kAOMEOBPT64,
       kAOMEOBPT128, kAOMEOBPT256, kAOMEOBPT512, kAOMEOBPT1024, kAOMEOBExtra,
       kAOMCoeffBaseEOB, kAOMCoeffBase, kAOMCoeffBaseRange, kAOMDCSign}));

  predictor_.reset(new (WP2Allocable::nothrow) SymbolCounter(recorder_));
  WP2_CHECK_ALLOC_OK(predictor_ != nullptr);
  WP2_CHECK_STATUS(
      predictor_->Allocate({kSymbolModeY, kSymbolModeUV, kSymbolModeA}));
  return WP2_STATUS_OK;
}

WP2Status Counters::CopyFrom(const Counters& other,
                             const SymbolRecorder& recorder) {
  if (&recorder != recorder_) {
    WP2_CHECK_STATUS(Init(recorder));
  } else {
    // Init() was already called.
    assert(transform_ != nullptr && residuals_ != nullptr &&
           residuals_aom_ != nullptr && predictor_ != nullptr);
  }
  WP2_CHECK_STATUS(transform_->CopyFrom(*other.transform_));
  WP2_CHECK_STATUS(residuals_->CopyFrom(*other.residuals_));
  WP2_CHECK_STATUS(residuals_aom_->CopyFrom(*other.residuals_aom_));
  WP2_CHECK_STATUS(predictor_->CopyFrom(*other.predictor_));
  return WP2_STATUS_OK;
}

//------------------------------------------------------------------------------

const int16_t* ContextCache::GetLarge(const CodedBlock& cb, Channel channel,
                                      bool fill_in, bool extend_right) {
  if (!is_computed_[channel][fill_in][extend_right]) {
    cb.GetContext(channel, fill_in, extend_right,
                  large_[channel][fill_in][extend_right]);
    is_computed_[channel][fill_in][extend_right] = true;
  }
  return large_[channel][fill_in][extend_right];
}

void ContextCache::Reset() {
  for (Channel channel : {kYChannel, kUChannel, kVChannel, kAChannel}) {
    for (bool fill_in : {false, true}) {
      for (bool extend_right : {false, true}) {
        is_computed_[channel][fill_in][extend_right] = false;
      }
    }
  }
}

//------------------------------------------------------------------------------

void CodedBlock::SetRange(int16_t yuv_min, int16_t yuv_max) {
  yuv_min_ = yuv_min;
  yuv_max_ = yuv_max;
}

void CodedBlock::SetDim(const Block& block, const FrontMgrBase& mgr) {
  blk_ = block;
  left_occupancy_ = mgr.GetOccupancy((int)blk_.x() - 1);
  right_occupancy_ = mgr.GetOccupancy(blk_.x() + blk_.w());
  top_context_extent_ = mgr.GetRightContextExtent(blk_);
}

void CodedBlock::SetDimDefault(const Block& block, bool full_left_ctx) {
  blk_ = block;
  left_occupancy_ =
      blk_.y() + (blk_.x() > 0 ? (full_left_ctx ? blk_.h() : 1) : 0);
  right_occupancy_ = blk_.y();
  top_context_extent_ = 0u;
}

void CodedBlock::SetXY(const FrontMgrBase& mgr, uint32_t x, uint32_t y) {
  blk_.SetXY(x, y);
  left_occupancy_ = mgr.GetOccupancy((int)blk_.x() - 1);
  right_occupancy_ = mgr.GetOccupancy(blk_.x() + blk_.w());
  top_context_extent_ = mgr.GetRightContextExtent(blk_);
}

void CodedBlock::SetSrcInput(const YUVPlane& in) {
  WP2_ASSERT_STATUS(in_.SetView(in, blk_.rect_pix()));
}

void CodedBlock::SetContextInput(const YUVPlane& in,
                                 ContextCache* const context_cache) {
  WP2_ASSERT_STATUS(predict_from_.SetView(in, blk_.rect_pix()));
  if (context_cache != nullptr) context_cache_ = context_cache;
  ResetContextCache();
}

void CodedBlock::SetReconstructedOutput(YUVPlane* const out) {
  WP2_ASSERT_STATUS(out_.SetView(*out, blk_.rect_pix()));
}

//------------------------------------------------------------------------------

void CodedBlock::ExtractFrom(const YUVPlane& yuv, Channel channel) const {
  Plane16 src, dst;
  WP2_ASSERT_STATUS(src.SetView(yuv.GetChannel(channel), blk_.rect_pix()));
  WP2_ASSERT_STATUS(dst.SetView(out_.GetChannel(channel),
                                {0, 0, blk_.w_pix(), blk_.h_pix()}));
  WP2_ASSERT_STATUS(dst.Copy(src, /*resize_if_needed=*/false));
}

//------------------------------------------------------------------------------

void CodedBlock::GetOccupancy(int8_t* const left, int8_t* const right,
                              int8_t* const top_context_extent) const {
  *left = std::max(-1, (int)(left_occupancy_ - blk_.y()));
  *right = std::max(-1, (int)(right_occupancy_ - blk_.y()));
  if (top_context_extent != nullptr) {
    *top_context_extent = top_context_extent_;
  }
}

void CodedBlock::SetContextCache(ContextCache* const context_cache) {
  context_cache_ = context_cache;
}

void CodedBlock::ResetContextCache() const {
  assert(context_cache_ != nullptr);
  context_cache_->Reset();
}

const int16_t* CodedBlock::GetContext(Channel channel, bool fill_in,
                                      bool extend_right) const {
  assert(context_cache_ != nullptr);
  return context_cache_->GetLarge(*this, channel, fill_in, extend_right);
}

// Macro replicating AtUnclamped for speed-up.
#define CONTEXT_SRC_AT_UNCLAMPED(X, Y) \
  *(context_src_data + (X) + (Y) * context_src_step)

void CodedBlock::GetContext(Channel channel, bool fill_in, bool extend_right,
                            int16_t context[]) const {
  const uint32_t w = blk_.w_pix(), h = blk_.h_pix();
  const Plane16& context_src = predict_from_.GetChannel(channel);
  // absolute position within dst_ picture
  const uint32_t y_pix = blk_.y_pix();
  // implement all the neighbouring sample logic:
  const bool have_top = (y_pix > 0);
  int8_t y_left_occ, y_right_occ, right_extent;
  GetOccupancy(&y_left_occ, &y_right_occ, &right_extent);
  // Number of available context pixels on the left and right sides.
  const int left_context_size =
      (std::max(0, (int)y_left_occ) * kMinBlockSizePix);
  const bool have_top_left = have_top && (-1 < y_left_occ);
  const bool have_top_right = have_top && (-1 < y_right_occ);

  // collect the context samples:
  size_t m = 0;
  // this is the value we reuse in case of non-availability
  int16_t last = fill_in ? 0 : kMissing;
  // left (bh samples)
  const int16_t* const context_src_data =
      context_src.IsEmpty() ? nullptr : &context_src.AtUnclamped(0, 0);
  const int context_src_step = context_src.Step();
  if (left_context_size > 0) {
    if (fill_in) last = CONTEXT_SRC_AT_UNCLAMPED(-1, left_context_size - 1);
    int j = h - 1;
    for (; j >= left_context_size; --j) context[m++] = last;
    for (; j >= 0; --j) context[m++] = CONTEXT_SRC_AT_UNCLAMPED(-1, j);
    if (fill_in) last = context[m - 1];
  } else {
    if (fill_in) {
      last =
          have_top ? CONTEXT_SRC_AT_UNCLAMPED(have_top_left ? -1 : 0, -1) : 0;
    }
    for (uint32_t j = 0; j < h; ++j) context[m++] = last;
  }
  // top-left
  if (have_top_left) {
    context[m++] = CONTEXT_SRC_AT_UNCLAMPED(-1, -1);
    if (fill_in) last = context[m - 1];
  } else {
    context[m++] = last;
  }
  // top (bw samples)
  if (have_top) {
    for (uint32_t i = 0; i < w; ++i) {
      context[m++] = CONTEXT_SRC_AT_UNCLAMPED(i, -1);
    }
    if (fill_in) last = context[m - 1];
  } else {
    if (fill_in) {
      last = (left_context_size > 0) ? CONTEXT_SRC_AT_UNCLAMPED(-1, 0) : 0;
    }
    for (uint32_t i = 0; i < w; ++i) context[m++] = last;
  }

  if (extend_right) {
    const uint32_t context_size = ContextWithTRSize(w, h);
    right_extent *= kMinBlockSizePix;
    int j = 0;
    for (; j < right_extent && m < context_size; ++j) {
      context[m++] = CONTEXT_SRC_AT_UNCLAMPED(w + j, -1);
    }
    if (fill_in) last = context[m - 1];
    for (; m < context_size; ++j) context[m++] = last;
  } else {
    // top-right
    if (have_top_right) {
      context[m++] = CONTEXT_SRC_AT_UNCLAMPED(w, -1);
    } else {
      context[m++] = last;
    }
    // right (bh samples)
    uint32_t j = 0;
    const uint32_t right_context_size =
        std::min(h, std::max(0, (int)y_right_occ) * kMinBlockSizePix);
    for (; j < right_context_size; ++j) {
      context[m++] = CONTEXT_SRC_AT_UNCLAMPED(w, j);
    }
    if (fill_in) last = context[m - 1];
    for (; j < h; ++j) context[m++] = last;
    assert(m == ContextSize(w, h));
  }
  const uint32_t kSizeThreshold = 48;  // context size above which to filter
  const int32_t kDiffThreshold = 96;   // slope threshold below which to filter
  if (m > kSizeThreshold) {   // block is large enough?
    int32_t prev = context[0], cur = context[0];
    uint32_t prev_diff = 0;
    for (uint32_t i = 0; i < m - 1; ++i) {
      const int32_t next = context[i + 1];
      // Smaller or un-aligned blocks in the context can lead to edge smearing
      // during prediction (not fixed by the deblocking filter because it
      // happens inside the blocks). Hence low-pass smooth the context samples
      // that don't seem to belong to a real edge.
      const uint32_t next_diff = std::abs(cur - next);
      if (prev_diff + next_diff < kDiffThreshold) {
        context[i] = (6 * cur + prev + next + 4) >> 3;
      }
      prev_diff = next_diff;
      prev = cur;
      cur = next;
    }
    const int32_t next = context[m - 1];
    if (prev_diff < kDiffThreshold) {
      context[m - 1] = (6 * cur + prev + next + 4) >> 3;
    }
  }
}

#undef CONTEXT_SRC_AT_UNCLAMPED

//------------------------------------------------------------------------------

static int CodingParamsIndex(Channel channel) {
  if (channel == kYChannel) return 0;
  if (channel == kUChannel || channel == kVChannel) return 1;
  if (channel == kAChannel) return 2;
  assert(false);
  return -1;
}

const CodedBlock::CodingParams& CodedBlock::GetCodingParams(
    Channel channel) const {
  return params_[CodingParamsIndex(channel)];
}

CodedBlock::CodingParams* CodedBlock::GetCodingParams(Channel channel) {
  return &params_[CodingParamsIndex(channel)];
}

TransformPair CodedBlock::GetImplicitTf(Channel channel) const {
  // TODO(maryla) check if it is still needed.
  if (channel == kAChannel) return kDctDct;
  TransformPair tf;
  if (GetCodingParams(channel).pred->OverridesTransform(channel, &tf)) {
    return tf;
  }
  return kUnknownTf;
}

//------------------------------------------------------------------------------

uint32_t CodedBlock::GetDisto(Channel channel,
                              const Rectangle& tile_rect) const {
  const Plane16& p_in = in_.GetChannel(channel);
  const Plane16& p_out = out_.GetChannel(channel);
  return WP2SumSquaredErrorBlock(p_in.Row(0), p_in.Step(), p_out.Row(0),
                                 p_out.Step(), visible_w_pix(tile_rect),
                                 visible_h_pix(tile_rect));
}

uint32_t CodedBlock::GetNumTransforms(Channel channel) const {
  return GetNumTransformsInBlock(dim(), GetCodingParams(channel).split_tf);
}

uint32_t CodedBlock::NumCoeffsPerTransform(Channel channel) const {
  return kNumBlocks[GetSplitSize(blk_.dim(),
                                 GetCodingParams(channel).split_tf)] *
         (((channel == kUChannel || channel == kVChannel) && is420_) ? 4 : 16);
}

bool CodedBlock::IsFirstCoeffDC(Channel channel) const {
  return (GetCodingParams(channel).tf_x() != kIdentity &&
          GetCodingParams(channel).tf_y() != kIdentity);
}

bool CodedBlock::HasCoeffs(Channel channel) const {
  for (uint32_t tf_i = 0; tf_i < GetNumTransforms(channel); ++tf_i) {
    if (num_coeffs_[channel][tf_i] != 0) return true;
  }
  return false;
}

uint32_t CodedBlock::visible_w_pix(const Rectangle& tile_rect) const {
  assert(x_pix() < tile_rect.width);
  return std::min(w_pix(), tile_rect.width - x_pix());
}
uint32_t CodedBlock::visible_h_pix(const Rectangle& tile_rect) const {
  assert(y_pix() < tile_rect.height);
  return std::min(h_pix(), tile_rect.height - y_pix());
}

//------------------------------------------------------------------------------

WP2Status CodedBlock::ResidualRate(const BlockContext& context, Channel channel,
                                   uint32_t num_channels,
                                   SymbolCounter* const counter,
                                   SymbolCounter* const aom_counter,
                                   float* rate) const {
  const bool is_uv = (channel == kUChannel || channel == kVChannel);
  // TODO(vrabaud) Using AOM context for U/V is worse for now...
  if (context.use_aom() && !is_uv) {
    WP2_CHECK_STATUS(ResidualWriter::GetRateAOM(*this, channel, *context.aom(),
                                                aom_counter, rate));
  } else {
    // Whether to use the more accurate but a lot more expensive rate.
    // Even this one is not perfect because the cost depends on global
    // distributions which are not known at this time.
    static constexpr bool kUseExpensiveRate = false;
    *rate = 0.f;
    for (uint32_t tf_i = 0; tf_i < GetNumTransforms(channel); ++tf_i) {
      float tf_rate;
      if (kUseExpensiveRate) {
        WP2_CHECK_STATUS(ResidualWriter::GetRate(
            channel, num_channels, tdim(channel), coeffs_[channel][tf_i],
            num_coeffs_[channel][tf_i], IsFirstCoeffDC(channel), counter,
            &tf_rate));
      } else {
        WP2_CHECK_STATUS(ResidualWriter::GetPseudoRate(
            channel, num_channels, tdim(channel), coeffs_[channel][tf_i],
            num_coeffs_[channel][tf_i], IsFirstCoeffDC(channel), counter,
            &tf_rate));
      }
      *rate += tf_rate;
    }
  }
  return WP2_STATUS_OK;
}

// Quantizes 'coeffs', then dequantizes, detransforms them and computes
// distortion compared to 'res'.
static uint32_t RoundTripDisto(CodedBlock* const cb, const QuantMtx& quant,
                               Channel channel, TrfSize tdim,
                               const int16_t* const res,
                               const int32_t* const coeffs, bool reduced,
                               int16_t* const quantized_coeffs,
                               uint32_t* const num_coeffs) {
  int32_t dequantized_coeffs[kMaxBlockSizePix2];
  quant.Quantize(coeffs, tdim, cb->IsFirstCoeffDC(channel), quantized_coeffs,
                 num_coeffs);
  quant.Dequantize(quantized_coeffs, *num_coeffs, tdim, dequantized_coeffs);
  const CodedBlock::CodingParams& params = *cb->GetCodingParams(channel);
  WP2InvTransform2D(dequantized_coeffs, params.tf_x(), params.tf_y(),
                    cb->w_pix(), cb->h_pix(), dequantized_coeffs, reduced);
  // Copy of dequantized_coeffs as int16 to compute error.
  int16_t dequantized_coeffs16[kMaxBlockSizePix2];
  std::copy(dequantized_coeffs, dequantized_coeffs + cb->w_pix() * cb->h_pix(),
            dequantized_coeffs16);
  return WP2SumSquaredErrorBlock(res, cb->w_pix(), dequantized_coeffs16,
                                 cb->w_pix(), cb->w_pix(), cb->h_pix());
}

// Computes the rate-distortion score for the U and V channels.
// 'res_u' and 'res_v' should be in frequency space. 'orig_u' and 'orig_v' are
// the original spatial residuals to compare against for distortion.
static WP2Status UVScore(const EncoderConfig& config, uint32_t tile_pos_x,
                         uint32_t tile_pos_y, bool has_alpha,
                         const QuantMtx& quant_u, const QuantMtx& quant_v,
                         const int16_t orig_u[kMaxBlockSizePix2],
                         const int16_t orig_v[kMaxBlockSizePix2],
                         const int32_t res_u[kMaxBlockSizePix2],
                         const int32_t res_v[kMaxBlockSizePix2], TrfSize dim,
                         bool reduced, SymbolCounter* const counter,
                         CodedBlock* const cb, float* score) {
  const float lambda_u = quant_u.lambda2 * 2.f;  // Another magic constant.
  const float lambda_v = quant_v.lambda2 * 2.f;

  uint32_t tmp_num_coeffs_u, tmp_num_coeffs_v;
  int16_t tmp_coeffs_u[kMaxBlockSizePix2];
  int16_t tmp_coeffs_v[kMaxBlockSizePix2];
  const uint32_t disto =
      RoundTripDisto(cb, quant_u, kUChannel, dim, orig_u, res_u, reduced,
                     tmp_coeffs_u, &tmp_num_coeffs_u) +
      RoundTripDisto(cb, quant_v, kVChannel, dim, orig_v, res_v, reduced,
                     tmp_coeffs_v, &tmp_num_coeffs_v);
  // early out if there's ~only DC
  if (!reduced && tmp_num_coeffs_u <= 1 && tmp_num_coeffs_v <= 1) {
    cb->Store420Scores(config, tile_pos_x, tile_pos_y, lambda_u, lambda_v,
                       reduced, disto, 0, 0);
    *score = 0;
    return WP2_STATUS_OK;
  }
  const uint32_t num_channels = (has_alpha ? 4 : 3);
  float rate_u, rate_v;
  WP2_CHECK_STATUS(ResidualWriter::GetPseudoRate(
      kUChannel, num_channels, dim, tmp_coeffs_u, tmp_num_coeffs_u,
      cb->IsFirstCoeffDC(kUChannel), counter, &rate_u));
  WP2_CHECK_STATUS(ResidualWriter::GetPseudoRate(
      kVChannel, num_channels, dim, tmp_coeffs_v, tmp_num_coeffs_v,
      cb->IsFirstCoeffDC(kVChannel), counter, &rate_v));
  *score = disto + lambda_u * rate_u + lambda_v * rate_v;
  cb->Store420Scores(config, tile_pos_x, tile_pos_y, lambda_u, lambda_v,
                     reduced, disto, rate_u, rate_v);
  return WP2_STATUS_OK;
}

WP2Status CodedBlock::DecideChromaSubsampling(
    const EncoderConfig& config, uint32_t tile_pos_x, uint32_t tile_pos_y,
    bool has_alpha, const QuantMtx& quant_u, const QuantMtx& quant_v,
    Counters* const counters) {
  const BlockSize dim = blk_.dim();

  const TrfSize full_dim = kFullDim[dim];
  const TrfSize half_dim = kHalfDim[dim];
  const uint32_t w_pix = blk_.w_pix();
  const uint32_t h_pix = blk_.h_pix();
  SubtractBlockFunc sub_block = WP2::SubtractBlock[TrfLog2[w_pix] - 2];

  const CodingParams* const u_params = GetCodingParams(kUChannel);
  const CodingParams* const v_params = GetCodingParams(kVChannel);

  int32_t res_u[kMaxBlockSizePix2], res_v[kMaxBlockSizePix2];
  {
    int16_t prediction[kMaxBlockSizePix2];
    u_params->pred->Predict(*this, kUChannel, prediction, w_pix);
    sub_block(in_.U.Row(0), in_.U.Step(), prediction, /*preds_step=*/w_pix,
              res_u, /*dst_step=*/w_pix, h_pix);
    v_params->pred->Predict(*this, kVChannel, prediction, w_pix);
    sub_block(in_.V.Row(0), in_.V.Step(), prediction, /*preds_step=*/w_pix,
              res_v, /*dst_step=*/w_pix, h_pix);
  }
  int16_t orig_u[kMaxBlockSizePix2], orig_v[kMaxBlockSizePix2];
  std::copy(res_u, res_u + w_pix * h_pix, orig_u);
  std::copy(res_v, res_v + w_pix * h_pix, orig_v);

  // transform at full resolution
  WP2Transform2D(res_u, u_params->tf_x(), u_params->tf_y(), w_pix, h_pix,
                 res_u);
  WP2Transform2D(res_v, v_params->tf_x(), v_params->tf_y(), w_pix, h_pix,
                 res_v);

  float score444;
  WP2_CHECK_STATUS(UVScore(config, tile_pos_x, tile_pos_y, has_alpha, quant_u,
                           quant_v, orig_u, orig_v, res_u, res_v, full_dim,
                           /*reduced=*/false, counters->residuals(), this,
                           &score444));
  if (score444 == 0) {
    is420_ = false;
    WP2_CHECK_STATUS(Store420Decision(config, tile_pos_x, tile_pos_y,
                                      Debug420Decision::k444EarlyExit));
    return WP2_STATUS_OK;
  }

  // inspect half-resolution now
  WP2ReduceCoeffs(res_u, w_pix, h_pix, res_u);
  WP2ReduceCoeffs(res_v, w_pix, h_pix, res_v);
  float score420;
  WP2_CHECK_STATUS(UVScore(config, tile_pos_x, tile_pos_y, has_alpha, quant_u,
                           quant_v, orig_u, orig_v, res_u, res_v, half_dim,
                           /*reduced=*/true, counters->residuals(), this,
                           &score420));

  is420_ = (score420 < score444);
  WP2_CHECK_STATUS(Store420Decision(
      config, tile_pos_x, tile_pos_y,
      is420_ ? Debug420Decision::k420 : Debug420Decision::k444));
  return WP2_STATUS_OK;
}

//------------------------------------------------------------------------------

void CodedBlock::PredictBlock(Channel channel, int16_t* out,
                              uint32_t step) const {
  const CodingParams& params = GetCodingParams(channel);
  params.pred->Predict(*this, channel, out, step);
}

void CodedBlock::GetResiduals(Channel channel, const Plane16& prediction,
                              BlockCoeffs32* const res) const {
  const Plane16& in = in_.GetChannel(channel);
  const BlockSize split_size =
      GetSplitSize(dim(), GetCodingParams(channel).split_tf);
  const uint32_t split_w = BlockWidthPix(split_size);
  const uint32_t split_h = BlockHeightPix(split_size);
  SubtractBlockFunc sub_block = WP2::SubtractBlock[TrfLog2[split_w] - 2];
  uint32_t tf_i = 0;
  for (uint32_t split_y = 0; split_y < h_pix(); split_y += split_h) {
    for (uint32_t split_x = 0; split_x < w_pix(); split_x += split_w) {
      sub_block(&in.At(split_x, split_y), in.Step(),
                &prediction.At(split_x, split_y), prediction.Step(),
                (*res)[tf_i], split_w, split_h);

#if defined(WP2_BITTRACE)
      if (original_res_ != nullptr) {
        memcpy((*original_res_)[channel][tf_i], &(*res)[tf_i][0],
               split_w * split_h * sizeof((*res)[tf_i][0]));
      }
#endif
      ++tf_i;
    }
  }
}

void CodedBlock::TransformAndReconstruct(const QuantMtx& quant, Channel channel,
                                         bool reduced_transform,
                                         const Plane16& prediction,
                                         BlockCoeffs32* const res) {
  const CodingParams* const params = GetCodingParams(channel);

  const BlockSize split_size =
      GetSplitSize(dim(), GetCodingParams(channel)->split_tf);
  const uint32_t split_w = BlockWidthPix(split_size);
  const uint32_t split_h = BlockHeightPix(split_size);
  for (uint32_t tf_i = 0; tf_i < GetNumTransforms(channel); ++tf_i) {
    int32_t tf_res[kMaxBlockSizePix2];
    WP2Transform2D((*res)[tf_i], params->tf_x(), params->tf_y(), split_w,
                   split_h, tf_res, reduced_transform);

    if (channel == kUChannel || channel == kVChannel) {
      assert(params->tf == kDctDct && !params->split_tf);
      // Note DC is scaled by NumPix(dim) by WP2Transform2D. We don't
      // normalize it when propagating the error, even though the same value
      // will have a bigger impact on a smaller block than on a large one.
      // This is intentional: since we have fewer pixels in a smaller block,
      // we shift them more. But it would also be possible to normalize by a
      // factor of NumPix(dim) (multiply dc_error_ by this and devide
      // dc_error_next by the same factor). This would mean we apply the same
      // color shift regardless of the block size.
      tf_res[0] += dc_error_[channel];
      dc_error_next_[channel] = quant.DCError(tf_res[0], tdim(channel));
    }
    quant.Quantize(tf_res, tdim(channel), IsFirstCoeffDC(channel),
                   coeffs_[channel][tf_i], &num_coeffs_[channel][tf_i]);
    if (quant.qstats != nullptr && channel == kYChannel) {
      quant.qstats->Record(tdim(channel), tf_res, coeffs_[channel][tf_i]);
    }
    if (channel == kYChannel && mtx_set_ != nullptr) {
      mtx_set_->DecideUseRndMtx(this);  // TODO(skal): later: U/V channel?
    }
  }
  Dequantize(quant, channel, res);
  Reconstruct(channel, reduced_transform, res, prediction);
}

void CodedBlock::Reconstruct(Channel channel, bool reduced_transform,
                             BlockCoeffs32* const res,
                             const Plane16& prediction) const {
  const CodingParams& params = GetCodingParams(channel);
  const BlockSize split_size = GetSplitSize(dim(), params.split_tf);
  const uint32_t split_w = BlockWidthPix(split_size);
  const uint32_t split_h = BlockHeightPix(split_size);
  for (uint32_t tf_i = 0; tf_i < GetNumTransforms(channel); ++tf_i) {
    WP2InvTransform2D((*res)[tf_i], params.tf_x(), params.tf_y(), split_w,
                      split_h, (*res)[tf_i], reduced_transform);
  }

  Plane16* const out = &out_.GetChannel(channel);
  Plane16 pred_view;
  if (prediction.IsEmpty()) {
    WP2_ASSERT_STATUS(pred_view.SetView(*out));
    PredictBlock(channel, pred_view.Row(0), pred_view.Step());
  } else {
    WP2_ASSERT_STATUS(pred_view.SetView(prediction));
  }

  AddBlockFunc add_block = WP2::AddBlock[TrfLog2[split_w] - 2];
  const int32_t min = (channel == kAChannel) ? 0 : yuv_min_;
  const int32_t max = (channel == kAChannel) ? kAlphaMax : yuv_max_;
  for (uint32_t tf_i = 0, split_y = 0; split_y < h_pix(); split_y += split_h) {
    for (uint32_t split_x = 0; split_x < w_pix(); split_x += split_w, ++tf_i) {
      add_block(&pred_view.At(split_x, split_y), pred_view.Step(), (*res)[tf_i],
                split_w, min, max, &out->At(split_x, split_y), out->Step(),
                split_h);
    }
  }
}

void CodedBlock::Quantize(const QuantMtx& quant, Channel channel,
                          bool reduced_transform) {
  BlockCoeffs32 res;
  Plane16* const out = &out_.GetChannel(channel);
  PredictBlock(channel, out->Row(0), out->Step());
  GetResiduals(channel, /*prediction=*/*out, &res);
  TransformAndReconstruct(quant, channel, reduced_transform,
                          /*prediction=*/*out, &res);
}

void CodedBlock::Dequantize(const QuantMtx& quant, Channel channel,
                            BlockCoeffs32* const res) {
  const BlockSize size =
      GetSplitSize(dim(), GetCodingParams(channel)->split_tf);
  for (uint32_t tf_i = 0; tf_i < GetNumTransforms(channel); ++tf_i) {
    int16_t* src = coeffs_[channel][tf_i];
    int16_t tmp[kMaxBlockSizePix2];
    if (use_mtx_ && channel == kYChannel) {
      assert(mtx_set_->IsOk(mtx_[channel], size));
      std::copy(src, src + NumCoeffsPerTransform(channel), tmp);
      mtx_set_->AddMatrix(mtx_[channel], size, tmp);
      src = tmp;
    }
    quant.Dequantize(src, num_coeffs_[channel][tf_i], tdim(channel),
                     (*res)[tf_i]);
  }
}

//------------------------------------------------------------------------------

WP2Status CodedBlock::FindBestTf(const EncoderConfig& config,
                                 const Rectangle& tile_rect,
                                 const QuantMtx& quant,
                                 const BlockContext& context, Channel channel,
                                 uint32_t num_channels, bool reduced,
                                 const Vector<TransformPair>& transforms,
                                 Counters* const counters) {
  assert(!transforms.empty());
  CodingParams* const params = GetCodingParams(channel);

  params->tf = GetImplicitTf(channel);
  if (params->tf != kUnknownTf) return WP2_STATUS_OK;
  params->tf = transforms.front();
  if (transforms.size() == 1) return WP2_STATUS_OK;

  TransformPair best_tf = transforms.front();
  float best_score = 0.f;
  const float lambda =
      lambda_mult_ * ((channel == kYChannel) ? quant.lambda1 : quant.lambda2);
  for (const TransformPair tf : transforms) {
    params->tf = tf;
    Quantize(quant, channel, reduced);

    const uint32_t distortion = GetDisto(channel, tile_rect);
    float tf_rate, res_rate;
    WP2_CHECK_STATUS(TransformRate(channel, transforms.size(),
                                   counters->transform(), &tf_rate));
    WP2_CHECK_STATUS(ResidualRate(context, channel, num_channels,
                                  counters->residuals(),
                                  counters->residuals_aom(), &res_rate));
    const float current_score = distortion + lambda * (res_rate + tf_rate);
    const bool is_best =
        (tf == transforms.front() || current_score < best_score);
    if (is_best) {
      best_tf = tf;
      best_score = current_score;
    }
    StorePredictionScore(config, tile_rect, channel, *params->pred, tf,
                         distortion, lambda, res_rate, /*pred_rate=*/0, tf_rate,
                         current_score, is_best);
  }

  params->tf = best_tf;
  return WP2_STATUS_OK;
}

void CodedBlock::FindBestPred(const Predictors& preds, Channel channel,
                              bool reduced) {
  assert(preds.size() > 1);  // Otherwise this function should be skipped.
  // kVChannel is tested along kUChannel.
  assert(channel == kYChannel || channel == kUChannel || channel == kAChannel);
  const Plane16& in = in_.GetChannel(channel);
  Plane16& out = out_.GetChannel(channel);
  CodingParams* const params = GetCodingParams(channel);

  uint64_t best_score = 0;
  const Predictor* best_pred = nullptr;
  for (const Predictor* const p : preds) {
    params->pred = p;
    p->Predict(*this, channel, out.Row(0), out.Step());
    uint64_t score =
        WP2SumSquaredErrorBlock(in.Row(0), in.Step(), out.Row(0), out.Step(),
                                blk_.w_pix(), blk_.h_pix());
    if (channel == kUChannel) {
      p->Predict(*this, kVChannel, out_.V.Row(0), out_.V.Step());
      score +=
          WP2SumSquaredErrorBlock(in_.V.Row(0), in_.V.Step(), out_.V.Row(0),
                                  out_.V.Step(), blk_.w_pix(), blk_.h_pix());
    }
    if (best_pred == nullptr || score < best_score) {
      best_score = score;
      best_pred = p;
    }
  }
  params->pred = best_pred;
  pred_scores_[channel] = best_score;
  if (channel == kUChannel) pred_scores_[kVChannel] = pred_scores_[kUChannel];
}

WP2Status CodedBlock::PredictorRate(
    const YModePredictor* const ymode_predictor, Channel channel,
    SymbolCounter* const counter, float* const rate) const {
  counter->Clear();
  ANSEncCounter enc;
  if (channel == kYChannel || channel == kAChannel) {
    SyntaxWriter::WriteYAPredictors(*this, channel, ymode_predictor,
                                    /*update_ymodes=*/false, counter, &enc);
  } else {
    SyntaxWriter::WriteUVPredictors(*this, channel, counter, &enc);
  }
  *rate = enc.GetCost();
  return WP2_STATUS_OK;
}

WP2Status CodedBlock::TransformRate(Channel channel, uint32_t num_transforms,
                                    SymbolCounter* const counter,
                                    float* const rate) const {
  ANSEncCounter enc;
  counter->Clear();
  SyntaxWriter::WriteSplitTransform(*this, channel, counter, &enc);
  SyntaxWriter::WriteHasCoeffs(*this, channel, counter, &enc);
  SyntaxWriter::WriteTransform(*this, channel, counter, &enc);
  *rate = enc.GetCost();
  return WP2_STATUS_OK;
}

// Get residual and predictor rates.
WP2Status CodedBlock::GetRates(const Rectangle& tile_rect,
                               const BlockContext& context,
                               const QuantMtx& quant, Channel channel,
                               uint32_t num_channels, float tf_rate,
                               Counters* const counters, uint32_t* const dist,
                               float* const res_rate, float* const pred_rate,
                               float* const score) const {
  *dist = GetDisto(channel, tile_rect);

  const bool is_uv = (channel == kUChannel || channel == kVChannel);
  WP2_CHECK_STATUS(ResidualRate(context, channel, num_channels,
                                counters->residuals(),
                                counters->residuals_aom(), res_rate));
  WP2_CHECK_STATUS(PredictorRate(is_uv ? nullptr : &context.ymodes(), channel,
                                 counters->predictor(), pred_rate));

  const float rate = *res_rate + *pred_rate + tf_rate;
  *score =
      *dist + lambda_mult_ * (is_uv ? quant.lambda2 : quant.lambda1) * rate;

  return WP2_STATUS_OK;
}

WP2Status CodedBlock::FindBestPredTf(
    const EncoderConfig& config, const Rectangle& tile_rect,
    const Predictors& preds, const Segment& segment,
    const BlockContext& context, Channel channel, uint32_t num_channels,
    bool reduced, const Vector<TransformPair>& transforms,
    Counters* const counters) {
  assert(channel == kYChannel || channel == kUChannel || channel == kAChannel);
  if (channel == kUChannel) {
    is420_ = reduced;  // Final is420 will be decided later but is needed now.
    assert(transforms.empty());  // Empty because UV tfs are defined by preds.
  } else {
    assert(!transforms.empty());
  }

  float best_score = std::numeric_limits<float>::max();
  TransformPair best_tf = kUnknownTf;
  const Predictor* best_pred = nullptr;
  WP2_CHECK_STATUS(preds.FindBest(config, tile_rect, num_channels, context,
                                  segment, reduced, transforms, counters, this,
                                  &best_score, &best_pred, &best_tf));
  SetPredictor(channel, best_pred);
  pred_scores_[channel] = best_score;
  if (channel == kUChannel) {
    pred_scores_[kVChannel] = pred_scores_[kUChannel];
    assert(GetCodingParams(channel)->tf == best_tf);
  } else {
    GetCodingParams(channel)->tf = best_tf;
  }
  return WP2_STATUS_OK;
}

WP2Status CodedBlock::FindBestSplitTf(const EncoderConfig& config,
                                      const Rectangle& tile_rect,
                                      const QuantMtx& quant,
                                      const BlockContext& context,
                                      Channel channel, uint32_t num_channels,
                                      bool reduced, Counters* const counters) {
  CodingParams* const params = GetCodingParams(channel);
  if (channel != kYChannel || GetSplitSize(dim(), /*split=*/true) == dim()) {
    assert(!params->split_tf);
    return WP2_STATUS_OK;
  }

  bool best_split_tf = false;  // avoids -Werror=maybe-uninitialized
  float best_score = 0.f;
  const float lambda =
      lambda_mult_ * ((channel == kYChannel) ? quant.lambda1 : quant.lambda2);
  for (bool split_tf : {false, true}) {
    params->split_tf = split_tf;
    Quantize(quant, channel, reduced);

    if (params->split_tf && !HasCoeffs(channel)) {
      continue;  // Split-tf with no coeff is forbidden.
    }

    const uint32_t distortion = GetDisto(channel, tile_rect);
    float tf_rate;
    float res_rate;
    WP2_CHECK_STATUS(TransformRate(channel, /*num_transforms=*/1,
                                   counters->transform(), &tf_rate));
    WP2_CHECK_STATUS(ResidualRate(context, channel, num_channels,
                                  counters->residuals(),
                                  counters->residuals_aom(), &res_rate));
    const float current_score = distortion + lambda * (tf_rate + res_rate);
    const bool is_best = (!split_tf || current_score < best_score);
    if (is_best) {
      best_split_tf = split_tf;
      best_score = current_score;
    }
    StorePredictionScore(config, tile_rect, channel, *params->pred, params->tf,
                         distortion, lambda, res_rate, /*pred_rate=*/0, tf_rate,
                         current_score, is_best);
  }

  params->split_tf = best_split_tf;
  return WP2_STATUS_OK;
}

WP2Status CodedBlock::OptimizeModesLuma(
    const EncoderConfig& config, const Rectangle& tile_rect, bool has_alpha,
    const Predictors& preds, const Segment& segment,
    const BlockContext& context, const Vector<TransformPair>& transforms,
    const Vector<TransformPair>& transforms_subset, Counters* const counters) {
  assert(preds.size() > 1);
  for (const Predictor* const pred : preds) {
    pred->ComputeParams(this, kYChannel);
  }

  const uint32_t num_channels = (has_alpha ? 4 : 3);
  GetCodingParams(kYChannel)->split_tf = false;
  bool found_best_tf = false;
  if (y_context_is_constant_) {
    // If the context is constant, all luma predictors will predict the same
    // values (plus or minus 1, because of rounding errors), so we just force
    // the predictor to 0 and we don't write it to the bitstream. We can't do
    // the same for UV and alpha because the chroma from luma predictor doesn't
    // have this property, but we could send a single bit saying whether the
    // predictor is cfl or "other".
    SetLumaUniformPredictor(preds);
    pred_scores_[kYChannel] = 1.f;
  } else if (SetForcedPredictor(preds, config, tile_rect.x, tile_rect.y,
                                kYChannel, this)) {
    // Predictor is set by EncoderConfig::info (debug).
  } else if (config.speed == 0) {
    // Select the first predictor to save time.
    SetPredictor(kYChannel, preds[0]);
    pred_scores_[kYChannel] = 1.f;
  } else if (config.speed == 1) {
    // Quickly select the best predictor and then choose the best transform.
    FindBestPred(preds, kYChannel, /*reduced=*/false);
  } else {
    // Select the best predictor+transform pair and then choose the best
    // transform if there are more to try.
    WP2_CHECK_STATUS(FindBestPredTf(config, tile_rect, preds, segment, context,
                                    kYChannel, num_channels, /*reduced=*/false,
                                    transforms_subset, counters));
    if (transforms.size() == transforms_subset.size()) found_best_tf = true;
  }
  if (!found_best_tf) {
    // TODO(yguyon): Do not recompute scores of 'transforms_subset'
    WP2_CHECK_STATUS(FindBestTf(config, tile_rect, segment.quant_y_, context,
                                kYChannel, num_channels, /*reduced=*/false,
                                transforms, counters));
  }
  WP2_CHECK_STATUS(FindBestSplitTf(config, tile_rect, segment.quant_y_, context,
                                   kYChannel, num_channels, /*reduced=*/false,
                                   counters));

  // TODO(yguyon): Remember the params during the last Quantize() call and omit
  //               this one if they are equivalent.
  Quantize(segment.quant_y_, kYChannel, /*reduced_transform=*/false);
  return WP2_STATUS_OK;
}

WP2Status CodedBlock::OptimizeModesChroma(
    const EncoderConfig& config, const Rectangle& tile_rect, bool has_alpha,
    const FrontMgrBase& mgr, const UVPredictors& preds, const Segment& segment,
    const BlockContext& context, ChromaSubsampling chroma_subsampling,
    DCDiffusionMap* const dc_error_u, DCDiffusionMap* const dc_error_v,
    Counters* const counters) {
  const uint32_t diffusion =
      DCDiffusionMap::GetDiffusion(config.error_diffusion);
  for (const Predictor* const pred : preds) {
    pred->ComputeParams(this, kUChannel);
    pred->ComputeParams(this, kVChannel);
  }
  // Decide if we reduce chroma for finding the best predictors. The
  // logic is based on what is done a few lines after.
  const bool reduced =
      ((chroma_subsampling == ChromaSubsampling::kSingleBlock &&
        config.uv_mode == EncoderConfig::UVMode420) ||
       chroma_subsampling == ChromaSubsampling::k420);

  if (SetForcedPredictor(preds, config, tile_rect.x, tile_rect.y, kUChannel,
                         this)) {
    // Predictor is set by EncoderConfig::info (debug).
  } else if (config.speed == 0 || preds.size() == 1) {
    // Select the first predictor to save time.
    SetUVPredictor(preds[0]);
    pred_scores_[kUChannel] = pred_scores_[kVChannel] = 1.f;
  } else if (config.speed == 1) {
    // Quickly select the best predictor.
    FindBestPred(preds, kUChannel, /*reduced=*/false);
  } else {
    const uint32_t num_channels = (has_alpha ? 4 : 3);
    // Empty because UV transforms are defined by predictors.
    const Vector<TransformPair> transforms;
    // Spend time to select the best predictor.
    WP2_CHECK_STATUS(FindBestPredTf(config, tile_rect, preds, segment, context,
                                    kUChannel, num_channels, reduced,
                                    transforms, counters));
  }

  if (diffusion > 0) {  // Add error diffusion
    dc_error_[kUChannel] = dc_error_u->Get(blk(), diffusion);
    dc_error_[kVChannel] = dc_error_v->Get(blk(), diffusion);
  }

  // TODO(skal): decide based on coded luma too?
  if (reduced) {
    is420_ = true;
  } else if ((chroma_subsampling == ChromaSubsampling::kSingleBlock &&
              config.uv_mode == EncoderConfig::UVMode444) ||
             chroma_subsampling == ChromaSubsampling::k444) {
    is420_ = false;
  } else {
    assert(chroma_subsampling == ChromaSubsampling::kSingleBlock ||
           chroma_subsampling == ChromaSubsampling::kAdaptive);
    WP2_CHECK_STATUS(DecideChromaSubsampling(config, tile_rect.x, tile_rect.y,
                                             has_alpha, segment.quant_u_,
                                             segment.quant_v_, counters));
  }
  Quantize(segment.quant_u_, kUChannel, is420_);
  Quantize(segment.quant_v_, kVChannel, is420_);
  if (diffusion > 0) {
    dc_error_u->Store(mgr, blk(), dc_error_next_[kUChannel]);
    dc_error_v->Store(mgr, blk(), dc_error_next_[kVChannel]);
  }
  return WP2_STATUS_OK;
}

//------------------------------------------------------------------------------

bool CodedBlock::ContextIsConstant(Channel channel) const {
  assert(context_cache_ != nullptr);
  const int16_t* context = context_cache_->GetLarge(
      *this, channel, /*fill_in=*/true, /*extend_right=*/false);
  const uint32_t context_size = ContextSize(blk_.w_pix(), blk_.h_pix());
  for (uint32_t i = 1; i < context_size; ++i) {
    if (context[i] != context[0]) {
      return false;
    }
  }
  return true;
}

bool CodedBlock::HasLossyAlpha() const {
  return (alpha_mode_ == kBlockAlphaLossy);
}

//------------------------------------------------------------------------------

void CodedBlock::SetUVPredictor(const Predictor* const pred) {
  CodedBlock::CodingParams* const params = GetCodingParams(kUChannel);
  params->pred = pred;
  if (pred->OverridesTransform(kUChannel, &params->tf)) {
    assert(params->tf != kUnknownTf);
  } else {
    assert(false);  // The transform should always be overridden for U and V.
  }
}

void CodedBlock::SetPredictor(Channel channel, const Predictor* const pred) {
  if (channel == kUChannel || channel == kVChannel) {
    SetUVPredictor(pred);
  } else {
    CodingParams* const params = GetCodingParams(channel);
    params->pred = pred;
  }
}

void CodedBlock::SetLumaUniformPredictor(const Predictors& preds) {
  SetPredictor(kYChannel, preds[0]);
}

//------------------------------------------------------------------------------

}  // namespace WP2
