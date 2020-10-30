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
// WP2 lossy encoding of the alpha plane.
//
// Author: Maryla (maryla@google.com)

#include "src/enc/lossless/losslessi_enc.h"
#include "src/enc/wp2_enc_i.h"
#include "src/wp2/format_constants.h"

namespace WP2 {

WP2Status AlphaWriter::Init(const EncoderConfig& config,
                            const GlobalParams& gparams,
                            const BlockContext& context, const YUVPlane& yuv,
                            const Rectangle& tile_rect) {
  assert(gparams.has_alpha_);

  tile_rect_ = tile_rect;
  context_ = &context;
  config_ = config;
  gparams_ = &gparams;
  lossless_enc_.AddDebugPrefix("Alpha");

  alpha_mode_ =
      (gparams.maybe_use_lossy_alpha_) ? kAlphaModeLossy : kAlphaModeLossless;

  if (alpha_mode_ == kAlphaModeLossy) {
    WP2_CHECK_STATUS(mode_predictor_.Init(tile_rect_.width, tile_rect_.height));
  } else {
    // Convert alpha plane to a grayscale image (not padded).
    WP2_CHECK_STATUS(alpha_.Resize(tile_rect_.width, tile_rect_.height));
    for (uint32_t y = 0; y < tile_rect_.height; ++y) {
      const int16_t* const src_row = yuv.A.Row(y);
      uint8_t* const dst = (uint8_t*)alpha_.GetRow(y);
      for (uint32_t x = 0; x < tile_rect_.width; ++x) {
        const uint8_t alpha = Clamp<int16_t>(src_row[x], 0, kAlphaMax);
        dst[4 * x + 0] = kAlphaMax;
        dst[4 * x + 1] = dst[4 * x + 2] = dst[4 * x + 3] = alpha;
      }
    }
    ArgbBuffer* buffer = &alpha_;
    ArgbBuffer pre_processed_buffer;
    if (config.alpha_quality > kMaxLossyQuality &&
        config.alpha_quality < kMaxQuality) {
      WP2_CHECK_STATUS(PreprocessNearLossless(alpha_, config, /*is_alpha=*/true,
                                              &pre_processed_buffer));
      buffer = &pre_processed_buffer;
    }

    // Encode the whole alpha losslessly.
    WP2_CHECK_ALLOC_OK(lossless_encode_info_.bits_per_pixel.resize(
        tile_rect_.width * tile_rect_.height));
    WP2_CHECK_ALLOC_OK(
        lossless_encode_info_.line_tokens.resize(tile_rect_.height));
    WP2_CHECK_STATUS(WP2L::EncodeImage(config_, *buffer, /*has_alpha=*/true,
                                       &lossless_enc_, &lossless_encode_info_));
    assert(lossless_encode_info_.line_tokens[alpha_.height - 1] ==
           lossless_enc_.NumTokens());
  }

  return WP2_STATUS_OK;
}

WP2Status AlphaWriter::CopyFrom(const AlphaWriter& other,
                                const BlockContext& context) {
  gparams_ = other.gparams_;
  tile_rect_ = other.tile_rect_;
  config_ = other.config_;
  gparams_ = other.gparams_;
  alpha_mode_ = other.alpha_mode_;

  if (alpha_mode_ == kAlphaModeLossless) {
    WP2_CHECK_OK(!other.alpha_.IsView(), WP2_STATUS_INVALID_PARAMETER);
    WP2_CHECK_STATUS(alpha_.CopyFrom(other.alpha_));
    WP2_CHECK_STATUS(lossless_enc_.Clone(other.lossless_enc_));
    WP2_CHECK_ALLOC_OK(
        lossless_encode_info_.CopyFrom(other.lossless_encode_info_));
    next_line_ = other.next_line_;
  } else {
    WP2_CHECK_STATUS(mode_predictor_.CopyFrom(other.mode_predictor_));
  }

  WP2_CHECK_OK(&context != other.context_, WP2_STATUS_INVALID_PARAMETER);
  context_ = &context;
  return WP2_STATUS_OK;
}

WP2Status AlphaWriter::WriteHeader(uint32_t num_coeffs_max,
                                   ANSEncBase* const enc) {
  ANSDebugPrefix prefix(enc, "Alpha");
  if (gparams_->maybe_use_lossy_alpha_) {
    enc->PutRValue(alpha_mode_, kAlphaModeNum, "alpha_mode");
  } else {
    assert(alpha_mode_ == kAlphaModeLossless);
  }
  return WP2_STATUS_OK;
}

WP2Status AlphaWriter::DecideAlpha(CodedBlock* const cb,
                                   const ResidualWriter& residual_writer,
                                   Counters* const counters) {
  if (alpha_mode_ == kAlphaModeLossless) {
    cb->alpha_mode_ = kBlockAlphaLossless;
  } else {
    WP2_CHECK_STATUS(ProcessLossy(
        cb, residual_writer, lossless_encode_info_.bits_per_pixel, counters));
  }

  return WP2_STATUS_OK;
}

WP2Status AlphaWriter::Record(const CodedBlock& cb) { return WP2_STATUS_OK; }

void AlphaWriter::WriteBlockBeforeCoeffs(const CodedBlock& cb,
                                         SymbolManager* const sm,
                                         ANSEncBase* const enc) {
  if (alpha_mode_ == kAlphaModeLossy) {
    mode_predictor_.Write(cb, enc, sm);
  }
}

WP2Status AlphaWriter::ResetRecord() { return WP2_STATUS_OK; }

WP2Status AlphaWriter::Write(const CodedBlock& cb, ANSEnc* const enc) {
  ANSDebugPrefix prefix(enc, "Alpha");
  if (alpha_mode_ == kAlphaModeLossless) {
    WP2_CHECK_STATUS(WriteLossless(cb, enc));
  }
  return enc->GetStatus();
}

WP2Status AlphaWriter::WriteLossless(const CodedBlock& cb, ANSEnc* const enc) {
  const uint32_t max_line = std::min(cb.y_pix() + cb.h_pix(), alpha_.height);
  if (next_line_ < max_line) {
    // Write as many new lines as we need for this block.
    const uint32_t lines_to_write = max_line - next_line_;

    const uint32_t start_token =
        (next_line_ > 0) ? lossless_encode_info_.line_tokens[next_line_ - 1]
                         : 0;
    const uint32_t n_tokens =
        lossless_encode_info_.line_tokens[next_line_ + lines_to_write - 1] -
        start_token;
    if (n_tokens > 0) {
      WP2_CHECK_STATUS(enc->AppendTokens(lossless_enc_, start_token, n_tokens));
    }

    next_line_ = max_line;
  }

  return WP2_STATUS_OK;
}

BlockAlphaMode AlphaWriter::GetBlockAlphaMode(const CodedBlock& cb,
                                              const Plane16& alpha) const {
  bool all_opaque = true;
  bool all_transp = true;
  for (uint32_t y = 0, y_offset = cb.y_pix();
       y < cb.h_pix() && y_offset < tile_rect_.height; ++y, ++y_offset) {
    const int16_t* const row = alpha.Row(y);
    for (uint32_t x = 0, x_offset = cb.x_pix();
         x < cb.w_pix() && x_offset < tile_rect_.width; ++x, ++x_offset) {
      const int16_t a = row[x];
      all_opaque = all_opaque && (a == kAlphaMax);
      all_transp = all_transp && (a == 0);
      if (!all_opaque && !all_transp) return kBlockAlphaLossy;
    }
  }
  assert(all_opaque || all_transp);
  return all_opaque ? kBlockAlphaFullOpaque : kBlockAlphaFullTransp;
}

WP2Status AlphaWriter::ProcessLossy(CodedBlock* const cb,
                                    const ResidualWriter& residual_writer,
                                    const Vector_f& bits_per_pixel,
                                    Counters* const counters) {
  cb->alpha_mode_ = GetBlockAlphaMode(*cb, cb->in_.A);
  if (cb->alpha_mode_ != kBlockAlphaLossy) {
    assert(cb->alpha_mode_ == kBlockAlphaFullOpaque ||
           cb->alpha_mode_ == kBlockAlphaFullTransp);
    cb->out_.A.Fill(cb->alpha_mode_ == kBlockAlphaFullTransp ? 0 : kAlphaMax);
    return WP2_STATUS_OK;
  }

  constexpr uint32_t kNumChannels = 4;  // Otherwise we would not be here.
  const bool forced_predictor = SetForcedPredictor(
      gparams_->a_preds_, config_, tile_rect_.x, tile_rect_.y, kAChannel, cb);

  const APredictors& preds = gparams_->a_preds_;
  const Segment& segment = gparams_->a_segments_[0];
  cb->GetCodingParams(kAChannel)->tf = kDctDct;

  if (forced_predictor) {
    // Predictor is set by EncoderConfig::info (debug).
  } else if (config_.speed == 1) {
    // Select the first predictor to save time.
    cb->SetPredictor(kAChannel, preds[0]);
    cb->pred_scores_[kAChannel] = 1.f;
  } else if (config_.speed == 2) {
    // Quickly select the best predictor.
    cb->FindBestPred(preds, kAChannel, /*reduced=*/false);
  } else {
    Vector<TransformPair> transforms;
    WP2_CHECK_ALLOC_OK(transforms.push_back(kDctDct));
    // Spend more time to select the best predictor.
    WP2_CHECK_STATUS(cb->FindBestPredTf(
        config_, tile_rect_, preds, segment, *context_, kAChannel, kNumChannels,
        /*reduced=*/false, transforms, counters));
  }
  // Quantizes and stores reconstructed pixels in cb->out_.
  cb->Quantize(segment.quant_a_, kAChannel, /*reduced_transform=*/false);

  // Update alpha_mode in case the reconstructed pixels are fully opaque or
  // fully transparent, even if the input pixels are not.
  cb->alpha_mode_ = GetBlockAlphaMode(*cb, cb->out_.A);
  if (cb->alpha_mode_ != kBlockAlphaLossy) return WP2_STATUS_OK;

  for (uint32_t tf_i = 0; tf_i < cb->GetNumTransforms(kAChannel); ++tf_i) {
    WP2_CHECK_STATUS(residual_writer.FindBestEncodingMethod(
        cb->tdim(kAChannel), cb->coeffs_[kAChannel][tf_i],
        cb->num_coeffs_[kAChannel][tf_i], cb->IsFirstCoeffDC(kAChannel),
        kAChannel, kNumChannels, counters->residuals(),
        &cb->method_[kAChannel][tf_i]));
  }

  return WP2_STATUS_OK;
}

}  // namespace WP2
