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
// WP2 lossy decoding of syntax.
//
// Author: Skal (pascal.massimino@gmail.com)

#include "src/common/constants.h"
#include "src/common/lossy/block_size.h"
#include "src/dec/tile_dec.h"
#include "src/dec/wp2_dec_i.h"
#include "src/utils/ans_utils.h"
#include "src/utils/front_mgr.h"
#include "src/utils/utils.h"
#include "src/wp2/decode.h"

namespace WP2 {

//------------------------------------------------------------------------------

SyntaxReader::SyntaxReader(ANSDec* const dec, const Rectangle& rect)
    : dec_(dec), width_(rect.width), height_(rect.height) {}

WP2Status SyntaxReader::LoadDictionaries() {
  const uint32_t max_num_blocks = SizeBlocks(width_) * SizeBlocks(height_);
  if (gparams_->explicit_segment_ids_) {
    WP2_CHECK_STATUS(
        sr_.ReadHeader(num_blocks_, kSymbolSegmentId, "segment_id"));
  }
  if (gparams_->use_rnd_mtx_) {
    WP2_CHECK_STATUS(sr_.ReadHeader(num_blocks_, kSymbolRndMtx0, "rnd_mtx"));
    WP2_CHECK_STATUS(
        sr_.ReadHeader(num_blocks_, kSymbolUseRandomMatrix, "use_rnd_mtx"));
  }
  WP2_CHECK_STATUS(
      sr_.ReadHeader(num_blocks_, kSymbolSplitTransform, "split_tf"));
  if (chroma_subsampling_ == ChromaSubsampling::kSingleBlock ||
      chroma_subsampling_ == ChromaSubsampling::kAdaptive) {
    WP2_CHECK_STATUS(sr_.ReadHeader(num_blocks_, kSymbolYuv420, "yuv420"));
  }
  WP2_CHECK_STATUS(sr_.ReadHeader(num_blocks_, kSymbolTransform, "transform"));
  sr_.SetRange(kSymbolModeY, /*min=*/0,
               /*max=*/gparams_->y_preds_.GetMaxMode());
  WP2_CHECK_STATUS(sr_.ReadHeader(num_blocks_, kSymbolModeY, "mode"));
  sr_.SetRange(kSymbolModeUV, /*min=*/0,
               /*max=*/gparams_->uv_preds_.GetMaxMode());
  WP2_CHECK_STATUS(sr_.ReadHeader(num_blocks_, kSymbolModeUV, "mode"));

  WP2_CHECK_STATUS(sr_.ReadHeader(2 * num_blocks_, kSymbolCflSlope, "cfl"));

  for (uint32_t i = 0; i < GetNumUniqueBounds(gparams_->partition_set_); ++i) {
    WP2_CHECK_STATUS(sr_.ReadHeader(/*cluster=*/i, num_blocks_,
                                    kSymbolBlockSize, "block_size"));
  }

  const uint32_t num_coeffs_max =
      max_num_blocks * kMinBlockSizePix * kMinBlockSizePix;
  const uint32_t num_coeffs_max_uv =
      (chroma_subsampling_ == ChromaSubsampling::k420) ? num_coeffs_max / 4
                                                       : num_coeffs_max;

  WP2_CHECK_STATUS(
      segment_ids_.ReadHeader(dec_, (uint32_t)gparams_->segments_.size(),
                              gparams_->explicit_segment_ids_, width_));

  if (gparams_->has_alpha_) {
    WP2_CHECK_STATUS(alpha_reader_->ReadHeader(*gparams_));

    if (alpha_reader_->GetAlphaMode() == kAlphaModeLossy) {
      WP2_CHECK_STATUS(sr_.ReadHeader(num_blocks_, kSymbolBlockAlphaMode,
                                      "block_alpha_mode"));
    }
    if (alpha_reader_->GetAlphaMode() != kAlphaModeLossless) {
      sr_.SetRange(kSymbolModeA, /*min=*/0, /*max=*/kAPredModeNum - 1);
      WP2_CHECK_STATUS(sr_.ReadHeader(num_blocks_, kSymbolModeA, "mode"));
    }
  }

  WP2_CHECK_STATUS(residual_reader_.ReadHeader(
      &sr_, num_coeffs_max, num_coeffs_max_uv, num_transforms_,
      gparams_->maybe_use_lossy_alpha_,
      gparams_->has_alpha_ &&
          (alpha_reader_->GetAlphaMode() != kAlphaModeLossless)));

  return WP2_STATUS_OK;
}

//------------------------------------------------------------------------------

WP2Status SyntaxReader::ReadPredModes(CodedBlock* const cb) {
  ANSDebugPrefix prefix(dec_, "pred_modes");
  cb->y_context_is_constant_ = cb->ContextIsConstant(kYChannel);
  for (Channel channel : {kYChannel, kAChannel}) {
    ANSDebugPrefix prefix2(dec_, (channel == kYChannel) ? "y" : "a");
    if (channel == kAChannel && !cb->HasLossyAlpha()) continue;
    if (channel == kYChannel && cb->y_context_is_constant_) {
      cb->SetLumaUniformPredictor(gparams_->y_preds_);
    } else {
      CodedBlock::CodingParams* const params = cb->GetCodingParams(channel);
      const Predictors& preds = gparams_->predictors(channel);
      uint32_t mode;
      if (channel == kYChannel) {
        mode = context_.ymodes().Read(*cb, preds, &sr_);
      } else {
        assert(channel == kAChannel);
        mode = sr_.Read(kSymbolModeA, "pred");
      }
      preds.GetFirstWithMode(mode)->ReadParams(cb, channel, &sr_, dec_);
      params->pred = preds.GetWithMode(mode, params->pred_sub_mode);
    }
  }
  const UVPredictors& uv_preds = gparams_->uv_preds_;
  const uint8_t uv_mode = sr_.Read(kSymbolModeUV, "uv");
  const Predictor* const p = uv_preds.GetFirstWithMode(uv_mode);
  p->ReadParams(cb, kUChannel, &sr_, dec_);
  p->ReadParams(cb, kVChannel, &sr_, dec_);
  const CodedBlock::CodingParams& params = *cb->GetCodingParams(kUChannel);
  cb->SetUVPredictor(uv_preds.GetWithMode(uv_mode, params.pred_sub_mode));
  return WP2_STATUS_OK;
}

//------------------------------------------------------------------------------

WP2Status SyntaxReader::ReadHeader(const BitstreamFeatures& features,
                                   const GlobalParams& gparams,
                                   const Tile& tile) {
  // TODO(skal): use differential coding to update gparam locally
  gparams_ = &gparams;

  if (gparams_->has_alpha_) {
    alpha_reader_.reset(new (WP2Allocable::nothrow) AlphaReader(dec_, tile));
    WP2_CHECK_ALLOC_OK(alpha_reader_ != nullptr);
    WP2_CHECK_STATUS(alpha_reader_->Allocate());
  }

  ANSDebugPrefix prefix(dec_, "GlobalHeader");

  use_aom_coeffs_ = dec_->ReadBool("use_aom_coeffs");
  WP2_CHECK_STATUS(context_.Init(use_aom_coeffs_, width_, height_));
  residual_reader_.Init(use_aom_coeffs_);

  WP2_CHECK_STATUS(symbols_info_.InitLossy(
      (uint32_t)gparams_->segments_.size(), gparams.partition_set_,
      gparams_->maybe_use_lossy_alpha_, features.quality_hint,
      use_aom_coeffs_));
  const uint32_t num_channels = (gparams_->maybe_use_lossy_alpha_ ? 4 : 3);
  for (Channel c : {kYChannel, kUChannel, kVChannel, kAChannel}) {
    if (c == kAChannel && !gparams_->maybe_use_lossy_alpha_) continue;
    WP2_CHECK_STATUS(symbols_info_.SetMinMax(
        ResidualIO::GetCluster(c, num_channels), kSymbolDC, /*min=*/0,
        /*max=*/gparams_->GetMaxDCRange(c) - 1));
  }
  symbols_info_.SetClusters(kSymbolDC, symbols_info_.NumClusters(kSymbolDC));
  WP2_CHECK_STATUS(sr_.Init(symbols_info_, dec_));
  WP2_CHECK_STATUS(sr_.Allocate());

  const uint32_t max_num_blocks = SizeBlocks(width_) * SizeBlocks(height_);
  const uint32_t min_num_blocks = std::max(1u, max_num_blocks / kMaxBlockSize2);
  num_blocks_ =
      ReadLargeRange(min_num_blocks, max_num_blocks, dec_, "num_blocks");
  num_transforms_ = std::min(num_blocks_ * 4, max_num_blocks);

  if (num_blocks_ == 1) {
    chroma_subsampling_ = ChromaSubsampling::kSingleBlock;
  } else {
    chroma_subsampling_ = (ChromaSubsampling)dec_->ReadRValue(
        (int)ChromaSubsampling::kSingleBlock, "chroma_subsampling");
  }

  WP2_CHECK_STATUS(LoadDictionaries());

  return dec_->GetStatus();
}

//------------------------------------------------------------------------------

void SyntaxReader::ReadSplitTransform(Channel channel, CodedBlock* const cb) {
  cb->GetCodingParams(channel)->split_tf =
      (channel == kYChannel &&
       GetSplitSize(cb->dim(), /*split=*/true) != cb->dim() &&
       (bool)sr_.Read(kSymbolSplitTransform, "split_tf"));
}

void SyntaxReader::ReadHasCoeffs(Channel channel, CodedBlock* const cb) {
  bool block_has_coeffs = false;
  for (uint32_t tf_i = 0; tf_i < cb->GetNumTransforms(channel); ++tf_i) {
    // TODO(yguyon): If '!is420_' was signaled, deduce 'has_coeffs'.
    // TODO(yguyon): If 'has_lossy_alpha_' was signaled, deduce 'has_coeffs'.
    bool tf_has_coeffs;
    if (!block_has_coeffs && cb->GetCodingParams(channel)->split_tf &&
        tf_i + 1 == cb->GetNumTransforms(channel)) {
      // At least one transform has at least one non-zero coeff, otherwise
      // 'split_tf' would be pointless.
      tf_has_coeffs = true;
    } else {
      const uint32_t cluster = (uint32_t)channel;
      tf_has_coeffs = (bool)sr_.Read(cluster, kSymbolHasCoeffs, "has_coeffs");
    }
    // 1 means "at least one coeff is not zero", 0 means "all coeffs are 0".
    cb->num_coeffs_[channel][tf_i] = tf_has_coeffs ? 1 : 0;
    if (tf_has_coeffs) block_has_coeffs = true;
  }
}

void SyntaxReader::ReadTransform(Channel channel, CodedBlock* const cb) {
  CodedBlock::CodingParams* params = cb->GetCodingParams(channel);
  const TransformPair implicit_tf = cb->GetImplicitTf(channel);
  if (implicit_tf != kUnknownTf) {
    // The transform is already known, no need to signal it.
    params->tf = implicit_tf;
  } else if (cb->HasCoeffs(channel)) {
    params->tf = (TransformPair)sr_.Read(kSymbolTransform, "transform");
  } else {
    // 0 coeff, transform does not matter.
    // TODO(maryla): we could have a noop transform to avoid doing any
    //               computation.
    params->tf = kIdentityIdentity;
  }
}

void SyntaxReader::ReadIs420(CodedBlock* const cb) {
  dec_->PushBitTracesCustomPrefix("is_420");
  if (chroma_subsampling_ == ChromaSubsampling::kSingleBlock ||
      chroma_subsampling_ == ChromaSubsampling::kAdaptive) {
    cb->is420_ = (bool)sr_.Read(kSymbolYuv420, "yuv420");
  } else if (chroma_subsampling_ == ChromaSubsampling::k420) {
    cb->is420_ = true;
  } else if (chroma_subsampling_ == ChromaSubsampling::k444) {
    cb->is420_ = false;
  } else {
    assert(false);
  }
  dec_->PopBitTracesCustomPrefix("is_420");
}

WP2Status SyntaxReader::GetBlock(CodedBlock* const cb, YUVPlane* const out,
                                 BlockInfo* const info) {
  // Read some block data.
  dec_->AddDebugPrefix("BlockHeader");

  // Initialize the views now that 'blk_.dim' is known.
  cb->SetContextInput(*out);  // Predict from other reconstructed pixels.
  cb->ResetContextCache();
  cb->SetReconstructedOutput(out);

  dec_->PushBitTracesCustomPrefix("segment_id");
  segment_ids_.ReadId(&sr_, cb);  // Look at surrounding availability.
  dec_->PopBitTracesCustomPrefix("segment_id");

  if (gparams_->has_alpha_) {
    alpha_reader_->GetBlockHeader(&sr_, cb);
  } else {
    cb->alpha_mode_ = kBlockAlphaLossless;
  }

  // Decode the predictor used and thus get access to the implicit tf, if any.
  dec_->PushBitTracesCustomPrefix("mode");
  WP2_CHECK_STATUS(ReadPredModes(cb));
  dec_->PopBitTracesCustomPrefix("mode");

  // Bound the number of coeffs to decode for chroma.
  ReadIs420(cb);

  // If there is at least one coeff to read and if the transform is not
  // implicit, decode the transform. The coeff method may differ depending on
  // the transform used.
  dec_->PushBitTracesCustomPrefix("transforms");
  for (Channel channel : {kYChannel, kUChannel, kVChannel, kAChannel}) {
    if (channel == kAChannel && !cb->HasLossyAlpha()) continue;
    ReadSplitTransform(channel, cb);
    ReadHasCoeffs(channel, cb);
    ReadTransform(channel, cb);
  }
  dec_->PopBitTracesCustomPrefix("transforms");

  // Signal the method with which coeffs were encoded.
  if (!use_aom_coeffs_) {
    ANSDebugPrefix prefix(dec_, "coeff_method");
    const Segment& segment = gparams_->segments_[cb->id_];
    WP2_CHECK_STATUS(segment.ReadEncodingMethods(&sr_, cb));
  }

  // Read random matrix.
  // TODO(skal): check rnd_mtx behavior when no coeff and/or identity transform
  dec_->PushBitTracesCustomPrefix("rnd_mtx");
  if (gparams_->use_rnd_mtx_) {
    cb->use_mtx_ = (bool)sr_.Read(kSymbolUseRandomMatrix, "use_rnd_mtx");
    if (cb->use_mtx_) {
      cb->mtx_[kYChannel] = (uint8_t)sr_.Read(kSymbolRndMtx0, "rnd_mtx_y");
      // TODO(skal): remove this bitstream possibility for error.
      WP2_CHECK_OK(gparams_->mtx_set_.IsOk(cb->mtx_[kYChannel], cb->dim()),
                   WP2_STATUS_BITSTREAM_ERROR);
    }
  }
  dec_->PopBitTracesCustomPrefix("rnd_mtx");

  dec_->PopDebugPrefix();

  // Read coefficients.
  for (Channel c : {kYChannel, kUChannel, kVChannel, kAChannel}) {
    if (c == kAChannel && !cb->HasLossyAlpha()) continue;
    dec_->PushBitTracesCustomPrefix(kCoeffsStr[c]);
    WP2_CHECK_STATUS(
        residual_reader_.ReadCoeffs(c, dec_, &sr_, cb, context_.aom(), info));
    dec_->PopBitTracesCustomPrefix(kCoeffsStr[c]);
  }

  // alpha
  if (gparams_->has_alpha_) {
    dec_->PushBitTracesCustomPrefix("alpha");
    WP2_CHECK_STATUS(alpha_reader_->GetBlock(cb));
    dec_->PopBitTracesCustomPrefix("alpha");
  }

  return dec_->GetStatus();
}

BlockSize SyntaxReader::GetBlockSize(const FrontMgrNxNBase& mgr) {
  ANSDebugPrefix prefix(dec_, "BlockHeader");
  return ReadBlockSize(mgr, &sr_);
}

//------------------------------------------------------------------------------

}  // namespace WP2
