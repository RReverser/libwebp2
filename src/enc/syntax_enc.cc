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
// WP2 lossy encoding.
//

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <memory>

#include "src/enc/wp2_enc_i.h"
#include "src/utils/utils.h"
#include "src/wp2/format_constants.h"
#include "src/common/lossy/transforms.h"

namespace WP2 {

WP2Status SyntaxWriter::Init(ANSDictionaries* const dicts,
                             const EncoderConfig& config,
                             const GlobalParams& gparams, const YUVPlane& yuv,
                             ChromaSubsampling chroma_subsampling,
                             const Rectangle& tile_rect, uint32_t num_blocks,
                             bool use_aom_coeffs) {
  config_ = &config;
  dicts_ = dicts;
  recorded_blocks_ = 0;
  tile_rect_ = tile_rect;
  gparams_ = &gparams;
  assert(gparams_->IsOk());
  use_aom_coeffs_ = use_aom_coeffs;
  chroma_subsampling_ = chroma_subsampling;

  num_blocks_ = num_blocks;
  // Number of transforms is unknown yet, use a higher bound.
  num_transforms_ = std::min(
      num_blocks * 4, SizeBlocks(yuv.GetWidth()) * SizeBlocks(yuv.GetHeight()));

  WP2_CHECK_STATUS(
      residual_writer_.Init(use_aom_coeffs_, gparams_->maybe_use_lossy_alpha_));

  const uint32_t num_segments =
      std::min((uint32_t)gparams_->segments_.size(),
               GetMaxNumSegments(gparams_->explicit_segment_ids_,
                                 GetQualityHint(config_->quality),
                                 config_->partition_set));
  // Deal with segment ids.
  WP2_CHECK_STATUS(segment_ids_.InitWrite(
      num_segments, gparams_->explicit_segment_ids_, tile_rect.width));
  // Deal with context.
  WP2_CHECK_STATUS(
      context_.Init(use_aom_coeffs_, tile_rect.width, tile_rect.height));

  if (gparams_->has_alpha_) {
    WP2_CHECK_STATUS(
        alpha_writer_.Init(*config_, gparams, context_, yuv, tile_rect));
  }
  WP2_CHECK_STATUS(InitPass());

  if (kDebugPrintRate) {
    WP2_CHECK_ALLOC_OK(residual_rate_.resize(num_blocks));
  }

  return WP2_STATUS_OK;
}

WP2Status SyntaxWriter::CopyFrom(const SyntaxWriter& other,
                                 ANSDictionaries* const dicts) {
  config_ = other.config_;

  WP2_CHECK_OK(dicts != nullptr, WP2_STATUS_NULL_PARAMETER);
  WP2_CHECK_OK(dicts != other.dicts_, WP2_STATUS_INVALID_PARAMETER);
  dicts_ = dicts;

  WP2_CHECK_STATUS(context_.CopyFrom(other.context_));
  if (other.gparams_->has_alpha_) {
    WP2_CHECK_STATUS(alpha_writer_.CopyFrom(other.alpha_writer_, context_));
  }

  num_blocks_ = other.num_blocks_;
  num_transforms_ = other.num_transforms_;
  recorded_blocks_ = other.recorded_blocks_;
  tile_rect_ = other.tile_rect_;
  chroma_subsampling_ = other.chroma_subsampling_;
  gparams_ = other.gparams_;
  use_aom_coeffs_ = other.use_aom_coeffs_;

  WP2_CHECK_STATUS(segment_ids_.CopyFrom(other.segment_ids_));
  WP2_CHECK_STATUS(
      symbol_writer_.CopyFrom(other.symbol_writer_, *dicts_, *other.dicts_));
  WP2_CHECK_STATUS(symbol_recorder_.CopyFrom(other.symbol_recorder_));
  WP2_CHECK_STATUS(counters_.CopyFrom(other.counters_, symbol_recorder_));
  residual_writer_.CopyFrom(other.residual_writer_);

  WP2_CHECK_ALLOC_OK(residual_rate_.copy_from(other.residual_rate_));
  return WP2_STATUS_OK;
}

WP2Status SyntaxWriter::InitPass() {
  SymbolsInfo symbols_info;
  WP2_CHECK_STATUS(symbols_info.InitLossy(
      gparams_->segments_.size(), gparams_->partition_set_,
      gparams_->maybe_use_lossy_alpha_, GetQualityHint(config_->quality),
      use_aom_coeffs_));
  symbols_info.SetMinMax(kSymbolModeY, /*min=*/0,
                         /*max=*/gparams_->y_preds_.GetMaxMode());
  if (gparams_->maybe_use_lossy_alpha_) {
    symbols_info.SetMinMax(kSymbolModeA, /*min=*/0,
                           /*max=*/gparams_->a_preds_.GetMaxMode());
  }
  symbols_info.SetMinMax(kSymbolModeUV, /*min=*/0,
                         /*max=*/gparams_->uv_preds_.GetMaxMode());
  const uint32_t num_channels = (gparams_->maybe_use_lossy_alpha_ ? 4 : 3);
  for (Channel c : {kYChannel, kUChannel, kVChannel, kAChannel}) {
    if (c == kAChannel && !gparams_->maybe_use_lossy_alpha_) continue;
    WP2_CHECK_STATUS(symbols_info.SetMinMax(
        ResidualIO::GetCluster(c, num_channels), kSymbolDC, /*min=*/0,
        /*max=*/gparams_->GetMaxDCRange(c) - 1));
  }
  symbols_info.SetClusters(kSymbolDC, symbols_info.NumClusters(kSymbolDC));

  WP2_CHECK_STATUS(symbol_writer_.Init(symbols_info));
  WP2_CHECK_STATUS(symbol_writer_.Allocate());
  if (pass_number_ == 0) {
    WP2_CHECK_STATUS(symbol_recorder_.Allocate(
        symbols_info, /*num_records=*/num_blocks_ * kMaxBlockSize2));
  } else {
    WP2_CHECK_STATUS(symbol_recorder_.MakeBackup());
  }
  ++pass_number_;
  WP2_CHECK_STATUS(counters_.Init(symbol_recorder_));
  context_.Reset();
  WP2_CHECK_STATUS(ResetRecord());
  return WP2_STATUS_OK;
}

WP2Status SyntaxWriter::WriteHeader(ANSEncBase* const enc) {
  const uint32_t max_num_blocks =
      SizeBlocks(tile_rect_.width) * SizeBlocks(tile_rect_.height);
  // The minimum number of blocks is not accurate but it is not worth it to have
  // a more complex solution (neglictible pessimization).
  const uint32_t min_num_blocks = std::max(1u, max_num_blocks / kMaxBlockSize2);
  // This is also used for partitioning with an upper bound of 'num_blocks_'.
  assert((num_blocks_ == max_num_blocks) || (recorded_blocks_ == num_blocks_));
  ANSDebugPrefix prefix(enc, "GlobalHeader");

  // TODO(skal): differential-write of gparams_

  enc->PutBool(use_aom_coeffs_, "use_aom_coeffs");

  // Write dictionaries.
  PutLargeRange(num_blocks_, min_num_blocks, max_num_blocks, enc, "num_blocks");

  // If there's only one block, we don't need this global value, we'll signal
  // the subsampling mode for that block later.
  if (chroma_subsampling_ != ChromaSubsampling::kSingleBlock) {
    enc->PutRValue((int)chroma_subsampling_,
                   (int)ChromaSubsampling::kSingleBlock, "chroma_subsampling");
  }

  if (gparams_->explicit_segment_ids_) {
    WP2_CHECK_STATUS(symbol_writer_.WriteHeader(num_blocks_, kSymbolSegmentId,
                                                symbol_recorder_, "segment_id",
                                                enc, dicts_));
  }

  if (gparams_->use_rnd_mtx_) {
    WP2_CHECK_STATUS(symbol_writer_.WriteHeader(
        num_blocks_, kSymbolRndMtx0, symbol_recorder_, "rnd_mtx", enc, dicts_));
    WP2_CHECK_STATUS(symbol_writer_.WriteHeader(
        num_blocks_, kSymbolUseRandomMatrix, symbol_recorder_, "use_rnd_mtx",
        enc, dicts_));
  }
  WP2_CHECK_STATUS(
      symbol_writer_.WriteHeader(num_blocks_, kSymbolSplitTransform,
                                 symbol_recorder_, "split_tf", enc, dicts_));
  if (chroma_subsampling_ == ChromaSubsampling::kSingleBlock ||
      chroma_subsampling_ == ChromaSubsampling::kAdaptive) {
    WP2_CHECK_STATUS(symbol_writer_.WriteHeader(
        num_blocks_, kSymbolYuv420, symbol_recorder_, "yuv420", enc, dicts_));
  }
  WP2_CHECK_STATUS(symbol_writer_.WriteHeader(num_blocks_, kSymbolTransform,
                                              symbol_recorder_, "transform",
                                              enc, dicts_));

  WP2_CHECK_STATUS(symbol_writer_.WriteHeader(
      num_blocks_, kSymbolModeY, symbol_recorder_, "mode", enc, dicts_));
  WP2_CHECK_STATUS(symbol_writer_.WriteHeader(
      num_blocks_, kSymbolModeUV, symbol_recorder_, "mode", enc, dicts_));

  // TODO(maryla): not necessary if the SignalingCflPredictor isn't used.
  WP2_CHECK_STATUS(symbol_writer_.WriteHeader(
      2 * num_blocks_, kSymbolCflSlope, symbol_recorder_, "cfl", enc, dicts_));

  for (uint32_t i = 0; i < GetNumUniqueBounds(gparams_->partition_set_); ++i) {
    WP2_CHECK_STATUS(symbol_writer_.WriteHeader(
        /*cluster=*/i, num_blocks_, kSymbolBlockSize, symbol_recorder_,
        "block_size", enc, dicts_));
  }

  const uint32_t num_coeffs_max =
      max_num_blocks * kMinBlockSizePix * kMinBlockSizePix;
  const uint32_t num_coeffs_max_uv =
      (chroma_subsampling_ == ChromaSubsampling::k420) ? num_coeffs_max / 4
                                                       : num_coeffs_max;
  // We should always have a UV block.
  assert(num_coeffs_max_uv > 0);

  WP2_CHECK_STATUS(segment_ids_.WriteHeader(enc));

  if (gparams_->has_alpha_) {
    WP2_CHECK_STATUS(alpha_writer_.WriteHeader(num_coeffs_max, enc));

    if (alpha_writer_.GetAlphaMode() == kAlphaModeLossy) {
      WP2_CHECK_STATUS(symbol_writer_.WriteHeader(
          num_blocks_, kSymbolBlockAlphaMode, symbol_recorder_,
          "block_alpha_mode", enc, dicts_));
    }
    if (alpha_writer_.GetAlphaMode() != kAlphaModeLossless) {
      WP2_CHECK_STATUS(symbol_writer_.WriteHeader(
          num_blocks_, kSymbolModeA, symbol_recorder_, "mode", enc, dicts_));
    }
  }

  WP2_CHECK_STATUS(residual_writer_.WriteHeader(
      num_coeffs_max, num_coeffs_max_uv, num_transforms_,
      gparams_->has_alpha_ &&
          (alpha_writer_.GetAlphaMode() != kAlphaModeLossless),
      symbol_recorder_, dicts_, enc, &symbol_writer_));

  return enc->GetStatus();
}

void SyntaxWriter::WriteBlockBeforeCoeffs(const CodedBlock& cb,
                                          bool update_ymodes,
                                          SymbolManager* const sm,
                                          ANSEncBase* const enc) {
  segment_ids_.WriteId(cb, sm, enc);

  if (gparams_->has_alpha_) {
    alpha_writer_.WriteBlockBeforeCoeffs(cb, sm, enc);
  } else if (!gparams_->maybe_use_lossy_alpha_) {
    assert(!cb.HasLossyAlpha());
  }

  for (Channel channel : {kYChannel, kAChannel, kUChannel, kVChannel}) {
    if (channel == kAChannel && !cb.HasLossyAlpha()) continue;
    if (channel == kYChannel || channel == kAChannel) {
      WriteYAPredictors(cb, channel, &context_.ymodes(), update_ymodes, sm,
                        enc);
    } else {
      WriteUVPredictors(cb, channel, sm, enc);
    }
  }

  if (chroma_subsampling_ == ChromaSubsampling::kSingleBlock ||
      chroma_subsampling_ == ChromaSubsampling::kAdaptive) {
    sm->Process(kSymbolYuv420, cb.is420_, "yuv420", enc);
  }

  for (Channel channel : {kYChannel, kUChannel, kVChannel, kAChannel}) {
    if (channel == kAChannel && !cb.HasLossyAlpha()) continue;
    WriteSplitTransform(cb, channel, sm, enc);
    WriteHasCoeffs(cb, channel, sm, enc);
    WriteTransform(cb, channel, sm, enc);
  }

  if (!use_aom_coeffs_) {
    ANSDebugPrefix prefix(enc, "coeff_method");
    gparams_->segments_[cb.id_].StoreEncodingMethods(cb, sm, enc);
  }

  if (cb.mtx_set_ != nullptr) {
    sm->Process(kSymbolUseRandomMatrix, cb.use_mtx_, "use_rnd_mtx", enc);
    if (cb.use_mtx_) {
      sm->Process(kSymbolRndMtx0, cb.mtx_[kYChannel], "rnd_mtx_y", enc);
    }
  }
}

void SyntaxWriter::WriteYAPredictors(
    const CodedBlock& cb, Channel channel,
    const YModePredictor* const ymode_predictor, bool update_ymodes,
    SymbolManager* const sm, ANSEncBase* const enc) {
  ANSDebugPrefix prefix(enc, "pred_modes");
  const CodedBlock::CodingParams& params = cb.GetCodingParams(channel);
  ANSDebugPrefix prefix2(enc, (channel == kYChannel) ? "y" : "a");
  if (channel == kYChannel && cb.y_context_is_constant_) return;
  if (channel == kYChannel) {
    assert(ymode_predictor != nullptr);
    const_cast<YModePredictor*>(ymode_predictor)
        ->Write(cb, *params.pred, update_ymodes, enc, sm);
  } else {
    assert(channel == kAChannel);
    sm->Process(kSymbolModeA, params.pred->mode(), "pred", enc);
  }
  params.pred->WriteParams(cb, channel, sm, enc);
}

void SyntaxWriter::WriteUVPredictors(const CodedBlock& cb, Channel channel,
                                     SymbolManager* const sm,
                                     ANSEncBase* const enc) {
  ANSDebugPrefix prefix(enc, "pred_modes");
  const CodedBlock::CodingParams& params = cb.GetCodingParams(channel);
  if (channel == kUChannel) {
    sm->Process(kSymbolModeUV, params.pred->mode(), "uv", enc);
  }
  params.pred->WriteParams(cb, channel, sm, enc);
}

void SyntaxWriter::WriteSplitTransform(const CodedBlock& cb, Channel channel,
                                       SymbolManager* const sm,
                                       ANSEncBase* const enc) {
  if (channel == kYChannel &&
      GetSplitSize(cb.dim(), /*split=*/true) != cb.dim()) {
    sm->Process(kSymbolSplitTransform, cb.GetCodingParams(channel).split_tf,
                "split_tf", enc);
  } else {
    assert(!cb.GetCodingParams(channel).split_tf);
  }
}

void SyntaxWriter::WriteHasCoeffs(const CodedBlock& cb, Channel channel,
                                  SymbolManager* const sm,
                                  ANSEncBase* const enc) {
  bool block_has_coeffs = false;
  for (uint32_t tf_i = 0; tf_i < cb.GetNumTransforms(channel); ++tf_i) {
    const bool tf_has_coeffs = (cb.num_coeffs_[channel][tf_i] != 0);
    if (!block_has_coeffs && cb.GetCodingParams(channel).split_tf &&
        tf_i + 1 == cb.GetNumTransforms(channel)) {
      // At least one transform has at least one non-zero coeff, otherwise
      // 'split_tf' would be pointless.
      assert(tf_has_coeffs);
    } else {
      const uint32_t cluster = (uint32_t)channel;
      sm->Process(cluster, kSymbolHasCoeffs, tf_has_coeffs ? 1 : 0,
                  "has_coeffs", enc);
    }
    if (tf_has_coeffs) block_has_coeffs = true;
  }
}

void SyntaxWriter::WriteTransform(const CodedBlock& cb, Channel channel,
                                  SymbolManager* const sm,
                                  ANSEncBase* const enc) {
  const CodedBlock::CodingParams& params = cb.GetCodingParams(channel);
  const TransformPair implicit_tf = cb.GetImplicitTf(channel);
  if (implicit_tf != kUnknownTf) {
    // The transform is already known.
    assert(params.tf == implicit_tf);
  } else if (cb.HasCoeffs(channel)) {
    // TODO(maryla): try more ways to signal the transform which might be
    //               correlated to block size.
    sm->Process(kSymbolTransform, params.tf, "transform", enc);
  } else {
    // No need to signal the transform as there is no coeff.
  }
}

void SyntaxWriter::RecordBlockHeader(const CodedBlock& cb) {
  ANSEncNoop enc;
  WriteBlockBeforeCoeffs(cb, /*update_ymodes=*/true, &symbol_recorder_, &enc);
}

WP2Status SyntaxWriter::WriteBlocks(const Vector<CodedBlock>& cblocks,
                                    FrontMgrNxNBase* const mgr,
                                    ANSEnc* const enc) {
  assert(cblocks.size() == num_blocks_);
  // Note: we could call segment_ids_.InitMap(bw, bh) again here.
  mgr->Clear();
  WP2_CHECK_STATUS(
      residual_writer_.Init(use_aom_coeffs_, gparams_->maybe_use_lossy_alpha_));
  context_.Reset();
  for (uint16_t ind : mgr->SizeIndices()) {
    {
      ANSDebugPrefix prefix(enc, "BlockHeader");
      const CodedBlock& cb = cblocks[ind];
      WriteBlockSize(*mgr, cb.dim(), &symbol_writer_, enc);
      Block block_tmp;
      WP2_CHECK_ALLOC_OK(mgr->UseSize(cb.dim(), ind, &block_tmp));
      assert(block_tmp == cb.blk());
    }
    while (true) {
      FrontMgrNxNBase::Info info;
      if (!mgr->UseFinal(&info)) break;
      const CodedBlock& cb = cblocks[info.ind];
      WP2_CHECK_STATUS(WriteBlock(cb, info.ind, enc));
    }
  }
  return enc->GetStatus();
}

//------------------------------------------------------------------------------

WP2Status SyntaxWriter::FindBestEncodingMethods(CodedBlock* const cb) {
  const uint32_t num_channels = (gparams_->maybe_use_lossy_alpha_ ? 4 : 3);
  for (auto c : {kYChannel, kUChannel, kVChannel, kAChannel}) {
    if (c == kAChannel && !cb->HasLossyAlpha()) continue;
    for (uint32_t tf_i = 0; tf_i < cb->GetNumTransforms(c); ++tf_i) {
      if (use_aom_coeffs_) {
        ResidualWriter::SetGeometry(cb->num_coeffs_[c][tf_i],
                                    &cb->method_[c][tf_i]);
      } else {
        WP2_CHECK_STATUS(residual_writer_.FindBestEncodingMethod(
            cb->tdim(c), cb->coeffs_[c][tf_i], cb->num_coeffs_[c][tf_i],
            cb->IsFirstCoeffDC(c), c, num_channels, counters_.residuals(),
            &cb->method_[c][tf_i]));
      }
    }
  }
  return WP2_STATUS_OK;
}

WP2Status SyntaxWriter::DecideAlpha(CodedBlock* const cb) {
  return alpha_writer_.DecideAlpha(cb, residual_writer_, &counters_);
}

WP2Status SyntaxWriter::RecordAlpha(const CodedBlock& cb) {
  return alpha_writer_.Record(cb);
}

WP2Status SyntaxWriter::SetInitialSegmentIds() {
  const Rectangle padded_tile_rect = {tile_rect_.x, tile_rect_.y,
                                      Pad(tile_rect_.width, kPredWidth),
                                      Pad(tile_rect_.height, kPredWidth)};
  CodedBlock cb;
  cb.SetDimDefault(Block(0, 0, BLK_4x4));
  FrontMgr4x4 mgr;
  WP2_CHECK_STATUS(mgr.Init(config_->partition_set, config_->partition_snapping,
                            padded_tile_rect.width, padded_tile_rect.height));
  for (uint32_t y = 0; y < SizeBlocks(padded_tile_rect.height); ++y) {
    for (uint32_t x = 0; x < SizeBlocks(padded_tile_rect.width); ++x) {
      cb.SetXY(mgr, x, y);
      cb.id_ = AssignSegmentId(*config_, *gparams_, padded_tile_rect, cb.blk());
      segment_ids_.InitInitialSegmentId(cb, cb.id_);
    }
  }
  return WP2_STATUS_OK;
}

//------------------------------------------------------------------------------

WP2Status SyntaxWriter::Record(const CodedBlock& cb) {
  assert(recorded_blocks_ < num_blocks_);

  RecordBlockHeader(cb);
  const uint32_t num_channels = (gparams_->maybe_use_lossy_alpha_ ? 4 : 3);
  for (Channel channel : {kYChannel, kUChannel, kVChannel, kAChannel}) {
    if (channel == kAChannel && !cb.HasLossyAlpha()) continue;

    if (kDebugPrintRate) {
      float rate;
      WP2_CHECK_STATUS(cb.ResidualRate(context_, channel, num_channels,
                                       counters_.residuals(),
                                       counters_.residuals_aom(), &rate));
      residual_rate_[recorded_blocks_][channel] = rate;
    }

    WP2_CHECK_STATUS(residual_writer_.RecordCoeffs(
        cb, channel, &symbol_recorder_, context_.aom()));
  }

  ++recorded_blocks_;

  return WP2_STATUS_OK;
}

WP2Status SyntaxWriter::RecordSize(const FrontMgrNxNBase& mgr, BlockSize dim) {
  ANSEncNoop enc;
  WriteBlockSize(mgr, dim, &symbol_recorder_, &enc);
  return WP2_STATUS_OK;
}

WP2Status SyntaxWriter::ResetRecord() {
  recorded_blocks_ = 0;
  symbol_recorder_.ResetRecord(/*reset_backup=*/false);
  if (gparams_->has_alpha_) {
    WP2_CHECK_STATUS(alpha_writer_.ResetRecord());
  }
  return WP2_STATUS_OK;
}

WP2Status SyntaxWriter::WriteBlock(const CodedBlock& cb, uint32_t block_index,
                                   ANSEnc* const enc) {
  {
    ANSDebugPrefix prefix(enc, "BlockHeader");
    WriteBlockBeforeCoeffs(cb, /*update_ymodes=*/true, &symbol_writer_, enc);
  }

  float previous_cost;
  if (kDebugPrintRate) {
    previous_cost = enc->GetCost();
  }

  for (Channel c : {kYChannel, kUChannel, kVChannel, kAChannel}) {
    if (c == kAChannel && !cb.HasLossyAlpha()) continue;

    WP2_CHECK_STATUS(residual_writer_.WriteCoeffs(cb, c, enc, &symbol_writer_,
                                                  context_.aom()));
    if (kDebugPrintRate) {
      const float cost = enc->GetCost();
      const float cost_diff = cost - previous_cost;
      previous_cost = cost;

      const float rate = residual_rate_[block_index][c];
      if (cost_diff != 0 || rate != 0) {
        fprintf(stderr,
                "residual rate block %d %d channel %d\n%.2f vs real %.2f, "
                "%+.2f%%\n",
                cb.x(), cb.y(), c, rate, cost_diff,
                (rate / cost_diff - 1.f) * 100.f);
      }
    }
  }

  if (gparams_->has_alpha_) {
    WP2_CHECK_STATUS(alpha_writer_.Write(cb, enc));
  }
  if (cb.HasLossyAlpha()) assert(gparams_->has_alpha_);

#if defined(WP2_ENC_DEC_DEEP_MATCH)
  PutRawPixels(cb.blk(), cb.out_, enc);  // Debug
#endif
  return WP2_STATUS_OK;
}

}  // namespace WP2
