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
// Block position/size scoring functions.
//
// Author: Yannis Guyon (yguyon@google.com)

#include "src/enc/partition_score_func.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <numeric>

#include "src/common/integral.h"
#include "src/dsp/dsp.h"
#include "src/dsp/math.h"
#include "src/enc/analysis.h"
#include "src/enc/partitioner.h"
#include "src/enc/wp2_enc_i.h"
#include "src/wp2/format_constants.h"

namespace WP2 {

struct BaseSlope {
  float base, slope;
};

// Maps 'config.quality' from [0, 95] to [min, (min + slope)].
static float MapQuality(const EncoderConfig& config, float min, float slope) {
  const float x = 1.f * config.quality / kMaxLossyQuality;
  return std::max(0.f, min + x * slope);
}

//------------------------------------------------------------------------------

WP2Status PartitionScoreFunc::Init(const EncoderConfig& config,
                                   const Rectangle& tile_rect,
                                   const YUVPlane& yuv,
                                   const GlobalParams& gparams) {
  config_ = &config;
  gparams_ = &gparams;
  tile_rect_ = tile_rect;
  src_ = &yuv;
  num_block_cols_ = SizeBlocks(yuv.Y.w_);
  num_block_rows_ = SizeBlocks(yuv.Y.h_);
  WP2_CHECK_STATUS(ClearVDebug());
  return WP2_STATUS_OK;
}

WP2Status PartitionScoreFunc::Use(const Block& block) { return WP2_STATUS_OK; }

//------------------------------------------------------------------------------

WP2Status BlockScoreFunc::Init(const EncoderConfig& config,
                               const Rectangle& tile_rect, const YUVPlane& yuv,
                               const GlobalParams& gparams) {
  WP2EncDspInit();
  WP2_CHECK_STATUS(PartitionScoreFunc::Init(config, tile_rect, yuv, gparams));

  const ChromaSubsampling chroma_subsampling =
      DecideChromaSubsampling(*config_, /*more_than_one_block=*/true);
  const bool use_aom_coeffs = DecideAOMCoeffs(*config_, tile_rect_);
  WP2_CHECK_STATUS(DecideTransforms(config, &transforms_, &transforms_subset_));

  // Store the reconstructed pixels of the temporary and final blocks.
  WP2_CHECK_STATUS(
      buffer_.Resize(src_->Y.w_, src_->Y.h_, /*pad=*/1, src_->HasAlpha()));
  WP2_CHECK_STATUS(front_mgr_.Init(config_->partition_set,
                                   config_->partition_snapping,
                                   tile_rect_.width, tile_rect_.height));

  // Initialize the instances recording only final blocks.
  WP2_CHECK_STATUS(syntax_writer_.Init(
      &dicts_, *config_, *gparams_, yuv, chroma_subsampling, tile_rect,
      num_block_cols_ * num_block_rows_, use_aom_coeffs));
  WP2_CHECK_STATUS(syntax_writer_.SetInitialSegmentIds());
  WP2_CHECK_STATUS(syntax_writer_.InitPass());

  if (DCDiffusionMap::GetDiffusion(config_->error_diffusion) > 0) {
    WP2_CHECK_STATUS(dc_error_u_.Init(tile_rect_.width));
    WP2_CHECK_STATUS(dc_error_v_.Init(tile_rect_.width));
  }

  // Initialize the cache.
  WP2_CHECK_STATUS(context_.Init(use_aom_coeffs, yuv.Y.w_, yuv.Y.h_));

  return WP2_STATUS_OK;
}

WP2Status BlockScoreFunc::GetScore(const Block& block, float* const score) {
  return GetScore(&block, /*num_blocks=*/1, score);
}

WP2Status BlockScoreFunc::GetScore(const Block blocks[4], uint32_t num_blocks,
                                   float* const score) {
  // Copy the final state into the temporary scratch one.
  WP2_CHECK_STATUS(tmp_dicts_.CopyFrom(dicts_));
  WP2_CHECK_STATUS(tmp_syntax_writer_.CopyFrom(syntax_writer_, &tmp_dicts_));
  if (DCDiffusionMap::GetDiffusion(config_->error_diffusion) > 0) {
    WP2_CHECK_STATUS(tmp_dc_error_u_.CopyFrom(dc_error_u_));
    WP2_CHECK_STATUS(tmp_dc_error_v_.CopyFrom(dc_error_v_));
  }

  float total_rate = 0.f, total_disto = 0.f, rate[4], disto[4];
  assert(num_blocks >= 1 && num_blocks <= 4);
  for (uint32_t i = 0; i < num_blocks; ++i) {
    const Block& block = blocks[i];

    // Encode the block.
    Block tmp_blk;
    if (!front_mgr_.SetNextBlockPosition(block.x(), block.y())) assert(false);
    assert(front_mgr_.TryGetNextBlock(block.dim(), &tmp_blk));
    assert(block == tmp_blk);
    tmp_cb_.SetDim(block, front_mgr_);
    WP2_CHECK_STATUS(EncodeBlock(front_mgr_, &tmp_cb_, &tmp_syntax_writer_,
                                 &tmp_dc_error_u_, &tmp_dc_error_v_, &buffer_));

    // Write the bits.
    ANSEnc enc;
    WP2_CHECK_STATUS(tmp_syntax_writer_.WriteHeader(&enc));
    const float header_rate = enc.GetCost(tmp_dicts_);
    WP2_CHECK_STATUS(
        WriteBlock(front_mgr_, tmp_cb_, &tmp_syntax_writer_, &enc));
    const uint32_t num_pixels = block.rect_pix().GetArea();
    // Exclude the header to prevent early decision from impacting later blocks.
    rate[i] = (enc.GetCost(tmp_dicts_) - header_rate) / num_pixels;
    assert(rate[i] >= 0.f);
    total_rate += rate[i];

    // Compute the distortion per pixel.
    disto[i] = 0.f;
    constexpr float disto_scale[] = {0.4f, 0.2f, 0.2f, 0.2f};
    for (Channel c : {kYChannel, kUChannel, kVChannel, kAChannel}) {
      if (c == kAChannel && !gparams_->has_alpha_) continue;
      disto[i] += disto_scale[c] * tmp_cb_.GetDisto(c, tile_rect_);
    }
    disto[i] /= num_pixels;  // Per pixel is not necessary, just nicer debug.
    total_disto += disto[i];

    // Register the blocks in the 'front_mgr_' except the last one.
    if (i + 1 < num_blocks) {
      if (!front_mgr_.UseSize(block.dim(), /*ind=*/0, &tmp_blk)) assert(false);
      assert(block == tmp_blk);
      front_mgr_.Use(block);
    }
  }

  // Unregister the blocks in the 'front_mgr_' except the last one.
  for (uint32_t i = num_blocks - 1; i-- > 0;) {
    front_mgr_.UndoUse(blocks[i]);
    front_mgr_.UndoUseSize(blocks[i]);
  }

  total_rate /= num_blocks;  // Average per pixel values.
  total_disto /= num_blocks;

  // Estimate a score from the written bits and the distortion.
  const float lambda = MapQuality(*config_, 10.f, -9.f);  // Empirical.
  constexpr float kNiceScale = 0.01f;  // Has no impact on the result.
  *score = 1.0f / (1.0f + kNiceScale * (lambda * total_rate + total_disto));
  RegisterScoreForVDebug(blocks, num_blocks, rate, disto, total_rate,
                         total_disto, *score);
  return WP2_STATUS_OK;
}

WP2Status BlockScoreFunc::Use(const Block& block) {
  if (!front_mgr_.SetNextBlockPosition(block.x(), block.y())) assert(false);
  tmp_cb_.SetDim(block, front_mgr_);
  // Write the final pixels for future context and rate computation.
  WP2_CHECK_STATUS(EncodeBlock(front_mgr_, &tmp_cb_, &syntax_writer_,
                               &dc_error_u_, &dc_error_v_, &buffer_));
  WP2_CHECK_ALLOC_OK(front_mgr_.UseSize(tmp_cb_.dim(),
                                        /*ind=*/0, /*block=*/nullptr));
  front_mgr_.Use(tmp_cb_.blk());
  return WP2_STATUS_OK;
}

WP2Status BlockScoreFunc::FindBestBlockParams(const FrontMgrNxNBase& front_mgr,
                                              const BlockContext& block_context,
                                              SyntaxWriter* const writer,
                                              DCDiffusionMap* const dc_error_u,
                                              DCDiffusionMap* const dc_error_v,
                                              CodedBlock* const cb) const {
  const Rectangle padded_tile_rect = {tile_rect_.x, tile_rect_.y,
                                      src_->GetWidth(), src_->GetHeight()};
  cb->id_ = AssignSegmentId(*config_, *gparams_, padded_tile_rect, cb->blk());
  const Segment& segment = gparams_->segments_.at(cb->id_);
  cb->mtx_set_ = gparams_->use_rnd_mtx_ ? &gparams_->mtx_set_ : nullptr;

  cb->ResetContextCache();
  cb->y_context_is_constant_ = cb->ContextIsConstant(kYChannel);

  WP2_CHECK_STATUS(cb->OptimizeModesLuma(
      *config_, tile_rect_, gparams_->maybe_use_lossy_alpha_,
      gparams_->y_preds_, segment, writer->context(), transforms_,
      transforms_subset_, writer->counters()));

  WP2_CHECK_STATUS(cb->OptimizeModesChroma(
      *config_, tile_rect_, gparams_->maybe_use_lossy_alpha_, front_mgr,
      gparams_->uv_preds_, segment, writer->context(),
      writer->chroma_subsampling(), dc_error_u, dc_error_v,
      writer->counters()));

  if (gparams_->has_alpha_) {
    WP2_CHECK_STATUS(writer->DecideAlpha(cb));
    if (!cb->HasLossyAlpha()) {
      // Consider no loss by copying original samples to the buffer.
      WP2_CHECK_STATUS(cb->out_.A.Copy(cb->in_.A, /*resize_if_needed=*/false));
    }
  }
  return WP2_STATUS_OK;
}

WP2Status BlockScoreFunc::EncodeBlock(const FrontMgrNxNBase& front_mgr,
                                      CodedBlock* const cb,
                                      SyntaxWriter* const syntax_writer,
                                      DCDiffusionMap* const dc_error_u,
                                      DCDiffusionMap* const dc_error_v,
                                      YUVPlane* const buffer) const {
  cb->SetRange(gparams_->transf_.GetYUVMin(), gparams_->transf_.GetYUVMax());
  cb->SetSrcInput(*src_);
  ContextCache pred_context;
  cb->SetContextInput(buffer_, &pred_context);
  cb->SetReconstructedOutput(buffer);

  // This is the slowest part: finding the best transform, predictor etc.
  WP2_CHECK_STATUS(FindBestBlockParams(front_mgr, context_, syntax_writer,
                                       dc_error_u, dc_error_v, cb));
  // CodedBlock::Quantize() should be called already.
  WP2_CHECK_STATUS(syntax_writer->FindBestEncodingMethods(cb));
  WP2_CHECK_STATUS(syntax_writer->RecordSize(front_mgr, cb->dim()));
  WP2_CHECK_STATUS(syntax_writer->Record(*cb));
  if (gparams_->has_alpha_) {
    WP2_CHECK_STATUS(syntax_writer->RecordAlpha(*cb));
  }
  return WP2_STATUS_OK;
}

WP2Status BlockScoreFunc::WriteBlock(const FrontMgrNxNBase& front_mgr,
                                     const CodedBlock& cb,
                                     SyntaxWriter* const syntax_writer,
                                     ANSEnc* const enc) const {
  assert(front_mgr.GetMaxFittingBlock().x() == cb.x() &&
         front_mgr.GetMaxFittingBlock().y() == cb.y());
  assert(front_mgr.GetMaxPossibleBlock().rect().Contains(cb.blk().rect()));
  {
    ANSDebugPrefix prefix(enc, "BlockHeader");
    WriteBlockSize(front_mgr, cb.dim(), syntax_writer->symbol_writer(), enc);
  }
  WP2_CHECK_STATUS(syntax_writer->WriteBlock(cb, /*block_index=*/0, enc));
  return WP2_STATUS_OK;
}

//------------------------------------------------------------------------------

AreaScoreFunc::AreaScoreFunc(uint32_t area_width, uint32_t area_height)
    : area_width_(area_width),
      area_height_(area_height),
      area_front_mgr_(area_width, area_height),
      comp_(kMaxTileSize / kMinBlockSizePix, area_width / kMinBlockSizePix,
            area_height / kMinBlockSizePix) {}

WP2Status AreaScoreFunc::Init(const EncoderConfig& config,
                              const Rectangle& tile_rect, const YUVPlane& yuv,
                              const GlobalParams& gparams) {
  WP2_CHECK_STATUS(PartitionScoreFunc::Init(config, tile_rect, yuv, gparams));

  // This scoring function needs a strict block layout.
  const BlockSize max_block_size =
      GetSmallestBounds(config.partition_set, BLK_32x32);
  WP2_CHECK_OK(config.partition_snapping, WP2_STATUS_INVALID_CONFIGURATION);
  WP2_CHECK_OK(BlockWidthPix(max_block_size) <= area_width_ &&
                   BlockHeightPix(max_block_size) <= area_height_,
               WP2_STATUS_INVALID_CONFIGURATION);

  // Matches TileEncoder::LossyEncode() behavior.
  const ChromaSubsampling chroma_subsampling =
      DecideChromaSubsampling(*config_, /*more_than_one_block=*/true);
  const bool use_aom_coeffs = DecideAOMCoeffs(*config_, tile_rect_);
  WP2_CHECK_STATUS(DecideTransforms(config, &transforms_, &transforms_subset_));

  // Retrieve the default partition (multipass) to have a reference.
  {
    MultiScoreFunc score_func;
    MultiPassPartitioner partitioner(&score_func);

    EncoderConfig cfg_no_dbg = config;
    cfg_no_dbg.info = nullptr;  // Do not use or output any debugging data.
    WP2_CHECK_STATUS(score_func.Init(cfg_no_dbg, tile_rect_, *src_, *gparams_));
    WP2_CHECK_STATUS(
        partitioner.Init(cfg_no_dbg, *src_, tile_rect_, &score_func));
    // TODO(yguyon): Add forced blocks from 'config' to 'default_partition_'
    WP2_CHECK_STATUS(partitioner.GetBestPartition(&default_partition_));
    // Sort the blocks to find those in a given area faster.
    std::sort(default_partition_.begin(), default_partition_.end(), comp_);
  }

  // Initialize the instances recording only final blocks.
  WP2_CHECK_STATUS(syntax_writer_.Init(
      &dicts_, *config_, *gparams_, yuv, chroma_subsampling, tile_rect,
      num_block_cols_ * num_block_rows_, use_aom_coeffs));
  WP2_CHECK_STATUS(syntax_writer_.SetInitialSegmentIds());
  WP2_CHECK_STATUS(syntax_writer_.InitPass());

  WP2_CHECK_STATUS(context_.Init(use_aom_coeffs, yuv.Y.w_, yuv.Y.h_));

  if (DCDiffusionMap::GetDiffusion(config_->error_diffusion) > 0) {
    WP2_CHECK_STATUS(dc_error_u_.Init(tile_rect_.width));
    WP2_CHECK_STATUS(dc_error_v_.Init(tile_rect_.width));
  }

  // Store the reconstructed pixels of the current area partitioning and of
  // the final selected blocks.
  WP2_CHECK_STATUS(
      buffer_.Resize(src_->Y.w_, src_->Y.h_, /*pad=*/1, src_->HasAlpha()));

  // Setup the first area as the top-left corner.
  WP2_CHECK_STATUS(area_front_mgr_.Init(config_->partition_set,
                                        config_->partition_snapping,
                                        tile_rect_.width, tile_rect_.height));
  WP2_CHECK_STATUS(BeginArea(/*area_x=*/0, /*area_y=*/0));  // x,y within tile
  return WP2_STATUS_OK;
}

// Returns a score in [0:1] representing how much a 'disto/rate' pair is good
// compared to a reference, where 0 is discarded, 1 is way better than 'ref' and
// 1/3 is equivalent to 'ref'.
static float ComputeScore(const EncoderConfig& config, float ref_disto,
                          float ref_rate, float disto, float rate) {
  if (disto > ref_disto || rate > ref_rate) return 0.f;  // No mercy.
  if (ref_disto > 0.f) {
    disto /= ref_disto;  // Normalize.
  } else  {
    if (disto > 0.f) return 0.f;  // Discard, no way of assessing 'disto'/'ref'.
    disto = 1.f;  // Both values are 0 so set to 1 to signal it is equivalent.
  }
  if (ref_rate > 0.f) {
    rate /= ref_rate;  // Normalize.
  } else  {
    if (rate > 0.f) return 0.f;  // Discard, no way of assessing 'rate'/'ref'.
    rate = 1.f;  // Both values are 0 so set to 1 to signal it is equivalent.
  }
  const float delta = Clamp(config.quality / kMaxLossyQuality, 0.01f, 0.99f);
  return 1.f / (1.f + delta * disto + (1.f - delta) * rate);
}

WP2Status AreaScoreFunc::GetScore(const Block&, float* const) {
  return WP2_STATUS_UNSUPPORTED_FEATURE;  // Unused.
}

WP2Status AreaScoreFunc::GetAreaDefaultScore(VectorNoCtor<Block>* const blocks,
                                             float* const score) {
  Vector<CodedBlock> partition;
  WP2_CHECK_STATUS(GetAreaDefaultPartition(&partition));
  assert(!partition.empty());
  WP2_CHECK_STATUS(GetDistoRate(&partition, &default_disto_, &default_rate_));
  // Compute the score with itself as reference, to be easily comparable.
  *score = ComputeScore(*config_, default_disto_, default_rate_,
                        default_disto_, default_rate_);
  for (const CodedBlock& cb : partition) {
    WP2_CHECK_ALLOC_OK(blocks->push_back(cb.blk()));
  }

  RegisterScoreForVDebug(BLK_LAST, partition, *score, default_disto_,
                         default_rate_);
  return WP2_STATUS_OK;
}

WP2Status AreaScoreFunc::GetAreaGridScore(BlockSize block_size,
                                          VectorNoCtor<Block>* const blocks,
                                          float* const score) {
  Vector<CodedBlock> partition;
  // Fill 'partition' with as many 'block_size' as possible.
  // Fill the remaining space with blocks as big as possible.
  FrontMgrArea front_mgr(area_width_, area_height_);
  WP2_CHECK_STATUS(front_mgr.CopyFrom(area_front_mgr_));
  uint32_t num_block_units = 0;
  while (!front_mgr.Done()) {
    Block block;
    if (!front_mgr.TryGetNextBlock(block_size, &block)) {
      block = front_mgr.GetMaxFittingBlock();
    }
    if (!area_.Contains(block.x_pix(), block.y_pix())) break;
    WP2_CHECK_ALLOC_OK(front_mgr.UseSize(block.dim(), 0, nullptr));
    front_mgr.Use(block);

    WP2_CHECK_ALLOC_OK(partition.resize(partition.size() + 1));
    partition.back().SetDimDefault(block);
    WP2_CHECK_ALLOC_OK(blocks->push_back(block));
    num_block_units += block.rect().GetArea();
  }
  assert(num_block_units == SizeBlocks(area_.width) * SizeBlocks(area_.height));

  assert(default_disto_ >= 0.f && default_rate_ >= 0.f);
  float disto = 0.f, rate = 0.f;
  WP2_CHECK_STATUS(GetDistoRate(&partition, &disto, &rate));
  *score = ComputeScore(*config_, default_disto_, default_rate_, disto, rate);

  RegisterScoreForVDebug(block_size, partition, *score, disto, rate);
  return WP2_STATUS_OK;
}

WP2Status AreaScoreFunc::GetDistoRate(Vector<CodedBlock>* const area_blocks,
                                      float* const disto,
                                      float* const rate) const {
  *rate = *disto = 0.f;

  ANSDictionaries dicts;
  WP2_CHECK_STATUS(dicts.CopyFrom(dicts_));
  SyntaxWriter syntax_writer;
  WP2_CHECK_STATUS(syntax_writer.CopyFrom(syntax_writer_, &dicts));

  DCDiffusionMap dc_error_u, dc_error_v;
  if (DCDiffusionMap::GetDiffusion(config_->error_diffusion) > 0) {
    WP2_CHECK_STATUS(dc_error_u.CopyFrom(dc_error_u_));
    WP2_CHECK_STATUS(dc_error_v.CopyFrom(dc_error_v_));
  }

  // Encode all 'area_blocks' using previously finished areas and blocks in this
  // 'area_' as prediction context (stored in 'buffer_').
  {
    FrontMgrArea front_mgr(area_width_, area_height_);
    WP2_CHECK_STATUS(front_mgr.CopyFrom(area_front_mgr_));
    for (CodedBlock& cb : *area_blocks) {
      WP2_CHECK_OK(!front_mgr.Done(), WP2_STATUS_INVALID_PARAMETER);
      assert(cb.x() == front_mgr.GetMaxPossibleBlock().x() &&
             cb.y() == front_mgr.GetMaxPossibleBlock().y() &&
             front_mgr.GetMaxPossibleBlock().rect().Contains(cb.blk().rect()));
      cb.SetDim(cb.blk(), front_mgr);
      WP2_CHECK_STATUS(EncodeBlock(front_mgr, &cb, &syntax_writer, &dc_error_u,
                                   &dc_error_v, &buffer_));
      WP2_CHECK_ALLOC_OK(front_mgr.UseSize(cb.dim(),
                                           /*ind=*/0, /*block=*/nullptr));
      front_mgr.Use(cb.blk());
    }
  }

  // Now estimate the bits necessary to encode all blocks in 'area_'.
  {
    // Also write the headers so that symbols are correctly set up.
    ANSEnc enc;
    WP2_CHECK_STATUS(syntax_writer.WriteHeader(&enc));
    const float header_rate = enc.GetCost(dicts);
    FrontMgrArea front_mgr(area_width_, area_height_);
    WP2_CHECK_STATUS(front_mgr.CopyFrom(area_front_mgr_));

    for (const CodedBlock& cb : *area_blocks) {
      WP2_CHECK_STATUS(WriteBlock(front_mgr, cb, &syntax_writer, &enc));
      WP2_CHECK_ALLOC_OK(front_mgr.UseSize(cb.dim(),
                                           /*ind=*/0, /*block=*/nullptr));
      front_mgr.Use(cb.blk());
    }
    // Exclude the header to prevent early decisions from impacting later areas.
    *rate = (enc.GetCost(dicts) - header_rate) / area_.GetArea();
    assert(*rate >= 0.f);
  }

  // Measure the distortion of this tested partition of 'area_'.
  {
    constexpr float disto_scale[] = {0.4f, 0.2f, 0.2f, 0.2f};
    for (Channel c : {kYChannel, kUChannel, kVChannel, kAChannel}) {
      if (c == kAChannel && !gparams_->has_alpha_) continue;
      Plane16 src_area_view, buffer_area_view;
      WP2_CHECK_STATUS(src_area_view.SetView(src_->GetChannel(c), area_));
      WP2_CHECK_STATUS(buffer_area_view.SetView(buffer_.GetChannel(c), area_));
      *disto +=
          disto_scale[c] * WP2SumSquaredErrorBlock(
                               src_area_view.Row(0), src_area_view.Step(),
                               buffer_area_view.Row(0), buffer_area_view.Step(),
                               area_.width, area_.height);
    }
    *disto /= area_.GetArea();
  }
  return WP2_STATUS_OK;
}

WP2Status AreaScoreFunc::BeginArea(uint32_t area_x, uint32_t area_y) {
  // Make sure areas are done in order.
  assert((area_x == 0 && area_y == 0) ||
         (area_x == area_.x + area_width_ && area_y == area_.y) ||
         (area_x == 0 && area_y == area_.y + area_height_));
  area_ = Rectangle(area_x, area_y, area_width_, area_height_)
              .ClipWith({0, 0, tile_rect_.width, tile_rect_.height});  // no pad
  // Assuming all areas are done in order row by row, top and left contexts
  // outside this 'area_' are available, if any.

  // Next available block should match the current 'area_'.
  assert(area_front_mgr_.GetMaxFittingBlock().x_pix() == area_.x &&
         area_front_mgr_.GetMaxFittingBlock().y_pix() == area_.y);

  return WP2_STATUS_OK;
}

WP2Status AreaScoreFunc::Use(const Block& block) {
  CodedBlock cb;
  cb.SetDim(block, area_front_mgr_);
  assert(cb.blk() == block);
  // Write the final pixels for future context and rate computation.
  WP2_CHECK_STATUS(EncodeBlock(area_front_mgr_, &cb, &syntax_writer_,
                               &dc_error_u_, &dc_error_v_, &buffer_));
  WP2_CHECK_ALLOC_OK(area_front_mgr_.UseSize(cb.dim(), /*ind=*/0,
                                             /*block=*/nullptr));
  area_front_mgr_.Use(cb.blk());
  if (area_front_mgr_.Done()) {
    // All areas are complete.
    area_ = {};
    area_front_mgr_.Clear();
  } else {
    const Block max_block = area_front_mgr_.GetMaxPossibleBlock();
    if (!area_.Contains(max_block.x_pix(), max_block.y_pix())) {
      // This 'area_' is complete. Prepare the next one.
      uint32_t area_x = area_.x + area_width_, area_y = area_.y;  // Next col.
      if (area_x >= tile_rect_.width) {                           // Next row.
        area_x = 0;
        area_y += area_height_;
      }
      assert(area_x < tile_rect_.width && area_y < tile_rect_.height);
      WP2_CHECK_STATUS(BeginArea(area_x, area_y));
    }
  }
  default_disto_ = default_rate_ = -1;
  return WP2_STATUS_OK;
}

WP2Status AreaScoreFunc::GetAreaDefaultPartition(
    Vector<CodedBlock>* const area_blocks) const {
  // Find the blocks in 'default_partition_' belonging to the current 'area_'.
  const Block area_pos(area_.x / kMinBlockSizePix, area_.y / kMinBlockSizePix,
                       BLK_32x32);
  VectorNoCtor<Block>::const_iterator block_it = std::lower_bound(
      default_partition_.begin(), default_partition_.end(), area_pos, comp_);

  uint32_t num_block_units = 0;
  for (; block_it != default_partition_.end(); ++block_it) {
    if (!area_.Contains(block_it->x_pix(), block_it->y_pix())) break;
    WP2_CHECK_ALLOC_OK(area_blocks->resize(area_blocks->size() + 1));
    area_blocks->back().SetDimDefault(*block_it);
    num_block_units += block_it->rect().GetArea();
  }
  // If snapping is not enabled or if the 'area_' dimensions do not match it,
  // the default partition cannot be used as is.
  WP2_CHECK_OK(
      num_block_units == SizeBlocks(area_.width) * SizeBlocks(area_.height),
      WP2_STATUS_INVALID_CONFIGURATION);
  return WP2_STATUS_OK;
}

//------------------------------------------------------------------------------

WP2Status SubAreaScoreFunc::GetScore(const Block& block, float* const score) {
  if (default_block_.dim() == BLK_LAST) {
    // Get the default partitioning and score of the remaining non-final blocks.
    // This is done once for each block position within each area.
    Vector<CodedBlock> default_partition;
    WP2_CHECK_STATUS(GetAreaRemainingDefaultPartition(&default_partition));
    assert(!default_partition.empty() &&
           default_partition.front().dim() != BLK_LAST);
    // TODO(yguyon): Also include the 'area_used_blocks_' into the disto/rate
    WP2_CHECK_STATUS(GetDistoRate(&default_partition, &default_block_disto_,
                                  &default_block_rate_));
    const float default_score =
        ComputeScore(*config_, default_block_disto_, default_block_rate_,
                     default_block_disto_, default_block_rate_);
    RegisterScoreForVDebug(default_partition.front().blk(), default_partition,
                           default_score, default_block_disto_,
                           default_block_rate_);
    default_block_ = default_partition.front().blk();
  }
  // 'default_block_' now contains the size of the first block among the
  // remaining non-final ones given by the default partitioning.
  assert(default_block_.x() == block.x() && default_block_.y() == block.y());

  float disto = 0.f, rate = 0.f;
  if (block == default_block_) {
    disto = default_block_disto_;
    rate = default_block_rate_;
    *score = ComputeScore(*config_, default_block_disto_, default_block_rate_,
                          disto, rate);
  } else {
    // 'block' here is the currently evaluated size for a given position.
    // 'area_used_blocks_' represent the final blocks previously encoded and
    // recorded. 'area_remaining_blocks' exist to fill the 'area_' partition and
    // compare the same surface by rate-distortion with the default partition.
    Vector<CodedBlock> area_remaining_blocks;
    WP2_CHECK_ALLOC_OK(area_remaining_blocks.resize(1));
    area_remaining_blocks.back().SetDimDefault(block);  // Force it.
    WP2_CHECK_STATUS(GetAreaRemainingDefaultPartition(&area_remaining_blocks));

    WP2_CHECK_STATUS(GetDistoRate(&area_remaining_blocks, &disto, &rate));
    *score = ComputeScore(*config_, default_block_disto_, default_block_rate_,
                          disto, rate);
    RegisterScoreForVDebug(block, area_remaining_blocks, *score, disto, rate);
  }
  return WP2_STATUS_OK;
}

WP2Status SubAreaScoreFunc::BeginArea(uint32_t area_x, uint32_t area_y) {
  WP2_CHECK_STATUS(AreaScoreFunc::BeginArea(area_x, area_y));
  area_used_blocks_.clear();
  return WP2_STATUS_OK;
}

WP2Status SubAreaScoreFunc::Use(const Block& block) {
  WP2_CHECK_ALLOC_OK(area_used_blocks_.resize(area_used_blocks_.size() + 1));
  area_used_blocks_.back().SetDim(block, area_front_mgr_);
  WP2_CHECK_STATUS(AreaScoreFunc::Use(block));
  default_block_ = Block();  // Reset.
  default_block_disto_ = default_block_rate_ = 0.f;
  return WP2_STATUS_OK;
}

WP2Status SubAreaScoreFunc::GetAreaRemainingDefaultPartition(
    Vector<CodedBlock>* const area_remaining_blocks) const {
  // It would be faster to reuse the 'AreaScoreFunc::default_partition_' (or
  // part of it) but it is not always compatible with the already selected
  // 'area_used_blocks_' so for simplicity it is always recomputed for each
  // block size of each block position.
  MultiScoreFunc score_func;
  MultiPassPartitioner partitioner(&score_func);

  EncoderConfig config = *config_;
  config.info = nullptr;  // Remove any debugging input/output.
  WP2_CHECK_STATUS(score_func.Init(config, tile_rect_, *src_, *gparams_));
  WP2_CHECK_STATUS(partitioner.Init(config, *src_, tile_rect_, &score_func));

  // 'blocks' will first contain all irrelevant blocks that will be forced into
  // the 'partitioner' to extract only the interesting ones.
  VectorNoCtor<Block> blocks;
  // Force all blocks external to the current 'area_' (their size and count do
  // not matter).
  {
    FrontMgrArea front_mgr(area_width_, area_height_);
    WP2_CHECK_STATUS(front_mgr.Init(config_->partition_set,
                                    config_->partition_snapping,
                                    tile_rect_.width, tile_rect_.height));
    uint32_t surface_kept = 0;
    while (!front_mgr.Done()) {
      const Block max_block = front_mgr.GetMaxFittingBlock();
      WP2_CHECK_ALLOC_OK(front_mgr.UseSize(max_block.dim(), 0, nullptr));
      front_mgr.Use(max_block);
      if (!area_.Contains(max_block.x_pix(), max_block.y_pix())) {
        WP2_CHECK_ALLOC_OK(blocks.push_back(max_block));
      } else {
        surface_kept += max_block.rect().GetArea();
      }
    }
    WP2_CHECK_OK(
        surface_kept == SizeBlocks(area_.width) * SizeBlocks(area_.height),
        WP2_STATUS_INVALID_PARAMETER);
  }
  // Force all already selected final blocks.
  for (const CodedBlock& cb : area_used_blocks_) {
    WP2_CHECK_ALLOC_OK(blocks.push_back(cb.blk()));
  }
  // Force the block under review, if any.
  for (const CodedBlock& cb : *area_remaining_blocks) {
    WP2_CHECK_ALLOC_OK(blocks.push_back(cb.blk()));
  }
  const size_t num_default_blocks_to_ignore = blocks.size();
  const size_t num_already_added_area_blocks = area_remaining_blocks->size();

  // Get the default partitioning of the remaining empty spaces into 'blocks'.
  WP2_CHECK_STATUS(partitioner.GetBestPartition(&blocks));
  // Sort in lexico order only the new blocks.
  std::sort(blocks.begin() + num_default_blocks_to_ignore, blocks.end());

  // Copy from Block struct to CodedBlock class into 'area_remaining_blocks'.
  WP2_CHECK_ALLOC_OK(area_remaining_blocks->resize(
      area_remaining_blocks->size() +
      (blocks.size() - num_default_blocks_to_ignore)));
  for (size_t i = num_default_blocks_to_ignore,
              j = num_already_added_area_blocks;
       i < blocks.size(); ++i, ++j) {
    area_remaining_blocks->at(j).SetDimDefault(blocks[i]);
  }
  return WP2_STATUS_OK;
}

//------------------------------------------------------------------------------

static constexpr uint32_t kMaxCertainty = 3;
static constexpr int32_t kKernelSize = 5, kMaxDist = kKernelSize / 2;

static void GetMinMax(const Plane16& src, Plane16& min, Plane16& max) {
  assert(min.w_ == max.w_ && min.h_ == max.h_);
  const uint32_t step = src.Step();
  const uint32_t num_blocks = min.w_;
  const int16_t* row = (const int16_t*)src.Row(0);
  for (uint32_t y = 0; y < min.h_; ++y) {
    int16_t* const min_row = (int16_t*)min.Row(y);
    int16_t* const max_row = (int16_t*)max.Row(y);
    WP2::GetBlockMinMax(row, step, num_blocks, min_row, max_row);
    row += kMinBlockSizePix * step;
  }
}

// Range of original YUV values in a kKernelSize window around each pixel.
// The higher the value, the more heterogenous the area is.
static void GetSpread(const Plane16& src, Plane16* const dst) {
  const int32_t w = src.w_, h = src.h_;
  int16_t* p_dst = dst->Row(0);
  const uint32_t dst_step = dst->Step();
  const int16_t* p_src = src.Row(0);
  const uint32_t src_step = src.Step();
  for (int32_t y = 0; y < h; ++y) {
    int16_t min, max;
    if (y < kMaxDist || y >= h - kMaxDist) {
      const int32_t min_sub_y = std::max(-kMaxDist, -y);
      const int32_t max_sub_y = std::min(kMaxDist, h - 1 - y);
      const int16_t* const row = &src.At(0, y + min_sub_y);
      for (int32_t x = 0; x < w; ++x) {
        const int32_t min_sub_x = std::max(-kMaxDist, -x);
        const int32_t max_sub_x = std::min(kMaxDist, w - 1 - x);
        GetBlockMinMaxGeneric(row + x + min_sub_x, src_step,
                              max_sub_x + 1 - min_sub_x,
                              max_sub_y + 1 - min_sub_y, &min, &max);
        p_dst[x] = (int16_t)ClampToSigned(max - min, kMaxYuvBits + 1u);
      }
    } else {
      const int16_t* const row = p_src - src_step * kMaxDist;
      for (int32_t x = 0; x < w; ++x) {
        if (x < kMaxDist || x - kMaxDist + 8 > w) {
          const int32_t min_sub_x = std::max(-kMaxDist, -x);
          const int32_t max_sub_x = std::min(kMaxDist, w - 1 - x);
          GetBlockMinMaxGeneric(row + x + min_sub_x, src_step,
                                max_sub_x + 1 - min_sub_x, kKernelSize,
                                &min, &max);
        } else {
          GetBlockMinMax_5x5(row + x - kMaxDist, src_step, &min, &max);
        }
        p_dst[x] = (int16_t)ClampToSigned(max - min, kMaxYuvBits + 1u);
      }
    }
    p_dst += dst_step;
    p_src += src_step;
  }
}

constexpr float MultiScoreFunc::kMinScore;
constexpr int MultiScoreFunc::kMinSpeedForGoodQuantDCT;

WP2Status MultiScoreFunc::Init(const EncoderConfig& config,
                               const Rectangle& tile_rect, const YUVPlane& yuv,
                               const GlobalParams& gparams) {
  DrctFilterInit();
  ScoreDspInit();

  WP2_CHECK_STATUS(PartitionScoreFunc::Init(config, tile_rect, yuv, gparams));
  yuv_range_ =
      (float)(gparams.transf_.GetYUVMax() - gparams.transf_.GetYUVMin());
  a_range_ratio_ = yuv_range_ / kAlphaMax;

  // Cache the 'min_' and 'max_' luma/chroma values per kMinBlockSize square.
  WP2_CHECK_STATUS(min_.Resize(SizeBlocks(src_->Y.w_), SizeBlocks(src_->Y.h_)));
  WP2_CHECK_STATUS(max_.Resize(min_.Y.w_, min_.Y.h_));
  for (Channel channel : {kYChannel, kUChannel, kVChannel}) {
    const Plane16& src_plane = src_->GetChannel(channel);
    assert(src_plane.w_ % kMinBlockSizePix == 0 &&
           src_plane.h_ % kMinBlockSizePix == 0);
    GetMinMax(src_plane, min_.GetChannel(channel), max_.GetChannel(channel));
  }

  if (config_->speed >= kMinSpeedForGoodQuantDCT) {
    // Image processing.
    WP2_CHECK_STATUS(spread_.Resize(src_->GetWidth(), src_->GetHeight(),
                                    /*pad=*/1, src_->HasAlpha()));
    for (Channel channel : {kYChannel, kUChannel, kVChannel, kAChannel}) {
      if (channel == kAChannel && !src_->HasAlpha()) continue;
      GetSpread(src_->GetChannel(channel), &spread_.GetChannel(channel));
    }
  }

  // Per-block standard deviation.
  WP2_CHECK_STATUS(
      stddev_.Allocate(num_block_cols_, num_block_rows_, kMinBlockSizePix));
  stddev_.AddValues(*src_);

  if (src_->HasAlpha()) {
    WP2_CHECK_STATUS(
        a_stddev_.Allocate(num_block_cols_, num_block_rows_, kMinBlockSizePix));
    a_stddev_.AddValues(src_->A);
  }

  // Per-block luma general direction.
  // In might give better results to compute that for each NxN block instead of
  // aggregating pre-computed 4x4 ones but it is probably too expensive.
  WP2_CHECK_ALLOC_OK(direction_.resize(num_block_cols_ * num_block_rows_));
  WP2_CHECK_ALLOC_OK(
      direction_certainty_.resize(num_block_cols_ * num_block_rows_));

  const uint32_t bitdepth = gparams.transf_.GetYUVPrecisionBits() + 1;
  for (uint32_t y = 0; y < num_block_rows_; ++y) {
    const int16_t* x_ptr = &src_->Y.At(/*x=*/0, y * kMinBlockSizePix);
    for (uint32_t i = y * num_block_cols_; i < (y + 1) * num_block_cols_;
         ++i, x_ptr += kMinBlockSizePix) {
      uint32_t variance;
      CdefDirection4x4(x_ptr, src_->Y.Step(), bitdepth, &direction_[i],
                       &variance);
      direction_certainty_[i] = Clamp(variance >> 4, 0u, kMaxCertainty);
      // TODO(yguyon): Also compute, store, use 8x8 direction for bigger blocks
    }
  }

  WP2_CHECK_STATUS(DrawVDebug());
  return WP2_STATUS_OK;
}

//------------------------------------------------------------------------------

WP2Status MultiScoreFunc::GetScore(const Block& block, float* const score) {
  float value = 0.f, threshold = 1.f;  // Passing if 'value <= threshold'.
  switch (pass_) {
    case Pass::LumaAlphaGradient:
      value = GetLumaAlphaGradient(block);
      threshold = GetLumaAlphaGradientThreshold(block);
      break;
    case Pass::NarrowStdDev:
      value = GetStdDevRange(block);
      threshold = GetStdDevRangeThreshold(block);
      break;
    case Pass::GoodQuantDCT:
      value = GetQuantDCT(block);
      threshold = GetQuantDCTThreshold(block);
      break;
    case Pass::Direction:
      value = GetDirection(block);
      threshold = GetDirectionThreshold(block);
      break;
    case Pass::Any:
      value = 0.f;
      threshold = 1.f;
      break;
    default:
      assert(false);
  }

  // Convert to "higher score is better" in [0:1], 0.5 being the threshold.
  if (value <= threshold) {
    *score = 1.001f - 0.5f * value / threshold;  // Will pass.
  } else {
    *score = 0.499f * threshold / value;  // Will not pass.
  }
  RegisterScoreForVDebug(block, *score);
  return WP2_STATUS_OK;
}

//------------------------------------------------------------------------------

void MultiScoreFunc::GetCoeffs(Channel channel, const Block& block,
                               int32_t coeffs[kMaxBlockSizePix2]) const {
  const Plane16& plane = src_->GetChannel(channel);
  // Copy luma pixels.
  int32_t* dst_row = coeffs;
  const int16_t* src_row = &plane.At(block.x_pix(), block.y_pix());
  for (uint32_t y = 0; y < block.h_pix(); ++y) {
    for (uint32_t x = 0; x < block.w_pix(); ++x) {
      dst_row[x] = src_row[x];
    }
    dst_row += block.w_pix();
    src_row += plane.Step();
  }
}

void MultiScoreFunc::QuantizeCoeffs(Channel channel, const Block& block,
                                    int32_t coeffs[kMaxBlockSizePix2]) const {
  GetCoeffs(channel, block, coeffs);

  // Transform them.
  WP2Transform2D(coeffs, kDct, kDct, block.w_pix(), block.h_pix(), coeffs,
                 /*reduced=*/false);

  // Find the segment of the block.
  const Rectangle padded_tile_rect = {tile_rect_.x, tile_rect_.y,
                                      src_->GetWidth(), src_->GetHeight()};
  const uint8_t segment_id =
      AssignSegmentId(*config_, *gparams_, padded_tile_rect, block);
  const QuantMtx& quant_mtx =
      gparams_->segments_[segment_id].GetQuant(channel);

  // Quantize and dequantize coefficients.
  int16_t quantized_coeffs[kMaxBlockSizePix2];
  const TrfSize tdim = GetTransform(block.dim());
  uint32_t num_coeffs;
  quant_mtx.Quantize(coeffs, tdim, /*first_is_dc=*/true, quantized_coeffs,
                     &num_coeffs);
  quant_mtx.Dequantize(quantized_coeffs, num_coeffs, tdim, coeffs);

  // Inverse transform them back and compare the result.
  WP2InvTransform2D(coeffs, kDct, kDct, block.w_pix(), block.h_pix(), coeffs,
                    /*reduced=*/false);
}

void MultiScoreFunc::QuantizeCoeffs(Channel channel, const Block& block,
                                    BlockSize sub_block_size,
                                    int32_t coeffs[kMaxBlockSizePix2],
                                    int32_t* const max_range) const {
  if (max_range != nullptr) *max_range = 0;
  Block sub_block(block.x(), block.y(), sub_block_size);
  while (sub_block.y() + sub_block.h() <= block.y() + block.h()) {
    int32_t sub_coeffs[kMaxBlockSizePix2];
    QuantizeCoeffs(channel, sub_block, sub_coeffs);
    int32_t min = sub_coeffs[0], max = sub_coeffs[0];
    const uint32_t y = sub_block.y_pix() - block.y_pix();
    const uint32_t x = sub_block.x_pix() - block.x_pix();
    for (uint32_t sub_y = 0; sub_y < sub_block.h_pix(); ++sub_y) {
      for (uint32_t sub_x = 0; sub_x < sub_block.w_pix(); ++sub_x) {
        const int32_t sub_coeff = sub_coeffs[sub_y * sub_block.w_pix() + sub_x];
        min = std::min(min, sub_coeff);
        max = std::max(max, sub_coeff);
        coeffs[(y + sub_y) * block.w_pix() + (x + sub_x)] = sub_coeff;
      }
    }
    if (max_range != nullptr) *max_range = std::max(*max_range, max - min);
    sub_block.SetXY(sub_block.x() + sub_block.w(), sub_block.y());
    if (sub_block.x() + sub_block.w() > block.x() + block.w()) {
      sub_block.SetXY(block.x(), sub_block.y() + sub_block.h());
    }
  }
}

//------------------------------------------------------------------------------

// Returns the maximum difference between 'src' and its predicted gradient.
// 'step' goes from a 'src' line to the next.
static int32_t GetGradientDiff(const int16_t* src, int32_t step,
                               int32_t w, int32_t h) {
  assert(w >= 4 && h >= 4);
  const int32_t max_x = w - 1, max_y = h - 1;

  // Average the corners (division by 3 is done at the very end).
  const int32_t top_left = src[0] + src[1] + src[step];
  const int32_t bottom_left =
      src[(max_y - 1) * step] + src[max_y * step] + src[max_y * step + 1];
  const int32_t top_right = src[max_x - 1] + src[max_x] + src[step + max_x];
  const int32_t bottom_right = src[(max_y - 1) * step + max_x] +
                               src[max_y * step + max_x - 1] +
                               src[max_y * step + max_x];

  // Create a gradient by bidimensional interpolation and compare with 'src'.
  int32_t max_diff = 0;
  for (int32_t y = 0; y <= max_y; ++y) {
    const int32_t left = top_left * (max_y - y) + bottom_left * y;
    const int32_t right = top_right * (max_y - y) + bottom_right * y;
    for (int32_t x = 0; x <= max_x; ++x) {
      const int32_t gradient_pixel =
          DivRound(left * (max_x - x) + right * x, 3 * max_x * max_y);
      max_diff = std::max(max_diff, std::abs(src[x] - gradient_pixel));
    }
    src += step;
  }
  return max_diff;
}

float MultiScoreFunc::GetLumaAlphaGradient(const Block& block) const {
  // Empirically chosen values.
  const float kDiffScale[] = {1.f, 1.f, 1.f, 0.25f / kAlphaMax * yuv_range_};
  float max_diff = 0.f;
  // kUChannel and kVChannel do not bring valuable partition decision-making
  // here so skip them for speed.
  for (Channel c : {kYChannel, kAChannel}) {
    if (c == kAChannel && !src_->HasAlpha()) continue;
    const int32_t diff = GetGradientDiff(
        &src_->GetChannel(c).At(block.x_pix(), block.y_pix()),
        src_->GetChannel(c).Step(), block.w_pix(), block.h_pix());
    max_diff = std::max(max_diff, diff * kDiffScale[c]);
  }
  return max_diff;
}

float MultiScoreFunc::GetLumaAlphaGradientThreshold(const Block& block) const {
  // The threshold is tighter for bigger blocks at higher qualities.
  // Medium blocks are ignored except at high qualities.
  constexpr BaseSlope kBaseSlope[] = {
      {0.f, 0.f},        //  4x4, unused
      {0.f, 0.f},        //  8x4, unused
      {-0.6f, 0.7f},     //  8x8  16x4
      {-0.06f, 0.08f},   // 16x8
      {-0.03f, 0.045f},  // 16x16 32x8
      {0.04f, -0.035f},  // 32x16
      {0.05f, -0.045f},  // 32x32
  };
  STATIC_ASSERT_ARRAY_SIZE(kBaseSlope, WP2Log2Ceil_k(kMaxBlockSize2) + 1);
  const uint32_t index = (uint32_t)WP2Log2Floor(block.rect().GetArea());
  assert((1u << index) == block.rect().GetArea());
  return yuv_range_ *
         MapQuality(*config_, kBaseSlope[index].base, kBaseSlope[index].slope);
}

//------------------------------------------------------------------------------

static float StdDevRange(const Block& block, const Integral& variance) {
  // Consider the standard deviation of the whole block as sub-blocks could
  // be coherent within themselves but not with other sub-blocks.
  const uint8_t overall_variance = variance.StdDevUint8(
      block.x(), block.y(), block.x() + block.w(), block.y() + block.h());
  uint8_t min = overall_variance, max = overall_variance;

  // Now check the sub-blocks.
  for (uint32_t sub_y = block.y(); sub_y < block.y() + block.h(); ++sub_y) {
    for (uint32_t sub_x = block.x(); sub_x < block.x() + block.w(); ++sub_x) {
      const uint8_t variance_tmp =
          variance.StdDevUint8(sub_x, sub_y, sub_x + 1, sub_y + 1);
      if (variance_tmp < min) {
        min = variance_tmp;
      } else if (variance_tmp > max) {
        max = variance_tmp;
      }
    }
  }
  return (max - min) / 255.f;
}

float MultiScoreFunc::GetStdDevRange(const Block& block) const {
  float range = StdDevRange(block, stddev_);
  if (!a_stddev_.empty()) {
    const float a_range =
        StdDevRange(block, a_stddev_) * a_range_ratio_ * a_range_ratio_;
    range = std::max(range, a_range);
  }

  return range;
}

float MultiScoreFunc::GetStdDevRangeThreshold(const Block& block) const {
  // The higher the quality, the narrower the standard deviation range needs to
  // be for a block to be accepted. Accept bigger blocks during partitioning at
  // low qualities, and seek smaller blocks at high qualities.
  return MapQuality(*config_, 0.50f, -0.38f);
}

//------------------------------------------------------------------------------

// Returns the average difference between the original luma coefficients and the
// quantized ones, weighted per pixel by the 'spread_' in order to give more
// importance to flat areas (penalize distant ripples more than noise on edges).
float MultiScoreFunc::GetQuantDCT(const Block& block, Channel channel) const {
  int32_t coeffs[kMaxBlockSizePix2];
  QuantizeCoeffs(channel, block, coeffs);

  float avg_diff = 0.f;
  int32_t* dst_row = coeffs;
  const int16_t* src_row =
      &src_->GetChannel(channel).At(block.x_pix(), block.y_pix());
  const int16_t* spread_row =
      &spread_.GetChannel(channel).At(block.x_pix(), block.y_pix());
  for (uint32_t y = 0; y < block.h_pix(); ++y) {
    for (uint32_t x = 0; x < block.w_pix(); ++x) {
      const int32_t diff = std::abs(dst_row[x] - src_row[x]);
      avg_diff += diff / Clamp(spread_row[x] / 20.f, 0.1f, 10.f);
    }
    dst_row += block.w_pix();
    src_row += src_->GetChannel(channel).Step();
    spread_row += spread_.GetChannel(channel).Step();
  }
  avg_diff /= (block.w_pix() * block.h_pix());
  return avg_diff;
}

float MultiScoreFunc::GetQuantDCT(const Block& block) const {
  // Take chroma and alpha into account just enough to discard bad layouts (for
  // example an entirely black image with alpha patterns).
  // TODO(yguyon): Check ADST too
  constexpr float kScale[] = {1.f, 0.1f, 0.1f, 0.1f};
  return std::max({GetQuantDCT(block, kYChannel) * kScale[kYChannel],
                   GetQuantDCT(block, kUChannel) * kScale[kUChannel],
                   GetQuantDCT(block, kVChannel) * kScale[kVChannel],
                   src_->HasAlpha()
                       ? GetQuantDCT(block, kAChannel) * kScale[kAChannel]
                       : 0.f});
}

float MultiScoreFunc::GetQuantDCTThreshold(const Block& block) const {
  // This metric is only useful at low qualities, for large blocks that would
  // still look fine with an aggressive quantization.
  return MapQuality(*config_, 4.00f, -4.75f);
}

//------------------------------------------------------------------------------

// Returns a low score for a 'block' that appears to have an obvious and uniform
// orientation.
float MultiScoreFunc::GetDirection(const Block& block) const {
  const uint32_t stride = num_block_cols_;
  const uint32_t* direction =
      direction_.data() + block.y() * stride + block.x();
  const uint32_t* certainty =
      direction_certainty_.data() + block.y() * stride + block.x();

  uint32_t weight[kDrctFltNumDirs] = {0};
  for (uint32_t sub_y = 0; sub_y < block.h(); ++sub_y) {
    for (uint32_t sub_x = 0; sub_x < block.w(); ++sub_x) {
      assert(direction[sub_x] < kDrctFltNumDirs);
      weight[direction[sub_x]] += std::min(certainty[sub_x], kMaxCertainty);
    }
    direction += stride;
    certainty += stride;
  }
  const uint32_t heaviest =
      std::max_element(weight, weight + kDrctFltNumDirs) - weight;
  const uint32_t max_weight = block.w() * block.h() * kMaxCertainty;
  assert(weight[heaviest] <= max_weight);
  const uint32_t previous_direction =
      (heaviest - 1 + kDrctFltNumDirs) % kDrctFltNumDirs;
  const uint32_t next_direction = (heaviest + 1) % kDrctFltNumDirs;
  const uint32_t weight_with_close_directions =
      weight[heaviest] +
      (weight[previous_direction] + weight[next_direction]) / 4;
  const float direction_score =
      1.f - Clamp(weight_with_close_directions / (float)max_weight, 0.f, 1.f);

  return direction_score;
}

float MultiScoreFunc::GetDirectionThreshold(const Block& block) const {
  return MapQuality(*config_, 0.15f, -0.0475f);
}

//------------------------------------------------------------------------------

WP2Status TileScoreFunc::Init(const EncoderConfig& config,
                              const Rectangle& tile_rect, const YUVPlane& yuv,
                              const GlobalParams& gparams) {
  WP2EncDspInit();
  WP2_CHECK_STATUS(PartitionScoreFunc::Init(config, tile_rect, yuv, gparams));
  local_gparams_.features_ = &features_map_;
  WP2_CHECK_STATUS(GlobalAnalysis(ArgbBuffer(), yuv, gparams.transf_,
                                  config, &local_gparams_));
  WP2_CHECK_STATUS(InitForEncode());

  // Initialize the best score with the partition containing only the forced
  // blocks.
  const Rectangle padded_tile_rect = {tile_rect_.x, tile_rect_.y,
                                      yuv.GetWidth(), yuv.GetHeight()};
  WP2_CHECK_STATUS(AddForcedBlocks(config, padded_tile_rect, &blocks_));
  WP2_CHECK_STATUS(TryEncode(blocks_, &best_score_));
  RegisterScoreForVDebug("starting", {}, best_score_);
  cached_best_score_ = 0.f;
  blocks_.clear();
  return WP2_STATUS_OK;
}

WP2Status TileScoreFunc::GetScore(const Block& block, float* const score) {
  WP2_CHECK_ALLOC_OK(blocks_.push_back(block));
  WP2_CHECK_STATUS(TryEncode(blocks_, score));
  if (*score > cached_best_score_) cached_best_score_ = *score;
  if (*score > best_score_) RegisterScoreForVDebug("new best", block, *score);
  blocks_.pop_back();
  return WP2_STATUS_OK;
}

WP2Status TileScoreFunc::InitForEncode() {
  WP2_CHECK_ALLOC_OK(blocks_.reserve((tile_rect_.width / kMaxBlockSizePix) *
                                     (tile_rect_.height / kMaxBlockSizePix)));

  enc_tiles_layout_.num_tiles_x = enc_tiles_layout_.num_tiles_y = 1;
  enc_tiles_layout_.tile_width = tile_rect_.width;
  enc_tiles_layout_.tile_height = tile_rect_.height;
  WP2_CHECK_ALLOC_OK(enc_tiles_layout_.tiles.resize(1));
  enc_tiles_layout_.first_unassigned_tile_index = 0;

  enc_tiles_layout_.tiles.front().rect = tile_rect_;
  enc_tiles_layout_.tiles.front().rgb_input.Deallocate();  // This is lossy.
  assert(!src_->IsEmpty());
  WP2_CHECK_STATUS(enc_tiles_layout_.tiles.front().yuv_input.SetView(*src_));
  tmp_config_ = *config_;
  tmp_config_.partition_method = sub_partition_method_;
  tmp_config_.info = nullptr;
  tile_encoder_.config_ = &tmp_config_;
  tile_encoder_.use_lossless_ = (tmp_config_.quality > kMaxLossyQuality);
  tile_encoder_.tiles_layout_ = &enc_tiles_layout_;
  WP2_CHECK_STATUS(tile_encoder_.AssignNextTile());

  // Recursion is too dangerous here. It's potentially creating
  // (kMaxTileSize/kMinBlockSizePix)^2 = a lot of recursive encoding contexts.
  assert(sub_partition_method_ != AUTO_PARTITIONING &&
         sub_partition_method_ != TILE_ENCODE_PARTITIONING);

  dec_config_.thread_level = 0;

  WP2_CHECK_STATUS(decompressed_yuv_.Copy(*src_, /*resize_if_needed=*/true));
  // Needed for API compliance. The pixels will not be accessed.
  WP2_CHECK_STATUS(
      decompressed_argb_.Resize(tile_rect_.width, tile_rect_.height));

  // A BitstreamFeatures instance is needed by LossyDecode(). Make up one.
  MemoryWriter writer;
  WP2_CHECK_STATUS(
      EncodeHeader(tmp_config_, tile_rect_.width, tile_rect_.height,
                   src_->HasAlpha(), /*is_anim=*/false, kInfiniteLoop,
                   /*background_color=*/{}, /*preview_color=*/{},
                   /*has_icc=*/false, /*has_trailing_data=*/false, &writer));
  WP2_CHECK_STATUS(features_.Read(writer.mem_, writer.size_));
  return WP2_STATUS_OK;
}

WP2Status TileScoreFunc::TryEncode(const VectorNoCtor<Block>& blocks,
                                   float* const score) {
  ANSEnc& enc = enc_tiles_layout_.tiles.front().enc;
  enc.WipeOut();
  enc_tiles_layout_.gparams = &local_gparams_;
  // Encode the whole tile with the forced 'blocks'.
  WP2_CHECK_STATUS(tile_encoder_.LossyEncode(blocks, &enc));
  WP2_CHECK_STATUS(enc.Assemble());

  // Reset the unique tile to a fresh state.
  const uint32_t width = decompressed_argb_.width;
  const uint32_t height = decompressed_argb_.height;
  const uint32_t tile_width = TileWidth(FinalTileShape(*config_), width);
  const uint32_t tile_height =
      TileHeight(FinalTileShape(*config_), /*image_width=*/width);
  WP2_CHECK_STATUS(GetTilesLayout(width, height, tile_width, tile_height,
                                  &decompressed_argb_, &decompressed_yuv_,
                                  &tiles_layout_));
  assert(tiles_layout_.tiles.size() == 1 &&
         enc_tiles_layout_.tiles.size() == 1);

  // Plug ANSEnc output to ANSDec input.
  Tile* const tile = &tiles_layout_.tiles.front();
  tile->chunk_size_is_known = true;
  tile->chunk_size = enc.BufferSize();
  tiles_layout_.gparams = &local_gparams_;
  tile->private_input = ExternalDataSource(enc.Buffer(), enc.BufferSize());
  tile->input = &tile->private_input;

  // Decode to 'decompressed_argb_'.
  ANSDec dec(tile->input);
  WP2_CHECK_STATUS(LossyDecode(features_, dec_config_, &tiles_layout_,
                               /*progress=*/nullptr, &dec, tile));

  // Compare the pixels of the non-padded area only.
  YUVPlane original_view, decompressed_view;
  WP2_CHECK_STATUS(original_view.SetView(*src_, {0, 0, width, height}));
  WP2_CHECK_STATUS(
      decompressed_view.SetView(decompressed_yuv_, {0, 0, width, height}));
  WP2_CHECK_STATUS(decompressed_view.GetDistortion(
      original_view, kMaxYuvBits + 1, PSNR, distortion_));

  // Compute a score based on distortion and the number of bits per pixel.
  const float ssim = distortion_[4];
  const float bpp =
      std::max(1u, (uint32_t)enc.BufferSize()) * 8.f / (width * height);
  const float lambda = MapQuality(*config_, 7.00f, -2.85f);
  *score = ssim - lambda * bpp;
  return WP2_STATUS_OK;
}

WP2Status TileScoreFunc::Use(const Block& block) {
  WP2_CHECK_ALLOC_OK(blocks_.push_back(block));
  if (cached_best_score_ > best_score_) best_score_ = cached_best_score_;
  cached_best_score_ = 0.f;
  return WP2_STATUS_OK;
}

//------------------------------------------------------------------------------

WP2Status FixedSizeScoreFunc::GetScore(const Block& block, float* const score) {
  *score = (block.dim() == size_) ? 1.f : 0.f;
  return WP2_STATUS_OK;
}

//------------------------------------------------------------------------------

}  // namespace WP2
