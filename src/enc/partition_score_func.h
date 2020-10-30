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

#ifndef WP2_ENC_PARTITION_SCORE_FUNC_H_
#define WP2_ENC_PARTITION_SCORE_FUNC_H_

#include "src/common/integral.h"
#include "src/common/lossy/block.h"
#include "src/common/lossy/block_size.h"
#include "src/common/lossy/block_size_io.h"
#include "src/common/lossy/predictor.h"
#include "src/common/symbols.h"
#include "src/dec/wp2_dec_i.h"
#include "src/enc/symbols_enc.h"
#include "src/enc/wp2_enc_i.h"
#include "src/utils/csp.h"
#include "src/utils/plane.h"
#include "src/utils/vector.h"
#include "src/wp2/encode.h"

namespace WP2 {

//------------------------------------------------------------------------------

// Used by Partitioner to evaluate a block layout in a tile.
class PartitionScoreFunc : public WP2Allocable {
 public:
  virtual ~PartitionScoreFunc() = default;

  // Configures the partition settings and allocates memory.
  // 'tile_rect' is in pixels, not padded.
  virtual WP2Status Init(const EncoderConfig& config,
                         const Rectangle& tile_rect, const YUVPlane& yuv,
                         const GlobalParams& gparams);

  // For a given block position and size, returns a score (higher means better).
  virtual WP2Status GetScore(const Block& block, float* const score) = 0;

  // Called every time a block is definitively chosen.
  virtual WP2Status Use(const Block& block);

 protected:
  const EncoderConfig* config_ = nullptr;
  const GlobalParams* gparams_ = nullptr;
  Rectangle tile_rect_;            // In pixels, not padded.
  const YUVPlane* src_ = nullptr;  // Original image (tile).
  uint32_t num_block_cols_ = 0;  // Maximum number of blocks on the horizontal
  uint32_t num_block_rows_ = 0;  // and vertical axis.

  // Visual debug.
  WP2Status ClearVDebug() const;
  bool RegisterScoreForVDebug(const Block& block, float score,
                              bool ending_new_line = true) const;
};

//------------------------------------------------------------------------------

// Based on block encoding. Can only be used with a TopLeftBestPartitioner.
class BlockScoreFunc : public PartitionScoreFunc {
 public:
  WP2Status Init(const EncoderConfig& config, const Rectangle& tile_rect,
                 const YUVPlane& yuv, const GlobalParams& gparams) override;

  WP2Status GetScore(const Block& block, float* const score) override;

  // Returns the score for encoding several 'blocks' (at most 4).
  // Can be used for evaluating a block against its sub-blocks, for example.
  WP2Status GetScore(const Block blocks[4], uint32_t num_blocks,
                     float* const score);

  WP2Status Use(const Block& block) override;

 protected:
  // Tries some or all transform/predictor pairs and keeps the best one in 'cb'.
  WP2Status FindBestBlockParams(const FrontMgrNxNBase& front_mgr,
                                const BlockContext& block_context,
                                SyntaxWriter* const syntax_writer,
                                DCDiffusionMap* const dc_error_u,
                                DCDiffusionMap* const dc_error_v,
                                CodedBlock* const cb) const;
  // Finds the best encoding parameters then encodes and records 'cb'.
  // Advances the 'front_mgr'.
  WP2Status EncodeBlock(const FrontMgrNxNBase& front_mgr, CodedBlock* const cb,
                        SyntaxWriter* const syntax_writer,
                        DCDiffusionMap* const dc_error_u,
                        DCDiffusionMap* const dc_error_v,
                        YUVPlane* const buffer) const;
  // Writes 'cb' to 'enc'. Advances the 'front_mgr'.
  WP2Status WriteBlock(const FrontMgrNxNBase& front_mgr, const CodedBlock& cb,
                       SyntaxWriter* const syntax_writer,
                       ANSEnc* const enc) const;

  // Defined by 'config_' and not modified outside Init().
  Vector<TransformPair> transforms_, transforms_subset_;

  // Instances for recording and estimating the bit-rate/disto per block.
  mutable YUVPlane buffer_;  // Reconstructed blocks (final + under test).
  FrontMgrCustom front_mgr_;
  ANSDictionaries dicts_;
  SyntaxWriter syntax_writer_;
  DCDiffusionMap dc_error_u_, dc_error_v_;

  // Temporary instances.
  CodedBlock tmp_cb_;
  ANSDictionaries tmp_dicts_;
  SyntaxWriter tmp_syntax_writer_;
  DCDiffusionMap tmp_dc_error_u_, tmp_dc_error_v_;
  BlockContext context_;

 private:
  bool RegisterScoreForVDebug(const Block blocks[], uint32_t num_blocks,
                              const float rate[4], const float disto[4],
                              float total_rate, float total_disto, float score);
};

//------------------------------------------------------------------------------

// Encodes all the blocks in an area to estimate a rate-distortion score.
// Can only be used with a AreaRDOptPartitioner.
class AreaScoreFunc : public BlockScoreFunc {
 public:
  // Dimensions in pixels of the areas.
  AreaScoreFunc(uint32_t area_width, uint32_t area_height);
  WP2Status Init(const EncoderConfig& config, const Rectangle& tile_rect,
                 const YUVPlane& yuv, const GlobalParams& gparams) override;
  const Rectangle& GetCurrentArea() const { return area_; }

  // Unused.
  WP2Status GetScore(const Block& block, float* const score) override;

  // Returns the 'blocks' and 'score' of encoding the 'area_' with the default
  // partitioning (multipass).
  WP2Status GetAreaDefaultScore(VectorNoCtor<Block>* const blocks,
                                float* const score);
  // Returns the 'score' of encoding the 'area_' with all 'blocks' having the
  // same 'block_size'.
  WP2Status GetAreaGridScore(BlockSize block_size,
                             VectorNoCtor<Block>* const blocks,
                             float* const score);

  // Sets a block as final. After choosing the partition of the current 'area_',
  // call this for each block in it.
  WP2Status Use(const Block& block) override;

 protected:
  // Returns the distortion and the bitstream size of encoding 'area_blocks'.
  WP2Status GetDistoRate(Vector<CodedBlock>* const area_blocks,
                         float* const disto, float* const rate) const;

  // Setup a new zone to focus on for partitioning.
  // 'area_x' and 'area_y' are coordinates in pixels within the tile.
  virtual WP2Status BeginArea(uint32_t area_x, uint32_t area_y);

  // Returns the 'area_blocks' as if they were selected by MultiScoreFunc.
  WP2Status GetAreaDefaultPartition(
      Vector<CodedBlock>* const area_blocks) const;

  // Defined by 'config_' and not modified outside Init().
  VectorNoCtor<Block> default_partition_;  // Generated by MultiScoreFunc.

  // Objects for the current zone.
  const uint32_t area_width_, area_height_;  // In pixels. Can be less on edges.
  Rectangle area_;  // Boundaries in pixels (no pad) of the region to consider
                    // for the RD-opt. Only use the surroundings as context for
                    // prediction. Next area begins when this area is complete.
  FrontMgrArea area_front_mgr_;    // Used to browse blocks area by area.
  const FrontMgrArea::Comp comp_;  // Used to group blocks per area.

  // Cached default partitioning stats of the current 'area_' for comparison.
  float default_disto_ = -1.f, default_rate_ = -1.f;

 private:
  // Visual debug.
  bool RegisterScoreForVDebug(BlockSize grid_size,
                              const Vector<CodedBlock>& area_blocks,
                              float score, float disto, float rate) const;
};

//------------------------------------------------------------------------------

// Similar to AreaScoreFunc but the default partitioning is computed and
// compared for each block inside every area instead of once per area.
// Can only be used with a SubAreaRDOptPartitioner.
class SubAreaScoreFunc : public AreaScoreFunc {
 public:
  SubAreaScoreFunc(uint32_t area_width, uint32_t area_height)
      : AreaScoreFunc(area_width, area_height) {}

  WP2Status GetScore(const Block& block, float* const score) override;
  WP2Status Use(const Block& block) override;

 protected:
  WP2Status BeginArea(uint32_t area_x, uint32_t area_y) override;

  WP2Status GetAreaRemainingDefaultPartition(
      Vector<CodedBlock>* const area_remaining_blocks) const;

  Vector<CodedBlock> area_used_blocks_;  // Final encoded blocks in this 'area_'

  // Cached default partitioning stats of the current block for comparison.
  Block default_block_;  // BLK_LAST if undefined.
  float default_block_disto_ = 0.f, default_block_rate_ = 0.f;

 private:
  // Visual debug.
  bool RegisterScoreForVDebug(const Block& block,
                              const Vector<CodedBlock>& area_remaining_blocks,
                              float score, float disto, float rate) const;
};

//------------------------------------------------------------------------------

// Uses several metrics to select block sizes.
class MultiScoreFunc : public PartitionScoreFunc {
 public:
  static constexpr float kMinScore = 0.5f;  // Blocks below that are discarded.
  static constexpr int kMinSpeedForGoodQuantDCT = 4;

  enum class Pass {  // The following passes select blocks with:
    LumaAlphaGradient,  // gradient-like original luma and alpha values
    NarrowStdDev,  // a narrow range of orig luma/alpha std dev per 4x4 block
    GoodQuantDCT,  // low loss when luma->DCT->quantized->dequantized->luma
    Direction,         // a common direction with a high enough certainty
    Any,               // no criterion (used for filling what's left)
  };

  WP2Status Init(const EncoderConfig& config, const Rectangle& tile_rect,
                 const YUVPlane& yuv, const GlobalParams& gparams) override;

  // Sets the 'pass_' of the next GetScore().
  void SetPass(Pass pass) { pass_ = pass; }
  // Returns a value/threshold ratio of the current 'pass_' in ]0:1].
  // Higher means better, scores above 0.5 are selected.
  WP2Status GetScore(const Block& block, float* const score) override;

 protected:
  // Helper functions.
  void GetCoeffs(Channel channel, const Block& block,
                 int32_t coeffs[kMaxBlockSizePix2]) const;
  void QuantizeCoeffs(Channel channel, const Block& block,
                      int32_t coeffs[kMaxBlockSizePix2]) const;
  void QuantizeCoeffs(Channel channel, const Block& block,
                      BlockSize sub_block_size,
                      int32_t coeffs[kMaxBlockSizePix2],
                      int32_t* const max_range = nullptr) const;
  float GetQuantDCT(const Block& block, Channel channel) const;

  // Passes. Values below threshold are selected.
  float GetLumaAlphaGradient(const Block& block) const;
  float GetLumaAlphaGradientThreshold(const Block& block) const;

  float GetStdDevRange(const Block& block) const;
  float GetStdDevRangeThreshold(const Block& block) const;

  float GetQuantDCT(const Block& block) const;
  float GetQuantDCTThreshold(const Block& block) const;

  float GetDirection(const Block& block) const;
  float GetDirectionThreshold(const Block& block) const;

  float yuv_range_ = 0.f;  // Range of allowed original luma values.
  float a_range_ratio_;    // Range of yuv values over range of alpha values.

  // Results of image processing done on the original planes.
  YUVPlane min_, max_;  // Minimum and maximum original values per block.
  YUVPlane spread_;     // Range of surrounding values.
  Integral stddev_;     // Contains quantized luma std dev of the block grid.
  Integral a_stddev_;   // Contains std dev for the alpha plane (if any).

  // Pixel orientation of the block grid. TODO(yguyon): Use less memory
  Vector_u32 direction_;            // [0:7] maps from +45 to -112.5 degrees.
  Vector_u32 direction_certainty_;  // [0:3] the higher the more confident.

  Pass pass_ = Pass::Any;

  // Visual debug.
  void DrawValueThresholdColor(uint32_t cell_width, uint32_t cell_height,
                               const Block& block,
                               ArgbBuffer* const debug_output) const;
  WP2Status DrawLossyLuma(const Block& block, const Rectangle& block_rect,
                          bool draw_side, ArgbBuffer* const debug_output) const;
  WP2Status DrawSelection(const Block& block, const Rectangle& block_rect,
                          ArgbBuffer* const debug_output) const;
  WP2Status DrawVDebug() const;
};

//------------------------------------------------------------------------------

// Based on tile encoding then decoding. Warning: very slow.
class TileScoreFunc : public PartitionScoreFunc {
 public:
  explicit TileScoreFunc(
      PartitionMethod sub_partition_method = MULTIPASS_PARTITIONING)
      : sub_partition_method_(sub_partition_method) {}

  WP2Status Init(const EncoderConfig& config, const Rectangle& tile_rect,
                 const YUVPlane& yuv, const GlobalParams& gparams) override;

  // Encodes the whole tile with the forced 'blocks' and returns a score based
  // on its distortion and size in bytes.
  WP2Status TryEncode(const VectorNoCtor<Block>& blocks, float* const score);

  // Temporarily forces the 'block' for a TryEncode() call.
  WP2Status GetScore(const Block& block, float* const score) override;

  // Forces the 'block' for all future TryEncode() calls.
  WP2Status Use(const Block& block) override;

 protected:
  WP2Status InitForEncode();

  const PartitionMethod sub_partition_method_;
  EncoderConfig tmp_config_;
  VectorNoCtor<Block> blocks_;
  EncTilesLayout enc_tiles_layout_;
  TileEncoder tile_encoder_;
  GlobalParams local_gparams_;
  ArgbBuffer argb_;
  FeatureMap features_map_;
  float best_score_ = 0.f;         // Best score before last call to Use().
  float cached_best_score_ = 0.f;  // Best score since last call to Use().

  YUVPlane decompressed_yuv_;
  ArgbBuffer decompressed_argb_;
  BitstreamFeatures features_;
  DecoderConfig dec_config_;
  TilesLayout tiles_layout_;
  float distortion_[5];

  // Visual debug.
  bool RegisterScoreForVDebug(const char label[], const Block& block,
                              float score) const;
};

//------------------------------------------------------------------------------

// Gives a high score to a predetermined block size and 0 to everything else.
class FixedSizeScoreFunc : public PartitionScoreFunc {
 public:
  explicit FixedSizeScoreFunc(BlockSize size) : size_(size) {}

  WP2Status GetScore(const Block& block, float* const score) override;

 protected:
  BlockSize size_;
};

//------------------------------------------------------------------------------

}  // namespace WP2

#endif  // WP2_ENC_PARTITION_SCORE_FUNC_H_
