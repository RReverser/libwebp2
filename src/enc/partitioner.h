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
// Tool for finding the best block layout.
//
// Author: Yannis Guyon (yguyon@google.com)

#ifndef WP2_ENC_PARTITIONER_H_
#define WP2_ENC_PARTITIONER_H_

#include "src/common/lossy/block.h"
#include "src/common/lossy/block_size.h"
#include "src/common/lossy/predictor.h"
#include "src/common/symbols.h"
#include "src/enc/partition_score_func.h"
#include "src/enc/symbols_enc.h"
#include "src/utils/csp.h"
#include "src/utils/plane.h"
#include "src/utils/vector.h"
#include "src/wp2/encode.h"

namespace WP2 {

//------------------------------------------------------------------------------

// Verifies and adds the blocks concerning this 'tile' from
// 'EncoderInfo::force_partition' to 'out'.
WP2Status AddForcedBlocks(const EncoderConfig& config,
                          const Rectangle& padded_tile_rect,
                          VectorNoCtor<Block>* const out);

//------------------------------------------------------------------------------

// Used to find the best block layout in a tile.
class Partitioner : public WP2Allocable {
 public:
  virtual ~Partitioner() = default;

  // Configures the partition settings and allocates memory.
  virtual WP2Status Init(const EncoderConfig& config, const YUVPlane& yuv,
                         const Rectangle& tile_rect,
                         PartitionScoreFunc* const score_func);

  // Computes and returns the best layout. Input 'blocks' are considered forced.
  virtual WP2Status GetBestPartition(VectorNoCtor<Block>* const blocks) = 0;

 protected:
  // Gets or sets areas as occupied.
  bool IsOccupied(const Block& block) const;
  void Occupy(const Block& block);

  // Returns true if the 'block' fits within 'max_block', is snapped etc.
  bool IsBlockValid(const Block& block, const Block& max_block,
                    size_t num_forced_blocks) const;

  // Verifies that the 'forced_blocks' don't overlap each other.
  // 'max_num_blocks_left' will be decremented by the sum of their areas.
  WP2Status RegisterForcedBlocks(const VectorNoCtor<Block>& forced_blocks,
                                 uint32_t* const max_num_blocks_left);

  const EncoderConfig* config_ = nullptr;
  Rectangle tile_rect_;            // In pixels, not padded.
  const YUVPlane* src_ = nullptr;  // Original image, padded.
  uint32_t num_block_cols_ = 0;    // Maximum number of blocks on the horizontal
  uint32_t num_block_rows_ = 0;    // and vertical axis.
  PartitionScoreFunc* score_func_ = nullptr;

  VectorNoCtor<bool> occupancy_;  // For a given pos, true if occupied.

  // Visual debug: display the order in which the 'block' was chosen and during
  // which pass.
  void RegisterOrderForVDebug(uint32_t pass, uint32_t pass_index,
                              const Block& block, uint32_t block_index,
                              uint32_t max_num_blocks) const;
};

//------------------------------------------------------------------------------

// For each block (starting with the top-left one), takes the best size.
class TopLeftBestPartitioner : public Partitioner {
 public:
  WP2Status GetBestPartition(VectorNoCtor<Block>* const blocks) override;

 private:
  // Returns the 'best_block' size within 'max_block'.
  WP2Status GetBestSize(const Block& max_block, size_t num_forced_blocks,
                        Block* const best_block);
};

//------------------------------------------------------------------------------

// Chooses the best sub-partition in each area.
// To be used only with AreaScoreFunc.
class AreaRDOptPartitioner : public Partitioner {
 public:
  AreaRDOptPartitioner(uint32_t area_width, uint32_t area_height)  // In pixels.
      : area_width_(area_width), area_height_(area_height) {}
  WP2Status GetBestPartition(VectorNoCtor<Block>* const blocks) override;

 private:
  virtual WP2Status GetAreaBestPartition(const Rectangle& area,
                                         VectorNoCtor<Block>* const blocks,
                                         uint32_t* const max_num_blocks_left);

  const uint32_t area_width_, area_height_;  // In pixels.
};

//------------------------------------------------------------------------------

// Same as AreaRDOptPartitioner but encodes the whole area for each size of each
// block inside it. To be used only with SubAreaScoreFunc.
class SubAreaRDOptPartitioner : public AreaRDOptPartitioner {
 public:
  SubAreaRDOptPartitioner(uint32_t area_width, uint32_t area_height)
      : AreaRDOptPartitioner(area_width, area_height),
        front_mgr_(area_width, area_height) {}
  WP2Status GetBestPartition(VectorNoCtor<Block>* const blocks) override;

 private:
  WP2Status GetAreaBestPartition(const Rectangle& area,
                                 VectorNoCtor<Block>* const blocks,
                                 uint32_t* const max_num_blocks_left) override;
  FrontMgrArea front_mgr_;
};

//------------------------------------------------------------------------------

// For each block, tries some split patterns and keeps the best one then
// recurses on sub-blocks.
class SplitRecursePartitioner : public Partitioner {
 public:
  WP2Status Init(const EncoderConfig& config, const YUVPlane& yuv,
                 const Rectangle& tile_rect,
                 PartitionScoreFunc* const score_func) override;
  WP2Status GetBestPartition(VectorNoCtor<Block>* const blocks) override;

 private:
  // Returns the maximum BlockSize fitting in 'max_block' and respecting all the
  // constraints (not overlapping a forced block or a previously-decided split
  // pattern, belonging to the partition set etc.).
  BlockSize GetMaxFittingBlockSize(const Block& max_block,
                                   uint32_t num_forced_blocks) const;

  // Contains one way of splitting a block into smaller ones.
  struct SplitPattern {
    uint32_t num_blocks;       // Number of blocks in this split.
    BlockSize block_sizes[4];  // Top-left,top-right,bottom-left,bottom-right
                               // (or subset in same order if less than 4).
  };
  // Returns all possible SplitPatterns given a 'set' and a 'size'.
  static void GetSplitPatterns(PartitionSet set, BlockSize size,
                               Vector<SplitPattern>* const p);
  // Returns the 'blocks' given a 'max_block' (size and position) and a
  // 'split_pattern'.
  static void GetSplitPatternBlocks(const Block& max_block,
                                    const SplitPattern& split_pattern,
                                    Block blocks[4]);

  // Contains whether a SplitPattern was decided for a given block position.
  struct SplitInfo {
    BlockSize size;  // BLK_LAST if nothing decided yet.
    bool recursive;  // False if the decision is definitive
                     // or true if smaller blocks may be tried.
    bool operator==(const SplitInfo& other) const {
      return size == other.size && recursive == other.recursive;
    }
    bool operator!=(const SplitInfo& other) const { return !operator==(other); }
  };
  static constexpr SplitInfo kNoSplitInfo = {BLK_LAST, false};  // No decision.

  // Gets the split pattern decision for a given block position.
  const SplitInfo& GetSplitInfo(uint32_t x, uint32_t y) const;
  // Returns true if the split pattern decision is 'split_info' in all the given
  // 'block' area.
  bool IsSplitInfo(const Block& block, const SplitInfo& split_info) const;
  // Sets the split pattern decision for a given 'block' area.
  void SetSplitInfo(const Block& block, bool recursive);

  // Contains whether a SplitPattern was decided or not for each block position.
  Plane<SplitInfo> split_info_;
};

//------------------------------------------------------------------------------

// Tries all possible block layouts. Warning: super slow.
class ExhaustivePartitioner : public Partitioner {
 public:
  // This partitioner needs a specific PartitionScoreFunc.
  explicit ExhaustivePartitioner(TileScoreFunc* const score_func)
      : tile_score_func_(score_func) {}

  WP2Status Init(const EncoderConfig& config, const YUVPlane& yuv,
                 const Rectangle& tile_rect,
                 PartitionScoreFunc* const score_func) override;

  WP2Status GetBestPartition(VectorNoCtor<Block>* const blocks) override;

 private:
  // Finds the next valid partition. If 'found', returns its 'score'.
  WP2Status GetNextPartitionScore(bool* const found, float* const score);

  // Pops blocks from the 'partition_' until a new path can be explored.
  // Returns false if none found or returns true and the 'next_block_size'.
  bool FindNewBranch(BlockSize* const next_block_size);

  TileScoreFunc* const tile_score_func_;  // Specific 'score_func_'.

  // For a given block size, stores the next one in the current partition set.
  std::array<BlockSize, BLK_LAST + 1> next_block_sizes_;

  VectorNoCtor<Block> forced_blocks_;
  FrontMgrLexico front_mgr_;       // To easily browse the partitions.
  VectorNoCtor<Block> partition_;  // Matches the layout of the 'front_mgr_'.

  // VDebug
  void RegisterScoreForVDebug(float best_partition_score,
                              size_t best_partition_size,
                              size_t num_iterations) const;
};

//------------------------------------------------------------------------------

// Uses several metrics to select block sizes.
class MultiPassPartitioner : public Partitioner {
 public:
  // This partitioner needs a specific MultiScoreFunc.
  explicit MultiPassPartitioner(MultiScoreFunc* const score_func)
      : multi_score_func_(score_func) {}

  WP2Status GetBestPartition(VectorNoCtor<Block>* const blocks) override;

 private:
  enum class Grid { Snapped, NonSnapped, All };

  // Selects blocks of dimension 'block_size', good enough according to the
  // given 'pass' and 'restrict_to' some alignment.
  WP2Status GetBlocks(MultiScoreFunc::Pass pass, BlockSize block_size,
                      Grid restrict_to, uint32_t* const max_num_blocks_left,
                      VectorNoCtor<Block>* const out);

  WP2Status SelectBlock(MultiScoreFunc::Pass pass, const Block& block,
                        uint32_t* const max_num_blocks_left,
                        VectorNoCtor<Block>* const out);

  MultiScoreFunc* const multi_score_func_;  // Specific 'score_func_'.

  // VDebug
  uint32_t pass_index_ = 0;

  void RegisterPassForVDebug(MultiScoreFunc::Pass pass, BlockSize block_size,
                             uint32_t num_chosen_blocks) const;
};

//------------------------------------------------------------------------------

}  // namespace WP2

#endif  // WP2_ENC_PARTITIONER_H_
