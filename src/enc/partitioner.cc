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

#include "src/enc/partitioner.h"

#include <algorithm>
#include <numeric>

#include "src/common/lossy/block_size.h"
#include "src/common/lossy/block_size_io.h"
#include "src/enc/analysis.h"
#include "src/utils/utils.h"
#include "src/wp2/base.h"

namespace WP2 {

//------------------------------------------------------------------------------

WP2Status Partitioner::Init(const EncoderConfig& config, const YUVPlane& yuv,
                            const Rectangle& tile_rect,
                            PartitionScoreFunc* const score_func) {
  config_ = &config;
  tile_rect_ = tile_rect;
  src_ = &yuv;
  num_block_cols_ = SizeBlocks(yuv.Y.w_);
  num_block_rows_ = SizeBlocks(yuv.Y.h_);
  score_func_ = score_func;

  WP2_CHECK_ALLOC_OK(occupancy_.resize(num_block_cols_ * num_block_rows_));
  return WP2_STATUS_OK;
}

bool Partitioner::IsOccupied(const Block& block) const {
  const uint32_t stride = num_block_cols_;
  const bool* occupancy = occupancy_.data() + block.y() * stride + block.x();
  for (uint32_t y = 0; y < block.h(); ++y) {
    for (uint32_t x = 0; x < block.w(); ++x) {
      if (occupancy[x]) return true;
    }
    occupancy += num_block_cols_;
  }
  return false;
}

void Partitioner::Occupy(const Block& block) {
  const uint32_t stride = num_block_cols_;
  bool* occupancy = occupancy_.data() + block.y() * stride + block.x();
  for (uint32_t y = 0; y < block.h(); ++y) {
    for (uint32_t x = 0; x < block.w(); ++x) occupancy[x] = true;
    occupancy += num_block_cols_;
  }
}

bool Partitioner::IsBlockValid(const Block& block,
                                          const Block& max_block,
                                          size_t num_forced_blocks) const {
  if (block.w() > max_block.w()) return false;
  if (block.h() > max_block.h()) return false;
  // Thanks to the FrontMgrLexico, this 'block' can only be rejected
  // because of other forced blocks or forced snapping.
  if (config_->partition_snapping && !block.IsSnapped()) return false;
  if (num_forced_blocks > 0) return !IsOccupied(block);
  assert(!IsOccupied(block));
  return true;
}

//------------------------------------------------------------------------------

WP2Status AddForcedBlocks(const EncoderConfig& config,
                          const Rectangle& padded_tile_rect,
                          VectorNoCtor<Block>* const out) {
  if (config.info == nullptr) return WP2_STATUS_OK;

  for (const Rectangle& rect_px : config.info->force_partition) {
    if (rect_px.x < padded_tile_rect.x || rect_px.y < padded_tile_rect.y) {
      continue;  // Ignore this block, it's outside of this tile.
    }
    const Rectangle local_rect = {rect_px.x - padded_tile_rect.x,
                                  rect_px.y - padded_tile_rect.y,
                                  rect_px.width, rect_px.height};
    if (local_rect.x >= padded_tile_rect.width ||
        local_rect.y >= padded_tile_rect.height) {
      continue;  // Ignore this block, it's outside of this tile.
    }
    WP2_CHECK_OK((local_rect.x + local_rect.width <= padded_tile_rect.width &&
                  local_rect.y + local_rect.height <= padded_tile_rect.height),
                 WP2_STATUS_INVALID_CONFIGURATION);

    WP2::BlockSize matching_size = BLK_LAST;
    for (const BlockSize* size = GetBlockSizes(config.partition_set);
         *size != BLK_LAST; ++size) {
      if (local_rect.width == WP2::BlockWidthPix(*size) &&
          local_rect.height == WP2::BlockHeightPix(*size)) {
        matching_size = *size;
        break;
      }
    }
    WP2_CHECK_OK(matching_size != BLK_LAST, WP2_STATUS_INVALID_CONFIGURATION);

    // Create a Block that matches the rectangle.
    Block block(local_rect.x / kMinBlockSizePix,
                local_rect.y / kMinBlockSizePix, matching_size);
    assert(local_rect.width == block.w_pix() &&
           local_rect.height == block.h_pix());
    // Verify that position coordinates are aligned on the block grid (multiples
    // of kMinBlockSizePix). Width and height are checked and asserted above.
    WP2_CHECK_OK(local_rect.x == block.x_pix() && local_rect.y == block.y_pix(),
                 WP2_STATUS_INVALID_CONFIGURATION);
    // No IsSnapped() check, trust the user's intention on that.

    WP2_CHECK_ALLOC_OK(out->push_back(block));
  }
  return WP2_STATUS_OK;
}

WP2Status Partitioner::RegisterForcedBlocks(
    const VectorNoCtor<Block>& forced_blocks,
    uint32_t* const max_num_blocks_left) {
  for (const Block& block : forced_blocks) {
    WP2_CHECK_OK(!IsOccupied(block), WP2_STATUS_INVALID_CONFIGURATION);
    Occupy(block);
    assert(*max_num_blocks_left >= block.w() * block.h());
    *max_num_blocks_left -= block.w() * block.h();
  }
  return WP2_STATUS_OK;
}

//------------------------------------------------------------------------------

WP2Status TopLeftBestPartitioner::GetBestPartition(
    VectorNoCtor<Block>* const blocks) {
  std::fill(occupancy_.begin(), occupancy_.end(), false);
  const uint32_t max_num_blocks = num_block_cols_ * num_block_rows_;
  uint32_t max_num_blocks_left = max_num_blocks;

  // Begin with the user-defined blocks, if any (debug).
  WP2_CHECK_STATUS(RegisterForcedBlocks(*blocks, &max_num_blocks_left));

  // Sort forced blocks by position for binary search below.
  VectorNoCtor<Block> forced_blocks;
  WP2_CHECK_ALLOC_OK(forced_blocks.resize(blocks->size()));
  std::copy(blocks->begin(), blocks->end(), forced_blocks.begin());
  std::sort(forced_blocks.begin(), forced_blocks.end());

  // Use a front manager to easily browse the blocks.
  FrontMgrLexico front_mgr;
  WP2_CHECK_STATUS(front_mgr.Init(config_->partition_set,
                                  config_->partition_snapping,
                                  src_->Y.w_, src_->Y.h_));

  while (!front_mgr.Done()) {
    // TODO(yguyon): This is using a top-most first order. Find an order that
    //               maximizes top and left context instead.
    const Block max_block = front_mgr.GetMaxFittingBlock();

    // Binary search to find if 'max_block' is already one of the forced blocks.
    const auto forced_block_it =
        std::lower_bound(forced_blocks.begin(), forced_blocks.end(), max_block);

    Block best_block;
    if (forced_block_it != forced_blocks.end() &&
        forced_block_it->x() == max_block.x() &&
        forced_block_it->y() == max_block.y()) {
      best_block = *forced_block_it;
    } else {
      assert(max_num_blocks_left > 0);
      WP2_CHECK_STATUS(
          GetBestSize(max_block, forced_blocks.size(), &best_block));
      RegisterOrderForVDebug(0, 0, best_block, blocks->size(), max_num_blocks);
      WP2_CHECK_ALLOC_OK(blocks->push_back(best_block));

      // The following instructions are done only for assertions/debug
      // because the front manager should handle everything, except for
      // the forced blocks but these are checked above.
      assert((Occupy(best_block), true));
      assert(max_num_blocks_left >= best_block.w() * best_block.h());
      max_num_blocks_left -= best_block.w() * best_block.h();
    }
    assert(best_block.dim() != BLK_LAST);
    WP2_CHECK_ALLOC_OK(front_mgr.UseSize(best_block.dim(),
                                         /*ind=*/0u, nullptr));
    front_mgr.Use(best_block);
    WP2_CHECK_STATUS(score_func_->Use(best_block));
  }
  assert(max_num_blocks_left == 0);
  return WP2_STATUS_OK;
}

WP2Status TopLeftBestPartitioner::GetBestSize(const Block& max_block,
                                              size_t num_forced_blocks,
                                              Block* const best_block) {
  const BlockSize* const sizes = GetBlockSizes(config_->partition_set);
  float best_score = 0.f;
  // Test each size and keep the best score.
  for (const BlockSize* size = sizes; *size != BLK_LAST; ++size) {
    const Block block = Block(max_block.x(), max_block.y(), *size);

    if (IsBlockValid(block, max_block, num_forced_blocks)) {
      float score;
      WP2_CHECK_STATUS(score_func_->GetScore(block, &score));

      const bool first_score = (size == sizes);
      // Increasing sizes: keep the largest one if same score (rare).
      if (first_score || score >= best_score) {
        best_score = score;
        *best_block = block;
      }
    }
  }
  return WP2_STATUS_OK;
}

//------------------------------------------------------------------------------

WP2Status AreaRDOptPartitioner::GetBestPartition(
    VectorNoCtor<Block>* const blocks) {
  std::fill(occupancy_.begin(), occupancy_.end(), false);
  uint32_t max_num_blocks_left = num_block_cols_ * num_block_rows_;

  // TODO(yguyon): Handle forced blocks
  WP2_CHECK_OK(blocks->empty(), WP2_STATUS_UNSUPPORTED_FEATURE);

  for (uint32_t area_y = 0; area_y < src_->Y.h_; area_y += area_height_) {
    for (uint32_t area_x = 0; area_x < src_->Y.w_; area_x += area_width_) {
      const Rectangle area =
          Rectangle(area_x, area_y, area_width_, area_height_)
              .ClipWith({0, 0, src_->Y.w_, src_->Y.h_});  // Padded.

      WP2_CHECK_STATUS(
          GetAreaBestPartition(area, blocks, &max_num_blocks_left));
    }
  }
  assert(max_num_blocks_left == 0);
  return WP2_STATUS_OK;
}

WP2Status AreaRDOptPartitioner::GetAreaBestPartition(
    const Rectangle& area, VectorNoCtor<Block>* const blocks,
    uint32_t* const max_num_blocks_left) {
  AreaScoreFunc* const score_func =
      reinterpret_cast<AreaScoreFunc*>(score_func_);
  assert(area.x == score_func->GetCurrentArea().x &&
         area.y == score_func->GetCurrentArea().y);
  VectorNoCtor<Block> area_blocks, area_best_blocks;

  // Start with the default partition as a reference score.
  float best_score;
  WP2_CHECK_STATUS(
      score_func->GetAreaDefaultScore(&area_best_blocks, &best_score));

  // Test a uniform block grid for each block size and keep the best score.
  for (const BlockSize* size = GetBlockSizes(config_->partition_set);
       *size != BLK_LAST; ++size) {
    if (BlockWidthPix(*size) <= area.width &&
        BlockHeightPix(*size) <= area.height) {
      area_blocks.clear();
      float score;
      WP2_CHECK_STATUS(
          score_func->GetAreaGridScore(*size, &area_blocks, &score));
      if (score > best_score) {
        best_score = score;
        WP2_CHECK_ALLOC_OK(area_best_blocks.resize(area_blocks.size()));
        std::copy(area_blocks.begin(), area_blocks.end(),
                  area_best_blocks.begin());
      }
    }
  }

  // Mark each block of the best partition as final.
  for (const Block& block : area_best_blocks) {
    assert(*max_num_blocks_left > 0);
    *max_num_blocks_left -= block.rect().GetArea();
    WP2_CHECK_ALLOC_OK(blocks->push_back(block));
    WP2_CHECK_STATUS(score_func_->Use(block));
  }
  return WP2_STATUS_OK;
}

//------------------------------------------------------------------------------

WP2Status SubAreaRDOptPartitioner::GetBestPartition(
    VectorNoCtor<Block>* const blocks) {
  WP2_CHECK_STATUS(front_mgr_.Init(config_->partition_set,
                                   config_->partition_snapping,
                                   tile_rect_.width, tile_rect_.height));
  return AreaRDOptPartitioner::GetBestPartition(blocks);
}

WP2Status SubAreaRDOptPartitioner::GetAreaBestPartition(
    const Rectangle& area, VectorNoCtor<Block>* const blocks,
    uint32_t* const max_num_blocks_left) {
  const uint32_t max_num_blocks = num_block_cols_ * num_block_rows_;
  assert(area.x == front_mgr_.GetMaxPossibleBlock().x_pix() &&
         area.y == front_mgr_.GetMaxPossibleBlock().y_pix());

  while (!front_mgr_.Done()) {
    const Block max_block = front_mgr_.GetMaxPossibleBlock();
    // Exit the loop when the 'area' is complete.
    if (!area.Contains(max_block.x_pix(), max_block.y_pix())) break;

    assert(*max_num_blocks_left > 0);
    Block best_block;
    float best_score = std::numeric_limits<float>::lowest();
    const BlockSize* const sizes = GetBlockSizes(config_->partition_set);
    for (const BlockSize* size = sizes; *size != BLK_LAST; ++size) {
      const Block block = Block(max_block.x(), max_block.y(), *size);
      if (block.w() <= max_block.w() && block.h() <= max_block.h()) {
        float score;
        WP2_CHECK_STATUS(score_func_->GetScore(block, &score));

        if (score >= best_score) {
          best_score = score;
          best_block = block;
        }
      }
    }
    const uint32_t block_index = max_num_blocks - *max_num_blocks_left;
    RegisterOrderForVDebug(0, 0, best_block, block_index, max_num_blocks);
    *max_num_blocks_left -= best_block.rect().GetArea();
    WP2_CHECK_ALLOC_OK(blocks->push_back(best_block));
    WP2_CHECK_ALLOC_OK(front_mgr_.UseSize(best_block.dim(), 0u, nullptr));
    front_mgr_.Use(best_block);
    WP2_CHECK_STATUS(score_func_->Use(best_block));
  }
  return WP2_STATUS_OK;
}

//------------------------------------------------------------------------------

constexpr SplitRecursePartitioner::SplitInfo
    SplitRecursePartitioner::kNoSplitInfo;

WP2Status SplitRecursePartitioner::Init(const EncoderConfig& config,
                                        const YUVPlane& yuv,
                                        const Rectangle& tile_rect,
                                        PartitionScoreFunc* const score_func) {
  WP2_CHECK_STATUS(Partitioner::Init(config, yuv, tile_rect, score_func));
  WP2_CHECK_STATUS(split_info_.Resize(num_block_cols_, num_block_rows_));
  return WP2_STATUS_OK;
}

WP2Status SplitRecursePartitioner::GetBestPartition(
    VectorNoCtor<Block>* const blocks) {
  std::fill(occupancy_.begin(), occupancy_.end(), false);
  split_info_.Fill(kNoSplitInfo);
  const uint32_t max_num_blocks = num_block_cols_ * num_block_rows_;
  uint32_t max_num_blocks_left = max_num_blocks;

  // Begin with the user-defined blocks, if any (debug).
  WP2_CHECK_STATUS(RegisterForcedBlocks(*blocks, &max_num_blocks_left));

  // Sort forced blocks by position for binary search below.
  VectorNoCtor<Block> forced_blocks;
  WP2_CHECK_ALLOC_OK(forced_blocks.resize(blocks->size()));
  std::copy(blocks->begin(), blocks->end(), forced_blocks.begin());
  std::sort(forced_blocks.begin(), forced_blocks.end());

  // Use a front manager to easily browse the blocks.
  FrontMgrLexico front_mgr;
  WP2_CHECK_STATUS(front_mgr.Init(config_->partition_set,
                                  config_->partition_snapping, src_->Y.w_,
                                  src_->Y.h_));
  Vector<SplitPattern> patterns;
  WP2_CHECK_ALLOC_OK(patterns.reserve(10));  // maximum 10 ways to split a block
  BlockScoreFunc* const score_func = (BlockScoreFunc*)score_func_;

  while (!front_mgr.Done()) {
    const Block max_block = front_mgr.GetMaxPossibleBlock();

    // Binary search to find if 'max_block' is already one of the forced blocks.
    const auto forced_block_it =
        std::lower_bound(forced_blocks.begin(), forced_blocks.end(), max_block);

    if (forced_block_it != forced_blocks.end() &&
        forced_block_it->x() == max_block.x() &&
        forced_block_it->y() == max_block.y()) {
      const Block block = *forced_block_it;
      assert(block.dim() != BLK_LAST);
      WP2_CHECK_ALLOC_OK(front_mgr.UseSize(block.dim(), /*ind=*/0u, nullptr));
      front_mgr.Use(block);
      WP2_CHECK_STATUS(score_func_->Use(block));
    } else {
      assert(max_num_blocks_left > 0);

      // Retrieve the previous decision to split this block, if any.
      const SplitInfo& split_info = GetSplitInfo(max_block.x(), max_block.y());
      if (split_info == kNoSplitInfo || split_info.recursive) {
        // Nothing decided or recursive allowed: try more split patterns.
        Block block;
        if (split_info == kNoSplitInfo) {
          block =
              Block(max_block.x(), max_block.y(),
                    GetMaxFittingBlockSize(max_block, forced_blocks.size()));
        } else {
          block = Block(max_block.x(), max_block.y(), split_info.size);
        }
        assert(IsBlockValid(block, max_block, forced_blocks.size()));
        assert(IsSplitInfo(block, split_info));

        // Get all split patterns and keep the one with the best score.
        patterns.clear();
        GetSplitPatterns(config_->partition_set, block.dim(), &patterns);
        SplitPattern best_pattern;
        if (patterns.size() == 1) {
          best_pattern = patterns.front();
        } else {
          best_pattern = {/*num_blocks=*/0, {BLK_LAST}};
          float best_score = 0.f;
          for (const SplitPattern& pattern : patterns) {
            Block pattern_blocks[4];
            GetSplitPatternBlocks(block, pattern, pattern_blocks);

            float score;
            WP2_CHECK_STATUS(score_func->GetScore(pattern_blocks,
                                                  pattern.num_blocks, &score));
            if (best_pattern.num_blocks == 0 || score > best_score) {
              best_pattern = pattern;
              best_score = score;
            }
          }
          assert(best_pattern.num_blocks > 0);
        }

        // Remember the decision.
        Block best_pattern_blocks[4];
        GetSplitPatternBlocks(block, best_pattern, best_pattern_blocks);
        for (uint32_t i = 0; i < best_pattern.num_blocks; ++i) {
          const bool recursive = (best_pattern.num_blocks > 1 &&
                                  best_pattern_blocks[i].dim() != BLK_4x4);
          SetSplitInfo(best_pattern_blocks[i], recursive);
          // TODO(yguyon): Recurse in fewer situations to save CPU.
        }
      } else {
        // This 'block' was decided as is with no possibility of recursion.
        const Block block =
            Block(max_block.x(), max_block.y(), split_info.size);
        assert(IsBlockValid(block, max_block, forced_blocks.size()));
        assert(IsSplitInfo(block, split_info));
        RegisterOrderForVDebug(0, 0, block, blocks->size(), max_num_blocks);
        WP2_CHECK_ALLOC_OK(blocks->push_back(block));
        WP2_CHECK_ALLOC_OK(front_mgr.UseSize(block.dim(), /*ind=*/0u, nullptr));
        front_mgr.Use(block);
        WP2_CHECK_STATUS(score_func_->Use(block));

        // The following instructions are done only for assertions/debug
        // because the front manager should handle everything, except for
        // the forced blocks but these are checked above.
        assert((Occupy(block), true));
        assert(max_num_blocks_left >= block.w() * block.h());
        max_num_blocks_left -= block.w() * block.h();
        assert(block.dim() != BLK_LAST);
      }
    }
  }
  assert(max_num_blocks_left == 0);
  return WP2_STATUS_OK;
}

BlockSize SplitRecursePartitioner::GetMaxFittingBlockSize(
    const Block& max_block, uint32_t num_forced_blocks) const {
  BlockSize best_size = BLK_LAST;
  const BlockSize* const sizes = GetBlockSizes(config_->partition_set);
  for (const BlockSize* size = sizes; *size != BLK_LAST; ++size) {
    const Block block = Block(max_block.x(), max_block.y(), *size);

    if ((best_size == BLK_LAST ||
         BlockWidth[*size] * BlockHeight[*size] >
             BlockWidth[best_size] * BlockHeight[best_size]) &&
        IsBlockValid(block, max_block, num_forced_blocks) &&
        IsSplitInfo(block, kNoSplitInfo)) {
      best_size = block.dim();
    }
  }
  assert(best_size != BLK_LAST);
  return best_size;
}

static BlockSize GetExactBlockSize(PartitionSet set, uint32_t w, uint32_t h) {
  const BlockSize s = GetBlockSize(std::max(1u, w), std::max(1u, h));
  if (IsBlockSizeInSet(set, s) && BlockWidth[s] == w && BlockHeight[s] == h) {
    return s;
  }
  return BLK_LAST;
}

void SplitRecursePartitioner::GetSplitPatterns(PartitionSet set, BlockSize size,
                                               Vector<SplitPattern>* const p) {
  assert(p->empty() && p->capacity() >= 10);
  if (IsBlockSizeInSet(set, size)) p->push_back_no_resize({1, {size}});

  const uint32_t width = BlockWidth[size], height = BlockHeight[size];
  // The BlockSize variable names below are based on a 4x4 block reference.
  const BlockSize b2x2 = GetExactBlockSize(set, width / 2, height / 2);
  const BlockSize b4x1 = GetExactBlockSize(set, width, height / 4);
  const BlockSize b1x4 = GetExactBlockSize(set, width / 4, height);
  const BlockSize b4x2 = GetExactBlockSize(set, width, height / 2);
  const BlockSize b2x4 = GetExactBlockSize(set, width / 2, height);

  if (b2x2 != BLK_LAST) p->push_back_no_resize({4, {b2x2, b2x2, b2x2, b2x2}});
  if (b4x1 != BLK_LAST) p->push_back_no_resize({4, {b4x1, b4x1, b4x1, b4x1}});
  if (b1x4 != BLK_LAST) p->push_back_no_resize({4, {b1x4, b1x4, b1x4, b1x4}});
  if (b4x2 != BLK_LAST) {
    p->push_back_no_resize({2, {b4x2, b4x2}});
    if (b2x2 != BLK_LAST) p->push_back_no_resize({3, {b4x2, b2x2, b2x2}});
    if (b2x2 != BLK_LAST) p->push_back_no_resize({3, {b2x2, b2x2, b4x2}});
  }
  if (b2x4 != BLK_LAST) {
    p->push_back_no_resize({2, {b2x4, b2x4}});
    if (b2x2 != BLK_LAST) p->push_back_no_resize({3, {b2x4, b2x2, b2x2}});
    if (b2x2 != BLK_LAST) p->push_back_no_resize({3, {b2x2, b2x4, b2x2}});
  }
  assert(!p->empty());
}

void SplitRecursePartitioner::GetSplitPatternBlocks(const Block& block,
                                                    const SplitPattern& pattern,
                                                    Block blocks[4]) {
  assert(pattern.num_blocks >= 1 && pattern.num_blocks <= 4);
  // Treat this special case differently. Others can be done row by row.
  const bool is_2x4_2x2_2x2 =
      (pattern.num_blocks == 3 &&
       BlockHeight[pattern.block_sizes[0]] == block.h());
  uint32_t num_blocks = 0, x = block.x(), y = block.y();
  while (num_blocks < pattern.num_blocks) {
    if (x == block.x() + block.w()) {  // End of a row.
      // Start a new row aligned under the first block, unless 'is_2x4_2x2_2x2'.
      x = blocks[is_2x4_2x2_2x2 ? 1 : 0].x();
      y += blocks[is_2x4_2x2_2x2 ? 1 : 0].h();
    }
    assert(x < block.x() + block.w() && y < block.y() + block.h());
    blocks[num_blocks] = Block(x, y, pattern.block_sizes[num_blocks]);
    x += blocks[num_blocks].w();
    ++num_blocks;
  }
}

const SplitRecursePartitioner::SplitInfo& SplitRecursePartitioner::GetSplitInfo(
    uint32_t x, uint32_t y) const {
  return split_info_.At(x, y);
}
bool SplitRecursePartitioner::IsSplitInfo(const Block& block,
                                          const SplitInfo& split_info) const {
  const SplitInfo* row = &split_info_.At(block.x(), block.y());
  for (uint32_t y = 0; y < block.h(); ++y, row += split_info_.Step()) {
    for (uint32_t x = 0; x < block.w(); ++x) {
      if (row[x] != split_info) return false;
    }
  }
  return true;
}
void SplitRecursePartitioner::SetSplitInfo(const Block& block, bool recursive) {
  SplitInfo* row = &split_info_.At(block.x(), block.y());
  for (uint32_t y = 0; y < block.h(); ++y, row += split_info_.Step()) {
    for (uint32_t x = 0; x < block.w(); ++x) row[x] = {block.dim(), recursive};
  }
}

//------------------------------------------------------------------------------

static constexpr uint32_t kBlkStart = BLK_LAST;

WP2Status ExhaustivePartitioner::Init(const EncoderConfig& config,
                                      const YUVPlane& yuv,
                                      const Rectangle& tile_rect,
                                      PartitionScoreFunc* const score_func) {
  WP2_CHECK_STATUS(Partitioner::Init(config, yuv, tile_rect, score_func));

  const BlockSize* size = GetBlockSizes(config_->partition_set);
  next_block_sizes_[kBlkStart] = *size;  // First one.
  for (; *size != BLK_LAST; ++size) next_block_sizes_[*size] = *(size + 1);

  WP2_CHECK_STATUS(front_mgr_.Init(config_->partition_set,
                                   config_->partition_snapping,
                                   src_->Y.w_, src_->Y.h_));
  return WP2_STATUS_OK;
}

WP2Status ExhaustivePartitioner::GetBestPartition(
    VectorNoCtor<Block>* const blocks) {
  std::fill(occupancy_.begin(), occupancy_.end(), false);

  WP2_CHECK_ALLOC_OK(forced_blocks_.resize(blocks->size()));
  std::copy(blocks->begin(), blocks->end(), forced_blocks_.begin());
  // Only occupy the forced blocks, the others are determined by 'front_mgr_'.
  for (const Block& block : forced_blocks_) Occupy(block);
  // Sort the 'forced_blocks_' by position for binary search later on.
  std::sort(forced_blocks_.begin(), forced_blocks_.end());

  float best_score = 0.f;
  VectorNoCtor<Block>* const best_partition = blocks;
  best_partition->clear();
  front_mgr_.Clear();
  partition_.clear();

  size_t num_iterations = 0;
  while (true) {
    bool found;
    float score;
    WP2_CHECK_STATUS(GetNextPartitionScore(&found, &score));
    if (!found) break;

    // Make sure the partition fills exactly the tile.
    assert(std::accumulate(partition_.begin(), partition_.end(), 0u,
                           [](uint32_t sum, const Block& block) {
                             return sum + (BlockWidth[block.dim()] *
                                           BlockHeight[block.dim()]);
                           }) == num_block_cols_ * num_block_rows_);

    if (best_partition->empty() || score > best_score) {
      best_score = score;
      WP2_CHECK_ALLOC_OK(best_partition->resize(partition_.size()));
      std::copy(partition_.begin(), partition_.end(), best_partition->begin());
    }
    ++num_iterations;
  }

  RegisterScoreForVDebug(best_score, best_partition->size(), num_iterations);
  return WP2_STATUS_OK;
}

WP2Status ExhaustivePartitioner::GetNextPartitionScore(bool* const found,
                                                       float* const score) {
  // Begin the tree search.
  BlockSize next_block_size = next_block_sizes_[kBlkStart];
  bool is_first_branch = partition_.empty();
  // Find a new branch in the search tree if possible.
  while (is_first_branch || FindNewBranch(&next_block_size)) {
    is_first_branch = false;
    // Try to go down the new branch.
    do {
      Block block;
      VectorNoCtor<Block>::iterator forced_block_it;
      // Check if the 'next_block_size' can fit, then if the 'block' is snapped
      // or does not need to, then if the 'block' is forced or intersecting one.
      if (front_mgr_.TryGetNextBlock(next_block_size, &block) &&
          (!config_->partition_snapping || block.IsSnapped()) &&
          (forced_blocks_.empty() || !IsOccupied(block) ||
           ((forced_block_it = std::lower_bound(forced_blocks_.begin(),
                                                forced_blocks_.end(), block)) !=
                forced_blocks_.end() &&
            *forced_block_it == block))) {
        if (forced_blocks_.empty()) assert(!IsOccupied(block));

        // 'block' is valid.
        WP2_CHECK_ALLOC_OK(front_mgr_.UseSize(next_block_size,
                                              /*ind=*/0u, &block));
        front_mgr_.Use(block);
        WP2_CHECK_ALLOC_OK(partition_.push_back(block));

        if (front_mgr_.Done()) {
          // 'partition_' is complete. Get its score.
          WP2_CHECK_STATUS(tile_score_func_->TryEncode(partition_, score));
          *found = true;
          return WP2_STATUS_OK;
        } else {
          // 'partition_' is incomplete. Go deeper.
          next_block_size = next_block_sizes_[kBlkStart];
        }
      } else {
        // 'block' is invalid. Try the next size or leave this branch.
        next_block_size = next_block_sizes_[next_block_size];
      }
    } while (next_block_size != BLK_LAST);
  }

  *found = false;
  return WP2_STATUS_OK;
}

bool ExhaustivePartitioner::FindNewBranch(BlockSize* const next_block_size) {
  while (!partition_.empty()) {
    *next_block_size = next_block_sizes_[partition_.back().dim()];
    front_mgr_.UndoUse(partition_.back());
    front_mgr_.UndoUseSize(partition_.back());
    partition_.pop_back();
    if (*next_block_size != BLK_LAST) return true;  // Found new branch.
  }
  return false;  // Browsed the whole tree.
}

//------------------------------------------------------------------------------

using Pass = MultiScoreFunc::Pass;

WP2Status MultiPassPartitioner::GetBestPartition(
    VectorNoCtor<Block>* const blocks) {
  std::fill(occupancy_.begin(), occupancy_.end(), false);
  uint32_t max_num_blocks_left = num_block_cols_ * num_block_rows_;
  WP2_CHECK_STATUS(RegisterForcedBlocks(*blocks, &max_num_blocks_left));

  // Run several passes to select matching blocks with a given size.
  // GetBlocks() will skip any pass with a BlockSize not belonging to the
  // current PartitionSet.

  // Start with easy-to-spot flat blocks (narrow range of values in the
  // original luma and alpha planes).
  for (BlockSize block_size : {BLK_32x32, BLK_32x16, BLK_16x32, BLK_16x16,
                               BLK_16x8, BLK_8x16, BLK_8x8}) {
    WP2_CHECK_STATUS(GetBlocks(Pass::LumaAlphaGradient, block_size,
                               Grid::Snapped, &max_num_blocks_left, blocks));
  }
  if (config_->speed >= MultiScoreFunc::kMinSpeedForGoodQuantDCT) {
    // Blocks having a DCT that still gives a good result when quantized.
    WP2_CHECK_STATUS(GetBlocks(Pass::GoodQuantDCT, BLK_32x32, Grid::Snapped,
                               &max_num_blocks_left, blocks));
  }

  // Pick a few not-too-small blocks that are following a luma edge.
  for (const BlockSize* block_size =
           GetBlockSizes(ALL_RECTS) + GetNumBlockSizes(ALL_RECTS) - 1u;
       *block_size != BLK_4x4; --block_size) {
    if (BlockWidth[*block_size] * BlockHeight[*block_size] < 4) continue;
    WP2_CHECK_STATUS(GetBlocks(Pass::Direction, *block_size, Grid::All,
                               &max_num_blocks_left, blocks));
  }

  // Choose smaller and smaller blocks based on variance.
  constexpr PartitionSet variance_partition_set = SMALL_RECTS;
  assert(*GetBlockSizes(variance_partition_set) == BLK_4x4);
  for (const BlockSize* allowed_block_size =
           GetBlockSizes(variance_partition_set) +
           GetNumBlockSizes(variance_partition_set) - 1u;
       *allowed_block_size != BLK_4x4; --allowed_block_size) {
    WP2_CHECK_STATUS(GetBlocks(Pass::NarrowStdDev, *allowed_block_size,
                               Grid::All, &max_num_blocks_left, blocks));
  }

  // Merge all remaining small blocks that are following a luma edge.
  for (const BlockSize* block_size =
           GetBlockSizes(ALL_RECTS) + GetNumBlockSizes(ALL_RECTS) - 1u;
       *block_size != BLK_4x4; --block_size) {
    if (BlockWidth[*block_size] * BlockHeight[*block_size] >= 4) continue;
    WP2_CHECK_STATUS(GetBlocks(Pass::Direction, *block_size, Grid::All,
                               &max_num_blocks_left, blocks));
  }

  // Fill the remaining empty areas.
  WP2_CHECK_STATUS(
      GetBlocks(Pass::Any, BLK_4x4, Grid::All, &max_num_blocks_left, blocks));
  assert(max_num_blocks_left == 0);
  return WP2_STATUS_OK;
}

WP2Status MultiPassPartitioner::GetBlocks(Pass pass,
                                          BlockSize block_size,
                                          Grid restrict_to,
                                          uint32_t* const max_num_blocks_left,
                                          VectorNoCtor<Block>* const out) {
  if (GetFittingBlockSize(config_->partition_set, block_size) != block_size) {
    return WP2_STATUS_OK;
  }

  // Only keep snapped blocks if the config forces it.
  if (config_->partition_snapping) {
    if (restrict_to == Grid::NonSnapped) return WP2_STATUS_OK;
    restrict_to = Grid::Snapped;
  }

  const uint32_t num_blocks_before = out->size();
  Block block(0, 0, block_size);

  // Make sure the 'block_size' is valid and there is at least one remaining.
  if (block.w() > num_block_cols_ || block.h() > num_block_rows_ ||
      block.rect().GetArea() > *max_num_blocks_left) {
    RegisterPassForVDebug(pass, block_size, out->size() - num_blocks_before);
    ++pass_index_;
    return WP2_STATUS_OK;
  }

  struct BlockScore {
    Block blk;
    float score;
  };
  VectorNoCtor<BlockScore> block_scores;
  multi_score_func_->SetPass(pass);

  uint32_t x = 0, y = 0;
  const uint32_t max_x = num_block_cols_ - block.w();
  const uint32_t max_y = num_block_rows_ - block.h();
  const bool* occupancy = occupancy_.data();

  // Browse the grid of selectable blocks.
  while (y <= max_y) {
    while (x <= max_x) {
      if (!occupancy[x]) {  // Quickly skip top-left-occupied blocks.
        block.SetXY(x, y);
        if (restrict_to == Grid::Snapped) assert(block.IsSnapped());
        // Assume snapped blocks were already tried if 'restrict_to==NonSnapped'
        const bool tried_block_earlier =
            (restrict_to == Grid::NonSnapped && block.IsSnapped());

        if (!tried_block_earlier && !IsOccupied(block)) {
          float score;
          WP2_CHECK_STATUS(score_func_->GetScore(block, &score));
          if (score >= MultiScoreFunc::kMinScore) {
            if (restrict_to == Grid::Snapped) {
              // Immediately select any passing block.
              WP2_CHECK_STATUS(
                  SelectBlock(pass, block, max_num_blocks_left, out));
            } else {
              // Keep the scores of the passing blocks in memory.
              WP2_CHECK_ALLOC_OK(block_scores.push_back({block, score}));
            }
          }
        }
      }

      x += (restrict_to == Grid::Snapped) ? block.w() : 1;
    }
    const uint32_t num_incr_rows =
        (restrict_to == Grid::Snapped) ? block.h() : 1;
    y += num_incr_rows;
    occupancy += num_incr_rows * num_block_cols_;
    x = 0;
  }

  if (restrict_to != Grid::Snapped) {
    // Select blocks with best (=lowest) score first.
    // stable_sort() is used for its determinism.
    std::stable_sort(block_scores.begin(), block_scores.end(),
                     [](const BlockScore& l, const BlockScore& r) {
                       return l.score < r.score;
                     });

    for (const BlockScore& bs : block_scores) {
      if (IsOccupied(bs.blk)) continue;  // 'occupancy_' might have changed.
      WP2_CHECK_STATUS(SelectBlock(pass, bs.blk, max_num_blocks_left, out));
    }
  }

  RegisterPassForVDebug(pass, block_size, out->size() - num_blocks_before);
  ++pass_index_;
  return WP2_STATUS_OK;
}

WP2Status MultiPassPartitioner::SelectBlock(MultiScoreFunc::Pass pass,
                                            const Block& block,
                                            uint32_t* const max_num_blocks_left,
                                            VectorNoCtor<Block>* const out) {
  RegisterOrderForVDebug((uint32_t)pass, pass_index_, block, out->size(),
                         num_block_cols_ * num_block_rows_);
  Occupy(block);
  const uint32_t block_area = block.w() * block.h();
  assert(*max_num_blocks_left >= block_area);
  *max_num_blocks_left -= block_area;
  WP2_CHECK_ALLOC_OK(out->push_back(block));
  WP2_CHECK_STATUS(score_func_->Use(block));
  return WP2_STATUS_OK;
}

//------------------------------------------------------------------------------

}  // namespace WP2
