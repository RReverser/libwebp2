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
//   Block sizes and partitions
//
// Author: Skal (pascal.massimino@gmail.com)

#ifndef WP2_COMMON_LOSSY_BLOCK_SIZE_IO_H_
#define WP2_COMMON_LOSSY_BLOCK_SIZE_IO_H_

#include "src/dec/symbols_dec.h"
#include "src/enc/symbols_enc.h"
#include "src/utils/ans.h"
#include "src/utils/front_mgr.h"

namespace WP2 {

//------------------------------------------------------------------------------

class PartitionSetStats {
 public:
  static constexpr uint8_t kInvalidSlot = BLK_LAST;
  PartitionSetStats(PartitionSet partition_set, const BlockSize blocks[]);

 public:
  const PartitionSet partition_set_;
  uint8_t num_blocks_[BLK_LAST];              // Per bounding size.
  BlockSize blocks_[BLK_LAST][BLK_LAST + 1];  // Per bounding size, size array.
  uint8_t slots_[BLK_LAST][BLK_LAST];         // Per bounding size, index array.
  BlockSize smallest_bounds_[BLK_LAST];       // Per bounding size.
  uint32_t num_unique_bounds_;                // Among 'smallest_bounds_'.
  BlockSize unique_bounds_[BLK_LAST];         // At most 'BLK_LAST' elements.
  uint32_t unique_bounds_index_[BLK_LAST];    // Each is < 'num_unique_bounds_'.
  uint32_t used_trfs_;                        // Bitset for used transforms.
};
static_assert(TRF_LAST < 32u, "can't use 32b bitset for transforms used.");

extern const PartitionSetStats kPartSetsStats[];

//------------------------------------------------------------------------------

// Predefined partitioning sets of blocks for a given bounding BlockSize.
const BlockSize* GetBlockSizes(PartitionSet set, BlockSize bounds = BLK_32x32);
uint32_t GetNumBlockSizes(PartitionSet set, BlockSize bounds = BLK_32x32);
bool TrfIsUsed(PartitionSet set, TrfSize tdim);

// Returns one of the biggest BlockSizes from the 'set' fitting in 'bounds'.
BlockSize GetFittingBlockSize(PartitionSet set, BlockSize bounds);

// Returns true if 'size' belongs in 'set'.
bool IsBlockSizeInSet(PartitionSet set, BlockSize size);

// For given 'bounds', contains the smallest BlockSize at least as wide and
// tall as each BlockSize belonging to the 'set' fitting into these 'bounds'.
// Same as GetFittingBlockSize() unless the 'set' is sparse (meaning some
// width/height combinations are missing up to the biggest square).
BlockSize GetSmallestBounds(PartitionSet set, BlockSize bounds);

// Returns the number of bounds that lead to different subsets of BlockSizes.
uint32_t GetNumUniqueBounds(PartitionSet set);
// Returns the 'index'th among all unique smallest bounds.
BlockSize GetUniqueBounds(PartitionSet set, uint32_t index);

//------------------------------------------------------------------------------

// Writes 'dim'.
void WriteBlockSize(const FrontMgrNxNBase& mgr, BlockSize dim,
                    SymbolManager* const sm, ANSEncBase* const enc);
// Reads, updates and returns 'dim'.
BlockSize ReadBlockSize(const FrontMgrNxNBase& mgr, SymbolReader* const sr);

//------------------------------------------------------------------------------

}  // namespace WP2

#endif  // WP2_COMMON_LOSSY_BLOCK_SIZE_IO_H_
