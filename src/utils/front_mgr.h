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
// Front manager: common struct for encoder/decoder, handling the
// front line of unprocessed blocks.
//
// Authors: Skal (pascal.massimino@gmail.com)

#ifndef WP2_UTILS_FRONT_MGR_H_
#define WP2_UTILS_FRONT_MGR_H_

#include <array>

#include "src/common/lossy/block_size.h"
#include "src/utils/vector.h"
#include "src/wp2/base.h"

namespace WP2 {

// Base class for a simple front manager that only fills the wave front
// and keeps track of the filling height of the columns.
class FrontMgrBase {
 public:
  // 'width' and 'height' are the dimensions of the current canvas in pixels.
  FrontMgrBase() = default;
  virtual ~FrontMgrBase() = default;
  FrontMgrBase(const FrontMgrBase&) = delete;

  // Allocates storage. Returns false in case of memory error.
  // Must be called first (once) before using the other methods.
  virtual WP2Status Init(PartitionSet partition_set, bool snapped,
                         uint32_t width, uint32_t height);

  // Clears the front manager and set it to an initial state.
  virtual void Clear();

  // Returns the y position of the first free sub-block at column position 'x'
  // (hence, returns 0 if column is empty). Note that parameter 'x' can be out
  // of bound, in which case it returns 0.
  // 'x' and the returned value are in kMinBlockSizePix units.
  uint32_t GetOccupancy(int x) const;

  // Returns the extent of the top context to the right of the block top in
  // kMinBlockSizePix units.
  uint32_t GetRightContextExtent(const Block& blk) const;

  // Returns true if the block at position x, y (in kMinBlockSizePix units)
  // has already been marked as 'used'.
  bool IsDone(uint32_t x, uint32_t y) const;

  // Returns true if all the surface is covered.
  bool Done() const { return (left_ == 0); }

  PartitionSet GetPartitionSet() const { return partition_set_; }

 protected:
  WP2Status CopyInternal(const FrontMgrBase& other);

  // Fills the 'occupancy' for a given block.
  void UseInternal(const Block& block, Vector_u32* const occupancy) const;
  void UndoUseInternal(const Block& block, Vector_u32* const occupancy) const;

  PartitionSet partition_set_ = NUM_PARTITION_SETS;  // Undefined.
  bool snapped_ = false;
  Vector_u32 occupancy_;
  uint32_t width_ = 0, height_ = 0;    // In pixels.
  uint32_t bwidth_ = 0, bheight_ = 0;  // In kMinBlockSizePix units.
  uint32_t left_ = 0;
};

//------------------------------------------------------------------------------

// Simple front manager that only fills 4x4 blocks.
// It is only useful to get the context of 4x4 blocks.
class FrontMgr4x4 : public FrontMgrBase {
 public:
  using FrontMgrBase::FrontMgrBase;

  WP2Status CopyFrom(const FrontMgr4x4& other);
  void Clear() override;

  // Use the next available block in lexicographic order.
  void Use();

 private:
  uint32_t bx_ = 0, by_ = 0;
};

//------------------------------------------------------------------------------

// Base class managing the reading/writing of NxN blocks in a given order.
// There are two things concerning a block: its size and everything else.
// Those two can be managed independently: the size of a block can depend on
// other block sizes according to different logics while the content is a given.
// There is a size order, and a block order.
// TODO(vrabaud): split this file in enc/dec to get lighter APIs.
class FrontMgrNxNBase : public FrontMgrBase {
 public:
  // Information of a stacked block. 'ind' is an index in whatever external
  // structure you want to refer to.
  struct Info {
    Block block;
    uint32_t ind;
  };

  // Allocates storage. Returns false in case of memory error.
  // Must be called first (once) before using the other methods.
  WP2Status Init(PartitionSet partition_set, bool snapped,
                 uint32_t width, uint32_t height) override;

  void Clear() override;

  // Returns all the underlying blocks of the tile.
  const VectorNoCtor<Block>& Blocks() const { return blocks_; }
  VectorNoCtor<Block>& Blocks() { return blocks_; }

  // Returns the indices of the blocks in Blocks() (that are sorted in block
  // order) so that they are in size order.
  const Vector_u16& SizeIndices() const;

  // Sorts the Blocks() in block order.
  virtual WP2Status Sort() = 0;

  // Block-order related functions.

  // Marks a given block as occupied.
  virtual void Use(const Block& block) = 0;

  // Returns whether an element is considered to be computable, i.e. its context
  // is considered filled at its maximum.
  virtual bool IsFinal(const WP2::Block& block) const = 0;

  // Size-order related functions.

  // Returns the biggest possible next block to use in size order.
  // This may not be a BlockSize allowed by the partition set.
  Block GetMaxPossibleBlock() const { return next_; }
  // Returns the biggest allowed next block to use in size order.
  Block GetMaxFittingBlock() const;

  // Uses a block in the size-order stack if it is final. Returns whether a
  // final block was popped.
  bool UseFinal(Info* const info = nullptr);

  // For a valid 'size', returns true and the next block (pos) in size order.
  // Assumes 'size' is part of the 'partition_set_'.
  virtual WP2_NO_DISCARD
  bool TryGetNextBlock(BlockSize size, Block* const block) const = 0;

  // Writes a block of a given size 'dim' in size order. The full block is
  // returned in 'block' if not nullptr.
  // 'ind' is some index to an external structure that contains information
  // about the block. It can be set to whatever if unused and it will belong to
  // an Info struct when popped in useFinal.
  // Returns false in case of error.
  virtual WP2_NO_DISCARD
  bool UseSize(BlockSize dim, uint32_t ind, Block* const block) = 0;

 protected:
  WP2Status CopyInternal(const FrontMgrNxNBase& other);

  // Finds the coordinates of the next block in lexicographic order.
  void FindNextLexico(uint32_t* const bx, uint32_t* const by) const;
  // Same as above but we know the last block that got added.
  void FindNextLexico(const Block& added_last, uint32_t* const bx,
                      uint32_t* const by) const;
  // Updates the coordinates of the next block to (x, y), and find the
  // dimensions that fit.
  void UpdateNextBlock(uint32_t x, uint32_t y);

  // Stack keeping track of the written sizes that need to be popped when the
  // rest of the block is written.
  VectorNoCtor<Info> size_stack_;
  Vector_u32 occupancy_size_;
  VectorNoCtor<Block> blocks_;
  Vector_u16 size_indices_;
  // Next max possible block.
  Block next_;
};

//------------------------------------------------------------------------------
// Chooses between a lexicographic or maximum top and left context front
// manager. For now, we only use one or the other (hence the typedef) but both
// could be used independently.

// Lexicographic front manager.
// Front manager for which the size and block order are the same.
class FrontMgrLexico : public FrontMgrNxNBase {
 public:
  void Clear() override;
  WP2Status CopyFrom(const FrontMgrLexico& other);
  WP2_NO_DISCARD
  bool TryGetNextBlock(BlockSize size, Block* const block) const override;
  void Use(const Block& block) override;
  virtual void UndoUse(const Block& block);
  WP2_NO_DISCARD
  bool UseSize(BlockSize dim, uint32_t ind, Block* const block) override;
  virtual void UndoUseSize(const Block& block);
  WP2Status Sort() override;
  bool IsFinal(const WP2::Block& block) const override;

 private:
  // Returns whether the 'Blocks()' are in the order they will be written.
  WP2Status CheckBlockList() const;
};

//------------------------------------------------------------------------------

// Front manager that maximizes the top and left context before writing the next
// block in size order.
class FrontMgrMax : public FrontMgrNxNBase {
 public:
  void Clear() override;
  WP2Status CopyFrom(const FrontMgrMax& other);
  WP2_NO_DISCARD
  bool TryGetNextBlock(BlockSize size, Block* const block) const override;
  void Use(const Block& block) override;
  WP2_NO_DISCARD
  bool UseSize(BlockSize dim, uint32_t ind, Block* const block) override;
  WP2Status Sort() override;
  bool IsFinal(const WP2::Block& block) const override;

 private:
  // Returns true if a block is considered final in the size occupancy.
  bool IsFinalSize(const WP2::Block& block) const;
  // Updates the 'next_' block.
  void UpdateNextBlock();
  // Returns whether the 'Blocks()' are in the order they will be written.
  WP2Status CheckBlockList() const;
};

// This is the FrontMgr class that should be used by default.
// TODO(vrabaud): Make sure the size order does not need anything outside the 32
//                lines of the cache for FrontMgrMax.
typedef FrontMgrMax FrontMgrDefault;

//------------------------------------------------------------------------------

// Generates all blocks lexicographically in the first WxH area, then all blocks
// in the second WxH area etc. (W and H being typically kMaxBlockSizePix).
class FrontMgrArea : public FrontMgrNxNBase {
 public:
  // BinaryPredicate used to sort blocks so that they are grouped by area.
  class Comp {
   public:
    // In kMinBlockSizePix.
    Comp(uint32_t tile_width, uint32_t area_width, uint32_t area_height);

    // Returns true if 'a' belongs to a lexicographically earlier area than 'b'
    // or if 'a' comes before 'b' in the same area.
    bool operator()(const Block& a, const Block& b) const;

   private:
    // Returns the unique index of the area whose 'block' belongs to.
    uint32_t GetAreaIndex(const Block& block) const;

    const uint32_t area_step_, area_width_, area_height_;
  };

  FrontMgrArea(uint32_t area_width, uint32_t area_height)  // In pixels.
      : area_width_(area_width / kMinBlockSizePix),
        area_height_(area_height / kMinBlockSizePix) {
    assert(area_width_ * kMinBlockSizePix == area_width);
    assert(area_height_ * kMinBlockSizePix == area_height);
  }
  // Only 'snapped=true' is allowed.
  WP2Status Init(PartitionSet partition_set, bool snapped,
                 uint32_t width, uint32_t height) override;

  void Clear() override;
  WP2Status CopyFrom(const FrontMgrArea& other);
  WP2_NO_DISCARD
  bool TryGetNextBlock(BlockSize size, Block* const block) const override;
  void Use(const Block& block) override;
  WP2_NO_DISCARD
  bool UseSize(BlockSize dim, uint32_t ind, Block* const block) override;
  WP2Status Sort() override;
  bool IsFinal(const WP2::Block& block) const override;

 protected:
  // Returns the next available block in the current area, or the first block in
  // the next area if the current area is full.
  void FindNextInArea(const Vector_u32& occupancy,
                      uint32_t* const bx, uint32_t* const by) const;

  const uint32_t area_width_, area_height_;  // In kMinBlockSizePix units.
};

//------------------------------------------------------------------------------

// FrontMgr for arbitrary block order. Each block must be manually set before
// any call to TryGetNextBlock(), UseSize() etc. and still has to respect the
// constraint that any block above it must be already used.
// Reuses the FrontMgrLexico implementation for Sort(), Use() etc.
class FrontMgrCustom : public FrontMgrLexico {
 public:
  WP2Status Init(PartitionSet partition_set, bool snapped,
                 uint32_t width, uint32_t height) override;

  void Clear() override;
  // Sets 'x, y' (in kMinBlockSizePix units) as the next block position.
  // Returns false if not possible.
  WP2_NO_DISCARD
  bool SetNextBlockPosition(uint32_t x, uint32_t y);
  WP2_NO_DISCARD
  bool TryGetNextBlock(BlockSize size, Block* const block) const override;
  WP2_NO_DISCARD
  bool UseSize(BlockSize dim, uint32_t ind, Block* const block) override;
  void UndoUseSize(const Block& block) override;

 protected:
  // If false, TryGetNextBlock() and UseSize() will return false.
  bool next_is_set_ = false;
};

//------------------------------------------------------------------------------

}  // namespace WP2

#endif  // WP2_UTILS_FRONT_MGR_H_
