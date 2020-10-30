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
// Author: Vincent Rabaud (rabaud@google.com)
//

#ifndef WP2_ENC_LOSSLESS_BACKWARD_REFERENCES_H_
#define WP2_ENC_LOSSLESS_BACKWARD_REFERENCES_H_

#include <cassert>
#include <cstdlib>
#include <memory>

#include "src/common/lossless/color_cache.h"
#include "src/common/symbols.h"
#include "src/utils/vector.h"
#include "src/wp2/base.h"
#include "src/wp2/format_constants.h"

namespace WP2L {

// -----------------------------------------------------------------------------
// PixelMode

// Class detailing what information a pixel contains.
class PixelMode {
 public:
  PixelMode() = default;

  static inline PixelMode CreateCopy(uint32_t distance, uint16_t len) {
    return {kSymbolTypeCopy, len, distance};
  }

  static inline PixelMode CreateCacheIdx(uint32_t idx, uint32_t range) {
    assert(idx < (1 << kMaxCacheBits));
    return {kSymbolTypeCacheIdx, 1, idx, range};
  }

  static inline PixelMode CreateLiteral(const int16_t* const argb) {
    return {kSymbolTypeLiteral, 1, argb};
  }

  inline SymbolType GetMode() const { return mode_; }

  // component 0,1,2,3 in order A,R,G,B.
  inline int16_t GetLiteral(int component) const {
    assert(mode_ == kSymbolTypeLiteral);
    return argb_or_distance_.argb[component];
  }

  inline uint32_t GetLength() const { return len_; }

  inline uint32_t GetCacheIdx() const {
    assert(mode_ == kSymbolTypeCacheIdx);
    assert(argb_or_distance_.cache_idx < (1U << kMaxCacheBits));
    return argb_or_distance_.cache_idx;
  }

  inline uint32_t GetCacheIdxRange() const {
    assert(mode_ == kSymbolTypeCacheIdx);
    return range_;
  }

  inline uint32_t GetDistance() const {
    assert(mode_ == kSymbolTypeCopy);
    return argb_or_distance_.distance;
  }

  inline void SetDistance(uint32_t distance) {
    assert(mode_ == kSymbolTypeCopy);
    argb_or_distance_.distance = distance;
  }

 private:
  PixelMode(SymbolType mode, uint16_t len, uint32_t distance,
            uint32_t range = 0)
      : mode_(mode), len_(len), range_(range) {
    argb_or_distance_.distance = distance;
  }
  PixelMode(SymbolType mode, uint16_t len, const int16_t* const argb)
      : mode_(mode), len_(len), range_(0) {
    for (uint32_t i = 0; i < 4; ++i) {
      argb_or_distance_.argb[i] = argb[i];
    }
  }

  SymbolType mode_;
  uint16_t len_;
  union {
    int16_t argb[4];
    uint32_t distance;
    uint32_t cache_idx;
  } argb_or_distance_;  // Or index in color cache.
  uint32_t range_;  // Range of 'argb_or_distance_' (currently used for cache).
};

// -----------------------------------------------------------------------------
// WP2LHashChain

class HashChain {
 public:
  HashChain() {
    static_assert(kMaxLengthBits + kWindowSizeBits <= 8 * sizeof(uint64_t),
                  "Wrong bits");
  }
  // Must be called first, to set size.
  // This is the maximum size of the hash_chain that can be constructed.
  // Typically this is the pixel count (width x height) for a given image.
  WP2_NO_DISCARD
  bool Allocate(uint32_t size);
  // Pre-compute the best matches for argb.
  WP2Status Fill(int speed, const int16_t* const argb, uint32_t width,
                 uint32_t height);

  inline uint32_t FindOffset(const uint32_t base_position) const {
    return offset_length_[base_position] >> kMaxLengthBits;
  }

  inline uint32_t FindLength(const uint32_t base_position) const {
    return offset_length_[base_position] & ((1U << kMaxLengthBits) - 1);
  }

  inline void FindCopy(uint32_t base_position, uint32_t* const offset_ptr,
                       uint32_t* const length_ptr) const {
    *offset_ptr = FindOffset(base_position);
    *length_ptr = FindLength(base_position);
  }

  // The kWindowSizeBits most significant bits contain the offset at which the
  // best match is found. The lower kMaxLengthBits bits contain the length of
  // the match.
  WP2::Vector_u64 offset_length_;

  // We want the max value to be attainable and stored in kMaxLengthBits bits.
  static constexpr uint8_t kMaxLengthBits = 1 + WP2Log2Floor_k(kMaxLZ77Length);

 private:
  static constexpr uint8_t kHashBits = 18;
  static constexpr uint32_t kHashSize = (1u << kHashBits);

  static constexpr uint8_t kWindowSizeBits =
      1 + WP2Log2Floor_k(kMaxLZ77Distance);
};

// -----------------------------------------------------------------------------
// WP2LBackwardRefs (block-based backward-references storage)

// maximum number of reference blocks the image will be segmented into
#define MAX_REFS_BLOCK_PER_IMAGE 16

struct PixOrCopyBlock {
  std::shared_ptr<PixOrCopyBlock> next;  // next block (or nullptr)
  WP2::VectorNoCtor<PixelMode> modes;    // pixel modes
  int size;                              // currently used size
};
class RefsCursor;
class BackwardRefsPool;

// Container for blocks chain
class BackwardRefs {
 public:
  BackwardRefs()
      : block_size_(0),
        refs_(nullptr),
        tail_(nullptr),
        free_blocks_(nullptr),
        last_block_(nullptr),
        pool_(nullptr) {}
  ~BackwardRefs() { Reset(); }
  WP2Status CopyFrom(const BackwardRefs& refs);
  // Releases memory for backward references.
  void Reset();
  void Clear();
  WP2Status CursorAdd(const PixelMode& v);
  // Returns true if the path covers "num_pixels" pixels.
  bool IsValid(uint32_t num_pixels) const;

  friend RefsCursor;  // To access refs_
  friend BackwardRefsPool;  // To access Init

 private:
  // Initialize the object. 'block_size' is the common block size to store
  // references (typically, width * height / MAX_REFS_BLOCK_PER_IMAGE).
  void Init(uint32_t block_size, BackwardRefsPool* const pool);
  WP2Status AddBlock();
  static void FreeFromPool(BackwardRefs* const refs);

  // minimum block size for backward references
  static constexpr uint32_t kMinBlockSize = 256;

  int block_size_;                         // common block-size
  std::shared_ptr<PixOrCopyBlock> refs_;   // list of currently used blocks
  std::shared_ptr<PixOrCopyBlock>* tail_;  // for list recycling
  std::shared_ptr<PixOrCopyBlock> free_blocks_;  // free-list
  // used for adding new refs (internal)
  std::shared_ptr<PixOrCopyBlock> last_block_;
  BackwardRefsPool* pool_;  // pool to which it belongs
};

// Class holding several allocated BackwardRefs.
// As we compare several BackwardRefs but only keep the best one, this
// class is useful not to keep reallocating BackwardRefs.
class BackwardRefsPool {
 public:
  // Initialize the BackwardRefs.
  void Init(uint32_t num_pixels) {
    // We round the block size up, so we're guaranteed to have
    // at most MAX_REFS_BLOCK_PER_IMAGE blocks used:
    const uint32_t block_size = (num_pixels - 1) / MAX_REFS_BLOCK_PER_IMAGE + 1;
    for (size_t i = 0; i < kNumBackwardRefs; ++i) {
      refs_[i].Reset();
      refs_[i].Init(block_size, this);
      is_used_[i] = false;
    }
  }
  typedef std::unique_ptr<BackwardRefs, void (*)(BackwardRefs*)> RefsPtr;
  static RefsPtr GetEmptyBackwardRefs() {
    return RefsPtr(nullptr, BackwardRefs::FreeFromPool);
  }
  RefsPtr GetFreeBackwardRefs() {
    // Send back the first BackwardRefs that is not used.
    for (size_t i = 0; i < kNumBackwardRefs; ++i) {
      if (!is_used_[i]) {
        is_used_[i] = true;
        return RefsPtr(&refs_[i], BackwardRefs::FreeFromPool);
      }
    }
    assert(false);
    return RefsPtr(nullptr, BackwardRefs::FreeFromPool);
  }
  // Reset the unused backward refs.
  void Reset() {
    for (size_t i = 0; i < kNumBackwardRefs; ++i) {
      if (!is_used_[i]) refs_[i].Reset();
    }
  }
  friend BackwardRefs;

 private:
  void Release(BackwardRefs* const refs) {
    for (size_t i = 0; i < kNumBackwardRefs; ++i) {
      if (&refs_[i] == refs) {
        assert(is_used_[i]);
        is_used_[i] = false;
        return;
      }
    }
    assert(false);
  }
  // 5 is enough in the code for now: 1 as the chosen backward ref between
  // iterations of the cruncher, 1 best per iteration, and 2 tmps in addition in
  // GetBackwardReferences
  static constexpr size_t kNumBackwardRefs = 5;
  BackwardRefs refs_[kNumBackwardRefs];
  bool is_used_[kNumBackwardRefs];
};

// Cursor for iterating on references content
class RefsCursor {
 public:
  // Positions the cursor at the beginning of the references list.
  explicit RefsCursor(const BackwardRefs& refs);
  // Returns true if cursor is pointing at a valid position.
  inline bool Ok() const { return (cur_pos_ != nullptr); }
  inline void Next() {
    assert(Ok());
    if (++cur_pos_ == last_pos_) NextBlock();
  }

  PixelMode* cur_pos_;  // current position
 private:
  // Move to next block of references.
  void NextBlock();

  std::shared_ptr<PixOrCopyBlock> cur_block_;  // current block in the refs list
  const PixelMode* last_pos_;  // sentinel for switching to next block
};

// -----------------------------------------------------------------------------
// Main entry points

enum LZ77Type {
  kLZ77Standard = 1,
  kLZ77RLE = 2,
  // LZ77 where matches are forced to happen within a given distance cost.
  kLZ77Box = 4,
  kLZ77None = 8
};

// Evaluates best possible backward references for specified speed.
// The input cache_bits to 'GetBackwardReferences' sets the maximum cache
// bits to use (passing 0 implies disabling the local color cache).
// The optimal cache bits is evaluated and set for the *cache_bits parameter.
// The return value is the pointer to the best of the two backward refs viz,
// refs[0] or refs[1].
WP2Status GetBackwardReferences(
    uint32_t width, uint32_t height, const int16_t* const argb, int speed,
    int lz77_types_to_try, uint32_t cache_bits_max, const HashChain& hash_chain,
    const LosslessSymbolsInfo& symbols_info, CacheConfig* const cache_config,
    BackwardRefsPool* const ref_pool, BackwardRefsPool::RefsPtr* const refs);
}  // namespace WP2L

#endif  // WP2_ENC_LOSSLESS_BACKWARD_REFERENCES_H_
