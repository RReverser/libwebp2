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
// Author: Jyrki Alakuijala (jyrki@google.com)
//

#include "src/enc/lossless/backward_references_enc.h"

#include "src/common/lossless/color_cache.h"
#include "src/dsp/math.h"
#include "src/enc/lossless/histogram_enc.h"
#include "src/enc/lossless/losslessi_enc.h"
#include "src/enc/symbols_enc.h"

namespace WP2L {

// Minimum number of pixels for which it is cheaper to encode a
// distance + length instead of each pixel as a literal.
#define MIN_LENGTH 4

// -----------------------------------------------------------------------------

static const uint8_t plane_to_code_lut[128] = {
  96,   73,  55,  39,  23,  13,   5,  1, 255, 255, 255, 255, 255, 255, 255, 255,
  101,  78,  58,  42,  26,  16,   8,  2,   0,   3,  9,   17,  27,  43,  59,  79,
  102,  86,  62,  46,  32,  20,  10,  6,   4,   7,  11,  21,  33,  47,  63,  87,
  105,  90,  70,  52,  37,  28,  18,  14, 12,  15,  19,  29,  38,  53,  71,  91,
  110,  99,  82,  66,  48,  35,  30,  24, 22,  25,  31,  36,  49,  67,  83, 100,
  115, 108,  94,  76,  64,  50,  44,  40, 34,  41,  45,  51,  65,  77,  95, 109,
  118, 113, 103,  92,  80,  68,  60,  56, 54,  57,  61,  69,  81,  93, 104, 114,
  119, 116, 111, 106,  97,  88,  84,  74, 72,  75,  85,  89,  98, 107, 112, 117
};

uint32_t DistanceToPlaneCode(uint32_t width, uint32_t dist) {
  const uint32_t yoffset = dist / width;
  const uint32_t xoffset = dist - yoffset * width;
  if (xoffset <= 8 && yoffset < 8) {
    return plane_to_code_lut[yoffset * 16 + 8 - xoffset] + 1;
  } else if (xoffset + 8 > width && yoffset < 7) {
    return plane_to_code_lut[(yoffset + 1) * 16 + 8 + (width - xoffset)] + 1;
  }
  return dist + kCodeToPlaneCodes;
}

// Returns the exact index where array1 and array2 are different. For an index
// inferior or equal to best_len_match, the return value just has to be strictly
// inferior to best_len_match. The current behavior is to return 0 if this index
// is best_len_match, and the index itself otherwise.
// If no two elements are the same, it returns max_limit.
static inline int FindMatchLength(const int16_t* const array1,
                                  const int16_t* const array2,
                                  uint32_t best_len_match, uint32_t max_limit) {
  // Before 'expensive' linear match, check if the two arrays match at the
  // current best length index.
  if (!ColorEq(&array1[4 * best_len_match], &array2[4 * best_len_match])) {
    return 0;
  }

  const auto m = std::mismatch(array1, array1 + 4 * max_limit, array2);
  return (m.first - array1) / 4;
}

// -----------------------------------------------------------------------------
//  BackwardRefs

void BackwardRefs::Clear() {
  if (tail_ != nullptr) {
    *tail_ = free_blocks_;  // recycle all blocks at once
  }
  free_blocks_ = refs_;
  tail_ = &refs_;
  last_block_ = nullptr;
  refs_ = nullptr;
}

WP2Status BackwardRefs::CopyFrom(const BackwardRefs& refs) {
  Clear();
  const PixOrCopyBlock* refs_in = refs.refs_.get();
  // Copy whole blocks over.
  while (refs_in != nullptr) {
    WP2_CHECK_STATUS(AddBlock());
    std::copy(refs_in->modes.data(), refs_in->modes.data() + refs_in->size,
              last_block_->modes.data());
    last_block_->size = refs_in->size;
    refs_in = refs_in->next.get();
  }
  return WP2_STATUS_OK;
}

void BackwardRefs::Reset() {
  Clear();
  free_blocks_.reset();
}

void BackwardRefs::Init(uint32_t block_size, BackwardRefsPool* const pool) {
  refs_ = nullptr;
  tail_ = nullptr;
  free_blocks_ = nullptr;
  last_block_ = nullptr;
  tail_ = &refs_;
  block_size_ = std::max(block_size, +kMinBlockSize);
  pool_ = pool;
}

RefsCursor::RefsCursor(const BackwardRefs& refs) {
  cur_block_ = refs.refs_;
  if (refs.refs_ != nullptr) {
    cur_pos_ = cur_block_->modes.data();
    last_pos_ = cur_pos_ + cur_block_->size;
  } else {
    cur_pos_ = nullptr;
    last_pos_ = nullptr;
  }
}

void RefsCursor::NextBlock() {
  auto b = cur_block_->next;
  cur_pos_ = (b == nullptr) ? nullptr : b->modes.data();
  last_pos_ = (b == nullptr) ? nullptr : b->modes.data() + b->size;
  cur_block_ = b;
}

// Create a new block, either from the free list or allocated
WP2Status BackwardRefs::AddBlock() {
  std::shared_ptr<PixOrCopyBlock> b = free_blocks_;
  if (b == nullptr) {  // allocate new memory chunk
    b = std::make_shared<PixOrCopyBlock>();
    if (!b->modes.resize(block_size_)) return WP2_STATUS_OUT_OF_MEMORY;
  } else {  // recycle from free-list
    free_blocks_ = b->next;
  }
  *tail_ = b;
  tail_ = &b->next;
  last_block_ = b;
  b->next = nullptr;
  b->size = 0;
  return WP2_STATUS_OK;
}

void BackwardRefs::FreeFromPool(BackwardRefs* const refs) {
  if (refs != nullptr && refs->pool_ != nullptr) refs->pool_->Release(refs);
}

WP2Status BackwardRefs::CursorAdd(const PixelMode& v) {
  auto b = last_block_;
  if (b == nullptr || b->size == block_size_) {
    WP2_CHECK_STATUS(AddBlock());
    b = last_block_;
  }
  b->modes[b->size++] = v;
  return WP2_STATUS_OK;
}

bool BackwardRefs::IsValid(uint32_t num_pixels) const {
  // Check whether the cursor goes over the whole image.
  RefsCursor c(*this);
  uint32_t tot = 0;
  while (c.Ok()) {
    tot += c.cur_pos_->GetLength();
    c.Next();
  }
  return (tot == num_pixels);
}

// -----------------------------------------------------------------------------
// Hash chains

bool HashChain::Allocate(uint32_t size) {
  assert(size > 0);
  if (!offset_length_.resize(size)) return false;

  return true;
}

// -----------------------------------------------------------------------------

#define HASH_MULTIPLIER_HI (0xc6a4a793ULL)
#define HASH_MULTIPLIER_LO (0x5bd1e996ULL)

static inline uint32_t GetPixPairHash64(const int16_t argb[8],
                                        const uint8_t hash_bits) {
  const uint32_t a = ((uint16_t)argb[0] << 24) | ((uint16_t)argb[1] << 16) |
                     ((uint16_t)argb[2] <<  8) | ((uint16_t)argb[3] <<  0);
  const uint32_t b = ((uint16_t)argb[4] << 24) | ((uint16_t)argb[5] << 16) |
                     ((uint16_t)argb[6] <<  8) | ((uint16_t)argb[7] <<  0);
  uint32_t key;
  key = (b * HASH_MULTIPLIER_HI) & 0xffffffffu;
  key += (a * HASH_MULTIPLIER_LO) & 0xffffffffu;
  key = key >> (32 - hash_bits);
  return key;
}

// Returns the maximum number of hash chain lookups to do for a
// given compression speed. Return value in range [45, 126].
// TODO(vrabaud) re-adjust those values and have a different logic (the number
// of iteration does not necessarily matter as much as the number of lines)
static int GetMaxItersForSpeed(int speed) {
  return (speed == 3) ? 8 + 20 * 20 / 128 :
         (speed == 5) ? 8 + 75 * 75 / 128 : 45 + speed * speed;
}

static uint32_t GetWindowSizeForHashChain(int speed, uint32_t width,
                                          const uint32_t window_size) {
  const uint32_t max_window_size = (speed >= 7) ? window_size
                                 : (speed >= 5) ? (width << 8)
                                 : (speed >= 3) ? (width << 6)
                                 : (width << 4);
  assert(width > 0);
  return std::min(max_window_size, window_size);
}

WP2Status HashChain::Fill(int speed, const int16_t* const argb, uint32_t width,
                          uint32_t height) {
  const int size = width * height;
  const int iter_max = GetMaxItersForSpeed(speed);
  const uint32_t window_size =
      GetWindowSizeForHashChain(speed, width, kMaxLZ77Distance);
  int pos;
  // Temporarily use the offset_length_ as a hash chain.
  auto chain = reinterpret_cast<int32_t*>(offset_length_.data());
  assert(size > 0);
  assert(!offset_length_.empty());

  if (size <= 2) {
    offset_length_[0] = offset_length_[size - 1] = 0;
    return WP2_STATUS_OK;
  }

  WP2::Vector_s32 hash_to_first_index;
  WP2_CHECK_ALLOC_OK(hash_to_first_index.resize(kHashSize));

  // Set the int32_t array to -1.
  memset(hash_to_first_index.data(), 0xff,
         kHashSize * sizeof(*hash_to_first_index.data()));
  // Fill the chain linking pixels with the same hash.
  bool argb_comp = ColorEq(&argb[0], &argb[4]);
  for (pos = 0; pos < size - 2;) {
    uint32_t hash_code;
    const bool argb_comp_next =
        ColorEq(&argb[4 * (pos + 1)], &argb[4 * (pos + 2)]);
    if (argb_comp && argb_comp_next) {
      // Consecutive pixels with the same color will share the same hash.
      // We therefore use a different hash: the color and its repetition
      // length.
      int16_t tmp[8] = {0};
      uint32_t len = 1;
      ColorCopy(&argb[4 * pos], &tmp[0]);
      // Figure out how far the pixels are the same.
      // The last pixel has a different 64 bit hash, as its next pixel does
      // not have the same color, so we just need to get to the last pixel equal
      // to its follower.
      while (pos + (int)len + 2 < size &&
             ColorEq(&argb[4 * (pos + len + 2)], &argb[4 * pos])) {
        ++len;
      }
      if (len > kMaxLZ77Length) {
        // Skip the pixels that match for distance=1 and length>kMaxLZ77Length
        // because they are linked to their predecessor and we automatically
        // check that in the main for loop below. Skipping means setting no
        // predecessor in the chain, hence -1.
        memset(chain + pos, 0xff, (len - kMaxLZ77Length) * sizeof(*chain));
        pos += len - kMaxLZ77Length;
        len = kMaxLZ77Length;
      }
      // Process the rest of the hash chain.
      while (len) {
        tmp[7] = len--;
        hash_code = GetPixPairHash64(tmp, kHashBits);
        chain[pos] = hash_to_first_index[hash_code];
        hash_to_first_index[hash_code] = pos++;
      }
      argb_comp = false;
    } else {
      // Just move one pixel forward.
      hash_code = GetPixPairHash64(&argb[4 * pos], kHashBits);
      chain[pos] = hash_to_first_index[hash_code];
      hash_to_first_index[hash_code] = pos++;
      argb_comp = argb_comp_next;
    }
  }
  // Process the penultimate pixel.
  chain[pos] =
      hash_to_first_index[GetPixPairHash64(&argb[4 * pos], kHashBits)];

  hash_to_first_index.reset();

  // Find the best match interval at each pixel, defined by an offset to the
  // pixel and a length. The right-most pixel cannot match anything to the right
  // (hence a best length of 0) and the left-most pixel nothing to the left
  // (hence an offset of 0).
  assert(size > 2);
  offset_length_[0] = offset_length_[size - 1] = 0;
  for (uint32_t base_position = size - 2; base_position > 0;) {
    const uint32_t max_len = std::min(size - 1 - base_position, kMaxLZ77Length);
    const int16_t* const argb_start = &argb[4 * base_position];
    int iter = iter_max;
    uint32_t best_length = 0;
    uint32_t best_distance = 0;
    int16_t best_argb[4];
    const int min_pos =
        (base_position > window_size) ? base_position - window_size : 0;
    const uint32_t length_max = std::min(max_len, 256u);
    uint32_t max_base_position;

    pos = chain[base_position];
    if (speed > 0) {
      uint32_t curr_length;
      // Heuristic: use the comparison with the above line as an initialization.
      if (base_position >= width) {
        curr_length = FindMatchLength(argb_start - 4 * width, argb_start,
                                      best_length, max_len);
        if (curr_length > best_length) {
          best_length = curr_length;
          best_distance = width;
        }
        --iter;
      }
      // Heuristic: compare to the previous pixel.
      curr_length =
          FindMatchLength(argb_start - 4, argb_start, best_length, max_len);
      if (curr_length > best_length) {
        best_length = curr_length;
        best_distance = 1;
      }
      --iter;
      // Skip the for loop if we already have the maximum.
      if (best_length == kMaxLZ77Length) pos = min_pos - 1;
    }
    ColorCopy(&argb_start[4 * best_length], best_argb);

    for (; pos >= min_pos && --iter; pos = chain[pos]) {
      uint32_t curr_length;
      assert(base_position > (uint32_t)pos);

      if (!ColorEq(&argb[4 * (pos + best_length)], best_argb)) continue;

      const auto start = &argb[4 * pos];
      const auto p = std::mismatch(start, start + 4 * max_len, argb_start);
      curr_length = (p.first - start) / 4;
      if (best_length < curr_length) {
        best_length = curr_length;
        best_distance = base_position - pos;
        ColorCopy(&argb_start[4 * best_length], best_argb);
        // Stop if we have reached a good enough length.
        if (best_length >= length_max) break;
      }
    }
    // We have the best match but in case the two intervals continue matching
    // to the left, we have the best matches for the left-extended pixels.
    max_base_position = base_position;
    while (true) {
      assert(best_length <= kMaxLZ77Length);
      assert(best_distance <= kMaxLZ77Distance);
      offset_length_[base_position] =
          ((uint64_t)best_distance << kMaxLengthBits) | (uint64_t)best_length;
      --base_position;
      // Stop if we don't have a match or if we are out of bounds.
      if (best_distance == 0 || base_position == 0) break;
      // Stop if we cannot extend the matching intervals to the left.
      if (base_position < best_distance ||
          !ColorEq(&argb[4 * (base_position - best_distance)],
                   &argb[4 * base_position])) {
        break;
      }
      // Stop if we are matching at its limit because there could be a closer
      // matching interval with the same maximum length. Then again, if the
      // matching interval is as close as possible (best_distance == 1), we will
      // never find anything better so let's continue.
      if (best_length == kMaxLZ77Length && best_distance != 1 &&
          base_position + kMaxLZ77Length < max_base_position) {
        break;
      }
      if (best_length < kMaxLZ77Length) {
        ++best_length;
        max_base_position = base_position;
      }
    }
  }
  return WP2_STATUS_OK;
}

static inline WP2Status AddSingleLiteral(const int16_t* const pixel,
                                         bool use_color_cache,
                                         ColorCache* const color_cache,
                                         BackwardRefs* const refs) {
  PixelMode v;
  if (use_color_cache) {
    uint32_t cache_index;
    const uint32_t index_range = color_cache->IndexRange();
    if (color_cache->Insert(pixel, &cache_index)) {
      v = PixelMode::CreateLiteral(pixel);
    } else {
      v = PixelMode::CreateCacheIdx(cache_index, index_range);
    }
  } else {
    v = PixelMode::CreateLiteral(pixel);
  }
  WP2_CHECK_STATUS(refs->CursorAdd(v));
  return WP2_STATUS_OK;
}

static WP2Status BackwardReferencesRle(uint32_t width, uint32_t height,
                                       const int16_t* const argb,
                                       const CacheConfig& cache_config,
                                       BackwardRefs* const refs) {
  const uint32_t num_pixels = width * height;

  const bool use_color_cache = (cache_config.type != CacheType::kNone);
  ColorCachePtr color_cache;
  WP2_CHECK_STATUS(color_cache.Init(cache_config));

  refs->Clear();
  // Add first pixel as literal.
  WP2_CHECK_STATUS(
      AddSingleLiteral(&argb[0], use_color_cache, color_cache.get(), refs));
  for (uint32_t i = 1; i < num_pixels;) {
    const uint32_t max_len = std::min(num_pixels - i, kMaxLZ77Length);
    const uint32_t rle_len =
        FindMatchLength(&argb[4 * i], &argb[4 * (i - 1)], 0, max_len);
    const uint32_t prev_row_len =
        (i < width)
            ? 0
            : FindMatchLength(&argb[4 * i], &argb[4 * (i - width)], 0, max_len);
    if (rle_len >= prev_row_len && rle_len >= MIN_LENGTH) {
      WP2_CHECK_STATUS(refs->CursorAdd(PixelMode::CreateCopy(1, rle_len)));
      // We don't need to update the color cache here since it is always the
      // same pixel being copied, and that does not change the color cache
      // state.
      i += rle_len;
    } else if (prev_row_len >= MIN_LENGTH) {
      WP2_CHECK_STATUS(
          refs->CursorAdd(PixelMode::CreateCopy(width, prev_row_len)));
      if (use_color_cache) {
        for (uint32_t k = 0; k < prev_row_len; ++k) {
          color_cache->Insert(&argb[4 * (i + k)], /*index_ptr=*/nullptr);
        }
      }
      i += prev_row_len;
    } else {
      WP2_CHECK_STATUS(AddSingleLiteral(&argb[4 * i], use_color_cache,
                                        color_cache.get(), refs));
      i++;
    }
  }
  return WP2_STATUS_OK;
}

static WP2Status BackwardReferencesNone(uint32_t num_pixels,
                                        const int16_t* const argb,
                                        const CacheConfig& cache_config,
                                        BackwardRefs* const refs) {
  const bool use_color_cache = (cache_config.type != CacheType::kNone);
  ColorCachePtr color_cache;
  WP2_CHECK_STATUS(color_cache.Init(cache_config));

  refs->Clear();
  for (uint32_t i = 0; i < num_pixels; ++i) {
    WP2_CHECK_STATUS(AddSingleLiteral(&argb[4 * i], use_color_cache,
                                      color_cache.get(), refs));
  }
  return WP2_STATUS_OK;
}

static WP2Status BackwardReferencesLz77(uint32_t width, uint32_t height,
                                        const int16_t* const argb,
                                        const CacheConfig& cache_config,
                                        const HashChain& hash_chain,
                                        BackwardRefs* const refs) {
  const uint32_t num_pixels = width * height;

  const bool use_color_cache = (cache_config.type != CacheType::kNone);
  ColorCachePtr color_cache;
  WP2_CHECK_STATUS(color_cache.Init(cache_config));

  refs->Clear();
  for (uint32_t i = 0, i_last_check = 0; i < num_pixels;) {
    // Alternative#1: Code the pixels starting at 'i' using backward reference.
    uint32_t offset = 0, len = 0;
    hash_chain.FindCopy(i, &offset, &len);
    if (len >= MIN_LENGTH) {
      const int len_ini = len;
      uint32_t max_reach = 0;
      const uint32_t j_max =
          (i + len_ini >= num_pixels) ? num_pixels - 1 : i + len_ini;
      // Only start from what we have not checked already.
      i_last_check = (i > i_last_check) ? i : i_last_check;
      // We know the best match for the current pixel but we try to find the
      // best matches for the current pixel AND the next one combined.
      // The naive method would use the intervals:
      // [i,i+len) + [i+len, length of best match at i+len)
      // while we check if we can use:
      // [i,j) (where j<=i+len) + [j, length of best match at j)
      for (uint32_t j = i_last_check + 1; j <= j_max; ++j) {
        const uint32_t len_j = hash_chain.FindLength(j);
        const uint32_t reach =
            j + (len_j >= MIN_LENGTH ? len_j : 1);  // 1 for single literal.
        if (reach > max_reach) {
          len = j - i;
          max_reach = reach;
          if (max_reach >= num_pixels) break;
        }
      }
    } else {
      len = 1;
    }
    // Go with literal or backward reference.
    assert(len > 0);
    if (len == 1) {
      WP2_CHECK_STATUS(AddSingleLiteral(&argb[4 * i], use_color_cache,
                                        color_cache.get(), refs));
    } else {
      WP2_CHECK_STATUS(refs->CursorAdd(PixelMode::CreateCopy(offset, len)));
      if (use_color_cache) {
        for (uint32_t j = i; j < i + len; ++j)
          color_cache->Insert(&argb[4 * j], /*index_ptr=*/nullptr);
      }
    }
    i += len;
  }

  return WP2_STATUS_OK;
}

// Compute an LZ77 by forcing matches to happen within a given distance cost.
// We therefore limit the algorithm to the lowest 32 values in the PlaneCode
// definition.
#define WINDOW_OFFSETS_SIZE_MAX 32
static WP2Status BackwardReferencesLz77Box(uint32_t width, uint32_t height,
                                           const int16_t* const argb,
                                           const CacheConfig& cache_config,
                                           const HashChain& hash_chain_best,
                                           HashChain* hash_chain,
                                           BackwardRefs* const refs) {
  const uint32_t num_pixels = width * height;
  int window_offsets[WINDOW_OFFSETS_SIZE_MAX] = {0};
  int window_offsets_new[WINDOW_OFFSETS_SIZE_MAX] = {0};
  uint32_t window_offsets_size = 0;
  int window_offsets_new_size = 0;
  uint32_t best_offset_prev = 0, best_length_prev = 0;
  WP2::Vector_u16 counts_ini;
  WP2_CHECK_ALLOC_OK(counts_ini.resize(width * height));

  // counts[i] counts how many times a pixel is repeated starting at position i.
  {
    int i = (int)num_pixels - 2;
    uint16_t* counts = counts_ini.data() + i;
    counts[1] = 1;
    for (; i >= 0; --i, --counts) {
      if (ColorEq(&argb[4 * i], &argb[4 * (i + 1)])) {
        // Max out the counts to kMaxLength.
        counts[0] = counts[1] + (counts[1] != kMaxLZ77Length);
      } else {
        counts[0] = 1;
      }
    }
  }

  // Figure out the window offsets around a pixel. They are stored in a
  // spiraling order around the pixel as defined by DistanceToPlaneCode.
  {
    int x, y;
    for (y = 0; y <= 6; ++y) {
      for (x = -6; x <= 6; ++x) {
        const int offset = y * width + x;
        int plane_code;
        // Ignore offsets that bring us after the pixel.
        if (offset <= 0) continue;
        plane_code = DistanceToPlaneCode(width, offset) - 1;
        if (plane_code >= WINDOW_OFFSETS_SIZE_MAX) continue;
        window_offsets[plane_code] = offset;
      }
    }
    // For narrow images, not all plane codes are reached, so remove those.
    for (const auto& offset : window_offsets) {
      if (offset == 0) continue;
      window_offsets[window_offsets_size++] = offset;
    }
    // Given a pixel P, find the offsets that reach pixels unreachable from P-1
    // with any of the offsets in window_offsets[].
    for (uint32_t i = 0; i < window_offsets_size; ++i) {
      bool is_reachable = false;
      for (uint32_t j = 0; j < window_offsets_size && !is_reachable; ++j) {
        is_reachable |= (window_offsets[i] == window_offsets[j] + 1);
      }
      if (!is_reachable) {
        window_offsets_new[window_offsets_new_size] = window_offsets[i];
        ++window_offsets_new_size;
      }
    }
  }

  hash_chain->offset_length_[0] = 0;
  for (uint32_t i = 1; i < num_pixels; ++i) {
    uint32_t best_length = hash_chain_best.FindLength(i);
    int best_offset;
    int do_compute = 1;

    if (best_length >= kMaxLZ77Length) {
      // Do not recompute the best match if we already have a maximal one in the
      // window.
      best_offset = hash_chain_best.FindOffset(i);
      for (uint32_t ind = 0; ind < window_offsets_size; ++ind) {
        if (best_offset == window_offsets[ind]) {
          do_compute = 0;
          break;
        }
      }
    }
    if (do_compute) {
      // Figure out if we should use the offset/length from the previous pixel
      // as an initial guess and therefore only inspect the offsets in
      // window_offsets_new[].
      const bool use_prev =
          (best_length_prev > 1) && (best_length_prev < kMaxLZ77Length);
      const uint32_t num_ind =
          use_prev ? window_offsets_new_size : window_offsets_size;
      best_length = use_prev ? best_length_prev - 1 : 0;
      best_offset = use_prev ? best_offset_prev : 0;
      // Find the longest match in a window around the pixel.
      for (uint32_t ind = 0; ind < num_ind; ++ind) {
        uint32_t curr_length = 0;
        uint32_t j = i;
        int j_offset =
            use_prev ? i - window_offsets_new[ind] : i - window_offsets[ind];
        if (j_offset < 0 || !ColorEq(&argb[4 * j_offset], &argb[4 * i])) {
          continue;
        }
        // The longest match is the sum of how many times each pixel is
        // repeated.
        do {
          const int counts_j_offset = counts_ini[j_offset];
          const int counts_j = counts_ini[j];
          if (counts_j_offset != counts_j) {
            curr_length +=
                (counts_j_offset < counts_j) ? counts_j_offset : counts_j;
            break;
          }
          // The same color is repeated counts_pos times at j_offset and j.
          curr_length += counts_j_offset;
          j_offset += counts_j_offset;
          j += counts_j_offset;
        } while (curr_length <= kMaxLZ77Length && j < num_pixels &&
                 ColorEq(&argb[4 * j_offset], &argb[4 * j]));
        if (best_length < curr_length) {
          best_offset =
              use_prev ? window_offsets_new[ind] : window_offsets[ind];
          if (curr_length >= kMaxLZ77Length) {
            best_length = kMaxLZ77Length;
            break;
          } else {
            best_length = curr_length;
          }
        }
      }
    }

    assert(i + best_length <= num_pixels);
    assert(best_length <= kMaxLZ77Length);
    if (best_length <= MIN_LENGTH) {
      hash_chain->offset_length_[i] = 0;
      best_offset_prev = 0;
      best_length_prev = 0;
    } else {
      hash_chain->offset_length_[i] =
          ((uint64_t)best_offset << hash_chain_best.kMaxLengthBits) |
          (uint64_t)best_length;
      best_offset_prev = best_offset;
      best_length_prev = best_length;
    }
  }
  hash_chain->offset_length_[0] = 0;

  return BackwardReferencesLz77(width, height, argb, cache_config, *hash_chain,
                                refs);
}

// -----------------------------------------------------------------------------

WP2Status BackwardRefsWithLocalCache(const int16_t* const argb,
                                     const CacheConfig& cache_config,
                                     BackwardRefs* const refs);

// Evaluate optimal cache bits for the local color cache.
// The input *best_cache_bits sets the maximum cache bits to use (passing 0
// implies disabling the local color cache). The local color cache is also
// disabled for the lower (<= 3) speed.
static WP2Status CalculateBestCacheConfig(
    uint32_t width, const int16_t* const argb, uint32_t num_pixels, int speed,
    const BackwardRefs& refs_in, BackwardRefsPool* const ref_pool,
    uint32_t cache_bits_max, const LosslessSymbolsInfo& symbols_info,
    CacheConfig* const cache_config) {
  if (speed <= 3 || cache_bits_max == 0) {
    cache_config->type = CacheType::kNone;
    cache_config->cache_bits = 0;
    return WP2_STATUS_OK;
  }

  assert(cache_bits_max > 0 && cache_bits_max <= kMaxCacheBits);

  // Allocate data.
  float cost_min = std::numeric_limits<float>::max();
  auto refs_tmp = ref_pool->GetFreeBackwardRefs();
  for (uint32_t i = 0; i <= cache_bits_max; ++i) {
    // Compute the symbol statistics.
    LosslessSymbolsInfo info;
    WP2_CHECK_STATUS(info.CopyFrom(symbols_info));

    // TODO(maryla): try other types.
    const CacheConfig cache_config_tmp = {
        (i > 0) ? CacheType::kHash : CacheType::kNone, i};

    // Update the refs to deal with a cache of size 'i' bits.
    info.SetCacheRange(GetColorCacheRange(cache_config_tmp));
    if (i > 0) {
      WP2_CHECK_STATUS(refs_tmp->CopyFrom(refs_in));
      WP2_CHECK_STATUS(
          BackwardRefsWithLocalCache(argb, cache_config_tmp, refs_tmp.get()));
    }

    const BackwardRefs& refs = (i == 0) ? refs_in : *refs_tmp;

    WP2::ANSEncNoop enc_noop;
    WP2::SymbolRecorder symbol_recorder;
    WP2_CHECK_STATUS(symbol_recorder.Allocate(info, /*num_records=*/0));
    WP2_CHECK_STATUS(StorePixels(width, refs, &enc_noop, &symbol_recorder));

    // Get the header cost in bits.
    WP2::ANSEncCounter enc_counter;
    WP2::SymbolWriter symbol_writer;
    WP2_CHECK_STATUS(symbol_writer.Init(info));
    WP2_CHECK_STATUS(symbol_writer.Allocate());
    WP2::ANSDictionaries dicts;
    WP2_CHECK_STATUS(WriteHeaders(symbol_recorder, info, num_pixels,
                                  &enc_counter, &dicts, &symbol_writer));
    const float cost_header = enc_counter.GetCost();

    // Get the storage cost by going through a proper write.
    // TODO(vrabaud) this is slow: we cannot initialize a SymbolCounter from the
    //               SymbolWriter because stats are lost.
    WP2::ANSEnc enc;
    WP2_CHECK_STATUS(StorePixels(width, refs, &enc, &symbol_writer));
    const float cost_main = enc.GetCost(dicts);
    const float cost = cost_header + cost_main;

    if (cost < cost_min) {
      *cache_config = cache_config_tmp;
      cost_min = cost;
    }
  }

  return WP2_STATUS_OK;
}

// Update (in-place) backward references for specified cache_bits.
WP2Status BackwardRefsWithLocalCache(const int16_t* const argb,
                                     const CacheConfig& cache_config,
                                     BackwardRefs* const refs) {
  ColorCachePtr color_cache;
  WP2_CHECK_STATUS(color_cache.Init(cache_config));

  uint32_t pixel_index = 0;
  RefsCursor c(*refs);
  while (c.Ok()) {
    PixelMode* const v = c.cur_pos_;
    switch (v->GetMode()) {
      case kSymbolTypeLiteral: {
        uint32_t ix;
        if (color_cache->Contains(&argb[4 * pixel_index], &ix)) {
          // color_cache contains argb_literal
          *v = PixelMode::CreateCacheIdx(ix, color_cache->IndexRange());
        }
        color_cache->Insert(&argb[4 * pixel_index], /*index_ptr=*/nullptr);
        ++pixel_index;
        break;
      }
      case kSymbolTypeCopy: {
        for (uint32_t k = 0; k < v->GetLength(); ++k) {
          color_cache->Insert(&argb[4 * pixel_index++], /*index_ptr=*/nullptr);
        }
        break;
      }
      default:
        // refs was created without local cache, so it can not have cache
        // indexes.
        assert(false);
    }
    c.Next();
  }
  return WP2_STATUS_OK;
}

static WP2Status GetBackwardReferencesLowEffort(
    uint32_t width, uint32_t height, const int16_t* const argb,
    const HashChain& hash_chain, BackwardRefsPool* const ref_pool,
    BackwardRefsPool::RefsPtr* const refs) {
  refs->reset();
  auto refs_lz77 = ref_pool->GetFreeBackwardRefs();
  WP2_CHECK_STATUS(BackwardReferencesLz77(
      width, height, argb, {CacheType::kNone, 0}, hash_chain, &*refs_lz77));
  *refs = std::move(refs_lz77);
  return WP2_STATUS_OK;
}

// Gets the overall cost of storing the refs: symbol headers and pixels.
static WP2Status GetStorageCost(uint32_t width, const BackwardRefs& refs,
                                const LosslessSymbolsInfo& symbols_info,
                                uint32_t num_pixels,
                                float* const cost) {
  // Get the symbols stats.
  WP2::SymbolRecorder recorder;
  WP2::ANSEncNoop enc_noop;
  WP2_CHECK_STATUS(recorder.Allocate(symbols_info, /*num_records=*/0));
  WP2_CHECK_STATUS(StorePixels(width, refs, &enc_noop, &recorder));

  // Compute how much it would take to write the headers and the image.
  WP2::ANSEncCounter enc_counter;
  WP2::ANSDictionaries dicts;
  WP2::SymbolWriter symbol_writer;
  WP2_CHECK_STATUS(WriteHeaders(recorder, symbols_info, num_pixels,
                                &enc_counter, &dicts, &symbol_writer));
  WP2_CHECK_STATUS(StorePixels(width, refs, &enc_counter, &symbol_writer));
  *cost = enc_counter.GetCost();
  return WP2_STATUS_OK;
}

extern WP2Status BackwardReferencesTraceBackwards(
    uint32_t width, uint32_t height, const int16_t* const argb,
    const LosslessSymbolsInfo& symbols_info, const HashChain& hash_chain,
    const CacheConfig& cache_config, const BackwardRefs& refs_src,
    BackwardRefs* const refs_dst, bool* const check_histogram);
static WP2Status GetBackwardReferencesImpl(
    uint32_t width, uint32_t height, const int16_t* const argb, int speed,
    int lz77_types_to_try, uint32_t cache_bits_max, const HashChain& hash_chain,
    const LosslessSymbolsInfo& symbols_info, CacheConfig* const cache_config,
    BackwardRefsPool* ref_pool, BackwardRefsPool::RefsPtr* const refs_best) {
  LosslessSymbolsInfo symbols_info_tmp;
  WP2_CHECK_STATUS(symbols_info_tmp.CopyFrom(symbols_info));
  const uint32_t num_pixels = width * height;
  float bit_cost_best = std::numeric_limits<float>::max();
  HashChain hash_chain_box;

  auto best = ref_pool->GetFreeBackwardRefs();

  // TODO(vrabaud): Assess the utility of kLZ77RLE and kLZ77Box now that we use
  // ANS. Those methods were good to reduce the entropy of the distances but we
  // might only use LZ77 now that costs are better estimated.

  // Always try not to use LZ77. This is useful for one-color images for
  // example.
  lz77_types_to_try |= kLZ77None;

  const CacheConfig kNoCache{CacheType::kNone, 0};

  for (int lz77_type = 1; lz77_types_to_try;
       lz77_types_to_try &= ~lz77_type, lz77_type <<= 1) {
    if ((lz77_types_to_try & lz77_type) == 0) continue;
    auto refs_tmp = ref_pool->GetFreeBackwardRefs();
    BackwardRefs* const refs_tmp_ptr = refs_tmp.get();
    switch (lz77_type) {
      case kLZ77RLE:
        WP2_CHECK_STATUS(
            BackwardReferencesRle(width, height, argb, kNoCache, refs_tmp_ptr));
        break;
      case kLZ77Standard:
        // Compute LZ77 with no cache (0 bits), as the ideal LZ77 with a color
        // cache is not that different in practice.
        WP2_CHECK_STATUS(BackwardReferencesLz77(width, height, argb, kNoCache,
                                                hash_chain, refs_tmp_ptr));
        break;
      case kLZ77Box:
        WP2_CHECK_ALLOC_OK(hash_chain_box.Allocate(num_pixels));
        WP2_CHECK_STATUS(
            BackwardReferencesLz77Box(width, height, argb, kNoCache, hash_chain,
                                      &hash_chain_box, refs_tmp_ptr));
        break;
      case kLZ77None:
        WP2_CHECK_STATUS(
            BackwardReferencesNone(num_pixels, argb, kNoCache, refs_tmp_ptr));
        break;
      default:
        assert(false);
    }
    assert(refs_tmp->IsValid(num_pixels));

    // Next, try with a color cache and update the references.
    symbols_info_tmp.SetNumClusters(/*num_clusters=*/1);
    CacheConfig cache_config_tmp;
    WP2_CHECK_STATUS(CalculateBestCacheConfig(
        width, argb, num_pixels, speed, *refs_tmp, ref_pool, cache_bits_max,
        symbols_info_tmp, &cache_config_tmp));
    if (cache_config_tmp.type != CacheType::kNone) {
      symbols_info_tmp.SetCacheRange(GetColorCacheRange(cache_config_tmp));
      WP2_CHECK_STATUS(
          BackwardRefsWithLocalCache(argb, cache_config_tmp, refs_tmp_ptr));
    } else {
      symbols_info_tmp.SetCacheRange(0);
    }

    float bit_cost;
    WP2_CHECK_STATUS(GetStorageCost(width, *refs_tmp, symbols_info_tmp,
                                    num_pixels, &bit_cost));

    // Improve on simple LZ77 but only for high speed (TraceBackwards is
    // costly).
    if ((lz77_type == kLZ77Standard || lz77_type == kLZ77Box) && speed >= 3) {
      const HashChain& hash_chain_tmp =
          (lz77_type == kLZ77Standard) ? hash_chain : hash_chain_box;
      auto refs_tmp_trace = ref_pool->GetFreeBackwardRefs();

      bool check_histogram;
      WP2_CHECK_STATUS(BackwardReferencesTraceBackwards(
          width, height, argb, symbols_info_tmp, hash_chain_tmp,
          cache_config_tmp, *refs_tmp, &*refs_tmp_trace, &check_histogram));
      if (check_histogram) {
        assert(refs_tmp_trace->IsValid(num_pixels));
        float bit_cost_trace;
        WP2_CHECK_STATUS(GetStorageCost(width, *refs_tmp_trace,
                                        symbols_info_tmp, num_pixels,
                                        &bit_cost_trace));
        if (bit_cost_trace < bit_cost) {
          bit_cost = bit_cost_trace;
          refs_tmp.swap(refs_tmp_trace);
          //          *cache_config = cache_config_tmp;
        }
      }
    }

    // Keep the best backward references.
    if (bit_cost < bit_cost_best) {
      best.swap(refs_tmp);
      bit_cost_best = bit_cost;
      *cache_config = cache_config_tmp;
    }
  }

  *refs_best = std::move(best);

  return WP2_STATUS_OK;
}

WP2Status GetBackwardReferences(
    uint32_t width, uint32_t height, const int16_t* const argb, int speed,
    int lz77_types_to_try, uint32_t cache_bits_max, const HashChain& hash_chain,
    const LosslessSymbolsInfo& symbols_info, CacheConfig* const cache_config,
    BackwardRefsPool* const ref_pool, BackwardRefsPool::RefsPtr* const refs) {
  if (speed == 0) {
    cache_config->type = CacheType::kNone;
    cache_config->cache_bits = 0;
    return GetBackwardReferencesLowEffort(width, height, argb, hash_chain,
                                          ref_pool, refs);
  } else {
    return GetBackwardReferencesImpl(
        width, height, argb, speed, lz77_types_to_try, cache_bits_max,
        hash_chain, symbols_info, cache_config, ref_pool, refs);
  }
}

}  // namespace WP2L
