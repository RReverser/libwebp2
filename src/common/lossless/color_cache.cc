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
// Color Cache for WebP Lossless
//
// Author: Vincent Rabaud (vrabaud@google.com)

#include "src/common/lossless/color_cache.h"

#include <algorithm>

#include "src/dsp/math.h"

namespace WP2L {

WP2Status ColorCachePtr::Init(CacheConfig config) {
  switch (config.type) {
    case CacheType::kHash: {
      HashColorCache* hash_cache = new (WP2Allocable::nothrow) HashColorCache();
      WP2_CHECK_ALLOC_OK(hash_cache != nullptr);
      cache = hash_cache;
      WP2_CHECK_ALLOC_OK(hash_cache->Allocate(config.cache_bits));
      return WP2_STATUS_OK;
    }
    case CacheType::kFifo: {
      FifoColorCache* fifo_cache = new (WP2Allocable::nothrow) FifoColorCache();
      WP2_CHECK_ALLOC_OK(fifo_cache != nullptr);
      cache = fifo_cache;
      WP2_CHECK_ALLOC_OK(fifo_cache->Allocate(config.cache_bits));
      return WP2_STATUS_OK;
    }
    case CacheType::kNone:
      return WP2_STATUS_OK;
  }
  assert(false);
  return WP2_STATUS_OK;
}

uint32_t GetColorCacheRange(CacheConfig config) {
  switch (config.type) {
    case CacheType::kFifo:
    case CacheType::kHash:
      return 1 << config.cache_bits;
    case CacheType::kNone:
      return 0;
  }
  assert(false);
  return 0;
}

//------------------------------------------------------------------------------
// ColorCacheMap.

WP2Status ColorCacheMap::Allocate(size_t color_num) {
  hash_bits_ = 1 + WP2Log2Floor(color_num);
  WP2_CHECK_ALLOC_OK(colors_.resize(1 << hash_bits_));
  WP2_CHECK_ALLOC_OK(is_present_.resize(colors_.size()));
  Reset();
  return WP2_STATUS_OK;
}

void ColorCacheMap::Reset() {
  // Make sure all the Color.is_present is false.
  memset(is_present_.data(), 0, is_present_.size() * sizeof(bool));
}

bool ColorCacheMap::Insert(const int16_t color[4], uint32_t value) {
  uint32_t tag = HashPix(color, hash_bits_);
  const uint32_t start_tag = tag;
  const size_t size = colors_.size();
  while (true) {
    Color* const c = &colors_[tag];
    bool& is_present = is_present_[tag];
    if (!is_present) {
      ColorCopy(color, &c->color[0]);
      is_present = true;
      c->value = value;
      return true;
    } else if (ColorEq(c->color, color)) {
      c->value = value;
      return false;  // The color is already there.
    } else {
      // Some other color sits here, so do linear conflict resolution.
      ++tag;
      tag %= size;
      if (tag == start_tag) {
        // We looped around, the map is full! Maybe you need a bigger map?
        assert(false);
      }
    }
  }
}

void ColorCacheMap::GetPalette(int16_t* const palette) const {
  for (uint32_t i = 0; i < colors_.size(); ++i) {
    if (is_present_[i]) {
      // For palettes, we assume the color's 'value' is its insertion index.
      ColorCopy(colors_[i].color, &palette[4 * colors_[i].value]);
    }
  }
}

bool ColorCacheMap::Contains(const int16_t color[4],
                             uint32_t* const value) const {
  uint32_t tag = HashPix(color, hash_bits_);
  const uint32_t start_tag = tag;
  const size_t size = colors_.size();
  while (true) {
    const Color& c = colors_[tag];
    if (!is_present_[tag]) {
      return false;
    } else if (ColorEq(c.color, color)) {
      *value = c.value;
      return true;  // The color is already there.
    } else {
      // Some other color sits here, so do linear conflict resolution.
      ++tag;
      tag %= size;
      if (tag == start_tag) return false;
    }
  }
}

//------------------------------------------------------------------------------
// HashColorCache.

void HashColorCache::Lookup(uint32_t index, int16_t* const dest) const {
  assert((index >> hash_bits_) == 0u);
  ColorCopy(&colors_[4 * index], dest);
}

bool HashColorCache::Insert(const int16_t color[4], uint32_t* const index_ptr) {
  const uint32_t tag = HashPix(color, hash_bits_);
  if (index_ptr != nullptr) *index_ptr = tag;
  if (ColorEq(&colors_[4 * tag], color)) {
    return false;
  } else {
    ColorCopy(color, &colors_[4 * tag]);
    return true;
  }
}

bool HashColorCache::Contains(const int16_t color[4],
                              uint32_t* const index_ptr) const {
  const uint32_t index = HashPix(color, hash_bits_);
  if (ColorEq(color, &colors_[4 * index])) {
    *index_ptr = index;
    return true;
  } else {
    return false;
  }
}

bool HashColorCache::Allocate(uint32_t hash_bits) {
  const uint64_t hash_size = (1 << hash_bits);
  assert(hash_bits > 0);
  if (!colors_.resize(4 * hash_size)) return false;
  hash_bits_ = hash_bits;
  Reset();
  return true;
}

uint32_t HashColorCache::IndexRange() const { return colors_.size() / 4; }

void HashColorCache::Reset() {
  // TODO(yguyon): Create a kUndefined, fill the vector with that and verify
  //               that no kUndefined color is read from the cache.
  std::fill(colors_.begin(), colors_.end(), 0xffu);
}

//------------------------------------------------------------------------------
// FifoColorCache.

void FifoColorCache::Lookup(uint32_t index, int16_t* const dest) const {
  assert((index >> cache_bits_) == 0u);
  ColorCopy(&colors_[4 * index], dest);
}

bool FifoColorCache::Insert(const int16_t color[4], uint32_t* const index_ptr) {
  // Naive (slow) implementation. TODO(maryla): optimize.
  const uint32_t max_colors = colors_.size() / 4;
  uint32_t current_index = max_colors - 1;
  const bool already_present = Contains(color, &current_index);
  if (index_ptr != nullptr) *index_ptr = current_index;
  // Shift colors back by one.
  for (uint32_t i = current_index; i > 0; --i) {
    ColorCopy(&colors_[4 * (i - 1)], &colors_[4 * i]);
  }
  // Insert the new color at the front.
  ColorCopy(color, &colors_[4 * 0]);
  if (!already_present && num_colors_ < max_colors) {
    ++num_colors_;
  }
  return !already_present;
}

bool FifoColorCache::Contains(const int16_t color[4],
                              uint32_t* const index) const {
  // Naive (slow) implementation. TODO(maryla): optimize.
  for (uint32_t i = 0; i < num_colors_; ++i) {
    if (ColorEq(color, &colors_[4 * i])) {
      *index = i;
      return true;
    }
  }
  return false;
}

bool FifoColorCache::Allocate(uint32_t cache_bits) {
  const uint64_t cache_size = (1 << cache_bits);
  assert(cache_bits > 0);
  if (!colors_.resize(4 * cache_size)) return false;
  cache_bits_ = cache_bits;
  Reset();
  return true;
}

uint32_t FifoColorCache::IndexRange() const { return num_colors_; }

void FifoColorCache::Reset() {
  num_colors_ = 0;
  std::fill(colors_.begin(), colors_.end(), 0xffu);
}

//------------------------------------------------------------------------------
// HybridColorCache.

void HybridColorCache::Init(ColorCache* const cache1,
                            ColorCache* const cache2) {
  cache1_ = cache1;
  cache2_ = cache2;
}

void HybridColorCache::Lookup(uint32_t index, int16_t* const dest) const {
  const uint32_t range1 = cache1_->IndexRange();
  if (index < range1) {
    cache1_->Lookup(index, dest);
  } else {
    cache2_->Lookup(index - range1, dest);
  }
}

bool HybridColorCache::Insert(const int16_t color[4],
                              uint32_t* const index_ptr) {
  uint32_t index1, index2;
  const bool already_present1 = !cache1_->Insert(color, &index1);
  const bool already_present2 = !cache2_->Insert(color, &index2);
  if (index_ptr != nullptr) {
    if (already_present1) {
      *index_ptr = index1;
    } else if (already_present2) {
      *index_ptr = index2 + cache1_->IndexRange();
    }
  }
  return !(already_present1 || already_present2);
}

bool HybridColorCache::Contains(const int16_t color[4],
                                uint32_t* const index_ptr) const {
  if (cache1_->Contains(color, index_ptr)) return true;
  if (cache2_->Contains(color, index_ptr)) {
    *index_ptr = *index_ptr + cache1_->IndexRange();
    return true;
  }
  return false;
}

uint32_t HybridColorCache::IndexRange() const {
  return cache1_->IndexRange() + cache2_->IndexRange();
}

void HybridColorCache::Reset() {
  cache1_->Reset();
  cache2_->Reset();
}

//------------------------------------------------------------------------------
// MoveToFrontCache.

WP2Status MoveToFrontCache::Init(bool enabled, uint32_t num_colors) {
  enabled_ = enabled;
  num_colors_ = num_colors;
  if (!enabled_) return WP2_STATUS_OK;
  WP2_CHECK_ALLOC_OK(colors_.resize(num_colors_));
  WP2_CHECK_ALLOC_OK(indices_.resize(num_colors_));
  Reset();
  return WP2_STATUS_OK;
}

void MoveToFrontCache::Reset() {
  if (!enabled_) return;
  for (uint16_t i = 0; i < num_colors_; ++i) {
    colors_[i] = indices_[i] = i;
  }
}

void MoveToFrontCache::MoveToFront(uint16_t color) {
  assert(color < num_colors_);
  if (!enabled_) return;
  const uint16_t idx = GetIndex(color);
  for (uint16_t i = idx; i > 0; --i) {
    colors_[i] = colors_[i - 1];
    indices_[colors_[i]] = i;
  }
  colors_[0] = color;
  indices_[color] = 0;
}

uint16_t MoveToFrontCache::GetIndex(uint16_t color) const {
  assert(color < num_colors_);
  if (!enabled_) return color;
  assert(colors_[indices_[color]] == color);
  return indices_[color];
}

uint16_t MoveToFrontCache::GetColor(uint16_t idx) const {
  assert(idx < num_colors_);
  if (!enabled_) return idx;
  assert(indices_[colors_[idx]] == idx);
  return colors_[idx];
}

}  // namespace WP2L
