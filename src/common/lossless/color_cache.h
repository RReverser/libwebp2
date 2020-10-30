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
// Authors: Vincent Rabaud (vrabaud@google.com)

#ifndef WP2_UTILS_LOSSLESS_COLOR_CACHE_H_
#define WP2_UTILS_LOSSLESS_COLOR_CACHE_H_

#include <cassert>
#include <cstring>

#include "src/dsp/math.h"
#include "src/utils/vector.h"
#include "src/wp2/base.h"

namespace WP2L {

// Returns true if the two colors are equal.
template <typename T>
static inline bool ColorEq(const T* const color1, const T* const color2) {
  return (std::memcmp(color1, color2, 4 * sizeof(T)) == 0);
}
template <typename T>
static inline void ColorCopy(const T* const color1, T* const color2) {
  std::memcpy(color2, color1, 4 * sizeof(T));
}

// Define a n_bits hash value of a color.
template <typename T>
static inline uint32_t HashPix(const T* const color, int n_bits) {
  // TODO(vrabaud) Figure out a common hash for all 16-bit ARGB colors.
  // Note that masking with 0xffffffffu is for preventing an
  // 'unsigned int overflow' warning. Doesn't impact the compiled code.
  const uint64_t kHashMul = 0x1e35a7bdull;
  const uint32_t color32 =
      (color[0] << 24) | (color[1] << 16) | (color[2] << 8) | (color[3]);
  return (uint32_t)(((color32 * kHashMul) & 0xffffffffu) >> (32 - n_bits));
}

// Color cache that behaves like a map with (key,value), where the key is a
// color and the value whatever uint32_t the user wants.
class ColorCacheMap {
 public:
  WP2Status Allocate(size_t color_num_max);
  void Reset();

  // Inserts the tuple (color, value). Returns true if the tuple was actually
  // inserted, false if it was already present.
  bool Insert(const int16_t color[4], uint32_t value);

  // Returns whether the cache contains the value, and also gives back the index
  // accordingly.
  bool Contains(const int16_t color[4], uint32_t* const value) const;

  // Copy the palette into a buffer.
  void GetPalette(int16_t* const palette) const;

 private:
  struct Color {
    int16_t color[4];
    uint16_t value;
  };

  WP2::VectorNoCtor<Color> colors_;
  WP2::VectorNoCtor<bool> is_present_;
  int hash_bits_ = 0;
};

enum class CacheType {
  kNone,
  kHash,
  kFifo,
};

struct CacheConfig {
  CacheType type;
  uint32_t cache_bits;
};

// A generic color cache.
class ColorCache : public WP2Allocable {
 public:
  virtual ~ColorCache() = default;

  virtual void Lookup(uint32_t index, int16_t* const dest) const = 0;

  // Inserts the color in the cache. Returns true if the color was actually
  // inserted, false if it was already present. If false, sets index_ptr to the
  // index for this color *prior to inserting*.
  virtual bool Insert(const int16_t color[4], uint32_t* const index_ptr) = 0;

  // Returns whether the cache contains the value, and also gives back the index
  // accordingly.
  virtual bool Contains(const int16_t color[4],
                        uint32_t* const index_ptr) const = 0;

  // Returns the range of a color index given the current state of the cache
  // (e.g. number of colors in the cache).
  virtual uint32_t IndexRange() const = 0;

  virtual void Reset() = 0;
};

// Struct that own a color cache pointer.
struct ColorCachePtr {
  // Instantiates and inits the color cache based on the given config.
  WP2Status Init(CacheConfig config);
  ~ColorCachePtr() { delete cache; }
  ColorCache* operator->() { return cache; }
  ColorCache* get() const { return cache; }
  ColorCache* cache = nullptr;
};

// Max range of the cache index for the given cache config.
uint32_t GetColorCacheRange(CacheConfig config);

// Very fast hash-based color cache. A given color always has the same index,
// and only gets evicted when a different color with the same hash is inserted.
class HashColorCache : public ColorCache {
 public:
  // Initializes the color cache with 'hash_bits' bits for the hashes.
  // Returns false in case of memory error.
  WP2_NO_DISCARD bool Allocate(uint32_t hash_bits);

  void Lookup(uint32_t index, int16_t* const dest) const override;
  bool Insert(const int16_t color[4], uint32_t* const index_ptr) override;
  bool Contains(const int16_t color[4],
                uint32_t* const index_ptr) const override;
  uint32_t IndexRange() const override;
  void Reset() override;

 private:
  WP2::Vector_s16 colors_;  // color entries
  int hash_bits_ = 0;
};

class FifoColorCache : public ColorCache {
 public:
  // Initializes the color cache with a size of 1<<cache_bits.
  // Returns false in case of memory error.
  WP2_NO_DISCARD bool Allocate(uint32_t cache_bits);

  void Lookup(uint32_t index, int16_t* const dest) const override;
  bool Insert(const int16_t color[4], uint32_t* const index_ptr) override;
  bool Contains(const int16_t color[4],
                uint32_t* const index_ptr) const override;
  uint32_t IndexRange() const override;
  void Reset() override;

 private:
  WP2::Vector_s16 colors_;  // color entries
  uint32_t cache_bits_ = 0;
  uint32_t num_colors_ = 0;
};

// Color cache that mixes two other caches.
class HybridColorCache : public ColorCache {
 public:
  void Init(ColorCache* const cache1, ColorCache* const cache2);

  void Lookup(uint32_t index, int16_t* const dest) const override;
  bool Insert(const int16_t color[4], uint32_t* const index_ptr) override;
  bool Contains(const int16_t color[4],
                uint32_t* const index_ptr) const override;
  uint32_t IndexRange() const override;
  void Reset() override;

 private:
  ColorCache* cache1_;
  ColorCache* cache2_;
};

// A cache for indexed colors (a color is just an uint16_t, so this is
// basically a remapping). Colors start in increasing order, and any color can
// be moved to the front, with all other colors shifted to the back.
class MoveToFrontCache {
 public:
  // enabled: if false, this cache never actually moves colors
  WP2Status Init(bool enabled, uint32_t num_colors);

  // Resets colors to the initial sorted order.
  void Reset();

  // Moves the given color to index 0. Colors that were in front of this one get
  // shifted towards the back.
  void MoveToFront(uint16_t color);

  // Returns the index of this color in the cache.
  uint16_t GetIndex(uint16_t color) const;
  // Returns the color at this index.
  uint16_t GetColor(uint16_t idx) const;

  uint16_t num_colors() const { return num_colors_; }
  bool enabled() const { return enabled_; }

 private:
  bool enabled_;
  uint16_t num_colors_;
  WP2::Vector_u16 colors_;
  WP2::Vector_u16 indices_;
};

}  // namespace WP2L

#endif  // WP2_UTILS_LOSSLESS_COLOR_CACHE_H_
