// Copyright 2020 Google LLC
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//      http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.
// -----------------------------------------------------------------------------
//
// Test color caches.
//
// Author: Maryla (maryla@google.com)

#include "src/common/lossless/color_cache.h"
#include "src/common/symbols.h"
#include "src/utils/random.h"
#include "tests/include/helpers.h"

namespace WP2L {
namespace {

TEST(HashColorCache, Simple) {
  HashColorCache cache;
  ASSERT_TRUE(cache.Allocate(/*hash_bits=*/2));

  const int16_t kBlack[] = {255, 0, 0, 0};
  const int16_t kGreen[] = {255, 0, 225, 0};
  const int16_t kRed[] = {255, 255, 0, 0};

  uint32_t index;
  EXPECT_FALSE(cache.Contains(kBlack, &index));
  EXPECT_FALSE(cache.Contains(kGreen, &index));
  EXPECT_FALSE(cache.Contains(kRed, &index));

  // Insert black.
  EXPECT_TRUE(cache.Insert(kBlack, &index));
  const uint32_t black_index = index;
  int16_t color[4];
  cache.Lookup(black_index, color);
  EXPECT_THAT(color, testing::ElementsAreArray(kBlack));

  EXPECT_TRUE(cache.Contains(kBlack, &index));
  EXPECT_FALSE(cache.Contains(kGreen, &index));
  EXPECT_FALSE(cache.Contains(kRed, &index));

  // Insert green.
  EXPECT_TRUE(cache.Insert(kGreen, &index));
  const uint32_t green_index = index;
  const bool black_was_removed = (green_index == black_index);
  cache.Lookup(green_index, color);
  EXPECT_THAT(color, testing::ElementsAreArray(kGreen));

  EXPECT_TRUE(cache.Contains(kGreen, &index));
  EXPECT_EQ(cache.Contains(kBlack, &index), black_index != green_index);
  if (!black_was_removed) {
    cache.Lookup(black_index, color);
    EXPECT_THAT(color, testing::ElementsAreArray(kBlack));
  }

  // Re insert black, check the index hasn't changed.
  EXPECT_EQ(cache.Insert(kBlack, &index), black_was_removed);
  EXPECT_EQ(black_index, index);

  // Insert red.
  EXPECT_TRUE(cache.Insert(kRed, &index));
  const uint32_t red_index = index;
  cache.Lookup(red_index, color);
  EXPECT_THAT(color, testing::ElementsAreArray(kRed));
  EXPECT_EQ(cache.Contains(kBlack, &index), red_index != black_index);
  EXPECT_EQ(cache.Contains(kGreen, &index), red_index != green_index);
}

// Tests that the cache API is consistent: e.g. a color that is present
// can be retrieved, a color that was just inserted is present, etc.
static void TestColorCacheAPI(uint32_t num_colors, ColorCache* const cache,
                              WP2::UniformIntDistribution* const random) {
  const uint32_t num_pixels = num_colors * 10;
  std::vector<int16_t> colors(num_colors * 4);
  for (uint32_t i = 0; i < colors.size(); ++i) {
    colors[i] = random->Get(0, 255);
  }

  std::vector<bool> color_seen(num_colors, false);
  std::vector<uint32_t> color_indices(num_pixels);
  const uint32_t kNotInCache = std::numeric_limits<uint32_t>::max();
  std::vector<uint32_t> cache_indices(num_pixels, kNotInCache);

  // Simulates encoding stage.
  for (uint32_t i = 0; i < num_pixels; ++i) {
    const uint32_t color_index = random->Get<uint32_t>(0, num_colors - 1);
    color_indices[i] = color_index;
    const int16_t* const color = &colors[color_index * 4];

    // If we've never seen this color, it shouldn't be in the cache
    // (but it's not necessarily in the cache even if we've seen it before).
    const bool present = cache->Contains(color, &cache_indices[i]);
    if (!color_seen[color_index]) {
      EXPECT_FALSE(present);
    }
    color_seen[color_index] = true;

    int16_t lookup[4];
    if (present) {
      cache->Lookup(cache_indices[i], lookup);
      EXPECT_THAT(lookup,
                  testing::ElementsAre(color[0], color[1], color[2], color[3]));
    }

    uint32_t insert_index;
    const bool inserted = cache->Insert(color, &insert_index);
    EXPECT_NE(inserted, present);
    if (present) {
      EXPECT_EQ(insert_index, cache_indices[i]);
    }

    // A color that was just inserted is in the cache.
    uint32_t new_index;
    EXPECT_TRUE(cache->Contains(color, &new_index));
    cache->Lookup(new_index, lookup);
    EXPECT_THAT(lookup,
                testing::ElementsAre(color[0], color[1], color[2], color[3]));
  }

  cache->Reset();

  // Simulates decoding stage.
  for (uint32_t i = 0; i < num_pixels; ++i) {
    const int16_t* const expected_color = &colors[color_indices[i] * 4];
    if (cache_indices[i] == kNotInCache) {
      uint32_t cache_index;
      EXPECT_FALSE(cache->Contains(expected_color, &cache_index));
    } else {
      uint32_t cache_index;
      EXPECT_TRUE(cache->Contains(expected_color, &cache_index));

      int16_t lookup[4];
      cache->Lookup(cache_indices[i], lookup);
      EXPECT_THAT(lookup,
                  testing::ElementsAre(expected_color[0], expected_color[1],
                                       expected_color[2], expected_color[3]));
    }
    cache->Insert(expected_color, nullptr);
  }
}

TEST(HashColorCache, Random) {
  for (uint32_t bits = 1; bits < 5; ++bits) {
    WP2::UniformIntDistribution random;
    HashColorCache cache;
    ASSERT_TRUE(cache.Allocate(bits));

    const uint32_t num_colors = random.Get(1 << (bits - 1), 1 << (bits + 1));
    TestColorCacheAPI(num_colors, &cache, &random);
  }
}

TEST(FifoColorCache, Simple) {
  FifoColorCache cache;
  ASSERT_TRUE(cache.Allocate(/*cache_bits=*/2));

  int16_t kColors[5 * 4] = {0, 0, 0, 0, 1, 1, 1, 1, 2, 2,
                            2, 2, 3, 3, 3, 3, 4, 4, 4, 4};

  uint32_t index = 0;
  EXPECT_FALSE(cache.Contains(&kColors[4 * 0], &index));
  EXPECT_TRUE(cache.Insert(&kColors[4 * 0], &index));
  EXPECT_EQ(cache.IndexRange(), 1u);

  EXPECT_TRUE(cache.Insert(&kColors[4 * 1], &index));
  EXPECT_TRUE(cache.Contains(&kColors[4 * 1], &index));
  EXPECT_EQ(index, 0u);
  EXPECT_TRUE(cache.Contains(&kColors[4 * 0], &index));
  EXPECT_EQ(index, 1u);
  EXPECT_EQ(cache.IndexRange(), 2u);

  EXPECT_FALSE(cache.Insert(&kColors[4 * 0], &index));
  EXPECT_EQ(index, 1u);
  EXPECT_TRUE(cache.Contains(&kColors[4 * 0], &index));
  EXPECT_EQ(index, 0u);
  EXPECT_TRUE(cache.Contains(&kColors[4 * 1], &index));
  EXPECT_EQ(index, 1u);
  EXPECT_EQ(cache.IndexRange(), 2u);

  EXPECT_TRUE(cache.Insert(&kColors[4 * 2], &index));
  EXPECT_TRUE(cache.Contains(&kColors[4 * 2], &index));
  EXPECT_EQ(index, 0u);
  EXPECT_TRUE(cache.Contains(&kColors[4 * 0], &index));
  EXPECT_EQ(index, 1u);
  EXPECT_TRUE(cache.Contains(&kColors[4 * 1], &index));
  EXPECT_EQ(index, 2u);
  EXPECT_EQ(cache.IndexRange(), 3u);

  EXPECT_TRUE(cache.Insert(&kColors[4 * 3], &index));
  EXPECT_TRUE(cache.Contains(&kColors[4 * 3], &index));
  EXPECT_EQ(index, 0u);
  EXPECT_TRUE(cache.Contains(&kColors[4 * 2], &index));
  EXPECT_EQ(index, 1u);
  EXPECT_TRUE(cache.Contains(&kColors[4 * 0], &index));
  EXPECT_EQ(index, 2u);
  EXPECT_TRUE(cache.Contains(&kColors[4 * 1], &index));
  EXPECT_EQ(index, 3u);
  EXPECT_EQ(cache.IndexRange(), 4u);

  EXPECT_FALSE(cache.Insert(&kColors[4 * 1], &index));
  EXPECT_EQ(index, 3u);
  EXPECT_TRUE(cache.Contains(&kColors[4 * 1], &index));
  EXPECT_EQ(index, 0u);
  EXPECT_TRUE(cache.Contains(&kColors[4 * 3], &index));
  EXPECT_EQ(index, 1u);
  EXPECT_TRUE(cache.Contains(&kColors[4 * 2], &index));
  EXPECT_EQ(index, 2u);
  EXPECT_TRUE(cache.Contains(&kColors[4 * 0], &index));
  EXPECT_EQ(index, 3u);
  EXPECT_EQ(cache.IndexRange(), 4u);

  EXPECT_TRUE(cache.Insert(&kColors[4 * 4], &index));
  EXPECT_TRUE(cache.Contains(&kColors[4 * 4], &index));
  EXPECT_EQ(index, 0u);
  EXPECT_TRUE(cache.Contains(&kColors[4 * 1], &index));
  EXPECT_EQ(index, 1u);
  EXPECT_TRUE(cache.Contains(&kColors[4 * 3], &index));
  EXPECT_EQ(index, 2u);
  EXPECT_TRUE(cache.Contains(&kColors[4 * 2], &index));
  EXPECT_EQ(index, 3u);
  EXPECT_FALSE(cache.Contains(&kColors[4 * 0], &index));
  EXPECT_EQ(cache.IndexRange(), 4u);
}

TEST(FifoColorCache, Random) {
  for (uint32_t bits = 1; bits < 5; ++bits) {
    WP2::UniformIntDistribution random;
    FifoColorCache cache;
    ASSERT_TRUE(cache.Allocate(bits));

    const uint32_t num_colors = random.Get(1 << (bits - 1), 1 << (bits + 1));
    TestColorCacheAPI(num_colors, &cache, &random);
  }
}

TEST(HybridColorCache, Random) {
  FifoColorCache fifo_cache;
  const uint32_t fifo_bits = 2;
  ASSERT_TRUE(fifo_cache.Allocate(fifo_bits));
  HashColorCache hash_cache;
  const uint32_t hash_bits = 8;
  ASSERT_TRUE(hash_cache.Allocate(hash_bits));
  const uint32_t bits = fifo_bits + hash_bits;

  HybridColorCache hybrid_cache;
  hybrid_cache.Init(&fifo_cache, &hash_cache);

  WP2::UniformIntDistribution random;
  const uint32_t num_colors = random.Get(1 << (bits - 1), 1 << (bits + 1));
  TestColorCacheAPI(num_colors, &hybrid_cache, &random);
}

TEST(MoveToFrontCache, Simple) {
  MoveToFrontCache cache;
  ASSERT_WP2_OK(cache.Init(/*enabled=*/true, /*num_colors=*/4));

  EXPECT_EQ(cache.GetIndex(0), 0);
  EXPECT_EQ(cache.GetIndex(1), 1);
  EXPECT_EQ(cache.GetIndex(2), 2);
  EXPECT_EQ(cache.GetIndex(3), 3);
  EXPECT_EQ(cache.GetColor(0), 0);
  EXPECT_EQ(cache.GetColor(1), 1);
  EXPECT_EQ(cache.GetColor(2), 2);
  EXPECT_EQ(cache.GetColor(3), 3);

  cache.MoveToFront(2);
  EXPECT_EQ(cache.GetIndex(2), 0);
  EXPECT_EQ(cache.GetIndex(0), 1);
  EXPECT_EQ(cache.GetIndex(1), 2);
  EXPECT_EQ(cache.GetIndex(3), 3);
  EXPECT_EQ(cache.GetColor(0), 2);
  EXPECT_EQ(cache.GetColor(1), 0);
  EXPECT_EQ(cache.GetColor(2), 1);
  EXPECT_EQ(cache.GetColor(3), 3);

  for (uint32_t i = 0; i < 2; ++i) {
    cache.MoveToFront(3);
    EXPECT_EQ(cache.GetIndex(3), 0);
    EXPECT_EQ(cache.GetIndex(2), 1);
    EXPECT_EQ(cache.GetIndex(0), 2);
    EXPECT_EQ(cache.GetIndex(1), 3);
    EXPECT_EQ(cache.GetColor(0), 3);
    EXPECT_EQ(cache.GetColor(1), 2);
    EXPECT_EQ(cache.GetColor(2), 0);
    EXPECT_EQ(cache.GetColor(3), 1);
  }

  cache.MoveToFront(0);
  EXPECT_EQ(cache.GetIndex(0), 0);
  EXPECT_EQ(cache.GetIndex(3), 1);
  EXPECT_EQ(cache.GetIndex(2), 2);
  EXPECT_EQ(cache.GetIndex(1), 3);
  EXPECT_EQ(cache.GetColor(0), 0);
  EXPECT_EQ(cache.GetColor(1), 3);
  EXPECT_EQ(cache.GetColor(2), 2);
  EXPECT_EQ(cache.GetColor(3), 1);
}

TEST(MoveToFrontCache, RandomTest) {
  constexpr uint32_t kNumColors = 100;

  for (bool enabled : {true, false}) {
    std::vector<uint16_t> colors(kNumColors);
    std::iota(colors.begin(), colors.end(), 0);
    WP2::Shuffle(colors.begin(), colors.end(), /*seed=*/0);

    MoveToFrontCache cache;
    ASSERT_WP2_OK(cache.Init(enabled, kNumColors));

    for (uint32_t i = 0; i < kNumColors; ++i) {
      cache.MoveToFront(colors[kNumColors - i - 1]);
    }

    for (uint32_t i = 0; i < kNumColors; ++i) {
      if (enabled) {
        EXPECT_EQ(cache.GetIndex(colors[i]), i);
        EXPECT_EQ(cache.GetColor(i), colors[i]);
      } else {
        EXPECT_EQ(cache.GetIndex(i), i);
        EXPECT_EQ(cache.GetColor(i), i);
      }
    }
  }
}

}  // namespace
}  // namespace WP2L
