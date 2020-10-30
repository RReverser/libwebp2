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

// BlockSize tests.

#include "include/helpers.h"
#include "src/common/lossy/block_size.h"
#include "src/common/lossy/block_size_io.h"

namespace WP2 {
namespace {

//------------------------------------------------------------------------------

TEST(BlockSizeTest, GetBlockSize) {
  for (BlockSize block_size : kAllBlockSizes) {
    ASSERT_EQ(GetBlockSize(BlockWidth[block_size], BlockHeight[block_size]),
              block_size);
  }
}

//------------------------------------------------------------------------------

TEST(BlockSizeTest, PartitionSetStats) {
  for (uint32_t i = 0; i < NUM_PARTITION_SETS; ++i) {
    const PartitionSet set = (PartitionSet)i;
    const PartitionSetStats& stats = kPartSetsStats[set];
    ASSERT_EQ(stats.partition_set_, set);
    ASSERT_GE(stats.num_unique_bounds_, stats.num_blocks_[BLK_32x32]);
    uint32_t num_unique_bounds = 0;

    for (BlockSize bounds : kAllBlockSizes) {
      const BlockSize smallest_bounds = stats.smallest_bounds_[bounds];
      ASSERT_LE(BlockWidth[smallest_bounds], BlockWidth[bounds]);
      ASSERT_LE(BlockHeight[smallest_bounds], BlockHeight[bounds]);
      ASSERT_EQ(stats.blocks_[bounds][stats.num_blocks_[bounds]], BLK_LAST);
      ASSERT_GE(stats.num_blocks_[bounds], 1);
      if (bounds == BLK_4x4) {
        ASSERT_EQ(stats.num_blocks_[bounds], 1);
      }
      ASSERT_EQ(stats.num_blocks_[bounds], stats.num_blocks_[smallest_bounds]);
      for (uint32_t b = 0; b < stats.num_blocks_[bounds]; ++b) {
        const BlockSize block_size = stats.blocks_[bounds][b];
        ASSERT_LE(BlockWidth[block_size], BlockWidth[smallest_bounds]);
        ASSERT_LE(BlockHeight[block_size], BlockHeight[smallest_bounds]);
        ASSERT_EQ(stats.blocks_[bounds][b], stats.blocks_[smallest_bounds][b]);
        ASSERT_EQ(stats.slots_[bounds][block_size], b);
      }
      if (stats.slots_[bounds][bounds] != PartitionSetStats::kInvalidSlot) {
        ASSERT_EQ(stats.blocks_[bounds][stats.num_blocks_[bounds] - 1], bounds);
        ASSERT_EQ(stats.blocks_[bounds][stats.slots_[bounds][bounds]], bounds);
      }
      if (smallest_bounds == bounds) {
        ++num_unique_bounds;
      } else {
        // No more than one hop.
        ASSERT_EQ(stats.smallest_bounds_[smallest_bounds], smallest_bounds);
      }
      ASSERT_LT(stats.unique_bounds_index_[bounds], stats.num_unique_bounds_);
    }
    ASSERT_EQ(num_unique_bounds, stats.num_unique_bounds_);
  }
}

//------------------------------------------------------------------------------

class PartitionSetTest : public ::testing::TestWithParam<
                             std::tuple<std::string, float, PartitionSet>> {};

TEST_P(PartitionSetTest, Simple) {
  const std::string& file_name = std::get<0>(GetParam());
  EncoderConfig config = EncoderConfig::kDefault;
  config.quality = std::get<1>(GetParam());
  config.partition_set = std::get<2>(GetParam());

  ASSERT_WP2_OK(testing::EncodeDecodeCompare(file_name, config));
}

INSTANTIATE_TEST_SUITE_P(
    PartitionSetTestInstantiation, PartitionSetTest,
    ::testing::Combine(::testing::Values("source1_64x48.png"),
                       /*quality=*/::testing::Values(0, 75, kMaxLossyQuality),
                       ::testing::Range(SMALL_SQUARES, NUM_PARTITION_SETS)));

//------------------------------------------------------------------------------

}  // namespace
}  // namespace WP2
