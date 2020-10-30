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

// testing FrontMgr class

#include <algorithm>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <set>
#include <vector>

#include "include/helpers.h"
#include "src/utils/front_mgr.h"
#include "src/wp2/format_constants.h"

namespace WP2 {
namespace {

//------------------------------------------------------------------------------

void CreateData(uint32_t w, uint32_t h, bool snapped,
                UniformIntDistribution* const gen,
                VectorNoCtor<Block>* const blocks, bool dump = false) {
  testing::CreatePartition(w, h, snapped, gen, blocks);
  std::vector<uint8_t> occupancy(w * h, 0x00);
  uint8_t color = 0xff;
  for (const Block& block : *blocks) {
    const uint32_t px = block.x_pix();
    const uint32_t py = block.y_pix();
    const uint32_t dx = block.w_pix();
    const uint32_t dy = block.h_pix();

    // draw a rectangle, check we're not overlapping
    for (uint32_t y = py; y < std::min(py + dy, h); ++y) {
      for (uint32_t x = px; x < std::min(px + dx, w); ++x) {
        EXPECT_EQ(occupancy[y * w + x], 0x00) << "Overlap check failed.";
        // Draw top and left borders with a different color for each block.
        occupancy[y * w + x] = (x == px || y == py) ? color : 0x01;
      }
    }
    color += 17;
    color |= 0x80;  // Make sure it is not zero and bright enough.
  }
  // check we've covered all the area
  ASSERT_TRUE(std::find(occupancy.begin(), occupancy.end(), 0x00) ==
              occupancy.end())
      << "All the area should have been covered!";

  if (dump) {
    FILE* const f = fopen("/tmp/dump.ppm", "wb");
    fprintf(f, "P5\n%d %d\n255\n", w, h);
    for (uint8_t v : occupancy) fprintf(f, "%c", v);
    fclose(f);
  }
}

//------------------------------------------------------------------------------

enum FrontMgrType { kLexico, kMax, kArea };
FrontMgrNxNBase* GetMgr(FrontMgrType front_mgr_type,
                        FrontMgrLexico* const mgr_lex,
                        FrontMgrMax* const mgr_max,
                        FrontMgrArea* const mgr_area) {
  if (front_mgr_type == kLexico) return mgr_lex;
  if (front_mgr_type == kMax) return mgr_max;
  if (front_mgr_type == kArea) return mgr_area;
  return nullptr;
}
bool GetSnapped(FrontMgrType type) { return (type == kArea); }

constexpr uint32_t kAreaSize = kMaxBlockSizePix;

class TestFrontMgr : public ::testing::TestWithParam<
                         std::tuple<uint32_t, uint32_t, FrontMgrType>> {};

TEST_P(TestFrontMgr, LeftContext) {
  const uint32_t seed = std::get<0>(GetParam());
  const uint32_t num_iterations = std::get<1>(GetParam());
  const FrontMgrType front_mgr_type = std::get<2>(GetParam());

  UniformIntDistribution gen(seed);
  const uint32_t kRangeMax = 100;
  for (size_t n = 0; n < num_iterations; ++n) {
    const uint32_t w = gen.Get(1u, kRangeMax);
    const uint32_t h = gen.Get(1u, kRangeMax);

    FrontMgrLexico mgr_lex;
    FrontMgrMax mgr_max;
    FrontMgrArea mgr_area(kAreaSize, kAreaSize);

    FrontMgrNxNBase* const mgr =
        GetMgr(front_mgr_type, &mgr_lex, &mgr_max, &mgr_area);
    CreateData(w, h, GetSnapped(front_mgr_type), &gen, &mgr->Blocks());

    ASSERT_WP2_OK(mgr->Init(ALL_RECTS, GetSnapped(front_mgr_type), w, h));
    ASSERT_WP2_OK(mgr->Sort());
    ASSERT_EQ(mgr->Blocks().size(), mgr->SizeIndices().size());
    // Check all indices are used only once.
    ASSERT_EQ(
        mgr->Blocks().size(),
        std::set<uint16_t>(mgr->SizeIndices().begin(), mgr->SizeIndices().end())
            .size());
    // Go over all indices and make sure it works.
    mgr->Clear();
    for (uint16_t i : mgr->SizeIndices()) {
      const Block& block = mgr->Blocks()[i];

      const Block max_block = mgr->GetMaxPossibleBlock();
      ASSERT_GE(block.x(), max_block.x());
      ASSERT_GE(block.y(), max_block.y());
      ASSERT_LE(block.x() + block.w(), max_block.x() + max_block.w());
      ASSERT_LE(block.y() + block.h(), max_block.y() + max_block.h());

      Block block_tmp;
      ASSERT_TRUE(mgr->TryGetNextBlock(block.dim(), &block_tmp));
      ASSERT_EQ(block, block_tmp);
      ASSERT_FALSE(mgr->Done());
      ASSERT_TRUE(mgr->UseSize(block.dim(), /*ind=*/0u, &block_tmp));
      ASSERT_EQ(block, block_tmp);
      while (mgr->UseFinal()) {
      }
    }
    ASSERT_TRUE(mgr->Done());
  }
}

// Make sure the deep copy is perfect.
TEST_P(TestFrontMgr, CopyFrom) {
  const uint32_t seed = std::get<0>(GetParam());
  const uint32_t num_iterations = std::get<1>(GetParam());
  const FrontMgrType front_mgr_type = std::get<2>(GetParam());

  UniformIntDistribution gen(seed);
  const uint32_t kRangeMax = 100;

  FrontMgrLexico mgr_lex_clone;
  FrontMgrMax mgr_max_clone;
  FrontMgrArea mgr_area_clone(kAreaSize, kAreaSize);
  FrontMgrNxNBase* const mgr_clone =
      GetMgr(front_mgr_type, &mgr_lex_clone, &mgr_max_clone, &mgr_area_clone);

  for (size_t n = 0; n < num_iterations; ++n) {
    const uint32_t w = gen.Get(1u, kRangeMax), h = gen.Get(1u, kRangeMax);

    FrontMgrLexico mgr_lex;
    FrontMgrMax mgr_max;
    FrontMgrArea mgr_area(kAreaSize, kAreaSize);
    FrontMgrNxNBase* const mgr =
        GetMgr(front_mgr_type, &mgr_lex, &mgr_max, &mgr_area);

    ASSERT_WP2_OK(mgr->Init(ALL_RECTS, GetSnapped(front_mgr_type), w, h));

    CreateData(w, h, GetSnapped(front_mgr_type), &gen, &mgr->Blocks());
    ASSERT_WP2_OK(mgr->Sort());
    mgr->Clear();  // Sort() messes up the FrontMgr state so clear it.

    int32_t num_blocks_before_cloning =
        gen.Get<int32_t>(1, mgr->Blocks().size());

    for (uint16_t i : mgr->SizeIndices()) {
      if (num_blocks_before_cloning == 0) {
        ASSERT_WP2_OK((front_mgr_type == kLexico)
                          ? mgr_lex_clone.CopyFrom(mgr_lex)
                          : (front_mgr_type == kMax)
                                ? mgr_max_clone.CopyFrom(mgr_max)
                                : mgr_area_clone.CopyFrom(mgr_area));
      }

      const Block& block = mgr->Blocks()[i];
      Block block_tmp;
      ASSERT_TRUE(mgr->TryGetNextBlock(block.dim(), &block_tmp));

      if (num_blocks_before_cloning <= 0) {
        // The clone was copied. Verify that it returns the same values as the
        // original 'map' and apply the same modifications to it from now on.
        const Block& block_clone = mgr_clone->Blocks()[i];
        Block block_tmp_clone;
        ASSERT_TRUE(
            mgr_clone->TryGetNextBlock(block_clone.dim(), &block_tmp_clone));

        ASSERT_EQ(block, block_clone);
        ASSERT_EQ(block_tmp, block_tmp_clone);
        ASSERT_EQ(mgr->Done(), mgr_clone->Done());

        ASSERT_TRUE(mgr_clone->UseSize(
            block_clone.dim(), /*ind=*/0u, &block_tmp_clone));
        while (mgr_clone->UseFinal()) {
        }
      }

      ASSERT_TRUE(mgr->UseSize(block.dim(), /*ind=*/0u, &block_tmp));
      while (mgr->UseFinal()) {
      }
      --num_blocks_before_cloning;
    }
    ASSERT_TRUE(mgr->Done());
    ASSERT_TRUE(mgr_clone->Done());
  }
}

INSTANTIATE_TEST_SUITE_P(
    TestFrontMgrInstanciation, TestFrontMgr,
    ::testing::Combine(
        /*seed=*/::testing::Values(37, 42),
        /*num_iterations=*/::testing::Values(7),
        /*front_mgr_type=*/::testing::Values(kLexico, kMax, kArea)));

//------------------------------------------------------------------------------

TEST(TestFrontMgr, Area) {
  FrontMgrArea mgr(kAreaSize, kAreaSize);
  ASSERT_WP2_OK(mgr.Init(ALL_RECTS, /*snapped=*/true, 39, 64));
  const Block blocks[] = {{0, 0, BLK_16x16}, {4, 0, BLK_16x16},  // First area
                          {0, 4, BLK_16x16}, {4, 4, BLK_16x16},

                          {8, 0, BLK_8x4},   {8, 1, BLK_8x4},  // Second area
                          {8, 2, BLK_8x8},   {8, 4, BLK_8x16},

                          {0, 8, BLK_32x32},  // Third area

                          {8, 8, BLK_8x32}};  // Fourth area

  for (const Block& block : blocks) {
    ASSERT_FALSE(mgr.Done());
    Block next_block;
    ASSERT_TRUE(mgr.UseSize(block.dim(), /*ind=*/0, &next_block));
    ASSERT_EQ(block, next_block);
    mgr.Use(block);
  }
  ASSERT_TRUE(mgr.Done());
}

TEST(TestFrontMgr, Comp) {
  const FrontMgrArea::Comp comp(kMaxTileSize / kMinBlockSizePix, kMaxBlockSize,
                                kMaxBlockSize);
  std::vector<Block> blocks = {{4, 0, BLK_16x8},  {8, 0, BLK_16x8},
                               {12, 0, BLK_16x4}, {0, 4, BLK_16x4},
                               {0, 0, BLK_16x4},  {12, 4, BLK_16x4}};
  std::sort(blocks.begin(), blocks.end(), comp);
  // Expected block order:   aaaabbbbddddeeee
  //                         ccccbbbbddddffff
  EXPECT_THAT(blocks, ::testing::ElementsAre(
                          Block{0, 0, BLK_16x4}, Block{4, 0, BLK_16x8},
                          Block{0, 4, BLK_16x4}, Block{8, 0, BLK_16x8},
                          Block{12, 0, BLK_16x4}, Block{12, 4, BLK_16x4}));
  // Expects e then f.
  auto it = std::lower_bound(blocks.begin(), blocks.end(),
                             Block{12, 0, BLK_32x32}, comp);
  EXPECT_EQ(*it, Block(12, 0, BLK_16x4));
  EXPECT_EQ(*(++it), Block(12, 4, BLK_16x4));
}

//------------------------------------------------------------------------------

}  // namespace
}  // namespace WP2
