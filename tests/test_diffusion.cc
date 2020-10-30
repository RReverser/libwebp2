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

#include "include/helpers.h"
#include "src/common/lossy/context.h"
#include "src/enc/wp2_enc_i.h"

namespace WP2 {
namespace {

//------------------------------------------------------------------------------

// Block layout: eleven blocks (where a = 10 and b = 11).
// 1111222233333333
// 1111222233333333
// 4444222255555555
// 4444222255555555
// 4444667755555555
// 4444667755555555
// 4444897755555555
// 4444ab7755555555
const Block blk1 = {0, 0, BLK_16x16};
const Block blk2 = {4, 0, BLK_16x32};
const Block blk3 = {8, 0, BLK_32x16};
const Block blk4 = {0, 4, BLK_16x32};
const Block blk5 = {8, 4, BLK_32x32};
const Block blk6 = {4, 8, BLK_8x8};
const Block blk7 = {6, 8, BLK_8x16};
const Block blk8 = {4, 10, BLK_4x4};
const Block blk9 = {5, 10, BLK_4x4};
const Block blk10 = {4, 11, BLK_4x4};
const Block blk11 = {5, 11, BLK_4x4};

const uint32_t kWidth = 64, kHeight = 64;  // Dimensions in pixels.
const uint32_t kMaxStrength = 255;

//------------------------------------------------------------------------------

TEST(TestDCDiffusionMap, Test) {
  DCDiffusionMap map;
  ASSERT_WP2_OK(map.Init(kWidth));

  EXPECT_EQ(map.Get(Block(10, 30, BLK_8x8), kMaxStrength), 0);

  FrontMgrDefault mgr;
  ASSERT_WP2_OK(mgr.Init(ALL_RECTS, /*snapped=*/false, kWidth, kHeight));

  EXPECT_EQ(map.Get(Block(0, 0, BLK_8x8), kMaxStrength), 0);
  map.Store(mgr, blk1, 10);
  mgr.Use(blk1);

  // Half of blk1 error
  EXPECT_EQ(map.Get(blk2, kMaxStrength), std::round(10.f / 2.f));
  map.Store(mgr, blk2, 20);
  mgr.Use(blk2);

  // 1/4th of blk2 error
  EXPECT_EQ(map.Get(blk3, kMaxStrength), std::round(20.f / 4.f));
  EXPECT_EQ(map.Get(blk3, kMaxStrength / 2 + 1), std::round(20.f / 4.f / 2.f));
  EXPECT_EQ(map.Get(blk3, 0), 0);
  map.Store(mgr, blk3, 30);
  mgr.Use(blk3);

  // Half of blk1 + 1/4th of blk2 error
  EXPECT_EQ(map.Get(blk4, kMaxStrength), std::round(10.f / 2.f + 20.f / 4.f));
  map.Store(mgr, blk4, 40);
  mgr.Use(blk4);

  // etc.
  EXPECT_EQ(map.Get(blk5, kMaxStrength), std::round(20.f / 4.f + 30.f));
  map.Store(mgr, blk5, 50);
  mgr.Use(blk5);

  EXPECT_EQ(map.Get(blk6, kMaxStrength), std::round(40.f / 4.f + 20.f / 8.f));
  map.Store(mgr, blk6, 60);
  mgr.Use(blk6);

  EXPECT_EQ(map.Get(blk7, kMaxStrength),
            std::round(20.f / 8.f + 50.f / 3.f + 60.f / 2.f));
  map.Store(mgr, blk7, 70);
  mgr.Use(blk7);

  EXPECT_EQ(map.Get(blk8, kMaxStrength), std::round(40.f / 8.f + 60.f / 4.f));
  map.Store(mgr, blk8, 80);
  mgr.Use(blk8);

  EXPECT_EQ(map.Get(blk9, kMaxStrength),
            std::round(60.f / 4.f + 70.f / 4.f + 80.f / 2.f));
  map.Store(mgr, blk9, 90);
  mgr.Use(blk9);

  EXPECT_EQ(map.Get(blk10, kMaxStrength), std::round(40.f / 8.f + 80.f / 2.f));
  map.Store(mgr, blk10, 100);
  mgr.Use(blk10);

  EXPECT_EQ(map.Get(blk11, kMaxStrength),
            std::round(70.f / 4.f + 90.f + 100.f / 2.f));
  map.Store(mgr, blk11, 110);
  mgr.Use(blk11);

  map.Clear();
}

//------------------------------------------------------------------------------

// Same test as above but a clone is created in the middle.
TEST(TestDCDiffusionMap, CopyFrom) {
  const Block kBlkBeforeCloning[] = {blk1, blk2, blk3, blk4, blk5};
  const Block kBlkAfterCloning[] = {blk6, blk7, blk8, blk9, blk10, blk11};

  DCDiffusionMap map_clone;
  FrontMgrDefault mgr_clone;
  int16_t error_clone;

  std::vector<int16_t> map_results;
  {
    DCDiffusionMap map;
    ASSERT_WP2_OK(map.Init(kWidth));
    FrontMgrDefault mgr;
    ASSERT_WP2_OK(mgr.Init(ALL_RECTS, /*snapped=*/false, kWidth, kHeight));

    int16_t error = 10;
    for (const Block& blk : kBlkBeforeCloning) {
      map.Store(mgr, blk, error);
      mgr.Use(blk);
      error += 10;
    }

    ASSERT_WP2_OK(map_clone.CopyFrom(map));
    ASSERT_WP2_OK(mgr_clone.CopyFrom(mgr));
    error_clone = error;

    for (const Block& blk : kBlkAfterCloning) {
      map_results.push_back(map.Get(blk, kMaxStrength));
      map.Store(mgr, blk, error);
      mgr.Use(blk);
      error += 10;
    }
  }
  // Verify that the 'map_results' match with the 'map_clone' now that the
  // original 'map' was deleted.

  uint32_t i = 0;
  for (const Block& blk : kBlkAfterCloning) {
    EXPECT_EQ(map_results[i++], map_clone.Get(blk, kMaxStrength)) << i;
    map_clone.Store(mgr_clone, blk, error_clone);
    mgr_clone.Use(blk);
    error_clone += 10;
  }
}

//------------------------------------------------------------------------------

}  // namespace
}  // namespace WP2
