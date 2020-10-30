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

// RectangleGroup test.

#include "include/helpers.h"
#include "src/enc/anim/anim_enc.h"

namespace WP2 {
namespace {

//------------------------------------------------------------------------------

bool Equal(const Rectangle& a, const Rectangle& b) {
  return a.x == b.x && a.y == b.y && a.width == b.width && a.height == b.height;
}

constexpr uint32_t kMaxNumRectangles = 4;

//------------------------------------------------------------------------------

TEST(RectangleGroupTest, Square) {
  RectangleGroup rectangle_group(kMaxNumRectangles);
  rectangle_group.AddPoint(0, 0);
  rectangle_group.AddPoint(1, 0);
  rectangle_group.AddPoint(0, 1);
  rectangle_group.AddPoint(1, 1);

  ASSERT_EQ(rectangle_group.GetNumRectangles(), 1u);
  EXPECT_TRUE(Equal(rectangle_group.GetRectangle(0), {0, 0, 2, 2}));
}

//------------------------------------------------------------------------------

TEST(RectangleGroupTest, Points) {
  RectangleGroup rectangle_group(kMaxNumRectangles);
  rectangle_group.AddPoint(0, 0);
  rectangle_group.AddPoint(2, 0);

  ASSERT_EQ(rectangle_group.GetNumRectangles(), 2u);
  EXPECT_TRUE(Equal(rectangle_group.GetRectangle(0), {0, 0, 1, 1}));
  EXPECT_TRUE(Equal(rectangle_group.GetRectangle(1), {2, 0, 1, 1}));
}

//------------------------------------------------------------------------------

TEST(RectangleGroupTest, Cross) {
  RectangleGroup rectangle_group(kMaxNumRectangles);
  for (uint32_t x = 0; x < 10; ++x) rectangle_group.AddPoint(x, 5);
  for (uint32_t y = 0; y < 10; ++y) rectangle_group.AddPoint(5, y);
  EXPECT_EQ(rectangle_group.GetCumulativeArea(), 10u + 9u);
}

//------------------------------------------------------------------------------

TEST(RectangleGroupTest, ScatteredLines) {
  RectangleGroup rectangle_group(kMaxNumRectangles);
  uint32_t num_points = 0;
  for (uint32_t x : {0u, 1u, 2u, 3u, 5u, 7u, 8u, 1000u}) {
    rectangle_group.AddPoint(x, 0);
    ++num_points;
  }
  ASSERT_EQ(rectangle_group.GetNumRectangles(), 4u);
  EXPECT_EQ(rectangle_group.GetCumulativeArea(), num_points);
}

//------------------------------------------------------------------------------

TEST(RectangleGroupTest, ScatteredRectangles) {
  RectangleGroup rectangle_group(kMaxNumRectangles);
  uint32_t num_points = 0;
  for (uint32_t x : {101u, 100u, 20u, 21u, 0u, 99u, 102u, 1u, 1000u}) {
    for (uint32_t y : {0u, 1u, 2u}) {
      rectangle_group.AddPoint(x, y);
      ++num_points;
    }
  }
  ASSERT_EQ(rectangle_group.GetNumRectangles(), 4u);
  EXPECT_EQ(rectangle_group.GetCumulativeArea(), num_points);
}

//------------------------------------------------------------------------------

}  // namespace
}  // namespace WP2
