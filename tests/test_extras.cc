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

// Test functions in extras/extras.h and src/utils/random.h

#include <random>

#include "extras/extras.h"
#include "include/helpers.h"
#include "src/utils/random.h"
#include "src/common/preview/preview.h"   // TODO(skal): make better include

namespace WP2 {
namespace {

//------------------------------------------------------------------------------

class RescaleTest : public ::testing::TestWithParam<int> {};

TEST_P(RescaleTest, ArgbBuffer) {
  static constexpr uint32_t min_size = 1, max_size = 301;
  UniformIntDistribution random(/*seed=*/(size_t)GetParam());  // Deterministic.

  // Fill ArgbBuffer randomly and rescale it.
  ArgbBuffer b;
  ASSERT_WP2_OK(
      b.Resize(random.Get(min_size, max_size), random.Get(min_size, max_size)));
  uint8_t* dst = (uint8_t*)b.GetRow(0);
  for (uint32_t j = 0; j < b.height; ++j) {
    for (uint32_t i = 0; i < 4 * b.width; ++i) {
      dst[i] = random.Get<uint8_t>(0x00u, 0xffu);
    }
    dst += b.stride;
  }
  ASSERT_WP2_OK(ArgbBufferRescale(&b, random.Get(min_size, max_size),
                                  random.Get(min_size, max_size)));
}

INSTANTIATE_TEST_SUITE_P(RescaleTestInstantiation, RescaleTest,
                         ::testing::Range(1, 20));

//------------------------------------------------------------------------------
// Make sure results of random number generation functions are always the same
// for a given seed, no matter the compiler or environment.
// Deterministic preview encoding relies on it.

TEST(RandomTest, Deterministic) {
  UniformIntDistribution random(/*seed=*/123);
  EXPECT_EQ(random.Get<uint8_t>(0x00u, 0xffu), 0xFEu);
  EXPECT_EQ(random.Get<uint32_t>(456u, 789000u), 200010u);
  EXPECT_EQ(random.Get<int16_t>(-5, 5), 3);
  EXPECT_FALSE(random.FlipACoin());

  std::vector<int> integers(10);
  std::iota(integers.begin(), integers.end(), 0);
  Shuffle(integers.begin(), integers.end(), /*seed=*/321);
  EXPECT_EQ(integers, std::vector<int>({ 5, 2, 7, 1, 0, 9, 6, 4, 3, 8 }));
}

//------------------------------------------------------------------------------

TEST(PreviewTest, TestBase64) {
  PreviewData p;
  std::string s;
  EXPECT_EQ(PreviewFromBase64("dasdascacas", &p), WP2_STATUS_BITSTREAM_ERROR);
  char b64[] = "G1XE1MCBS+nBgDcL1mnWis3SkHw=";
  EXPECT_WP2_OK(PreviewFromBase64(b64, &p));
  b64[4] = '\0';
  EXPECT_EQ(PreviewFromBase64(b64, &p), WP2_STATUS_BITSTREAM_ERROR);
}

TEST(PreviewTest, TestText) {
  PreviewData p;
  std::string in =
    "12 7 3 1\n"       // grid_width grid_height palette_size use_noise
    "8 31 1 1\n"       // colormap entries
    "16 31 7 1\n"
    "45 31 9 0\n"
    "a..........c\n"   // grid content
    ".........a..\n"
    "...........a\n"
    "............\n"
    "....b.....b.\n"
    ".a.........c\n"
    "b.....c....a\n";
  EXPECT_WP2_OK(PreviewFromText(in, &p));
  std::string out = PreviewToText(p);
  EXPECT_EQ(out, in);
  EXPECT_EQ(p.grid_width_, 12u);
  EXPECT_EQ(p.grid_height_, 7u);
  EXPECT_EQ(p.palette_.size(), 3u);
  EXPECT_EQ(p.vertices_.size(), 7u);
  EXPECT_EQ(p.corners_[0], 0);
  EXPECT_EQ(p.corners_[1], 2);
  EXPECT_EQ(p.corners_[2], 1);
  EXPECT_EQ(p.corners_[3], 0);
  EXPECT_EQ(p.use_noise_, true);

  const std::string s_out = PreviewToBase64(p);
  EXPECT_TRUE(!s_out.empty());
  EXPECT_WP2_OK(PreviewFromBase64(s_out, &p));
  printf("%s", PrintPreview(p, /* reduced= */false).c_str());
  EXPECT_EQ(p.grid_width_, 12u);
  EXPECT_EQ(p.grid_height_, 7u);
  EXPECT_EQ(p.palette_.size(), 3u);
  EXPECT_EQ(p.vertices_.size(), 7u);
  EXPECT_EQ(p.corners_[0], 0);
  EXPECT_EQ(p.corners_[1], 2);
  EXPECT_EQ(p.corners_[2], 1);
  EXPECT_EQ(p.corners_[3], 0);
  EXPECT_EQ(p.use_noise_, true);
}

//------------------------------------------------------------------------------

}  // namespace
}  // namespace WP2
