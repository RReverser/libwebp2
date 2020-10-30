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

// Test example utils.

#include <string>
#include <vector>
#include <random>
#include <deque>

#include "examples/example_utils.h"
#include "include/helpers.h"
#include "src/utils/random.h"

namespace {

//------------------------------------------------------------------------------

TEST(PathTest, Simple) {
  ASSERT_TRUE(WP2IsDirectory(WP2::testing::GetTestDataPath("")));
  ASSERT_FALSE(WP2IsDirectory(WP2::testing::GetTestDataPath("source1.png")));

#if defined(_WIN32)
  ASSERT_EQ(WP2JoinPath("path", "file.png"), "path\\file.png");
#else
  ASSERT_EQ(WP2JoinPath("path", "file.png"), "path/file.png");
#endif
  ASSERT_EQ(WP2GetFileName("path/file.png"), "file.png");
  ASSERT_EQ(WP2RemoveFileExtension("file.png"), "file");
}

//------------------------------------------------------------------------------

std::string RandomString(WP2::UniformIntDistribution* const r, size_t len,
                         int density = 0) {
  std::string s(len + 1, 0);
  for (auto& c : s) {
    c = (r->Get(0, 100) < density) ? ' ' : (char)r->Get('a', 'z');
  }
  s[len] = '\0';
  return s;
}

TEST(ProgramOptionTest, Display) {
  ProgramOptions opt;
  WP2::UniformIntDistribution rnd(21565);
  opt.Add("");
  opt.Add("test", "");
  std::deque<std::string> flags;  // permanent storage for flag strings
  for (uint32_t i = 0; i < 200; ++i) {
    const int density = rnd.Get(0, 100);
    flags.push_back(RandomString(&rnd, rnd.Get(0, 50), 10).data());
    const std::string desc = RandomString(&rnd, rnd.Get(0, 200), density);
    opt.Add(&flags.back()[0], desc);
    opt.Add("a", RandomString(&rnd, 79, density));
    opt.Add("b", RandomString(&rnd, 80, density));
    opt.Add("c", RandomString(&rnd, 81, density));
    opt.Add("glah", RandomString(&rnd, 100, 101));
  }
  opt.Add("test", "one  two   three    four");
  opt.Add(
      "a",
      "hgjfkaloirhgjfkaloirhgjfkaloirhgjfkaloirhgjfkaloirhgjfkaloirhgjfkaloir");
  opt.Add("a", "h");
  opt.AddSystemOptionSection();
  opt.Print();
}

//------------------------------------------------------------------------------

static void CheckAndReset(float* const values, uint32_t num, uint32_t max,
                          bool* const error, bool error_value) {
  for (uint32_t n = num; n < max; ++n) EXPECT_EQ(values[n], -1.f) << " #" << n;
  for (uint32_t n = 0; n < max; ++n) values[n] = -1.f;
  EXPECT_EQ(*error, error_value);
  *error = false;
}

TEST(ArgParser, Floats) {
  constexpr uint32_t kNumValues = 4;
  bool error = false;
  float values[kNumValues];
  CheckAndReset(values, kNumValues, kNumValues, &error, false);

  // Some corner cases
  EXPECT_EQ(ExUtilGetFloats(nullptr, values, kNumValues, &error, ' '), 0u);
  CheckAndReset(values, 0, kNumValues, &error, true);

  EXPECT_EQ(ExUtilGetFloats("", values, kNumValues, &error, ' '), 0u);
  CheckAndReset(values, 0, kNumValues, &error, false);

  EXPECT_EQ(ExUtilGetFloats("", values, kNumValues, nullptr, ' '), 0u);
  CheckAndReset(values, 0, kNumValues, &error, false);

  EXPECT_EQ(ExUtilGetFloats(" ", values, kNumValues, nullptr, ' '), 0u);
  CheckAndReset(values, 0, kNumValues, &error, false);

  EXPECT_EQ(ExUtilGetFloats(" ", values, kNumValues, nullptr, ' '), 0u);
  CheckAndReset(values, 0, kNumValues, &error, false);

  // some proper parsing
  EXPECT_EQ(ExUtilGetFloats("-2.54", values, kNumValues, &error, ' '), 1u);
  EXPECT_EQ(values[0], -2.54f);
  CheckAndReset(values, 1u, kNumValues, &error, false);

  EXPECT_EQ(ExUtilGetFloats("-2.54,312", values, kNumValues, &error, ','), 2u);
  EXPECT_EQ(values[0], -2.54f);
  EXPECT_EQ(values[1], 312.f);
  CheckAndReset(values, 2u, kNumValues, &error, false);

  EXPECT_EQ(ExUtilGetFloats(" -2.54,3, -2., 4",
                            values, kNumValues, &error, ','), 4u);
  EXPECT_EQ(values[0], -2.54f);
  EXPECT_EQ(values[1], 3.f);
  EXPECT_EQ(values[2], -2.f);
  EXPECT_EQ(values[3], 4.f);
  CheckAndReset(values, 4u, kNumValues, &error, false);

  EXPECT_EQ(ExUtilGetFloats(" -2.4| 3.2 |-2.1 4|",
                            values, kNumValues, &error, '|'), 3u);
  EXPECT_EQ(values[0], -2.4f);
  EXPECT_EQ(values[1], 3.2f);
  EXPECT_EQ(values[2], -2.1f);
  CheckAndReset(values, 4u, kNumValues, &error, false);

  EXPECT_EQ(ExUtilGetFloats("0,-1,2,-3,4,-5,6,-7,1,1,1,1,1,1,1,1,1",
                            values, kNumValues, &error, ','), kNumValues);
  EXPECT_EQ(values[0], 0.f);
  EXPECT_EQ(values[1], -1.f);
  EXPECT_EQ(values[2], 2.f);
  EXPECT_EQ(values[3], -3.f);
  CheckAndReset(values, kNumValues, kNumValues, &error, false);
}

//------------------------------------------------------------------------------

}  // namespace
