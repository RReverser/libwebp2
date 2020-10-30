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

// Color precision test. Ensures loss is within expected bounds and that
// conversion Argb <-> YCoCg is stable at some point.
// Warning: These tests can be easily optimized out. Compile with Debug and -O0

#include <algorithm>

#include "include/helpers.h"
#include "src/common/color_precision.h"

namespace WP2 {
namespace {

//------------------------------------------------------------------------------

TEST(ColorPrecisionTest, Argb32b_Argb38b) {
  for (uint32_t v = 0; v < (1u << 8); ++v) {
    const Argb32b original{(uint8_t)v, (uint8_t)v, (uint8_t)v, (uint8_t)v};
    const Argb32b result = ToArgb32b(ToArgb38b(original));
    ASSERT_EQ(original.a, result.a);
    ASSERT_EQ(original.r, result.r);
    ASSERT_EQ(original.g, result.g);
    ASSERT_EQ(original.b, result.b);
  }

  for (uint32_t v = 0; v < (1u << 10); ++v) {
    const Argb38b original{(uint8_t)v, (uint8_t)v, (uint8_t)v, (uint8_t)v};
    const Argb38b result = ToArgb38b(ToArgb32b(original));
    ASSERT_LE(original.a, result.a);
    static constexpr int32_t kMaxLoss = (1u << 2) - 1u;
    ASSERT_LE(std::abs((int32_t)original.r - result.r), kMaxLoss);
    ASSERT_LE(std::abs((int32_t)original.g - result.g), kMaxLoss);
    ASSERT_LE(std::abs((int32_t)original.b - result.b), kMaxLoss);
  }
}

TEST(ColorPrecisionTest, Argb32b_RGB12b) {
  for (uint32_t v = 0; v < (1u << 4); ++v) {
    const RGB12b original{(uint8_t)v, (uint8_t)v, (uint8_t)v};
    const RGB12b result = ToRGB12b(ToArgb32b(original));
    ASSERT_EQ(original.r, result.r);
    ASSERT_EQ(original.g, result.g);
    ASSERT_EQ(original.b, result.b);
  }

  for (uint32_t v = 0; v < (1u << 8); ++v) {
    const Argb32b original{255u, (uint8_t)v, (uint8_t)v, (uint8_t)v};
    const Argb32b result = ToArgb32b(ToRGB12b(original));
    ASSERT_LE(original.a, result.a);
    static constexpr int32_t kMaxLoss = (1u << 4) - 1u;
    ASSERT_LE(std::abs((int32_t)original.r - result.r), kMaxLoss);
    ASSERT_LE(std::abs((int32_t)original.g - result.g), kMaxLoss);
    ASSERT_LE(std::abs((int32_t)original.b - result.b), kMaxLoss);
  }
}

//------------------------------------------------------------------------------

TEST(ColorPrecisionTest, Argb32b_AYCoCg19b) {
  uint32_t total = 0;
  for (uint32_t r = 0; r < (1u << 8); ++r) {
    for (uint32_t g = 0; g < (1u << 8); ++g) {
      for (uint32_t b = 0; b < (1u << 8); ++b) {
        const Argb32b original{255u, (uint8_t)r, (uint8_t)g, (uint8_t)b};
        const Argb32b result = ToArgb32b(ToAYCoCg19b(original));
        const uint32_t err_a = std::abs((int32_t)original.a - result.a);
        const uint32_t err_r = std::abs((int32_t)original.r - result.r);
        const uint32_t err_g = std::abs((int32_t)original.g - result.g);
        const uint32_t err_b = std::abs((int32_t)original.b - result.b);
        static constexpr uint32_t kMaxLossAlpha = (1u << 7) - 1u;
        static constexpr uint32_t kMaxLoss = (1u << 3);
        ASSERT_LE(err_a, kMaxLossAlpha);
        ASSERT_LE(err_r, kMaxLoss);
        ASSERT_LE(err_g, kMaxLoss);
        ASSERT_LE(err_b, kMaxLoss);
        total += err_r + err_g + err_b;
        const Argb32b stable_value = ToArgb32b(ToAYCoCg19b(result));
        const Argb32b stable_result = ToArgb32b(ToAYCoCg19b(stable_value));
        ASSERT_EQ(stable_value.r, stable_result.r);
        ASSERT_EQ(stable_value.g, stable_result.g);
        ASSERT_EQ(stable_value.b, stable_result.b);
      }
    }
  }
  printf("Total error: %.3f\n", total / (256 * 256 * 256.f));
}

//------------------------------------------------------------------------------

TEST(ColorPrecisionTest, Argb32b_UInt8) {
  uint8_t argb[4];
  for (uint32_t v = 0; v < (1u << 8); ++v) {
    const Argb32b original{(uint8_t)v, (uint8_t)v, (uint8_t)v, (uint8_t)v};
    ToUInt8(original, argb);
    const Argb32b result = ToArgb32b(argb);
    ASSERT_EQ(original.a, result.a);
    ASSERT_EQ(original.r, result.r);
    ASSERT_EQ(original.g, result.g);
    ASSERT_EQ(original.b, result.b);
  }
}

//------------------------------------------------------------------------------

TEST(ColorPrecisionTest, CheckPremultiplied) {
  EXPECT_TRUE(CheckPremultiplied(Argb32b{255, 255, 255, 255}));
  EXPECT_TRUE(CheckPremultiplied(Argb32b{255, 12, 255, 58}));
  EXPECT_TRUE(CheckPremultiplied(Argb32b{159, 45, 159, 78}));
  EXPECT_TRUE(CheckPremultiplied(Argb32b{89, 45, 23, 59}));
  EXPECT_TRUE(CheckPremultiplied(Argb32b{0, 0, 0, 0}));
  EXPECT_TRUE(CheckPremultiplied(ToArgb38b(Argb32b{255, 255, 255, 255})));
  EXPECT_TRUE(CheckPremultiplied(ToArgb38b(Argb32b{255, 12, 255, 58})));
  EXPECT_TRUE(CheckPremultiplied(ToArgb38b(Argb32b{159, 45, 159, 78})));
  EXPECT_TRUE(CheckPremultiplied(ToArgb38b(Argb32b{89, 45, 23, 59})));
  EXPECT_TRUE(CheckPremultiplied(ToArgb38b(Argb32b{0, 0, 0, 0})));

  EXPECT_FALSE(CheckPremultiplied(Argb32b{236, 24, 247, 12}));
  EXPECT_FALSE(CheckPremultiplied(Argb32b{56, 48, 57, 57}));
  EXPECT_FALSE(CheckPremultiplied(Argb32b{0, 24, 5, 0}));
  EXPECT_FALSE(CheckPremultiplied(ToArgb38b(Argb32b{236, 24, 247, 12})));
  EXPECT_FALSE(CheckPremultiplied(ToArgb38b(Argb32b{56, 48, 57, 57})));
  EXPECT_FALSE(CheckPremultiplied(ToArgb38b(Argb32b{0, 24, 5, 0})));
}

}  // namespace
}  // namespace WP2
