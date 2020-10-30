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

// Test some DSP functions in dsp/dsp.h

#include <cstdio>
#include <string>

#include "examples/example_utils.h"
#include "include/helpers.h"
#include "src/dsp/math.h"
#include "src/dsp/dsp.h"
#include "src/utils/plane.h"
#include "src/utils/random.h"

namespace WP2 {
namespace {

class DspTest : public ::testing::TestWithParam<WP2CPUInfo> {
  void SetUp() override {
    WP2DspReset();
    WP2GetCPUInfo = GetParam();
  }
};

TEST_P(DspTest, SSETest16s) {
  WP2PSNRInit();
  testing::PrecalculatedRandom<65536u, int16_t> rnd(-1024, 1024);
  uint64_t sum = 0;
  constexpr uint32_t max_h = 128;
  constexpr uint32_t max_w = 128;
  constexpr uint32_t stride1 = max_w + 53, stride2 = max_w + 39;
  for (uint32_t N = 0; N < 5000; ++N) {
    const uint32_t h = 1 + rnd.GetUnsigned(max_w);
    for (uint32_t w = 4; w <= max_w; w += 4) {
      int16_t tmp1[stride1 * max_h], tmp2[stride2 * max_h];
      rnd.Fill(tmp1, w, h, stride1);
      rnd.Fill(tmp2, w, h, stride2);
      sum = ((sum << 2) + (sum << 3))
          ^ WP2SumSquaredErrorBlock(tmp1, stride1, tmp2, stride2, w, h);
    }
  }
  EXPECT_EQ(sum, 2548820040083235724llu);
}

TEST_P(DspTest, SSETest8u) {
  WP2PSNRInit();
  testing::PrecalculatedRandom<65536u, uint8_t> rnd(0, 256);
  uint64_t sum = 0;
  constexpr uint32_t max_w = 128;
  for (uint32_t N = 0; N < 3000000; ++N) {
    const uint32_t w = 1 + rnd.GetUnsigned(max_w);
    uint8_t tmp1[max_w], tmp2[max_w];
    rnd.Fill(tmp1, w);
    rnd.Fill(tmp2, w);
    sum = ((sum >> 3) + (sum << 2))
        ^ WP2SumSquaredError8u(tmp1, tmp2, w);
  }
  EXPECT_EQ(sum, 16994877026985094642llu);
}

TEST_P(DspTest, SSETest4x8u) {
  WP2PSNRInit();
  testing::PrecalculatedRandom<65536u, uint8_t> rnd(0, 256);
  constexpr uint32_t max_w = 128;
  uint64_t sums[4] = { 0 };
  for (uint32_t N = 0; N < 2000000; ++N) {
    const uint32_t w = 1 + rnd.GetUnsigned(max_w);
    uint8_t tmp1[4 * max_w], tmp2[4 * max_w];
    rnd.Fill(tmp1, 4 * w);
    rnd.Fill(tmp2, 4 * w);
    WP2SumSquaredError4x8u(tmp1, tmp2, w, sums);
  }
  const uint64_t sum = sums[0] ^ (sums[1] + sums[2] - sums[3]);
  EXPECT_EQ(sum, 211731994llu);
}

//------------------------------------------------------------------------------

TEST_P(DspTest, ANS) {
  ANSInit();
  const uint16_t kCDFBase[] = {0, 1, 2, 3, 4, 5, 6, 7};
  const uint16_t kCDFVar[] = {0,     16377, 16377, 16377,
                              16377, 16377, 16377, 16377};
  uint16_t cumul[] = {0, 2340, 4681, 7021, 9362, 11702, 14043, 16384};
  ANSUpdateCDF(/*n=*/8, kCDFBase, kCDFVar, /*mult=*/12345, cumul);
  EXPECT_THAT(cumul, ::testing::ElementsAre(0, 4984, 6884, 8783, 10684, 12583,
                                            14483, 16384));
}

// Special case where 'mult' is 65535 meaning factor is 1.
TEST_P(DspTest, ANSCopy) {
  ANSInit();
  const uint16_t kCDFBase[] = {0, 1, 2, 3, 4, 5, 6, 7};
  const uint16_t kCDFVar[] = {0,     16377, 16377, 16377,
                              16377, 16377, 16377, 16377};
  uint16_t cumul[] = {0, 2340, 4681, 7021, 9362, 11702, 14043, 16384};
  ANSUpdateCDF(/*n=*/7, kCDFBase, kCDFVar, /*mult=*/65535, cumul);
  const uint16_t cumul_to_check[] = {cumul[0], cumul[1], cumul[2], cumul[3],
                                     cumul[4], cumul[5], cumul[6]};
  EXPECT_THAT(cumul_to_check, ::testing::ElementsAre(0, 16378, 16379, 16380,
                                                     16381, 16382, 16383));
}

//------------------------------------------------------------------------------

TEST_P(DspTest, Alpha) {
  WP2AlphaInit();
  testing::PrecalculatedRandom<65536u, uint8_t> rnd8(0, 256);
  testing::PrecalculatedRandom<65536u, int16_t> rnd16(-32768, 32768);
  const uint32_t kMaxLen = 500;
  std::vector<uint8_t> tmp8(kMaxLen, 0);
  std::vector<int16_t> tmp16(kMaxLen, 0);
  // simple 8b tests first
  for (uint32_t value8 = 0; value8 < 256; ++value8) {
    for (auto& v : tmp8) v = value8 ^ 255;
    EXPECT_FALSE(WP2HasValue8b(tmp8.data(), kMaxLen, value8));
    EXPECT_TRUE(WP2HasOtherValue8b(tmp8.data(), kMaxLen, value8));
    EXPECT_TRUE(WP2HasOtherValue8b32b(tmp8.data(), kMaxLen / 4, value8));

    tmp8[rnd8.GetUnsigned(kMaxLen)] = value8;
    EXPECT_TRUE(WP2HasValue8b(tmp8.data(), kMaxLen, value8));
    EXPECT_TRUE(WP2HasOtherValue8b(tmp8.data(), kMaxLen, value8));

    tmp8[rnd8.GetUnsigned(kMaxLen) & ~3] = value8;
    EXPECT_TRUE(WP2HasOtherValue8b32b(tmp8.data(), kMaxLen / 4, value8));

    for (auto& v : tmp8) v = value8;
    EXPECT_TRUE(WP2HasValue8b(tmp8.data(), kMaxLen, value8));
    EXPECT_FALSE(WP2HasOtherValue8b(tmp8.data(), kMaxLen, value8));
    EXPECT_FALSE(WP2HasOtherValue8b32b(tmp8.data(), kMaxLen / 4, value8));
  }
  // also 8b-in-a-32b-ARGB-sample
  for (uint32_t value8 = 0; value8 < 256; ++value8) {
    for (auto& v : tmp8) v = value8 ^ 255;
    EXPECT_TRUE(WP2HasOtherValue8b32b(tmp8.data(), kMaxLen / 4, value8));

    for (uint32_t c : { 1, 2, 3 }) {
      for (uint32_t i = 0; i < kMaxLen / 4; ++i) tmp8[4 * i + c] = value8;
      EXPECT_TRUE(WP2HasOtherValue8b32b(tmp8.data(), kMaxLen / 4, value8));
    }
    tmp8[4 * rnd8.GetUnsigned(kMaxLen / 4) + 0] = value8;
    EXPECT_TRUE(WP2HasOtherValue8b32b(tmp8.data(), kMaxLen / 4, value8));

    for (uint32_t i = 0; i < kMaxLen; ++i) {
      tmp8[i] = ((i % 4) == 0) ? value8 : value8 ^ 0xff;
    }
    EXPECT_FALSE(WP2HasOtherValue8b32b(tmp8.data(), kMaxLen / 4, value8));

    tmp8[4 * rnd8.GetUnsigned(kMaxLen / 4) + 0] = value8 ^ 0xff;
    EXPECT_TRUE(WP2HasOtherValue8b32b(tmp8.data(), kMaxLen / 4, value8));
  }
  // simple 16b tests
  for (int32_t value16 = -32768; value16 < 32768; ++value16) {
    for (auto& v : tmp16) v = value16 ^ 255;
    EXPECT_FALSE(WP2HasValue16b(tmp16.data(), kMaxLen, value16));
    EXPECT_TRUE(WP2HasOtherValue16b(tmp16.data(), kMaxLen, value16));

    tmp16[rnd16.GetUnsigned(kMaxLen)] = value16;
    EXPECT_TRUE(WP2HasValue16b(tmp16.data(), kMaxLen, value16));
    EXPECT_TRUE(WP2HasOtherValue16b(tmp16.data(), kMaxLen, value16));

    for (auto& v : tmp16) v = value16;
    EXPECT_TRUE(WP2HasValue16b(tmp16.data(), kMaxLen, value16));
    EXPECT_FALSE(WP2HasOtherValue16b(tmp16.data(), kMaxLen, value16));
  }

  // then random tests
  for (uint32_t N = 0; N < 50000; ++N) {
    const uint8_t value8 = (uint8_t)(N % 256u);
    const int16_t value16 = (int16_t)(N % 65536u) - 32768;
    const uint32_t len = 100 + rnd8.GetUnsigned(kMaxLen - 100);
    rnd8.Fill(tmp8.data(), len);
    rnd16.Fill(tmp16.data(), len);
    for (auto& v : tmp8) v = (v == value8) ? value8 ^ 0xff : v;
    for (auto& v : tmp16) v = (v == value16) ? value8 ^ 0xffff : v;
    EXPECT_FALSE(WP2HasValue8b(tmp8.data(), len, value8));
    EXPECT_FALSE(WP2HasValue16b(tmp16.data(), len, value16));
    EXPECT_TRUE(WP2HasOtherValue8b(tmp8.data(), len, value8));
    EXPECT_TRUE(WP2HasOtherValue16b(tmp16.data(), len, value16));
  }
}

//------------------------------------------------------------------------------

// Creates a 's x s' pattern with some 'v' elements aligned in direction of each
// of 8 angles from 45 to -112.5 degrees.
template <int s>
void Fill8DirPatterns(int16_t patterns[kDrctFltNumDirs][s][s], int v) {
  for (int i = 0; i < s; ++i) {
    patterns[0][/*y=*/i][/*x=*/s - 1 - i] = v;              // 45 degrees
    patterns[1][/*y=*/s / 4 + i / 2][/*x=*/s - 1 - i] = v;  // 22.5
    patterns[2][/*y=*/s / 2][/*x=*/i] = v;                  // horizontal
    patterns[3][/*y=*/s / 4 + i / 2][/*x=*/i] = v;          // -22.5
    patterns[4][/*y=*/i][/*x=*/i] = v;                      // -45
    patterns[5][/*y=*/i][/*x=*/s / 4 + i / 2] = v;          // -67.5
    patterns[6][/*y=*/i][/*x=*/s / 2] = v;                  // vertical
    patterns[7][/*y=*/s - 1 - i][/*x=*/s / 4 + i / 2] = v;  // -112.5
  }
}

TEST_P(DspTest, CdefDirection4x4) {
  DrctFilterInit();
  const uint32_t kNumPrecisionBits = 10;
  int16_t patterns[kDrctFltNumDirs + 1][4][4] = {{{0}}};    // direction, y, x
  Fill8DirPatterns<4>(patterns, (1 << (kNumPrecisionBits - 1)) - 1);

  for (uint32_t d = 0; d <= kDrctFltNumDirs; ++d) {
    uint32_t direction, certainty;
    CdefDirection4x4(&patterns[d][0][0], /*step=*/4, kNumPrecisionBits,
                     &direction, &certainty);
    if (d < kDrctFltNumDirs) {
      ASSERT_EQ(direction, d);
      EXPECT_GE(certainty, 63u);  // Maximum confidence.
    } else {
      EXPECT_EQ(certainty, 0u);  // Expect uniform pixels.
    }
  }
}

TEST_P(DspTest, CdefDirection8x8) {
  DrctFilterInit();
  const uint32_t kNumPrecisionBits = 10;
  int16_t patterns[kDrctFltNumDirs + 1][8][8] = {{{0}}};    // direction, y, x
  Fill8DirPatterns<8>(patterns, (1 << (kNumPrecisionBits - 1)) - 1);

  for (uint32_t d = 0; d <= kDrctFltNumDirs; ++d) {
    uint32_t direction, certainty;
    CdefDirection8x8(&patterns[d][0][0], /*step=*/8, kNumPrecisionBits,
                     &direction, &certainty);
    if (d < kDrctFltNumDirs) {
      ASSERT_EQ(direction, d);
      EXPECT_GE(certainty, 63u);  // Maximum confidence.
    } else {
      EXPECT_EQ(certainty, 0u);  // Expect uniform pixels.
    }
  }
}

//------------------------------------------------------------------------------

INSTANTIATE_TEST_SUITE_P(DspTestInstantiation, DspTest,
    ::testing::ValuesIn(testing::kWP2CpuInfos));

}  // namespace
}  // namespace WP2
