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

// Test speed of some DSP functions

#include <memory>

#include "examples/example_utils.h"
#include "include/helpers.h"
#include "src/utils/plane.h"
#include "src/utils/random.h"

namespace WP2 {
namespace {

class SpeedTestCpuInfo : public ::testing::TestWithParam<WP2CPUInfo> {
  void SetUp() override {
    WP2DspReset();
    WP2GetCPUInfo = GetParam();
    WP2DecDspInit();
    WP2EncDspInit();
    WP2SSIMInit();
    WP2PSNRInit();
  }
};

//------------------------------------------------------------------------------

TEST_P(SpeedTestCpuInfo, TestAddSub) {
  constexpr uint32_t kMaxLen = 32 * 32;
  int16_t src[kMaxLen], pred[kMaxLen], dst[kMaxLen];
  int32_t res[kMaxLen];
  constexpr uint32_t kNumTests = 100;
  constexpr uint32_t kNumLoops = 10000;
  PseudoRNG rng;
  for (uint32_t t = 0; t < kNumTests; ++t) {
    const uint32_t len = (rng.GetUnsigned(kMaxLen - 1) + 1) & ~3;
    const int32_t min = -(1 + rng.GetUnsigned(1023u));
    const int32_t max =  (1 + rng.GetUnsigned(1023u));
    testing::PrecalculatedRandom<kMaxLen * 3, int16_t> rnd16(-2048, 2048);
    for (uint32_t l = 0; l < kNumLoops; ++l) {
      rnd16.Fill(src, len);
      rnd16.Fill(pred, len);
      WP2::SubtractRow(src, pred, res, len);
      WP2::AddRow(pred, res, min, max, dst, len);
    }
    uint32_t nb_errors = 0;
    for (uint32_t i = 0; i < len; ++i) {
      if (src[i] <= min) {
        nb_errors += (dst[i] != min);
      } else if (src[i] >= max) {
        nb_errors += (dst[i] != max);
      } else {
        nb_errors += (src[i] != dst[i]);
      }
    }
    EXPECT_EQ(nb_errors, 0u) << " rate: " << (100 * nb_errors / len) << "%";
  }
}

TEST_P(SpeedTestCpuInfo, TestAddSubBlock) {
  constexpr uint32_t kMaxLen = 32 * 32;
  int16_t src[kMaxLen], pred[kMaxLen], dst[kMaxLen];
  int32_t res[kMaxLen];
  constexpr uint32_t kNumTests = 1000;
  constexpr uint32_t kNumLoops = 200;
  PseudoRNG rng;
  for (uint32_t t = 0; t < kNumTests; ++t) {
    const uint32_t height = 1 + rng.GetUnsigned(32 - 1);
    const int32_t min = -(1 + rng.GetUnsigned(1023u));
    const int32_t max =  (1 + rng.GetUnsigned(1023u));
    testing::PrecalculatedRandom<kMaxLen * 3, int16_t> rnd16(-2048, 2048);
    for (uint32_t l = 0; l < kNumLoops; ++l) {
      rnd16.Fill(src, kMaxLen);
      rnd16.Fill(pred, kMaxLen);
      for (uint32_t idx : {3, 2, 1, 0}) {   // reverse order, for large to small
        const uint32_t len = (4 << idx);
        WP2::SubtractBlock[idx](src, 32, pred, 32, res, len, height);
        WP2::AddBlock[idx](pred, 32, res, len, min, max, dst, 32, height);
      }
    }
    uint32_t nb_errors = 0;
    for (uint32_t i = 0; i < height * 32; ++i) {
      if (src[i] <= min) {
        nb_errors += (dst[i] != min);
      } else if (src[i] >= max) {
        nb_errors += (dst[i] != max);
      } else {
        nb_errors += (src[i] != dst[i]);
      }
    }
    EXPECT_EQ(nb_errors, 0u) << " rate: " << (100 * nb_errors / kMaxLen) << "%";
  }
}

//------------------------------------------------------------------------------

template<typename T, int BitDepth>
void DoTestSSIM(int min_range, int max_range, double result) {
  constexpr uint32_t kMaxDim = 300;
  constexpr uint32_t kNumTests = 15;
  constexpr uint32_t kNumLoops = 7;
  PseudoRNG rng;
  double total_score = 0.;
  for (uint32_t t = 0; t < kNumTests; ++t) {
    const uint32_t width  = 1 + rng.GetUnsigned(kMaxDim);
    const uint32_t height = 1 + rng.GetUnsigned(kMaxDim);
    testing::PrecalculatedRandom<kMaxDim * 7, T> rnd(min_range, max_range);
    double score;
    Plane<T> plane1, plane2, plane3;
    ASSERT_WP2_OK(plane1.Resize(width, height));
    EXPECT_EQ(plane1.GetSSIM(plane2, BitDepth, &score),
              WP2_STATUS_BAD_DIMENSION);
    ASSERT_WP2_OK(plane2.Resize(width, height));
    ASSERT_WP2_OK(plane3.Resize(width, height));
    plane3.Fill(0);
    for (uint32_t y = 0; y < height; ++y) {
      T* const row1 = plane1.Row(y);
      T* const row2 = plane2.Row(y);
      rnd.Fill(row1, width);
      for (uint32_t x = 0; x < width; ++x) row2[x] = row1[x] ^ 63;
    }
    for (uint32_t l = 0; l < kNumLoops; ++l) {
      ASSERT_WP2_OK(plane1.GetSSIM(plane1, BitDepth, &score));
      EXPECT_EQ(score, 1. * width * height);
      ASSERT_WP2_OK(plane1.GetSSIM(plane2, BitDepth, &score));
      EXPECT_LT(score, 0.99999 * width * height);
      ASSERT_WP2_OK(plane1.GetSSIM(plane3, BitDepth, &score));
    }
    total_score += score;
  }
  EXPECT_NEAR(total_score, result, 1e-16);
}

TEST_P(SpeedTestCpuInfo, TestSSIM_16b) {
  DoTestSSIM<int16_t, 10>(-1024, 1024, 0.094818372359986899);
}

TEST_P(SpeedTestCpuInfo, TestSSIM_8b) {
  DoTestSSIM<uint8_t, 8>(0, 256, 5.2010943029822574);
}

//------------------------------------------------------------------------------

template<typename T>
void DoTestPSNR(int min_range, int max_range, uint64_t result) {
  constexpr uint32_t kMaxDim = 500;
  constexpr uint32_t kNumTests = 100;
  constexpr uint32_t kNumLoops = 50;
  PseudoRNG rng;
  uint64_t total_score = 0;
  for (uint32_t t = 0; t < kNumTests; ++t) {
    const uint32_t width  = 1 + rng.GetUnsigned(kMaxDim);
    const uint32_t height = 1 + rng.GetUnsigned(kMaxDim);
    testing::PrecalculatedRandom<kMaxDim * 7, T> rnd(min_range, max_range);
    uint64_t score;
    Plane<T> plane1, plane2, plane3;
    ASSERT_WP2_OK(plane1.Resize(width, height));
    EXPECT_EQ(plane1.GetSSE(plane2, &score), WP2_STATUS_BAD_DIMENSION);
    ASSERT_WP2_OK(plane2.Resize(width, height));
    ASSERT_WP2_OK(plane3.Resize(width, height));
    plane3.Fill(0);
    for (uint32_t y = 0; y < height; ++y) {
      T* const row1 = plane1.Row(y);
      T* const row2 = plane2.Row(y);
      rnd.Fill(row1, width);
      for (uint32_t x = 0; x < width; ++x) row2[x] = row1[x] ^ 1;
    }
    for (uint32_t l = 0; l < kNumLoops; ++l) {
      ASSERT_WP2_OK(plane1.GetSSE(plane1, &score));
      EXPECT_EQ(score, 0u);
      ASSERT_WP2_OK(plane1.GetSSE(plane2, &score));
      EXPECT_EQ(score, width * height);  // diff is '1' for each pixel
      ASSERT_WP2_OK(plane1.GetSSE(plane3, &score));
    }
    total_score += score;
  }
  EXPECT_EQ(total_score, result);
}

TEST_P(SpeedTestCpuInfo, TestPSNR_16b) {
  DoTestPSNR<int16_t>(-2048, 2048, 64311554390702llu);
}

TEST_P(SpeedTestCpuInfo, TestPSNR_8b) {
  DoTestPSNR<uint8_t>(0, 256, 147478357240llu);
}

//------------------------------------------------------------------------------

INSTANTIATE_TEST_SUITE_P(SpeedTestInstantiation, SpeedTestCpuInfo,
                         ::testing::ValuesIn(testing::kWP2CpuInfos));

}  // namespace
}  // namespace WP2
