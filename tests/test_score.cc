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

// tests for scoring functions

#include "examples/example_utils.h"
#include "include/helpers.h"
#include "src/dsp/dsp.h"
#include "src/dsp/math.h"
#include "src/utils/random.h"

namespace WP2 {
namespace {

class ScoreTestCpuInfo : public ::testing::TestWithParam<WP2CPUInfo> {
  void SetUp() override {
    WP2DspReset();
    WP2GetCPUInfo = GetParam();
    WP2MathInit();
    WP2::ScoreDspInit();
  }
};

TEST_P(ScoreTestCpuInfo, TestMinMax) {
  UniformIntDistribution gen(76234);
  constexpr uint32_t step = kMinBlockSizePix + 7;
  constexpr uint32_t kNumBlks = 25;
  int16_t data0[kMinBlockSizePix * step * kNumBlks];
  for (uint32_t y = 0; y < kMinBlockSizePix; ++y) {
    for (uint32_t x = 0; x < kMinBlockSizePix * kNumBlks; ++x) {
      data0[x + y * step] = gen.Get(-1024, 1024);
    }
  }

  constexpr uint32_t kNumTests = 100;
  uint32_t CRC = 4124;
  for (uint32_t range = 1; range < 1024; ++range) {
    int16_t data[kMinBlockSizePix * step * kNumBlks];
    for (uint32_t y = 0; y < kMinBlockSizePix; ++y) {
      for (uint32_t x = 0; x < kMinBlockSizePix * kNumBlks; ++x) {
        data[x + y * step] = data0[x + y * step] % range;
      }
    }
    int16_t min[kNumBlks], max[kNumBlks];
    for (uint32_t ind = 0; ind < kNumTests; ++ind) {
      WP2::GetBlockMinMax(data, step, kNumBlks, min, max);
    }
    for (auto& v : min) CRC = testing::SillyCRC(CRC, v);
    for (auto& v : max) CRC = testing::SillyCRC(CRC, v);
  }
  EXPECT_EQ(CRC, 3632092995u);
}

TEST_P(ScoreTestCpuInfo, TestMinMax_5x5) {
  UniformIntDistribution gen(76234);
  const uint32_t kW = 47, kH = 53;
  constexpr uint32_t kStep = kW + 7;
  int16_t data0[kW * kH];
  for (auto& v : data0) v = gen.Get(-1024, 1024);

  constexpr uint32_t kNumTests = 1000;
  uint32_t CRC = 4315;
  int16_t data[kH * kStep];
  uint64_t num_errors = 0;
  for (uint32_t range = 1; range < 1024; ++range) {
    for (uint32_t y = 0; y < kH; ++y) {
      for (uint32_t x = 0; x < kW; ++x) {
        data[x + y * kStep] = data0[x + y * kW] % range;
      }
    }
    int16_t min, max;
    for (uint32_t ind = 0; ind < kNumTests; ++ind) {
      const uint32_t off = (ind % (kW - 8)) + kStep * ((ind / kW) % (kH - 5));
      WP2::GetBlockMinMax_5x5(data + off, kStep, &min, &max);
      CRC = testing::SillyCRC(CRC, min);
      CRC = testing::SillyCRC(CRC, max);
    }
    int16_t min2, max2;
    WP2::GetBlockMinMax_5x5(data, kStep, &min, &max);
    GetBlockMinMaxGeneric(data, kStep, 5, 5, &min2, &max2);
    num_errors += (min != min2);
    num_errors += (max != max2);
  }
  EXPECT_EQ(CRC, 2416416445u);
  EXPECT_EQ(num_errors, 0u);
}

INSTANTIATE_TEST_SUITE_P(ScoreTestCpuInfoInstantiation, ScoreTestCpuInfo,
                         ::testing::ValuesIn(testing::kWP2CpuInfos));

}  // namespace
}  // namespace WP2
