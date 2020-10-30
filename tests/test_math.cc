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

// Tests for math functions and random generators

#include <iostream>
#include <random>
#include <limits>
#include <vector>

#include "include/helpers.h"
#include "src/dsp/math.h"
#include "src/utils/random.h"

namespace WP2 {
namespace {

TEST(Math, InnerProduct) {
  WP2MathInit();
  int16_t a[100];
  int16_t b[100];
  UniformIntDistribution random(/*seed=*/125);
  for (size_t i = 0; i < 1000; ++i) {
    // Choose values small enough to fit residuals.
    for (int16_t& j : a) j = random.Get<int16_t>(0, 999);
    for (int16_t& j : b) j = random.Get<int16_t>(0, 999);
    // Try for different sizes of the vectors.
    for (size_t j = 1; j < 100; ++j) {
      const int sum = WP2InnerProduct(a, b, j);
      EXPECT_EQ(sum, std::inner_product(a, a + j, b, 0));
    }
  }
}

TEST(Math, PseudoRandom) {
  PseudoRNG rng;
  for (int32_t amp : { 1, 4, 16, 65536, 0x7fffffff }) {
    for (uint32_t n = 0; n < 10000; ++n) {
      const int32_t v = rng.GetSigned(amp);
      EXPECT_LT(v, amp);
      EXPECT_GT(v, -amp);
    }
  }
}

TEST(Math, PseudoRandomStats) {
  PseudoRNG rng;
  const int32_t kAmp = 100;
  const int32_t kCount = 1000000;
  int32_t histo[2 * kAmp - 1] = { 0 };
  const int32_t num_samples = (2 * kAmp - 1) * kCount;
  for (uint32_t n = 0; n < num_samples; ++n) {
    ++histo[kAmp - 1 + rng.GetSigned(kAmp)];
  }
  const int32_t tolerance = sqrt(num_samples / 25);
  float max_err = 0.;
  for (const auto h : histo) {
    EXPECT_LT(std::abs(h - kCount), tolerance);
    const float err = 100.f * (h - kCount) / kCount;
    max_err = std::max(max_err, err);
  }
  printf("Max error: %.3f%%\n", max_err);
}

TEST(Math, PseudoRandomStats2) {
  PseudoRNG rng;
  for (int32_t amp : {77, 101, 259, 1001}) {
    constexpr int32_t kCount = 200000;
    std::vector<int32_t> histo(amp, 0);
    const int32_t num_samples = amp * kCount;
    for (int32_t n = 0; n < num_samples; ++n) {
      ++histo[rng.GetUnsigned(amp)];
    }
    const int32_t tolerance = sqrt(num_samples / 11);
    float max_err = 0.;
    for (int32_t i = 0; i < amp; ++i) {
      EXPECT_LT(std::abs(histo[i] - kCount), tolerance)
          << "# " << i << ";  Amp: " << amp;
      const float err = 100.f * (histo[i] - kCount) / kCount;
      max_err = std::max(max_err, err);
    }
    printf("Amp:%u, Max error: %.3f%%\n", amp, max_err);
  }
}

TEST(Math, Random) {
  UniformIntDistribution dist;
  for (int32_t amp : { 1, 4, 16, 65536, 0x7fffffff }) {
    for (uint32_t n = 0; n < 10000; ++n) {
      const int32_t v = dist.Get(-(amp - 1), amp - 1);
      EXPECT_LT(v, amp);
      EXPECT_GT(v, -amp);
    }
  }
}

TEST(Math, RandomStats) {
  UniformIntDistribution dist;
  const int32_t kAmp = 100;
  const int32_t kCount = 100000;
  int32_t histo[2 * kAmp - 1] = { 0 };
  const int32_t num_samples = (2 * kAmp - 1) * kCount;
  for (uint32_t n = 0; n < num_samples; ++n) {
    ++histo[kAmp - 1 + dist.Get(-kAmp + 1, kAmp - 1)];
  }
  const int32_t tolerance = sqrt(num_samples / 18);
  float max_err = 0.;
  for (const auto h : histo) {
    EXPECT_LT(std::abs(h - kCount), tolerance);
    const float err = 100.f * (h - kCount) / kCount;
    max_err = std::max(max_err, err);
  }
  printf("Max error: %.3f%%\n", max_err);
  EXPECT_LT(max_err, 0.9f);
}

TEST(Math, Log2) {
  WP2MathInit();
  const uint32_t kMax = 1u << 17;  // kApproxLogWithCorrectionMax is 1 << 16
  const uint32_t kIter = 100;
  double err = 0., sum = 0.;
  for (uint32_t iter = 0; iter < kIter; ++iter) {
    for (uint32_t v = 0; v < kMax; ++v) {
      const float log2_v = WP2Log2(v);
      err += std::abs(log2_v - (v > 0 ? std::log2(v) : 0));
      sum += log2_v;
    }
  }
  err /= kIter * kMax * sum;
  printf("log2 error: %.3e\n", err);
  EXPECT_LT(err, 9.5e-14);
}

TEST(Math, SLog2) {
  WP2MathInit();
  const uint32_t kMax = 1u << 17;  // kApproxLogWithCorrectionMax is 1 << 16
  const uint32_t kIter = 100;
  double err = 0., sum = 0.;
  for (uint32_t iter = 0; iter < kIter; ++iter) {
    for (uint32_t v = 0; v < kMax; ++v) {
      const float slog2_v = WP2SLog2(v);
      err += std::abs(slog2_v - (v > 0 ? v * std::log2(v) : 0));
      sum += slog2_v;
    }
  }
  err /= kIter * kMax * sum;
  printf("slog2 error: %.3e\n", err);
  EXPECT_LT(err, 1.0e-14);
}

TEST(Math, Log2Floor) {
  ASSERT_EQ(WP2Log2Floor_k(0), 0u);
  ASSERT_EQ(WP2Log2Floor_k(1), 0u);
  for (uint32_t v = 2; v < 4096; ++v) {
    const uint32_t log2 = WP2Log2Floor_k(v);
    ASSERT_GE(v, 1u << log2);
    ASSERT_LT(v, 1u << (log2 + 1u));
  }
  ASSERT_EQ(WP2Log2Floor_k(255), 7u);
  ASSERT_EQ(WP2Log2Floor_k(256), 8u);
  ASSERT_EQ(WP2Log2Floor_k(257), 8u);
  ASSERT_EQ(WP2Log2Floor_k(0x7FFFFFFFFFFFFFFEull), 62u);
  ASSERT_EQ(WP2Log2Floor_k(0x7FFFFFFFFFFFFFFFull), 62u);
  ASSERT_EQ(WP2Log2Floor_k(0x8000000000000000ull), 63u);
  ASSERT_EQ(WP2Log2Floor_k(0xFFFFFFFFFFFFFFFEull), 63u);
  ASSERT_EQ(WP2Log2Floor_k(0xFFFFFFFFFFFFFFFFull), 63u);
}

TEST(Math, CtzClz) {
  UniformIntDistribution random(/*seed=*/421523);
  for (int i = 0; i < 32; ++i) {
    for (uint32_t n = 0; n < 1000; ++n) {
      const uint32_t v0 = random.Get<uint16_t>(0u, 0xffff)
         | (random.Get<uint16_t>(0u, 0xffff) << 16);
      // clear 'i' lower bits, and force ith bit to 1
      const uint32_t v1 = ((v0 >> i) | 1) << i;
      ASSERT_EQ(WP2Ctz(v1), i) << " bit : " << i << " value: " << v1;
      // same with high bits
      const uint32_t v2 = ((v0 << i) | (1u << 31)) >> i;
      ASSERT_EQ(WP2Log2Floor(v2), 31 - i) << " bit : " << i << " value: " << v2;
    }
    ASSERT_EQ(WP2Ctz(1 << i), i);
    ASSERT_EQ(WP2Log2Floor(1u << i), i);
  }
}

TEST(Math, Log2Ceil) {
  ASSERT_EQ(WP2Log2Ceil_k(0u), 0u);
  ASSERT_EQ(WP2Log2Ceil_k(1u), 0u);
  for (uint32_t v = 2; v < 4096; ++v) {
    const uint32_t log2 = WP2Log2Ceil_k(v);
    ASSERT_LE(v, 1u << log2);
    ASSERT_GT(v, 1u << (log2 - 1u));
  }
  ASSERT_EQ(WP2Log2Ceil_k(255), 8u);
  ASSERT_EQ(WP2Log2Ceil_k(256), 8u);
  ASSERT_EQ(WP2Log2Ceil_k(257), 9u);
  ASSERT_EQ(WP2Log2Ceil_k(0x7FFFFFFFFFFFFFFEull), 63u);
  ASSERT_EQ(WP2Log2Ceil_k(0x7FFFFFFFFFFFFFFFull), 63u);
  ASSERT_EQ(WP2Log2Ceil_k(0x8000000000000000ull), 63u);
  ASSERT_EQ(WP2Log2Ceil_k(0x8000000000000001ull), 64u);
  ASSERT_EQ(WP2Log2Ceil_k(0xFFFFFFFFFFFFFFFEull), 64u);
  ASSERT_EQ(WP2Log2Ceil_k(0xFFFFFFFFFFFFFFFFull), 64u);
}

TEST(Math, SqrtFloor) {
  for (uint32_t v = 0; v < 10000; ++v) {
    const uint32_t r = SqrtFloor(v);
    ASSERT_LE(r * r, v);
    ASSERT_GT((r + 1) * (r + 1), v);
  }
}

}  // namespace
}  // namespace WP2
