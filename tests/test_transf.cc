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

// Test transforms.

#include <cmath>
#include <cstdio>
#include <cstring>
#include <random>
#include <tuple>

#include "examples/example_utils.h"
#include "include/helpers.h"
#include "src/common/lossy/transforms.h"
#include "src/dsp/dsp.h"
#include "src/utils/random.h"
#include "src/utils/stats.h"
#include "src/utils/utils.h"
#include "src/wp2/base.h"
#include "extras/extras.h"

namespace WP2 {
namespace {

constexpr WP2TransformType kTransformTypes[] = {kDct, kAdst, kHadamard,
                                                kIdentity};
STATIC_ASSERT_ARRAY_SIZE(kTransformTypes, kNumTransforms);

constexpr uint32_t kTransformSizes[] = {2, 4, 8, 16, 32};

//------------------------------------------------------------------------------

class TransposeTest : public ::testing::TestWithParam<WP2CPUInfo> {
  void SetUp() override {
    WP2DspReset();
    WP2GetCPUInfo = GetParam();
  }
};

TEST_P(TransposeTest, Transpose) {
  WP2TransformInit();

  for (uint32_t size_x : kTransformSizes) {
    for (uint32_t size_y : kTransformSizes) {
      SCOPED_TRACE(SPrintf("size_x: %d size_y: %d", size_x, size_y));
      const uint32_t num_coeffs = size_x * size_y;
      int32_t in[32 * 32], tmp[32 * 32], ref[32 * 32];
      for (uint32_t i = 0; i < num_coeffs; ++i) ref[i] = i;
      for (uint32_t i = num_coeffs; i < 32 * 32; ++i) ref[i] = -1;
      memcpy(in, ref, sizeof(ref));

      // Verify transposing.
      WP2Transpose(in, size_x, size_y, tmp);
      WP2Transpose(tmp, size_y, size_x, in);
      ASSERT_EQ(memcmp(in, ref, sizeof(ref[0]) * num_coeffs), 0);
      for (uint32_t i = num_coeffs; i < 32 * 32; ++i) ASSERT_EQ(in[i], -1);
    }
  }
}

INSTANTIATE_TEST_SUITE_P(TransposeTestInstantiation, TransposeTest,
    ::testing::ValuesIn(testing::kWP2CpuInfos));

//------------------------------------------------------------------------------

class RandomTransform2DTest : public ::testing::TestWithParam<size_t> {
 public:
  void BuildAutoRegression(int32_t dst[], uint32_t size_x, uint32_t size_y,
                           UniformIntDistribution* const random) {
    int32_t first = 0;
    const int32_t hx = size_x / 2, hy = size_y / 2;   // noise amplitude
    for (uint32_t y = 0; y < size_y; ++y) {
      int32_t last = first;
      for (uint32_t x = 0; x < size_x; ++x) {
        last += random->Get<int32_t>(-hx, hx);
        if (x & y & 1) {   // average with above
          last = (last + dst[x + (y - 1) * size_x] + 1) / 2;
        }
        dst[x + y * size_x] = last;
        if (x == 0) first = last + random->Get<int32_t>(-hy, hy);
      }
    }
  }
};

TEST_P(RandomTransform2DTest, All) {
  const int N = 50;
  const double tolerance = 0.60;

  UniformIntDistribution random(/*seed=*/GetParam());  // Deterministic.

  WP2TransformInit();

  // Used for coverage of WP2::Stats (not thread-safe).
  WP2::Stats<int> stats("simple-test", "%d", true);
  for (WP2TransformType tf_x : kTransformTypes) {
    for (WP2TransformType tf_y : kTransformTypes) {
      for (uint32_t size_x : kTransformSizes) {
        for (uint32_t size_y : kTransformSizes) {
          for (bool reduced : {false, true}) {
            if (reduced && (size_x == 2 || size_y == 2)) continue;
            for (bool inplace : {false, true}) {
              const uint32_t num_coeffs = size_x * size_y;
              uint32_t max_err = 0;
              SCOPED_TRACE(SPrintf(
                  "tf_x: %s size_x: %d tf_y: %s size_y: %d r:%d",
                  WP2TransformNames[tf_x], size_x,
                  WP2TransformNames[tf_y], size_y, reduced));
              for (int n = 0; n < N; ++n) {
                int32_t dst[32 * 32], ref[32 * 32] = {0};
                if (!reduced) {
                  for (auto& v : ref) {
                    v = random.Get<int32_t>(-1024, 1023);
                    stats.Add(v / 128);
                  }
                  // A pathological case.
                  if (size_x == 2 && size_y == 2 && n == 0) {
                    ref[0] = -966;
                    ref[1] = 978;
                    ref[2] = -394;
                    ref[3] = -370;
                  }
                } else {
                  // for reduced transform, the input needs to be quite smooth
                  // so we use a controlled auto-regression to fill the input
                  BuildAutoRegression(ref, size_x, size_y, &random);
                }

                // Verify transforms round-trip.
                if (inplace) {
                  std::memcpy(dst, ref, sizeof(ref));
                  WP2Transform2D(dst, tf_x, tf_y, size_x, size_y, dst, reduced);
                  WP2InvTransform2D(dst, tf_x, tf_y, size_x, size_y, dst,
                                    reduced);
                } else {
                  int32_t tmp[32 * 32];
                  WP2Transform2D(ref, tf_x, tf_y, size_x, size_y, tmp, reduced);
                  std::memset(dst, 0, sizeof(dst));
                  WP2InvTransform2D(tmp, tf_x, tf_y, size_x, size_y, dst,
                                    reduced);
                }

                // no spill
                for (size_t i = num_coeffs; i < 32 * 32; ++i) {
                  ASSERT_EQ(dst[i], inplace ? ref[i] : 0);
                }

                // Check error tolerance.
                uint32_t err = 0.;
                for (size_t i = 0; i < num_coeffs; ++i) {
                  err += std::abs(dst[i] - ref[i]);
                }
                double threshold =
                    tolerance * num_coeffs * (log2(size_x) + log2(size_y));
                // the dct-2 is less precise overall, so we need some headroom
                if (size_x == 2 || size_y == 2) {
                  threshold *= 2.1;
                }
                // reduced transform is less precise
                if (reduced) {
                  threshold *= 7.0;
                }
                // identity is a special case
                if (!reduced && tf_x == kIdentity && tf_y == kIdentity) {
                  threshold = 0;
                }
                if (err > threshold) {
                  printf("actual:   ");
                  for (size_t i = 0; i < num_coeffs; ++i) printf("%d ", dst[i]);
                  printf("\n");
                  printf("expected: ");
                  for (size_t i = 0; i < num_coeffs; ++i) printf("%d ", ref[i]);
                  printf("\n");
                }
                ASSERT_LE(err, threshold)
                    << "Precision error: " << size_x << ", " << size_y;
                max_err = std::max(err, max_err);
              }
            }
          }
        }
      }
    }
  }
}

// Quick and partial test for WP2ReduceCoeffs(), which only
// works with kDct and kAdst mostly.
TEST_P(RandomTransform2DTest, Compare) {
  const int kNumIterations = 25;
  const double tolerance = 2.6;

  UniformIntDistribution random(/*seed=*/GetParam());  // Deterministic.

  WP2TransformInit();
  for (WP2TransformType tf_x : {kDct, kAdst}) {
    for (WP2TransformType tf_y : {kDct, kAdst}) {
      for (uint32_t size_x : {4, 8, 16}) {
        for (uint32_t size_y : {4, 8, 16}) {
          const uint32_t tx = size_x / 2, ty = size_y / 2;
          const uint32_t thresh = tolerance * tx * ty * (log2(tx) + log2(ty));
          uint32_t max_err = 0;
          for (int n = 0; n < kNumIterations; ++n) {
            int32_t in1[32 * 32], in2[32 * 32], out1[32 * 32], out2[32 * 32];
            BuildAutoRegression(in1, size_x, size_y, &random);
            memcpy(in2, in1, sizeof(in1));
            // transform with both methods
            WP2Transform2D(in1, tf_x, tf_y, size_x, size_y, out1, true);
            WP2Transform2D(in2, tf_x, tf_y, size_x, size_y, out2, false);
            WP2ReduceCoeffs(out2, size_x, size_y, out2);
            // compare with inplace result
            WP2Transform2D(in1, tf_x, tf_y, size_x, size_y, in1, true);
            WP2Transform2D(in2, tf_x, tf_y, size_x, size_y, in2, false);
            WP2ReduceCoeffs(in2, size_x, size_y, in2);
            for (uint32_t i = 0; i < tx * ty / 4; ++i) {
              ASSERT_EQ(in1[i], out1[i]);
            }
            for (uint32_t i = 0; i < tx * ty; ++i) {
              ASSERT_EQ(in2[i], out2[i]);
            }
            // check error tolerance
            uint32_t err = 0.;
            for (uint32_t i = 0; i < tx * ty; ++i) {
              err += std::abs(in1[i] - in2[i]);
            }
            if (err > thresh) {
              printf("normal | crop:\n");
              for (uint32_t y = 0; y < ty; ++y) {
                for (uint32_t x = 0; x < tx; ++x) {
                  printf("%3d ", in1[x + tx * y]);
                }
                printf("    ");
                for (uint32_t x = 0; x < tx; ++x) {
                  printf("%3d ", in2[x + tx * y]);
                }
                printf("\n");
              }
              printf("-------------\n");
            }
            EXPECT_LE(err, thresh)
              << "Precision error: " << size_x << ", " << size_y;
            max_err = std::max(err, max_err);
          }
        }
      }
    }
  }
}

// Verifies that inplace and src/dst WP2Transform2D() versions give the same
// results.
TEST_P(RandomTransform2DTest, Inplace) {
  const int N = 50;

  UniformIntDistribution random(/*seed=*/GetParam() + 123);

  WP2TransformInit();

  for (WP2TransformType tf_x : kTransformTypes) {
    for (WP2TransformType tf_y : kTransformTypes) {
      for (uint32_t size_x : kTransformSizes) {
        for (uint32_t size_y : kTransformSizes) {
          for (bool reduced : {false, true}) {
            if (reduced && (size_x == 2 || size_y == 2)) continue;
            const uint32_t num_coeffs = size_x * size_y;
            const uint32_t num_tf_coeffs = num_coeffs / (reduced ? 4 : 1);
            for (int n = 0; n < N; ++n) {
              int32_t src[32 * 32] = {0};
              for (int32_t& v : src) v = random.Get<int32_t>(-1024, 1023);

              int32_t inplace_tf[32 * 32] = {0}, dst_tf[32 * 32] = {0};
              std::memcpy(inplace_tf, src, sizeof(src));
              WP2Transform2D(inplace_tf, tf_x, tf_y, size_x, size_y, inplace_tf,
                             reduced);
              WP2Transform2D(src, tf_x, tf_y, size_x, size_y, dst_tf, reduced);
              for (size_t i = 0; i < num_tf_coeffs; ++i) {
                ASSERT_EQ(inplace_tf[i], dst_tf[i]) << i;
              }

              int32_t inplace_invtf[32 * 32] = {0}, dst_invtf[32 * 32] = {0};
              std::memcpy(inplace_invtf, src, sizeof(src));
              WP2InvTransform2D(inplace_invtf, tf_x, tf_y, size_x, size_y,
                                inplace_invtf, reduced);
              WP2InvTransform2D(src, tf_x, tf_y, size_x, size_y, dst_invtf,
                                reduced);
              for (size_t i = 0; i < num_coeffs; ++i) {
                ASSERT_EQ(inplace_invtf[i], dst_invtf[i]) << i;
              }
            }
          }
        }
      }
    }
  }
}

// Using different seeds might trigger some failures but as long as the errors
// are close enough to the thresholds it's fine.
INSTANTIATE_TEST_SUITE_P(RandomTransform2DTestInstantiation,
                         RandomTransform2DTest,
                         /*seed=*/::testing::Values(1, 2, 3, 4));

TEST(Transform2DTest, Zero) {
  WP2TransformInit();

  for (bool reduced : {true, false}) {
    for (WP2TransformType tf_x : kTransformTypes) {
      for (WP2TransformType tf_y : kTransformTypes) {
        for (uint32_t size_x : kTransformSizes) {
          for (uint32_t size_y : kTransformSizes) {
            if (reduced && (size_x == 2 || size_y == 2)) continue;
            SCOPED_TRACE(SPrintf("tf_x: %s size_x: %d tf_y: %s size_y: %d r:%d",
                                 WP2TransformNames[tf_x], size_x,
                                 WP2TransformNames[tf_y], size_y, reduced));
            // Check that if all coeffs are zero then the transformed coeffs are
            // also zero, and we still get back zero after inverse transform.
            int32_t in[32 * 32] = {0};
            WP2Transform2D(in, tf_x, tf_y, size_x, size_y, in, reduced);
            for (auto& i : in) ASSERT_EQ(i, 0);

            WP2InvTransform2D(in, tf_x, tf_y, size_x, size_y, in, reduced);
            for (auto& i : in) ASSERT_EQ(i, 0);
          }
        }
      }
    }
  }
}

//------------------------------------------------------------------------------

class RandomTransformTest
    : public ::testing::TestWithParam<std::tuple<WP2TransformType, uint32_t>> {
};

TEST_P(RandomTransformTest, TestTransformPrecision) {
  const WP2TransformType type = std::get<0>(GetParam());
  const uint32_t size = std::get<1>(GetParam());
  const int N = 20000;
  const float tolerance = (type == kHadamard) ? 0.89f :
                          (type == kAdst) ? 0.30f :
                          (type == kIdentity) ? 0.00f :
                          0.15f;

  UniformIntDistribution random(/*seed=*/(size_t)size +
                                1000 * (size_t)type);  // Deterministic.

  WP2TransformInit();

  int32_t in[32] = {0}, out[32] = {0}, check[32] = {0};
  const uint32_t S = size;
  double err[32] = {0.};
  for (int n = 0; n < N; ++n) {
    for (size_t i = 0; i < S; ++i) {
      in[i] = random.Get<int32_t>(-1024, 1023);
      out[i] = random.Get<int32_t>(0, 2147483647);
      check[i] = random.Get<int32_t>(0, 2147483647);
    }
    WP2Transform(in, type, size, out);
    WP2InvTransform(out, type, size, check);

    for (size_t i = 0; i < S; ++i) err[i] += abs(in[i] - check[i]);
    for (size_t i = S; i < 32; ++i) {  // Check there's no spill.
      ASSERT_EQ(in[i], 0);
      ASSERT_EQ(out[i], 0);
      ASSERT_EQ(check[i], 0);
    }
  }
  const float tolerance_max = tolerance * log2(S);
  double max_err = 0.;
  for (size_t i = 0; i < S; ++i) {
    err[i] /= N;
    printf("%.2lf ", err[i]);
    max_err = std::max(max_err, err[i]);
    ASSERT_LE(err[i], tolerance_max);
  }
  printf("  type=%d size=%d max-err=%.3lf tolerance=%.3f\n",
         type, S, max_err, tolerance_max);
}

INSTANTIATE_TEST_SUITE_P(
    RandomTransformTestInstantiation, RandomTransformTest,
    ::testing::Combine(::testing::ValuesIn(kTransformTypes),
                       ::testing::ValuesIn(kTransformSizes)));

//------------------------------------------------------------------------------

TEST(TransformClass, All) {
  EXPECT_EQ(GetTransformClass(kDctDct), TransformClass::kTwoD);
  EXPECT_EQ(GetTransformClass(kDctAdst), TransformClass::kTwoD);
  EXPECT_EQ(GetTransformClass(kIdentityIdentity), TransformClass::kTwoD);
  EXPECT_EQ(GetTransformClass(kIdentityDct), TransformClass::kVertical);
  EXPECT_EQ(GetTransformClass(kDctIdentity), TransformClass::kHorizontal);
}

//------------------------------------------------------------------------------
// This is not really a test; this is used to generate the representation of
// each transform coefficient into an image.

class PrintTransformTest
    : public ::testing::TestWithParam<
          std::tuple<WP2TransformType, WP2TransformType, uint32_t>> {};

TEST_P(PrintTransformTest, Simple) {
  WP2TransformInit();
  const WP2TransformType tf_x = std::get<0>(GetParam());
  const WP2TransformType tf_y = std::get<1>(GetParam());
  const uint32_t size = std::get<2>(GetParam());
  const uint32_t scale = 16;

  ArgbBuffer image;
  ASSERT_WP2_OK(image.Resize(size * size, size * size));
  const uint32_t bpp = WP2FormatBpp(image.format);
  for (uint32_t y = 0; y < size; ++y) {
    for (uint32_t x = 0; x < size; ++x) {
      int32_t coeffs[32 * 32] = {0};
      coeffs[y * size + x] = 125;
      WP2InvTransform2D(coeffs, tf_x, tf_y, size, size, coeffs);

      uint8_t* row = (uint8_t*)image.GetRow(y * size) + (x * size) * bpp;
      for (uint32_t sub_y = 0; sub_y < size; ++sub_y) {
        for (uint32_t sub_x = 0; sub_x < size; ++sub_x) {
          const int32_t color = 128 + coeffs[sub_y * size + sub_x];
          assert(color >= 0 && color <= 255);
          row[sub_x * bpp] = 255;        // Alpha
          row[sub_x * bpp + 1] = color;  // Grey
          row[sub_x * bpp + 2] = color;
          row[sub_x * bpp + 3] = color;
        }
        row += image.stride;
      }
    }
  }

  ASSERT_WP2_OK(
      ArgbBufferRescale(&image, image.width * scale, image.height * scale));
  const Argb32b red = {0xFF, 0xFF, 0x00, 0x00};
  for (uint32_t i = 1; i < size; ++i) {
    const uint32_t row_or_col = i * size * scale;
    image.Fill({0, row_or_col, image.width, 1}, red);
    image.Fill({row_or_col, 0, 1, image.height}, red);
  }

  char file_path[100];  // Change output folder here.
  std::snprintf(file_path, sizeof(file_path) / sizeof(file_path[0]),
                "/tmp/transform_%s_%s_%ux%u.png", WP2TransformNames[tf_x],
                WP2TransformNames[tf_y], size, size);
  testing::DumpImage(image, file_path);
}

// "DISABLED_" keeps it from running.
INSTANTIATE_TEST_SUITE_P(
    DISABLED_PrintTransformTestInstantiation, PrintTransformTest,
    ::testing::Combine(::testing::ValuesIn(kTransformTypes),
                       ::testing::ValuesIn(kTransformTypes),
                       ::testing::Values(8)));

//------------------------------------------------------------------------------

// Test which transforms give a constant result for DC only coeffs.
TEST(TransformProperty, DCOnlyImpliesConstant) {
  WP2TransformInit();

  int32_t coeffs[kMaxBlockSizePix2];
  for (int32_t i = -10; i < 10; ++i) {
    for (uint32_t size_x : kTransformSizes) {
      for (uint32_t size_y : kTransformSizes) {
        for (uint32_t tf_x = 0; tf_x < (uint32_t)kNumTransforms; ++tf_x) {
          for (uint32_t tf_y = 0; tf_y < (uint32_t)kNumTransforms; ++tf_y) {
            std::fill(coeffs, coeffs + size_x * size_y, 0);
            coeffs[0] = i;
            WP2InvTransform2D(coeffs, (WP2TransformType)tf_x,
                              (WP2TransformType)tf_y, size_x, size_y, coeffs,
                              /*reduced=*/false);
            const bool is_constant =
                std::all_of(coeffs, coeffs + size_x * size_y,
                            [&coeffs](int32_t v) { return (v == coeffs[0]); });
            if (i == 0) {
              EXPECT_TRUE(coeffs[0] == 0 && is_constant);
              continue;
            }
            EXPECT_TRUE(tf_x != kDct || tf_y != kDct || is_constant)
                << "DC = " << i << " w=" << size_x << " h=" << size_y;
          }
        }
      }
    }
  }
}

//------------------------------------------------------------------------------

}  // namespace
}  // namespace WP2
