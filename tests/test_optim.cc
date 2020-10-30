// Copyright 2018 Google LLC
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
//
// test for various Wiener-filtering optimization
//
// Author: Skal (pascal.massimino@gmail.com)

#include <cmath>
#include <cstdio>
#include <cstring>
#include <string>

#include "imageio/image_dec.h"
#include "imageio/image_enc.h"
#include "imageio/imageio_util.h"
#include "include/helpers.h"
#include "src/utils/plane.h"
#include "src/utils/wiener.h"
#include "src/wp2/base.h"

#include "examples/example_utils.h"  // must be last to avoid CPUInfo conflict

namespace WP2 {
namespace {

int16_t GetRand(int min, int max) {
  return (min < max) ? (int16_t)(min + (random() % (max - min))) : (int16_t)0;
}

void FindBestUpsampling(const Plane16& src, const Plane16& ref, float scale,
                        int16_t result[4]) {
  WienerOptimizer<3, 1> opt;
  for (uint32_t y = 0; y < src.h_ - 1; ++y) {
    for (uint32_t x = 0; x < src.w_ - 1; ++x) {
      const int16_t a = src.At(x + 0, y + 0);
      const int16_t b = src.At(x + 1, y + 0);
      const int16_t c = src.At(x + 0, y + 1);
      const int16_t d = src.At(x + 1, y + 1);
      const int16_t A = ref.At(2 * x + 0, 2 * y + 0);
      const int16_t B = ref.At(2 * x + 1, 2 * y + 0);
      const int16_t C = ref.At(2 * x + 0, 2 * y + 1);
      const int16_t D = ref.At(2 * x + 1, 2 * y + 1);
      // a . . b
      // . A B .
      // . C D .
      // c . . d
      // Note: we symmetrize the coefficients #1 and #2
      const int16_t s0[] = {a, (int16_t)(b + c), d}, t0[] = {A};
      const int16_t s1[] = {b, (int16_t)(d + a), c}, t1[] = {B};
      const int16_t s2[] = {c, (int16_t)(a + d), b}, t2[] = {C};
      const int16_t s3[] = {d, (int16_t)(c + b), a}, t3[] = {D};
      opt.AddSample(s0, t0);
      opt.AddSample(s1, t1);
      opt.AddSample(s2, t2);
      opt.AddSample(s3, t3);
    }
  }

  float tmp[1][4];
  EXPECT_TRUE(opt.Optimize(tmp)) << "Solver error";

  for (int i = 0; i < 3; ++i) {
    result[i] = (int16_t)std::lround(scale * tmp[0][i]);
  }
  EXPECT_LE(fabs((double)tmp[0][3]), 2.) << "Average is too large!";
  result[3] = 0;  // force mean to zero, as should be
}

void Upsample(const Plane16& uv, int16_t w[4], Plane16* const dst) {
  const uint32_t sum = (uint32_t)(w[0] + 2 * w[1] + w[2]);
  const double norm = (sum != 0) ? 1. / sum : 0.;
  const int16_t filter[2 * 2][4] = {{w[0], w[1], w[1], w[2]},
                                    {w[1], w[0], w[2], w[1]},
                                    {w[1], w[2], w[0], w[1]},
                                    {w[2], w[1], w[1], w[0]}};
  for (uint32_t y = 0; y < dst->h_; ++y) {
    for (uint32_t x = 0; x < dst->w_; ++x) {
      const uint32_t slot = (x & 1) + (y & 1) * 2;
      const uint32_t i = x >> 1;
      const uint32_t j = y >> 1;
      const int16_t a = uv.At(i + 0, j + 0);
      const int16_t b = uv.AtClamped(i + 1, j + 0);
      const int16_t c = uv.AtClamped(i + 0, j + 1);
      const int16_t d = uv.AtClamped(i + 1, j + 1);
      const int16_t* const W = &filter[slot][0];
      dst->At(x, y) = (int16_t)std::lround(
          (W[0] * a + W[1] * b + W[2] * c + W[3] * d) * norm + w[3]);
    }
  }
}

////////////////////////////////////////////////////////////////////////////////
// test for recovering pure gradients (w/ noise)

template <int kNum> void DoTestGradient() {
  for (int n = 0; n < 100; ++n) {
    const int16_t max_range = (int16_t)sqrt(32768. / (kNum + 1));
    const int16_t range = GetRand(2, max_range);
    EXPECT_LT(range * range * (kNum + 1), 32768)
         << "range overflows 16b signed calculations!";

    int16_t alpha[kNum + 1];
    for (size_t i = 0; i < kNum; ++i) alpha[i] = GetRand(-range, range);

    // average component:
    alpha[kNum] = GetRand(-range * range, range * range);

    const int16_t noise_range = GetRand(0, range / 4);  // add some noise
    const size_t num_samples = (size_t)(100 + GetRand(0, 500));

    WienerOptimizer<kNum, 1> opt;
    for (size_t m = 0; m < num_samples; ++m) {
      int16_t context[kNum];
      for (size_t i = 0; i < kNum; ++i) context[i] = GetRand(-range, range);
      int16_t value[1] = {0};
      // gradient component
      for (size_t i = 0; i < kNum; ++i) value[0] += context[i] * alpha[i];
      // DC component
      value[0] += alpha[kNum];
      // add noise
      value[0] += GetRand(-noise_range, noise_range);
      opt.AddSample(context, value);
    }

    float out[1][kNum + 1];
    EXPECT_TRUE(opt.Optimize(out)) << "Optimization shouldn't fail.";
    for (int i = 0; i <= kNum; ++i) {
      const double error = fabs(1. * alpha[i] - out[0][i]);
      EXPECT_LE(error, 4) << "optimization result outside tolerance";
    }
  }
}

TEST(OptimTest, GradientTest) {
  ASSERT_NO_FATAL_FAILURE(WP2::DoTestGradient<2>());
  ASSERT_NO_FATAL_FAILURE(WP2::DoTestGradient<3>());
  ASSERT_NO_FATAL_FAILURE(WP2::DoTestGradient<5>());
  ASSERT_NO_FATAL_FAILURE(WP2::DoTestGradient<10>());
}

////////////////////////////////////////////////////////////////////////////////

class ImageTest :
    public ::testing::TestWithParam<std::tuple<std::string, bool>> {};

TEST_P(ImageTest, ImageOptimTest) {
  bool debug = false;
  bool use_default = std::get<1>(GetParam());
  double scale = 16000.;  // ~14 bits

  EXPECT_TRUE(WP2CheckVersion());

  ArgbBuffer pic;
  size_t size = 0;
  const std::string name = testing::GetTestDataPath(std::get<0>(GetParam()));
  ASSERT_WP2_OK(ReadImage(name.c_str(), &pic, &size));

  YUVPlane yuv;
  CSPTransform transf;
  ASSERT_WP2_OK(yuv.Import(pic, pic.HasTransparency(), transf,
                           /*resize_if_needed=*/true));

  YUVPlane ref;
  ASSERT_WP2_OK(ref.Copy(yuv, /*resize_if_needed=*/true));
  ASSERT_WP2_OK(yuv.Downsample());
  int16_t Uw[4], Vw[4];
  if (!use_default) {
    FindBestUpsampling(yuv.U, ref.U, (float)scale, Uw);
    FindBestUpsampling(yuv.V, ref.V, (float)scale, Vw);
  } else {
    Uw[0] = Vw[0] = (int16_t)std::lround(scale * 4. / 9.);
    Uw[1] = Vw[1] = (int16_t)std::lround(scale * 2. / 9.);
    Uw[2] = Vw[2] = (int16_t)std::lround(scale * 1. / 9.);
    Uw[3] = Vw[3] = 0;
  }
  if (debug) {
    for (const auto& w : Uw) printf("[%d]", w);
    printf("\n");
    for (const auto& w : Vw) printf("[%d]", w);
    printf("\n");
  }
  EXPECT_GT(Uw[0], Uw[1]) << "strange coefficient results!";
  EXPECT_GT(Uw[1], Uw[2]) << "strange coefficient results!";
  EXPECT_GT(Vw[0], Vw[1]) << "strange coefficient results!";
  EXPECT_GT(Vw[1], Vw[2]) << "strange coefficient results!";

  // apply upsampling
  Upsample(yuv.U, Uw, &ref.U);
  Upsample(yuv.V, Vw, &ref.V);

  ArgbBuffer pic2;
  ASSERT_WP2_OK(ref.Export(transf, /*resize_if_needed=*/true, &pic2));
  if (debug) {
    const std::string temp_file = testing::GetTempDataPath("tmp2.png");
    ASSERT_TRUE(SaveImage(pic2, temp_file.c_str(), /*overwrite=*/true))
        << "Error while saving temp file " << temp_file;
    printf("saved %s\n", temp_file.c_str());
  }
  float disto[5];
  ASSERT_WP2_OK(pic2.GetDistortion(pic, PSNR, disto));
  if (debug) printf("PSNR: %f\n", disto[4]);
  EXPECT_GE(disto[4], 30.f);
}

INSTANTIATE_TEST_SUITE_P(ImageTestInstantiation, ImageTest,
  ::testing::Combine(::testing::Values("source0.ppm", "source1_64x48.png"),
                     ::testing::Values(true, false) /* use_default */));

////////////////////////////////////////////////////////////////////////////////

}  // namespace
}  // namespace WP2
