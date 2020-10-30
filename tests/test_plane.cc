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
// test for Plane
//
// Author: Skal (pascal.massimino@gmail.com)

#include <cstdio>
#include <cstdlib>
#include <string>

#include "imageio/image_dec.h"
#include "include/helpers.h"
#include "include/helpers_filter.h"
#include "src/utils/csp.h"
#include "src/utils/plane.h"
#include "src/wp2/base.h"

namespace WP2 {
namespace {

//------------------------------------------------------------------------------

template <class T>
void Fill(Plane<T>* const p) {
  for (uint32_t y = 0; y < p->h_; ++y) {
    for (uint32_t x = 0; x < p->w_; ++x) {
      p->At(x, y) = (x * 3 + y * 2 + 1) & 0xffffu;
    }
  }
}

//------------------------------------------------------------------------------

class PlaneTest : public ::testing::TestWithParam<std::string> {};

TEST_P(PlaneTest, Simple) {
  const std::string& file_name = GetParam();
  const tap_t up_taps[5] = {0, -3, 16, 5, -2};
  const SamplingTaps up_filter = {SamplingTaps::Taps(2, up_taps + 2),
                                  SamplingTaps::Taps(2, up_taps + 2)};
  const tap_t down_taps[2 * 2 + 1] = {1, 2, 10, 2, 1};
  const SamplingTaps down_filter{SamplingTaps::Taps(2, down_taps + 2),
                                 SamplingTaps::Taps(2, down_taps + 2)};

  ArgbBuffer pic;
  size_t size = 0;
  ASSERT_WP2_OK(
      ReadImage(testing::GetTestDataPath(file_name).c_str(), &pic, &size));

  ArgbBuffer pic2;
  YUVPlane yuv, yuv2;
  CSPTransform transf;
  const uint32_t bit_depth = transf.GetYUVPrecisionBits() + 1;
  ASSERT_TRUE(yuv.IsEmpty());
  ASSERT_WP2_OK(yuv.Import(pic, pic.HasTransparency(), transf,
                           /*resize_if_needed=*/true));
  ASSERT_TRUE(!yuv.IsEmpty());
  if (pic.HasTransparency()) {
    ASSERT_TRUE(yuv.HasAlpha());
  }
  ASSERT_WP2_OK(yuv2.Copy(yuv, /*resize_if_needed=*/true));
  ASSERT_EQ(yuv.IsEmpty(), yuv2.IsEmpty());
  ASSERT_EQ(yuv.IsView(), yuv2.IsView());
  // Expect equality.
  EXPECT_TRUE(
      testing::Compare(yuv, yuv2, transf.GetYUVPrecisionBits() + 1, file_name));

  ASSERT_WP2_OK(yuv.Downsample(down_filter));
  // Should fail because downsampled.
  EXPECT_EQ(yuv.Export(transf, /*resize_if_needed=*/true, &pic2),
            WP2_STATUS_INVALID_COLORSPACE);
  ASSERT_WP2_OK(yuv.Upsample(up_filter));

  float disto[5];
  ASSERT_WP2_OK(yuv2.GetDistortion(yuv, bit_depth, PSNR, disto));
  ASSERT_GT(disto[4], 35.f);
  ASSERT_LT(disto[4], 60.f);
  ASSERT_WP2_OK(yuv2.GetDistortion(yuv2, bit_depth, PSNR, disto));
  ASSERT_EQ(disto[4], 99.f);

  ASSERT_WP2_OK(yuv.Export(transf, /*resize_if_needed=*/true, &pic2));
  ASSERT_WP2_OK(pic2.GetDistortion(pic, PSNR, disto));
  ASSERT_GT(disto[4], 35.f);

  // some test view
  ASSERT_WP2_OK(yuv.SetView(yuv, {0, 0, pic.width, pic.height}));
  ASSERT_FALSE(yuv.Y.IsView());
  Fill(&yuv.Y);

  Plane<int16_t> p;
  ASSERT_WP2_OK(p.Resize(3, 6));
  Fill(&p);
  ASSERT_WP2_OK(p.SetView(p, {1, 1, 2, 5}));
  Fill(&p);
  ASSERT_WP2_OK(p.Resize(0, 0));
  Fill(&p);

  int16_t tmp[32 * 32] = {0};
  ASSERT_WP2_OK(p.SetView(tmp, /*width=*/7, /*height=*/7, /*step=*/32));
  Fill(&p);
  for (int i = 0; i < 32 * 32; ++i) {
    const int x = i % 32;
    const int y = i / 32;
    ASSERT_TRUE(tmp[x + y * 32] != 0 || x >= 7 || y >= 7);
  }
}

TEST_P(PlaneTest, TestSpeedSample) {
  const std::string& file_name = testing::GetTestDataPath(GetParam());
  const tap_t up_taps[5] = {0, -3, 16, 5, -2};
  const SamplingTaps up_filter = {SamplingTaps::Taps(2, up_taps + 2),
                                  SamplingTaps::Taps(2, up_taps + 2)};
  const tap_t down_taps[2 * 2 + 1] = {1, 2, 10, 2, 1};
  const SamplingTaps down_filter{SamplingTaps::Taps(2, down_taps + 2),
                                 SamplingTaps::Taps(2, down_taps + 2)};

  ArgbBuffer pic;
  ASSERT_WP2_OK(ReadImage(file_name.c_str(), &pic, nullptr));

  YUVPlane yuv, yuv2;
  CSPTransform transf;
  ASSERT_TRUE(yuv.IsEmpty());
  ASSERT_WP2_OK(yuv.Import(pic, pic.HasTransparency(), transf,
                           /*resize_if_needed=*/true));
  ASSERT_TRUE(!yuv.IsEmpty());

  ASSERT_WP2_OK(yuv2.Copy(yuv, /*resize_if_needed=*/true));
  constexpr uint32_t kNumTests = 8;
  for (uint32_t n = 0; n < kNumTests; ++n) {
    ASSERT_WP2_OK(yuv.Downsample(down_filter));
    ASSERT_WP2_OK(yuv.Upsample(up_filter));
  }
  float disto[5];
  const uint32_t bit_depth = transf.GetYUVPrecisionBits() + 1;
  ASSERT_WP2_OK(yuv2.GetDistortion(yuv, bit_depth, PSNR, disto));
  ASSERT_GT(disto[4], 32.f);
  ASSERT_LT(disto[4], 45.f);
}

INSTANTIATE_TEST_SUITE_P(PlaneTestInstantiation, PlaneTest,
                         ::testing::Values("source0.ppm", "source1_64x48.png",
                                           "source1_1x48.png",
                                           "source1_64x1.png"));

//------------------------------------------------------------------------------

TEST(PlaneTest, GetDistortion) {
  ArgbBuffer pic;
  ASSERT_WP2_OK(
      ReadImage(testing::GetTestDataPath("source1_64x48.png").c_str(), &pic));

  YUVPlane original, comparison;
  CSPTransform transf;
  const uint32_t bit_depth = transf.GetYUVPrecisionBits() + 1;
  ASSERT_WP2_OK(original.Import(pic, pic.HasTransparency(), transf,
                                /*resize_if_needed=*/true));
  ASSERT_WP2_OK(comparison.Copy(original, /*resize_if_needed=*/true));

  float disto[5];
  ASSERT_WP2_OK(comparison.GetDistortion(original, bit_depth, PSNR, disto));
  for (float distortion : disto) {
    ASSERT_EQ(distortion, 99.f);
  }
  for (Channel channel : {kYChannel, kUChannel, kVChannel, kAChannel}) {
    testing::Noise(transf.GetYUVMin(), transf.GetYUVMax(), /*seed=*/channel,
                   /*strength=*/5, &comparison.GetChannel(channel));
    ASSERT_WP2_OK(comparison.GetDistortion(original, bit_depth, PSNR, disto));
    const bool modified_channel[] = {channel >= kAChannel, channel >= kYChannel,
                                     channel >= kUChannel, channel >= kVChannel,
                                     /*all=*/true};
    for (uint32_t i = 0; i < 5; ++i) {
      if (modified_channel[i]) {
        EXPECT_GT(disto[i], 35.f);
        EXPECT_LT(disto[i], 60.f);
      } else {
        EXPECT_EQ(disto[i], 99.f);
      }
    }
  }
}

//------------------------------------------------------------------------------

}  // namespace
}  // namespace WP2
