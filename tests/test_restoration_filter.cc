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

// Test the directional filter.

#include <array>
#include <tuple>

#include "include/helpers.h"
#include "include/helpers_filter.h"
#include "src/dec/filters/restoration_filter.h"
#include "src/dsp/math.h"
#include "src/utils/utils.h"
#include "src/wp2/format_constants.h"

using ::testing::ElementsAre;

namespace WP2 {
namespace {

//------------------------------------------------------------------------------

YUVPlane CreateNoisyPixels(uint32_t width, uint32_t height,
                           int16_t min_value, int16_t max_value) {
  YUVPlane pixels;
  WP2_ASSERT_STATUS(pixels.Resize(width, height));
  pixels.Fill(Ayuv38b{kAlphaMax, min_value, 0, max_value});
  for (Channel channel : {kYChannel, kUChannel, kVChannel, kAChannel}) {
    testing::Noise(min_value, max_value, /*seed=*/channel, 5,
                   &pixels.GetChannel(channel));
  }
  return pixels;
}

// Applies the restoration filter on 'pixels' with 'filter_strength'.
void Filter(uint32_t num_precision_bits, int16_t min_value, int16_t max_value,
            const int8_t tap_weights[kWieFltTapDist], uint32_t quality_hint,
            YUVPlane* const pixels) {
  const Vector<Segment> segments =
      testing::CreateSegments(quality_hint, min_value, max_value);
  const double residual_density = 1. * quality_hint / kMaxLossyQualityHint;
  const double min_bpp = 4. * quality_hint / kMaxLossyQualityHint;
  FilterBlockMap blocks;
  testing::CreateBlocks(/*tile_x=*/0, /*tile_y=*/0, segments, residual_density,
                        min_bpp, num_precision_bits, pixels, &blocks);

  RstrFltParams params(pixels->Y.w_, pixels->Y.h_);
  params.SetAll(tap_weights);
  DecoderConfig config = DecoderConfig::kDefault;
  config.enable_restoration_filter = true;
  RestorationFilter restoration_filter(config, blocks, params);
  ASSERT_WP2_OK(restoration_filter.Allocate());

  uint32_t num_deblocked_rows = 0;
  for (uint32_t num_decoded_yuv_rows = 0; num_decoded_yuv_rows <= pixels->Y.h_;
       ++num_decoded_yuv_rows) {
    const uint32_t n = restoration_filter.Enhances(num_decoded_yuv_rows);
    ASSERT_GE(n, num_deblocked_rows);
    num_deblocked_rows = n;
  }
  EXPECT_EQ(num_deblocked_rows, pixels->Y.h_);
}

//------------------------------------------------------------------------------

constexpr int8_t kIdentityTapWeights[kWieFltTapDist]{0, 0, 0};
constexpr int8_t kDefaultTapWeights[kWieFltTapDist]{3, -7, 15};

class RestorationFilterTest : public ::testing::TestWithParam<
                                  std::tuple<uint32_t, uint32_t, uint32_t>> {};

// Tests that no pixel is modified if each plane has one unique color.
TEST_P(RestorationFilterTest, FilteringPlain) {
  const uint32_t num_precision_bits = std::get<0>(GetParam());
  const int16_t min_value = testing::GetMinValueSigned(num_precision_bits);
  const int16_t max_value = testing::GetMaxValueSigned(num_precision_bits);
  const uint32_t width = std::get<1>(GetParam());
  const uint32_t height = std::get<2>(GetParam());

  YUVPlane pixels =
      testing::CreatePlainPixels(width, height, min_value, max_value);
  YUVPlane reference;
  ASSERT_WP2_OK(reference.Copy(pixels, /*resize_if_needed=*/true));

  ASSERT_NO_FATAL_FAILURE(Filter(num_precision_bits, min_value, max_value,
                                 kDefaultTapWeights,
                                 /*quality_hint=*/0, &pixels));

  EXPECT_TRUE(testing::AreEqual(pixels.Y, reference.Y));
  EXPECT_TRUE(testing::AreEqual(pixels.U, reference.U));
  EXPECT_TRUE(testing::AreEqual(pixels.V, reference.V));
}

// Tests that pixels are filtered.
TEST_P(RestorationFilterTest, MaxFilteringBlocked) {
  const uint32_t num_precision_bits = std::get<0>(GetParam());
  const int16_t min_value = testing::GetMinValueSigned(num_precision_bits);
  const int16_t max_value = testing::GetMaxValueSigned(num_precision_bits);
  const uint32_t width = std::get<1>(GetParam());
  const uint32_t height = std::get<2>(GetParam());

  YUVPlane pixels = CreateNoisyPixels(width, height, min_value, max_value);
  YUVPlane reference;
  ASSERT_WP2_OK(reference.Copy(pixels, /*resize_if_needed=*/true));

  ASSERT_NO_FATAL_FAILURE(Filter(num_precision_bits, min_value, max_value,
                                 kDefaultTapWeights,
                                 /*quality_hint=*/0, &pixels));

  EXPECT_FALSE(testing::AreEqual(pixels.Y, reference.Y));
  EXPECT_FALSE(testing::AreEqual(pixels.U, reference.U));
  EXPECT_FALSE(testing::AreEqual(pixels.V, reference.V));
}

// Tests that no pixel is modified if it's the highest quality.
TEST_P(RestorationFilterTest, NoFiltering) {
  const uint32_t num_precision_bits = std::get<0>(GetParam());
  const int16_t min_value = testing::GetMinValueSigned(num_precision_bits);
  const int16_t max_value = testing::GetMaxValueSigned(num_precision_bits);
  const uint32_t width = std::get<1>(GetParam());
  const uint32_t height = std::get<2>(GetParam());

  YUVPlane pixels = CreateNoisyPixels(width, height, min_value, max_value);
  YUVPlane reference;
  ASSERT_WP2_OK(reference.Copy(pixels, /*resize_if_needed=*/true));

  ASSERT_NO_FATAL_FAILURE(
      Filter(num_precision_bits, min_value, max_value, kDefaultTapWeights,
             /*quality_hint=*/kMaxLossyQualityHint, &pixels));

  EXPECT_TRUE(testing::AreEqual(pixels.Y, reference.Y));
  EXPECT_TRUE(testing::AreEqual(pixels.U, reference.U));
  EXPECT_TRUE(testing::AreEqual(pixels.V, reference.V));
}

// Tests that no pixel is modified if the filter tap weights are the identity.
TEST_P(RestorationFilterTest, Identity) {
  const uint32_t num_precision_bits = std::get<0>(GetParam());
  const int16_t min_value = testing::GetMinValueSigned(num_precision_bits);
  const int16_t max_value = testing::GetMaxValueSigned(num_precision_bits);
  const uint32_t width = std::get<1>(GetParam());
  const uint32_t height = std::get<2>(GetParam());

  YUVPlane pixels = CreateNoisyPixels(width, height, min_value, max_value);
  YUVPlane reference;
  ASSERT_WP2_OK(reference.Copy(pixels, /*resize_if_needed=*/true));

  ASSERT_NO_FATAL_FAILURE(Filter(num_precision_bits, min_value, max_value,
                                 kIdentityTapWeights,
                                 /*quality_hint=*/0, &pixels));

  EXPECT_TRUE(testing::AreEqual(pixels.Y, reference.Y));
  EXPECT_TRUE(testing::AreEqual(pixels.U, reference.U));
  EXPECT_TRUE(testing::AreEqual(pixels.V, reference.V));
}

INSTANTIATE_TEST_SUITE_P(
    RestorationFilterTestInstantiation, RestorationFilterTest,
    ::testing::Combine(
        /*num_precision_bits=*/::testing::Values(8u, 10u),
        /*width=*/::testing::Values(kWieFltWidth / 3),
        /*height=*/::testing::Values(3 * kWieFltHeight - 1u)));

// Disabled to reduce the number of test cases.
// Can still be run with flag --test_arg=--gunit_also_run_disabled_tests
INSTANTIATE_TEST_SUITE_P(
    DISABLED_MoreRestorationFilterTestInstantiation, RestorationFilterTest,
    ::testing::Combine(
        /*num_precision_bits=*/::testing::Values(8u, 10u, 12u),
        /*width=*/::testing::Values(kWieFltWidth / 2, kWieFltWidth),
        /*height=*/
        ::testing::Values(kWieFltHeight / 5, kWieFltHeight,
                          3 * kWieFltHeight + 7u)));

//------------------------------------------------------------------------------

// Tests that PSNR is good with only restoration filtering.
TEST(RestorationFilterTest, Simple) {
  EncoderConfig encoder_config = EncoderConfig::kDefault;
  encoder_config.quality = 40;  // Expect compression artifacts.

  DecoderConfig decoder_config = DecoderConfig::kDefault;
  decoder_config.enable_deblocking_filter = false;
  decoder_config.enable_directional_filter = false;
  decoder_config.enable_restoration_filter = true;

  ASSERT_WP2_OK(testing::EncodeDecodeCompare("source4.webp", encoder_config,
                                             decoder_config));
}

// Tests that PSNR is better with the restoration filter than without.
TEST(RestorationFilterTest, BetterPSNR) {
  EncoderConfig encoder_config = EncoderConfig::kDefault;
  encoder_config.quality = 50;  // Expect compression artifacts.

  DecoderConfig decoder_config_no_filter = DecoderConfig::kDefault;
  decoder_config_no_filter.enable_restoration_filter = false;

  DecoderConfig decoder_config_filter = decoder_config_no_filter;
  decoder_config_filter.enable_restoration_filter = true;

  float disto_diff[5];
  ASSERT_WP2_OK(testing::GetDistortionDiff("source1.png", {1, 0, 300, 200},
                                           encoder_config,
                                           decoder_config_no_filter,
                                           decoder_config_filter, disto_diff));
  EXPECT_GT(disto_diff[4], 0.f);
}

//------------------------------------------------------------------------------

TEST(WienerTest, WienerHalfToFullWgts) {
  std::vector<int32_t> h{1, 2, 3, 0, 0, 0, 0};
  WienerHalfToFullWgts(h.data(), h.data());

  EXPECT_THAT(h, ElementsAre(h[0], h[1], h[2], 128 - 2 * (h[0] + h[1] + h[2]),
                             h[2], h[1], h[0]));
}

//------------------------------------------------------------------------------
// Test some 3x3 patterns, filtered by kernels of 3 horizontal and 3 vertical
// taps.

WP2Status FilterCanvas(const std::array<int32_t, 3>& half_tap_weights_h,
                       const std::array<int32_t, 3>& half_tap_weights_v,
                       uint32_t width, uint32_t height,
                       const std::vector<int16_t>& canvas,
                       std::vector<int16_t>* const filtered) {
  WienerFilterInit();
  int32_t tap_weights_h[kWieFltNumTaps];
  int32_t tap_weights_v[kWieFltNumTaps];
  WienerHalfToFullWgts(half_tap_weights_h.data(), tap_weights_h);
  WienerHalfToFullWgts(half_tap_weights_v.data(), tap_weights_v);
  std::vector<uint8_t> strength_map(canvas.size(), 63);
  *filtered = canvas;
  return WienerFilter(width, height, nullptr, 0, nullptr, 0, nullptr, 0,
                      nullptr, 0, 0, 0, 0, 0, tap_weights_h, tap_weights_v,
                      strength_map.data(), width, /*num_precision_bits=*/8,
                      width, filtered->data());
}

// clang-format off
TEST(WienerTest, ExpectOnes_3x3) {
  std::vector<int16_t> filtered;
  ASSERT_WP2_OK(
      FilterCanvas(/*half_tap_weights_h=*/{0, 0, 0},  // = 0.0, 1.0, 0.0
                   /*half_tap_weights_v=*/{   0,
                                              0,
                                              0},
                   3, 3, /*canvas=*/{1, 2, 3,
                                     4, 5, 6,
                                     6, 6, 6}, &filtered));
      EXPECT_THAT(filtered, ElementsAre(1, 2, 3,
                  4, 5, 6,
                  6, 6, 6));  // Expect no change.
}

TEST(WienerTest, ExpectZeros_3x3) {
  const int32_t h = 1 << (kWieFltNumBitsTapWgts - 2);
  std::vector<int16_t> filtered;
  ASSERT_WP2_OK(
      FilterCanvas(/*half_tap_weights_h=*/{0, 0, h},  // = 0.5, 0.0, 0.5
                   /*half_tap_weights_v=*/{   0,
                                              0,
                                              h},
                   3, 3, /*canvas=*/{0, 1, 0,
                                     1, 1, 1,
                                     0, 1, 0}, &filtered));
      EXPECT_THAT(filtered, ElementsAre(0, 0, 0,
                  0, 0, 0,
                  0, 0, 0));  // Expect all values erased.
}

TEST(WienerTest, ExpectBlurH_3x3) {
  const int32_t p = 1 << (kWieFltNumBitsTapWgts - 3);
  std::vector<int16_t> filtered;
  ASSERT_WP2_OK(
      FilterCanvas(/*half_tap_weights_h=*/{0, 0, p},  // = 0.25, 0.5, 0.25
                   /*half_tap_weights_v=*/{   0,
                                              0,
                                              0},
                   3, 3, /*canvas=*/{9, 0, 9,
                                     0, 9, 0,
                                     9, 9, 9}, &filtered));
      EXPECT_THAT(filtered, ElementsAre(6, 4, 6,
                  2, 4, 2,
                  9, 9, 9));  // Horizontal blur.
}

TEST(WienerTest, ExpectBlurV_3x3) {
  const int32_t p = 1 << (kWieFltNumBitsTapWgts - 3);
  std::vector<int16_t> filtered;
  ASSERT_WP2_OK(
      FilterCanvas(/*half_tap_weights_h=*/{0, 0, 0},
                   /*half_tap_weights_v=*/{   0,
                                              0,
                                              p},
                   3, 3, /*canvas=*/{9, 0, 9,
                                     0, 9, 0,
                                     9, 9, 9}, &filtered));
      EXPECT_THAT(filtered, ElementsAre(6, 2, 6,
                  4, 6, 4,
                  6, 9, 6));  // Vertical blur.
}

TEST(WienerTest, ExpectBlur_2x2) {
  const int32_t p = 1 << (kWieFltNumBitsTapWgts - 3);
  std::vector<int16_t> filtered;
  ASSERT_WP2_OK(
      FilterCanvas(/*half_tap_weights_h=*/{p / 8, p / 4, p / 2},
                   /*half_tap_weights_v=*/{       p / 8,
                                                  p / 4,
                                                  p / 2},
                   2, 2, /*canvas=*/{9, 0,
                                     0, 9}, &filtered));
      EXPECT_THAT(filtered, ElementsAre(5, 3,
                  3, 5));  // Blur.
}

TEST(WienerTest, ExpectBlur_8x1) {
  const int32_t p = 1 << (kWieFltNumBitsTapWgts - 3);
  std::vector<int16_t> filtered;
  ASSERT_WP2_OK(
      FilterCanvas(/*half_tap_weights_h=*/{p / 8, p / 4, p / 2},
                   /*half_tap_weights_v=*/{       p / 8,
                                                  p / 4,
                                                  p / 2},
                   8, 1, /*canvas=*/{9, 0, 9, 0, 9, 0, 9, 0}, &filtered));
      EXPECT_THAT(filtered, ElementsAre(7, 3, 6, 2, 6, 2, 5, 1));  // Blur.
}

TEST(WienerTest, ExpectSharpenV_3x5) {
  const int32_t p = 1 << (kWieFltNumBitsTapWgts - 3);
  std::vector<int16_t> filtered;
  ASSERT_WP2_OK(
      FilterCanvas(/*half_tap_weights_h=*/{0, 0, 0},
                   /*half_tap_weights_v=*/{  -p / 8,
                                             -p / 4,
                                             -p / 2},
                   3, 5, /*canvas=*/{3,  3,  3,
                                     0,  3,  3,
                                     3,  3,  3,
                                     3,  3, 10,
                                     3,  3,  3}, &filtered));
      EXPECT_THAT(filtered, ElementsAre( 3,  3,  2,
                  -2,  3,  2,
                   3,  3,  2,
                   3,  3, 13,
                   3,  3,  2));  // Vertical sharpening.
}
// clang-format on

//------------------------------------------------------------------------------

}  // namespace
}  // namespace WP2
