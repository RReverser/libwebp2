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

#include <cmath>
#include <tuple>

#include "include/helpers.h"
#include "include/helpers_filter.h"
#include "src/dec/filters/directional_filter.h"
#include "src/dsp/math.h"
#include "src/utils/utils.h"
#include "src/wp2/format_constants.h"

namespace WP2 {
namespace {

//------------------------------------------------------------------------------

const double kPi = std::acos(-1);

// Fills the 'plane' with a sinusoidal pattern oriented at 'angle'.
void Wave(int16_t min_value, int16_t max_value, double angle, double wavelength,
          Plane16* const plane) {
  angle += kPi * 0.5;  // 'angle' is for pattern orientation, not gradient.
  const double cosinus = std::cos(angle), sinus = std::sin(angle);
  wavelength *= kPi * 2.;  // For cos(x/wavelength) step.
  for (uint32_t y = 0; y < plane->h_; ++y) {
    for (uint32_t x = 0; x < plane->w_; ++x) {
      const double oriented_x = x * cosinus + y * sinus;
      const double amplitude =
          std::cos(oriented_x / wavelength) * 0.5 * (max_value - min_value);
      plane->At(x, y) =
          (int16_t)std::lround(amplitude + (min_value + max_value) * 0.5);
      assert(plane->At(x, y) >= min_value && plane->At(x, y) <= max_value);
    }
  }
}

YUVPlane CreateNoisyWavyPixels(uint32_t width, uint32_t height,
                               int16_t min_value, int16_t max_value) {
  YUVPlane pixels;
  WP2_ASSERT_STATUS(pixels.Resize(width, height));
  Wave(min_value, max_value, kPi * 0.0, pixels.Y.w_, &pixels.Y);  // Horizontal.
  Wave(min_value, max_value, kPi * 0.5, pixels.U.h_, &pixels.U);  // Vertical.
  Wave(min_value, max_value, kPi * 2.1, 10., &pixels.V);          // Something.
  for (Channel channel : {kYChannel, kUChannel, kVChannel, kAChannel}) {
    testing::Noise(min_value, max_value, /*seed=*/channel, /*strength=*/5,
                   &pixels.GetChannel(channel));
  }
  return pixels;
}

// Applies the directional filter on 'pixels' with 'filter_strength'.
void Filter(uint32_t num_precision_bits, int16_t min_value, int16_t max_value,
            uint32_t quality_hint, YUVPlane* const pixels) {
  const Vector<Segment> segments =
      testing::CreateSegments(quality_hint, min_value, max_value);
  const double residual_density = 1. * quality_hint / kMaxLossyQualityHint;
  const double min_bpp = 4. * quality_hint / kMaxLossyQualityHint;
  FilterBlockMap blocks;
  testing::CreateBlocks(/*tile_x=*/0, /*tile_y=*/0, segments, residual_density,
                        min_bpp, num_precision_bits, pixels, &blocks);

  DecoderConfig config = DecoderConfig::kDefault;
  config.enable_directional_filter = true;
  BitstreamFeatures features;
  features.quality_hint = quality_hint;
  DirectionalFilter directional_filter(config, features, blocks);
  ASSERT_WP2_OK(directional_filter.Allocate());

  uint32_t num_deblocked_rows = 0;
  for (uint32_t num_decoded_yuv_rows = 0; num_decoded_yuv_rows <= pixels->Y.h_;
       ++num_decoded_yuv_rows) {
    const uint32_t n = directional_filter.Smooth(num_decoded_yuv_rows);
    ASSERT_GE(n, num_deblocked_rows);
    num_deblocked_rows = n;
  }
  ASSERT_EQ(num_deblocked_rows, pixels->Y.h_);
}

//------------------------------------------------------------------------------

class DirectionalFilterTest : public ::testing::TestWithParam<
                                  std::tuple<uint32_t, uint32_t, uint32_t>> {};

// Tests that no pixel is modified if each plane has one unique color.
TEST_P(DirectionalFilterTest, FilteringPlain) {
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
                                 /*quality_hint=*/0, &pixels));

  EXPECT_TRUE(testing::AreEqual(pixels.Y, reference.Y));
  EXPECT_TRUE(testing::AreEqual(pixels.U, reference.U));
  EXPECT_TRUE(testing::AreEqual(pixels.V, reference.V));
}

// Tests that pixels are filtered.
TEST_P(DirectionalFilterTest, MaxFilteringBlocked) {
  const uint32_t num_precision_bits = std::get<0>(GetParam());
  const int16_t min_value = testing::GetMinValueSigned(num_precision_bits);
  const int16_t max_value = testing::GetMaxValueSigned(num_precision_bits);
  const uint32_t width = std::get<1>(GetParam());
  const uint32_t height = std::get<2>(GetParam());

  YUVPlane pixels = CreateNoisyWavyPixels(width, height, min_value, max_value);
  YUVPlane reference;
  ASSERT_WP2_OK(reference.Copy(pixels, /*resize_if_needed=*/true));

  ASSERT_NO_FATAL_FAILURE(Filter(num_precision_bits, min_value, max_value,
                                 /*quality_hint=*/0, &pixels));

  EXPECT_FALSE(testing::AreEqual(pixels.Y, reference.Y));
  EXPECT_FALSE(testing::AreEqual(pixels.U, reference.U));
  EXPECT_FALSE(testing::AreEqual(pixels.V, reference.V));
}

// Tests that no pixel is modified if it's the highest quality.
TEST_P(DirectionalFilterTest, NoFiltering) {
  const uint32_t num_precision_bits = std::get<0>(GetParam());
  const int16_t min_value = testing::GetMinValueSigned(num_precision_bits);
  const int16_t max_value = testing::GetMaxValueSigned(num_precision_bits);
  const uint32_t width = std::get<1>(GetParam());
  const uint32_t height = std::get<2>(GetParam());

  YUVPlane pixels = CreateNoisyWavyPixels(width, height, min_value, max_value);
  YUVPlane reference;
  ASSERT_WP2_OK(reference.Copy(pixels, /*resize_if_needed=*/true));

  ASSERT_NO_FATAL_FAILURE(Filter(num_precision_bits, min_value, max_value,
                                 /*quality_hint=*/kMaxLossyQualityHint,
                                 &pixels));

  EXPECT_TRUE(testing::AreEqual(pixels.Y, reference.Y));
  EXPECT_TRUE(testing::AreEqual(pixels.U, reference.U));
  EXPECT_TRUE(testing::AreEqual(pixels.V, reference.V));
}

INSTANTIATE_TEST_SUITE_P(
    DirectionalFilterTestInstantiation, DirectionalFilterTest,
    ::testing::Combine(
        /*num_precision_bits=*/::testing::Values(8u, 10u, 12u),
        /*width=*/::testing::Values(kDrctFltSize, 3 * kDrctFltSize - 1u),
        /*height=*/::testing::Values(kDrctFltSize, 3 * kDrctFltSize - 1u)));

//------------------------------------------------------------------------------

// Tests that PSNR is good even if there is only directional filtering.
TEST(DirectionalFilterTest, Simple) {
  EncoderConfig encoder_config = EncoderConfig::kDefault;
  encoder_config.quality = 40;  // Expect compression artifacts.

  DecoderConfig decoder_config = DecoderConfig::kDefault;
  decoder_config.enable_deblocking_filter = false;
  decoder_config.enable_directional_filter = true;
  decoder_config.enable_restoration_filter = false;

  ASSERT_WP2_OK(testing::EncodeDecodeCompare("source4.webp", encoder_config,
                                             decoder_config));
}

// Tests that PSNR is better with the directional filter than without.
TEST(DirectionalFilterTest, BetterPSNR) {
  EncoderConfig encoder_config = EncoderConfig::kDefault;
  encoder_config.quality = 50;  // Expect compression artifacts.

  DecoderConfig decoder_config_no_filter = DecoderConfig::kDefault;
  decoder_config_no_filter.enable_directional_filter = false;

  DecoderConfig decoder_config_filter = decoder_config_no_filter;
  decoder_config_filter.enable_directional_filter = true;

  float disto_diff[5];
  ASSERT_WP2_OK(
      testing::GetDistortionDiff("source1.png", {1, 42, 300, 66},
                                 encoder_config, decoder_config_no_filter,
                                 decoder_config_filter, disto_diff));
  EXPECT_GT(disto_diff[4], 0.f);
}

//------------------------------------------------------------------------------

}  // namespace
}  // namespace WP2
