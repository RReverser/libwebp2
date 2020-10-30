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

// Test the alpha filter.

#include "examples/example_utils.h"
#include "include/helpers.h"
#include "include/helpers_filter.h"
#include "src/common/integral.h"
#include "src/dec/filters/alpha_filter.h"

namespace WP2 {
namespace {

constexpr uint32_t kNumBits = 8;
constexpr int32_t kYUVMin = -(1 << (kNumBits - 1));
constexpr int32_t kYUVMax = (1 << (kNumBits - 1)) - 1;

void RegisterBlock(uint16_t x, uint16_t y, BlockSize dim, bool has_lossy_alpha,
                   FilterBlockMap* const blocks) {
  testing::RegisterBlock(x, y, dim, /*segment_id=*/0, /*res_den=*/0.5f,
                         /*min_bpp=*/1.f, has_lossy_alpha, blocks);
}

void Filter(const FilterBlockMap& blocks, uint32_t h) {
  GlobalParams gparams;
  ASSERT_TRUE(gparams.a_segments_.resize(1));
  gparams.a_segments_[0].quality_factor_ = 30;
  gparams.has_alpha_ = true;
  gparams.maybe_use_lossy_alpha_ = true;
  gparams.enable_alpha_filter_ = true;
  AlphaFilter alpha_filter(DecoderConfig::kDefault, gparams, blocks);
  ASSERT_WP2_OK(alpha_filter.Allocate());
  ASSERT_TRUE(IsAlphaFilterEnabled(DecoderConfig::kDefault, gparams));

  uint32_t num_filtered_rows = 0;
  for (uint32_t num_decoded_yuv_rows = 0; num_decoded_yuv_rows <= h;
       ++num_decoded_yuv_rows) {
    const uint32_t n = alpha_filter.Smooth(num_decoded_yuv_rows);
    ASSERT_GE(n, num_filtered_rows);
    num_filtered_rows = n;
  }
  ASSERT_EQ(num_filtered_rows, h);
}

//------------------------------------------------------------------------------

TEST(AlphaFilterTest, DoesntFilterLossless) {
  YUVPlane pixels;
  ASSERT_WP2_OK(pixels.Resize(32, 16, /*pad=*/1, /*has_alpha=*/true));
  pixels.A.Fill(0);
  testing::Noise(0, kAlphaMax, /*seed=*/0, /*strength=*/kAlphaMax / 2,
                 &pixels.A);
  YUVPlane reference;
  ASSERT_WP2_OK(reference.Copy(pixels, /*resize_if_needed=*/true));

  // ______________________
  // |         | 8x8 |8x16|
  // |  16x16  |_____|    | all blocks lossless
  // |         | 8x8 |    |
  // | ________|_____|____|
  FilterBlockMap blocks;
  blocks.Init(/*tile_rect=*/{0, 0, 32, 16}, kNumBits, kYUVMin, kYUVMax,
              &pixels);
  ASSERT_WP2_OK(blocks.Allocate());
  RegisterBlock(0, 0, BLK_16x16, /*has_lossy_alpha=*/false, &blocks);
  RegisterBlock(4, 0, BLK_8x8, /*has_lossy_alpha=*/false, &blocks);
  RegisterBlock(6, 0, BLK_8x16, /*has_lossy_alpha=*/false, &blocks);
  RegisterBlock(4, 2, BLK_8x8, /*has_lossy_alpha=*/false, &blocks);

  ASSERT_NO_FATAL_FAILURE(Filter(blocks, pixels.A.h_));

  EXPECT_TRUE(testing::AreEqual(pixels.A, reference.A));
}

TEST(AlphaFilterTest, DoesFilterLossyBlocks) {
  YUVPlane pixels;
  const uint32_t width = 32;
  const uint32_t height = 16;
  ASSERT_WP2_OK(pixels.Resize(width, height, /*pad=*/1, /*has_alpha=*/true));
  pixels.A.Fill(0);
  testing::Noise(0, kAlphaMax, /*seed=*/0, /*strength=*/kAlphaMax / 2,
                 &pixels.A);
  YUVPlane reference;
  ASSERT_WP2_OK(reference.Copy(pixels, /*resize_if_needed=*/true));

  // ______________________
  // |         | 8x8 |8x16|
  // |  16x16  |_____|    |
  // |  lossy  | 8x8 |    |
  // | ________|lossy|____|
  FilterBlockMap blocks;
  blocks.Init(/*tile_rect=*/{0, 0, width, height}, kNumBits, kYUVMin, kYUVMax,
              &pixels);
  ASSERT_WP2_OK(blocks.Allocate());
  RegisterBlock(0, 0, BLK_16x16, /*has_lossy_alpha=*/true, &blocks);
  RegisterBlock(4, 0, BLK_8x8, /*has_lossy_alpha=*/false, &blocks);
  RegisterBlock(6, 0, BLK_8x16, /*has_lossy_alpha=*/false, &blocks);
  RegisterBlock(4, 2, BLK_8x8, /*has_lossy_alpha=*/true, &blocks);

  ASSERT_NO_FATAL_FAILURE(Filter(blocks, pixels.A.h_));

  // Lossless blocks should be unchanged, lossy blocks should be changed.
  YUVPlane ref_view, pixels_view;
  ASSERT_WP2_OK(ref_view.SetView(reference, {0, 0, 16, 16}));
  ASSERT_WP2_OK(pixels_view.SetView(pixels, {0, 0, 16, 16}));
  EXPECT_FALSE(testing::AreEqual(pixels_view.A, ref_view.A));
  ASSERT_WP2_OK(ref_view.SetView(reference, {16, 0, 8, 8}));
  ASSERT_WP2_OK(pixels_view.SetView(pixels, {16, 0, 8, 8}));
  EXPECT_TRUE(testing::AreEqual(pixels_view.A, ref_view.A));
  ASSERT_WP2_OK(ref_view.SetView(reference, {24, 0, 8, 16}));
  ASSERT_WP2_OK(pixels_view.SetView(pixels, {24, 0, 8, 16}));
  EXPECT_TRUE(testing::AreEqual(pixels_view.A, ref_view.A));
  ASSERT_WP2_OK(ref_view.SetView(reference, {16, 8, 8, 8}));
  ASSERT_WP2_OK(pixels_view.SetView(pixels, {16, 8, 8, 8}));
  EXPECT_FALSE(testing::AreEqual(pixels_view.A, ref_view.A));

  const uint32_t block_size = 4;
  Integral ref_integral;
  ASSERT_WP2_OK(ref_integral.Allocate(8, block_size, block_size));
  ref_integral.AddValues(reference.A);
  Integral pixels_integral;
  ASSERT_WP2_OK(pixels_integral.Allocate(8, block_size, block_size));
  pixels_integral.AddValues(pixels.A);
  // Alpha filter is a blurring filter so it should lower the variance.
  EXPECT_LT(
      pixels_integral.StdDev(0, 0, width / block_size, height / block_size),
      ref_integral.StdDev(0, 0, width / block_size, height / block_size));
}

TEST(AlphaFilterTest, PreservesHardEdges) {
  YUVPlane pixels;
  const uint32_t width = 32;
  const uint32_t height = 32;
  ASSERT_WP2_OK(pixels.Resize(width, height, /*pad=*/1, /*has_alpha=*/true));
  // Fill half with 0, half with max. This is a very hard edge and should be
  // preserved.
  pixels.A.Fill({0, 0, width / 2, height}, /*value=*/0);
  pixels.A.Fill({width / 2, 0, width / 2, height}, /*value=*/kAlphaMax);
  YUVPlane reference;
  ASSERT_WP2_OK(reference.Copy(pixels, /*resize_if_needed=*/true));

  FilterBlockMap blocks;
  blocks.Init(/*tile_rect=*/{0, 0, width, height}, kNumBits, kYUVMin, kYUVMax,
              &pixels);
  ASSERT_WP2_OK(blocks.Allocate());
  RegisterBlock(0, 0, BLK_32x32, /*has_lossy_alpha=*/true, &blocks);

  ASSERT_NO_FATAL_FAILURE(Filter(blocks, pixels.A.h_));
  EXPECT_TRUE(testing::AreEqual(pixels.A, reference.A));
}

TEST(AlphaFilterTest, BlursSoftEdges) {
  YUVPlane pixels;
  const uint32_t width = 32;
  const uint32_t height = 32;
  ASSERT_WP2_OK(pixels.Resize(width, height, /*pad=*/1, /*has_alpha=*/true));
  // Fill left half with a value, right half with a slightly large one.
  // This is a soft edge which should get blurred.
  const uint32_t left_value = 25;
  const uint32_t right_value = left_value + 10;
  pixels.A.Fill({0, 0, width / 2, height}, left_value);
  pixels.A.Fill({width / 2, 0, width / 2, height}, right_value);
  YUVPlane reference;
  ASSERT_WP2_OK(reference.Copy(pixels, /*resize_if_needed=*/true));

  FilterBlockMap blocks;
  blocks.Init(/*tile_rect=*/{0, 0, width, height}, kNumBits, kYUVMin, kYUVMax,
              &pixels);
  ASSERT_WP2_OK(blocks.Allocate());
  RegisterBlock(0, 0, BLK_32x32, /*has_lossy_alpha=*/true, &blocks);

  ASSERT_NO_FATAL_FAILURE(Filter(blocks, pixels.A.h_));
  EXPECT_FALSE(testing::AreEqual(pixels.A, reference.A));

  // Expect horizontally increasing values.
  for (uint32_t x = 1; x < pixels.A.w_; ++x) {
    SCOPED_TRACE(SPrintf("x %d", x));
    const int16_t v = pixels.A.At(x, 0);
    EXPECT_GE(v, pixels.A.At(x - 1, 0));
    for (uint32_t y = 0; y < pixels.A.h_; ++y) {
      SCOPED_TRACE(SPrintf("y %d", y));
      EXPECT_EQ(pixels.A.At(x, y), v);
    }
  }
}

TEST(AlphaFilterTest, BetterPSNR) {
  EncoderConfig encoder_config_no_filter = EncoderConfig::kDefault;
  encoder_config_no_filter.alpha_quality = 50;  // Expect compression artifacts.
  encoder_config_no_filter.enable_alpha_filter = false;

  EncoderConfig encoder_config_filter = encoder_config_no_filter;
  encoder_config_filter.enable_alpha_filter = true;

  DecoderConfig decoder_config = DecoderConfig::kDefault;

  float disto_diff[5];
  ASSERT_WP2_OK(testing::GetDistortionDiff(
      "source1.png", {25, 27, 200, 208}, encoder_config_no_filter,
      encoder_config_filter, decoder_config, disto_diff));
  // PSNR should be better with the alpha filter than without.
  EXPECT_GT(disto_diff[4], 0.f);
}

}  // namespace
}  // namespace WP2
