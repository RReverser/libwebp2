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
//
// Test the deblocking filter.

#include <tuple>
#include <limits>

#include "include/helpers.h"
#include "include/helpers_filter.h"
#include "src/dec/filters/deblocking_filter.h"
#include "src/wp2/format_constants.h"

namespace WP2 {
namespace {

//------------------------------------------------------------------------------

// Constants used for the tests below.
constexpr uint32_t kMaxNumBlocksX = 8;
constexpr uint32_t kMaxNumBlocksY = 8;
constexpr uint32_t kMaxPixelsWidth = kMaxNumBlocksX * kMinBlockSizePix;
constexpr uint32_t kMaxPixelsHeight = kMaxNumBlocksY * kMinBlockSizePix;
// Minus pixels to exercise cropped blocks (bottom/right).
constexpr uint32_t kPixelsWidth = kMaxPixelsWidth - 1u;
constexpr uint32_t kPixelsHeight = kMaxPixelsHeight - (kMinBlockSizePix - 1u);

// Creates areas with close values, otherwise the filter will not be applied.
// If with_alpha is true, adds a fully opaque (255) alpha plane.
YUVPlane CreateBlockedPixels(int16_t min_value, int16_t max_value,
                             bool with_alpha = false) {
  YUVPlane pixels;
  WP2_ASSERT_STATUS(
      pixels.Resize(kPixelsWidth, kPixelsHeight, /*pad=*/1, with_alpha));
  const uint32_t half_x = kMaxPixelsWidth / 2;
  const uint32_t half_y = kMaxPixelsHeight / 2;
  // Testing horizontal edge.
  pixels.Y.Fill({0, 0, half_x, kMaxPixelsHeight}, 42);
  pixels.Y.Fill({half_x, 0, half_x, kMaxPixelsHeight}, 38);
  // Testing vertical edge.
  pixels.U.Fill({0, 0, kMaxPixelsWidth, half_y}, min_value);
  pixels.U.Fill({0, half_y, kMaxPixelsWidth, half_y}, (int16_t)(min_value + 4));
  // Testing both.
  pixels.V.Fill({0, 0, half_x, half_y}, (int16_t)(max_value - 7));
  pixels.V.Fill({half_x, 0, half_x, half_y}, max_value);
  pixels.V.Fill({0, half_y, half_x, half_y}, (int16_t)(max_value - 12));
  pixels.V.Fill({half_x, half_y, half_x, half_y}, (int16_t)(max_value - 5));

  if (with_alpha) pixels.A.Fill(kAlphaMax);

  return pixels;
}

FilterBlockMap CreateAlignedBlocks(const Vector<Segment>& segments,
                                   double res_den, double min_bpp,
                                   uint32_t num_precision_bits,
                                   YUVPlane* const pixels) {
  assert((pixels->Y.w_ + kMinBlockSizePix - 1u) / kMinBlockSizePix == 8u);
  assert((pixels->Y.h_ + kMinBlockSizePix - 1u) / kMinBlockSizePix == 8u);
  FilterBlockMap blocks;
  const int32_t min = -(1 << (num_precision_bits - 1)), max = -min - 1;
  blocks.Init(/*tile_rect=*/{0, 0, pixels->Y.w_, pixels->Y.h_},
              num_precision_bits, min, max, pixels);
  EXPECT_WP2_OK(blocks.Allocate());
  const uint32_t s = segments.size();
  uint32_t id = 0;  // Segment id
  testing::RegisterBlock(0, 0, BLK_16x8, id++ % s, res_den, min_bpp,
                         /*with_lossy_alpha=*/true, &blocks);
  testing::RegisterBlock(4, 0, BLK_16x8, id++ % s, res_den, min_bpp,
                         /*with_lossy_alpha=*/true, &blocks);
  testing::RegisterBlock(0, 2, BLK_8x16, id++ % s, res_den, min_bpp,
                         /*with_lossy_alpha=*/true, &blocks);
  testing::RegisterBlock(2, 2, BLK_16x16, id++ % s, res_den, min_bpp,
                         /*with_lossy_alpha=*/true, &blocks);
  testing::RegisterBlock(6, 2, BLK_8x4, id++ % s, res_den, min_bpp,
                         /*with_lossy_alpha=*/true, &blocks);
  testing::RegisterBlock(6, 3, BLK_4x4, id++ % s, res_den, min_bpp,
                         /*with_lossy_alpha=*/false, &blocks);
  testing::RegisterBlock(7, 3, BLK_4x4, id++ % s, res_den, min_bpp,
                         /*with_lossy_alpha=*/true, &blocks);
  testing::RegisterBlock(6, 4, BLK_8x16, id++ % s, res_den, min_bpp,
                         /*with_lossy_alpha=*/true, &blocks);
  testing::RegisterBlock(0, 6, BLK_8x8, id++ % s, res_den, min_bpp,
                         /*with_lossy_alpha=*/true, &blocks);
  testing::RegisterBlock(2, 6, BLK_16x8, id++ % s, res_den, min_bpp,
                         /*with_lossy_alpha=*/false, &blocks);
  return blocks;
}

// Applies the deblocking filter on 'pixels' with 'filter_strength'.
void Deblock(uint32_t num_precision_bits, int16_t min_value, int16_t max_value,
             uint32_t quality_hint, YUVPlane* const pixels) {
  const Vector<Segment> segments =
      testing::CreateSegments(quality_hint, min_value, max_value);
  const double residual_density = 1. * quality_hint / kMaxLossyQualityHint;
  const double min_bpp = 4. * quality_hint / kMaxLossyQualityHint;
  // The created block edges should match at least some of the boundaries of
  // the areas defined in CreateBlockedPixels().
  const FilterBlockMap blocks = CreateAlignedBlocks(
      segments, residual_density, min_bpp, num_precision_bits, pixels);

  DecoderConfig config = DecoderConfig::kDefault;
  config.enable_deblocking_filter = true;
  BitstreamFeatures features;
  features.quality_hint = quality_hint;
  features.tile_width = features.tile_height = 256;
  DeblockingFilter deblocking_filter(config, features, blocks);
  ASSERT_WP2_OK(deblocking_filter.Allocate());

  uint32_t num_deblocked_rows = 0;
  for (uint32_t num_decoded_yuv_rows = 0; num_decoded_yuv_rows <= pixels->Y.h_;
       ++num_decoded_yuv_rows) {
    const uint32_t n = deblocking_filter.Deblock(num_decoded_yuv_rows);
    ASSERT_GE(n, num_deblocked_rows);
    num_deblocked_rows = n;
  }
  ASSERT_EQ(num_deblocked_rows, pixels->Y.h_);
}

//------------------------------------------------------------------------------

// Param: num_precision_bits
class DeblockingFilterTest : public ::testing::TestWithParam<uint32_t> {};

// Tests that no pixel is modified if each plane has one unique color.
TEST_P(DeblockingFilterTest, FilteringPlain) {
  const uint32_t num_precision_bits = GetParam();
  const int16_t min_value = testing::GetMinValueSigned(num_precision_bits);
  const int16_t max_value = testing::GetMaxValueSigned(num_precision_bits);

  YUVPlane pixels = testing::CreatePlainPixels(kPixelsWidth, kPixelsHeight,
                                               min_value, max_value);
  YUVPlane reference;
  ASSERT_WP2_OK(reference.Copy(pixels, /*resize_if_needed=*/true));

  Deblock(num_precision_bits, min_value, max_value, /*quality_hint=*/0,
          &pixels);

  EXPECT_TRUE(testing::AreEqual(pixels.Y, reference.Y));
  EXPECT_TRUE(testing::AreEqual(pixels.U, reference.U));
  EXPECT_TRUE(testing::AreEqual(pixels.V, reference.V));
}

// Tests that some pixels are deblocked for horizontal and vertical edges.
TEST_P(DeblockingFilterTest, MaxFilteringBlocked) {
  const uint32_t num_precision_bits = GetParam();
  const int16_t min_value = testing::GetMinValueSigned(num_precision_bits);
  const int16_t max_value = testing::GetMaxValueSigned(num_precision_bits);

  YUVPlane pixels =
      CreateBlockedPixels(min_value, max_value, /*with_alpha=*/false);
  YUVPlane reference;
  ASSERT_WP2_OK(reference.Copy(pixels, /*resize_if_needed=*/true));

  Deblock(num_precision_bits, min_value, max_value, /*quality_hint=*/0,
          &pixels);

  EXPECT_FALSE(testing::AreEqual(pixels.Y, reference.Y));
  EXPECT_FALSE(testing::AreEqual(pixels.U, reference.U));
  EXPECT_FALSE(testing::AreEqual(pixels.V, reference.V));
}

// Tests that some pixels are deblocked for horizontal and vertical edges
// when the alpha plane is constant.
TEST_P(DeblockingFilterTest, MaxFilteringBlockedWithPlainAlpha) {
  const uint32_t num_precision_bits = GetParam();
  const int16_t min_value = testing::GetMinValueSigned(num_precision_bits);
  const int16_t max_value = testing::GetMaxValueSigned(num_precision_bits);

  YUVPlane pixels =
      CreateBlockedPixels(min_value, max_value, /*with_alpha=*/true);
  YUVPlane reference;
  ASSERT_WP2_OK(reference.Copy(pixels, /*resize_if_needed=*/true));

  Deblock(num_precision_bits, min_value, max_value, /*quality_hint=*/0,
          &pixels);

  EXPECT_FALSE(testing::AreEqual(pixels.Y, reference.Y));
  EXPECT_FALSE(testing::AreEqual(pixels.U, reference.U));
  EXPECT_FALSE(testing::AreEqual(pixels.V, reference.V));
  EXPECT_TRUE(testing::AreEqual(pixels.A, reference.A));
}

// Tests that some pixels are deblocked for horizontal and vertical edges
// when the alpha plane has low variance. Alpha plane should also be deblocked.
TEST_P(DeblockingFilterTest, MaxFilteringBlockedWithLowVarianceAlpha) {
  const uint32_t num_precision_bits = GetParam();
  const int16_t min_value = testing::GetMinValueSigned(num_precision_bits);
  const int16_t max_value = testing::GetMaxValueSigned(num_precision_bits);

  YUVPlane pixels =
      CreateBlockedPixels(min_value, max_value, /*with_alpha=*/true);
  // The rest is 255.
  pixels.A.Fill({0, 0, kMaxPixelsWidth / 2, kMaxPixelsWidth / 2}, 250);
  YUVPlane reference;
  ASSERT_WP2_OK(reference.Copy(pixels, /*resize_if_needed=*/true));

  Deblock(num_precision_bits, min_value, max_value, /*quality_hint=*/0,
          &pixels);

  EXPECT_FALSE(testing::AreEqual(pixels.Y, reference.Y));
  EXPECT_FALSE(testing::AreEqual(pixels.U, reference.U));
  EXPECT_FALSE(testing::AreEqual(pixels.V, reference.V));
  EXPECT_FALSE(testing::AreEqual(pixels.A, reference.A));
}

// Tests that pixels are not deblocked when alpha plane has high variance.
TEST_P(DeblockingFilterTest, MaxFilteringBlockedWithHighVarianceAlpha) {
  const uint32_t num_precision_bits = GetParam();
  const int16_t min_value = testing::GetMinValueSigned(num_precision_bits);
  const int16_t max_value = testing::GetMaxValueSigned(num_precision_bits);

  YUVPlane pixels =
      CreateBlockedPixels(min_value, max_value, /*with_alpha=*/true);
  const uint32_t half_x = kMaxPixelsWidth / 2;
  const uint32_t half_y = kMaxPixelsHeight / 2;
  // Testing horizontal edge.
  pixels.A.Fill({0, 0, half_x, half_y}, 50);
  pixels.A.Fill({half_x, half_y, half_x, half_y}, 50);

  YUVPlane reference;
  ASSERT_WP2_OK(reference.Copy(pixels, /*resize_if_needed=*/true));

  Deblock(num_precision_bits, min_value, max_value, /*quality_hint=*/0,
          &pixels);

  EXPECT_TRUE(testing::AreEqual(pixels.Y, reference.Y));
  EXPECT_TRUE(testing::AreEqual(pixels.U, reference.U));
  EXPECT_TRUE(testing::AreEqual(pixels.V, reference.V));
  EXPECT_TRUE(testing::AreEqual(pixels.A, reference.A));
}

// Tests that no pixel is modified if it's the highest quality.
TEST_P(DeblockingFilterTest, NoFiltering) {
  const uint32_t num_precision_bits = GetParam();
  const int16_t min_value = testing::GetMinValueSigned(num_precision_bits);
  const int16_t max_value = testing::GetMaxValueSigned(num_precision_bits);

  YUVPlane pixels = CreateBlockedPixels(min_value, max_value);
  YUVPlane reference;
  ASSERT_WP2_OK(reference.Copy(pixels, /*resize_if_needed=*/true));

  Deblock(num_precision_bits, min_value, max_value,
          /*quality_hint=*/kMaxLossyQualityHint, &pixels);

  EXPECT_TRUE(testing::AreEqual(pixels.Y, reference.Y));
  EXPECT_TRUE(testing::AreEqual(pixels.U, reference.U));
  EXPECT_TRUE(testing::AreEqual(pixels.V, reference.V));
}

INSTANTIATE_TEST_SUITE_P(
    DeblockingFilterTestInstantiation, DeblockingFilterTest,
    /*num_precision_bits=*/::testing::Values(8u, 10u, 12u));

//------------------------------------------------------------------------------

// Tests that PSNR is good even if there is only deblocking.
TEST(DeblockingFilterTest, Simple) {
  EncoderConfig encoder_config = EncoderConfig::kDefault;
  encoder_config.quality = 40;  // Expect compression artifacts.

  DecoderConfig decoder_config = DecoderConfig::kDefault;
  decoder_config.enable_deblocking_filter = true;
  decoder_config.enable_directional_filter = false;
  decoder_config.enable_restoration_filter = false;

  ASSERT_WP2_OK(testing::EncodeDecodeCompare("source4.webp", encoder_config,
                                             decoder_config));
}

// Tests that PSNR is better with the deblocking filter than without.
TEST(DeblockingFilterTest, BetterPSNR) {
  EncoderConfig encoder_config = EncoderConfig::kDefault;
  encoder_config.quality = 50;  // Expect compression artifacts.

  DecoderConfig decoder_config_no_filter = DecoderConfig::kDefault;
  decoder_config_no_filter.enable_deblocking_filter = false;

  DecoderConfig decoder_config_filter = decoder_config_no_filter;
  decoder_config_filter.enable_deblocking_filter = true;

  float disto_diff[5];
  ASSERT_WP2_OK(
      testing::GetDistortionDiff("source1.png", {25, 1, 150, 200},
                                 encoder_config, decoder_config_no_filter,
                                 decoder_config_filter, disto_diff));
  EXPECT_GT(disto_diff[4], 0.f);
}

//------------------------------------------------------------------------------

TEST(HelpersFilterTest, MinMaxValueSigned) {
  EXPECT_EQ(testing::GetMinValueSigned(/*num_precision_bits=*/8),
            std::numeric_limits<int8_t>::min());
  EXPECT_EQ(testing::GetMaxValueSigned(/*num_precision_bits=*/8),
            std::numeric_limits<int8_t>::max());

  EXPECT_EQ(testing::GetMinValueSigned(/*num_precision_bits=*/10), -512);
  EXPECT_EQ(testing::GetMaxValueSigned(/*num_precision_bits=*/10), 511);

  EXPECT_EQ(testing::GetMinValueSigned(/*num_precision_bits=*/16),
            std::numeric_limits<int16_t>::min());
  EXPECT_EQ(testing::GetMaxValueSigned(/*num_precision_bits=*/16),
            std::numeric_limits<int16_t>::max());
}

//------------------------------------------------------------------------------

class FilterTestCpuInfo : public ::testing::TestWithParam<WP2CPUInfo> {
  void SetUp() override {
    WP2DspReset();
    WP2GetCPUInfo = GetParam();
    DblkFilterInit();
  }
};

TEST_P(FilterTestCpuInfo, TestMeasureLength) {
  UniformIntDistribution gen(52452);
  const uint32_t kNumTests = 7000;
  const uint32_t kStep = 13;
  int16_t A[kStep * 16], *B = A + kStep * 8;
  uint32_t crc = 323452;
  for (uint32_t test = 0; test < kNumTests; ++test) {
    for (auto& a : A) a = gen.Get(-1024, 1024);
    for (uint32_t half = 1; half <= kDblkMaxHalf; ++half) {
      for (uint32_t s = 0; s <= kDblkMaxSharpness; ++s) {
        for (uint32_t step = 0; step <= kStep; ++step) {
          const uint32_t len = MeasureFlatLength(s, half, B, s);
          crc = testing::SillyCRC(crc, len);
        }
      }
    }
  }
  EXPECT_EQ(2013140036u, crc);
}

INSTANTIATE_TEST_SUITE_P(FilterTestCpuInfoInstantiation, FilterTestCpuInfo,
                         ::testing::ValuesIn(testing::kWP2CpuInfos));

//------------------------------------------------------------------------------

}  // namespace
}  // namespace WP2
