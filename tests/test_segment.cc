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

// Tests for histograms and segments.

#include <cstdint>

#include "include/helpers.h"
#include "src/common/constants.h"
#include "src/common/lossy/block.h"
#include "src/enc/wp2_enc_i.h"
#include "src/utils/plane.h"
#include "src/wp2/encode.h"

namespace WP2 {
namespace {

using ::testing::ElementsAre;

TEST(ClusterHistogramTest, Simple) {
  const uint32_t histogram[10] = {1, 10, 2, 0, 1, 5, 6, 50, 3, 10};
  uint32_t result[10];

  uint32_t max_clusters = 1;
  uint32_t num_clusters;
  uint32_t center_positions[kMaxNumSegments];
  num_clusters = ClusterHistogram(histogram, 10, max_clusters, result,
                                  center_positions);
  EXPECT_THAT(result, ElementsAre(0, 0, 0, 0, 0, 0, 0, 0, 0, 0));
  EXPECT_EQ(max_clusters, num_clusters);
  EXPECT_EQ(6u, center_positions[0]);

  max_clusters = 2;
  num_clusters = ClusterHistogram(histogram, 10, max_clusters, result,
                                  center_positions);
  EXPECT_THAT(result, ElementsAre(0, 0, 0, 0, 0, 1, 1, 1, 1, 1));
  EXPECT_EQ(max_clusters, num_clusters);
  EXPECT_EQ(1u, center_positions[0]);
  EXPECT_EQ(7u, center_positions[1]);

  max_clusters = 3;
  num_clusters = ClusterHistogram(histogram, 10, max_clusters, result,
                                  center_positions);
  EXPECT_THAT(result, ElementsAre(0, 0, 0, 0, 1, 1, 1, 2, 2, 2));
  EXPECT_EQ(max_clusters, num_clusters);
  EXPECT_EQ(1u, center_positions[0]);
  EXPECT_EQ(5u, center_positions[1]);
  EXPECT_EQ(7u, center_positions[2]);

  max_clusters = 4;
  num_clusters = ClusterHistogram(histogram, 10, max_clusters, result,
                                  center_positions);
  EXPECT_THAT(result, ElementsAre(0, 0, 0, 1, 1, 1, 2, 2, 3, 3));
  EXPECT_EQ(max_clusters, num_clusters);
  EXPECT_EQ(1u, center_positions[0]);
  EXPECT_EQ(4u, center_positions[1]);
  EXPECT_EQ(6u, center_positions[2]);
  EXPECT_EQ(8u, center_positions[3]);
}

TEST(ClusterHistogramTest, FewNonZeroValues) {
  const uint32_t histogram[10] = {1, 0, 0, 0, 0, 0, 0, 2, 5, 6};
  uint32_t result[10];

  uint32_t max_clusters = 2;
  uint32_t num_clusters;
  uint32_t center_positions[kMaxNumSegments];
  num_clusters = ClusterHistogram(histogram, 10, max_clusters, result,
                                  center_positions);
  EXPECT_THAT(result, ElementsAre(0, 0, 0, 0, 0, 1, 1, 1, 1, 1));
  EXPECT_EQ(max_clusters, num_clusters);
  EXPECT_EQ(0u, center_positions[0]);
  EXPECT_EQ(8u, center_positions[1]);

  max_clusters = 3;
  num_clusters = ClusterHistogram(histogram, 10, max_clusters, result,
                                  center_positions);
  EXPECT_THAT(result, ElementsAre(0, 0, 0, 0, 1, 1, 1, 1, 1, 2));
  EXPECT_EQ(max_clusters, num_clusters);
  EXPECT_EQ(0u, center_positions[0]);
  EXPECT_EQ(7u, center_positions[1]);
  EXPECT_EQ(9u, center_positions[2]);

  max_clusters = 4;
  num_clusters = ClusterHistogram(histogram, 10, max_clusters, result,
                                  center_positions);
  EXPECT_THAT(result, ElementsAre(0, 0, 0, 0, 1, 1, 1, 1, 1, 2));
  EXPECT_EQ(max_clusters - 1, num_clusters);   // one less cluster than expected
  EXPECT_EQ(0u, center_positions[0]);
  EXPECT_EQ(7u, center_positions[1]);
  EXPECT_EQ(9u, center_positions[2]);
}

Vector<CodedBlock> SplitInBlocks(const YUVPlane& yuv) {
  Vector<CodedBlock> cblocks;
  const uint32_t bw = SizeBlocks(yuv.Y.w_);
  const uint32_t bh = SizeBlocks(yuv.Y.h_);
  EXPECT_TRUE(cblocks.resize(bw * bh));
  uint32_t n = 0;
  for (uint32_t y = 0; y < bh; ++y) {
    for (uint32_t x = 0; x < bw; ++x) {
      if (x + 4 <= bw && y + 4 <= bh) {
        CodedBlock* const cb = &cblocks[n++];
        cb->SetDimDefault({x, y, BLK_4x4});
      }
    }
  }
  EXPECT_TRUE(cblocks.resize(n));
  return cblocks;
}

TEST(TestSegment, Simple) {
  // A uniform 10x10 image with a square of a different color in the bottom
  // right quarter.
  YUVPlane yuv;
  ASSERT_WP2_OK(yuv.Resize(10, 10));
  yuv.Fill(Ayuv38b{kAlphaMax, 0, 0, 0});
  yuv.Y.Fill({5, 5, 5, 5}, 100);

  EncoderConfig config;
  config.segments = 3;
  Vector<CodedBlock> cblocks = SplitInBlocks(yuv);
  GlobalParams gparams;
  FeatureMap map;
  gparams.features_ = &map;
  ASSERT_WP2_OK(FindSegments(yuv, config, &gparams));
  EXPECT_LE(gparams.segments_.size(), 3u);
}

TEST(TestSegment, OneSegment) {
  YUVPlane yuv;
  ASSERT_WP2_OK(yuv.Resize(10, 10));
  yuv.Fill(Ayuv38b{kAlphaMax, 0, 0, 0});

  EncoderConfig config;
  config.segments = 1;
  GlobalParams gparams;
  FeatureMap map;
  gparams.features_ = &map;
  ASSERT_WP2_OK(FindSegments(yuv, config, &gparams));
  EXPECT_EQ(1u, gparams.segments_.size());
}

TEST(TestSegment, TinyImage) {
  YUVPlane yuv;
  ASSERT_WP2_OK(yuv.Resize(5, 2));
  yuv.Fill(Ayuv38b{kAlphaMax, 0, 0, 0});

  EncoderConfig config;
  config.segments = 2;
  GlobalParams gparams;
  FeatureMap map;
  gparams.features_ = &map;
  ASSERT_WP2_OK(FindSegments(yuv, config, &gparams));
  // Number of segments was reduced to 1.
  EXPECT_EQ(1u, gparams.segments_.size());
}

TEST(TestSegment, Quantization) {
  Segment segment;
  segment.SetYUVBounds(-500, 500);
  segment.SetQuality(GetQualityHint(60.f), kQualityToQFactor(60.f),
                     kNeutralMultiplier / 2, kNeutralMultiplier * 2);
  const uint16_t lambda = 100;
  segment.FinalizeQuant(/*lambda=*/lambda);
  const uint32_t quantized_range1 = segment.quant_y_.GetMaxAbsResDC(TRF_16x16);
  const uint32_t u_quantized_range1 =
      segment.quant_u_.GetMaxAbsResDC(TRF_16x16);
  const uint32_t v_quantized_range1 =
      segment.quant_v_.GetMaxAbsResDC(TRF_16x16);
  segment.SetYUVBounds(-200, 200);
  segment.FinalizeQuant(lambda);
  const uint32_t quantized_range2 = segment.quant_y_.GetMaxAbsResDC(TRF_16x16);
  segment.SetYUVBounds(-1000, 1000);
  segment.FinalizeQuant(lambda);
  const uint32_t quantized_range3 = segment.quant_y_.GetMaxAbsResDC(TRF_16x16);
  // V is more quantized than U.
  EXPECT_LT(v_quantized_range1, u_quantized_range1);
  // Check the range of quantized values is more or less the same (modulo
  // rounding errors) regardless of the yuv bounds, as long as the yuv range
  // is centered.
  EXPECT_NEAR(quantized_range1, quantized_range2, 2);
  EXPECT_NEAR(quantized_range2, quantized_range3, 2);
  // Here the range of YUV values is not centered so the max quantized value
  // should be twice as high.
  segment.SetYUVBounds(0, 1000);
  segment.FinalizeQuant(lambda);
  const uint32_t not_centered_range =
      segment.quant_y_.GetMaxAbsResDC(TRF_16x16);
  EXPECT_NEAR(not_centered_range, quantized_range1 * 2, 1);

  segment.SetQuality(GetQualityHint(20.f), kQualityToQFactor(20.f),
                     kDefaultQuantMultiplier, kDefaultQuantMultiplier);
  segment.FinalizeQuant(lambda);
  const uint32_t lower_quality_factor_range =
      segment.quant_y_.GetMaxAbsResDC(TRF_16x16);
  // Counter intuitively, a lower "quality" setting actually means we quantize
  // less (better actual quality), and therefore the range gets larger.
  EXPECT_GT(lower_quality_factor_range, quantized_range3);
}

TEST(TestSegment, UVQuantMultiplier) {
  Segment segment;
  segment.SetYUVBounds(-500, 500);
  segment.SetQuality(GetQualityHint(50.f), kQualityToQFactor(50.f),
                     /*u_quant_multiplier=*/kNeutralMultiplier,
                     /*v_quant_multiplier=*/kNeutralMultiplier);
  const uint16_t lambda = 10;
  segment.FinalizeQuant(lambda);
  EXPECT_EQ(segment.quant_steps_[kUChannel][3],
            segment.quant_steps_[kVChannel][3]);

  segment.SetQuality(GetQualityHint(50.f), kQualityToQFactor(50.f),
                     /*u_quant_multiplier=*/1,
                     /*v_quant_multiplier=*/(1 << kQuantMultiplierBits));
  segment.FinalizeQuant(lambda);
  EXPECT_LE(segment.quant_steps_[kUChannel][3],
            segment.quant_steps_[kVChannel][3]);
}

TEST(TestSegment, IsMergeableWith) {
  Segment segment1;
  segment1.SetYUVBounds(-512, 512);
  segment1.SetQuality(GetQualityHint(42.f), kQualityToQFactor(42.f),
                      kDefaultQuantMultiplier, kDefaultQuantMultiplier);
  segment1.FinalizeQuant(100);

  Segment segment2;
  segment2.SetYUVBounds(-512, 512);
  segment2.SetQuality(GetQualityHint(42.f), kQualityToQFactor(42.f),
                      kDefaultQuantMultiplier, kDefaultQuantMultiplier);
  segment2.FinalizeQuant(/*lambda=*/100);

  // Segments are the same.
  EXPECT_TRUE(segment1.IsMergeableWith(segment2));
  EXPECT_TRUE(segment2.IsMergeableWith(segment1));

  // Change lambda, segments are still mergeable (we don't take lambda into
  // account right now).
  segment2.FinalizeQuant(/*lambda=*/200);
  EXPECT_TRUE(segment1.IsMergeableWith(segment2));

  // Set a different quality. Segments are now different.
  segment2.SetQuality(GetQualityHint(43.f), kQualityToQFactor(43.f),
                      kDefaultQuantMultiplier, kDefaultQuantMultiplier);
  EXPECT_FALSE(segment1.IsMergeableWith(segment2));
  EXPECT_FALSE(segment2.IsMergeableWith(segment1));

  // Revert quality and change the bounds. Segments are also different.
  segment2.SetQuality(GetQualityHint(42.f), kQualityToQFactor(42.f),
                      kDefaultQuantMultiplier, kDefaultQuantMultiplier);
  segment2.SetYUVBounds(-1024, 1024);
  EXPECT_FALSE(segment1.IsMergeableWith(segment2));
}

TEST(TestSegment, MinQuant) {
  Segment segment;
  // Largest possible input range.
  segment.SetYUVBounds(-(1 << kMaxYuvBits), (1 << kMaxYuvBits) - 1);
  // Lowest possible quantization.
  segment.SetQuality(kMaxLossyQualityHint, /*quality_factor=*/0,
                     /*u_quant_multiplier=*/0,
                     /*v_quant_multiplier=*/0);
  segment.FinalizeQuant(/*lambda=*/100);
  // Make sure the max dc value is STRICTLY lower than kMaxDcValue which means
  // it wasn't capped. If this test fails, it probably means kMaxDcValue is too
  // small.
  EXPECT_LT(segment.GetMaxAbsDC(kYChannel), kMaxDcValue);
  EXPECT_LT(segment.GetMaxAbsDC(kUChannel), kMaxDcValue);
  EXPECT_LT(segment.GetMaxAbsDC(kVChannel), kMaxDcValue);
  EXPECT_LT(segment.GetMaxAbsDC(kAChannel), kMaxDcValue);
}

}  // namespace
}  // namespace WP2
