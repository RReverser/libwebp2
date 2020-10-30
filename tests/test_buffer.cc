// Copyright 2019 Google LLC
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
// Test for ArgbBuffer
//
// Author: Skal (pascal.massimino@gmail.com)

#include <array>
#include <string>
#include <tuple>
#include <vector>

#include "examples/example_utils.h"
#include "include/helpers.h"
#include "src/wp2/decode.h"

namespace WP2 {
namespace {

//------------------------------------------------------------------------------
// Test that some functions clear metadata, some don't.

TEST(BufferTest, Metadata) {
  std::array<uint8_t, 128> metadata{0};
  std::array<uint8_t, 64 * 4> pixels8{0};
  std::array<uint16_t, 64 * 4> pixels16{0};
  const uint8_t color8[] = {255, 1, 2, 3};
  const uint16_t color16[] = {255, 1, 2, 3};
  for (uint32_t format_ind = 0; format_ind < (uint32_t)WP2_FORMAT_NUM;
       ++format_ind) {
    const WP2SampleFormat format = (WP2SampleFormat)format_ind;
    ArgbBuffer buffer(format);

    ASSERT_WP2_OK(
        buffer.metadata.exif.CopyFrom(metadata.data(), metadata.size()));
    ASSERT_FALSE(buffer.metadata.IsEmpty());
    // SetExternal() should clear the metadata.
    if (WP2FormatBpc(format) == 1) {
      ASSERT_WP2_OK(buffer.SetExternal(8, 8, pixels8.data(), 8 * 4));
    } else {
      ASSERT_WP2_OK(
          buffer.SetExternal(8, 8, pixels16.data(), 8 * 4 * sizeof(uint16_t)));
    }
    ASSERT_TRUE(buffer.metadata.IsEmpty());

    ASSERT_WP2_OK(
        buffer.metadata.iccp.CopyFrom(metadata.data(), metadata.size()));
    ASSERT_FALSE(buffer.metadata.IsEmpty());
    ASSERT_WP2_OK(buffer.SetView(buffer, {1, 1, 1, 1}));
    ASSERT_TRUE(buffer.IsView());
    ASSERT_FALSE(
        buffer.metadata.IsEmpty());  // Kept metadata because same buffer

    ArgbBuffer other_buffer(format);
    ASSERT_WP2_OK(other_buffer.Resize(2, 2));
    ASSERT_WP2_OK(
        other_buffer.metadata.iccp.CopyFrom(metadata.data(), metadata.size()));
    ASSERT_WP2_OK(buffer.SetView(other_buffer, {1, 1, 1, 1}));
    ASSERT_EQ(other_buffer.SetView(buffer), WP2_STATUS_INVALID_PARAMETER);
    ASSERT_TRUE(buffer.metadata.IsEmpty());

    // Resize() should clear the buffer anew, including metadata.
    ASSERT_WP2_OK(
        buffer.metadata.xmp.CopyFrom(metadata.data(), metadata.size()));
    ASSERT_FALSE(buffer.metadata.IsEmpty());
    ASSERT_EQ(buffer.Resize(8, 8),
              WP2_STATUS_BAD_DIMENSION);  // Not exactly same
    ASSERT_WP2_OK(buffer.Resize(buffer.width, buffer.height));
    ASSERT_TRUE(buffer.IsView());
    ASSERT_TRUE(buffer.metadata.IsEmpty());

    ASSERT_WP2_OK(
        buffer.metadata.exif.CopyFrom(metadata.data(), metadata.size()));
    ASSERT_FALSE(buffer.metadata.IsEmpty());
    buffer.Deallocate();
    ASSERT_TRUE(buffer.metadata.IsEmpty());

    ASSERT_WP2_OK(
        buffer.metadata.exif.CopyFrom(metadata.data(), metadata.size()));
    ASSERT_FALSE(buffer.metadata.IsEmpty());
    ASSERT_WP2_OK(buffer.Resize(16, 16));
    ASSERT_TRUE(buffer.metadata.IsEmpty());

    if (WP2FormatBpc(format) == 1) {
      other_buffer.Fill({0, 0, 2, 2}, color8);
    } else {
      other_buffer.Fill({0, 0, 2, 2}, color16);
    }
    ASSERT_WP2_OK(
        buffer.metadata.exif.CopyFrom(metadata.data(), metadata.size()));
    ASSERT_FALSE(buffer.metadata.IsEmpty());
    ASSERT_WP2_OK(buffer.CopyFrom(other_buffer));
    ASSERT_TRUE(buffer.metadata.IsEmpty());  // Clears and doesn't copy metadata

    ArgbBuffer other_format(
        (WP2SampleFormat)((format_ind == 0) ? 1 : format_ind - 1));
    ASSERT_WP2_OK(other_format.Resize(2, 2));
    ASSERT_WP2_OK(
        other_format.metadata.iccp.CopyFrom(metadata.data(), metadata.size()));
    ASSERT_WP2_OK(
        buffer.metadata.iccp.CopyFrom(metadata.data(), metadata.size()));
    ASSERT_FALSE(other_format.metadata.IsEmpty());
    if (buffer.format == WP2_Argb_32 &&
        WP2FormatBpc(other_format.format) == 1) {
      ASSERT_WP2_OK(other_format.ConvertFrom(buffer));
      ASSERT_TRUE(other_format.metadata.IsEmpty());  // Doesn't copy metadata
    }

    ASSERT_FALSE(buffer.metadata.IsEmpty());
    buffer.Deallocate();
    ASSERT_TRUE(buffer.metadata.IsEmpty());
    ASSERT_FALSE(other_buffer.metadata.IsEmpty());
    ASSERT_WP2_OK(buffer.Swap(&other_buffer));
    ASSERT_FALSE(buffer.metadata.IsEmpty());
    ASSERT_TRUE(other_buffer.metadata.IsEmpty());
  }
}

//------------------------------------------------------------------------------

class BufferTest
    : public ::testing::TestWithParam<std::tuple<std::string, bool, bool>> {
};

WP2Status Decode(const uint8_t* data, size_t data_size,
                 ArgbBuffer* output, bool use_decoder) {
  if (use_decoder) {
    ArrayDecoder decoder(data, data_size, DecoderConfig::kDefault, output);
    const bool decoded_frame = decoder.ReadFrame();
    return decoder.Failed() ? decoder.GetStatus()
           : decoded_frame  ? WP2_STATUS_OK
                            : WP2_STATUS_NOT_ENOUGH_DATA;
  }
  return Decode(data, data_size, output);
}

TEST_P(BufferTest, API) {
  const std::string filename = std::get<0>(GetParam());
  const float quality = std::get<1>(GetParam()) ? 100.f : 75.f;
  const int speed = 2;
  const bool use_decoder = std::get<2>(GetParam());
  ArgbBuffer src;
  MemoryWriter data;
  ASSERT_WP2_OK(testing::CompressImage(filename, &data, &src, quality, speed));
  EXPECT_FALSE(src.IsView());
  ASSERT_WP2_OK(src.SetView(src, {1, 1, src.width / 2, src.height / 2}));
  EXPECT_FALSE(src.IsView());  // We still own the memory.

  ArgbBuffer output1;
  ASSERT_WP2_OK(Decode(data.mem_, data.size_, &output1, use_decoder));

  ArgbBuffer output2;
  // Try a too small buffer
  ASSERT_WP2_OK(
      output2.SetView(output1, {0, 0, output1.width / 2, output1.height / 2}));
  ASSERT_EQ(Decode(data.mem_, data.size_, &output2, use_decoder),
            WP2_STATUS_BAD_DIMENSION);
  EXPECT_FALSE(testing::Compare(output1, output2, filename));
  // Try a tightly fit one
  ASSERT_WP2_OK(output2.SetView(output1));
  EXPECT_TRUE(output2.IsView());
  ASSERT_WP2_OK(Decode(data.mem_, data.size_, &output2, use_decoder));
  EXPECT_TRUE(output2.IsView());  // should decode directly in the view
  EXPECT_TRUE(testing::Compare(output1, output2, filename));

  // And then a large internal one
  output2.Deallocate();
  ASSERT_WP2_OK(output2.Resize(output1.width * 2, output1.height * 3));
  ASSERT_WP2_OK(Decode(data.mem_, data.size_, &output2, use_decoder));
  EXPECT_TRUE(testing::Compare(output1, output2, filename));

  // And a too large external one
  ASSERT_WP2_OK(output1.Resize(output1.width * 2, output1.height * 3));
  ASSERT_WP2_OK(output2.SetView(output1));
  ASSERT_EQ(Decode(data.mem_, data.size_, &output2, use_decoder),
            WP2_STATUS_BAD_DIMENSION);
}

INSTANTIATE_TEST_SUITE_P(
    BufferTestInstantiation1, BufferTest,
    ::testing::Combine(
        ::testing::Values("source1.png"),
        ::testing::Values(true, false) /* lossless */,
        ::testing::Values(true, false) /* use_decoder */));

//------------------------------------------------------------------------------

TEST(BufferTest, Composite) {
  ArgbBuffer buffer(WP2_Argb_32);
  ASSERT_WP2_OK(buffer.Resize(4, 1));
  buffer.Fill({0, 0, 1, 1}, Argb32b({255, 0, 30, 60}));
  buffer.Fill({1, 0, 1, 1}, Argb32b({200, 120, 30, 60}));
  buffer.Fill({2, 0, 1, 1}, Argb32b({30, 20, 25, 0}));
  buffer.Fill({3, 0, 1, 1}, Argb32b({0, 0, 0, 0}));

  EXPECT_TRUE(buffer.HasTransparency());

  ASSERT_WP2_OK(buffer.CompositeOver({127, 13, 57, 116}));

  uint8_t* row = (uint8_t*)buffer.GetRow(0);

  EXPECT_TRUE(buffer.HasTransparency());

  EXPECT_EQ(row[4 * 0 + 0], 255);
  EXPECT_EQ(row[4 * 0 + 1], 0);
  EXPECT_EQ(row[4 * 0 + 2], 30);
  EXPECT_EQ(row[4 * 0 + 3], 60);

  EXPECT_EQ(row[4 * 1 + 0], 227);
  EXPECT_EQ(row[4 * 1 + 1], 123);
  EXPECT_EQ(row[4 * 1 + 2], 42);
  EXPECT_EQ(row[4 * 1 + 3], 85);

  EXPECT_EQ(row[4 * 2 + 0], 142);
  EXPECT_EQ(row[4 * 2 + 1], 31);
  EXPECT_EQ(row[4 * 2 + 2], 75);
  EXPECT_EQ(row[4 * 2 + 3], 102);

  EXPECT_EQ(row[4 * 3 + 0], 127);
  EXPECT_EQ(row[4 * 3 + 1], 13);
  EXPECT_EQ(row[4 * 3 + 2], 57);
  EXPECT_EQ(row[4 * 3 + 3], 116);

  ASSERT_WP2_OK(buffer.CompositeOver({255, 0, 0, 0}));

  EXPECT_FALSE(buffer.HasTransparency());

  EXPECT_EQ(row[4 * 0 + 0], 255);
  EXPECT_EQ(row[4 * 0 + 1], 0);
  EXPECT_EQ(row[4 * 0 + 2], 30);
  EXPECT_EQ(row[4 * 0 + 3], 60);

  EXPECT_EQ(row[4 * 1 + 0], 255);
  EXPECT_EQ(row[4 * 1 + 1], 123);
  EXPECT_EQ(row[4 * 1 + 2], 42);
  EXPECT_EQ(row[4 * 1 + 3], 85);

  EXPECT_EQ(row[4 * 2 + 0], 255);
  EXPECT_EQ(row[4 * 2 + 1], 31);
  EXPECT_EQ(row[4 * 2 + 2], 75);
  EXPECT_EQ(row[4 * 2 + 3], 102);

  EXPECT_EQ(row[4 * 3 + 0], 255);
  EXPECT_EQ(row[4 * 3 + 1], 13);
  EXPECT_EQ(row[4 * 3 + 2], 57);
  EXPECT_EQ(row[4 * 3 + 3], 116);
}

TEST(BufferTest, CompositeOnBuffer) {
  ArgbBuffer buffer(WP2_Argb_32);
  ASSERT_WP2_OK(buffer.Resize(4, 2));
  buffer.Fill({0, 0, 1, 2}, Argb32b({255, 0, 30, 60}));
  buffer.Fill({1, 0, 1, 2}, Argb32b({200, 120, 30, 60}));
  buffer.Fill({2, 0, 1, 2}, Argb32b({30, 20, 25, 0}));
  buffer.Fill({3, 0, 1, 2}, Argb32b({0, 0, 0, 0}));
  ArgbBuffer buffer2(WP2_Argb_32);
  ASSERT_WP2_OK(buffer2.CopyFrom(buffer));

  ArgbBuffer bg_buffer(WP2_Argb_32);
  ASSERT_WP2_OK(bg_buffer.Resize(4, 2));
  bg_buffer.Fill({0, 0, 4, 1}, Argb32b({127, 13, 57, 116}));
  bg_buffer.Fill({0, 1, 4, 1}, Argb32b({255, 0, 0, 0}));
  ArgbBuffer bg_buffer2(WP2_Argb_32);
  ASSERT_WP2_OK(bg_buffer2.CopyFrom(bg_buffer));

  ASSERT_WP2_OK(buffer.CompositeOver(bg_buffer));
  ASSERT_WP2_OK(bg_buffer2.CompositeUnder(buffer2));

  ASSERT_TRUE(testing::Compare(buffer, bg_buffer2, "buffer"));

  uint8_t* row = (uint8_t*)buffer.GetRow(0);

  EXPECT_EQ(row[4 * 0 + 0], 255);
  EXPECT_EQ(row[4 * 0 + 1], 0);
  EXPECT_EQ(row[4 * 0 + 2], 30);
  EXPECT_EQ(row[4 * 0 + 3], 60);

  EXPECT_EQ(row[4 * 1 + 0], 227);
  EXPECT_EQ(row[4 * 1 + 1], 123);
  EXPECT_EQ(row[4 * 1 + 2], 42);
  EXPECT_EQ(row[4 * 1 + 3], 85);

  EXPECT_EQ(row[4 * 2 + 0], 142);
  EXPECT_EQ(row[4 * 2 + 1], 31);
  EXPECT_EQ(row[4 * 2 + 2], 75);
  EXPECT_EQ(row[4 * 2 + 3], 102);

  EXPECT_EQ(row[4 * 3 + 0], 127);
  EXPECT_EQ(row[4 * 3 + 1], 13);
  EXPECT_EQ(row[4 * 3 + 2], 57);
  EXPECT_EQ(row[4 * 3 + 3], 116);

  row = (uint8_t*)buffer.GetRow(1);

  EXPECT_EQ(row[4 * 0 + 0], 255);
  EXPECT_EQ(row[4 * 0 + 1], 0);
  EXPECT_EQ(row[4 * 0 + 2], 30);
  EXPECT_EQ(row[4 * 0 + 3], 60);

  EXPECT_EQ(row[4 * 1 + 0], 255);
  EXPECT_EQ(row[4 * 1 + 1], 120);
  EXPECT_EQ(row[4 * 1 + 2], 30);
  EXPECT_EQ(row[4 * 1 + 3], 60);

  EXPECT_EQ(row[4 * 2 + 0], 255);
  EXPECT_EQ(row[4 * 2 + 1], 20);
  EXPECT_EQ(row[4 * 2 + 2], 25);
  EXPECT_EQ(row[4 * 2 + 3], 0);

  EXPECT_EQ(row[4 * 3 + 0], 255);
  EXPECT_EQ(row[4 * 3 + 1], 0);
  EXPECT_EQ(row[4 * 3 + 2], 0);
  EXPECT_EQ(row[4 * 3 + 3], 0);
}

void FillWithNoise(size_t seed, ArgbBuffer* const buffer) {
  assert(buffer->format == WP2_Argb_32);
  const uint32_t bpp = WP2FormatBpp(buffer->format);
  UniformIntDistribution random(seed);  // Deterministic.
  for (uint32_t y = 0; y < buffer->height; ++y) {
    uint8_t* const row = (uint8_t*)buffer->GetRow(y);
    for (uint32_t x = 0; x < buffer->width; ++x) {
      row[x * bpp] = random.Get<uint8_t>(0, 255);  // Alpha
      for (uint32_t c = 1; c < bpp; ++c) {
        // Premultiplied
        row[x * bpp + c] = random.Get<uint8_t>(0, row[x * bpp]);
      }
    }
  }
}

TEST(BufferTest, CompositeYuv) {
  ArgbBuffer foreground(WP2_Argb_32);
  const uint32_t width = 100;
  const uint32_t height = 50;
  ASSERT_WP2_OK(foreground.Resize(width, height));
  FillWithNoise(/*seed=*/1, &foreground);

  ArgbBuffer background(WP2_Argb_32);
  ASSERT_WP2_OK(background.Resize(width, height));
  FillWithNoise(/*seed=*/2, &background);

  CSPTransform csp;
  static constexpr int16_t kYuvToRgb[] = {
      1188, 0, 1638, 1188, -399, -829, 1188, 2058, 0,
  };
  static constexpr int16_t kRgbAvg[] = {-100, 50, 110};
  ASSERT_TRUE(csp.Init(kYuvToRgb, kRgbAvg));
  YUVPlane yuv_foreground;
  ASSERT_WP2_OK(yuv_foreground.Import(foreground, /*import_alpha=*/true, csp,
                                      /*resize_if_needed=*/true));
  YUVPlane yuv_background;
  ASSERT_WP2_OK(yuv_background.Import(background, /*import_alpha=*/true, csp,
                                      /*resize_if_needed=*/true));

  // Composite yuv foreground on yuv background.
  YUVPlane yuv_composite_over;
  ASSERT_WP2_OK(
      yuv_composite_over.Copy(yuv_foreground, /*resize_if_needed=*/true));
  ASSERT_WP2_OK(yuv_composite_over.CompositeOver(yuv_background, csp));

  // Composite yuv background under yuv foreground.
  YUVPlane yuv_composite_under;
  ASSERT_WP2_OK(
      yuv_composite_under.Copy(yuv_background, /*resize_if_needed=*/true));
  ASSERT_WP2_OK(yuv_composite_under.CompositeUnder(yuv_foreground, csp));

  // Composite in rgb then convert to yuv.
  ASSERT_WP2_OK(foreground.CompositeOver(background));
  YUVPlane yuv_imported_frop_rgb;
  ASSERT_WP2_OK(yuv_imported_frop_rgb.Import(foreground, /*import_alpha=*/true,
                                             csp,
                                             /*resize_if_needed=*/true));

  // All three methods should give the same result.
  // CompositeOVer/CompositeUnder.
  EXPECT_TRUE(testing::Compare(yuv_composite_over, yuv_composite_under, 16,
                               "test", 99));
  // Composite in RGB vs YUV. We allow for some rounding error.
  for (uint32_t y = 0; y < height; ++y) {
    for (uint32_t x = 0; x < width; ++x) {
      SCOPED_TRACE(SPrintf("x %d y %d", x, y));
      EXPECT_NEAR(yuv_composite_over.Y.At(x, y),
                  yuv_imported_frop_rgb.Y.At(x, y), 3);
      EXPECT_NEAR(yuv_composite_over.U.At(x, y),
                  yuv_imported_frop_rgb.U.At(x, y), 3);
      EXPECT_NEAR(yuv_composite_over.V.At(x, y),
                  yuv_imported_frop_rgb.V.At(x, y), 3);
      EXPECT_EQ(yuv_composite_over.A.At(x, y),
                yuv_imported_frop_rgb.A.At(x, y));
    }
  }
}

//------------------------------------------------------------------------------

TEST(BufferTest, GetDistortion_BlackSquare) {
  // Transparent background.
  const uint8_t transparent_8[] = {0, 0, 0, 0};
  const uint16_t transparent_16[] = {0, 0, 0, 0};
  // Black square.
  const uint8_t black_8[] = {255, 0, 0, 0};
  const uint16_t black_16[] = {255, 0, 0, 0};
  // Transparent black square.
  const uint8_t black_transp_8[] = {250, 0, 0, 0};
  const uint16_t black_transp_16[] = {250, 0, 0, 0};
  for (WP2SampleFormat format : {WP2_Argb_32, WP2_Argb_38}) {
    ArgbBuffer orig(format);
    ASSERT_WP2_OK(orig.Resize(100, 100));
    // Opaque black square on a transparent background.
    if (WP2FormatBpc(format) == 1) {
      orig.Fill(/*window=*/{0, 0, 100, 100}, transparent_8);
      orig.Fill(/*window=*/{25, 25, 50, 50}, black_8);
    } else {
      orig.Fill(/*window=*/{0, 0, 100, 100}, transparent_16);
      orig.Fill(/*window=*/{25, 25, 50, 50}, black_16);
    }

    ArgbBuffer compressed(format);
    ASSERT_WP2_OK(compressed.Resize(100, 100));
    // Slightly transparent black square on a transparent background.
    if (WP2FormatBpc(format) == 1) {
      compressed.Fill(/*window=*/{0, 0, 100, 100}, transparent_8);
      compressed.Fill(/*window=*/{25, 25, 50, 50}, black_transp_8);
    } else {
      compressed.Fill(/*window=*/{0, 0, 100, 100}, transparent_16);
      compressed.Fill(/*window=*/{25, 25, 50, 50}, black_transp_16);
    }

    for (uint32_t metric = 0; metric < NUM_METRIC_TYPES; ++metric) {
      // Not yet implemented.
      if (format == WP2_Argb_38 && (MetricType)metric != PSNR) continue;
      SCOPED_TRACE(SPrintf("metric %s", kWP2MetricNames[metric]));
      float result[5];
      ASSERT_WP2_OK(compressed.GetDistortion(orig, (MetricType)metric, result));
      EXPECT_LT(result[0], 50.f);  // a
      EXPECT_EQ(result[1], 99.f);  // r (or y for PSNRHVS)
      EXPECT_EQ(result[2], 99.f);  // g (or u for PSNRHVS)
      EXPECT_EQ(result[3], 99.f);  // b (or v for PSNRHVS)
      if (metric == MetricType::PSNR_YUV) {
        EXPECT_LT(result[4], 60.f);  // all
      } else {
        EXPECT_LT(result[4], 50.f);  // all
      }
      // Not yet implemented.
      if (format == WP2_Argb_38) continue;
      ASSERT_WP2_OK(compressed.GetDistortionBlackOrWhiteBackground(
          orig, (MetricType)metric, result));

      EXPECT_LT(result[0], 50.f);  // a
      if (metric == MetricType::PSNRHVS || metric == MetricType::PSNR_YUV ||
          metric == MetricType::SSIM_YUV) {
        if (metric == MetricType::PSNR_YUV) {
          EXPECT_EQ(result[1], 99.f);  // y
        } else {
          EXPECT_LT(result[1], 40.f);  // y
        }
        EXPECT_EQ(result[2], 99.f);  // u
        EXPECT_EQ(result[3], 99.f);  // v
      } else {
        EXPECT_LT(result[1], 50.f);  // r
        EXPECT_LT(result[2], 50.f);  // g
        EXPECT_LT(result[3], 50.f);  // b
      }
      EXPECT_LT(result[4], 50.f);  // all
    }
  }
}

TEST(BufferTest, GetDistortion_WhiteSquare) {
  ArgbBuffer orig(WP2_Argb_32);
  ASSERT_WP2_OK(orig.Resize(100, 100));
  // Opaque white square on a transparent background.
  orig.Fill({0, 0, 100, 100}, Argb32b({0, 0, 0, 0}));
  orig.Fill({25, 25, 50, 50}, Argb32b({255, 255, 255, 255}));

  ArgbBuffer compressed(WP2_Argb_32);
  ASSERT_WP2_OK(compressed.Resize(100, 100));
  // Slightly transparent white square on a transparent background.
  compressed.Fill({0, 0, 100, 100}, Argb32b({0, 0, 0, 0}));
  compressed.Fill({25, 25, 50, 50}, Argb32b({250, 250, 250, 250}));

  for (uint32_t metric = 0; metric < NUM_METRIC_TYPES; ++metric) {
    SCOPED_TRACE(SPrintf("metric %s", kWP2MetricNames[metric]));
    float result[5];
    ASSERT_WP2_OK(compressed.GetDistortion(orig, (MetricType)metric, result));

    EXPECT_LT(result[0], 50.f);  // a
    if (metric == MetricType::PSNRHVS || metric == MetricType::PSNR_YUV ||
        metric == MetricType::SSIM_YUV) {
      EXPECT_LT(result[1], 50.f);  // y
      EXPECT_EQ(result[2], 99.f);  // u
      EXPECT_EQ(result[3], 99.f);  // v
    } else {
      EXPECT_LT(result[1], 50.f);  // r
      EXPECT_LT(result[2], 50.f);  // g
      EXPECT_LT(result[3], 50.f);  // b
    }
    EXPECT_LT(result[4], 50.f);  // all

    ASSERT_WP2_OK(compressed.GetDistortionBlackOrWhiteBackground(
        orig, (MetricType)metric, result));

    EXPECT_LT(result[0], 50.f);  // a
    if (metric == MetricType::PSNRHVS || metric == MetricType::PSNR_YUV ||
        metric == MetricType::SSIM_YUV) {
      EXPECT_LT(result[1], 50.f);  // y
      EXPECT_EQ(result[2], 99.f);  // u
      EXPECT_EQ(result[3], 99.f);  // v
    } else {
      EXPECT_LT(result[1], 50.f);  // r
      EXPECT_LT(result[2], 50.f);  // g
      EXPECT_LT(result[3], 50.f);  // b
    }
    EXPECT_LT(result[4], 50.f);  // all
  }
}

//------------------------------------------------------------------------------

TEST(RectangleTest, Contains) {
  Rectangle r1{5, 0, 10, 12};
  EXPECT_TRUE(r1.Contains(5, 0));
  EXPECT_TRUE(r1.Contains(6, 0));
  EXPECT_TRUE(r1.Contains(8, 11));
  EXPECT_TRUE(r1.Contains(14, 11));
  EXPECT_FALSE(r1.Contains(4, 0));
  EXPECT_FALSE(r1.Contains(15, 2));
  EXPECT_FALSE(r1.Contains(7, 12));
  EXPECT_FALSE(r1.Contains(15, 12));

  EXPECT_TRUE(r1.Contains(r1));
  EXPECT_TRUE(r1.Contains({5, 0, 0, 0}));
  EXPECT_TRUE(r1.Contains({6, 2, 3, 7}));
  EXPECT_TRUE(r1.Contains({6, 2, 9, 7}));
  EXPECT_FALSE(r1.Contains({6, 2, 10, 7}));
}

//------------------------------------------------------------------------------

}  // namespace
}  // namespace WP2
