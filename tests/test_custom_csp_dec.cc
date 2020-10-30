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

// Decoding to custom color space test.

#include <string>
#include <tuple>

#include "extras/ccsp_imageio.h"
#include "extras/extras.h"
#include "imageio/image_dec.h"
#include "include/helpers.h"
#include "src/dsp/math.h"
#include "src/utils/csp.h"
#include "src/utils/plane.h"
#include "src/wp2/decode.h"
#include "src/wp2/encode.h"

namespace WP2 {
namespace {

//------------------------------------------------------------------------------

// Matrices and associated shifts from extras.h
enum class Matrix { kRGB8, kRGB10, kYCoCg, kYpUV };
constexpr const int16_t* const kRGBToCCSPMatrices[] = {
    kRGBToRGB8Matrix, kRGBToRGB10Matrix, kRGBToYCoCgMatrix, kRGBToYpUVMatrix};
constexpr const uint32_t kRGBToCCSPShifts[] = {
    kRGBToRGB8Shift, kRGBToRGB10Shift, kRGBToYCoCgShift, kRGBToYpUVShift};

class DecodeCustomCspTest
    : public ::testing::TestWithParam<
          std::tuple<std::string, bool, float, Matrix, bool, bool, bool, bool,
                     bool, bool>> {
  void SetUp() override { WP2CSPConverterInit(); }
};

//------------------------------------------------------------------------------

TEST_P(DecodeCustomCspTest, Comparison) {
  const std::string& src_file_name = std::get<0>(GetParam());
  const bool keep_alpha = std::get<1>(GetParam());
  const float quality = std::get<2>(GetParam());
  const Matrix matrix = std::get<3>(GetParam());
  const int16_t* const rgb_to_ccsp_matrix = kRGBToCCSPMatrices[(int)matrix];
  const uint32_t rgb_to_ccsp_shift = kRGBToCCSPShifts[(int)matrix];
  const bool pad_y = std::get<4>(GetParam());
  const bool pad_u = std::get<5>(GetParam());
  const bool pad_v = std::get<6>(GetParam());
  const bool pad_a = std::get<7>(GetParam());
  const bool get_alpha = std::get<8>(GetParam());
  const bool get_metadata = std::get<9>(GetParam());

  ArgbBuffer src;
  ASSERT_WP2_OK(
      ReadImage(testing::GetTestDataPath(src_file_name).c_str(), &src));
  if (!keep_alpha) {
    ASSERT_WP2_OK(src.CompositeOver(Argb32b{255, 0, 0, 0}));
  }
  const bool is_opaque = !src.HasTransparency();

  EncoderConfig encoder_config = EncoderConfig::kDefault;
  encoder_config.quality = quality;
  MemoryWriter data;
  ASSERT_WP2_OK(Encode(src, &data, encoder_config));

  // Decode through the regular pipeline to have a reference output.
  ArgbBuffer ref;
  ASSERT_WP2_OK(Decode(data.mem_, data.size_, &ref));
  // No point in going further if the reference does not pass.
  ASSERT_TRUE(testing::Compare(src, ref, src_file_name,
                               testing::GetExpectedDistortion(quality)));

  YUVPlane custom_output;  // Might contain RGB.
  Metadata metadata;
  const uint32_t padded_width = Pad(src.width, kPredWidth);
  const uint32_t padded_height = Pad(src.height, kPredWidth);
  ASSERT_WP2_OK(custom_output.Y.Resize(pad_y ? padded_width : src.width,
                                       pad_y ? padded_height : src.height));
  ASSERT_WP2_OK(custom_output.U.Resize(pad_u ? padded_width : src.width,
                                       pad_u ? padded_height : src.height));
  ASSERT_WP2_OK(custom_output.V.Resize(pad_v ? padded_width : src.width,
                                       pad_v ? padded_height : src.height));
  if (get_alpha) {
    ASSERT_WP2_OK(custom_output.A.Resize(pad_a ? padded_width : src.width,
                                        pad_a ? padded_height : src.height));
  }

  ASSERT_WP2_OK(Decode(
      data.mem_, data.size_, rgb_to_ccsp_matrix, rgb_to_ccsp_shift,
      custom_output.Y.Row(0), custom_output.Y.Step(), custom_output.Y.Size(),
      custom_output.U.Row(0), custom_output.U.Step(), custom_output.U.Size(),
      custom_output.V.Row(0), custom_output.V.Step(), custom_output.V.Size(),
      custom_output.A.IsEmpty() ? nullptr : custom_output.A.Row(0),
      custom_output.A.Step(), custom_output.A.Size(),
      get_metadata ? &metadata : nullptr));

  if (get_metadata) {
    EXPECT_TRUE(testing::HasSameData(metadata.exif, ref.metadata.exif));
    EXPECT_TRUE(testing::HasSameData(metadata.xmp, ref.metadata.xmp));
    EXPECT_TRUE(testing::HasSameData(metadata.iccp, ref.metadata.iccp));
  }

  switch (matrix) {
    case Matrix::kRGB8:
    case Matrix::kRGB10: {
      // Try the YUV pipeline with a custom YUV which is actually RGB, for
      // easier comparison.
      const YUVPlane& rgb = custom_output;
      // Change the container but keep the data as is, except kRGB10 shift.
      ArgbBuffer output;
      ASSERT_WP2_OK(output.Resize(ref.width, ref.height));
      for (uint32_t y = 0; y < ref.height; ++y) {
        for (uint32_t x = 0; x < ref.width; ++x) {
          const uint32_t px = x * WP2FormatBpp(output.format);
          uint8_t* const pixel = &((uint8_t*)output.GetRow(y))[px];
          pixel[0] = rgb.A.IsEmpty() ? ((uint8_t*)ref.GetRow(y))[px + 0]
                                     : (uint8_t)rgb.A.At(x, y);
          constexpr int16_t min = 0, max = 255;
          if (matrix == Matrix::kRGB8) {
            pixel[1] = (uint8_t)Clamp(rgb.Y.At(x, y), min, max);
            pixel[2] = (uint8_t)Clamp(rgb.U.At(x, y), min, max);
            pixel[3] = (uint8_t)Clamp(rgb.V.At(x, y), min, max);
          } else {
            pixel[1] =
                (uint8_t)Clamp(RightShiftRound(rgb.Y.At(x, y), 2), min, max);
            pixel[2] =
                (uint8_t)Clamp(RightShiftRound(rgb.U.At(x, y), 2), min, max);
            pixel[3] =
                (uint8_t)Clamp(RightShiftRound(rgb.V.At(x, y), 2), min, max);
          }
        }
      }
      // Compare the source input with the custom RGB output.
      EXPECT_TRUE(testing::Compare(src, output, src_file_name,
                                   testing::GetExpectedDistortion(quality)));
      // Compare the reference output with the custom RGB output.
      // For kRGB10 there might be a slight difference due to rounding.
      // For transparent images, the lack of clamping by alpha of
      // CSPTransform::YUVToCustom() can lead to noticeable differences.
      const float expected_disto =
          is_opaque ? ((matrix == Matrix::kRGB10) ? 45.f : 99.f) : 35.f;
      EXPECT_TRUE(testing::Compare(ref, output, src_file_name, expected_disto));
      break;
    }
    case Matrix::kYCoCg:
    case Matrix::kYpUV: {
      // Actual luma/chroma output.
      CSPTransform csp_transform;  // Match Decode() behavior.
      constexpr int16_t kRgbAvg[3] = {0, 0, 0};
      if (matrix == Matrix::kYCoCg) {
        //  1  1 -1      with 10b fixed-point precision (x 1<<10)
        //  1  0  1      so that the inverse will be 14b and thus
        //  1 -1 -1      after the shift of 12b, 8+2b remain
        constexpr int16_t kYCoCgMatrix[] = {1024, 1024, -1024, 1024, 0,
                                            1024, 1024, -1024, -1024};
        ASSERT_TRUE(csp_transform.Init(kYCoCgMatrix, kRgbAvg));
      } else {
        //  1.0  0.00000  1.13983      with 9b fixed-point precision (x 1<<9)
        //  1.0 -0.39465 -0.58060      so that the inverse will be 15b and thus
        //  1.0  2.03211  0.00000      after the shift of 12b, 8+3b remain
        constexpr int16_t kYpUVMatrix[] = {512,  0,   584,  512, -202,
                                           -297, 512, 1040, 0};
        ASSERT_TRUE(csp_transform.Init(kYpUVMatrix, kRgbAvg));
      }
      YUVPlane ref_yuv;
      ASSERT_WP2_OK(ref_yuv.Import(ref, get_alpha, csp_transform,
                                   /*resize_if_needed=*/true));
      // SetView() plane by plane to avoid IsSubsampled() checks.
      YUVPlane non_padded_yuv;
      const Rectangle rect = {0, 0, ref.width, ref.height};
      ASSERT_WP2_OK(non_padded_yuv.Y.SetView(custom_output.Y, rect));
      ASSERT_WP2_OK(non_padded_yuv.U.SetView(custom_output.U, rect));
      ASSERT_WP2_OK(non_padded_yuv.V.SetView(custom_output.V, rect));
      if (!custom_output.A.IsEmpty()) {
        ASSERT_WP2_OK(non_padded_yuv.A.SetView(custom_output.A, rect));
      }
      const uint32_t bit_depth = csp_transform.GetYUVPrecisionBits() + 1;
      // Some error is expected for lossy qualities because YCoCg is no longer
      // "lossless" when it has been tampered with during compression (it can
      // output RGB values that are not exact multiples of 1<<kMtxShift).
      // For transparent images, the lack of clamping by alpha of
      // CSPTransform::YUVToCustom() can lead to noticeable differences. There
      // might also be some slight rounding discrepancies.
      const float expected_disto =
          (quality <= kMaxLossyQuality) ? (is_opaque ? 55.f : 45.f) : 85.f;
      EXPECT_TRUE(
          testing::Compare(ref_yuv, non_padded_yuv, bit_depth, src_file_name,
                           expected_disto)) << "quality : " << quality;
      break;
    }
  }
}

//------------------------------------------------------------------------------

INSTANTIATE_TEST_SUITE_P(
    DecodeCustomCspTestInstantiation, DecodeCustomCspTest,
    ::testing::Combine(::testing::Values("source1_1x1.png", "source1_64x48.png",
                                         "test_exif_xmp.webp"),
                       ::testing::Values(false) /* keep_alpha */,
                       ::testing::Values(5.f, 100.f) /* quality */,
                       ::testing::Values(Matrix::kRGB8, Matrix::kRGB10,
                                         Matrix::kYCoCg, Matrix::kYpUV),
                       ::testing::Values(false, true) /* get_alpha */,
                       ::testing::Values(false, true) /* pad_y */,
                       ::testing::Values(true) /* pad_u */,
                       ::testing::Values(true) /* pad_v */,
                       ::testing::Values(true) /* pad_a */,
                       ::testing::Values(false) /* get_metadata */));

INSTANTIATE_TEST_SUITE_P(
    DecodeCustomCspTestInstantiationPadding, DecodeCustomCspTest,
    ::testing::Combine(::testing::Values("source1_1x1.png"),
                       ::testing::Values(false, true) /* keep_alpha */,
                       ::testing::Values(5.f, 100.f) /* quality */,
                       ::testing::Values(Matrix::kRGB8),
                       ::testing::Values(false, true) /* get_alpha */,
                       ::testing::Values(false, true) /* pad_y */,
                       ::testing::Values(false, true) /* pad_u */,
                       ::testing::Values(false, true) /* pad_v */,
                       ::testing::Values(false, true) /* pad_a */,
                       ::testing::Values(false) /* get_metadata */));

INSTANTIATE_TEST_SUITE_P(
    DecodeCustomCspTestInstantiationMetadata, DecodeCustomCspTest,
    ::testing::Combine(::testing::Values("source1_1x1.png",
                                         "test_exif_xmp.webp"),
                       ::testing::Values(false, true) /* keep_alpha */,
                       ::testing::Values(50.f) /* quality */,
                       ::testing::Values(Matrix::kYCoCg),
                       ::testing::Values(true) /* get_alpha */,
                       ::testing::Values(false) /* pad_y */,
                       ::testing::Values(true) /* pad_u */,
                       ::testing::Values(false) /* pad_v */,
                       ::testing::Values(true) /* pad_a */,
                       ::testing::Values(true) /* get_metadata */));

//------------------------------------------------------------------------------

static const struct BufferParam {
  const char* const label;
  const char* const src_file_name;
  uint32_t r_step, r_buffer_size;
  uint32_t g_step, g_buffer_size;
  uint32_t b_step, b_buffer_size;
  uint32_t a_step, a_buffer_size;
  WP2Status expected_status;
} kBufferParams[] = {
    { "0-sized buffers are nullptr.",
      "source1_1x1.png", 0, 0, 0, 0, 0, 0, 0, 0, WP2_STATUS_NULL_PARAMETER },
    { "allowed steps.",
      "source1_1x1.png", 2, 4, 2, 4, 1, 2, 3, 6, WP2_STATUS_OK },
    { "Odd buffer sizes are allowed.",
      "source1_1x1.png", 2, 5, 2, 5, 2, 7, 4, 9, WP2_STATUS_OK },
    { "buffer size are ok (reference).",
      "source1_1x1.png", 8, 16, 8, 16, 10, 20, 8, 16, WP2_STATUS_OK },
    { "Smaller C0-buffer size than stride.", "source1_1x1.png",
      8, 15, 8, 16, 10, 20, 8, 16, WP2_STATUS_BAD_DIMENSION },
    { "Smaller C1-buffer size than stride.", "source1_1x1.png",
      8, 16, 8, 15, 10, 20, 8, 16, WP2_STATUS_BAD_DIMENSION },
    { "Smaller C2-buffer size than stride.", "source1_1x1.png",
      8, 16, 8, 16, 10, 19, 8, 16, WP2_STATUS_BAD_DIMENSION },
    { "Smaller A-buffer size than stride.", "source1_1x1.png",
      8, 16, 8, 16, 10, 20, 8, 15, WP2_STATUS_BAD_DIMENSION },
    { "Buffer is too small.",
      "source1_32x32.png", 32, 2048, 32, 2046, 32, 2048, 0, 0,
      WP2_STATUS_BAD_DIMENSION },
    { "Buffer is fine.",
      "source1_32x32.png", 32, 2048, 32, 2048, 32, 2048, 0, 0, WP2_STATUS_OK },
    { "Buffer is fine.",
      "source1_32x32.png", 33, 4000, 35, 3000, 34, 2887, 0, 0, WP2_STATUS_OK }
};

class StrideBufferSizeTest : public ::testing::TestWithParam<BufferParam> {};

TEST_P(StrideBufferSizeTest, Combination) {
  const BufferParam& p = GetParam();

  MemoryWriter data;
  ASSERT_WP2_OK(testing::CompressImage(p.src_file_name, &data));

  Data r, g, b, a;
  ASSERT_WP2_OK(r.Resize(p.r_buffer_size, /*keep_bytes=*/false));
  ASSERT_WP2_OK(g.Resize(p.g_buffer_size, /*keep_bytes=*/false));
  ASSERT_WP2_OK(b.Resize(p.b_buffer_size, /*keep_bytes=*/false));
  ASSERT_WP2_OK(a.Resize(p.a_buffer_size, /*keep_bytes=*/false));
  int16_t* const r_buffer = reinterpret_cast<int16_t*>(r.bytes);
  int16_t* const g_buffer = reinterpret_cast<int16_t*>(g.bytes);
  int16_t* const b_buffer = reinterpret_cast<int16_t*>(b.bytes);
  int16_t* const a_buffer = reinterpret_cast<int16_t*>(a.bytes);

  ASSERT_EQ(Decode(data.mem_, data.size_, kRGBToRGB8Matrix, kRGBToRGB8Shift,
                   r_buffer, p.r_step, p.r_buffer_size,
                   g_buffer, p.g_step, p.g_buffer_size,
                   b_buffer, p.b_step, p.b_buffer_size,
                   a_buffer, p.a_step, p.a_buffer_size),
            p.expected_status) << p.label;
}

INSTANTIATE_TEST_SUITE_P(
    StrideBufferSizeTestInstantiation, StrideBufferSizeTest,
    ::testing::ValuesIn(kBufferParams));

//------------------------------------------------------------------------------

class YUVPlaneTest : public ::testing::TestWithParam<std::string> {};

TEST_P(YUVPlaneTest, LosslessUpsampling) {
  const std::string& file_name = GetParam();
  const std::string file_path = testing::GetTestDataPath(file_name);

  YUVPlane original;
  CSPMtx ccsp_to_rgb = {};
  ASSERT_WP2_OK(ReadImage(file_path.c_str(), &original, &ccsp_to_rgb));
  // Ignoring the 'ccsp_to_rgb_matrix', the color space is not tested here.
  ASSERT_WP2_OK(original.Downsample());

  uint32_t bit_depth;
  ASSERT_WP2_OK(ReadBitDepth(file_path.c_str(), &bit_depth));

  YUVPlane upsampled;
  // Upsampling then downsampling with no interpolation is lossless.
  ASSERT_WP2_OK(upsampled.UpsampleFrom(original, SamplingTaps::kUpNearest));
  ASSERT_WP2_OK(upsampled.Downsample(SamplingTaps::kDownAvg));
  EXPECT_TRUE(testing::Compare(original, upsampled, bit_depth, file_name));

  // Lossy upsampling and downsampling.
  ASSERT_WP2_OK(upsampled.UpsampleFrom(original, SamplingTaps::kUpSmooth));
  ASSERT_WP2_OK(upsampled.Downsample(SamplingTaps::kDownSharp));
  EXPECT_TRUE(
      testing::Compare(original, upsampled, bit_depth, file_name, 30.f));
}

TEST_P(YUVPlaneTest, NotTooLossyDownsampling) {
  const std::string& file_name = GetParam();
  const std::string file_path = testing::GetTestDataPath(file_name);

  YUVPlane original;
  CSPMtx ccsp_to_rgb = {};
  ASSERT_WP2_OK(ReadImage(file_path.c_str(), &original, &ccsp_to_rgb));
  // Ignoring 'ccsp_to_rgb' because the color space is not tested here.

  uint32_t bit_depth;
  ASSERT_WP2_OK(ReadBitDepth(file_path.c_str(), &bit_depth));

  YUVPlane downsampled;
  // ouch
  ASSERT_WP2_OK(downsampled.DownsampleFrom(original, SamplingTaps::kDownAvg));
  ASSERT_WP2_OK(downsampled.Upsample(SamplingTaps::kUpNearest));
  EXPECT_TRUE(
      testing::Compare(original, downsampled, bit_depth, file_name, 28.f));

  // save the furniture
  ASSERT_WP2_OK(downsampled.DownsampleFrom(original, SamplingTaps::kDownSharp));
  ASSERT_WP2_OK(downsampled.Upsample(SamplingTaps::kUpSmooth));
  EXPECT_TRUE(
      testing::Compare(original, downsampled, bit_depth, file_name, 29.f));

  float disto[5];
  ASSERT_WP2_OK(downsampled.GetDistortion(original, bit_depth, PSNR, disto));
  EXPECT_EQ(disto[0], 99.f);  // Luma and alpha are not impacted by
  EXPECT_EQ(disto[1], 99.f);  // downsampling, only chroma.
}

INSTANTIATE_TEST_SUITE_P(YUVPlaneTestInstantiation, YUVPlaneTest,
                         ::testing::Values("source0.ppm", "source1_64x48.png",
                                           "source1_1x48.png",
                                           "source1_64x1.png",
                                           "ccsp/source3_C444p12.y4m"));

//------------------------------------------------------------------------------

}  // namespace
}  // namespace WP2
