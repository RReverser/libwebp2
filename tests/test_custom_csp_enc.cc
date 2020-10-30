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

// Encoding from custom color space test.

#include <string>
#include <tuple>

#include "extras/ccsp_imageio.h"
#include "extras/extras.h"
#include "imageio/image_dec.h"
#include "include/helpers.h"
#include "src/dsp/math.h"
#include "src/utils/plane.h"
#include "src/utils/vector.h"
#include "src/wp2/decode.h"
#include "src/wp2/encode.h"

namespace WP2 {
namespace {

//------------------------------------------------------------------------------

// Matrices and associated shifts from extras.h
enum class Matrix { kRGB8, kRGB10, kYCoCg, kYpUV };
constexpr const int16_t* const kCCSPToRGBMatrices[] = {
    kRGB8ToRGBMatrix, kRGB10ToRGBMatrix, kYCoCgToRGBMatrix, kYpUVToRGBMatrix};
constexpr const uint32_t kCCSPToRGBShifts[] = {
    kRGB8ToRGBShift, kRGB10ToRGBShift, kYCoCgToRGBShift, kYpUVToRGBShift};

class EncodeCustomCspTest
    : public ::testing::TestWithParam<std::tuple<
          std::string, float, Matrix, bool, bool, bool, bool, bool, bool>> {
  void SetUp() override { WP2CSPConverterInit(); }
};

//------------------------------------------------------------------------------

TEST_P(EncodeCustomCspTest, Comparison) {
  const std::string& src_file_name = std::get<0>(GetParam());
  const float quality = std::get<1>(GetParam());
  const Matrix matrix = std::get<2>(GetParam());
  const int16_t* const ccsp_to_rgb_matrix = kCCSPToRGBMatrices[(int)matrix];
  const uint32_t ccsp_to_rgb_shift = kCCSPToRGBShifts[(int)matrix];
  const bool pad_y = std::get<3>(GetParam());
  const bool pad_u = std::get<4>(GetParam());
  const bool pad_v = std::get<5>(GetParam());
  const bool pad_a = std::get<6>(GetParam());
  const bool encode_alpha = std::get<7>(GetParam());
  const bool encode_metadata = std::get<8>(GetParam());

  ArgbBuffer src;
  MemoryWriter ref_data;
  ASSERT_WP2_OK(
      testing::CompressImage(src_file_name, &ref_data, &src, quality));

  YUVPlane custom_input;  // Might contain RGB.
  const uint32_t padded_width = Pad(src.width, kPredWidth);
  const uint32_t padded_height = Pad(src.height, kPredWidth);
  ASSERT_WP2_OK(custom_input.Y.Resize(pad_y ? padded_width : src.width,
                                      pad_y ? padded_height : src.height));
  ASSERT_WP2_OK(custom_input.U.Resize(pad_u ? padded_width : src.width,
                                      pad_u ? padded_height : src.height));
  ASSERT_WP2_OK(custom_input.V.Resize(pad_v ? padded_width : src.width,
                                      pad_v ? padded_height : src.height));
  if (encode_alpha) {
    ASSERT_WP2_OK(custom_input.A.Resize(pad_a ? padded_width : src.width,
                                        pad_a ? padded_height : src.height));
  }

  CSPTransform csp_transform;  // Match Encode() behavior.
  constexpr int16_t kRgbAvg[3] = {0, 0, 0};
  if (matrix == Matrix::kRGB8 || matrix == Matrix::kRGB10) {
    const int16_t v =
        1u << (CSPTransform::kMtxShift - ((matrix == Matrix::kRGB8) ? 0 : 2));
    const int16_t kRGBToRGBMatrix[9] = {v, 0, 0, 0, v, 0, 0, 0, v};  // Identity
    ASSERT_TRUE(csp_transform.Init(kRGBToRGBMatrix, kRgbAvg));
  } else if (matrix == Matrix::kYCoCg) {
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

  for (uint32_t y = 0; y < padded_height; ++y) {
    for (uint32_t x = 0; x < padded_width; ++x) {
      if (x < src.width && y < src.height) {
        const uint32_t px = x * WP2FormatBpp(src.format);
        const uint8_t* const pixel = &((uint8_t*)src.GetRow(y))[px];
        if (encode_alpha) custom_input.A.At(x, y) = pixel[0];
        csp_transform.ToYUV(pixel[1], pixel[2], pixel[3],
                            &custom_input.Y.At(x, y), &custom_input.U.At(x, y),
                            &custom_input.V.At(x, y));
      }
    }
  }
  for (Channel c : {kYChannel, kUChannel, kVChannel, kAChannel}) {
    if (custom_input.GetChannel(c).IsEmpty()) continue;
    ASSERT_WP2_OK(custom_input.GetChannel(c).FillPad(src.width, src.height));
  }

  MemoryWriter data;
  EncoderConfig config = EncoderConfig::kDefault;
  config.quality = quality;
  Metadata empty_metadata;
  ASSERT_WP2_OK(Encode(
      src.width, src.height,
      custom_input.Y.Row(0), custom_input.Y.Step(),
      custom_input.U.Row(0), custom_input.U.Step(),
      custom_input.V.Row(0), custom_input.V.Step(),
      custom_input.HasAlpha() ? custom_input.A.Row(0) : nullptr,
      custom_input.A.Step(), ccsp_to_rgb_matrix, ccsp_to_rgb_shift, &data,
      config, encode_metadata ? src.metadata : empty_metadata));

  // Decode the reference and the custom bitstreams.
  ArgbBuffer ref, dst;
  ASSERT_WP2_OK(Decode(ref_data.mem_, ref_data.size_, &ref));
  ASSERT_WP2_OK(Decode(data.mem_, data.size_, &dst));

  ASSERT_TRUE(testing::Compare(src, ref, src_file_name,
                               testing::GetExpectedDistortion(quality)));
  ASSERT_TRUE(testing::Compare(src, dst, src_file_name,
                               testing::GetExpectedDistortion(quality)));
  ASSERT_TRUE(testing::Compare(
      ref, dst, src_file_name,
      (matrix == Matrix::kRGB8 || matrix == Matrix::kRGB10) ? 99.f : 50.f));

  if (encode_metadata) {
    EXPECT_TRUE(testing::HasSameData(dst.metadata.exif, ref.metadata.exif));
    EXPECT_TRUE(testing::HasSameData(dst.metadata.xmp, ref.metadata.xmp));
    EXPECT_TRUE(testing::HasSameData(dst.metadata.iccp, ref.metadata.iccp));
  }
}

//------------------------------------------------------------------------------

INSTANTIATE_TEST_SUITE_P(
    EncodeCustomCspTestInstantiation, EncodeCustomCspTest,
    ::testing::Combine(::testing::Values("source1_1x1.png", "source1_64x48.png",
                                         "test_exif_xmp.webp"),
                       ::testing::Values(0.f, 100.f) /* quality */,
                       ::testing::Values(Matrix::kRGB8, Matrix::kRGB10,
                                         Matrix::kYCoCg, Matrix::kYpUV),
                       ::testing::Values(false, true) /* encode_alpha */,
                       ::testing::Values(true) /* pad_y */,
                       ::testing::Values(true) /* pad_u */,
                       ::testing::Values(true) /* pad_v */,
                       ::testing::Values(true) /* pad_a */,
                       ::testing::Values(false) /* encode_metadata */));

INSTANTIATE_TEST_SUITE_P(
    EncodeCustomCspTestInstantiationPadding, EncodeCustomCspTest,
    ::testing::Combine(::testing::Values("source1_1x1.png"),
                       ::testing::Values(0.f, 100.f) /* quality */,
                       ::testing::Values(Matrix::kRGB8),
                       ::testing::Values(false, true) /* encode_alpha */,
                       ::testing::Values(false, true) /* pad_y */,
                       ::testing::Values(false, true) /* pad_u */,
                       ::testing::Values(false, true) /* pad_v */,
                       ::testing::Values(false, true) /* pad_a */,
                       ::testing::Values(false) /* encode_metadata */));

INSTANTIATE_TEST_SUITE_P(
    EncodeCustomCspTestInstantiationMetadata, EncodeCustomCspTest,
    ::testing::Combine(::testing::Values("source1_1x1.png",
                                         "test_exif_xmp.webp"),
                       ::testing::Values(50.f) /* quality */,
                       ::testing::Values(Matrix::kRGB8),
                       ::testing::Values(true) /* encode_alpha */,
                       ::testing::Values(false) /* pad_y */,
                       ::testing::Values(true) /* pad_u */,
                       ::testing::Values(false) /* pad_v */,
                       ::testing::Values(true) /* pad_a */,
                       ::testing::Values(true) /* encode_metadata */));

//------------------------------------------------------------------------------
// Encode the bytes read from a y4m file (YCbCr).

class ReadCustomCspTest
    : public ::testing::TestWithParam<std::tuple<const char*, float, bool>> {};

TEST_P(ReadCustomCspTest, Y4M) {
  const std::string file_name = std::get<0>(GetParam());
  const float quality = std::get<1>(GetParam());
  const bool use_animation_encoder = std::get<2>(GetParam());

  YUVPlane ccsp;
  CSPMtx ccsp_to_rgb = {};
  Metadata metadata;
  ASSERT_WP2_OK(ReadImage(testing::GetTestDataPath(file_name).c_str(), &ccsp,
                          &ccsp_to_rgb, &metadata));
  // Limit canvas size to avoid a timeout.
  ASSERT_WP2_OK(ccsp.SetView(ccsp, {0, 0, std::min(ccsp.GetWidth(), 255u),
                                    std::min(ccsp.GetHeight(), 254u)}));

  ArgbBuffer converted;
  ASSERT_WP2_OK(
      ccsp.Export(ccsp_to_rgb, /*resize_if_needed=*/true, &converted));

  EncoderConfig config = EncoderConfig::kDefault;
  config.quality = quality;
  MemoryWriter data;
  if (use_animation_encoder) {
    ASSERT_WP2_OK(
        Encode(ccsp.GetWidth(), ccsp.GetHeight(),
               ccsp.Y.Row(0), ccsp.Y.Step(),
               ccsp.U.Row(0), ccsp.U.Step(),
               ccsp.V.Row(0), ccsp.V.Step(),
               ccsp.HasAlpha() ? ccsp.A.Row(0) : nullptr, ccsp.A.Step(),
               ccsp_to_rgb.mtx(), ccsp_to_rgb.shift, &data, config, metadata));
  } else {
    AnimationEncoder encoder;
    for (uint32_t duration_ms : {1, 2}) {
      ASSERT_WP2_OK(encoder.AddFrame(
          ccsp.GetWidth(), ccsp.GetHeight(),
          ccsp.Y.Row(0), ccsp.Y.Step(),
          ccsp.U.Row(0), ccsp.U.Step(),
          ccsp.V.Row(0), ccsp.V.Step(),
          ccsp.HasAlpha() ? ccsp.A.Row(0) : nullptr, ccsp.A.Step(),
          ccsp_to_rgb.mtx(), ccsp_to_rgb.shift, duration_ms));
    }
    ASSERT_WP2_OK(encoder.Encode(&data, config, kInfiniteLoop));
  }

  ArgbBuffer decoded;
  ASSERT_WP2_OK(Decode(data.mem_, data.size_, &decoded));

  EXPECT_TRUE(testing::Compare(converted, decoded, file_name,
                               testing::GetExpectedDistortion(quality)));
}

INSTANTIATE_TEST_SUITE_P(
    ReadCustomCspTestInstantiation, ReadCustomCspTest,
    ::testing::Combine(
        ::testing::Values("source1_64x48.png", "ccsp/source3_C444p8.y4m",
                          "ccsp/source3_C444p12.y4m"),
        ::testing::Values(20.f, 100.f), /* quality */
        ::testing::Values(false, true) /* use_animation_encoder */));

//------------------------------------------------------------------------------

static struct BufferParam {
  const char* label;
  const char* file_name;
  uint32_t steps[4];
  WP2Status expected_status;
} const kBufferParams[] = {
  { "0-sized buffers are nullptr.",
    "source1_1x1.png", {0, 0, 0, 0}, WP2_STATUS_NULL_PARAMETER },
  { "0-sized buffers are nullptr.",
    "source1_1x1.png", {0, 2, 2, 2}, WP2_STATUS_NULL_PARAMETER },
  { "0-sized buffers are nullptr.",
    "source1_1x1.png", {2, 0, 2, 2}, WP2_STATUS_NULL_PARAMETER },
  { "0-sized buffers are nullptr.",
    "source1_1x1.png", {2, 2, 0, 2}, WP2_STATUS_NULL_PARAMETER },
  { "OK strides.",
    "source1_1x1.png", {1, 1, 1, 2}, WP2_STATUS_OK },
  { "R-Stride is too small.",
    "source1_32x32.png", {31, 32, 32, 0}, WP2_STATUS_BAD_DIMENSION },
  { "G-Stride is too small.",
    "source1_32x32.png", {32, 31, 32, 0}, WP2_STATUS_BAD_DIMENSION },
  { "B-Stride is too small.",
    "source1_32x32.png", {32, 32, 31, 0}, WP2_STATUS_BAD_DIMENSION },
  { "Strides are fine I.",
    "source1_32x32.png", {32, 32, 32, 0}, WP2_STATUS_OK },
  { "Strides are fine II.",
    "source1_32x32.png", {33, 35, 36, 0}, WP2_STATUS_OK }
};

class StrideTest : public ::testing::TestWithParam<BufferParam> {};

TEST_P(StrideTest, Combination) {
  const BufferParam& p = GetParam();

  ArgbBuffer src;
  ASSERT_WP2_OK(ReadImage(testing::GetTestDataPath(p.file_name).c_str(), &src));

  Vector_s16 rgba[4];
  for (uint32_t c : {0, 1, 2, 3}) {
    if (p.steps[c] == 0) continue;
    ASSERT_TRUE(rgba[c].resize(p.steps[c] * src.height));
    const uint32_t max_width = std::min(src.width, p.steps[c]);
    for (uint32_t y = 0; y < src.height; ++y) {
      const uint8_t* const src_row = (const uint8_t*)src.GetRow(y);
      int16_t* const dst_row = &rgba[c][p.steps[c] * y];
      for (uint32_t x = 0; x < max_width; ++x) {
        dst_row[x] = src_row[x * WP2FormatBpp(src.format) + c];
      }
    }
    // Pixels outside of the image are undefined on purpose for the sanitizers.
  }

  MemoryWriter writer;
  ASSERT_EQ(Encode(src.width, src.height,
                   rgba[0].data(), p.steps[0], rgba[1].data(), p.steps[1],
                   rgba[2].data(), p.steps[2], rgba[3].data(), p.steps[3],
                   kRGB8ToRGBMatrix, kRGB8ToRGBShift, &writer),
            p.expected_status) << p.label;
}

INSTANTIATE_TEST_SUITE_P(
    StrideTestInstantiation, StrideTest, ::testing::ValuesIn(kBufferParams));

//------------------------------------------------------------------------------

}  // namespace
}  // namespace WP2
