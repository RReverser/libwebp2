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

// Encoding from and to custom color space test.

#include <string>
#include <tuple>

#include "extras/extras.h"
#include "imageio/image_dec.h"
#include "include/helpers.h"
#include "src/utils/plane.h"
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
constexpr const int16_t* const kRGBToCCSPMatrices[] = {
    kRGBToRGB8Matrix, kRGBToRGB10Matrix, kRGBToYCoCgMatrix, kRGBToYpUVMatrix};
constexpr const uint32_t kRGBToCCSPShifts[] = {
    kRGBToRGB8Shift, kRGBToRGB10Shift, kRGBToYCoCgShift, kRGBToYpUVShift};

class EncodeDecodeCustomCspTest
    : public ::testing::TestWithParam<
          std::tuple<std::string, float, Matrix, bool>> {};

//------------------------------------------------------------------------------

TEST_P(EncodeDecodeCustomCspTest, Comparison) {
  const std::string& file_name = std::get<0>(GetParam());
  const float quality = std::get<1>(GetParam());
  const Matrix matrix = std::get<2>(GetParam());
  const int16_t* const ccsp_to_rgb_matrix = kCCSPToRGBMatrices[(int)matrix];
  const uint32_t ccsp_to_rgb_shift = kCCSPToRGBShifts[(int)matrix];
  const int16_t* const rgb_to_ccsp_matrix = kRGBToCCSPMatrices[(int)matrix];
  const uint32_t rgb_to_ccsp_shift = kRGBToCCSPShifts[(int)matrix];
  const bool encode_alpha = std::get<3>(GetParam());

  ArgbBuffer src;
  ASSERT_WP2_OK(ReadImage(testing::GetTestDataPath(file_name).c_str(), &src));

  YUVPlane custom_input, custom_output;  // Might contain RGB.

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

  ASSERT_WP2_OK(custom_input.Import(src, encode_alpha, csp_transform,
                                    /*resize_if_needed=*/true, /*pad=*/1));
  ASSERT_WP2_OK(
      custom_output.Resize(src.width, src.height, /*pad=*/1, encode_alpha));

  MemoryWriter data;
  EncoderConfig config = EncoderConfig::kDefault;
  config.quality = quality;
  ASSERT_WP2_OK(Encode(
      src.width, src.height,
      custom_input.Y.Row(0), custom_input.Y.Step(),
      custom_input.U.Row(0), custom_input.U.Step(),
      custom_input.V.Row(0), custom_input.V.Step(),
      custom_input.HasAlpha() ? custom_input.A.Row(0) : nullptr,
      custom_input.A.Step(), ccsp_to_rgb_matrix, ccsp_to_rgb_shift, &data,
      config));

  ASSERT_WP2_OK(Decode(
      data.mem_, data.size_, rgb_to_ccsp_matrix, rgb_to_ccsp_shift,
      custom_output.Y.Row(0), custom_output.Y.Step(), custom_output.Y.Size(),
      custom_output.U.Row(0), custom_output.U.Step(), custom_output.U.Size(),
      custom_output.V.Row(0), custom_output.V.Step(), custom_output.V.Size(),
      custom_output.HasAlpha() ? custom_output.A.Row(0) : nullptr,
      custom_output.A.Step(), custom_output.A.Size()));

  // Even if the color space is custom, the actual source is RGB so there is no
  // loss due to conversion as there is no real precision beyond 8 bits. However
  // there might be some slight rounding differences.
  const uint32_t bit_depth = csp_transform.GetYUVPrecisionBits() + 1;
  const float expected_distortion =
      std::min(testing::GetExpectedDistortion(quality), 80.f);
  ASSERT_TRUE(testing::Compare(custom_input, custom_output, bit_depth,
                               file_name, expected_distortion));
}

//------------------------------------------------------------------------------

INSTANTIATE_TEST_SUITE_P(
    EncodeDecodeCustomCspTestAlpha, EncodeDecodeCustomCspTest,
    ::testing::Combine(::testing::Values("source1_1x1.png", "alpha_ramp.png"),
                       ::testing::Values(0.f, 100.f) /* quality */,
                       ::testing::Values(Matrix::kRGB8, Matrix::kRGB10,
                                         Matrix::kYCoCg, Matrix::kYpUV),
                       ::testing::Values(true) /* encode_alpha */));

INSTANTIATE_TEST_SUITE_P(
    EncodeDecodeCustomCspTestOpaque, EncodeDecodeCustomCspTest,
    ::testing::Combine(::testing::Values("source3_222x167.jpg"),
                       ::testing::Values(50.f, 100.f) /* quality */,
                       ::testing::Values(Matrix::kRGB8, Matrix::kRGB10,
                                         Matrix::kYCoCg, Matrix::kYpUV),
                       ::testing::Values(false) /* encode_alpha */));

//------------------------------------------------------------------------------

}  // namespace
}  // namespace WP2
