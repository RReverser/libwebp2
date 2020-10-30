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
// Test CSPTransform::Import() and Export() with different CspTypes.

#include <cstdio>
#include <string>
#include <tuple>

#include "examples/example_utils.h"
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

class CspTest : public ::testing::TestWithParam<Csp> {};

TEST_P(CspTest, Ranges) {
  const Csp csp_type = GetParam();

  CSPTransform transf;
  ASSERT_WP2_OK(transf.Init(csp_type));

  int32_t min[3] = {0}, max[3] = {0};
  for (uint32_t i = 0; i < 3; ++i) {
    for (uint32_t j = 0; j < 3; ++j) {
      const int16_t rgb_min = 0 - transf.GetRgbAverage()[j];
      const int16_t rgb_max = kRgbMax - transf.GetRgbAverage()[j];

      const int32_t coeff = transf.GetRgbToYuvMatrix()[i * 3 + j];
      min[i] += coeff * ((coeff > 0) ? rgb_min : rgb_max);
      max[i] += coeff * ((coeff > 0) ? rgb_max : rgb_min);
    }
  }
  const float norm = 1.f / (1 << 12u);
  for (uint32_t i = 0; i < 3; ++i) {
    min[i] = std::round(min[i] * norm);
    max[i] = std::round(max[i] * norm);
  }
  transf.Print();
  printf("======== Min / Max ==== \n");
  for (uint32_t i = 0; i < 3; ++i) {
    const int16_t min_value = min[i];
    const int16_t max_value = max[i];
    printf("[%d %d]", min_value, max_value);
    // Check min/max values fit within YUV precision bits.
    // TODO(maryla): these are higher bounds for now but it would be nice to
    // track actual min/max values.
    EXPECT_LE(std::abs(min_value), 1 << transf.GetYUVPrecisionBits());
    EXPECT_LT(std::abs(max_value), 1 << transf.GetYUVPrecisionBits());
  }
  // Check we couldn't be using fewer bits.
  const int32_t global_min = std::min({min[0], min[1], min[2]});
  const int32_t global_max = std::max({max[0], max[1], max[2]});
  EXPECT_EQ(global_min, transf.GetYUVMin());
  EXPECT_EQ(global_max, transf.GetYUVMax());
  EXPECT_GT(std::abs(global_min), 1 << (transf.GetYUVPrecisionBits() - 1));
  EXPECT_GT(std::abs(global_max), 1 << (transf.GetYUVPrecisionBits() - 1));

  printf("\nNorms:\n  {");
  for (uint32_t i = 0; i < 3; ++i) {
    printf("%.3f, ", (float)(max[i] - min[i]) / kRgbMax);
  }
  printf("}\n");

  uint32_t kExpectedYUVBits[kNumCspTypes] = {/* YCoCg */ 9, /* YCbCr */ 9,
                                             /* Custom, not tested */ 0,
                                             /* YIQ */ 9};
  EXPECT_EQ(transf.GetYUVPrecisionBits(), kExpectedYUVBits[(uint32_t)csp_type]);
}

TEST_P(CspTest, Precision) {
  const Csp csp_type = GetParam();
  // YCoCg is the only lossless conversion from and to RGB.
  const int16_t error_margin = (csp_type == Csp::kYCoCg) ? 0 : 1;

  CSPTransform transf;
  ASSERT_WP2_OK(transf.Init(csp_type));
  double error_sum = 0., error_avg = 0.;
  uint32_t num_values = 0;

  // Test some Argb combinations. It is easier than testing YUV combinations
  // that lead to valid Argb.
  constexpr int32_t kStep = 7;
  for (int16_t r = 0; r < (int16_t)kRgbMax + kStep; r += kStep) {
    for (int16_t g = 0; g < (int16_t)kRgbMax + kStep; g += kStep) {
      for (int16_t b = 0; b < (int16_t)kRgbMax + kStep; b += kStep) {
        // Make sure valid extremes are tested.
        const int16_t in_r = std::min(r, (int16_t)kRgbMax);
        const int16_t in_g = std::min(g, (int16_t)kRgbMax);
        const int16_t in_b = std::min(b, (int16_t)kRgbMax);

        int16_t in_y, in_u, in_v;
        transf.ToYUV(in_r, in_g, in_b, &in_y, &in_u, &in_v);

        int16_t out_r, out_g, out_b;
        transf.ToRGB(in_y, in_u, in_v, &out_r, &out_g, &out_b);
        error_sum += out_r - in_r + out_g - in_g + out_b - in_b;
        error_avg += std::abs(out_r - in_r) + std::abs(out_g - in_g) +
                     std::abs(out_b - in_b);
        num_values += 3;

        ASSERT_NEAR(out_r, in_r, error_margin);
        ASSERT_NEAR(out_g, in_g, error_margin);
        ASSERT_NEAR(out_b, in_b, error_margin);
      }
    }
  }

  error_sum /= num_values;
  error_avg /= num_values;
  EXPECT_LT(std::abs(error_sum), 0.0001);  // This should be almost 0.
  EXPECT_LT(error_avg, 0.05);  // Less than 5% of values have an error of +-1.
}

INSTANTIATE_TEST_SUITE_P(CspTestInstantiation, CspTest,
                         ::testing::Values(Csp::kYCoCg, Csp::kYCbCr,
                                           Csp::kYIQ));

//------------------------------------------------------------------------------

class CspImageTest
    : public ::testing::TestWithParam<std::tuple<std::string, Csp>> {};

TEST_P(CspImageTest, Simple) {
  const std::string& src_file_name = std::get<0>(GetParam());
  const Csp csp_type = std::get<1>(GetParam());

  ArgbBuffer src;
  ASSERT_WP2_OK(
      ReadImage(testing::GetTestDataPath(src_file_name).c_str(), &src));

  CSPTransform transf;
  ASSERT_WP2_OK(transf.Init(csp_type, src));
  transf.Print();

  YUVPlane yuv;
  ASSERT_WP2_OK(yuv.Import(src, src.HasTransparency(), transf,
                           /*resize_if_needed=*/true, /*pad=*/kPredWidth));

  ArgbBuffer dst;
  ASSERT_WP2_OK(yuv.Export(transf, /*resize_if_needed=*/true, &dst));

  float disto[5];
  ASSERT_WP2_OK(dst.GetDistortion(src, PSNR, disto));

  // Alpha is not supported, total is not either.
  const float min_in_disto = std::min({disto[1], disto[2], disto[3]});
  const float min_expected_disto = (csp_type == Csp::kYCoCg) ? 99.f : 45.f;
  EXPECT_GE(min_in_disto, min_expected_disto)
      << std::endl
      << "Distortion A " << disto[0] << ", R " << disto[1] << ", G " << disto[2]
      << ", B " << disto[3] << ", total " << disto[4] << " for file "
      << src_file_name << ", csp " << (uint32_t)csp_type;
}

INSTANTIATE_TEST_SUITE_P(
    CspImageTestInstantiation, CspImageTest,
    ::testing::Combine(::testing::Values("source1_64x48.png"),
                       ::testing::Values(Csp::kYCoCg, Csp::kYCbCr, Csp::kCustom,
                                         Csp::kYIQ)));

// This one takes a while to run so it is disabled.
// Can still be run with flag --test_arg=--gunit_also_run_disabled_tests
INSTANTIATE_TEST_SUITE_P(
    DISABLED_CspImageTestInstantiation, CspImageTest,
    ::testing::Combine(::testing::Values("source0.pgm", "source0.ppm",
                                         "source4.webp", "test_exif_xmp.webp"),
                       ::testing::Values(Csp::kYCoCg, Csp::kYCbCr, Csp::kCustom,
                                         Csp::kYIQ)));

//------------------------------------------------------------------------------
// This test exercizes a value underflow/overflow in ToYUV().
// RGB are multiplied with RGBtoYUV matrix and after rounding, Y is outside the
// maximum allowed values of +/-1023 (10 bits + sign). RGBtoYUV matrix float
// computation and rounding to int is probably the root cause.

class CustomCspTest : public ::testing::TestWithParam<const char*> {};

TEST_P(CustomCspTest, RoundingError) {
  const char* const file_name = GetParam();
  EncoderConfig config = EncoderConfig::kDefault;
  config.csp_type = Csp::kCustom;
  ASSERT_WP2_OK(testing::EncodeDecodeCompare(file_name, config));
}

INSTANTIATE_TEST_SUITE_P(CustomCspTestInstantiation, CustomCspTest,
                         ::testing::Values("source1_64x48.png"));

// This one takes a while to run so it is disabled.
// Can still be run with flag --test_arg=--gunit_also_run_disabled_tests
INSTANTIATE_TEST_SUITE_P(DISABLED_CustomCspTestInstantiation, CustomCspTest,
                         ::testing::Values("source0.pgm", "source0.ppm",
                                           "source1.png", "source3.jpg",
                                           "source4.webp"));

//------------------------------------------------------------------------------
// Make sure premultiplied color components are at most equal to alpha.

class DecodeArgbTest
    : public ::testing::TestWithParam<std::tuple<std::string, float>> {};

TEST_P(DecodeArgbTest, Valid) {
  const std::string& src_file_name = std::get<0>(GetParam());
  const float quality = std::get<1>(GetParam());

  ArgbBuffer src(WP2_Argb_32), dst(WP2_Argb_32);
  MemoryWriter data;
  ASSERT_WP2_OK(testing::CompressImage(src_file_name, &data, &src, quality));
  ASSERT_WP2_OK(Decode(data.mem_, data.size_, &dst));
  ASSERT_TRUE(testing::Compare(src, dst, src_file_name,
                               testing::GetExpectedDistortion(quality)));

  for (const ArgbBuffer* const buffer : {&src, &dst}) {
    for (uint32_t y = 0; y < buffer->height; ++y) {
      const uint8_t* pixel = (const uint8_t*)buffer->GetRow(y);
      for (uint32_t x = 0; x < buffer->width; ++x) {
        ASSERT_GE(pixel[0], pixel[1]);
        ASSERT_GE(pixel[0], pixel[2]);
        ASSERT_GE(pixel[0], pixel[3]);
        pixel += 4;
      }
    }
  }
}

INSTANTIATE_TEST_SUITE_P(
    DecodeArgbTestInstantiation, DecodeArgbTest,
    ::testing::Combine(::testing::Values("source1_4x4.png"),
                       ::testing::Values(0.f, 100.f) /* quality */));

//------------------------------------------------------------------------------
// Tests CSPTransform::MakeEigenMatrixEncodable() results.

struct MtxValidity {
  double error;                  // Maximum offset among all elements.
  std::array<double, 9> matrix;  // Elements.
};

class CspTestCustomCsp
    : public ::testing::TestWithParam<MtxValidity> {};

TEST_P(CspTestCustomCsp, EigenMtx) {
  const MtxValidity& param = GetParam();

  std::array<int16_t, 9> fixed_point_matrix;
  for (uint32_t i = 0; i < 9; ++i) {
    fixed_point_matrix[i] =
        (int16_t)std::lround(param.matrix[i] * (1 << CSPTransform::kMtxShift));
  }

  int32_t error;
  std::array<int16_t, 9> encodable_matrix;
  ASSERT_WP2_OK(CSPTransform::MakeEigenMatrixEncodable(
      fixed_point_matrix.data(), encodable_matrix.data(), &error));

  EXPECT_NEAR(error / (double)(1 << CSPTransform::kMtxShift), param.error,
              0.0001);
}

INSTANTIATE_TEST_SUITE_P(
    CspTestCustomCspInstantiation, CspTestCustomCsp,
    ::testing::Values(  // error, matrix
        MtxValidity{0.,
                    {0., 0., 0.,  // Empty matrix
                     0., 0., 0.,  //
                     0., 0., 0.}},
        MtxValidity{0.,
                    {1.,  0.,  0.,  // Identity matrix
                     0., -1.,  0.,  //
                     0.,  0., -1.}},
        MtxValidity{0.9988,
                    {1.,  0., 0.,  // Bad matrix
                     0., -1., 0.,  //
                     0.,  0., 0.}},
        MtxValidity{0.1399,
                    {1.,  0.,   0.,  // Bad matrix
                     0., -0.99, 0.,  //
                     0.,  0.,  -1.}},
        MtxValidity{0.,
                    {0.263184, 0.321289, 0.199707,   // Valid output from
                     0.27002, 0.0107422, -0.373291,  // NormalizeFromInput()
                     0.264893, -0.330078, 0.181885}},
        MtxValidity{0.526367,
                    {-0.263184, 0.321289, 0.199707,  // Same but flipped
                     0.27002, 0.0107422, -0.373291,  // first sign.
                     0.264893, -0.330078, 0.181885}},
        MtxValidity{0.,
                    {0.264893, 0.374023, 0.115967,     // Valid
                     0.269531, -0.0727539, -0.381592,  //
                     0.28418, -0.279785, 0.253906}},
        MtxValidity{0.,
                    {0.287354, 0.217773, 0.22168,      // Valid
                     0.261719, -0.00708008, -0.33252,  //
                     0.167236, -0.362793, 0.139404}},
        MtxValidity{0.,
                    {0.333740, 0.265136, 0.364501,   // Valid
                     0.329101, 0.166503, -0.422363,  //
                     0.307861, -0.465087, 0.056152}},
        MtxValidity{0.014648,
                    {0.285400, 0.240234, 0.331543,   // Valid output but too
                     0.288330, -0.399170, 0.029297,  // imprecise to encode.
                     0.282227, 0.164795, -0.365234}}));

//------------------------------------------------------------------------------

}  // namespace
}  // namespace WP2
