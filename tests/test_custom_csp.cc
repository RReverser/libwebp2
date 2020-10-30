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

// Custom color space conversion test.

#include <tuple>
#include <vector>

#include "extras/ccsp_imageio.h"
#include "extras/extras.h"
#include "imageio/image_dec.h"
#include "include/helpers.h"
#include "src/dsp/math.h"
#include "src/utils/csp.h"
#include "src/utils/random.h"

namespace WP2 {
namespace {

//------------------------------------------------------------------------------

// Matrices and associated shifts from extras.h
enum class Matrix { kRGB8, kRGB10, kYCoCg, kYpUV };
constexpr CSPMtx kRGBToCCSPMatrices[] = {{kRGBToRGB8Matrix, kRGBToRGB8Shift},
                                         {kRGBToRGB10Matrix, kRGBToRGB10Shift},
                                         {kRGBToYCoCgMatrix, kRGBToYCoCgShift},
                                         {kRGBToYpUVMatrix, kRGBToYpUVShift}};
constexpr CSPMtx kCCSPToRGBShifts[] = {{kRGB8ToRGBMatrix, kRGB8ToRGBShift},
                                       {kRGB10ToRGBMatrix, kRGB10ToRGBShift},
                                       {kYCoCgToRGBMatrix, kYCoCgToRGBShift},
                                       {kYpUVToRGBMatrix, kYpUVToRGBShift}};

constexpr const uint32_t kCCSPBitDepth[] = {8, 10, 10, 11};

//------------------------------------------------------------------------------

class CspMatrixTest : public ::testing::TestWithParam<std::tuple<Csp, Matrix>> {
  void SetUp() override { WP2CSPConverterInit(); }
};

TEST_P(CspMatrixTest, All) {
  const Csp csp_type = std::get<0>(GetParam());
  const Matrix matrix = std::get<1>(GetParam());
  const CSPMtx& rgb_to_ccsp = kRGBToCCSPMatrices[(int)matrix];
  const CSPMtx& ccsp_to_rgb = kCCSPToRGBShifts[(int)matrix];

  // YCoCg is the only lossless conversion from and to RGB. For the others, the
  // 'CSPTransform::im_' is not exactly the inverse matrix of 'm_' and thus
  // having more precision bits for the custom color space is not enough to keep
  // the exact same channel values because there is a temporary RGB conversion.
  const int16_t error_margin = (csp_type == Csp::kYCoCg) ? 0 : 1;

  CSPTransform transf;
  ASSERT_WP2_OK(transf.Init(csp_type));

  // Test most Argb combinations. It is easier than testing YUV combinations
  // that lead to valid Argb.
  constexpr int32_t kStep = 7;
  constexpr int16_t in_a = kAlphaMax;
  for (int16_t r = 0; r < in_a + kStep; r += kStep) {
    for (int16_t g = 0; g < in_a + kStep; g += kStep) {
      for (int16_t b = 0; b < in_a + kStep; b += kStep) {
        // Make sure valid extremes are tested.
        const int16_t in_r = std::min(r, in_a);
        const int16_t in_g = std::min(g, in_a);
        const int16_t in_b = std::min(b, in_a);

        int16_t in_y, in_u, in_v;
        transf.ToYUV(in_r, in_g, in_b, &in_y, &in_u, &in_v);

        int16_t c0, c1, c2;
        const uint32_t width = 1, height = 1, step = 1;
        ASSERT_WP2_OK(transf.YUVToCustom(width, height, &in_y, step, &in_u,
                                         step, &in_v, step, rgb_to_ccsp, &c0,
                                         step, &c1, step, &c2, step));
        if (matrix == Matrix::kRGB8 || matrix == Matrix::kRGB10) {
          ASSERT_NEAR(RightShiftRound(c0, ccsp_to_rgb.shift), in_r,
                      error_margin);
          ASSERT_NEAR(RightShiftRound(c1, ccsp_to_rgb.shift), in_g,
                      error_margin);
          ASSERT_NEAR(RightShiftRound(c2, ccsp_to_rgb.shift), in_b,
                      error_margin);
        } else {
          ASSERT_LE(std::abs(c0), (1 << 11) - 1);  // Reasonable range.
          ASSERT_LE(std::abs(c1), (1 << 11) - 1);
          ASSERT_LE(std::abs(c2), (1 << 11) - 1);
        }

        int16_t out_y, out_u, out_v;
        ASSERT_WP2_OK(transf.CustomToYUV(width, height, &c0, step, &c1, step,
                                         &c2, step, ccsp_to_rgb, &out_y, step,
                                         &out_u, step, &out_v, step));

        ASSERT_NEAR(out_y, in_y, error_margin);
        ASSERT_NEAR(out_u, in_u, error_margin);
        ASSERT_NEAR(out_v, in_v, error_margin);
      }
    }
  }
}

INSTANTIATE_TEST_SUITE_P(
    CspMatrixTestInstantiation, CspMatrixTest,
    ::testing::Combine(::testing::Values(Csp::kYCoCg, Csp::kYCbCr, Csp::kYIQ),
                       ::testing::Values(Matrix::kRGB8, Matrix::kRGB10,
                                         Matrix::kYCoCg, Matrix::kYpUV)));

//------------------------------------------------------------------------------

class BackAndForthTest
    : public ::testing::TestWithParam<std::tuple<const char*, Matrix>> {
  void SetUp() override { WP2CSPConverterInit(); }
};

TEST_P(BackAndForthTest, Comparison) {
  const char* const file_name = std::get<0>(GetParam());
  const Matrix matrix = std::get<1>(GetParam());
  const CSPMtx& rgb_to_ccsp = kRGBToCCSPMatrices[(int)matrix];
  const CSPMtx& ccsp_to_rgb = kCCSPToRGBShifts[(int)matrix];
  const uint32_t ccsp_bit_depth = kCCSPBitDepth[(int)matrix];

  YUVPlane original;
  CSPMtx file_to_rgb = {};
  ASSERT_WP2_OK(ReadImage(testing::GetTestDataPath(file_name).c_str(),
                          &original, &file_to_rgb));
  if (original.IsDownsampled()) {
    ASSERT_WP2_OK(original.Upsample());
  }

  YUVPlane rgb;
  ASSERT_WP2_OK(rgb.Copy(original, /*resize_if_needed=*/true));
  const uint32_t kept_precision = std::min(file_to_rgb.shift, 1u);
  ASSERT_WP2_OK(
      rgb.Apply(file_to_rgb.mtx(), file_to_rgb.shift - kept_precision));

  YUVPlane converted;
  ASSERT_WP2_OK(converted.Copy(rgb, /*resize_if_needed=*/true));
  ASSERT_WP2_OK(converted.Apply(rgb_to_ccsp.mtx(), rgb_to_ccsp.shift));

  {
    YUVPlane bypass_rgb;
    ASSERT_WP2_OK(bypass_rgb.Copy(original, /*resize_if_needed=*/true));
    int32_t file_to_ccsp_matrix[9];
    Multiply(rgb_to_ccsp.mtx(), file_to_rgb.mtx(), file_to_ccsp_matrix);
    ASSERT_WP2_OK((bypass_rgb.Apply<int32_t, int64_t>(
        file_to_ccsp_matrix,
        /*shift=*/file_to_rgb.shift - kept_precision + rgb_to_ccsp.shift)));
    // The loss can only come from a high precision input matrix that forces a
    // right shift or a clamping of the samples (int16_t) during the RGB step.
    const float expected_distortion =
        (file_to_rgb.shift - kept_precision == 0) ? 99.f : 50.f;
    EXPECT_TRUE(testing::Compare(converted, bypass_rgb, ccsp_bit_depth,
                                 file_name, expected_distortion));
  }

  {
    YUVPlane back_to_rgb;
    ASSERT_WP2_OK(back_to_rgb.Copy(converted, /*resize_if_needed=*/true));
    ASSERT_WP2_OK(back_to_rgb.Apply(ccsp_to_rgb.mtx(), ccsp_to_rgb.shift));
    // Almost no loss thanks to the 'kept_precision' above.
    EXPECT_TRUE(testing::Compare(rgb, back_to_rgb, ccsp_bit_depth, file_name,
                                 /*expected_distortion=*/90.f));
  }
}

INSTANTIATE_TEST_SUITE_P(
    BackAndForthTestInstantiation, BackAndForthTest,
    ::testing::Combine(::testing::Values("alpha_ramp.pam", "source1_1x1.png",
                                         "source1_64x48.png", "source3.jpg",
                                         "ccsp/source3_C420p8.y4m",
                                         "ccsp/source3_C420p10.y4m",
                                         "ccsp/source3_C420p12.y4m",
                                         "ccsp/source3_C444p8.y4m",
                                         "ccsp/source3_C444p10.y4m",
                                         "ccsp/source3_C444p12.y4m"),
                       ::testing::Values(Matrix::kRGB8, Matrix::kRGB10,
                                         Matrix::kYCoCg, Matrix::kYpUV)));

//------------------------------------------------------------------------------

TEST(CspCustomTest, Simple) {
  const std::string& file_name = "composited.png";

  // Make sure Csp::kCustom is going to be used and not the default Csp::kYCoCg.
  ArgbBuffer img;
  ASSERT_WP2_OK(ReadImage(testing::GetTestDataPath(file_name).c_str(), &img));
  CSPTransform transf;
  ASSERT_WP2_OK(transf.Init(Csp::kCustom, img));
  ASSERT_EQ(transf.GetType(), Csp::kCustom);

  EncoderConfig config = EncoderConfig::kDefault;
  config.csp_type = Csp::kCustom;
  ASSERT_WP2_OK(testing::EncodeDecodeCompare(file_name, config));
}

//------------------------------------------------------------------------------

class SpeedTest : public ::testing::TestWithParam<WP2CPUInfo> {
  void SetUp() override {
    WP2DspReset();
    WP2GetCPUInfo = GetParam();
    WP2CSPConverterInit();
  }

 protected:
  static constexpr uint32_t kIter = 100u;
  static constexpr uint32_t kNumTest = 10000u / kIter;
  static constexpr uint32_t kMaxLen = 16384u;
  static constexpr int32_t kMaxAvg = (1u << 10) - 1;
  static constexpr int32_t kMaxYUV = (1u << 10) - 1;
  static constexpr int32_t kMaxMtx = 1u << 14;
  std::vector<int16_t> yuv_;
  int16_t avg_[3], mtx_[9];

  void Init(UniformIntDistribution* const random) {
    yuv_.resize(3 * kMaxLen);
    for (auto& v : yuv_) {
      // from time to time, inject some extreme values
      const int rnd = random->Get(0, 30);
      v = (rnd == 0) ?  kMaxYUV
        : (rnd == 1) ? -kMaxYUV
        : random->Get<int16_t>(-kMaxYUV, kMaxYUV);
    }
  }
  void FillMtx(UniformIntDistribution* const random) {
    for (auto& i : avg_) i = random->Get(-kMaxAvg, kMaxAvg);
    for (auto& i : mtx_) i = random->Get(-kMaxMtx, kMaxMtx);
  }

  uint32_t DoTest(WP2YUVToArgbFunc func, UniformIntDistribution* const random) {
    Init(random);
    uint32_t crc = 0;
    std::vector<uint8_t> out(4 * kMaxLen, 0);
    for (uint32_t test = 0; test < kNumTest; ++test) {
      const uint32_t len = (test == 0) ? kMaxLen : random->Get(1u, kMaxLen);
      const int16_t* const y = &yuv_[random->Get(0u, kMaxLen - len)];
      const int16_t* const u = &yuv_[random->Get(0u, kMaxLen - len)];
      const int16_t* const v = &yuv_[random->Get(0u, kMaxLen - len)];
      const int16_t* const a = &yuv_[random->Get(0u, kMaxLen - len)];
      FillMtx(random);
      for (uint32_t iter = 0; iter < kIter; ++iter) {
        func(y, u, v, a, avg_, mtx_, &out[0], len);
      }
      for (auto i : out) crc = testing::SillyCRC(crc, i);
    }
    return crc;
  }
};

TEST_P(SpeedTest, YUVToCustom) {
  UniformIntDistribution random(431412);
  Init(&random);
  std::vector<int16_t> out(3 * kMaxLen, 0);
  uint32_t crc = 0;
  for (uint32_t test = 0; test < kNumTest; ++test) {
    const uint32_t len = (test == 0) ? kMaxLen : random.Get(1u, kMaxLen);
    const int16_t* const y = &yuv_[random.Get(0u, kMaxLen - len)];
    const int16_t* const u = &yuv_[random.Get(0u, kMaxLen - len)];
    const int16_t* const v = &yuv_[random.Get(0u, kMaxLen - len)];
    FillMtx(&random);
    for (uint32_t iter = 0; iter < kIter; ++iter) {
      WP2YUVToCustom(y, u, v, avg_, mtx_,
                     &out[0 * kMaxLen], &out[1 * kMaxLen], &out[2 * kMaxLen],
                     len);
    }
    for (auto i : out) crc = testing::SillyCRC(crc, i);
  }
  EXPECT_EQ(crc, 1929850794u);
}

TEST_P(SpeedTest, YUVToXRGB) {
  UniformIntDistribution random(1532);
  EXPECT_EQ(DoTest(WP2YUVToXRGB, &random), 806577452u);
}

TEST_P(SpeedTest, YUVToARGB) {
  UniformIntDistribution random(8724);
  EXPECT_EQ(DoTest(WP2YUVToARGB, &random), 3144549024u);
}

TEST_P(SpeedTest, YUVToArgb) {
  UniformIntDistribution random(48623);
  EXPECT_EQ(DoTest(WP2YUVToArgb, &random), 3017869648u);
}

INSTANTIATE_TEST_SUITE_P(SpeedTestInstantiation, SpeedTest,
    ::testing::ValuesIn(testing::kWP2CpuInfos));

//------------------------------------------------------------------------------

}  // namespace
}  // namespace WP2
