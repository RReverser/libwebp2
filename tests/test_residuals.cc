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

// Test ResidualWriter/ResidualReader.

#include "include/helpers.h"
#include "src/common/lossy/block.h"
#include "src/enc/wp2_enc_i.h"
#include "src/wp2/format_constants.h"

namespace WP2 {
namespace {

//------------------------------------------------------------------------------

class ResidualTest : public ::testing::TestWithParam<
                         std::tuple<Channel, TrfSize, bool, bool>> {};

TEST_P(ResidualTest, Simple) {
  const Channel channel = std::get<0>(GetParam());
  const uint32_t num_channels = (uint32_t)channel + 1;
  const TrfSize tf_size = std::get<1>(GetParam());
  const bool use_pseudo_rate = std::get<2>(GetParam());
  const bool first_coeff_is_dc = std::get<3>(GetParam());

  WP2MathInit();
  SymbolsInfo symbols_info;
  ASSERT_WP2_OK(symbols_info.InitLossy(/*num_segments=*/1, ALL_RECTS,
                                       /*has_alpha=*/(num_channels >= 4),
                                       GetQualityHint(/*quality=*/75.f),
                                       /*use_aom_coeffs=*/false));
  SymbolRecorder recorder;
  ASSERT_WP2_OK(recorder.Allocate(symbols_info, /*num_records=*/0));
  Counters counters;
  ASSERT_WP2_OK(counters.Init(recorder));

  constexpr uint32_t kNumTests = 5;
  constexpr uint32_t kMaxNumCoeffs = kMaxBlockSizePix2;
  // Increasingly costly coeffs.
  constexpr int16_t kResiduals[kNumTests][kMaxNumCoeffs] = {
      {0},
      {-3},
      {42},
      {5, 1, -1},
      {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15}};

  float previous_rate;
  for (uint32_t i = 0; i < kNumTests; ++i) {
    const uint32_t num_coeffs =
        kMaxNumCoeffs -
        std::count(kResiduals[i], kResiduals[i] + kMaxNumCoeffs, 0);

    float rate;
    if (use_pseudo_rate) {
      EncodingMethod encoding_method;
      ASSERT_WP2_OK(ResidualWriter::GetRate(
          channel, num_channels, tf_size, kResiduals[i], num_coeffs,
          first_coeff_is_dc, counters.residuals(), &rate, &encoding_method));
      if (num_coeffs == 0) {
        EXPECT_EQ(encoding_method, EncodingMethod::kAllZero);
      } else if (num_coeffs == 1 && first_coeff_is_dc) {
        EXPECT_EQ(encoding_method, EncodingMethod::kDCOnly);
      } else {
        EXPECT_THAT(encoding_method,
                    ::testing::AnyOf(::testing::Eq(EncodingMethod::kMethod0),
                                     ::testing::Eq(EncodingMethod::kMethod1)));
      }
    } else {
      ASSERT_WP2_OK(ResidualWriter::GetPseudoRate(
          channel, num_channels, tf_size, kResiduals[i], num_coeffs,
          first_coeff_is_dc, counters.residuals(), &rate));
    }

    if (num_coeffs == 0) {
      EXPECT_EQ(rate, 0);
    } else {
      EXPECT_GT(rate, previous_rate);
    }
    previous_rate = rate;
  }
}

INSTANTIATE_TEST_SUITE_P(
    ResidualTestInstantiation, ResidualTest,
    ::testing::Combine(::testing::Values(kYChannel, kUChannel, kVChannel,
                                         kAChannel),
                       ::testing::Values(TRF_4x4, TRF_8x32), ::testing::Bool(),
                       ::testing::Bool()));

//------------------------------------------------------------------------------

}  // namespace
}  // namespace WP2
