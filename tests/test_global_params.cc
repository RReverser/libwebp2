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

#include "include/helpers.h"
#include "src/common/global_params.h"
#include "src/enc/wp2_enc_i.h"

using ::testing::AssertionResult;

namespace WP2 {
namespace {

WP2Status WriteAndRead(const GlobalParams& params_to_write,
                       GlobalParams* read_params) {
  ANSEnc enc;
  constexpr int kQualityHint = 8;
  constexpr bool kImageHasAlpha = true;
  WP2_CHECK_STATUS(params_to_write.Write(kQualityHint, kImageHasAlpha, &enc));
  WP2_CHECK_STATUS(enc.Assemble());
  ExternalDataSource data_source(enc.Buffer(), enc.BufferSize());
  ANSDec dec(&data_source);
  WP2_CHECK_STATUS(read_params->Read(kQualityHint, kImageHasAlpha, &dec));
  WP2_CHECK_OK(read_params->IsOk(), WP2_STATUS_BAD_READ);
  return WP2_STATUS_OK;
}

TEST(GlobalParmas, ReadWriteUVQUantMultiplier) {
  GlobalParams params;
  ASSERT_TRUE(params.segments_.resize(3));
  const uint32_t lambda = 1000;
  params.u_quant_multiplier_ = kDefaultQuantMultiplier;
  params.v_quant_multiplier_ = kDefaultQuantMultiplier + 2;
  for (Segment& s : params.segments_) {
    s.SetYUVBounds(-512, 512);
    s.SetQuality(kMaxLossyQualityHint / 2, 10, params.u_quant_multiplier_,
                 params.v_quant_multiplier_);
    s.FinalizeQuant(lambda);
  }
  ASSERT_TRUE(params.IsOk());

  GlobalParams read_params;
  ASSERT_WP2_OK(WriteAndRead(params, &read_params));
  EXPECT_EQ(read_params.u_quant_multiplier_, kDefaultQuantMultiplier);
  EXPECT_EQ(read_params.v_quant_multiplier_, kDefaultQuantMultiplier + 2);

  params.u_quant_multiplier_ = kDefaultQuantMultiplier;
  params.v_quant_multiplier_ = kDefaultQuantMultiplier - 3;
  ASSERT_WP2_OK(WriteAndRead(params, &read_params));
  EXPECT_EQ(read_params.u_quant_multiplier_, kDefaultQuantMultiplier);
  EXPECT_EQ(read_params.v_quant_multiplier_, kDefaultQuantMultiplier - 3);

  params.u_quant_multiplier_ = kDefaultQuantMultiplier - 2;
  params.v_quant_multiplier_ = kDefaultQuantMultiplier + 1;
  ASSERT_WP2_OK(WriteAndRead(params, &read_params));
  EXPECT_EQ(read_params.u_quant_multiplier_, kDefaultQuantMultiplier - 2);
  EXPECT_EQ(read_params.v_quant_multiplier_, kDefaultQuantMultiplier + 1);
}

TEST(GlobalParmas, FixedPredictors) {
  GlobalParams params;
  ASSERT_TRUE(params.segments_.resize(3));
  ASSERT_TRUE(params.a_segments_.resize(1));
  params.has_alpha_ = params.maybe_use_lossy_alpha_ = true;
  // Init twice (should not matter).
  ASSERT_WP2_OK(params.InitFixedPredictors());
  ASSERT_WP2_OK(params.InitFixedPredictors());
  EXPECT_EQ(params.y_preds_.size(), kYPredNum);
  EXPECT_EQ(params.y_preds_.GetMaxMode(), kYPredModeNum - 1);
  EXPECT_EQ(params.uv_preds_.size(), kUVPredNum);
  EXPECT_EQ(params.a_preds_.size(), kAPredNum);
  EXPECT_EQ(params.a_preds_.GetMaxMode(), kAPredModeNum - 1);

  GlobalParams read_params;
  ASSERT_WP2_OK(WriteAndRead(params, &read_params));
  EXPECT_EQ(read_params.y_preds_.size(), kYPredNum);
  EXPECT_EQ(read_params.y_preds_.GetMaxMode(), kYPredModeNum - 1);
  EXPECT_EQ(read_params.uv_preds_.size(), kUVPredNum);
  EXPECT_EQ(read_params.a_preds_.size(), kAPredNum);
  EXPECT_EQ(read_params.a_preds_.GetMaxMode(), kAPredModeNum - 1);
}

}  // namespace
}  // namespace WP2
