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

// Partition scoring function test.

#include <string>

#include "imageio/image_dec.h"
#include "include/helpers.h"
#include "src/common/global_params.h"
#include "src/enc/partition_score_func.h"
#include "src/utils/csp.h"
#include "src/utils/plane.h"
#include "src/wp2/decode.h"
#include "src/wp2/encode.h"

namespace WP2 {
namespace {

//------------------------------------------------------------------------------

TEST(PartitionScoreFuncTest, BlockScoreFunc) {
  ArgbBuffer rgb;
  ASSERT_WP2_OK(
      ReadImage(testing::GetTestDataPath("source1_32x32.png").c_str(), &rgb));

  const EncoderConfig& config = EncoderConfig::kDefault;
  CSPTransform transf;
  YUVPlane yuv;
  ASSERT_WP2_OK(yuv.Import(rgb, /*import_alpha=*/rgb.HasTransparency(), transf,
                           /*resize_if_needed=*/true, /*pad=*/kPredWidth));
  GlobalParams gparams;
  FeatureMap features;
  gparams.features_ = &features;
  ASSERT_WP2_OK(GlobalAnalysis(rgb, yuv, transf, config, &gparams));
  ASSERT_WP2_OK(gparams.InitRndMtxSet());
  ASSERT_WP2_OK(gparams.Init());

  const Rectangle tile_rect = {0, 0, yuv.GetWidth(), yuv.GetHeight()};
  const Block kBlocks[] = {{0, 0, BLK_4x4}, {1, 0, BLK_4x8}, {0, 1, BLK_4x4}};
  float score;

  BlockScoreFunc score_func;
  ASSERT_WP2_OK(score_func.Init(config, tile_rect, yuv, gparams));

  // Test block by block.
  for (const Block& block : kBlocks) {
    if (block.y() != 0) continue;  // Blocks above must be used.
    ASSERT_WP2_OK(score_func.GetScore(block, &score));
    EXPECT_GT(score, 0);
  }

  // Test a group of blocks.
  ASSERT_WP2_OK(score_func.GetScore(
      kBlocks, sizeof(kBlocks) / sizeof(kBlocks[0]), &score));
  EXPECT_GT(score, 0);
}

//------------------------------------------------------------------------------

}  // namespace
}  // namespace WP2
