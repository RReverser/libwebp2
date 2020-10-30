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

// Test SyntaxWriter class.

#include "imageio/image_dec.h"
#include "include/helpers.h"
#include "src/common/global_params.h"
#include "src/enc/wp2_enc_i.h"
#include "src/utils/ans.h"
#include "src/utils/vector.h"
#include "src/wp2/base.h"
#include "src/wp2/encode.h"

namespace WP2 {
namespace {

//------------------------------------------------------------------------------
// Verifies SyntaxWriter::CopyFrom() creates a valid deep clone.

TEST(TestSyntaxWriter, CopyFrom) {
  WP2DspReset();
  WP2QuantizeInit();
  WP2MathInit();
  ANSInit();

  // Settings.
  const EncoderConfig& config = EncoderConfig::kDefault;
  const Rectangle tile_rect = {0, 0, 64, 64};
  const std::array<BlockSize, 4> kBlocks = {BLK_8x16, BLK_4x4, BLK_4x4,
                                            BLK_16x8};
  const uint32_t kCopyAtBlockIndex = kBlocks.size() / 2;

  // Load the sample image and convert it to luma/chroma.
  ArgbBuffer rgb;
  ASSERT_WP2_OK(
      ReadImage(testing::GetTestDataPath("source1.png").c_str(), &rgb));
  ASSERT_WP2_OK(rgb.SetView(rgb, tile_rect));
  CSPTransform csp_transform;
  YUVPlane yuv, out;
  ASSERT_WP2_OK(yuv.Import(rgb, rgb.HasTransparency(), csp_transform,
                           /*resize_if_needed=*/true));
  ASSERT_WP2_OK(
      out.Resize(yuv.GetWidth(), yuv.GetHeight(), /*pad=*/1, yuv.HasAlpha()));

  // Create and initialize mockups.
  GlobalParams gparams;
  FeatureMap features;
  gparams.features_ = &features;
  ASSERT_WP2_OK(GlobalAnalysis(rgb, yuv, csp_transform, config, &gparams));
  ContextCache context_cache;

  // Declare a SyntaxWriter instance and a clone; the latter outlives the former
  // to verify there is no hidden dependency left.
  SyntaxWriter syntax_writer_copy;
  ANSDictionaries dicts_copy;

  // Instances that will be needed after 'syntax_writer' is deleted.
  Vector<CodedBlock> cblocks;
  FrontMgrLexico mgr;
  float cost, cost_dicts;
  uint32_t num_tokens;
  Vector_u8 buffer;

  {
    SyntaxWriter syntax_writer;
    ANSDictionaries dicts;

    ASSERT_WP2_OK(syntax_writer.Init(
        &dicts, config, gparams, yuv, ChromaSubsampling::k420, tile_rect,
        (uint32_t)kBlocks.size(), /*use_aom_coeffs=*/true));

    // Encode a bunch of blocks.
    ASSERT_WP2_OK(mgr.Init(config.partition_set, /*snapped=*/false,
                           tile_rect.width, tile_rect.height));
    for (uint32_t i = 0; i < kBlocks.size(); ++i) {
      Block block;
      ASSERT_TRUE(mgr.TryGetNextBlock(kBlocks[i], &block));

      ASSERT_TRUE(cblocks.resize(cblocks.size() + 1));
      CodedBlock& cb = cblocks.back();
      cb.is420_ = false;
      cb.id_ = AssignSegmentId(config, gparams, tile_rect, block);
      cb.SetRange(gparams.transf_.GetYUVMin(), gparams.transf_.GetYUVMax());
      cb.SetDim(block, mgr);
      cb.SetPredictor(kYChannel, gparams.y_preds_[0]);
      cb.SetUVPredictor(gparams.uv_preds_[0]);
      if (!gparams.a_preds_.empty()) {
        cb.SetPredictor(kAChannel, gparams.a_preds_[0]);
      }
      cb.SetSrcInput(yuv);
      cb.SetContextInput(out, &context_cache);
      cb.SetReconstructedOutput(&out);
      cb.y_context_is_constant_ = false;

      cb.Quantize(gparams.segments_[cb.id_].quant_y_, kYChannel, false);
      cb.Quantize(gparams.segments_[cb.id_].quant_u_, kUChannel, cb.is420_);
      cb.Quantize(gparams.segments_[cb.id_].quant_v_, kVChannel, cb.is420_);
      if (!gparams.a_preds_.empty()) {
        cb.Quantize(gparams.segments_[cb.id_].quant_a_, kAChannel, false);
      }

      // All blocks before 'kCopyAtBlockIndex' are only recorded by
      // 'syntax_writer', all blocks after are also recorded by the clone.
      if (i == kCopyAtBlockIndex) {
        ASSERT_WP2_OK(dicts_copy.CopyFrom(dicts));
        ASSERT_WP2_OK(syntax_writer_copy.CopyFrom(syntax_writer, &dicts_copy));
      }

      ASSERT_WP2_OK(syntax_writer.FindBestEncodingMethods(&cb));
      ASSERT_WP2_OK(syntax_writer.Record(cb));
      ASSERT_WP2_OK(syntax_writer.RecordSize(mgr, cb.dim()));

      if (i >= kCopyAtBlockIndex) {
        ASSERT_WP2_OK(syntax_writer_copy.FindBestEncodingMethods(&cb));
        ASSERT_WP2_OK(syntax_writer_copy.Record(cb));
        ASSERT_WP2_OK(syntax_writer_copy.RecordSize(mgr, cb.dim()));
      }

      ASSERT_TRUE(mgr.UseSize(block.dim(),
                              (uint32_t)cblocks.size() - 1, &block));
      mgr.Use(block);
    }

    ANSEnc enc;
    mgr.Clear();
    ASSERT_WP2_OK(syntax_writer.WriteHeader(&enc));
    ASSERT_WP2_OK(syntax_writer.WriteBlocks(cblocks, &mgr, &enc));

    cost = enc.GetCost();
    cost_dicts = enc.GetCost(dicts);
    num_tokens = enc.NumTokens();
    ASSERT_WP2_OK(enc.Assemble());
    ASSERT_TRUE(buffer.resize(enc.BufferSize()));
    std::copy(enc.Buffer(), enc.Buffer() + enc.BufferSize(), buffer.data());
  }

  // 'syntax_writer' does not exist anymore, compare the clone output now.
  ANSEnc enc_copy;
  mgr.Clear();
  ASSERT_WP2_OK(syntax_writer_copy.WriteHeader(&enc_copy));
  ASSERT_WP2_OK(syntax_writer_copy.WriteBlocks(cblocks, &mgr, &enc_copy));

  ASSERT_EQ(cost, enc_copy.GetCost());
  ASSERT_EQ(cost_dicts, enc_copy.GetCost(dicts_copy));
  ASSERT_EQ(num_tokens, enc_copy.NumTokens());
  ASSERT_WP2_OK(enc_copy.Assemble());
  ASSERT_EQ(buffer.size(), enc_copy.BufferSize());
  ASSERT_TRUE(std::equal(buffer.data(), buffer.data() + buffer.size(),
                         enc_copy.Buffer()));
}

//------------------------------------------------------------------------------

}  // namespace
}  // namespace WP2
