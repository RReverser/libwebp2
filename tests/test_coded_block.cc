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

// Test CodedBlock.

#include "include/helpers.h"
#include "src/common/lossy/block.h"
#include "src/common/lossy/block_size.h"
#include "src/common/lossy/transforms.h"
#include "src/enc/wp2_enc_i.h"
#include "src/utils/plane.h"
#include "src/utils/random.h"

namespace WP2 {
namespace {

using ::testing::ElementsAre;

class ConstantPredictor : public Predictor {
 public:
  explicit ConstantPredictor(int16_t c) : c_(c) {}

  void Predict(const CodedBlock& cb, Channel channel, int16_t output[],
               uint32_t step) const override {
    for (uint32_t j = 0; j < kPredHeight; ++j) {
      for (uint32_t i = 0; i < kPredWidth; ++i) output[i] = c_;
      output += step;
    }
  }

  std::string GetName() const override { return ""; }
  std::string GetFakePredStr() const override { return ""; }
  std::string GetPredStr(const CodedBlock& cb, Channel channel) const override {
    return "";
  }
  void Draw(const WP2::Rectangle& rect,
            ArgbBuffer* const debug_output) const override {}

 private:
  int16_t c_;
};

struct FrontMgrLexicoTest {
  typedef FrontMgrLexico Type;
  static const uint32_t kId = 0;
};

struct FrontMgrMaxTest {
  typedef FrontMgrMax Type;
  static const uint32_t kId = 1;
};

typedef ::testing::Types<FrontMgrLexicoTest, FrontMgrMaxTest> FrontMgrs;

template <class T> class CodedBlockTest : public ::testing::Test {
  void SetUp() override {
    WP2DspReset();
    WP2QuantizeInit();
  }
};

TYPED_TEST_SUITE(CodedBlockTest, FrontMgrs);

TYPED_TEST(CodedBlockTest, FindBestModes) {
  WP2::ANSInit();
  WP2MathInit();
  WP2PSNRInit();
  WP2TransformInit();
  WP2::PredictionInit();

  constexpr uint32_t kHasAlpha = false;
  constexpr uint32_t kNumChannels = (kHasAlpha ? 4 : 3);
  const uint32_t height = kPredHeight * 2;
  const uint32_t width = kPredWidth * 2;
  const Rectangle tile_rect = {0, 0, width, height};
  typename TypeParam::Type mgr;
  EXPECT_WP2_OK(mgr.Init(ALL_RECTS, /*snapped=*/false, width, height));

  const int16_t v1 = 10;
  const int16_t v2 = 11;

  YUVPlane in;
  ASSERT_WP2_OK(in.Resize(width, height));
  in.Y.Fill(v1);
  in.Y.Fill({0, 0, kPredWidth, kPredHeight}, v2);
  YUVPlane out;
  ASSERT_WP2_OK(out.Resize(width, height));

  CodedBlock cb;
  cb.id_ = 0;
  cb.SetRange(/*yuv_min=*/-512, /*yuv_max=*/511);
  cb.SetDim(/*block=*/{/*x=*/0, /*y=*/0, /*dim=*/BLK_8x8}, mgr);
  cb.SetSrcInput(in);
  // No context is available so no SetContextInput().
  cb.SetReconstructedOutput(&out);
  ContextCache context_cache;
  cb.SetContextCache(&context_cache);
  Vector_u8 fake_map;
  ASSERT_TRUE(fake_map.resize(width * height));
  for (auto& m : fake_map) m = 0;

  YPredictors predictors;
  ASSERT_WP2_OK(predictors.Fill(/*min_value=*/-512, /*max_value=*/511));

  Segment segment;
  ASSERT_WP2_OK(segment.AllocateForEncoder());
  const uint16_t quants[kNumQuantZones] = {128, 128, 128, 128};
  segment.quant_y_.Init(/*max_residual=*/1024, /*q_scale=*/1, ALL_RECTS,
                        quants);
  // Force lambda value to make test less brittle.
  segment.quant_y_.SetLambda(1.f);

  SymbolsInfo symbols_info;
  ASSERT_WP2_OK(
      symbols_info.InitLossy(/*num_segments=*/1, ALL_RECTS, kHasAlpha,
                             GetQualityHint(EncoderConfig::kDefault.quality),
                             /*use_aom_coeffs=*/false));
  SymbolRecorder recorder;
  symbols_info.SetInfo(kSymbolModeY, /*min=*/0, /*max=*/kAPredModeNum - 1,
                       /*num_clusters=*/1, SymbolsInfo::StorageMethod::kAuto);
  ASSERT_WP2_OK(recorder.Allocate(symbols_info, /*num_records=*/0));
  Counters counters;
  ASSERT_WP2_OK(counters.Init(recorder));

  Vector<TransformPair> transforms;
  ASSERT_TRUE(transforms.push_back(kDctDct));

  BlockContext context;
  WP2_ASSERT_STATUS(context.Init(/*use_aom=*/false, width, height));
  ASSERT_WP2_OK(cb.FindBestPredTf(
      EncoderConfig::kDefault, tile_rect, predictors, segment, context,
      kYChannel, kNumChannels, /*reduced=*/false, transforms, &counters));
}

// Turns two coordinates into a number.
static int16_t Pix(uint32_t x, uint32_t y) {
  assert(y < 100);
  return x * 100 + y;
}

template <class T>
static void SetCoords4x4(uint32_t num_blocks, uint32_t x, uint32_t y,
                         T* const mgr, YUVPlane* const yuv,
                         CodedBlock* const cb) {
  for (uint32_t i = 0; i < num_blocks; ++i) {
    Block blk;
    ASSERT_TRUE(mgr->UseSize(BLK_4x4, /*ind=*/1u, &blk));
    mgr->Use(blk);
  }
  cb->SetDim(/*block=*/{x, y, /*dim=*/BLK_4x4}, *mgr);
  cb->SetSrcInput(*yuv);
  cb->SetContextInput(*yuv);
  cb->SetReconstructedOutput(yuv);
}

TYPED_TEST(CodedBlockTest, GetContext_FirstBlock) {
  const uint32_t kSize = 16;
  YUVPlane yuv;
  ASSERT_WP2_OK(yuv.Resize(kSize, kSize));
  for (uint32_t j = 0; j < kSize; ++j) {
    for (uint32_t i = 0; i < kSize; ++i) yuv.Y.At(i, j) = Pix(i, j);
  }
  yuv.U.Fill(0);
  yuv.V.Fill(0);

  typename TypeParam::Type mgr;
  EXPECT_WP2_OK(mgr.Init(ALL_RECTS, /*snapped=*/false, kSize, kSize));

  CodedBlock cb;
  cb.SetRange(/*yuv_min=*/-512, /*yuv_max=*/512);
  cb.SetDim(/*block=*/{/*x=*/0, /*y=*/0, /*dim=*/BLK_16x8}, mgr);
  cb.SetSrcInput(yuv);
  ContextCache context_cache;
  cb.SetContextInput(yuv, &context_cache);
  cb.SetReconstructedOutput(&yuv);

  const int16_t missing = CodedBlock::kMissing;

  for (bool extend_right : {false, true}) {
    int16_t truncated_context[34];
    const int16_t* full_context =
        cb.GetContext(kYChannel, /*fill_in=*/false, extend_right);
    std::copy(full_context, full_context + 34, truncated_context);
    EXPECT_THAT(truncated_context,
                ElementsAre(
                    // Left.
                    missing, missing, missing, missing, missing, missing,
                    missing, missing,
                    // Top left.
                    missing,
                    // Top.
                    missing, missing, missing, missing, missing, missing,
                    missing, missing, missing, missing, missing, missing,
                    missing, missing, missing, missing,
                    // Top right.
                    missing,
                    // Right.
                    missing, missing, missing, missing, missing, missing,
                    missing, missing));
    full_context = cb.GetContext(kYChannel, /*fill_in=*/true, extend_right);
    std::copy(full_context, full_context + 34, truncated_context);
    EXPECT_THAT(
        truncated_context,
        ElementsAre(0, 0, 0, 0, 0, 0, 0, 0,  // Left.
                    0,                       // Top left.
                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  // Top.
                    0,                         // Top right.
                    0, 0, 0, 0, 0, 0, 0, 0));  // Right.
  }

  // This is the first block of the image. The first (top left) sub-block has
  // therefore zero context.
  int16_t context[kContextSize];
  SetCoords4x4(/*num_blocks=*/0, /*x=*/0, /*y=*/0, &mgr, &yuv, &cb);
  const int16_t* context_ini =
      cb.GetContext(kYChannel, /*fill_in=*/false, /*extend_right=*/false);
  std::copy(context_ini, context_ini + kContextSize, context);
  EXPECT_THAT(context,
              ElementsAre(missing, missing, missing, missing,    // Left.
                          missing,                               // Top left.
                          missing, missing, missing, missing,    // Top.
                          missing,                               // Top right.
                          missing, missing, missing, missing));  // Right.
  context_ini =
      cb.GetContext(kYChannel, /*fill_in=*/true, /*extend_right=*/false);
  std::copy(context_ini, context_ini + kContextSize, context);
  EXPECT_THAT(context, ElementsAre(0, 0, 0, 0,    // Left.
                                   0,             // Top left.
                                   0, 0, 0, 0,    // Top.
                                   0,             // Top right.
                                   0, 0, 0, 0));  // Right.

  // Top middle sub block. Only left context is available.
  SetCoords4x4(/*num_blocks=*/1, /*x=*/1, /*y=*/0, &mgr, &yuv, &cb);
  context_ini =
      cb.GetContext(kYChannel, /*fill_in=*/false, /*extend_right=*/false);
  std::copy(context_ini, context_ini + kContextSize, context);
  EXPECT_THAT(context,
              ElementsAre(Pix(3, 3), Pix(3, 2), Pix(3, 1), Pix(3, 0),  // Left.
                          missing,                               // Top left.
                          missing, missing, missing, missing,    // Top.
                          missing,                               // Top right.
                          missing, missing, missing, missing));  // Right.
  context_ini =
      cb.GetContext(kYChannel, /*fill_in=*/true, /*extend_right=*/false);
  std::copy(context_ini, context_ini + kContextSize, context);
  EXPECT_THAT(
      context,
      ElementsAre(Pix(3, 3), Pix(3, 2), Pix(3, 1), Pix(3, 0),    // Left.
                  Pix(3, 0),                                     // Top left.
                  Pix(3, 0), Pix(3, 0), Pix(3, 0), Pix(3, 0),    // Top.
                  Pix(3, 0),                                     // Top right.
                  Pix(3, 0), Pix(3, 0), Pix(3, 0), Pix(3, 0)));  // Right.

  // Bottom left sub block. Only top context is available.
  SetCoords4x4(/*num_blocks=*/3, /*x=*/0, /*y=*/1, &mgr, &yuv, &cb);
  context_ini =
      cb.GetContext(kYChannel, /*fill_in=*/false, /*extend_right=*/false);
  std::copy(context_ini, context_ini + kContextSize, context);
  EXPECT_THAT(context,
              ElementsAre(missing, missing, missing, missing,  // Left.
                          missing,                             // Top left.
                          Pix(0, 3), Pix(1, 3), Pix(2, 3), Pix(3, 3),  // Top.
                          Pix(4, 3),                             // Top right.
                          missing, missing, missing, missing));  // Right.
  context_ini =
      cb.GetContext(kYChannel, /*fill_in=*/true, /*extend_right=*/false);
  std::copy(context_ini, context_ini + kContextSize, context);
  EXPECT_THAT(
      context,
      ElementsAre(Pix(0, 3), Pix(0, 3), Pix(0, 3), Pix(0, 3),    // Left.
                  Pix(0, 3),                                     // Top left.
                  Pix(0, 3), Pix(1, 3), Pix(2, 3), Pix(3, 3),    // Top.
                  Pix(4, 3),                                     // Top right.
                  Pix(4, 3), Pix(4, 3), Pix(4, 3), Pix(4, 3)));  // Right.

  // Bottom middle sub block. Left and top contexts are available.
  SetCoords4x4(/*num_blocks=*/1, /*x=*/1, /*y=*/1, &mgr, &yuv, &cb);
  context_ini =
      cb.GetContext(kYChannel, /*fill_in=*/false, /*extend_right=*/false);
  std::copy(context_ini, context_ini + kContextSize, context);
  EXPECT_THAT(context,
              ElementsAre(Pix(3, 7), Pix(3, 6), Pix(3, 5), Pix(3, 4),  // Left.
                          Pix(3, 3),  // Top left.
                          Pix(4, 3), Pix(5, 3), Pix(6, 3), Pix(7, 3),  // Top.
                          Pix(8, 3),                             // Top right.
                          missing, missing, missing, missing));  // Right.
  context_ini =
      cb.GetContext(kYChannel, /*fill_in=*/true, /*extend_right=*/false);
  std::copy(context_ini, context_ini + kContextSize, context);
  EXPECT_THAT(
      context,
      ElementsAre(Pix(3, 7), Pix(3, 6), Pix(3, 5), Pix(3, 4),    // Left.
                  Pix(3, 3),                                     // Top left.
                  Pix(4, 3), Pix(5, 3), Pix(6, 3), Pix(7, 3),    // Top.
                  Pix(8, 3),                                     // Top right.
                  Pix(8, 3), Pix(8, 3), Pix(8, 3), Pix(8, 3)));  // Right.
}

TYPED_TEST(CodedBlockTest, GetContext_HasRight) {
  const uint32_t kSize = 16;
  YUVPlane yuv;
  ASSERT_WP2_OK(yuv.Resize(kSize, kSize));
  for (uint32_t j = 0; j < kSize; ++j) {
    for (uint32_t i = 0; i < kSize; ++i) yuv.Y.At(i, j) = Pix(i, j);
  }

  // The image has 6 blocks, laid out like this:
  // 1 2 2 3
  // 1 4 4 3
  // 5 4 4 6
  // The block being tested is block number 4. Blocks 5 and 6 have not been
  // seen yet so their data is unavailable for context.
  typename TypeParam::Type mgr;
  EXPECT_WP2_OK(mgr.Init(ALL_RECTS, /*snapped=*/false, kSize, kSize));
  Block blk;
  ASSERT_TRUE(mgr.UseSize(BLK_4x8, /*ind=*/0u, &blk));
  mgr.Use(blk);
  ASSERT_TRUE(mgr.UseSize(BLK_8x4, /*ind=*/1u, &blk));
  mgr.Use(blk);

  // Test an imaginary 4x4 block at (0, 1).
  CodedBlock cb;
  cb.SetRange(/*yuv_min=*/-512, /*yuv_max=*/512);
  int8_t left, right, top_context_extent;
  cb.SetDim(/*block=*/{/*x=*/0, /*y=*/1, /*dim=*/BLK_4x4}, mgr);
  cb.GetOccupancy(&left, &right, &top_context_extent);
  EXPECT_EQ(left, -1);
  EXPECT_EQ(right, 0);
  EXPECT_EQ(top_context_extent, 2);

  ASSERT_TRUE(mgr.UseSize(BLK_4x8, /*ind=*/2u, &blk));
  mgr.Use(blk);

  cb.SetDim(/*block=*/{/*x=*/1, /*y=*/1, /*dim=*/BLK_8x8}, mgr);
  cb.GetOccupancy(&left, &right, &top_context_extent);
  EXPECT_EQ(left, 1);
  EXPECT_EQ(right, 1);
  EXPECT_EQ(top_context_extent, 1);

  cb.SetSrcInput(yuv);
  ContextCache context_cache;
  cb.SetContextInput(yuv, &context_cache);

  const int16_t missing = CodedBlock::kMissing;

  // Check full context.
  const int16_t* full_context =
      cb.GetContext(kYChannel, /*fill_in=*/false, /*extend_right=*/false);
  int16_t truncated_context[26];
  std::copy(full_context, full_context + 26, truncated_context);
  EXPECT_THAT(truncated_context,
              ElementsAre(
                  // Left.
                  missing, missing, missing, missing, Pix(3, 7), Pix(3, 6),
                  Pix(3, 5), Pix(3, 4),
                  // Top left
                  Pix(3, 3),
                  // Top
                  Pix(4, 3), Pix(5, 3), Pix(6, 3), Pix(7, 3), Pix(8, 3),
                  Pix(9, 3), Pix(10, 3), Pix(11, 3),
                  // Top rightS
                  Pix(12, 3),
                  // Right
                  Pix(12, 4), Pix(12, 5), Pix(12, 6), Pix(12, 7), missing,
                  missing, missing, missing));
  full_context =
      cb.GetContext(kYChannel, /*fill_in=*/false, /*extend_right=*/true);
  std::copy(full_context, full_context + 26, truncated_context);
  EXPECT_THAT(truncated_context,
              ElementsAre(
                  // Left.
                  missing, missing, missing, missing, Pix(3, 7), Pix(3, 6),
                  Pix(3, 5), Pix(3, 4),
                  // Top left
                  Pix(3, 3),
                  // Top
                  Pix(4, 3), Pix(5, 3), Pix(6, 3), Pix(7, 3), Pix(8, 3),
                  Pix(9, 3), Pix(10, 3), Pix(11, 3),
                  // Extended top right
                  Pix(12, 3), Pix(13, 3), Pix(14, 3), Pix(15, 3),
                  // Further
                  missing, missing, missing, missing, missing));
  full_context =
      cb.GetContext(kYChannel, /*fill_in=*/true, /*extend_right=*/false);
  std::copy(full_context, full_context + 26, truncated_context);
  EXPECT_THAT(truncated_context,
              ElementsAre(
                  // Left.
                  Pix(3, 7), Pix(3, 7), Pix(3, 7), Pix(3, 7), Pix(3, 7),
                  Pix(3, 6), Pix(3, 5), Pix(3, 4),
                  // Top left
                  Pix(3, 3),
                  // Top
                  Pix(4, 3), Pix(5, 3), Pix(6, 3), Pix(7, 3), Pix(8, 3),
                  Pix(9, 3), Pix(10, 3), Pix(11, 3),
                  // Top right
                  Pix(12, 3),
                  // Right
                  Pix(12, 4), Pix(12, 5), Pix(12, 6), Pix(12, 7), Pix(12, 7),
                  Pix(12, 7), Pix(12, 7), Pix(12, 7)));
  full_context =
      cb.GetContext(kYChannel, /*fill_in=*/true, /*extend_right=*/true);
  std::copy(full_context, full_context + 26, truncated_context);
  EXPECT_THAT(truncated_context,
              ElementsAre(
                  // Left.
                  Pix(3, 7), Pix(3, 7), Pix(3, 7), Pix(3, 7), Pix(3, 7),
                  Pix(3, 6), Pix(3, 5), Pix(3, 4),
                  // Top left
                  Pix(3, 3),
                  // Top
                  Pix(4, 3), Pix(5, 3), Pix(6, 3), Pix(7, 3), Pix(8, 3),
                  Pix(9, 3), Pix(10, 3), Pix(11, 3),
                  // Extended top right
                  Pix(12, 3), Pix(13, 3), Pix(14, 3), Pix(15, 3),
                  // Further
                  Pix(15, 3), Pix(15, 3), Pix(15, 3), Pix(15, 3), Pix(15, 3)));

  // Top left block. Left and top contexts are available
  // (from previously decoded blocks 1 and 2).
  int16_t context[kContextSize];
  SetCoords4x4(/*num_blocks=*/0, /*x=*/1, /*y=*/1, &mgr, &yuv, &cb);
  const int16_t* context_ini =
      cb.GetContext(kYChannel, /*fill_in=*/false, /*extend_right=*/false);
  std::copy(context_ini, context_ini + kContextSize, context);
  EXPECT_THAT(context,
              ElementsAre(Pix(3, 7), Pix(3, 6), Pix(3, 5), Pix(3, 4),  // Left.
                          Pix(3, 3),  // Top left.
                          Pix(4, 3), Pix(5, 3), Pix(6, 3), Pix(7, 3),  // Top.
                          Pix(8, 3),                             // Top right.
                          missing, missing, missing, missing));  // Right.
  context_ini =
      cb.GetContext(kYChannel, /*fill_in=*/true, /*extend_right=*/false);
  std::copy(context_ini, context_ini + kContextSize, context);
  EXPECT_THAT(
      context,
      ElementsAre(Pix(3, 7), Pix(3, 6), Pix(3, 5), Pix(3, 4),    // Left.
                  Pix(3, 3),                                     // Top left.
                  Pix(4, 3), Pix(5, 3), Pix(6, 3), Pix(7, 3),    // Top.
                  Pix(8, 3),                                     // Top right.
                  Pix(8, 3), Pix(8, 3), Pix(8, 3), Pix(8, 3)));  // Right.

  // Top right sub block. Full context (left, top and right) is available.
  SetCoords4x4(/*num_blocks=*/1, /*x=*/2, /*y=*/1, &mgr, &yuv, &cb);
  context_ini =
      cb.GetContext(kYChannel, /*fill_in=*/false, /*extend_right=*/false);
  std::copy(context_ini, context_ini + kContextSize, context);
  if (TypeParam::kId == 0) {
    EXPECT_THAT(
        context,
        ElementsAre(Pix(7, 7), Pix(7, 6), Pix(7, 5), Pix(7, 4),    // Left.
                    Pix(7, 3),                                     // Top left.
                    Pix(8, 3), Pix(9, 3), Pix(10, 3), Pix(11, 3),  // Top.
                    Pix(12, 3),                                    // Top right.
                    Pix(12, 4), Pix(12, 5), Pix(12, 6), Pix(12, 7)));  // Right.
    context_ini =
        cb.GetContext(kYChannel, /*fill_in=*/true, /*extend_right=*/false);
    std::copy(context_ini, context_ini + kContextSize, context);
    EXPECT_THAT(
        context,
        ElementsAre(Pix(7, 7), Pix(7, 6), Pix(7, 5), Pix(7, 4),    // Left.
                    Pix(7, 3),                                     // Top left.
                    Pix(8, 3), Pix(9, 3), Pix(10, 3), Pix(11, 3),  // Top.
                    Pix(12, 3),                                    // Top right.
                    Pix(12, 4), Pix(12, 5), Pix(12, 6), Pix(12, 7)));  // Right.
  } else {
    assert(TypeParam::kId == 1);
    EXPECT_THAT(
        context,
        ElementsAre(missing, missing, missing, missing,            // Left.
                    Pix(7, 3),                                     // Top left.
                    Pix(8, 3), Pix(9, 3), Pix(10, 3), Pix(11, 3),  // Top.
                    Pix(12, 3),                                    // Top right.
                    Pix(12, 4), Pix(12, 5), Pix(12, 6), Pix(12, 7)));  // Right.
    context_ini =
        cb.GetContext(kYChannel, /*fill_in=*/true, /*extend_right=*/false);
    std::copy(context_ini, context_ini + kContextSize, context);
    EXPECT_THAT(
        context,
        ElementsAre(Pix(7, 3), Pix(7, 3), Pix(7, 3), Pix(7, 3),    // Left.
                    Pix(7, 3),                                     // Top left.
                    Pix(8, 3), Pix(9, 3), Pix(10, 3), Pix(11, 3),  // Top.
                    Pix(12, 3),                                    // Top right.
                    Pix(12, 4), Pix(12, 5), Pix(12, 6), Pix(12, 7)));  // Right.
  }

  // Bottom left sub block. Only top context is available.
  SetCoords4x4(/*num_blocks=*/1, /*x=*/1, /*y=*/2, &mgr, &yuv, &cb);
  context_ini =
      cb.GetContext(kYChannel, /*fill_in=*/false, /*extend_right=*/false);
  std::copy(context_ini, context_ini + kContextSize, context);
  EXPECT_THAT(context,
              ElementsAre(missing, missing, missing, missing,  // Left.
                          Pix(3, 7),                           // Top left.
                          Pix(4, 7), Pix(5, 7), Pix(6, 7), Pix(7, 7),  // Top.
                          Pix(8, 7),                             // Top right.
                          missing, missing, missing, missing));  // Right.
  context_ini =
      cb.GetContext(kYChannel, /*fill_in=*/true, /*extend_right=*/false);
  std::copy(context_ini, context_ini + kContextSize, context);
  EXPECT_THAT(
      context,
      ElementsAre(Pix(3, 7), Pix(3, 7), Pix(3, 7), Pix(3, 7),    // Left.
                  Pix(3, 7),                                     // Top left.
                  Pix(4, 7), Pix(5, 7), Pix(6, 7), Pix(7, 7),    // Top.
                  Pix(8, 7),                                     // Top right.
                  Pix(8, 7), Pix(8, 7), Pix(8, 7), Pix(8, 7)));  // Right.

  // Bottom right sub block. Left and top contexts are available.
  SetCoords4x4(/*num_blocks=*/2, /*x=*/2, /*y=*/2, &mgr, &yuv, &cb);
  context_ini =
      cb.GetContext(kYChannel, /*fill_in=*/false, /*extend_right=*/false);
  std::copy(context_ini, context_ini + kContextSize, context);
  EXPECT_THAT(
      context,
      ElementsAre(Pix(7, 11), Pix(7, 10), Pix(7, 9), Pix(7, 8),  // Left.
                  Pix(7, 7),                                     // Top left.
                  Pix(8, 7), Pix(9, 7), Pix(10, 7), Pix(11, 7),  // Top.
                  Pix(12, 7),                                    // Top right.
                  missing, missing, missing, missing));          // Right.
  context_ini =
      cb.GetContext(kYChannel, /*fill_in=*/true, /*extend_right=*/false);
  std::copy(context_ini, context_ini + kContextSize, context);
  EXPECT_THAT(
      context,
      ElementsAre(Pix(7, 11), Pix(7, 10), Pix(7, 9), Pix(7, 8),  // Left.
                  Pix(7, 7),                                     // Top left.
                  Pix(8, 7), Pix(9, 7), Pix(10, 7), Pix(11, 7),  // Top.
                  Pix(12, 7),                                    // Top right.
                  Pix(12, 7), Pix(12, 7), Pix(12, 7), Pix(12, 7)));  // Right.
}

TYPED_TEST(CodedBlockTest, GetContext_Large) {
  // The grid is 16 by 16 and we only add 8x8 blocks. We will number them:
  //  0 1
  //  2 3
  const uint32_t kSize = 16;
  YUVPlane yuv;
  ASSERT_WP2_OK(yuv.Resize(kSize, kSize));
  for (uint32_t j = 0; j < kSize; ++j) {
    for (uint32_t i = 0; i < kSize; ++i) yuv.Y.At(i, j) = Pix(i, j);
  }

  typename TypeParam::Type mgr;
  EXPECT_WP2_OK(mgr.Init(ALL_RECTS, /*snapped=*/false, kSize, kSize));
  // Add block 0.
  ASSERT_TRUE(mgr.UseSize(BLK_8x8, /*ind=*/0u, /*block=*/nullptr));

  // Create block 2.
  CodedBlock cb;
  cb.SetRange(/*yuv_min=*/-512, /*yuv_max=*/512);
  cb.SetDim(/*block=*/{/*x=*/0, /*y=*/2, /*dim=*/BLK_8x8}, mgr);
  cb.SetSrcInput(yuv);
  ContextCache context_cache;
  cb.SetContextInput(yuv, &context_cache);
  int8_t left, right, top_context_extent;
  cb.GetOccupancy(&left, &right, &top_context_extent);
  EXPECT_EQ(left, -1);
  EXPECT_EQ(right, -1);
  EXPECT_EQ(top_context_extent, 0);

  const int16_t missing = CodedBlock::kMissing;

  // The context of block 2 is only the bottom of block 0.
  const int16_t* full_context;
  int16_t context[8 + 10 + 8];
  for (bool extend_right : {false, true}) {
    full_context = cb.GetContext(kYChannel, /*fill_in=*/false, extend_right);
    std::copy(full_context, full_context + 26, context);
    EXPECT_THAT(context, ElementsAre(
                             // Left.
                             missing, missing, missing, missing, missing,
                             missing, missing, missing, missing,
                             // Top.
                             Pix(0, 7), Pix(1, 7), Pix(2, 7), Pix(3, 7),
                             Pix(4, 7), Pix(5, 7), Pix(6, 7), Pix(7, 7),
                             // Right.
                             missing, missing, missing, missing, missing,
                             missing, missing, missing, missing));
    full_context = cb.GetContext(kYChannel, /*fill_in=*/true, extend_right);
    std::copy(full_context, full_context + 26, context);
    EXPECT_THAT(context,
                ElementsAre(
                    // Left.
                    Pix(0, 7), Pix(0, 7), Pix(0, 7), Pix(0, 7), Pix(0, 7),
                    Pix(0, 7), Pix(0, 7), Pix(0, 7), Pix(0, 7),
                    // Top.
                    Pix(0, 7), Pix(1, 7), Pix(2, 7), Pix(3, 7), Pix(4, 7),
                    Pix(5, 7), Pix(6, 7), Pix(7, 7),
                    // Right.
                    Pix(7, 7), Pix(7, 7), Pix(7, 7), Pix(7, 7), Pix(7, 7),
                    Pix(7, 7), Pix(7, 7), Pix(7, 7), Pix(7, 7)));
  }

  // Add block 1 and redefine block 2.
  Block blk;
  ASSERT_TRUE(mgr.UseSize(BLK_8x8, /*ind=*/1u, &blk));
  mgr.Use(blk);
  cb.SetDim({/*x=*/0, /*y=*/2, /*dim=*/BLK_8x8}, mgr);
  cb.ResetContextCache();
  // The context of block 2 now has a top right context.
  full_context =
      cb.GetContext(kYChannel, /*fill_in=*/false, /*extend_right=*/false);
  std::copy(full_context, full_context + 26, context);
  EXPECT_THAT(context, ElementsAre(
                           // Left.
                           missing, missing, missing, missing, missing, missing,
                           missing, missing, missing,
                           // Top.
                           Pix(0, 7), Pix(1, 7), Pix(2, 7), Pix(3, 7),
                           Pix(4, 7), Pix(5, 7), Pix(6, 7), Pix(7, 7),
                           // Right.
                           Pix(8, 7), missing, missing, missing, missing,
                           missing, missing, missing, missing));
  full_context =
      cb.GetContext(kYChannel, /*fill_in=*/false, /*extend_right=*/true);
  std::copy(full_context, full_context + 26, context);
  EXPECT_THAT(context,
              ElementsAre(
                  // Left.
                  missing, missing, missing, missing, missing, missing, missing,
                  missing, missing,
                  // Top.
                  Pix(0, 7), Pix(1, 7), Pix(2, 7), Pix(3, 7), Pix(4, 7),
                  Pix(5, 7), Pix(6, 7), Pix(7, 7),
                  // Top right.
                  Pix(8, 7), Pix(9, 7), Pix(10, 7), Pix(11, 7), Pix(12, 7),
                  Pix(13, 7), Pix(14, 7), Pix(15, 7), missing));
  full_context =
      cb.GetContext(kYChannel, /*fill_in=*/true, /*extend_right=*/false);
  std::copy(full_context, full_context + 26, context);
  EXPECT_THAT(context,
              ElementsAre(
                  // Left.
                  Pix(0, 7), Pix(0, 7), Pix(0, 7), Pix(0, 7), Pix(0, 7),
                  Pix(0, 7), Pix(0, 7), Pix(0, 7), Pix(0, 7),
                  // Top.
                  Pix(0, 7), Pix(1, 7), Pix(2, 7), Pix(3, 7), Pix(4, 7),
                  Pix(5, 7), Pix(6, 7), Pix(7, 7),
                  // Right.
                  Pix(8, 7), Pix(8, 7), Pix(8, 7), Pix(8, 7), Pix(8, 7),
                  Pix(8, 7), Pix(8, 7), Pix(8, 7), Pix(8, 7)));
  full_context =
      cb.GetContext(kYChannel, /*fill_in=*/true, /*extend_right=*/true);
  std::copy(full_context, full_context + 26, context);
  EXPECT_THAT(context,
              ElementsAre(
                  // Left.
                  Pix(0, 7), Pix(0, 7), Pix(0, 7), Pix(0, 7), Pix(0, 7),
                  Pix(0, 7), Pix(0, 7), Pix(0, 7), Pix(0, 7),
                  // Top.
                  Pix(0, 7), Pix(1, 7), Pix(2, 7), Pix(3, 7), Pix(4, 7),
                  Pix(5, 7), Pix(6, 7), Pix(7, 7),
                  // Top right.
                  Pix(8, 7), Pix(9, 7), Pix(10, 7), Pix(11, 7), Pix(12, 7),
                  Pix(13, 7), Pix(14, 7), Pix(15, 7), Pix(15, 7)));
}

TYPED_TEST(CodedBlockTest, ContextIsConstant) {
  const uint32_t kSize = 32;
  YUVPlane yuv;
  ASSERT_WP2_OK(yuv.Resize(kSize, kSize));
  // Y plane is constant.
  yuv.Y.Fill(42);
  // U plane has a bunch of different values.
  for (uint32_t j = 0; j < kSize; ++j) {
    for (uint32_t i = 0; i < kSize; ++i) yuv.U.At(i, j) = Pix(i, j);
  }

  typename TypeParam::Type mgr;
  EXPECT_WP2_OK(mgr.Init(ALL_RECTS, /*snapped=*/false, kSize, kSize));
  Block blk;
  ASSERT_TRUE(mgr.UseSize(BLK_32x8, /*ind=*/0u, &blk));
  mgr.Use(blk);
  ASSERT_TRUE(mgr.UseSize(BLK_4x4, /*ind=*/1u, &blk));
  mgr.Use(blk);

  CodedBlock cb;
  cb.SetRange(/*yuv_min=*/-512, /*yuv_max=*/512);
  cb.SetDim(/*block=*/{/*x=*/1, /*y=*/2, /*dim=*/BLK_16x8}, mgr);
  cb.SetSrcInput(yuv);
  ContextCache context_cache;
  cb.SetContextInput(yuv, &context_cache);

  EXPECT_TRUE(cb.ContextIsConstant(kYChannel));
  EXPECT_FALSE(cb.ContextIsConstant(kUChannel));

  // Change a pixel inside the block.
  yuv.Y.At(4, 8) = 666;
  EXPECT_TRUE(cb.ContextIsConstant(kYChannel));
  // Change a pixel outside the block (in the context).
  yuv.Y.At(3, 7) = 666;
  cb.ResetContextCache();
  EXPECT_FALSE(cb.ContextIsConstant(kYChannel));
}

//------------------------------------------------------------------------------

// Make sure that the computed distortion does not vary based on partitioning.
TEST(GetDistoTest, BlockSizeDoesNotMatter) {
  WP2PSNRInit();
  constexpr uint32_t kNumSubX = 4, kNumSubY = 2;
  CodedBlock cb;
  cb.SetDimDefault(Block(0, 0, BLK_32x32));
  const uint32_t sub_w = cb.w() / kNumSubX, sub_h = cb.h() / kNumSubY;

  YUVPlane original, decoded;
  ASSERT_WP2_OK(
      original.Resize(cb.w_pix(), cb.h_pix(), /*pad=*/1, /*has_alpha=*/true));
  original.Fill(Ayuv38b{200, 321, 123, -10});
  ASSERT_WP2_OK(decoded.Copy(original, /*resize_if_needed=*/true));

  UniformIntDistribution random(/*seed==*/666);
  for (Channel channel : {kYChannel, kUChannel, kVChannel, kAChannel}) {
    for (uint32_t y = 0; y < original.GetHeight(); ++y) {
      for (uint32_t x = 0; x < original.GetWidth(); ++x) {
        decoded.GetChannel(channel).At(x, y) += random.Get(-37, 50);
      }
    }

    uint32_t total_sub_disto = 0;
    for (uint32_t y = 0; y < kNumSubY; ++y) {
      for (uint32_t x = 0; x < kNumSubX; ++x) {
        CodedBlock sub_cb;
        sub_cb.SetDimDefault(
            Block(x * sub_w, y * sub_h, GetBlockSize(sub_w, sub_h)));
        sub_cb.SetSrcInput(original);  // Must be called after SetDimDefault().
        sub_cb.SetReconstructedOutput(&decoded);
        total_sub_disto += sub_cb.GetDisto(channel, cb.blk().rect_pix());
      }
    }
    cb.SetSrcInput(original);
    cb.SetReconstructedOutput(&decoded);
    EXPECT_EQ(total_sub_disto, cb.GetDisto(channel, cb.blk().rect_pix()));
  }
}

}  // namespace
}  // namespace WP2
