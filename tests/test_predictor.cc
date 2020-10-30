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

// Unit tests for predictors.

#include "examples/example_utils.h"
#include "include/helpers.h"
#include "src/common/lossy/block.h"
#include "src/common/lossy/predictor.h"
#include "src/dsp/math.h"
#include "src/utils/random.h"

namespace WP2 {
namespace {

using ::testing::ElementsAre;
using ::testing::Ge;

class CodedBlockForTest : public CodedBlock {
 public:
  void SetContext(Channel channel, std::vector<int16_t> context) {
    std::copy(context.begin(), context.begin() + kContextSize,
              context_[channel]);
  }

  const int16_t* GetContext(Channel channel, bool, bool) const override {
    return context_[channel];
  }

 private:
  int16_t context_[3][kContextSize];
};

// Test that given a context where all pixels have the same value, all
// regular predictors predict that value for all pixels (except for rounding
// errors). The only predictors that don't have this property are the zero
// predictor (always predicts zero) and CflPredictor. They are not tested here.

WP2Status init(int16_t min, int16_t max, YPredictors* const predictors) {
  return predictors->Fill(min, max);
}

template <class V>
void TestConstantContext(const std::vector<BlockSize>& block_sizes) {
  ContextCache context_cache;
  for (BlockSize size : block_sizes) {
    CodedBlock cb;
    cb.SetDimDefault({/*x=*/1, /*y=*/1, size});

    const int16_t kValue = 42;

    YUVPlane plane;
    ASSERT_WP2_OK(
        plane.Resize(cb.x_pix() + cb.w_pix(), cb.y_pix() + cb.h_pix()));
    plane.Y.Fill(kValue);
    cb.SetContextInput(plane, &context_cache);

    V predictors;
    ASSERT_WP2_OK(init(-512, 512, &predictors));
    for (const Predictor* p : predictors) {
      int16_t prediction[kMaxBlockSizePix2];
      p->Predict(cb, kYChannel, prediction, cb.w_pix());
      for (uint32_t i = 0; i < cb.w_pix() * cb.h_pix(); ++i) {
        EXPECT_NEAR(kValue, prediction[i], 1);
      }
    }
  }
}

TEST(AllPredictors, ConstantContext) {
  WP2MathInit();
  WP2::PredictionInit();

  TestConstantContext<YPredictors>({BLK_8x8, BLK_16x4, BLK_16x32});
}

template <class V>
void TestPredictionRange(const std::vector<BlockSize>& block_sizes) {
  UniformIntDistribution gen(53251);

  // Since this is a random test, we run it multiple times.
  ContextCache context_cache;
  constexpr uint32_t kNumIter = 5;
  for (uint32_t iter = 0; iter < kNumIter; ++iter) {
    for (BlockSize size : block_sizes) {
      CodedBlock cb;
      cb.SetDimDefault({/*x=*/1, /*y=*/1, size});

      const int16_t max = 1 << gen.Get(2, 10);
      const int16_t min = -(1 << gen.Get(2, 10));

      V predictors;
      ASSERT_WP2_OK(init(min, max, &predictors));

      YUVPlane plane;
      ASSERT_WP2_OK(
          plane.Resize(cb.x_pix() + cb.w_pix(), cb.y_pix() + cb.h_pix()));
      for (uint32_t j = 0; j < plane.Y.h_; ++j) {
        for (uint32_t i = 0; i < plane.Y.w_; ++i) {
          plane.Y.At(i, j) = gen.Get(min, max);
        }
      }
      cb.SetContextInput(plane, &context_cache);

      for (const Predictor* p : predictors) {
        int16_t prediction[kMaxBlockSizePix2];
        p->Predict(cb, kYChannel, prediction, cb.w_pix());

        for (uint32_t i = 0; i < cb.w_pix() * cb.h_pix(); ++i) {
          const int16_t v = prediction[i];
          EXPECT_GE(v, min) << "Failed for predictor " << p->GetName()
                            << " block size " << size;
          EXPECT_LE(v, max) << "Failed for predictor " << p->GetName()
                            << " block size " << size;
        }
      }
    }
  }
}

// Tests that predicted values are within the given min/max bounds.
TEST(AllPredictors, ValueRange) {
  WP2MathInit();
  WP2::PredictionInit();

  TestPredictionRange<YPredictors>({BLK_8x8, BLK_16x4, BLK_16x32});
}

TEST(ImplicitCflPredictor, SimpleLinearRelationship) {
  WP2MathInit();
  WP2::PredictionInit();

  CflPredictor p(/*min_value=*/-512, /*max_value=*/511);

  CodedBlockForTest cb;
  cb.SetDimDefault({/*x=*/0, /*y=*/0, /*dim=*/BLK_4x4});
  // Set up context so that U = 2 * Y + 50
  // Missing values should be ignored.
  const int16_t missing = CodedBlock::kMissing;
  cb.SetContext(kYChannel, {100, 100, 100, 100, 100,       // Left + top left
                            200, 200, 200, 200, 200,       // Top + top right
                            100, 100, missing, missing});  // Right
  cb.SetContext(kUChannel, {250, 250, 250, 250, 250,       // Left + top left
                            450, 450, 450, 450, 450,       // Top + top right
                            250, 250, missing, missing});  // Right

  YUVPlane plane;
  ASSERT_WP2_OK(plane.Resize(4, 4));
  plane.Y.Fill(50);
  ContextCache context_cache;
  cb.SetContextInput(plane, &context_cache);

  int16_t prediction[kPredSize];
  p.Predict(cb, kUChannel, prediction, kPredWidth);
  for (uint32_t i = 0; i < kPredSize; ++i) {
    // 2 * 50 + 50 = 150
    EXPECT_EQ(150, prediction[i]);
  }

  plane.Y.Fill(300);
  p.Predict(cb, kUChannel, prediction, kPredWidth);
  EXPECT_EQ(511, prediction[0]);  // Value is clipped.

  // Try with more extreme values to test accuracy of the prediction.
  cb.SetContext(kYChannel, {2, 2, 2, 2, 2,                 // Left + top left
                            4, 4, 4, 4, 4,                 // Top + top right
                            2, 2, missing, missing});      // Right
  cb.SetContext(kUChannel, {250, 250, 250, 250, 250,       // Left + top left
                            450, 450, 450, 450, 450,       // Top + top right
                            250, 250, missing, missing});  // Right
  plane.Y.Fill(4);
  cb.SetContextInput(plane);
  p.Predict(cb, kUChannel, prediction, kPredWidth);
  for (uint32_t i = 0; i < kPredSize; ++i) {
    // 100 * 4 + 50 = 450
    EXPECT_EQ(450, prediction[i]);
  }
}

TEST(ImplicitCflPredictor, AllMissing) {
  WP2MathInit();
  WP2::PredictionInit();

  CflPredictor p(/*min_value=*/-512, /*max_value=*/511);

  CodedBlockForTest cb;
  cb.SetDimDefault({/*x=*/0, /*y=*/0, /*dim=*/BLK_4x4});
  // No context available.
  const int16_t missing = CodedBlock::kMissing;
  cb.SetContext(
      kYChannel,
      {missing, missing, missing, missing, missing,  // Left + top left
       missing, missing, missing, missing, missing,  // Top + top right
       missing, missing, missing, missing});         // Right
  cb.SetContext(
      kUChannel,
      {missing, missing, missing, missing, missing,  // Left + top left
       missing, missing, missing, missing, missing,  // Top + top right
       missing, missing, missing, missing});         // Right

  YUVPlane plane;
  ASSERT_WP2_OK(plane.Resize(4, 4));
  plane.Y.Fill(50);
  ContextCache context_cache;
  cb.SetContextInput(plane, &context_cache);

  // Without context, the prediction should be 0.
  int16_t prediction[kPredSize];
  p.Predict(cb, kUChannel, prediction, kPredWidth);
  for (uint32_t i = 0; i < kPredSize; ++i) {
    EXPECT_EQ(0, prediction[i]);
  }
}

TEST(ImplicitCflPredictor, WithErrors) {
  WP2MathInit();
  WP2::PredictionInit();

  CflPredictor p(/*min_value=*/-512, /*max_value=*/511);

  CodedBlockForTest cb;
  cb.SetDimDefault({/*x=*/0, /*y=*/0, /*dim=*/BLK_4x4});
  // Set up context so that U is more or less equal to 2 * Y + 50, with some
  // errors.
  cb.SetContext(kYChannel, {100, 100, 100, 100, 100,  // Left + top left
                            200, 200, 200, 200, 200,  // Top + top right
                            100, 100, 100, 100});     // Right
  cb.SetContext(kUChannel,
                {250 + 10, 250 - 5, 250, 250 + 5, 250,        // Left + top left
                 450 + 11, 450 + 15, 450 - 10, 450 - 2, 450,  // Top + top right
                 250, 250 + 5, 250 - 3, 250 + 4});            // Right

  YUVPlane plane;
  ASSERT_WP2_OK(plane.Resize(4, 4));
  plane.Y.Fill(50);
  ContextCache context_cache;
  cb.SetContextInput(plane, &context_cache);

  int16_t prediction[kPredSize];
  p.Predict(cb, kUChannel, prediction, kPredWidth);
  for (uint32_t i = 0; i < kPredSize; ++i) {
    // 2 * 50 + 50 = 150, with some error.
    EXPECT_NEAR(150, prediction[i], 2);
  }
}

// Tests that if the context is constant, the prediction will be constant.
TEST(LargePredictors, UniformTest) {
  WP2::PredictionInit();
  UniformIntDistribution gen(53251);
  const uint32_t step = kMaxBlockSizePix + 5;
  const uint32_t kNumTests = 5;
  for (BlockSize dim : kAllBlockSizes) {
    const uint32_t w = BlockWidthPix(dim), h = BlockHeightPix(dim);
    const uint32_t ctx_size = ContextSize(w, h);
    int16_t ctx[kMaxContextSize];
    const int16_t max = 1 << WP2::kMaxYuvBits;
    const int16_t min = -max;
    for (uint32_t n = 0; n < kNumTests; ++n) {
      const int16_t r = gen.Get(min, max);
      for (uint32_t i = 0; i < ctx_size; ++i) ctx[i] = r;
      for (BasePredictor d :
           {BPRED_DC, BPRED_DC_L, BPRED_DC_T, BPRED_SMOOTH, BPRED_SMOOTH_H,
            BPRED_SMOOTH_V, BPRED_TM}) {
        int16_t prediction[kMaxBlockSizePix * step];
        BasePredictors[d](ctx, w, h, min, max, prediction, step);
        for (uint32_t j = 0; j < h; ++j) {
          for (uint32_t i = 0; i < w; ++i) {
            EXPECT_EQ(prediction[i + j * step], r);
          }
        }
      }
    }
    for (uint32_t ind = 0u; true; ++ind) {
      const float angle_deg = round(10.f * (12.5f + 12.5f / 3.f * ind)) / 10.f;
      if (angle_deg > 220.f) break;
      const int16_t r = gen.Get(min, max);
      for (uint32_t i = 0; i < ctx_size; ++i) ctx[i] = r;
      int16_t prediction[kMaxBlockSizePix * step];
      BaseAnglePredictor(angle_deg, ctx, w, h, prediction, step, min, max);
      for (uint32_t j = 0; j < h; ++j) {
        for (uint32_t i = 0; i < w; ++i) {
          ASSERT_NEAR(prediction[i + j * step], r, 1)
              << " at i=" << i << " j=" << j;
        }
      }
    }
    for (float strength : {0.0f, 0.1f, 0.5f, 1.0f}) {
      const int16_t r = gen.Get(min, max);
      for (uint32_t i = 0; i < ctx_size; ++i) ctx[i] = r;
      int16_t prediction[kMaxBlockSizePix * step];
      LargeWeightTable table;
      PrecomputeLargeWeightTable(strength, table);
      BaseFusePredictor(table, ctx, w, h, prediction, step, min, max);
      for (uint32_t j = 0; j < h; ++j) {
        for (uint32_t i = 0; i < w; ++i) {
          ASSERT_NEAR(prediction[i + j * step], r, 1)
              << " at i=" << i << " j=" << j;
        }
      }
    }
  }
}

TEST(LargePredictors, TestSimpleAngle) {
  WP2::PredictionInit();
  UniformIntDistribution gen(53251);
  const uint32_t step = kMaxBlockSizePix + 5;
  constexpr int16_t kMagic = 0x7324;
  int16_t prediction[kMaxBlockSizePix * step];
  for (uint32_t ind = 0; ind < kNumDirectionalAngles; ++ind) {
    for (BlockSize dim : kAllBlockSizes) {
      const uint32_t w = BlockWidthPix(dim), h = BlockHeightPix(dim);
      const uint32_t w_idx = TrfLog2[w];
      const uint32_t h_idx = TrfLog2[h];
      const uint32_t ctx_size =
          (ind < (int)kAngle_45) ? ContextWithTRSize(w, h) : ContextSize(w, h);
      int16_t ctx[kMaxContextSize];
      const int16_t max = 1 << WP2::kMaxYuvBits;
      const int16_t min = -max;
      const int16_t r = gen.Get(min, max);
      for (uint32_t i = 0; i < ctx_size; ++i) ctx[i] = r;
      for (auto& p : prediction) p = kMagic;
      SimpleAnglePredictor(ind, ctx, w_idx, h_idx, prediction, step);
      for (uint32_t j = 0; j < h; ++j) {
        for (uint32_t i = 0; i < w; ++i) {
          EXPECT_NEAR(prediction[i + j * step], r, 1)
              << " at i=" << i << " j=" << j << " index=" << ind;
        }
        for (uint32_t i = w; i < step; ++i) {
          EXPECT_EQ(prediction[i + j * step], kMagic);
        }
      }
    }
  }
}

TEST(LargePredictors, TestSimpleAnglePredictorSpeed) {
  WP2::PredictionInit();
  UniformIntDistribution gen(76234);
  constexpr int16_t max = 1 << WP2::kMaxYuvBits;
  constexpr int16_t min = -max;
  int16_t ctx[kMaxContextSize];
  for (auto& v : ctx) v = gen.Get(min, max);
  const uint32_t step = kMaxBlockSizePix + 7;
  int16_t prediction[kMaxBlockSizePix * step];
  constexpr uint32_t kNumTests = 1000;
  for (uint32_t ind = 0; ind < kNumDirectionalAngles; ++ind) {
    uint32_t CRC = 0;
    for (BlockSize dim : kAllBlockSizes) {
      const uint32_t w = BlockWidthPix(dim), h = BlockHeightPix(dim);
      const uint32_t w_idx = TrfLog2[w];
      const uint32_t h_idx = TrfLog2[h];
      for (uint32_t N = 0; N < kNumTests; ++N) {
        SimpleAnglePredictor(ind, ctx, w_idx, h_idx, prediction, step);
      }
      for (uint32_t j = 0; j < h; ++j) {
        for (uint32_t i = 0; i < w; ++i) {
          CRC = testing::SillyCRC(CRC, prediction[i + j * step]);
        }
      }
    }
    // It's a bit extreme to hardcode a CRC for each angle_idx, but it
    // will hopefully make the debugging easier for some particular angles
    // when writing SIMD code.
    constexpr uint32_t kCRCs[kNumDirectionalAngles] = {
        0xfc0766af, 0x24076af3, 0x7fa2a433, 0xf7faebf3, 0xe41835a5, 0x08683a98,
        0x63676ed7, 0xba68a9ac, 0x56daf1bd, 0x40b882b4, 0x592d6725, 0x6b74c1fe,
        0x2dee5147, 0x579abdb2, 0xcb1606a4, 0x407c836f, 0xb10777c1, 0x1564dfd5,
        0x62c17ea6, 0x1f50d632, 0xeeb93e3c, 0x42bf14c0, 0x80f43826, 0x4b67e1d2,
        0xa8d8ef27, 0xc067f6a2, 0x6b81ad22, 0xad969dc5, 0x43cfa60c, 0x28809de1,
        0xc272024e, 0xef4d3d45, 0xee4f256b, 0xe8d34800, 0xfc31e532, 0x9566491e,
        0x5e130153, 0x9006678f, 0x203b30c0, 0x9100f6ca, 0xd378e42f, 0x0f6dc9a5,
        0x0377ddd9, 0xebbe7dfb, 0x1dde713c, 0x210c9655, 0x901aac59, 0x04e2e256,
        0xb5b9fc5f, 0x31a52495, 0x5182df7d, 0xbe18889f, 0x841ce1fb, 0x09de7afe,
        0xcc13dc00, 0xfc14f13d
    };
    EXPECT_EQ(CRC, kCRCs[ind]) << "Angle idx : " << ind;
  }
}

TEST(LargePredictors, AlphaPredictors) {
  WP2::PredictionInit();
  APredictors preds;
  CSPTransform csp;
  csp.InitYCoCg();
  ASSERT_WP2_OK(preds.Fill(/*min_value=*/0, /*max_value=*/kAlphaMax, &csp));

  CodedBlock cb;
  cb.SetDimDefault({/*x=*/0, /*y=*/0, /*dim=*/BLK_4x4});
  ArgbBuffer rgb;
  ASSERT_WP2_OK(rgb.Resize(4, 4));
  rgb.Fill(/*window=*/{0, 0, 2, 2}, Argb32b{128, 32, 100, 54});
  rgb.Fill(/*window=*/{2, 0, 2, 2}, Argb32b{254, 230, 50, 0});
  rgb.Fill(/*window=*/{0, 2, 2, 2}, Argb32b{255, 255, 0, 0});
  rgb.Fill(/*window=*/{2, 2, 2, 2}, Argb32b{0, 0, 0, 0});
  YUVPlane yuv;
  ASSERT_WP2_OK(
      yuv.Import(rgb, /*import_alpha=*/false, csp, /*resize_if_needed=*/true));
  ContextCache context_cache;
  cb.SetContextInput(yuv, &context_cache);

  for (const Predictor* const p : preds) {
    SCOPED_TRACE(SPrintf("predictor %s", p->GetName().c_str()));
    int16_t prediction[16];
    p->Predict(cb, kAChannel, prediction, cb.w_pix());
    EXPECT_THAT(prediction,
                ElementsAre(Ge(100), Ge(100), Ge(230), Ge(230),  // Row 1.
                            Ge(100), Ge(100), Ge(230), Ge(230),  // Row 2.
                            Ge(255), Ge(255), Ge(0), Ge(0),      // Row 3.
                            Ge(255), Ge(255), Ge(0), Ge(0)));    // Row 4.
  }
}

// not exactly a test, but it's good to have a way to easily generate the tables
TEST(LargePredictors, PrintSimpleAngleTable) {
  constexpr uint32_t kAnglePredPrecision = 15;
  printf("static constexpr int32_t kCoTanTable[] = { /* %dbit precision */\n  ",
         kAnglePredPrecision);
  for (uint32_t i = 0; i < kNumDirectionalAngles; ++i) {
    const double alpha = M_PI / 180. * (45.f + 22.5f * ((int)i - 3) / 7.);
    printf("%d,", (i == 45) ? 0 :
           (int32_t)(1. / std::tan(alpha) * (1u << kAnglePredPrecision)));
    printf("%s", ((i & 7) == 7) ? "\n  " : " ");
  }
  printf("};\n");

  printf("static constexpr int32_t kTanTable[] = { /* %dbit precision */\n  ",
         kAnglePredPrecision);
  for (uint32_t i = 0; i < kNumDirectionalAngles; ++i) {
    const double alpha = M_PI / 180. * (45.f + 22.5f * ((int)i - 3) / 7.);
    printf("%d,", (i == 17) ? 0 :
           (int32_t)(std::tan(alpha) * (1u << kAnglePredPrecision)));
    printf("%s", ((i & 7) == 7) ? "\n  " : " ");
  }
  printf("};\n");
}

}  // namespace
}  // namespace WP2
