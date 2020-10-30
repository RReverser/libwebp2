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

// Test quantization matrix.

#include <memory>

#include "examples/example_utils.h"
#include "include/helpers.h"
#include "src/common/lossy/block_size.h"
#include "src/enc/wp2_enc_i.h"

namespace WP2 {
namespace {

class QuantMtxTest : public ::testing::Test {
  void SetUp() override {
    WP2DspReset();
    WP2TransformInit();
    WP2QuantizeInit();
  }
};

TEST_F(QuantMtxTest, RoundTrip) {
  const uint32_t kMaxResidual = 150;
  for (TrfSize trf : kAllTrfSizes) {
    const uint32_t num_coeffs = kNumCoeffs[trf];
    const uint32_t scale = num_coeffs >> 2;
    SCOPED_TRACE(SPrintf("trf: %s (scale %d)", kTDimNames[trf], scale));
    for (const uint16_t quant : {128, 493, 1499, 3021, 4523, 5000}) {
      SCOPED_TRACE(SPrintf("quant: %d", quant));
      QuantMtx mtx;
      ASSERT_WP2_OK(mtx.Allocate());
      const uint16_t quants[kNumQuantZones] = {quant, quant, quant, quant};
      mtx.Init(kMaxResidual, /*q_scale=*/1, ALL_RECTS, quants);
      int32_t residuals[kMaxBlockSizePix2];
      for (uint32_t i = 0; i < num_coeffs; ++i) {
        residuals[i] = i * (kMaxResidual / num_coeffs) * scale;
        if (i % 2 == 0) residuals[i] *= -1;
      }
      for (bool first_coeff_is_dc : {true, false}) {
        int16_t coeffs[kMaxBlockSizePix2];
        uint32_t num_coeffs_unused;
        const uint32_t error = mtx.RoundTrip(residuals, trf, first_coeff_is_dc,
                                             coeffs, &num_coeffs_unused);
        const float root_mean_square = std::sqrt(error / (float)num_coeffs);
        EXPECT_LT(root_mean_square, (1 + 2 * quant) * scale * 0.25f);
      }
    }
  }
}

TEST_F(QuantMtxTest, DCError) {
  for (TrfSize trf : kAllTrfSizes) {
    const uint32_t scale = kNumCoeffs[trf] >> 2;
    SCOPED_TRACE(SPrintf("trf: %s (scale %d)", kTDimNames[trf], scale));
    for (uint16_t quant : {128, 200, 500, 2001, 4000, 5000}) {
      SCOPED_TRACE(SPrintf("quant: %d", quant));
      QuantMtx mtx;
      ASSERT_WP2_OK(mtx.Allocate());
      const uint16_t quants[kNumQuantZones] = {quant, quant, quant, quant};
      mtx.Init(/*max_residual=*/150, /*q_scale=*/1, ALL_RECTS, quants);
      for (int32_t dc : {1000 / quant, -1000 / quant}) {
        const int32_t dc_error = mtx.DCError(dc, trf);
        EXPECT_LT(std::abs(dc_error), (1 + 2 * quant) * scale * 0.03f);
      }
    }
  }
}

TEST_F(QuantMtxTest, DCErrorRoundTrip) {
  UniformIntDistribution random(/*seed=*/123);

  const uint32_t kMaxResidual = 15000;
  for (TrfSize trf : kAllTrfSizes) {
    const uint32_t num_coeffs = kNumCoeffs[trf];
    const uint32_t scale = num_coeffs >> 2;
    SCOPED_TRACE(SPrintf("trf: %s (scale %d)", kTDimNames[trf], scale));
    for (uint16_t quant : {128, 200, 507, 2001, 4000, 5000}) {
      SCOPED_TRACE(SPrintf("quant: %d", quant));
      QuantMtx mtx;
      ASSERT_WP2_OK(mtx.Allocate());
      const uint16_t quants[kNumQuantZones] = {quant, quant, quant, quant};
      mtx.Init(kMaxResidual, /*q_scale=*/1, ALL_RECTS, quants);
      int32_t residuals[kMaxBlockSizePix2] = {0};
      int32_t range = kMaxResidual / num_coeffs * scale;
      residuals[0] = random.Get(-range, range);

      for (bool first_coeff_is_dc : {true, false}) {
        int16_t quantized[kMaxBlockSizePix2];
        uint32_t num_coeffs_unused;
        mtx.Quantize(residuals, trf, first_coeff_is_dc, quantized,
                     &num_coeffs_unused);
        int32_t dequantized[kMaxBlockSizePix2];
        mtx.Dequantize(quantized, num_coeffs, trf, dequantized);
        const int16_t dc_error = residuals[0] - dequantized[0];
        EXPECT_EQ(dc_error, mtx.DCError(residuals[0], trf));
      }
    }
  }
}

//------------------------------------------------------------------------------

class QuantMtxTestCpuInfo : public ::testing::TestWithParam<WP2CPUInfo> {
  void SetUp() override {
    WP2DspReset();
    WP2GetCPUInfo = GetParam();
    WP2TransformInit();
    WP2QuantizeInit();
  }
};

TEST_P(QuantMtxTestCpuInfo, VerifyDcScaling) {
  constexpr WP2TransformType kTrfs[] = { kDct, kAdst };
  for (const BlockSize block_size : kAllBlockSizes) {
    for (bool reduced : {false, false}) {
      const uint32_t W = BlockWidthPix(block_size);
      const uint32_t H = BlockHeightPix(block_size);
      const float expected_scale = sqrt(W * H) / (reduced ? 2. : 1.);
      for (auto& trf_x : kTrfs) {
        for (auto& trf_y : kTrfs) {
          const float ratio0 = (trf_x == kAdst && trf_y == kAdst) ? 1.2f : 1.f;
          for (int32_t dc = 50; dc < 1024; ++dc) {
            int32_t tmp[kMaxBlockSizePix2];
              for (auto& res : tmp) res = dc;
              WP2Transform2D(tmp, trf_x, trf_y, W, H, tmp, reduced);
              const float observed_scale = 1.f * tmp[0] / dc;
              EXPECT_NEAR(expected_scale / observed_scale, ratio0, 0.6)
                   << kDimNames[block_size] << "  DC=" << dc
                   << " transf=" << WP2TransformNames[trf_x]
                   << " x " << WP2TransformNames[trf_y];
          }
        }
      }
    }
  }
}

TEST_P(QuantMtxTestCpuInfo, DcQuantizedDynamicRange) {
  constexpr uint32_t kMaxCoeffs = kMaxBlockSizePix2;
  int32_t residuals[kMaxCoeffs];

  constexpr uint32_t kMaxResidual = 130;
  QuantMtx mtx;
  ASSERT_WP2_OK(mtx.Allocate());

  for (const BlockSize block_size : kAllBlockSizes) {
    std::fill(residuals, residuals + kMaxCoeffs, kMaxResidual);
    const uint32_t W = BlockWidthPix(block_size);
    const uint32_t H = BlockHeightPix(block_size);
    constexpr WP2TransformType kTrfs[] = { kDct, kAdst };
    for (auto& trf_x : kTrfs) {
      for (auto& trf_y : kTrfs) {
        WP2Transform2D(residuals, trf_x, trf_y, W, H,
                       residuals, /*reduced=*/false);
        SCOPED_TRACE(SPrintf("block size: %s  DC: %d  ratio=%.1f",
                             kDimNames[block_size], residuals[0], sqrt(W * H)));
        for (uint16_t quant : {129u, 149u, 230u, 1003u, 2017u,
                               4021u, 5722u, (uint32_t)kMaxQuantStep}) {
          SCOPED_TRACE(SPrintf("quant: %d", quant));
          const uint16_t quants[kNumQuantZones] = {quant, quant, quant, quant};
          mtx.Init(kMaxResidual, /*q_scale=*/1, ALL_RECTS, quants);
          int16_t coeffs[kMaxCoeffs];
          const TrfSize tdim = GetTransform(block_size);
          uint32_t num_coeffs;
          mtx.Quantize(residuals, tdim, /*first_is_dc=*/true, coeffs,
                       &num_coeffs);
          EXPECT_LE(num_coeffs, kNumCoeffs[tdim]);

          // Check DC is equal or near equal (because of rounding errors)
          // to its max range.
          if (trf_x == kDct && trf_y == kDct) {
            EXPECT_NEAR(std::abs(coeffs[0]), mtx.GetMaxAbsResDC(tdim), 2);
          }
          EXPECT_LE(std::abs(coeffs[0]), mtx.GetMaxAbsResDC(tdim));
        }
      }
    }
  }
}

TEST_P(QuantMtxTestCpuInfo, DequantizedDynamicRange) {
  constexpr uint32_t kMaxCoeffs = kMaxBlockSizePix2;
  int16_t coeffs[kMaxCoeffs];
  for (auto& v : coeffs) {
    v = kMaxCoeffValue;
    if (v % 2 == 0) {
      v *= -1;
    }
  }

  for (const TrfSize tdim : kAllTrfSizes) {
    const uint32_t num_coeffs = kNumCoeffs[tdim];
    SCOPED_TRACE(SPrintf("transform size: %d", tdim));

    for (uint16_t quant : {512u, 1200u, 4310u, (uint32_t)kMaxQuantStep}) {
      SCOPED_TRACE(SPrintf("quant: %d", quant));
      QuantMtx mtx;
      const uint16_t quants[kNumQuantZones] = { quant, quant, quant, quant };
      mtx.Init(/*max_residual=*/(1 << kMaxYuvBits), /*q_scale=*/1, ALL_RECTS,
               quants);

      int32_t residuals[kMaxCoeffs];
      mtx.Dequantize(coeffs, kNumCoeffs[tdim], tdim, residuals);

      const int32_t kMax = (1 << 18) * num_coeffs;
      for (uint32_t i = 0; i < num_coeffs; ++i) {
        EXPECT_LE(residuals[i], kMax);
        EXPECT_GE(residuals[i], -kMax);
      }
    }
  }
}

TEST_P(QuantMtxTestCpuInfo, TestQuantSpeed) {
  constexpr uint32_t kMaxCoeffs = kMaxBlockSizePix2;

  constexpr uint32_t kMaxResidual = 130;
  static constexpr uint32_t kBufferSize = 30000;
  int32_t buf[kBufferSize + kMaxBlockSizePix2];
  for (uint32_t i = 0; i < kBufferSize + kMaxBlockSizePix2; ++i) {
    buf[i] = (int32_t)(i % (2 * kMaxResidual + 1)) - kMaxResidual;
  }
  uint32_t bias[2 * kMaxCoeffs], iq[2 * kMaxCoeffs];
  int16_t dq[2 * kMaxCoeffs];
  for (uint32_t i = 0; i < 2 * kMaxCoeffs; ++i) bias[i] = i % (1 << WP2QBits);
  for (uint32_t i = 0; i < 2 * kMaxCoeffs; ++i) {
    const uint32_t Q = (i * 17 + 4) % kMaxQuantStep;
    iq[i] = (1u << WP2QBits) / (Q + 1);
    dq[i] = Q;
  }

  uint32_t crc1 = 335u, crc2 = 64221u;
  const uint32_t kLoop = 50;
  const uint32_t kNumTest = 32 * 32 * kLoop;

  for (const BlockSize block_size : kAllBlockSizes) {
    int32_t out0[kMaxCoeffs];
    int16_t out1[kMaxCoeffs];
    int32_t out2[kMaxCoeffs];
    for (uint32_t num = 1; num <= kNumTest; num += kLoop) {
      const TrfSize tdim = GetTransform(block_size);
      const uint32_t bsize = kNumCoeffs[tdim];
      uint32_t len = bsize % num;
      for (uint32_t i = 0; i < len; ++i) out0[i] = buf[(num % kBufferSize) + i];
      for (uint32_t i = len; i < bsize; ++i) out0[i] = 0;
      for (uint32_t k = 0; k < kLoop; ++k) {   // timing loop
        WP2Quantize(&iq[num % kMaxCoeffs], &bias[num % kMaxCoeffs],
                    out0, out1, bsize);
        WP2Dequantize(out1, &dq[num % kMaxCoeffs], out2, len, bsize);
      }
      for (uint32_t i = len; i < bsize; ++i) {
        EXPECT_EQ(out0[i], 0);
        EXPECT_EQ(out1[i], 0);
        EXPECT_EQ(out2[i], 0);
      }
      for (uint32_t i = 0; i < len; ++i) {
        crc1 = testing::SillyCRC(crc1, out1[i]);
        crc2 = testing::SillyCRC(crc2, out2[i]);
      }
    }
  }
  EXPECT_EQ(crc1, 4142269684u);
  EXPECT_EQ(crc2, 3318598408u);
}

INSTANTIATE_TEST_SUITE_P(QuantMtxTest1Instantiation, QuantMtxTestCpuInfo,
                         ::testing::ValuesIn(testing::kWP2CpuInfos));

//------------------------------------------------------------------------------

typedef ::testing::Types<FrontMgrLexico, FrontMgrMax> FrontMgrs;

template <class T>
class QuantMtxTestFrontMgrs : public ::testing::Test {
  void SetUp() override {
    WP2DspReset();
    WP2QuantizeInit();
    WP2MathInit();
    ANSInit();
  }
};

TYPED_TEST_SUITE(QuantMtxTestFrontMgrs, FrontMgrs);

// Checks that the quantizer outputs coeffs that can be successfully encoded
// by the syntax writer.
TYPED_TEST(QuantMtxTestFrontMgrs, PlaysWellWithSyntaxWriter) {
  // Unfortunately there's quite a bit of setup/boilerplate to make the
  // SyntaxWriter work.
  EncoderConfig config;
  GlobalParams params;
  config.partition_set = params.partition_set_ = ALL_RECTS;
  config.partition_snapping = params.partition_snapping_ = true;

  ANSDictionaries dicts;
  constexpr bool kUseAOMCoeffs = false, kHasAlpha = false;
  constexpr uint32_t kNumChannels = (kHasAlpha ? 4 : 3);
  SymbolsInfo symbols_info;
  ASSERT_WP2_OK(symbols_info.InitLossy(
      /*num_segments=*/1, ALL_RECTS, /*has_alpha=*/false,
      /*quality_hint=*/GetQualityHint(/*quality=*/75.f), kUseAOMCoeffs));
  ResidualWriter residual_writer;
  ASSERT_WP2_OK(residual_writer.Init(kUseAOMCoeffs, kHasAlpha));
  ANSEnc enc;
  std::unique_ptr<WP2::SymbolWriter>
      sw(new (WP2Allocable::nothrow) WP2::SymbolWriter);
  ASSERT_TRUE(sw != nullptr);

  ASSERT_WP2_OK(sw->Init(symbols_info));
  ASSERT_WP2_OK(sw->Allocate());

  const int16_t yuv_min = -512;
  const int16_t yuv_max = 511;
  ASSERT_WP2_OK(params.uv_preds_.Fill(yuv_min, yuv_max));
  ASSERT_WP2_OK(params.y_preds_.Fill(yuv_min, yuv_max));

  ASSERT_TRUE(params.segments_.resize(1));

  constexpr uint32_t kMaxCoeffs = kMaxBlockSizePix2;
  const uint32_t tf_i = 0;  // One single transform per block.
  int32_t residuals[kMaxCoeffs];
  ArgbBuffer buffer;

  // Try all block sizes.
  for (const BlockSize block_size : kAllBlockSizes) {
    SCOPED_TRACE(SPrintf("block size: %d", block_size));
    Vector<CodedBlock> cblocks;
    EXPECT_TRUE(cblocks.resize(1));
    CodedBlock &cb = cblocks[0];
    cb.is420_ = false;
    cb.id_ = 0;  // Segment id.

    const Rectangle rect = {0, 0, BlockWidthPix(block_size),
                            BlockHeightPix(block_size)};

    TypeParam mgr;
    ASSERT_WP2_OK(
        mgr.Init(ALL_RECTS, /*snapped=*/false, rect.width, rect.height));
    EXPECT_TRUE(mgr.Blocks().resize(1));
    cb.SetRange(yuv_min, yuv_max);
    cb.SetDim(/*block=*/{/*x=*/0, /*y=*/0, /*dim=*/block_size}, mgr);
    cb.y_context_is_constant_ = false;
    cb.SetPredictor(kYChannel, params.y_preds_[BPRED_DC_L]);
    cb.SetUVPredictor(params.uv_preds_[0]);

    ASSERT_WP2_OK(buffer.Resize(rect.width, rect.height));
    buffer.Fill({0xFFu, 0xFFu, 0xFFu, 0xFFu});
    YUVPlane yuv;
    ASSERT_WP2_OK(yuv.Import(buffer, buffer.HasTransparency(), params.transf_,
                             /*resize_if_needed=*/true));

    const TrfSize tdim = GetTransform(block_size);
    const uint32_t num_coeffs = kNumCoeffs[tdim];
    // Set the raw coeffs to the min/max possible values.
    for (uint32_t i = 0; i < kMaxCoeffs; ++i) {
      residuals[i] = (1 << kMaxYuvBits) * sqrt(num_coeffs);
      if (i % 2 == 0) {
        residuals[i] *= -1;
      }
    }

    // Try different quantization levels.
    for (uint16_t quant : {511, 1366, 2041, 3100, 5000}) {
      SCOPED_TRACE(SPrintf("quant: %d", quant));
      Vector<Segment>& segments = params.segments_;
      for (auto& s : segments) ASSERT_WP2_OK(s.AllocateForEncoder());

      // Quantize!
      const uint16_t quants[kNumQuantZones] = { quant, quant, quant, quant };
      segments[0].quant_y_.Init(/*max_residual=*/(1 << kMaxYuvBits),
                                /*q_scale=*/1, ALL_RECTS, quants);
      segments[0].quant_u_.Init(/*max_residual=*/(1 << kMaxYuvBits),
                                /*q_scale=*/1, ALL_RECTS, quants);
      segments[0].quant_v_.Init(/*max_residual=*/(1 << kMaxYuvBits),
                                /*q_scale=*/1, ALL_RECTS, quants);

      for (bool first_coeff_is_dc : {true, false}) {
        segments[0].quant_y_.Quantize(residuals, tdim, first_coeff_is_dc,
                                      cb.coeffs_[kYChannel][tf_i],
                                      &cb.num_coeffs_[kYChannel][tf_i]);
        segments[0].quant_u_.Quantize(residuals, tdim, first_coeff_is_dc,
                                      cb.coeffs_[kUChannel][tf_i],
                                      &cb.num_coeffs_[kUChannel][tf_i]);
        segments[0].quant_v_.Quantize(residuals, tdim, first_coeff_is_dc,
                                      cb.coeffs_[kVChannel][tf_i],
                                      &cb.num_coeffs_[kVChannel][tf_i]);

        SymbolRecorder recorder;
        ASSERT_WP2_OK(recorder.Allocate(symbols_info, /*num_records=*/0));
        Counters counters;
        ASSERT_WP2_OK(counters.Init(recorder));

        // Check the coeffs can be successfully written by SyntaxWriter.
        SyntaxWriter syntax_writer;
        ASSERT_WP2_OK(syntax_writer.Init(
            &dicts, config, params, yuv, ChromaSubsampling::k420, rect,
            (uint32_t)cblocks.size(), kUseAOMCoeffs));

        // Just check that we can get the pseudo rate without crashing.
        // WARNING! We must call this *after* syntax_writer.Init() because
        // syntax_writer.Init() will populate symbols_info (kSymbolDC) with
        // the max DC values from segment::GetMaxAbsDC
        const TrfSize dim = GetTransform(cb.dim());
        float pseudo_rate;
        ASSERT_WP2_OK(ResidualWriter::GetPseudoRate(
            kYChannel, kNumChannels, dim, cb.coeffs_[kYChannel][tf_i],
            cb.num_coeffs_[kYChannel][tf_i], cb.IsFirstCoeffDC(kYChannel),
            counters.residuals(), &pseudo_rate));

        ASSERT_WP2_OK(syntax_writer.FindBestEncodingMethods(&cb));
        ASSERT_WP2_OK(syntax_writer.Record(cb));
        ASSERT_WP2_OK(syntax_writer.RecordSize(mgr, cb.dim()));
        ASSERT_WP2_OK(syntax_writer.WriteHeader(&enc));
        ASSERT_WP2_OK(syntax_writer.WriteBlocks(cblocks, &mgr, &enc));
      }
    }
  }
}

}  // namespace
}  // namespace WP2
