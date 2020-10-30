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

#include <memory>

#include "include/helpers.h"
#include "src/common/lossy/block.h"
#include "src/enc/wp2_enc_i.h"
#include "src/dec/residuals_dec_aom.h"
#include "src/enc/residuals_enc_aom.h"

namespace WP2 {
namespace {

//------------------------------------------------------------------------------

using AOMTuple = std::tuple<BlockSize, uint32_t>;
class ResidualWriterTest
    : public ::testing::TestWithParam<std::tuple<AOMTuple, float>> {};

TEST_P(ResidualWriterTest, AOM) {
  WP2MathInit();
  ANSInit();
  UniformIntDistribution random(/*seed=*/0);
  FrontMgrLexico mgr;
  EXPECT_WP2_OK(
      mgr.Init(ALL_RECTS, /*snapped=*/false, /*width=*/10, /*height=*/10));
  libgav1::AOMResidualWriter writer;
  libgav1::AOMResidualReader reader;
  libgav1::AOMContext aom_context_enc, aom_context_dec;
  const uint32_t quality_hint = GetQualityHint(std::get<1>(GetParam()));
  const AOMTuple aom_tuple = std::get<0>(GetParam());
  const uint32_t num_iters = std::get<1>(aom_tuple);
  const BlockSize dim = std::get<0>(aom_tuple);
  SymbolsInfo info;
  ASSERT_WP2_OK(info.InitLossy(/*num_segments=*/1, ALL_RECTS,
                               /*has_alpha=*/false, quality_hint,
                               /*use_aom_coeffs=*/true));
  // Create some data.
  const uint32_t bw = BlockWidth[dim];
  const uint32_t bh = BlockHeight[dim];
  const uint32_t bw_pix = BlockWidthPix(dim);
  const uint32_t bh_pix = BlockHeightPix(dim);
  const uint32_t num_elem = bw_pix * bh_pix;
  const uint32_t tf_i = 0;  // One single transform per block.

  CodedBlock cb;
  Block block(/*x=*/0, /*y=*/0, dim);
  cb.SetDim(block, mgr);
  cb.is420_ = false;
  writer.Init();
  reader.Init();
  EXPECT_WP2_OK(aom_context_enc.Init(bw, bh));
  EXPECT_WP2_OK(aom_context_dec.Init(bw, bh));
  for (bool do_reset : {false, true}) {
    std::unique_ptr<WP2::SymbolWriter>
      sw(new (WP2Allocable::nothrow) WP2::SymbolWriter);
    ASSERT_TRUE(sw != nullptr);
    ASSERT_WP2_OK(sw->Init(info));
    ASSERT_WP2_OK(sw->Allocate());
    SymbolRecorder recorder;
    ASSERT_WP2_OK(recorder.Allocate(info, /*num_records=*/num_iters));
    for (uint32_t num = 0; num < num_iters; ++num) {
      for (Channel channel : {kYChannel, kUChannel, kVChannel}) {
        for (bool fill_all_zero : {false, true}) {
          // Create some data.
          int16_t* const coeffs = cb.coeffs_[channel][tf_i];
          if (fill_all_zero) {
            std::fill(coeffs, coeffs + num_elem, 0);
          } else {
            for (uint32_t y = 0, i = 0; y < bh_pix; ++y) {
              for (uint32_t x = 0; x < bw_pix; ++x, ++i) {
                if (i == 0) {
                  // DC.
                  coeffs[0] = random.Get(-20, 20) *
                              std::exp(-4. * random.Get(0, 20) / 20);
                } else {
                  coeffs[i] = (4. - 4.f * x / (bw_pix - 1)) *
                              (4. - 4.f * y / (bh_pix - 1)) *
                              std::exp(-4. * random.Get(0, 20) / 20);
                }
                if (random.FlipACoin()) coeffs[i] = -coeffs[i];
              }
            }
          }
          // Record data to get stats.
          if (do_reset) {
            recorder.ResetRecord(/*reset_backup=*/true);
            writer.Init();
            aom_context_enc.Reset();
          }
          ANSEncNoop enc_noop;
          EXPECT_WP2_OK(
              writer.WriteCoeffs(cb.x_pix(), cb.y_pix(), cb.dim(), kDctDct,
                                 /*first_is_dc=*/true, cb.is420_, coeffs,
                                 channel, /*do_update=*/false, &recorder,
                                 &enc_noop, &aom_context_enc));

          // Write the data.
          ANSEnc enc;
          if (do_reset) {
            writer.Init();
            aom_context_enc.Reset();
          }
          // Write the headers.
          const uint32_t max_nnz = 32 * 32;
          ANSDictionaries dicts;
          // Write the headers.
          for (uint32_t s = 0; s < info.Size(); ++s) {
            for (uint32_t c = 0; c < info.NumClusters(s); ++c) {
              ASSERT_WP2_OK(sw->WriteHeader(c, max_nnz, s, recorder, "counts",
                                            &enc, &dicts));
            }
          }
          // Write the coefficients.
          EXPECT_WP2_OK(
              writer.WriteCoeffs(cb.x_pix(), cb.y_pix(), cb.dim(), kDctDct,
                                 /*first_is_dc=*/true, cb.is420_, coeffs,
                                 channel, /*do_update=*/true, sw.get(),
                                 &enc, &aom_context_enc));

          // Read the data.
          EXPECT_WP2_OK(enc.Assemble());
          const size_t size = enc.BufferSize();
          const uint8_t* buf = enc.Buffer();
          // Initialize decoder and symbol reader.
          ExternalDataSource data_source(buf, size);
          ANSDec dec(&data_source);
          SymbolReader sr;
          ASSERT_WP2_OK(sr.Init(info, &dec));
          ASSERT_WP2_OK(sr.Allocate());
          for (uint32_t s = 0; s < info.Size(); ++s) {
            for (uint32_t c = 0; c < info.NumClusters(s); ++c) {
              ASSERT_WP2_OK(sr.ReadHeader(c, max_nnz, s, "counts"));
            }
          }

          if (do_reset) {
            reader.Init();
            aom_context_dec.Reset();
          }
          CodedBlock cb_out;
          cb_out.SetDim(block, mgr);
          cb_out.is420_ = cb.is420_;
          // Rand might generate zeros, checking 'fill_all_zero' is not enough.
          const bool all_zero =
              (std::count(coeffs, coeffs + num_elem, 0) == num_elem);
          // This is usually decoded before by ReadHasCoeffs().
          cb_out.num_coeffs_[channel][tf_i] = (all_zero ? 0 : num_elem);
          std::fill(cb_out.coeffs_[channel][tf_i],
                    cb_out.coeffs_[channel][tf_i] + num_elem, 0);
          const int max_num_coeffs =
              all_zero ? 0 : WP2::kNumCoeffs[cb.tdim(channel)];
          EXPECT_WP2_OK(reader.ReadCoeffs(
              channel, cb.x_pix(), cb.y_pix(), cb.dim(), cb.is420_,
              cb_out.GetCodingParams(channel)->tf,
              cb_out.IsFirstCoeffDC(channel), max_num_coeffs, &dec,
              &sr, &aom_context_dec, cb_out.coeffs_[channel][tf_i],
              &cb_out.num_coeffs_[channel][tf_i]));
          for (uint32_t i = 0; i < num_elem; ++i) {
            EXPECT_EQ(cb.coeffs_[channel][tf_i][i],
                      cb_out.coeffs_[channel][tf_i][i])
                << "pos: " << i << "  is_all_zero:" << fill_all_zero;
          }
        }
      }
    }
  }
}

INSTANTIATE_TEST_SUITE_P(
    ResidualWriterTestAllSizes, ResidualWriterTest,
    ::testing::Combine(::testing::Values(AOMTuple(BLK_4x4, 4),
                                         AOMTuple(BLK_8x8, 3),
                                         AOMTuple(BLK_8x4, 3),
                                         AOMTuple(BLK_4x8, 3),
                                         AOMTuple(BLK_16x16, 2),
                                         AOMTuple(BLK_16x8, 3),
                                         AOMTuple(BLK_8x16, 3),
                                         AOMTuple(BLK_16x32, 2),
                                         AOMTuple(BLK_32x16, 2),
                                         AOMTuple(BLK_32x32, 2)),
                       ::testing::Values(0.f, 40.f, 75.f)));

//------------------------------------------------------------------------------

}  // namespace
}  // namespace WP2
