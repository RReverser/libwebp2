// Copyright 2019 Google LLC
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//     https://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.
// -----------------------------------------------------------------------------
//
// WP2 lossy decoding.
//
// Author: Skal (pascal.massimino@gmail.com)

#include <algorithm>
#if !defined(WP2_BITTRACE)
#include <cstdio>
#endif
#include <utility>

#include "src/common/constants.h"
#include "src/common/filters/rstr_flt_params.h"
#include "src/common/global_params.h"
#include "src/common/lossy/residuals_aom.h"
#include "src/dec/filters/block_map_filter.h"
#include "src/dec/filters/intertile_filter.h"
#include "src/dec/filters/intratile_filter.h"
#include "src/dec/wp2_dec_i.h"
#include "src/dsp/math.h"
#include "src/utils/ans_utils.h"
#include "src/utils/front_mgr.h"
#include "src/utils/utils.h"
#include "src/wp2/decode.h"

namespace WP2 {

//------------------------------------------------------------------------------

#if defined(WP2_BITTRACE)

// Threshold above which ANSDec::GetBitCount() not matching the sum of bit
// traces is an issue not linked to floating point precision.
static constexpr double kDoubleSumTolerance = 0.0001;

// Returns the absolute normalized difference of 'a' and 'b' in [0:2].
static double NormDiff(double a, double b, double reduce_diff_by = 0) {
  const double diff = std::abs(a - b);
  if (diff <= reduce_diff_by) return 0.;
  return (diff - reduce_diff_by) / std::max(std::abs(a), std::abs(b));
}

#endif  // defined(WP2_BITTRACE)

//------------------------------------------------------------------------------

WP2Status LossyDecode(const BitstreamFeatures& features,
                      const DecoderConfig& config,
                      TilesLayout* const tiles_layout,
                      ProgressWatcher* const progress, ANSDec* const dec,
                      Tile* const tile) {
  WP2TransformInit();
  WP2QuantizeInit();
  PredictionInit();

  assert(tile != nullptr);
  tile->num_decoded_rows = 0;

  dec->AddDebugPrefix("GlobalHeader");
  const bool is_neural =
      (tiles_layout->gparams->type_ == GlobalParams::GP_BOTH &&
       dec->ReadBool("is_neural"));
  dec->PopDebugPrefix();
  if (is_neural) {
    assert(!tile->rgb_output.IsEmpty());
    tile->output_is_yuv = false;
    return NeuralDecode(&tile->rgb_output, dec, features);
  }
  assert(!tile->yuv_output.IsEmpty());
  tile->output_is_yuv = true;

  const GlobalParams& gparams = *tiles_layout->gparams;
  SyntaxReader reader(dec, tile->rect);
  WP2_CHECK_STATUS(reader.ReadHeader(features, gparams, *tile));

  tile->block_map.Init(tile->rect, gparams.transf_.GetYUVPrecisionBits() + 1u,
                       gparams.transf_.GetYUVMin(), gparams.transf_.GetYUVMax(),
                       &tile->yuv_output);
  if (IsFilterBlockMapNeeded(config, features, gparams, *tiles_layout)) {
    WP2_CHECK_STATUS(tile->block_map.Allocate());
  }

  RstrFltParams rstr_params(tile->rect.width, tile->rect.height);
  WP2_CHECK_STATUS(ReadRestorationFilterParams(dec, &rstr_params));
  IntratileFilter filter(config, features, gparams, tile->block_map,
                         rstr_params);
  WP2_CHECK_STATUS(filter.Allocate());

  // VisualDebug
  Vector<CodedBlock> all_blocks;
  ArgbBuffer debug_output;  // Tile view of 'config.info->debug_output'
  Plane16 vd_plane;         // Visual debug plane.
  Plane16 vd_block_view;    // View of vd_plane for the current block.
  if (VDMatch(config, "")) {
    WP2_CHECK_STATUS(
        debug_output.SetView(config.info->debug_output, tile->rect));
    if (VDMatch(config, "prediction/raw") || VDMatch(config, "residuals") ||
        VDMatch(config, "coeff-method") || VDMatch(config, "bits-per-pixel") ||
        VDMatch(config, "chroma-from-luma") || VDMatch(config, "a/lossless") ||
        VDMatch(config, "a/lossy") || VDMatch(config, "a/is_lossy")) {
      WP2_CHECK_STATUS(
          vd_plane.Resize(tile->yuv_output.Y.w_, tile->yuv_output.Y.h_));
      // Areas that are not overwritten will show as black in the debug UI.
      // This is only relevant for the alpha prediction/residual/modes where
      // we only show lossy blocks.
      vd_plane.Fill(-512);
    }
  }
  if (VDMatch(config, "bits-per-pixel") ||
      (config.info != nullptr && config.info->store_blocks)) {
#if !defined(WP2_BITTRACE)
    static bool warning_printed = false;
    if (!warning_printed) {
      printf("Warning! VDEBUG_BITS_PER_PIXEL or store_blocks won't work ");
      printf("without WP2_BITTRACE compile flag!\n");
      warning_printed = true;
    }
    if (VDMatch(config, "bits-per-pixel")) vd_plane.Clear();
#endif
  }

#if defined(WP2_BITTRACE)
  double blocks_bit_traces_sum = dec->GetBitCount();  // Start with header.
  // Mapping from block index to block bit count.
  VectorNoCtor<double> size_bits;
#endif

  FrontMgrDefault mgr;
  WP2_CHECK_STATUS(mgr.Init(reader.GetPartitionSet(), reader.GetPartitionSnap(),
                            tile->rect.width, tile->rect.height));

  while ((dec->GetStatus() == WP2_STATUS_OK) && !mgr.Done()) {
    // Read dimensions until we get a full context.
    FrontMgrNxNBase::Info info;
    while (true) {
      if (mgr.UseFinal(&info)) break;
      // Read block sizes until a final one can be used.
#if defined(WP2_BITTRACE)
      const double prev_bit_count = dec->GetBitCount();
#endif
      const BlockSize dim = reader.GetBlockSize(mgr);
      Block block;
#if defined(WP2_BITTRACE)
      const uint32_t block_ind = size_bits.size();
      WP2_CHECK_ALLOC_OK(mgr.UseSize(dim, block_ind, &block));
      const double num_bits = dec->GetBitCount() - prev_bit_count;
      WP2_CHECK_ALLOC_OK(size_bits.push_back(num_bits));
#else
      WP2_CHECK_ALLOC_OK(mgr.UseSize(dim, /*ind=*/0u, &block));
#endif
    }

    // When reaching the end of a line, it is certain that all rows above
    // current block are decoded.
    if (info.block.x_pix() + info.block.w_pix() == tile->yuv_output.Y.w_) {
      const uint32_t num_decoded_yuv_rows = info.block.y_pix();
      tile->num_decoded_rows = filter.Apply(num_decoded_yuv_rows);
    }

    CodedBlock cb;
    cb.SetRange(gparams.transf_.GetYUVMin(), gparams.transf_.GetYUVMax());
    cb.SetDim(info.block, mgr);
    ContextCache context_cache;
    cb.SetContextCache(&context_cache);
    assert(cb.y_pix() + cb.h_pix() <= tile->yuv_output.Y.h_);
    cb.mtx_set_ = &gparams.mtx_set_;

    if (!vd_plane.IsEmpty()) {
      WP2_CHECK_STATUS(vd_block_view.SetView(vd_plane, cb.blk().rect_pix()));
    }

    // Estimation of the bits-per-pixel for the FilterBlockMap.
    const size_t num_used_bytes_before = dec->GetNumUsedBytes();

    // Read the block's syntax and set the output view.
#if defined(WP2_BITTRACE)
    dec->ClearBitTracesCustom();
    // Register the "block_size" bit trace during the decoding of the associated
    // block.
    dec->GetBitTracesCustom()["block_size"].bits = size_bits[info.ind];
    dec->GetBitTracesCustom()["block_size"].num_occurrences = 1;

    BlockInfo block_info;
    const double last_bit_pos = dec->GetBitCount();
    WP2_CHECK_STATUS(reader.GetBlock(&cb, &tile->yuv_output, &block_info));
#else
    WP2_CHECK_STATUS(reader.GetBlock(&cb, &tile->yuv_output, /*info=*/nullptr));
#endif

    // reconstruct
    const Segment& segment = gparams.segments_[cb.id_];
    for (Channel c : {kYChannel, kUChannel, kVChannel}) {
      const QuantMtx& quant = segment.GetQuant(c);
      BlockCoeffs32 res;
      cb.Dequantize(quant, c, &res);
      cb.Reconstruct(c, (c == kYChannel) ? false : cb.is420_, &res);
    }
    if (gparams.has_alpha_) {
      reader.alpha_reader().Reconstruct(&cb);

      if (VDMatch(config, "a/lossless")) {
        reader.alpha_reader().ReconstructLossless(&cb, &vd_block_view);
      } else if (VDMatch(config, "a/lossy") && cb.HasLossyAlpha()) {
        WP2_CHECK_STATUS(
            vd_block_view.Copy(cb.out_.A, /*resize_if_needed=*/false));
      } else if (VDMatch(config, "a/is_lossy")) {
        vd_block_view.Fill(cb.HasLossyAlpha() ? 255 : 0);
      }
    }

#if defined(WP2_ENC_DEC_DEEP_MATCH)
    reader.ReadAndCompareRawPixels(cb.blk(), cb.out_, dec);  // Debug
#endif

    // VisualDebug
    if (VDMatch(config, "prediction/raw")) {
      const Channel c = VDChannel(config);
      if (c != kAChannel || cb.HasLossyAlpha()) {
        cb.StorePredictionModes(config, tile->rect, c, gparams.predictors(c),
                                &vd_plane, nullptr);
      }
    } else if (VDMatch(config, "chroma-from-luma/prediction")) {
      const Channel channel = VDChannel(config);
      const CodedBlock::CodingParams saved = *cb.GetCodingParams(channel);
      assert(gparams.uv_preds_[0]->DependsOnLuma());
      cb.SetUVPredictor(gparams.uv_preds_[0]);
      cb.StorePredictionModes(config, tile->rect, channel,
                              gparams.predictors(channel), &vd_plane, nullptr);
      *cb.GetCodingParams(channel) = saved;
    } else if (VDMatch(config, "chroma-from-luma/slope") ||
               VDMatch(config, "chroma-from-luma/intercept")) {
      const bool selected =
          VDSelected(tile->rect.x, tile->rect.y, cb.blk().rect_pix(), config);
      std::string* const debug_str =
          selected ? &config.info->selection_info : nullptr;

      if (VDMatch(config, "slope")) {
        cb.StoreCflSlope(VDChannel(config), gparams.transf_.GetYUVMin(),
                         gparams.transf_.GetYUVMax(), &vd_plane, debug_str);
      } else {
        cb.StoreCflIntercept(VDChannel(config), gparams.transf_.GetYUVMin(),
                             gparams.transf_.GetYUVMax(), &vd_plane, debug_str);
      }
    } else if (VDMatch(config, "transform")) {
      cb.StoreTransform(config, tile->rect.x, tile->rect.y, &debug_output);
    } else if (VDMatch(config, "residuals")) {
      const Channel c = VDChannel(config);
      if (c != kAChannel) {
        cb.StoreResiduals(config, tile->rect.x, tile->rect.y,
                          segment.GetQuant(c), c, &vd_plane);
      } else if (cb.HasLossyAlpha()) {
        cb.StoreResiduals(config, tile->rect.x, tile->rect.y,
                          gparams.a_segments_[0].quant_a_, kAChannel,
                          &vd_plane);
      }
    } else if (VDMatch(config, "prediction/modes")) {
      const Channel c = VDChannel(config);
      if (c != kAChannel || cb.HasLossyAlpha()) {
        cb.StorePredictionModes(config, tile->rect, c, gparams.predictors(c),
                                nullptr, &debug_output);
      }
    } else if (VDMatch(config, "original") && !VDMatch(config, "histogram")) {
      cb.AppendOriginalPixels(config, tile->rect.x, tile->rect.y,
                              gparams.transf_, &debug_output);
    } else if (VDMatch(config, "compressed")) {
      cb.AppendCompressedPixels(config, tile->rect.x, tile->rect.y,
                                &debug_output);
    } else if (VDMatch(config, "coeff-method")) {
      cb.StoreCoeffMethod(config, &vd_plane);
    } else if (VDMatch(config, "bits-per-pixel")) {
      reader.StoreBitCost(config, tile->rect.x, tile->rect.y, cb.blk(),
                          &vd_plane);
    }

#if defined(WP2_BITTRACE)
    const double num_bits =
        dec->GetBitCount() - last_bit_pos + size_bits[info.ind];
    double bt_sum = 0;
    for (const auto& bt : dec->GetBitTracesCustom()) bt_sum += bt.second.bits;
    // Check that all reads were registered as custom bit traces.
    if (NormDiff(bt_sum, num_bits) > kDoubleSumTolerance) assert(false);
    blocks_bit_traces_sum += bt_sum;
    if (config.info != nullptr && config.info->store_blocks) {
      cb.ToBlockInfo(&block_info);
      block_info.rect.x += tile->rect.x;
      block_info.rect.y += tile->rect.y;
      block_info.bits = num_bits;
      block_info.bit_traces.clear();
      block_info.bit_traces.insert(dec->GetBitTracesCustom().begin(),
                                   dec->GetBitTracesCustom().end());
      WP2_CHECK_STATUS(tiles_layout->assignment_lock.Acquire());
      config.info->blocks.push_back(block_info);
      tiles_layout->assignment_lock.Release();
    }
    if (config.info != nullptr && config.info->bits_per_pixel != nullptr) {
      config.info->bits_per_pixel->Fill(
          {tile->rect.x + cb.x_pix(), tile->rect.y + cb.y_pix(), cb.w_pix(),
           cb.h_pix()},
          num_bits / (cb.w_pix() * cb.h_pix()));
    }
#endif

    const uint32_t bw = cb.w();
    const uint32_t bh = cb.h();
    const uint32_t actual_block_width_px =
        std::min(bw * kMinBlockSizePix, tile->rect.width - cb.x_pix());
    const uint32_t actual_block_height_px =
        std::min(bh * kMinBlockSizePix, tile->rect.height - cb.y_pix());

    size_t min_num_used_bytes = ANSDec::GetMinNumUsedBytesDiff(
        num_used_bytes_before, dec->GetNumUsedBytes());
    if (gparams.has_alpha_) {
      // Don't include lossless alpha. TODO(maryla): remove lossy alpha as well?
      min_num_used_bytes =
          SafeSub(min_num_used_bytes, cb.alpha_lossless_bytes_);
    }
    tile->block_map.RegisterBlock(cb, (uint32_t)min_num_used_bytes);

    if (progress != nullptr) {
      WP2_CHECK_STATUS(
          progress->AdvanceBy(actual_block_width_px * actual_block_height_px));
    }

    if (VDMatch(config, "blocks")) {
      WP2_CHECK_ALLOC_OK(all_blocks.push_back(std::move(cb)));
    }
  }

  WP2_CHECK_STATUS(dec->GetStatus());
  // Not using the exact whole chunk is an issue.
  WP2_CHECK_OK(dec->GetNumUsedBytes() == tile->chunk_size,
               WP2_STATUS_BITSTREAM_ERROR);

#if defined(WP2_BITTRACE)
  // Check that custom bit traces sum up to the expected bit count.
  if (NormDiff(blocks_bit_traces_sum, dec->GetBitCount()) >
      kDoubleSumTolerance) {
    assert(false);
  }
  // Verify that the estimated bit count match the tile size.
  if (NormDiff(dec->GetBitCount(), tile->chunk_size * 8.,
               /*reduce_diff_by=*/kANSPaddingCost) >
      ANSDec::kBitCountAccuracy) {
    assert(false);
  }
#endif

  tile->num_decoded_rows = filter.Apply(tile->rect.height);
  assert(tile->num_decoded_rows == tile->rect.height);

  // VisualDebug
  if (VDMatch(config, "original") || VDMatch(config, "compressed")) {
    ApplyVDebugBeforeAfter(config, gparams.transf_, *tile, &debug_output);
  } else if (VDMatch(config, "blocks")) {
    YUVPlane non_padded;
    WP2_CHECK_STATUS(non_padded.SetView(
        tile->yuv_output, {0, 0, tile->rect.width, tile->rect.height}));
    WP2_CHECK_STATUS(non_padded.Export(
        gparams.transf_, /*resize_if_needed=*/false, &debug_output));
    for (const CodedBlock& cb : all_blocks) {
      WP2_CHECK_STATUS(
          cb.Draw(config, tile->rect.x, tile->rect.y, gparams, &debug_output));
    }
  } else if (VDMatch(config, "filter-block-map")) {
    tile->block_map.ApplyVDebug(config, &debug_output);
  } else if (!vd_plane.IsEmpty()) {
    if (VDMatch(config, "a")) {
      WP2_CHECK_STATUS(vd_plane.ToGray(&debug_output, /*num_bits=*/8,
                                       /*shift=*/VDMatch(config, "residuals")));
    } else {
      WP2_CHECK_STATUS(vd_plane.ToGray(&debug_output, /*num_bits=*/10));
    }
  }

  // Only clear data that is no longer valid. The rest might be used by the
  // IntertileFilter.
  tile->block_map.pixels_ = nullptr;
  return WP2_STATUS_OK;
}

//------------------------------------------------------------------------------
// Empty neural decoding

#if !defined(WP2_EXPERIMENTAL)
WP2Status NeuralDecode(ArgbBuffer* const picture, ANSDec* const dec,
                       const BitstreamFeatures& features) {
  (void)picture;
  (void)dec;
  (void)features;
  return WP2_STATUS_UNSUPPORTED_FEATURE;
}
#endif  // !WP2_EXPERIMENTAL

}  // namespace WP2
