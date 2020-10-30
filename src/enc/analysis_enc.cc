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
//  Source analysis
//
// Author: Skal (pascal.massimino@gmail.com)
//

#include "src/dsp/dsp.h"
#include "src/enc/analysis.h"
#include "src/enc/wp2_enc_i.h"
#include "src/utils/vector.h"

// #define PRINT_STATS   // define, perform stats on predictors

#if defined(PRINT_STATS)
#include "src/utils/stats.h"
#endif

namespace WP2 {

WP2Status SourceStats::AnalyzeSource(const YUVPlane& yuv,
                                     const PredictorVector& preds) {
  WP2TransformInit();
  num_block_predictors = preds.size();
  bw = SizeBlocks(yuv.Y.w_);
  bh = SizeBlocks(yuv.Y.h_);
  WP2_CHECK_ALLOC_OK(stats.resize(bw * bh));

  const Plane16& luma = yuv.GetChannel(kYChannel);
  CodedBlock cb;
  ContextCache context_cache;
  int16_t values[kPredSize];
  int32_t res[kPredSize];
  auto* stat = &stats[0];
  for (uint32_t y = 0; y < bh; ++y) {
    for (uint32_t x = 0; x < bw; ++x, ++stat) {
      WP2_CHECK_ALLOC_OK(stat->scores.resize(num_block_predictors));
      cb.SetDimDefault({x, y, BLK_4x4});
      cb.SetSrcInput(yuv);
      cb.SetContextInput(yuv, &context_cache);  // Predict from source pixels.
      uint32_t p = 0;
      for (const auto& P : preds) {
        P->Predict(cb, kYChannel, values, kPredWidth);
        for (size_t j = 0; j < kPredHeight; ++j) {
          const int16_t* const src =
              &luma.At(x * kPredWidth, y * kPredHeight + j);
          for (size_t i = 0; i < kPredWidth; ++i) {
            const int16_t v = values[i + kPredWidth * j];
            res[i + j * kPredWidth] = src[i] - v;
          }
        }
        // TODO(skal): special scoring function
        WP2Transform2D(res, kDct, kDct, 4, 4, res);
        score_t score = 0;
        for (uint32_t i = 0; i < 16; ++i) {
          // shift down to reduce noise sensitivity
          score += (res[i] * res[i]) >> 12;
        }
        stat->scores[p] = score;
        if (p == 0 || score < stat->best_score) {
          stat->best_score = score;
          stat->best_mode = p;
        }
        ++p;
      }
    }
  }
  return WP2_STATUS_OK;
}

WP2Status SourceStats::ExtractToMap(Vector_u8* const map) const {
#if defined(PRINT_STATS)
  static Stats<uint8_t> stats0("old", "%1u", true), stats1("new", "%1u", true);
  for (uint32_t i = 0; i < bw * bh; ++i) {
    stats0.Add((*map)[i]);
    stats1.Add(stats[i].best_mode);
  }
#endif

  WP2_CHECK_ALLOC_OK(map->resize(bw * bh));
  for (uint32_t i = 0; i < bw * bh; ++i) {
    (*map)[i] = stats[i].best_mode;
  }
  return WP2_STATUS_OK;
}

//------------------------------------------------------------------------------

static void AssignGrainLevel(const EncoderConfig& config,
                             Vector<Segment>* const segments) {
  bool store_grain = config.store_grain;
  if (store_grain) {
    bool have_grain = false;
    for (const Segment& s : *segments) {
      have_grain |= s.grain_.IsUsed();
      if (have_grain) break;
    }
    if (!have_grain) store_grain = false;
  }
  if (!store_grain) {
    for (Segment& s : *segments) s.grain_.Reset();
  }
}

#if 0   // unused for now
// Check if some segments are identical (same quantization) and reduce their
// number accordingly. This usually happens at very low or very high quality.
static WP2Status SimplifySegments(Vector<Segment>* const segments) {
  uint32_t size = 1;
  for (uint32_t i = 1; i < segments->size(); ++i) {
    bool found = false;
    for (uint32_t k = 0; k < size; ++k) {
      found = (*segments)[i].IsMergeableWith((*segments)[k]);
      if (found) break;
    }
    if (!found) {
      if (i > size) {
        // WP2_CHECK_STATUS((*segments)[size].CopyFrom((*segments)[i]));
        std::swap((*segments)[size], (*segments)[i]);
      }
      ++size;
    }
  }
  WP2_CHECK_ALLOC_OK(segments->resize(size));
  return WP2_STATUS_OK;
}
#endif

//------------------------------------------------------------------------------

GlobalParams::Type DecideGlobalParamsType(const EncoderConfig& config) {
  // Neural compression is considered as GP_BOTH for now because it is lossy
  // compression outputting RGB.
  if (config.use_neural_compression) return GlobalParams::GP_BOTH;
  // TODO(skal): mixed lossy / lossless case.
  if (config.quality <= kMaxLossyQuality) return GlobalParams::GP_LOSSY;
  return GlobalParams::GP_LOSSLESS;
}

WP2Status GlobalAnalysis(const ArgbBuffer& rgb, const YUVPlane& yuv,
                         const CSPTransform& transf,
                         const EncoderConfig& config,
                         GlobalParams* const gparams) {
  WP2_CHECK_OK(gparams != nullptr, WP2_STATUS_NULL_PARAMETER);
  WP2_CHECK_OK(!rgb.IsEmpty() || !yuv.IsEmpty(), WP2_STATUS_INVALID_PARAMETER);

  gparams->type_ = DecideGlobalParamsType(config);

  if (gparams->type_ != GlobalParams::GP_LOSSY) {
    WP2_CHECK_OK(!rgb.IsEmpty(), WP2_STATUS_INVALID_PARAMETER);
    gparams->has_alpha_ = rgb.HasTransparency();
    // TODO(skal): extract global params from lossless code
  }
  if (gparams->type_ != GlobalParams::GP_LOSSLESS) {
    WP2_CHECK_OK(!yuv.IsEmpty(), WP2_STATUS_INVALID_PARAMETER);
    // If CSPTransform::Init() failed with kCustom, it felt back to another Csp.
    WP2_CHECK_OK(
        transf.GetType() == config.csp_type || config.csp_type == Csp::kCustom,
        WP2_STATUS_INVALID_PARAMETER);

    WP2MathInit();
    WP2TransformInit();
    PredictionInit();

    gparams->partition_set_ = config.partition_set;
    gparams->partition_snapping_ = config.partition_snapping;
    gparams->explicit_segment_ids_ =
        (config.segment_id_mode == WP2::EncoderConfig::SEGMENT_ID_AUTO) ?
            (config.quality > 15.f) :
        (config.segment_id_mode == WP2::EncoderConfig::SEGMENT_ID_EXPLICIT);
    gparams->use_rnd_mtx_ = config.use_random_matrix;
    gparams->transf_ = transf;
    gparams->u_quant_multiplier_ = config.u_quant_multiplier;
    gparams->v_quant_multiplier_ = config.v_quant_multiplier;

    // perform perceptual analysis to extract segments
    WP2_CHECK_STATUS(FindSegments(yuv, config, gparams));
    WP2_CHECK_STATUS(gparams->AssignQuantizations(config));
    WP2_CHECK_STATUS(gparams->AssignAlphaQuantizations(yuv, config));
    WP2_CHECK_STATUS(AssignGrainParams(yuv, config, gparams));
    AssignGrainLevel(config, &gparams->segments_);

    // TODO(skal): not working yet, we need to remap cluster_map[]
    // WP2_CHECK_STATUS(SimplifySegments(&gparams->segments_));

    if (config.info != nullptr) {
      for (const auto& forced_segment : config.info->force_segment) {
        WP2_CHECK_OK(forced_segment.segment_id < gparams->segments_.size(),
                     WP2_STATUS_INVALID_CONFIGURATION);
      }
    }

    // No decision for large Y and U/V predictors.
    WP2_CHECK_STATUS(gparams->InitFixedPredictors());
  }
  return WP2_STATUS_OK;
}

}  // namespace WP2
