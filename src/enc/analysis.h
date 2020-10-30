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
//   Analysis of source frame: partitioning, segmentation, ...
//
// Author: Skal (pascal.massimino@gmail.com)

#ifndef WP2_SRC_ENC_ANALYSIS_H_
#define WP2_SRC_ENC_ANALYSIS_H_

#include "src/common/lossy/block.h"
#include "src/utils/plane.h"
#include "src/utils/vector.h"
#include "src/wp2/encode.h"

namespace WP2 {

//------------------------------------------------------------------------------

class FeatureMap {
 public:
  static constexpr uint32_t kHistogramSize =
      kMinBlockSizePix * kMinBlockSizePix;

  Vector_u8 smap_;    // score map or segment-id map
  Vector_u8 cplx_;    // map of block's complexity (0 = easy, 255 = hard)
  struct NoiseStr {
    uint16_t level[3];
    uint16_t cut[3];
  };
  VectorNoCtor<NoiseStr> noise_map_;
  uint32_t step_;    // step for smap_ and noise_map_
  uint32_t cluster_map_[kHistogramSize];   // map from score to segment id
};

// Returns whether to have lossless and/or lossy tiles depending on 'config'.
GlobalParams::Type DecideGlobalParamsType(const EncoderConfig& config);

// This function will populate *gparams and *gparams->features_, taking 'config'
// into account. It must be called prior to calling the tile encoding
// (LossyEncode(), LosslessEncode(), etc.).
WP2Status GlobalAnalysis(const ArgbBuffer& rgb, const YUVPlane& yuv,
                         const CSPTransform& transf,
                         const EncoderConfig& config,
                         GlobalParams* const gparams);

//------------------------------------------------------------------------------

class FrontMgr4x4;

// Returns a block layout. 'tile_rect' is in pixels, not padded.
// Input 'blocks' are considered forced.
WP2Status ExtractBlockPartition(const EncoderConfig& config,
                                const GlobalParams& gparams,
                                const YUVPlane& yuv, const Rectangle& tile_rect,
                                VectorNoCtor<Block>* const blocks);

// Finds a segmentation of the source in clusters with similar-features.
// gparam's 'segments_' vector is created but only their susceptibility
// characteristics (variance, risk_class, ...) are filled.
// The max number of segments is supplied by 'config' (fewer segments might
// be used if appropriate). The global feature map 'gparams->features_' is
// filled with the best segment ids and scores.
WP2Status FindSegments(const YUVPlane& yuv, const EncoderConfig& config,
                       GlobalParams* const gparams);

// Using gparams.features_, assign noise/grain level to gparams->segments_
WP2Status AssignGrainParams(const YUVPlane& yuv, const EncoderConfig& config,
                            GlobalParams* const gparams);

// Returns the segment id of the 'block' based on the 'gparams'.
// The 'segment_score' can also be retrieved if not null.
uint8_t AssignSegmentId(const EncoderConfig& config,
                        const GlobalParams& gparams,
                        const Rectangle& padded_tile, const Block& block,
                        uint32_t* const segment_score = nullptr);
// Assigns a segment to each of the 'blocks' (however blocks smaller than 8x8
// will be later assigned a predicted segment-id instead of this one).
WP2Status AssignSegmentIds(const EncoderConfig& config,
                           const GlobalParams& gparams,
                           const Rectangle& padded_tile,
                           Vector<CodedBlock>* const cblocks);

// Segmentation debug information for a given rectangle.
void SegmentationGridVDebug(const EncoderConfig& config, const Rectangle& rect,
                            uint32_t tile_pos_x, uint32_t tile_pos_y,
                            float value);
// Segmentation debug information for a given block.
void SegmentationBlockVDebug(const EncoderConfig& config, const CodedBlock& cb,
                             uint32_t tile_pos_x, uint32_t tile_pos_y,
                             float score, const Segment& segment);

// Spreads the values of the 'channel' of 'src' into at most 'histogram_size'
// buckets. Returns the number of buckets used.
uint32_t GetHistogram(const ArgbBuffer& src, uint32_t channel,
                      uint32_t max_num_buckets, uint32_t histogram[],
                      uint8_t* const min_value = nullptr,
                      uint8_t* const max_value = nullptr);

// Clusters the given 'histogram' using the k-means algorithm.
// 'clusters[]' must be 'histogram_size' long, and will contain the cluster id
// for each histogram bucket. Ids are in [0:num_clusters[, num_clusters being
// returned and is at most 'max_clusters'.
// 'centers[]' is the position of each cluster center.
uint32_t ClusterHistogram(const uint32_t histogram[], uint32_t histogram_size,
                          uint32_t max_clusters, uint32_t clusters[],
                          uint32_t centers[]);

struct SourceStats {
  typedef uint64_t score_t;
  struct BlockStat {
    VectorNoCtor<score_t> scores;
    uint8_t best_mode;
    score_t best_score;
  };
  uint32_t num_block_predictors;   // number of predictors per block
  uint32_t bw, bh;                 // dimension of image in block units
  Vector<BlockStat> stats;

  WP2Status AnalyzeSource(const YUVPlane& yuv,
                          const PredictorVector& preds);

  WP2Status ExtractToMap(Vector_u8* const map) const;
};

//------------------------------------------------------------------------------
// Prediction block size: kPredWidth * kPredHeight

// Applies forced predictors from the EncoderConfig (if any). Returns true if
// this block's predictor was forced.
bool SetForcedPredictor(const Predictors& preds, const EncoderConfig& config,
                        uint32_t tile_pos_x, uint32_t tile_pos_y,
                        Channel channel, CodedBlock* const cb);

//------------------------------------------------------------------------------
// Near-lossless

// Allocate 'out_buffer' and pre-process 'in_buffer' near-losslessly, based
// on config's quality setting.
WP2Status PreprocessNearLossless(const ArgbBuffer& in_buffer,
                                 const EncoderConfig& config, bool is_alpha,
                                 ArgbBuffer* const out_buffer);

//------------------------------------------------------------------------------

}  // namespace WP2

#endif  /* WP2_SRC_ENC_ANALYSIS_H_ */
