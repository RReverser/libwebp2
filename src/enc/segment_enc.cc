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
//  Analysis and segmentation of frame
//
// Author: Skal (pascal.massimino@gmail.com)
//

#include <algorithm>
#include <cmath>

#include "src/common/integral.h"
#include "src/common/lossy/context.h"
#include "src/common/lossy/residuals.h"
#include "src/enc/analysis.h"
#include "src/enc/wp2_enc_i.h"

// #define LOG_NOISE   // uncomment to print noise measurements

namespace WP2 {

//------------------------------------------------------------------------------
// Segmentation

uint32_t GetHistogram(const ArgbBuffer& src, uint32_t channel,
                      uint32_t max_num_buckets, uint32_t histogram[],
                      uint8_t* const min_value, uint8_t* const max_value) {
  assert(max_num_buckets > 0);
  const uint32_t num_channels = WP2FormatBpp(src.format);
  const uint8_t* row = (const uint8_t*)src.GetRow(0) + channel;
  uint8_t min = row[0], max = row[0];
  for (uint32_t y = 0; y < src.height; ++y) {
    for (uint32_t x = 0; x < src.width; ++x) {
      const uint8_t value = row[x * num_channels];
      min = std::min(min, value);
      max = std::max(max, value);
    }
    row += src.stride;
  }
  if (min_value != nullptr) *min_value = min;
  if (max_value != nullptr) *max_value = max;

  std::fill(histogram, histogram + max_num_buckets, 0u);
  const uint32_t range = std::max(1, max - min);
  max_num_buckets = std::min(range, max_num_buckets);
  const uint32_t max_bucket = max_num_buckets - 1;

  row = (const uint8_t*)src.GetRow(0) + channel;
  for (uint32_t y = 0; y < src.height; ++y) {
    for (uint32_t x = 0; x < src.width; ++x) {
      const uint8_t value = row[x * num_channels];
      const uint32_t bucket =
          DivRound((uint32_t)(value - min) * max_bucket, range);
      ++histogram[bucket];
    }
    row += src.stride;
  }
  return max_num_buckets;
}

uint32_t ClusterHistogram(const uint32_t histogram[], uint32_t histogram_size,
                          uint32_t max_clusters, uint32_t clusters[],
                          uint32_t centers[kMaxNumSegments]) {
  assert(max_clusters <= kMaxNumSegments);
  max_clusters = std::min(max_clusters, histogram_size);
  if (max_clusters == 0) return 0;

  uint32_t total_sum = 0;
  for (uint32_t i = 0; i < histogram_size; ++i) total_sum += histogram[i];
  max_clusters = std::min(max_clusters, total_sum);
  if (max_clusters == 0) {
    // corner case
    for (uint32_t c = 0; c < histogram_size; ++c) clusters[c] = 0;
    centers[0] = histogram_size / 2;
    return 1;
  }

  // initial center positions, split by roughly equal area
  uint32_t last_sum = 0, cur_sum = 0;
  uint32_t last_wsum = 0, cur_wsum = 0;  // wsum = weighted-sum
  uint32_t num_clusters = 0;
  for (uint32_t i = 0; num_clusters < max_clusters && i < histogram_size; ++i) {
    cur_sum += histogram[i];
    cur_wsum += histogram[i] * i;
    if (cur_sum > last_sum &&
        cur_sum * max_clusters > total_sum * num_clusters) {
      const uint32_t avg = (cur_wsum - last_wsum) / (cur_sum - last_sum);
      if (num_clusters == 0 || centers[num_clusters - 1] != avg) {
        centers[num_clusters++] = avg;
      }
      last_sum = cur_sum;
      last_wsum = cur_wsum;
    }
  }

  const uint32_t kMaxIter = 3;
  for (uint32_t iter = 0; iter < kMaxIter; ++iter) {
    uint32_t wsum[kMaxNumSegments] = {0}, sum[kMaxNumSegments] = {0};
    // brute-force loop to assign centers
    for (uint32_t c = 0; c < histogram_size; ++c) {
      uint32_t best_dist = 0;
      uint32_t best_center = 0;
      for (uint32_t i = 0; i < num_clusters; ++i) {
        const uint32_t d = std::abs((int)centers[i] - (int)c);
        if (i == 0 || d < best_dist) {
          best_dist = d;
          best_center = i;
        }
      }
      clusters[c] = best_center;
      wsum[best_center] += histogram[c] * c;
      sum[best_center] += histogram[c];
    }
    if (iter + 1 == kMaxIter) break;

    // move centers closer to average
    uint32_t k = 0;
    bool changed = false;
    for (uint32_t i = 0; i < num_clusters; ++i) {
      if (sum[i] > 0) {
        const uint32_t new_center = wsum[i] / sum[i];
        if (new_center != centers[k]) {
          centers[k] = new_center;
          changed = true;
        }
        ++k;
      }
    }
    num_clusters = k;
    if (!changed) break;
  }
  return num_clusters;
}

WP2Status FindSegments(const YUVPlane& yuv, const EncoderConfig& config,
                       GlobalParams* const gparams) {
  assert(gparams != nullptr && gparams->features_ != nullptr);

  const uint32_t bw = SizeBlocks(yuv.Y.w_);
  const uint32_t bh = SizeBlocks(yuv.Y.h_);

  // Create a histogram of the variance values of all blocks.
  FeatureMap* const map = gparams->features_;
  map->step_ = bw;
  WP2_CHECK_ALLOC_OK(map->smap_.resize(bw * bh));
  WP2_CHECK_ALLOC_OK(map->noise_map_.resize(bw * bh));
  const bool store_cplx = config.tune_perceptual;
  if (store_cplx) WP2_CHECK_ALLOC_OK(map->cplx_.resize(bw * bh));
  uint32_t histogram[FeatureMap::kHistogramSize] = {0};

  Integral integral;
  WP2_CHECK_STATUS(integral.Allocate(bw, bh, kMinBlockSizePix));
  integral.AddValues(yuv);

  WP2TransformInit();
  const float factor = (kMaxLossyQuality - config.quality) / kMaxLossyQuality;
  const uint32_t cut_off[2 /* 0=Y, 1=U/V */] = {1u + (uint32_t)(factor * 100.),
                                                1u + (uint32_t)(factor * 100.)};
  const uint16_t* const zigzag = ResidualIterator::GetZigzag(TRF_4x4);
  for (uint32_t y = 0; y < bh; ++y) {
    const uint32_t y_pix = y * kMinBlockSizePix;
    const uint32_t max_j = std::min(y_pix + kMinBlockSizePix, yuv.Y.h_);

    for (uint32_t x = 0; x < bw; ++x) {
      const uint32_t idx = x + y * bw;
      float variance = integral.StdDev(x, y, x + 1, y + 1);
      if (store_cplx) {
        map->cplx_[idx] = (uint8_t)std::min(255.f, 40.f * variance);
      }
      if (variance != 0) variance = std::log(variance);

      // Convert variance to a "risk" factor. In general, lower variance means
      // higher risk, because there is less detail and if you quantize too much
      // all the detail will be gone. However, very low variance is also low
      // risk because the region is very smooth. So the function nmapping from
      // variance to risk is shaped like this: /\.

      // Variance that has the most risk. The risk decreases the farther away
      // the variance is from this value.
      const float kHighestRiskVariance = 2.f;

      const float kMax = 7.f;
      const float kCutoff = kMax - kHighestRiskVariance;
      float risk = kMax - Clamp(variance, 0.f, kMax);
      if (risk > kCutoff) {
        risk = kCutoff * (1 - (risk - kCutoff) / (kMax - kCutoff));
      }
      const uint8_t score =
          std::round(risk / kCutoff * (FeatureMap::kHistogramSize - 1));
      assert(score < FeatureMap::kHistogramSize);
      map->smap_[idx] = (uint8_t)score;
      ++histogram[score];

      // Noise estimation.
      const uint32_t x_pix = x * kMinBlockSizePix;
      const uint32_t max_i = std::min(x_pix + kMinBlockSizePix, yuv.Y.w_);
      int32_t res[kMinBlockSizePix * kMinBlockSizePix] = {0};
      for (const auto channel : {kYChannel, kUChannel, kVChannel}) {
        const Plane16& p = yuv.GetChannel(channel);
        for (uint32_t k = 0, j = y_pix; j < max_j; ++j) {
          for (uint32_t i = x_pix; i < max_i; ++i) res[k++] = p.At(i, j);
        }
        WP2Transform2D(res, WP2TransformType::kDct, WP2TransformType::kDct,
                       kMinBlockSizePix, kMinBlockSizePix, res,
                       /*reduced=*/false);
        uint32_t cnt = 0;
        int16_t last = 0;
        for (uint32_t k = 1; k < kMinBlockSizePix * kMinBlockSizePix; ++k) {
          const uint32_t l = zigzag[k];
          const uint32_t v = std::abs(res[l]);
          if (v > cut_off[channel != kYChannel]) {
            last = l;
            cnt = 0;  // restart counter
          } else {
            ++cnt;
          }
        }
        map->noise_map_[idx].level[channel] = cnt;
        const uint32_t X = last % kMinBlockSizePix;
        const uint32_t Y = last / kMinBlockSizePix;
        map->noise_map_[idx].cut[channel] = X + Y;
      }

      if (VDMatch(config, "encoder/segmentation")) {
        const Rectangle rect = {x * kMinBlockSizePix, y * kMinBlockSizePix,
                                kMinBlockSizePix, kMinBlockSizePix};
        const uint32_t tile_x = 0, tile_y = 0;
        if (VDMatch(config, "variance")) {
          SegmentationGridVDebug(config, rect, tile_x, tile_y, variance);
        } else if (VDMatch(config, "score")) {
          SegmentationGridVDebug(config, rect, tile_x, tile_y, score);
        }
      }
    }
  }

  uint32_t num_non_empty_buckets = 0;
  for (uint32_t i = 0; i < FeatureMap::kHistogramSize; ++i) {
    if (histogram[i] != 0) ++num_non_empty_buckets;
  }

  // Cluster the histogram.
  const uint32_t max_num_segments =
      std::max(1u, std::min({(uint32_t)config.segments, num_non_empty_buckets,
                             GetMaxNumSegments(gparams->explicit_segment_ids_,
                                               GetQualityHint(config.quality),
                                               config.partition_set)}));

  // Cluster id for each variance value.
  uint32_t centers[kMaxNumSegments];
  const uint32_t num_segments =
      ClusterHistogram(histogram, FeatureMap::kHistogramSize, max_num_segments,
                       map->cluster_map_, centers);
  assert(num_segments > 0 && num_segments <= max_num_segments);
  Vector<Segment>& segments = gparams->segments_;  // shortcut
  WP2_CHECK_ALLOC_OK(segments.resize(num_segments));
  for (uint32_t i = 0; i < num_segments; ++i) {
    segments[i].risk_ = centers[i] / (float)FeatureMap::kHistogramSize;
    segments[i].risk_class_ = i;
    segments[i].avg_lambda_ = Segment::kDefaultAvgLambda;
  }
  if (config.tune_perceptual) {
    uint32_t counts[kMaxNumSegments] = { 0 };
    uint32_t avg_cplx[kMaxNumSegments] = { 0 };
    for (uint32_t idx = 0; idx < map->smap_.size(); ++idx) {
      const uint32_t id = map->cluster_map_[map->smap_[idx]];
      counts[id] += 1;
      avg_cplx[id] += map->cplx_[idx];
    }
    for (uint32_t i = 0; i < num_segments; ++i) {
      if (counts[i] > 0) {
        segments[i].avg_lambda_ = 1.f * avg_cplx[i] / counts[i];
      }
    }
  }
  return WP2_STATUS_OK;
}

WP2Status AssignGrainParams(const YUVPlane& yuv, const EncoderConfig& config,
                            GlobalParams* const gparams) {
  const FeatureMap* const map = gparams->features_;
  Vector<Segment>* const segments = &gparams->segments_;   // shortcut
  if (!config.store_grain) {
    for (auto& s : *segments) s.use_grain_ = false;
    return WP2_STATUS_OK;
  }

  // Noise/Grain level records
  uint64_t noise[3][kMaxNumSegments] = {{0}, {0}};
  uint64_t cut[3][kMaxNumSegments] = {{0}, {0}};
  uint64_t noise_cnt[kMaxNumSegments] = {0};
  for (uint32_t idx = 0; idx < map->smap_.size(); ++idx) {
    const uint32_t id = map->cluster_map_[map->smap_[idx]];
    for (const auto channel : {kYChannel, kUChannel, kVChannel}) {
      noise[channel][id] += map->noise_map_[idx].level[channel];
      cut[channel][id] += map->noise_map_[idx].cut[channel];
    }
    noise_cnt[id] += 1;
  }

  const uint32_t num_segments = segments->size();
#if defined(LOG_NOISE)
  for (auto c : {kYChannel, kUChannel, kVChannel}) {
    for (uint32_t i = 0; i < num_segments; ++i) {
      printf("Noise #%d: %.2f %%\n",
             i, 100.f * noise[c][i] / noise_cnt[i] / 16.);
    }
    printf("--- [channel:%d] --- \n", c);
  }
  printf("=========================\n");
#endif
  for (uint32_t s = 0; s < num_segments; ++s) {
    Segment* const seg = &(*segments)[s];
    uint32_t level = 0;
    if (noise_cnt[s] > 0) {
      for (auto c : {kYChannel, kUChannel, kVChannel}) {
        level += noise[c][s];
      }
      level /= noise_cnt[s];
      seg->grain_.y_ = 15 * noise[0][s] / noise_cnt[s] / 16;
      seg->grain_.uv_ = 15 * (noise[1][s] + noise[2][s]) / noise_cnt[s] / 32;
      const uint32_t cut_y = cut[0][s] / noise_cnt[s];
      const uint32_t cut_uv = (cut[1][s] + cut[2][s]) / noise_cnt[s] / 2;
      seg->grain_.cut_y_ = std::min(cut_y, 7u);
      seg->grain_.cut_uv_ = std::min(cut_uv, 7u);
    }
    seg->use_grain_ = (level > 0);
  }
  return WP2_STATUS_OK;
}

uint8_t AssignSegmentId(const EncoderConfig& config,
                        const GlobalParams& gparams,
                        const Rectangle& padded_tile, const Block& block,
                        uint32_t* const segment_score) {
  if (config.info != nullptr && gparams.explicit_segment_ids_) {
    for (const auto& forced_segment : config.info->force_segment) {
      if (block.x_pix() + padded_tile.x == forced_segment.x &&
          block.y_pix() + padded_tile.y == forced_segment.y) {
        if (segment_score != nullptr) {
          *segment_score = FeatureMap::kHistogramSize;
        }
        return forced_segment.segment_id;
      }
    }
  }

  if (block.IsSmall() || !gparams.explicit_segment_ids_) {
    return SegmentIdPredictor::GetIdFromSize(
        (uint32_t)gparams.segments_.size(), block.dim());
  }

  const FeatureMap& map = *gparams.features_;
  const uint32_t pos_x = padded_tile.x / kMinBlockSizePix;
  const uint32_t pos_y = padded_tile.y / kMinBlockSizePix;
  const uint8_t* const smap = &map.smap_[pos_x + map.step_ * pos_y];

  const uint32_t x0 = block.x(), x1 = block.x() + block.w();
  const uint32_t y0 = block.y(), y1 = block.y() + block.h();
  assert(x1 <= padded_tile.width && y1 <= padded_tile.height);

  // Strategy for choosing the score for a whole block.
  enum ScoringStrategy {
    kTopLeft,
    kMostFrequent,
    kAverage,
  } method = kAverage;

  uint32_t score = 0;
  switch (method) {
    case kTopLeft: {
      // Simple (bad) solution: just choose the top-left one.
      score = smap[x0 + y0 * map.step_];
      break;
    }
    case kMostFrequent: {
      // we use the most frequent cluster-id ("distribution mode")
      uint32_t counts[FeatureMap::kHistogramSize] = {0};
      uint32_t max = 0;
      uint32_t max_id = 0;
      for (uint32_t y = y0; y < y1; ++y) {
        for (uint32_t x = x0; x < x1; ++x) {
          const uint32_t id = smap[x + y * map.step_];
          assert(id < FeatureMap::kHistogramSize);
          ++counts[id];
          if (counts[id] > max) {
            max = counts[id];
            max_id = id;
          }
        }
      }
      score = max_id;
      break;
    }
    case kAverage: {
      // Even simpler method: take the average (ids are mapped proportional to
      // complexity, so it makes sense to take the average).
      score = 0;
      for (uint32_t y = y0; y < y1; ++y) {
        for (uint32_t x = x0; x < x1; ++x) {
          const uint32_t id = smap[x + y * map.step_];
          score += id;
          assert(id < FeatureMap::kHistogramSize);
        }
      }
      const uint32_t area = (x1 - x0) * (y1 - y0);
      score /= area;
      break;
    }
  }
  // Earlier, we stored the score in 'map'. Now we map it to a cluster id.
  assert(score < FeatureMap::kHistogramSize);
  if (segment_score != nullptr) *segment_score = score;
  return map.cluster_map_[score];
}

// block.id_ must be assigned prior to calling this function
static float AssignLambda(const EncoderConfig& config,
                          const GlobalParams& gparams,
                          const Rectangle& padded_tile,
                          const CodedBlock& cb) {
  if (!config.tune_perceptual) return Segment::kDefaultAvgLambda;

  const FeatureMap& map = *gparams.features_;
  const uint32_t pos_x = padded_tile.x / kMinBlockSizePix;
  const uint32_t pos_y = padded_tile.y / kMinBlockSizePix;

  const uint32_t x0 = cb.blk().x(), w = cb.blk().w();
  const uint32_t y0 = cb.blk().y(), h = cb.blk().h();
  assert(x0 + w <= padded_tile.width && y0 + h <= padded_tile.height);

  uint32_t avg_cplx = 0;
  const uint8_t* cplx_map = &map.cplx_[pos_x + map.step_ * (pos_y + y0)];
  for (uint32_t y = y0; y < y0 + h; ++y) {
    for (uint32_t x = x0; x < x0 + w; ++x) avg_cplx += cplx_map[x];
    cplx_map += map.step_;
  }
  const float avg_lambda = 1.f * avg_cplx / (w * h);
  const float seg_lambda = gparams.segments_[cb.id_].avg_lambda_;
  if (seg_lambda < 0.00001f) return Segment::kDefaultAvgLambda;

  // The complexity multiplier is estimated as:
  //   1 + multiplier * (local_cplx - average_cplx) / average_cplx
  // and we restrict to a reasonable [0.8, 1.4] range to avoid crazy cases.
  // The multiplier linearly depends on the other perceptual parameter 'sns'.
  const float multiplier = 10.f * config.sns / 100.f;
  return Clamp(1.0f + multiplier * (avg_lambda - seg_lambda) / seg_lambda,
               0.8f, 1.4f);
}

WP2Status AssignSegmentIds(const EncoderConfig& config,
                           const GlobalParams& gparams,
                           const Rectangle& padded_tile,
                           Vector<CodedBlock>* const cblocks) {
  if (VDMatch(config, "encoder/segmentation/id")) {
    const FeatureMap& map = *gparams.features_;
    const uint32_t tile_x = padded_tile.x / kMinBlockSizePix;
    const uint32_t tile_y = padded_tile.y / kMinBlockSizePix;
    for (uint32_t y = 0; y < SizeBlocks(padded_tile.height); ++y) {
      for (uint32_t x = 0; x < SizeBlocks(padded_tile.width); ++x) {
        const uint32_t score =
            gparams.features_->smap_[(tile_y + y) * map.step_ + (tile_x + x)];
        CodedBlock cb;
        cb.SetDimDefault(Block(x, y, BLK_4x4));
        cb.id_ = map.cluster_map_[score];
        SegmentationBlockVDebug(config, cb, padded_tile.x, padded_tile.y, score,
                                gparams.segments_[cb.id_]);
      }
    }
  }

  for (CodedBlock& cb : *cblocks) {
    uint32_t score;
    cb.id_ = AssignSegmentId(config, gparams, padded_tile, cb.blk(), &score);
    cb.lambda_mult_ = AssignLambda(config, gparams, padded_tile, cb);
    if (VDMatch(config, "encoder/segmentation/blocks")) {
      SegmentationBlockVDebug(config, cb, padded_tile.x, padded_tile.y, score,
                              gparams.segments_[cb.id_]);
    }
    assert(cb.id_ < gparams.segments_.size());
  }
  return WP2_STATUS_OK;
}

//------------------------------------------------------------------------------

}  // namespace WP2
