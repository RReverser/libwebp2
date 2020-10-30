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
//  Predictor clustering
//
// Author: Skal (pascal.massimino@gmail.com)
//

#include <cstdio>
#include <limits>

#include "src/common/lossy/block.h"
#include "src/enc/analysis.h"
#include "src/utils/front_mgr.h"
#include "src/utils/plane.h"

namespace WP2 {

//------------------------------------------------------------------------------
// Linear predictors

static constexpr uint8_t kUnassigned = 0xff;

// 4x4 block element for stats
class ClusterElement {
 public:
  void Init(const YUVPlane* const p, Channel channel, uint32_t x, uint32_t y) {
    p_ = p;
    channel_ = channel;
    x_ = x;
    y_ = y;
    id_ = kUnassigned;
    dist_ = 0;
  }

  void Predict(const Predictor& p, const FrontMgr4x4& mgr, uint32_t step,
               int16_t prediction[]) const;
  uint32_t ComputeDistance(const int16_t prediction[]) const;
  void GetContext(int16_t context[]) const;
  void GetValues(int16_t values[]) const;

  template <typename TPredVec>  // Can be a PredictorVector or equivalent
  void AssignBestPredictor(const FrontMgr4x4& mgr, const TPredVec& preds) {
    uint8_t best_id = 0;
    uint32_t best_score = ~0u;
    for (size_t id = 0; id < preds.size(); ++id) {
      int16_t values[kPredSize];
      Predict(*preds[id], mgr, kPredWidth, values);
      const uint32_t score = ComputeDistance(values);
      if (score < best_score) {
        best_score = score;
        best_id = id;
        if (best_score == 0) {
          break;
        }
      }
    }
    dist_ = best_score;
    id_ = best_id;
  }

  mutable ContextCache context_cache_;
  const YUVPlane* p_;
  Channel channel_;
  uint32_t x_, y_;
  uint8_t id_;
  uint32_t dist_;

 private:
  const Plane16& GetPlane() const { return p_->GetChannel(channel_); }
};

void ClusterElement::Predict(const Predictor& p, const FrontMgr4x4& mgr,
                             uint32_t step, int16_t prediction[]) const {
  // Build a fake CodedBlock to represent this ClusterElement.
  CodedBlock cb;
  // yuv_min and yuv_max are used for reconstruction, which we don't care about
  // here.
  cb.SetRange(std::numeric_limits<int16_t>::min(),
              std::numeric_limits<int16_t>::max());
  cb.SetDim(Block(x_ / kMinBlockSizePix, y_ / kMinBlockSizePix, BLK_4x4), mgr);
  cb.SetContextInput(*p_, &context_cache_);
  p.Predict(cb, channel_, prediction, step);
}

// Retrieve the context (clipped) surrounding the block at x,y of size bw x bh
// (in pixel units). 'p' is the full plane (and not the sub-view), with extra
// hidden pixel boundary.
void ClusterElement::GetContext(int16_t context[/*bw + 2*(bh+1)*/]) const {
  size_t m = 0;
  // left + top_left
  const Plane16& p = GetPlane();
  for (int j = kPredHeight - 1; j >= -1; --j) {
    context[m++] = p.AtClamped((int)x_ - 1, (int)y_ + j);
  }
  // top + top-right
  for (uint32_t i = 0; i <= kPredWidth; ++i) {
    context[m++] = p.AtClamped((int)(x_ + i), (int)y_ - 1);
  }
  // right
  for (uint32_t j = 0; j < kPredHeight; ++j) {
    context[m++] = p.AtClamped(x_ + kPredWidth, y_ + j);
  }
}

void ClusterElement::GetValues(int16_t values[]) const {
  for (uint32_t k = 0, j = 0; j < kPredHeight; ++j) {
    for (uint32_t i = 0; i < kPredWidth; ++i) {
      values[k++] = GetPlane().At(x_ + i, y_ + j);
    }
  }
}

uint32_t ClusterElement::ComputeDistance(const int16_t prediction[]) const {
  const Plane16& p = GetPlane();
  // TODO(skal): better scoring
  uint32_t dist = 0;
  for (uint32_t k = 0, y = 0; y < kPredHeight; ++y) {
    for (uint32_t x = 0; x < kPredWidth; ++x, ++k) {
      const int16_t v = p.At(x_ + x, y_ + y);
      dist += (prediction[k] - v) * (prediction[k] - v);
      // dist += -log2(1. - fabs(prediction[k] - v) / 512.);
    }
  }
  return dist;
}

//------------------------------------------------------------------------------

bool SetForcedPredictor(const Predictors& preds, const EncoderConfig& config,
                        uint32_t tile_pos_x, uint32_t tile_pos_y,
                        Channel channel, CodedBlock* const cb) {
  EncoderInfo* const info = config.info;
  if (info == nullptr || info->force_predictor.empty()) return false;

  for (const EncoderInfo::ForcedPredictor& forced : info->force_predictor) {
    if (forced.channel != channel || forced.x != (cb->x_pix() + tile_pos_x) ||
        forced.y != (cb->y_pix() + tile_pos_y)) {
      continue;
    }

    assert(!preds.empty());
    const Predictor* const pred =
        preds[std::min(forced.predictor_id, (uint8_t)(preds.size() - 1))];
    cb->SetPredictor(channel, pred);
    cb->pred_scores_[(channel == kVChannel) ? kUChannel : channel] = 1.f;
    return true;
  }

  return false;
}

//------------------------------------------------------------------------------

}  // namespace WP2
