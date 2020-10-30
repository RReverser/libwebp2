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
//   Contexted predictions and syntax I/O
//
// Author: Skal (pascal.massimino@gmail.com)

#include "src/common/lossy/context.h"

#include <algorithm>
#include <functional>   // for std::greater<>

#include "src/common/lossy/block.h"
#include "src/common/lossy/predictor.h"
#include "src/dec/symbols_dec.h"
#include "src/enc/symbols_enc.h"

namespace WP2 {

constexpr uint32_t kMapHeight = kMaxBlockSize + 1;

//------------------------------------------------------------------------------
// BlockContextPredictor
//------------------------------------------------------------------------------

WP2Status BlockContextPredictor::Init(uint8_t symbol, uint32_t symbol_range,
                                      uint32_t width) {
  WP2_CHECK_OK(symbol_range <= kMaxRange, WP2_STATUS_INVALID_PARAMETER);
  symbol_ = symbol;
  range_ = symbol_range;
  WP2_CHECK_ALLOC_OK(map_.resize(SizeBlocks(width) * kMapHeight));
  ctxt_.Init(&map_, SizeBlocks(width), kMapHeight);
  return WP2_STATUS_OK;
}

void BlockContextPredictor::SetValue(const CodedBlock& cb, uint8_t value) {
  ctxt_.Set(cb.blk(), value);
}

WP2Status BlockContextPredictor::CopyFrom(const BlockContextPredictor& other) {
  symbol_ = other.symbol_;
  range_ = other.range_;
  WP2_CHECK_ALLOC_OK(map_.copy_from(other.map_));
  const uint32_t width = (uint32_t)map_.size() / kMapHeight;
  assert(width * kMapHeight == other.map_.size());
  ctxt_.Init(&map_, width, kMapHeight);
  // Nothing else to copy in 'ctxt_' as it only points to 'map_'.
  return WP2_STATUS_OK;
}

void BlockContextPredictor::WriteValue(const CodedBlock& cb, uint8_t value,
                                       SymbolManager* const sm,
                                       ANSEncBase* const enc, WP2_OPT_LABEL) {
  sm->Process(symbol_, GetSlot(cb, value), label, enc);
}

uint8_t BlockContextPredictor::ReadValue(const CodedBlock& cb,
                                         SymbolReader* const sr,
                                         WP2_OPT_LABEL) {
  const uint32_t slot = sr->Read(symbol_, label);
  uint8_t context[kMaxNumSegments];
  GetContext(cb, context);
  return context[slot];
}

void BlockContextPredictor::GetContext(const CodedBlock& cb,
                                       uint8_t* const context) const {
  uint8_t context_tmp[kMaxContextSize];
  int8_t y_left_occ, y_right_occ;
  cb.GetOccupancy(&y_left_occ, &y_right_occ, /*top_context_extent=*/nullptr);
  const uint32_t size =
      ctxt_.GetContext(cb.blk(), y_left_occ, y_right_occ, context_tmp);
  // we use the lower 8bit to store the id, the upper 8 to store the counts
  uint16_t counts[kMaxNumSegments] = {0};
  for (uint32_t i = 0; i < size; ++i) counts[context_tmp[i]] += (1 << 8);
  for (uint32_t v = 0; v < range_; ++v) counts[v] |= v;
  std::sort(counts, counts + range_, std::greater<int>());
  for (uint32_t i = 0; i < range_; ++i) {
    context[i] = counts[i] & 0xff;  // retrieve the id
  }
}

uint32_t BlockContextPredictor::GetSlot(const CodedBlock& cb,
                                        uint8_t id) const {
  uint8_t id_ctxt[kMaxRange];
  GetContext(cb, id_ctxt);
  for (uint32_t slot = 0; slot < range_; ++slot) {
    if (id == id_ctxt[slot]) return slot;
  }
  assert(0);  // shouldn't be reached
  return 0;
}

//------------------------------------------------------------------------------
// SegmentIdPredictor
//------------------------------------------------------------------------------

void SegmentIdPredictor::InitInitialSegmentId(const CodedBlock& cb,
                                              uint32_t id) {
  predictor_.SetValue(cb, id);
}

WP2Status SegmentIdPredictor::CopyFrom(const SegmentIdPredictor& other) {
  num_segments_ = other.num_segments_;
  explicit_segment_ids_ = other.explicit_segment_ids_;
  WP2_CHECK_STATUS(predictor_.CopyFrom(other.predictor_));
  return WP2_STATUS_OK;
}

//------------------------------------------------------------------------------
// Writing

WP2Status SegmentIdPredictor::InitWrite(uint32_t num_segments,
                                        bool explicit_segment_ids,
                                        uint32_t width) {
  WP2_CHECK_STATUS(predictor_.Init(kSymbolSegmentId, num_segments, width));
  num_segments_ = num_segments;
  explicit_segment_ids_ = explicit_segment_ids;
  return WP2_STATUS_OK;
}

WP2Status SegmentIdPredictor::WriteHeader(ANSEncBase* const enc) {
  return WP2_STATUS_OK;
}

void SegmentIdPredictor::WriteId(const CodedBlock& cb, SymbolManager* const sm,
                                 ANSEncBase* const enc) {
  if (explicit_segment_ids_ && !cb.blk().IsSmall()) {
    predictor_.WriteValue(cb, cb.id_, sm, enc, "segment_id");
  } else {
    assert(cb.id_ == GetIdFromSize(num_segments_, cb.dim()));
  }
  // The following is still needed, even if RecordId() did it already: the
  // map in segment_id_ctx_ is only partially stored.
  predictor_.SetValue(cb, cb.id_);  // Fill the segment ids.
}

//------------------------------------------------------------------------------
// Reading

WP2Status SegmentIdPredictor::ReadHeader(ANSDec* const dec,
                                         uint32_t num_segments,
                                         bool explicit_segment_ids,
                                         uint32_t width) {
  (void)dec;
  num_segments_ = num_segments;
  explicit_segment_ids_ = explicit_segment_ids;
  WP2_CHECK_STATUS(predictor_.Init(kSymbolSegmentId, num_segments, width));
  // For now, we consider the statistics as being the same for the whole image.
  return WP2_STATUS_OK;
}

void SegmentIdPredictor::ReadId(SymbolReader* const sr, CodedBlock* const cb) {
  if (explicit_segment_ids_ && !cb->blk().IsSmall()) {
    cb->id_ = predictor_.ReadValue(*cb, sr, "segment_id");
  } else {
    cb->id_ = GetIdFromSize(num_segments_, cb->dim());
  }
  assert(cb->id_ < num_segments_);
  predictor_.SetValue(*cb, cb->id_);  // Fill the segment ids.
}

uint8_t SegmentIdPredictor::GetIdFromSize(uint32_t num_segments,
                                          BlockSize block_size) {
  // Map [4x4:32x32] to [0:3].
  const uint8_t segment_id = (uint8_t)WP2Log2Floor(
      (BlockWidth[block_size] + BlockHeight[block_size]) / 2);
  // Scale it to the actual number of segments.
  return DivRound<uint8_t>(segment_id * (num_segments - 1), 3);
}

//------------------------------------------------------------------------------
// DCPredictor
//------------------------------------------------------------------------------

WP2Status DCPredictor::InitMap(uint32_t width, uint32_t) {
  const uint32_t map_height = kMaxBlockSize + 1;
  WP2_CHECK_ALLOC_OK(map_.resize(SizeBlocks(width) * map_height));
  dc_ctxt_.Init(&map_, SizeBlocks(width), map_height);
  return WP2_STATUS_OK;
}

//------------------------------------------------------------------------------
// Writing

WP2Status DCPredictor::InitWrite(ANSDictionaries* const dicts,
                                 uint32_t prefix_bits) {
  prefix_bits_ = prefix_bits;
  return WP2_STATUS_OK;
}

uint16_t DCPredictor::GetPrefix(int16_t dc) const {
  assert(dc > 0 && dc < (1 << prefix_bits_));
  uint16_t n = 0;
  while ((dc >> n) != 1u) ++n;
  return n;
}

WP2Status DCPredictor::WriteHeader(uint32_t cluster, uint32_t max_count,
                                   const SymbolRecorder& recorder,
                                   ANSEncBase* const enc,
                                   ANSDictionaries* const dicts,
                                   SymbolWriter* const sw) {
  // For now, we consider the statistics as being the same for the whole image.
  assert(cluster == 0);
  WP2_CHECK_STATUS(sw->WriteHeader(cluster, max_count, kSymbolDCDict, recorder,
                                   "dc_dict", enc, dicts));
  return WP2_STATUS_OK;
}

void DCPredictor::Record(const CodedBlock& cb, const DCValues& DCs,
                         SymbolRecorder* const recorder) {
  ANSEncNoop enc;
  Write(cb, DCs, &enc, recorder);
}

void DCPredictor::Write(const CodedBlock& cb, const DCValues& DCs,
                        ANSEncBase* const enc, SymbolManager* const sm) {
  DCValues dc_preds;
  GetDCPrediction(cb, &dc_preds);
  for (uint32_t c = 0; c < 3; ++c) {
    int16_t v = DCs.DCs[c] - dc_preds.DCs[c];
    if (enc->PutABit(v == 0, &is_zero_, "C0-dc_zero")) continue;
    const bool sign = (v < 0);
    v = sign ? -v : v;
    const int16_t n = GetPrefix(v);
    sm->Process(kSymbolDCDict, n, "C0-dc_prefix", enc);
    if (n > 0) {
      enc->PutUValue(v & ((1 << n) - 1), n, "C0-dc_abs_val");
    }
    enc->PutBool(sign, "C0-dc_sign");
  }
  dc_ctxt_.Set(cb.blk(), DCs);
}

//------------------------------------------------------------------------------
// Reading

WP2Status DCPredictor::InitRead(uint32_t cluster, SymbolReader* const sr,
                                uint32_t prefix_bits, uint32_t max_count) {
  prefix_bits_ = prefix_bits;
  assert(cluster == 0);
  return sr->ReadHeader(cluster, max_count, kSymbolDCDict, "dc_dict");
}

void DCPredictor::Read(ANSDec* const dec, SymbolReader* const sr,
                       const CodedBlock& cb, DCValues* const DCs) {
  GetDCPrediction(cb, DCs);
  for (uint32_t c = 0; c < 3; ++c) {
    if (!dec->ReadABit(&is_zero_, "C0-dc_zero")) {
      const uint16_t n = sr->Read(kSymbolDCDict, "C0-dc_prefix");
      int16_t v = 1 << n;
      if (n > 0) {
        v |= dec->ReadUValue(n, "C0-dc_abs_val");
      }
      if (dec->ReadBool("C0-dc_sign")) v = -v;
      DCs->DCs[c] += v;
    }
  }
  dc_ctxt_.Set(cb.blk(), *DCs);
}

//------------------------------------------------------------------------------

void DCPredictor::GetDCPrediction(const CodedBlock& cb,
                                  DCValues* const values) {
  DCValues context[kMaxContextSize];
  int8_t y_left_occ, y_right_occ;
  cb.GetOccupancy(&y_left_occ, &y_right_occ, /*top_context_extent=*/nullptr);
  const uint32_t size =
      dc_ctxt_.GetContext(cb.blk(), y_left_occ, y_right_occ, context);
  // Get Average DC values.
  // TODO(skal): try median value as predictor too. Or else.
  values->DCs[0] = 0;
  values->DCs[1] = 0;
  values->DCs[2] = 0;
  for (size_t i = 0; i < size; ++i) {
    values->DCs[0] += context[i].DCs[0];
    values->DCs[1] += context[i].DCs[1];
    values->DCs[2] += context[i].DCs[2];
  }
  values->DCs[0] /= size;
  values->DCs[1] /= size;
  values->DCs[2] /= size;
}

//------------------------------------------------------------------------------
// LargeModePredictor

WP2Status YModePredictor::InitMap(uint32_t width, uint32_t) {
  const uint32_t map_width = SizeBlocks(width);
  WP2_CHECK_ALLOC_OK(map_.resize(map_width * kMapHeight));
  predictor_ctxt_.Init(&map_, map_width, kMapHeight);
  Reset();
  return WP2_STATUS_OK;
}

WP2Status YModePredictor::CopyFrom(const YModePredictor& other) {
  WP2_CHECK_ALLOC_OK(map_.copy_from(other.map_));
  const uint32_t map_width = (uint32_t)map_.size() / kMapHeight;
  assert(map_width * kMapHeight == other.map_.size());
  predictor_ctxt_.Init(&map_, map_width, kMapHeight);
  // Nothing else to copy in 'predictor_ctxt_' as it only points to 'map_'.
  return WP2_STATUS_OK;
}

void YModePredictor::Reset() { std::fill(map_.begin(), map_.end(), -1); }

static uint32_t GetCluster(const ContextImpl<int8_t>& ctxt,
                           const CodedBlock& cb) {
  const int above_mode = (cb.y() > 0) ? ctxt.Get(cb.x(), cb.y() - 1) : 0;
  const int left_mode = (cb.x() > 0) ? ctxt.Get(cb.x() - 1, cb.y()) : 0;
  return (above_mode == left_mode && above_mode >= 0)
             ? above_mode
             : (uint32_t)Predictor::ModeContext::kNum;
}

void YModePredictor::Write(const CodedBlock& cb, const Predictor& pred,
                           bool do_update, ANSEncBase* const enc,
                           SymbolManager* const sw) {
  const uint32_t cluster = GetCluster(predictor_ctxt_, cb);
  sw->Process(cluster, kSymbolModeY, pred.mode(), "predictor", enc);
  if (do_update) {
    const Predictor::ModeContext ctxt = pred.mode_context();
    predictor_ctxt_.Set(cb.blk(), (int)ctxt);
  }
}

uint32_t YModePredictor::Read(const CodedBlock& cb, const Predictors& preds,
                              SymbolReader* const sr) {
  const uint32_t cluster = GetCluster(predictor_ctxt_, cb);
  const uint32_t mode = sr->Read(cluster, kSymbolModeY, "predictor");
  const Predictor::ModeContext ctxt =
      preds.GetFirstWithMode(mode)->mode_context();
  predictor_ctxt_.Set(cb.blk(), (int)ctxt);
  return mode;
}

//------------------------------------------------------------------------------
// AlphaModePredictor
//------------------------------------------------------------------------------

WP2Status AlphaModePredictor::Init(uint32_t width, uint32_t height) {
  return predictor_.Init(kSymbolBlockAlphaMode, /*symbol_range=*/3, width);
}

WP2Status AlphaModePredictor::CopyFrom(const AlphaModePredictor& other) {
  return predictor_.CopyFrom(other.predictor_);
}

void AlphaModePredictor::Update(const CodedBlock& cb) {
  predictor_.SetValue(cb, (uint8_t)cb.alpha_mode_);
}

void AlphaModePredictor::Write(const CodedBlock& cb, ANSEncBase* const enc,
                               SymbolManager* const sw) {
  predictor_.WriteValue(cb, (uint8_t)cb.alpha_mode_, sw, enc, "alpha_mode");
  Update(cb);
}

BlockAlphaMode AlphaModePredictor::Read(const CodedBlock& cb,
                                        SymbolReader* const sr) {
  const uint32_t mode = predictor_.ReadValue(cb, sr, "alpha_mode");
  predictor_.SetValue(cb, (int)mode);
  return (BlockAlphaMode)mode;
}

//------------------------------------------------------------------------------
// DCDiffusionMap

WP2Status DCDiffusionMap::Init(uint32_t width) {
  const uint32_t bw = SizeBlocks(width);
  for (uint32_t i = 0; i < kErrorHeight; ++i) {
    WP2_CHECK_ALLOC_OK(error_[i].resize(bw));
  }
  Clear();
  return WP2_STATUS_OK;
}

WP2Status DCDiffusionMap::CopyFrom(const DCDiffusionMap& other) {
  for (uint32_t i = 0; i < kErrorHeight; ++i) {
    WP2_CHECK_ALLOC_OK(error_[i].copy_from(other.error_[i]));
  }
  max_row_ = other.max_row_;
  return WP2_STATUS_OK;
}

void DCDiffusionMap::Clear() {
  max_row_ = 0;
  for (auto& l : error_[0]) l = 0.f;
}

void DCDiffusionMap::Store(const FrontMgrBase& mgr, const Block& blk,
                           int16_t error) {
  const uint32_t block_max_row = blk.y() + blk.h();
  if (block_max_row > max_row_) {
    for (uint32_t i = max_row_ + 1; i <= block_max_row; ++i) {
      for (auto& l : error_[i % kErrorHeight]) l = 0.f;
    }
    max_row_ = block_max_row;
  }
  const uint32_t left_occupancy =
      (blk.x() > 0) ? std::min(mgr.GetOccupancy(blk.x() - 1), block_max_row)
                    : block_max_row;
  const uint32_t width = error_[0].size();
  const uint32_t right_occupancy =
      (blk.x() + blk.w() < width)
          ? std::min(mgr.GetOccupancy(blk.x() + blk.w()), block_max_row)
          : block_max_row;

  const uint32_t border_size = blk.w() + (block_max_row - left_occupancy) +
                               (block_max_row - right_occupancy);
  const float error_per_block = (float)error / border_size;
  for (uint32_t y = left_occupancy; y < block_max_row; ++y) {
    error_[y % kErrorHeight][blk.x() - 1] += error_per_block;
  }
  for (uint32_t y = right_occupancy; y < block_max_row; ++y) {
    error_[y % kErrorHeight][blk.x() + blk.w()] += error_per_block;
  }
  for (uint32_t x = blk.x(); x < blk.x() + blk.w(); ++x) {
    error_[block_max_row % kErrorHeight][x] += error_per_block;
  }
}

int16_t DCDiffusionMap::Get(const Block& blk, uint32_t strength) const {
  float new_error = 0;
  for (uint32_t y = blk.y(); y < (blk.y() + blk.h()) && y <= max_row_; ++y) {
    new_error += error_[y % kErrorHeight][blk.x()];
    if (blk.w() > 1) {
      new_error += error_[y % kErrorHeight][blk.x() + blk.w() - 1];
    }
  }
  for (uint32_t x = blk.x() + 1; x < (blk.x() + blk.w()) - 1; ++x) {
    new_error += error_[blk.y() % kErrorHeight][x];
  }
  return std::round(new_error * strength / 255.f);
}

}   // namespace WP2
