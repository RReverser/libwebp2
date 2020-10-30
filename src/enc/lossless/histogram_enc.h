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
// Author: Vincent Rabaud (vrabaud@google.com)
//
// Models the histograms of literal and distance codes.

#ifndef WP2_ENC_LOSSLESS_HISTOGRAM_ENC_H_
#define WP2_ENC_LOSSLESS_HISTOGRAM_ENC_H_

#include <algorithm>
#include <cstring>

#include "src/common/symbols.h"
#include "src/enc/lossless/backward_references_enc.h"
#include "src/utils/quantizer.h"
#include "src/wp2/format_constants.h"

namespace WP2L {

// Class that abstracts the storage of counts for all the existing symbols.
template <typename T>
class CountsBuffer {
 public:
  CountsBuffer() noexcept { Init(); }
  ~CountsBuffer() { Reset(); }
  void Init() noexcept {
    buffer_ = nullptr;
    has_own_buffer_ = false;
  }
  // Allocate the buffer and set it to 0.
  WP2_NO_DISCARD
  bool Allocate(const LosslessSymbolsInfo* const symbols_info) {
    Reset();
    symbols_info_ = symbols_info;
    buffer_ = (T*)WP2Calloc(symbols_info_->MaxRangeSum(), sizeof(*buffer_));
    if (buffer_ == nullptr) return false;
    has_own_buffer_ = true;
    InitPointers();
    return true;
  }
  // Reset the counts to 0.
  void Clear() {
    memset(buffer_, 0, sizeof(*buffer_) * symbols_info_->MaxRangeSum());
  }

  void CopyTo(CountsBuffer<T>* const dst) const {
    memcpy(dst->buffer_, buffer_,
           sizeof(*buffer_) * symbols_info_->MaxRangeSum());
    dst->InitPointers();
  }

  void SetBuffer(T* buffer, const LosslessSymbolsInfo& symbols_info) {
    symbols_info_ = &symbols_info;
    Reset();
    assert(buffer != nullptr);
    buffer_ = buffer;
    InitPointers();
  }

  T* counts_[kSymbolNum];
  const LosslessSymbolsInfo* symbols_info_;

 private:
  void Reset() {
    if (has_own_buffer_) {
      WP2Free(buffer_);
      buffer_ = nullptr;
      has_own_buffer_ = false;
    }
  }
  // Set the pointers for the counts of each histogram to the right spots
  // in the overall buffer_.
  void InitPointers() {
    counts_[0] = buffer_;
    for (size_t i = 1; i < kSymbolNum; ++i) {
      counts_[i] = counts_[i - 1] + symbols_info_->GetMaxRange((Symbol)(i - 1));
    }
  }
  bool has_own_buffer_;
  T* buffer_;  // The buffer contains all the symbol histograms one after
               // another.
};

// A simple container for histograms of data.
class Histogram : public CountsBuffer<uint32_t> {
 public:
  Histogram() noexcept : CountsBuffer<uint32_t>() { Init(); }
  void CopyTo(Histogram* const dst) const {
    super::CopyTo(dst);
    memcpy(dst->trivial_symbols_, trivial_symbols_, sizeof(trivial_symbols_));
    memcpy(dst->nonzeros_, nonzeros_, sizeof(nonzeros_));
    dst->bit_cost_ = bit_cost_;
    memcpy(dst->costs_, costs_, sizeof(costs_));
  }
  void Init() noexcept {
    super::Init();
    memset(trivial_symbols_, 0, sizeof(trivial_symbols_));
    memset(nonzeros_, 0, sizeof(nonzeros_));
    bit_cost_ = 0.;
    memset(costs_, 0, sizeof(costs_));
  }
  void Clear() {
    super::Clear();
    memset(trivial_symbols_, 0, sizeof(trivial_symbols_));
    memset(nonzeros_, 0, sizeof(nonzeros_));
    memset(costs_, 0, sizeof(costs_));
  }

  uint32_t trivial_symbols_[kSymbolNum];
  uint32_t nonzeros_[kSymbolNum];
  double bit_cost_;           // cached value of bit cost.
  double costs_[kSymbolNum];  // cached values of entropy costs.

 private:
  typedef CountsBuffer<uint32_t> super;
};

// Collection of histograms with fixed capacity, allocated as one
// big memory chunk.
class HistogramSet {
 public:
  // Allocate an array of pointer to histograms, allocated and initialized
  // using the number of bits of 'symbols_info'.
  WP2_NO_DISCARD
  bool Allocate(uint32_t size, const LosslessSymbolsInfo& symbols_info);

  Histogram* tmp_histo_;
  // histograms sharing the same buffer
  WP2::VectorNoCtor<Histogram*> histograms_;
  const LosslessSymbolsInfo* symbols_info_ = nullptr;

 private:
  WP2::Vector<Histogram> histo_buffer_;  // internal buffer of histograms
  WP2::Vector_u32 buffer_;               // common buffer to all histograms
};

// Fills the histograms initialized in HistogramSet from the backward
// references. If there is only one histogram, the parameters 'width',
// 'height'and 'histo_bits' do not matter.
WP2Status BuildHistograms(uint32_t width, uint32_t height, uint32_t histo_bits,
                          const LosslessSymbolsInfo& symbols_info,
                          const BackwardRefs& refs,
                          HistogramSet* const histogram_set);

// Builds the histogram image.
WP2Status GetHistoImageSymbols(uint32_t width, uint32_t height,
                               const BackwardRefs& refs, int speed,
                               int histogram_bits,
                               const LosslessSymbolsInfo& symbols_info,
                               HistogramSet* const image_histo,
                               uint16_t* const histogram_symbols);

// Class comparing histograms for the sake of merging them.
// It contains several buffers that speed-up the computation of the different
// entropies.
class PopulationAnalyzer {
 public:
  explicit PopulationAnalyzer(uint32_t image_size);

  WP2Status Allocate(const LosslessSymbolsInfo& symbols_info) {
    const size_t histogram_size_max_ = symbols_info.GetMaxRange();
    WP2_CHECK_ALLOC_OK(histogram_.resize(histogram_size_max_));
    WP2_CHECK_ALLOC_OK(mapping_.resize(histogram_size_max_));
    WP2_CHECK_ALLOC_OK(sum_tmp_.resize(histogram_size_max_));
    WP2_CHECK_STATUS(quantizer_.Allocate(histogram_size_max_));
    return WP2_STATUS_OK;
  }
  void PopulationCost(const uint32_t* const population, uint32_t range,
                      uint32_t* const nonzeros, uint32_t* const trivial_sym,
                      double* cost);

  void GetCombinedEntropy(const uint32_t* const X, const uint32_t* const Y,
                          uint32_t X_nonzeros, uint32_t Y_nonzeros,
                          uint32_t X_trivial, uint32_t Y_trivial,
                          uint32_t length, double* cost);

  // Various histogram combine/cost-eval functions

  int GetCombinedHistogramEntropy(const Histogram& a, const Histogram& b,
                                  double cost_threshold, double* cost);

  // Performs out = a + b, computing the cost C(a+b) - C(a) - C(b) while
  // comparing to the threshold value 'cost_threshold'. The score returned is
  //  Score = C(a+b) - C(a) - C(b), where C(a) + C(b) is known and fixed.
  // Since the previous score passed is 'cost_threshold', we only need to
  // compare the partial cost against 'cost_threshold + C(a) + C(b)' to possibly
  // bail-out early.
  double HistogramAddEval(const Histogram& a, const Histogram& b,
                          Histogram* const out, double cost_threshold);

  // Same as HistogramAddEval(), except that the resulting histogram
  // is not stored. Only the cost C(a+b) - C(a) is evaluated. We omit
  // the term C(b) which is constant over all the evaluations.
  double HistogramAddThresh(const Histogram& a, const Histogram& b,
                            double cost_threshold);

  // Update the individual/overall bit costs of an histogram.
  void UpdateHistogramCost(Histogram* const h);

 private:
  WP2::Vector_u32 sum_tmp_;

  // Variables used in most computations.
  WP2::Vector_u32 histogram_;
  WP2::Vector_u16 mapping_;
  WP2::Quantizer quantizer_;
  uint32_t image_size_;
};

// Pair of histograms used for histogram merges.
typedef struct {
  uint16_t idx1;
  uint16_t idx2;
  // Cost difference between the merged histogram cost and the individual costs
  // for idx1 and idx2.
  double cost_diff;
  // Cost of the histogram merging the histograms at idx1 and idx2.
  double cost_combo;
} HistogramPair;

// Queue storing pairs of symbol histograms and maintaining the best one (with
// minimal cost) at the front. It is not a traditional queue because when the
// front is popped, the two histograms are merged and all pairs using one of the
// two has to recompute its cost.
// TODO(vrabaud) Keep track of the analyzed pairs to avoid duplicates and
//               re-analyzing pairs that give worse results.
class HistoQueue {
 public:
  bool empty() const { return queue_.empty(); }
  const HistogramPair& front() const { return queue_.front(); }
  bool full() const { return queue_.size() == queue_.capacity(); }
  // Returns the list of valid remaining valid histograms (the ones that have
  // not been merged).
  const WP2::VectorNoCtor<uint16_t>& ValidIndices() const { return indices_; }

  // Initializes the queue to be used with a given set of histograms using a
  // given 'analyzer'
  WP2Status Init(WP2::VectorNoCtor<Histogram*>* const histograms,
                 PopulationAnalyzer* const analyzer);
  // Resizes the queue: the bigger the size the more flexibility, but the
  // slower it gets.
  WP2Status Resize(uint32_t size);

  // Merges the histograms from the front (idx1 and idx2), invalidate idx2 and
  // also checks if any pair that dealt with idx1 or idx2 is still worth it (by
  // updating the cost or removing the pair).
  void MergeFront();

  // Creates a pair from indices "idx1" and "idx2" provided its cost
  // is inferior to "threshold", a negative entropy.
  // It returns the cost of the pair, or 0. if it superior to threshold.
  WP2Status Push(uint16_t idx1, uint16_t idx2, double threshold,
                 double* const cost_diff = nullptr);

  // Removes histograms that are not part of the valid indices. Also empties the
  // queue as elements become invalid. The queue needs to be resized to be used
  // again.
  WP2Status RemoveInvalidHistograms();

 private:
  // Checks whether a pair should be updated as head or not.
  void MaybeUpdateHead(HistogramPair* const pair);

  WP2::VectorNoCtor<HistogramPair> queue_;
  // List of valid indices (idx2 of pairs that have been merged are removed).
  WP2::VectorNoCtor<uint16_t> indices_;
  PopulationAnalyzer* analyzer_;
  WP2::VectorNoCtor<Histogram*>* histograms_;
};

}  // namespace WP2L

#endif  // WP2_ENC_LOSSLESS_HISTOGRAM_ENC_H_
