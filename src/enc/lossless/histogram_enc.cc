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
// Author: Jyrki Alakuijala (jyrki@google.com)
//

#include "src/enc/lossless/histogram_enc.h"

#include <functional>
#include <limits>
#include <numeric>

#ifdef HAVE_CONFIG_H
#include "src/webp/config.h"
#endif

#include "src/common/color_precision.h"
#include "src/dsp/lossless/lossless.h"
#include "src/dsp/lossless/lossless_common.h"
#include "src/enc/lossless/losslessi_enc.h"
#include "src/enc/symbols_enc.h"
#include "src/utils/random.h"
#include "src/utils/vector.h"

namespace WP2L {

// Maximum number of histograms allowed in greedy combining algorithm.
#define MAX_HISTO_GREEDY 100

//------------------------------------------------------------------------------

static inline void CopySubHistogram(const Histogram& a, const size_t s,
                                    Histogram* const out) {
  const size_t size =
      sizeof(out->counts_[s][0]) * a.symbols_info_->GetMaxRange(s);
  memcpy(out->counts_[s], a.counts_[s], size);
  out->trivial_symbols_[s] = a.trivial_symbols_[s];
  out->nonzeros_[s] = a.nonzeros_[s];
}

// Add two histograms, taking into account the number of non-zero elements for
// speed-up.
static void HistogramAdd(const Histogram& a, const Histogram& b,
                         Histogram* const out) {
  assert(a.symbols_info_ == b.symbols_info_);
  if (&b != out) {
    for (size_t s = 0; s < kSymbolNum; ++s) {
      if (a.nonzeros_[s] == 0) {
        CopySubHistogram(b, s, out);
      } else if (b.nonzeros_[s] == 0) {
        CopySubHistogram(a, s, out);
      } else if (a.nonzeros_[s] == 1 && b.nonzeros_[s] == 1 &&
                 a.trivial_symbols_[s] == b.trivial_symbols_[s]) {
        const auto& trivial_symbol = a.trivial_symbols_[s];
        const size_t size =
            sizeof(out->counts_[s][0]) * a.symbols_info_->GetMaxRange(s);
        memset(out->counts_[s], 0, size);
        out->counts_[s][trivial_symbol] =
            a.counts_[s][trivial_symbol] + b.counts_[s][trivial_symbol];
        out->trivial_symbols_[s] = trivial_symbol;
        out->nonzeros_[s] = 1;
      } else {
        out->trivial_symbols_[s] = 0;
        WP2LBufferAdd(a.counts_[s], b.counts_[s],
                      a.symbols_info_->GetMaxRange(s), out->counts_[s],
                      &out->nonzeros_[s]);
        assert(out->nonzeros_[s] > 1);
      }
    }
  } else {
    for (size_t s = 0; s < kSymbolNum; ++s) {
      if (a.nonzeros_[s] == 0) {
        // out is b and no modification is added.
      } else if (b.nonzeros_[s] == 0) {
        CopySubHistogram(a, s, out);
      } else if (a.nonzeros_[s] == 1 && b.nonzeros_[s] == 1 &&
                 a.trivial_symbols_[s] == b.trivial_symbols_[s]) {
        // a and b are basically the same.
        out->counts_[s][a.trivial_symbols_[s]] +=
            a.counts_[s][a.trivial_symbols_[s]];
      } else {
        out->trivial_symbols_[s] = 0;
        WP2LBufferAdd(a.counts_[s], b.counts_[s],
                      a.symbols_info_->GetMaxRange(s), out->counts_[s],
                      &out->nonzeros_[s]);
        assert(out->nonzeros_[s] > 1);
      }
    }
  }
}

bool HistogramSet::Allocate(uint32_t size,
                            const LosslessSymbolsInfo& symbols_info) {
  symbols_info_ = &symbols_info;
  const uint32_t cache_num = symbols_info.MaxRangeSum();
  if (!histograms_.resize(size) || !histo_buffer_.resize(size + 1) ||
      !buffer_.resize((size + 1) * cache_num)) {
    return false;
  }
  std::fill(buffer_.begin(), buffer_.end(), 0);
  uint32_t* memory = buffer_.data();
  for (size_t i = 0; i < size; ++i) {
    histograms_[i] = &histo_buffer_[i];
    histo_buffer_[i].Init();
    histo_buffer_[i].SetBuffer(memory, symbols_info);
    memory += cache_num;
  }
  tmp_histo_ = &histo_buffer_[size];
  tmp_histo_->Init();
  tmp_histo_->SetBuffer(memory, symbols_info);
  return true;
}

// -----------------------------------------------------------------------------
// Entropy-related functions.

PopulationAnalyzer::PopulationAnalyzer(uint32_t image_size)
    : image_size_(image_size) {
  WP2LEncDspInit();
}

void PopulationAnalyzer::PopulationCost(const uint32_t* const population,
                                        uint32_t range,
                                        uint32_t* const nonzeros,
                                        uint32_t* const trivial_sym,
                                        double* cost) {
  WP2::Quantizer::Config* config;

  uint32_t size = 0;
  for (uint32_t k = 0; k < range; ++k) {
    if (population[k] == 0) continue;
    histogram_[size] = population[k];
    mapping_[size] = k;
    ++size;
  }
  if (nonzeros != nullptr) *nonzeros = size;
  // TODO(vrabaud) replace the cost logic below with a proper call to
  //               SymbolWriter::WriteHeader.
  if (size == 0) {
    if (cost != nullptr) *cost = WP2Log2(3);
    return;
  }
  if (size == 1) {
    if (trivial_sym != nullptr) *trivial_sym = mapping_[0];
    // TODO(vrabaud) do in the Quantizer.
    if (cost != nullptr) *cost = WP2Log2(3) + WP2Log2(range);
    return;
  }

  // TODO(vrabaud) use the speed parameter instead of imposing 0.
  quantizer_.Quantize(histogram_.data(), mapping_.data(), size, range,
                      image_size_, /*speed=*/0, &config);
  if (cost != nullptr) {
    *cost = WP2Log2(3) +
            WP2Log2(std::min(image_size_, range) - 1) +
            config->cost_;
  }
}

void PopulationAnalyzer::GetCombinedEntropy(
    const uint32_t* const X, const uint32_t* const Y, uint32_t X_nonzeros,
    uint32_t Y_nonzeros, uint32_t X_trivial, uint32_t Y_trivial,
    uint32_t length, double* cost) {
  // If one of the histograms is basically unused:
  if (X_nonzeros == 0) {
    PopulationCost(Y, length, nullptr, nullptr, cost);
    return;
  }
  if (Y_nonzeros == 0) {
    PopulationCost(X, length, nullptr, nullptr, cost);
    return;
  }
  // If they both contain the same element.
  if (X_nonzeros == 1 && Y_nonzeros == 1 && X_trivial == Y_trivial) {
    *cost = WP2Log2(3) + WP2Log2(length);
    return;
  }

  WP2LBufferAdd(X, Y, length, sum_tmp_.data(), nullptr);
  PopulationCost(sum_tmp_.data(), length, nullptr, nullptr, cost);
}

// -----------------------------------------------------------------------------
// Various histogram combine/cost-eval functions

int PopulationAnalyzer::GetCombinedHistogramEntropy(const Histogram& a,
                                                    const Histogram& b,
                                                    double cost_threshold,
                                                    double* cost) {
  assert(a.symbols_info_ == b.symbols_info_);

  // Deal with the raw bits of length and distance.
  for (const auto& sym : {kSymbolLen, kSymbolDist}) {
    *cost += WP2LExtraCostCombined(a.counts_[sym], b.counts_[sym],
                                   a.symbols_info_->GetMaxRange(sym));
    if (*cost > cost_threshold) return 0;
  }

  bool is_used[kSymbolTypeNum];
  for (size_t i = 0; i < kSymbolTypeNum; ++i) {
    is_used[i] = (a.counts_[kSymbolType][i] + b.counts_[kSymbolType][i] > 0);
  }

  // Deal with all the symbols.
  for (size_t s = 0; s < kSymbolNum; ++s) {
    if (a.symbols_info_->IsSymbolUseless((Symbol)s, is_used) &&
        b.symbols_info_->IsSymbolUseless((Symbol)s, is_used)) {
      continue;
    }
    double tmp;
    GetCombinedEntropy(a.counts_[s], b.counts_[s], a.nonzeros_[s],
                       b.nonzeros_[s], a.trivial_symbols_[s],
                       b.trivial_symbols_[s], a.symbols_info_->GetMaxRange(s),
                       &tmp);
    *cost += tmp;
    if (*cost > cost_threshold) return 0;
  }
  return 1;
}

// Performs out = a + b, computing the cost C(a+b) - C(a) - C(b) while comparing
// to the threshold value 'cost_threshold'. The score returned is
//  Score = C(a+b) - C(a) - C(b), where C(a) + C(b) is known and fixed.
// Since the previous score passed is 'cost_threshold', we only need to compare
// the partial cost against 'cost_threshold + C(a) + C(b)' to possibly bail-out
// early.
double PopulationAnalyzer::HistogramAddEval(const Histogram& a,
                                            const Histogram& b,
                                            Histogram* const out,
                                            double cost_threshold) {
  double cost = 0;
  const double sum_cost = a.bit_cost_ + b.bit_cost_;
  cost_threshold += sum_cost;

  if (GetCombinedHistogramEntropy(a, b, cost_threshold, &cost)) {
    HistogramAdd(a, b, out);
    out->bit_cost_ = cost;
    assert(a.symbols_info_ == b.symbols_info_);
  }

  return cost - sum_cost;
}

// Same as HistogramAddEval(), except that the resulting histogram
// is not stored. Only the cost C(a+b) - C(a) is evaluated. We omit
// the term C(b) which is constant over all the evaluations.
double PopulationAnalyzer::HistogramAddThresh(const Histogram& a,
                                              const Histogram& b,
                                              double cost_threshold) {
  double cost = -a.bit_cost_;
  GetCombinedHistogramEntropy(a, b, cost_threshold, &cost);
  return cost;
}

// -----------------------------------------------------------------------------

void PopulationAnalyzer::UpdateHistogramCost(Histogram* const h) {
  h->bit_cost_ = 0.;
  bool is_used[kSymbolTypeNum];
  bool has_at_least_one_symbol = false;
  for (uint32_t i = 0; i < kSymbolTypeNum; ++i) {
    is_used[i] = (h->counts_[kSymbolType][i] > 0);
    has_at_least_one_symbol |= is_used[i];
  }
  if (!has_at_least_one_symbol) {
    // If the histogram is completely empty, it contains no information so no
    // cost.
    std::fill(h->costs_, h->costs_ + kSymbolNum, 0.);
    return;
  }
  for (size_t s = 0; s < kSymbolNum; ++s) {
    const Symbol sym = (Symbol)s;
    if (h->symbols_info_->IsSymbolUseless(sym, is_used)) {
      h->costs_[s] = 0.;
    } else {
      PopulationCost(h->counts_[s], h->symbols_info_->GetMaxRange(sym),
                     &h->nonzeros_[s], &h->trivial_symbols_[s],
                     &(h->costs_[s]));
    }
  }
  for (const auto& sym : {kSymbolLen, kSymbolDist}) {
    h->costs_[sym] +=
        WP2LExtraCost(h->counts_[sym], h->symbols_info_->GetMaxRange(sym));
  }
  // Update the total bit_cost.
  for (const auto& c : h->costs_) h->bit_cost_ += c;
}

WP2Status BuildHistograms(uint32_t width, uint32_t height, uint32_t histo_bits,
                          const LosslessSymbolsInfo& symbols_info,
                          const BackwardRefs& refs,
                          HistogramSet* const histogram_set) {
  // Count the symbols using a SymbolRecorder.
  WP2::SymbolRecorder recorder;
  WP2::ANSEncNoop enc_noop;
  LosslessSymbolsInfo symbols_info_counter;
  WP2_CHECK_STATUS(symbols_info_counter.CopyFrom(symbols_info));
  const uint32_t num_histos = histogram_set->histograms_.size();
  symbols_info_counter.SetNumClusters(num_histos);
  WP2_CHECK_STATUS(recorder.Allocate(symbols_info_counter, /*num_records=*/0));
  WP2_CHECK_STATUS(StorePixels(width, height, histo_bits, refs,
                               /*histogram_symbols=*/nullptr, &enc_noop,
                               &recorder));

  // Copy the dictionaries to the HistogramSet.
  Histogram** const histograms = histogram_set->histograms_.data();
  for (uint32_t i = 0; i < num_histos; ++i) {
    for (uint32_t s = 0; s < kSymbolNum; ++s) {
      if (symbols_info_counter.Range(i, s) == 0) continue;
      const WP2::Vector_u32& counts = recorder.GetRecordedDict(i, s).Counts();
      std::copy(counts.begin(), counts.end(), histograms[i]->counts_[s]);
    }
  }
  return WP2_STATUS_OK;
}

// Copies the histograms and computes its bit_cost.
static WP2Status HistogramSetCopy(const HistogramSet& in,
                                  const LosslessSymbolsInfo& symbols_info,
                                  HistogramSet* const out) {
  uint32_t size = 0;
  for (Histogram* const histo : in.histograms_) {
    // Skip this histogram if it is not even used once.
    if (histo->bit_cost_ != 0) ++size;
  }
  WP2_CHECK_ALLOC_OK(out->Allocate(size, symbols_info));

  size = 0;
  for (Histogram* const histo : in.histograms_) {
    // Skip this histogram if it is not even used once.
    if (histo->bit_cost_ == 0) continue;
    // Copy histograms from 'in' to 'out'.
    histo->CopyTo(out->histograms_[size++]);
  }
  return WP2_STATUS_OK;
}

// Class that merges histograms together if they have common entropies for some
// of their symbols.
class BinMerger {
 public:
  BinMerger(uint32_t histo_size, uint32_t speed) : speed_(speed) {
    combine_cost_factor_ = 0.16;
    if (speed < 8) {
      if (histo_size > 256) combine_cost_factor_ /= 2.;
      if (histo_size > 512) combine_cost_factor_ /= 2.;
      if (histo_size > 1024) combine_cost_factor_ /= 2.;
      if (speed < 5) combine_cost_factor_ /= 2.;
    }
  }
  // Merges empty histograms together.
  static WP2Status RemoveEmptyHistograms(HistogramSet* const image_histo) {
    uint32_t new_size = 0;
    auto& histograms = image_histo->histograms_;
    for (uint32_t i = 0; i < histograms.size(); ++i) {
      if (histograms[i]->bit_cost_ != 0.) {
        std::swap(histograms[i], histograms[new_size]);
        ++new_size;
      }
    }
    WP2_CHECK_ALLOC_OK(histograms.resize(new_size));
    return WP2_STATUS_OK;
  }
  // Partition histograms to different entropy bins for three dominant (literal,
  // red and blue) symbol costs and compute the histogram aggregate bit_cost.
  bool AnalyzeEntropyBin(const HistogramSet& image_histo,
                         uint16_t* const bin_map) {
    // Don't attempt linear bin-partition heuristic for
    // histograms of small sizes (as bin_map will be very sparse) and
    // maximum speed 9 (to preserve the compression gains at that level).
    const uint32_t entropy_combine_num_bins =
        speed_ == 0 ? kNumPartitions : kBinSizeTot;
    const bool entropy_combine =
        (image_histo.histograms_.size() > entropy_combine_num_bins * 2) &&
        (speed_ < 9);
    if (!entropy_combine) return false;

    std::fill(max_costs_, max_costs_ + kSymbolNum, 0);
    std::fill(min_costs_, min_costs_ + kSymbolNum,
              std::numeric_limits<double>::max());

    // Analyze the entropy cost bounds.
    for (const Histogram* h : image_histo.histograms_) {
      for (uint32_t s = 0; s < h->symbols_info_->Size(); ++s) {
        if (h->costs_[s] < min_costs_[s]) {
          min_costs_[s] = h->costs_[s];
        } else if (h->costs_[s] > max_costs_[s]) {
          max_costs_[s] = h->costs_[s];
        }
      }
    }

    // Find the widest cost disparities. Variance and coefficient of variation
    // are actually worse. PCA should be tried.
    for (uint32_t i = 0; i < image_histo.symbols_info_->Size(); ++i) {
      disparities_[i].first = max_costs_[i] - min_costs_[i];
      disparities_[i].second = i;
    }
    std::sort(disparities_, disparities_ + kSymbolNum,
              std::greater<std::pair<double, uint32_t>>());

    // bin-hash histograms on three of the dominant (literal, red and blue)
    // symbol costs and store the resulting bin_id for each histogram.
    uint32_t i = 0;
    for (const Histogram* h : image_histo.histograms_) {
      bin_map[i++] = GetHistoBinIndex(*h);
    }
    return true;
  }

  // Compact image_histo[] by merging some histograms with same bin_id together
  // if it's advantageous.
  WP2Status CombineEntropyBin(const WP2::Vector_u16& bin_map,
                              PopulationAnalyzer* const analyzer,
                              HistogramSet* const image_histo) const {
    Histogram** const histograms = image_histo->histograms_.data();
    // Work in-place: processed histograms are put at the beginning of
    // image_histo[]. At the end, we just have to truncate the array.
    struct {
      int16_t first;  // position of the histogram that accumulates all
                      // histograms with the same bin_id
      uint16_t num_combine_failures;  // number of combine failures per bin_id
    } bin_info[kBinSizeTot];

    const uint32_t num_bins = (speed_ == 0) ? kNumPartitions : kBinSizeTot;
    assert(num_bins <= kBinSizeTot);
    for (uint32_t idx = 0; idx < num_bins; ++idx) {
      bin_info[idx].first = -1;
      bin_info[idx].num_combine_failures = 0;
    }

    uint32_t size = 0;
    for (uint32_t idx = 0; idx < image_histo->histograms_.size(); ++idx) {
      const int bin_id = bin_map[idx];
      const int first = bin_info[bin_id].first;
      assert(size <= idx);
      if (first == -1) {
        // just move histogram #idx to its final position
        histograms[size] = histograms[idx];
        bin_info[bin_id].first = size++;
      } else if (speed_ == 0) {
        HistogramAdd(*histograms[idx], *histograms[first], histograms[first]);
      } else {
        // try to merge #idx into #first (both share the same bin_id)
        const double bit_cost = histograms[idx]->bit_cost_;
        const double bit_cost_thresh = -bit_cost * combine_cost_factor_;
        const double curr_cost_diff = analyzer->HistogramAddEval(
            *histograms[first], *histograms[idx], image_histo->tmp_histo_,
            bit_cost_thresh);
        if (curr_cost_diff < bit_cost_thresh) {
          std::swap(image_histo->tmp_histo_, histograms[first]);
        } else {
          histograms[size++] = histograms[idx];
        }
      }
    }
    WP2_CHECK_ALLOC_OK(image_histo->histograms_.resize(size));
    if (speed_ == 0) {
      // for low_effort case, update the final cost when everything is merged
      for (uint32_t idx = 0; idx < size; ++idx) {
        analyzer->UpdateHistogramCost(histograms[idx]);
      }
    }
    return WP2_STATUS_OK;
  }

 private:
  uint32_t GetHistoBinIndex(const Histogram& h) const {
    uint32_t bin_id = 0;
    // Carve the space depending on the three widest disparities.
    for (uint32_t i = 0; i < (speed_ > 0 ? 3u : 1u); ++i) {
      const uint32_t s = disparities_[i].second;
      const double range = disparities_[i].first;
      if (range == 0.) break;
      // Fix potential rounding issues.
      const uint32_t new_bin_id = std::min(
          (uint32_t)(kNumPartitions * (h.costs_[s] - min_costs_[s]) / range),
          kNumPartitions - 1);
      bin_id = bin_id * kNumPartitions + new_bin_id;
      assert(bin_id < kBinSizeTot);
    }
    return bin_id;
  }

  // Number of partitions for the dominant symbols.
  static constexpr uint32_t kNumPartitions = 4;
  // The size of the bin-hash corresponding to the three dominant costs.
  static constexpr uint32_t kBinSizeTot =
      (kNumPartitions * kNumPartitions * kNumPartitions);

  double min_costs_[kSymbolNum];
  double max_costs_[kSymbolNum];
  std::pair<double, uint32_t> disparities_[kSymbolNum];

  double combine_cost_factor_;
  const uint32_t speed_;
};

// -----------------------------------------------------------------------------
// Histogram pairs priority queue

WP2Status HistoQueue::Init(WP2::VectorNoCtor<Histogram*>* const histograms,
                           PopulationAnalyzer* const analyzer) {
  histograms_ = histograms;
  analyzer_ = analyzer;
  return WP2_STATUS_OK;
}

WP2Status HistoQueue::Resize(uint32_t size) {
  WP2_CHECK_ALLOC_OK(queue_.reserve(size));
  WP2_CHECK_ALLOC_OK(indices_.resize(histograms_->size()));
  std::iota(indices_.begin(), indices_.end(), 0u);
  return WP2_STATUS_OK;
}

void HistoQueue::MergeFront() {
  const uint16_t idx1 = queue_.front().idx1;
  const uint16_t idx2 = queue_.front().idx2;
  assert(idx1 < idx2);
  // Merge idx2 into idx1.
  HistogramAdd(*(*histograms_)[idx2], *(*histograms_)[idx1],
               (*histograms_)[idx1]);
  (*histograms_)[idx1]->bit_cost_ = queue_.front().cost_combo;

  for (size_t i = 0; i < queue_.size();) {
    HistogramPair* const p = &queue_[i];
    if (p->idx1 == idx1 || p->idx2 == idx1 || p->idx1 == idx2 ||
        p->idx2 == idx2) {
      // Pop the pair from the list.
      *p = queue_.back();
      queue_.pop_back();
      continue;
    }
    MaybeUpdateHead(p);
    ++i;
  }
  // Remove idx2 from valid indices.
  const auto iter = std::lower_bound(indices_.begin(), indices_.end(), idx2);
  assert(iter != indices_.end() && *iter == idx2);
  std::copy(iter + 1, indices_.end(), iter);
  indices_.pop_back();
}

WP2Status HistoQueue::Push(uint16_t idx1, uint16_t idx2, double threshold,
                           double* const cost_diff) {
  if (idx1 > idx2) std::swap(idx1, idx2);
  HistogramPair pair;
  pair.idx1 = idx1;
  pair.idx2 = idx2;
  const Histogram* const h1 = (*histograms_)[idx1];
  const Histogram* const h2 = (*histograms_)[idx2];
  const double sum_cost = h1->bit_cost_ + h2->bit_cost_;
  pair.cost_combo = 0.;
  analyzer_->GetCombinedHistogramEntropy(*h1, *h2, sum_cost + threshold,
                                         &pair.cost_combo);
  pair.cost_diff = pair.cost_combo - sum_cost;

  // Do not even consider the pair if it does not improve the entropy.
  if (cost_diff != nullptr) *cost_diff = 0;
  if (pair.cost_diff >= threshold) return WP2_STATUS_OK;

  // We cannot add more elements than the capacity.
  assert(queue_.size() < queue_.capacity());
  WP2_CHECK_ALLOC_OK(queue_.push_back(pair, /*resize_if_needed=*/false));
  MaybeUpdateHead(&queue_.back());

  if (cost_diff != nullptr) *cost_diff = pair.cost_diff;
  return WP2_STATUS_OK;
}

WP2Status HistoQueue::RemoveInvalidHistograms() {
  uint32_t size = 0;
  for (uint32_t i : indices_) {
    std::swap((*histograms_)[i], (*histograms_)[size]);
    ++size;
  }
  WP2_CHECK_ALLOC_OK(histograms_->resize(size));
  queue_.clear();
  return WP2_STATUS_OK;
}

void HistoQueue::MaybeUpdateHead(HistogramPair* const pair) {
  assert(pair >= queue_.data() && pair <= &queue_.back());
  assert(!queue_.empty());
  if (pair->cost_diff < queue_.front().cost_diff) {
    // Replace the best pair.
    const HistogramPair tmp = queue_.front();
    queue_.front() = *pair;
    *pair = tmp;
  }
}

// -----------------------------------------------------------------------------

// Combines histograms by continuously choosing the one with the highest cost
// reduction.
// Returns true on success, false otherwise.
static WP2Status HistogramCombineGreedy(HistogramSet* const image_histo,
                                        PopulationAnalyzer* const analyzer) {
  const uint32_t image_histo_size = image_histo->histograms_.size();

  // Priority queue of histogram pairs.
  HistoQueue queue;
  WP2_CHECK_STATUS(queue.Init(&image_histo->histograms_, analyzer));
  WP2_CHECK_STATUS(queue.Resize((image_histo_size - 1) * image_histo_size / 2));

  for (uint32_t i = 0; i < image_histo_size; ++i) {
    for (uint32_t j = i + 1; j < image_histo_size; ++j) {
      WP2_CHECK_STATUS(queue.Push(i, j, 0.));
    }
  }

  while (queue.ValidIndices().size() > 1 && !queue.empty()) {
    const uint32_t idx1 = queue.front().idx1;
    queue.MergeFront();

     // Push new pairs formed with combined histogram to the queue.
    for (uint32_t i : queue.ValidIndices()) {
      if (i != idx1) {
        WP2_CHECK_STATUS(queue.Push(idx1, i, /*threshold=*/0.));
      }
    }
  }
  // Remove histograms that have been merged.
  WP2_CHECK_STATUS(queue.RemoveInvalidHistograms());

  return WP2_STATUS_OK;
}

// Perform histogram aggregation using a stochastic approach.
// 'do_greedy' is set to 1 if a greedy approach needs to be performed
// afterwards, 0 otherwise.
// Returns true on success, false otherwise.
static WP2Status HistogramCombineStochastic(HistogramSet* const image_histo,
                                            uint32_t min_cluster_size,
                                            float min_diff,
                                            PopulationAnalyzer* const analyzer,
                                            bool* const do_greedy) {
  WP2::UniformIntDistribution random;
  const uint32_t outer_iters = image_histo->histograms_.size();
  const uint32_t num_tries_no_success = outer_iters / 2;
  // Priority queue of histogram pairs. Its size impacts the quality of the
  // compression and the speed: the smaller the faster but the worse for the
  // compression.
  HistoQueue queue;
  constexpr uint32_t kHistoQueueSize = 9;

  WP2_CHECK_STATUS(queue.Init(&image_histo->histograms_, analyzer));
  WP2_CHECK_STATUS(queue.Resize(kHistoQueueSize));

  // Collapse similar histograms in 'image_histo'.
  ++min_cluster_size;
  uint32_t valid_size = queue.ValidIndices().size();
  for (uint32_t iter = 0, tries_with_no_success = 0;
       iter < outer_iters && valid_size >= min_cluster_size &&
       ++tries_with_no_success < num_tries_no_success;
       ++iter) {
    double best_cost = queue.empty() ? 0. : queue.front().cost_diff;
    // valid_size / 2 was chosen empirically. Less means faster but worse
    // compression.
    const uint32_t num_tries = valid_size / 2;
    for (uint32_t j = 0; j < num_tries; ++j) {
      double curr_cost;
      // Choose two different histograms at random and try to combine them.
      uint32_t idx1 = random.Get(0u, valid_size - 1);
      uint32_t idx2 = random.Get(0u, valid_size - 2);
      if (idx2 >= idx1) ++idx2;
      idx1 = queue.ValidIndices()[idx1];
      idx2 = queue.ValidIndices()[idx2];
      if (idx1 > idx2) std::swap(idx1, idx2);

      // Calculate cost reduction on combination.
      WP2_CHECK_STATUS(queue.Push(idx1, idx2, best_cost, &curr_cost));
      if (curr_cost < min_diff) {  // found a better pair?
        best_cost = curr_cost;
        // Empty the queue if we reached full capacity.
        if (queue.full()) break;
      }
    }
    if (queue.empty()) continue;

    // Merge the two best histograms.
    queue.MergeFront();
    valid_size = queue.ValidIndices().size();

    tries_with_no_success = 0;
  }
  WP2_CHECK_STATUS(queue.RemoveInvalidHistograms());
  *do_greedy = (valid_size <= min_cluster_size);
  return WP2_STATUS_OK;
}

// -----------------------------------------------------------------------------
// Histogram refinement

// Find the best 'out' histogram for each of the 'in' histograms.
// Note: we assume that out[]->bit_cost_ is already up-to-date.
static WP2Status HistogramRemap(const HistogramSet& in, HistogramSet* const out,
                                PopulationAnalyzer* const analyzer,
                                uint16_t* const mapping) {
  Histogram* const* const in_histo = in.histograms_.data();
  Histogram* const* const out_histo = out->histograms_.data();
  const uint32_t in_size = in.histograms_.size();
  const uint32_t out_size = out->histograms_.size();
  if (out_size > 1) {
    for (size_t i = 0; i < in_size; ++i) {
      int best_out = 0;
      double best_bits = std::numeric_limits<double>::max();
      for (uint32_t k = 0; k < out_size; ++k) {
        const double cur_bits = analyzer->HistogramAddThresh(
            *out_histo[k], *in_histo[i], best_bits);
        if (k == 0 || cur_bits < best_bits) {
          best_bits = cur_bits;
          best_out = k;
        }
      }
      mapping[i] = best_out;
    }
  } else {
    assert(out_size == 1);
    for (size_t i = 0; i < in_size; ++i) {
      mapping[i] = 0;
    }
  }
  // Re-write histogram_symbols so that it is a proper label image with values
  // increasing.
  // TODO(vrabaud) Which index should be given to empty histograms ? The most
  //               common one or the same as the previous block to favor LZ77 ?
  const uint32_t num_histos = in.histograms_.size();
  WP2::Vector_u32 mapping_tmp;
  WP2_CHECK_ALLOC_OK(mapping_tmp.resize(num_histos));
  std::fill(mapping_tmp.begin(), mapping_tmp.end(), num_histos);
  uint32_t mapping_max = 0;
  for (uint32_t i = 0; i < num_histos; ++i) {
    if (mapping_tmp[mapping[i]] == num_histos) {
      mapping_tmp[mapping[i]] = mapping_max++;
    }
    mapping[i] = mapping_tmp[mapping[i]];
  }

  // Recompute each out based on raw and symbols.
  for (size_t i = 0; i < out_size; ++i) {
    out_histo[i]->Clear();
  }

  for (size_t i = 0; i < in_size; ++i) {
    const int idx = mapping[i];
    HistogramAdd(*in_histo[i], *out_histo[idx], out_histo[idx]);
  }
  for (size_t i = 0; i < out_size; ++i) {
    analyzer->UpdateHistogramCost(out_histo[i]);
  }
  return WP2_STATUS_OK;
}

// Number of histograms for which we try to merge them all into one to see if
// it is more efficient.
// TODO(vrabaud) Make this parameter speed dependent.
const int kAllHistogramMergingNum = 10;

WP2Status GetHistoImageSymbols(uint32_t width, uint32_t height,
                               const BackwardRefs& refs, int speed,
                               int histogram_bits,
                               const LosslessSymbolsInfo& symbols_info,
                               HistogramSet* const image_histo,
                               uint16_t* const histogram_symbols) {
  const uint32_t histo_xsize =
      histogram_bits ? WP2LSubSampleSize(width, histogram_bits) : 1;
  const uint32_t histo_ysize =
      histogram_bits ? WP2LSubSampleSize(height, histogram_bits) : 1;
  const uint32_t image_histo_raw_size = histo_xsize * histo_ysize;
  HistogramSet orig_histo;
  WP2_CHECK_ALLOC_OK(orig_histo.Allocate(image_histo_raw_size, symbols_info));

  // Construct the histograms from backward references.
  WP2_CHECK_STATUS(BuildHistograms(width, height, histogram_bits, symbols_info,
                                   refs, &orig_histo));
  // Copies the histograms and computes its bit_cost.
  PopulationAnalyzer analyzer(image_histo_raw_size);
  WP2_CHECK_STATUS(analyzer.Allocate(symbols_info));
  for (Histogram* const h : orig_histo.histograms_) {
    analyzer.UpdateHistogramCost(h);
  }
  WP2_CHECK_STATUS(HistogramSetCopy(orig_histo, symbols_info, image_histo));

  WP2::Vector_u16 bin_map;
  BinMerger merger(orig_histo.histograms_.size(), speed);
  WP2_CHECK_STATUS(merger.RemoveEmptyHistograms(image_histo));
  WP2_CHECK_ALLOC_OK(bin_map.resize(image_histo->histograms_.size()));
  const bool entropy_combine =
      merger.AnalyzeEntropyBin(*image_histo, bin_map.data());
  if (entropy_combine) {
    // Collapse histograms with similar entropy.
    WP2_CHECK_STATUS(merger.CombineEntropyBin(bin_map, &analyzer, image_histo));
  }

  // Don't combine the histograms using stochastic and greedy heuristics for
  // low-effort compression mode.
  if (speed > 0 || !entropy_combine) {
    const float x = (speed == 5) ? 75.f / 100.f : speed / 9.f;
    // cubic ramp between 1 and MAX_HISTO_GREEDY:
    const auto threshold_size = (int)(1 + (x * x * x) * (MAX_HISTO_GREEDY - 1));
    bool do_greedy;
    WP2_CHECK_STATUS(HistogramCombineStochastic(image_histo, threshold_size,
                                                0.f, &analyzer, &do_greedy));
    if (do_greedy) {
      WP2_CHECK_STATUS(HistogramCombineGreedy(image_histo, &analyzer));
    }
  }

  // Sometimes, combining histograms provides worse entropy. But if we
  // combine all the histograms, we could have the overall extra entropy be
  // less than the cost of adding the entropy image.
  if (image_histo->histograms_.size() > 1 &&
      image_histo->histograms_.size() < kAllHistogramMergingNum &&
      image_histo->histograms_.size() < (orig_histo.histograms_.size() - 1)) {
    WP2_CHECK_STATUS(
        HistogramRemap(orig_histo, image_histo, &analyzer, histogram_symbols));

    // Cost of all the histograms and the entropy image.
    size_t cost_all;
    // Figure out the entropy of the symbols image.
    {
      WP2::ANSEnc enc;
      HashChain hash_chain;
      WP2_CHECK_ALLOC_OK(hash_chain.Allocate(image_histo_raw_size));
      BackwardRefsPool ref_pool;
      ref_pool.Init(image_histo_raw_size);

      WP2::Vector_s16 labels;
      WP2_CHECK_ALLOC_OK(labels.resize(4 * image_histo_raw_size));

      for (size_t i = 0; i < image_histo_raw_size; ++i) {
        labels[4 * i + 0] = 0;
        labels[4 * i + 1] = 0;
        labels[4 * i + 2] = histogram_symbols[i];
        labels[4 * i + 3] = 0;
      }
      WP2::ANSDictionaries dicts;
      LosslessSymbolsInfo symbols_info_tmp;
      // 'has_alpha' is false as it is a label image.
      symbols_info_tmp.Init(/*has_alpha=*/false, symbols_info.SampleFormat());
      symbols_info_tmp.SetAsLabelImage();
      WP2_CHECK_STATUS(EncodeImageNoClusters(
          &enc, &dicts, labels.data(), &hash_chain, &ref_pool, histo_ysize,
          histo_xsize, symbols_info_tmp, speed));
      cost_all = enc.GetCost(dicts);
    }

    // Compute the cost if we merge all histograms and update cost_all with the
    // cost of each histogram.
    Histogram* const tmp_histo = image_histo->tmp_histo_;
    image_histo->histograms_[0]->CopyTo(tmp_histo);
    cost_all += image_histo->histograms_[0]->bit_cost_;
    for (uint32_t i = 1; i < image_histo->histograms_.size(); ++i) {
      cost_all += image_histo->histograms_[i]->bit_cost_;
      HistogramAdd(*image_histo->histograms_[i], *tmp_histo, tmp_histo);
    }
    analyzer.UpdateHistogramCost(tmp_histo);

    if (tmp_histo->bit_cost_ < cost_all) {
      // Use that merged histogram.
      WP2_CHECK_ALLOC_OK(image_histo->histograms_.resize(1));
      tmp_histo->CopyTo(image_histo->histograms_[0]);
    }
  }

  // In case we have too many histograms, merge the least bad pairs.
  while (image_histo->histograms_.size() > kMaxHistogramImageSize) {
    bool do_greedy;
    WP2_CHECK_STATUS(HistogramCombineStochastic(
        image_histo, kMaxHistogramImageSize, std::numeric_limits<float>::max(),
        &analyzer, &do_greedy));
  }

  WP2_CHECK_STATUS(HistogramRemap(orig_histo, image_histo, &analyzer,
                                  histogram_symbols));

  return WP2_STATUS_OK;
}

}  // namespace WP2L
