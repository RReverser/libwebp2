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
// Improves a given set of backward references by analyzing its bit cost.
// The algorithm is similar to the Zopfli compression algorithm but tailored to
// images.
//
// Author: Vincent Rabaud (vrabaud@google.com)
//

#include <cassert>

#include "src/common/lossless/color_cache.h"
#include "src/dsp/lossless/lossless_common.h"
#include "src/enc/lossless/backward_references_enc.h"
#include "src/enc/lossless/histogram_enc.h"
#include "src/enc/lossless/losslessi_enc.h"
#include "src/enc/symbols_enc.h"

namespace WP2L {

static void ConvertPopulationCountTableToBitEstimates(
    int num_symbols, const uint32_t population_counts[], double output[]) {
  uint32_t sum = 0;
  int nonzeros = 0;
  int i;
  for (i = 0; i < num_symbols; ++i) {
    sum += population_counts[i];
    if (population_counts[i] > 0) {
      ++nonzeros;
    }
  }
  if (nonzeros <= 1) {
    memset(output, 0, num_symbols * sizeof(*output));
  } else {
    const double logsum = WP2Log2(sum);
    for (i = 0; i < num_symbols; ++i) {
      output[i] = logsum - WP2Log2(population_counts[i]);
    }
  }
}

// Class storing the cost of each symbol in bits ,based on its count.
class CostModel : public CountsBuffer<double> {
 public:
  // 'refs' must not have its distance modified by DistanceToPlaneCode.
  WP2Status Build(uint32_t width, const LosslessSymbolsInfo& info,
                  const BackwardRefs& refs) {
    width_ = width;
    HistogramSet histo;
    WP2_CHECK_ALLOC_OK(histo.Allocate(/*size=*/1, info));

    WP2_CHECK_STATUS(BuildHistograms(width, /*height=*/0, /*histo_bits=*/0,
                                     info, refs, &histo));
    has_alpha_ = LosslessSymbolsInfo::HasAlpha(info);
    is_label_image_ = LosslessSymbolsInfo::IsLabelImage(info);

    for (size_t i = 0; i < kSymbolNum; ++i) {
      ConvertPopulationCountTableToBitEstimates(
          histo.symbols_info_->GetMaxRange(i), histo.histograms_[0]->counts_[i],
          counts_[i]);
    }
    return WP2_STATUS_OK;
  }

  inline double GetLiteralCost(const int16_t* const v) const {
    return counts_[kSymbolType][kSymbolTypeLiteral] +
           (has_alpha_ ? counts_[kSymbolA][MakeIndex(v[0])] : 0.) +
           (is_label_image_ ? 0.
                            : counts_[kSymbolR][MakeIndex(v[1])] +
                                  counts_[kSymbolB][MakeIndex(v[3])]) +
           counts_[kSymbolG][MakeIndex(v[2])];
  }

  inline double GetCacheCost(uint32_t idx) const {
    return counts_[kSymbolType][kSymbolTypeCacheIdx] +
           counts_[kSymbolCache][idx];
  }

  inline double GetLengthCost(uint32_t length) const {
    int code, extra_bits;
    WP2LPrefixEncodeBits(length, &code, &extra_bits);
    return counts_[kSymbolType][kSymbolTypeCopy] +
           counts_[kSymbolLen][code] + extra_bits;
  }

  inline double GetDistanceCost(uint32_t distance) const {
    distance = DistanceToPlaneCode(width_, distance);
    int code, extra_bits;
    WP2LPrefixEncodeBits(distance, &code, &extra_bits);
    return counts_[kSymbolDist][code] + extra_bits;
  }

 private:
  uint32_t width_;
  bool has_alpha_;
  bool is_label_image_;
};

// Describes how a segment of pixels of length 'len' is represented (copy ...)
struct Segment {
  uint16_t len;
  SymbolType mode;
};

static inline void AddSingleLiteralWithCostModel(
    const int16_t* const argb, ColorCache* const color_cache,
    const CostModel& cost_model, int idx, int use_color_cache, float prev_cost,
    float* const cost, Segment* const segments) {
  double cost_val = prev_cost;
  const int16_t* const color = &argb[4 * idx];
  uint32_t tag;
  if (use_color_cache && color_cache->Contains(color, &tag)) {
    // use_color_cache is true and color_cache contains color
    const double mul0 = 0.68;
    cost_val += cost_model.GetCacheCost(tag) * mul0;
    if (cost[idx] > cost_val) {
      cost[idx] = (float)cost_val;
      // Only one is inserted.
      segments[idx].len = 1;
      segments[idx].mode = kSymbolTypeCacheIdx;
    }
  } else {
    const double mul1 = 0.82;
    if (use_color_cache) color_cache->Insert(color, /*index_ptr=*/nullptr);
    cost_val += cost_model.GetLiteralCost(color) * mul1;
    if (cost[idx] > cost_val) {
      cost[idx] = (float)cost_val;
      // Only one is inserted.
      segments[idx].len = 1;
      segments[idx].mode = kSymbolTypeLiteral;
    }
  }
}

// -----------------------------------------------------------------------------
// CostManager and interval handling

// Empirical value to avoid high memory consumption but good for performance.
#define COST_CACHE_INTERVAL_SIZE_MAX 500

// To perform backward reference every pixel at 'index' is considered and
// the cost for the kMaxLZ77Length following pixels computed. Those following
// pixels at index 'index' + k (k from 0 to kMaxLZ77Length) have a cost of:
//     cost_ = distance cost at index + GetLengthCost(cost_model, k)
// and the minimum value is kept. GetLengthCost(cost_model, k) is cached in an
// array of size kMaxLZ77Length.
// Instead of performing kMaxLZ77Length comparisons per pixel, we keep track of
// the minimal values using intervals of constant cost.
// An interval is defined by the 'index' of the pixel that generated it and is
// only useful in a range of indices from 'start' to 'end' (exclusive), i.e. it
// contains the minimum value for pixels between 'start' and 'end'. Intervals
// are stored in a linked list and ordered by 'start'. When a new interval has a
// better value, old intervals are split or removed. There are therefore no
// overlapping intervals.
struct CostInterval : WP2Allocable {
  float cost;
  int start;
  int end;
  int index;
  CostInterval* previous;
  CostInterval* next;
};

// The GetLengthCost(cost_model, k) are cached in a CostCacheInterval.
typedef struct {
  double cost;
  int start;
  int end;  // Exclusive.
} CostCacheInterval;

// This structure is in charge of managing intervals and costs.
// It caches the different CostCacheInterval, caches the different
// GetLengthCost(cost_model, k) in cost_cache_ and the CostInterval's (whose
// count_ is limited by COST_CACHE_INTERVAL_SIZE_MAX).
#define COST_MANAGER_MAX_FREE_LIST 10
class CostManager {
 public:
  CostManager() = default;
  ~CostManager() { Reset(); }
  bool Init(Segment* const segments, uint32_t num_pixels,
            const CostModel& cost_model);
  void Reset();

  CostInterval* head_ = nullptr;
  int count_ = 0;  // The number of stored intervals.
  WP2::VectorNoCtor<CostCacheInterval> cache_intervals_;
  size_t cache_intervals_size_ = 0;
  // Contains the GetLengthCost(cost_model, k). Size max: kMaxLZ77Length.
  WP2::VectorNoCtor<double> cost_cache_;
  WP2::Vector_f costs_;
  Segment* segments_ = nullptr;
  // Most of the time, we only need few intervals -> use a free-list, to avoid
  // fragmentation with small allocs in most common cases.
  CostInterval intervals_[COST_MANAGER_MAX_FREE_LIST];
  CostInterval* free_intervals_ = nullptr;
  // These are regularly malloc'd remains. This list can't grow larger than than
  // size COST_CACHE_INTERVAL_SIZE_MAX - COST_MANAGER_MAX_FREE_LIST, note.
  CostInterval* recycled_intervals_ = nullptr;
};

static void CostIntervalAddToFreeList(CostManager* const manager,
                                      CostInterval* const interval) {
  interval->next = manager->free_intervals_;
  manager->free_intervals_ = interval;
}

static int CostIntervalIsInFreeList(const CostManager* const manager,
                                    const CostInterval* const interval) {
  return (interval >= &manager->intervals_[0] &&
          interval <= &manager->intervals_[COST_MANAGER_MAX_FREE_LIST - 1]);
}

static void CostManagerInitFreeList(CostManager* const manager) {
  int i;
  manager->free_intervals_ = NULL;
  for (i = 0; i < COST_MANAGER_MAX_FREE_LIST; ++i) {
    CostIntervalAddToFreeList(manager, &manager->intervals_[i]);
  }
}

static void DeleteIntervalList(CostManager* const manager,
                               const CostInterval* interval) {
  while (interval != NULL) {
    const CostInterval* const next = interval->next;
    if (!CostIntervalIsInFreeList(manager, interval)) {
      WP2Free((void*)interval);
    }  // else: do nothing
    interval = next;
  }
}

void CostManager::Reset() {
  costs_.reset();
  cache_intervals_.reset();

  // Clear the interval lists.
  DeleteIntervalList(this, head_);
  head_ = nullptr;
  DeleteIntervalList(this, recycled_intervals_);
  recycled_intervals_ = nullptr;

  // Reset pointers, count_ and cache_intervals_size_.
  count_ = 0;
  cache_intervals_size_ = 0;
  cost_cache_.reset();
  segments_ = nullptr;
  memset(intervals_, 0, sizeof(intervals_));
  free_intervals_ = nullptr;
  CostManagerInitFreeList(this);
}

bool CostManager::Init(Segment* const segments, uint32_t num_pixels,
                       const CostModel& cost_model) {
  Reset();
  segments_ = segments;

  const uint32_t cost_cache_size = std::min(kMaxLZ77Length, num_pixels);
  if (!cost_cache_.resize(cost_cache_size)) return false;

  // Fill in the cost_cache_.
  cache_intervals_size_ = 1;
  cost_cache_[0] = cost_model.GetLengthCost(0);
  for (uint32_t i = 1; i < cost_cache_size; ++i) {
    cost_cache_[i] = cost_model.GetLengthCost(i);
    // Get the number of bound intervals.
    if (cost_cache_[i] != cost_cache_[i - 1]) {
      ++cache_intervals_size_;
    }
  }

  // With the current cost model, we usually have below 20 intervals.
  // The worst case scenario with a cost model would be if every length has a
  // different cost, hence kMaxLZ77Length but that is impossible with the
  // current implementation that spirals around a pixel.
  assert(cache_intervals_size_ <= kMaxLZ77Length);
  if (!cache_intervals_.resize(cache_intervals_size_)) {
    Reset();
    return false;
  }

  // Fill in the cache_intervals_.
  {
    CostCacheInterval* cur = cache_intervals_.data();

    // Consecutive values in cost_cache_ are compared and if a big enough
    // difference is found, a new interval is created and bounded.
    cur->start = 0;
    cur->end = 1;
    cur->cost = cost_cache_[0];
    for (uint32_t i = 1; i < cost_cache_size; ++i) {
      const double cost_val = cost_cache_[i];
      if (cost_val != cur->cost) {
        ++cur;
        // Initialize an interval.
        cur->start = i;
        cur->cost = cost_val;
      }
      cur->end = i + 1;
    }
  }

  if (!costs_.resize(num_pixels)) {
    Reset();
    return false;
  }
  // Set the initial costs_ high for every pixel as we will keep the minimum.
  std::fill(costs_.begin(), costs_.end(), 1e38f);

  return true;
}

// Given the cost and the position that define an interval, update the cost at
// pixel 'i' if it is smaller than the previously computed value.
static inline void UpdateCost(CostManager* const manager, int i, int position,
                              float cost) {
  const uint32_t k = i - position;
  assert(i >= position && k < kMaxLZ77Length);

  if (manager->costs_[i] > cost) {
    manager->costs_[i] = cost;
    // Force the mode to be a copy.
    manager->segments_[i].len = (uint16_t)(k + 1);
    manager->segments_[i].mode = kSymbolTypeCopy;
  }
}

// Given the cost and the position that define an interval, update the cost for
// all the pixels between 'start' and 'end' excluded.
static inline void UpdateCostPerInterval(CostManager* const manager,
                                              int start, int end, int position,
                                              float cost) {
  int i;
  for (i = start; i < end; ++i) UpdateCost(manager, i, position, cost);
}

// Given two intervals, make 'prev' be the previous one of 'next' in 'manager'.
static inline void ConnectIntervals(CostManager* const manager,
                                         CostInterval* const prev,
                                         CostInterval* const next) {
  if (prev != NULL) {
    prev->next = next;
  } else {
    manager->head_ = next;
  }

  if (next != NULL) next->previous = prev;
}

// Pop an interval in the manager.
static inline void PopInterval(CostManager* const manager,
                                    CostInterval* const interval) {
  if (interval == NULL) return;

  ConnectIntervals(manager, interval->previous, interval->next);
  if (CostIntervalIsInFreeList(manager, interval)) {
    CostIntervalAddToFreeList(manager, interval);
  } else {  // recycle regularly malloc'd intervals too
    interval->next = manager->recycled_intervals_;
    manager->recycled_intervals_ = interval;
  }
  --manager->count_;
  assert(manager->count_ >= 0);
}

// Update the cost at index i by going over all the stored intervals that
// overlap with i.
// If 'do_clean_intervals' is set to something different than 0, intervals that
// end before 'i' will be popped.
static inline void UpdateCostAtIndex(CostManager* const manager, int i,
                                          int do_clean_intervals) {
  CostInterval* current = manager->head_;

  while (current != NULL && current->start <= i) {
    CostInterval* const next = current->next;
    if (current->end <= i) {
      if (do_clean_intervals) {
        // We have an outdated interval, remove it.
        PopInterval(manager, current);
      }
    } else {
      UpdateCost(manager, i, current->index, current->cost);
    }
    current = next;
  }
}

// Given a current orphan interval and its previous interval, before
// it was orphaned (which can be NULL), set it at the right place in the list
// of intervals using the start_ ordering and the previous interval as a hint.
static inline void PositionOrphanInterval(CostManager* const manager,
                                               CostInterval* const current,
                                               CostInterval* previous) {
  assert(current != NULL);

  if (previous == NULL) previous = manager->head_;
  while (previous != NULL && current->start < previous->start) {
    previous = previous->previous;
  }
  while (previous != NULL && previous->next != NULL &&
         previous->next->start < current->start) {
    previous = previous->next;
  }

  if (previous != NULL) {
    ConnectIntervals(manager, current, previous->next);
  } else {
    ConnectIntervals(manager, current, manager->head_);
  }
  ConnectIntervals(manager, previous, current);
}

// Insert an interval in the list contained in the manager by starting at
// interval_in as a hint. The intervals are sorted by start_ value.
static inline WP2Status InsertInterval(CostManager* const manager,
                                       CostInterval* const interval_in,
                                       float cost, int position, int start,
                                       int end) {
  CostInterval* interval_new;

  if (start >= end) return WP2_STATUS_OK;
  if (manager->count_ >= COST_CACHE_INTERVAL_SIZE_MAX) {
    // Serialize the interval if we cannot store it.
    UpdateCostPerInterval(manager, start, end, position, cost);
    return WP2_STATUS_OK;
  }
  if (manager->free_intervals_ != NULL) {
    interval_new = manager->free_intervals_;
    manager->free_intervals_ = interval_new->next;
  } else if (manager->recycled_intervals_ != NULL) {
    interval_new = manager->recycled_intervals_;
    manager->recycled_intervals_ = interval_new->next;
  } else {  // malloc for good
    interval_new = new (WP2Allocable::nothrow) CostInterval();
    WP2_CHECK_ALLOC_OK(interval_new != nullptr);
  }

  interval_new->cost = cost;
  interval_new->index = position;
  interval_new->start = start;
  interval_new->end = end;
  PositionOrphanInterval(manager, interval_new, interval_in);

  ++manager->count_;
  return WP2_STATUS_OK;
}

// Given a new cost interval defined by its start at position, its length value
// and distance_cost, add its contributions to the previous intervals and costs.
// If handling the interval or one of its subintervals becomes to heavy, its
// contribution is added to the costs right away.
static inline WP2Status PushInterval(CostManager* const manager,
                                     double distance_cost, int position,
                                     int len) {
  size_t i;
  CostInterval* interval = manager->head_;
  CostInterval* interval_next;
  const CostCacheInterval* const cost_cache_intervals =
      manager->cache_intervals_.data();
  // If the interval is small enough, no need to deal with the heavy
  // interval logic, just serialize it right away. This constant is empirical.
  const int kSkipDistance = 10;

  if (len < kSkipDistance) {
    int j;
    for (j = position; j < position + len; ++j) {
      const uint32_t k = j - position;
      float cost_tmp;
      assert(j >= position && k < kMaxLZ77Length);
      cost_tmp = (float)(distance_cost + manager->cost_cache_[k]);

      if (manager->costs_[j] > cost_tmp) {
        manager->costs_[j] = cost_tmp;
        // Force the mode to be a copy.
        manager->segments_[j].len = (uint16_t)(k + 1);
        manager->segments_[j].mode = kSymbolTypeCopy;
      }
    }
    return WP2_STATUS_OK;
  }

  for (i = 0; i < manager->cache_intervals_size_ &&
              cost_cache_intervals[i].start < len;
       ++i) {
    // Define the intersection of the ith interval with the new one.
    int start = position + cost_cache_intervals[i].start;
    const int end =
        position +
        (cost_cache_intervals[i].end > len ? len : cost_cache_intervals[i].end);
    const float cost = (float)(distance_cost + cost_cache_intervals[i].cost);

    for (; interval != NULL && interval->start < end;
         interval = interval_next) {
      interval_next = interval->next;

      // Make sure we have some overlap
      if (start >= interval->end) continue;

      if (cost >= interval->cost) {
        // When intervals are represented, the lower, the better.
        // [**********************************************************[
        // start                                                    end
        //                   [----------------------------------[
        //                   interval->start_       interval->end_
        // If we are worse than what we already have, add whatever we have so
        // far up to interval.
        const int start_new = interval->end;
        WP2_CHECK_STATUS(InsertInterval(manager, interval, cost, position,
                                        start, interval->start));
        start = start_new;
        if (start >= end) break;
        continue;
      }

      if (start <= interval->start) {
        if (interval->end <= end) {
          //                   [----------------------------------[
          //                   interval->start_       interval->end_
          // [**************************************************************[
          // start                                                        end
          // We can safely remove the old interval as it is fully included.
          PopInterval(manager, interval);
        } else {
          //              [------------------------------------[
          //              interval->start_        interval->end_
          // [*****************************[
          // start                       end
          interval->start = end;
          break;
        }
      } else {
        if (end < interval->end) {
          // [--------------------------------------------------------------[
          // interval->start_                                  interval->end_
          //                     [*****************************[
          //                     start                       end
          // We have to split the old interval as it fully contains the new one.
          const int end_original = interval->end;
          interval->end = start;
          WP2_CHECK_STATUS(InsertInterval(manager, interval, interval->cost,
                                          interval->index, end, end_original));
          interval = interval->next;
          break;
        } else {
          // [------------------------------------[
          // interval->start_        interval->end_
          //                     [*****************************[
          //                     start                       end
          interval->end = start;
        }
      }
    }
    // Insert the remaining interval from start to end.
    WP2_CHECK_STATUS(
        InsertInterval(manager, interval, cost, position, start, end));
  }
  return WP2_STATUS_OK;
}

// Optimizes the BackwardReferences based on the distance code.
static WP2Status BackwardReferencesHashChainDistanceOnly(
    uint32_t width, uint32_t height, const int16_t* const argb,
    const LosslessSymbolsInfo& symbols_info, const HashChain& hash_chain,
    const CacheConfig& cache_config, const BackwardRefs& refs,
    Segment* const segments) {
  const uint32_t num_pixels = width * height;

  CostModel cost_model;
  uint32_t offset_prev = -1, len_prev = -1;
  double offset_cost = -1;
  int first_offset_is_constant = -1;  // initialized with 'impossible' value
  uint32_t reach = 0;

  WP2_CHECK_ALLOC_OK(cost_model.Allocate(&symbols_info));
  const bool use_color_cache = (cache_config.type != CacheType::kNone);
  ColorCachePtr color_cache;
  WP2_CHECK_STATUS(color_cache.Init(cache_config));

  WP2_CHECK_STATUS(cost_model.Build(width, symbols_info, refs));

  CostManager cost_manager;
  WP2_CHECK_ALLOC_OK(cost_manager.Init(segments, num_pixels, cost_model));

  // We loop one pixel at a time, but store all currently best points to
  // non-processed locations from this point.
  segments[0].len = 1;
  segments[0].mode = kSymbolTypeLiteral;
  // Add first pixel as literal.
  AddSingleLiteralWithCostModel(argb, color_cache.get(), cost_model, 0,
                                use_color_cache, 0.f,
                                cost_manager.costs_.data(), segments);

  for (size_t i = 1; i < num_pixels; ++i) {
    const float prev_cost = cost_manager.costs_[i - 1];
    uint32_t offset, len;
    hash_chain.FindCopy(i, &offset, &len);

    // Try adding the pixel as a literal.
    AddSingleLiteralWithCostModel(argb, color_cache.get(), cost_model, i,
                                  use_color_cache, prev_cost,
                                  cost_manager.costs_.data(), segments);

    // If we are dealing with a non-literal.
    if (len >= 2) {
      if (offset != offset_prev) {
        offset_cost = cost_model.GetDistanceCost(offset);
        first_offset_is_constant = 1;
        WP2_CHECK_STATUS(
            PushInterval(&cost_manager, prev_cost + offset_cost, i, len));
      } else {
        assert(offset_cost >= 0);
        assert(len_prev >= 0);
        assert(first_offset_is_constant == 0 || first_offset_is_constant == 1);
        // Instead of considering all contributions from a pixel i by calling:
        //         PushInterval(cost_manager, prev_cost + offset_cost, i, len);
        // we optimize these contributions in case offset_cost stays the same
        // for consecutive pixels. This describes a set of pixels similar to a
        // previous set (e.g. constant color regions).
        if (first_offset_is_constant) {
          reach = i - 1 + len_prev - 1;
          first_offset_is_constant = 0;
        }

        if (i + len - 1 > reach) {
          // We can only be go further with the same offset if the previous
          // length was maxed, hence len_prev == len == kMaxLZ77Length.
          // TODO(vrabaud): bump i to the end right away (insert cache and
          //                update cost).
          // TODO(vrabaud): check if one of the points in between does not have
          //                a lower cost.
          // Already consider the pixel at "reach" to add intervals that are
          // better than whatever we add.
          uint32_t offset_j, len_j = 0;
          uint32_t j;
          assert(len == kMaxLZ77Length || len == num_pixels - i);
          // Figure out the last consecutive pixel within [i, reach + 1] with
          // the same offset.
          for (j = i; j <= reach; ++j) {
            hash_chain.FindCopy(j + 1, &offset_j, &len_j);
            if (offset_j != offset) {
              hash_chain.FindCopy(j, &offset_j, &len_j);
              break;
            }
          }
          // Update the cost at j - 1 and j.
          UpdateCostAtIndex(&cost_manager, j - 1, 0);
          UpdateCostAtIndex(&cost_manager, j, 0);

          WP2_CHECK_STATUS(
              PushInterval(&cost_manager,
                           cost_manager.costs_[j - 1] + offset_cost, j, len_j));
          reach = j + len_j - 1;
        }
      }
    }

    UpdateCostAtIndex(&cost_manager, i, 1);
    offset_prev = offset;
    len_prev = len;
  }

  return WP2_STATUS_OK;
}

// We pack the path at the end of *segments and return
// a pointer to this part of the array. Example:
// segments = [1x2xx3x2] => packed [1x2x1232], chosen_path = [1232]
static void TraceBackwards(Segment* const segments, size_t segments_size,
                           Segment** const chosen_path,
                           size_t* const chosen_path_size) {
  Segment* path = segments + segments_size;
  Segment* cur = segments + segments_size - 1;
  while (cur >= segments) {
    --path;
    *path = *cur;
    cur -= cur->len;
  }
  *chosen_path = path;
  *chosen_path_size = (size_t)(segments + segments_size - path);
}

static WP2Status BackwardReferencesHashChainFollowChosenPath(
    const int16_t* const argb, const CacheConfig& cache_config,
    const Segment* const chosen_path, uint32_t chosen_path_size,
    const HashChain& hash_chain, BackwardRefs* const refs,
    bool* const check_histogram) {
  const bool use_color_cache = (cache_config.type != CacheType::kNone);
  ColorCachePtr color_cache;
  WP2_CHECK_STATUS(color_cache.Init(cache_config));

  refs->Clear();
  for (uint32_t ix = 0, i = 0; ix < chosen_path_size; ++ix) {
    const Segment& seg = chosen_path[ix];
    if (seg.mode == kSymbolTypeCopy) {
      // TODO(vrabaud) do not go through here when seg.len_ == 1
      const uint32_t len = seg.len;
      const uint32_t offset = hash_chain.FindOffset(i);
      WP2_CHECK_STATUS(refs->CursorAdd(PixelMode::CreateCopy(offset, len)));
      if (use_color_cache) {
        for (uint32_t k = 0; k < len; ++k) {
          color_cache->Insert(&argb[4 * (i + k)], /*index_ptr=*/nullptr);
        }
      }
      i += len;
    } else if (seg.mode == kSymbolTypeCacheIdx) {
      uint32_t tag;
      if (!color_cache->Contains(&argb[4 * i], &tag)) {
        *check_histogram = false;
        break;
      }
      WP2_CHECK_STATUS(refs->CursorAdd(
          PixelMode::CreateCacheIdx(tag, color_cache->IndexRange())));
      color_cache->Insert(&argb[4 * i], /*index_ptr=*/nullptr);
      ++i;
    } else {
      assert(seg.mode == kSymbolTypeLiteral);
      if (use_color_cache) {
        color_cache->Insert(&argb[4 * i], /*index_ptr=*/nullptr);
      }
      WP2_CHECK_STATUS(refs->CursorAdd(PixelMode::CreateLiteral(&argb[4 * i])));
      ++i;
    }
  }
  return WP2_STATUS_OK;
}

// Returns 1 on success.
extern WP2Status BackwardReferencesTraceBackwards(
    uint32_t width, uint32_t height, const int16_t* const argb,
    const LosslessSymbolsInfo& symbols_info, const HashChain& hash_chain,
    const CacheConfig& cache_config, const BackwardRefs& refs_src,
    BackwardRefs* const refs_dst, bool* const check_histogram);
WP2Status BackwardReferencesTraceBackwards(
    uint32_t width, uint32_t height, const int16_t* const argb,
    const LosslessSymbolsInfo& symbols_info, const HashChain& hash_chain,
    const CacheConfig& cache_config, const BackwardRefs& refs_src,
    BackwardRefs* const refs_dst, bool* const check_histogram) {
  const uint32_t num_pixels = width * height;
  WP2::VectorNoCtor<Segment> segments;
  *check_histogram = true;

  WP2_CHECK_ALLOC_OK(segments.resize(num_pixels));

  WP2_CHECK_STATUS(BackwardReferencesHashChainDistanceOnly(
      width, height, argb, symbols_info, hash_chain, cache_config, refs_src,
      segments.data()));
  Segment* chosen_path = nullptr;
  size_t chosen_path_size = 0;
  TraceBackwards(segments.data(), num_pixels, &chosen_path, &chosen_path_size);
  WP2_CHECK_STATUS(BackwardReferencesHashChainFollowChosenPath(
      argb, cache_config, chosen_path, chosen_path_size, hash_chain, refs_dst,
      check_histogram));
  return WP2_STATUS_OK;
}

}  // namespace WP2L
