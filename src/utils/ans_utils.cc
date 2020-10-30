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
// Assymetric Numeral System helper functions.
//
// Author: Skal (pascal.massimino@gmail.com)

#include "src/utils/ans_utils.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>

namespace WP2 {

//------------------------------------------------------------------------------
// Header bit I/O helpers.

uint32_t HeaderDec::ReadBitsInternal(int nb_bits) {
  assert(nb_bits <= 24);
  uint32_t value = 0;
  if (Ok()) {
    if (bit_pos_ + nb_bits > max_pos_) {
      error_ = true;
    } else {
      size_t pos = bit_pos_ >> 3;
      const int shift = bit_pos_ & 7;
      bit_pos_ += nb_bits;
      for (int n = 0; n < shift + nb_bits; n += 8) {
        value |= (uint32_t)buf_[pos++] << n;
      }
      value >>= shift;
      value &= (1 << nb_bits) - 1;
    }
  }
  return value;
}

void HeaderEnc::PutBitsInternal(uint32_t value, int nb_bits) {
  assert(nb_bits <= 24);
  assert(value < (1u << nb_bits));

  const int shift = bit_pos_ & 7;
  size_t pos = bit_pos_ >> 3;
  bit_pos_ += nb_bits;
  if (bit_pos_ > max_pos_) {
    error_ = true;
    return;
  }
  value <<= shift;
  for (int n = 0; n < nb_bits + shift; n += 8) {
    buf_[pos++] |= value & 0xff;
    value >>= 8;
  }
}

//------------------------------------------------------------------------------
// ANS quantizations.

void ANSCountsQuantizeHuffman(uint32_t size,
                              const uint32_t* const counts,
                              uint32_t* const out,
                              uint32_t max_bits,
                              float* const cost) {
  const uint32_t current_max_value = *std::max_element(counts, counts + size);
  const uint32_t max_value = std::min(1u << (max_bits - 1), current_max_value);

  uint32_t sum = 0, sum_q = 0;
  const float norm = 1. * max_value / current_max_value;
  if (cost != nullptr) *cost = 0;
  for (size_t c = 0; c < size; ++c) {
    if (counts[c] > 0) {
      const float scaled = std::max(counts[c] * norm, 1.f);
      const uint32_t num_bits = (uint32_t)(std::log2(scaled) + .5);
      const uint32_t huffmanized = 1u << num_bits;
      out[c] = huffmanized;
      if (cost != nullptr) {
        *cost -= counts[c] * num_bits;
        sum += counts[c];
        sum_q += huffmanized;
      }
    } else {
      out[c] = 0;
    }
  }

  if (cost != nullptr) {
    *cost += sum * log2(sum_q);
  }
}

bool ANSCountsQuantize(bool do_expand, uint32_t max_freq, uint32_t size,
                       uint32_t* const counts, float* const cost) {
  return ANSCountsQuantize(do_expand, max_freq, size, counts, counts, cost);
}

bool ANSCountsQuantize(bool do_expand, uint32_t max_freq, uint32_t size,
                       const uint32_t* const counts, uint32_t* const out,
                       float* const cost) {
  assert(max_freq > 0);
  if (size == 0) return false;
  const uint32_t max_value = *std::max_element(counts, counts + size);
  if (max_value == 0 || max_freq == max_value ||
      (!do_expand && max_value <= max_freq)) {
    if (out != counts) memcpy(out, counts, size * sizeof(uint32_t));
    if (cost != nullptr) *cost = ANSCountsQuantizedCost(counts, counts, size);
    return false;
  }
  uint32_t sum = 0, sum_q = 0;
  if (cost != nullptr) *cost = 0;
  const float norm = 1. * max_freq / max_value;
  for (size_t c = 0; c < size; ++c) {
    if (counts[c] > 0) {
      const int new_count = (int)(counts[c] * norm + .5);
      out[c] = (new_count < 1) ? 1 : new_count;
      if (cost != nullptr) {
        *cost -= counts[c] * log2(out[c]);
        sum += counts[c];
        sum_q += out[c];
      }
    } else {
      out[c] = 0;
    }
  }
  if (cost != nullptr) *cost += sum * log2(sum_q);
  return true;
}

// The bit cost of an element "i" is:
// log2(1/proba[i]) = log2(sum(distribution)/distribution[i])
float ANSCountsQuantizedCost(const uint32_t* counts,
                             const uint32_t* distribution, uint32_t size) {
  uint32_t sum = 0, sum_q = 0;
  float cost = 0.f;
  for (size_t i = 0; i < size; ++i) {
    if (distribution[i] == 0) continue;
    cost -= counts[i] * log2(distribution[i]);
    sum += counts[i];
    sum_q += distribution[i];
  }
  if (sum_q > 0) cost += sum * log2(sum_q);
  // Remove rounding errors when only one count is non-zero.
  if (cost < 0.f) cost = 0.f;
  return cost;
}

//------------------------------------------------------------------------------
// ANS vector storage.

// Compute the bit cost when merging two stats.
static inline uint32_t
ComputeCost(const OptimizeArrayStorageStat* const s1,
            const OptimizeArrayStorageStat* const s2) {
  const int bit_diff = (s1->n_bits - s2->n_bits);
  int total_diff = bit_diff * ((bit_diff > 0) ? s2->count : -s1->count);
  // If there was only one value, by merging we lose the optimization that we
  // don't need to store the high order bit.
  // TODO(maryla): This is approximate. In the case where n_bits == n_bits_max
  // (when values are stored with an RValue instead of a UValue), we potentially
  // win a lot more than one bit.
  if (s1->count == 1) ++total_diff;
  if (s2->count == 1) ++total_diff;
  return total_diff;
}

// 'overhead' is the size in bits of: size in bits to store the number of
// elements and their size in bits.
// The list of stats is simplified to represent the same number of elements but
// with a smaller overall bit cost.
void OptimizeArrayStorage(OptimizeArrayStorageStat* const stats,
                          size_t* const size_in, float overhead) {
// TODO(vrabaud) include a speed parameter to choose between the two methods.
// or at least, choose the slower one if size is small (e.g. < 5).
// The first method gives a speedup of an order of magnitude for a compression
// hit of less than 0.1%
  size_t size = *size_in;
#if OPTIMIZE_ARRAY_STORAGE
  while (size > 1) {
    auto s2 = stats + 1, s_end = stats + size;
    // Find the first possible merge, if any.
    while (s2 < s_end && ComputeCost(s2 - 1, s2) >= overhead) ++s2;

    auto s1 = s2 - 1;
    for (; s2 < s_end; ++s2) {
      if (ComputeCost(s1, s2) < overhead) {
        // Merge.
        s1->count = s2->count + s1->count;
        s1->n_bits = std::max(s2->n_bits, s1->n_bits);
      } else {
        ++s1;
        *s1 = *s2;
      }
    }
    const size_t size_new = s1 - stats + 1;
    if (size_new == size) break;
    size = size_new;
  }
#else
  // TODO(vrabaud) overhead actually decreases as the number of elements becomes
  // less and less: take that into account.
  // Fill the costs once and for all.
  for (size_t i = 1; i < size; ++i) {
    stats[i].cost_ = ComputeCost(&stats[i - 1], &stats[i]);
  }

  size_t start = 1, end = 1;
  while (start < size) {
    // Find the first element that could be optimized.
    while (start < size && stats[start].cost_ >= overhead) {
      ++start;
    }
    if (start == size) break;

    // Find the first following element that cannot be optimized (or the end).
    if (end <= start) {
      end = start + 1;
      while (end < size && stats[end].cost_ < overhead) {
        ++end;
      }
    }

    do {
      size_t i = start, i_max = i;
      // Find the minimal cost.
      uint32_t max_cost = stats[i].cost_;
      for (++i; i < end; ++i) {
        if (stats[i].cost_ < max_cost) {
          max_cost = stats[i].cost_;
          i_max = i;
        }
      }
      if (max_cost < overhead) break;
      // Only continue if we were able to have a good enough cost.
      auto max_s2 = stats + i_max, max_s1 = max_s2 - 1;
      max_s1->count_ = max_s1->count_ + max_s2->count_;
      max_s1->n_bits_ = std::max(max_s1->n_bits_, max_s2->n_bits_);
      if (max_s1 != stats) max_s1->cost_ = ComputeCost(max_s1 - 1, max_s1);
      // Remove max_s2.
      memmove(stats + i_max, stats + i_max + 1,
              (size - i_max - 1) * sizeof(OptimizeArrayStorageStat));
      --size;
      --end;
      if (i_max != size) {
        stats[i_max].cost_ = ComputeCost(stats + i_max - 1, stats + i_max);
      }
    } while (start < end);
  }
#endif
  *size_in = size;
}


////////////////////////////////////////////////////////////////////////////////
// ANSEncCounter

float ANSEncCounter::GetCost() const { return cost_ + symbol_cost_; }

float ANSEncCounter::GetCost(const ANSDictionaries&dicts) const {
  float cost = cost_;
  for (const auto& d : dicts) {
    if (d != nullptr) cost += d->Cost();
  }
  return cost;
}

uint32_t ANSEncCounter::PutBit(uint32_t bit, uint32_t proba, WP2_OPT_LABEL) {
  cost_ += PROBA_BITS - WP2Log2(bit ? PROBA_MAX - proba : proba);
  return bit;
}

uint32_t ANSEncCounter::PutABit(uint32_t bit, ANSBinSymbol* const stats,
                                WP2_OPT_LABEL) {
  ANSEncCounter::PutBit(bit, stats->Proba(), label);   // will update cost_
  stats->Update(bit);
  return bit;
}

uint32_t ANSEncCounter::PutSymbol(uint32_t symbol, const ANSDictionary& dict,
                                  WP2_OPT_LABEL) {
  symbol_cost_ += dict.SymbolCost(symbol);
  return symbol;
}

uint32_t ANSEncCounter::PutSymbol(uint32_t symbol,
                                  const ANSAdaptiveSymbol& asym,
                                  WP2_OPT_LABEL) {
  cost_ += asym.GetCost(symbol);
  return symbol;
}

uint32_t ANSEncCounter::PutSymbol(uint32_t symbol, uint32_t max_symbol,
                                  const ANSDictionary& dict, WP2_OPT_LABEL) {
  symbol_cost_ += dict.SymbolCost(symbol, max_symbol);
  return symbol;
}

uint32_t ANSEncCounter::PutSymbol(uint32_t symbol, uint32_t max_symbol,
                                  const ANSAdaptiveSymbol& asym,
                                  WP2_OPT_LABEL) {
  cost_ += asym.GetCost(symbol, max_symbol);
  return symbol;
}

uint32_t ANSEncCounter::PutUValue(uint32_t value, uint32_t bits,
                                  WP2_OPT_LABEL) {
  cost_ += bits;
  return value;
}

uint32_t ANSEncCounter::PutRValue(uint32_t value, uint32_t range,
                                  WP2_OPT_LABEL) {
  cost_ += WP2Log2(range);
  return value;
}

WP2Status ANSEncCounter::Append(const ANSEncCounter& enc) {
  cost_ += enc.cost_;
  return WP2_STATUS_OK;
}

void ANSEncCounter::Reset() {
  cost_ = 0;
  symbol_cost_ = 0;
}

//------------------------------------------------------------------------------
// Mapping utils.

enum class MappingMethod {
  kRange = 0,
  kStoreVector,
  kAdaptativeBit,
  kBit,
};
static constexpr uint32_t kNumMapping = 4;

WP2Status LoadMapping(ANSDec* const dec, uint32_t size, uint32_t range,
                      Vector_u16& mapping) {
  assert(size <= range);
  WP2_CHECK_ALLOC_OK(mapping.resize(size));
  if (size == 0) return WP2_STATUS_OK;
  ANSDebugPrefix prefix(dec, "mapping");
  // TODO(vrabaud) Have it be a ReadBit with learned probability.
  if (size == range || dec->ReadBool("is_series")) {
    std::iota(mapping.begin(), mapping.end(), 0);
    return WP2_STATUS_OK;
  }

  const auto method = (MappingMethod)dec->ReadRValue(kNumMapping, "method");
  switch (method) {
    case MappingMethod::kRange: {
      // Read the runs, alternating between the beginning and end ones.
      uint32_t i = 0, j = size;
      for (uint32_t num_left = size; num_left > 0; ++i) {
        if (i == 0) {
          mapping[i] = dec->ReadRange(0, range - num_left, "run_beg");
        } else {
          mapping[i] = dec->ReadRange(mapping[i - 1] + 1, mapping[j] - num_left,
                                      "run_beg");
        }
        --num_left;
        if (num_left > 0) {
          --j;
          mapping[j] = dec->ReadRange(mapping[i] + num_left,
                                      (i == 0) ? range - 1 : mapping[j + 1] - 1,
                                      "run_end");
          --num_left;
          // early out to avoid no-op calls to ReadRange()
          if (mapping[i] + num_left + 1 == mapping[j]) break;
        }
      }
      while (++i < j) mapping[i] = mapping[i - 1] + 1;
      break;
    }
    case MappingMethod::kAdaptativeBit: {
      uint32_t num_zeros = range - size;
      for (uint32_t i = 0, j = 0; i < size; ++j) {
        const uint32_t num_left = range - j;
        if (num_zeros == 0 ||
            dec->ReadBit(num_zeros * PROBA_MAX / num_left, "is_used")) {
          mapping[i] = j;
          ++i;
        } else {
          --num_zeros;
        }
      }
      break;
    }
    case MappingMethod::kBit: {
      uint32_t num_zeros = range - size;
      for (uint32_t i = 0, j = 0; i < size; ++j) {
        if (num_zeros == 0 || dec->ReadBool("is_used")) {
          mapping[i] = j;
          ++i;
        } else {
          --num_zeros;
        }
      }
      break;
    }
    case MappingMethod::kStoreVector: {
      size_t k = 0;
      uint32_t n_bits_prev = 0;

      // 'range_left' is the range of the remaining runs, hence +1.
      uint32_t range_left = range - size + 1;
      uint32_t size_left = size;
      while (size_left > 0 && range_left > 1) {
        const uint32_t n_bits_max = FindLastSet(range_left - 1);
        // Figure out the number of bits.
        uint32_t n_bits;
        // The number of bits is in [0:n_bits_max] for the first set or when the
        // previous value was out of the current range.
        if (k == 0 || n_bits_prev > n_bits_max) {
          n_bits = dec->ReadRValue(n_bits_max + 1, "n_bits");
        } else {
           // If we know n_bits_prev is in [0:n_bits_max], we can save one bit
           // as we know the value has to be different from before.
          n_bits = dec->ReadRValue(n_bits_max, "n_bits");
          if (n_bits >= n_bits_prev) ++n_bits;
        }
        n_bits_prev = n_bits;

        // Figure out the number of elements with that bit depth.
        const uint32_t n = dec->ReadRange(1, size_left, "count");
        if (n_bits == 0) {
          // For a depth of 0 bits, we know it is a 0.
          for (size_t i = 0; i < n; ++i) {
            mapping[k] = (k == 0) ? 0 : mapping[k - 1] + 1;
           ++k;
          }
        } else if (n == 1) {
          // If we only have one number, we didn't store the high order
          // bit as it's always 1, otherwise it would fit in (n_bits_ - 1).
          uint32_t delta;
          if (n_bits == 1) {
            // And if that number is only on one bit to begin with, we didn't
            // store anything.
            delta = 1u;
          } else {
            delta = 1u << (n_bits - 1);
            assert(range_left >= delta);
            // After this value, the range will be at  most range_left - delta.
            // If it is smaller than the number of bits we use, we can just use
            // a range.
            if (range_left - delta < delta) {
              delta |= dec->ReadRValue(range_left - delta, "value");
            } else {
              delta |= dec->ReadUValue(n_bits - 1, "value");
            }
          }
          mapping[k] = (k == 0) ? delta : mapping[k - 1] + delta + 1;
          ++k;
          assert(range_left >= delta);
          range_left -= delta;
        } else {
          for (size_t i = 0; i < n; ++i) {
            uint32_t value;
            if (range_left < (1u << n_bits)) {
              value = dec->ReadRValue(range_left, "value");
            } else {
              value = dec->ReadUValue(n_bits, "value");
            }
            mapping[k] = (k == 0) ? value : mapping[k - 1] + value + 1;
            ++k;
            assert(range_left >= value);
            range_left -= value;
            // Stop here if everything else is in [0, 1) hence is 0.
            if (range_left == 1) break;
          }
        }
        assert(size_left >= n);
        size_left -= n;
      }

      // If we have not filled 'mapping' yet, it means we only have 0's left
      // for the runs.
      assert(k > 0);
      std::iota(&mapping[k], &mapping[size], mapping[k - 1] + 1);
      break;
    }
  }
  return WP2_STATUS_OK;
}

// Returns the effective reduced range of mapping values that are not trivial.
static uint32_t GetEffectiveRange(const uint16_t* const mapping,
                                  uint32_t size, uint32_t range) {
  uint32_t j = 0;
  for (uint32_t i = 0; i < size; ++i) {
    j = mapping[i] + 1;
    if (j > range - size + i) return range - size + i;
  }
  return j;
}

// Helper function that is a copy/paste of StoreVector with a few tweaked
// values.
static void StoreVectorForMapping(const uint16_t* const mapping, size_t size,
                                  uint32_t range,
                                  const OptimizeArrayStorageStat* const stats,
                                  uint32_t stats_size, ANSEncBase* const enc) {
  assert(size <= range);
  uint32_t n_bits_prev = 0;
  size_t size_left = size;
  uint32_t range_left = range - size + 1;
  for (uint32_t ind = 0, i = 0; i < stats_size; ++i) {
    assert(range_left > 1);
    assert(ind < size);
    // We need n_bits_max to store what is left.
    const uint32_t n_bits_max = FindLastSet(range_left - 1);
    const OptimizeArrayStorageStat& s = stats[i];
    // The number of bits is in [0:n_bits_max] for the first set or when the
    // previous value was out of the current range.
    if (ind == 0 || n_bits_prev > n_bits_max) {
      enc->PutRValue(s.n_bits, n_bits_max + 1, "n_bits");
    } else {
      // If we know n_bits_prev is in [0:n_bits_max], we can save one bit as we
      // know the value has to be different from before.
      // If we store n_bits0 for a set, the following n_bits1 will be stored
      // as:
      //   n_bits1       if n_bits1 < n_bits0,
      //   n_bits1 - 1   otherwise.
      if (s.n_bits < n_bits_prev) {
        enc->PutRValue(s.n_bits, n_bits_max, "n_bits");
      } else {
        enc->PutRValue(s.n_bits - 1, n_bits_max, "n_bits");
      }
    }
    n_bits_prev = s.n_bits;
    // Store the number of values.
    enc->PutRange(s.count, 1, size_left, "count");

    if (s.n_bits == 0) {
      ind += s.count;
    } else if (s.count == 1) {
      // If we only have one number, we don't need to store the high order
      // bit as it's always 1, otherwise it would fit in (n_bits_ - 1).
      // And if that number is only on one bit to begin with, we don't need
      // to store anything.
      const uint16_t val =
          (ind == 0) ? mapping[ind] : mapping[ind] - mapping[ind - 1] - 1;
      if (s.n_bits > 1) {
        const uint32_t high_order_bit = 1u << (s.n_bits - 1);
        const uint32_t value_to_code = val ^ high_order_bit;
        // After this value, the range will be at
        // most range_left - high_order_bit. If it is smaller than the
        // number of bits we use, we can just use a range.
        if ((uint32_t)(range_left - high_order_bit) < high_order_bit) {
          enc->PutRValue(value_to_code, range_left - high_order_bit, "value");
        } else {
          enc->PutUValue(value_to_code, s.n_bits - 1, "value");
        }
      }
      range_left -= val;
      ++ind;
    } else {
      // Store the values.
      for (size_t j = 0; j < s.count; ++j) {
        const uint16_t val =
            (ind == 0) ? mapping[ind] : mapping[ind] - mapping[ind - 1] - 1;
        if (range_left < (1u << s.n_bits)) {
          enc->PutRValue(val, range_left, "value");
        } else {
          enc->PutUValue(val, s.n_bits, "value");
        }
        range_left -= val;
        ++ind;
      }
    }
    size_left -= s.count;
  }
}

// Helper function storing the mapping with one of the different methods.
static void StoreMappingHelper(const uint16_t* const mapping, size_t size,
                               uint32_t range,
                               const OptimizeArrayStorageStat* const stats,
                               size_t stats_size, MappingMethod method,
                               ANSEncBase* const enc) {
  switch (method) {
    case MappingMethod::kRange: {
      // Store everything as a range by alternating between the beginning and
      // end runs to get 'num_left' as small as possible.
      // This method is usually good for few values in a big range.
      for (uint32_t i = 0, j = size, num_left = size; num_left > 0; ++i) {
        // Deal with the beginning runs.
        if (i == 0) {
          enc->PutRange(mapping[i], 0, range - num_left, "run_beg");
        } else {
          enc->PutRange(mapping[i], mapping[i - 1] + 1,
                        mapping[j] - num_left, "run_beg");
        }
        --num_left;
        // Deal with the end runs.
        if (num_left > 0) {
          --j;
          enc->PutRange(mapping[j], mapping[i] + num_left,
                        (i == 0) ? range - 1 : mapping[j + 1] - 1,
                        "run_end");
          --num_left;
          // early-out to avoid empty calls to PutRange()
          if (mapping[i] + num_left + 1 == mapping[j]) break;
        }
      }
      break;
    }
    case MappingMethod::kStoreVector: {
      // Store eveything as a list of tuples like in StoreVector.
      // This method is usually good for many values in a big range.
      StoreVectorForMapping(mapping, size, range, stats, stats_size, enc);
      break;
    }
    case MappingMethod::kAdaptativeBit: {
      // Store everything as a succession of probability bits indicating whether
      // the index is used.
      // This method is usually good for many values in a small range.
      uint32_t num_zeros = range - size;
      // Once there are no more zeros, it will cost nothing to store the bit.
      for (uint32_t i = 0, j = 0; i < size; ++j) {
        const uint32_t num_left = range - j;
        if (enc->PutBit(j == mapping[i],
                        num_zeros * PROBA_MAX / num_left, "is_used")) {
          ++i;
        } else {
          // early out to avoid calling PutBit() with proba 0
          if (--num_zeros == 0) break;
        }
      }
      break;
    }
    case MappingMethod::kBit: {
      // Store everything as a succession of bits indicating whether the index
      // is used.
      // This method is usually good for few values in a small range, and
      // located at the beginning.
      uint32_t num_zeros = range - size;
      // We keep going until we cannot deduce what is left (i.e. as long as
      // there is a 1 or 0).
      for (uint32_t i = 0, j = 0; i < size && num_zeros > 0; ++j) {
        if (enc->PutBool(j == mapping[i], "is_used")) {
          ++i;
        } else {
          --num_zeros;
        }
      }
      break;
    }
  }
}

// stats[] stores number of bits, and numbers of consecutive values with
// that bit depth.
static uint32_t CollectOptimalStats(
    const uint16_t* const mapping, size_t size, uint32_t range,
    OptimizeArrayStorageStat* const stats) {
  size_t stats_size = 0;
  uint16_t val_max = 0, val_min = 0;
  uint32_t range_left = range - size + 1;
  uint16_t prev_val = 0;
  for (size_t i = 0; i < size; ++i) {
    assert(mapping[i] < range);
    const uint16_t val = mapping[i] - prev_val;
    prev_val = mapping[i] + 1;
    // Check if the highest set bit is the same (faster than calling
    // FindLastSet).
    if (val < val_max && val >= val_min) {
      ++stats[stats_size - 1].count;
    } else {
      stats[stats_size].count = 1;
      stats[stats_size].n_bits = FindLastSet(val);
      val_max = (1u << stats[stats_size].n_bits);
      val_min = (val_max >> 1);
      ++stats_size;
    }
    // Stop here as everything else is in [0, 1) hence is 0.
    range_left -= val;
    if (range_left == 1) break;
  }
  // Merge the pairs to give an optimal cost.
  const uint8_t n_bits_max = FindLastSet(range_left - 1);
  assert(n_bits_max <= kANSMaxRangeBits);
  OptimizeArrayStorage(stats, &stats_size,
                       WP2Log2Fast(n_bits_max) + WP2Log2(size + 1));
  return stats_size;
}

float StoreMapping(const uint16_t* const mapping, size_t size, uint32_t range,
                   OptimizeArrayStorageStat* const stats,
                   ANSEncBase* const enc) {
  if (size == 0 || size == range) return 0.f;
  // Check whether all is consecutive, starting at 0.
  const bool is_series = (mapping[0] == 0 && mapping[size - 1] == size - 1);
  ANSDebugPrefix prefix(enc, "mapping");
  if (enc != nullptr) {
    enc->PutBool(is_series, "is_series");
  }
  if (is_series) return 1.f;

  WP2MathInit();  // Initialization for WP2Log2.

  // Find the optimal stats[] sequence for the kStoreVector case
  const size_t stats_size = CollectOptimalStats(mapping, size, range, stats);

  // Find most cost-efficient mapping method and finalize bitstream
  // if enc is non-null. Since kBit cost is easy to compute, use it as
  // initial value.
  float cost_min = GetEffectiveRange(mapping, size, range);
  MappingMethod best_method = MappingMethod::kBit;
  for (MappingMethod method : { MappingMethod::kRange,
                                MappingMethod::kStoreVector,
                                MappingMethod::kAdaptativeBit }) {
    ANSEncCounter counter;
    StoreMappingHelper(mapping, size, range, stats, stats_size,
                       method, &counter);
    const float cost = counter.GetCost();
    if (cost < cost_min) {
      cost_min = cost;
      best_method = method;
    }
  }

  if (enc != nullptr) {
    enc->PutRValue((uint32_t)best_method, kNumMapping, "method");
    StoreMappingHelper(mapping, size, range, stats, stats_size, best_method,
                       enc);
  }

  // Add the bit for 'is_series'.
  return 1.f + cost_min;
}

//------------------------------------------------------------------------------

uint32_t PutLargeRange(uint32_t v, uint32_t min, uint32_t max,
                       ANSEncBase* const enc, WP2_OPT_LABEL) {
  assert(v >= min);
  assert(v <= max);
  const uint32_t range = max - min + 1;
  if (range <= kANSMaxRange) return enc->PutRange(v, min, max, label);

  v = v - min;
  const uint32_t num_intervals = DivCeil(range, kANSMaxRange);
  const uint32_t interval = v / kANSMaxRange;
  if (range % kANSMaxRange == 0) {
    enc->PutRValue(v % kANSMaxRange, kANSMaxRange, label);
    enc->PutRValue(v / kANSMaxRange, num_intervals, label);
  } else {
    const bool is_last_interval = (interval == (num_intervals - 1));
    enc->PutBit(is_last_interval,
                PROBA_MAX - (range % kANSMaxRange) * PROBA_MAX / range, label);
    if (is_last_interval) {
      enc->PutRValue(v % kANSMaxRange, range % kANSMaxRange, label);
    } else {
      enc->PutRValue(v % kANSMaxRange, kANSMaxRange, label);
      enc->PutRValue(v / kANSMaxRange, num_intervals - 1, label);
    }
  }
  return min + v;
}

uint32_t ReadLargeRange(uint32_t min, uint32_t max, ANSDec* const dec,
                        WP2_OPT_LABEL) {
  const uint32_t range = max - min + 1;
  if (range <= kANSMaxRange) return dec->ReadRange(min, max, label);

  const uint32_t num_intervals = DivCeil(range, kANSMaxRange);
  if (range % kANSMaxRange == 0) {
    uint32_t v = dec->ReadRValue(kANSMaxRange, label);
    v += kANSMaxRange * dec->ReadRValue(num_intervals, label);
    return min + v;
  } else {
    const bool is_last_interval = dec->ReadBit(
        PROBA_MAX - (range % kANSMaxRange) * PROBA_MAX / range, label);
    if (is_last_interval) {
      return min + dec->ReadRValue(range % kANSMaxRange, label) +
             kANSMaxRange * (num_intervals - 1);
    } else {
      uint32_t v = dec->ReadRValue(kANSMaxRange, label);
      v += kANSMaxRange * dec->ReadRValue(num_intervals - 1, label);
      return min + v;
    }
  }
}

}  // namespace WP2
