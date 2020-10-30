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
//  Encoding of ANS symbols.
//
// Author: Vincent Rabaud (vrabaud@google.com)

#include "src/enc/symbols_enc.h"

#include <algorithm>

#include "src/dsp/math.h"
#include "src/utils/ans.h"

namespace WP2 {

WP2Status SymbolRecorder::MakeBackup() {
  // Save the stats to the previous pass dicts.
  for (uint32_t i = 0; i < dicts_.size(); ++i) {
    WP2_CHECK_STATUS(dicts_previous_pass_[i]->CopyFrom(*dicts_[i]));
  }
  return WP2_STATUS_OK;
}

void SymbolRecorder::ResetRecord(bool reset_backup) {
  // Save the stats to the previous pass dicts.
  if (reset_backup) {
    for (auto& d : dicts_previous_pass_) d->ResetCounts();
  }
  for (ANSDictionary* const d : dicts_) d->ResetCounts();

  uint32_t i = 0;
  for (uint32_t sym = 0; sym < symbols_info_.Size(); ++sym) {
    if (symbols_info_.Method(sym) != SymbolsInfo::StorageMethod::kAdaptiveBit) {
      continue;
    }
    for (uint32_t cluster = 0; cluster < symbols_info_.NumClusters(sym);
         ++cluster) {
      a_bits_[i] = ANSBinSymbol(symbols_info_.StartingProbaP0(cluster, sym),
                                symbols_info_.StartingProbaP1(cluster, sym));
      ++i;
    }
  }
  assert(i == a_bits_.size());

  for (ANSAdaptiveSymbol& sym : a_symbols_) {
    sym.InitFromUniform(sym.NumSymbols());
  }
  for (Vector_u8& v : values_) v.clear();
}

WP2Status SymbolRecorder::CopyFrom(const SymbolRecorder& other) {
  WP2_CHECK_STATUS(symbols_info_.CopyFrom(other.symbols_info_));
  index_ = other.index_;
  values_index_ = other.values_index_;
  WP2_CHECK_STATUS(dicts_.CopyFrom(other.dicts_));
  WP2_CHECK_STATUS(dicts_previous_pass_.CopyFrom(other.dicts_previous_pass_));
  WP2_CHECK_STATUS(a_bits_.CopyFrom(other.a_bits_));
  WP2_CHECK_STATUS(a_symbols_.CopyFrom(other.a_symbols_));
  WP2_CHECK_ALLOC_OK(values_.resize(other.values_.size()));
  for (uint32_t i = 0; i < values_.size(); ++i) {
    values_[i].clear();
    // Each sub vector is allocated only once with a fixed capacity.
    WP2_CHECK_ALLOC_OK(values_[i].reserve(other.values_[i].capacity()));
    WP2_CHECK_ALLOC_OK(values_[i].copy_from(other.values_[i]));
  }
  return WP2_STATUS_OK;
}

void SymbolRecorder::DeepClear() {
  index_ = {};
  dicts_.DeepClear();
  dicts_previous_pass_.DeepClear();
  a_bits_.clear();
  a_symbols_.clear();
  values_.clear();
}

WP2Status SymbolRecorder::Allocate(const SymbolsInfo& symbols_info,
                                   uint32_t num_records) {
  WP2_CHECK_STATUS(symbols_info_.CopyFrom(symbols_info));
  DeepClear();

  // Allocate dictionaries.
  for (uint32_t sym = 0; sym < symbols_info_.Size(); ++sym) {
    const uint32_t num_clusters = symbols_info_.NumClusters(sym);
    switch (symbols_info_.Method(sym)) {
      case SymbolsInfo::StorageMethod::kAuto: {
        index_[sym] = dicts_.size();
        for (uint32_t cluster = 0; cluster < num_clusters; ++cluster) {
          const uint32_t range = symbols_info_.Range(cluster, sym);
          if (range > 0) {
            WP2_CHECK_STATUS(dicts_.Add(range));
            WP2_CHECK_STATUS(dicts_previous_pass_.Add(range));
          }
        }
        break;
      }
      case SymbolsInfo::StorageMethod::kAdaptiveBit: {
        index_[sym] = a_bits_.size();
        for (uint32_t cluster = 0; cluster < num_clusters; ++cluster) {
          const ANSBinSymbol s(symbols_info_.StartingProbaP0(cluster, sym),
                               symbols_info_.StartingProbaP1(cluster, sym));
          WP2_CHECK_ALLOC_OK(a_bits_.push_back(s));
        }
        break;
      }
      case SymbolsInfo::StorageMethod::kAdaptiveSym: {
        index_[sym] = a_symbols_.size();
        for (uint32_t cluster = 0; cluster < num_clusters; ++cluster) {
          const uint32_t range = symbols_info_.Range(cluster, sym);
          if (range > 0) WP2_CHECK_ALLOC_OK(a_symbols_.Add(range));
        }
        break;
      }
      case SymbolsInfo::StorageMethod::kAdaptiveWithAutoSpeed: {
        // Must fit in a uint8_t
        assert(symbols_info_.GetMaxRange(sym) < (1 << 8));
        index_[sym] = a_symbols_.size();
        values_index_[sym] = values_.size();
        WP2_CHECK_ALLOC_OK(values_.resize(values_.size() + num_clusters));
        for (uint32_t cluster = 0; cluster < num_clusters; ++cluster) {
          const uint32_t range = symbols_info_.Range(cluster, sym);
          if (range > 0) {
            WP2_CHECK_ALLOC_OK(a_symbols_.Add(range));
            WP2_CHECK_ALLOC_OK(
                values_[values_index_[sym] + cluster].reserve(num_records));
          }
        }
        break;
      }
      case SymbolsInfo::StorageMethod::kUnused: {
        break;
      }
    }
  }
  return WP2_STATUS_OK;
}

const ANSDictionary& SymbolRecorder::GetRecordedDict(uint32_t cluster,
                                                     uint32_t sym) const {
  assert(symbols_info_.Range(cluster, sym) != 0);
  assert(symbols_info_.Method(sym) == SymbolsInfo::StorageMethod::kAuto);
  return *dicts_[index_[sym] + cluster];
}

const ANSDictionary& SymbolRecorder::GetDictPreviousPass(uint32_t cluster,
                                                         uint32_t sym) const {
  assert(symbols_info_.Range(cluster, sym) != 0);
  assert(symbols_info_.Method(sym) == SymbolsInfo::StorageMethod::kAuto);
  return *dicts_previous_pass_[index_[sym] + cluster];
}

ANSDictionary& SymbolRecorder::GetDictPreviousPass(uint32_t cluster,
                                                   uint32_t sym) {
  assert(symbols_info_.Range(cluster, sym) != 0);
  assert(symbols_info_.Method(sym) == SymbolsInfo::StorageMethod::kAuto);
  return *dicts_previous_pass_[index_[sym] + cluster];
}

ANSDictionary* SymbolRecorder::GetDict(uint32_t cluster, uint32_t sym) {
  assert(symbols_info_.Range(cluster, sym) != 0);
  assert(symbols_info_.Method(sym) == SymbolsInfo::StorageMethod::kAuto);
  return dicts_[index_[sym] + cluster];
}

const ANSBinSymbol& SymbolRecorder::GetABit(uint32_t cluster,
                                            uint32_t sym) const {
  assert(symbols_info_.Method(sym) == SymbolsInfo::StorageMethod::kAdaptiveBit);
  return a_bits_[index_[sym] + cluster];
}

ANSBinSymbol* SymbolRecorder::GetABit(uint32_t cluster, uint32_t sym) {
  assert(symbols_info_.Method(sym) == SymbolsInfo::StorageMethod::kAdaptiveBit);
  return &a_bits_[index_[sym] + cluster];
}

const ANSAdaptiveSymbol& SymbolRecorder::GetASymbol(uint32_t cluster,
                                                    uint32_t sym) const {
  assert(symbols_info_.Range(cluster, sym) != 0);
  assert(symbols_info_.Method(sym) == SymbolsInfo::StorageMethod::kAdaptiveSym);
  return a_symbols_[index_[sym] + cluster];
}

ANSAdaptiveSymbol* SymbolRecorder::GetASymbol(uint32_t cluster, uint32_t sym) {
  assert(symbols_info_.Range(cluster, sym) != 0);
  assert(symbols_info_.Method(sym) == SymbolsInfo::StorageMethod::kAdaptiveSym);
  return &a_symbols_[index_[sym] + cluster];
}

const ANSAdaptiveSymbol& SymbolRecorder::GetASymbolWithSpeed(
    uint32_t cluster, uint32_t sym) const {
  assert(symbols_info_.Range(cluster, sym) != 0);
  assert(symbols_info_.Method(sym) ==
         SymbolsInfo::StorageMethod::kAdaptiveWithAutoSpeed);
  return a_symbols_[index_[sym] + cluster];
}

const Vector_u8& SymbolRecorder::GetRecordedValues(uint32_t cluster,
                                                   uint32_t sym) const {
  assert(symbols_info_.Method(sym) ==
         SymbolsInfo::StorageMethod::kAdaptiveWithAutoSpeed);
  return values_[values_index_[sym] + cluster];
}

float SymbolRecorder::GetCost(uint32_t cluster, uint32_t sym,
                              int32_t value) const {
  switch (symbols_info_.Method(sym)) {
    case SymbolsInfo::StorageMethod::kAuto:
      return GetRecordedDict(cluster, sym)
          .SymbolCost(value - symbols_info().Min(cluster, sym));
    // TODO(yguyon): case SymbolsInfo::StorageMethod::kAdaptiveBit:
    case SymbolsInfo::StorageMethod::kAdaptiveSym:
      return GetASymbol(cluster, sym)
          .GetCost(value - symbols_info().Min(cluster, sym));
    case SymbolsInfo::StorageMethod::kAdaptiveWithAutoSpeed: {
      assert(symbols_info_.Range(cluster, sym) > 0);
      // TODO(yguyon): Could be more precise (use GetRecordedValues()).
      return WP2Log2(symbols_info_.Range(cluster, sym));
    }
    default:
      assert(false);
      return 0.f;
  }
}

int32_t SymbolRecorder::ProcessInternal(uint32_t cluster, uint32_t sym,
                                        int32_t value, bool use_max_value,
                                        uint32_t max_value, WP2_OPT_LABEL,
                                        ANSEncBase* const enc, double* const) {
  assert((uint32_t)std::abs(value) < symbols_info_.Range(cluster, sym));
  assert(!use_max_value || (uint32_t)std::abs(value) <= max_value);
  switch (symbols_info_.Method(sym)) {
    case SymbolsInfo::StorageMethod::kAuto:
      GetDict(cluster, sym)
          ->RecordSymbol(value - symbols_info().Min(cluster, sym));
      break;
    case SymbolsInfo::StorageMethod::kAdaptiveBit:
      GetABit(cluster, sym)->Update(value);
      break;
    case SymbolsInfo::StorageMethod::kAdaptiveSym:
      GetASymbol(cluster, sym)
          ->Update(value - symbols_info().Min(cluster, sym));
      break;
    case SymbolsInfo::StorageMethod::kAdaptiveWithAutoSpeed:
      a_symbols_[index_[sym] + cluster].Update(
          value - symbols_info().Min(cluster, sym));
      if (!values_[values_index_[sym] + cluster].push_back(
              value - symbols_info().Min(cluster, sym),
              /*resize_if_needed=*/false)) {
        // Fail if we get more values than expected.
        assert(false);
      }
      break;
    case SymbolsInfo::StorageMethod::kUnused:
      // Even though the object can be created with unused symbols (e.g. alpha
      // symbols), those should never be used.
      assert(false);
      break;
  }
  return value;
}

//------------------------------------------------------------------------------

constexpr uint32_t SymbolCounter::kInvalidIndex;

WP2Status SymbolCounter::Allocate(std::initializer_list<uint32_t> syms) {
  std::array<bool, kSymbolNumMax> is_used;
  is_used.fill(false);
  for (uint32_t sym : syms) is_used[sym] = true;

  num_a_bits_ = 0u;
  num_a_symbols_ = 0u;
  indices_.fill(kInvalidIndex);
  for (uint32_t sym = 0; sym < kSymbolNumMax; ++sym) {
    if (!is_used[sym]) continue;
    // In some partitioning (e.g. BlockScoreFunc), symbols_info is
    // actually empty, so choose the maximum possible.
    const uint32_t num_clusters = (symbols_info_.Size() == 0)
                                      ? kMaxLossyClusters
                                      : symbols_info_.NumClusters(sym);
    switch (symbols_info_.Method(sym)) {
      case SymbolsInfo::StorageMethod::kAdaptiveBit:
        indices_[sym] = num_a_bits_;
        num_a_bits_ += num_clusters;
        break;
      case SymbolsInfo::StorageMethod::kAdaptiveSym:
      case SymbolsInfo::StorageMethod::kAdaptiveWithAutoSpeed:
        indices_[sym] = num_a_symbols_;
        num_a_symbols_ += num_clusters;
        break;
      case SymbolsInfo::StorageMethod::kUnused:
        break;
      default:
        // In practice, it should be assert(false) but some symbols might have
        // a different method during tests.
        break;
    }
  }
  WP2_CHECK_ALLOC_OK(a_bits_.resize(num_a_bits_));
  WP2_CHECK_ALLOC_OK(a_symbols_.resize(num_a_symbols_));
  // If the following asserts ever break, just increase the size of the arrays.
  // For now, the sizes are tuned to our current usage.
  assert(num_a_bits_ <= a_bit_initialized_.size());
  assert(num_a_symbols_ <= a_symbol_initialized_.size());
  Clear();
  return WP2_STATUS_OK;
}

WP2Status SymbolCounter::CopyFrom(const SymbolCounter& other) {
  // Nothing to copy in base classes SymbolManager and WP2Allocable.
  indices_ = other.indices_;
  num_a_bits_ = other.num_a_bits_;
  a_bit_initialized_ = other.a_bit_initialized_;
  WP2_CHECK_STATUS(a_bits_.CopyFrom(other.a_bits_));
  num_a_symbols_ = other.num_a_symbols_;
  a_symbol_initialized_ = other.a_symbol_initialized_;
  WP2_CHECK_STATUS(a_symbols_.CopyFrom(other.a_symbols_));
  return WP2_STATUS_OK;
}

void SymbolCounter::Clear() {
  // Somehow, memset is faster than std::fill.
  memset(&a_bit_initialized_[0], 0, num_a_bits_ * sizeof(true));
  memset(&a_symbol_initialized_[0], 0, num_a_symbols_ * sizeof(true));
}

ANSBinSymbol* SymbolCounter::GetABit(uint32_t cluster, uint32_t sym) {
  assert(symbols_info_.Method(sym) ==
         SymbolsInfo::StorageMethod::kAdaptiveBit);
  assert(indices_[sym] != kInvalidIndex);
  const uint32_t ind = indices_[sym] + cluster;
  if (!a_bit_initialized_[ind]) {
    a_bits_[ind] = recorder_->GetABit(cluster, sym);
    a_bit_initialized_[ind] = true;
  }
  return &a_bits_[ind];
}

ANSAdaptiveSymbol* SymbolCounter::GetASymbol(uint32_t cluster, uint32_t sym) {
  assert(symbols_info_.Range(cluster, sym) != 0);
  const SymbolsInfo::StorageMethod method = symbols_info_.Method(sym);
  assert(method == SymbolsInfo::StorageMethod::kAdaptiveSym ||
         method == SymbolsInfo::StorageMethod::kAdaptiveWithAutoSpeed);
  assert(indices_[sym] != kInvalidIndex);
  const uint32_t ind = indices_[sym] + cluster;
  if (!a_symbol_initialized_[ind]) {
    if (method == SymbolsInfo::StorageMethod::kAdaptiveSym) {
      a_symbols_[ind] = recorder_->GetASymbol(cluster, sym);
    } else {
      // Adaptive symbol with auto speed. We approximate it as regular adaptive
      // symbol with the default speed.
      // TODO(maryla): not sure this is the best. Maybe just assume a range?
      a_symbols_[ind] = recorder_->GetASymbolWithSpeed(cluster, sym);
    }
    a_symbol_initialized_[ind] = true;
  }
  return &a_symbols_[ind];
}

int32_t SymbolCounter::ProcessInternal(uint32_t cluster, uint32_t sym,
                                       int32_t value_in, bool use_max_value,
                                       uint32_t max_value, WP2_OPT_LABEL,
                                       ANSEncBase* const enc, double* const) {
  assert((uint32_t)std::abs(value_in) < symbols_info_.Range(cluster, sym));
  assert(!use_max_value || (uint32_t)std::abs(value_in) <= max_value);
  switch (symbols_info_.Method(sym)) {
    case SymbolsInfo::StorageMethod::kAuto: {
      const ANSDictionary& current_dict =
          recorder_->GetRecordedDict(cluster, sym);
      const ANSDictionary& previous_dict =
          recorder_->GetDictPreviousPass(cluster, sym);

      // Whether the recorder has accumulated enough stats (by some arbitrary
      // criteria).
      const bool enough_stats =
          (current_dict.Total() > 50 ||
           current_dict.Total() > 10 * symbols_info_.Range(cluster, sym));
      const ANSDictionary* const dict =
          enough_stats ? &current_dict : &previous_dict;
      const uint32_t value = value_in - symbols_info().Min(cluster, sym);
      const bool symbol_present = (dict->Counts()[value] > 0);
      if (symbol_present) {
        // TODO(maryla): average the two dicts instead of sharply transitioning
        // between the two?
        if (use_max_value) {
          enc->PutSymbol(value, max_value, *dict, label);
        } else {
          enc->PutSymbol(value, *dict, label);
        }
      } else if (enough_stats || previous_dict.Total() > 0) {
        // If the symbol is not present but we have enough data, assume a proba
        // of 1/Total.
        // Emulate the cost of log2(1/Total) by using a range.
        enc->PutRValue(0, dict->Total(), label);
      } else {
        // Otherwise we assume all values have the same probability, i.e. it's
        // a range.
        uint32_t range = symbols_info_.Range(cluster, sym);
        if (use_max_value) range = std::min(range, max_value + 1);
        if (symbols_info_.ProbablyGeometric(sym)) {
          // Very approximate. Favors lower values. TODO(maryla): find a better
          // approximation of golomb cost.
          range = std::min(range, (value == 0) ? 2 : ((uint32_t)value << 2));
        }
        enc->PutRValue(value, range, label);
      }
      break;
    }
    case SymbolsInfo::StorageMethod::kAdaptiveBit:
      enc->PutABit(value_in, GetABit(cluster, sym), label);
      break;
    case SymbolsInfo::StorageMethod::kAdaptiveSym:
    case SymbolsInfo::StorageMethod::kAdaptiveWithAutoSpeed:
      enc->PutASymbol(value_in, GetASymbol(cluster, sym), label);
      break;
    case SymbolsInfo::StorageMethod::kUnused:
      // Even though the object can be created with unused symbols (e.g. alpha
      // symbols), those should never be used.
      assert(false);
      break;
  }
  return value_in;
}

//------------------------------------------------------------------------------

WP2Status SymbolWriter::Allocate() {
  WP2_CHECK_STATUS(SymbolIO<SymbolWriterStatExtra>::Allocate());
  const uint32_t symbol_range_max = symbols_info_.GetMaxRange();
  WP2_CHECK_ALLOC_OK(histogram_.resize(symbol_range_max));
  WP2_CHECK_ALLOC_OK(mapping_.resize(symbol_range_max));
  WP2_CHECK_ALLOC_OK(infos_tmp_.resize(symbol_range_max));
  // The following buffers are also used when storing recursive histograms
  // (storage type HuffmanANS). Assuming a symbol appears at most once per
  // pixel, the range of a huffmanized histogram is log2(max pixels).
  constexpr uint32_t kMaxHistoValueBits = WP2Log2Ceil_k(kMaxTilePixels);
  const uint32_t histo_max_range =
      std::max(symbol_range_max, kMaxHistoValueBits);
  WP2_CHECK_ALLOC_OK(stats_buffer_.resize(histo_max_range));
  WP2_CHECK_ALLOC_OK(counts_tmp_.reserve(histo_max_range));

  const Golomb golomb1(symbol_range_max - 1, 1);
  const Golomb golomb2(symbol_range_max - 1, 2);
  const uint32_t golomb_range_max =
      1 + std::max(golomb1.prefix, golomb2.prefix);
  WP2_CHECK_ALLOC_OK(histogram_golomb_.resize(golomb_range_max));
  WP2_CHECK_ALLOC_OK(mapping_golomb_.resize(golomb_range_max));
  WP2_CHECK_STATUS(quantizer_.Allocate(symbol_range_max));
  WP2_CHECK_STATUS(quantizer_golomb_.Allocate(golomb_range_max));
  return WP2_STATUS_OK;
}

WP2Status SymbolWriter::CopyFrom(const SymbolWriter& other,
                                 const ANSDictionaries& original_dicts,
                                 const ANSDictionaries& copied_dicts) {
  // Nothing to copy in base classes SymbolManager and WP2Allocable.

  // Deep copy base class SymbolIO<SymbolWriterStatExtra>. This is easier to do
  // it here rather than in SymbolIO because the template is known.
  WP2_CHECK_STATUS(symbols_info_.CopyFrom(other.symbols_info_));
  WP2_CHECK_STATUS(Allocate());  // Needed for 'stats_start_', Stat::mappings

  WP2_CHECK_STATUS(a_bits_.CopyFrom(other.a_bits_));
  WP2_CHECK_STATUS(a_symbols_.CopyFrom(other.a_symbols_));
  // 'all_stats_' cannot be copied as is, it contains references.
  assert(all_stats_.size() == other.all_stats_.size());
  for (uint32_t i = 0; i < all_stats_.size(); ++i) {
    Stat& stat = all_stats_[i];
    const Stat& other_stat = other.all_stats_[i];

    stat.type = other_stat.type;
    stat.use_mapping = other_stat.use_mapping;
    // Keep the Stat::mappings that was assigned during Allocate().
    stat.range = other_stat.range;
    stat.param = other_stat.param;

    if (stat.type != SymbolWriter::Stat::Type::kUnknown) {
      SymbolWriterStatExtra& extra = stat.extra;
      const SymbolWriterStatExtra& other_extra = other_stat.extra;

      extra.dict = copied_dicts.GetEquivalent(original_dicts, other_extra.dict);
      assert((extra.dict != nullptr) == (other_extra.dict != nullptr));
      extra.mapping_size = other_extra.mapping_size;
    }
  }
  // 'mappings_buffer_' cannot be reallocated, it is referenced in Allocate().
  assert(mappings_buffer_.size() == other.mappings_buffer_.size());
  std::copy(other.mappings_buffer_.begin(), other.mappings_buffer_.end(),
            mappings_buffer_.begin());

  // Deep copy remaining direct members of SymbolWriter.
  WP2_CHECK_ALLOC_OK(infos_tmp_.copy_from(other.infos_tmp_));
  WP2_CHECK_ALLOC_OK(counts_tmp_.copy_from(other.counts_tmp_));
  WP2_CHECK_ALLOC_OK(histogram_.copy_from(other.histogram_));
  WP2_CHECK_ALLOC_OK(mapping_.copy_from(other.mapping_));
  WP2_CHECK_ALLOC_OK(histogram_golomb_.copy_from(other.histogram_golomb_));
  WP2_CHECK_ALLOC_OK(mapping_golomb_.copy_from(other.mapping_golomb_));
  // No need to copy Quantizer instances. They only contain data members to
  // avoid reallocations. TODO(yguyon): Check if more could be skipped
  WP2_CHECK_ALLOC_OK(stats_buffer_.copy_from(other.stats_buffer_));
  return WP2_STATUS_OK;
}

//------------------------------------------------------------------------------

// Function writing the histogram to the encoder.
WP2Status SymbolWriter::WriteHistogram(const Quantizer::Config& config,
                                       uint32_t symbol_range,
                                       uint32_t max_count,
                                       ANSEncBase* const enc) {
  // Store the mapping first.
  const auto& param = config.param_;
  if (config.h_.nnz < symbol_range) {
    enc->PutBool(param.is_sparse_, "is_sparse");
    if (param.is_sparse_) {
      // The diff between two consecutive elements is at most symbol_range-nbr
      // of symbols.
      assert(config.h_.nnz <= stats_buffer_.size());
      StoreMapping(config.h_.mapping, config.h_.nnz, symbol_range,
                   stats_buffer_.data(), enc);
    } else {
      enc->PutRange(config.size_to_write_, config.h_.nnz, symbol_range,
                    "histogram_size");
    }
  } else {
    // We use the full range so no need to do anything.
    assert(param.is_sparse_);
  }
  if (param.is_sparse_) assert(config.size_to_write_ == config.h_.nnz);
  assert(config.size_to_write_ >= config.h_.nnz);

  enc->PutRValue(param.type_, 3, "coding_type");
  // Save the probabilities.
  const uint32_t max_count_bits =
      std::min(kMaxFreqBits, 1u + (uint32_t)WP2Log2Floor(max_count));
  for (uint32_t i = 0; i < config.h_.nnz; ++i) {
    assert(config.h_.counts[i] < (1u << (max_count_bits + 1)));
  }
  if (param.type_ == Quantizer::HuffmanANS) {
    // First, save the following probabilities.
    assert(config.next_ != nullptr);
    const auto& config_next = *config.next_;
    enc->PutRange(config_next.h_.nnz, 1, 1 + kMaxFreqBits, "size");
    WP2_CHECK_STATUS(
        WriteHistogram(config_next, 1 + kMaxFreqBits,
        config.size_to_write_, enc));
    // Convert to an histogram.
    const uint32_t h_size = config_next.h_.mapping[config_next.h_.nnz - 1] + 1;
    ANSDictionary dict;
    WP2_CHECK_STATUS(dict.Init(h_size));
    for (size_t k = 0; k < config_next.h_.nnz; ++k) {
      dict.RecordSymbol(config_next.h_.mapping[k], config_next.h_.counts[k]);
    }
    WP2_CHECK_STATUS(dict.ToCodingTable());

    uint32_t nnz_left = config.h_.nnz;
    const uint32_t* const k_end =
        config.histogram_to_write_ + config.size_to_write_;
    for (const uint32_t* k = config.histogram_to_write_;
         nnz_left > 0 && k < k_end; ++k) {
      if (!param.is_sparse_ && k + nnz_left == k_end) {
        // If we only have non-zeros left, we can store counts - 1.
        assert(*k > 0);
        enc->PutSymbol(*k - 1, dict, "counts");
        --nnz_left;
      } else {
        enc->PutSymbol(*k, dict, "counts");
        if (*k != 0) --nnz_left;
      }
    }
  } else {
    assert(param.max_freq_bits_ >= 1);
    if (param.type_ == Quantizer::Huffman) {
      enc->PutRange(param.max_freq_bits_, 1, 1 + WP2Log2Floor(max_count_bits),
                    "max_freq_bits");
    } else {
      enc->PutRange(param.max_freq_bits_, 1, max_count_bits, "max_freq_bits");
    }
    assert(config.size_to_write_ <= stats_buffer_.size());
    if (param.is_sparse_) {
      StoreVector(config.histogram_to_write_, config.size_to_write_,
                  (1 << param.max_freq_bits_) - 1, stats_buffer_.data(), enc);
    } else {
      StoreVectorNnz(config.histogram_to_write_, config.size_to_write_,
                     config.h_.nnz, (1 << param.max_freq_bits_) - 1,
                     stats_buffer_.data(), enc);
    }
  }
  return WP2_STATUS_OK;
}

//------------------------------------------------------------------------------

static bool ANSSymbolInfoCmp(const ANSSymbolInfo& i1, const ANSSymbolInfo& i2) {
  return (i1.symbol < i2.symbol);
}

void SymbolWriter::AddTrivial(uint32_t cluster, uint32_t sym, int32_t value) {
  Stat* const stat = GetStats(cluster, sym);

  // No dictionary for trivial symbols.
  stat->type = Stat::Type::kTrivial;
  stat->param.trivial_value = value;
  stat->extra.dict = nullptr;
  stat->extra.mapping_size = 0;
}

void SymbolWriter::AddRange(uint32_t cluster, uint32_t sym,
                            const uint16_t* const mapping, uint32_t size,
                            uint16_t range) {
  Stat* const stat = GetStats(cluster, sym);

  stat->type = Stat::Type::kRange;
  stat->use_mapping = (mapping != nullptr);
  stat->extra.dict = nullptr;
  stat->extra.mapping_size = 0;
  if (stat->use_mapping) {
    stat->range = size;
    stat->extra.mapping_size =
        1 + *std::max_element(&mapping[0], &mapping[0] + size);
    std::fill(stat->mappings, stat->mappings + stat->extra.mapping_size,
              Stat::kInvalidMapping);
    for (uint32_t i = 0; i < size; ++i) {
      stat->mappings[mapping[i]] = i;
    }
  } else {
    stat->range = range - size;
  }
}

WP2Status SymbolWriter::AddDict(uint32_t cluster, uint32_t sym,
                                const uint32_t* const counts,
                                const uint32_t* const quantized_counts,
                                const uint16_t* const mapping, size_t size,
                                ANSDictionaries* const dicts) {
  assert(mapping != nullptr);
  assert(quantized_counts != nullptr);

  Stat* const stat = GetStats(cluster, sym);

  stat->type = Stat::Type::kDict;
  stat->use_mapping = true;

  // Sort with the most frequent symbols first.
  for (uint32_t i = 0; i < size; ++i) {
    infos_tmp_[i].freq = quantized_counts[i];
    assert(infos_tmp_[i].freq == quantized_counts[i]); /* no overflow */
    infos_tmp_[i].symbol = mapping[i];
  }
  std::sort(&infos_tmp_[0], &infos_tmp_[0] + size, ANSSymbolInfoCmp);
  // Figure out the mappings.
  stat->extra.mapping_size = 0;
  for (uint32_t k = 0; k < size; ++k) {
    if (infos_tmp_[k].symbol >= stat->extra.mapping_size) {
      stat->extra.mapping_size = 1 + infos_tmp_[k].symbol;
    }
  }
  std::fill(stat->mappings, stat->mappings + stat->extra.mapping_size,
            Stat::kInvalidMapping);
  WP2_CHECK_ALLOC_OK(counts_tmp_.resize(size));
  for (uint32_t k = 0; k < size; ++k) {
    stat->mappings[infos_tmp_[k].symbol] = k;
    counts_tmp_[k] = infos_tmp_[k].freq;
  }
  WP2_CHECK_STATUS(dicts->Add(size));
  auto* const d = dicts->back();
  // In case the counts were already quantized.
  for (uint32_t k = 0; k < size; ++k) {
    d->RecordSymbol(stat->mappings[mapping[k]], counts[k]);
  }
  WP2_CHECK_STATUS(d->SetQuantizedCounts(counts_tmp_));
  WP2_CHECK_STATUS(d->ToCodingTable());
  stat->extra.dict = d;
  return WP2_STATUS_OK;
}

WP2Status SymbolWriter::AddGolomb(uint32_t cluster, uint32_t sym,
                                  const uint32_t* const counts,
                                  const uint32_t* const quantized_counts,
                                  const uint16_t* const mapping, uint32_t size,
                                  uint32_t prefix_size,
                                  ANSDictionaries* const dicts) {
  WP2_CHECK_STATUS(
      AddDict(cluster, sym, counts, quantized_counts, mapping, size, dicts));
  SetGolombStat(cluster, sym, prefix_size);
  return WP2_STATUS_OK;
}

uint32_t SymbolWriter::FindBiggestMappingIndex(const Stat& stat,
                                               uint32_t max_index) {
  const uint32_t index_ini =
      std::min(max_index, (uint32_t)(stat.extra.mapping_size - 1));
  uint32_t index = index_ini;
  while (index > 0 && stat.mappings[index] == Stat::kInvalidMapping) {
    --index;
  }
  if (index == 0 && stat.mappings[index] == Stat::kInvalidMapping) {
    index = index_ini;
    while (index < stat.extra.mapping_size &&
           stat.mappings[index] == Stat::kInvalidMapping) {
      ++index;
    }
    assert(index < stat.extra.mapping_size);
  }
  assert(stat.mappings[index] != Stat::kInvalidMapping);
  return index;
}

int32_t SymbolWriter::ProcessInternal(uint32_t cluster, uint32_t sym,
                                      int32_t value_in, bool use_max_value,
                                      uint32_t max_value, WP2_OPT_LABEL,
                                      ANSEncBase* const enc,
                                      double* const cost) {
  if (use_max_value) assert(value_in <= (int32_t)max_value);
  const Stat& stat = *GetStats(cluster, sym);
  // Only write a symbol if it is not trivial.
  if (stat.type == Stat::Type::kTrivial) {
    assert(value_in == stat.param.trivial_value);
    return value_in;
  }
  ANSDebugPrefix debug_prefix(enc, label);
  const uint32_t value = (value_in < 0) ? -value_in : value_in;
  switch (stat.type) {
    case Stat::Type::kRange: {
      uint32_t range;
      if (!use_max_value) {
        range = stat.range;
      } else if (stat.use_mapping) {
        // Find the biggest range <= max_value that is valid.
        const uint32_t max_index = FindBiggestMappingIndex(stat, max_value);
        range = stat.mappings[max_index] + 1;
        assert(range <= stat.range);
      } else {
        range = std::min(stat.range, (uint16_t)(max_value + 1));
      }
      enc->PutRValue(stat.use_mapping ? stat.mappings[value] : value, range,
                     "range");
      if (cost != nullptr) {
        *cost += std::log2(range);
      }
      break;
    }
    case Stat::Type::kDict: {
      if (use_max_value) {
        // Get the info of the maximally usable symbol.
        max_value =
            std::min(max_value, (uint32_t)(stat.extra.mapping_size - 1));
        while (max_value > 0 &&
               stat.mappings[max_value] == Stat::kInvalidMapping) {
          --max_value;
        }
        assert(stat.mappings[max_value] != Stat::kInvalidMapping);
        enc->PutSymbol(stat.mappings[value], stat.mappings[max_value],
                       *stat.extra.dict, "dict");
        if (cost != nullptr) {
          *cost += stat.extra.dict->SymbolCost(stat.mappings[value],
                                               stat.mappings[max_value]);
        }
      } else {
        enc->PutSymbol(stat.mappings[value], *stat.extra.dict, "dict");
        if (cost != nullptr) {
          *cost += stat.extra.dict->SymbolCost(stat.mappings[value]);
        }
      }
      break;
    }
    case Stat::Type::kGolomb: {
      // TODO(vrabaud) restrict when use_max_value is set.
      const uint32_t prefix_size = stat.param.golomb.prefix_size;
      const Golomb golomb(value, prefix_size);
      bool use_range;
      uint32_t range = 0;
      if (use_max_value) {
        const Golomb golomb_max(max_value, prefix_size);
        const uint32_t max_prefix =
            FindBiggestMappingIndex(stat, golomb_max.prefix);
        enc->PutSymbol(stat.mappings[golomb.prefix], stat.mappings[max_prefix],
                       *stat.extra.dict, "golomb");
        use_range = (golomb.prefix == max_prefix);
        if (use_range) {
          // Get (the biggest value with a prefix of 'max_prefix') + 1.
          range = std::min(
              max_value + 1,
              Golomb::Merge(max_prefix, prefix_size, 0) +
                  (1 << Golomb::NumExtraBits(max_prefix, prefix_size)));
        }
        if (cost != nullptr) {
          *cost += stat.extra.dict->SymbolCost(stat.mappings[golomb.prefix],
                                               stat.mappings[max_prefix]);
        }
      } else {
        enc->PutSymbol(stat.mappings[golomb.prefix], *stat.extra.dict,
                       "golomb");
        // If the last bit interval is truncated, go with ranges.
        use_range = (golomb.prefix + 1 == stat.param.golomb.range);
        if (use_range) range = stat.range;
        if (cost != nullptr) {
          *cost += stat.extra.dict->SymbolCost(stat.mappings[golomb.prefix]);
        }
      }
      if (golomb.extra_bits_num > 0) {
        if (use_range) {
          assert(Golomb::Merge(golomb.prefix, prefix_size, 0) < range);
          range -= Golomb::Merge(golomb.prefix, prefix_size, 0);
          enc->PutRValue(golomb.extra_bits_value, range, "extra_bits_value");
          if (cost != nullptr) {
            *cost += std::log2(range);
          }
        } else {
          enc->PutUValue(golomb.extra_bits_value, golomb.extra_bits_num,
                         "extra_bits_value");
          if (cost != nullptr) {
            *cost += golomb.extra_bits_num;
          }
        }
      }
      break;
    }
    case Stat::Type::kAdaptiveBit: {
      if (!use_max_value || max_value > 0) {
        enc->PutBit(value, a_bits_[stat.param.a_bit_index].Proba(), "a_bit");
        if (cost != nullptr) {
          *cost += a_bits_[stat.param.a_bit_index].GetCost(value);
        }
      }
      a_bits_[stat.param.a_bit_index].Update(value);
      break;
    }
    case Stat::Type::kAdaptiveSymbol: {
      if (use_max_value) {
        enc->PutSymbol(value, max_value, a_symbols_[stat.param.a_symbol_index],
                       "a_symbol");
        if (cost != nullptr) {
          *cost +=
              a_symbols_[stat.param.a_symbol_index].GetCost(value, max_value);
        }
      } else {
        enc->PutSymbol(value, a_symbols_[stat.param.a_symbol_index],
                       "a_symbol");
        if (cost != nullptr) {
          *cost += a_symbols_[stat.param.a_symbol_index].GetCost(value);
        }
      }
      // TODO(vrabaud) take the max_value into account when updating.
      a_symbols_[stat.param.a_symbol_index].Update(value);
      break;
    }
    default:
      // Did you forget to call WriteHeader for this symbol?
      assert(false);
      break;
  }
  // Store the sign.
  if (value_in != 0 && symbols_info().Min(cluster, sym) < 0) {
    enc->PutBool(value_in < 0, "is_negative");
    if (cost != nullptr) *cost += 1.;
  }
  return value_in;
}

void SymbolWriter::ComputeGolombCost(uint32_t range, uint32_t max_nnz,
                                     uint32_t size, uint32_t prefix_size,
                                     float* cost,
                                     Quantizer::Config** const config_golomb,
                                     uint32_t* const size_golomb) {
  // Analyze whether a Golomb encoding would work.
  const Golomb golomb(range - 1, prefix_size);
  const uint32_t range_golomb = golomb.prefix + 1;
  assert(range_golomb <= histogram_golomb_.size());
  std::fill(&histogram_golomb_[0], &histogram_golomb_[0] + range_golomb, 0);
  *cost = 0.f;
  *size_golomb = 0;
  for (uint32_t i = 0; i < size; ++i) {
    const Golomb golomb_i(mapping_[i], prefix_size);
    if (*size_golomb == 0 ||
        golomb_i.prefix != mapping_golomb_[*size_golomb - 1]) {
      mapping_golomb_[(*size_golomb)++] = golomb_i.prefix;
      assert(*size_golomb - 1 < mapping_golomb_.size());
    }
    histogram_golomb_[*size_golomb - 1] += histogram_[i];
    // Add the bit cost.
    *cost += golomb_i.extra_bits_num * histogram_[i];
  }
  if (*size_golomb == 1) {
    // TODO(vrabaud) Allow for a unique value (probably requires recursion).
    *cost = std::numeric_limits<float>::max();
    return;
  }
  // Quantize the histogram to get the cost of using a dictionary.
  quantizer_golomb_.Quantize(histogram_golomb_.data(), mapping_golomb_.data(),
                             *size_golomb, range_golomb, max_nnz,
                             /*speed=*/9, config_golomb);
  *cost += (*config_golomb)->cost_;
  // Do not use Golomb if there is only one value for now.
  const uint32_t nnz_range_golomb = std::min(max_nnz, range_golomb);
  // Add the prefix_size_m1 and size_m2 cost.
  *cost += 1.f + WP2Log2(nnz_range_golomb - 1);
}

WP2Status SymbolWriter::WriteHeader(uint32_t cluster, uint32_t max_nnz,
                                    uint32_t sym,
                                    const SymbolRecorder& recorder,
                                    WP2_OPT_LABEL, ANSEncBase* const enc,
                                    ANSDictionaries* const dicts) {
  if (symbols_info_.Method(sym) == SymbolsInfo::StorageMethod::kUnused) {
    return WP2_STATUS_OK;
  }
  const uint32_t range = symbols_info_.Range(cluster, sym);
  if (range < 2) {
    assert(range > 0);
    AddTrivial(cluster, sym, 0);
    return WP2_STATUS_OK;
  }
  // We only allow negative values if we have kAuto.
  assert(symbols_info_.Min(cluster, sym) >= 0 ||
         symbols_info_.Method(sym) == SymbolsInfo::StorageMethod::kAuto);

  switch (symbols_info_.Method(sym)) {
    case SymbolsInfo::StorageMethod::kAdaptiveBit:
      WP2_CHECK_STATUS(AddAdaptiveBit(
          cluster, sym, symbols_info_.StartingProbaP0(cluster, sym),
          symbols_info_.StartingProbaP1(cluster, sym)));
      break;
    case SymbolsInfo::StorageMethod::kAdaptiveSym: {
      WP2_CHECK_STATUS(AddAdaptiveSymbol(cluster, sym,
                                         ANSAdaptiveSymbol::Method::kAOM,
                                         kANSAProbaInvalidSpeed));
      break;
    }
    case SymbolsInfo::StorageMethod::kAdaptiveWithAutoSpeed: {
      if (max_nnz <= 1) {
        // Adaptation method/speed don't matter if we only use this symbol once.
        return AddAdaptiveSymbol(cluster, sym,
                                 ANSAdaptiveSymbol::Method::kConstant, 0);
      }
      ANSDebugPrefix prefix(enc, label);
      ANSAdaptiveSymbol s;
      s.SetAdaptationSpeed(ANSAdaptiveSymbol::Method::kAOM);
      s.InitFromUniform(symbols_info_.Range(cluster, sym));
      ANSAdaptiveSymbol::Method method;
      uint32_t speed;
      s.FindBestAdaptationSpeed(recorder.GetRecordedValues(cluster, sym),
                                &method, &speed);
      enc->PutRValue((uint32_t)method,
                     (uint32_t)ANSAdaptiveSymbol::Method::kNum,
                     "adaptation_method");
      if (method == ANSAdaptiveSymbol::Method::kConstant) {
        enc->PutRValue(speed, kNumAdaptationSpeeds, "adaptation_speed");
        speed = kAdaptationSpeeds[speed];
      } else {
        speed = kANSAProbaInvalidSpeed;
      }
      WP2_CHECK_STATUS(AddAdaptiveSymbol(cluster, sym, method, speed));
      return enc->GetStatus();
    }
    case SymbolsInfo::StorageMethod::kAuto: {
      return WriteHeader(cluster, max_nnz, sym,
                         recorder.GetRecordedDict(cluster, sym).Counts().data(),
                         label, enc, dicts);
    }
    case SymbolsInfo::StorageMethod::kUnused:
      // Dealt with above.
      assert(false);
      break;
  }
  return WP2_STATUS_OK;
}

WP2Status SymbolWriter::WriteHeader(uint32_t max_nnz, uint32_t sym,
                                    const SymbolRecorder& recorder,
                                    WP2_OPT_LABEL, ANSEncBase* const enc,
                                    ANSDictionaries* const dicts) {
  for (uint32_t cluster = 0; cluster < symbols_info_.NumClusters(sym);
       ++cluster) {
    WP2_CHECK_STATUS(
        WriteHeader(cluster, max_nnz, sym, recorder, label, enc, dicts));
  }
  return WP2_STATUS_OK;
}

WP2Status SymbolWriter::WriteHeader(uint32_t cluster, uint32_t max_nnz,
                                    uint32_t sym, const uint32_t* const counts,
                                    WP2_OPT_LABEL, ANSEncBase* const enc,
                                    ANSDictionaries* const dicts) {
  // Figure out how many symbols are used and the sizes in bits of the
  // probabilities and difference of consecutive symbols.
  uint32_t size = 0;
  const uint32_t range = symbols_info_.Range(cluster, sym);
  uint32_t counts_total = 0;
  const int32_t min = symbols_info_.Min(cluster, sym);
  const int32_t max = symbols_info_.Max(cluster, sym);
  const bool can_be_negative = (min < 0);
  assert(min == 0 || min == -max);
  int32_t trivial_value = std::numeric_limits<int32_t>::max();
  uint32_t num_nnz_count = 0;
  for (uint32_t k = 0; k < (can_be_negative ? max + 1 : range); ++k) {
    // If the values can be negative, we clump to gether the statistics for x
    // and -x (the symbols we use for now do have symmetric distributions).
    const uint32_t count =
        can_be_negative
            ? (k == 0) ? counts[max] : counts[max + k] + counts[max - k]
            : counts[k];
    if (count == 0) continue;
    if (can_be_negative && k != 0) {
      for (uint32_t c : {-k, k}) {
        if (counts[max + c] != 0) {
          ++num_nnz_count;
          trivial_value = c;
        }
      }
    } else {
      ++num_nnz_count;
      trivial_value = k;
    }
    histogram_[size] = count;
    counts_total += count;
    mapping_[size] = k;
    ++size;
  }

  // Deal with the trivial cases.
  ANSDebugPrefix prefix(enc, label);
  const uint32_t nnz_range = std::min(max_nnz, range);
  assert(size <= nnz_range);
  if (num_nnz_count <= 1) {
    if (size == 0) {
      enc->PutRValue(kSymbolCountZero,
                     nnz_range == 1 ? kSymbolCountLast - 1 : kSymbolCountLast,
                     "scount");
      AddTrivial(cluster, sym, 0);
    } else {
      enc->PutRValue(kSymbolCountOne,
                     nnz_range == 1 ? kSymbolCountLast - 1 : kSymbolCountLast,
                     "scount");
      enc->PutRange(trivial_value, min, max, "symbol");
      AddTrivial(cluster, sym, trivial_value);
    }
    return WP2_STATUS_OK;
  }
  enc->PutRValue(kSymbolCountMoreThanOne, kSymbolCountLast, "scount");

  // Figure out the cost of using ranges.
  // Cost of mapping + cost of storing each symbol as a range + cost of size_m2.
  const float cost_range = StoreMapping(mapping_.data(), size, range,
                                        stats_buffer_.data(), nullptr) +
                           counts_total * WP2Log2(size) +
                           WP2Log2(nnz_range - 1);

  // Compute the best cost of Golomb encoding.
  float cost_golomb, cost_golomb1, cost_golomb2;
  Quantizer::Config* config_golomb = nullptr;
  uint32_t size_golomb1, size_golomb2;
  ComputeGolombCost(range, max_nnz, size, 1, &cost_golomb1, &config_golomb,
                    &size_golomb1);
  ComputeGolombCost(range, max_nnz, size, 2, &cost_golomb2, &config_golomb,
                    &size_golomb2);
  uint32_t golomb_prefix_size;
  uint32_t size_golomb;
  if (cost_golomb1 < cost_golomb2) {
    // Re-initialize the different Golomb parameters.
    ComputeGolombCost(range, max_nnz, size, 1, &cost_golomb1, &config_golomb,
                      &size_golomb1);
    cost_golomb = cost_golomb1;
    size_golomb = size_golomb1;
    golomb_prefix_size = 1;
  } else {
    cost_golomb = cost_golomb2;
    size_golomb = size_golomb2;
    golomb_prefix_size = 2;
  }

  // If we have too many symbols for dictionaries, go with ranges.
  Quantizer::Config* config;
  float cost_dict;
  if (size < ANS_MAX_SYMBOLS) {
    // Quantize the histogram to get the cost of using a dictionary.
    quantizer_.Quantize(histogram_.data(), mapping_.data(), size, range,
                        max_nnz, /*speed=*/9, &config);
    cost_dict = config->cost_;
    cost_dict += WP2Log2(nnz_range);  // Add the size_m1 cost.
  } else {
    cost_dict = std::numeric_limits<float>::max();
  }

  // Find the best storage method.
  Stat::Type type;
  if (cost_golomb < cost_dict && cost_golomb < cost_range) {
    type = Stat::Type::kGolomb;
  } else if (cost_range < cost_dict) {
    type = Stat::Type::kRange;
  } else {
    type = Stat::Type::kDict;
  }

  // TODO(vrabaud) For very small images, always choose range to save bits on
  // this symbol.
  enc->PutRValue(type == Stat::Type::kRange ? 0 :
                 type == Stat::Type::kDict ? 1 : 2,
                 3, "type");
  if (type == Stat::Type::kGolomb) {
    const Golomb golomb(range - 1, golomb_prefix_size);
    const uint32_t range_golomb = golomb.prefix + 1;
    const uint32_t nnz_range_golomb = std::min(max_nnz, range_golomb);
    enc->PutRange(golomb_prefix_size, 1, 2, "prefix_size_m1");
    enc->PutRange(size_golomb, 2, nnz_range_golomb, "size_m2");
    // Store the dictionary.
    WP2_CHECK_STATUS(WriteHistogram(*config_golomb, range_golomb,
                                    max_nnz, enc));
    // Create needed dictionaries.
    WP2_CHECK_STATUS(AddGolomb(
        cluster, sym, histogram_golomb_.data(), config_golomb->h_.counts,
        mapping_golomb_.data(), size_golomb, golomb_prefix_size, dicts));
  } else if (type == Stat::Type::kDict) {
    // The number of symbols is bounded by the number of possible symbols
    // and the size of the image.
    enc->PutRange(size, 2, nnz_range, "size_m2");
    // Store the dictionary.
    WP2_CHECK_STATUS(WriteHistogram(*config, range, max_nnz, enc));
    // Create needed dictionaries.
    WP2_CHECK_STATUS(AddDict(cluster, sym, histogram_.data(),
                             config->h_.counts, mapping_.data(), size,
                             dicts));
  } else {  // kRange
    enc->PutRange(size, 1, nnz_range, "size_m1");
    // Store the mapping.
    StoreMapping(mapping_.data(), size, range, stats_buffer_.data(), enc);
    // Initialize the SymbolWriter.
    AddRange(cluster, sym, mapping_.data(), size, range);
  }
  return WP2_STATUS_OK;
}

void SymbolWriter::GetPotentialUsage(uint32_t cluster, uint32_t sym,
                                     bool is_maybe_used[], size_t size) const {
  const Stat& stat = *GetStats(cluster, sym);
  switch (stat.type) {
    case (Stat::Type::kTrivial):
      // Nothing is used but the one value.
      std::fill(is_maybe_used, is_maybe_used + size, false);
      is_maybe_used[stat.param.trivial_value] = true;
      break;
    case (Stat::Type::kRange):
      // We have no idea of what is used or not.
      std::fill(is_maybe_used, is_maybe_used + size, true);
      break;
    case (Stat::Type::kDict): {
      // Go over the stats to figure out what is used or not.
      if (stat.use_mapping) {
        std::fill(is_maybe_used, is_maybe_used + size, false);
      }
      const Vector_u32& v = stat.extra.dict->Counts();
      if (stat.use_mapping) {
        for (size_t i = 0; i < stat.extra.mapping_size; ++i) {
          if (stat.mappings[i] != Stat::kInvalidMapping) {
            is_maybe_used[i] = (v[stat.mappings[i]] > 0);
          }
        }
      } else {
        for (size_t i = 0; i < stat.extra.dict->MaxSymbol(); ++i) {
          is_maybe_used[i] = (v[i] > 0);
        }
      }
      break;
    }
    case Stat::Type::kGolomb: {
      // Go over the stats to figure out what is used or not.
      if (stat.use_mapping) {
        std::fill(is_maybe_used, is_maybe_used + size, false);
      }
      const Vector_u32& v = stat.extra.dict->Counts();
      const uint32_t prefix_size = stat.param.golomb.prefix_size;
      if (stat.use_mapping) {
        for (size_t i = 0; i < stat.extra.mapping_size; ++i) {
          const auto m = stat.mappings[i];
          if (m != Stat::kInvalidMapping && v[m] > 0) {
            const uint32_t extra_bits_num =
                Golomb::NumExtraBits(i, prefix_size);
            const uint32_t i1 = Golomb::Merge(i, prefix_size, 0);
            assert(i1 <= size);
            const size_t i2 =
                Golomb::Merge(i, prefix_size, (1 << extra_bits_num) - 1);
            std::fill(is_maybe_used + i1,
                      is_maybe_used + std::min(i2 + 1, size), true);
          }
        }
      } else {
        for (size_t i = 0; i < stat.extra.dict->MaxSymbol(); ++i) {
          if (v[i] == 0) continue;
          const uint32_t extra_bits_num = Golomb::NumExtraBits(i, prefix_size);
          const uint32_t i1 = Golomb::Merge(i, prefix_size, 0);
          assert(i1 <= size);
          const size_t i2 =
              Golomb::Merge(i, prefix_size, (1 << extra_bits_num) - 1);
          std::fill(is_maybe_used + i1, is_maybe_used + std::min(i2 + 1, size),
                    true);
        }
      }
      break;
    }
    default:
      assert(false);
  }
}

}  // namespace WP2
