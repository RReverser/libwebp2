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
//  Decoding of ANS symbols.
//
// Author: Vincent Rabaud (vrabaud@google.com)

#include "src/dec/symbols_dec.h"

#include <algorithm>
#include <numeric>

#include "src/dsp/math.h"
#include "src/utils/ans_utils.h"

namespace WP2 {

// 'nnz' is the number of non-zero values in the histogram.
// 'max_count' is the maximum value of a count in the histogram.
static WP2Status ReadHistogram(uint32_t nnz, uint32_t symbol_range,
                               uint32_t max_count, ANSDec* const dec,
                               VectorNoCtor<ANSSymbolInfo>& infos) {
  // Read the mappings first.
  Vector_u16 mapping;
  size_t histo_size = nnz;
  bool is_sparse;
  if (nnz < symbol_range) {
    is_sparse = dec->ReadBool("is_sparse");
    if (is_sparse) {
      WP2_CHECK_STATUS(LoadMapping(dec, nnz, symbol_range, mapping));
    } else {
      histo_size = dec->ReadRange(nnz, symbol_range, "histogram_size");
    }
  } else {
    // If we use the whole range, actually force sparse usage. This seems
    // counter-intuitive but it's actually due to the optimization done for
    // sparse data: 1 is subtracted to all counts[] values.
    is_sparse = true;
  }

  Vector_u16 counts;
  WP2_CHECK_ALLOC_OK(counts.resize(histo_size));

  // 'type' is 0 for raw probabilities, 1 for Huffman-represented
  // probabilities and 2 for Huffman-represented probabilities compressed
  // recursively.
  const uint8_t type = dec->ReadRValue(3, "coding_type");
  const uint32_t max_count_bits =
      std::min(kMaxFreqBits, 1u + (uint32_t)WP2Log2Floor(max_count));
  if (type == 2) {  // HuffmanANS
    // Read the probabilities we depend upon.
    VectorNoCtor<ANSSymbolInfo> infos_sub;
    const size_t size_sub = dec->ReadRange(1, 1 + kMaxFreqBits, "size");
    WP2_CHECK_STATUS(
        ReadHistogram(size_sub, 1 + kMaxFreqBits, histo_size, dec, infos_sub));
    // Deduce the counts.
    uint16_t max_symbol = 0;
    for (const auto& i : infos_sub) {
      max_symbol = std::max(max_symbol, i.symbol);
    }
    // Count the symbols.
    Vector_u32 counts_sub;
    WP2_CHECK_ALLOC_OK(counts_sub.resize(max_symbol + 1));
    std::fill(counts_sub.begin(), counts_sub.end(), 0);
    for (const auto& i : infos_sub) counts_sub[i.symbol] = i.freq;
    // Build the spread table.
    VectorNoCtor<ANSSymbolInfo> codes;
    WP2_CHECK_STATUS(ANSCountsToSpreadTable(
        counts_sub.data(), counts_sub.size(), ANS_LOG_TAB_SIZE, codes));
    // Read the symbols that make up the final stats.
    uint32_t nnz_left = nnz;
    for (size_t i = 0; i < histo_size; ++i) {
      if (nnz_left > 0) {
        counts[i] = dec->ReadSymbol(codes.data(), ANS_LOG_TAB_SIZE, "counts");
      } else {
        counts[i] = 0;  // Only zeros left.
      }
      // If we only have non-zeros left (or in sparse mode), counts were
      // stored - 1.
      if (is_sparse || (i + nnz_left == histo_size)) ++counts[i];
      if (counts[i] != 0) --nnz_left;
    }
  } else {
    // Read the probabilities.
    uint32_t max_freq_bits;
    if (type == 1) {  // Huffman
      max_freq_bits =
          dec->ReadRange(1, 1 + WP2Log2Floor(max_count_bits), "max_freq_bits");
    } else {  // Raw
      max_freq_bits = dec->ReadRange(1, max_count_bits, "max_freq_bits");
    }
    if (is_sparse) {
      ReadVector(dec, (1 << max_freq_bits) - 1, counts);
      for (auto& c : counts) c += 1;
    } else {
      ReadVectorNnz(dec, nnz, (1 << max_freq_bits) - 1, counts);
    }
  }

  // Store everything in the info.
  WP2_CHECK_ALLOC_OK(infos.resize(nnz));
  const bool do_huffman = (type == 1 || type == 2);
  for (size_t ind = 0, k = 0; k < counts.size(); ++k) {
    if (counts[k] > 0) {
      assert(ind < nnz);  // We cannot have more non-zeros than needed.
      infos[ind].symbol = mapping.empty() ? k : mapping[k];
      infos[ind].freq = do_huffman ? (1 << (counts[k] - 1)) : counts[k];
      ++ind;
    }
  }
  return WP2_STATUS_OK;
}

//------------------------------------------------------------------------------

static bool ANSSymbolInfoCmp(const ANSSymbolInfo& i1, const ANSSymbolInfo& i2) {
  return (i1.symbol < i2.symbol);
}

WP2Status SymbolReader::Init(const SymbolsInfo& symbols_info,
                             ANSDec* const dec) {
  dec_ = dec;
  // Get an upper bound on the number of symbols that might use a dictionary.
  uint32_t num_auto = 0u;
  for (uint32_t s = 0; s < symbols_info.Size(); ++s) {
    if (symbols_info.Method(s) == SymbolsInfo::StorageMethod::kAuto) {
      num_auto += symbols_info.NumClusters(s);
    }
  }
  WP2_CHECK_ALLOC_OK(codes_.reserve(num_auto));
  return SymbolIO<StatExtra>::Init(symbols_info);
}

void SymbolReader::AddTrivial(uint32_t cluster, uint32_t sym, int32_t value) {
  Stat* const stat = GetStats(cluster, sym);
  stat->type = Stat::Type::kTrivial;
  stat->param.trivial_value = value;
}

void SymbolReader::AddRange(uint32_t cluster, uint32_t sym,
                            const VectorNoCtor<ANSSymbolInfo>* const infos,
                            uint16_t max_range) {
  Stat* const stat = GetStats(cluster, sym);

  stat->type = Stat::Type::kRange;
  if (infos != nullptr) {
    stat->use_mapping = true;
    for (size_t k = 0; k < infos->size(); ++k) {
      stat->mappings[k] = (*infos)[k].symbol;
    }
    stat->range = infos->size();
  } else {
    stat->use_mapping = false;
    stat->range = max_range;
  }
}

WP2Status SymbolReader::AddDict(uint32_t cluster, uint32_t sym,
                                VectorNoCtor<ANSSymbolInfo>* const infos) {
  Stat* const stat = GetStats(cluster, sym);

  assert(infos->size() > 1);
  stat->type = Stat::Type::kDict;
  stat->use_mapping = true;
  // Make sure the symbols are ordered like in the encoding.
  std::sort(infos->begin(), infos->end(), ANSSymbolInfoCmp);
  // Create the mappings.
  Vector_u32 counts;
  WP2_CHECK_ALLOC_OK(counts.resize(infos->size()));
  for (size_t k = 0; k < counts.size(); ++k) {
    stat->mappings[k] = (*infos)[k].symbol;
    counts[k] = (*infos)[k].freq;
  }
  assert(infos->size() <= ANS_MAX_SYMBOLS);
  // Convert the histogram to a spread table.

  WP2_CHECK_ALLOC_OK(codes_.resize(codes_.size() + 1));
  stat->extra.log2_tab_size = ANS_LOG_TAB_SIZE;
  WP2_CHECK_STATUS(ANSCountsToSpreadTable(
      &counts[0], counts.size(), stat->extra.log2_tab_size,
      codes_.back()));   // will allocate codes_.back()
  stat->extra.codes = codes_.back().data();
  return WP2_STATUS_OK;
}

WP2Status SymbolReader::AddGolomb(uint32_t cluster, uint32_t sym,
                                  VectorNoCtor<WP2::ANSSymbolInfo>* const infos,
                                  uint32_t prefix_size) {
  WP2_CHECK_STATUS(AddDict(cluster, sym, infos));
  SetGolombStat(cluster, sym, prefix_size);

  return WP2_STATUS_OK;
}

int32_t SymbolReader::Read(uint32_t cluster, uint32_t sym, WP2_OPT_LABEL,
                           double* const cost) {
  int32_t value;
  const WP2Status status = ReadInternal(cluster, sym, /*use_max_value=*/false,
                                        /*max_value=*/0, label, &value, cost);
  (void)status;
  assert(status == WP2_STATUS_OK);
  return value;
}

WP2Status SymbolReader::TryRead(uint32_t cluster, uint32_t sym,
                                uint32_t max_value, WP2_OPT_LABEL,
                                int32_t* const value, double* const cost) {
  return ReadInternal(cluster, sym, /*use_max_value=*/true, max_value, label,
                      value, cost);
}

uint32_t SymbolReader::FindSymbol(const Stat& stat, uint32_t max_value) {
  if (stat.mappings[stat.extra.codes[0].symbol] > max_value) {
    // If the first symbol is already bigger than the max_value, stop here.
    return 0;
  }
  uint32_t i_sup = (1u << stat.extra.log2_tab_size) - 1;
  i_sup -= stat.extra.codes[i_sup].offset;
  if (stat.mappings[stat.extra.codes[i_sup].symbol] <= max_value) {
    // Stop if the last interval is actually the right one.
    return i_sup;
  }

  // Use a binary search to find i_inf such that:
  // stat.mappings[symbol at i_inf] <= max_value
  // and stat.mappings[next symbol after the one at i_inf] > max_value
  uint32_t i_inf = 0;
  while (true) {
    const uint32_t i_inf_next = i_inf + stat.extra.codes[i_inf].freq;
    // Stop if we reached a good interval.
    if (stat.mappings[stat.extra.codes[i_inf_next].symbol] > max_value) {
      break;
    }
    i_inf = i_inf_next;
    uint32_t i_mid = (i_inf + i_sup) / 2;
    const ANSSymbolInfo& info_mid = stat.extra.codes[i_mid];
    i_mid -= info_mid.offset;
    if (stat.mappings[info_mid.symbol] > max_value) {
      i_sup = i_mid;
    } else {
      i_inf = i_mid;
    }
  }

  // Make sure the interval is valid: check its symbol.
  const ANSSymbolInfo& info_inf = stat.extra.codes[i_inf];
  assert(stat.mappings[info_inf.symbol] <= max_value);
  // Check if we are the last interval.
  if (i_inf + info_inf.freq == (1u << stat.extra.log2_tab_size)) return i_inf;
  // If not, check the symbol of the next interval.
  const ANSSymbolInfo& info_inf_next = stat.extra.codes[i_inf + info_inf.freq];
  if (stat.mappings[info_inf_next.symbol] <= max_value) assert(false);
  return i_inf;
}

static ANSSymbolInfo GetInfo(const ANSSymbolInfo* const codes, uint32_t value) {
  uint32_t i = 0;
  while (codes[i].symbol < value &&
         i + codes[i].freq < (1 << ANS_LOG_TAB_SIZE)) {
    i += codes[i].freq;
  }
  ANSSymbolInfo info;
  info.offset = i;
  info.freq = codes[i].freq;
  info.symbol = codes[i].symbol;
  return info;
}

// Returns the index of 'mapped' in the 'mappings' array. If 'mapped' is not
// in the array, the index of the first value below it is returned.
static uint16_t Unmap(const uint16_t* const mappings, uint32_t mapping_size,
                      uint16_t mapped) {
  if (mappings[0] > mapped) return 0;
  uint32_t i_sup = mapping_size - 1;
  if (mappings[i_sup] <= mapped) return i_sup;

  // Use a binary search to find i_inf such that:
  // mappings[i_inf] <= max_value
  // and mappings[i_inf + 1] > max_value
  uint32_t i_inf = 0;
  while (true) {
    // Stop if we reached a good interval.
    if (mappings[i_inf + 1] > mapped) {
      break;
    }
    ++i_inf;
    uint32_t i_mid = (i_inf + i_sup) / 2;
    if (mappings[i_mid] > mapped) {
      i_sup = i_mid;
    } else {
      i_inf = i_mid;
    }
  }
  assert(mappings[i_inf] <= mapped);
  assert(i_inf == (mapping_size - 1) || mappings[i_inf + 1] > mapped);
  return i_inf;
}

// Returns the frequency value for symbol "value", scaled to take into account
// that the max possible value is "max_value". Both "value" and "max_value"
// should be raw values (not mapped).
static uint32_t ScaleFreq(const ANSSymbolInfo* codes, uint32_t log2_tab_size,
                          uint32_t value, uint32_t max_value) {
  const uint32_t tab_size = (1u << log2_tab_size);
  const ANSSymbolInfo& info = GetInfo(codes, value);
  const ANSSymbolInfo& info_max = GetInfo(codes, max_value);
  const uint32_t offset =
      info.offset * tab_size / (info_max.offset + info_max.freq);
  const uint32_t offset_next =
      (info.offset + info.freq) * tab_size / (info_max.offset + info_max.freq);
  return offset_next - offset;
}

WP2Status SymbolReader::ReadInternal(uint32_t cluster, uint32_t sym,
                                     bool use_max_value, uint32_t max_value,
                                     WP2_OPT_LABEL, int32_t* const value,
                                     double* const cost) {
  // Consider reading a symbol as a single read occurrence.
  dec_->PushBitTracesCustomPrefix(label, /*merge_until_pop=*/true);

  const Stat& stat = *GetStats(cluster, sym);
  if (use_max_value && stat.use_mapping) {
    // Make sure at least one value in the mapping is inferior to the max_value.
    // We check the minimal/first one.
    // TODO(vrabaud) could be optimized if at encoding time, every time we hit
    //               the minimal stat.mapping values, they are max_value. We
    //               would then remove them from stat.mapping values and let
    //               max_value be used instead.
    WP2_CHECK_OK(stat.mappings[0] <= max_value, WP2_STATUS_BITSTREAM_ERROR);
  }
  ANSDebugPrefix debug_prefix(dec_, label);
  switch (stat.type) {
    case Stat::Type::kTrivial:
      *value = stat.param.trivial_value;
      if (use_max_value) {
        // TODO(vrabaud) could be optimized if at encoding time we only have one
        //               value and hit max_value.
        WP2_CHECK_OK(*value <= (int32_t)max_value, WP2_STATUS_BITSTREAM_ERROR);
      }
      break;
    case Stat::Type::kRange: {
      uint32_t range;
      if (!use_max_value) {
        range = stat.range;
      } else if (stat.use_mapping) {
        // Find the biggest range such that mapping[range] <= max_value.
        range = std::upper_bound(stat.mappings, stat.mappings + stat.range,
                                 max_value) -
                stat.mappings;
        assert(range == stat.range || stat.mappings[range] > max_value);
        assert(stat.mappings[range - 1] <= max_value);
      } else {
        range = std::min(stat.range, (uint16_t)(max_value + 1));
      }
      *value = dec_->ReadRValue(range, "range");
      if (cost != nullptr) *cost += std::log2(range);
      if (stat.use_mapping) *value = stat.mappings[*value];
      break;
    }
    case Stat::Type::kDict: {
      uint32_t raw_value;
      if (use_max_value) {
        const uint32_t max_index = FindSymbol(stat, max_value);
        raw_value = dec_->ReadSymbol(stat.extra.codes, stat.extra.log2_tab_size,
                                     max_index, "dict");
        if (cost != nullptr) {
          const uint32_t mapping_size =
              stat.extra.codes[(1 << stat.extra.log2_tab_size) - 1].symbol + 1;
          const uint32_t scaled_freq =
              ScaleFreq(stat.extra.codes, stat.extra.log2_tab_size, raw_value,
                        Unmap(stat.mappings, mapping_size, max_value));
          *cost -= std::log2(scaled_freq) - stat.extra.log2_tab_size;
        }
      } else {
        raw_value = dec_->ReadSymbol(&stat.extra.codes[0],
                                     stat.extra.log2_tab_size, "dict");
        if (cost != nullptr) {
          const uint16_t freq = GetInfo(stat.extra.codes, raw_value).freq;
          *cost -= std::log2(freq) - stat.extra.log2_tab_size;
        }
      }
      assert(stat.use_mapping);
      *value = stat.mappings[raw_value];
      break;
    }
    case Stat::Type::kGolomb: {
      const uint32_t prefix_size = stat.param.golomb.prefix_size;
      uint32_t prefix;
      bool use_range;
      uint32_t range = 0;
      if (use_max_value) {
        const Golomb golomb_max(max_value, prefix_size);
        const uint32_t i_inf = FindSymbol(stat, golomb_max.prefix);
        const uint32_t raw_prefix = dec_->ReadSymbol(
            &stat.extra.codes[0], stat.extra.log2_tab_size, i_inf, "golomb");
        prefix =
            std::min(golomb_max.prefix, (uint32_t)stat.mappings[raw_prefix]);
        const uint32_t max_prefix =
            std::min(golomb_max.prefix,
                     (uint32_t)stat.mappings[stat.extra.codes[i_inf].symbol]);
        // Use ranges if we are at the last interval.
        use_range = (prefix == max_prefix);
        if (use_range) {
          // Get (the biggest value with a prefix of 'max_prefix') + 1.
          range = std::min(
              max_value + 1,
              Golomb::Merge(max_prefix, prefix_size, 0) +
                  (1 << Golomb::NumExtraBits(max_prefix, prefix_size)));
        }
        if (cost != nullptr) {
          const uint32_t scaled_freq =
              ScaleFreq(stat.extra.codes, stat.extra.log2_tab_size, raw_prefix,
                        golomb_max.prefix);
          *cost -= std::log2(scaled_freq) - stat.extra.log2_tab_size;
        }
      } else {
        const uint32_t raw_prefix = dec_->ReadSymbol(
            &stat.extra.codes[0], stat.extra.log2_tab_size, "golomb");
        prefix = stat.mappings[raw_prefix];
        // Use ranges if we are at the last interval.
        use_range = (prefix + 1 == stat.param.golomb.range);
        if (use_range) range = stat.range;
        if (cost != nullptr) {
          const uint16_t freq = GetInfo(stat.extra.codes, raw_prefix).freq;
          *cost -= std::log2(freq) - stat.extra.log2_tab_size;
        }
      }
      const uint32_t extra_bits = Golomb::NumExtraBits(prefix, prefix_size);
      uint32_t extra_bits_value;
      if (extra_bits == 0) {
        extra_bits_value = 0;
      } else {
        if (use_range) {
          assert(Golomb::Merge(prefix, prefix_size, 0) < range);
          range -= Golomb::Merge(prefix, prefix_size, 0);
          extra_bits_value = dec_->ReadRValue(range, "extra_bits_value");
          if (cost != nullptr) {
            *cost += std::log2(range);
          }
        } else {
          extra_bits_value = dec_->ReadUValue(extra_bits, "extra_bits_value");
          if (cost != nullptr) {
            *cost += extra_bits;
          }
        }
      }
      *value = Golomb::Merge(prefix, prefix_size, extra_bits_value);
      break;
    }
    case Stat::Type::kAdaptiveBit: {
      if (use_max_value && max_value == 0) {
        *value = 0;
      } else {
        *value =
            dec_->ReadBit(a_bits_[stat.param.a_bit_index].Proba(), "a_bit");
        if (cost != nullptr) {
          *cost += a_bits_[stat.param.a_bit_index].GetCost(*value);
        }
      }
      a_bits_[stat.param.a_bit_index].Update(*value);
      break;
    }
    case Stat::Type::kAdaptiveSymbol: {
      ANSAdaptiveSymbol& asym = a_symbols_[stat.param.a_symbol_index];
      if (use_max_value) {
        *value = dec_->ReadSymbol(asym, max_value, "a_symbol");
        if (cost != nullptr) {
          *cost += asym.GetCost(*value, max_value);
        }
      } else {
        *value = dec_->ReadSymbol(asym, "a_symbol");
        if (cost != nullptr) {
          *cost += asym.GetCost(*value);
        }
      }
      asym.Update(*value);
      break;
    }
    default:
      // Did you forget to call ReadHeader for this symbol?
      assert(false);
      *value = 0;
      dec_->PopBitTracesCustomPrefix(label);
      return WP2_STATUS_BITSTREAM_ERROR;
  }
  if (*value != 0 && stat.type != Stat::Type::kTrivial &&
      symbols_info_.Min(cluster, sym) < 0) {
    // TODO(vrabaud) implement other methods to store negative values.
    if (dec_->ReadBool("is_negative")) *value = -(*value);
    if (cost != nullptr) *cost += 1.;
  }
  if (use_max_value) assert((uint32_t)std::abs(*value) <= max_value);
  dec_->PopBitTracesCustomPrefix(label);
  return WP2_STATUS_OK;
}

WP2Status SymbolReader::ReadHeader(uint32_t max_nnz, uint32_t sym,
                                   WP2_OPT_LABEL) {
  for (uint32_t cluster = 0; cluster < symbols_info_.NumClusters(sym);
       ++cluster) {
    WP2_CHECK_STATUS(ReadHeader(cluster, max_nnz, sym, label));
  }
  return WP2_STATUS_OK;
}

WP2Status SymbolReader::ReadHeader(uint32_t cluster, uint32_t max_nnz,
                                   uint32_t sym, WP2_OPT_LABEL) {
  if (symbols_info_.Method(sym) == SymbolsInfo::StorageMethod::kUnused) {
    return WP2_STATUS_OK;
  }
  const uint32_t range = symbols_info_.Range(cluster, sym);
  if (range < 2) {
    assert(range > 0);
    AddTrivial(cluster, sym, 0);
    return WP2_STATUS_OK;
  }

  ANSDebugPrefix prefix(dec_, label);

  switch (symbols_info_.Method(sym)) {
    case SymbolsInfo::StorageMethod::kAuto:
      // This main case is handled after.
      break;
    case SymbolsInfo::StorageMethod::kAdaptiveBit:
      WP2_CHECK_STATUS(AddAdaptiveBit(
          cluster, sym, symbols_info_.StartingProbaP0(cluster, sym),
          symbols_info_.StartingProbaP1(cluster, sym)));
      return WP2_STATUS_OK;
    case SymbolsInfo::StorageMethod::kAdaptiveSym: {
      WP2_CHECK_STATUS(AddAdaptiveSymbol(cluster, sym,
                                         ANSAdaptiveSymbol::Method::kAOM,
                                         kANSAProbaInvalidSpeed));
      return WP2_STATUS_OK;
    }
    case SymbolsInfo::StorageMethod::kAdaptiveWithAutoSpeed: {
      if (max_nnz <= 1) {
        // Adaptation method/speed don't matter if we only use this symbol once.
        return AddAdaptiveSymbol(cluster, sym,
                                 ANSAdaptiveSymbol::Method::kConstant, 0);
      }
      const auto method = (ANSAdaptiveSymbol::Method)dec_->ReadRValue(
          (uint32_t)ANSAdaptiveSymbol::Method::kNum, "adaptation_method");
      uint32_t speed;
      if (method == ANSAdaptiveSymbol::Method::kConstant) {
        speed = kAdaptationSpeeds[dec_->ReadRValue(kNumAdaptationSpeeds,
                                                   "adaptation_speed")];
      } else {
        speed = kANSAProbaInvalidSpeed;
      }
      WP2_CHECK_STATUS(AddAdaptiveSymbol(cluster, sym, method, speed));
      return WP2_STATUS_OK;
    }
    case SymbolsInfo::StorageMethod::kUnused:
      break;
  }

  const uint32_t nnz_range = std::min(max_nnz, range);
  const SymbolCount c = (SymbolCount)dec_->ReadRValue(
      nnz_range == 1 ? kSymbolCountLast - 1 : kSymbolCountLast, "scount");

  // Deal with trivial cases.
  if (c == kSymbolCountZero) {
    AddTrivial(cluster, sym, 0);
    return WP2_STATUS_OK;
  }
  if (c == kSymbolCountOne) {
    const int32_t value =
        dec_->ReadRange(symbols_info_.Min(cluster, sym),
                        symbols_info_.Max(cluster, sym), "symbol");
    AddTrivial(cluster, sym, value);
    return WP2_STATUS_OK;
  }

  const uint32_t type = dec_->ReadRValue(3, "type");
  if (type == 2) {         // kGolomb
    const uint32_t golomb_prefix_size = dec_->ReadRange(1, 2, "prefix_size_m1");
    const Golomb golomb(range - 1, golomb_prefix_size);
    const uint32_t range_golomb = golomb.prefix + 1;
    const uint32_t nnz =
        dec_->ReadRange(2, std::min(max_nnz, range_golomb), "size_m2");
    WP2_CHECK_STATUS(
        ReadHistogram(nnz, range_golomb, max_nnz, dec_, infos_));
    WP2_CHECK_STATUS(AddGolomb(cluster, sym, &infos_, golomb_prefix_size));
  } else if (type == 1) {  // kDict
    const uint32_t nnz = dec_->ReadRange(2, nnz_range, "size_m2");
    // TODO(vrabaud) Solve this differently.
    WP2_CHECK_OK(nnz < ANS_MAX_SYMBOLS, WP2_STATUS_BITSTREAM_ERROR);
    WP2_CHECK_STATUS(ReadHistogram(nnz, range, max_nnz, dec_, infos_));
    WP2_CHECK_STATUS(AddDict(cluster, sym, &infos_));
  } else {                 // kRange
    const uint32_t nnz = dec_->ReadRange(1, nnz_range, "size_m1");
    WP2_CHECK_STATUS(LoadMapping(dec_, nnz, range, mapping_));
    // Store everything in the info.
    WP2_CHECK_ALLOC_OK(infos_.resize(nnz));
    for (size_t k = 0; k < infos_.size(); ++k) {
      infos_[k].symbol = mapping_[k];
      infos_[k].freq = 1;
    }
    AddRange(cluster, sym, &infos_, range - nnz);
  }
  return WP2_STATUS_OK;
}

void SymbolReader::GetPotentialUsage(uint32_t cluster, uint32_t sym,
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
    case (Stat::Type::kDict):
      // Go over the stats to figure out what is used or not.
      if (stat.use_mapping) {
        std::fill(is_maybe_used, is_maybe_used + size, false);
      }
      for (uint32_t k = 0; k < (1u << stat.extra.log2_tab_size);
           k += stat.extra.codes[k].freq) {
        is_maybe_used[stat.use_mapping
                          ? stat.mappings[stat.extra.codes[k].symbol]
                          : stat.extra.codes[k].symbol] = true;
      }
      break;
    case Stat::Type::kGolomb: {
      // Go over the stats to figure out what is used or not.
      if (stat.use_mapping) {
        std::fill(is_maybe_used, is_maybe_used + size, false);
      }
      const uint32_t prefix_size = stat.param.golomb.prefix_size;
      for (uint32_t k = 0; k < (1u << stat.extra.log2_tab_size);
           k += stat.extra.codes[k].freq) {
        const uint32_t m = stat.use_mapping
                               ? stat.mappings[stat.extra.codes[k].symbol]
                               : stat.extra.codes[k].symbol;
        const uint32_t extra_bits_num = Golomb::NumExtraBits(m, prefix_size);
        const uint32_t m1 = Golomb::Merge(m, prefix_size, 0);
        const size_t
            m2 = Golomb::Merge(m, prefix_size, (1 << extra_bits_num) - 1);
        std::fill(is_maybe_used + m1, is_maybe_used + std::min(m2 + 1, size),
                  true);
      }
      break;
    }
    default:
      assert(false);
      break;
  }
}

}  // namespace WP2
