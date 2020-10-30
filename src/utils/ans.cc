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
// Assymetric Numeral System coder
//
// Author: Skal (pascal.massimino@gmail.com)

#include "src/utils/ans.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <memory>
#include <numeric>
#include <utility>

#include "src/dsp/dsp.h"
#include "src/dsp/math.h"
#include "src/utils/data_source.h"
#include "src/utils/utils.h"

namespace WP2 {

ANSSymbolInfo ANSSymbolInfo::Rescale(const ANSSymbolInfo& max,
                                     uint32_t tab_size) const {
  ANSSymbolInfo info;
  const uint32_t max_offset = max.offset + max.freq;
  const uint32_t offset_next = (offset + freq) * tab_size / max_offset;
  info.offset                =  offset         * tab_size / max_offset;
  // Compute the frequency with the next offset to not leave any gap.
  info.freq = offset_next - info.offset;
  info.symbol = symbol;
  assert(offset < max.offset || offset_next == tab_size);
  return info;
}

//------------------------------------------------------------------------------
// ANSBinSymbol

ANSBinSymbol::ANSBinSymbol(uint32_t p0, uint32_t p1) {
  p0_ = p0;
  sum_ = p0 + p1;
  p_ = (sum_ > 0) ? (p0_ << PROBA_BITS) / sum_ : PROBA_MAX / 2;
}

void ANSBinSymbol::UpdateCount(uint32_t bit) {
  assert(sum_ < kMaxSum);
  p0_ += !bit;
  ++sum_;
  p_ = (p0_ << PROBA_BITS) / sum_;   // cache the current proba
}

//------------------------------------------------------------------------------
// update-table (hardcoded for APROBA_MAX_SYMBOL == 16)

// For mixing distributions while preserve norms, see:
//    https://fgiesen.wordpress.com/2015/02/20/mixing-discrete-probability-distributions/    NOLINT

void ANSAdaptiveSymbol::Update(uint32_t sym) {
  const uint32_t mult = adapt_factor_;
  if (method_ == Method::kAOM && adapt_factor_ > kMinAOMAdaptSpeed) {
    adapt_factor_ -= adapt_factor_ >> kAOMAdaptSpeedShift;
  }
  // TODO(skal): SSE or table-free for plain-C
  assert(sym < max_symbol_);
  const uint16_t* const cdf_var = &cdf_var_[APROBA_MAX_SYMBOL - sym - 1];
  ANSUpdateCDF(max_symbol_, cdf_base_.data(), cdf_var, mult, cumul_.data());
  assert(cumul_[0] == 0);
  assert(cumul_[max_symbol_ - 1] < APROBA_MAX);
  assert(max_symbol_ == APROBA_MAX_SYMBOL || cumul_[max_symbol_] == APROBA_MAX);
}

//------------------------------------------------------------------------------

static void GenerateVariableCDF(uint32_t nnz, uint16_t* const table) {
  for (uint32_t i = 0; i < APROBA_MAX_SYMBOL; ++i) table[i] = 0;
  for (uint32_t i = APROBA_MAX_SYMBOL; i < 2 * APROBA_MAX_SYMBOL - 1; ++i) {
    table[i] = APROBA_MAX - nnz;
  }
}

void ANSAdaptiveSymbol::InitCDF() {
  uint32_t nnz = 0;
  cumul_[max_symbol_] = APROBA_MAX;
  for (uint32_t s = 0; s < max_symbol_; ++s) {
    cdf_base_[s] = nnz;
    const bool is_used = (cumul_[s + 1] > cumul_[s]);
    if (is_used) nnz += 1;
  }
  for (uint32_t s = max_symbol_; s < APROBA_MAX_SYMBOL; ++s) cdf_base_[s] = nnz;

  // TODO(skal): could be static tables, indexed by 'nnz'.
  GenerateVariableCDF(nnz, &cdf_var_[0]);

  assert(max_symbol_ >= 1);
  if (max_symbol_ + 1 <= APROBA_MAX_SYMBOL) {
    // It is useful to fill those values for the later search within cumul_.
    std::fill(&cumul_[max_symbol_ + 1], &cumul_[APROBA_MAX_SYMBOL] + 1,
              APROBA_MAX + 1);
  }
}

void ANSAdaptiveSymbol::InitFromUniform(uint32_t max_symbol) {
  assert(max_symbol > 0 && max_symbol <= APROBA_MAX_SYMBOL);
  max_symbol_ = max_symbol;
  for (uint32_t i = 0; i < max_symbol_; ++i) {
    cumul_[i] = (uint16_t)(i * APROBA_MAX / max_symbol_);
  }
  InitCDF();
}

WP2Status ANSAdaptiveSymbol::InitFromCounts(const uint32_t counts[],
                                            uint32_t max_symbol) {
  assert(max_symbol <= APROBA_MAX_SYMBOL);
  max_symbol_ = 0;
  uint32_t pdf[APROBA_MAX_SYMBOL];

  std::copy(counts, counts + max_symbol, pdf);
  WP2_CHECK_STATUS(ANSNormalizeCounts(pdf, max_symbol, APROBA_MAX));
  counts = pdf;

  cumul_[0] = 0;
  for (uint32_t i = 0; i < max_symbol; ++i) {
    if (counts[i] > 0) max_symbol_ = i + 1;
    cumul_[i + 1] = cumul_[i] + (uint16_t)counts[i];
  }
  assert(cumul_[max_symbol_] == APROBA_MAX);
  InitCDF();
  return WP2_STATUS_OK;
}

WP2Status ANSAdaptiveSymbol::InitFromCDF(const uint16_t* const cdf,
                                         uint32_t max_symbol,
                                         uint32_t max_proba) {
  assert(max_symbol > 0 && max_symbol <= APROBA_MAX_SYMBOL);
  assert(max_proba > 0);
  assert(cdf[0] == 0);
  max_symbol_ = max_symbol;
  uint32_t counts[APROBA_MAX_SYMBOL];
  for (uint32_t i = 0; i < max_symbol - 1; ++i) {
    assert(cdf[i + 1] >= cdf[i] && cdf[i + 1] <= max_proba);
    counts[i] = cdf[i + 1] - cdf[i];
  }
  counts[max_symbol - 1] = max_proba - cdf[max_symbol - 1];
  WP2_CHECK_STATUS(InitFromCounts(counts, max_symbol));
  return WP2_STATUS_OK;
}

void ANSAdaptiveSymbol::CopyFrom(const ANSAdaptiveSymbol& other) {
  max_symbol_ = other.max_symbol_;
  adapt_factor_ = other.adapt_factor_;
  std::copy(std::begin(other.cumul_), std::end(other.cumul_),
            std::begin(cumul_));
  std::copy(std::begin(other.cdf_base_), std::end(other.cdf_base_),
            std::begin(cdf_base_));
  std::copy(std::begin(other.cdf_var_), std::end(other.cdf_var_),
            std::begin(cdf_var_));
}

void ANSAdaptiveSymbol::SetAdaptationSpeed(Method method, uint32_t speed) {
  method_ = method;
  if (method_ == Method::kAOM) {
    assert(speed == kANSAProbaInvalidSpeed);
    adapt_factor_ = kMaxAOMAdaptSpeed;
  } else {
    assert(speed != kANSAProbaInvalidSpeed);
    adapt_factor_ = std::min(speed, 0xffffu);
  }
}

void ANSAdaptiveSymbol::Print(bool use_percents) const {
  printf(" -- adaptive [max_symbol:%d] (speed factor=%u): ",
         max_symbol_, adapt_factor_);
  if (use_percents) {
    for (uint32_t i = 0; i < APROBA_MAX_SYMBOL; ++i) {
      printf("[%5.2f%%]", 100.f * GetProba(i));
    }
  } else {
    for (uint32_t i = 0; i <= APROBA_MAX_SYMBOL; ++i) {
      printf("[0x%.4x]", cumul_[i]);
    }
  }
  printf("\n");
}

void ANSAdaptiveSymbol::PrintProbas(uint32_t norm) const {
  if (norm == 0) norm = APROBA_MAX;
  uint32_t pdf[APROBA_MAX_SYMBOL];
  for (uint32_t i = 0; i < APROBA_MAX_SYMBOL; ++i) {
    pdf[i] = cumul_[i + 1] - cumul_[i];
  }
  if (norm != APROBA_MAX) {
    (void)ANSNormalizeCounts(pdf, APROBA_MAX_SYMBOL, APROBA_MAX);
  }
  printf(" { ");
  for (uint32_t i = 0; i < APROBA_MAX_SYMBOL; ++i) {
    if (i > 0) printf(", ");
    printf("%d", pdf[i]);
  }
  printf(" },\n");
}

float ANSAdaptiveSymbol::ComputeDistance(const ANSAdaptiveSymbol& from) const {
  assert(from.max_symbol_ <= max_symbol_);
  float score = 0.f;
  if (max_symbol_ > 0) {
    for (uint32_t i = 1; i <= from.max_symbol_; ++i) {
      const int delta0 = cumul_[i] - cumul_[i - 1];
      const int delta1 = from.cumul_[i] - from.cumul_[i - 1];
      score += abs(delta0 - delta1);
    }
    score /= (APROBA_MAX * from.max_symbol_);
  }
  return score;
}

float ANSAdaptiveSymbol::ComputeNonUniformDistance() const {
  uint32_t score = 0;
  if (max_symbol_ > 0) {
    const float average = 1.f * APROBA_MAX / max_symbol_;
    for (uint32_t i = 1; i <= max_symbol_; ++i) {
      const int delta = cumul_[i] - cumul_[i - 1];
      score += (delta - average) * (delta - average);
    }
  }
  return score / (float)(APROBA_MAX * APROBA_MAX);
}

void ANSAdaptiveSymbol::FindBestAdaptationSpeed(
    const Vector_u8& syms, Method* const method, uint32_t* const speed_idx,
    uint32_t* const header_bits) const {
  if (header_bits != nullptr) *header_bits = 0;
  *speed_idx = 0;
  if (syms.empty()) {
    // Choose AOM by default if there are no elements as it does not require
    // storing a speed.
    *method = Method::kAOM;
    return;
  }

  *method = Method::kConstant;
  float best_cost;
  for (uint32_t idx = 0; idx < kNumAdaptationSpeeds; ++idx) {
    const uint32_t speed = kAdaptationSpeeds[idx];
    const float cost = ScoreSequence(Method::kConstant, speed, syms);
    if (idx == 0 || cost < best_cost) {
      best_cost = cost;
      *speed_idx = idx;
    }
    // Add the cost of including the speed.
    best_cost += WP2Log2Fast(kNumAdaptationSpeeds);
  }
  // Try out AOM's method.
  const float cost_aom =
      ScoreSequence(Method::kAOM, kANSAProbaInvalidSpeed, syms);
  if (cost_aom < best_cost) {
    best_cost = cost_aom;
    *method = Method::kAOM;
    *speed_idx = 0;
  }

  if (header_bits != nullptr) {
    uint32_t counts[APROBA_MAX_SYMBOL] = { 0 };
    for (uint32_t s : syms) ++counts[s];
    float cost = WP2SLog2(syms.size());
    for (uint32_t c : counts) cost -= WP2SLog2(c);
    // fixed probabilities seems better
    if (cost + 2.f < best_cost) {    // +2.f = headroom for inaccuracy
      *header_bits = (uint32_t)(best_cost - cost);
    }
  }
}

float ANSAdaptiveSymbol::ScoreSequence(Method method, uint32_t speed,
                                       const Vector_u8& syms) const {
  ANSAdaptiveSymbol repro;
  repro.InitFromUniform(max_symbol_);
  repro.SetAdaptationSpeed(method, speed);
  float cost = 0.f;
  for (uint32_t s : syms) {
    cost += repro.GetCost(s);
    repro.Update(s);
  }
  return cost;
}

//------------------------------------------------------------------------------
// Enc/Dec tokens.

// Internal token marking. There are 4 variants of token, for binary symbols,
// uniformly drawn values (range or bits) and arbitrary symbols.
typedef uint32_t ANSToken;

enum ANSTokenType {
  ANSTokenTypeBit  = (1u << 29),
  ANSTokenTypeU    = (2u << 29),
  ANSTokenTypeR    = (3u << 29),
  ANSTokenTypeSym  = (4u << 29),
  ANSTokenTypeASym = (5u << 29)
};
inline ANSTokenType ANSTokenGetType(const ANSToken& t) {
  return (ANSTokenType)(t & (7u << 29));
}

// Class generating ANS tokens and extracting info from them.
// A token contains a type stored in the first 3 bits, then a first and a
// second value.
template <ANSTokenType T, uint32_t SecondBits, bool DoCheckSecondBits = true>
class ANSTokenInfo {
 public:
  static constexpr uint32_t kBitsFirst =
      (32 - 3 - SecondBits < IO_BITS) ? (32 - 3 - SecondBits) : IO_BITS;
  static constexpr uint32_t kBitsSecond = SecondBits;
  // We limit values to be at most IO_BITS bits.
  // TODO(skal): check if we can actually use more.
  static_assert(3 + kBitsFirst + kBitsSecond <= 32 && kBitsFirst > 0 &&
                    kBitsFirst <= IO_BITS,
                "wrong token bits");
  static_assert(!DoCheckSecondBits ||
                    (kBitsSecond > 0 && kBitsSecond <= IO_BITS),
                "wrong token bits");
  static uint32_t First(const ANSToken& t) {
    return (t >> kBitsSecond) & ((1u << kBitsFirst) - 1);
  }
  static uint32_t Second(const ANSToken& t) {
    return t & ((1u << kBitsSecond) - 1);
  }
  static ANSToken Init(const uint32_t& f, const uint32_t& s) {
    assert(f < (1u << kBitsFirst));
    assert(s < (1u << kBitsSecond));
    return (ANSToken)(T | (f << kBitsSecond) | s);
  }
};

// We need PROBA_BITS + 1, to be able to store PROBA_MAX.
typedef ANSTokenInfo<ANSTokenTypeBit, PROBA_BITS + 1, false> ANSTokenInfoBit;
// Dictionary and adaptive symbols have the same implementation, just a
// different name.
typedef ANSTokenInfo<ANSTokenTypeSym, APROBA_BITS> ANSTokenInfoSym;
typedef ANSTokenInfo<ANSTokenTypeASym, APROBA_BITS> ANSTokenInfoASym;

// We need enough bits to represent the number of bits with which the symbol
// will be represented. 5 is enough as (1 << 5) bits is enough to represent
// the symbol.
typedef ANSTokenInfo<ANSTokenTypeU, 5> ANSTokenInfoU;
static_assert(kANSMaxUniformBits < (1u << ANSTokenInfoU::kBitsSecond),
              "Too many bits for maximum uniform value.");
static_assert(3 + kANSMaxUniformBits + ANSTokenInfoU::kBitsSecond <= 32,
              "Too many bits for maximum uniform value.");
// kANSMaxRangeBits bits for range.
typedef ANSTokenInfo<ANSTokenTypeR, kANSMaxRangeBits> ANSTokenInfoR;
static_assert(3 + 2 * ANSTokenInfoR::kBitsSecond <= 32u,
              "Too many bits for range.");

// Specializations.
template <>
ANSToken ANSTokenInfoBit::Init(const uint32_t& bit, const uint32_t& proba) {
  assert(proba >= 0 && proba <= PROBA_MAX);
  return (ANSToken)(ANSTokenTypeBit | (!!bit << kBitsSecond) | proba);
}
template <>
ANSToken ANSTokenInfoBit::First(const uint32_t& t) {
  // Optimized computation.
  return !!(t & (1u << kBitsSecond));
}
template <>
ANSToken ANSTokenInfoR::Init(const uint32_t& v, const uint32_t& R) {
  assert(v < (1u << kBitsFirst));
  assert(R > 0 && R <= (1u << kBitsSecond));
  // The range is always > 0 so we can remove 1.
  return (ANSToken)(ANSTokenTypeR | (v << kBitsSecond) | (R - 1));
}
template <>
ANSToken ANSTokenInfoR::Second(const uint32_t& t) {
  return 1 + (t & ((1u << kBitsSecond) - 1));
}

// Adaptive symbol token. We store the offset+freq ANSSymbolInfo.
template <>
ANSToken ANSTokenInfoASym::Init(const uint32_t& s, const uint32_t& f) {
  assert(s < (1u << kBitsFirst));
  assert(f > 0 && f <= (1u << kBitsSecond));
  return (ANSToken)(ANSTokenTypeASym | (s << kBitsSecond) | (f - 1));
}
template <>
ANSToken ANSTokenInfoASym::Second(const uint32_t& f) {
  return 1 + (f & ((1u << kBitsSecond) - 1));
}

// Dictionary symbol: same as adaptive symbol.
template <>
ANSToken ANSTokenInfoSym::Init(const uint32_t& s, const uint32_t& f) {
  assert(s < (1u << kBitsFirst));
  assert(f > 0 && f <= (1u << kBitsSecond));
  return (ANSToken)(ANSTokenTypeSym | (s << kBitsSecond) | (f - 1));
}
template <>
ANSToken ANSTokenInfoSym::Second(const uint32_t& f) {
  return 1 + (f & ((1u << kBitsSecond) - 1));
}
//------------------------------------------------------------------------------
// Decoder

#if defined(WP2_BITTRACE)
constexpr double ANSDec::kBitCountAccuracy;
#endif

uint32_t ANSDec::ReadNextWord(uint32_t s) {
  s <<= IO_BITS;
  const uint8_t* data;
  if (data_source_->TryReadNext(IO_BYTES, &data)) {
#if IO_BITS == 8
    s |= data[0];
#else
    s |= (data[0] << 8) | data[1];
#endif
  } else {
    status_ = WP2_STATUS_BITSTREAM_ERROR;
  }
  return s;
}

uint32_t ANSDec::ReadBitInternal(uint32_t p0) {
  const int q0 = PROBA_MAX - p0;
  if (p0 == 0) {
    return 1;
  } else if (q0 == 0) {
    return 0;
  } else {
    const uint32_t xfrac = state_ & PROBA_MASK;
    const uint32_t bit = (xfrac >= p0);
    if (!bit) {
      state_ = p0 * (state_ >> PROBA_BITS) + xfrac;
    } else {
      state_ = q0 * (state_ >> PROBA_BITS) + xfrac - p0;
    }
    if (state_ < IO_LIMIT_LO) state_ = ReadNextWord(state_);
#if defined(WP2_BITTRACE)
    cur_pos_ += PROBA_BITS - WP2Log2(bit ? q0 : p0);
#endif
    return bit;
  }
}

uint32_t ANSDec::ReadABitInternal(ANSBinSymbol* const stats) {
  const uint32_t bit = ReadBitInternal(stats->Proba());
  return stats->Update(bit);
}

uint32_t ANSDec::ReadSymbolInternal(const ANSSymbolInfo codes[],
                                    uint32_t log2_tab_size,
                                    uint32_t max_index) {
  const uint32_t tab_size = (1u << log2_tab_size);
  const uint32_t res = state_ & (tab_size - 1);
  // Deduce the scaled freq/offset.
  ANSSymbolInfo s;
  // A good start is the scaled index. It is almost always the right value,
  // except for very small frequencies, because of rounding errors.
  const uint32_t max_pos = max_index + codes[max_index].freq;
  uint32_t i = res * max_pos / tab_size;
  i -= codes[i].offset;
  uint32_t offset_next;
  while (true) {
    offset_next = (i + codes[i].freq) * tab_size / max_pos;
    if (res < offset_next) break;
    i += codes[i].freq;
  }
  s.offset = i * tab_size / max_pos;
  s.freq = offset_next - s.offset;
  s.symbol = codes[i].symbol;
  assert(res >= s.offset && res < offset_next);
  state_ = s.freq * (state_ >> log2_tab_size) + (res - s.offset);
  if (state_ < IO_LIMIT_LO) state_ = ReadNextWord(state_);
#if defined(WP2_BITTRACE)
  cur_pos_ += log2_tab_size - WP2Log2(s.freq);
#endif
  return s.symbol;
}

uint32_t ANSDec::ReadSymbolInternal(const ANSAdaptiveSymbol& asym,
                                    uint32_t max_index) {
  const uint32_t res = state_ & ((1u << ANS_LOG_TAB_SIZE) - 1);
  // Deduce the scaled freq/offset.
  // TODO(maryla): optimize so we don't do a linear scan from he first symbol.
  ANSSymbolInfo s;
  s.offset = s.freq = s.symbol = 0;
  const ANSSymbolInfo info_max = asym.GetInfo(max_index);
  for (uint32_t i = 0u; i <= max_index; ++i) {
    s = asym.GetInfo(i).Rescale(info_max);
    if (res < (uint32_t)s.offset + s.freq) break;
  }
  assert(res >= s.offset);
  state_ = s.freq * (state_ >> ANS_LOG_TAB_SIZE) + (res - s.offset);
  if (state_ < IO_LIMIT_LO) state_ = ReadNextWord(state_);
#if defined(WP2_BITTRACE)
  cur_pos_ += ANS_LOG_TAB_SIZE - WP2Log2(s.freq);
#endif
  return s.symbol;
}

uint32_t ANSDec::ReadSymbolInternal(const ANSSymbolInfo codes[],
                                    uint32_t log2_tab_size) {
  const uint32_t res = state_ & ((1u << log2_tab_size) - 1);
  const ANSSymbolInfo s = codes[res];
  state_ = s.freq * (state_ >> log2_tab_size) + s.offset;
  if (state_ < IO_LIMIT_LO) state_ = ReadNextWord(state_);
#if defined(WP2_BITTRACE)
  cur_pos_ += log2_tab_size - WP2Log2(s.freq);
#endif
  return s.symbol;
}

uint32_t ANSDec::ReadSymbolInternal(const ANSAdaptiveSymbol& asym) {
  const uint32_t proba = state_ & (APROBA_MAX - 1);
  const ANSSymbolInfo s = asym.GetSymbol(proba);
  state_ = s.freq * (state_ >> APROBA_BITS) + (proba - s.offset);
  if (state_ < IO_LIMIT_LO) state_ = ReadNextWord(state_);
#if defined(WP2_BITTRACE)
  cur_pos_ += ANSAdaptiveSymbol::GetFreqCost(s.freq);
#endif
  return s.symbol;
}

uint32_t ANSDec::ReadUValueInternal(uint32_t bits) {
  assert(bits <= kANSMaxUniformBits);
  const uint32_t value = state_ & ((1u << bits) - 1);
  state_ >>= bits;
  if (state_ < IO_LIMIT_LO) state_ = ReadNextWord(state_);
#if defined(WP2_BITTRACE)
  cur_pos_ += bits;
#endif
  return value;
}

uint32_t ANSDec::ReadRValueInternal(uint32_t range) {
  // Deal with the power of 2 case.
  if ((range & (range - 1)) == 0) {
    return ReadUValueInternal((uint32_t)WP2Log2Floor(range));
  }
  assert(range <= (1u << ANSTokenInfoR::kBitsSecond));
  // Warning! Unlike other symbols, for RValue we read the NextWord()
  // *before* the extraction operations, therefore making the state_
  // potentially go larger than 32b.
  uint32_t s;
  if (state_ < IO_LIMIT_LO * range) {
    uint64_t state = (uint64_t)state_;
    state = (state << IO_LIMIT_BITS) | ReadNextWord(0);
    s = state % range;
    state_ = (uint32_t)(state / range);
  } else {  // only need 32b operations in this case
    s = state_ % range;
    state_ /= range;
  }
#if defined(WP2_BITTRACE)
  cur_pos_ += WP2Log2(range);
#endif
  return s;
}

void ANSDec::Init(DataSource* const data_source) {
  status_ = WP2_STATUS_OK;
  data_source_ = data_source;
  state_ = 0;
  if (data_source == nullptr) {
    status_ = WP2_STATUS_NULL_PARAMETER;
  } else {
    state_ = ReadNextWord(state_);
    state_ = ReadNextWord(state_);
#if IO_BITS == 8
    state_ = ReadNextWord(state_);
#endif
  }
  if (status_ == WP2_STATUS_OK && state_ < IO_LIMIT_LO) {
    status_ = WP2_STATUS_BITSTREAM_ERROR;
  }
}

//------------------------------------------------------------------------------
// Bit-counters

#if defined(WP2_BITTRACE)
void ANSDec::BitTrace(uint32_t value, const char* const type,
                      const char label[]) {
  if (!is_padding_counted_) {
    const double incr = kANSPaddingCost;
    counters_["ANS_padding"].bits += incr;
    cur_pos_ += incr;
    last_pos_ += incr;
    counters_["ANS_padding"].num_occurrences += 1;
    is_padding_counted_ = true;
  }
  const std::string key = debug_prefix_ + std::string(label);
  counters_[key].bits += cur_pos_ - last_pos_;
  ++counters_[key].num_occurrences;
  ++counters_[key].histo[value];
  if (std::strcmp(type, "bit") == 0) {
    counters_[key].type = LabelStats::Type::Bit;
  } else if (std::strcmp(type, "abit") == 0) {
    counters_[key].type = LabelStats::Type::ABit;
  } else if (std::strcmp(type, "rvalue") == 0) {
    counters_[key].type = LabelStats::Type::RValue;
  } else if (std::strcmp(type, "uvalue") == 0) {
    counters_[key].type = LabelStats::Type::UValue;
  } else if (std::strcmp(type, "asymbol") == 0) {
    counters_[key].type = LabelStats::Type::ASymbol;
  } else if (std::strcmp(type, "symbol") == 0) {
    counters_[key].type = LabelStats::Type::Symbol;
  } else {
    assert(false);
  }

  // Add stats for the custom counter (even if no custom key, not to miss any).
  if (merge_custom_until_pop_) {
    assert(!counters_custom_key_.empty());
    counters_custom_[counters_custom_key_].bits += cur_pos_ - last_pos_;
  } else {
    const std::string counters_custom_key =
        counters_custom_key_.empty() ? key
                                     : (counters_custom_key_ + "/" + label);
    counters_custom_[counters_custom_key].bits += cur_pos_ - last_pos_;
    ++counters_custom_[counters_custom_key].num_occurrences;
  }
  last_pos_ = cur_pos_;
}

void ANSDec::PushBitTracesCustomPrefix(const char* const suffix_of_prefix,
                                       bool merge_until_pop) {
  assert(std::strlen(suffix_of_prefix) > 0);
  if (!counters_custom_key_.empty()) counters_custom_key_ += "/";
  counters_custom_key_ += suffix_of_prefix;

  assert(!merge_custom_until_pop_);
  if (merge_until_pop) {
    merge_custom_until_pop_ = true;
    // Increment by one now, merge next ones.
    ++counters_custom_[counters_custom_key_].num_occurrences;
  }
}

void ANSDec::PopBitTracesCustomPrefix(const char* const suffix_of_prefix) {
  if (counters_custom_key_ == suffix_of_prefix) {
    counters_custom_key_.clear();
  } else {
    const size_t suffix_length = std::strlen(suffix_of_prefix);
    assert(counters_custom_key_.size() > suffix_length);
    const size_t new_length = counters_custom_key_.size() - suffix_length - 1;
    assert(counters_custom_key_[new_length] == '/' &&
           counters_custom_key_.substr(new_length + 1) == suffix_of_prefix);
    counters_custom_key_.resize(new_length);
  }
  merge_custom_until_pop_ = false;
}
#endif

size_t ANSDec::GetNumUsedBytes() const {
  return data_source_->GetNumDiscardedBytes() + data_source_->GetNumReadBytes();
}

size_t ANSDec::GetMinNumUsedBytesDiff(size_t num_used_bytes_before,
                                      size_t num_used_bytes_after) {
  // Precision is +-IO_BYTES.
  return SafeSub(num_used_bytes_after, num_used_bytes_before + IO_BYTES);
}

//------------------------------------------------------------------------------
// Dictionary management

WP2Status ANSDictionary::CopyFrom(const ANSDictionary& dict) {
  max_symbol_ = dict.max_symbol_;
  total_ = dict.total_;
  total_quantized_ = dict.total_quantized_;
  WP2_CHECK_ALLOC_OK(counts_.copy_from(dict.counts_));
  WP2_CHECK_ALLOC_OK(quantized_counts_.copy_from(dict.quantized_counts_));
  WP2_CHECK_ALLOC_OK(infos_.copy_from(dict.infos_));
  log2_tab_size_ = dict.log2_tab_size_;
  return WP2_STATUS_OK;
}

WP2Status ANSDictionary::Init(uint32_t max_symbol) {
  // Note: max_symbol can be > ANS_MAX_SYMBOLS at this point
  // (then, ToCodingTable() will not be called on this dictionary)
  WP2_CHECK_OK(max_symbol > 0, WP2_STATUS_INVALID_PARAMETER);

  log2_tab_size_ = ANS_LOG_TAB_SIZE;

  Clear();
  if (max_symbol > 0) {
    WP2_CHECK_ALLOC_OK(infos_.resize(max_symbol));
    WP2_CHECK_ALLOC_OK(counts_.resize(max_symbol));
    std::fill(counts_.begin(), counts_.end(), 0);
    max_symbol_ = max_symbol;
  }
  return WP2_STATUS_OK;
}

void ANSDictionary::ResetCounts() {
  total_ = 0u;
  total_quantized_ = 0u;
  assert(counts_.size() == max_symbol_);
  std::fill(counts_.begin(), counts_.end(), 0);
  quantized_counts_.clear();
}

void ANSDictionary::Clear() {
  counts_.clear();
  quantized_counts_.clear();
  infos_.clear();
  max_symbol_ = 0;
  ResetCounts();
  // don't reset log2_tab_size_.
}

void ANSDictionary::RecordSymbol(uint32_t symbol, uint32_t count) {
  assert(symbol < max_symbol_);
  if (count == 0) return;
  counts_[symbol] += count;
  total_ += count;
}

float ANSDictionary::SymbolCost(uint32_t symbol) const {
  const uint32_t total = IsQuantized() ? total_quantized_ : total_;
  const float c = IsQuantized() ? quantized_counts_[symbol] : counts_[symbol];
  if (total == 0 || c == 0) return 0.f;
  return WP2Log2(total) - WP2Log2(c);
}

float ANSDictionary::SymbolCost(uint32_t symbol, uint32_t max_symbol) const {
  assert(symbol <= max_symbol);
  const uint32_t total_full = IsQuantized() ? total_quantized_ : total_;
  const Vector_u32& counts = IsQuantized() ? quantized_counts_ : counts_;
  const float c = counts[symbol];
  if (total_full == 0 || c == 0) return 0.f;
  const uint32_t total =
      std::accumulate(&counts[0], &counts[max_symbol + 1], 0u);
  assert(total != 0 && total <= total_full);
  return WP2Log2(total) - WP2Log2(c);
}

float ANSDictionary::Cost() const {
  // TODO(vrabaud) cache the following results if the dictionary has not
  // changed.
  float cost = 0.;
  for (uint32_t i = 0; i < max_symbol_; ++i) {
    // TODO(skal): WP2Log2(total_) is computed each time. Optimize.
    cost += counts_[i] * SymbolCost(i);
  }
  return cost;
}

static void CountsToInfos(const Vector_u32& counts,
                          VectorNoCtor<ANSSymbolInfo>& infos) {
  assert(counts.size() <= infos.size());
  for (uint32_t total = 0, s = 0; s < counts.size(); ++s) {
    const uint32_t freq = counts[s];
    infos[s].freq = freq;
    infos[s].offset = total;
    infos[s].symbol = s;  // not really needed
    total += freq;
  }
}

WP2Status ANSCountsToSpreadTable(uint32_t counts[], uint32_t max_symbol,
                                 uint32_t log2_tab_size,
                                 VectorNoCtor<ANSSymbolInfo>& codes) {
  if (max_symbol == 0) return WP2_STATUS_OK;
  assert(log2_tab_size > 0);
  const uint32_t tab_size = (1u << log2_tab_size);

  WP2_CHECK_STATUS(ANSNormalizeCounts(counts, max_symbol, tab_size));
  WP2_CHECK_ALLOC_OK(codes.resize(tab_size));
  uint32_t pos = 0;
  for (uint16_t s = 0; s < max_symbol; ++s) {
    const uint16_t freq = counts[s];
    for (uint16_t k = 0; k < freq; ++k, ++pos) {
      codes[pos].freq = freq;
      codes[pos].offset = k;
      codes[pos].symbol = s;
    }
  }
  return WP2_STATUS_OK;
}

WP2Status ANSDictionary::ToCodingTable() {
  WP2_CHECK_OK(max_symbol_ <= ANS_MAX_SYMBOLS, WP2_STATUS_INVALID_PARAMETER);
  WP2_CHECK_OK(total_ > 0, WP2_STATUS_INVALID_PARAMETER);

  for (uint32_t i = 0; IsQuantized() && i < max_symbol_; ++i) {
    // If a symbol is used, make sure its quantized count is non-zero.
    assert((counts_[i] == 0) == (quantized_counts_[i] == 0u));
  }
  // total_ has been updated all along with RecordSymbol
  assert(total_ == std::accumulate(counts_.begin(), counts_.end(), 0u));
  if (IsQuantized()) {  // total_quantized_ is updated by SetQuantizedCounts()
    assert(total_quantized_ ==
           std::accumulate(quantized_counts_.begin(),
                           quantized_counts_.end(), 0u));
  }
  const uint32_t tab_size = 1u << log2_tab_size_;
  if (total_quantized_ != tab_size) {
    Vector_u32 counts;
    WP2_CHECK_ALLOC_OK(
        counts.copy_from(IsQuantized() ? quantized_counts_ : counts_));
    WP2_CHECK_STATUS(ANSNormalizeCounts(&counts[0], max_symbol_, tab_size));
    CountsToInfos(counts, infos_);
  } else {
    CountsToInfos(quantized_counts_, infos_);
  }
  return WP2_STATUS_OK;
}

WP2Status ANSDictionary::SetQuantizedCounts(
    const Vector_u32& quantized_counts) {
  WP2_CHECK_OK(quantized_counts.size() == max_symbol_,
               WP2_STATUS_INVALID_PARAMETER);
  WP2_CHECK_ALLOC_OK(quantized_counts_.copy_from(quantized_counts));
  total_quantized_ = std::accumulate(quantized_counts.begin(),
                                     quantized_counts.end(), 0u);
  return WP2_STATUS_OK;
}

//------------------------------------------------------------------------------
// ANSEncDictionaries

WP2Status ANSDictionaries::CopyFrom(const ANSDictionaries& dicts) {
  DeepClear();
  // Copy the dictionaries.
  WP2_CHECK_ALLOC_OK(reserve(dicts.size()));
  for (const ANSDictionary* dict : dicts) {
    std::unique_ptr<ANSDictionary>
        s(new (WP2Allocable::nothrow) ANSDictionary());
    WP2_CHECK_ALLOC_OK(s != nullptr);
    WP2_CHECK_STATUS(s->CopyFrom(*dict));
    WP2_CHECK_ALLOC_OK(push_back(s.get(), /*resize_if_needed=*/false));
    s.release();
  }
  return WP2_STATUS_OK;
}

ANSDictionary* ANSDictionaries::GetEquivalent(
    const ANSDictionaries& original,
    const ANSDictionary* const dict_to_find) const {
  if (size() != original.size()) return nullptr;
  const auto it = std::find(original.begin(), original.end(), dict_to_find);
  if (it == end()) return nullptr;
  ANSDictionary* const equivalent = at(it - original.begin());
  assert(equivalent->SeemsEquivalentTo(*dict_to_find));
  return equivalent;
}

WP2Status ANSDictionaries::Add(uint32_t max_symbol) {
  std::unique_ptr<ANSDictionary>
        s(new (WP2Allocable::nothrow) ANSDictionary());
  WP2_CHECK_ALLOC_OK(s != nullptr);
  WP2_CHECK_STATUS(s->Init(max_symbol));
  WP2_CHECK_ALLOC_OK(push_back(s.get()));
  s.release();
  return WP2_STATUS_OK;
}

WP2Status ANSDictionaries::AppendAndClear(ANSDictionaries* const in) {
  if (in->empty()) return WP2_STATUS_OK;

  WP2_CHECK_ALLOC_OK(reserve(size() + in->size()));
  for (const auto& d : *in) {
    WP2_CHECK_ALLOC_OK(push_back(d, /*resize_if_needed=*/false));
  }
  // Leave ownership to the current ANSDictionaries.
  in->clear();
  return WP2_STATUS_OK;
}

WP2Status ANSDictionaries::ToCodingTable() {
  for (auto& d : *this) WP2_CHECK_STATUS(d->ToCodingTable());
  return WP2_STATUS_OK;
}

//------------------------------------------------------------------------------
// ANSAdaptiveBits

WP2Status ANSAdaptiveBits::CopyFrom(const ANSAdaptiveBits& dicts) {
  clear();
  WP2_CHECK_ALLOC_OK(reserve(dicts.size()));
  for (const ANSBinSymbol& other : dicts) {
    WP2_CHECK_ALLOC_OK(push_back(other, /*resize_if_needed=*/false));
  }
  return WP2_STATUS_OK;
}

//------------------------------------------------------------------------------
// ANSAdaptiveSymbols

bool ANSAdaptiveSymbols::Add(uint32_t max_symbol, bool resize_if_needed) {
  assert(max_symbol <= APROBA_MAX_SYMBOL);
  ANSAdaptiveSymbol s;
  s.InitFromUniform(max_symbol);
  s.SetAdaptationSpeed(ANSAdaptiveSymbol::Method::kAOM);
  return push_back(s, resize_if_needed);
}

WP2Status ANSAdaptiveSymbols::CopyFrom(const ANSAdaptiveSymbols& dicts) {
  clear();
  WP2_CHECK_ALLOC_OK(reserve(dicts.size()));
  for (const ANSAdaptiveSymbol& other : dicts) {
    WP2_CHECK_ALLOC_OK(push_back(other, /*resize_if_needed=*/false));
  }
  return WP2_STATUS_OK;
}

//------------------------------------------------------------------------------
// Encoder

// internal structure for storing tokens in linked pages.
const uint32_t kAnsMaxTokenPerPage = 4096;
struct ANSTokenPage : public WP2Allocable {
  const ANSTokenPage* prev;  // previous -finished- pages
  ANSToken tokens[kAnsMaxTokenPerPage];
  int nb_tokens;
};

// Reset the bitstream buffer.
void ANSEnc::ResetBuffer() {
  buffer_.clear();
  buffer_pos_ = 0;
  state_ = ANS_SIGNATURE;
  assert(state_ >= IO_LIMIT_LO);
}

ANSEnc::ANSEnc() {
  tokens_ = nullptr;
  list_ = nullptr;
  status_ = WP2_STATUS_OK;
  ResetBuffer();
  Reset(true);
}

ANSEnc::~ANSEnc() { WipeOut(); }

ANSEnc::ANSEnc(ANSEnc&& other) noexcept { TrivialMoveCtor(this, &other); }

WP2Status ANSEnc::Clone(const ANSEnc& e) {
  WipeOut();
#if defined(WP2_BITTRACE) || defined(WP2_TRACE) || defined(WP2_ENC_DEC_MATCH)
  debug_prefix_ = e.debug_prefix_;
#endif
  status_ = e.status_;
  WP2_CHECK_STATUS(status_);
  state_ = e.state_;
  buffer_pos_ = e.buffer_pos_;
  WP2_CHECK_ALLOC_OK(buffer_.copy_from(e.buffer_));
  // Copy the tokens.
  const ANSTokenPage* page_e = e.tokens_;
  while (page_e != nullptr) {
    if (!NewPage()) {
      WipeOut();
      return WP2_STATUS_OUT_OF_MEMORY;
    }
    memcpy(tokens_->tokens, page_e->tokens, sizeof(page_e->tokens));
    tokens_->nb_tokens = page_e->nb_tokens;
    page_e = page_e->prev;
  }
  nb_tokens_ = e.nb_tokens_;
  nb_asymbols_ = e.nb_asymbols_;
  nb_uvalues_ = e.nb_uvalues_;
  nb_bits_ = e.nb_bits_;
  return WP2_STATUS_OK;
}

//------------------------------------------------------------------------------
// emit symbol/bit/uniform-value

bool ANSEnc::ReallocBuffer() {
  const uint32_t cur_size = BufferSize();
  const uint32_t new_size = std::max(4096u, 2u * cur_size);
  Vector_u8 new_buffer;
  if (!new_buffer.resize(new_size)) {
    status_ = WP2_STATUS_OUT_OF_MEMORY;
    return false;
  }
  if (cur_size > 0) {
    memcpy(&new_buffer[new_size - cur_size], &buffer_[buffer_pos_], cur_size);
  }
  buffer_pos_ = new_size - cur_size;
  swap(buffer_, new_buffer);
  return true;
}

uint32_t ANSEnc::EmitWord(uint32_t s) {
  if (buffer_pos_ < IO_BYTES) {
    if (!ReallocBuffer()) return (s >> IO_BITS);   // error
  }
  buffer_[--buffer_pos_] = (s >> 0) & 0xff;
#if IO_BITS == 16
  buffer_[--buffer_pos_] = (s >> 8) & 0xff;
#endif
  return (s >> IO_BITS);
}

void ANSEnc::EmitBit(uint32_t bit, uint32_t p0) {
  const uint32_t q0 = PROBA_MAX - p0;
  if (p0 == 0) {
    assert(bit != 0);
  } else if (q0 == 0) {
    assert(bit == 0);
  } else {
    uint32_t s = state_;
    if ((s >> (IO_LIMIT_BITS + IO_BITS - PROBA_BITS)) >= (bit ? q0 : p0)) {
      s = EmitWord(s);
      assert(s < IO_LIMIT_LO);
    }
    if (bit) {
      s = ((s / q0) << PROBA_BITS) + (s % q0) + p0;
    } else {
      s = ((s / p0) << PROBA_BITS) + (s % p0);
    }
    state_ = s;
  }
  assert(state_ >= IO_LIMIT_LO);
}

void ANSEnc::EmitASymbol(uint32_t offset, uint32_t freq) {
  uint32_t s = state_;
  assert(freq > 0);   // coding a symbol that was never seen ?!
  assert(freq <= APROBA_MAX);
  if ((s >> (IO_LIMIT_BITS + IO_BITS - APROBA_BITS)) >= freq) {
    s = EmitWord(s);
  }
  state_ = ((s / freq) << APROBA_BITS) + (s % freq) + offset;
  assert(state_ >= IO_LIMIT_LO);
}

void ANSEnc::EmitUValue(uint32_t value, int bits) {
  // No need to do anything for 0 bits.
  // We cannot go through the code normally as the standard says: if the value
  // of the right operand is negative or is greater or equal to the number of
  // bits in the promoted left operand, the behavior is undefined.
  if (bits == 0) return;
  uint32_t s = state_;
  if (s >> (IO_LIMIT_BITS + IO_BITS - bits) > 0) {
    s = EmitWord(s);
  }
  state_ = (s << bits) + value;
  assert(state_ >= IO_LIMIT_LO);
}

void ANSEnc::EmitRValue(uint32_t value, uint32_t range) {
  // If the range is a power of 2, go with bits: they are faster and more
  // accurate.
  if ((range & (range - 1)) == 0) {
    EmitUValue(value, WP2Log2Floor(range));
    return;
  }
  assert(state_ >= IO_LIMIT_LO);
  // Warning! Unlike other symbols, for RValue we call EmitWord() first,
  // and use 64b for temporary state_.
  uint64_t state = state_ * (uint64_t)range + value;
  if (state >= IO_LIMIT_HI) {
    EmitWord((uint32_t)state);   // flush the lower bits of 'state'
    state >>= IO_LIMIT_BITS;     // rotate-in the higher bits
  }
  assert(state >= IO_LIMIT_LO && state < IO_LIMIT_HI);
  state_ = (uint32_t)state;
}

//------------------------------------------------------------------------------

WP2Status ANSEnc::FinishStream() {
  state_ = EmitWord(state_);
  state_ = EmitWord(state_);
#if IO_BITS == 8
  state_ = EmitWord(state_);
#endif
  return status_;
}

WP2Status ANSEnc::EmitOnlyBits() {
  const ANSTokenPage* page = tokens_;
  while (page != nullptr) {
    int n = page->nb_tokens;
    while (--n >= 0) {
      const ANSToken tok = page->tokens[n];
      const uint32_t p0 = ANSTokenInfoBit::Second(tok);
      const int bit = ANSTokenInfoBit::First(tok);
      assert(ANSTokenGetType(tok) == ANSTokenTypeBit);
      EmitBit(bit, p0);
    }
    page = page->prev;
  }
  return FinishStream();
}

WP2Status ANSEnc::EmitAll() {
  const ANSTokenPage* page = tokens_;
  while (page != nullptr) {
    int n = page->nb_tokens;
    while (--n >= 0) {
      const ANSToken tok = page->tokens[n];
      switch (ANSTokenGetType(tok)) {
        case (ANSTokenTypeBit): {
          const uint32_t p0 = ANSTokenInfoBit::Second(tok);
          const int bit = ANSTokenInfoBit::First(tok);
          EmitBit(bit, p0);
          break;
        }
        case (ANSTokenTypeU):
          EmitUValue(ANSTokenInfoU::First(tok), ANSTokenInfoU::Second(tok));
          break;
        case (ANSTokenTypeR):
          EmitRValue(ANSTokenInfoR::First(tok), ANSTokenInfoR::Second(tok));
          break;
        case (ANSTokenTypeASym): {
          const uint32_t offset = ANSTokenInfoASym::First(tok);
          const uint32_t freq = ANSTokenInfoASym::Second(tok);
          EmitASymbol(offset, freq);
          break;
        }
        case (ANSTokenTypeSym): {
          const uint32_t offset = ANSTokenInfoSym::First(tok);
          const uint32_t freq = ANSTokenInfoSym::Second(tok);
          EmitASymbol(offset, freq);
          break;
        }
      }
    }
    page = page->prev;
  }
  return FinishStream();
}

//------------------------------------------------------------------------------

// Analyze counts[] and renormalize with Squeaky Wheel fix, so that
// the total is rescaled to be equal to 'tab_size' exactly.
// Returns 0 in case of failure. Otherwise, returns the final size.
WP2Status ANSNormalizeCounts(uint32_t counts[], uint32_t max_symbol,
                             uint32_t tab_size) {
  if (max_symbol == 0) return WP2_STATUS_OK;   // corner case, but ok.

  uint64_t total = 0;
  for (uint32_t s = 0; s < max_symbol; ++s) {
    total += counts[s];
  }
  WP2_CHECK_OK(total <= 0xffffffffu, WP2_STATUS_INVALID_PARAMETER);  // overflow
  if (total == tab_size) return WP2_STATUS_OK;   // already ok.
  if (total == 0) {   // force a sane state
    counts[0] = tab_size;
    return WP2_STATUS_OK;
  }

  Vector_u32 keys;
  WP2_CHECK_ALLOC_OK(keys.resize(max_symbol));

  const float norm = 1.f * tab_size / total;
  const float kKeyNorm = (float)(1u << 24) / ANS_MAX_SYMBOLS;
  int miss = tab_size, non_zero = 0;
  for (uint32_t s = 0; s < max_symbol; ++s) {
    if (counts[s] > 0) {
      const float target = norm * counts[s];
      counts[s] = (uint32_t)(target + .5);  // round
      if (counts[s] == 0) counts[s] = 1;
      miss -= counts[s];
      const uint32_t error = (uint32_t)fabsf(kKeyNorm * (target - counts[s]));
      keys[non_zero++] = (error * ANS_MAX_SYMBOLS) + s;
    }
  }
  if (miss == 0) return WP2_STATUS_OK;

  if (miss > 0) {
    if (miss < non_zero) {
      std::nth_element(&keys[0], &keys[non_zero - miss], &keys[non_zero]);
    }
    for (int n = (miss < non_zero) ? non_zero - miss : 0; n < non_zero; ++n) {
      ++counts[keys[n] % ANS_MAX_SYMBOLS];
    }
  } else if (non_zero <= (int)tab_size) {
    const uint32_t cap_count = (1u << 23) - 1;  // to avoid overflow
    // Overflow case. We need to decrease some counts, but need extra care
    // to not make any counts[] go to zero. So we just loop and shave off
    // the largest elements greater than 2 until we're good. It's guaranteed
    // to terminate.
    non_zero = 0;
    for (uint32_t s = 0; s < max_symbol; ++s) {
      if (counts[s] > 1) {
        const uint32_t c = (counts[s] > cap_count) ? cap_count : counts[s];
        keys[non_zero++] = (c * ANS_MAX_SYMBOLS) + s;
      }
    }
    assert(non_zero > 0);
    miss = -miss;
    if (miss < non_zero) {
      std::nth_element(&keys[0], &keys[non_zero - miss], &keys[non_zero]);
    }
    {
      int to_fix = miss;
      while (to_fix > 0) {
        for (int n = (miss < non_zero) ? non_zero - miss : 0;
             n < non_zero && to_fix > 0; ++n) {
          const uint32_t idx = keys[n] % ANS_MAX_SYMBOLS;
          if (counts[idx] > 1) {
            --counts[idx];
            --to_fix;
          }
        }
      }
    }
  } else {
    assert(0);  // impossible: more symbols than final table size expected.
    return WP2_STATUS_INVALID_PARAMETER;
  }
  return WP2_STATUS_OK;
}

//------------------------------------------------------------------------------
// Main call, that will gather stats and code the buffer out.

// Returns the cost of a token that is not a ANSTokenTypeSym.
static inline float GetTokenCost(const ANSToken& tok) {
  switch (ANSTokenGetType(tok)) {
    case ANSTokenTypeBit: {
      const uint32_t p0 = ANSTokenInfoBit::Second(tok);
      if (p0 > 0 && p0 < PROBA_MAX) {
        const bool bit = ANSTokenInfoBit::First(tok);
        return PROBA_BITS - WP2Log2(bit ? PROBA_MAX - p0 : p0);
      }
      break;
    }
    case ANSTokenTypeASym:
      return ANSAdaptiveSymbol::GetFreqCost(ANSTokenInfoASym::Second(tok));
    case ANSTokenTypeSym:
      // This has to be handled by the calling function.
      assert(false);
      break;
    case ANSTokenTypeU:
      return ANSTokenInfoU::Second(tok);
    case ANSTokenTypeR:
      return WP2Log2(ANSTokenInfoR::Second(tok));
    default:
      assert(false);
  }
  return 0.f;
}

float ANSEnc::GetCost() const {
  float cost = 0.f;
  for (const ANSTokenPage* page = tokens_; page != nullptr; page = page->prev) {
    int n = page->nb_tokens;
    while (--n >= 0) {
      const ANSToken& tok = page->tokens[n];
      switch (ANSTokenGetType(tok)) {
        case ANSTokenTypeSym:
          cost += ANSAdaptiveSymbol::GetFreqCost(ANSTokenInfoSym::Second(tok));
          break;
        default:
          cost += GetTokenCost(tok);
          break;
      }
    }
  }
  return cost;
}

float ANSEnc::GetCost(const ANSDictionaries& dicts) const {
  float cost = 0.f;
  for (const ANSTokenPage* page = tokens_; page != nullptr; page = page->prev) {
    int n = page->nb_tokens;
    while (--n >= 0) {
      const ANSToken& tok = page->tokens[n];
      switch (ANSTokenGetType(tok)) {
        case ANSTokenTypeSym:
          // This cost is handled by the dictionaries below. Skip them for now.
          break;
        default:
          cost += GetTokenCost(tok);
          break;
      }
    }
  }
  // Add the symbol costs.
  for (const auto* d : dicts) cost += d->Cost();
  return cost;
}

float ANSEncBase::GetCostFull() const { return GetCost() + kANSPaddingCost; }

float ANSEncBase::GetCostFull(const ANSDictionaries& dicts) const {
  return GetCost(dicts) + kANSPaddingCost;
}

WP2Status ANSEnc::Assemble() {
  ResetBuffer();
  if (nb_uvalues_ == 0 && nb_rvalues_ == 0 && nb_asymbols_ == 0) {
    WP2_CHECK_STATUS(EmitOnlyBits());     // only bits
  } else {
    WP2_CHECK_STATUS(EmitAll());   // bits + uvalues + rvalues + asym
  }
  return status_;
}

void ANSEnc::WipeOut() {
  ResetBuffer();
  Clear();
}

WP2Status ANSEnc::Append(const ANSEnc& enc) {
  return AppendTokens(enc, 0, enc.NumTokens());
}

WP2Status ANSEnc::AppendTokens(const ANSEnc& enc, uint32_t start_token,
                               uint32_t n_tokens) {
  // Figure out the pages in order
  size_t n_pages = 0;
  for (const ANSTokenPage* page = enc.tokens_; page != nullptr;
       page = page->prev) {
    ++n_pages;
  }
  VectorNoCtor<const ANSTokenPage*> pages;
  if (!pages.resize(n_pages)) return WP2_STATUS_OUT_OF_MEMORY;
  size_t idx = n_pages;
  for (const ANSTokenPage* page = enc.tokens_; page != nullptr;
       page = page->prev) {
    pages[--idx] = page;
  }

  uint32_t i_token = 0;
  for (const auto page : pages) {
    if (i_token + page->nb_tokens <= start_token) {
      // Skip this page completely.
      i_token += page->nb_tokens;
      continue;
    } else if (i_token >= start_token + n_tokens) {
      break;
    }
    const uint32_t skip = (start_token > i_token) ? start_token - i_token : 0;
    const uint32_t end_token =
        std::min((uint32_t)page->nb_tokens, start_token + n_tokens - i_token);
    i_token += skip;
    for (uint32_t n = skip; n < end_token; ++n, ++i_token) {
      const ANSToken tok = page->tokens[n];
      switch (ANSTokenGetType(tok)) {
        case ANSTokenTypeBit:
          PutBitInternal(ANSTokenInfoBit::First(tok),
                         ANSTokenInfoBit::Second(tok));
          break;
        case ANSTokenTypeU:
          PutUValueInternal(ANSTokenInfoU::First(tok),
                            ANSTokenInfoU::Second(tok));
          break;
        case ANSTokenTypeR:
          PutRValueInternal(ANSTokenInfoR::First(tok),
                            ANSTokenInfoR::Second(tok));
          break;
        case ANSTokenTypeASym: {
          ANSSymbolInfo info;
          info.offset = ANSTokenInfoASym::First(tok);
          info.freq = ANSTokenInfoASym::Second(tok);
          PutSymbolInternal(info, /*is_adaptive=*/true);
          break;
        }
        case ANSTokenTypeSym: {
          ANSSymbolInfo info;
          info.offset = ANSTokenInfoSym::First(tok);
          info.freq = ANSTokenInfoSym::Second(tok);
          PutSymbolInternal(info, /*is_adaptive=*/false);
          break;
        }
        default:
          assert(0);
          return WP2_STATUS_BAD_WRITE;
      }
    }
  }

  return WP2_STATUS_OK;
}

uint32_t ANSEnc::PutBitInternal(uint32_t bit, uint32_t proba) {
  assert((proba == 0) ? (bit == 1) : (proba != PROBA_MAX || bit == 0));
  EnqueueToken(ANSTokenInfoBit::Init(bit, proba));
  ++nb_bits_;
  return bit;
}

// Append an adaptive binary symbol 'bit', updating 'stats' afterward.
uint32_t ANSEnc::PutABitInternal(uint32_t bit, ANSBinSymbol* const stats) {
  PutBitInternal(bit, stats->Proba());
  return stats->Update(bit);
}

// Append a symbol 'symbol' to the message. Probabilities are optimally
// evaluated at the end of the message, when calling ANSEncAssemble().
// 'symbol' should be in range [0, ANS_MAX_SYMBOLS).
// Returns the symbol.
void ANSEnc::PutSymbolInternal(const ANSSymbolInfo& info, bool is_adaptive) {
  if (is_adaptive) {
    EnqueueToken(ANSTokenInfoASym::Init(info.offset, info.freq));
  } else {
    EnqueueToken(ANSTokenInfoSym::Init(info.offset, info.freq));
  }
  ++nb_asymbols_;
}

// Append a 'value' uniformly distributed in the range [0..1 << bits).
// Current constants impose 'bits' to fit in 8 bits.
// Returns the value, for convenience.
uint32_t ANSEnc::PutUValueInternal(uint32_t value, uint32_t bits) {
  EnqueueToken(ANSTokenInfoU::Init(value, bits));
  assert(bits <= kANSMaxUniformBits);
  assert(value < (1u << bits));
  ++nb_uvalues_;
  return value;
}

// Append a 'value' uniformly distributed in [0..range).
// Current constants impose 'range' to fit in kANSMaxRangeBits bits.
// Returns the value, for convenience.
uint32_t ANSEnc::PutRValueInternal(uint32_t value, uint32_t range) {
  EnqueueToken(ANSTokenInfoR::Init(value, range));
  assert(value < range);
  ++nb_rvalues_;
  return value;
}

//------------------------------------------------------------------------------

void ANSEnc::Reset(bool delete_pages) {
  ANSTokenPage* page = tokens_;
  while (page != nullptr) {
    ANSTokenPage* const next = (ANSTokenPage*)page->prev;
    if (delete_pages) {
      delete page;
    } else {
      page->prev = list_;
      list_ = page;
    }
    page = next;
  }
  tokens_ = nullptr;
  nb_tokens_ = 0;
  nb_asymbols_ = 0;
  nb_uvalues_ = 0;
  nb_rvalues_ = 0;
  nb_bits_ = 0;
}

// Reset the bistream, the message and free old pages.
void ANSEnc::Clear() {
  Reset(true);
  state_ = ANS_SIGNATURE;
  while (list_ != nullptr) {
    ANSTokenPage* const next = (ANSTokenPage*)list_->prev;
    delete list_;
    list_ = next;
  }
}

//------------------------------------------------------------------------------
// Token pages management

bool ANSEnc::NewPage() {
  ANSTokenPage* page = list_;
  if (page != nullptr) {
    list_ = (ANSTokenPage*)page->prev;
  } else {
    page = new (WP2Allocable::nothrow) ANSTokenPage;
    if (page == nullptr) {
      status_ = WP2_STATUS_OUT_OF_MEMORY;
      return false;
    }
  }
  page->prev = tokens_;
  page->nb_tokens = 0;
  tokens_ = page;
  return true;
}

void ANSEnc::EnqueueToken(ANSToken tok) {
  ANSTokenPage* page = tokens_;
  if (page == nullptr || page->nb_tokens == kAnsMaxTokenPerPage) {
    if (!NewPage()) return;
    page = tokens_;
  }
  page->tokens[page->nb_tokens++] = tok;
  ++nb_tokens_;
}

uint32_t ANSEnc::NumTokens() const { return nb_tokens_; }

//------------------------------------------------------------------------------

void swap(ANSEnc& e1, ANSEnc& e2) {
  using std::swap;
#if defined(WP2_BITTRACE) || defined(WP2_TRACE) || defined(WP2_ENC_DEC_MATCH)
  swap(e1.debug_prefix_, e2.debug_prefix_);
#endif
  swap(e1.status_, e2.status_);
  swap(e1.state_, e2.state_);
  swap(e1.buffer_, e2.buffer_);
  swap(e1.buffer_pos_, e2.buffer_pos_);
  swap(e1.tokens_, e2.tokens_);
  swap(e1.list_, e2.list_);
  swap(e1.nb_tokens_, e2.nb_tokens_);
  swap(e1.nb_asymbols_, e2.nb_asymbols_);
  swap(e1.nb_uvalues_, e2.nb_uvalues_);
  swap(e1.nb_bits_, e2.nb_bits_);
}

}   // namespace WP2
