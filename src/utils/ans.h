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

#ifndef WP2_UTILS_ANS_H_
#define WP2_UTILS_ANS_H_

#include <algorithm>
#include <array>
#include <cassert>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <map>
#include <memory>
#if defined(WP2_BITTRACE) || defined(WP2_TRACE) || defined(WP2_ENC_DEC_MATCH)
#include <string>

#include "src/wp2/debug.h"
#endif

#ifdef HAVE_CONFIG_H
#include "src/wp2/config.h"
#endif

#include "src/utils/data_source.h"
#include "src/utils/utils.h"
#include "src/utils/vector.h"
#include "src/dsp/dsp.h"
#include "src/dsp/math.h"

#define WP2_OPT_LABEL const char label[]

#if defined(WP2_BITTRACE)
// For each label, store the amount of bits, and the number of occurrences.
// (can be useful whether WP2_BITTRACE is defined or not)
typedef std::map<const std::string, WP2::LabelStats> WP2BitCounts;
#endif

// When uncommenting the following, the 'label' string used when inputting
// anything in the ANS encoder is also inserted, as a hash or as the whole
// string if WP2_ENC_DEC_DEEP_MATCH is defined. When decoding the ANS
// stream, it is asserted that the same string is used to make sure that the
// same information is read in the same order.
// #define WP2_ENC_DEC_MATCH

namespace WP2 {

#define ANS_LOG_TAB_SIZE 14
#define ANS_TAB_SIZE (1u << ANS_LOG_TAB_SIZE)
#define ANS_TAB_MASK (ANS_TAB_SIZE - 1)
#define ANS_MAX_SYMBOLS_BITS 10
#define ANS_MAX_SYMBOLS (1 << ANS_MAX_SYMBOLS_BITS)

#define PROBA_BITS 16
// bits have a probability within [0, PROBA_MAX]
#define PROBA_MAX (1u << PROBA_BITS)
#define PROBA_MASK (PROBA_MAX - 1)

#define IO_BITS 16    // 8 or 16 only
#define IO_BYTES 2    // IO_BITS = 8 << IO_BYTES
#define IO_MASK ((1u << IO_BITS) - 1u)

#define IO_LIMIT_BITS 16
#define IO_LIMIT_LO (1ull << IO_LIMIT_BITS)
#define IO_LIMIT_LO_MASK (IO_LIMIT_LO - 1)
#define IO_LIMIT_HI (1ull << (IO_LIMIT_BITS + IO_BITS))
#define ANS_SIGNATURE (0xf3 * IO_LIMIT_LO)    // Initial state, used as CRC.

// Max number (inclusive) of bits for the range used during range-value coding.
static constexpr uint32_t kANSMaxRangeBits = 14;
static constexpr uint32_t kANSMaxRange = (1u << kANSMaxRangeBits) - 1;
// Max number (inclusive) of bits to use to describe a uniform value.
static constexpr uint32_t kANSMaxUniformBits = IO_BITS;

// Adaptive symbol with a small dictionary
// APROBA_BITS must be a little less precision than PROBA_BITS in order to fit
// the range in the ANSTokenInfoASym as 14+14 bits.
#define APROBA_BITS 14
#define APROBA_MAX (1 << APROBA_BITS)
#define APROBA_ADAPT_SPEED 4096   // must be strictly less than 32768
static constexpr uint32_t kANSAProbaInvalidSpeed = (1 << PROBA_BITS);
#define APROBA_MAX_SYMBOL 16u     // for now, hardcoded for 4,8 or 16.

#if IO_BITS == 8
static constexpr uint32_t kANSPaddingCost = 3 * 8;
#else
static_assert(IO_BITS == 16, "IO_BITS should be == 16");
static constexpr uint32_t kANSPaddingCost = 4 * 8;
#endif

// TODO(vrabaud): make sure the ANS read/write functions are optimized/inlined.
//                A solution could be to get Dec and Enc to be friends and have
//                the Token logic be part of one of them.

//------------------------------------------------------------------------------
// Decoding

struct ANSSymbolInfo {  // Symbol information.
  uint16_t offset;      // offset to the first ANSSymbolInfo of that symbol
  uint16_t freq;        // frequency in [0, PROBA_MAX]
  uint16_t symbol;
  ANSSymbolInfo Rescale(const ANSSymbolInfo& max_symbol,
                        uint32_t tab_size = ANS_TAB_SIZE) const;
};
typedef VectorNoCtor<ANSSymbolInfo> ANSCodes;

// Struct for recording binary event and updating probability.
// TODO(maryla): merge this into ANSAdaptiveSymbol to make them simpler to use?
class ANSBinSymbol {
 public:
  explicit ANSBinSymbol(uint32_t p0 = 1, uint32_t p1 = 1);
  ANSBinSymbol(ANSBinSymbol&&) = default;
  ANSBinSymbol(const ANSBinSymbol& other) = default;
  ANSBinSymbol& operator=(const ANSBinSymbol&) = default;
  // update observation and return 'bit'
  inline uint32_t Update(uint32_t bit) {
    if (sum_ < kMaxSum) UpdateCount(bit);
    return bit;
  }
  // get the cached probability for symbol '1' out of the current state
  constexpr uint32_t Proba() const { return p_; }

  // cost of storing the given 'bit'.
  float GetCost(uint32_t bit) const {
    const uint32_t p = Proba();
    if (p == 0 || p == PROBA_MAX) return 0;
    return -WP2Log2(bit ? PROBA_MAX - p : p) + PROBA_BITS;
  }

 private:
  // we are dealing with quite stationary sources, so updating the proba
  // past kMaxSum events is usually irrelevant.
  static constexpr uint32_t kMaxSum = 256u;
  uint16_t p0_, sum_;
  uint32_t p_;

  void UpdateCount(uint32_t bit);
};

// Generate a 'flat' spread table from a set of frequencies. The counts[]
// distributions is re-normalized such that its sum is equal to 1<<log2_tab_size
// and codes will be resized to '1<<log2_tab_size' if needed.
// Returns ok if there are no ill-defined counts.
WP2Status ANSCountsToSpreadTable(uint32_t counts[], uint32_t max_symbol,
                                 uint32_t log2_tab_size,
                                 VectorNoCtor<ANSSymbolInfo>& codes);

// Analyze counts[] and renormalize with the total is rescaled to be equal
// to 'tab_size' exactly.
WP2Status ANSNormalizeCounts(uint32_t counts[], uint32_t max_symbol,
                             uint32_t tab_size);

// Stores and adapts small APROBA_MAX_SYMBOL-dictionary
class ANSAdaptiveSymbol {
 public:
  enum class Method {
    // Adaptation is constant with a speed in kAdaptationSpeeds.
    kConstant,
    // Adaptation is exactly the AOM one, with a rate depending on how many
    // symbols got already written.
    kAOM,
    kNum
  };

  ANSAdaptiveSymbol() : method_(Method::kNum) {}

  // Initializes with a uniform distribution.
  void InitFromUniform(uint32_t max_symbol);
  // Initializes from an un-normalized pdf. Norm is APROBA_MAX.
  // Returns false if the distribution can't be normalized.
  WP2Status InitFromCounts(const uint32_t counts[], uint32_t max_symbol);
  // Initializes from a cdf, normalization to APROBA_MAX is performed if the CDF
  // is normalized to max_proba.
  WP2Status InitFromCDF(const uint16_t* const cdf, uint32_t max_symbol,
                        uint32_t max_proba = APROBA_MAX);
  void CopyFrom(const ANSAdaptiveSymbol& other);

  // Adaptation speed goes from 0x3fffu (fast adaptation) to 0 (no adaptation).
  // Special case is speed 0xffffu (instant adaptation)
  void SetAdaptationSpeed(Method method,
                          uint32_t speed = kANSAProbaInvalidSpeed);

  // proba is in [0, APROBA_MAX) range
  ANSSymbolInfo GetSymbol(uint32_t proba) const {
    uint32_t s = 0;
    // TODO(yguyon): Move SSE/C code to dsp/
#if defined(WP2_USE_SSE)
    const __m128i A0 = _mm_loadu_si128((const __m128i*)&cumul_[0]);
    const __m128i A1 = _mm_loadu_si128((const __m128i*)&cumul_[8]);
    const __m128i B = _mm_set1_epi16(proba + 1);
    const __m128i C0 = _mm_cmplt_epi16(A0, B);  // we'd need mm_cmple_epi16 !!
    const __m128i C1 = _mm_cmplt_epi16(A1, B);
    const __m128i D = _mm_packs_epi16(C0, C1);
    const uint32_t bits = _mm_movemask_epi8(D);
    s = WP2Log2Floor(bits);
#else
#if 0
    for (; s < max_symbol_ - 1; ++s) {
      if (proba < cumul_[s + 1]) break;
    }
#else  // dichotomy
#if (APROBA_MAX_SYMBOL > 8)
    if (proba >= cumul_[8]) s += 8;
#endif
#if (APROBA_MAX_SYMBOL > 4)
    if (proba >= cumul_[s + 4]) s += 4;
#endif
    if (proba >= cumul_[s + 2]) s += 2;
    if (proba >= cumul_[s + 1]) s += 1;
#endif
#endif  // !WP2_USE_SSE
    return GetInfo(s);
  }
  // Convenience representation of ranges as ANSSymbolInfo.
  // 'sym' must be in [0, max_symbol_) range.
  inline ANSSymbolInfo GetInfo(uint32_t sym) const {
    assert(sym < max_symbol_);
    ANSSymbolInfo info;
    info.symbol = (uint16_t)sym;
    info.offset = cumul_[sym];
    info.freq = cumul_[sym + 1] - info.offset;
    return info;
  }

  // Updates the cumulative distribution after coding symbol 'sym',
  // preserving the norm equal to APROBA_MAX.
  void Update(uint32_t sym);

  // Returns proba = (cumul[sym + 1] - cumul[sym]) / APROBA_MAX
  float GetProba(uint32_t sym) const {
    assert(sym < max_symbol_);
    const uint32_t w = cumul_[sym + 1] - cumul_[sym];
    return (1.f / APROBA_MAX) * w;
  }
  float GetCost(uint32_t sym) const {
    return kCostTable[cumul_[sym + 1] - cumul_[sym]];
  }
  float GetCost(uint32_t sym, uint32_t max_symbol) const {
    assert(sym <= max_symbol);
    return kCostTable[cumul_[sym + 1] - cumul_[sym]]
         - kCostTable[cumul_[max_symbol + 1] - cumul_[0]];
  }

  // debug
  void Print(bool use_percents = true) const;
  void PrintProbas(uint32_t norm = 0) const;

  // Computes the 'distance' between two cdf. Can be used to
  // monitor convergence to final stationary distribution.
  float ComputeDistance(const ANSAdaptiveSymbol& from) const;
  // Distance (in [0..1]) based on variance.
  float ComputeNonUniformDistance() const;

  // Will return the optimal speed index to use to reach the final distribution.
  // The index refers to kAdaptationSpeeds[].
  // If header_bits is not nullptr, and in case a fixed distribution is better
  // than an adaptive one (ie.: adaptation speed looks too fast), will return
  // the number of header-bits that can be spent to code the fixed distribution.
  // If fixed-proba is not advantageous, 0 is returned as header_bits.
  void FindBestAdaptationSpeed(const Vector_u8& syms, Method* const method,
                               uint32_t* const speed,
                               uint32_t* const header_bits = nullptr) const;
  // Evaluates the cost of the sequence of symbols for the given
  // adaptation speed.
  float ScoreSequence(Method method, uint32_t speed,
                      const Vector_u8& syms) const;

  const uint16_t* GetCumul() const { return cumul_.data(); }
  uint32_t NumSymbols() const { return max_symbol_; }

  static inline float GetFreqCost(uint32_t freq) { return kCostTable[freq]; }

 protected:
  Method method_;
  uint32_t max_symbol_;  // At most APROBA_MAX_SYMBOL.
  uint32_t adapt_factor_ = APROBA_ADAPT_SPEED;
  // cumulative frequency. cumul[0] is always 0, cumul[last] is always
  // APROBA_MAX.
  std::array<uint16_t, APROBA_MAX_SYMBOL + 1> cumul_;

  void InitCDF();  // fills cdf_base_[] and cdf_var_[] according to cumul_[]
  std::array<uint16_t, APROBA_MAX_SYMBOL> cdf_base_;
  std::array<uint16_t, APROBA_MAX_SYMBOL * 2 - 1> cdf_var_;

  // precalc for -log2(GetProba())
  //  TODO(skal): Costs should be normalized to int16_t, not float.
  static const float kCostTable[APROBA_MAX + 1];
};

#if defined(WP2_ENC_DEC_MATCH)
#if defined(WP2_ENC_DEC_DEEP_MATCH)
template <typename... Ts>
inline std::string ANSString(const std::string& debug_prefix,
                             const char label[], Ts... extra) {
  std::string str = debug_prefix;
  str += label;
  const uint16_t extra_arr[] = {(const uint16_t)extra...};
  for (uint16_t e : extra_arr) str += ", " + std::to_string(e);
  return str;
}
#else
// djb2 string hash
template <typename... Ts>
inline uint16_t ANSStringHash(const std::string& debug_prefix,
                              const char label[], Ts... extra) {
  uint16_t hash = 5381;
  for (char c : debug_prefix) hash = ((hash << 5) + hash) + c;
  for (int i = 0; label[i] != 0; ++i) hash = ((hash << 5) + hash) + label[i];
  const uint16_t extra_arr[] = {(const uint16_t)extra...};
  for (uint16_t e : extra_arr) hash = ((hash << 5) + hash) + e;
  return hash;
}
#endif  // defined(WP2_ENC_DEC_DEEP_MATCH)
#endif  // defined(WP2_ENC_DEC_MATCH)

// Class for holding the decoding state.
class ANSDec {
 public:
  explicit ANSDec(DataSource* const data_source) { Init(data_source); }
  ANSDec(const ANSDec&) = delete;

  // initializes a new ANSDec object.
  void Init(DataSource* const data_source);

  // Decodes a symbol, according to the spread table 'codes'.
  uint32_t ReadSymbol(const ANSSymbolInfo codes[], uint32_t log2_tab_size,
                      WP2_OPT_LABEL) {
    const uint32_t symbol = ReadSymbolInternal(codes, log2_tab_size);
    Trace("%s: symbol=%u", label, symbol);
    BitTrace(symbol, "symbol", label);
    return symbol;
  }

  // Decodes a symbol, according to the info.
  // 'max_index' is the index in 'codes' for which we know for sure the read
  // value is not strictly superior.
  uint32_t ReadSymbol(const ANSSymbolInfo codes[], uint32_t log2_tab_size,
                      uint32_t max_index, WP2_OPT_LABEL) {
    const uint32_t symbol = ReadSymbolInternal(codes, log2_tab_size, max_index);
    Trace("%s: symbol=%u", label, symbol);
    BitTrace(symbol, "symbol", label);
    return symbol;
  }

  // Decodes a symbol from small adaptive dictionary.
  // 'asym' probability is adapted.
  uint32_t ReadASymbol(ANSAdaptiveSymbol* const asym, WP2_OPT_LABEL) {
    const uint32_t symbol = ReadSymbol(*asym, label);
    asym->Update(symbol);
    return symbol;
  }
  // Decodes a symbol from small adaptive dictionary.
  // 'asym' probability is *NOT* adapted.
  uint32_t ReadSymbol(const ANSAdaptiveSymbol& asym, WP2_OPT_LABEL) {
    const uint32_t symbol = ReadSymbolInternal(asym);
    Trace("%s: symbol=%u", label, symbol);
    BitTrace(symbol, "asymbol", label);
    return symbol;
  }
  uint32_t ReadSymbol(const ANSAdaptiveSymbol& asym, uint32_t max_index,
                      WP2_OPT_LABEL) {
    const uint32_t symbol = ReadSymbolInternal(asym, max_index);
    Trace("%s: symbol=%u", label, symbol);
    BitTrace(symbol, "asymbol", label);
    return symbol;
  }

  // Decodes a binary symbol with probability 'proba' in range [0, PROBA_MAX].
  uint32_t ReadBit(uint32_t proba, WP2_OPT_LABEL) {
    const uint32_t bit = ReadBitInternal(proba);
    Trace("%s: bit=%u proba=%p", label, bit, proba);
    BitTrace(bit, "bit", label);
    return bit;
  }

  // Decodes an adaptive binary symbol with statistics 'stats'.
  // 'stats' is updated upon return.
  uint32_t ReadABit(ANSBinSymbol* const stats, WP2_OPT_LABEL) {
    const uint32_t bit = ReadABitInternal(stats);
    Trace("%s: abit=%u", label, bit);
    BitTrace(bit, "abit", label);
    return bit;
  }

  // Decodes a bool with a fifty-fifty probability.
  inline bool ReadBool(WP2_OPT_LABEL) { return (ReadUValue(1, label) != 0); }

  // Decodes a uniform value known to be in range [0..1 << bits), with 'bits'
  // in [0, IO_BITS] range.
  // If 'bits' == 0, it returns 0.
  uint32_t ReadUValue(uint32_t bits, WP2_OPT_LABEL) {
    const uint32_t value = ReadUValueInternal(bits);
    Trace("%s: value=0x%x bits=%u", label, value, bits);
    BitTrace(value, "uvalue", label);
    return value;
  }

  // Same as ReadUValue() with signed values in [ -2^(bits-1) .. 2^(bits-1) ).
  inline int32_t ReadSUValue(uint32_t bits, WP2_OPT_LABEL) {
    return (int32_t)ReadUValue(bits, label) - ((1 << (bits)) >> 1);
  }

  // Decodes a uniform value known to be in [0..range), an interval fitting in
  // kANSMaxRangeBits bits.
  uint32_t ReadRValue(uint32_t range, WP2_OPT_LABEL) {
    const uint32_t value = ReadRValueInternal(range);
    Trace("%s: value=0x%x range=0x%x", label, value, range);
    BitTrace(value, "rvalue", label);
    return value;
  }

  // Decodes a uniform value known to be in [min..max], an interval fitting in
  // kANSMaxRangeBits bits.
  inline int32_t ReadRange(int32_t min, int32_t max, WP2_OPT_LABEL) {
    assert(min <= max);
    return ReadRValue(max - min + 1, label) + min;
  }

  // Returns true if no error occurred.
  WP2Status GetStatus() const { return status_; }

  // Add something to the prefix that will appear when adding bits.
  void AddDebugPrefix(const char prefix[]) {
#if defined(WP2_BITTRACE) || defined(WP2_TRACE) || defined(WP2_ENC_DEC_MATCH)
    debug_prefix_ += prefix;
#if defined(WP2_BITTRACE)
    counters_[kPrefixStr + debug_prefix_].bits = 0.f;
    ++counters_[kPrefixStr + debug_prefix_].num_occurrences;
#endif
    debug_prefix_ += "/";
#else
    (void)prefix;
#endif
  }
  // Pop the latest addition to the prefix.
  void PopDebugPrefix() {
#if defined(WP2_BITTRACE) || defined(WP2_TRACE) || defined(WP2_ENC_DEC_MATCH)
    assert(!debug_prefix_.empty());
    int i = debug_prefix_.size() - 2;
    while (i >= 0 && debug_prefix_[i] != '/') --i;
    debug_prefix_.erase(i + 1);
#endif
  }

 private:
  // Logs something and makes sure the label is the same.
  template<typename ...Ts>
  inline void Trace(const char* format, const char label[], Ts...extra) {
    (void)sizeof...(extra);
    (void)label;
    WP2Trace(format, debug_prefix_, label, extra...);
#if defined(WP2_ENC_DEC_MATCH)
#if defined(WP2_ENC_DEC_DEEP_MATCH)
    const std::string str = ANSString(debug_prefix_, label, extra...);
    std::string str_read(ReadUValueInternal(8), '.');
    for (char& c : str_read) c = (char)ReadUValueInternal(8);
    if (str_read != str) assert(false);
#else
    const uint16_t hash = ANSStringHash(debug_prefix_, label, extra...);
    const uint16_t hash_read = ReadUValueInternal(16);
    if (hash_read != hash) assert(false);
#endif  // defined(WP2_ENC_DEC_DEEP_MATCH)
#endif  // defined(WP2_ENC_DEC_MATCH)
  }

  // internal ReadXXX versions without label[]
  uint32_t ReadSymbolInternal(const ANSSymbolInfo codes[],
                              uint32_t log2_tab_size);
  uint32_t ReadSymbolInternal(const ANSSymbolInfo codes[],
                              uint32_t log2_tab_size, uint32_t max_index);
  uint32_t ReadSymbolInternal(const ANSAdaptiveSymbol& asym);
  uint32_t ReadSymbolInternal(const ANSAdaptiveSymbol& asym,
                              uint32_t max_index);
  uint32_t ReadRValueInternal(uint32_t range);
  uint32_t ReadUValueInternal(uint32_t bits);
  uint32_t ReadABitInternal(ANSBinSymbol* const stats);
  uint32_t ReadBitInternal(uint32_t p0);

  // Shift 'state' and read-in the lower bits. Return the new value.
  // Typically called as: 'state_ = ReadNextWord(state_);'
  inline uint32_t ReadNextWord(uint32_t state);

  DataSource* data_source_;
  uint32_t state_;
  WP2Status status_;
#if defined(WP2_BITTRACE) || defined(WP2_TRACE) || defined(WP2_ENC_DEC_MATCH)
  std::string debug_prefix_;
#endif

 public:
  // Returns the number of bytes read from the bitstream.
  size_t GetNumUsedBytes() const;
  // Returns a lower bound of the number of bytes completely decoded between
  // 'num_bytes_before' and 'num_bytes_after'.
  static size_t GetMinNumUsedBytesDiff(size_t num_used_bytes_before,
                                       size_t num_used_bytes_after);

 protected:
#if defined(WP2_BITTRACE)
  double last_pos_ = 0;
  double cur_pos_ = 0;
  WP2BitCounts counters_;
  // Next reads will be registered in 'counters_custom_' with this key.
  std::string counters_custom_key_;
  // Custom counter on which the user has write access (e.g. to
  // initialize/clean).
  WP2BitCounts counters_custom_;
  // If true, merge custom bit traces until next Pop().
  bool merge_custom_until_pop_ = false;
  // Tells whether the padding has been added to counters_ yet.
  bool is_padding_counted_ = false;
  void BitTrace(uint32_t value, const char* const type, WP2_OPT_LABEL);

 public:
  // Precision of GetBitCount() compared to the actual used bits is at most
  // kANSPaddingCost or 0.16%.
  static constexpr double kBitCountAccuracy = 0.0016;

  // This string will be prepended to debug prefixes when added to the ANSDec
  // to keep track of their call number. It should be a string the user
  // will not use in a prefix.
  static constexpr const char* kPrefixStr = "____BITTRACE____";
  const WP2BitCounts& GetBitTraces() const { return counters_; }
  WP2BitCounts& GetBitTracesCustom() { return counters_custom_; }
  void ClearBitTracesCustom() { counters_custom_.clear(); }
  const std::string& GetBitTracesCustomKey() const {
    return counters_custom_key_;
  }
  void PushBitTracesCustomPrefix(const char* const suffix_of_prefix,
                                 bool merge_until_pop = false);
  void PopBitTracesCustomPrefix(const char* const suffix_of_prefix);
  double GetBitCount() const { return cur_pos_; }
#else
  static inline void BitTrace(uint32_t value, const char* const type,
                              WP2_OPT_LABEL) {
    (void)value;
    (void)type;
    (void)label;
  }

 public:
  void PushBitTracesCustomPrefix(const char* const, bool = false) {}
  void PopBitTracesCustomPrefix(const char* const) {}
  double GetBitCount() const { return 0.f; }
#endif
};

//------------------------------------------------------------------------------
// Encoding

class ANSDictionary : public WP2Allocable {
 public:
  ~ANSDictionary() noexcept { Clear(); }
  // Allocates and copy data from 'dict'.
  WP2Status CopyFrom(const ANSDictionary& dict);

  // Initializes dictionary and allocate memory of holding records
  // of 'max_symbol'. The previous memory (if any) is released.
  WP2Status Init(uint32_t max_symbol);

  // Deallocates memory and reset the dictionary structure.
  void Clear();

  // Records use of a given symbol, 'count' times.
  void RecordSymbol(uint32_t symbol, uint32_t count = 1);
  // Number of recorded symbols.
  uint32_t Total() const { return total_; }

  // Returns the current cost in bits for the symbol. Symbol is not recorded.
  float SymbolCost(uint32_t symbol) const;
  float SymbolCost(uint32_t symbol, uint32_t max_symbol) const;
  // Returns the total cost for all symbols.
  float Cost() const;

  // Processes a closed dictionary to populate infos_[] prior to coding.
  // Fails if there are too many symbols (max_symbol_ > ANS_MAX_SYMBOLS).
  WP2Status ToCodingTable();

  // Resets counts (both regular and quantized).
  void ResetCounts();
  // Returns the current counts vector.
  const Vector_u32& Counts() const { return counts_; }

  // Sets counts to a specific value.
  // When recording a symbol afterwards, the quantized counts will not be
  // updated while the normal counts will be.
  WP2Status SetQuantizedCounts(const Vector_u32& quantized_counts);
  // Returns true if quantized_counts have been set.
  bool IsQuantized() const { return !quantized_counts_.empty(); }

  uint32_t MaxSymbol() const { return max_symbol_; }
  const ANSSymbolInfo& GetInfo(uint32_t symbol) const { return infos_[symbol]; }

  uint32_t Log2TabSize() const { return log2_tab_size_; }

  // Returns whether this instance contains mostly the same data as 'other' and
  // thus is probably a clone.
  bool SeemsEquivalentTo(const ANSDictionary& other) const {
    return (max_symbol_ == other.max_symbol_ &&
            total_ == other.total_ &&
            total_quantized_ == other.total_quantized_ &&
            counts_.size() == other.counts_.size() &&
            quantized_counts_.size() == other.quantized_counts_.size() &&
            log2_tab_size_ == other.log2_tab_size_ &&
            infos_.size() == other.infos_.size());
  }

 private:
  uint32_t max_symbol_;
  // Sum of non quantized counts.
  uint32_t total_;
  // Sum of quantized counts.
  uint32_t total_quantized_;
  // counts of symbols.
  Vector_u32 counts_;
  // quantized counts (non-empty if some quantized counts have been set)
  Vector_u32 quantized_counts_;
  // the following params must be set before calling ANSEnc::Assemble()
  uint32_t log2_tab_size_;  // final expect total of counts[], used for coding
  // (0="default to ANS_LOG_TAB_SIZE")
  // Must be in [0,16] range.

  VectorNoCtor<ANSSymbolInfo> infos_;    // used during final coding
};

//------------------------------------------------------------------------------
// Decorated container classes for holding pointers

// Class containing several ANSDictionary pointers.
class ANSDictionaries : public VectorNoCtor<ANSDictionary*> {
 public:
  ~ANSDictionaries() { DeepClear(); }

  // Frees everything (deletes the added pointers).
  void DeepClear() {
    for (const ANSDictionary* const d : *this) delete d;
    clear();
  }

  // Copies all dictionaries from 'dicts'.
  WP2Status CopyFrom(const ANSDictionaries& dicts);

  // Considering this ANSDictionaries to be an intact deep clone of 'original',
  // returns the equivalent of the 'dict_to_find' or null if none.
  ANSDictionary* GetEquivalent(const ANSDictionaries& original,
                               const ANSDictionary* const dict_to_find) const;

  // The dictionary is placed last, and can be immediately access using Back().
  WP2Status Add(uint32_t max_symbol);

  // Appends the dictionaries of 'in' into the invoking object.
  // As we are dealing with pointers, 'in' loses ownership and is cleared.
  WP2Status AppendAndClear(ANSDictionaries* const in);

  // Finalize the dictionaries.
  WP2Status ToCodingTable();
};

// Class containing several ANSBinSymbol pointers.
class ANSAdaptiveBits : public VectorNoCtor<ANSBinSymbol> {
 public:
  WP2Status CopyFrom(const ANSAdaptiveBits& dicts);
};

// Class containing several ANSAdaptiveSymbol pointers.
class ANSAdaptiveSymbols : public VectorNoCtor<ANSAdaptiveSymbol> {
 public:
  WP2_NO_DISCARD bool Add(uint32_t max_symbol, bool resize_if_needed = true);
  WP2Status CopyFrom(const ANSAdaptiveSymbols& dicts);
};

//------------------------------------------------------------------------------
// Encoding structure for storing message and generating bitstream.

// Base class for ANS encoders.

class ANSEncBase {
 public:
  virtual ~ANSEncBase() = default;

  // Appends a binary symbol 'bit' with probability 'proba' (of 'bit' being
  // zero) to the message. The probability is in range [0, PROBA_MAX].
  // Returns the value of 'bit'.
  virtual uint32_t PutBit(uint32_t bit, uint32_t proba, WP2_OPT_LABEL) = 0;

  // Appends an adaptive binary symbol 'bit', updating 'stats' afterward.
  virtual uint32_t PutABit(uint32_t bit, ANSBinSymbol* const stats,
                           WP2_OPT_LABEL) = 0;

  // Appends a bool with a fifty-fifty probability. Takes one bit.
  inline bool PutBool(bool value, WP2_OPT_LABEL) {
    PutUValue(value ? 1 : 0, 1, label);
    return value;
  }

  // Appends a symbol 'symbol' to the message.
  // 'symbol' should be in range [0, ANS_MAX_SYMBOLS).
  // Returns the symbol.
  virtual uint32_t PutSymbol(uint32_t symbol, const ANSDictionary& dict,
                             WP2_OPT_LABEL) = 0;

  // Appends a symbol 'symbol' to the message using an adaptive small
  // dictionary. 'asym' is *NOT* updated after coding. 'symbol' should be in
  // range [0, APROBA_MAX_SYMBOL). Returns the symbol.
  virtual uint32_t PutSymbol(uint32_t symbol, const ANSAdaptiveSymbol& asym,
                             WP2_OPT_LABEL) = 0;

  // Variants with a max value.
  virtual uint32_t PutSymbol(uint32_t symbol, uint32_t max_symbol,
                             const ANSDictionary& dict, WP2_OPT_LABEL) = 0;
  virtual uint32_t PutSymbol(uint32_t symbol, uint32_t max_symbol,
                             const ANSAdaptiveSymbol& asym, WP2_OPT_LABEL) = 0;

  // Appends a symbol 'symbol' to the message using an adaptive small
  // dictionary. 'asym' is updated after coding. 'symbol' should be in
  // range [0, APROBA_MAX_SYMBOL). Returns the symbol.
  uint32_t PutASymbol(uint32_t symbol, ANSAdaptiveSymbol* const asym,
                      WP2_OPT_LABEL) {
    PutSymbol(symbol, *asym, label);
    asym->Update(symbol);
    return symbol;
  }

  // Appends a 'value' uniformly distributed in the range [0..1 << bits),
  // with 'bits' in range [0, IO_BITS]. Returns the value.
  // If 'bits' == 0, a value of 0 is stored.
  virtual uint32_t PutUValue(uint32_t value, uint32_t bits, WP2_OPT_LABEL) = 0;

  // Same as PutUValue() with signed values in [ -2^(bits-1) .. 2^(bits-1) ).
  inline int32_t PutSUValue(int32_t value, uint32_t bits, WP2_OPT_LABEL) {
    const int32_t half_range = ((1 << (bits)) >> 1);
    assert(bits == 0 || (-half_range <= value && value < half_range));
    PutUValue((uint32_t)(value + half_range), bits, label);
    return value;
  }

  // Appends a 'value' uniformly distributed in [0..range), an interval fitting
  // in kANSMaxRangeBits bits. Returns the value, for convenience.
  virtual uint32_t PutRValue(uint32_t value, uint32_t range, WP2_OPT_LABEL) = 0;

  // Appends a 'value' uniformly distributed in [min..max], an interval fitting
  // in kANSMaxRangeBits bits. Returns the value, for convenience.
  inline int32_t PutRange(int32_t value, int32_t min, int32_t max,
                          WP2_OPT_LABEL) {
    assert(min <= value && value <= max);
    assert((max - min + 1) <= (int32_t)kANSMaxRange);
    PutRValue(value - min, max - min + 1, label);
    return value;
  }

  // Returns the estimated cost of the final stream in bits.
  virtual float GetCost() const = 0;
  // Returns the estimated cost of the final stream in bits.
  // This overload is a speed optimization: the cost for symbols will be
  // computed per dictionary, instead of per symbol.
  // There can be a small difference with GetCost() which does a cumulative
  // sum over floats.
  virtual float GetCost(const ANSDictionaries& dicts) const = 0;
  // Same as GetCost but also adds the cost of the few extra padding words
  float GetCostFull() const;
  float GetCostFull(const ANSDictionaries& dicts) const;

  // Returns the number of tokens (bits, symbols, uvalues, rvalues...) that have
  // been buffered in this encoder.
  virtual uint32_t NumTokens() const = 0;

  virtual WP2Status GetStatus() const = 0;

  // Adds something to the prefix that will appear when adding bits.
  void AddDebugPrefix(const char prefix[]) {
#if defined(WP2_BITTRACE) || defined(WP2_TRACE) || defined(WP2_ENC_DEC_MATCH)
    debug_prefix_ += prefix;
    debug_prefix_ += "/";
#else
    (void)prefix;
#endif
  }
  // Pops the latest addition to the prefix.
  void PopDebugPrefix() {
#if defined(WP2_BITTRACE) || defined(WP2_TRACE) || defined(WP2_ENC_DEC_MATCH)
    assert(!debug_prefix_.empty());
    int i = debug_prefix_.size() - 2;
    while (i >= 0 && debug_prefix_[i] != '/') --i;
    debug_prefix_.erase(i + 1);
#endif
  }
  // Copies the prefix from another ANSEnc.
  void CopyDebugPrefix(const ANSEncBase& from) {
#if defined(WP2_BITTRACE) || defined(WP2_TRACE) || defined(WP2_ENC_DEC_MATCH)
    debug_prefix_ = from.debug_prefix_;
#endif
    (void)from;
  }

 protected:
#if defined(WP2_BITTRACE) || defined(WP2_TRACE) || defined(WP2_ENC_DEC_MATCH)
  std::string debug_prefix_;
#endif
};

struct ANSTokenPage;  // Forward declaration.

class ANSEnc : public ANSEncBase {
 public:
  ANSEnc();
  ~ANSEnc() override;
  ANSEnc(ANSEnc&&) noexcept;
  ANSEnc(const ANSEnc&) = delete;
  ANSEnc& operator=(const ANSEnc&) = delete;

  // Deep copies the given ANSEnc into the current ANSEnc.
  WP2Status Clone(const ANSEnc& e);

  // Deallocates the coded buffer generated by Assemble().
  // Will not deallocate the token pages / stats / dictionaries.
  void ResetBuffer();

  // Generates a bitstream from the current state.
  // The bitstream is stored internally and can be accessed using
  // Buffer() and BufferSize(). If any symbol was used, 'dicts' must not be
  // nullptr.
  WP2Status Assemble();

  // Reclaims all memory, reset the encoder.
  // (= reset the bistream, the message and free old pages).
  void Clear();

  // Resets the stored message (but not the bitstream, if any).
  // If delete_pages is true, past message pages are free'd. Otherwise, they are
  // marked for re-use in a future message.
  void Reset(bool delete_pages);

  // Frees write-buffer memory (in case of error, for instance)
  void WipeOut();

  // Appends the tokens of 'enc' into the invoking object.
  WP2Status Append(const ANSEnc& enc);

  // Appends tokens from 'enc' into the invoking object, starting from
  // 'start_token' (0 indexed), and stopping after the given number of tokens.
  WP2Status AppendTokens(const ANSEnc& enc, uint32_t start_token,
                         uint32_t n_tokens);

  // Retrieves the coded buffer and size (after the last call to Assemble()).
  // This buffer is only valid until the next call to Assemble() or WipeOut().
  const uint8_t* Buffer() const { return buffer_.data() + buffer_pos_; }
  size_t BufferSize() const { return (buffer_.size() - buffer_pos_); }

 public:   // overridden methods
  uint32_t PutBit(uint32_t bit, uint32_t proba, WP2_OPT_LABEL) override {
    const uint32_t res = PutBitInternal(bit, proba);
    Trace("%s: bit=%d proba=%p", label, bit, proba);
    return res;
  }
  uint32_t PutABit(uint32_t bit, ANSBinSymbol* const stats,
                   WP2_OPT_LABEL) override {
    const uint32_t res = PutABitInternal(bit, stats);
    Trace("%s: abit=%d", label, bit);
    return res;
  }
  uint32_t PutSymbol(uint32_t symbol, const ANSDictionary& dict,
                     WP2_OPT_LABEL) override {
    // TODO(vrabaud) Have it work with different log2_table.
    assert(dict.Log2TabSize() == ANS_LOG_TAB_SIZE);
    PutSymbolInternal(dict.GetInfo(symbol), /*is_adaptive=*/false);
    Trace("%s: symbol=%u", label, symbol);
    return symbol;
  }
  uint32_t PutSymbol(uint32_t symbol, const ANSAdaptiveSymbol& asym,
                     WP2_OPT_LABEL) override {
    PutSymbolInternal(asym.GetInfo(symbol), /*is_adaptive=*/true);
    Trace("%s: symbol=%u", label, symbol);
    return symbol;
  }
  uint32_t PutSymbol(uint32_t symbol, uint32_t max_symbol,
                     const ANSDictionary& dict, WP2_OPT_LABEL) override {
    return PutSymbolInternal(symbol, max_symbol, dict, label);
  }
  uint32_t PutSymbol(uint32_t symbol, uint32_t max_symbol,
                     const ANSAdaptiveSymbol& asym, WP2_OPT_LABEL) override {
    return PutSymbolInternal(symbol, max_symbol, asym, label);
  }
  uint32_t PutUValue(uint32_t value, uint32_t bits, WP2_OPT_LABEL) override {
    const uint32_t res = PutUValueInternal(value, bits);
    Trace("%s: value=0x%x bits=%u", label, value, bits);
    return res;
  }
  uint32_t PutRValue(uint32_t value, uint32_t range, WP2_OPT_LABEL) override {
    const uint32_t res = PutRValueInternal(value, range);
    Trace("%s: value=0x%x range=0x%x", label, value, range);
    return res;
  }
  float GetCost() const override;
  float GetCost(const ANSDictionaries& dicts) const override;

  WP2Status GetStatus() const override { return status_; }

  uint32_t NumTokens() const override;

 private:
  friend void swap(ANSEnc& e1, ANSEnc& e2);

  // internal PutXXX versions without label[]
  uint32_t PutBitInternal(uint32_t bit, uint32_t proba);
  uint32_t PutABitInternal(uint32_t bit, ANSBinSymbol* const stats);
  // is_adaptive indicates whether the symbol is adaptive or not.
  void PutSymbolInternal(const ANSSymbolInfo& info, bool is_adaptive);
  uint32_t PutUValueInternal(uint32_t value, uint32_t bits);
  uint32_t PutRValueInternal(uint32_t value, uint32_t range);
  template <class T>    // T must have a ::GetInfo(symbol)
  uint32_t PutSymbolInternal(uint32_t symbol, uint32_t max_symbol,
                             const T& dict, WP2_OPT_LABEL) {
    assert(symbol <= max_symbol);
    const ANSSymbolInfo info_max = dict.GetInfo(max_symbol);
    const ANSSymbolInfo info = dict.GetInfo(symbol).Rescale(info_max);
    PutSymbolInternal(info, /*is_adaptive=*/false);
    Trace("%s: symbol=%u", label, info.symbol);
    return info.symbol;
  }

  template<typename ...Ts>
  inline void Trace(const char* format, const char label[], Ts...extra) {
    (void)sizeof...(extra);
    (void)label;
    WP2Trace(format, debug_prefix_, label, extra...);
#if defined(WP2_ENC_DEC_MATCH)
#if defined(WP2_ENC_DEC_DEEP_MATCH)
    const std::string str = ANSString(debug_prefix_, label, extra...);
    PutUValueInternal(str.size(), 8);
    for (char c : str) PutUValueInternal(c, 8);
#else
    const uint16_t hash = ANSStringHash(debug_prefix_, label, extra...);
    if (GetStatus() == WP2_STATUS_OK) PutUValueInternal(hash, 16);
#endif  // defined(WP2_ENC_DEC_DEEP_MATCH)
#endif  // defined(WP2_ENC_DEC_MATCH)
  }

  void EnqueueToken(uint32_t tok);
  bool NewPage();

  // to be called last
  WP2Status FinishStream();
  WP2Status EmitOnlyBits();
  WP2Status EmitAll();

  // emits the lower bits of 's' and returns the new value
  inline uint32_t EmitWord(uint32_t s);

  void EmitBit(uint32_t bit, uint32_t p0);
  void EmitASymbol(uint32_t offset, uint32_t freq);
  void EmitUValue(uint32_t value, int bits);
  void EmitRValue(uint32_t value, uint32_t range);

  bool ReallocBuffer();

  WP2Status status_;
  uint32_t state_;         // state for the coder

  Vector_u8 buffer_;       // output coded stream
  size_t buffer_pos_;      // position from the *end* of buffer_

  // Stored message, composed of "tokens" which can be bits, uvalues or symbols.
  // Pointer to the most recent (last) page of tokens. Use tokens_->prev_ to
  // access previous pages.
  ANSTokenPage* tokens_;
  // Pointer to old allocated pages, which can be reused next time we need a
  // new page (see the delete_pages parameter of Reset()).
  ANSTokenPage* list_;
  // misc stats:
  int nb_tokens_;          // total number of tokens in the message
  int nb_asymbols_;        // number of adaptive symbols stored
  int nb_uvalues_;         // number of uniform values stored
  int nb_rvalues_;         // number of range-based values stored
  int nb_bits_;            // number of binary symbols stored
};

void swap(ANSEnc& e1, ANSEnc&e2);

//------------------------------------------------------------------------------

}   // namespace WP2

#endif   // WP2_UTILS_ANS_H_
