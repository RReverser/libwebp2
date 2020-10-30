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
// ANS symbol encoding.
//
// Author: Vincent Rabaud (vrabaud@google.com)

#ifndef WP2_ENC_SYMBOLS_ENC_H_
#define WP2_ENC_SYMBOLS_ENC_H_

#include <array>
#include <initializer_list>
#include <limits>

#include "src/common/lossy/residuals.h"
#include "src/common/symbols.h"
#include "src/utils/ans.h"
#include "src/utils/quantizer.h"
#include "src/utils/vector.h"

namespace WP2 {

// Interface for a class that processes symbols (to record/write/count them...)
class SymbolManager {
 public:
  virtual ~SymbolManager() = default;

  // Processes the symbol of type 'sym' and value 'value' in the 'cluster'th
  // cluster. Returns the value.
  int32_t Process(uint32_t cluster, uint32_t sym, int32_t value, WP2_OPT_LABEL,
                  ANSEncBase* const enc) {
    return ProcessInternal(cluster, sym, value, /*use_max_value=*/false,
                           /*max_value=*/0, label, enc, /*cost=*/nullptr);
  }

  int32_t Process(uint32_t sym, int32_t value, WP2_OPT_LABEL,
                  ANSEncBase* const enc) {
    return Process(/*cluster=*/0, sym, value, label, enc);
  }

  // Does the same as above except that it caps the value to be in
  // [0, max_value]. This overload can be costly on the decoder side as it
  // requires a binary search in [0, max_value]. This function is usually to
  // be used at the end of an interval where we know we can cap values because
  // they cannot be out of it.
  int32_t Process(uint32_t cluster, uint32_t sym, int32_t value,
                  uint32_t max_value, WP2_OPT_LABEL, ANSEncBase* const enc) {
    return ProcessInternal(cluster, sym, value, /*use_max_value=*/true,
                           max_value, label, enc, /*cost=*/nullptr);
  }

  // Same as process, but if 'cost' is not nullptr, it is set to the cost of
  // storing this symbol.
  void ProcessWithCost(uint32_t cluster, uint32_t sym, int32_t value,
                       WP2_OPT_LABEL, ANSEncBase* const enc,
                       double* const cost) {
    ProcessInternal(cluster, sym, value, /*use_max_value=*/false,
                    /*max_value=*/0, label, enc, cost);
  }
  void ProcessWithCost(uint32_t cluster, uint32_t sym, int32_t value,
                       uint32_t max_value, WP2_OPT_LABEL, ANSEncBase* const enc,
                       double* const cost) {
    ProcessInternal(cluster, sym, value, /*use_max_value=*/true, max_value,
                    label, enc, cost);
  }

  virtual const SymbolsInfo& symbols_info() const = 0;

 protected:
  virtual int32_t ProcessInternal(uint32_t cluster, uint32_t sym, int32_t value,
                                  bool use_max_value, uint32_t max_value,
                                  WP2_OPT_LABEL, ANSEncBase* const enc,
                                  double* const cost) = 0;
};

// Class for counting the usage of symbols.
class SymbolRecorder : public SymbolManager {
 public:
  const SymbolsInfo& symbols_info() const override { return symbols_info_; }

  // Initializes and allocates memory. 'num_records' is an upper bound on the
  // number of records for a symbol with kAdaptiveWithAutoSpeed, so that storage
  // can be pre-allocated. If more symbols will need to be recorded, code will
  // fail.
  WP2Status Allocate(const SymbolsInfo& symbols_info, uint32_t num_records);
  WP2Status CopyFrom(const SymbolRecorder& other);
  void DeepClear();

  // Copies the dictionaries to a backup: useful for multi-pass.
  WP2Status MakeBackup();

  // Resets stats.
  void ResetRecord(bool reset_backup);

  const ANSDictionary& GetRecordedDict(uint32_t cluster, uint32_t sym) const;
  const Vector_u8& GetRecordedValues(uint32_t cluster, uint32_t sym) const;
  const ANSBinSymbol& GetABit(uint32_t cluster, uint32_t sym) const;
  const ANSAdaptiveSymbol& GetASymbol(uint32_t cluster, uint32_t sym) const;
  // Returns the current value updated with everything stored so far at default
  // speed.
  const ANSAdaptiveSymbol& GetASymbolWithSpeed(uint32_t cluster,
                                               uint32_t sym) const;
  float GetCost(uint32_t cluster, uint32_t sym, int32_t value) const;

  // Returns the dictionary of recorded values from the previous pass (before
  // the last call to ResetRecord()).
  const ANSDictionary& GetDictPreviousPass(uint32_t cluster,
                                           uint32_t sym) const;
  // Returns the dictionary of recorded values from the previous pass. Can be
  // used to manually set initial statistics.
  ANSDictionary& GetDictPreviousPass(uint32_t cluster, uint32_t sym);

 protected:
  int32_t ProcessInternal(uint32_t cluster, uint32_t sym, int32_t value,
                          bool use_max_value, uint32_t max_value, WP2_OPT_LABEL,
                          ANSEncBase* const enc, double* const) override;

 private:
  ANSDictionary* GetDict(uint32_t cluster, uint32_t sym);
  ANSBinSymbol* GetABit(uint32_t cluster, uint32_t sym);
  ANSAdaptiveSymbol* GetASymbol(uint32_t cluster, uint32_t sym);

  SymbolsInfo symbols_info_;
  // For each symbol 'sym', the statistics for all clusters start at
  // index_[sym], either in dicts_, a_bits_, a_symbols_ or values_ depending on
  // the type.
  std::array<uint32_t, kSymbolNumMax> index_;
  // Index for ASymbol in the values_ vector.
  std::array<uint32_t, kSymbolNumMax> values_index_;
  // Dictionaries, for StorageMethod::kAuto
  ANSDictionaries dicts_;
  ANSDictionaries dicts_previous_pass_;
  // Adaptive bits, for StorageMethod::kAdaptive with range == 2
  ANSAdaptiveBits a_bits_;
  // Adaptive symbols, for StorageMethod::kAdaptive with range != 2
  ANSAdaptiveSymbols a_symbols_;
  // Full values, for StorageMethod::kAdaptiveWithAutoSpeed
  Vector<Vector_u8> values_;
};

// Class for using with an ANSEncCounter to count space taken. It's meant to be
// used with an ANSEncCounter. Then call GetCost() on the counter to
// know the size.
class SymbolCounter : public SymbolManager, public WP2Allocable {
 public:
  static constexpr uint32_t kMaxLossyClusters =
      kResidualClusters * ResidualIO::kNumSectors;

  explicit SymbolCounter(const SymbolRecorder* const recorder)
      : recorder_(recorder), symbols_info_(recorder_->symbols_info()) {}

  const SymbolsInfo& symbols_info() const override { return symbols_info_; }

  // Allocates memory. Callers must decide in advance how many adaptive bits/
  // symbols to allocate. This number must take into account the number of
  // clusters (num distinct symbol types * num clusters). If an attempt is made
  // to call Process() with more cluster/symbol types than were allocated, an
  // assertion will fail.
  WP2Status Allocate(std::initializer_list<uint32_t> syms);
  WP2Status CopyFrom(const SymbolCounter& other);
  void Clear();

 protected:
  int32_t ProcessInternal(uint32_t cluster, uint32_t sym, int32_t value_in,
                          bool use_max_value, uint32_t max_value, WP2_OPT_LABEL,
                          ANSEncBase* const enc, double* const) override;

 private:
  ANSBinSymbol* GetABit(uint32_t cluster, uint32_t sym);
  ANSAdaptiveSymbol* GetASymbol(uint32_t cluster, uint32_t sym);

  const SymbolRecorder* const recorder_;
  const SymbolsInfo& symbols_info_;
  static constexpr uint32_t kInvalidIndex =
      std::numeric_limits<uint32_t>::max();
  // For each symbol, the index in a_bits_/a_symbols_ where its storage for each
  // cluster starts.
  std::array<uint32_t, kSymbolNumMax> indices_;
  // For ABit/ASymbol, store how many can used at most in num_*_, if they are
  // initialized in *_initialized_ and all of them in *_.
  // TODO(vrabaud) Check the dimensions once we remove the clusters in the
  //               residual code.
  // Adaptive bits, for StorageMethod::kAdaptive with range == 2
  uint32_t num_a_bits_;
  std::array<bool, 3094> a_bit_initialized_;
  ANSAdaptiveBits a_bits_;
  // Adaptive symbols, for StorageMethod::kAdaptive with range != 2
  uint32_t num_a_symbols_;
  std::array<bool, 892> a_symbol_initialized_;
  ANSAdaptiveSymbols a_symbols_;
};

// Class simplifying the writing of symbols.
// When writing a symbol to the bitstream:
//  - it does not have to be written if it is trivial (it is always the same
//  value)
//  - it is written in a certain dictionary depending on its type
// WriteHeader() must be called for each symbol type before it can be
// written with Process().
struct SymbolWriterStatExtra {
  const ANSDictionary* dict;
  size_t mapping_size;
};
class SymbolWriter : public SymbolManager,
                     public SymbolIO<SymbolWriterStatExtra>,
                     public WP2Allocable {
 public:
  const SymbolsInfo& symbols_info() const override { return symbols_info_; }

  // Allocates memory. Should be called after calling SymbolIO::Init().
  WP2Status Allocate() override;

  // Deep copy. 'copied_dicts' must be a copy of 'original_dicts'.
  WP2Status CopyFrom(const SymbolWriter& other,
                     const ANSDictionaries& original_dicts,
                     const ANSDictionaries& copied_dicts);

  // For the 'cluster'th cluster of symbol 'sym', sets the writing mode to
  // trivial.
  void AddTrivial(uint32_t cluster, uint32_t sym, int32_t value);
  // For the 'cluster'th cluster of symbol 'sym', (given a maximum value of
  // non-zero values of 'max_nnz', which can be the number of pixels for
  // examples) of statistics 'counts' (first element being the count for the
  // minimum value of the symbol), decides how to store it
  // (trivial/range/dictionary) and writes the meta-information used to later be
  // able to decode the symbols.
  WP2Status WriteHeader(uint32_t cluster, uint32_t max_nnz, uint32_t sym,
                        const uint32_t* const counts, WP2_OPT_LABEL,
                        ANSEncBase* const enc, ANSDictionaries* const dicts);
  // Same as above but uses counts from the provided SymbolRecorder.
  WP2Status WriteHeader(uint32_t cluster, uint32_t max_nnz, uint32_t sym,
                        const SymbolRecorder& syntax_recorder, WP2_OPT_LABEL,
                        ANSEncBase* const enc, ANSDictionaries* const dicts);
  // Same as above but writes headers for all clusters.
  WP2Status WriteHeader(uint32_t max_nnz, uint32_t sym,
                        const SymbolRecorder& syntax_recorder, WP2_OPT_LABEL,
                        ANSEncBase* const enc, ANSDictionaries* const dicts);

  // Fills the array 'is_maybe_used' with true when the symbol might be used,
  // false if we are sure it is not used.
  void GetPotentialUsage(uint32_t cluster, uint32_t sym, bool is_maybe_used[],
                         size_t size) const override;

 protected:
  // For the 'cluster'th cluster of symbol 'sym', sets the writing mode to
  // a range.
  // 'mapping' is used in case we don't use the raw values but a mapping.
  void AddRange(uint32_t cluster, uint32_t sym, const uint16_t* const mapping,
                uint32_t size, uint16_t max_range);
  // For the 'cluster'th cluster of symbol 'sym', sets the writing mode to
  // a dictionary of statistics represented by 'counts', 'quantized_counts' and
  // 'mapping'.
  WP2Status AddDict(uint32_t cluster, uint32_t sym,
                    const uint32_t* const counts,
                    const uint32_t* const quantized_counts,
                    const uint16_t* const mapping, size_t size,
                    ANSDictionaries* const dicts);
  // For the 'cluster'th cluster of symbol 'sym', sets the writing mode to
  // Golomb of prefix size 'prefix_size' using the statistics represented by
  // 'counts', 'quantized_counts' and 'mapping'.
  WP2Status AddGolomb(uint32_t cluster, uint32_t sym,
                      const uint32_t* const counts,
                      const uint32_t* const quantized_counts,
                      const uint16_t* const mapping, uint32_t size,
                      uint32_t prefix_size, ANSDictionaries* const dicts);

 private:
  // Writes the symbol of type 'sym' in the 'cluster'th cluster. 'WriteHeader'
  // needs to have been called before so that the SymbolWriter knows how to
  // store 'sym'. If 'cost' is not nullptr, it's set to the cost in bits of
  // storing 'sym'.
  int32_t ProcessInternal(uint32_t cluster, uint32_t sym, int32_t value_in,
                          bool use_max_value, uint32_t max_value, WP2_OPT_LABEL,
                          ANSEncBase* const enc, double* const cost) override;

  // Writes a histogram as defined by the quantizer.
  WP2Status WriteHistogram(const Quantizer::Config& config,
                           uint32_t symbol_range,
                           uint32_t max_count, ANSEncBase* const enc);
  // Computes the cost of storing a histogram as Golomb.
  void ComputeGolombCost(uint32_t range, uint32_t max_nnz, uint32_t size,
                         uint32_t prefix_size, float* cost,
                         Quantizer::Config** const config_golomb,
                         uint32_t* const size_golomb);

  // Returns the biggest index such that mapping[index] is valid and index <=
  // max_value. If no such value exists, returns the first valid index above
  // max_value.
  static uint32_t FindBiggestMappingIndex(const Stat& stat, uint32_t max_index);
  // Tmp variables to avoid multiple mallocs.
  ANSCodes infos_tmp_;
  Vector_u32 counts_tmp_;
  // Tmp variables for updating the SymbolWriter.
  Vector_u32 histogram_;
  Vector_u16 mapping_;
  Vector_u32 histogram_golomb_;
  Vector_u16 mapping_golomb_;
  Quantizer quantizer_;
  Quantizer quantizer_golomb_;
  VectorNoCtor<OptimizeArrayStorageStat> stats_buffer_;
};

}  // namespace WP2

#endif /* WP2_ENC_SYMBOLS_ENC_H_ */
