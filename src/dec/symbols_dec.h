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
// ANS symbol decoding.
//
// Author: Vincent Rabaud (vrabaud@google.com)

#ifndef WP2_DEC_SYMBOLS_DEC_H_
#define WP2_DEC_SYMBOLS_DEC_H_

#include "src/common/symbols.h"
#include "src/utils/ans.h"

namespace WP2 {

// Class simplifying the reading of symbols.
// A symbol can be:
// - trivial if it is always the same value
// - from  a specific dictionary depending on its type
struct StatExtra {
  const ANSSymbolInfo* codes;
  uint32_t log2_tab_size;
};

class SymbolReader : public SymbolIO<StatExtra> {
 public:
  WP2Status Init(const SymbolsInfo& symbols_info, ANSDec* const dec);

  const SymbolsInfo& symbols_info() const { return symbols_info_; }

  // Sets the range for a particular symbol.
  inline void SetRange(uint32_t sym, int32_t min, int32_t max) {
    symbols_info_.SetMinMax(sym, min, max);
  }

  // For the 'cluster'th cluster of symbol 'sym', sets the mode to trivial value
  // with value 'value'.
  void AddTrivial(uint32_t cluster, uint32_t sym, int32_t value);
  // For the 'cluster'th cluster of symbol 'sym' (given a maximum value of
  // non-zero values of 'max_nnz', which can be the number of pixels for
  // examples) reads its header to later decide how to read it.
  WP2Status ReadHeader(uint32_t cluster, uint32_t max_nnz, uint32_t sym,
                       WP2_OPT_LABEL);
  // Same as above but reads headers for all clusters.
  WP2Status ReadHeader(uint32_t max_nnz, uint32_t sym, WP2_OPT_LABEL);
  // Reads the symbol of type 'sym' in the 'cluster'th cluster.
  // If 'cost' is not nullptr, the cost in bits of storing 'sym' is ADDED to
  // 'cost'.
  int32_t Read(uint32_t cluster, uint32_t sym, WP2_OPT_LABEL,
               double* const cost = nullptr);
  int32_t Read(uint32_t sym, WP2_OPT_LABEL, double* const cost = nullptr) {
    return Read(/*cluster=*/0, sym, label, cost);
  }
  // Same as above except we know that the read value cannot be strictly
  // superior to 'max_value'.
  WP2Status TryRead(uint32_t cluster, uint32_t sym, uint32_t max_value,
                    WP2_OPT_LABEL, int32_t* const value,
                    double* const cost = nullptr);
  // Fills the array 'is_maybe_used' with true when the symbol might be used,
  // false if we are sure it is not used.
  void GetPotentialUsage(uint32_t cluster, uint32_t sym, bool is_maybe_used[],
                         size_t size) const override;

 protected:
  // For the 'cluster'th cluster of symbol 'sym', sets the mode to range value,
  // in [0:'max_range'[
  // 'infos' is only useful for the mapping.
  void AddRange(uint32_t cluster, uint32_t sym,
                const VectorNoCtor<ANSSymbolInfo>* const infos,
                uint16_t max_range);
  // For the 'cluster'th cluster of symbol 'sym', sets the mode to dictionary
  // value, using info from 'infos'. 'infos' is modified in this call.
  WP2Status AddDict(uint32_t cluster, uint32_t sym,
                    VectorNoCtor<ANSSymbolInfo>* const infos);
  // For the 'cluster'th cluster of symbol 'sym', sets the mode to Golomb
  // value, using info from 'infos'. 'infos' is modified in this call.
  // 'prefix_size' is the number of bits in the Golomb prefix.
  WP2Status AddGolomb(uint32_t cluster, uint32_t sym,
                      VectorNoCtor<ANSSymbolInfo>* const infos,
                      uint32_t prefix_size);

 private:
  // Generic implementation called by the different Read specializations above.
  WP2Status ReadInternal(uint32_t cluster, uint32_t sym, bool use_max_value,
                         uint32_t max_value, WP2_OPT_LABEL,
                         int32_t* const value, double* const cost = nullptr);
  // Finds the biggest symbol 'sym' such that stat.mappings[sym] <= max_value.
  uint32_t FindSymbol(const Stat& stat, uint32_t max_value);

  ANSDec* dec_;
  // Tmp variables to use during the update.
  Vector_u16 mapping_;
  // Codes used for fast probability check for symbols relying on dictionaries.
  typedef VectorNoCtor<ANSSymbolInfo> Codes;
  Vector<Codes> codes_;
  VectorNoCtor<ANSSymbolInfo> infos_;
};

}  // namespace WP2

#endif /* WP2_DEC_SYMBOLS_DEC_H_ */
