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
// Common code dealing with symbols.
//
// Author: Vincent Rabaud (vrabaud@google.com)

#ifndef WP2_COMMON_SYMBOLS_H_
#define WP2_COMMON_SYMBOLS_H_

#include <algorithm>
#include <array>
#include <cassert>
#include <limits>

#include "src/common/constants.h"
#include "src/common/lossy/block_size.h"
#include "src/utils/ans.h"
#include "src/utils/utils.h"
#include "src/utils/vector.h"
#include "src/wp2/base.h"

//------------------------------------------------------------------------------

namespace WP2 {

// Maximum number of symbol types to use.
static constexpr uint32_t kSymbolNumMax = 41;
// Number of residual clusters, not including sectors.
// Y/UV/A * Method0/Method1 * transform.
static constexpr uint32_t kResidualClusters = 3 * kNumResidualStats * TRF_LAST;

// Class containing info about the size (in bits) of the different symbols
// (colors, cache index and LZ77 distance/length).
// WARNING! This is a *big* object. Don't use on the stack!
class SymbolsInfo {
 public:
  // Bound to use when a symbol is unused.
  static constexpr int32_t kInvalidBound = std::numeric_limits<int32_t>::max();

  enum class StorageMethod {
    // Automatically choose between dictionary, range, golomb.
    kAuto,
    // Adaptive bit using counts of 0/1 so far for probabilities.
    kAdaptiveBit,
    // Adaptive symbol (that can be of size 2, but uses the adaptive symbol
    // method).
    kAdaptiveSym,
    // Adaptive symbol, with automatically chosen adaptation speed (signaled
    // in the header). Only allowed for symbols that occur at most once per
    // block and whose values fit in a uint8_t, because we keep all symbols
    // in memory which might quickly add up.
    kAdaptiveWithAutoSpeed,
    kUnused
  };

  SymbolsInfo() = default;
  SymbolsInfo(const SymbolsInfo& other) = delete;
  SymbolsInfo& operator=(const SymbolsInfo& info) = delete;

  // TODO(maryla): CopyFrom is fairly expensive (this is a big object) and I
  //               think we currently use it in some places where we don't
  //               actually need a copy.
  WP2Status CopyFrom(const SymbolsInfo& other);

  WP2Status InitLossy(uint32_t num_segments, PartitionSet partition_set,
                      bool has_alpha, uint32_t quality_hint,
                      bool use_aom_coeffs);

  bool operator==(const SymbolsInfo& info) const;  // only for tests

  void SetInfo(uint32_t sym, int32_t min, int32_t max, uint32_t num_clusters,
               StorageMethod method, bool probably_geometric = false);

  // Sets symbol range for all clusters. Indicates that symbol values are in
  // [min; max]. For now, only min=0 or min=-max is supported.
  // Overwrites any values previously set through
  // SetMinMax(sym,min,max) or SetMinMax(cluster,min,max).
  // The convention min==max==kInvalidBound means the symbol will never be used.
  void SetMinMax(uint32_t sym, int32_t min, int32_t max);
  // Sets symbol range for a specific cluster. Overrides any value previously
  // set through SetMinMax(sym,min,max) for that cluster.
  WP2Status SetMinMax(uint32_t cluster, uint32_t sym, int32_t min, int32_t max);

  void SetClusters(uint32_t sym, uint32_t num_clusters);

  WP2Status SetStartingProba(uint32_t cluster, uint32_t sym,
                             uint32_t p0, uint32_t p1);

  uint32_t StartingProbaP0(uint32_t cluster, uint32_t sym) const {
    assert(Method(sym) == StorageMethod::kAdaptiveBit &&
           Range(cluster, sym) == 2);
    return GetClusterInfo(cluster, sym).p0;
  }

  uint32_t StartingProbaP1(uint32_t cluster, uint32_t sym) const {
    assert(Method(sym) == StorageMethod::kAdaptiveBit &&
           Range(cluster, sym) == 2);
    return GetClusterInfo(cluster, sym).p1;
  }

  WP2Status SetInitialCDF(const uint16_t* const cdf, uint16_t max_proba,
                          uint32_t cluster, uint32_t sym) {
    assert(IsAdaptiveMethod(sym));
    ClusterInfo* const info = GetOrMakeClusterInfo(cluster, sym);
    WP2_CHECK_ALLOC_OK(info != nullptr);
    info->cdf = cdf;
    info->max_proba = max_proba;
    return WP2_STATUS_OK;
  }

  void InitialCDF(uint32_t cluster, uint32_t sym, const uint16_t** const cdf,
                  uint16_t* const max_proba) const {
    assert(IsAdaptiveMethod(sym));
    *cdf = GetClusterInfo(cluster, sym).cdf;
    *max_proba = GetClusterInfo(cluster, sym).max_proba;
  }

  // Range of the symbol of bounds [min,max]: max-min+1.
  inline uint32_t Range(uint32_t cluster, uint32_t sym) const {
    const ClusterInfo& info_cluster = GetClusterInfo(cluster, sym);
    if (info_cluster.min != kInvalidBound &&
        info_cluster.max != kInvalidBound) {
      return info_cluster.max - info_cluster.min + 1;
    }
    const SymbolInfo& info_sym = GetSymbolInfo(sym);
    if (info_sym.min != kInvalidBound && info_sym.max != kInvalidBound) {
      return info_sym.max - info_sym.min + 1;
    }
    return 0u;
  }

  inline int32_t Min(uint32_t cluster, uint32_t sym) const {
    const ClusterInfo& info_cluster = GetClusterInfo(cluster, sym);
    if (info_cluster.min != kInvalidBound &&
        info_cluster.max != kInvalidBound) {
      return info_cluster.min;
    }
    return GetSymbolInfo(sym).min;
  }

  inline int32_t Max(uint32_t cluster, uint32_t sym) const {
    const ClusterInfo& info_cluster = GetClusterInfo(cluster, sym);
    if (info_cluster.min != kInvalidBound &&
        info_cluster.max != kInvalidBound) {
      return info_cluster.max;
    }
    return GetSymbolInfo(sym).max;
  }

  inline StorageMethod Method(const uint32_t sym) const {
    return GetSymbolInfo(sym).storage_method;
  }

  inline bool IsAdaptiveMethod(const uint32_t sym) const {
    const StorageMethod method = Method(sym);
    return (method == StorageMethod::kAdaptiveSym ||
            method == StorageMethod::kAdaptiveWithAutoSpeed);
  }

  // Whether values for this symbol are likely to follow a geometric
  // distribution. This information can be used when estimating the cost
  // of storing this symbol.
  inline bool ProbablyGeometric(const uint32_t sym) const {
    return GetSymbolInfo(sym).probably_geometric;
  }

  inline uint32_t NumClusters(const uint32_t sym) const {
    return GetSymbolInfo(sym).num_clusters;
  }

  // Gets the sum of the max ranges of all the symbols.
  // TODO(maryla): not sure this is needed, it seems it's usually called for
  // symbols that have 1 cluster, in which case it's equivalent to RangeSum().
  uint32_t MaxRangeSum() const {
    uint32_t sum = 0;
    for (uint32_t sym = 0; sym < num_symbols_; ++sym) {
      sum += GetMaxRange(sym);
    }
    return sum;
  }
  // Gets the sum of the ranges of all the symbols/clusters.
  uint32_t RangeSum() const {
    uint32_t sum = 0;
    for (uint32_t sym = 0; sym < num_symbols_; ++sym) {
      for (uint32_t c = 0; c < NumClusters(sym); ++c) {
        sum += Range(c, sym);
      }
    }
    return sum;
  }
  // Gets the maximum range among all the symbols and clusters.
  uint32_t GetMaxRange() const {
    const auto iter =
        std::max_element(symbols_info_.begin(), symbols_info_.end(),
                         [](const SymbolInfo& lhs, const SymbolInfo& rhs) {
                           return (lhs.max - lhs.min < rhs.max - rhs.min);
                         });
    return iter->max - iter->min + 1;
  }
  // Gets the maximum range for a given symbol among all clusters.
  uint32_t GetMaxRange(uint32_t sym) const {
    uint32_t max_range = 0;
    for (uint32_t c = 0; c < NumClusters(sym); ++c) {
      max_range = std::max(max_range, Range(c, sym));
    }
    return max_range;
  }

  uint32_t Size() const { return num_symbols_; }

 private:
  struct ClusterInfo {
    // Overrides SymbolInfo.min if min <= max.
    int min = kInvalidBound;
    int max = kInvalidBound;
    // Starting probabilities for StorageMethod::kAdaptiveBit.
    uint32_t p0 = 1, p1 = 1;
    // Pointer to a reference cdf to use.
    const uint16_t* cdf = nullptr;
    // Max proba in the CDF.
    uint16_t max_proba;

    // not a real comparison operator, just used for tests
    bool operator==(const ClusterInfo& rhs) const;
  };

  struct SymbolInfo {
    StorageMethod storage_method = StorageMethod::kAuto;
    int min = kInvalidBound;
    int max = kInvalidBound;
    // Index of start of per-cluster info in cluster_info_.
    int cluster_info_index = kNoPerClusterInfo;
    uint32_t num_clusters = 0;
    // Whether values are likely to follow a geometric distribution.
    bool probably_geometric = false;

    static constexpr int kNoPerClusterInfo = -1;

    bool MinMaxIsValid(int new_min, int new_max) const;

    // not a real comparison operator, just used for tests
    bool operator==(const SymbolInfo& rhs) const;
  };

  inline SymbolInfo* GetSymbolInfo(uint32_t sym) {
    assert(sym < kSymbolNumMax);
    return &symbols_info_[sym];
  }

  inline const SymbolInfo& GetSymbolInfo(uint32_t sym) const {
    assert(sym < num_symbols_);
    return symbols_info_[sym];
  }

  // Returns the ClusterInfo for the given cluster/sym, initializing the
  // per-cluster info storage if needed.
  ClusterInfo* GetOrMakeClusterInfo(uint32_t cluster, uint32_t sym);

  const ClusterInfo& GetClusterInfo(uint32_t cluster, uint32_t sym) const;

  static const ClusterInfo kDefaultClusterInfo;

  uint32_t num_symbols_ = 0;
  std::array<SymbolInfo, kSymbolNumMax> symbols_info_;

  Vector<ClusterInfo> cluster_info_;
  uint32_t next_cluster_info_index_ = 0;
};

//------------------------------------------------------------------------------

// Split the symbol 'v' into a 'prefix' (corresponding to the first or first
// two bits depending on 'prefix_size'), the number of extra bits,
// and the value contained in those extra bits.
struct Golomb {
  Golomb(uint32_t v, uint32_t prefix_size) {
    assert(prefix_size == 1 || prefix_size == 2);
    if (prefix_size == 2) {
      if (v < 4) {
        prefix = v;
        extra_bits_num = 0;
      } else {
        const uint32_t highest_bit_index = WP2Log2Floor(v);
        const uint32_t second_highest_bit = (v >> (highest_bit_index - 1)) & 1;
        extra_bits_num = highest_bit_index - 1;
        // Store highest bit position and whether the next highest bit is 1 into
        // an int.
        prefix = 2 * highest_bit_index + second_highest_bit;
      }
    } else {
      if (v == 0) {
        prefix = 0;
        extra_bits_num = 0;
      } else {
        const uint16_t n = WP2Log2Floor(v);
        extra_bits_num = n;
        prefix = n + 1;
      }
    }
    extra_bits_value = v & ((1 << extra_bits_num) - 1);
  }
  // Given a 'prefix', returns the number of extra bits required to complete
  // it.
  static uint32_t NumExtraBits(uint32_t prefix, uint32_t prefix_size) {
    assert(prefix_size == 1 || prefix_size == 2);
    if (prefix_size == 2) {
      if (prefix < 4) return 0;
      return (prefix - 2) >> 1;
    } else {
      if (prefix <= 1) return 0;
      return prefix - 1;
    }
  }
  // Given a 'prefix' and 'extra_bits_value', combine the two to reconstruct
  // the final value.
  static uint32_t Merge(uint32_t prefix, uint32_t prefix_size,
                        uint32_t extra_bits_value) {
    assert(prefix_size == 1 || prefix_size == 2);
    if (prefix_size == 2) {
      if (prefix < 3) return prefix;
      const uint32_t extra_bits_num = (prefix - 2) >> 1;
      const uint32_t offset = (2 + (prefix & 1)) << extra_bits_num;
      return offset + extra_bits_value;
    } else {
      if (prefix < 2) return prefix;
      return (1 << (prefix - 1)) + extra_bits_value;
    }
  }

  uint32_t prefix;
  uint32_t extra_bits_num;
  uint32_t extra_bits_value;
};

// Read/write a value using Golomb encoding.
void WriteGolomb(uint32_t value, uint32_t min, uint32_t max,
                 uint32_t prefix_size, ANSEncBase* const enc, WP2_OPT_LABEL);
uint32_t ReadGolomb(uint32_t min, uint32_t max, uint32_t prefix_size,
                    ANSDec* const dec, WP2_OPT_LABEL);

//------------------------------------------------------------------------------

// Class common to decoder to read/write statistics in an image.
// It is made common to have a similar behavior when analyzing tiles.
template <class EXTRA>
class SymbolIO {
 public:
  virtual ~SymbolIO() = default;
  // Sets basic information about symbols.
  WP2Status Init(const SymbolsInfo& symbols_info) {
    WP2_CHECK_STATUS(symbols_info_.CopyFrom(symbols_info));
    return WP2_STATUS_OK;
  }

  // Allocates data after ranges/cluster sizes have been specified.
  virtual WP2Status Allocate() {
    const uint32_t size = symbols_info_.RangeSum();
    uint32_t size_contexts = 0;
    for (uint32_t sym = 0; sym < symbols_info_.Size(); ++sym) {
      size_contexts += symbols_info_.NumClusters(sym);
    }
    WP2_CHECK_ALLOC_OK(mappings_buffer_.resize(size));
    WP2_CHECK_ALLOC_OK(all_stats_.resize(size_contexts));
    // Prepare the mappings to point to the right spot in the buffer.
    auto* mapping_buffer = mappings_buffer_.data();
    for (uint32_t sym = 0, ind = 0; sym < symbols_info_.Size(); ++sym) {
      stats_start_[sym] = &all_stats_[ind];
      for (uint32_t c = 0; c < symbols_info_.NumClusters(sym); ++c) {
        all_stats_[ind].type = Stat::Type::kUnknown;
        all_stats_[ind].use_mapping = false;
        all_stats_[ind].range = 0;
        all_stats_[ind].mappings = mapping_buffer;
        mapping_buffer += symbols_info_.Range(c, sym);
        ++ind;
      }
    }
    return WP2_STATUS_OK;
  }
  // Fills the array 'is_maybe_used' with true when the symbol might be used,
  // false if we are sure it is not used.
  virtual void GetPotentialUsage(uint32_t cluster, uint32_t sym,
                                 bool is_maybe_used[], size_t size) const = 0;

 protected:
  // The stats of a symbol.
  struct Stat {
    static const uint16_t kInvalidMapping;
    enum class Type {
      // We have no information on this symbol.
      kUnknown,
      // Symbol always has the same value or is unused.
      kTrivial,
      // Symbol is stored as a range (all values assumed to have equal
      // probability).
      kRange,
      // Symbol is stored as a dictionary.
      kDict,
      // Symbol is stored with a Golomb encoding with a prefix of 1 or 2 bits.
      kGolomb,
      // Symbol is stored as an adaptive bit.
      kAdaptiveBit,
      // Symbol is stored as an adaptive symbol.
      kAdaptiveSymbol,
    };
    Type type;
    // Whether to use 'mappings' to map raw stored values to symbol values.
    bool use_mapping;
    uint16_t* mappings;
    // For type kRange, the symbol is stored as a value in [0;range[
    uint16_t range;
    union Param {
      // For type kGolomb.
      struct {
        // Range of the prefix.
        uint16_t range;
        // Size in bits of the Golomb prefix (1 or 2 for now).
        uint8_t prefix_size;
      } golomb;
      // For type kTrivial, there is only one value.
      int trivial_value;
      // For type kAdaptive where range is 2.
      // Index of this adaptive bit a_bits_.
      uint32_t a_bit_index;
      // For type kAdaptive where range is != 2.
      // Index of the adaptive symbol in a_symbols_.
      uint32_t a_symbol_index;
    };
    Param param;
    EXTRA extra;
  };
  // Gets the statistics for symbol 'sym' in cluster 'cluster'.
  inline Stat* GetStats(uint32_t cluster, uint32_t sym) {
    assert(sym < symbols_info_.Size() &&
           cluster < symbols_info_.NumClusters(sym));
    return &stats_start_[sym][cluster];
  }
  inline const Stat* GetStats(uint32_t cluster, uint32_t sym) const {
    assert(sym < symbols_info_.Size() &&
           cluster < symbols_info_.NumClusters(sym));
    return &stats_start_[sym][cluster];
  }

  // Helper function to set the stats of symbol 'sym' at 'cluster' to be proper
  // Golomb.
  void SetGolombStat(uint32_t cluster, uint32_t sym, uint32_t prefix_size) {
    Stat* const stat = GetStats(cluster, sym);
    stat->type = Stat::Type::kGolomb;
    stat->range = symbols_info_.Range(cluster, sym);
    const Golomb g = Golomb(stat->range - 1, prefix_size);
    stat->param.golomb.range = g.prefix + 1;
    stat->param.golomb.prefix_size = prefix_size;
  }

  // Sets up the given symbol as an adaptive bit.
  WP2Status AddAdaptiveBit(uint32_t cluster, uint32_t sym,
                           uint32_t p0 = 1, uint32_t p1 = 1) {
    Stat* const stat = GetStats(cluster, sym);
    stat->type = Stat::Type::kAdaptiveBit;
    stat->use_mapping = false;
    stat->range = 2;
    stat->param.a_bit_index = a_bits_.size();
    WP2_CHECK_ALLOC_OK(a_bits_.push_back(ANSBinSymbol(p0, p1)));
    return WP2_STATUS_OK;
  }

  // Sets up the given symbol as an adaptive symbol.
  WP2Status AddAdaptiveSymbol(uint32_t cluster, uint32_t sym,
                              ANSAdaptiveSymbol::Method method,
                              uint32_t speed) {
    Stat* const stat = GetStats(cluster, sym);
    stat->type = Stat::Type::kAdaptiveSymbol;
    stat->use_mapping = false;
    const uint32_t range = symbols_info_.Range(cluster, sym);
    stat->range = range;
    stat->param.a_symbol_index = a_symbols_.size();
    WP2_CHECK_ALLOC_OK(a_symbols_.Add(range));
    const uint16_t* cdf;
    uint16_t max_proba;
    symbols_info_.InitialCDF(cluster, sym, &cdf, &max_proba);
    auto& asym = a_symbols_[stat->param.a_symbol_index];
    if (cdf == nullptr) {
      // If no CDF has been set, assume a uniform distribution.
      asym.InitFromUniform(range);
    } else {
      WP2_CHECK_STATUS(asym.InitFromCDF(cdf, range, max_proba));
    }
    asym.SetAdaptationSpeed(method, speed);
    return WP2_STATUS_OK;
  }

  SymbolsInfo symbols_info_;

  ANSAdaptiveBits a_bits_;
  ANSAdaptiveSymbols a_symbols_;

  // For each symbol 'sym', the statistics for all clusters start at
  // stats_start_[sym].
  Stat* stats_start_[kSymbolNumMax];
  VectorNoCtor<Stat> all_stats_;
  Vector_u16 mappings_buffer_;
};

template <class T>
const uint16_t SymbolIO<T>::Stat::kInvalidMapping =
    std::numeric_limits<uint16_t>::max();

}  // namespace WP2

//------------------------------------------------------------------------------
// Info about the symbols.

namespace WP2 {

enum Symbol {
  kSymbolModeY = 0,  // prediction modes for Y
  kSymbolModeUV,     // prediction modes for U/V
  kSymbolModeA,      // prediction mode for A
  kSymbolSegmentId,
  kSymbolRndMtx0,
  // Symbol for the block transform (same in X and Y).
  kSymbolTransform,
  // Symbol for whether the block contains several transforms or not.
  kSymbolSplitTransform,
  // Symbol for BlockAlphaMode for each block.
  kSymbolBlockAlphaMode,
  // Symbol for whether the block uses 420 (downscaled UV).
  kSymbolYuv420,
  // Symbol for whether the block uses random matrices.
  kSymbolUseRandomMatrix,
  // Symbol for whether at least one coeff is not zero.
  kSymbolHasCoeffs,
  // Symbol for the residual encoding method for Y, U and V.
  kSymbolEncodingMethod,
  // Symbol for the residual encoding method for A.
  kSymbolEncodingMethodA,
  // Slope ('a' in 'chroma = a * luma + b') for chroma-from-luma prediction.
  kSymbolCflSlope,
  // Block size. One cluster per bounding BlockSize.
  kSymbolBlockSize,

  // ------------------
  // Residual encoding
  // ------------------
  // Number of zero residuals before the enc of blocK
  kSymbolResNumZeros,
  // Symbols for small residuals.
  kSymbolBits0,
  // Symbols for the prefix of big residuals.
  kSymbolBits1,
  // Symbols for DC residuals.
  kSymbolDC,
  // Symbol for the prefix of DC coefficients.
  kSymbolDCDict,

  // Whether we use 1 or 2 boundaries for the coefficients.
  kSymbolResidualUseBounds,
  // Whether the first used boundary is X.
  kSymbolResidualBound1IsX,
  // Whether we store the second boundary for the coefficients.
  kSymbolResidualUseBound2,
  // Whether the residual value is 0.
  kSymbolResidualIsZero,
  // Whether the residual value is 1.
  kSymbolResidualIsOne,
  // Whether the residual value is 2.
  kSymbolResidualIsTwo,
  // Whether all residuals after this one are 0.
  kSymbolResidualEndOfBlock,
  // Whether all residuals after this one are 1 or -1.
  kSymbolResidualHasOnlyOnesLeft,

  // ---------------------
  // AOM residual encoding
  // ---------------------
  kAOMAllZero,
  kAOMEOBPT16,
  kAOMEOBPT32,
  kAOMEOBPT64,
  kAOMEOBPT128,
  kAOMEOBPT256,
  kAOMEOBPT512,
  kAOMEOBPT1024,
  kAOMEOBExtra,
  kAOMCoeffBaseEOB,
  kAOMCoeffBase,
  kAOMCoeffBaseRange,
  kAOMDCSign,
  kSymbolNum
};

static_assert(kSymbolNum <= kSymbolNumMax, "kSymbolNumMax too small");

constexpr Symbol kSymbolsForCoeffs[] = {
    kSymbolResidualIsZero, kSymbolResidualIsOne, kSymbolResidualIsTwo,
    kSymbolResidualEndOfBlock, kSymbolResidualHasOnlyOnesLeft};
constexpr Symbol kSymbolsForCoeffsPerTf[] = {
    kSymbolDC, kSymbolResidualUseBounds, kSymbolResidualBound1IsX,
    kSymbolResidualUseBound2};
constexpr Symbol kSymbolsForAOMCoeffs[] = {
    kAOMAllZero,  kAOMEOBPT16,      kAOMEOBPT32,   kAOMEOBPT64,
    kAOMEOBPT128, kAOMEOBPT256,     kAOMEOBPT512,  kAOMEOBPT1024,
    kAOMEOBExtra, kAOMCoeffBaseEOB, kAOMCoeffBase, kAOMCoeffBaseRange,
    kAOMDCSign};

enum SymbolCount {
  kSymbolCountZero = 0,
  kSymbolCountOne,
  kSymbolCountMoreThanOne,
  kSymbolCountLast
};

enum AlphaMode {
  kAlphaModeLossless = 0,
  kAlphaModeLossy,
  kAlphaModeNum
};

enum BlockAlphaMode {
  kBlockAlphaFullTransp,  // All pixels have an alpha of 0.
  kBlockAlphaFullOpaque,  // All pixels have an alpha of kAlphaMax.
  kBlockAlphaLossy,       // Block alpha is encoded lossily.
  kBlockAlphaLossless,    // Block alpha is encoded losslessly.
};

enum class ChromaSubsampling {
  k420,
  k444,
  kAdaptive,     // Can differ from one block to the next.
  kSingleBlock,  // Do not signal anything before decoding the only block.
};

}  // namespace WP2

namespace WP2L {

enum Symbol {
  kSymbolType = 0,  // one of the SymbolType
  kSymbolA,
  kSymbolR,
  kSymbolG,
  kSymbolB,
  kSymbolCache,
  kSymbolDist,
  kSymbolLen,
  kSymbolNum
};
static const char* const kSymbolNames[] = {"type", "a",     "r",    "g",
                                           "b",    "cache", "dist", "len"};
STATIC_ASSERT_ARRAY_SIZE(kSymbolNames, kSymbolNum);

// Group 4 symbols.
enum SymbolGroup4 {
  kSymbolG4Type,
  kSymbolG4HorizontalDist,
  kSymbolG4VerticalDist,
  // Boolean to signal if the next color is different from the expected one
  // (as deduced from previous row).
  kSymbolG4ColorChange,
  // Index of the new color if kSymbolG4ColorChange was true.
  kSymbolG4NewColor,

  kSymbolG4Num
};
static const char* const kSymbolGroup4Names[] = {"type", "hdist", "vdist",
                                                 "color_change", "new_color"};
STATIC_ASSERT_ARRAY_SIZE(kSymbolGroup4Names, kSymbolG4Num);

// Initializes a SymbolsInfo to be used with Group 4 algorithm.
void InitGroup4(uint32_t width, uint32_t num_colors,
                WP2::SymbolsInfo* const info);

static_assert(kSymbolNum <= WP2::kSymbolNumMax, "kSymbolNumMax too small");

enum SymbolType {
  kSymbolTypeLiteral,   // literal value
  kSymbolTypeCopy,      // copy from (distance, length)
  kSymbolTypeCacheIdx,  // index in the color cache
  kSymbolTypeNum
};

// Class containing info about the size (in bits) of the different symbols
// (colors, cache index and LZ77 distance/length).
class LosslessSymbolsInfo : public WP2::SymbolsInfo {
 public:
  LosslessSymbolsInfo();
  void Init(bool has_alpha, WP2SampleFormat format, uint32_t num_clusters = 1);

  virtual ~LosslessSymbolsInfo() = default;

  WP2Status CopyFrom(const LosslessSymbolsInfo& other);

  // Changes the number of clusters for all symbols.
  void SetNumClusters(uint32_t num_clusters);

  void SetCacheRange(uint32_t cache_range);
  // Set the symbols A,R,B as useless as all the information will be stored in
  // G channel.
  void SetAsLabelImage();
  // Set the symbol A as useless and R,G,B as ranging with the cross color image
  // bounds kRedToBlueMax, kGreenToBlueMax, kGreenToRedMax
  void SetAsCrossColorImage();
  // TODO(vrabaud) ChannelBits is now used for wrap-around only but there are
  //               hard coded values that need to be removed for a full
  //               clean-up.
  uint32_t ChannelBits() const { return 8; }
  WP2SampleFormat SampleFormat() const { return sample_format_; }
  // Returns whether a symbol is useless to encode given the counts of its type.
  // E.g., if there is no literal, the A,R,G,B symbols are useless.
  bool IsSymbolUseless(const Symbol sym,
                       const bool is_maybe_used[kSymbolTypeNum]) const;
  // Returns whether alpha is used and whether it is a label image (ARB unused).
  static bool HasAlpha(const SymbolsInfo& info);
  static bool IsLabelImage(const SymbolsInfo& info);

 private:
  WP2SampleFormat sample_format_;
};

}  // namespace WP2L

//------------------------------------------------------------------------------

#endif  // WP2_COMMON_SYMBOLS_H_
