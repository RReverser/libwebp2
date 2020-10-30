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

#include "src/common/symbols.h"

#include "src/common/color_precision.h"
#include "src/common/lossy/block.h"
#include "src/common/lossy/block_size_io.h"
#include "src/common/lossy/predictor.h"
#include "src/common/lossy/transforms.h"
#include "src/dsp/lossless/lossless_common.h"

//------------------------------------------------------------------------------

namespace WP2 {

const SymbolsInfo::ClusterInfo SymbolsInfo::kDefaultClusterInfo = {};

bool SymbolsInfo::operator==(const SymbolsInfo& info) const {
  if (num_symbols_ != info.num_symbols_) return false;
  for (uint32_t sym = 0; sym < num_symbols_; ++sym) {
    if (!(GetSymbolInfo(sym) == info.GetSymbolInfo(sym))) return false;
    for (uint32_t cluster = 0; cluster < NumClusters(sym); ++cluster) {
      if (!(GetClusterInfo(cluster, sym) ==
            info.GetClusterInfo(cluster, sym))) {
        return false;
      }
    }
  }
  return true;
}

void SymbolsInfo::SetInfo(uint32_t sym, int32_t min, int32_t max,
                          uint32_t num_clusters, StorageMethod method,
                          bool probably_geometric) {
  assert(sym < kSymbolNumMax);
  // Used as a guarantee that Init() was called.
  assert(method == StorageMethod::kUnused || num_clusters > 0);
  num_symbols_ = std::max(num_symbols_, sym + 1);
  SymbolInfo* const info = GetSymbolInfo(sym);
  info->storage_method = method;
  info->num_clusters = num_clusters;
  info->probably_geometric = probably_geometric;
  SetMinMax(sym, min, max);
}

bool SymbolsInfo::SymbolInfo::MinMaxIsValid(int32_t new_min,
                                            int32_t new_max) const {
  if (new_min == kInvalidBound && new_max == kInvalidBound) {
    return true;
  }
  assert(new_min <= new_max);
  // For now, only support positive and symbols equiprobable with respect to
  // sign.
  assert(new_min == 0 || new_min == -new_max);
  const int32_t new_range = new_max - new_min + 1;
  assert(new_range >= 0 && new_range < (1 << kANSMaxRangeBits));
  switch (storage_method) {
    case SymbolsInfo::StorageMethod::kAuto:
      return ((uint32_t)new_range <= kANSMaxRange);
    case SymbolsInfo::StorageMethod::kAdaptiveBit:
      return (new_range == 2);
    case SymbolsInfo::StorageMethod::kAdaptiveSym:
    case SymbolsInfo::StorageMethod::kAdaptiveWithAutoSpeed:
      return (new_range <= (int32_t)APROBA_MAX_SYMBOL);
    case SymbolsInfo::StorageMethod::kUnused:
      return true;
  }
  assert(false);
  return true;
}

bool SymbolsInfo::SymbolInfo::operator==(
    const SymbolsInfo::SymbolInfo& rhs) const {
  return num_clusters == rhs.num_clusters &&
         storage_method == rhs.storage_method && min == rhs.min &&
         max == rhs.max;
}
bool SymbolsInfo::ClusterInfo::operator==(
    const SymbolsInfo::ClusterInfo& rhs) const {
  return (p0 == rhs.p0 && p1 == rhs.p1 && min == rhs.min && max == rhs.max &&
          cdf == rhs.cdf);
}

WP2Status SymbolsInfo::CopyFrom(const SymbolsInfo& other) {
  num_symbols_ = other.num_symbols_;
  symbols_info_ = other.symbols_info_;
  WP2_CHECK_ALLOC_OK(cluster_info_.resize(other.cluster_info_.size()));
  std::copy(other.cluster_info_.begin(), other.cluster_info_.end(),
            cluster_info_.begin());
  next_cluster_info_index_ = other.next_cluster_info_index_;
  return WP2_STATUS_OK;
}

WP2Status SymbolsInfo::SetMinMax(uint32_t cluster, uint32_t sym, int32_t min,
                                 int32_t max) {
  assert(sym < kSymbolNumMax);
  assert(cluster < GetSymbolInfo(sym)->num_clusters);
  assert(GetSymbolInfo(sym)->MinMaxIsValid(min, max));
  ClusterInfo* const info = GetOrMakeClusterInfo(cluster, sym);
  WP2_CHECK_ALLOC_OK(info != nullptr);
  info->min = min;
  info->max = max;
  return WP2_STATUS_OK;
}

void SymbolsInfo::SetMinMax(uint32_t sym, int32_t min, int32_t max) {
  SymbolInfo* const info = GetSymbolInfo(sym);
  // Check Init() was called.
  assert(info->storage_method == StorageMethod::kUnused ||
         info->num_clusters > 0);
  assert(info->MinMaxIsValid(min, max));
  info->min = min;
  info->max = max;
  for (uint32_t c = 0; c < info->num_clusters; ++c) {
    // Should not set a global range after setting per-cluster range.
    assert(GetClusterInfo(c, sym).min == kInvalidBound &&
           GetClusterInfo(c, sym).max == kInvalidBound);
  }
}

void SymbolsInfo::SetClusters(uint32_t sym, uint32_t num_clusters) {
  assert(sym < kSymbolNumMax);
  // Check Init() was called.
  assert(Method(sym) == StorageMethod::kUnused ||
         GetSymbolInfo(sym)->num_clusters > 0);
  assert(num_clusters > 0);
  GetSymbolInfo(sym)->num_clusters = num_clusters;
}

WP2Status SymbolsInfo::SetStartingProba(uint32_t cluster, uint32_t sym,
                                        uint32_t p0, uint32_t p1) {
  assert(p0 > 0);
  assert(p1 > 0);
  assert(Method(sym) == StorageMethod::kAdaptiveBit);
  ClusterInfo* const info = GetOrMakeClusterInfo(cluster, sym);
  WP2_CHECK_ALLOC_OK(info != nullptr);
  info->p0 = p0;
  info->p1 = p1;
  return WP2_STATUS_OK;
}

SymbolsInfo::ClusterInfo* SymbolsInfo::GetOrMakeClusterInfo(uint32_t cluster,
                                                            uint32_t sym) {
  SymbolInfo* const info = GetSymbolInfo(sym);
  assert(cluster < info->num_clusters);
  if (info->cluster_info_index == SymbolInfo::kNoPerClusterInfo) {
    info->cluster_info_index = next_cluster_info_index_;
    next_cluster_info_index_ += info->num_clusters;
    if (!cluster_info_.resize(next_cluster_info_index_)) {
      return nullptr;
    }
  }
  return &cluster_info_[info->cluster_info_index + cluster];
}

const SymbolsInfo::ClusterInfo& SymbolsInfo::GetClusterInfo(
    uint32_t cluster, uint32_t sym) const {
  assert(cluster < NumClusters(sym));
  const SymbolInfo& info = GetSymbolInfo(sym);
  if (info.cluster_info_index == SymbolInfo::kNoPerClusterInfo) {
    return kDefaultClusterInfo;
  }
  return cluster_info_[info.cluster_info_index + cluster];
}

//------------------------------------------------------------------------------

WP2Status SymbolsInfo::InitLossy(uint32_t num_segments,
                                 PartitionSet partition_set, bool has_alpha,
                                 uint32_t quality_hint, bool use_aom_coeffs) {
  const uint32_t num_channels = (has_alpha ? 4 : 3);

  // kNum is used as an escape code when the upper/left mode do not match.
  SetInfo(kSymbolModeY, /*min=*/0, /*max=*/kYPredModeNum - 1,
          (uint32_t)Predictor::ModeContext::kNum + 1,
          StorageMethod::kAdaptiveSym);
  SetInfo(kSymbolModeUV, /*min=*/0, /*max=*/kUVPredModeNum - 1,
          /*num_clusters=*/1, StorageMethod::kAdaptiveWithAutoSpeed);
  SetInfo(kSymbolSegmentId, /*min=*/0, /*max=*/num_segments - 1,
          /*num_clusters=*/1, StorageMethod::kAdaptiveWithAutoSpeed,
          /*probably_geometric=*/true);
  SetInfo(kSymbolResNumZeros, /*min=*/0, /*max=*/kResidualCons0Max,
          /*num_clusters=*/kNumResidualStats * num_channels,
          StorageMethod::kAuto);
  SetInfo(kSymbolBits0, /*min=*/0, /*max=*/kResidual1Max,
          /*num_clusters=*/kNumResidualStats * num_channels,
          StorageMethod::kAuto);
  SetInfo(kSymbolBits1, /*min=*/0, /*max=*/kMaxCoeffBits,
          /*num_clusters=*/kNumResidualStats * num_channels,
          StorageMethod::kAuto);
  // 1 bit is added by the fact that we subtract the predicted value from
  // the actual value, and 1 bit is added by the 2D transform.
  // Finally, one bit is added for the sign.
  SetInfo(kSymbolDC, /*min=*/0, /*max=*/2 * kMaxDcValue,
          /*num_clusters=*/num_channels, StorageMethod::kAuto,
          /*probably_geometric=*/true);
  SetInfo(kSymbolRndMtx0, /*min=*/0, /*max=*/kMaxNumRndMtx - 1,
          /*num_clusters=*/1, StorageMethod::kAuto);

  SetInfo(kSymbolTransform, /*min=*/0, /*max=*/kNumTransformPairs - 1,
          /*num_clusters=*/1, StorageMethod::kAdaptiveWithAutoSpeed);
  SetInfo(kSymbolBlockAlphaMode, /*min=*/0, /*max=*/2,
          /*num_clusters=*/1, StorageMethod::kAuto);
  SetInfo(kSymbolSplitTransform, /*min=*/0, /*max=*/1, /*num_clusters=*/1,
          StorageMethod::kAdaptiveWithAutoSpeed);
  SetInfo(kSymbolYuv420, /*min=*/0, /*max=*/1, /*num_clusters=*/1,
          StorageMethod::kAdaptiveWithAutoSpeed);
  SetInfo(kSymbolUseRandomMatrix, /*min=*/0, /*max=*/1,
          /*num_clusters=*/1, StorageMethod::kAdaptiveBit);

  // Minus 1 bit because it doesn't include the sign.
  SetInfo(kSymbolCflSlope, /*min=*/0,
          /*max=*/(1 << (SignalingCflPredictor::kAResBits - 1)) - 1,
          /*num_clusters=*/1, StorageMethod::kAuto,
          /*probably_geometric=*/true);

  SetInfo(kSymbolHasCoeffs, /*min=*/0, /*max=*/1, /*num_clusters=*/num_channels,
          StorageMethod::kAdaptiveWithAutoSpeed);
  if (has_alpha) {
    SetInfo(kSymbolModeA, /*min=*/0, /*max=*/kAPredModeNum - 1,
            /*num_clusters=*/1, StorageMethod::kAdaptiveSym);
    SetInfo(kSymbolEncodingMethodA, /*min=*/0,
            /*max=*/(int)EncodingMethod::kNumMethod - 1, /*num_clusters=*/1,
            StorageMethod::kAdaptiveWithAutoSpeed);
  } else {
    SetInfo(kSymbolModeA, /*min=*/SymbolsInfo::kInvalidBound,
            /*max=*/SymbolsInfo::kInvalidBound, /*num_clusters=*/1,
            StorageMethod::kUnused);
    SetInfo(kSymbolEncodingMethodA, /*min=*/SymbolsInfo::kInvalidBound,
            /*max=*/SymbolsInfo::kInvalidBound, /*num_clusters=*/1,
            StorageMethod::kUnused);
  }
  SetInfo(kSymbolEncodingMethod, /*min=*/0,
          /*max=*/(int32_t)EncodingMethod::kNumMethod - 1,
          /*num_clusters=*/3 * num_segments,
          StorageMethod::kAdaptiveWithAutoSpeed);

  const uint32_t max_num_sizes = GetNumBlockSizes(partition_set, BLK_32x32);
  SetInfo(kSymbolBlockSize, /*min=*/0, /*max=*/max_num_sizes-1,
          /*num_clusters=*/GetNumUniqueBounds(partition_set),
          (max_num_sizes > APROBA_MAX_SYMBOL) ? StorageMethod::kAuto
                                              : StorageMethod::kAdaptiveSym);
  for (uint32_t i = 0; i < GetNumUniqueBounds(partition_set); ++i) {
    const BlockSize bounds = GetUniqueBounds(partition_set, i);
    const uint32_t num_sizes = GetNumBlockSizes(partition_set, bounds);
    assert(num_sizes <= max_num_sizes);
    WP2_CHECK_STATUS(SetMinMax(/*cluster=*/i, kSymbolBlockSize,
                               /*min=*/0, /*max=*/num_sizes - 1));
  }

  SetInfo(kSymbolResidualUseBounds, /*min=*/0, /*max=*/1,
          /*num_clusters=*/kResidualClusters, StorageMethod::kAdaptiveBit);
  SetInfo(kSymbolResidualBound1IsX, /*min=*/0, /*max=*/1,
          /*num_clusters=*/kResidualClusters, StorageMethod::kAdaptiveBit);
  SetInfo(kSymbolResidualUseBound2, /*min=*/0, /*max=*/1,
          /*num_clusters=*/kResidualClusters, StorageMethod::kAdaptiveBit);

  // Non-AOM coeffs are always used for UV so they need to be initialized.
  for (uint32_t sym : kSymbolsForCoeffs) {
    // Y/UV, per residual method, per block size, per sector.
    const uint32_t num_clusters = kResidualClusters * ResidualIO::kNumSectors;
    SetInfo(sym, /*min=*/0, /*max=*/1, num_clusters,
            StorageMethod::kAdaptiveBit);
  }

  // Initialize coefficient probabilities.
  WP2_CHECK_STATUS(ResidualIO::InitializeInfo(has_alpha, quality_hint,
                                              use_aom_coeffs, this));
  return WP2_STATUS_OK;
}

//------------------------------------------------------------------------------

void WriteGolomb(uint32_t value, uint32_t min, uint32_t max,
                 uint32_t prefix_size, ANSEncBase* const enc, WP2_OPT_LABEL) {
  ANSDebugPrefix ans_prefix(enc, label);

  const Golomb golomb(value - min, prefix_size);
  uint32_t range = max - min + 1;
  const Golomb max_golomb(range - 1, prefix_size);

  enc->PutRValue(golomb.prefix, max_golomb.prefix + 1, "golomb");
  if (golomb.extra_bits_num > 0) {
    // If the last bit interval is truncated, go with ranges.
    if (golomb.prefix == max_golomb.prefix) {
      assert(Golomb::Merge(golomb.prefix, prefix_size, 0) < range);
      range -= Golomb::Merge(golomb.prefix, prefix_size, 0);
      enc->PutRValue(golomb.extra_bits_value, range, "extra_bits_value");
    } else {
      enc->PutUValue(golomb.extra_bits_value, golomb.extra_bits_num,
                     "extra_bits_value");
    }
  }
}

uint32_t ReadGolomb(uint32_t min, uint32_t max, uint32_t prefix_size,
                    ANSDec* const dec, WP2_OPT_LABEL) {
  ANSDebugPrefix ans_prefix(dec, label);

  uint32_t range = max - min + 1;
  const Golomb max_golomb(range - 1, prefix_size);

  const uint32_t prefix = dec->ReadRValue(max_golomb.prefix + 1, "golomb");
  const uint32_t extra_bits = Golomb::NumExtraBits(prefix, prefix_size);
  uint32_t extra_bits_value;
  if (extra_bits == 0) {
    extra_bits_value = 0;
  } else {
    // Use ranges if we are at the last interval.
    if (prefix == max_golomb.prefix) {
      assert(Golomb::Merge(prefix, prefix_size, 0) < range);
      range -= Golomb::Merge(prefix, prefix_size, 0);
      extra_bits_value = dec->ReadRValue(range, "extra_bits_value");
    } else {
      extra_bits_value = dec->ReadUValue(extra_bits, "extra_bits_value");
    }
  }
  return min + Golomb::Merge(prefix, prefix_size, extra_bits_value);
}

}  // namespace WP2

//------------------------------------------------------------------------------

namespace WP2L {

void InitGroup4(uint32_t width, uint32_t num_colors,
                WP2::SymbolsInfo* const info) {
  constexpr WP2::SymbolsInfo::StorageMethod storage_method =
      WP2::SymbolsInfo::StorageMethod::kAuto;
  info->SetInfo(kSymbolG4Type, /*min=*/0, /*max=*/2, /*num_clusters=*/1,
                storage_method);
  info->SetInfo(kSymbolG4HorizontalDist, /*min=*/0, /*max=*/width,
                /*num_clusters=*/2, storage_method);
  info->SetInfo(kSymbolG4VerticalDist, /*min=*/0, /*max=*/2 * kGroup4Window,
                /*num_clusters=*/1, storage_method);
  if (num_colors > 2) {
    info->SetInfo(kSymbolG4ColorChange, /*min=*/0, /*max=*/1,
                  /*num_clusters=*/1, storage_method);
    // Minus 2 because the new color cannot be the same as the current color or
    // the expected color.
    info->SetInfo(kSymbolG4NewColor, /*min=*/0, /*max=*/num_colors - 3,
                  /*num_clusters=*/1, storage_method);
  }
}

LosslessSymbolsInfo::LosslessSymbolsInfo() : sample_format_(WP2_FORMAT_NUM) {
  for (uint32_t i = 0; i < kSymbolNum; ++i) {
    SetInfo(i, /*min=*/SymbolsInfo::kInvalidBound,
            /*max=*/SymbolsInfo::kInvalidBound, /*num_clusters=*/1,
            StorageMethod::kAuto);
  }
}

void LosslessSymbolsInfo::Init(bool has_alpha, WP2SampleFormat sample_format,
                               uint32_t num_clusters) {
  sample_format_ = sample_format;
  constexpr StorageMethod storage_method = StorageMethod::kAuto;
  int extra_bits, code;
  SetInfo(kSymbolType, /*min=*/0, /*max=*/kSymbolTypeNum - 1, num_clusters,
          storage_method);
  SetInfo(kSymbolA, /*min=*/0, /*max=*/WP2::FormatMax(sample_format, 0),
          num_clusters, storage_method);
  SetInfo(kSymbolR, /*min=*/0, /*max=*/WP2::FormatMax(sample_format, 1),
          num_clusters, storage_method);
  SetInfo(kSymbolG, /*min=*/0, /*max=*/WP2::FormatMax(sample_format, 2),
          num_clusters, storage_method);
  SetInfo(kSymbolB, /*min=*/0, /*max=*/WP2::FormatMax(sample_format, 3),
          num_clusters, storage_method);
  WP2LPrefixEncodeBits(kMaxLZ77Distance + kCodeToPlaneCodes, &code,
                       &extra_bits);
  SetInfo(kSymbolDist, /*min=*/0, /*max=*/code, num_clusters, storage_method);
  WP2LPrefixEncodeBits(kMaxLZ77Length, &code, &extra_bits);
  SetInfo(kSymbolLen, /*min=*/0, /*max=*/code, num_clusters, storage_method);

  SetCacheRange(0);
  if (!has_alpha) {
    SetInfo(kSymbolA, /*min=*/SymbolsInfo::kInvalidBound,
            /*max=*/SymbolsInfo::kInvalidBound, /*num_clusters=*/0,
            StorageMethod::kUnused);
  }
}

WP2Status LosslessSymbolsInfo::CopyFrom(const LosslessSymbolsInfo& other) {
  WP2_CHECK_STATUS(SymbolsInfo::CopyFrom(other));
  sample_format_ = other.sample_format_;
  return WP2_STATUS_OK;
}

void LosslessSymbolsInfo::SetNumClusters(uint32_t num_clusters) {
  for (uint32_t sym = 0; sym < Size(); ++sym) SetClusters(sym, num_clusters);
}

void LosslessSymbolsInfo::SetCacheRange(uint32_t cache_range) {
  const uint32_t num_clusters = NumClusters(kSymbolType);
  if (cache_range == 0) {
    // Do not consider kSymbolTypeCacheIdx if the cache is unused.
    SetInfo(kSymbolType, /*min=*/0, /*max=*/kSymbolTypeNum - 2, num_clusters,
            StorageMethod::kAuto);
    SetInfo(kSymbolCache, /*min=*/SymbolsInfo::kInvalidBound,
            /*max=*/SymbolsInfo::kInvalidBound, /*num_clusters=*/0,
            StorageMethod::kUnused);
  } else {
    SetInfo(kSymbolType, /*min=*/0, /*max=*/kSymbolTypeNum - 1, num_clusters,
            StorageMethod::kAuto);
    SetInfo(kSymbolCache, /*min=*/0, /*max=*/cache_range - 1, num_clusters,
            StorageMethod::kAuto);
  }
}

void LosslessSymbolsInfo::SetAsLabelImage() {
  SetInfo(kSymbolA, /*min=*/SymbolsInfo::kInvalidBound,
          /*max=*/SymbolsInfo::kInvalidBound, /*num_clusters=*/0,
          StorageMethod::kUnused);
  SetInfo(kSymbolR, /*min=*/SymbolsInfo::kInvalidBound,
          /*max=*/SymbolsInfo::kInvalidBound, /*num_clusters=*/0,
          StorageMethod::kUnused);
  SetInfo(kSymbolB, /*min=*/SymbolsInfo::kInvalidBound,
          /*max=*/SymbolsInfo::kInvalidBound, /*num_clusters=*/0,
          StorageMethod::kUnused);
}

void LosslessSymbolsInfo::SetAsCrossColorImage() {
  SetInfo(kSymbolA, /*min=*/SymbolsInfo::kInvalidBound,
          /*max=*/SymbolsInfo::kInvalidBound, /*num_clusters=*/0,
          StorageMethod::kUnused);
  SetMinMax(kSymbolR, /*min=*/-(int32_t)kRedToBlueMax, /*max=*/kRedToBlueMax);
  SetMinMax(kSymbolG, /*min=*/-(int32_t)kGreenToBlueMax,
            /*max=*/kGreenToBlueMax);
  SetMinMax(kSymbolB, /*min=*/-(int32_t)kGreenToRedMax, /*max=*/kGreenToRedMax);
}

bool LosslessSymbolsInfo::IsSymbolUseless(
    const Symbol sym, const bool is_maybe_used[kSymbolTypeNum]) const {
  if (Method(sym) == StorageMethod::kUnused) return true;
  switch (sym) {
    case kSymbolType:
      // We always have symbols so the type is always useful.
      return false;
    case kSymbolA:
    case kSymbolR:
    case kSymbolG:
    case kSymbolB:
      return !is_maybe_used[kSymbolTypeLiteral];
    case kSymbolCache:
      return !is_maybe_used[kSymbolTypeCacheIdx];
    case kSymbolDist:
    case kSymbolLen:
      return !is_maybe_used[kSymbolTypeCopy];
    default:
      assert(false);
      return false;
  }
}

bool LosslessSymbolsInfo::HasAlpha(const SymbolsInfo& info) {
  return (info.Method(kSymbolA) != StorageMethod::kUnused);
}

bool LosslessSymbolsInfo::IsLabelImage(const SymbolsInfo& info) {
  const bool is_label_image = (info.Method(kSymbolR) == StorageMethod::kUnused);
  if (is_label_image) {
    assert(info.Method(kSymbolA) == StorageMethod::kUnused &&
           info.Method(kSymbolB) == StorageMethod::kUnused);
  }
  return is_label_image;
}

}  // namespace WP2L

//------------------------------------------------------------------------------
