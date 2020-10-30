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
//   Common code for AOM residuals, branched from libgav1, with as little
//   changes as possible.
//
// Author: Vincent Rabaud (vrabaud@google.com)

#ifndef WP2_COMMON_LOSSY_RESIDUALS_AOM_H_
#define WP2_COMMON_LOSSY_RESIDUALS_AOM_H_

#include <array>

#include "src/common/lossy/aom/array_2d.h"
#include "src/common/lossy/residuals.h"
#include "src/common/symbols.h"
#include "src/wp2/base.h"

namespace WP2 {
namespace libgav1 {

//------------------------------------------------------------------------------
// From libgav1/src/dsp/constants.h

enum {
  kCdfMaxProbability = 32768,
  kMaxSegments = 8,
  kQuantizedCoefficientBufferPadding = 4,
};  // anonymous enum

enum Plane : uint8_t { kPlaneY, kPlaneU, kPlaneV, kPlaneA };
enum : uint8_t { kMaxPlanesMonochrome = kPlaneY + 1, kMaxPlanes = kPlaneA + 1 };

// The plane types, called luma and chroma in the spec.
// Matches PlaneType in cdfs.inc
// TODO(maryla): should there be one separate for alpha?
enum PlaneType : uint8_t { kPlaneTypeY, kPlaneTypeUV, kNumPlaneTypes };

enum BlockSize : uint8_t {
  kBlock4x4,
  kBlock4x8,
  kBlock4x16,
  kBlock8x4,
  kBlock8x8,
  kBlock8x16,
  kBlock8x32,
  kBlock16x4,
  kBlock16x8,
  kBlock16x16,
  kBlock16x32,
  kBlock16x64,
  kBlock32x8,
  kBlock32x16,
  kBlock32x32,
  kBlock32x64,
  kBlock64x16,
  kBlock64x32,
  kBlock64x64,
  kBlock64x128,
  kBlock128x64,
  kBlock128x128,
  kMaxBlockSizes,
  kBlockInvalid
};

enum TransformSize : uint8_t {
  kTransformSize4x4,
  kTransformSize4x8,
  kTransformSize4x16,
  kTransformSize8x4,
  kTransformSize8x8,
  kTransformSize8x16,
  kTransformSize8x32,
  kTransformSize16x4,
  kTransformSize16x8,
  kTransformSize16x16,
  kTransformSize16x32,
  kTransformSize16x64,
  kTransformSize32x8,
  kTransformSize32x16,
  kTransformSize32x32,
  kTransformSize32x64,
  kTransformSize64x16,
  kTransformSize64x32,
  kTransformSize64x64,
  kNumTransformSizes
};

extern const uint8_t kBlockWidthPixels[kMaxBlockSizes];

extern const uint8_t kBlockHeightPixels[kMaxBlockSizes];

extern const uint8_t kTransformWidth[kNumTransformSizes];

extern const uint8_t kTransformHeight[kNumTransformSizes];

extern const uint8_t kTransformWidth4x4[kNumTransformSizes];

extern const uint8_t kTransformHeight4x4[kNumTransformSizes];

extern const uint8_t kTransformWidthLog2[kNumTransformSizes];

// Replaces all occurrences of 64x* and *x64 with 32x* and *x32 respectively.
constexpr TransformSize kAdjustedTransformSize[kNumTransformSizes] = {
    kTransformSize4x4,   kTransformSize4x8,   kTransformSize4x16,
    kTransformSize8x4,   kTransformSize8x8,   kTransformSize8x16,
    kTransformSize8x32,  kTransformSize16x4,  kTransformSize16x8,
    kTransformSize16x16, kTransformSize16x32, kTransformSize16x32,
    kTransformSize32x8,  kTransformSize32x16, kTransformSize32x32,
    kTransformSize32x32, kTransformSize32x16, kTransformSize32x32,
    kTransformSize32x32};

//------------------------------------------------------------------------------
// From libgav1/src/utils/common.h

constexpr int DivideBy2(int n) { return n >> 1; }
constexpr int DivideBy4(int n) { return n >> 2; }
constexpr int DivideBy8(int n) { return n >> 3; }

// Convert |value| to unsigned before shifting to avoid undefined behavior with
// negative values.
inline int LeftShift(int value, int bits) {
  assert(bits >= 0);
  assert(value >= -(int64_t{1} << (31 - bits)));
  assert(value <= (int64_t{1} << (31 - bits)) - ((bits == 0) ? 1 : 0));
  return static_cast<int>(static_cast<uint32_t>(value) << bits);
}
inline int MultiplyBy2(int n) { return LeftShift(n, 1); }
inline int MultiplyBy4(int n) { return LeftShift(n, 2); }

constexpr PlaneType GetPlaneType(Plane plane) {
  // Alpha uses kPlaneTypeY.
  return (plane == kPlaneU || plane == kPlaneV) ? kPlaneTypeUV : kPlaneTypeY;
}

// Input: 1d array index |index|, which indexes into a 2d array of width
//     1 << |tx_width_log2|.
// Output: 1d array index which indexes into a 2d array of width
//     (1 << |tx_width_log2|) + kQuantizedCoefficientBufferPadding.
inline int PaddedIndex(int index, int tx_width_log2) {
  return index + MultiplyBy4(index >> tx_width_log2);
}

constexpr int kNumQuantizerBaseLevels = 2;

//------------------------------------------------------------------------------
// From libgav1/src/symbol_decoder_context.h

enum {
  kBooleanFieldCdfSize = 3,
  kCoefficientQuantizerContexts = 4,
  kNumSquareTransformSizes = 5,
  kAllZeroContexts = 13,
  kEobPtContexts = 2,
  kEobPt16SymbolCount = 5,
  kEobPt32SymbolCount = 6,
  kEobPt64SymbolCount = 7,
  kEobPt128SymbolCount = 8,
  kEobPt256SymbolCount = 9,
  kEobPt512SymbolCount = 10,
  kEobPt1024SymbolCount = 11,
  kEobExtraContexts = 9,
  kCoeffBaseEobContexts = 4,
  kCoeffBaseEobSymbolCount = 3,
  kCoeffBaseContexts = 42,
  kCoeffBaseSymbolCount = 4,
  kCoeffBaseRangeContexts = 21,
  kCoeffBaseRangeSymbolCount = 4,
  kDcSignContexts = 3,
};  // anonymous enum

struct SymbolDecoderContext {
  SymbolDecoderContext() = default;

  // 'quantizer_context' is 0 for highest qualities, and 3 for lowest qualities.
  void Initialize(int quantizer_context);

  uint16_t all_zero_cdf[kNumSquareTransformSizes][kAllZeroContexts]
                       [kBooleanFieldCdfSize];
  uint16_t eob_pt_16_cdf[kNumPlaneTypes][kEobPtContexts]
                        [kEobPt16SymbolCount + 1];
  uint16_t eob_pt_32_cdf[kNumPlaneTypes][kEobPtContexts]
                        [kEobPt32SymbolCount + 1];
  uint16_t eob_pt_64_cdf[kNumPlaneTypes][kEobPtContexts]
                        [kEobPt64SymbolCount + 1];
  uint16_t eob_pt_128_cdf[kNumPlaneTypes][kEobPtContexts]
                         [kEobPt128SymbolCount + 1];
  uint16_t eob_pt_256_cdf[kNumPlaneTypes][kEobPtContexts]
                         [kEobPt256SymbolCount + 1];
  uint16_t eob_pt_512_cdf[kNumPlaneTypes][kEobPt512SymbolCount + 1];
  uint16_t eob_pt_1024_cdf[kNumPlaneTypes][kEobPt1024SymbolCount + 1];
  uint16_t eob_extra_cdf[kNumSquareTransformSizes][kNumPlaneTypes]
                        [kEobExtraContexts][kBooleanFieldCdfSize];
  uint16_t coeff_base_eob_cdf[kNumSquareTransformSizes][kNumPlaneTypes]
                             [kCoeffBaseEobContexts]
                             [kCoeffBaseEobSymbolCount + 1];
  uint16_t coeff_base_cdf[kNumSquareTransformSizes][kNumPlaneTypes]
                         [kCoeffBaseContexts][kCoeffBaseSymbolCount + 1];
  uint16_t coeff_base_range_cdf[kNumSquareTransformSizes][kNumPlaneTypes]
                               [kCoeffBaseRangeContexts]
                               [kCoeffBaseRangeSymbolCount + 1];
  uint16_t dc_sign_cdf[kNumPlaneTypes][kDcSignContexts][kBooleanFieldCdfSize];
};

// This is computed as:
// min(transform_width_log2, 5) + min(transform_height_log2, 5) - 4.
constexpr uint8_t kEobMultiSizeLookup[kNumTransformSizes] = {
    0, 1, 2, 1, 2, 3, 4, 2, 3, 4, 5, 5, 4, 5, 6, 6, 5, 6, 6};

class AOMContext;

typedef int (*CoeffBaseContextFunc)(
    const uint8_t qlevels[], TransformSize tx_size,
    int adjusted_tx_width_log2, uint16_t pos);
typedef int (*CoeffBaseRangeContextFunc)(
    const uint8_t qlevels[], int adjusted_tx_width_log2, int pos);

// Class specialized in reading transform residuals.
class AOMResidualIO {
 public:
  // Returns the cluster for a given symbol, given certain other parameters.
  // The overload GetClusterAOM1, GetClusterAOM2, GetClusterAOM3 are meant
  // for symbols using 1, 2 or 3 arguments.
  template <Symbol sym>
  static uint32_t GetClusterAOM1(uint32_t arg0) {
    static_assert(sym == kAOMEOBPT512 || sym == kAOMEOBPT1024, "wrong symbol");
    return arg0;
  }
  template <Symbol sym>
  static uint32_t GetClusterAOM2(uint32_t arg0, uint32_t arg1) {
    static_assert(sym == kAOMAllZero || sym == kAOMEOBPT16 ||
                      sym == kAOMEOBPT32 || sym == kAOMEOBPT64 ||
                      sym == kAOMEOBPT128 || sym == kAOMEOBPT256 ||
                      sym == kAOMDCSign,
                  "wrong symbol");
    const std::array<uint32_t, 2> extent =
        (sym == kAOMAllZero)
            ? std::array<uint32_t, 2>{kNumSquareTransformSizes,
                                      kAllZeroContexts}
            : (sym == kAOMDCSign)
                  ? std::array<uint32_t, 2>{kNumPlaneTypes, kDcSignContexts}
                  : std::array<uint32_t, 2>{kNumPlaneTypes, kEobPtContexts};
    uint32_t cluster = arg1;
    uint32_t stride = extent[1];
    cluster += arg0 * stride;
    return cluster;
  }
  template <Symbol sym>
  static uint32_t GetClusterAOM3(uint32_t arg0, uint32_t arg1, uint32_t arg2) {
    static_assert(sym == kAOMEOBExtra || sym == kAOMCoeffBaseEOB ||
                      sym == kAOMCoeffBase || sym == kAOMCoeffBaseRange,
                  "wrong symbol");
    const std::array<uint32_t, 3> extent =
        (sym == kAOMEOBExtra)
            ? std::array<uint32_t, 3>{kNumSquareTransformSizes, kNumPlaneTypes,
                                      kEobExtraContexts}
            : (sym == kAOMCoeffBaseEOB)
                  ? std::array<uint32_t, 3>{kNumSquareTransformSizes,
                                            kNumPlaneTypes,
                                            kCoeffBaseEobContexts}
                  : (sym == kAOMCoeffBase)
                        ? std::array<uint32_t, 3>{kNumSquareTransformSizes,
                                                  kNumPlaneTypes,
                                                  kCoeffBaseContexts}
                        : std::array<uint32_t, 3>{kNumSquareTransformSizes,
                                                  kNumPlaneTypes,
                                                  kCoeffBaseRangeContexts};
    uint32_t cluster = arg2;
    uint32_t stride = extent[2];
    cluster += arg1 * stride;
    stride *= extent[1];
    cluster += arg0 * stride;
    return cluster;
  }

 public:
  // Lookup used to call the right variant of GetCoeffBase*() based on
  // the transform class (in [0..2])
  static CoeffBaseContextFunc
      GetCoeffsBaseContextFunc(uint32_t tx_class);
  static CoeffBaseRangeContextFunc
      GetCoeffBaseRangeContextFunc(uint32_t tx_class);

 protected:
  // The residual pointer is used to traverse the |residual_buffer_|. It is
  // used in two different ways.
  // If |split_parse_and_decode_| is true:
  //    The pointer points to the beginning of the |residual_buffer_| when the
  //    "parse" and "decode" steps begin. It is then moved forward tx_size in
  //    each iteration of the "parse" and the "decode" steps. In this case, the
  //    ResidualPtr variable passed into various functions starting from
  //    ProcessSuperBlock is used as an in/out parameter to keep track of the
  //    residual pointer.
  // If |split_parse_and_decode_| is false:
  //    The pointer is reset to the beginning of the |residual_buffer_| for
  //    every transform block.
  using ResidualPtr = uint8_t*;

  // |x4| and |y4| are the column and row positions of the 4x4 block. |w4| and
  // |h4| are the width and height in 4x4 units of |tx_size|.
  static int GetTransformAllZeroContext(const AOMContext& aom_context,
                                        Plane plane, TransformSize tx_size,
                                        int x4, int y4, int w4, int h4);

  static int GetCoeffBaseContextEob(TransformSize tx_size, int index);

  static int GetCoeffBaseRangeContextEob(int adjusted_tx_width_log2, int pos);
  static int GetDcSignContext(const AOMContext& aom_context, int x4, int y4,
                              int w4, int h4, Plane plane);
  static void SetEntropyContexts(int x4, int y4, int w4, int h4, Plane plane,
                                 uint8_t coefficient_level, int8_t dc_category,
                                 AOMContext* const aom_context);
};

// Class containing the updated context from one set of residuals to the next.
class AOMContext {
 public:
  WP2Status Init(uint32_t width, uint32_t height);
  WP2Status CopyFrom(const AOMContext& other);
  void Reset();

  uint32_t Rows4x4(Plane plane) const {
    return rows4x4_ >> (IsReduced(plane) ? 1 : 0);
  }
  uint32_t Columns4x4(Plane plane) const {
    return columns4x4_ >> (IsReduced(plane) ? 1 : 0);
  }
  bool Is420() const { return is420_; }
  void SetIs420(bool is420) { is420_ = is420; }
  bool IsReduced(const Plane& plane) const {
    return ((plane == kPlaneU || plane == kPlaneV) && is420_);
  }

  // These two arrays (|coefficient_levels_| and |dc_categories_|) are used to
  // store the entropy context. Their dimensions are as follows: First -
  // left/top; Second - plane; Third - row4x4 (if first dimension is
  // left)/column4x4 (if first dimension is top).
  //
  // This is equivalent to the LeftLevelContext and AboveLevelContext arrays in
  // the spec. In the spec, it stores values from 0 through 63 (inclusive). The
  // stored values are used to compute the left and top contexts in
  // GetTransformAllZeroContext. In that function, we only care about the
  // following values: 0, 1, 2, 3 and >= 4. So instead of clamping to 63, we
  // clamp to 4 (i.e.) all the values greater than 4 are stored as 4.
  std::array<Array2D<uint8_t>, 2> coefficient_levels_;
  // This is equivalent to the LeftDcContext and AboveDcContext arrays in the
  // spec. In the spec, it can store 3 possible values: 0, 1 and 2 (where 1
  // means the value is < 0, 2 means the value is > 0 and 0 means the value is
  // equal to 0).
  //
  // The stored values are used in two places:
  //  * GetTransformAllZeroContext: Here, we only care about whether the
  //  value is 0 or not (whether it is 1 or 2 is irrelevant).
  //  * GetDcSignContext: Here, we do the following computation: if the
  //  stored value is 1, we decrement a counter. If the stored value is 2
  //  we increment a counter.
  //
  // Based on this usage, we can simply replace 1 with -1 and 2 with 1 and
  // use that value to compute the counter.
  //
  // The usage on GetTransformAllZeroContext is unaffected since there we
  // only care about whether it is 0 or not.
  std::array<Array2D<int8_t>, 2> dc_categories_;

 private:
  uint32_t columns4x4_;
  uint32_t rows4x4_;
  bool is420_;
};

//------------------------------------------------------------------------------
// From libgav1/src/tile/tile.cc

// ith entry of this array is computed as:
// DivideBy2(TransformSizeToSquareTransformIndex(kTransformSizeSquareMin[i]) +
//           TransformSizeToSquareTransformIndex(kTransformSizeSquareMax[i]) +
//           1)
constexpr uint8_t kTransformSizeContext[kNumTransformSizes] = {
    0, 1, 1, 1, 1, 2, 2, 1, 2, 2, 3, 3, 2, 3, 3, 4, 3, 4, 4};

// Range above kNumQuantizerBaseLevels which the exponential golomb coding
// process is activated.
constexpr int kQuantizerCoefficientBaseRange = 12;
constexpr int kQuantizerCoefficientBaseRangeContextClamp =
    kQuantizerCoefficientBaseRange + kNumQuantizerBaseLevels + 1;
constexpr int kCoeffBaseRangeMaxIterations =
    kQuantizerCoefficientBaseRange / (kCoeffBaseRangeSymbolCount - 1);
constexpr int kEntropyContextLeft = 0;
constexpr int kEntropyContextTop = 1;

constexpr TransformSize kSizeWP2ToAOM[] = {
    kTransformSize4x4,   kTransformSize8x4,   kTransformSize16x4,
    kTransformSize4x8,   kTransformSize8x8,   kTransformSize16x8,
    kTransformSize32x8,  kTransformSize4x16,  kTransformSize8x16,
    kTransformSize16x16, kTransformSize32x16, kTransformSize8x32,
    kTransformSize16x32, kTransformSize32x32};
STATIC_ASSERT_ARRAY_SIZE(kSizeWP2ToAOM, BLK_LAST);
constexpr WP2::BlockSize kSizeAOMToWP2[] = {
    BLK_4x4,   BLK_4x8,  BLK_4x16,  BLK_8x4,   BLK_8x8,  BLK_8x16, BLK_8x32,
    BLK_16x4,  BLK_16x8, BLK_16x16, BLK_16x32, BLK_LAST, BLK_32x8, BLK_32x16,
    BLK_32x32, BLK_LAST, BLK_LAST,  BLK_LAST,  BLK_LAST};
STATIC_ASSERT_ARRAY_SIZE(kSizeAOMToWP2, kNumTransformSizes);

inline TransformSize ConvertSize(WP2::BlockSize size) {
  const TransformSize res = kSizeWP2ToAOM[(uint32_t)size];
  assert(res != kNumTransformSizes && kSizeAOMToWP2[res] == size);
  return res;
}

inline WP2::BlockSize ConvertSize(TransformSize tx_size) {
  const WP2::BlockSize res = kSizeAOMToWP2[(uint32_t)tx_size];
  assert(res != BLK_LAST && kSizeWP2ToAOM[res] == tx_size);
  return res;
}

inline Plane ConvertPlane(Channel channel) {
  switch (channel) {
    case kYChannel:
      return kPlaneY;
    case kUChannel:
      return kPlaneU;
    case kVChannel:
      return kPlaneV;
    case kAChannel:
      return kPlaneA;
  }
  assert(false);
  return kPlaneY;
}

inline Channel ConvertPlane(Plane p) {
  switch (p) {
    case kPlaneY:
      return kYChannel;
    case kPlaneU:
      return kUChannel;
    case kPlaneV:
      return kVChannel;
    case kPlaneA:
      return kAChannel;
  }
  assert(false);
  return kYChannel;
}

inline const uint16_t* UberZigZag(TransformSize tx_size, bool is420) {
  return ResidualIterator::GetZigzag(GetTransform(ConvertSize(tx_size), is420));
}

}  // namespace libgav1

}  // namespace WP2

#endif /* WP2_COMMON_LOSSY_RESIDUALS_AOM_H_ */
