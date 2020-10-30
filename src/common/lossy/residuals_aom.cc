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

#include "src/common/lossy/residuals_aom.h"

#include <numeric>

#include "src/common/lossy/aom/cdfs.inc"
#include "src/utils/utils.h"

namespace WP2 {
namespace libgav1 {

//------------------------------------------------------------------------------
// From libgav1/src/symbol_decoder_context.cc

#define CDF_COPY(source, destination) \
  memcpy(destination, source, sizeof(source))

void SymbolDecoderContext::Initialize(int quantizer_context) {
  CDF_COPY(kDefaultAllZeroCdf[quantizer_context], all_zero_cdf);
  CDF_COPY(kDefaultEobPt16Cdf[quantizer_context], eob_pt_16_cdf);
  CDF_COPY(kDefaultEobPt32Cdf[quantizer_context], eob_pt_32_cdf);
  CDF_COPY(kDefaultEobPt64Cdf[quantizer_context], eob_pt_64_cdf);
  CDF_COPY(kDefaultEobPt128Cdf[quantizer_context], eob_pt_128_cdf);
  CDF_COPY(kDefaultEobPt256Cdf[quantizer_context], eob_pt_256_cdf);
  CDF_COPY(kDefaultEobPt512Cdf[quantizer_context], eob_pt_512_cdf);
  CDF_COPY(kDefaultEobPt1024Cdf[quantizer_context], eob_pt_1024_cdf);
  CDF_COPY(kDefaultEobExtraCdf[quantizer_context], eob_extra_cdf);
  CDF_COPY(kDefaultCoeffBaseEobCdf[quantizer_context], coeff_base_eob_cdf);
  CDF_COPY(kDefaultCoeffBaseCdf[quantizer_context], coeff_base_cdf);
  CDF_COPY(kDefaultCoeffBaseRangeCdf[quantizer_context], coeff_base_range_cdf);
  CDF_COPY(kDefaultDcSignCdf[quantizer_context], dc_sign_cdf);
}

#undef CDF_COPY

//------------------------------------------------------------------------------
// From libgav1/src/utils/constants.cc

const uint8_t kBlockWidthPixels[kMaxBlockSizes] = {
    4,  4,  4,  8,  8,  8,  8,  16, 16, 16,  16,
    16, 32, 32, 32, 32, 64, 64, 64, 64, 128, 128};

const uint8_t kBlockHeightPixels[kMaxBlockSizes] = {
    4,  8, 16, 4,  8,  16, 32, 4,  8,   16, 32,
    64, 8, 16, 32, 64, 16, 32, 64, 128, 64, 128};

const uint8_t kTransformWidth[kNumTransformSizes] = {
    4, 4, 4, 8, 8, 8, 8, 16, 16, 16, 16, 16, 32, 32, 32, 32, 64, 64, 64};

const uint8_t kTransformHeight[kNumTransformSizes] = {
    4, 8, 16, 4, 8, 16, 32, 4, 8, 16, 32, 64, 8, 16, 32, 64, 16, 32, 64};

const uint8_t kTransformWidth4x4[kNumTransformSizes] = {
    1, 1, 1, 2, 2, 2, 2, 4, 4, 4, 4, 4, 8, 8, 8, 8, 16, 16, 16};

const uint8_t kTransformHeight4x4[kNumTransformSizes] = {
    1, 2, 4, 1, 2, 4, 8, 1, 2, 4, 8, 16, 2, 4, 8, 16, 4, 8, 16};

const uint8_t kTransformWidthLog2[kNumTransformSizes] = {
    2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4, 4, 5, 5, 5, 5, 6, 6, 6};

//------------------------------------------------------------------------------

namespace {

/* clang-format off */
constexpr uint8_t kCoeffBaseContextOffset[kNumTransformSizes][5][5] = {
    {{0,  1,  6,  6, 0},
     {1,  6,  6, 21, 0},
     {6,  6, 21, 21, 0},
     {6, 21, 21, 21, 0},
     {0,  0,  0,  0, 0}},    // 4x4

    {{ 0, 11, 11, 11, 0},
     {11, 11, 11, 11, 0},
     { 6,  6, 21, 21, 0},
     { 6, 21, 21, 21, 0},
     {21, 21, 21, 21, 0}},   // 4x8

    {{ 0, 11, 11, 11, 0},
     {11, 11, 11, 11, 0},
     { 6,  6, 21, 21, 0},
     { 6, 21, 21, 21, 0},
     {21, 21, 21, 21, 0}},   // 4x16

    {{ 0, 16,  6,  6, 21},
     {16, 16,  6, 21, 21},
     {16, 16, 21, 21, 21},
     {16, 16, 21, 21, 21},
     { 0,  0,  0,  0,  0}},  // 8x4

    {{ 0,  1,  6,  6, 21},
     { 1,  6,  6, 21, 21},
     { 6,  6, 21, 21, 21},
     { 6, 21, 21, 21, 21},
     {21, 21, 21, 21, 21}},  // 8x8

    {{ 0, 11, 11, 11, 11},
     {11, 11, 11, 11, 11},
     { 6,  6, 21, 21, 21},
     { 6, 21, 21, 21, 21},
     {21, 21, 21, 21, 21}},  // 8x16

    {{ 0, 11, 11, 11, 11},
     {11, 11, 11, 11, 11},
     { 6,  6, 21, 21, 21},
     { 6, 21, 21, 21, 21},
     {21, 21, 21, 21, 21}},  // 8x32

    {{ 0, 16,  6,  6, 21},
     {16, 16,  6, 21, 21},
     {16, 16, 21, 21, 21},
     {16, 16, 21, 21, 21},
     { 0,  0,  0,  0,  0}},  // 16x4

    {{ 0, 16,  6,  6, 21},
     {16, 16,  6, 21, 21},
     {16, 16, 21, 21, 21},
     {16, 16, 21, 21, 21},
     {16, 16, 21, 21, 21}},  // 16x8

    {{ 0,  1,  6,  6, 21},
     { 1,  6,  6, 21, 21},
     { 6,  6, 21, 21, 21},
     { 6, 21, 21, 21, 21},
     {21, 21, 21, 21, 21}},  // 16x16

    {{ 0, 11, 11, 11, 11},
     {11, 11, 11, 11, 11},
     { 6,  6, 21, 21, 21},
     { 6, 21, 21, 21, 21},
     {21, 21, 21, 21, 21}},  // 16x32

    {{ 0, 11, 11, 11, 11},
     {11, 11, 11, 11, 11},
     { 6,  6, 21, 21, 21},
     { 6, 21, 21, 21, 21},
     {21, 21, 21, 21, 21}},  // 16x64

    {{ 0, 16,  6,  6, 21},
     {16, 16,  6, 21, 21},
     {16, 16, 21, 21, 21},
     {16, 16, 21, 21, 21},
     {16, 16, 21, 21, 21}},  // 32x8

    {{ 0, 16,  6,  6, 21},
     {16, 16,  6, 21, 21},
     {16, 16, 21, 21, 21},
     {16, 16, 21, 21, 21},
     {16, 16, 21, 21, 21}},  // 32x16

    {{ 0,  1,  6,  6, 21},
     { 1,  6,  6, 21, 21},
     { 6,  6, 21, 21, 21},
     { 6, 21, 21, 21, 21},
     {21, 21, 21, 21, 21}},  // 32x32

    {{ 0, 11, 11, 11, 11},
     {11, 11, 11, 11, 11},
     { 6,  6, 21, 21, 21},
     { 6, 21, 21, 21, 21},
     {21, 21, 21, 21, 21}},  // 32x64

    {{ 0, 16,  6,  6, 21},
     {16, 16,  6, 21, 21},
     {16, 16, 21, 21, 21},
     {16, 16, 21, 21, 21},
     {16, 16, 21, 21, 21}},  // 64x16

    {{ 0, 16,  6,  6, 21},
     {16, 16,  6, 21, 21},
     {16, 16, 21, 21, 21},
     {16, 16, 21, 21, 21},
     {16, 16, 21, 21, 21}},  // 64x32

    {{ 0,  1,  6,  6, 21},
     { 1,  6,  6, 21, 21},
     { 6,  6, 21, 21, 21},
     { 6, 21, 21, 21, 21},
     {21, 21, 21, 21, 21}}   // 64x64
};
/* clang-format on */

constexpr uint8_t kCoeffBasePositionContextOffset[3] = {26, 31, 36};

// Returns the minimum of |length| or |max|-|start|. This is used to clamp array
// indices when accessing arrays whose bound is equal to |max|.
int GetNumElements(int length, int start, int max) {
  return std::min(length, max - start);
}

}  // namespace

static int PlaneCount() {
  return kMaxPlanes;
}

WP2Status AOMContext::Init(uint32_t width, uint32_t height) {
  is420_ = false;
  columns4x4_ = (width + 3) >> 2;
  rows4x4_ = (height + 3) >> 2;
  for (size_t i = 0; i < coefficient_levels_.size(); ++i) {
    const int contexts_per_plane =
        (i == kEntropyContextLeft) ? rows4x4_ : columns4x4_;
    WP2_CHECK_ALLOC_OK(
        coefficient_levels_[i].Reset(PlaneCount(), contexts_per_plane));
    WP2_CHECK_ALLOC_OK(
        dc_categories_[i].Reset(PlaneCount(), contexts_per_plane));
  }
  return WP2_STATUS_OK;
}

WP2Status AOMContext::CopyFrom(const AOMContext& other) {
  for (uint32_t i = 0; i < coefficient_levels_.size(); ++i) {
    WP2_CHECK_ALLOC_OK(
        coefficient_levels_[i].CopyFrom(other.coefficient_levels_[i]));
  }
  for (uint32_t i = 0; i < dc_categories_.size(); ++i) {
    WP2_CHECK_ALLOC_OK(dc_categories_[i].CopyFrom(other.dc_categories_[i]));
  }
  columns4x4_ = other.columns4x4_;
  rows4x4_ = other.rows4x4_;
  is420_ = other.is420_;
  return WP2_STATUS_OK;
}

void AOMContext::Reset() {
  for (size_t i = 0; i < coefficient_levels_.size(); ++i) {
    const int contexts_per_plane =
        (i == kEntropyContextLeft) ? rows4x4_ : columns4x4_;
    if (!coefficient_levels_[i].Reset(PlaneCount(), contexts_per_plane)) {
      assert(false);
    }
    if (!dc_categories_[i].Reset(PlaneCount(), contexts_per_plane)) {
      assert(false);
    }
  }
}

// Contrary to the AV1 spec, the transform size and the residual size have the
// same size.
int AOMResidualIO::GetTransformAllZeroContext(const AOMContext& aom_context,
                                              Plane plane,
                                              TransformSize tx_size, int x4,
                                              int y4, int w4, int h4) {
  const int max_x4x4 = aom_context.Columns4x4(plane);
  const int max_y4x4 = aom_context.Rows4x4(plane);

  int top = 0;
  int left = 0;
  const int num_top_elements = GetNumElements(w4, x4, max_x4x4);
  const int num_left_elements = GetNumElements(h4, y4, max_y4x4);
  if (plane == kPlaneY || plane == kPlaneA) return 0;
  const uint8_t* coefficient_levels =
      &aom_context.coefficient_levels_[kEntropyContextTop][plane][x4];
  const int8_t* dc_categories =
      &aom_context.dc_categories_[kEntropyContextTop][plane][x4];
  for (int i = 0; i < num_top_elements; ++i) {
    top |= coefficient_levels[i];
    top |= dc_categories[i];
  }
  coefficient_levels =
      &aom_context.coefficient_levels_[kEntropyContextLeft][plane][y4];
  dc_categories = &aom_context.dc_categories_[kEntropyContextLeft][plane][y4];
  for (int i = 0; i < num_left_elements; ++i) {
    left |= coefficient_levels[i];
    left |= dc_categories[i];
  }
  return static_cast<int>(top != 0) + static_cast<int>(left != 0) + 7;
}

// Section 8.3.2 in the spec, under coeff_base_eob.
int AOMResidualIO::GetCoeffBaseContextEob(TransformSize tx_size, int index) {
  if (index == 0) return 0;
  const TransformSize adjusted_tx_size = kAdjustedTransformSize[tx_size];
  const int tx_width_log2 = kTransformWidthLog2[adjusted_tx_size];
  const int tx_height = kTransformHeight[adjusted_tx_size];
  if (index <= DivideBy8(tx_height << tx_width_log2)) return 1;
  if (index <= DivideBy4(tx_height << tx_width_log2)) return 2;
  return 3;
}

// Section 8.3.2 in the spec, under coeff_base.
static int GetCoeffBaseContext2D(const uint8_t qlevels[],
                                 TransformSize tx_size,
                                 int adjusted_tx_width_log2,
                                 uint16_t pos) {
  if (pos == 0) return 0;
  const int tx_width = 1 << adjusted_tx_width_log2;
  const int padded_tx_width = tx_width + kQuantizedCoefficientBufferPadding;
  const uint8_t* const levels =
      &qlevels[PaddedIndex(pos, adjusted_tx_width_log2)];
  const int context = std::min(
      4, DivideBy2(1 + levels[1] + levels[2] +
                       levels[padded_tx_width] +
                       levels[padded_tx_width + 1] +
                       levels[2 * padded_tx_width]));
  const int row = pos >> adjusted_tx_width_log2;
  const int column = pos & (tx_width - 1);
  return context + kCoeffBaseContextOffset[tx_size][std::min(row, 4)]
                                          [std::min(column, 4)];
}

// Section 8.3.2 in the spec, under coeff_base.
static int GetCoeffBaseContextHorizontal(
    const uint8_t qlevels[], TransformSize /*tx_size*/,
    int adjusted_tx_width_log2, uint16_t pos) {
  const int tx_width = 1 << adjusted_tx_width_log2;
  const int padded_tx_width = tx_width + kQuantizedCoefficientBufferPadding;
  const uint8_t* const levels =
      &qlevels[PaddedIndex(pos, adjusted_tx_width_log2)];
  const int context = std::min(
      4, DivideBy2(1 + levels[1] + levels[2] + levels[3] + levels[4] +
                       levels[padded_tx_width]));
  const int index = pos & (tx_width - 1);
  return context + kCoeffBasePositionContextOffset[std::min(index, 2)];
}

// Section 8.3.2 in the spec, under coeff_base.
static int GetCoeffBaseContextVertical(
    const uint8_t qlevels[], TransformSize /*tx_size*/,
    int adjusted_tx_width_log2, uint16_t pos) {
  const int tx_width = 1 << adjusted_tx_width_log2;
  const int padded_tx_width = tx_width + kQuantizedCoefficientBufferPadding;
  const uint8_t* const levels =
      &qlevels[PaddedIndex(pos, adjusted_tx_width_log2)];
  const int context = std::min(
      4, DivideBy2(1 + levels[1] +
                       levels[1 * padded_tx_width] +
                       levels[2 * padded_tx_width] +
                       levels[3 * padded_tx_width] +
                       levels[4 * padded_tx_width]));
  const int index = pos >> adjusted_tx_width_log2;
  return context + kCoeffBasePositionContextOffset[std::min(index, 2)];
}

static constexpr CoeffBaseContextFunc kCoeffBaseContextFuncs[] = {
  GetCoeffBaseContext2D,
  GetCoeffBaseContextHorizontal,
  GetCoeffBaseContextVertical
};

CoeffBaseContextFunc
  AOMResidualIO::GetCoeffsBaseContextFunc(uint32_t tx_class) {
  return kCoeffBaseContextFuncs[tx_class];
}

// Section 8.3.2 in the spec, under coeff_br.
static int GetCoeffBaseRangeContext2D(
    const uint8_t qlevels[], int adjusted_tx_width_log2,
    int pos) {
  const uint8_t tx_width = 1 << adjusted_tx_width_log2;
  const int padded_tx_width = tx_width + kQuantizedCoefficientBufferPadding;
  const uint8_t* const levels =
      &qlevels[PaddedIndex(pos, adjusted_tx_width_log2)];
  const int context = std::min(
      6, DivideBy2(1 + levels[1] +
                       levels[padded_tx_width] +
                       levels[padded_tx_width + 1]));
  if (pos == 0) return context;
  const int row = pos >> adjusted_tx_width_log2;
  const int column = pos & (tx_width - 1);
  return context + (((row | column) < 2) ? 7 : 14);
}

// Section 8.3.2 in the spec, under coeff_br.
static int GetCoeffBaseRangeContextHorizontal(
    const uint8_t qlevels[], int adjusted_tx_width_log2, int pos) {
  const uint8_t tx_width = 1 << adjusted_tx_width_log2;
  const int padded_tx_width = tx_width + kQuantizedCoefficientBufferPadding;
  const uint8_t* const levels =
      &qlevels[PaddedIndex(pos, adjusted_tx_width_log2)];
  const int context = std::min(
      6, DivideBy2(1 + levels[1] + levels[2] + levels[padded_tx_width]));
  if (pos == 0) return context;
  const int column = pos & (tx_width - 1);
  return context + ((column == 0) ? 7 : 14);
}

// Section 8.3.2 in the spec, under coeff_br.
static int GetCoeffBaseRangeContextVertical(
    const uint8_t qlevels[], int adjusted_tx_width_log2, int pos) {
  const uint8_t tx_width = 1 << adjusted_tx_width_log2;
  const int padded_tx_width = tx_width + kQuantizedCoefficientBufferPadding;
  const uint8_t* const levels =
      &qlevels[PaddedIndex(pos, adjusted_tx_width_log2)];
  const int context = std::min(
      6, DivideBy2(1 + levels[1] +
                       levels[padded_tx_width] + levels[2 * padded_tx_width]));
  if (pos == 0) return context;
  const int row = pos >> adjusted_tx_width_log2;
  return context + ((row == 0) ? 7 : 14);
}

static constexpr CoeffBaseRangeContextFunc kCoeffBaseRangeContextFuncs[] = {
  GetCoeffBaseRangeContext2D,
  GetCoeffBaseRangeContextHorizontal,
  GetCoeffBaseRangeContextVertical
};

CoeffBaseRangeContextFunc
  AOMResidualIO::GetCoeffBaseRangeContextFunc(uint32_t tx_class) {
  return kCoeffBaseRangeContextFuncs[tx_class];
}

// Section 8.3.2 in the spec, under coeff_br. Optimized for end of block based
// on the fact that {0, 1}, {1, 0}, {1, 1}, {0, 2} and {2, 0} will all be 0 in
// the end of block case.
int AOMResidualIO::GetCoeffBaseRangeContextEob(int adjusted_tx_width_log2,
                                               int pos) {
  if (pos == 0) return 0;
  const uint8_t tx_width = 1 << adjusted_tx_width_log2;
  const int row = pos >> adjusted_tx_width_log2;
  const int column = pos & (tx_width - 1);
  return ((row | column) < 2) ? 7 : 14;
}

int AOMResidualIO::GetDcSignContext(const AOMContext& aom_context, int x4,
                                    int y4, int w4, int h4, Plane plane) {
  const int max_x4x4 = aom_context.Columns4x4(plane);
  const int8_t* dc_categories =
      &aom_context.dc_categories_[kEntropyContextTop][plane][x4];
  int dc_sign = std::accumulate(
      dc_categories, dc_categories + GetNumElements(w4, x4, max_x4x4), 0);
  const int max_y4x4 = aom_context.Rows4x4(plane);
  dc_categories = &aom_context.dc_categories_[kEntropyContextLeft][plane][y4];
  dc_sign = std::accumulate(
      dc_categories, dc_categories + GetNumElements(h4, y4, max_y4x4), dc_sign);
  // This return statement is equivalent to:
  //   if (dc_sign < 0) return 1;
  //   if (dc_sign > 0) return 2;
  //   return 0;
  return static_cast<int>(dc_sign < 0) +
         MultiplyBy2(static_cast<int>(dc_sign > 0));
}

void AOMResidualIO::SetEntropyContexts(int x4, int y4, int w4, int h4,
                                       Plane plane, uint8_t coefficient_level,
                                       int8_t dc_category,
                                       AOMContext* const aom_context) {
  const int max_x4x4 = aom_context->Columns4x4(plane);
  const int num_top_elements = GetNumElements(w4, x4, max_x4x4);
  assert(sizeof(aom_context->coefficient_levels_[0][0][0]) == 1);
  memset(&aom_context->coefficient_levels_[kEntropyContextTop][plane][x4],
         coefficient_level, num_top_elements);
  memset(&aom_context->dc_categories_[kEntropyContextTop][plane][x4],
         dc_category, num_top_elements);
  const int max_y4x4 = aom_context->Rows4x4(plane);
  const int num_left_elements = GetNumElements(h4, y4, max_y4x4);
  memset(&aom_context->coefficient_levels_[kEntropyContextLeft][plane][y4],
         coefficient_level, num_left_elements);
  memset(&aom_context->dc_categories_[kEntropyContextLeft][plane][y4],
         dc_category, num_left_elements);
}

}  // namespace libgav1
}  // namespace WP2
