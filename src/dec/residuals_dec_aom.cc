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
//   AOM residual decoding, branched from libgav1, with as little changes as
//   possible.
//
// Author: Vincent Rabaud (vrabaud@google.com)

#include "src/dec/residuals_dec_aom.h"

#include <algorithm>
#include <numeric>

#include "src/utils/utils.h"

namespace WP2 {
namespace libgav1 {

//------------------------------------------------------------------------------
// From libgav1/src/tile/tile.cc

template <bool is_dc_coefficient>
bool AOMResidualReader::ReadSignAndApplyDequantization(
    int adjusted_tx_width_log2, uint16_t dc_sign_cluster,
    int16_t* const coeffs,
    int8_t* const dc_category, int* const coefficient_level) {
  int32_t coeff = *coeffs;
  if (coeff == 0) return true;
  const bool sign = is_dc_coefficient ?
      sr_->Read(dc_sign_cluster, kAOMDCSign, "is_negative_dc") :
      dec_->ReadBool("is_negative");
  if (coeff > kNumQuantizerBaseLevels + kQuantizerCoefficientBaseRange) {
    int length = 0;
    bool golomb_length_bit = false;
    do {
      golomb_length_bit =
          static_cast<bool>(dec_->ReadBool("golomb_length_bit"));
      ++length;
      if (length > 20) {
        // LIBGAV1_DLOG(ERROR, "Invalid golomb_length %d", length);
        return false;
      }
    } while (!golomb_length_bit);
    int x = 1;
    for (int ii = length - 2; ii >= 0; --ii) {
      x = (x << 1) | dec_->ReadBool("more_golomb");
    }
    coeff += x - 1;
  }
  if (is_dc_coefficient && coeff > 0) {
    *dc_category = (sign != 0) ? -1 : 1;
  }
  // Clamp to avoid any integer overflow.
  coeff = std::min(coeff, is_dc_coefficient ? kMaxDcValue : kMaxCoeffValue);
  *coefficient_level += coeff;
  if (sign) coeff = -coeff;
  *coeffs = coeff;
  return true;
}

int AOMResidualReader::ReadCoeffBaseRange(int clamped_tx_size_context,
                                          int cdf_context, int plane_type) {
  int level = 0;
  for (int j = 0; j < kCoeffBaseRangeMaxIterations; ++j) {
    const int coeff_base_range =
        sr_->Read(GetClusterAOM3<kAOMCoeffBaseRange>(clamped_tx_size_context,
                                plane_type, cdf_context),
                  kAOMCoeffBaseRange, "coeff_base_range");
    level += coeff_base_range;
    if (coeff_base_range < (kCoeffBaseRangeSymbolCount - 1)) break;
  }
  return level;
}

WP2Status AOMResidualReader::ReadCoeffs(
    Channel channel, uint32_t x_pix, uint32_t y_pix, WP2::BlockSize dim,
    bool is420, TransformPair tx_type, bool first_is_dc, int max_num_coeffs,
    ANSDec* const dec, SymbolReader* const sr, AOMContext* const aom_context,
    int16_t* const coeffs, uint32_t* const num_coeffs) {
  // Convert WP2 data to libgav1 data.
  const bool half_res = (channel == kUChannel || channel == kVChannel) && is420;
  const TransformSize tx_size = ConvertSize(dim);
  const Plane plane = ConvertPlane(channel);
  dec_ = dec;
  sr_ = sr;

  aom_context->SetIs420(is420);
  const uint32_t subsampling = half_res ? 1u : 0u;

  const int num_coeffs_in_transform = ReadTransformCoefficients(
      plane, x_pix >> subsampling, y_pix >> subsampling, tx_size, tx_type,
      first_is_dc, max_num_coeffs, coeffs, aom_context);
  WP2_CHECK_OK(num_coeffs_in_transform >= 0, WP2_STATUS_BITSTREAM_ERROR);
  *num_coeffs = (uint32_t)num_coeffs_in_transform;
  return WP2_STATUS_OK;
}

int AOMResidualReader::ReadTransformCoefficients(
    Plane plane, int start_x, int start_y, TransformSize tx_size,
    TransformPair tx_type, bool first_is_dc, uint32_t max_num_coeffs,
    int16_t* const coeffs, AOMContext* const aom_context) {
  const int x4 = DivideBy4(start_x);
  const int y4 = DivideBy4(start_y);
  const int w4 = kTransformWidth4x4[tx_size];
  const int h4 = kTransformHeight4x4[tx_size];
  if (max_num_coeffs == 0) {
    SetEntropyContexts(x4, y4, w4, h4, plane, 0, 0, aom_context);
    // This is not used in this case, so it can be set to any value.
    return 0;
  }
  const int tx_size_context = kTransformSizeContext[tx_size];
  const int tx_width = kTransformWidth[tx_size];
  const int tx_height = kTransformHeight[tx_size];
  const int clamped_tx_width = std::min(tx_width, 32);
  const int clamped_tx_height = std::min(tx_height, 32);
  const int padded_tx_width =
      clamped_tx_width + kQuantizedCoefficientBufferPadding;
  const int padded_tx_height =
      clamped_tx_height + kQuantizedCoefficientBufferPadding;
  const int eob_multi_size = kEobMultiSizeLookup[tx_size];
  const PlaneType plane_type = GetPlaneType(plane);
  const uint32_t tx_class = (uint32_t)GetTransformClass(tx_type);
  int context = (tx_class != (uint32_t)TransformClass::kTwoD) ? 1 : 0;
  int16_t eob_pt_tmp;
  switch (eob_multi_size) {
    case 0:
      eob_pt_tmp = sr_->Read(GetClusterAOM2<kAOMEOBPT16>(plane_type, context),
                             kAOMEOBPT16, "eob_pt");
      break;
    case 1:
      eob_pt_tmp = sr_->Read(GetClusterAOM2<kAOMEOBPT32>(plane_type, context),
                             kAOMEOBPT32, "eob_pt");
      break;
    case 2:
      eob_pt_tmp = sr_->Read(GetClusterAOM2<kAOMEOBPT64>(plane_type, context),
                             kAOMEOBPT64, "eob_pt");
      break;
    case 3:
      eob_pt_tmp = sr_->Read(GetClusterAOM2<kAOMEOBPT128>(plane_type, context),
                             kAOMEOBPT128, "eob_pt");
      break;
    case 4:
      eob_pt_tmp = sr_->Read(GetClusterAOM2<kAOMEOBPT256>(plane_type, context),
                             kAOMEOBPT256, "eob_pt");
      break;
    case 5:
      eob_pt_tmp = sr_->Read(GetClusterAOM1<kAOMEOBPT512>(plane_type),
                             kAOMEOBPT512, "eob_pt");
      break;
    case 6:
    default:
      eob_pt_tmp = sr_->Read(GetClusterAOM1<kAOMEOBPT1024>(plane_type),
                             kAOMEOBPT1024, "eob_pt");
      break;
  }
  const int16_t eob_pt = 1 + eob_pt_tmp;
  int16_t eob = (eob_pt < 2) ? eob_pt : ((1 << (eob_pt - 2)) + 1);
  if (eob_pt >= 3) {
    context = eob_pt - 3;
    const bool eob_extra = sr_->Read(
        GetClusterAOM3<kAOMEOBExtra>(tx_size_context, plane_type, context),
        kAOMEOBExtra, "eob_extra");
    if (eob_extra) eob += 1 << (eob_pt - 3);
    for (int i = 1; i < eob_pt - 2; ++i) {
      assert(eob_pt - i >= 3);
      assert(eob_pt <= kEobPt1024SymbolCount);
      if (static_cast<bool>(dec_->ReadBool("eob_more"))) {
        eob += 1 << (eob_pt - i - 3);
      }
    }
  }

  // Only the first |padded_tx_width| * |padded_tx_height| values of |quantized|
  // will be used by this function and the functions to which it is passed into.
  // So we simply need to zero out those values before it is being used. If
  // |eob| == 1, then only the first index will be populated and used. So there
  // is no need to initialize this array in that case.
  uint8_t quantized[2048], quantized2[2048];
  if (eob > 1) {
    std::fill(quantized, quantized + padded_tx_width * padded_tx_height, 0u);
    std::fill(quantized2, quantized2 + padded_tx_width * padded_tx_height, 0u);
  }

  const uint16_t* scan = UberZigZag(tx_size, aom_context->IsReduced(plane));
  const TransformSize adjusted_tx_size = kAdjustedTransformSize[tx_size];
  const int adjusted_tx_width_log2 = kTransformWidthLog2[adjusted_tx_size];
  const int clamped_tx_size_context = std::min(tx_size_context, 3);
  // Read the last coefficient.
  {
    context = GetCoeffBaseContextEob(tx_size, eob - 1);
    if (eob - 1 >= (int)max_num_coeffs) {
      return -1;
    }
    const uint16_t pos = scan[eob - 1];
    const uint32_t cluster =
        GetClusterAOM3<kAOMCoeffBaseEOB>(tx_size_context, plane_type, context);
    int level = 1 + sr_->Read(cluster, kAOMCoeffBaseEOB, "level");
    if (level > kNumQuantizerBaseLevels) {
      level += ReadCoeffBaseRange(
          clamped_tx_size_context,
          GetCoeffBaseRangeContextEob(adjusted_tx_width_log2, pos), plane_type);
    }
    coeffs[pos] = level;
    quantized[PaddedIndex(pos, adjusted_tx_width_log2)] = std::min(level, 3);
    quantized2[PaddedIndex(pos, adjusted_tx_width_log2)] =
        std::min(level, kQuantizerCoefficientBaseRangeContextClamp);
  }

  auto get_coeff_base_context_func =
      AOMResidualIO::GetCoeffsBaseContextFunc(tx_class);
  auto get_coeff_base_range_context_func =
      AOMResidualIO::GetCoeffBaseRangeContextFunc(tx_class);

  // Read all the other coefficients.
  for (int i = eob - 2; i >= 0; --i) {
    const uint16_t pos = scan[i];

    context = get_coeff_base_context_func(quantized, tx_size,
                                          adjusted_tx_width_log2, pos);
    const uint32_t cluster =
        GetClusterAOM3<kAOMCoeffBase>(tx_size_context, plane_type, context);
    int level = sr_->Read(cluster, kAOMCoeffBase, "level");
    if (level > kNumQuantizerBaseLevels) {
      const uint32_t cdf_context = get_coeff_base_range_context_func(
          quantized2, adjusted_tx_width_log2, pos);
      level +=
          ReadCoeffBaseRange(clamped_tx_size_context, cdf_context, plane_type);
    }
    coeffs[pos] = level;
    quantized[PaddedIndex(pos, adjusted_tx_width_log2)] = std::min(level, 3);
    quantized2[PaddedIndex(pos, adjusted_tx_width_log2)] =
        std::min(level, kQuantizerCoefficientBaseRangeContextClamp);
  }
  int coefficient_level = 0;
  int8_t dc_category = 0;
  const uint32_t dc_sign_cluster =
      (quantized[0] != 0)
          ? GetClusterAOM2<kAOMDCSign>(plane_type,
                          GetDcSignContext(*aom_context, x4, y4, w4, h4, plane))
          : 0;
  assert(scan[0] == 0);
  if (first_is_dc &&
      !ReadSignAndApplyDequantization</*is_dc_coefficient=*/true>(
          adjusted_tx_width_log2, dc_sign_cluster, &coeffs[0], &dc_category,
          &coefficient_level)) {
    return -1;
  }
  uint32_t last_coeff = (first_is_dc && coeffs[0]) ? 1 : 0;
  for (int i = (first_is_dc ? 1 : 0); i < eob; ++i) {
    if (coeffs[scan[i]]) {
      if (!ReadSignAndApplyDequantization</*is_dc_coefficient=*/false>(
              adjusted_tx_width_log2, dc_sign_cluster,
              &coeffs[scan[i]], nullptr, &coefficient_level)) {
        return -1;
      }
      last_coeff = std::max(last_coeff, scan[i] + 1u);
    }
  }

  SetEntropyContexts(x4, y4, w4, h4, plane, std::min(4, coefficient_level),
                     dc_category, aom_context);
  return last_coeff;
}

}  // namespace libgav1
}  // namespace WP2
