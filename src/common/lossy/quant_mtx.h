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
//   Quantization matrices
//
// Author: Skal (pascal.massimino@gmail.com)

#ifndef WP2_COMMON_LOSSY_QUANT_MTX_H_
#define WP2_COMMON_LOSSY_QUANT_MTX_H_

#include "src/common/constants.h"
#include "src/common/lossy/block_size.h"
#include "src/wp2/base.h"
#include "src/utils/vector.h"

namespace WP2 {

class RndMtxSet;

//------------------------------------------------------------------------------
// Quantization

static constexpr uint32_t kNumQuantZones = 4;
static constexpr uint32_t kMinQuantStep = 128u;
static constexpr uint32_t kMaxQuantStep = 5800u;

class QuantStat {
 public:
  void Record(TrfSize tdim,
              int32_t* const coeffs,  /* unquantized coeffs */
              int16_t* const qcoeffs  /* quantized coeffs */);
  void Reset();
  // Returns true if one quant_steps[] value was changed at least.
  bool AdjustQuantSteps(uint16_t quant_steps[kNumQuantZones]) const;

  void Print() const;   // for debug

 private:
  void Record(TrfSize tdim, uint32_t idx, int32_t qvalue, int32_t value);
  // Whether any values are quantized to > 1
  bool is_large_[TRF_LAST][kMaxBlockSizePix2];
  // Count of values that are quantized to 1 (per quant zone).
  uint32_t count_[kNumQuantZones];
  // Sum of values that are quantized to 1 (per quant zone).
  float value_[kNumQuantZones];
  uint32_t last_idx_[TRF_LAST];
};

struct QuantMtx {
 public:
  // allocate the *encoding* data space
  WP2Status Allocate();

  // Initialize quant_steps[] for steps[].
  // max_residual is the maximum absolute input values to be quantized.
  // q_scale is the expected precision scale of the input Y/U/V coeffs
  // (typically 1<<yuv_bits with yuv_bits ~= 8-10bits, depending on the CSP).
  // Does not initialize for TrfSizes that are not part of the 'partition_set'.
  void Init(uint32_t max_residual, float q_scale, PartitionSet partition_set,
            const uint16_t steps[kNumQuantZones]);
  // Set the lambda magic constants (encoder), using quant_steps[].
  // Must be called after Init().
  void SetLambda(float lambda);

  // Quantizes or dequantizes.
  // If 'first_is_dc', 'coeffs[0]' will be treated differently.
  void Quantize(const int32_t residuals[], TrfSize dim, bool first_is_dc,
                int16_t coeffs[], uint32_t* const num_coeffs) const;
  void Dequantize(const int16_t coeffs[], uint32_t num_coeffs, TrfSize dim,
                  int32_t residuals[]) const;
  // Performs quantization and dequantization and returns sum of squared errors.
  // 'coeffs' is set to the quantized values.
  uint32_t RoundTrip(const int32_t residuals[], TrfSize dim, bool first_is_dc,
                     int16_t coeffs[], uint32_t* const num_coeffs) const;

  // Computes the quantization error on DC.
  // Note that DC, and hence the DC error, is scaled by sqrt(dim/2) by
  // WP2Transform2D. To compare DC error values for blocks of different sizes,
  // the value should be normalized.
  int16_t DCError(int32_t DC, TrfSize dim) const;

  uint16_t GetMaxAbsResDC(TrfSize tdim) const { return max_abs_res_dc[tdim]; }

  QuantStat* qstats = nullptr;
  void SetQuantStats(QuantStat* const stats) { qstats = stats; }

  void Print() const;  // debug

 public:
  // lambda multipliers
  float lambda1;   // used for Y mode-decision
  float lambda2;   // used for UV mode-decision

 protected:
  uint16_t quant_steps[kNumQuantZones] = {0};   // quant levels for each zones

  // quantization and dequantization factors for converting between transform
  // coefficients and coded levels:
  //   level[] = (coeff[] * iquant[] + bias[]) >> WP2QBits
  //   coeffs[] = level[] * dequant[]
  // quant[] incorporates the transform scaling constant sqrt(W * H / 4).
  const int16_t* dequants[TRF_LAST];   // 15bits (always > 0)
  const uint32_t* iquants[TRF_LAST];
  const uint32_t* bias[TRF_LAST];
  // Ranges of quantized residuals - 1, i.e. their max absolute values.
  const uint16_t* max_abs_res[TRF_LAST];
  uint16_t max_abs_res_dc[TRF_LAST];   // permanent storage for DC

  // data for all quant matrices
  static constexpr size_t kDataSize =
      8 +  // < should be 2 * 2, but we round up to 8 bytes for WP2Quantize()
      4 * 4 + 8 * 8 + 16 * 16 +  32 * 32 +
      2 * (2 * 4 + 2 * 8 + 4 * 8 + 4 * 16 + 8 * 16 + 8 * 32 + 16 * 32);
  int16_t dequant_data[kDataSize];
  Vector_u32 iquant_data;
  Vector_u32 bias_data;
  Vector_u16 range_data;
};

//------------------------------------------------------------------------------

}    // namespace WP2

#endif  /* WP2_COMMON_LOSSY_QUANT_MTX_H_ */
