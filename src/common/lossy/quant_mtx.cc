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
//   Coefficient Quantization
//
// Author: Skal (pascal.massimino@gmail.com)

#include <algorithm>
#include <cstdio>

#include "src/common/lossy/block.h"
#include "src/common/lossy/block_size_io.h"
#include "src/utils/utils.h"

namespace WP2 {

//------------------------------------------------------------------------------

// idx at which we use the next quant parameter (for both x and y directions)
static constexpr uint32_t kZones[kNumQuantZones] = { 1, 4, 16, 32 };

// Linearly interpolate values[] at point 'n' in interval [s, e]
// TODO(skal): use fixed-point precision
static uint32_t Interpolate(const uint16_t values[],
                            uint32_t s, uint32_t e, uint32_t n) {
  const float v = 1.f * (values[0] * (e - n) + values[1] * (n - s)) / (e - s);
  return (uint32_t)std::lrintf(v);
}

WP2Status QuantMtx::Allocate() {
  WP2_CHECK_ALLOC_OK(iquant_data.resize(kDataSize));
  WP2_CHECK_ALLOC_OK(bias_data.resize(kDataSize));
  WP2_CHECK_ALLOC_OK(range_data.resize(kDataSize));
  return WP2_STATUS_OK;
}

void QuantMtx::Init(uint32_t max_residual, float q_scale,
                    PartitionSet partition_set,
                    const uint16_t steps[kNumQuantZones]) {
  for (uint32_t i = 0; i < kNumQuantZones; ++i) {
    assert(steps[i] >= kMinQuantStep && steps[i] <= kMaxQuantStep);
    quant_steps[i] = steps[i];
  }

  // build an interpolated ramp for the quant-steps
  uint32_t qramp[kMaxBlockSizePix];
  {
    uint32_t s = 0, last = 1u;
    for (uint32_t z = 0; z + 1 < kNumQuantZones; ++z) {
      const uint32_t e = kZones[z];
      for (uint16_t n = s; n < e; ++n) {
        last = qramp[n] = Interpolate(steps + z, s, e, n);
      }
      s = e;
    }
    // last segment is flat
    for (; s < kMaxBlockSizePix; ++s) qramp[s] = last;
  }

  // prepare the quant matrices
  const bool for_decoding = iquant_data.empty();
  int16_t* dq_ptr = dequant_data;
  uint32_t* iq_ptr = iquant_data.data();
  uint32_t* b_ptr = bias_data.data();
  uint16_t* range_ptr = range_data.data();
  for (const TrfSize tdim : kAllTrfSizes) {
    // Check if tdim is used in the PartitionSet.
    const bool is_used = TrfIsUsed(partition_set, tdim);
    const uint32_t W = TrfWidth[tdim];
    const uint32_t H = TrfHeight[tdim];
    const uint32_t max_residual_value =
        (uint32_t)(max_residual * std::sqrt(W * H));
    const float norm = 256. / kMaxQuantStep / q_scale;
    for (uint32_t j = is_used ? 0 : H; j < H; ++j) {
      const uint32_t qy = qramp[j * 32 / H];
      for (uint32_t i = 0; i < W; ++i) {
        const bool is_dc = (i == 0 && j == 0);
        const uint32_t qx = qramp[i * 32 / W];
        const uint32_t q = 1u + (uint32_t)((qx + qy) * norm);
        assert(q <= (uint32_t)std::numeric_limits<int16_t>::max());
        const uint32_t qbias = q - 1;
        const uint32_t idx = i + j * W;
        assert(q <= (uint32_t)std::numeric_limits<int16_t>::max());
        const uint32_t IQ = ((1ull << WP2QBits) + qbias) / q;
        const uint32_t B = (1u << WP2QBits) / 3;
        if (!for_decoding) {
          iq_ptr[idx] = IQ;
          b_ptr[idx] = B;
        }
        if (!for_decoding || is_dc) {
          const uint32_t max_range = is_dc ? kMaxDcValue : kMaxCoeffValue;
          // for max_range calculation, one must use the actual quantization
          // formula, instead of the exact divide by q.
          // This quantization formulae is the same as the one used in dsp/,
          // only rounded up for safety.
          const uint16_t range = std::min(max_range,
              (max_residual_value * IQ + (1 << WP2QBits) - 1) >> WP2QBits);
          if (!for_decoding) range_ptr[idx] = range;
          if (is_dc) max_abs_res_dc[tdim] = range;
        }
        dq_ptr[idx] = (int16_t)q;
      }
    }
    // Because of the assumption of WP2Quantize(), we need to round up the size
    // to next multiple of 8.
    const uint32_t size8 = (W * H + 7) & ~7;
    dequants[tdim] = dq_ptr;
    iquants[tdim] = iq_ptr;
    bias[tdim] = b_ptr;
    max_abs_res[tdim] = range_ptr;
    dq_ptr += size8;
    if (!for_decoding) {
      iq_ptr += size8;
      b_ptr += size8;
      range_ptr += size8;
    }
  }
  assert(dq_ptr == &dequant_data[kDataSize]);
}

void QuantMtx::SetLambda(float lambda) {
  // this expression was fit over the sampling of several images encoded at
  // fixed quality (=quant_steps[1]) with varying lambda.
  const double lambda_f = 450. * pow(quant_steps[1] / 1019., 1.8);
  assert(lambda_f > 0.);
  lambda1 = lambda_f;
  lambda2 = 0.84 * lambda_f;
}

void QuantMtx::Quantize(const int32_t residuals[], TrfSize dim,
                        bool first_is_dc, int16_t coeffs[],
                        uint32_t* const num_coeffs) const {
  const uint32_t num_coeffs_full = kNumCoeffs[dim];
  WP2Quantize(iquants[dim], bias[dim], residuals, coeffs, num_coeffs_full);
  uint32_t last_nz = num_coeffs_full;
  while (last_nz > 0 && coeffs[last_nz - 1] == 0) --last_nz;
  *num_coeffs = last_nz;

#if defined(WP2_ERROR_TRACE) || defined(WP2_TRACE)
  for (uint32_t i = 0; i < num_coeffs_full; ++i) {
    const uint32_t v = std::abs(coeffs[i]);
    assert(v <= max_abs_res[dim][i]);
    if (i == 0 && v > kMaxDcValue) {
      // TODO(maryla): make sure clipping never occurs? Clipping is bad, it's
      // a form of unexpected quantization that happens at high quality (but
      // skal disagrees).
      fprintf(stderr, "WARNING!! Clipping occured. Clipped DC %d to %d\n", v,
              kMaxDcValue);
    } else if (i > 0 && v > kMaxCoeffValue) {
      fprintf(stderr, "WARNING!! Clipping occured. Clipped value %d to %d\n", v,
              kMaxCoeffValue);
    }
  }
#endif
  if (first_is_dc) {
    coeffs[0] = Clamp<int16_t>(coeffs[0], -kMaxDcValue, kMaxDcValue);
  }
  for (uint32_t i = (first_is_dc ? 1 : 0); i < last_nz; ++i) {
    coeffs[i] = Clamp<int16_t>(coeffs[i], -kMaxCoeffValue, kMaxCoeffValue);
  }
}

int16_t QuantMtx::DCError(int32_t DC, TrfSize dim) const {
  const uint32_t iq = iquants[dim][0];
  const uint32_t b = bias[dim][0];
  int32_t v = ((uint64_t)std::abs(DC) * iq + b) >> WP2QBits;
  v = std::min(v, kMaxDcValue) * dequants[dim][0] * (DC < 0 ? -1 : 1);
  return DC - v;
}

void QuantMtx::Dequantize(const int16_t coeffs[], uint32_t num_coeffs,
                          TrfSize dim, int32_t residuals[]) const {
  WP2Dequantize(coeffs, dequants[dim], residuals, num_coeffs, kNumCoeffs[dim]);
}

uint32_t QuantMtx::RoundTrip(const int32_t residuals[], TrfSize dim,
                             bool first_is_dc, int16_t coeffs[],
                             uint32_t* const num_coeffs) const {
  int32_t tmp[kMaxBlockSizePix2];
  Quantize(residuals, dim, first_is_dc, coeffs, num_coeffs);
  // TODO(skal): rnd_mtx ?
  WP2Dequantize(coeffs, dequants[dim], tmp, *num_coeffs, kNumCoeffs[dim]);

  uint32_t error = 0;
  for (uint32_t i = 0; i < kNumCoeffs[dim]; ++i) {
    const int32_t diff = (int32_t)residuals[i] - tmp[i];
    error += diff * diff;
  }
  return error / kNumCoeffs[dim];
}

void QuantMtx::Print() const {
  for (const TrfSize tdim : kAllTrfSizes) {
    printf(" --- Q-%s ---\n", kTDimNames[tdim]);
    const uint32_t W = TrfWidth[tdim];
    const uint32_t H = TrfHeight[tdim];
    for (uint32_t j = 0; j < H; ++j) {
      for (uint32_t i = 0; i < W; ++i) {
        const uint32_t idx = i + W * j;
        printf("%4u ", (uint32_t)((1ull << WP2QBits) / iquants[tdim][idx]));
      }
      printf("\n");
    }
  }
  printf("quant_steps[]: ");
  for (const auto& q : quant_steps) printf("[%d]", q);
  printf("\n");
}

//------------------------------------------------------------------------------

void QuantStat::Reset() {
  for (const TrfSize tdim : kAllTrfSizes) {
    const uint32_t num_coeffs = kNumCoeffs[tdim];
    last_idx_[tdim] = 0;
    for (uint32_t i = 0; i < num_coeffs; ++i) {
      is_large_[tdim][i] = false;
    }
  }
  for (uint32_t i = 0; i < kNumQuantZones; ++i) {
    count_[i] = 0;
    value_[i] = 0;
  }
}

static uint16_t IdxToZone(TrfSize tdim, uint32_t idx) {
  const uint32_t W = TrfWidth[tdim];
  const uint32_t x = idx % W;
  const uint32_t y = idx / W;
  uint32_t zi = 0, zj = 0;
  while (x > kZones[zi]) ++zi;
  while (y > kZones[zj]) ++zj;
  return std::min(zi, zj);
}

// These global constants are inverse of kDynNorm[] * q_scale.
// TODO(maryla): this makes no sense since q_scale is not a constant.
static constexpr float kZoneNorm[TRF_LAST] = {
  0,      0,      0,
  0, 2.1609, 3.1458, 4.5276,
  0, 3.1458, 4.5796,  6.5912, 9.4160,
     4.5276, 6.5912, 9.4864, 13.5520,
            9.4160, 13.5520, 19.36000
};
STATIC_ASSERT_ARRAY_SIZE(kZoneNorm, TRF_LAST);

void QuantStat::Record(TrfSize tdim,
                       uint32_t idx, int32_t qvalue, int32_t value) {
  assert(idx < kNumCoeffs[tdim]);
  qvalue = std::abs(qvalue);
  const bool large = (qvalue >= 2);
  is_large_[tdim][idx] |= large;
  if (qvalue > 0) {
    last_idx_[tdim] = std::max(last_idx_[tdim], idx + 1);
  }
  if (qvalue == 1) {
    const uint32_t z = IdxToZone(tdim, idx);
    ++count_[z];
    value_[z] += std::abs(value) / kZoneNorm[tdim] / 2.;
  } else {
    // TODO(skal): value /= qvalue ?
  }
}

void QuantStat::Record(TrfSize tdim,
                       int32_t* const coeffs, int16_t* const qcoeffs) {
  const uint32_t num_coeffs = kNumCoeffs[tdim];
  for (uint32_t i = 0; i < num_coeffs; ++i) {
    if (qcoeffs[i]) Record(tdim, i, qcoeffs[i], coeffs[i]);
  }
}

void QuantStat::Print() const {
  for (uint32_t z = 0; z < kNumQuantZones; ++z) {
    if (count_[z] > 0) {
      const uint32_t q = (uint32_t)(value_[z] / count_[z] + .5);
      printf("%3d ", q);
    } else {
      printf("  .  ");
    }
  }
  printf("\n");
}

bool QuantStat::AdjustQuantSteps(uint16_t quant_steps[kNumQuantZones]) const {
  bool changed = false;
  for (uint32_t z = 0; z < kNumQuantZones; ++z) {
    if (count_[z] > 0) {
      const uint16_t q = (uint32_t)(value_[z] / count_[z] + .5);
      changed |= (q != quant_steps[z]);
      quant_steps[z] = q;
    }
  }
  return changed;
}

//------------------------------------------------------------------------------

}    // namespace WP2
