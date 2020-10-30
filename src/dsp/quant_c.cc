// Copyright 2020 Google LLC
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
// Quantization
//
// Author: Skal (pascal.massimino@gmail.com)

#include <algorithm>
#include <cassert>

#include "src/dsp/dsp.h"

//------------------------------------------------------------------------------

namespace {

static inline int16_t Quantize(uint32_t v, uint32_t iq, uint32_t bias) {
  return (int16_t)((v * iq + bias) >> WP2QBits);
}

void Quantize_C(const uint32_t iq[], const uint32_t bias[],
                const int32_t res[], int16_t coeffs[], uint32_t len) {
  for (uint32_t i = 0; i < len; ++i) {
    coeffs[i] = (res[i] < 0) ? -Quantize(-res[i], iq[i], bias[i])
                             :  Quantize( res[i], iq[i], bias[i]);
  }
}

void Dequantize_C(const int16_t in[], const int16_t dequants[],
                  int32_t out[], uint32_t len , uint32_t max_len) {
  for (uint32_t i = 0; i < len; ++i) out[i] = in[i] * dequants[i];
  std::fill(out + len, out + max_len, 0);
}

}  // namespace

//------------------------------------------------------------------------------

void (*WP2Quantize)(const uint32_t iq[], const uint32_t bias[],
                    const int32_t res[], int16_t coeffs[],
                    uint32_t len) = nullptr;
void (*WP2Dequantize)(const int16_t in[], const int16_t dequants[],
                      int32_t out[], uint32_t len , uint32_t max_len) = nullptr;

static volatile WP2CPUInfo quantize_last_cpuinfo_used =
    (WP2CPUInfo)&quantize_last_cpuinfo_used;

extern void WP2QuantizeInitSSE();

WP2_TSAN_IGNORE_FUNCTION void WP2QuantizeInit() {
  if (quantize_last_cpuinfo_used == WP2GetCPUInfo) return;

  WP2Quantize = Quantize_C;
  WP2Dequantize = Dequantize_C;

  if (WP2GetCPUInfo != nullptr) {
#if defined(WP2_USE_SSE)
    if (WP2GetCPUInfo(kSSE)) WP2QuantizeInitSSE();
#endif
  }

  quantize_last_cpuinfo_used = WP2GetCPUInfo;
}
