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
// Grain generation
//
// Author: Skal (pascal.massimino@gmail.com)

#include "src/dsp/dsp.h"
#include "src/dsp/math.h"

namespace {

void GenerateGrain4x4_C(WP2::PseudoRNG* const rng,
                        uint32_t amp, uint32_t cut_off, int32_t dst[]) {
  for (uint32_t j = 0; j < 4; ++j) {
    for (uint32_t i = 0; i < 4; ++i) {
      if (i + j > cut_off && i + j < 3 * cut_off) {
        dst[i + j * 4] = rng->GetSigned(amp) / (int)(i + j - cut_off);
      } else {
        dst[i + j * 4] = 0;
      }
    }
  }
  WP2InvTransform2D(dst, kDct, kDct, 4, 4, dst);
}

void AddGrain4x4_C(int16_t* samples, size_t step, WP2::PseudoRNG* const rng,
                   uint32_t amp, uint32_t cut_off) {
  if (amp == 0) return;
  int32_t tmp[4 * 4];
  GenerateGrain4x4_C(rng, amp, cut_off, tmp);
  const int32_t max_value = (1 << WP2::kMaxYuvBits) - 1;
  for (uint32_t j = 0; j < 4; ++j) {
    for (uint32_t i = 0; i < 4; ++i) {
      samples[i] = WP2::Clamp(samples[i] + tmp[i + j * 4],
                              -max_value, max_value);
    }
    samples += step;
  }
}

}  // namespace

namespace WP2 {

//------------------------------------------------------------------------------

AddGrainF AddGrain4x4 = nullptr;
GenerateGrainF GenerateGrain4x4 = nullptr;

static volatile WP2CPUInfo grain_filter_last_cpuinfo_used =
    (WP2CPUInfo)&grain_filter_last_cpuinfo_used;

WP2_TSAN_IGNORE_FUNCTION void GrainFilterInit() {
  if (grain_filter_last_cpuinfo_used == WP2GetCPUInfo) return;

  WP2TransformInit();
  GenerateGrain4x4 = GenerateGrain4x4_C;
  AddGrain4x4 = AddGrain4x4_C;

  grain_filter_last_cpuinfo_used = WP2GetCPUInfo;
}

//------------------------------------------------------------------------------

}  // namespace WP2
