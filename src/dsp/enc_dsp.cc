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
// Speed-critical encoding functions, default plain-C implementations.
//
// Author: Skal (pascal.massimino@gmail.com)

#include "src/dsp/dsp.h"

//------------------------------------------------------------------------------

extern void WP2EncDspInitSSE();
extern void WP2EncDspInitNEON();

static volatile WP2CPUInfo enc_last_cpuinfo_used =
    (WP2CPUInfo)&enc_last_cpuinfo_used;

WP2_TSAN_IGNORE_FUNCTION void WP2EncDspInit() {
  if (enc_last_cpuinfo_used == WP2GetCPUInfo) return;

  WP2DecDspInit();
  WP2PSNRInit();   // used by GetDisto, etc.
  WP2AlphaInit();
  WP2QuantizeInit();

  // If defined, use CPUInfo() to overwrite some pointers with faster versions.
  if (WP2GetCPUInfo != nullptr) {
#if defined(WP2_USE_SSE)
    if (WP2GetCPUInfo(kSSE)) {
      // WP2EncDspInitSSE();
    }
#endif
#if defined(WP2_USE_NEON)
    if (WP2GetCPUInfo(kNEON)) {
      // WP2EncDspInitNEON();
    }
#endif
  }
  enc_last_cpuinfo_used = WP2GetCPUInfo;
}

//------------------------------------------------------------------------------
