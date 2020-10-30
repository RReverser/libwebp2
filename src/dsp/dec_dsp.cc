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
// Speed-critical decoding functions, default plain-C implementations.
//
// Author: Skal (pascal.massimino@gmail.com)

#include "src/dsp/dsp.h"

#if defined(WP2_BITTRACE)   // we'll need some WP2Log2() calls for bit-trace
#include "src/dsp/math.h"
#endif

//------------------------------------------------------------------------------

static volatile WP2CPUInfo dec_last_cpuinfo_used =
    (WP2CPUInfo)&dec_last_cpuinfo_used;

WP2_TSAN_IGNORE_FUNCTION void WP2DecDspInit() {
  if (dec_last_cpuinfo_used == WP2GetCPUInfo) return;

#if defined(WP2_BITTRACE)
  WP2MathInit();
#endif
  WP2TransformInit();
  WP2QuantizeInit();
  WP2::ANSInit();
  WP2::PredictionInit();

  dec_last_cpuinfo_used = WP2GetCPUInfo;
}

//------------------------------------------------------------------------------
