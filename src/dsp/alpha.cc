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
// alpha-related functions
//
// Author: Skal (pascal.massimino@gmail.com)

#include <cstdint>

#include "src/dsp/dsp.h"

namespace {

template<typename T> bool HasValue_C(const T* src, size_t len, T value) {
  for (size_t i = 0; i < len; ++i) {
    if (src[i] == value) return true;
  }
  return false;
}

template<typename T> bool HasOtherValue_C(const T* src, size_t len, T value) {
  for (size_t i = 0; i < len; ++i) {
    if (src[i] != value) return true;
  }
  return false;
}

bool HasOtherValue8b32b_C(const uint8_t* src, size_t len, uint8_t value) {
  for (size_t i = 0; i < len; ++i) {
    if (src[4 * i + 0] != value) return true;
  }
  return false;
}

}  // namespace

//------------------------------------------------------------------------------

extern void WP2AlphaInitNEON();
extern void WP2AlphaInitSSE();

bool (*WP2HasValue8b)(const uint8_t* src, size_t len, uint8_t value) = nullptr;
bool (*WP2HasValue16b)(const int16_t* src, size_t len, int16_t value) = nullptr;
bool (*WP2HasOtherValue8b)(
    const uint8_t* src, size_t len, uint8_t value) = nullptr;
bool (*WP2HasOtherValue16b)(
    const int16_t* src, size_t len, int16_t value) = nullptr;
bool (*WP2HasOtherValue8b32b)(
    const uint8_t* src, size_t len, uint8_t value) = nullptr;

static volatile WP2CPUInfo alpha_last_cpuinfo_used =
    (WP2CPUInfo)&alpha_last_cpuinfo_used;

WP2_TSAN_IGNORE_FUNCTION void WP2AlphaInit() {
  if (alpha_last_cpuinfo_used == WP2GetCPUInfo) return;

  WP2HasValue8b = HasValue_C<uint8_t>;
  WP2HasValue16b = HasValue_C<int16_t>;
  WP2HasOtherValue8b = HasOtherValue_C<uint8_t>;
  WP2HasOtherValue16b = HasOtherValue_C<int16_t>;
  WP2HasOtherValue8b32b = HasOtherValue8b32b_C;

  if (WP2GetCPUInfo != nullptr) {
#if defined(WP2_USE_SSE)
    if (WP2GetCPUInfo(kSSE)) {
      // TODO(skal): WP2AlphaInitSSE();
    }
#endif
#if defined(WP2_USE_NEON)
    if (WP2GetCPUInfo(kNEON)) {
      // TODO(skal): WP2AlphaInitNEON();
    }
#endif
  }

  alpha_last_cpuinfo_used = WP2GetCPUInfo;
}

//------------------------------------------------------------------------------
