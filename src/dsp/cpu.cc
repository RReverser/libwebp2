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
// CPU detection
//
// Author: James Zern (jzern@google.com)

#include "src/dsp/dsp.h"

#if defined(WP2_HAVE_NEON_RTCD)
#include <cstdio>
#include <cstring>
#endif

#if defined(WP2_ANDROID_NEON)
#include <cpu-features.h>
#endif

//------------------------------------------------------------------------------
// SSExx detection.
//

// apple/darwin gcc-4.0.1 defines __PIC__, but not __pic__ with -fPIC.
#if (defined(__pic__) || defined(__PIC__)) && defined(__i386__)
static inline void GetCPUInfo(int cpu_info[4], int info_type) {
  __asm__ volatile (
    "mov %%ebx, %%edi\n"
    "cpuid\n"
    "xchg %%edi, %%ebx\n"
    : "=a"(cpu_info[0]), "=D"(cpu_info[1]), "=c"(cpu_info[2]), "=d"(cpu_info[3])
    : "a"(info_type), "c"(0));
}
#elif defined(__x86_64__) && \
      (defined(__code_model_medium__) || defined(__code_model_large__)) && \
      defined(__PIC__)
static inline void GetCPUInfo(int cpu_info[4], int info_type) {
  __asm__ volatile (
    "xchg{q}\t{%%rbx}, %q1\n"
    "cpuid\n"
    "xchg{q}\t{%%rbx}, %q1\n"
    : "=a"(cpu_info[0]), "=&r"(cpu_info[1]), "=c"(cpu_info[2]),
      "=d"(cpu_info[3])
    : "a"(info_type), "c"(0));
}
#elif defined(__i386__) || defined(__x86_64__)
static inline void GetCPUInfo(int cpu_info[4], int info_type) {
  __asm__ volatile (
    "cpuid\n"
    : "=a"(cpu_info[0]), "=b"(cpu_info[1]), "=c"(cpu_info[2]), "=d"(cpu_info[3])
    : "a"(info_type), "c"(0));
}
#elif (defined(_M_X64) || defined(_M_IX86)) && \
      defined(_MSC_FULL_VER) && _MSC_FULL_VER >= 150030729  // >= VS2008 SP1
#include <intrin.h>
#define GetCPUInfo(info, type) __cpuidex(info, type, 0)  // set ecx=0
#elif defined(WP2_MSC_SSE)
#define GetCPUInfo __cpuid
#endif

// NaCl has no support for xgetbv or the raw opcode.
#if !defined(__native_client__) && (defined(__i386__) || defined(__x86_64__))
static inline uint64_t xgetbv() {
  const uint32_t ecx = 0;
  uint32_t eax, edx;
  // Use the raw opcode for xgetbv for compatibility with older toolchains.
  __asm__ volatile (
    ".byte 0x0f, 0x01, 0xd0\n"
    : "=a"(eax), "=d"(edx) : "c" (ecx));
  return ((uint64_t)edx << 32) | eax;
}
#elif (defined(_M_X64) || defined(_M_IX86)) && \
      defined(_MSC_FULL_VER) && _MSC_FULL_VER >= 160040219  // >= VS2010 SP1
#include <immintrin.h>
#define xgetbv() _xgetbv(0)
#elif defined(_MSC_VER) && defined(_M_IX86)
static inline uint64_t xgetbv() {
  uint32_t eax_, edx_;
  __asm {
    xor ecx, ecx  // ecx = 0
    // Use the raw opcode for xgetbv for compatibility with older toolchains.
    __asm _emit 0x0f __asm _emit 0x01 __asm _emit 0xd0
    mov eax_, eax
    mov edx_, edx
  }
  return ((uint64_t)edx_ << 32) | eax_;
}
#else
#define xgetbv() 0U  // no AVX for older x64 or unrecognized toolchains.
#endif

#if defined(__i386__) || defined(__x86_64__) || defined(WP2_MSC_SSE)

// helper function for run-time detection of slow SSSE3 platforms
static bool CheckSlowModel(int info) {
  // Table listing display models with longer latencies for the bsr instruction
  // (ie 2 cycles vs 10/16 cycles) and some SSSE3 instructions like pshufb.
  // Refer to Intel 64 and IA-32 Architectures Optimization Reference Manual.
  static const uint8_t kSlowModels[] = {
    0x37, 0x4a, 0x4d,  // Silvermont Microarchitecture
    0x1c, 0x26, 0x27   // Atom Microarchitecture
  };
  const uint32_t model = ((info & 0xf0000) >> 12) | ((info >> 4) & 0xf);
  const uint32_t family = (info >> 8) & 0xf;
  if (family == 0x06) {
    size_t i;
    for (i = 0; i < sizeof(kSlowModels) / sizeof(kSlowModels[0]); ++i) {
      if (model == kSlowModels[i]) return true;
    }
  }
  return false;
}

static bool x86CPUInfo(WP2CPUFeature feature) {
  int max_cpuid_value;
  int cpu_info[4];
  int is_intel = 0;

  // get the highest feature value cpuid supports
  GetCPUInfo(cpu_info, 0);
  max_cpuid_value = cpu_info[0];
  if (max_cpuid_value < 1) {
    return false;
  } else {
    const int VENDOR_ID_INTEL_EBX = 0x756e6547;  // uneG
    const int VENDOR_ID_INTEL_EDX = 0x49656e69;  // Ieni
    const int VENDOR_ID_INTEL_ECX = 0x6c65746e;  // letn
    is_intel = (cpu_info[1] == VENDOR_ID_INTEL_EBX &&
                cpu_info[2] == VENDOR_ID_INTEL_ECX &&
                cpu_info[3] == VENDOR_ID_INTEL_EDX);    // genuine Intel?
  }

  GetCPUInfo(cpu_info, 1);
  if (feature == kSSE2) {
    return !!(cpu_info[3] & (1 << 26));
  }
  if (feature == kSSE3) {
    return !!(cpu_info[2] & (1 << 0));
  }
  if (feature == kSlowSSSE3) {
    if (is_intel && (cpu_info[2] & (1 << 0))) {   // SSSE3?
      return CheckSlowModel(cpu_info[0]);
    }
    return false;
  }

  if (feature == kSSE) {   // both SSE4.1 and SSE4.2 are set
    return ((cpu_info[2] & (3 << 19)) == (3 << 19));
  }
  if (feature == kSSE4_1) {
    return !!(cpu_info[2] & (1 << 19));
  }
  if (feature == kSSE4_2) {
    return !!(cpu_info[2] & (1 << 20));
  }
  if (feature == kAVX) {
    // bits 27 (OSXSAVE) & 28 (256-bit AVX)
    if ((cpu_info[2] & 0x18000000) == 0x18000000) {
      // XMM state and YMM state enabled by the OS.
      return (xgetbv() & 0x6) == 0x6;
    }
  }
  if (feature == kAVX2) {
    if (x86CPUInfo(kAVX) && max_cpuid_value >= 7) {
      GetCPUInfo(cpu_info, 7);
      return !!(cpu_info[1] & (1 << 5));
    }
  }
  return false;
}
WP2CPUInfo WP2GetCPUInfo = x86CPUInfo;
#elif defined(WP2_ANDROID_NEON)  // NB: needs to be before generic NEON test.
static bool AndroidCPUInfo(WP2CPUFeature feature) {
  const AndroidCpuFamily cpu_family = android_getCpuFamily();
  const uint64_t cpu_features = android_getCpuFeatures();
  if (feature == kNEON) {
    return (cpu_family == ANDROID_CPU_FAMILY_ARM &&
            0 != (cpu_features & ANDROID_CPU_ARM_FEATURE_NEON));
  }
  return false;
}
WP2CPUInfo WP2GetCPUInfo = AndroidCPUInfo;
#elif defined(EMSCRIPTEN) // also needs to be before generic NEON test
// Use compile flags as an indicator of SIMD support instead of a runtime check.
static bool wasmCPUInfo(WP2CPUFeature feature) {
  switch (feature) {
    case kSSE2:
    case kSSE3:
    case kSlowSSSE3:
    case kSSE4_1:
    case kSSE4_2:
    case kSSE:
#ifdef WP2_USE_SSE
      return true;
#endif
      return false;
    case kAVX:
    case kAVX2:
#ifdef WP2_USE_AVX2
      return true;
#endif
      return false;
    case kNEON:
#ifdef WP2_USE_NEON
      return true;
#endif
      return false;
    default:
      return false;
  }
}
WP2CPUInfo WP2GetCPUInfo = wasmCPUInfo;
#elif defined(WP2_USE_NEON)
// define a dummy function to enable turning off NEON at runtime by setting
// WP2DecGetCPUInfo = NULL
static bool armCPUInfo(WP2CPUFeature feature) {
  if (feature != kNEON) return false;
#if defined(__linux__) && defined(WP2_HAVE_NEON_RTCD)
  {
    bool has_neon = false;
    char line[200];
    FILE* const cpuinfo = fopen("/proc/cpuinfo", "r");
    if (cpuinfo == NULL) return false;
    while (fgets(line, sizeof(line), cpuinfo)) {
      if (!strncmp(line, "Features", 8)) {
        if (strstr(line, " neon ") != NULL) {
          has_neon = true;
          break;
        }
      }
    }
    fclose(cpuinfo);
    return has_neon;
  }
#else
  return true;
#endif
}
WP2CPUInfo WP2GetCPUInfo = armCPUInfo;
#elif defined(WP2_USE_MIPS32) || defined(WP2_USE_MIPS_DSP_R2) || \
      defined(WP2_USE_MSA)
static bool mipsCPUInfo(WP2CPUFeature feature) {
  if ((feature == kMIPS32) || (feature == kMIPSdspR2) || (feature == kMSA)) {
    return true;
  } else {
    return false;
  }
}
WP2CPUInfo WP2GetCPUInfo = mipsCPUInfo;
#else
WP2CPUInfo WP2GetCPUInfo = nullptr;
#endif

// undocumented function, for testing proper calls to WP2xxxInit()
extern void WP2DspReset();

// this function is just here to offer a different address than GetCPUInfo,
// so that calling WP2xxxInit() really resets everything.
static bool FakeGetInfo(WP2CPUFeature feature) { return false; }

void WP2DspReset() {
  // first time, we save the original WP2GetCPUInfo, for later reset
  static WP2CPUInfo saved_cpu_info = nullptr;
  if (saved_cpu_info == nullptr) saved_cpu_info = WP2GetCPUInfo;

  // first we clear the local xxx_last_cpuinfo_used internal pointers
  WP2GetCPUInfo = FakeGetInfo;
  WP2EncDspInit();
  WP2DecDspInit();
  WP2ArgbConverterInit();
  WP2TransformInit();
  WP2::DblkFilterInit();
  WP2::DrctFilterInit();
  WP2::WienerFilterInit();
  WP2::GrainFilterInit();
  WP2::PredictionInit();
  WP2::ScoreDspInit();
  WP2SSIMInit();
  WP2PSNRInit();
  WP2AlphaInit();
  WP2QuantizeInit();
  WP2::ANSInit();

  // then, we reset the pointer-to-functions for good
  WP2GetCPUInfo = saved_cpu_info;

  for (auto& f : WP2ArgbConvertFrom) f = nullptr;
  for (auto& f : WP2ArgbConvertTo) f = nullptr;

  WP2YUVToCustom = nullptr;
  WP2YUVToXRGB = nullptr;
  WP2YUVToArgb = nullptr;
  WP2YUVToARGB = nullptr;

  WP2Transpose = nullptr;
  for (uint32_t i = 0; i < 5; ++i) {
    WP2InvDct[i] = WP2InvAdst[i] = WP2InvHadamard[i] = nullptr;
    WP2FwdDct[i] = WP2FwdAdst[i] = WP2FwdHadamard[i] = nullptr;
  }
  WP2SlowDct8x8 = nullptr;

  WP2::GetBlockMinMax = nullptr;
  WP2::GetBlockMinMax_5x5 = nullptr;

  WP2::DeblockLine = nullptr;
  WP2::WouldDeblockLine = nullptr;
  WP2::MeasureFlatLength = nullptr;

  WP2::CdefDirection4x4 = nullptr;
  WP2::CdefDirection8x8 = nullptr;
  WP2::CdefPad = nullptr;
  WP2::CdefFiltering = nullptr;

  WP2::WienerFilter = nullptr;

  WP2::AddGrain4x4 = nullptr;
  WP2::GenerateGrain4x4 = nullptr;

  for (auto& f : WP2::BasePredictors) f = nullptr;
  WP2::AnglePredInterpolate = nullptr;
  WP2::SubtractRow = nullptr;
  WP2::AddRow = nullptr;
  for (auto& f : WP2::SubtractBlock) f = nullptr;
  for (auto& f : WP2::AddBlock) f = nullptr;

  WP2SSIMGet4x8u = nullptr;
  WP2SSIMGet8u = nullptr;
  WP2SSIMGet16s = nullptr;

  WP2SSIMGetClipped4x8u = nullptr;
  WP2SSIMGetClipped8u = nullptr;
  WP2SSIMGetClipped16s = nullptr;

  WP2SumSquaredError8u = nullptr;
  WP2SumSquaredError16s = nullptr;
  WP2SumSquaredError4x8u = nullptr;

  WP2HasValue8b = nullptr;
  WP2HasValue16b = nullptr;
  WP2HasOtherValue8b = nullptr;
  WP2HasOtherValue16b = nullptr;
  WP2HasOtherValue8b32b = nullptr;

  WP2::ANSUpdateCDF = nullptr;

  WP2Quantize = nullptr;
  WP2Dequantize = nullptr;
}
