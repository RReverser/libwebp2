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
// Image transforms and color space conversion methods for lossless decoder.
//
// Authors: Vincent Rabaud (vrabaud@google.com)

#ifndef WP2_DSP_LOSSLESS_LOSSLESS_COMMON_H_
#define WP2_DSP_LOSSLESS_LOSSLESS_COMMON_H_

#include <cmath>

#include "src/common/lossless/color_cache.h"
#include "src/dsp/dsp.h"
#include "src/dsp/math.h"
#include "src/utils/utils.h"

//------------------------------------------------------------------------------
// Misc methods.

// Computes sampled size of 'size' when sampling using 'sampling bits'.
// Equivalent to: ceil(size / 2^sampling_bits)
static inline uint32_t WP2LSubSampleSize(uint32_t size,
                                         uint32_t sampling_bits) {
  return (size + (1 << sampling_bits) - 1) >> sampling_bits;
}

// -----------------------------------------------------------------------------
// PrefixEncode()

// Splitting of distance and length codes into prefixes and
// extra bits. The prefixes are encoded with an entropy code
// while the extra bits are stored just as normal bits.
static inline void WP2LPrefixEncodeBitsNoLUT(int distance, int* const code,
                                                  int* const extra_bits) {
  const int highest_bit = WP2Log2Floor(--distance);
  const int second_highest_bit = (distance >> (highest_bit - 1)) & 1;
  *extra_bits = highest_bit - 1;
  *code = 2 * highest_bit + second_highest_bit;
}

static inline void WP2LPrefixEncodeNoLUT(int distance, int* const code,
                                              int* const extra_bits,
                                              int* const extra_bits_value) {
  const int highest_bit = WP2Log2Floor(--distance);
  const int second_highest_bit = (distance >> (highest_bit - 1)) & 1;
  *extra_bits = highest_bit - 1;
  *extra_bits_value = distance & ((1 << *extra_bits) - 1);
  *code = 2 * highest_bit + second_highest_bit;
}

#define PREFIX_LOOKUP_IDX_MAX   512
typedef struct {
  int8_t code_;
  int8_t extra_bits_;
} WP2LPrefixCode;

// These tables are derived using WP2LPrefixEncodeNoLUT.
extern const WP2LPrefixCode kWP2PrefixEncodeCode[PREFIX_LOOKUP_IDX_MAX];
extern const uint8_t kWP2PrefixEncodeExtraBitsValue[PREFIX_LOOKUP_IDX_MAX];
static inline void WP2LPrefixEncodeBits(int distance, int* const code,
                                             int* const extra_bits) {
  if (distance < PREFIX_LOOKUP_IDX_MAX) {
    const WP2LPrefixCode prefix_code = kWP2PrefixEncodeCode[distance];
    *code = prefix_code.code_;
    *extra_bits = prefix_code.extra_bits_;
  } else {
    WP2LPrefixEncodeBitsNoLUT(distance, code, extra_bits);
  }
}

static inline void WP2LPrefixEncode(int distance, int* const code,
                                         int* const extra_bits,
                                         int* const extra_bits_value) {
  if (distance < PREFIX_LOOKUP_IDX_MAX) {
    const WP2LPrefixCode prefix_code = kWP2PrefixEncodeCode[distance];
    *code = prefix_code.code_;
    *extra_bits = prefix_code.extra_bits_;
    *extra_bits_value = kWP2PrefixEncodeExtraBitsValue[distance];
  } else {
    WP2LPrefixEncodeNoLUT(distance, code, extra_bits, extra_bits_value);
  }
}

// Sum of each component, mod 256.
static WP2_UBSAN_IGNORE_UNSIGNED_OVERFLOW inline void WP2LAddPixels(
    const uint16_t* const a, bool has_alpha, const uint16_t* const b,
    uint16_t mask, uint16_t* const out) {
  out[0] = has_alpha ? ((a[0] + b[0]) & mask) : WP2::kAlphaMax;
  for (size_t i = 1; i < 4; ++i) out[i] = (a[i] + b[i]) & mask;
}

// Difference of each component, mod 256.
static inline void WP2LSubPixels(const uint16_t* const a, bool has_alpha,
                                 const uint16_t* const b, uint16_t mask,
                                 uint16_t* const out) {
  out[0] = has_alpha ? ((a[0] - b[0]) & mask) : WP2::kAlphaMax;
  for (size_t i = 1; i < 4; ++i) out[i] = (a[i] - b[i]) & mask;
}

//------------------------------------------------------------------------------
// Transform-related functions use din both encoding and decoding.

// Macros used to create a batch predictor that iteratively uses a
// one-pixel predictor.

// The predictor is added to the output pixel (which
// is therefore considered as a residual) to get the final prediction.
#define GENERATE_PREDICTOR_ADD(PREDICTOR, PREDICTOR_ADD)                      \
  static void PREDICTOR_ADD(const uint16_t* in, bool has_alpha,               \
                            const uint16_t* upper, uint32_t num_pixels,       \
                            uint32_t num_bits, uint16_t* out) {               \
    uint16_t pred[4];                                                         \
    const uint16_t mask = (1u << num_bits) - 1;                               \
    assert(upper != nullptr);                                                 \
    for (const uint16_t* const out_end = out + 4 * num_pixels; out < out_end; \
         out += 4, in += 4, upper += 4) {                                     \
      (PREDICTOR)(&out[-4], upper, mask, pred);                               \
      WP2LAddPixels(in, has_alpha, pred, mask, out);                          \
    }                                                                         \
  }

// It subtracts the prediction from the input pixel and stores the residual
// in the output pixel.
#define GENERATE_PREDICTOR_SUB(PREDICTOR, PREDICTOR_SUB)                      \
  static void PREDICTOR_SUB(const uint16_t* in, bool has_alpha,               \
                            const uint16_t* upper, uint32_t num_pixels,       \
                            uint32_t num_bits, uint16_t* out) {               \
    uint16_t pred[4];                                                         \
    const uint16_t mask = (1u << num_bits) - 1;                               \
    assert(upper != nullptr);                                                 \
    for (const uint16_t* const out_end = out + 4 * num_pixels; out < out_end; \
         out += 4, in += 4, upper += 4) {                                     \
      (PREDICTOR)(&in[-4], upper, mask, pred);                                \
      WP2LSubPixels(in, has_alpha, pred, mask, out);                          \
    }                                                                         \
  }

//------------------------------------------------------------------------------
// Image transforms.

static inline void Average2(const uint16_t* const a0, const uint16_t* const a1,
                            uint16_t mask, uint16_t* const out) {
  for (size_t i = 0; i < 4; ++i) out[i] = ((a0[i] + a1[i]) >> 1) & mask;
}

static inline void Average3(const uint16_t* const a0, const uint16_t* const a1,
                            const uint16_t* const a2, uint16_t mask,
                            uint16_t* const out) {
  Average2(a0, a2, mask, out);
  Average2(out, a1, mask, out);
}

static inline void Average4(const uint16_t* const a0, const uint16_t* const a1,
                            const uint16_t* const a2, const uint16_t* const a3,
                            uint16_t mask, uint16_t* const out) {
  for (size_t i = 0; i < 4; ++i) {
    const uint16_t tmp0 = ((a0[i] + a1[i]) >> 1) & mask;
    const uint16_t tmp1 = ((a2[i] + a3[i]) >> 1) & mask;
    out[i] = ((tmp0 + tmp1) >> 1) & mask;
  }
}

static inline int AddSubtractComponentFull(int a, int b, int c, int max_value) {
  return WP2::Clamp(a + b - c, 0, max_value);
}

static inline void ClampedAddSubtractFull(const uint16_t* const c0,
                                          const uint16_t* const c1,
                                          const uint16_t* const c2,
                                          uint16_t mask,
                                          uint16_t* const out) {
  for (size_t i = 0; i < 4; ++i) {
    out[i] = AddSubtractComponentFull(c0[i], c1[i], c2[i], mask) & mask;
  }
}

static inline int AddSubtractComponentHalf(int a, int b, int max_value) {
  return WP2::Clamp(a + (a - b) / 2, 0, max_value);
}

static inline void ClampedAddSubtractHalf(const uint16_t* const c0,
                                          const uint16_t* const c1,
                                          const uint16_t* const c2,
                                          uint16_t mask, uint16_t* const out) {
  Average2(c0, c1, mask, out);
  for (size_t i = 0; i < 4; ++i) {
    out[i] = AddSubtractComponentHalf(out[i], c2[i], mask) & mask;
  }
}

// gcc-4.9 on ARM generates incorrect code in Select() when Sub3() is inlined.
#if defined(__arm__) && \
    (LOCAL_GCC_VERSION == 0x409 || LOCAL_GCC_VERSION == 0x408)
# define LOCAL_INLINE __attribute__ ((noinline))
#else
# define LOCAL_INLINE inline
#endif

static LOCAL_INLINE int Sub3(int a, int b, int c) {
  const int pb = b - c;
  const int pa = a - c;
  return std::abs(pb) - std::abs(pa);
}

#undef LOCAL_INLINE

static inline void Select(const uint16_t* const a, const uint16_t* const b,
                          const uint16_t* const c, uint16_t* const out) {
  int pa_minus_pb = 0;
  for (size_t i = 0; i < 4; ++i) pa_minus_pb += Sub3(a[i], b[i], c[i]);
  if (pa_minus_pb <= 0) {
    for (size_t i = 0; i < 4; ++i) out[i] = a[i];
  } else {
    for (size_t i = 0; i < 4; ++i) out[i] = b[i];
  }
}

#endif  // WP2_DSP_LOSSLESS_LOSSLESS_COMMON_H_
