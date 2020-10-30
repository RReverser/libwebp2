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
// Wiener filter based on AV1 implementation at
// https://code.videolan.org/videolan/dav1d/blob/master/src/
// Attenuates the noise.
//
// Author: Yannis Guyon (yguyon@google.com)

#include <algorithm>
#include <cassert>
#include <cstring>
#include <limits>

#include "src/dsp/dsp.h"
#include "src/dsp/math.h"
#include "src/utils/vector.h"

namespace WP2 {

//------------------------------------------------------------------------------

static inline void Copy(const int16_t* src, size_t src_step,
                        int32_t width, int32_t height,
                        size_t dst_step, int16_t* dst) {
  for (int32_t y = 0; y < height; ++y, dst += dst_step, src += src_step) {
    std::memcpy(dst, src, width * sizeof(src[0]));
  }
}

void WienerHalfToFullWgts(const int32_t half_tap_weights[kWieFltTapDist],
                          int32_t full_tap_weights[kWieFltNumTaps]) {
  // unit vector
  full_tap_weights[kWieFltTapDist] = 1 << (kWieFltNumBitsTapWgts - 1);
  for (uint32_t i = 0; i < kWieFltTapDist; ++i) {
    full_tap_weights[i] = half_tap_weights[i];
    full_tap_weights[kWieFltNumTaps - 1u - i] = half_tap_weights[i];  // symm.
    full_tap_weights[kWieFltTapDist] -= 2 * half_tap_weights[i];  // unit vector
  }
}

//------------------------------------------------------------------------------

// Fills the extended input buffer with unfiltered pixels, stretched to the
// sides for areas beyond available source.
// Area is sliced into corners, sides, and center areas as follows ->  ABBBC
// 'tmp[0]' and 'src[0]' represent the top-left E.                     DEEEF
// 'n_*' should be at most kWieFltTapDist.                             GHHHI
static void Pad(int32_t width, int32_t height,
                const int16_t* const src, size_t src_step,
                const int16_t* const left, size_t left_step,
                const int16_t* const right, size_t right_step,
                const int16_t* const top, size_t top_step,
                const int16_t* const bottom, size_t bottom_step,
                int32_t n_left, int32_t n_right,
                int32_t n_top, int32_t n_bottom,
                size_t tmp_step, int16_t* const tmp) {
  static constexpr int32_t gap = (int32_t)kWieFltTapDist;
  assert(width > 0 && height > 0);
  assert(n_left <= gap && n_right <= gap && n_top <= gap && n_bottom <= gap);
  int32_t x_start = -n_left, x_end = width + n_right;
  int32_t y_start = -n_top, y_end = height + n_bottom;

  // Copy pixels as far as available.
  if (n_top > 0) {  // [A]B[C]
    Copy(top, top_step, /*width=*/x_end - x_start, /*height=*/n_top,
         tmp_step, /*dst=*/tmp - gap * (int32_t)tmp_step + x_start);
  }
  if (n_left > 0) {  // D
    Copy(left, left_step, /*width=*/n_left, /*height=*/height,
         tmp_step, /*dst=*/tmp + x_start);
  }
  {  // E
    Copy(src, src_step, /*width=*/width, /*height=*/height,
         tmp_step, /*dst=*/tmp);
  }
  if (n_right> 0) {  // F
    Copy(right, right_step, /*width=*/n_right, /*height=*/height,
         tmp_step, /*dst=*/tmp + width);
  }
  if (n_bottom > 0) {  // [G]H[I]
    Copy(bottom, bottom_step, /*width=*/x_end - x_start, /*height=*/n_bottom,
         tmp_step, /*dst=*/tmp + height * (int32_t)tmp_step + x_start);
  }

  // Stretch [A]B[C] from [D]E[F].
  for (int32_t y = -gap; y < y_start; ++y) {
    Copy(/*src=*/tmp + y_start * (int)tmp_step + x_start, tmp_step,
         /*width=*/x_end - x_start, /*height=*/1,
         tmp_step, /*dst=*/tmp + y * (int)tmp_step + x_start);
  }
  // Stretch [G]H[I] from [D]E[F].
  for (int32_t y = height + gap - 1; y >= y_end; --y) {
    Copy(/*src=*/tmp + (y_end - 1) * (int)tmp_step + x_start, tmp_step,
         /*width=*/x_end - x_start, /*height=*/1,
         tmp_step, /*dst=*/tmp + y * (int)tmp_step + x_start);
  }
  // Stretch ADG from BEH.
  for (int32_t x = -gap; x < x_start; ++x) {
    Copy(/*src=*/tmp - gap * (int)tmp_step + x_start, tmp_step,
         /*width=*/1, /*height=*/gap + height + gap,
         tmp_step, /*dst=*/tmp - gap * (int)tmp_step + x);
  }
  // Stretch CFI from BEH.
  for (int32_t x = width + gap - 1; x >= x_end; --x) {
    Copy(/*src=*/tmp - gap * (int)tmp_step + x_end - 1, tmp_step,
         /*width=*/1, /*height=*/gap + height + gap,
         tmp_step, /*dst=*/tmp - gap * (int)tmp_step + x);
  }
}

//------------------------------------------------------------------------------

static inline int32_t Lerp(int32_t from, int32_t to, int32_t t,
                           uint32_t t_num_bits) {
  const int32_t t_max = 1 << t_num_bits;
  if (t == 0) return from;
  if (t == t_max) return to;
  assert(t >= 0 && t <= t_max);
  return RightShiftRound(from * (t_max - t) + to * t, t_num_bits);
}

static WP2Status WienerFilter_C(
    uint32_t width, uint32_t height, const int16_t* const left,
    size_t left_step, const int16_t* const right, size_t right_step,
    const int16_t* const top, size_t top_step, const int16_t* const bottom,
    size_t bottom_step, uint32_t n_left, uint32_t n_right, uint32_t n_top,
    uint32_t n_bottom, const int32_t tap_weights_h[kWieFltNumTaps],
    const int32_t tap_weights_v[kWieFltNumTaps],
    const uint8_t* const strength_map, size_t strength_step,
    uint32_t num_precision_bits, size_t dst_step, int16_t* const dst) {
  const uint32_t filter_width = std::min(width, kWieFltWidth);
  const uint32_t filter_height = std::min(height, kWieFltHeight);
  const uint32_t buf_width = filter_width + 2 * kWieFltTapDist;
  const uint32_t buf_height = filter_height + 2 * kWieFltTapDist;

  // TODO(maryla): allocating this buffer all the time is pretty slow.
  Vector_s16 tmp;
  WP2_CHECK_ALLOC_OK(tmp.resize(buf_width * buf_height));
  const size_t tmp_step = width + 2 * kWieFltTapDist;
  const uint32_t num_precision_bits_tmp = sizeof(tmp[0]) * 8;

  Pad((int32_t)width, (int32_t)height, dst, dst_step, left, left_step, right,
      right_step, top, top_step, bottom, bottom_step, (int32_t)n_left,
      (int32_t)n_right, (int32_t)n_top, (int32_t)n_bottom, tmp_step,
      &tmp[0] + kWieFltTapDist * tmp_step + kWieFltTapDist);

  // Temporary buffer containing horizontal sums.
  Vector_s16 hor;
  WP2_CHECK_ALLOC_OK(hor.resize(tmp.size()));

  // Filter horizontally.
  const uint32_t num_precision_bits_h =
      num_precision_bits + kWieFltNumBitsTapWgts - 1;  // Both include sign.
  const uint32_t shift_h = SafeSub(
      num_precision_bits_h + kWieFltNumBitsOverflow, num_precision_bits_tmp);
  const int32_t round_h = shift_h / 2;

  const int16_t* tmp_row = &tmp[0];
  int16_t* hor_row = &hor[0];
  for (uint32_t j = 0; j < height + 2 * kWieFltTapDist; j++) {
    for (uint32_t i = 0; i < width; i++) {
      int32_t sum = 0;
      for (uint32_t k = 0; k < kWieFltNumTaps; k++) {
        sum += tmp_row[i + k] * tap_weights_h[k];
      }

      sum = RightShift(sum + round_h, shift_h);
      assert(sum >= std::numeric_limits<int16_t>::min() &&
             sum <= std::numeric_limits<int16_t>::max());
      hor_row[i] = (int16_t)sum;
    }
    tmp_row += tmp_step;
    hor_row += tmp_step;
  }

  // Filter vertically.
  const uint32_t num_precision_bits_v =
      num_precision_bits_h + kWieFltNumBitsTapWgts - 1;  // Both include sign.
  const uint32_t shift_v =
      SafeSub(num_precision_bits_v - shift_h, num_precision_bits);
  const int32_t round_v = shift_v / 2;

  for (uint32_t i = 0; i < width; i++) {
    for (uint32_t j = 0; j < height; j++) {
      const uint8_t strength = strength_map[j * strength_step + i];
      if (strength > 0) {
        int32_t sum = 0;
        for (uint32_t k = 0; k < kWieFltNumTaps; k++) {
          sum += hor[(j + k) * tmp_step + i] * tap_weights_v[k];
        }

        int16_t* const pixel = &dst[j * dst_step + i];
        // Verify input is within 'num_precision_bits'.
        assert(*pixel == ClampToSigned(*pixel, num_precision_bits));

        const int32_t filtered_pixel = RightShift(sum + round_v, shift_v);
        // 'strength' is at most 63 so add 1 for an easier right shift.
        assert(strength <= 63);
        const int32_t final_pixel =
            Lerp(*pixel, filtered_pixel, strength + 1, 6);
        // Values can overflow of 'kWieFltNumBitsOverflow' hence the clamping.
        *pixel = (int16_t)ClampToSigned(final_pixel, num_precision_bits);
      }
    }
  }

  return WP2_STATUS_OK;
}

//------------------------------------------------------------------------------

WienerFilterF WienerFilter = nullptr;

static volatile WP2CPUInfo wiener_filter_last_cpuinfo_used =
    (WP2CPUInfo)&wiener_filter_last_cpuinfo_used;

WP2_TSAN_IGNORE_FUNCTION void WienerFilterInit() {
  if (wiener_filter_last_cpuinfo_used == WP2GetCPUInfo) return;

  WienerFilter = WienerFilter_C;

  wiener_filter_last_cpuinfo_used = WP2GetCPUInfo;
}

//------------------------------------------------------------------------------

}  // namespace WP2
