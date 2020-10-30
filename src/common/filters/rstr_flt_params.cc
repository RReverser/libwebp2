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
// Restoration filter parameters.
// Optimized during encoding and read from the bitstream during decoding.
//
// Author: Yannis Guyon (yguyon@google.com)

#include "src/common/filters/rstr_flt_params.h"

#include <algorithm>
#include <cassert>

namespace WP2 {

//------------------------------------------------------------------------------

RstrFltParams::RstrFltParams(uint32_t width, uint32_t height)
    : num_areas(DivCeil(width, kWieFltWidth) * DivCeil(height, kWieFltHeight)) {
  assert(num_areas <= kWieMaxNumAreas);
}

void RstrFltParams::SetAll(const int8_t values[kWieFltTapDist]) {
  for (Channel channel : {kYChannel, kUChannel, kVChannel}) {
    for (uint32_t area = 0; area < num_areas; ++area) {
      for (uint32_t h_or_v : {0, 1}) {
        std::copy(values, values + kWieFltTapDist,
                  half_tap_weights[channel][area][h_or_v]);
      }
    }
  }
}

//------------------------------------------------------------------------------

bool WienerVerifyHalfTapWeights(
    const int8_t half_tap_weights[kWieFltTapDist]) {
  for (uint32_t i = 0; i < kWieFltTapDist; ++i) {
    if (half_tap_weights[i] < GetMinWienerTapWeight(i) ||
        half_tap_weights[i] > GetMaxWienerTapWeight(i)) {
      return false;
    }
  }

  // Make sure during compilation that if it fits within 'kWieNumBitsPerTapWgt',
  // it won't take more than 'kWieFltNumBitsOverflow' for computation.
  static_assert(WP2Log2Ceil_k((1u << (kWieFltNumBitsTapWgts - 1)) +
                              4 * ((1u << (kWieNumBitsPerTapWgt[0] - 1)) +
                                   (1u << (kWieNumBitsPerTapWgt[1] - 1)) +
                                   (1u << (kWieNumBitsPerTapWgt[2] - 1)))) <=
                    kWieFltNumBitsTapWgts - 1 + kWieFltNumBitsOverflow,
                "Overflowing kWieNumBitsPerTapWgt");
  return true;
}

//------------------------------------------------------------------------------

}    // namespace WP2
