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
// WP2 lossless encoding.
//

#include "src/enc/lossless/losslessi_enc.h"
#include "src/enc/wp2_enc_i.h"

namespace WP2 {

// -----------------------------------------------------------------------------

WP2Status TileEncoder::LosslessEncode(ANSEnc* const enc) {
  ANSDictionaries dicts;
  // If R,G,B fit on 10 bits, make sure the values do not exceed 8 bits.
  // TODO(vrabaud) remove
  if (tile_->rgb_input.format == WP2_Argb_38) {
    for (uint32_t y = 0; y < tile_->rgb_input.height; ++y) {
      const uint16_t* const row = tile_->rgb_input.GetRow16(y);
      for (uint32_t x = 0; x < tile_->rgb_input.width; ++x) {
        for (uint32_t c = 1; c < 4; ++c) {
          WP2_CHECK_OK(row[4 * x + c] < (1 << 8),
                       WP2_STATUS_INVALID_COLORSPACE);
        }
      }
    }
  }
  WP2_CHECK_STATUS(WP2L::EncodeImage(*config_, tile_->rgb_input,
                                     tiles_layout_->gparams->has_alpha_, enc));
  return enc->GetStatus();
}

// -----------------------------------------------------------------------------

}    // namespace WP2
