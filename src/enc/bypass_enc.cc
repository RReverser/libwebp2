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
//  Bypass with raw pixel coding
//
// Author: Yannis Guyon (yguyon@google.com)

#include "src/dec/tile_dec.h"
#include "src/dsp/math.h"
#include "src/enc/wp2_enc_i.h"
#include "src/utils/utils.h"
#include "src/wp2/base.h"
#include "src/wp2/encode.h"
#include "src/wp2/format_constants.h"

namespace WP2 {

//------------------------------------------------------------------------------

WP2Status BypassTileEnc(const GlobalParams& gparams, EncTile* const tile,
                        Writer* const output) {
  const uint32_t max_num_bytes = GetTileMaxNumBytes(gparams, tile->rect);

  if (gparams.type_ == GlobalParams::GP_LOSSLESS) {
    // RGB input is expected.
    assert(!tile->rgb_input.IsEmpty());
    const uint8_t* row = (const uint8_t*)tile->rgb_input.GetRow(0);
    for (uint32_t y = 0; y < tile->rect.height; ++y) {
      if (gparams.has_alpha_) {
        // Write rows of Argb samples.
        WP2_CHECK_ALLOC_OK(output->Append((void*)row, tile->rect.width * 4));
      } else {
        // Write only rgb.
        for (uint32_t x = 0; x < tile->rect.width; ++x) {
          WP2_CHECK_ALLOC_OK(output->Append((void*)&row[x * 4 + 1], 3));
        }
      }
      row += tile->rgb_input.stride;
    }
    // TODO(yguyon): Simulate ANSDebugPrefix to count these raw bytes?
  } else {
    // YUV input is expected even though RGB may also be available.
    assert(!tile->yuv_input.IsEmpty());
    const uint32_t channel_bits = gparams.transf_.GetYUVPrecisionBits() + 1u;

    // Allocate, fill and write the 'data' to 'output'.
    // TODO(yguyon): Use less memory by outputting 8 lines by 8 lines?
    Data data;
    WP2_CHECK_STATUS(data.Resize(max_num_bytes, /*keep_bytes=*/false));
    // Use a HeaderEnc instance to easily pack sample bits (e.g. 10 by 10).
    HeaderEnc bit_packer(data.bytes, data.size);
    bit_packer.SetDebugPrefix("raw_px");

    for (uint32_t y = 0; y < tile->rect.height; ++y) {
      // Interleave channels to output samples line by line.
      for (Channel channel : {kYChannel, kUChannel, kVChannel}) {
        const int16_t* const row = tile->yuv_input.GetChannel(channel).Row(y);
        for (uint32_t x = 0; x < tile->rect.width; ++x) {
          bit_packer.PutSBits((int32_t)row[x], channel_bits, "yuv");
        }
      }
      if (gparams.has_alpha_) {
        const int16_t* const row = tile->yuv_input.A.Row(y);
        for (uint32_t x = 0; x < tile->rect.width; ++x) {
          assert(row[x] >= 0);
          bit_packer.PutBits((uint32_t)row[x], kAlphaBits, "alpha");
        }
      }
    }

    WP2_CHECK_OK(bit_packer.Ok(), WP2_STATUS_BAD_WRITE);
    assert(bit_packer.Used() == max_num_bytes);
    WP2_CHECK_ALLOC_OK(output->Append(data.bytes, data.size));
  }
  return WP2_STATUS_OK;
}

//------------------------------------------------------------------------------

}  // namespace WP2
