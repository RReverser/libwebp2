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
// WP2 lossless decoding.
//
// Author: Skal (pascal.massimino@gmail.com)

#include "src/dec/wp2_dec_i.h"

#include "src/dec/lossless/losslessi_dec.h"

namespace WP2 {

//------------------------------------------------------------------------------
// Lossless decoding.

WP2Status LosslessDecode(const BitstreamFeatures& features,
                         const DecoderConfig& config,
                         const GlobalParams& gparams,
                         ProgressWatcher* const progress, ANSDec* const dec,
                         Tile* const tile) {
  (void)features;
  assert(!tile->rgb_output.IsEmpty());
  tile->output_is_yuv = false;

  WP2L::Decoder dec_lossless;
  dec_lossless.Init(config, gparams, progress, dec, tile);

  WP2_CHECK_STATUS(dec_lossless.DecodeHeader());
  WP2_CHECK_STATUS(dec_lossless.DecodeImage());

  WP2_CHECK_STATUS(dec->GetStatus());
  // Not using whole chunk is an issue.
  WP2_CHECK_OK(dec->GetNumUsedBytes() == tile->chunk_size,
               WP2_STATUS_BITSTREAM_ERROR);

  return WP2_STATUS_OK;
}

}    // namespace WP2
