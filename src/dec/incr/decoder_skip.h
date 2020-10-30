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
// Functions used to skip frames and glimpse at future frame features.
//
// Author: Yannis Guyon (yguyon@google.com)

#ifndef WP2_DEC_INCR_DECODER_SKIP_H_
#define WP2_DEC_INCR_DECODER_SKIP_H_

#include "src/dec/incr/decoder_state.h"
#include "src/wp2/decode.h"

namespace WP2 {

//------------------------------------------------------------------------------

// Sets the glimpsed frame to the same bitstream position as the current frame.
void CopyGlimpseFromCurrent(Decoder::State* const state);

// Skips and reads the bitstream till the next frame features are known.
WP2Status GlimpseTillNextANMFChunk(const DecoderConfig& config,
                                   const BitstreamFeatures& features,
                                   Decoder::State* const state);

// Advances the glimpsed frame without modifying the decoding state.
WP2Status GlimpseAtNextFrame(const DecoderConfig& config,
                             const BitstreamFeatures& features,
                             Decoder::State* const state);

// Returns true if a usable frame is already fully present in the output buffer.
// In this case, execution should continue as if the last tile was just decoded.
bool FrameIsAlreadyInOutputBuffer(uint32_t tile_width, uint32_t tile_height,
                                  Decoder::State* const state);

// If possible, skips entire frames and discards the related bitstream part.
bool DiscardWholeFrames(Decoder::State* const state);

// If needed and possible, all remaining frames (known or unknown) are discarded
// and the bistream is advanced till the ending metadata chunks.
bool DiscardTillMetadata(Decoder::State* const state);

//------------------------------------------------------------------------------

}  // namespace WP2

#endif /* WP2_DEC_INCR_DECODER_SKIP_H_ */
