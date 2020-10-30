// Copyright 2019 Google LLC
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//      http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.
// -----------------------------------------------------------------------------
//
//  Fuzzing of animation wp2 encoding settings.
//
// Author: Yannis Guyon (yguyon@google.com)

#include "imageio/anim_image_dec.h"
#include "src/utils/orientation.h"
#include "src/wp2/decode.h"
#include "src/wp2/encode.h"
#include "tests/include/helpers.h"

namespace WP2 {
namespace testing {

int TestAnimEncConfig(const std::vector<ArgbBuffer>& original_frames,
                      const std::vector<uint32_t>& durations_ms,
                      uint32_t loop_count, const EncoderConfig& config,
                      bool compare_after_decoding, float expected_distortion) {
  // std::cerr << config << std::endl;  // Uncomment to print config
  if (!config.IsValid()) abort();

  // Add frames to the encoder.
  AnimationEncoder animation_encoder;
  for (size_t i = 0; i < original_frames.size(); ++i) {
    if (animation_encoder.AddFrame(original_frames[i], durations_ms[i]) !=
        WP2_STATUS_OK) {
      abort();
    }
  }

  // Encode with a reasonable config, without diving into detailed settings.
  MemoryWriter memory_writer;
  const WP2Status encoder_status =
      animation_encoder.Encode(&memory_writer, config, (uint8_t)loop_count);
  if (encoder_status != WP2_STATUS_OK) abort();

  // Verify that it decodes fine and that the result is close enough.
  ArgbBuffer decoded;
  ArgbBuffer decoded_unoriented;
  ArrayDecoder decoder(memory_writer.mem_, memory_writer.size_,
                       DecoderConfig::kDefault, &decoded);
  for (size_t i = 0; i < original_frames.size(); ++i) {
    uint32_t duration_ms;
    if (!decoder.ReadFrame(&duration_ms)) abort();
    assert(!decoder.Failed());
    if (duration_ms != durations_ms[i]) abort();

    if (compare_after_decoding) {
      // 'decoded' can't be modified as long as 'decoder' is used, so
      // 'decoded_unoriented' is compared instead.
      if (config.decoding_orientation == Orientation::kOriginal) {
        if (decoded_unoriented.SetView(decoded) != WP2_STATUS_OK) {
          abort();
        }
      } else if (RotateBuffer(
                     GetInverseOrientation(config.decoding_orientation),
                     decoded, &decoded_unoriented) != WP2_STATUS_OK) {
        abort();
      }

      // Tiny frames are not interesting to compare.
      if (original_frames[i].width * original_frames[i].height >= 64 &&
          !Compare(original_frames[i], decoded_unoriented, "anim",
                   expected_distortion)) {
        // Uncomment to save the original and decoded frames as PNG files.
        // WP2::testing::SaveBeforeAfter(
        //     original_frames[i], decoded_unoriented,
        //     "/tmp/comparison_" + std::to_string(i));
        abort();
      }
    }
  }
  if (decoder.GetStatus() != WP2_STATUS_OK) abort();
  return 0;
}

}  // namespace testing
}  // namespace WP2
