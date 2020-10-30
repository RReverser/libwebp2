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
//  Fuzzing of still image wp2 encoding settings.
//
// Author: Yannis Guyon (yguyon@google.com)

#include "src/utils/orientation.h"
#include "src/wp2/decode.h"
#include "src/wp2/encode.h"
#include "tests/fuzz/fuzz_utils.h"
#include "tests/include/helpers.h"

namespace WP2 {
namespace testing {

int TestImageEncConfig(const ArgbBuffer& original, const EncoderConfig& config,
                       bool compare_after_decoding, float expected_distortion) {
  // std::cerr << config << std::endl;  // Uncomment to print config
  if (!config.IsValid()) abort();

  MemoryWriter memory_writer;
  const WP2Status status = Encode(original, &memory_writer, config);
  if (status != WP2_STATUS_OK) abort();

  // Verify that it decodes fine and that the result is close enough.
  ArgbBuffer decoded;
  if (Decode(memory_writer.mem_, memory_writer.size_, &decoded) !=
      WP2_STATUS_OK) {
    abort();
  }
  if (compare_after_decoding) {
    if (RotateBuffer(GetInverseOrientation(config.decoding_orientation),
                     &decoded) != WP2_STATUS_OK) {
      abort();
    }
    // Tiny images are not interesting to compare.
    if (original.width * original.height >= 64 &&
        !Compare(original, decoded, "image", expected_distortion)) {
      // Uncomment to save the original and decoded images as PNG files.
      // WP2::testing::SaveBeforeAfter(original, decoded, "/tmp/comparison");
      abort();
    }
  }
  return 0;
}

int TestImageDecConfig(DataView encoded_data, const DecoderConfig& config,
                       bool expected_success) {
  ArgbBuffer decoded;
  const WP2Status status =
      Decode(encoded_data.bytes, encoded_data.size, &decoded);
  if (status != WP2_STATUS_OK && expected_success) abort();
  return 0;
}

}  // namespace testing
}  // namespace WP2
