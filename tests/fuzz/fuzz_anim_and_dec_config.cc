// Copyright 2020 Google LLC
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
//  Fuzzing of a wp2 animated bitstream and decoding settings.
//
// Author: Yannis Guyon (yguyon@google.com)

#include <vector>

#include "src/wp2/decode.h"
#include "src/wp2/format_constants.h"
#include "tests/fuzz/fuzz_utils.h"
#include "tests/include/helpers_incr.h"

//------------------------------------------------------------------------------

namespace WP2 {
namespace testing {
namespace {

class FuzzedParametersSetup : public FuzzedParameters {
 public:
  using FuzzedParameters::FuzzedParameters;

  IncrementalDecodingTestSetup ExtractSetup(size_t max_bitstream_size) {
    const size_t incr_size_step =
        ExtractUInt32(1, std::max((size_t)kHeaderMinSize, max_bitstream_size));
    const DecoderType decoder_type =
        (DecoderType)ExtractUInt32((uint32_t)DecoderType::kSwappingCustom);

    std::vector<DecoderAction> actions;
    actions.reserve(ExtractUInt32(0, 15));
    for (uint32_t i = 0; i < actions.capacity(); ++i) {
      const DecoderAction::Type type = (DecoderAction::Type)ExtractUInt32(
          (uint32_t)DecoderAction::Type::kSkip);
      const size_t bistream_position = ExtractUInt32(0, max_bitstream_size);
      switch (type) {
        case DecoderAction::Type::kRewind:
        case DecoderAction::Type::kRewindKeepBytes: {
          const bool keep_output_buffer = ExtractBool();
          actions.emplace_back(type, bistream_position,
                               (uint32_t)keep_output_buffer);
          break;
        }
        case DecoderAction::Type::kSkip: {
          const uint32_t num_frames_to_skip =
              ExtractBool() ? ExtractUInt32(0, 7) : kMaxNumFrames;
          actions.emplace_back(type, bistream_position, num_frames_to_skip);
          break;
        }
      }
    }
    return IncrementalDecodingTestSetup(incr_size_step, decoder_type,
                                        std::move(actions));
  }
};

}  // namespace
}  // namespace testing
}  // namespace WP2

//------------------------------------------------------------------------------

extern "C" int LLVMFuzzerTestOneInput(const uint8_t* data, size_t size) {
  // Input is considered random data and converted to parameters.
  WP2::testing::FuzzedParametersSetup params(data, size);

  // Extract a random decoding pattern (when to rewind, which frames to skip).
  const WP2::testing::IncrementalDecodingTestSetup setup =
      params.ExtractSetup(size);

  // Extract a random config.
  WP2::DecoderConfig config;
  params.ExtractConfig(&config);

  // Create a (possibly invalid) wp2 bitstream, starting with the wp2 signature.
  // It's best to extract other params before because it will try to use all
  // remaining bytes.
  WP2::Data encoded_data;
  if (!params.TryExtractWP2Bitstream(&encoded_data)) return 0;  // malloc fail.

  // Success is not expected but it shouldn't crash or return a weird error.
  const WP2Status status = WP2::testing::DecodeIncremental(
      config, {encoded_data.bytes, encoded_data.size}, setup);
  if (status != WP2_STATUS_OK &&
      status != WP2_STATUS_NOT_ENOUGH_DATA &&
      status != WP2_STATUS_BITSTREAM_ERROR &&
      status != WP2_STATUS_OUT_OF_MEMORY) {
    abort();
  }
  return 0;
}
