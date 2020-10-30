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

// Incremental decoding test function.

#ifndef WP2_TESTS_INCLUDE_HELPERS_INCR_H_
#define WP2_TESTS_INCLUDE_HELPERS_INCR_H_

#include <vector>

#include "src/utils/utils.h"
#include "src/wp2/base.h"
#include "src/wp2/decode.h"

namespace WP2 {
namespace testing {

enum class DecoderType { kArray, kStream, kUnstableCustom, kSwappingCustom };
struct DecoderAction {
  enum class Type { kRewind, kRewindKeepBytes, kSkip } type;
  size_t bistream_position;
  uint32_t value;  // Whether to keep the output buffer if 'kRewind',
                   // or number of frames to skip.

  DecoderAction(Type t, size_t p, uint32_t v)
      : type(t), bistream_position(p), value(v) {}
};

struct IncrementalDecodingTestSetup {
  size_t incr_size_step = 1;
  DecoderType decoder_type = DecoderType::kArray;
  std::vector<DecoderAction> actions;

  IncrementalDecodingTestSetup(size_t s, DecoderType t,
                               std::vector<DecoderAction>&& a)
      : incr_size_step(s), decoder_type(t), actions(a) {}
};

WP2Status DecodeIncremental(
    const DecoderConfig& config, DataView input,
    const IncrementalDecodingTestSetup& setup,
    std::vector<ArgbBuffer>* const decoded_frames = nullptr,
    std::vector<uint32_t>* const decoded_durations_ms = nullptr);

}  // namespace testing
}  // namespace WP2

#endif  // WP2_TESTS_INCLUDE_HELPERS_INCR_H_
