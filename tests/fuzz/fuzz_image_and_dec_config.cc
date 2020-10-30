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
//  Fuzzing of a wp2 bitstream and still image decoding settings.
//
// Author: Yannis Guyon (yguyon@google.com)

#include "src/wp2/decode.h"
#include "tests/fuzz/fuzz_utils.h"

extern "C" int LLVMFuzzerTestOneInput(const uint8_t* data, size_t size) {
  // Input is considered random data and converted to parameters.
  WP2::testing::FuzzedParameters params(data, size);

  // Extract a random config.
  WP2::DecoderConfig config;
  params.ExtractConfig(&config);

  // Create a (possibly invalid) wp2 bitstream, starting with the wp2 signature.
  // It's best to extract other params before because it will try to use all
  // remaining bytes.
  WP2::Data encoded_data;
  if (!params.TryExtractWP2Bitstream(&encoded_data)) return 0;  // malloc fail.

  // Success is not expected but it shouldn't crash.
  return WP2::testing::TestImageDecConfig(
      {encoded_data.bytes, encoded_data.size}, config,
      /*expected_success=*/false);
}
