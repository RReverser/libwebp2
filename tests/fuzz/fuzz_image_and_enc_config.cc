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
//  Fuzzing of raw pixels and still image wp2 encoding settings.
//
// Author: Yannis Guyon (yguyon@google.com)

#include "src/wp2/encode.h"
#include "tests/fuzz/fuzz_utils.h"
#include "tests/include/helpers.h"

extern "C" int LLVMFuzzerTestOneInput(const uint8_t* data, size_t size) {
  // Input is considered random data and converted to parameters.
  WP2::testing::FuzzedParameters params(data, size);

  // Extract a random config.
  // Encoding settings are limited to guarantee a reasonable test duration.
  const bool use_fast_config = params.ExtractBool();
  const bool use_advanced_config = params.ExtractBool();

  WP2::EncoderConfig config;
  params.ExtractSimpleConfig(/*max_speed=*/use_fast_config ? 2 : 5, &config);
  if (use_advanced_config) {
    params.ExtractAdvancedConfig(use_fast_config, &config);
  }

  // Create a randomly sized valid buffer with random pixels. It's best to
  // extract other params before because it will try to use all remaining bytes.
  WP2::ArgbBuffer original = params.ExtractArgbBuffer(
      /*max_num_pixels=*/use_fast_config ? 256u * 256u : 64u * 64u);
  if (original.IsEmpty()) return 0;  // No byte left to extract as raw pixels.

  const bool compare_after_decoding = !use_advanced_config;
  // Input can be random noise; output can get heavily distorted but it's fine.
  const float expected_distortion =
      WP2::testing::GetExpectedDistortion(config) - 10.f;
  return WP2::testing::TestImageEncConfig(
      original, config, compare_after_decoding, expected_distortion);
}
