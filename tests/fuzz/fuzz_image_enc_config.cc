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

#include <iostream>
#include <string>

#include "imageio/image_dec.h"
#include "src/wp2/encode.h"
#include "tests/fuzz/fuzz_utils.h"
#include "tests/include/helpers.h"

extern "C" int LLVMFuzzerTestOneInput(const uint8_t* data, size_t size) {
  // Input is considered random data and converted to parameters.
  WP2::testing::FuzzedParameters params(data, size);
  const char* file_name;
  const uint8_t* file_data;
  size_t file_size;
  if (!params.ExtractSource(&file_name, &file_data, &file_size)) abort();

  // Read a randomly chosen valid image file.
  WP2::ArgbBuffer original;
  const WP2Status read_status =
      WP2::ReadImage(file_data, file_size, &original, WP2::FileFormat::AUTO,
                     WP2::LogLevel::QUIET);
  if (read_status != WP2_STATUS_OK) {
    std::cerr << "Unable to read " << file_name << ": "
              << WP2GetStatusMessage(read_status) << std::endl;
    abort();
  }

  // Extract a random config.
  // Encoding settings are limited to guarantee a reasonable test duration.
  const bool use_fast_config = (original.width > 64 || original.height > 64);
  const bool use_advanced_config = params.ExtractBool();

  WP2::EncoderConfig config;
  params.ExtractSimpleConfig(/*max_speed=*/use_fast_config ? 1 : 4, &config);
  if (use_advanced_config) {
    params.ExtractAdvancedConfig(use_fast_config, &config);
  }

  const bool compare_after_decoding = !use_advanced_config;
  const float expected_distortion = WP2::testing::GetExpectedDistortion(config);
  return WP2::testing::TestImageEncConfig(
      original, config, compare_after_decoding, expected_distortion);
}
