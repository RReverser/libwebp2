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

#include <iostream>
#include <string>
#include <vector>

#include "imageio/anim_image_dec.h"
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

  // Read a randomly chosen valid image file (possibly an animation).
  std::vector<WP2::ArgbBuffer> original_frames;
  std::vector<uint32_t> durations_ms;
  uint32_t loop_count = 1;
  const WP2Status read_status = WP2::testing::ReadAnimation(
      file_data, file_size, &original_frames, &durations_ms, &loop_count);
  if (read_status != WP2_STATUS_OK) {
    std::cerr << "Unable to read " << file_name << ": "
              << WP2GetStatusMessage(read_status) << std::endl;
    abort();
  }
  if (loop_count == WP2::ImageReader::kInfiniteLoopCount ||
      loop_count > WP2::kMaxLoopCount) {
    loop_count = WP2::kInfiniteLoop;
  }

  // Encode up to ~1M pixels to avoid timeouts.
  const size_t kMaxNumPixels = (1u << 20);
  const size_t num_pixels_per_frame =
      original_frames.front().width * original_frames.front().height;
  const size_t max_num_frames = 1u + kMaxNumPixels / num_pixels_per_frame;
  if (original_frames.size() > max_num_frames) {
    original_frames.resize(max_num_frames);
    durations_ms.resize(max_num_frames);
  }

  // Make sure durations are valid for wp2 (to keep the same number of frames).
  for (size_t i = 0; i < original_frames.size(); ++i) {
    if (durations_ms[i] == WP2::ImageReader::kInfiniteDuration ||
        durations_ms[i] > WP2::kMaxFrameDurationMs) {
      durations_ms[i] = params.ExtractUInt32(1, WP2::kMaxFrameDurationMs);
    }
  }

  // Extract a random config.
  // Encoding settings are limited to guarantee a reasonable test duration.
  const bool use_fast_config = (original_frames.front().width > 128u ||
                                original_frames.front().height > 128u);
  const bool use_advanced_config = params.ExtractBool();

  WP2::EncoderConfig config;
  params.ExtractSimpleConfig(/*max_speed=*/use_fast_config ? 0 : 3, &config);
  if (use_advanced_config) {
    params.ExtractAdvancedConfig(use_fast_config, &config);
  }

  const bool compare_after_decoding = !use_advanced_config;
  const float expected_distortion = WP2::testing::GetExpectedDistortion(config);
  return WP2::testing::TestAnimEncConfig(
      original_frames, durations_ms, loop_count, config, compare_after_decoding,
      expected_distortion);
}
