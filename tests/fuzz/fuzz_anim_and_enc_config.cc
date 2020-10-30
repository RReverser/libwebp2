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
//  Fuzzing of raw pixels and animation wp2 encoding settings.
//
// Author: Yannis Guyon (yguyon@google.com)

#include <algorithm>
#include <vector>

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
  params.ExtractSimpleConfig(/*max_speed=*/use_fast_config ? 0 : 3, &config);
  if (use_advanced_config) {
    params.ExtractAdvancedConfig(use_fast_config, &config);
  }

  // Extract a random (maximum) number of frames and random durations.
  uint32_t num_frames_x = params.ExtractUInt32(1, use_fast_config ? 4 : 1);
  uint32_t num_frames_y = params.ExtractUInt32(1, use_fast_config ? 4 : 1);
  uint32_t num_frames = num_frames_x * num_frames_y;  // Up to 16
  std::vector<uint32_t> durations_ms(num_frames);
  for (uint32_t& duration_ms : durations_ms) {
    duration_ms = params.ExtractUInt32(1, WP2::kMaxFrameDurationMs);
  }

  const uint32_t loop_count =
      params.ExtractBool() ? WP2::kInfiniteLoop
                           : params.ExtractUInt32(1, WP2::kMaxLoopCount);

  // Create a randomly sized valid buffer with random pixels and use it as an
  // atlas that contains all frames. It's best to extract other params before
  // because it will try to use all remaining bytes.
  WP2::ArgbBuffer atlas = params.ExtractArgbBuffer(
      /*max_num_pixels=*/use_fast_config ? (128u * 128u * num_frames)
                                         : (64u * 64u * num_frames));
  if (atlas.IsEmpty()) return 0;  // No byte left to extract as raw pixels.

  num_frames_x = std::min(num_frames_x, atlas.width);
  num_frames_y = std::min(num_frames_y, atlas.height);
  num_frames = num_frames_x * num_frames_y;
  const uint32_t frame_width = atlas.width / num_frames_x;
  const uint32_t frame_height = atlas.height / num_frames_y;

  std::vector<WP2::ArgbBuffer> original_frames;
  original_frames.reserve(num_frames);
  for (uint32_t frame_x = 0; frame_x < num_frames_x; ++frame_x) {
    for (uint32_t frame_y = 0; frame_y < num_frames_y; ++frame_y) {
      original_frames.emplace_back(atlas.format);
      if (original_frames.back().SetView(
              atlas, {frame_x * frame_width, frame_y * frame_height,
                      frame_width, frame_height}) != WP2_STATUS_OK) {
        abort();
      }
    }
  }
  durations_ms.resize(num_frames);

  const bool compare_after_decoding = !use_advanced_config;
  // Input can be random noise; output can get heavily distorted but it's fine.
  const float expected_distortion =
      WP2::testing::GetExpectedDistortion(config) - 10.f;
  return WP2::testing::TestAnimEncConfig(
      original_frames, durations_ms, loop_count, config, compare_after_decoding,
      expected_distortion);
}
