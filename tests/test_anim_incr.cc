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

// Incremental animation decoding test.

#include <string>
#include <tuple>
#include <vector>

#include "imageio/image_dec.h"
#include "include/helpers.h"
#include "include/helpers_incr.h"
#include "src/wp2/decode.h"
#include "src/wp2/encode.h"
#include "src/wp2/format_constants.h"

namespace WP2 {
namespace {

//------------------------------------------------------------------------------

typedef std::tuple<std::vector<const char*>, std::vector<uint32_t>, float, int,
                   uint32_t, testing::IncrementalDecodingTestSetup>
    Param;

class AnimTestIncr : public ::testing::TestWithParam<Param> {};

TEST_P(AnimTestIncr, Simple) {
  const std::vector<const char*>& file_names = std::get<0>(GetParam());
  const std::vector<uint32_t>& durations_ms = std::get<1>(GetParam());
  const float quality = std::get<2>(GetParam());
  const int speed = std::get<3>(GetParam());
  const uint32_t thread_level = std::get<4>(GetParam());
  testing::IncrementalDecodingTestSetup setup = std::get<5>(GetParam());

  std::vector<ArgbBuffer> frames;
  MemoryWriter encoded_data;
  ASSERT_WP2_OK(testing::CompressAnimation(file_names, durations_ms,
                                           &encoded_data, &frames, quality,
                                           speed, thread_level));

  DecoderConfig config = DecoderConfig::kDefault;
  config.thread_level = thread_level;

  // Try incremental decoding with different steps.
  const DataView input = {encoded_data.mem_, encoded_data.size_};
  for (const size_t step : {(size_t)1, input.size / 20, input.size}) {
    setup.incr_size_step = step;
    std::vector<ArgbBuffer> decoded_frames;
    std::vector<uint32_t> decoded_durations_ms;
    ASSERT_WP2_OK(testing::DecodeIncremental(
        config, input, setup, &decoded_frames, &decoded_durations_ms));

    // Hard to know which frames to compare if some were skipped.
    if (setup.actions.empty() ||
        setup.actions.back().type != testing::DecoderAction::Type::kSkip) {
      ASSERT_EQ(decoded_frames.size(), frames.size());
      for (size_t i = 0; i < frames.size(); ++i) {
        EXPECT_TRUE(testing::Compare(decoded_frames[i], frames[i],
                                     file_names[i],
                                     testing::GetExpectedDistortion(quality)));
        EXPECT_EQ(decoded_durations_ms[i], durations_ms[i]);
      }
    }
  }
}

// Input images are expected to be of equal sizes inside a Param.
INSTANTIATE_TEST_SUITE_P(
    AnimTestIncrInstantiation, AnimTestIncr,
    ::testing::Values(
        Param({"source1_64x48.png"},
              /*frame durations (ms):*/ {1000},
              /*quality=*/0.0f, /*speed=*/4, /*thread_level=*/2,
              {1, testing::DecoderType::kArray, {}}),
        Param({"alpha_ramp.png", "alpha_ramp.webp"},
              /*frame durations (ms):*/ {50, 1},
              /*quality=*/100.0f, /*speed=*/0, /*thread_level=*/4,
              {1, testing::DecoderType::kUnstableCustom, {}}),
        Param({"source1_64x48.png"},
              /*frame durations (ms):*/ {1000},
              /*quality=*/40.0f, /*speed=*/0, /*thread_level=*/0,
              {1, testing::DecoderType::kSwappingCustom, {}})));

// This one takes a while to run so it is disabled.
// Can still be run with flag --test_arg=--gunit_also_run_disabled_tests
INSTANTIATE_TEST_SUITE_P(
    DISABLED_AnimTestIncrInstantiation, AnimTestIncr,
    ::testing::Values(
        Param({"alpha_ramp.png", "alpha_ramp.lossy.webp", "alpha_ramp.webp"},
              /*frame durations (ms):*/ {50, 70, 1},
              /*quality=*/100.0f, /*speed=*/2, /*thread_level=*/0,
              {1, testing::DecoderType::kUnstableCustom, {}}),
        Param({"source0.pgm", "source1.png", "source2.tiff", "source3.jpg"},
              /*frame durations (ms):*/ {kMaxFrameDurationMs, 1, 20, 3},
              /*quality=*/75.0f, /*speed=*/4, /*thread_level=*/8,
              {1, testing::DecoderType::kStream, {}}),
        Param({"source3.jpg", "source2.tiff", "source3.jpg"},
              /*frame durations (ms):*/ {123, 456, 789},
              /*quality=*/25.0f, /*speed=*/3, /*thread_level=*/5,
              {1, testing::DecoderType::kSwappingCustom, {}}),
        Param({"source0.pgm", "source3.jpg"},
              /*frame durations (ms):*/ {kMaxFrameDurationMs, 3},
              /*quality=*/75.0f, /*speed=*/2, /*thread_level=*/8,
              {1, testing::DecoderType::kStream, {}})));

INSTANTIATE_TEST_SUITE_P(
    AnimTestIncrRewind, AnimTestIncr,
    ::testing::Values(
        Param({"source1_64x48.png"},
              /*frame durations (ms):*/ {1000},
              /*quality=*/0.0f, /*speed=*/4, /*thread_level=*/2,
              {1,
               testing::DecoderType::kArray,
               {{testing::DecoderAction::Type::kRewind, 0, 0},
                {testing::DecoderAction::Type::kRewind, 100, 0}}})));

// This one takes a while to run so it is disabled.
// Can still be run with flag --test_arg=--gunit_also_run_disabled_tests
INSTANTIATE_TEST_SUITE_P(
    DISABLED_AnimTestIncrRewind, AnimTestIncr,
    ::testing::Values(
        Param({"source1_64x48.png"},
              /*frame durations (ms):*/ {1},
              /*quality=*/100.0f, /*speed=*/7, /*thread_level=*/0,
              {1,
               testing::DecoderType::kArray,
               {{testing::DecoderAction::Type::kRewind, 100000000, 0}}}),
        Param({"test_exif_xmp.webp"},
              /*frame durations (ms):*/ {1},
              /*quality=*/33.0f, /*speed=*/0, /*thread_level=*/1,
              {1,
               testing::DecoderType::kArray,
               {{testing::DecoderAction::Type::kSkip, 35, 1}}}),
        Param({"alpha_ramp.png", "alpha_ramp.lossy.webp", "alpha_ramp.webp"},
              /*frame durations (ms):*/ {50, 70, 1},
              /*quality=*/100.0f, /*speed=*/2, /*thread_level=*/0,
              {1,
               testing::DecoderType::kUnstableCustom,
               {{testing::DecoderAction::Type::kSkip, 355, 1},
                {testing::DecoderAction::Type::kRewindKeepBytes, 356, 0},
                {testing::DecoderAction::Type::kSkip, 357, 0}}}),
        Param({"source0.pgm", "source1.png", "source2.tiff", "source3.jpg"},
              /*frame durations (ms):*/ {kMaxFrameDurationMs, 1, 20, 3},
              /*quality=*/75.0f, /*speed=*/4, /*thread_level=*/8,
              {1,
               testing::DecoderType::kStream,
               {{testing::DecoderAction::Type::kSkip, 15, kMaxNumFrames + 1},
                {testing::DecoderAction::Type::kRewindKeepBytes, 12000, 1}}}),
        Param({"source0.pgm", "source1.png", "source2.tiff", "source3.jpg"},
              /*frame durations (ms):*/ {5689, 21121, 52, 3},
              /*quality=*/50.0f, /*speed=*/5, /*thread_level=*/4,
              {1,
               testing::DecoderType::kArray,
               {{testing::DecoderAction::Type::kSkip, 0, kMaxNumFrames},
                {testing::DecoderAction::Type::kRewind, 100000000, 0},
                {testing::DecoderAction::Type::kSkip, 0, 1},
                {testing::DecoderAction::Type::kRewind,
                 kHeaderMaxSize + kANMFHeaderSize + 1, 0},
                {testing::DecoderAction::Type::kSkip, 0, kMaxNumFrames},
                {testing::DecoderAction::Type::kRewind, 100000000, 0}}}),
        Param({"source3.jpg", "source2.tiff", "source3.jpg"},
              /*frame durations (ms):*/ {123, 456, 789},
              /*quality=*/25.0f, /*speed=*/3, /*thread_level=*/5,
              {1,
               testing::DecoderType::kSwappingCustom,
               {{testing::DecoderAction::Type::kSkip, 0, 2},
                {testing::DecoderAction::Type::kRewind, 100000000, 1}}})));

//------------------------------------------------------------------------------

}  // namespace
}  // namespace WP2
