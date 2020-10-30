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

// Progress hook test.

#include <array>
#include <functional>
#include <string>
#include <tuple>

#include "include/helpers.h"
#include "src/common/constants.h"
#include "src/wp2/decode.h"
#include "src/wp2/encode.h"

namespace WP2 {
namespace {

//------------------------------------------------------------------------------

class ProgressTester : public ProgressHook {
 public:
  bool OnUpdate(float progress) override {
    assert(progress >= 0.f && progress <= 1.f);
    assert(progress >= last_progress_);
    last_progress_ = progress;
    return (progress < fail_at_);
  }

  float last_progress_ = -1.f;
  float fail_at_ = 2.f;
};

class ProgressTest
    : public ::testing::TestWithParam<std::tuple<std::string, bool, uint32_t>> {
};

//------------------------------------------------------------------------------

TEST_P(ProgressTest, Simple) {
  const std::string src_file_name = std::get<0>(GetParam());
  const float quality = std::get<1>(GetParam()) ? 100.f : 75.f;
  const uint32_t thread_level = std::get<2>(GetParam());

  ArgbBuffer src;
  MemoryWriter data;
  ASSERT_WP2_OK(testing::CompressImage(src_file_name, &data, &src, quality));

  ProgressTester progress_tester;

  ArgbBuffer output;
  DecoderConfig config;
  config.thread_level = thread_level;
  config.progress_hook = &progress_tester;

  // Test that it starts at 0.
  for (size_t data_size = 0; data_size < kHeaderMinSize; ++data_size) {
    ASSERT_EQ(Decode(data.mem_, data_size, &output, config),
              WP2_STATUS_BITSTREAM_ERROR);
    ASSERT_EQ(progress_tester.last_progress_, 0.f);
    progress_tester.last_progress_ = -1.f;
  }

  // Test that it fails as expected.
  for (float fail_at :
       {0.f, kProgressDecHeader, 0.5f, 0.75f, 1.f - kProgressDecEnd, 0.9999f}) {
    progress_tester.fail_at_ = fail_at;
    ASSERT_EQ(Decode(data.mem_, data.size_, &output, config),
              WP2_STATUS_USER_ABORT);
    progress_tester.last_progress_ = -1.f;
  }
  progress_tester.fail_at_ = 2.f;

  // Test that it ends at 1.
  ASSERT_WP2_OK(Decode(data.mem_, data.size_, &output, config));
  ASSERT_EQ(progress_tester.last_progress_, 1.f);

  ASSERT_TRUE(testing::Compare(src, output, src_file_name,
                               testing::GetExpectedDistortion(quality)));
}

//------------------------------------------------------------------------------

INSTANTIATE_TEST_SUITE_P(
    ProgressTestInstantiationSingleTile, ProgressTest,
    ::testing::Combine(
        ::testing::Values("source1_1x1.png",
                          "source1_64x48.png") /* one-tile image */,
        ::testing::Values(true) /* lossless */,
        ::testing::Values(0) /* no multithread; just one tile anyway */));

INSTANTIATE_TEST_SUITE_P(
    ProgressTestInstantiationMultiTile, ProgressTest,
    ::testing::Combine(::testing::Values("source0.pgm") /* several tiles */,
                       ::testing::Values(false) /* lossy */,
                       ::testing::Values(4) /* extra threads */));

// This one takes a while to run so it is disabled.
// Can still be run with flag --test_arg=--gunit_also_run_disabled_tests
INSTANTIATE_TEST_SUITE_P(
    DISABLED_ProgressTestInstantiation, ProgressTest,
    ::testing::Combine(::testing::Values("source3.jpg", "alpha_ramp.png"),
                       ::testing::Values(true, false) /* lossless, lossy */,
                       ::testing::Values(0, 1, 64) /* extra threads */));

//------------------------------------------------------------------------------
// Test that WP2_STATUS_USER_ABORT is returned first, even with a failing
// bitstream.

TEST(ProgressTest, UserAbortBeforeAnyDecoding) {
  const size_t encoded_size = 64;
  std::string encoded_data(encoded_size, 'u');  // Junk

  DecoderConfig config = DecoderConfig::kDefault;
  ProgressTester progress_stopper;
  config.progress_hook = &progress_stopper;

  std::array<uint8_t, encoded_size> decoded_pixels{0};
  const uint32_t stride = 32;
  ArgbBuffer decoded_buffer;

  progress_stopper.fail_at_ = 0.f;
  ASSERT_EQ(DecodeArgb_32((const uint8_t*)encoded_data.data(), encoded_size,
                          decoded_pixels.data(), stride, decoded_pixels.size(),
                          config),
            WP2_STATUS_USER_ABORT);
  ASSERT_EQ(Decode(encoded_data, &decoded_buffer, config),
            WP2_STATUS_USER_ABORT);

  progress_stopper.fail_at_ = 0.001f;
  ASSERT_EQ(DecodeArgb_32((const uint8_t*)encoded_data.data(), encoded_size,
                          decoded_pixels.data(), stride, decoded_pixels.size(),
                          config),
            WP2_STATUS_BITSTREAM_ERROR);
  ASSERT_EQ(Decode(encoded_data, &decoded_buffer, config),
            WP2_STATUS_BITSTREAM_ERROR);
}

//------------------------------------------------------------------------------
// Test progress precision with a lot of pixels.

TEST(ProgressTest, Precision) {
  ArgbBuffer big_source;
  ASSERT_WP2_OK(big_source.Resize(2048, 2048));
  big_source.Fill({130, 129, 128, 127});

  MemoryWriter encoded_bitstream;
  EncoderConfig encoder_config = EncoderConfig::kDefault;
  encoder_config.speed = 0;
  ASSERT_WP2_OK(Encode(big_source, &encoded_bitstream, encoder_config));

  DecoderConfig decoder_config = DecoderConfig::kDefault;
  ProgressTester progress_tester;
  decoder_config.progress_hook = &progress_tester;

  ArgbBuffer decoded_buffer;
  ASSERT_WP2_OK(Decode(encoded_bitstream.mem_, encoded_bitstream.size_,
                       &decoded_buffer, decoder_config));
}

//------------------------------------------------------------------------------

}  // namespace
}  // namespace WP2
