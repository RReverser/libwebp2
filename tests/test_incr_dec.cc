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

// Incremental decoding test.

#include <string>
#include <tuple>

#include "include/helpers.h"
#include "include/helpers_incr.h"
#include "src/wp2/decode.h"

namespace WP2 {
namespace {

//------------------------------------------------------------------------------

class IncrementalDecoderTest
    : public ::testing::TestWithParam<
          std::tuple<std::string, bool, uint32_t,
                     DecoderConfig::IncrementalMode, testing::DecoderType>> {};

TEST_P(IncrementalDecoderTest, Simple) {
  const std::string src_file_name = std::get<0>(GetParam());
  const bool lossless = std::get<1>(GetParam());
  const float quality = lossless ? 100.f : 75.f;
  const uint32_t thread_level = std::get<2>(GetParam());
  const DecoderConfig::IncrementalMode incremental_mode =
      std::get<3>(GetParam());
  const testing::DecoderType decoder_type = std::get<4>(GetParam());

  ArgbBuffer src;
  MemoryWriter writer;
  ASSERT_WP2_OK(testing::CompressImage(src_file_name, &writer, &src, quality));

  DecoderConfig config = DecoderConfig::kDefault;
  config.thread_level = thread_level;
  config.incremental_mode = incremental_mode;

  ArgbBuffer reference_output(WP2_Argb_32);
  if (!lossless) {
    ASSERT_WP2_OK(Decode(writer.mem_, writer.size_, &reference_output, config));
  }

  // Try incremental decoding with different steps.
  for (const size_t step : {(size_t)1, writer.size_ / 20, writer.size_}) {
    std::vector<ArgbBuffer> decoded_frames;
    ASSERT_WP2_OK(
        testing::DecodeIncremental(config, {writer.mem_, writer.size_},
                                   {step, decoder_type, {}}, &decoded_frames));
    ASSERT_EQ(decoded_frames.size(), 1u);
    EXPECT_TRUE(testing::Compare(src, decoded_frames.front(), src_file_name,
                                 testing::GetExpectedDistortion(quality)));

    // Assert that the output of the incremental decoder is the same as Decode()
    // for lossy.
    if (!lossless) {
      EXPECT_TRUE(testing::Compare(decoded_frames.front(), reference_output,
                                   src_file_name));
    }
  }
}

//------------------------------------------------------------------------------

INSTANTIATE_TEST_SUITE_P(
    IncrementalDecoderTestInstantiationPartialMultiTile, IncrementalDecoderTest,
    ::testing::Combine(
        ::testing::Values("source0.pgm") /* several tiles */,
        ::testing::Values(false) /* lossy */,
        ::testing::Values(64) /* extra threads */,
        ::testing::Values(DecoderConfig::IncrementalMode::PARTIAL_TILE_CONTEXT),
        ::testing::Values(testing::DecoderType::kUnstableCustom)));

//------------------------------------------------------------------------------

// Tests that fully decoding then rewinding an image generates the same metadata
// position in the bitstream.
TEST(IncrementalDecoderTest, RewindAfterMetadata) {
  MemoryWriter writer;
  ASSERT_WP2_OK(testing::CompressImage("source1_64x48.png", &writer));

  testing::IncrementalDecodingTestSetup setup(
      writer.size_, testing::DecoderType::kArray,
      {{testing::DecoderAction::Type::kRewindKeepBytes, writer.size_ + 1,
        /*output_buffer_changed*/ 0}});
  ASSERT_WP2_OK(testing::DecodeIncremental(DecoderConfig::kDefault,
                                           {writer.mem_, writer.size_}, setup));
}

//------------------------------------------------------------------------------

}  // namespace
}  // namespace WP2
