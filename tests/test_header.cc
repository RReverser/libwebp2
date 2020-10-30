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

// Test different headers (~encoder configs)

#include "include/helpers.h"
#include "src/wp2/encode.h"
#include "src/wp2/decode.h"
#include "src/wp2/format_constants.h"

namespace WP2 {
namespace {

//------------------------------------------------------------------------------

TEST(HeaderTest, MinHeaderSize) {
  ArgbBuffer image;
  ASSERT_WP2_OK(image.Resize(1, 1));
  image.Fill(Argb32b{255, 0, 0, 0});

  MemoryWriter memory_writer;
  ASSERT_WP2_OK(Encode(image, &memory_writer));

  BitstreamFeatures features;
  ASSERT_WP2_OK(features.Read(memory_writer.mem_, memory_writer.size_));
  EXPECT_EQ(features.header_size, kHeaderMinSize);
}

TEST(HeaderTest, MaxHeaderSize) {
  const uint8_t num_frames = 2;
  const uint32_t duration_ms = 100u;
  AnimationEncoder encoder;
  for (uint8_t frame_index = 0; frame_index < num_frames; ++frame_index) {
    ArgbBuffer frame;
    ASSERT_WP2_OK(frame.Resize(1, 1));
    const Argb32b color{255, frame_index, frame_index, frame_index};
    frame.Fill(color);
    ASSERT_WP2_OK(encoder.AddFrame(frame, duration_ms));
  }

  EncoderConfig encoder_config = EncoderConfig::kDefault;
  encoder_config.transfer_function = WP2_TF_ITU_R_BT709;  // Any besides BT2020
  MemoryWriter memory_writer;
  ASSERT_WP2_OK(encoder.Encode(&memory_writer, encoder_config));

  BitstreamFeatures features;
  ASSERT_WP2_OK(features.Read(memory_writer.mem_, memory_writer.size_));
  // TODO(yguyon): Compare with exactly 'kHeaderMaxSize' when it's possible to
  //               reach it (with a 'kBackgroundCustom' color)
  EXPECT_EQ(features.header_size, kHeaderMaxSize - 5);
}

//------------------------------------------------------------------------------

TEST(HeaderTest, TestAlpha) {
  for (bool is_opaque : {true, false}) {
    ArgbBuffer image;
    ASSERT_WP2_OK(image.Resize(1, 1));
    image.Fill(Argb32b{(uint8_t)(is_opaque ? 255 : 128), 0, 0, 0});

    MemoryWriter memory_writer;
    ASSERT_WP2_OK(Encode(image, &memory_writer));

    BitstreamFeatures features;
    ASSERT_WP2_OK(features.Read(memory_writer.mem_, memory_writer.size_));
    EXPECT_EQ(features.is_opaque, is_opaque);
  }
}

//------------------------------------------------------------------------------

}  // namespace
}  // namespace WP2
