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

// Test ImageIo animations.

#include <cstring>

#include "imageio/anim_image_dec.h"
#include "imageio/image_dec.h"
#include "include/helpers.h"
#include "src/wp2/encode.h"

namespace WP2 {
namespace {

//------------------------------------------------------------------------------

// TODO(yguyon): test broken / exotic GIF files
TEST(ReadImageTest, Animation) {
  const std::string file_gif = testing::GetTestDataPath("animation0.gif");
  const std::string file_webp = testing::GetTestDataPath("animation0.webp");
  const float expected_distortion = 35.f;

  std::vector<ArgbBuffer> frames_gif;
  std::vector<uint32_t> durations_gif;
  uint32_t loop_count_gif;

  std::vector<ArgbBuffer> frames_webp;
  std::vector<uint32_t> durations_webp;
  uint32_t loop_count_webp;

  ASSERT_WP2_OK(testing::ReadAnimation(file_gif, &frames_gif, &durations_gif,
                                       &loop_count_gif));
  ASSERT_WP2_OK(testing::ReadAnimation(file_webp, &frames_webp, &durations_webp,
                                       &loop_count_webp));

  ASSERT_EQ(frames_gif.size(), frames_webp.size());
  ASSERT_EQ(frames_gif.size(), durations_gif.size());
  ASSERT_EQ(frames_webp.size(), durations_webp.size());
  EXPECT_EQ(loop_count_gif, loop_count_webp);
  for (size_t i = 0; i < frames_gif.size(); ++i) {
    EXPECT_EQ(durations_gif[i], durations_webp[i]);
    EXPECT_TRUE(testing::Compare(frames_gif[i], frames_webp[i], "animation0",
                                 expected_distortion));
  }
}

//------------------------------------------------------------------------------

Data CreateData(const char* const str) {
  Data data;
  WP2_ASSERT_STATUS(
      data.CopyFrom(reinterpret_cast<const uint8_t*>(str), std::strlen(str)));
  return data;
}

TEST(ReadImageTest, AnimationMetadata) {
  EncoderConfig config;
  config.tile_shape = TILE_SHAPE_SQUARE_128;
  // Create input with metadata.
  ArgbBuffer input;
  ASSERT_WP2_OK(input.Resize(128 + 1, 128 * 2 + 1));
  input.Fill(Argb32b{0xFFu, 0xFFu, 0xFFu, 0xFFu});

  const Data good_exif = CreateData("good exif");
  const Data good_iccp = CreateData("good iccp");
  const Data good_xmp = CreateData("good xmp");
  ASSERT_WP2_OK(input.metadata.exif.CopyFrom(good_exif.bytes, good_exif.size));
  ASSERT_WP2_OK(input.metadata.iccp.CopyFrom(good_iccp.bytes, good_iccp.size));
  ASSERT_WP2_OK(input.metadata.xmp.CopyFrom(good_xmp.bytes, good_xmp.size));

  // Encode.
  AnimationEncoder encoder;
  for (int i = 0; i < 3; ++i) {
    ASSERT_WP2_OK(encoder.AddFrame(input, /*duration_ms=*/1));
  }
  MemoryWriter output;
  ASSERT_WP2_OK(encoder.Encode(&output, config, kInfiniteLoop, input.metadata));

  // Decode only first frame.
  {
    ArgbBuffer decoded;
    ASSERT_WP2_OK(ReadImage(output.mem_, output.size_, &decoded));
    EXPECT_TRUE(testing::HasSameData(decoded.metadata.exif, good_exif));
    EXPECT_TRUE(testing::HasSameData(decoded.metadata.iccp, good_iccp));
    EXPECT_TRUE(testing::HasSameData(decoded.metadata.xmp, good_xmp));
  }

  // Decode and check metadata for all frames.
  {
    ArgbBuffer decoded;
    ImageReader image_reader(output.mem_, output.size_, &decoded);
    bool is_last_frame = true;
    do {
      ASSERT_WP2_OK(image_reader.ReadFrame(&is_last_frame));
      EXPECT_TRUE(testing::HasSameData(decoded.metadata.exif, good_exif));
      EXPECT_TRUE(testing::HasSameData(decoded.metadata.iccp, good_iccp));
      EXPECT_TRUE(testing::HasSameData(decoded.metadata.xmp, good_xmp));
    } while (!is_last_frame);
  }
}

//------------------------------------------------------------------------------

}  // namespace
}  // namespace WP2
