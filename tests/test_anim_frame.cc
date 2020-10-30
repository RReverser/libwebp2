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

// Tests for FrameEncoder.

#include <cstring>
#include <string>
#include <tuple>
#include <vector>

#include "examples/example_utils.h"
#include "imageio/image_dec.h"
#include "imageio/image_enc.h"
#include "include/helpers.h"
#include "src/common/color_precision.h"
#include "src/enc/anim/anim_enc.h"
#include "src/enc/wp2_enc_i.h"
#include "src/utils/orientation.h"
#include "src/utils/utils.h"
#include "src/wp2/decode.h"
#include "src/wp2/encode.h"

namespace WP2 {
namespace {

using WP2::testing::GetTestDataPath;

TEST(AnimFrameTest, ForceLossless) {
  std::vector<ArgbBuffer> buffers;
  ASSERT_WP2_OK(
      testing::ReadImages({GetTestDataPath("logo.png"),
                           GetTestDataPath("background.png")},
                           &buffers));

  Vector<Frame> frames;
  ASSERT_TRUE(frames.resize(2));
  Frame& first_frame = frames[0];
  first_frame.duration_ms = 10;
  first_frame.force_lossless = true;
  ASSERT_WP2_OK(first_frame.rgb_pixels.SetView(buffers[0]));

  Frame& second_frame = frames[1];
  second_frame.duration_ms = 20;
  ASSERT_WP2_OK(second_frame.rgb_pixels.SetView(buffers[1]));

  const bool has_alpha = !AllFramesAreOpaque(frames);

  EncoderConfig config;
  config.quality = 80;
  config.tile_shape = TILE_SHAPE_SQUARE_256;
  MemoryWriter encoded_data;

  const Argb38b background_color = kTransparentArgb38b;
  ASSERT_WP2_OK(
      EncodeHeader(config, first_frame.GetWidth(), first_frame.GetHeight(),
                   has_alpha, /*is_anim=*/(frames.size() > 1), /*loop_count=*/1,
                   background_color, RGB12b{10, 0, 30}, /*has_icc=*/false,
                   /*has_trailing_data=*/false, &encoded_data));

  FrameEncoder frame_encode(config, background_color, frames);
  ASSERT_WP2_OK(frame_encode.Encode(&encoded_data));

  // Decode.
  ArrayDecoder decoder(encoded_data.mem_, encoded_data.size_);

  // First frame should be lossless (PSNR 99).
  uint32_t duration_ms = 0;
  ASSERT_TRUE(decoder.ReadFrame(&duration_ms));
  EXPECT_EQ(duration_ms, first_frame.duration_ms);
  EXPECT_TRUE(testing::Compare(first_frame.rgb_pixels, decoder.GetPixels(),
                               "logo", 99.f));

  // We don't force lossless on the second frame so it will probably be lossy.
  ASSERT_TRUE(decoder.ReadFrame(&duration_ms));
  EXPECT_EQ(duration_ms, second_frame.duration_ms);
  EXPECT_TRUE(testing::Compare(second_frame.rgb_pixels, decoder.GetPixels(),
                               "background", 35.f));

  // Then ReadFrame() returns false when all data is decoded.
  EXPECT_FALSE(decoder.ReadFrame());
  ASSERT_WP2_OK(decoder.GetStatus());
}

// orientation, is_slow.
typedef std::tuple<Orientation, bool> Param;

class ParameterizedAnimFrameTest : public ::testing::TestWithParam<Param> {
 public:
  // Inits decoded_frame_ buffer.
  WP2Status InitDecodedFrame(uint32_t width, uint32_t height) {
    const Orientation orientation = std::get<0>(GetParam());
    const bool is_slow = std::get<1>(GetParam());

    if (orientation == Orientation::kOriginal ||
        orientation == Orientation::k180) {
      WP2_CHECK_STATUS(decoded_frame_internal_.Resize(width, height));
    } else {
      WP2_CHECK_STATUS(decoded_frame_internal_.Resize(height, width));
    }
    // Allows setting "is_slow".
    WP2_CHECK_STATUS(decoded_frame_.SetExternal(
        decoded_frame_internal_.width, decoded_frame_internal_.height,
        (uint8_t*)decoded_frame_internal_.GetRow(0),
        decoded_frame_internal_.stride, is_slow));
    return WP2_STATUS_OK;
  }

  ArgbBuffer decoded_frame_;

 protected:
  ArgbBuffer decoded_frame_internal_;
};

TEST_P(ParameterizedAnimFrameTest, Blend) {
  Orientation orientation = std::get<0>(GetParam());

  std::vector<ArgbBuffer> buffers;
  ASSERT_WP2_OK(testing::ReadImages(
      {GetTestDataPath("background.png"), GetTestDataPath("logo.png"),
       GetTestDataPath("composited.png")},
      &buffers));
  // Smaller view to speed up the test.
  for (ArgbBuffer& b : buffers) {
    ASSERT_WP2_OK(b.SetView(b, {0, 0, b.width / 2, b.height / 2}));
  }

  const Rectangle window = {0, 0, buffers[0].width, buffers[0].height};

  Vector<Subframe> subframes;
  ASSERT_TRUE(subframes.resize(2));
  Subframe& first_frame = subframes[0];
  first_frame.duration_ms = 10;
  first_frame.window = window;
  first_frame.dispose = true;
  ASSERT_WP2_OK(first_frame.rgb_pixels.SetView(buffers[0]));

  Subframe& second_frame = subframes[1];
  second_frame.duration_ms = 5;
  second_frame.blend = true;
  second_frame.window = window;
  ASSERT_WP2_OK(second_frame.rgb_pixels.SetView(buffers[1]));

  bool has_alpha = false;
  for (const Subframe& subframe : subframes) {
    has_alpha |= subframe.rgb_pixels.HasTransparency();
    if (has_alpha) break;
  }

  EncoderConfig config;
  config.quality = 80;
  config.tile_shape = TILE_SHAPE_SQUARE_256;
  config.decoding_orientation = orientation;
  MemoryWriter encoded_data;

  const Argb38b background_color = kTransparentArgb38b;
  ASSERT_WP2_OK(EncodeHeader(
      config, first_frame.rgb_pixels.width, first_frame.rgb_pixels.height,
      has_alpha, /*is_anim=*/(subframes.size() > 1), /*loop_count=*/1,
      background_color, RGB12b{10, 0, 30}, /*has_icc=*/false,
      /*has_trailing_data=*/false, &encoded_data));

  SubframeEncoder subframe_encoder(GetQualityHint(config.quality), has_alpha);
  ASSERT_WP2_OK(subframe_encoder.EncodeSubframe(config, first_frame,
                                                /*is_last=*/false,
                                                &encoded_data));
  ASSERT_WP2_OK(subframe_encoder.EncodeSubframe(config, second_frame,
                                                /*is_last=*/true,
                                                &encoded_data));
  // Decode.
  ASSERT_WP2_OK(InitDecodedFrame(window.width, window.height));
  ArrayDecoder decoder(encoded_data.mem_, encoded_data.size_,
                       DecoderConfig::kDefault, &decoded_frame_);

  // First frame (shows "background.png").
  uint32_t duration_ms;
  ASSERT_TRUE(decoder.ReadFrame(&duration_ms));
  EXPECT_EQ(duration_ms, first_frame.duration_ms);
  EXPECT_TRUE(testing::Compare(buffers[0], decoded_frame_, "background", 35.,
                               PSNR, config.decoding_orientation));

  // Second frame (shows "logo.png" composited on "background.png" which is
  // "composited.png").
  ASSERT_TRUE(decoder.ReadFrame(&duration_ms));
  EXPECT_EQ(duration_ms, second_frame.duration_ms);
  EXPECT_TRUE(testing::Compare(buffers[2], decoded_frame_, "composited", 35.,
                               PSNR, config.decoding_orientation));

  // Then ReadFrame() returns false when all data is decoded.
  EXPECT_FALSE(decoder.ReadFrame());
  ASSERT_WP2_OK(decoder.GetStatus());
}

TEST_P(ParameterizedAnimFrameTest, BlendWithPreframe) {
  Orientation orientation = std::get<0>(GetParam());

  std::vector<ArgbBuffer> buffers;
  ASSERT_WP2_OK(testing::ReadImages(
      {GetTestDataPath("background.png"), GetTestDataPath("logo.png"),
       GetTestDataPath("composited.png")},
      &buffers));
  // Smaller view to speed up the test.
  for (ArgbBuffer& b : buffers) {
    ASSERT_WP2_OK(b.SetView(
        b, {b.width / 3, b.height / 3, b.width / 3, b.height / 3}));
  }

  const Rectangle window = {0, 0, buffers[0].width, buffers[0].height};

  Vector<Subframe> subframes;
  ASSERT_TRUE(subframes.resize(2));
  Subframe& first_frame = subframes[0];
  first_frame.duration_ms = 0;  // Preframe.
  first_frame.window = window;
  first_frame.dispose = true;
  ASSERT_WP2_OK(first_frame.rgb_pixels.SetView(buffers[0]));

  Subframe& second_frame = subframes[1];
  second_frame.duration_ms = 5;
  second_frame.blend = true;
  second_frame.window = window;
  ASSERT_WP2_OK(second_frame.rgb_pixels.SetView(buffers[1]));

  bool has_alpha = false;
  for (const Subframe& subframe : subframes) {
    has_alpha |= subframe.rgb_pixels.HasTransparency();
    if (has_alpha) break;
  }

  EncoderConfig config;
  config.quality = 80;
  config.tile_shape = TILE_SHAPE_SQUARE_256;
  config.decoding_orientation = orientation;
  MemoryWriter encoded_data;

  const Argb38b background_color = kTransparentArgb38b;
  ASSERT_WP2_OK(EncodeHeader(
      config, first_frame.rgb_pixels.width, first_frame.rgb_pixels.height,
      has_alpha, /*is_anim=*/(subframes.size() > 1), /*loop_count=*/1,
      background_color, RGB12b{10, 0, 30}, /*has_icc=*/false,
      /*has_trailing_data=*/false, &encoded_data));

  SubframeEncoder subframe_encoder(GetQualityHint(config.quality), has_alpha);
  ASSERT_WP2_OK(subframe_encoder.EncodeSubframe(config, first_frame,
                                                /*is_last=*/false,
                                                &encoded_data));
  ASSERT_WP2_OK(subframe_encoder.EncodeSubframe(config, second_frame,
                                                /*is_last=*/true,
                                                &encoded_data));
  // Decode.
  ASSERT_WP2_OK(InitDecodedFrame(window.width, window.height));
  ArrayDecoder decoder(encoded_data.mem_, encoded_data.size_,
                       DecoderConfig::kDefault, &decoded_frame_);

  // Final frame (shows "logo.png" composited on "background.png" which is
  // "composited.png").
  uint32_t duration_ms;
  ASSERT_TRUE(decoder.ReadFrame(&duration_ms));
  EXPECT_EQ(duration_ms, second_frame.duration_ms);
  EXPECT_TRUE(decoder.TryGetFrameDecodedFeatures(0)->is_last);
  EXPECT_TRUE(testing::Compare(buffers[2], decoded_frame_, "composited", 35.,
                               PSNR, config.decoding_orientation));

  // Normal (non incremental) decode.
  ArgbBuffer decoded_first_frame;
  ASSERT_WP2_OK(Decode(encoded_data.mem_, encoded_data.size_,
                       &decoded_first_frame, DecoderConfig::kDefault));
  // We should get the same result as incremental decoding.
  EXPECT_TRUE(testing::Compare(decoded_frame_, decoded_first_frame,
                               "composited", 99.f));

  // Then ReadFrame() returns false when all data is decoded.
  EXPECT_FALSE(decoder.ReadFrame());
  ASSERT_WP2_OK(decoder.GetStatus());
}

static Argb32b RandomColor(UniformIntDistribution* const random) {
  const bool opaque = random->FlipACoin();
  const uint8_t alpha = opaque ? 255 : random->Get(0, 255);
  const float a = alpha / 255.f;
  return {alpha, (uint8_t)(a * random->Get(0, 255)),
          (uint8_t)(a * random->Get(0, 255)),
          (uint8_t)(a * random->Get(0, 255))};
}

static Rectangle RandomRect(uint32_t frame_width, uint32_t frame_height,
                            UniformIntDistribution* const random) {
  const uint32_t x = random->Get<uint32_t>(0, frame_width - 2);
  const uint32_t y = random->Get<uint32_t>(0, frame_height - 2);
  return {x, y, random->Get<uint32_t>(1, frame_width - x - 1),
          random->Get<uint32_t>(1, frame_height - y - 1)};
}

// Set (e.g. to /tmp) to debug RandomFrames test below.
static const char* kDebugOutputPath = "";

static WP2Status FillRandom(const std::string& debug_str,
                            const ArgbBuffer* const previous,
                            UniformIntDistribution* const random,
                            ArgbBuffer* const buffer) {
  if (previous != nullptr && random->FlipACoin()) {
    WP2_CHECK_STATUS(buffer->CopyFrom(*previous));
  } else if (random->FlipACoin()) {
    buffer->Fill(RandomColor(random));
  } else {
    buffer->Fill({0, 0, 0, 0});
  }

  const uint32_t width = buffer->width;
  const uint32_t height = buffer->height;
  const uint32_t num_rects = random->Get(1, 3);
  for (uint32_t j = 0; j < num_rects; ++j) {
    Rectangle window = RandomRect(width, height, random);
    buffer->Fill(window, random->FlipACoin() ? Argb32b{0, 0, 0, 0}
                                             : RandomColor(random));
  }

  if (std::strlen(kDebugOutputPath) > 0) {
    WP2_CHECK_STATUS(WP2::SaveImage(
        *buffer,
        SPrintf("%s/%s.png", kDebugOutputPath, debug_str.c_str()).c_str(),
        true));
  }
  return WP2_STATUS_OK;
}

static WP2Status SaveDebug(const ArgbBuffer& expected_buffer,
                           const ArgbBuffer& decoded_frame,
                           Orientation orientation, uint32_t test_id,
                           uint32_t final_frame_index) {
  if (std::strlen(kDebugOutputPath) > 0) {
    ArgbBuffer expected_rotated;
    WP2_CHECK_STATUS(expected_rotated.CopyFrom(expected_buffer));
    WP2_CHECK_STATUS(RotateBuffer(orientation, &expected_rotated));
    WP2_CHECK_STATUS(
        WP2::SaveImage(decoded_frame,
                       SPrintf("%s/test%d-%d-decoded.png", kDebugOutputPath,
                               test_id, final_frame_index)
                           .c_str(),
                       true));
    WP2_CHECK_STATUS(
        WP2::SaveImage(expected_rotated,
                       SPrintf("%s/test%d-%d-expected.png", kDebugOutputPath,
                               test_id, final_frame_index)
                           .c_str(),
                       true));
  }
  return WP2_STATUS_OK;
}

TEST_P(ParameterizedAnimFrameTest, RandomFrames) {
  const Orientation orientation = std::get<0>(GetParam());
  const bool is_slow = std::get<1>(GetParam());
  const int test_id = (int)orientation * 2 + is_slow;
  UniformIntDistribution random(test_id);

  Vector<Frame> frames;
  const int num_frames = random.Get(20, 40);
  ASSERT_TRUE(frames.resize(num_frames));

  const uint32_t width = random.Get(1, 50);
  const uint32_t height = random.Get(1, 50);

  for (uint32_t i = 0; i < frames.size(); ++i) {
    Frame& frame = frames[i];
    frame.duration_ms = random.Get(1, 100);
    frame.force_dispose = random.FlipACoin();
    frame.force_lossless = random.FlipACoin();
    ASSERT_WP2_OK(frame.rgb_pixels.Resize(width, height));
    const std::string debug_str =
        SPrintf("test%d-frame%d-duration%d-force_dispose%d", test_id, i,
                frame.duration_ms, frame.force_dispose);
    ASSERT_WP2_OK(
        FillRandom(debug_str,
                   /*previous=*/i > 0 ? &(frames[i - 1].rgb_pixels) : nullptr,
                   &random, &frame.rgb_pixels));

    printf("frame %d duration %d force_dispose %d\n", i, frame.duration_ms,
           frame.force_dispose);
  }
  const bool has_alpha = !AllFramesAreOpaque(frames);

  EncoderConfig config;
  config.quality = random.Get(0, 95);
  config.tile_shape = TILE_SHAPE_SQUARE_128;
  config.decoding_orientation = orientation;
  MemoryWriter encoded_data;

  const Argb32b background_color = RandomColor(&random);
  printf("Background_color: Argb32(%d, %d, %d, %d)\n", background_color.a,
         background_color.r, background_color.g, background_color.b);
  ASSERT_WP2_OK(EncodeHeader(
      config, width, height, has_alpha, /*is_anim=*/(frames.size() > 1),
      /*loop_count=*/1, ToArgb38b(background_color),
      /*preview_color=*/RGB12b{10, 0, 30}, /*has_icc=*/false,
      /*has_trailing_data=*/false, &encoded_data));

  FrameEncoder frame_encode(config, ToArgb38b(background_color), frames);
  ASSERT_WP2_OK(frame_encode.Encode(&encoded_data));

  // Decode.
  ASSERT_WP2_OK(InitDecodedFrame(width, height));
  ArrayDecoder decoder(encoded_data.mem_, encoded_data.size_,
                       DecoderConfig::kDefault, &decoded_frame_);

  for (uint32_t i = 0; i < frames.size(); ++i) {
    const Frame& frame = frames[i];

    uint32_t duration_ms = 0;
    ASSERT_TRUE(decoder.ReadFrame(&duration_ms));
    EXPECT_TRUE(testing::Compare(frame.rgb_pixels, decoded_frame_, "fake data",
                                 testing::GetExpectedDistortion(config.quality),
                                 PSNR, orientation))
        << "Error on final frame " << i;

    ASSERT_WP2_OK(
        SaveDebug(frame.rgb_pixels, decoded_frame_, orientation, test_id, i));
  }

  // Then ReadFrame() returns false when all data is decoded.
  EXPECT_FALSE(decoder.ReadFrame());
  ASSERT_WP2_OK(decoder.GetStatus());
}

TEST_P(ParameterizedAnimFrameTest, RandomSubframes) {
  const Orientation orientation = std::get<0>(GetParam());
  const bool is_slow = std::get<1>(GetParam());
  const int test_id = (int)orientation * 2 + is_slow;
  UniformIntDistribution random(test_id);

  Vector<Subframe> subframes;
  const uint32_t num_subframes = random.Get(20, 40);
  ASSERT_TRUE(subframes.resize(num_subframes));

  const uint32_t width = random.Get(1, 50);
  const uint32_t height = random.Get(1, 50);
  bool has_alpha = false;

  EncoderConfig config;
  config.quality = random.Get(50, 95);
  config.tile_shape = TILE_SHAPE_SQUARE_128;
  config.decoding_orientation = orientation;
  MemoryWriter encoded_data;

  const Argb32b background_color = RandomColor(&random);
  printf("Background_color: Argb32(%d, %d, %d, %d)\n", background_color.a,
         background_color.r, background_color.g, background_color.b);

  uint32_t num_consecutive_preframes = 0;
  for (uint32_t i = 0; i < subframes.size(); ++i) {
    Subframe& subframe = subframes[i];
    const bool is_preframe = (i == subframes.size() - 1 ||
                              num_consecutive_preframes >= kMaxNumPreframes)
                                 ? false
                                 : random.FlipACoin();
    num_consecutive_preframes = is_preframe ? num_consecutive_preframes + 1 : 0;
    subframe.duration_ms = is_preframe ? 0 : random.Get(1, 100);
    subframe.dispose =
        (i == 0) ? true
                 : (num_consecutive_preframes > 0) ? false : random.FlipACoin();
    subframe.blend = random.FlipACoin();
    ASSERT_WP2_OK(subframe.rgb_pixels.Resize(width, height));
    const std::string debug_str =
        SPrintf("test%d-frame%d-duration%d-force_dispose%d-blend%d", test_id, i,
                subframe.duration_ms, subframe.dispose, subframe.blend);
    ASSERT_WP2_OK(FillRandom(
        debug_str,
        /*previous=*/i > 0 ? &(subframes[i - 1].rgb_pixels) : nullptr, &random,
        &subframe.rgb_pixels));
    has_alpha |= subframe.rgb_pixels.HasTransparency();

    if (random.FlipACoin()) {
      subframe.window = {0, 0, width, height};
    } else {
      subframe.window = RandomRect(width, height, &random);
    }

    printf("frame %d duration %d dispose %d blend %d\n", i,
           subframe.duration_ms, subframe.dispose, subframe.blend);
  }

  ASSERT_WP2_OK(EncodeHeader(
      config, width, height, has_alpha, /*is_anim=*/(subframes.size() > 1),
      /*loop_count=*/1, ToArgb38b(background_color),
      /*preview_color=*/RGB12b{10, 0, 30}, /*has_icc=*/false,
      /*has_trailing_data=*/false, &encoded_data));

  SubframeEncoder subframe_encoder(GetQualityHint(config.quality), has_alpha);

  for (uint32_t i = 0; i < subframes.size(); ++i) {
    const bool is_last = (i == num_subframes - 1);
    ASSERT_WP2_OK(subframe_encoder.EncodeSubframe(config, subframes[i], is_last,
                                                  &encoded_data));
  }

  // Decode.
  ASSERT_WP2_OK(InitDecodedFrame(width, height));
  ArrayDecoder decoder(encoded_data.mem_, encoded_data.size_,
                       DecoderConfig::kDefault, &decoded_frame_);

  // Use the classic (non incremental) decoder to decode the first frame.
  ArgbBuffer decoded_first_frame;
  ASSERT_WP2_OK(Decode(encoded_data.mem_, encoded_data.size_,
                       &decoded_first_frame));

  ArgbBuffer expected_buffer;
  ASSERT_WP2_OK(expected_buffer.Resize(width, height));
  int final_frame_index = 0;
  for (uint32_t i = 0; i < subframes.size(); ++i) {
    const Subframe& subframe = subframes[i];
    if (i == 0 || subframe.dispose) expected_buffer.Fill(background_color);
    if (subframe.blend) {
      ASSERT_WP2_OK(
          expected_buffer.CompositeUnder(subframe.rgb_pixels, subframe.window));
    } else {
      ArgbBuffer view1, view2;
      ASSERT_WP2_OK(view1.SetView(expected_buffer, subframe.window));
      ASSERT_WP2_OK(view2.SetView(subframe.rgb_pixels, subframe.window));
      ASSERT_WP2_OK(view1.CopyFrom(view2));
    }

    if (subframe.duration_ms == 0) continue;

    uint32_t duration_ms = 0;
    ASSERT_TRUE(decoder.ReadFrame(&duration_ms));
    ASSERT_TRUE(decoded_frame_.IsView() &&
                decoded_frame_.GetRow(0) == decoded_frame_internal_.GetRow(0));
    EXPECT_TRUE(testing::Compare(expected_buffer, decoded_frame_, "fake data",
                                 testing::GetExpectedDistortion(config.quality),
                                 PSNR, orientation))
        << "Error on final frame " << final_frame_index;

    if (final_frame_index == 0) {
      // Compare output of classic decoder (which only decodes the first frame).
      EXPECT_TRUE(testing::Compare(decoded_frame_, decoded_first_frame,
                                   "fake data", 99.))
          << "Error on first frame";
    }

    ASSERT_WP2_OK(SaveDebug(expected_buffer, decoded_frame_, orientation,
                            test_id, final_frame_index));
    ++final_frame_index;
  }

  // Then ReadFrame() returns false when all data is decoded.
  EXPECT_FALSE(decoder.ReadFrame());
  ASSERT_WP2_OK(decoder.GetStatus());
}

INSTANTIATE_TEST_SUITE_P(
    ParameterizedAnimFrameTestInstantiation, ParameterizedAnimFrameTest,
    ::testing::Combine(::testing::Values(Orientation::kOriginal,
                                         Orientation::k90Clockwise,
                                         Orientation::k180),
                       ::testing::Values(true, false)));  // is_slow

}  // namespace
}  // namespace WP2
