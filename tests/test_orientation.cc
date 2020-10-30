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

// Rotation functions test. Input should be non-square images.

#include <algorithm>

#include "imageio/image_dec.h"
#include "include/helpers.h"
#include "src/utils/orientation.h"
#include "src/wp2/base.h"
#include "src/wp2/decode.h"
#include "src/wp2/encode.h"

namespace WP2 {
namespace {

//------------------------------------------------------------------------------

TEST(OrientationTest, ArgbBuffer) {
  const std::string& file_name = "source1_64x48.png";
  ArgbBuffer original;
  ArgbBuffer oriented;
  ASSERT_WP2_OK(
      ReadImage(testing::GetTestDataPath(file_name).c_str(), &original));
  ASSERT_WP2_OK(oriented.CopyFrom(original));

  ASSERT_WP2_OK(RotateBuffer(Orientation::kOriginal, &oriented));
  EXPECT_TRUE(testing::Compare(original, oriented, file_name));

  ASSERT_WP2_OK(RotateBuffer(Orientation::k90Clockwise, &oriented));
  ASSERT_WP2_OK(RotateBuffer(Orientation::k270Clockwise, &oriented));
  EXPECT_TRUE(testing::Compare(original, oriented, file_name));

  ASSERT_WP2_OK(RotateBuffer(Orientation::k180, &oriented));
  ASSERT_WP2_OK(RotateBuffer(Orientation::k270Clockwise, &oriented));
  ASSERT_WP2_OK(RotateBuffer(Orientation::k270Clockwise, &oriented));
  EXPECT_TRUE(testing::Compare(original, oriented, file_name));

  ASSERT_WP2_OK(RotateBuffer(Orientation::k90Clockwise, &oriented));
  ASSERT_WP2_OK(RotateBuffer(Orientation::k180, &oriented));
  ASSERT_WP2_OK(RotateBuffer(Orientation::k90Clockwise, &oriented));
  EXPECT_TRUE(testing::Compare(original, oriented, file_name));

  // Make sure a rotation happens.
  ASSERT_WP2_OK(RotateBuffer(Orientation::k90Clockwise, &oriented));
  EXPECT_FALSE(testing::Compare(original, oriented, file_name));

  ASSERT_WP2_OK(oriented.CopyFrom(original));
  ASSERT_WP2_OK(RotateBuffer(Orientation::k180, &oriented));
  EXPECT_FALSE(testing::Compare(original, oriented, file_name));

  ASSERT_WP2_OK(oriented.CopyFrom(original));
  ASSERT_WP2_OK(RotateBuffer(Orientation::k270Clockwise, &oriented));
  EXPECT_FALSE(testing::Compare(original, oriented, file_name));

  ArgbBuffer view;
  ASSERT_WP2_OK(view.SetView(oriented));
  ASSERT_WP2_OK(RotateBuffer(Orientation::kOriginal, &view));
  ASSERT_NE(RotateBuffer(Orientation::k90Clockwise, &view), WP2_STATUS_OK);
  ASSERT_WP2_OK(RotateBuffer(Orientation::k180, &view));
  ASSERT_NE(RotateBuffer(Orientation::k270Clockwise, &view), WP2_STATUS_OK);

  ASSERT_NE(RotateBuffer(Orientation::kOriginal, original, &view),
            WP2_STATUS_OK);
  ASSERT_WP2_OK(RotateBuffer(Orientation::k90Clockwise, original, &view));
  ASSERT_NE(RotateBuffer(Orientation::k180, original, &view), WP2_STATUS_OK);
  ASSERT_WP2_OK(RotateBuffer(Orientation::k270Clockwise, original, &view));
}

//------------------------------------------------------------------------------

TEST(OrientationTest, YUVPlane) {
  const std::string& file_name = "source1_64x48.png";
  ArgbBuffer rgb;
  ASSERT_WP2_OK(ReadImage(testing::GetTestDataPath(file_name).c_str(), &rgb));
  CSPTransform csp_transform;
  const uint32_t bit_depth = csp_transform.GetYUVPrecisionBits() + 1;
  YUVPlane original;
  ASSERT_WP2_OK(original.Import(rgb, rgb.HasTransparency(), csp_transform,
                                /*resize_if_needed=*/true));
  YUVPlane oriented;
  ASSERT_WP2_OK(oriented.Copy(original, /*resize_if_needed=*/true));

  ASSERT_WP2_OK(RotateBuffer(Orientation::kOriginal, &oriented));
  EXPECT_TRUE(testing::Compare(original, oriented, bit_depth, file_name));

  ASSERT_WP2_OK(RotateBuffer(Orientation::k90Clockwise, &oriented));
  ASSERT_WP2_OK(RotateBuffer(Orientation::k270Clockwise, &oriented));
  EXPECT_TRUE(testing::Compare(original, oriented, bit_depth, file_name));

  ASSERT_WP2_OK(RotateBuffer(Orientation::k180, &oriented));
  ASSERT_WP2_OK(RotateBuffer(Orientation::k270Clockwise, &oriented));
  ASSERT_WP2_OK(RotateBuffer(Orientation::k270Clockwise, &oriented));
  EXPECT_TRUE(testing::Compare(original, oriented, bit_depth, file_name));

  ASSERT_WP2_OK(RotateBuffer(Orientation::k90Clockwise, &oriented));
  ASSERT_WP2_OK(RotateBuffer(Orientation::k180, &oriented));
  ASSERT_WP2_OK(RotateBuffer(Orientation::k90Clockwise, &oriented));
  EXPECT_TRUE(testing::Compare(original, oriented, bit_depth, file_name));

  // Make sure a rotation happens.
  ASSERT_WP2_OK(RotateBuffer(Orientation::k90Clockwise, &oriented));
  EXPECT_FALSE(testing::Compare(original, oriented, bit_depth, file_name));

  ASSERT_WP2_OK(oriented.Copy(original, /*resize_if_needed=*/true));
  ASSERT_WP2_OK(RotateBuffer(Orientation::k180, &oriented));
  EXPECT_FALSE(testing::Compare(original, oriented, bit_depth, file_name));

  ASSERT_WP2_OK(oriented.Copy(original, /*resize_if_needed=*/true));
  ASSERT_WP2_OK(RotateBuffer(Orientation::k270Clockwise, &oriented));
  EXPECT_FALSE(testing::Compare(original, oriented, bit_depth, file_name));

  YUVPlane view;
  ASSERT_WP2_OK(view.SetView(
      oriented, {0, 0, oriented.GetWidth(), oriented.GetHeight()}));
  ASSERT_WP2_OK(RotateBuffer(Orientation::kOriginal, &view));
  ASSERT_NE(RotateBuffer(Orientation::k90Clockwise, &view), WP2_STATUS_OK);
  ASSERT_WP2_OK(RotateBuffer(Orientation::k180, &view));
  ASSERT_NE(RotateBuffer(Orientation::k270Clockwise, &view), WP2_STATUS_OK);

  ASSERT_NE(RotateBuffer(Orientation::kOriginal, original, &view),
            WP2_STATUS_OK);
  ASSERT_WP2_OK(RotateBuffer(Orientation::k90Clockwise, original, &view));
  ASSERT_NE(RotateBuffer(Orientation::k180, original, &view), WP2_STATUS_OK);
  ASSERT_WP2_OK(RotateBuffer(Orientation::k270Clockwise, original, &view));
}

//------------------------------------------------------------------------------

TEST(OrientationTest, EncodeDecode) {
  const std::string& file_name = "source1_64x48.png";
  ArgbBuffer original;
  ASSERT_WP2_OK(
      ReadImage(testing::GetTestDataPath(file_name).c_str(), &original));

  EncoderConfig config;
  for (Orientation orientation :
       {Orientation::kOriginal, Orientation::k90Clockwise, Orientation::k180,
        Orientation::k270Clockwise}) {
    config.decoding_orientation = orientation;

    for (float quality : {0.f, 50.f, 100.f}) {
      config.quality = quality;

      MemoryWriter memory_writer;
      EXPECT_WP2_OK(Encode(original, &memory_writer, config));

      ArgbBuffer decoded;
      EXPECT_WP2_OK(Decode(memory_writer.mem_, memory_writer.size_, &decoded));
      ASSERT_WP2_OK(RotateBuffer(GetInverseOrientation(orientation), &decoded));
      EXPECT_TRUE(testing::Compare(original, decoded, file_name,
                                   testing::GetExpectedDistortion(config)));
    }
  }
}

//------------------------------------------------------------------------------

TEST(OrientationTest, Incremental) {
  const std::string& file_name = "source1_64x48.png";
  const uint32_t incr_size_step = 10;
  ArgbBuffer original;
  ASSERT_WP2_OK(
      ReadImage(testing::GetTestDataPath(file_name).c_str(), &original));

  EncoderConfig config;
  for (Orientation orientation :
       {Orientation::kOriginal, Orientation::k90Clockwise, Orientation::k180,
        Orientation::k270Clockwise}) {
    config.decoding_orientation = orientation;

    for (float quality : {25.f, 100.f}) {
      config.quality = quality;

      MemoryWriter memory_writer;
      EXPECT_WP2_OK(Encode(original, &memory_writer, config));

      ArgbBuffer decoded;
      ArrayDecoder idec(DecoderConfig::kDefault, &decoded);
      Rectangle decoded_area(0, 0, 0, 0);
      Rectangle previous_decoded_area(0, 0, 0, 0);
      ArgbBuffer desoriented;

      size_t available_input_size = 0;
      while (available_input_size < memory_writer.size_) {
        available_input_size = std::min(available_input_size + incr_size_step,
                                        memory_writer.size_);
        idec.SetInput(memory_writer.mem_, available_input_size);
        idec.ReadFrame();
        ASSERT_FALSE(idec.Failed());

        decoded_area = idec.GetDecodedArea();
        if (orientation == Orientation::kOriginal ||
            orientation == Orientation::k270Clockwise) {
          EXPECT_EQ(decoded_area.x, 0u);
        } else {
          if (previous_decoded_area.width > 0) {
            EXPECT_LE(decoded_area.x, previous_decoded_area.x);
          }
          if (decoded.width > 0) {
            EXPECT_LT(decoded_area.x, decoded.width);
          }
        }
        if (orientation == Orientation::kOriginal ||
            orientation == Orientation::k90Clockwise) {
          EXPECT_EQ(decoded_area.y, 0u);
        } else {
          if (previous_decoded_area.height > 0) {
            EXPECT_LE(decoded_area.y, previous_decoded_area.y);
          }
          if (decoded.height > 0) {
            EXPECT_LT(decoded_area.y, decoded.height);
          }
        }
        EXPECT_GE(decoded_area.width, previous_decoded_area.width);
        EXPECT_LE(decoded_area.width, decoded.width);
        EXPECT_GE(decoded_area.height, previous_decoded_area.height);
        EXPECT_LE(decoded_area.height, decoded.height);

        previous_decoded_area = decoded_area;
      }
      while (idec.ReadFrame()) continue;  // Flush anything left (metadata...).
      ASSERT_WP2_OK(idec.GetStatus());
      EXPECT_EQ(decoded_area.x, 0u);
      EXPECT_EQ(decoded_area.y, 0u);
      EXPECT_EQ(decoded_area.width, decoded.width);
      EXPECT_EQ(decoded_area.height, decoded.height);

      ASSERT_WP2_OK(RotateBuffer(GetInverseOrientation(orientation), &decoded));
      EXPECT_TRUE(testing::Compare(original, decoded, file_name,
                                   testing::GetExpectedDistortion(config)));
    }
  }
}

//------------------------------------------------------------------------------

TEST(OrientationTest, Preview) {
  const std::string file_name = "test_exif_xmp.webp";
  const std::string file_path = testing::GetTestDataPath(file_name);

  ArgbBuffer original_image;
  ASSERT_WP2_OK(ReadImage(file_path.c_str(), &original_image));
  original_image.height -= 2;  // Make sure it's not square.

  EncoderConfig config;
  config.quality = 10.f;
  config.speed = 1;
  config.create_preview = true;

  ArgbBuffer preview_reference(original_image.format);
  {
    MemoryWriter memory_writer;
    EXPECT_WP2_OK(Encode(original_image, &memory_writer, config));

    ASSERT_WP2_OK(ExtractPreview(memory_writer.mem_, memory_writer.size_,
                                 &preview_reference));
  }

  for (Orientation orientation : {Orientation::k90Clockwise, Orientation::k180,
                                  Orientation::k270Clockwise}) {
    config.decoding_orientation = orientation;

    MemoryWriter memory_writer;
    EXPECT_WP2_OK(Encode(original_image, &memory_writer, config));

    ArgbBuffer preview(original_image.format);  // Use original size.
    ASSERT_WP2_OK(
        ExtractPreview(memory_writer.mem_, memory_writer.size_, &preview));

    ASSERT_WP2_OK(RotateBuffer(GetInverseOrientation(orientation), &preview));
    EXPECT_TRUE(testing::Compare(preview_reference, preview, file_name));

    // Also test custom sizes.
    for (uint32_t width : {1u, 2u, 128u}) {
      for (uint32_t height : {1u, 2u, 128u}) {
        ASSERT_WP2_OK(preview.Resize(width, height));
        ASSERT_WP2_OK(
            ExtractPreview(memory_writer.mem_, memory_writer.size_, &preview));
      }
    }
  }
}

//------------------------------------------------------------------------------

}  // namespace
}  // namespace WP2
