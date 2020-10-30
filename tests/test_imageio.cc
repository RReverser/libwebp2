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

// Test ImageIo.

#include <string>

#include "extras/ccsp_imageio.h"
#include "imageio/anim_image_dec.h"
#include "imageio/image_dec.h"
#include "imageio/image_enc.h"
#include "include/helpers.h"

namespace WP2 {
namespace {

//------------------------------------------------------------------------------

TEST(ReadImageTest, Metadata) {
  ArgbBuffer buffer;

  // Read image with metadata.
  ASSERT_WP2_OK(ReadImage(
      testing::GetTestDataPath("test_exif_xmp.webp").c_str(), &buffer));
  EXPECT_TRUE(!buffer.metadata.IsEmpty());
  // Read image without metadata.
  ASSERT_WP2_OK(ReadImage(testing::GetTestDataPath("source1_64x48.png").c_str(),
                          &buffer));
  // Expect nothing left from first image.
  EXPECT_TRUE(buffer.metadata.IsEmpty());
}

TEST(ReadImageTest, NoStdoutFileSize) {
  ArgbBuffer buffer;
  ASSERT_WP2_OK(buffer.Resize(1, 1));
  buffer.Fill(Argb32b{0, 0, 0, 0});

  size_t file_size;
  ASSERT_EQ(
      SaveImage(buffer, "-", /*overwrite=*/true, FileFormat::AUTO, &file_size),
      WP2_STATUS_INVALID_PARAMETER);
}

//------------------------------------------------------------------------------
// Test all formats that support read and lossless save.

class WriteImageTest : public ::testing::TestWithParam<const char*> {};

TEST_P(WriteImageTest, Argb) {
  const char* const file_name = GetParam();
  Data data;
  ASSERT_WP2_OK(
      IoUtilReadFile(testing::GetTestDataPath(file_name).c_str(), &data));
  const FileFormat format = GuessImageFormat(data.bytes, data.size);

  ArgbBuffer original;
  if (IsCustomColorSpace(format)) {
    EXPECT_EQ(ReadImage(data.bytes, data.size, &original),
              WP2_STATUS_UNSUPPORTED_FEATURE);
    return;
  }
  ASSERT_WP2_OK(ReadImage(data.bytes, data.size, &original));

  const std::string temp_file_path = testing::GetTempDataPath(file_name);
  std::remove(temp_file_path.c_str());  // Make sure it doesn't exist.
  size_t written_file_size, read_file_size;
  ASSERT_WP2_OK(SaveImage(original, temp_file_path.c_str(), /*overwrite=*/true,
                          format, &written_file_size));

  ArgbBuffer comparison;
  ASSERT_WP2_OK(
      ReadImage(temp_file_path.c_str(), &comparison, &read_file_size));
  EXPECT_EQ(written_file_size, read_file_size);

  // Reading a custom color space into an ArgbBuffer will generate loss.
  const bool is_ccsp = IsCustomColorSpace(GetFormatFromExtension(file_name));
  const float expected_distortion = is_ccsp ? 40.f : 99.f;
  EXPECT_TRUE(
      testing::Compare(original, comparison, file_name, expected_distortion));
}

TEST_P(WriteImageTest, Ccsp) {
  const char* const file_name = GetParam();
  Data data;
  ASSERT_WP2_OK(
      IoUtilReadFile(testing::GetTestDataPath(file_name).c_str(), &data));
  const FileFormat format = GuessImageFormat(data.bytes, data.size);

  YUVPlane original;
  CSPMtx ccsp_to_rgb = {};
  ASSERT_WP2_OK(ReadImage(data.bytes, data.size, &original, &ccsp_to_rgb));

  uint32_t file_num_bits;
  if (IsCustomColorSpace(format)) {
    ASSERT_GE(ccsp_to_rgb.shift, CSPTransform::kMtxShift);
    file_num_bits = CCSPImageReader::kMinBitDepth +
                    (ccsp_to_rgb.shift - CSPTransform::kMtxShift);
  } else {
    file_num_bits = 8;
  }

  uint32_t bit_depth;
  ASSERT_WP2_OK(ReadBitDepth(data.bytes, data.size, &bit_depth));
  ASSERT_EQ(file_num_bits, bit_depth);

  ASSERT_EQ(IsChromaSubsampled(format), original.IsDownsampled());
  YUVPlane upsampled;
  if (original.IsDownsampled()) {
    // Upscale with no interpolation to keep the original samples.
    ASSERT_WP2_OK(upsampled.UpsampleFrom(original, SamplingTaps::kUpNearest));
  } else {
    ASSERT_WP2_OK(upsampled.SetView(
        original, {0, 0, original.GetWidth(), original.GetHeight()}));
  }

  const std::string temp_file_path = testing::GetTempDataPath(file_name);
  std::remove(temp_file_path.c_str());  // Make sure it doesn't exist.
  ASSERT_WP2_OK(SaveImage(upsampled, ccsp_to_rgb,
                          file_num_bits, temp_file_path.c_str(),
                          /*overwrite=*/true, format, Metadata(),
                          SamplingTaps::kDownAvg));

  YUVPlane comparison;
  ASSERT_WP2_OK(ReadImage(temp_file_path.c_str(), &comparison, &ccsp_to_rgb));

  // Reading an Argb color space into a YUVPlane will not generate loss.
  EXPECT_TRUE(testing::Compare(original, comparison, file_num_bits, file_name));
}

INSTANTIATE_TEST_SUITE_P(
    WriteImageTestInstantiation, WriteImageTest,
    ::testing::Values("ccsp/source3_C420p8.y4m"));

//------------------------------------------------------------------------------

// Simple test to exercise a corner case (1px wide, non-multiple-of-2 tall).
TEST(Y4MTest, Nearest) {
  const char* const file_name = "ccsp/source3_C420p8.y4m";

  YUVPlane original;
  CSPMtx ccsp_to_rgb = {};
  ASSERT_WP2_OK(
      ReadImage(testing::GetTestDataPath(file_name).c_str(), &original,
                &ccsp_to_rgb));
  ASSERT_TRUE(original.IsDownsampled());
  // Crop.
  ASSERT_WP2_OK(original.Y.SetView(original.Y, {0, 0, 1, 41}));
  ASSERT_WP2_OK(original.U.SetView(original.U, {0, 0, 1, 21}));
  ASSERT_WP2_OK(original.V.SetView(original.V, {0, 0, 1, 21}));
  original.A.Clear();

  YUVPlane upsampled;
  ASSERT_WP2_OK(upsampled.UpsampleFrom(original, SamplingTaps::kUpNearest));
  ASSERT_WP2_OK(upsampled.Downsample(SamplingTaps::kDownAvg));

  uint32_t bit_depth;
  ASSERT_WP2_OK(
      ReadBitDepth(testing::GetTestDataPath(file_name).c_str(), &bit_depth));
  EXPECT_TRUE(testing::Compare(original, upsampled, bit_depth, file_name));
}

//------------------------------------------------------------------------------

TEST(WriteImageTest, NegativeArgbBuffer) {
  ArgbBuffer buffer;
  const std::string dst = testing::GetTempDataPath("some_path.some_extension");

  ASSERT_TRUE(buffer.IsEmpty());
  EXPECT_EQ(SaveImage(buffer, dst.c_str()), WP2_STATUS_INVALID_PARAMETER);

  ASSERT_WP2_OK(buffer.Resize(1, 1));
  buffer.Fill({0, 0, 1, 1}, Argb32b({0x00u, 0x00u, 0x00u, 0x00u}));
  ASSERT_FALSE(buffer.IsEmpty());
  ASSERT_EQ(SaveImage(buffer, nullptr), WP2_STATUS_NULL_PARAMETER);
}

TEST(WriteImageTest, NegativeYUVPlane) {
  YUVPlane buffer;
  const CSPMtx good_ccsp_to_rgb({1, 0, 0, 0, 1, 0, 0, 0, 1}, 0);
  const CSPMtx bad_ccsp_to_rgb({1, 0, 0, 0, 1, 0, 0, 0, 1}, 1000);
  const uint32_t file_num_bits = 8;
  const uint32_t bad_file_num_bits = 2;
  const std::string dst = testing::GetTempDataPath("some_path.some_extension");

  ASSERT_TRUE(buffer.IsEmpty());
  EXPECT_EQ(SaveImage(buffer, good_ccsp_to_rgb, file_num_bits, dst.c_str()),
            WP2_STATUS_INVALID_PARAMETER);

  ASSERT_WP2_OK(buffer.Resize(1, 1));
  buffer.Fill(Rectangle{0, 0, 1, 1}, Ayuv38b{0xffu, 0x00u, 0x00u, 0x00u});
  ASSERT_FALSE(buffer.IsEmpty());
  ASSERT_EQ(SaveImage(buffer, good_ccsp_to_rgb, file_num_bits, nullptr),
            WP2_STATUS_NULL_PARAMETER);

  const FileFormat format = FileFormat::Y4M_444;
  ASSERT_EQ(SaveImage(buffer, bad_ccsp_to_rgb, file_num_bits, dst.c_str(),
                      /*overwrite=*/true, format),
            WP2_STATUS_INVALID_PARAMETER);
  ASSERT_EQ(SaveImage(buffer, good_ccsp_to_rgb, bad_file_num_bits, dst.c_str(),
                      /*overwrite=*/true, format),
            WP2_STATUS_INVALID_PARAMETER);
}

//------------------------------------------------------------------------------

}  // namespace
}  // namespace WP2
