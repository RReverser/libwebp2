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

// Test example utils.

#include <cstdio>
#include <string>
#include <tuple>

#include "examples/example_utils.h"
#include "imageio/image_dec.h"
#include "imageio/image_enc.h"
#include "include/helpers.h"
#include "src/dec/preview/preview_dec.h"
#include "src/enc/preview/preview_enc.h"
#include "src/utils/data_source.h"
#include "src/wp2/base.h"
#include "src/wp2/decode.h"
#include "src/wp2/encode.h"

namespace WP2 {
namespace {

constexpr float kMinExpectedSimilarity = 13.f;

//------------------------------------------------------------------------------

void ShowOutput(const std::string& file_name, const ArgbBuffer& original_image,
                const ArgbBuffer& decompressed_preview, const float disto[5]) {
#ifdef WP2_TEST_PREVIEW_API_SAVE_OUTPUT
  std::string original = WP2RemoveFileExtension(file_name) + "_original.png";
  EXPECT_WP2_OK(SaveImage(original_image, ("/tmp/" + original).c_str(), true));
  std::string decompressed = WP2RemoveFileExtension(file_name) + "_preview.png";
  EXPECT_WP2_OK(
      SaveImage(decompressed_preview, ("/tmp/" + decompressed).c_str(), true));
  fprintf(stderr,
          "Distortion: A: %2.2f  R: %2.2f  G: %2.2f  B: %2.2f  Total: %2.2f\n",
          disto[0], disto[1], disto[2], disto[3], disto[4]);
#else
  (void)file_name;
  (void)original_image;
  (void)decompressed_preview;
  (void)disto;
#endif
}

//------------------------------------------------------------------------------

class PreviewTest
    : public ::testing::TestWithParam<
          std::tuple<const char*, float, int, PreviewConfig::AnalysisMethod>> {
};

TEST_P(PreviewTest, SimpleTest) {
  const std::string file_name = std::get<0>(GetParam());
  const std::string file_path = testing::GetTestDataPath(file_name);
  const float quality = std::get<1>(GetParam());
  const int speed = std::get<2>(GetParam());
  const PreviewConfig::AnalysisMethod analysis_method = std::get<3>(GetParam());

  ArgbBuffer original_image;
  ASSERT_WP2_OK(ReadImage(file_path.c_str(), &original_image));

  MemoryWriter writer;
  PreviewConfig config(quality, speed);
  config.analysis_method = analysis_method;
  ASSERT_WP2_OK(EncodePreview(original_image, config, &writer));

  ArgbBuffer decompressed_preview(original_image.format);
  ASSERT_WP2_OK(
      decompressed_preview.Resize(original_image.width, original_image.height));
  ASSERT_WP2_OK(
      DecodePreview(writer.mem_, writer.size_, &decompressed_preview));

  float disto[5];
  ASSERT_WP2_OK(
      decompressed_preview.GetDistortion(original_image, PSNR, disto));
  EXPECT_GE(disto[4], kMinExpectedSimilarity);

  ShowOutput(file_name, original_image, decompressed_preview, disto);
}

INSTANTIATE_TEST_SUITE_P(
    PreviewTestInstantiationFast, PreviewTest,
    ::testing::Combine(
        ::testing::Values("source1_64x48.png"),
        ::testing::Values(40.f),  // Quality
        ::testing::Values(3),     // Speed
        ::testing::Values(
            PreviewConfig::AnalysisMethod::kEdgeSelectionAndRepulsion)));

// This one might take ~10min so it is disabled.
// Can still be run with flag --test_arg=--gunit_also_run_disabled_tests
INSTANTIATE_TEST_SUITE_P(
    DISABLED_PreviewTestInstantiation, PreviewTest,
    ::testing::Combine(
        ::testing::Values("source0.ppm", "source1.png"),
        ::testing::Values(10.f, 50.f),  // Quality
        ::testing::Values(2, 5),        // Speed
        ::testing::Values(
            PreviewConfig::AnalysisMethod::kEdgeSelectionAndRepulsion,
            PreviewConfig::AnalysisMethod::kColorDiffMaximization,
            PreviewConfig::AnalysisMethod::kRandomEdgeSelection)));

//------------------------------------------------------------------------------

TEST(PreviewTest, NoOptim) {
  const std::string file_name = "alpha_ramp.png";
  const std::string file_path = testing::GetTestDataPath(file_name);

  ArgbBuffer original_image;
  ASSERT_WP2_OK(ReadImage(file_path.c_str(), &original_image));

  MemoryWriter writer;
  PreviewConfig config;
  config.num_vertices = 100;
  config.num_colors = 14;
  config.grid_size = 32;
  config.grid_density = 1.0;
  config.blur_radius = 4;
  config.num_iterations = 0;
  ASSERT_WP2_OK(EncodePreview(original_image, config, &writer));

  ArgbBuffer decompressed_preview(original_image.format);
  ASSERT_WP2_OK(
      decompressed_preview.Resize(original_image.width, original_image.height));
  ASSERT_WP2_OK(
      DecodePreview(writer.mem_, writer.size_, &decompressed_preview));

  float disto[5];
  ASSERT_WP2_OK(
      decompressed_preview.GetDistortion(original_image, PSNR, disto));
  EXPECT_GE(disto[4], kMinExpectedSimilarity);

  ShowOutput(file_name, original_image, decompressed_preview, disto);
}

//------------------------------------------------------------------------------

TEST(PreviewTest, EncodeDecode) {
  const std::string file_name = "test_exif_xmp.webp";
  const std::string file_path = testing::GetTestDataPath(file_name);

  ArgbBuffer original_image;
  ASSERT_WP2_OK(ReadImage(file_path.c_str(), &original_image));

  EncoderConfig config;
  config.quality = 40.f;
  config.speed = 2;
  config.create_preview = true;
  MemoryWriter memory_writer;
  EXPECT_WP2_OK(Encode(original_image, &memory_writer, config));

  ArgbBuffer decompressed_preview(original_image.format);
  ASSERT_WP2_OK(ExtractPreview(memory_writer.mem_, memory_writer.size_,
                               &decompressed_preview));

  float disto[5];
  ASSERT_WP2_OK(
      decompressed_preview.GetDistortion(original_image, PSNR, disto));
  EXPECT_GE(disto[4], kMinExpectedSimilarity);

  ShowOutput(file_name, original_image, decompressed_preview, disto);
}

//------------------------------------------------------------------------------

TEST(PreviewTest, CustomSize) {
  const std::string file_name = "test_exif_xmp.webp";
  const std::string file_path = testing::GetTestDataPath(file_name);

  ArgbBuffer original_image;
  ASSERT_WP2_OK(ReadImage(file_path.c_str(), &original_image));

  EncoderConfig config;
  config.quality = 20.f;
  config.speed = 1;
  config.create_preview = true;
  MemoryWriter memory_writer;
  EXPECT_WP2_OK(Encode(original_image, &memory_writer, config));

  ArgbBuffer decompressed_preview(original_image.format);
  for (uint32_t width : {1u, 2u, 256u}) {
    for (uint32_t height : {1u, 2u, 256u}) {
      ASSERT_WP2_OK(decompressed_preview.Resize(width, height));
      ASSERT_WP2_OK(ExtractPreview(memory_writer.mem_, memory_writer.size_,
                                   &decompressed_preview));
    }
  }
}

//------------------------------------------------------------------------------

// To get coverage on the preview skipping code.
TEST(PreviewTest, Skipped) {
  ArgbBuffer original_image;
  ASSERT_WP2_OK(ReadImage(testing::GetTestDataPath("source1_64x48.png").c_str(),
                          &original_image));

  EncoderConfig config;
  config.quality = 0.f;  // Lowest preview settings; it will be skipped anyway.
  config.speed = 0;
  config.create_preview = true;
  MemoryWriter memory_writer;
  EXPECT_WP2_OK(Encode(original_image, &memory_writer, config));

  ArgbBuffer ignored_output;
  EXPECT_WP2_OK(
      Decode(memory_writer.mem_, memory_writer.size_, &ignored_output));
}

// Test the 'binary preview' encoder option
TEST(PreviewTest, PreviewData) {
  ArgbBuffer original_image;
  ASSERT_WP2_OK(ReadImage(testing::GetTestDataPath("source1_64x48.png").c_str(),
                          &original_image));
  ArgbBuffer ignored_output;
  MemoryWriter writer;
  EncoderConfig config;
  config.preview_size = 123;
  config.preview = nullptr;
  EXPECT_EQ(Encode(original_image, &writer, config),
                   WP2_STATUS_INVALID_CONFIGURATION);

  const char kPreviewdata[] = "this is a nice thumbnail!";
  config.preview = (const uint8_t*)kPreviewdata;
  config.preview_size = strlen(kPreviewdata);
  writer.Reset();
  EXPECT_WP2_OK(Encode(original_image, &writer, config));
  EXPECT_WP2_OK(Decode(writer.mem_, writer.size_, &ignored_output));

  config.preview_size = 0;
  writer.Reset();
  EXPECT_WP2_OK(Encode(original_image, &writer, config));
  EXPECT_WP2_OK(Decode(writer.mem_, writer.size_, &ignored_output));
}

//------------------------------------------------------------------------------

}  // namespace
}  // namespace WP2
