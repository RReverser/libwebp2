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

// Image analysis/preprocessing tests.

#include <unordered_set>

#include "examples/example_utils.h"
#include "imageio/image_dec.h"
#include "include/helpers.h"
#include "src/enc/analysis.h"
#include "src/wp2/encode.h"

namespace WP2 {
namespace {

uint32_t CountColors(const ArgbBuffer& buffer) {
  std::unordered_set<uint32_t> colors;
  for (uint32_t y = 0; y < buffer.height; ++y) {
    const uint32_t* const row = (uint32_t*)buffer.GetRow(y);
    for (uint32_t x = 0; x < buffer.width; ++x) {
      colors.insert(row[x]);
    }
  }
  return colors.size();
}

TEST(PreprocessNearLosslessTest, Photo) {
  ArgbBuffer original;
  ASSERT_WP2_OK(
      ReadImage(testing::GetTestDataPath("background.png").c_str(), &original));

  // Encode losslessly.
  EncoderConfig encoder_config;
  encoder_config.quality = 100;  // lossless
  MemoryWriter memory_writer;
  ASSERT_WP2_OK(Encode(original, &memory_writer, encoder_config));
  uint32_t previous_size = memory_writer.size_;

  uint32_t previous_num_colors = CountColors(original);

  // Encode at various near lossless levels.
  for (uint32_t q : {99, 98, 97, 96}) {
    SCOPED_TRACE(SPrintf("Quality: %d", q));
    encoder_config.quality = q;
    ArgbBuffer preprocessed;
    ASSERT_WP2_OK(PreprocessNearLossless(original, encoder_config,
                                         /*is_alpha=*/false, &preprocessed));
    // Check we have fewer colors and the encoded image is smaller than at
    // higher quality levels. While this is not strictly guaranteed by the
    // algorithm, it's generally the case for typical images and is a desirable
    // property (kind of the whole point).
    const uint32_t num_colors = CountColors(preprocessed);
    EXPECT_LT(num_colors, previous_num_colors);
    previous_num_colors = num_colors;

    memory_writer.Reset();
    ASSERT_WP2_OK(Encode(preprocessed, &memory_writer, encoder_config));
    const uint32_t size = memory_writer.size_;
    EXPECT_LT(size, previous_size);
    previous_size = size;
  }
}

TEST(PreprocessNearLosslessTest, Mask) {
  ArgbBuffer original;
  ASSERT_WP2_OK(
      ReadImage(testing::GetTestDataPath("logo.png").c_str(), &original));
  // Convert image to a grayscale mask.
  ArgbBuffer mask;
  ASSERT_WP2_OK(mask.Resize(original.width, original.height));
  for (uint32_t y = 0; y < original.height; ++y) {
    const uint8_t* const src_row = (uint8_t*)original.GetRow(y);
    uint8_t* const dst = (uint8_t*)mask.GetRow(y);
    for (uint32_t x = 0; x < original.width; ++x) {
      dst[4 * x + 0] = kAlphaMax;
      dst[4 * x + 1] = dst[4 * x + 2] = dst[4 * x + 3] = src_row[4 * x + 0];
    }
  }

  // Encode losslessly.
  EncoderConfig encoder_config;
  encoder_config.alpha_quality = 100;  // lossless
  MemoryWriter memory_writer;
  ASSERT_WP2_OK(Encode(mask, &memory_writer, encoder_config));
  uint32_t previous_size = memory_writer.size_;

  uint32_t previous_num_colors = CountColors(mask);

  // Encode at various near lossless levels.
  for (uint32_t q : {99, 98, 97, 96}) {
    SCOPED_TRACE(SPrintf("Quality: %d", q));
    encoder_config.alpha_quality = q;
    ArgbBuffer preprocessed;
    ASSERT_WP2_OK(PreprocessNearLossless(mask, encoder_config,
                                         /*is_alpha=*/true, &preprocessed));
    // Check we have fewer colors and the encoded image is smaller than at
    // higher quality levels. While this is not strictly guaranteed by the
    // algorithm, it's generally the case for typical images and is a desirable
    // property (kind of the whole point).
    const uint32_t num_colors = CountColors(preprocessed);
    // Check there are fewer colors (or equal, a bit less strict here as there
    // are not many colors to begin with).
    EXPECT_LT(num_colors, previous_num_colors);
    previous_num_colors = num_colors;

    memory_writer.Reset();
    ASSERT_WP2_OK(Encode(preprocessed, &memory_writer, encoder_config));
    const uint32_t size = memory_writer.size_;
    // Would be nice if it was always smaller, but that's not the case... Allow
    // a 2% increase :(
    EXPECT_LT(size, previous_size * 1.02);
    previous_size = size;

    // Check that fully white and fully black pixels are unchanged.
    for (uint32_t y = 0; y < original.height; ++y) {
      const uint8_t* const mask_row = (uint8_t*)mask.GetRow(y);
      const uint8_t* const nll_row = (uint8_t*)preprocessed.GetRow(y);
      for (uint32_t x = 0; x < original.width; ++x) {
        if (mask_row[4 * x + 1] == 0 || mask_row[4 * x + 1] == 255) {
          EXPECT_EQ(nll_row[4 * x + 1], mask_row[4 * x + 1]);
        }
      }
    }
  }
}

}  // namespace
}  // namespace WP2
