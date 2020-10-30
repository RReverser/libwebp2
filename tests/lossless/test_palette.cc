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

// Test Palette class.

#include "src/dec/lossless/losslessi_dec.h"
#include "src/enc/lossless/losslessi_enc.h"
#include "src/enc/lossless/palette.h"
#include "src/wp2/decode.h"
#include "src/wp2/encode.h"
#include "tests/include/helpers.h"

namespace WP2 {
namespace {

// num_colors, seed
class PaletteTest
    : public ::testing::TestWithParam<std::tuple<uint32_t, uint32_t>> {};

TEST_P(PaletteTest, RoundTrip) {
  WP2LDspInit();
  const uint32_t num_colors = std::get<0>(GetParam());
  const uint32_t seed = std::get<1>(GetParam());

  UniformIntDistribution random(seed);

  // Generate random colors. Note we don't check for collisions (we might
  // randomly create the same color twice) but the odds are ridiculously small.
  Vector_u8 colors;
  ASSERT_TRUE(colors.resize(num_colors * 4));
  for (uint32_t i = 0; i < num_colors; ++i) {
    const uint32_t a = random.Get(0u, kAlphaMax);
    colors[4 * i + 0] = a;
    for (uint32_t c = 1; c < 4; ++c) colors[4 * i + c] = random.Get(0u, a);
  }

  // Generate random source buffer.
  const uint32_t kWidth = 100;
  const uint32_t kHeight = num_colors;
  ArgbBuffer buffer;
  ASSERT_WP2_OK(buffer.Resize(kWidth, kHeight));
  for (uint32_t y = 0; y < kHeight; ++y) {
    uint8_t* const row = (uint8_t*)buffer.GetRow(y);
    for (uint32_t x = 0; x < kWidth; ++x) {
      const uint32_t color_index = random.Get(0u, num_colors - 1);
      std::copy(&colors[4 * color_index], &colors[4 * color_index] + 4,
                &row[4 * x]);
    }
  }

  // Create palette.
  WP2L::Palette palette;
  const bool has_alpha = buffer.HasTransparency();
  palette.Init(buffer.format, has_alpha);
  ASSERT_WP2_OK(palette.AnalyzeAndCreate(buffer));
  // There could be fewer colors if we never chose some of them but it's pretty
  // unlikely so check for equality.
  EXPECT_EQ(palette.Size(), num_colors);

  // Write palette.
  ANSEnc enc;
  WP2::ANSDictionaries dicts;
  const WP2::EncoderConfig encoder_config;
  WP2L::Encoder encoder(encoder_config, buffer, has_alpha);
  ASSERT_WP2_OK(encoder.Allocate());
  ASSERT_WP2_OK(palette.FindBestMethod(&encoder, /*speed=*/5));
  ASSERT_WP2_OK(palette.Write(&enc, &dicts, &encoder));

  // Apply palette on buffer.
  Vector_u16 paletted_buffer;
  ASSERT_TRUE(paletted_buffer.resize(kWidth * kHeight * 4));
  ASSERT_WP2_OK(palette.Apply(buffer, paletted_buffer.data()));

  // Decode palette.
  WP2L::Decoder decoder;
  const WP2::DecoderConfig decoder_config;
  GlobalParams gparams;
  gparams.has_alpha_ = has_alpha;
  ASSERT_WP2_OK(enc.Assemble());
  ExternalDataSource data_source(enc.Buffer(), enc.BufferSize());
  ANSDec dec(&data_source);
  Tile tile;
  decoder.Init(decoder_config, gparams,
               /*progress=*/nullptr, &dec, &tile);
  WP2L::Transform transform;
  transform.type_ = WP2L::COLOR_INDEXING_TRANSFORM;
  transform.width_ = kWidth;
  transform.height_ = kHeight;
  ASSERT_WP2_OK(decoder.ReadPalette(&transform));
  EXPECT_EQ(transform.data_.size(), colors.size());
  for (uint32_t i = 0; i < 4 * num_colors; ++i) {
    SCOPED_TRACE(i);
    EXPECT_EQ(transform.data_[i], palette.GetColors()[i]);
  }

  // Reconstruct original buffer.
  Vector_u16 reconstructed;
  ASSERT_TRUE(reconstructed.resize(kWidth * kHeight * 4));
  WP2LInverseTransform(&transform, /*row_start=*/0, /*row_end=*/kHeight,
                       /*channel_bits=*/8, paletted_buffer.data(),
                       /*has_alpha=*/true, reconstructed.data());

  for (uint32_t y = 0; y < kHeight; ++y) {
    uint8_t* const reference_row = (uint8_t*)buffer.GetRow(y);
    uint16_t* const reconstructed_row = &reconstructed[y * kWidth * 4];
    for (uint32_t x = 0; x < kWidth; ++x) {
      for (uint32_t c = 0; c < 4; ++c) {
        ASSERT_EQ(reconstructed_row[4 * x + c],
                  (uint16_t)reference_row[4 * x + c]);
      }
    }
  }
}

INSTANTIATE_TEST_SUITE_P(PaletteTestInstantiation, PaletteTest,
                         ::testing::Combine(
                             // num_colors
                             ::testing::Values<uint32_t>(2, 5, 10, 30, 50, 100),
                             // seed
                             ::testing::Range<uint32_t>(0, 5)));

}  // namespace
}  // namespace WP2
