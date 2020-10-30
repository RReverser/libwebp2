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

// Test Group4 lossless algorithm for binary images.

#include "src/common/color_precision.h"
#include "src/wp2/decode.h"
#include "src/wp2/encode.h"
#include "tests/include/helpers.h"

namespace WP2 {
namespace {

static std::string DisplayMask(int dim, const ArgbBuffer &image) {
  std::string out = "****\n";
  std::vector<uint32_t> colors;
  for (int y = 0; y < dim; ++y) {
    const auto row = (const uint32_t *)image.GetRow(y);
    for (int x = 0; x < dim; ++x) {
      uint32_t i = 0;
      while (i < colors.size() && colors[i] != row[x]) ++i;
      if (i == colors.size()) colors.push_back(row[x]);
      out.append(std::to_string(i));
    }
    out.append("\n");
  }
  out.append("\n");
  return out;
}

//------------------------------------------------------------------------------

enum class Type { kCircle, kRandom };

class Group4Test : public ::testing::TestWithParam<std::tuple<Type, int, int>> {
};

static void RandomColor(UniformIntDistribution *const gen, uint8_t color[4]) {
  for (uint32_t c = 0; c < 4; ++c) {
    if (c == 0 && gen->FlipACoin()) {
      color[c] = WP2::kAlphaMax;
    } else {
      color[c] = gen->Get(0u, (c == 0) ? WP2::kAlphaMax : color[0]);
    }
  }
}

TEST_P(Group4Test, Pattern) {
  const Type type = std::get<0>(GetParam());
  const uint32_t seed = std::get<1>(GetParam());
  const int dim = std::get<2>(GetParam());
  WP2::UniformIntDistribution gen(seed);

  ArgbBuffer original;
  ASSERT_EQ(original.Resize(dim, dim), WP2_STATUS_OK);

  // Draw a pattern.
  if (type == Type::kRandom) {
    // Choose two colors.
    uint8_t colors[2][4];
    for (uint32_t i = 0; i < 2; ++i) {
      RandomColor(&gen, colors[i]);
    }

    for (int y = 0; y < dim; ++y) {
      uint8_t *const row = (uint8_t *)original.GetRow(y);
      for (int x = 0; x < dim; ++x) {
        const uint32_t color_index = gen.FlipACoin();
        std::copy(colors[color_index], colors[color_index] + 4, &row[4 * x]);
      }
    }
  } else {
    uint8_t color[4];
    RandomColor(&gen, color);
    original.Fill(WP2::ToArgb32b(color));

    const uint32_t num_circles = gen.Get(1, dim / 10 + 1);
    for (uint32_t i = 1; i < num_circles; ++i) {
      const int center_x = gen.Get(0, dim - 1);
      const int center_y = gen.Get(0, dim - 1);
      const int radius = gen.Get(0, dim);

      RandomColor(&gen, color);

      for (int y = 0; y < dim; ++y) {
        uint8_t *const row = (uint8_t *)original.GetRow(y);
        for (int x = 0; x < dim; ++x) {
          if ((x - center_x) * (x - center_x) +
                  (y - center_y) * (y - center_y) <
              radius * radius) {
            std::copy(color, color + 4, &row[4 * x]);
          }
        }
      }
    }
  }

  MemoryWriter memory_writer;
  EncoderConfig encoder_config;
  encoder_config.quality = 100;
  WP2_ASSERT_STATUS(Encode(original, &memory_writer, encoder_config));

  DecoderConfig decoder_config;
  ArgbBuffer decompressed;
  WP2_ASSERT_STATUS(Decode(memory_writer.mem_, memory_writer.size_,
                           &decompressed, decoder_config));

  EXPECT_TRUE(testing::Compare(original, decompressed, /*file_name=*/"pattern"))
      << DisplayMask(dim, original) << DisplayMask(dim, decompressed);
}  // namespace

INSTANTIATE_TEST_SUITE_P(
    Group4Test, Group4Test,
    ::testing::Combine(::testing::Values(Type::kCircle, Type::kRandom),
                       ::testing::Range(0, 100),
                       ::testing::Values(WP2L::kGroup4Window, 101)));

}  // namespace
}  // namespace WP2
