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

// Test the intertile filter.

#include <vector>

#include "include/helpers.h"
#include "include/helpers_filter.h"
#include "src/dec/filters/intertile_filter.h"
#include "src/utils/random.h"

namespace WP2 {
namespace {

//------------------------------------------------------------------------------

static constexpr uint32_t kTileSize = 128;

void Filter(bool enabled, uint32_t quality_hint,
            const TilesLayout& tiles_layout, uint32_t step,
            YUVPlane* const canvas) {
  DecoderConfig config = DecoderConfig::kDefault;
  config.enable_deblocking_filter = enabled;
  IntertileFilter intertile_filter;
  BitstreamFeatures features;
  features.quality_hint = quality_hint;
  features.tile_width = features.tile_height = kTileSize;
  GlobalParams gparams;
  ASSERT_WP2_OK(gparams.ApplyDecoderConfig(config));
  intertile_filter.Init(config, features, gparams, tiles_layout, canvas);
  const Tile& last_tile = tiles_layout.tiles.back();
  const uint32_t height = last_tile.rect.y + last_tile.rect.height;

  uint32_t num_deblocked_rows = 0;
  for (uint32_t num_rows = 0; num_rows <= height;
       num_rows += std::max(1u, std::min(step, height - num_rows))) {
    intertile_filter.Deblock(num_rows);
    const uint32_t num_new_deblocked_rows =
        intertile_filter.GetNumFilteredRows();
    ASSERT_LE(num_new_deblocked_rows, num_rows);
    ASSERT_GE(num_new_deblocked_rows, num_deblocked_rows);
    num_deblocked_rows = num_new_deblocked_rows;
  }
  ASSERT_EQ(num_deblocked_rows, height);
}

void CreateCanvasAndFilter(uint32_t width, uint32_t height,
                           int32_t noise_strength, bool enable_filter,
                           uint32_t quality_hint, double res_den,
                           double min_bpp, uint32_t step,
                           YUVPlane* const canvas, YUVPlane* const unfiltered) {
  ArgbBuffer argb_buffer;  // Unused, for API compliance only.
  TilesLayout tiles_layout;
  ASSERT_WP2_OK(
      canvas->Resize(width, height, /*pad=*/kPredWidth, /*has_alpha=*/true));
  ASSERT_WP2_OK(GetTilesLayout(width, height, kTileSize, kTileSize,
                               &argb_buffer, canvas, &tiles_layout));

  std::vector<Vector<Segment>> segments;
  std::vector<YUVPlane> yuvs;
  segments.reserve(tiles_layout.tiles.size());
  yuvs.reserve(tiles_layout.tiles.size());  // Plane16 has no safe copy ctor.

  for (Tile& tile : tiles_layout.tiles) {
    segments.emplace_back(testing::CreateSegments(quality_hint, -128, 127));
    yuvs.emplace_back();
    ASSERT_WP2_OK(
        yuvs.back().Resize(tile.yuv_output.Y.w_, tile.yuv_output.Y.h_));
    for (Channel c : {kYChannel, kUChannel, kVChannel}) {
      yuvs.back().GetChannel(c).Fill(0);
    }

    testing::CreateBlocks(tile.rect.x, tile.rect.y, segments.back(), res_den,
                          min_bpp, /*num_precision_bits=*/8, &yuvs.back(),
                          &tile.block_map);
  }

  CSPTransform csp;
  canvas->Fill(csp.ToYUV(Argb32b{220, 220, 0, 60}));
  for (Channel c : {kYChannel, kUChannel, kVChannel}) {
    testing::Noise(csp.GetYUVMin(), csp.GetYUVMax(), /*seed=*/12 + c,
                   /*strength=*/noise_strength, &canvas->GetChannel(c));
  }
  ASSERT_WP2_OK(unfiltered->Copy(*canvas, /*resize_if_needed=*/true));

  Filter(enable_filter, quality_hint, tiles_layout, step, canvas);
}

//------------------------------------------------------------------------------

class IntertileFilterTest : public ::testing::TestWithParam<uint32_t> {
 public:
  IntertileFilterTest() : step_(GetParam()) {}
  const uint32_t step_;
};

// Tests that no pixel is modified if each plane has one unique color.
TEST_P(IntertileFilterTest, Plain) {
  YUVPlane canvas, unfiltered;
  CreateCanvasAndFilter(/*width=*/3 * kTileSize - 1, /*height=*/kTileSize + 1,
                        /*noise_strength=*/0, /*enable_filter=*/true,
                        /*quality_hint=*/0, /*res_den=*/0., /*min_bpp=*/0.,
                        step_, &canvas, &unfiltered);
  ASSERT_TRUE(testing::AreEqual(canvas, unfiltered));
}

// Tests that no pixel is modified if the filter is disabled.
TEST_P(IntertileFilterTest, DisabledFilter) {
  YUVPlane canvas, unfiltered;
  CreateCanvasAndFilter(/*width=*/3 * kTileSize - 1, /*height=*/kTileSize + 1,
                        /*noise_strength=*/10, /*enable_filter=*/false,
                        /*quality_hint=*/0, /*res_den=*/0., /*min_bpp=*/0.,
                        step_, &canvas, &unfiltered);
  ASSERT_TRUE(testing::AreEqual(canvas, unfiltered));
}

// Tests that no pixel is modified if the quality is too high.
TEST_P(IntertileFilterTest, HighQuality) {
  YUVPlane canvas, unfiltered;
  CreateCanvasAndFilter(/*width=*/3 * kTileSize - 1, /*height=*/kTileSize + 1,
                        /*noise_strength=*/10, /*enable_filter=*/true,
                        /*quality_hint=*/kMaxLossyQualityHint, /*res_den=*/0.,
                        /*min_bpp=*/0., step_, &canvas, &unfiltered);
  ASSERT_TRUE(testing::AreEqual(canvas, unfiltered));
}

// Tests that no pixel is modified if the bpp are too high.
TEST_P(IntertileFilterTest, HighBpp) {
  YUVPlane canvas, unfiltered;
  CreateCanvasAndFilter(/*width=*/3 * kTileSize - 1, /*height=*/kTileSize + 1,
                        /*noise_strength=*/10, /*enable_filter=*/true,
                        /*quality_hint=*/0, /*res_den=*/0., /*min_bpp=*/2.,
                        step_, &canvas, &unfiltered);
  ASSERT_TRUE(testing::AreEqual(canvas, unfiltered));
}

// Tests that no pixel is modified if there is only one tile.
TEST_P(IntertileFilterTest, OneTile) {
  YUVPlane canvas, unfiltered;
  CreateCanvasAndFilter(/*width=*/kTileSize, /*height=*/10,
                        /*noise_strength=*/10, /*enable_filter=*/true,
                        /*quality_hint=*/0, /*res_den=*/0., /*min_bpp=*/0.,
                        step_, &canvas, &unfiltered);
  ASSERT_TRUE(testing::AreEqual(canvas, unfiltered));
}

// Tests that pixels are modified if there are several tiles.
TEST_P(IntertileFilterTest, TileRow) {
  YUVPlane canvas, unfiltered;
  CreateCanvasAndFilter(/*width=*/3 * kTileSize - 1, /*height=*/10,
                        /*noise_strength=*/10, /*enable_filter=*/true,
                        /*quality_hint=*/0, /*res_den=*/0., /*min_bpp=*/0.,
                        step_, &canvas, &unfiltered);
  ASSERT_FALSE(testing::AreEqual(canvas.Y, unfiltered.Y));
  ASSERT_FALSE(testing::AreEqual(canvas.U, unfiltered.U));
  ASSERT_FALSE(testing::AreEqual(canvas.V, unfiltered.V));
}
TEST_P(IntertileFilterTest, TileColumn) {
  YUVPlane canvas, unfiltered;
  CreateCanvasAndFilter(/*width=*/10, /*height=*/2 * kTileSize + 1,
                        /*noise_strength=*/10, /*enable_filter=*/true,
                        /*quality_hint=*/0, /*res_den=*/0., /*min_bpp=*/0.,
                        step_, &canvas, &unfiltered);
  ASSERT_FALSE(testing::AreEqual(canvas.Y, unfiltered.Y));
  ASSERT_FALSE(testing::AreEqual(canvas.U, unfiltered.U));
  ASSERT_FALSE(testing::AreEqual(canvas.V, unfiltered.V));
}

// Tests that pixels are modified in the same way independently of the step.
TEST_P(IntertileFilterTest, Deterministic) {
  std::array<uint32_t, 3> steps{step_, step_ + 1, step_ + 1000};
  std::array<YUVPlane, steps.size()> canvas, unfiltered;

  for (uint32_t i = 0; i < steps.size(); ++i) {
    CreateCanvasAndFilter(/*width=*/3 * kTileSize - 1,
                          /*height=*/kTileSize + 100, /*noise_strength=*/10,
                          /*enable_filter=*/true, /*quality_hint=*/0,
                          /*res_den=*/0., /*min_bpp=*/0.,
                          steps[i], &canvas[i], &unfiltered[i]);
    ASSERT_FALSE(testing::AreEqual(canvas[i].Y, unfiltered[i].Y));
    ASSERT_FALSE(testing::AreEqual(canvas[i].U, unfiltered[i].U));
    ASSERT_FALSE(testing::AreEqual(canvas[i].V, unfiltered[i].V));
    if (i > 0) {
      // The noise seed being the same, the same base image is expected.
      ASSERT_TRUE(testing::AreEqual(unfiltered[0], unfiltered[i]));
      // Same expected filtering too.
      ASSERT_TRUE(testing::AreEqual(canvas[0], canvas[i]));
    }
  }
}

INSTANTIATE_TEST_SUITE_P(IntertileFilterTestInstantiation, IntertileFilterTest,
                         ::testing::Values(1, kTileSize - 1, kTileSize,
                                           kTileSize + 1, kTileSize * 10));

//------------------------------------------------------------------------------

void CreatePlainTiles(uint32_t width, uint32_t height, uint8_t min_color,
                      uint8_t max_color, ArgbBuffer* const canvas) {
  UniformIntDistribution random(/*seed=*/12);
  ASSERT_WP2_OK(canvas->Resize(width, height));
  uint32_t num_tiles_x, num_tiles_y;
  GetNumTiles(canvas->width, canvas->height, kTileSize, kTileSize, &num_tiles_x,
              &num_tiles_y);
  for (uint32_t x = 0; x < num_tiles_x; ++x) {
    for (uint32_t y = 0; y < num_tiles_y; ++y) {
      const Argb32b color{255, random.Get(min_color, max_color),
                          random.Get(min_color, max_color),
                          random.Get(min_color, max_color)};
      const Rectangle window{
          x * kTileSize, y * kTileSize,
          std::min(canvas->width - x * kTileSize, kTileSize),
          std::min(canvas->height - y * kTileSize, kTileSize)};
      canvas->Fill(window, color);
    }
  }
}

// Tests that PSNR improves if noise is added per tile, not per pixel.
// It reduces the impact of the intratile filter.
TEST(IntertileFilterTest, Simple) {
  const uint32_t min = 150, max = 220, avg = (min + max) / 2;
  ArgbBuffer plain_planes, plain_tiles;
  CreatePlainTiles(2 * kTileSize + 1, 2 * kTileSize - 1, avg, avg,
                   &plain_planes);
  CreatePlainTiles(plain_planes.width, plain_planes.height, min, max,
                   &plain_tiles);

  EncoderConfig encoder_config = EncoderConfig::kDefault;
  encoder_config.quality = 0;  // Maximum filtering.

  DecoderConfig decoder_config = DecoderConfig::kDefault;
  decoder_config.enable_directional_filter = false;
  decoder_config.enable_restoration_filter = false;

  float psnr_filter_on[5], psnr_filter_off[5];
  for (bool filter : {true, false}) {
    decoder_config.enable_deblocking_filter = filter;

    MemoryWriter memory_writer;
    ASSERT_WP2_OK(Encode(plain_tiles, &memory_writer, encoder_config));

    ArgbBuffer decompressed;
    ASSERT_WP2_OK(Decode(memory_writer.mem_, memory_writer.size_, &decompressed,
                         decoder_config));

    ASSERT_WP2_OK(decompressed.GetDistortion(plain_planes, PSNR,
                                             filter ? psnr_filter_on
                                                    : psnr_filter_off));
  }
  // Compared to the plain image, the blurred tiles should be closer than the
  // intact plain tiles.
  EXPECT_GT(psnr_filter_on[4], psnr_filter_off[4]);
}

//------------------------------------------------------------------------------

}  // namespace
}  // namespace WP2
