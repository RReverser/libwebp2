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

// Integration tests for tile shapes.

#include <string>

#include "include/helpers.h"
#include "src/wp2/decode.h"
#include "src/wp2/encode.h"

namespace WP2 {
namespace {

//------------------------------------------------------------------------------

class TileShapeTest
    : public ::testing::TestWithParam<std::tuple<std::string, TileShape>> {};

TEST_P(TileShapeTest, Simple) {
  const std::string& file_name = std::get<0>(GetParam());
  EncoderConfig config = EncoderConfig::kDefault;
  config.tile_shape = (TileShape)std::get<1>(GetParam());

  ASSERT_WP2_OK(testing::EncodeDecodeCompare(file_name, config));
}

INSTANTIATE_TEST_SUITE_P(
    PredictorTestInstantiationAll, TileShapeTest,
    ::testing::Combine(::testing::Values("source1_64x48.png", "background.png",
                                         "large.png", "verylarge.png"),
                       ::testing::Values(TILE_SHAPE_SQUARE_128,
                                         TILE_SHAPE_SQUARE_256,
                                         TILE_SHAPE_SQUARE_512,
                                         TILE_SHAPE_WIDE)));

//------------------------------------------------------------------------------

}  // namespace
}  // namespace WP2
