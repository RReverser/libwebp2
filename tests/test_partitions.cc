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

// Partitioning test.

#include <string>

#include "include/helpers.h"
#include "src/enc/partition_score_func.h"
#include "src/wp2/decode.h"
#include "src/wp2/encode.h"
#include "src/wp2/format_constants.h"

namespace WP2 {
namespace {

//------------------------------------------------------------------------------

class PartitionsTest
    : public ::testing::TestWithParam<std::tuple<
          std::string, float, int, PartitionMethod, PartitionSet, bool>> {};

TEST_P(PartitionsTest, Simple) {
  const std::string& file_name = std::get<0>(GetParam());
  EncoderConfig config = EncoderConfig::kDefault;
  config.quality = std::get<1>(GetParam());
  config.speed = std::get<2>(GetParam());
  config.partition_method = std::get<3>(GetParam());
  config.partition_set = std::get<4>(GetParam());
  config.partition_snapping = std::get<5>(GetParam());

  ASSERT_WP2_OK(testing::EncodeDecodeCompare(file_name, config));
}

INSTANTIATE_TEST_SUITE_P(
    PartitionsTestInstantiation, PartitionsTest,
    ::testing::Combine(
        ::testing::Values("source1_64x48.png"),
        /*quality=*/::testing::Values(80),
        ::testing::Values(EncoderConfig::kDefault.speed),
        ::testing::Values(
            MULTIPASS_PARTITIONING, BLOCK_ENCODE_PARTITIONING,
            AREA_ENCODE_PARTITIONING
            // TILE_ENCODE and EXHAUSTIVE_PARTITIONING are too slow.
            ),
        ::testing::Values(ALL_RECTS),
        ::testing::Values(EncoderConfig::kDefault.partition_snapping)));

INSTANTIATE_TEST_SUITE_P(
    PartitionsTestInstantiationMultipass, PartitionsTest,
    ::testing::Combine(
        ::testing::Values("source1_64x48.png"),
        /*quality=*/::testing::Values(50),
        ::testing::Values(EncoderConfig::kDefault.speed),
        ::testing::Values(MULTIPASS_PARTITIONING),
        ::testing::Values(ALL_SQUARES),
        ::testing::Values(EncoderConfig::kDefault.partition_snapping)));

INSTANTIATE_TEST_SUITE_P(
    PartitionsTestInstantiationSlow, PartitionsTest,
    ::testing::Combine(
        ::testing::Values("source1_32x32.png"),
        /*quality=*/::testing::Values(50),
        ::testing::Values(EncoderConfig::kDefault.speed),
        ::testing::Values(TILE_ENCODE_PARTITIONING),
        ::testing::Values(ALL_RECTS),
        ::testing::Values(EncoderConfig::kDefault.partition_snapping)));

INSTANTIATE_TEST_SUITE_P(
    PartitionsTestInstantiationMegaSlow, PartitionsTest,
    ::testing::Combine(
        ::testing::Values("source1_4x4.png"),
        /*quality=*/::testing::Values(35),
        ::testing::Values(EncoderConfig::kDefault.speed),
        ::testing::Values(EXHAUSTIVE_PARTITIONING),
        ::testing::Values(ALL_RECTS),
        ::testing::Values(EncoderConfig::kDefault.partition_snapping)));

INSTANTIATE_TEST_SUITE_P(
    PartitionsTestInstantiationFixed, PartitionsTest,
    ::testing::Combine(
        ::testing::Values("source1_64x48.png"),
        /*quality=*/::testing::Values(30),
        ::testing::Values(EncoderConfig::kDefault.speed),
        ::testing::Range(ALL_4X4_PARTITIONING,
                         NUM_PARTITION_METHODS),
        ::testing::Values(ALL_SQUARES),
        ::testing::Values(EncoderConfig::kDefault.partition_snapping)));

INSTANTIATE_TEST_SUITE_P(
    PartitionsTestInstantiationRDOpt, PartitionsTest,
    ::testing::Combine(::testing::Values("source1_64x48.png"),
                       /*quality=*/::testing::Values(75),
                       ::testing::Values(EncoderConfig::kDefault.speed),
                       ::testing::Values(SPLIT_RECURSE_PARTITIONING,
                                         AREA_ENCODE_PARTITIONING,
                                         SUB_AREA_ENCODE_PARTITIONING),
                       ::testing::Values(ALL_RECTS),
                       /*partition_snapping=*/::testing::Values(true)));

// This one takes a while to run so it is disabled.
// Can still be run with flag --test_arg=--gunit_also_run_disabled_tests
INSTANTIATE_TEST_SUITE_P(
    DISABLED_PartitionsTestInstantiationAll, PartitionsTest,
    ::testing::Combine(
        ::testing::Values("alpha_ramp.png", "source3.jpeg"),
        /*quality=*/::testing::Values(0.f, 25.f, 50.f, 85.f),
        /*speed=*/::testing::Values(0, 3, 6, 9),
        ::testing::Range(MULTIPASS_PARTITIONING, NUM_PARTITION_METHODS),
        ::testing::Values(SMALL_SQUARES, MEDIUM_SQUARES, ALL_SQUARES,
                          SMALL_RECTS, ALL_RECTS, THICK_RECTS),
        /*partition_snapping=*/::testing::Values(true, false)));

//------------------------------------------------------------------------------

}  // namespace
}  // namespace WP2
