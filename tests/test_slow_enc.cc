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

// Maximum compression effort test.

#include <string>

#include "include/helpers.h"
#include "src/wp2/decode.h"
#include "src/wp2/encode.h"

namespace WP2 {
namespace {

//------------------------------------------------------------------------------

class EncodeDecodeTest : public ::testing::TestWithParam<
                             std::tuple<std::string, float, int, float>> {};

TEST_P(EncodeDecodeTest, Simple) {
  const std::string& file_name = std::get<0>(GetParam());
  EncoderConfig config = EncoderConfig::kDefault;
  config.quality = std::get<1>(GetParam());
  config.speed = std::get<2>(GetParam());
  config.sns = std::get<3>(GetParam());

  ASSERT_WP2_OK(testing::EncodeDecodeCompare(file_name, config));
}

INSTANTIATE_TEST_SUITE_P(
    EncodeDecodeTestInstantiationLosslessMaxSpeed, EncodeDecodeTest,
    ::testing::Combine(::testing::Values("source1_64x48.png"),
                       /*quality=*/::testing::Values(100.f),
                       /*speed=*/::testing::Values(9),
                       /*sns=*/::testing::Values(30.f)));

INSTANTIATE_TEST_SUITE_P(
    EncodeDecodeTestInstantiationLossyMaxSpeed, EncodeDecodeTest,
    ::testing::Combine(::testing::Values("source1_64x48.png"),
                       /*quality=*/::testing::Values(75.f),
                       /*speed=*/::testing::Values(9),
                       ::testing::Values(EncoderConfig::kDefault.sns)));

// This one takes a while to run so it is disabled.
// Can still be run with flag --test_arg=--gunit_also_run_disabled_tests
INSTANTIATE_TEST_SUITE_P(
    DISABLED_EncodeDecodeTestInstantiation, EncodeDecodeTest,
    ::testing::Combine(::testing::Values("alpha_ramp.png", "source3.jpg"),
                       /*quality=*/::testing::Values(0.f, 50.f, 98.f),
                       /*speed=*/::testing::Values(6, 7, 8, 9),
                       ::testing::Values(EncoderConfig::kDefault.sns)));

//------------------------------------------------------------------------------

}  // namespace
}  // namespace WP2
