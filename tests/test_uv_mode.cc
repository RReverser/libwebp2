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
#include "src/wp2/decode.h"
#include "src/wp2/encode.h"

namespace WP2 {
namespace {

//------------------------------------------------------------------------------

class UVModeTest : public ::testing::TestWithParam<
                       std::tuple<std::string, EncoderConfig::UVMode>> {};

TEST_P(UVModeTest, Simple) {
  const std::string& file_name = std::get<0>(GetParam());
  EncoderConfig config = EncoderConfig::kDefault;
  config.uv_mode = std::get<1>(GetParam());

  ASSERT_WP2_OK(testing::EncodeDecodeCompare(file_name, config));
}

INSTANTIATE_TEST_SUITE_P(
    UVModeTestInstantiation, UVModeTest,
    ::testing::Combine(::testing::Values("source1_64x48.png",
                                         "source3_222x167.jpg",
                                         "source1_1x1.png"),
                       ::testing::Range(EncoderConfig::UVModeAdapt,
                                        EncoderConfig::NumUVMode)));

INSTANTIATE_TEST_SUITE_P(
    UVModeTestInstantiationSlow, UVModeTest,
    ::testing::Combine(::testing::Values("test_exif_xmp.webp"),
                       ::testing::Values(EncoderConfig::UVMode420)));

//------------------------------------------------------------------------------

}  // namespace
}  // namespace WP2
