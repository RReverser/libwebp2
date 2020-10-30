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

// Simple encoding test.

#include "src/wp2/encode.h"

#include "imageio/image_dec.h"
#include "include/helpers.h"

namespace WP2 {
namespace {

TEST(EncodeTest, Simple) {
  ArgbBuffer src;
  ASSERT_WP2_OK(
      ReadImage(testing::GetTestDataPath("source1.png").c_str(), &src));

  MemoryWriter memory_writer;
  EXPECT_WP2_OK(Encode(src, &memory_writer));
}

}  // namespace
}  // namespace WP2
