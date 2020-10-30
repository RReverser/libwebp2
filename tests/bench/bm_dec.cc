// Copyright 2019 Google LLC
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
// -----------------------------------------------------------------------------
//
// Measures WP2 decoding time.

#include <cstddef>

#include "imageio/image_dec.h"
#include "src/wp2/decode.h"
#include "src/wp2/encode.h"
#include "tests/bench/helpers_bm.h"
#include "tests/include/helpers.h"

namespace WP2 {
namespace {

void BM_DecodeImage(benchmark::State& state) {
  const char* const file_name = testing::kFiles[state.range(0)];
  ArgbBuffer original;
  ASSERT_EQ(ReadImage(testing::GetTestDataPath(file_name).c_str(), &original),
            WP2_STATUS_OK);

  MemoryWriter data;
  ASSERT_EQ(Encode(original, &data), WP2_STATUS_OK);

  for (const auto& _ : state) {
    (void)_;  // The content of this loop is benchmarked.
    ArgbBuffer output;
    EXPECT_EQ(Decode(data.mem_, data.size_, &output), WP2_STATUS_OK);
  }
}
BENCHMARK(BM_DecodeImage)->DenseRange(0, testing::kNumFiles - 1);

}  // namespace
}  // namespace WP2
