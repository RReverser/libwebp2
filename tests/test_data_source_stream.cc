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
// -----------------------------------------------------------------------------

// DataSource test.

#include <array>
#include <limits>
#include <numeric>

#include "include/helpers.h"
#include "src/utils/data_source_stream.h"

namespace WP2 {
namespace {

//------------------------------------------------------------------------------

TEST(StreamDataSourceTest, Simple) {
  std::array<uint8_t, 256> bytes = {{0}};
  std::iota(std::begin(bytes), std::end(bytes), 0);

  StreamDataSource data_source;
  const uint8_t* data;

  ASSERT_EQ(data_source.GetNumNextBytes(), 0u);
  ASSERT_FALSE(data_source.TryGetNext(1, &data));

  ASSERT_WP2_OK(data_source.AppendAsExternal(&bytes.at(0), 1));
  ASSERT_EQ(data_source.GetNumNextBytes(), 1u);
  ASSERT_TRUE(data_source.TryReadNext(1, &data));
  ASSERT_EQ(data[0], bytes[0]);

  ASSERT_WP2_OK(data_source.AppendAsExternal(&bytes.at(1), 10));
  ASSERT_WP2_OK(data_source.AppendExternalToInternal());
  ASSERT_EQ(data_source.GetNumNextBytes(), 10u);
  ASSERT_TRUE(data_source.TryReadNext(10, &data));
  ASSERT_EQ(data[0], bytes[1]);

  data_source.Discard(data_source.GetNumReadBytes());

  ASSERT_WP2_OK(data_source.AppendAsExternal(&bytes.at(11), 1));
  ASSERT_WP2_OK(data_source.AppendAsExternal(&bytes.at(12), 1));
  ASSERT_WP2_OK(data_source.AppendAsExternal(&bytes.at(13), 1));
  ASSERT_WP2_OK(data_source.AppendAsExternal(&bytes.at(14), 1));
  ASSERT_EQ(data_source.GetNumNextBytes(), 4u);
  ASSERT_TRUE(data_source.TryReadNext(4, &data));
  ASSERT_EQ(data[0], bytes[11]);

  data_source.UnmarkAllReadBytes();

  ASSERT_WP2_OK(data_source.AppendAsExternal(&bytes.at(11), 196));
  ASSERT_TRUE(data_source.TryGetNext(200, &data));
  ASSERT_EQ(data[0], bytes[11]);
}

//------------------------------------------------------------------------------

TEST(StreamDataSourceTest, DiscardOverflow) {
  std::array<uint8_t, 30> bytes = {{0}};
  std::iota(std::begin(bytes), std::end(bytes), 0);

  StreamDataSource data_source;
  const uint8_t* data;

  // Append 10
  ASSERT_WP2_OK(data_source.AppendAsExternal(&bytes.at(0), 10));
  ASSERT_TRUE(data_source.TryReadNext(10, &data));

  // Discard 20 (more than previously added)
  data_source.Discard(20);
  ASSERT_FALSE(data_source.TryReadNext(10, &data));

  // Append 10 (the 10 left to discard should be disposed of after that)
  ASSERT_WP2_OK(data_source.AppendAsExternal(&bytes.at(10), 10));
  ASSERT_FALSE(data_source.TryReadNext(10, &data));

  // Append 10, verify the values
  ASSERT_WP2_OK(data_source.AppendAsExternal(&bytes.at(20), 10));
  ASSERT_TRUE(data_source.TryReadNext(10, &data));
  for (int i = 0; i < 10; ++i) {
    ASSERT_EQ(data[i], bytes[20 + i]);
  }
  data_source.UnmarkAllReadBytes();

  data_source.Discard(std::numeric_limits<size_t>::max());
  ASSERT_FALSE(data_source.TryReadNext(10, &data));
}

//------------------------------------------------------------------------------

}  // namespace
}  // namespace WP2
