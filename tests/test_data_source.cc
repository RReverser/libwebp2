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
#include "src/utils/data_source_context.h"
#include "src/utils/data_source.h"
#include "src/utils/random.h"

namespace WP2 {
namespace {

//------------------------------------------------------------------------------

TEST(DataSourceTest, Simple) {
  std::array<uint8_t, 256> bytes = {{0}};
  std::iota(std::begin(bytes), std::end(bytes), 0);

  ExternalDataSource data_source(&bytes.at(0), bytes.size());
  const uint8_t* data;

  ASSERT_EQ(data_source.GetNumNextBytes(), bytes.size());
  ASSERT_TRUE(data_source.TryGetNext(0, &data));
  ASSERT_TRUE(data_source.TryGetNext(1, &data));
  ASSERT_EQ(data[0], bytes[0]);
  ASSERT_TRUE(data_source.TryGetNext(bytes.size(), &data));
  ASSERT_EQ(data[0], bytes[0]);
  ASSERT_FALSE(data_source.TryGetNext(bytes.size() + 1, &data));
  ASSERT_EQ(data_source.GetNumNextBytes(), bytes.size());
  data_source.MarkNumBytesAsRead(bytes.size() - 1);
  ASSERT_EQ(data_source.GetNumNextBytes(), 1u);
  data_source.UnmarkAllReadBytes();

  ASSERT_TRUE(data_source.TryReadNext(bytes.size() - 1, &data));
  ASSERT_EQ(data_source.GetNumNextBytes(), 1u);
  data_source.Discard(data_source.GetNumReadBytes());
  data_source.UnmarkAllReadBytes();
  ASSERT_TRUE(data_source.TryReadNext(1, &data));
  ASSERT_EQ(data[0], bytes.back());
  ASSERT_EQ(data_source.GetNumNextBytes(), 0u);
  data_source.UnmarkAllReadBytes();
  data_source.Discard(0);
  data_source.Discard(5);
  ASSERT_EQ(data_source.GetNumNextBytes(), 0u);
  ASSERT_FALSE(data_source.TryReadNext(1, &data));
}

//------------------------------------------------------------------------------

TEST(DataSourceTest, LoopRead) {
  std::array<uint8_t, 1024> bytes = {{0}};
  std::iota(std::begin(bytes), std::end(bytes), 0);

  UniformIntDistribution random(/*seed=*/bytes.size());  // Deterministic.

  ExternalDataSource data_source(&bytes.at(0), bytes.size());

  size_t size_left = bytes.size();
  while (size_left > 0) {
    const size_t step = std::min(size_left, random.Get<size_t>(0u, 10u));
    const uint8_t* data;
    ASSERT_TRUE(data_source.TryReadNext(step, &data));
    if (step > 0) {
      ASSERT_EQ(data[0], bytes.at(bytes.size() - size_left));
    }
    size_left -= step;
  }
  ASSERT_EQ(data_source.GetNumNextBytes(), 0u);
}

//------------------------------------------------------------------------------

TEST(DataSourceTest, LoopDiscard) {
  std::array<uint8_t, 512> bytes = {{0}};
  std::iota(std::begin(bytes), std::end(bytes), 0);

  UniformIntDistribution random(/*seed=*/bytes.size());  // Deterministic.

  ExternalDataSource data_source(&bytes.at(0), bytes.size());

  size_t size_left = bytes.size();
  while (size_left > 0) {
    const size_t step = std::min(size_left, random.Get<size_t>(0u, 10u));
    const uint8_t* data;
    ASSERT_TRUE(data_source.TryGetNext(step, &data));
    if (step > 0) {
      ASSERT_EQ(data[0], bytes.at(bytes.size() - size_left));
    }
    data_source.Discard(step);
    size_left -= step;
  }
  ASSERT_EQ(data_source.GetNumNextBytes(), 0u);
}

//------------------------------------------------------------------------------

TEST(DataSourceTest, LoopUpdate) {
  std::array<uint8_t, 512> bytes = {{0}};
  std::iota(std::begin(bytes), std::end(bytes), 0);

  UniformIntDistribution random(/*seed=*/bytes.size());  // Deterministic.

  ExternalDataSource data_source;

  size_t total_read_bytes = 0;
  size_t total_discarded_bytes = 0;
  while (total_read_bytes < bytes.size() &&
         total_discarded_bytes < bytes.size()) {
    const size_t size_left =
        bytes.size() - std::min(bytes.size(), std::max(total_read_bytes,
                                                       total_discarded_bytes));
    const size_t step = std::min(size_left, random.Get<size_t>(0u, 10u));
    const uint8_t* data;
    data_source.Update(&bytes.at(0), bytes.size() - size_left + step);
    ASSERT_TRUE(data_source.TryGetNext(step, &data));
    if (step > 0) {
      ASSERT_EQ(data[0], bytes.at(bytes.size() - size_left));
    }
    const size_t num_read_bytes = random.Get<size_t>(0u, 10u);
    const size_t num_discarded_bytes = random.Get<size_t>(0u, 10u);
    data_source.MarkNumBytesAsRead(num_read_bytes);
    data_source.Discard(num_discarded_bytes);
    total_discarded_bytes += num_discarded_bytes;
    total_read_bytes =
        std::max(total_read_bytes + num_read_bytes, total_discarded_bytes);
  }
  ASSERT_EQ(data_source.GetNumNextBytes(), 0u);
}

//------------------------------------------------------------------------------

TEST(DataSourceTest, Update) {
  std::array<uint8_t, 128> bytes = {{0}};
  ExternalDataSource data_source(&bytes.at(0), 1);
  const uint8_t* data;

  ASSERT_TRUE(data_source.TryReadNext(1, &data));
  data_source.Discard(data_source.GetNumReadBytes());
  data_source.Update(&bytes.at(0), bytes.size());
  ASSERT_TRUE(data_source.TryReadNext(bytes.size() - 1, &data));
  ASSERT_EQ(data[0], bytes[1]);
}

//------------------------------------------------------------------------------

class MyDataSource : public DataSource {
 public:
  MyDataSource() {
    std::iota(std::begin(bytes_), std::end(bytes_), 0);  // Fill with 0, 1, 2...
    available_bytes_ = bytes_.data();
    num_available_bytes_ = 0;
  }

 private:
  bool Fetch(size_t num_requested_bytes) override {
    if ((num_discarded_bytes_ + num_read_bytes_ + num_requested_bytes) >
        bytes_.size()) {
      return false;  // Can't fetch enough data.
    }
    available_bytes_ = &bytes_[num_discarded_bytes_];
    num_available_bytes_ = num_read_bytes_ + num_requested_bytes;
    return true;
  }
  void OnDiscard(size_t num_bytes) override {}

 public:
  std::array<uint8_t, 128> bytes_;
};

TEST(DataSourceTest, Fetch) {
  MyDataSource data_source;
  const uint8_t* data;

  const size_t kStep1 = 27;
  const size_t kStep2 = 3;
  const size_t kStep3 = data_source.bytes_.size() - kStep1 - kStep2;
  data_source.Discard(0);
  ASSERT_TRUE(data_source.TryReadNext(kStep1, &data));
  ASSERT_EQ(data[0], data_source.bytes_[0]);

  ASSERT_TRUE(data_source.TryReadNext(kStep2 - 1u, &data));
  ASSERT_EQ(data[0], data_source.bytes_[kStep1]);
  data_source.Discard(data_source.GetNumReadBytes() + 1u);

  ASSERT_TRUE(data_source.TryReadNext(kStep3, &data));
  ASSERT_EQ(data[0], data_source.bytes_[kStep1 + kStep2]);

  ASSERT_FALSE(data_source.TryReadNext(1, &data));
  ASSERT_EQ(data_source.GetNumNextBytes(), 0u);
  data_source.Discard(std::numeric_limits<uint32_t>::max());
  ASSERT_EQ(data_source.GetNumNextBytes(), 0u);
  ASSERT_FALSE(data_source.TryReadNext(1, &data));
}

//------------------------------------------------------------------------------

#if defined(WP2_USE_CONTEXT_SWITCH) && (WP2_USE_CONTEXT_SWITCH > 0)

TEST(DataSourceTest, SuspendableDataSourceCoverage) {
  SuspendableDataSource data_source;
  data_source.SetContext(nullptr);
  data_source.Update(nullptr, 0);
  const uint8_t* data;
  ASSERT_FALSE(data_source.TryReadNext(1, &data));
  data_source.Reset();
}

#endif  // WP2_USE_CONTEXT_SWITCH

//------------------------------------------------------------------------------

}  // namespace
}  // namespace WP2
