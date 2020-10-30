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

// Test WP2Malloc() failure at different places.

#include <cstring>
#include <iostream>
#include <vector>

#include "include/helpers.h"
#include "src/utils/vector.h"
#include "src/wp2/decode.h"

//------------------------------------------------------------------------------

extern void WP2SetMallocFailAt(int num_malloc_before_failure);

namespace WP2 {
namespace {

// In order to generate only WP2_STATUS_OUT_OF_MEMORY failures, MemoryWriter is
// replaced by this class.
class NoAllocationFailureMemoryWriter : public Writer {
 public:
  void Reset() { buffer_.clear(); }
  bool Append(const void* data, size_t data_size) override {
    const uint8_t* const bytes = reinterpret_cast<const uint8_t*>(data);
    buffer_.insert(buffer_.end(), bytes, bytes + data_size);
    return true;
  }
  const uint8_t* Data() const {
    return reinterpret_cast<const uint8_t*>(buffer_.data());
  }
  size_t Size() const { return buffer_.size(); }

 private:
  std::vector<uint8_t> buffer_;
};

class MallocFailureTest
    : public ::testing::TestWithParam<
          std::tuple<std::string, float, int, uint32_t, int, int>> {};

//------------------------------------------------------------------------------

TEST_P(MallocFailureTest, Encoding) {
  const std::string& src_file_name = std::get<0>(GetParam());
  const float quality = std::get<1>(GetParam());
  const int speed = std::get<2>(GetParam());
  const int thread_level = (std::get<3>(GetParam()) == 0) ? 0 : 1;
  const int min_malloc_fail_at = std::get<4>(GetParam());
  const int max_malloc_fail_at = std::get<5>(GetParam());

  ArgbBuffer src;
  NoAllocationFailureMemoryWriter data;
  int malloc_fail_at = min_malloc_fail_at;
  WP2Status status = WP2_STATUS_INVALID_PARAMETER;
  while (malloc_fail_at <= max_malloc_fail_at) {
    src.Deallocate();
    data.Reset();
    WP2SetMallocFailAt(malloc_fail_at);
    status = testing::CompressImage(src_file_name, &data, &src, quality, speed,
                                    thread_level);
    if (status != WP2_STATUS_OK) {
      ASSERT_EQ(status, WP2_STATUS_OUT_OF_MEMORY)
          << "Allocation failure number " << malloc_fail_at << " is uncaught";
    } else {
      break;
    }
    ++malloc_fail_at;
  }
  std::cout << "Compressing " << src_file_name << " (quality " << quality
            << ", speed " << speed << ", mt " << thread_level << ")"
            << " took at least " << malloc_fail_at << " WP2Malloc() calls "
            << "(" << status << ")." << std::endl;

  if (status == WP2_STATUS_OK) {
    // Check a few more times that it was not a false positive.
    for (int i = 1; i <= 5; ++i) {
      WP2SetMallocFailAt(malloc_fail_at + i);
      ArgbBuffer another_src;  // Keep first OK output though.
      NoAllocationFailureMemoryWriter another_data;
      ASSERT_WP2_OK(testing::CompressImage(src_file_name, &another_data,
                                           &another_src, quality, speed,
                                           thread_level))
          << "False positives from " << malloc_fail_at << " to "
          << (malloc_fail_at + i - 1);
    }

    // Decode to see if it matches.
    WP2SetMallocFailAt(0);
    ArgbBuffer output;
    DecoderConfig config;
    config.thread_level = 8;
    ASSERT_WP2_OK(Decode(data.Data(), data.Size(), &output, config));
    ASSERT_TRUE(testing::Compare(src, output, src_file_name,
                                 testing::GetExpectedDistortion(quality)));
  }
}

//------------------------------------------------------------------------------

TEST_P(MallocFailureTest, Decoding) {
  const std::string& src_file_name = std::get<0>(GetParam());
  const float quality = std::get<1>(GetParam());
  const int speed = std::get<2>(GetParam());
  const int thread_level = std::get<3>(GetParam());
  const int min_malloc_fail_at = std::get<4>(GetParam());
  const int max_malloc_fail_at = std::get<5>(GetParam());

  // Compress without allocation failure.
  ArgbBuffer src;
  MemoryWriter data;
  WP2SetMallocFailAt(0);
  ASSERT_WP2_OK(testing::CompressImage(src_file_name, &data, &src, quality,
                                       speed, /*thread_level=*/1));

  ArgbBuffer output;
  DecoderConfig config;
  config.thread_level = (uint32_t)thread_level;
  int malloc_fail_at = min_malloc_fail_at;
  WP2Status status = WP2_STATUS_INVALID_PARAMETER;
  while (malloc_fail_at <= max_malloc_fail_at) {
    output.Deallocate();
    WP2SetMallocFailAt(malloc_fail_at);
    status = Decode(data.mem_, data.size_, &output, config);
    if (status != WP2_STATUS_OK) {
      ASSERT_EQ(status, WP2_STATUS_OUT_OF_MEMORY)
          << "Allocation failure number " << malloc_fail_at << " is uncaught";
    } else {
      break;
    }
    ++malloc_fail_at;
  }
  std::cout << "Decompressing " << src_file_name << " (quality " << quality
            << ", speed " << speed << ", mt " << thread_level << ")"
            << " took at least " << malloc_fail_at << " WP2Malloc() calls "
            << "(" << status << ")." << std::endl;

  if (status == WP2_STATUS_OK) {
    // Check a few more times that it was not a false positive.
    for (int i = 1; i <= 5; ++i) {
      WP2SetMallocFailAt(malloc_fail_at + i);
      ArgbBuffer another_output;  // Keep first OK output though.
      ASSERT_WP2_OK(Decode(data.mem_, data.size_, &another_output, config))
          << "False positives from " << malloc_fail_at << " to "
          << (malloc_fail_at + i - 1);
    }

    WP2SetMallocFailAt(0);
    ASSERT_TRUE(testing::Compare(src, output, src_file_name,
                                 testing::GetExpectedDistortion(quality)));
  }
}

//------------------------------------------------------------------------------

INSTANTIATE_TEST_SUITE_P(
    MallocFailureTestInstantiationLossy, MallocFailureTest,
    ::testing::Combine(::testing::Values("source1_64x48.png"),
                       ::testing::Values(75.f) /* quality */,
                       ::testing::Values(0) /* speed */,
                       ::testing::Values(8) /* multithreading */,
                       ::testing::Values(3) /* min_malloc_fail_at */,
                       ::testing::Values(4) /* max_malloc_fail_at */));

INSTANTIATE_TEST_SUITE_P(
    MallocFailureTestInstantiationLossless, MallocFailureTest,
    ::testing::Combine(::testing::Values("source1_64x48.png"),
                       ::testing::Values(100.f) /* quality */,
                       ::testing::Values(0) /* speed */,
                       ::testing::Values(8) /* multithreading */,
                       ::testing::Values(8) /* min_malloc_fail_at */,
                       ::testing::Values(9) /* max_malloc_fail_at */));

INSTANTIATE_TEST_SUITE_P(
    MallocFailureTestInstantiationSmall, MallocFailureTest,
    ::testing::Combine(::testing::Values("source1_1x1.png"),
                       ::testing::Values(0.f, 99.f) /* quality */,
                       ::testing::Values(5) /* speed */,
                       ::testing::Values(0) /* multithreading */,
                       ::testing::Values(19) /* min_malloc_fail_at */,
                       ::testing::Values(20) /* max_malloc_fail_at */));

// This one might take ages so it is disabled.
// Can still be run with flag --test_arg=--gunit_also_run_disabled_tests
INSTANTIATE_TEST_SUITE_P(
    DISABLED_MallocFailureTestInstantiation, MallocFailureTest,
    ::testing::Combine(::testing::Values("source1.png", "source1_1x1.png",
                                         "source1_1x48.png", "source1_64x1.png",
                                         "source1_64x48.png",
                                         "alpha_ramp.lossy.webp", "source3.jpg",
                                         "test_exif_xmp.webp"),
                       ::testing::Values(0.f, 50.f, 99.f, 100.f) /* quality */,
                       ::testing::Values(0, 5, 9) /* speed */,
                       ::testing::Values(0, 8) /* multithreading */,
                       ::testing::Values(1) /* min_malloc_fail_at */,
                       ::testing::Values(1000000) /* max_malloc_fail_at */));

//------------------------------------------------------------------------------

struct CtorStruct {
  CtorStruct() noexcept : value(0) { ++counter; }
  explicit CtorStruct(int v) : value(v) { ++counter; }
  CtorStruct(const CtorStruct& o) noexcept : value(o.value) { ++counter; }
  CtorStruct(CtorStruct&& o) noexcept : value(o.value) { ++counter; }
  ~CtorStruct() { --counter; }
  int value;
  static int counter;
};
int CtorStruct::counter;

TEST(MallocFailureTest, WP2VectorResize) {
  CtorStruct::counter = -1;
  const CtorStruct element(13);
  {
    Vector<CtorStruct> vector;

    // Successfully put at least one element.
    WP2SetMallocFailAt(2);
    size_t valid_size = 0;
    while (vector.resize(valid_size + 1)) {
      ASSERT_EQ(vector.size(), valid_size + 1);
      valid_size = vector.size();
      vector.back().value = element.value;
      ASSERT_EQ(CtorStruct::counter, (int)valid_size);
    }
    ASSERT_GT(valid_size, 0u);
    ASSERT_EQ(vector.size(), valid_size);

    // Make sure data stored before failure can be retrieved.
    WP2SetMallocFailAt(0);
    for (const CtorStruct& constructed : vector) {
      ASSERT_EQ(constructed.value, element.value);
    }
  }
  ASSERT_EQ(CtorStruct::counter, 0);
}

TEST(MallocFailureTest, WP2VectorPushBack) {
  CtorStruct::counter = -1;
  const CtorStruct element(13);
  {
    Vector<CtorStruct> vector;

    // Successfully put at least one element.
    WP2SetMallocFailAt(2);
    size_t valid_size = 0;
    while (vector.push_back(element, /*resize_if_needed=*/true)) {
      ASSERT_EQ(vector.size(), valid_size + 1);
      valid_size = vector.size();
      ASSERT_EQ(CtorStruct::counter, (int)valid_size);
    }
    ASSERT_GT(valid_size, 0u);
    ASSERT_EQ(vector.size(), valid_size);

    // Make sure data stored before failure can be retrieved.
    WP2SetMallocFailAt(0);
    for (const CtorStruct& constructed : vector) {
      ASSERT_EQ(constructed.value, element.value);
    }
  }
  ASSERT_EQ(CtorStruct::counter, 0);
}

TEST(MallocFailureTest, WP2VectorReservePushBack) {
  CtorStruct::counter = -1;
  const CtorStruct element(13);
  {
    Vector<CtorStruct> vector;

    // Successfully put at least one element.
    WP2SetMallocFailAt(2);
    size_t valid_size = 0;
    while (vector.reserve(valid_size + 1)) {
      ASSERT_EQ(vector.size(), valid_size);
      ASSERT_TRUE(vector.push_back(element, /*resize_if_needed=*/false));
      ASSERT_EQ(vector.size(), valid_size + 1);
      valid_size = vector.size();
      ASSERT_EQ(CtorStruct::counter, (int)valid_size);
    }
    ASSERT_FALSE(vector.push_back(element, /*resize_if_needed=*/false));
    ASSERT_GT(valid_size, 0u);
    ASSERT_EQ(vector.size(), valid_size);
    ASSERT_EQ(vector.capacity(), valid_size);

    // Make sure data stored before failure can be retrieved.
    WP2SetMallocFailAt(0);
    for (const CtorStruct& constructed : vector) {
      ASSERT_EQ(constructed.value, element.value);
    }
  }
  ASSERT_EQ(CtorStruct::counter, 0);
}

//------------------------------------------------------------------------------

}  // namespace
}  // namespace WP2
