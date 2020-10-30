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

// Test example utils.

#include <limits>

#include "imageio/imageio_util.h"
#include "include/helpers.h"

namespace WP2 {
namespace {

//------------------------------------------------------------------------------

TEST(ImageIoUtilTest, Simple) {
  std::string data;
  ASSERT_WP2_OK(IoUtilReadFile(
                testing::GetTestDataPath("source1_64x48.png").c_str(), &data));

  const std::string file_path = testing::GetTempDataPath("test_image_io_util");
  std::remove(file_path.c_str());  // Make sure it doesn't exist.
  ASSERT_WP2_OK(IoUtilWriteFile(data, file_path.c_str(), /*overwrite=*/false));
  ASSERT_EQ(IoUtilWriteFile(data, file_path.c_str(), /*overwrite=*/false),
            WP2_STATUS_BAD_WRITE);
  ASSERT_WP2_OK(IoUtilWriteFile(data, file_path.c_str(), /*overwrite=*/true));

  ASSERT_EQ(IoUtilWriteFile(data, "/tmp", /*overwrite=*/true),
            WP2_STATUS_INVALID_PARAMETER);
  EXPECT_EQ(
      IoUtilReadFile(testing::GetTestDataPath("missing_file").c_str(), &data),
      WP2_STATUS_INVALID_PARAMETER);
}

//------------------------------------------------------------------------------
// Exercize type cast check. It is expected that size_t and uint64_t have the
// same range and are the biggest integer types.

template <typename ResultType, typename ValueType, typename MultType>
bool CheckMultFitsIn(ValueType value, MultType mult,
                     ResultType* const result = nullptr) {
  return (MultFitsIn(value, mult, result) == WP2_STATUS_OK);
}

template <typename T>
void TestMult() {
  T r = 0;

  EXPECT_TRUE(CheckMultFitsIn((T)0, (T)0, &r) && r == 0);
  EXPECT_TRUE(CheckMultFitsIn<T>((T)0, (T)1));
  EXPECT_TRUE(CheckMultFitsIn((T)1, (T)0, &r) && r == 0);
  EXPECT_TRUE(CheckMultFitsIn((T)1, (T)1, &r) && r == 1);
  EXPECT_TRUE(CheckMultFitsIn((int8_t)1, (uint16_t)1, &r) && r == 1);

  static constexpr T kMaxValue = std::numeric_limits<T>::max();
  EXPECT_TRUE(CheckMultFitsIn(kMaxValue, 1, &r) && r == kMaxValue);
  EXPECT_TRUE(CheckMultFitsIn<T>(1, kMaxValue));
  EXPECT_TRUE(CheckMultFitsIn(2, kMaxValue / 2, &r) &&
              r == (kMaxValue / 2) * 2);

  EXPECT_TRUE(CheckMultFitsIn((uint64_t)kMaxValue, (size_t)1, &r) &&
              r == kMaxValue);
  EXPECT_TRUE(CheckMultFitsIn((size_t)kMaxValue, (int8_t)1, &r) &&
              r == kMaxValue);

  EXPECT_FALSE(CheckMultFitsIn(kMaxValue, 2, &r));
  EXPECT_FALSE(CheckMultFitsIn(2, kMaxValue, &r));
  EXPECT_FALSE(CheckMultFitsIn(kMaxValue, kMaxValue, &r));
  EXPECT_FALSE(CheckMultFitsIn((size_t)kMaxValue, (size_t)kMaxValue, &r));
  EXPECT_FALSE(CheckMultFitsIn((size_t)kMaxValue, (uint64_t)kMaxValue, &r));
  EXPECT_FALSE(CheckMultFitsIn((uint64_t)kMaxValue, (uint64_t)kMaxValue, &r));

  if ((uint64_t)kMaxValue < std::numeric_limits<uint64_t>::max()) {
    static constexpr uint64_t kTooBig = (uint64_t)kMaxValue + 1;
    EXPECT_FALSE(CheckMultFitsIn<T>(kTooBig, 1));
    EXPECT_FALSE(CheckMultFitsIn(1, kTooBig, &r));
    EXPECT_FALSE(CheckMultFitsIn(kTooBig, kTooBig, &r));
  }

  static constexpr T kSquareRoot = kMaxValue >> (sizeof(kMaxValue) * 4);
  EXPECT_TRUE(CheckMultFitsIn(kSquareRoot, kSquareRoot, &r) &&
              r == kSquareRoot * kSquareRoot);
  EXPECT_TRUE(CheckMultFitsIn((uint64_t)kSquareRoot, (size_t)kSquareRoot, &r) &&
              r == kSquareRoot * kSquareRoot);
  EXPECT_TRUE(CheckMultFitsIn((uint32_t)kSquareRoot, (size_t)kSquareRoot, &r) &&
              r == kSquareRoot * kSquareRoot);

  static constexpr uint64_t kTooBigSquareRoot = 1ull << (sizeof(r) * 4);
  EXPECT_TRUE(CheckMultFitsIn(kTooBigSquareRoot, 1, &r) &&
              r == kTooBigSquareRoot);
  EXPECT_FALSE(CheckMultFitsIn(kTooBigSquareRoot, kTooBigSquareRoot, &r));
  EXPECT_FALSE(CheckMultFitsIn((size_t)kTooBigSquareRoot,
                               (uint64_t)kTooBigSquareRoot, &r));
  EXPECT_FALSE(CheckMultFitsIn((int64_t)kTooBigSquareRoot,
                               (uint64_t)kTooBigSquareRoot, &r));

  static constexpr T kMinValue = std::numeric_limits<T>::min();
  if (kMinValue < 0) {
    EXPECT_FALSE(CheckMultFitsIn((T)-1, (T)0, &r));
    EXPECT_FALSE(CheckMultFitsIn(kMaxValue, (T)-1, &r));
    EXPECT_FALSE(CheckMultFitsIn((T)1, kMinValue, &r));
    EXPECT_FALSE(CheckMultFitsIn(kMinValue, kMinValue, &r));
  }
  EXPECT_FALSE(CheckMultFitsIn((int32_t)-1, (T)0, &r));
  EXPECT_FALSE(CheckMultFitsIn((T)0, std::numeric_limits<int64_t>::min(), &r));
  EXPECT_FALSE(CheckMultFitsIn(std::numeric_limits<int64_t>::min(),
                               std::numeric_limits<int64_t>::min(), &r));
  EXPECT_FALSE(CheckMultFitsIn(std::numeric_limits<int64_t>::max(),
                               std::numeric_limits<uint64_t>::max(), &r));
  EXPECT_FALSE(CheckMultFitsIn(std::numeric_limits<int64_t>::min(),
                               std::numeric_limits<uint64_t>::max(), &r));
  EXPECT_FALSE(CheckMultFitsIn(std::numeric_limits<uint64_t>::max(),
                               std::numeric_limits<int64_t>::min(), &r));
}

TEST(ImageIoUtilTest, Mult) {
  TestMult<int8_t>();
  TestMult<uint8_t>();
  TestMult<int16_t>();
  TestMult<uint16_t>();
  TestMult<int32_t>();
  TestMult<uint32_t>();
  TestMult<int64_t>();
  TestMult<uint64_t>();

  TestMult<char>();
  TestMult<unsigned char>();
  TestMult<int>();
  TestMult<unsigned int>();
  TestMult<long>();                // NOLINT
  TestMult<unsigned long>();       // NOLINT
  TestMult<long long>();           // NOLINT
  TestMult<unsigned long long>();  // NOLINT
  TestMult<size_t>();
}

//------------------------------------------------------------------------------

}  // namespace
}  // namespace WP2
