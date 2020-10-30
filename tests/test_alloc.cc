// Copyright 2017 Google LLC
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
// Test allocator and vector
//
// Author: Skal (pascal.massimino@gmail.com)

#include <cstdint>
#include <cstdio>

#include "src/utils/vector.h"
#include "include/helpers.h"

namespace WP2 {
namespace {

//------------------------------------------------------------------------------

template <class T, class U>
bool DoTest(U v0) {
  T vec1;
  EXPECT_TRUE(vec1.resize(30));  // still uninitialized
  EXPECT_EQ(vec1.size(), 30u);

  T vec2;
  EXPECT_TRUE(vec2.resize(37));
  for (auto& v : vec2) v = v0;

  T vec3;
  for (auto& v : vec2) {
    EXPECT_TRUE(vec3.push_back(v, /*resize_if_needed=*/true));
  }
  for (const auto& v : vec3) EXPECT_EQ((U)v, v0);

  size_t i = 0;
  for (auto& v : vec3) v = i++;
  vec2.clear();
  for (auto& v : vec3) {
    EXPECT_TRUE(vec2.push_back(v, /*resize_if_needed=*/true));
  }
  // Make sure we keep the same values when expanding.
  EXPECT_TRUE(vec2.resize(vec3.size() * 2));
  for (i = 0; i < vec3.size(); ++i) EXPECT_EQ(vec2[i], vec3[i]);
  // Make sure we keep the same values when shrinking.
  EXPECT_TRUE(vec2.resize(vec3.size() / 2));
  for (i = 0; i < vec2.size(); ++i) EXPECT_EQ(vec2[i], vec3[i]);

  for (size_t k = 1; k < 23; ++k) {
    if (!vec3.resize(k)) {
      printf("Resizing to size %d failed!\n", (int)k);
      return false;
    }
    EXPECT_EQ(vec3.size(), k);
  }
  return true;
}

// A test structure to explore non POD types.
struct TestStruct {
  explicit TestStruct(uint32_t value = 37) noexcept : i(value) {}
  TestStruct& operator=(uint32_t value) {
    v = 13.;
    i = value;
    return *this;
  }
  double v = -2;
  uint32_t i = 3;
  char c = 'c';
  bool operator==(const TestStruct& other) const {
    return i == other.i && c == other.c && v == other.v;
  }
  bool operator==(int value) const {
    // compare to objects constructed with operator=(int)
    return (int)i == value && c == 'c' && v == 13;
  }
  explicit operator uint32_t() const { return i; }
};

TEST(TestAlloc, TestVectors) {
  ASSERT_TRUE(DoTest<Vector_s8>(-3));
  ASSERT_TRUE(DoTest<Vector_u8>(24u));
  ASSERT_TRUE(DoTest<Vector_s16>(-31464));
  ASSERT_TRUE(DoTest<Vector_u16>(65030u));
  ASSERT_TRUE(DoTest<Vector_s32>(-2));
  ASSERT_TRUE(DoTest<Vector_u32>(35264u));
  ASSERT_TRUE(DoTest<Vector_s64>(-0x3254235423ull));
  ASSERT_TRUE(DoTest<Vector_u64>(0x54254235423ull));
  ASSERT_TRUE(DoTest<Vector<TestStruct>>(TestStruct(3)));
  ASSERT_TRUE(DoTest<Vector_f>(.4f));
}

//------------------------------------------------------------------------------

struct CtorStruct {
  static uint32_t counter;
  CtorStruct() noexcept { ++counter; }
  CtorStruct(const CtorStruct&) { ++counter; }
  CtorStruct(CtorStruct&&) noexcept { ++counter; }
};
struct DtorStruct : public CtorStruct {  // Cannot be used by VectorNoCtor.
  ~DtorStruct() { --counter; }
};
uint32_t CtorStruct::counter;

TEST(TestAlloc, WP2VectorCtor) {
  CtorStruct::counter = 0u;
  Vector<DtorStruct> vector;
  ASSERT_TRUE(vector.resize(5));
  ASSERT_EQ(CtorStruct::counter, 5u);

  {
    Vector<DtorStruct> other_vector;
    swap(vector, other_vector);
    ASSERT_EQ(other_vector.size(), 5u);
    ASSERT_EQ(CtorStruct::counter, 5u);
  }
  ASSERT_EQ(vector.size(), 0u);
  ASSERT_EQ(CtorStruct::counter, 0u);

  for (int i = 0; i < 5; ++i) {
    ASSERT_TRUE(vector.push_back(DtorStruct(), /*resize_if_needed=*/true));
  }
  ASSERT_EQ(CtorStruct::counter, 5u);

  ASSERT_TRUE(vector.shrink_to_fit());
  ASSERT_TRUE(vector.reserve(10));
  ASSERT_EQ(CtorStruct::counter, 5u);
  ASSERT_TRUE(vector.reserve(10));
  ASSERT_TRUE(vector.shrink_to_fit());
  ASSERT_EQ(CtorStruct::counter, 5u);

  vector.erase(vector.begin(), vector.begin() + 3);
  ASSERT_EQ(CtorStruct::counter, 2u);
  vector.erase(vector.begin());
  ASSERT_EQ(CtorStruct::counter, 1u);
  vector.clear();
  ASSERT_EQ(CtorStruct::counter, 0u);
  ASSERT_TRUE(vector.shrink_to_fit());
}

TEST(TestAlloc, VectorNoCtor) {
  CtorStruct::counter = 0u;
  VectorNoCtor<CtorStruct> vector;
  ASSERT_TRUE(vector.resize(5));
  ASSERT_EQ(CtorStruct::counter, 0u);
}

//------------------------------------------------------------------------------

}  // namespace
}  // namespace WP2
