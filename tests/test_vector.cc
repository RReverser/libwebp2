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
// Test custom vector API

#include "include/helpers.h"
#include "src/utils/vector.h"

namespace WP2 {
namespace {

const int kMagic = 38, kMagic2 = 113;

class Foo {
 public:
  Foo(int x = kMagic) : x_(x) {}       // NOLINT (explicitly no 'explicit')
  operator int() const { return x_; }  // NOLINT (explicitly no 'explicit')

 private:
  int x_ = kMagic;
};

TEST(VectorTest, NoCtor) {
  VectorNoCtor<int> v;
  EXPECT_TRUE(v.resize(100));
  Vector<int> w;
  EXPECT_TRUE(w.resize(100));

#if WP2_HAVE_MSAN
  // Use MemorySanitizer to check VectorNoCtor::resize() does not initialize
  // the memory while Vector::resize() does.
  //
  // __msan_test_shadow(const void *x, uptr size) returns the offset of the
  // first (at least partially) poisoned byte in the range, or -1 if the whole
  // range is good.
  for (size_t i = 0; i < 100; ++i) {
    EXPECT_EQ(__msan_test_shadow(&v[i], sizeof(int)), 0);
    EXPECT_EQ(__msan_test_shadow(&w[i], sizeof(int)), -1);
  }
#endif
}

template<typename T> void DoTestCtor(const T& v0) {
  Vector<T> V;
  EXPECT_TRUE(V.resize(100));
  for (const T& v : V) EXPECT_EQ(v, v0);
}

TEST(VectorTest, Constructor) {
  DoTestCtor(Foo());
  DoTestCtor(0);
  DoTestCtor(0u);
  DoTestCtor((uint16_t)0u);
  DoTestCtor((int16_t)0);
  DoTestCtor(0.f);
  DoTestCtor(0.);
}

template<typename T> void DoTestPushBack(const T& vector_type) {
  (void)vector_type;

  T v;
  EXPECT_TRUE(v.reserve(8));
  EXPECT_EQ(v.size(), 0u);
  EXPECT_EQ(v.capacity(), 8u);

  EXPECT_TRUE(v.push_back(kMagic));
  EXPECT_EQ(v.size(), 1u);
  EXPECT_EQ((int)v[0], kMagic);

  EXPECT_TRUE(v.push_back(kMagic2));
  EXPECT_EQ(v.size(), 2u);
  EXPECT_EQ((int)v[0], kMagic);
  EXPECT_EQ((int)v[1], kMagic2);

  EXPECT_TRUE(v.resize(1));   // down-size back to 1
  EXPECT_EQ(v.size(), 1u);
  EXPECT_EQ((int)v[0], kMagic);
}

TEST(VectorTest, PushBack) {
  DoTestPushBack(Vector_u8());
  DoTestPushBack(Vector_s16());
  DoTestPushBack(Vector_u32());
  DoTestPushBack(VectorNoCtor<float>());
  DoTestPushBack(Vector<Foo>());
}

// Copy constructor and assignment are deleted, but move constructor and
// assignment are OK.
template<typename T> void DoTestMove(const T& vector_type) {
  (void)vector_type;

  T v1;
  EXPECT_TRUE(v1.reserve(4));
  EXPECT_TRUE(v1.push_back(2));
  EXPECT_TRUE(v1.push_back(3));
  EXPECT_TRUE(v1.push_back(5));
  EXPECT_TRUE(v1.push_back(7));

  // Move constructor.
  T v2(std::move(v1));
  EXPECT_EQ(v2.size(), 4u);
  EXPECT_EQ((int)v2[0], 2);
  EXPECT_EQ((int)v2[1], 3);
  EXPECT_EQ((int)v2[2], 5);
  EXPECT_EQ((int)v2[3], 7);

  // Move constructor again.
  T v3 = std::move(v2);
  EXPECT_EQ(v3.size(), 4u);
  EXPECT_EQ((int)v3[0], 2);
  EXPECT_EQ((int)v3[1], 3);
  EXPECT_EQ((int)v3[2], 5);
  EXPECT_EQ((int)v3[3], 7);
}

TEST(VectorTest, Move) {
  DoTestMove(Vector<int>());
  DoTestMove(Vector_u32());
  DoTestMove(Vector_s16());
  DoTestMove(Vector<uint8_t>());
  DoTestMove(Vector<float>());
}

template<typename T> void DoTestErase(const T& vector_type) {
  (void)vector_type;

  T v;
  EXPECT_TRUE(v.reserve(4));
  EXPECT_TRUE(v.push_back(2));
  EXPECT_TRUE(v.push_back(3));
  EXPECT_TRUE(v.push_back(5));
  EXPECT_TRUE(v.push_back(7));

  EXPECT_EQ(v.size(), 4u);
  EXPECT_EQ((int)v[0], 2);
  EXPECT_EQ((int)v[1], 3);
  EXPECT_EQ((int)v[2], 5);
  EXPECT_EQ((int)v[3], 7);

  v.erase(v.begin());
  EXPECT_EQ(v.size(), 3u);
  EXPECT_EQ((int)v[0], 3);
  EXPECT_EQ((int)v[1], 5);
  EXPECT_EQ((int)v[2], 7);
}

TEST(VectorTest, Erase) {
  DoTestErase(Vector<int>());
  DoTestErase(Vector_u32());
  DoTestErase(Vector_s16());
  DoTestErase(Vector<uint8_t>());
  DoTestErase(Vector<Foo>());
}

template<typename T> void DoTestCopy(const T& vector_type) {
  (void)vector_type;

  T v, copy;
  uint32_t N = 61;
  ASSERT_TRUE(v.reserve(N));
  for (uint32_t i = 0; i < N; ++i) {
    int value = (i * 13 - 47) % 22;
    EXPECT_TRUE(v.push_back(value));
  }
  ASSERT_TRUE(copy.copy_from(v));
  for (uint32_t i = 0; i < v.size(); ++i) {
    EXPECT_EQ(v[i], copy[i]);
  }
  ASSERT_TRUE(v.resize(v.size() / 2));
  for (uint32_t i = 0; i < v.size(); ++i) {
    EXPECT_EQ(v[i], copy[i]);
  }
}

TEST(VectorTest, Copy) {
  DoTestCopy(Vector<int>());
  DoTestCopy(Vector_u32());
  DoTestCopy(Vector_s16());
  DoTestCopy(Vector<uint8_t>());
  DoTestCopy(Vector<Foo>());
}

}  // namespace
}  // namespace WP2
