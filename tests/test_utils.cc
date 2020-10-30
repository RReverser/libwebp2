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

// Tests functions in utils.h

#include <initializer_list>

#include "include/helpers.h"
#include "src/utils/utils.h"

namespace WP2 {
namespace {

//------------------------------------------------------------------------------

constexpr int null_ptr = -1;

bool TestStrnstr(const char* const str, size_t count, const char* const target,
                 int expected_offset) {
  return (strnstr(str, count, target) ==
          ((expected_offset == null_ptr) ? nullptr : (str + expected_offset)));
}

TEST(Utils, Strnstr) {
  for (uint32_t count : {0, 1, 2, 3, 4, 10000}) {
    ASSERT_TRUE(TestStrnstr("", count, "", 0));
    ASSERT_TRUE(TestStrnstr("abc", count, "", 0));
    ASSERT_TRUE(TestStrnstr("abc", count, "a", (count < 1) ? null_ptr : 0));
    ASSERT_TRUE(TestStrnstr("abc", count, "b", (count < 2) ? null_ptr : 1));
    ASSERT_TRUE(TestStrnstr("abc", count, "c", (count < 3) ? null_ptr : 2));
    ASSERT_TRUE(TestStrnstr("abc", count, "abc", (count < 3) ? null_ptr : 0));
    ASSERT_TRUE(TestStrnstr("abc", count, "ab", (count < 2) ? null_ptr : 0));
    ASSERT_TRUE(TestStrnstr("abc", count, "bc", (count < 3) ? null_ptr : 1));

    ASSERT_TRUE(TestStrnstr("\0abc", count, "abc", null_ptr));
    ASSERT_TRUE(TestStrnstr("abc", count, "\0abc", 0));
    ASSERT_TRUE(TestStrnstr("ababab\0c", count, "ababc", null_ptr));

    ASSERT_TRUE(
        TestStrnstr("abababc", count, "ababc", (count < 2 + 5) ? null_ptr : 2));
    ASSERT_TRUE(
        TestStrnstr("YUV4MPEG2 W900 H1100 F24:1 Ip A0:0 C420jpeg\nFRAME", count,
                    " C420", (count < 34 + 5) ? null_ptr : 34));

    ASSERT_TRUE(TestStrnstr(nullptr, count, "", null_ptr));
    ASSERT_TRUE(TestStrnstr(nullptr, count, "a", null_ptr));
    ASSERT_TRUE(TestStrnstr("", count, nullptr, null_ptr));
    ASSERT_TRUE(TestStrnstr("a", count, nullptr, null_ptr));
    ASSERT_TRUE(TestStrnstr(nullptr, count, nullptr, null_ptr));
  }
}

//------------------------------------------------------------------------------

class TrivialClass {
 public:
  TrivialClass() = default;
  TrivialClass(const TrivialClass&) = delete;
  TrivialClass(TrivialClass&& other) noexcept { TrivialMoveCtor(this, &other); }
  ~TrivialClass() { assert(some_pointer_ == nullptr); }
  uint32_t* some_pointer_ = nullptr;
};

TEST(Utils, TrivialMoveCtor) {
  uint32_t some_value = 2;
  TrivialClass moved_instance;
  ASSERT_EQ(moved_instance.some_pointer_, nullptr);
  moved_instance.some_pointer_ = &some_value;

  TrivialClass created_instance(std::move(moved_instance));
  ASSERT_EQ(created_instance.some_pointer_, &some_value);
  ASSERT_EQ(moved_instance.some_pointer_, nullptr);  // NOLINT
  created_instance.some_pointer_ = nullptr;
}

//------------------------------------------------------------------------------

}  // namespace
}  // namespace WP2
