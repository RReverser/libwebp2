// Copyright 2019 Google LLC
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//     https://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.
// -----------------------------------------------------------------------------
//
// Misc. common utility functions
//
// Authors: Skal (pascal.massimino@gmail.com)
//

#ifndef WP2_UTILS_UTILS_H_
#define WP2_UTILS_UTILS_H_

#ifdef HAVE_CONFIG_H
#include "src/wp2/config.h"
#endif

#include <cassert>
#include <cstddef>
#include <cstring>
#include <initializer_list>
#include <type_traits>
#include <utility>
#include <vector>
#if defined(WP2_BITTRACE) || defined(WP2_TRACE) || defined(WP2_ENC_DEC_MATCH)
#include <string>
#endif

#include "src/wp2/base.h"
#include "src/wp2/format_constants.h"

//------------------------------------------------------------------------------
// Tracing and printing

#if defined(WP2_TRACE)
// we'll supply a default impl if WP2_TRACE is not defined
extern void WP2Trace(const char* format, const std::string prefix, ...);
#else
#define WP2Trace(format, ...) (void)format
#endif  // WP2_TRACE

// Safe replacement for printf(fmt, ...) when fmt is not a literal
void WP2Print(const char* const fmt, ...);

// Safe replacement for sprintf(fmt, ...) when fmt is not a literal.
// Writes the data to str. Used only for debugging pruposes, hence the
// replacement with an empty macro otherwise
#if defined(WP2_BITTRACE) || defined(WP2_TRACE) || defined(WP2_ENC_DEC_MATCH)
void WP2SPrint(std::string* const str, const char* const fmt, ...);
#else
#define WP2SPrint(str, ...) (void)str
#endif

// handy error-logging macros
void WP2ErrorLog(const char* const file, int line, WP2Status s,
                 const char* const extra_message = nullptr);

// performs "if (status != WP2_STATUS_OK) return status;"
#define WP2_CHECK_STATUS(status)              \
  do {                                        \
    const WP2Status __S__ = (status);         \
    if (__S__ != WP2_STATUS_OK) {             \
      WP2ErrorLog(__FILE__, __LINE__, __S__); \
      return __S__;                           \
    }                                         \
  } while (0)

// performs "if (!cond) return error_status;"
#define WP2_CHECK_OK(cond, error_status)             \
  do {                                               \
    if (!(cond)) {                                   \
      WP2ErrorLog(__FILE__, __LINE__, error_status); \
      return error_status;                           \
    }                                                \
  } while (0)

// specialized for vector alloc checks
#define WP2_CHECK_ALLOC_OK(cond) WP2_CHECK_OK(cond, WP2_STATUS_OUT_OF_MEMORY)

// performs "if (status != WP2_STATUS_OK) assert(false);"
#define WP2_ASSERT_STATUS(status)             \
  do {                                        \
    const WP2Status __S__ = (status);         \
    if (__S__ != WP2_STATUS_OK) {             \
      WP2ErrorLog(__FILE__, __LINE__, __S__); \
    }                                         \
    assert(__S__ == WP2_STATUS_OK);           \
  } while (0)

// Can be overridden and called instead of WP2ErrorLog() if 'wp2_error_tracer'
// is not null. Beware of concurrent accesses during multithreading.
#if defined(WP2_TRACE) || defined(WP2_ERROR_TRACE)
class WP2ErrorTracer {
 public:
  virtual ~WP2ErrorTracer() = default;
  virtual void Log(const char* const file, int line, WP2Status s,
                   const char* const extra_message);
};
extern WP2ErrorTracer* wp2_error_tracer;
#endif  // defined(WP2_TRACE) || defined(WP2_ERROR_TRACE)

//------------------------------------------------------------------------------
// Memory allocation

// This is the maximum memory amount that libwp2 will ever try to allocate.
#ifndef WP2_MAX_ALLOCABLE_MEMORY
#if SIZE_MAX > (1ULL << 34)
#define WP2_MAX_ALLOCABLE_MEMORY (1ULL << 34)
#else
// For 32-bit targets keep this below INT_MAX to avoid valgrind warnings.
#define WP2_MAX_ALLOCABLE_MEMORY ((1ULL << 31) - (1 << 16))
#endif
#endif  // WP2_MAX_ALLOCABLE_MEMORY

//------------------------------------------------------------------------------
// Memory functions

// Size-checking safe malloc/calloc: verify that the requested size is not too
// large, or return nullptr.
// Note that these expect the second argument type to be 'size_t'
// in order to favor the "calloc(num_foo, sizeof(foo))" pattern.
WP2_EXTERN void* WP2Malloc(uint64_t num_elements, size_t element_size);
WP2_EXTERN void* WP2Calloc(uint64_t num_elements, size_t element_size);
// Don't forget to deallocate 'ptr' if WP2Realloc() returns nullptr.
WP2_EXTERN void* WP2Realloc(void* const ptr, uint64_t num_elements,
                            size_t element_size);

// Releases memory allocated by the functions above or by WP2Decode*().
WP2_EXTERN void WP2Free(void* ptr);

// Inherit to call WP2Malloc and WP2Free on 'new' and 'delete'.
struct WP2Allocable {
  // Tag similar to std::nothrow.
  struct NoThrow {
  } static constexpr nothrow = NoThrow();  // NOLINT (Camelcase)

  // No virtual dtor to not force inherited objects to have a virtual table.

  // Enforce usage of 'WP2Allocable::nothrow' to warn the user.
  void* operator new(std::size_t size) = delete;
  void* operator new(std::size_t size, const std::nothrow_t&) = delete;
  void* operator new[](std::size_t size) = delete;
  void* operator new[](std::size_t size, const std::nothrow_t&) = delete;
  void* operator new(std::size_t size, const NoThrow&) noexcept {
    return WP2Malloc(1, size);
  }
  void* operator new[](std::size_t size, const NoThrow&) noexcept {
    return WP2Malloc(1, size);
  }
  void operator delete(void* ptr) noexcept { WP2Free(ptr); }
  void operator delete[](void* ptr) noexcept { WP2Free(ptr); }
  // Only called when a ctor throws during 'new (WP2Allocable::nothrow)'.
  void operator delete(void* ptr, const NoThrow&) noexcept { WP2Free(ptr); }
  void operator delete[](void* ptr, const NoThrow&) noexcept { WP2Free(ptr); }
};

//------------------------------------------------------------------------------

namespace WP2 {

// Convenient pointer/length tuple. Does not own or copy anything.
struct DataView {
  bool IsEmpty() const { return (size == 0); }

  const uint8_t* bytes;
  size_t size;
};

// Short move constructor. Should only be used for plain-old-data final classes
// that do not refer to themselves (pointer to another local member etc.).
template <typename T>
void TrivialMoveCtor(T* const created_instance, T* const moved_instance) {
  std::memcpy((void*)created_instance, (const void*)moved_instance, sizeof(T));
  new(moved_instance) T();
}

template <typename T>
void SetArray(T* array, const std::initializer_list<T>& values) {
  for (const T& value : values) *(array++) = value;
}

#define STATIC_ASSERT_ARRAY_SIZE(array, size)                           \
  static_assert(                                                        \
      sizeof(array) /                                                   \
              sizeof(std::remove_all_extents<decltype(array)>::type) == \
          (size),                                                       \
      "Expected " #array " to contain " #size " elements");

// Same as std::strstr() but considers 'str[count]' to be the null-character.
const char* strnstr(const char* str, size_t count, const char* const target);

// Same as std::strncmp() but clamps 'count' to std::strlen(rhs).
int strlcmp(const char* lhs, const char* rhs, size_t count);

}  // namespace WP2

//------------------------------------------------------------------------------

#endif  /* WP2_UTILS_UTILS_H_ */
