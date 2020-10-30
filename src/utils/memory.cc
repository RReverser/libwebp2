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
// Memory allocator and containers.
//
// Author: Skal (pascal.massimino@gmail.com)

#include <cstdlib>

#include "src/utils/thread_utils.h"
#include "src/utils/utils.h"

// If PRINT_MEM_INFO is defined, extra info (like total memory used, number of
// alloc/free etc) is printed. For debugging/tuning purpose only (it's slow,
// and not multi-thread safe!).
// An interesting alternative is valgrind's 'massif' tool:
//    http://valgrind.org/docs/manual/ms-manual.html
// Here is an example command line:
/*    valgrind --tool=massif --massif-out-file=massif.out \
               --stacks=yes --alloc-fn=WP2Malloc --alloc-fn=WP2Calloc
      ms_print massif.out
*/
// In addition:
// * if PRINT_MEM_TRAFFIC is defined, all the details of the malloc/free cycles
//   are printed.
// * if MALLOC_FAIL_AT is defined, the global environment variable
//   $MALLOC_FAIL_AT is used to simulate a memory error when calloc or malloc
//   is called for the nth time. Example usage:
//   export MALLOC_FAIL_AT=50 && ./examples/cwp2 input.png
// * if MALLOC_LIMIT is defined, the global environment variable $MALLOC_LIMIT
//   sets the maximum amount of memory (in bytes) made available to libwp2.
//   This can be used to emulate environment with very limited memory.
//   Example: export MALLOC_LIMIT=64000000 && ./examples/dwp2 picture.wp2

// #define PRINT_MEM_INFO
// #define PRINT_MEM_TRAFFIC
// #define MALLOC_LIMIT
// #define MALLOC_FAIL_AT   // can be defined alone, without PRINT_MEM_INFO

// For experimental set-up, we always define MALLOC_FAIL_AT
#if defined(WP2_EXPERIMENTAL) && !defined(MALLOC_FAIL_AT)
#define MALLOC_FAIL_AT
#endif

extern void WP2SetMallocFailAt(int num_malloc_before_failure);
static int countdown_to_fail = 0;     // 0 = off

#if defined(PRINT_MEM_INFO) || defined(MALLOC_LIMIT) || defined(MALLOC_FAIL_AT)
// A regular ThreadLock cannot be used as its ctor would call WP2Malloc(),
// which uses this 'thread_lock', which calls WP2Malloc() in its ctor etc.
// std::mutex might be forbidden but here it's for debugging.
#include <mutex>  // NOLINT
#endif

//------------------------------------------------------------------------------
// Checked memory allocation

#if defined(PRINT_MEM_INFO) || defined(MALLOC_LIMIT)

#include <cstdio>
#include <unordered_map>

static int num_malloc_calls = 0;
static int num_free_calls = 0;

static std::unordered_map<void*, size_t> all_blocks;
static size_t total_mem = 0;
static size_t total_mem_allocated = 0;
static size_t high_water_mark = 0;
static size_t mem_limit = 0;

static bool exit_registered = false;
static std::mutex mutex;

#if defined(PRINT_MEM_INFO)
static void PrintMemInfo() {
  std::lock_guard<std::mutex> lock(mutex);
  fprintf(stderr, "\nMEMORY INFO:\n");
  fprintf(stderr, "num calls to: malloc = %4d\n", num_malloc_calls);
  fprintf(stderr, "              free   = %4d\n", num_free_calls);
  fprintf(stderr, "total_mem: %u\n", (uint32_t)total_mem);
  fprintf(stderr, "total_mem allocated: %u\n", (uint32_t)total_mem_allocated);
  fprintf(stderr, "high-water mark: %u\n", (uint32_t)high_water_mark);
  all_blocks.clear();
}
#endif

static void Increment(int* const v) {
  std::lock_guard<std::mutex> lock(mutex);
  if (!exit_registered) {
#if defined(MALLOC_FAIL_AT)
    const char* const malloc_fail_at_str = getenv("MALLOC_FAIL_AT");
    if (malloc_fail_at_str != nullptr) {
      countdown_to_fail = atoi(malloc_fail_at_str);
    }
#endif
#if defined(MALLOC_LIMIT)
    const char* const malloc_limit_str = getenv("MALLOC_LIMIT");
    if (malloc_limit_str != nullptr) mem_limit = atoi(malloc_limit_str);
#endif
    (void)countdown_to_fail;
    (void)mem_limit;
#if defined(PRINT_MEM_INFO)
    atexit(PrintMemInfo);
#endif
    exit_registered = true;
  }
  ++*v;
}

static void AddMem(void* ptr, size_t size) {
  if (ptr != nullptr) {
    std::lock_guard<std::mutex> lock(mutex);
    all_blocks[ptr] = size;
    total_mem += size;
    total_mem_allocated += size;
#if defined(PRINT_MEM_TRAFFIC)
#if defined(MALLOC_FAIL_AT)
    fprintf(stderr, "fail-count: %5d [mem=%u]\n",
            num_malloc_calls, (uint32_t)total_mem);
#else
    fprintf(stderr, "Mem: %u (+%u)\n", (uint32_t)total_mem, (uint32_t)size);
#endif
#endif  // defined(PRINT_MEM_TRAFFIC)
    if (total_mem > high_water_mark) high_water_mark = total_mem;
  }
}

static void SubMem(void* ptr) {
  if (ptr != nullptr) {
    std::lock_guard<std::mutex> lock(mutex);
    size_t size;
    if (all_blocks.find(ptr) == all_blocks.end()) {
      fprintf(stderr, "Invalid pointer free! (%p)\n", ptr);
      abort();
    }
    size = all_blocks[ptr];
    all_blocks.erase(ptr);
    total_mem -= size;
#if defined(PRINT_MEM_TRAFFIC)
    fprintf(stderr, "Mem: %u (-%u)\n", (uint32_t)total_mem, (uint32_t)size);
#endif
  }
}

#else

#if defined(MALLOC_FAIL_AT)   // very light version for this case

static bool exit_registered = false;
static int num_malloc_calls = 0;
static int num_free_calls = 0;
static std::mutex mutex;

static void Increment(int* const) {
  std::lock_guard<std::mutex> lock(mutex);
  if (!exit_registered) {
    const char* const malloc_fail_at_str = getenv("MALLOC_FAIL_AT");
    if (malloc_fail_at_str != nullptr) {
      countdown_to_fail = atoi(malloc_fail_at_str);
    }
    exit_registered = true;
  }
}

void WP2SetMallocFailAt(int num_malloc_before_failure) {
  std::lock_guard<std::mutex> lock(mutex);
  countdown_to_fail = num_malloc_before_failure;
  exit_registered = true;  // To prevent getenv("MALLOC_FAIL_AT") in Increment()
}

#else

#define Increment(v) do {} while (false)

void WP2SetMallocFailAt(int num_malloc_before_failure) {
  countdown_to_fail = num_malloc_before_failure;
}

#endif  // defined(MALLOC_FAIL_AT)

#define AddMem(p, s) do {} while (false)
#define SubMem(p)    do {} while (false)

#endif

// Returns false in case of overflow of nmemb * size.
static bool CheckSizeArgumentsOverflow(uint64_t nmemb, size_t size) {
  const uint64_t total_size = nmemb * size;
  if (nmemb == 0) return true;
  if ((uint64_t)size > WP2_MAX_ALLOCABLE_MEMORY / nmemb) return false;
  if (total_size != (size_t)total_size) return false;

#if defined(MALLOC_FAIL_AT) || defined(MALLOC_LIMIT)
  {
    std::lock_guard<std::mutex> lock(mutex);

#if defined(MALLOC_FAIL_AT)
    if (countdown_to_fail > 0 && --countdown_to_fail == 0) {
      return false;  // fake fail!
    }
#endif
#if defined(MALLOC_LIMIT)
    if (mem_limit > 0) {
      const uint64_t new_total_mem = (uint64_t)total_mem + total_size;
      if (new_total_mem != (size_t)new_total_mem || new_total_mem > mem_limit) {
        return false;  // fake fail!
      }
    }
#endif
  }
#else
  // we insert the countdown_to_fail mechanism without mutex, for testing
  if (countdown_to_fail > 0 && --countdown_to_fail == 0) {
    return false;  // fake fail!
  }
#endif  // defined(MALLOC_FAIL_AT) || defined(MALLOC_LIMIT)
  return true;
}

void* WP2Malloc(uint64_t nmemb, size_t size) {
  Increment(&num_malloc_calls);
  if (!CheckSizeArgumentsOverflow(nmemb, size)) return nullptr;
  assert(nmemb * size > 0);
  void* const ptr = malloc((size_t)(nmemb * size));
  AddMem(ptr, (size_t)(nmemb * size));
  return ptr;
}

void* WP2Calloc(uint64_t nmemb, size_t size) {
  Increment(&num_malloc_calls);
  if (!CheckSizeArgumentsOverflow(nmemb, size)) return nullptr;
  assert(nmemb * size > 0);
  void* const ptr = calloc((size_t)nmemb, size);
  AddMem(ptr, (size_t)(nmemb * size));
  return ptr;
}

void* WP2Realloc(void* const ptr, uint64_t nmemb, size_t size) {
  Increment(&num_malloc_calls);
  if (!CheckSizeArgumentsOverflow(nmemb, size)) return nullptr;
  assert(nmemb * size > 0);
  void* const new_ptr = realloc(ptr, (size_t)(nmemb * size));
  // Let the caller free 'ptr' when 'new_ptr' is null - like realloc().
  if (new_ptr != nullptr) SubMem(ptr);
  AddMem(new_ptr, (size_t)(nmemb * size));
  return new_ptr;
}

void WP2Free(void* ptr) {
  if (ptr != nullptr) {
    Increment(&num_free_calls);
    SubMem(ptr);
  }
  free(ptr);
}

constexpr WP2Allocable::NoThrow WP2Allocable::nothrow;

//------------------------------------------------------------------------------
