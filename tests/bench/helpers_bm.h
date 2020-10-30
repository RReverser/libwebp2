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
// Benchmark corpus and tools. Uses Google Benchmark library.

#ifndef WP2_TESTS_BENCH_HELPERS_BM_H_
#define WP2_TESTS_BENCH_HELPERS_BM_H_

#include <cstddef>

#include "benchmark/benchmark.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"

namespace WP2 {
namespace testing {

// The benchmarks print only the indices and not the names of the files, that's
// why the association is written here for convenience:
static const char* const kFiles[]{
    "alpha_ramp.lossy.webp",  // 0
    "alpha_ramp.pam",         // 1
    "alpha_ramp.png",         // 2
    "alpha_ramp.tiff",        // 3
    "alpha_ramp.webp",        // 4
    "animation0.gif",         // 5
    "animation0.webp",        // 6
    "animation1.webp",        // 7
    "source0.pgm",            // 8
    "source0.ppm",            // 9
    "source1.itl.png",        // 10
    "source1.png",            // 11
    "source1_1x1.png",        // 12
    "source1_1x48.png",       // 13
    "source1_32x32.png",      // 14
    "source1_4x4.png",        // 15
    "source1_64x1.png",       // 16
    "source1_64x48.png",      // 17
    "source2.tiff",           // 18
    "source3.jpg",            // 19
    "source4.ll.webp",        // 20
    "source4.webp",           // 21
    "test_exif_xmp.webp",     // 22
};

static constexpr size_t kNumFiles = sizeof(kFiles) / sizeof(kFiles[0]);

}  // namespace testing
}  // namespace WP2

#endif  // WP2_TESTS_BENCH_HELPERS_BM_H_
