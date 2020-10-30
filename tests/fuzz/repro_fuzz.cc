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
//  Reproducer binary for easier debugging of reproducer files.
//  Usage:
//    repro_fuzz [reproducer file...]
//
// Author: Yannis Guyon (yguyon@google.com)

#include <cstdint>
#include <iostream>
#include <vector>

#include "examples/example_utils.h"
#include "imageio/image_dec.h"
#include "src/wp2/base.h"

extern "C" int LLVMFuzzerTestOneInput(const uint8_t* data, size_t size);

int main(int argc, char* argv[]) {
  if (argc < 2) {
    std::cerr << "Missing reproducer file" << std::endl;
    return 1;
  }

  for (int i = 1; i < argc; ++i) {
    WP2::Data data;
    if (WP2::IoUtilReadFile(argv[i], &data) != WP2_STATUS_OK) abort();
    std::cout << "Reproducing '" << argv[i] << "' (" << data.size << " bytes)"
              << std::endl;

    const double start_time = GetStopwatchTime();
    LLVMFuzzerTestOneInput(data.bytes, data.size);
    const double end_time = GetStopwatchTime();
    const double time_spent = end_time - start_time;
    std::cout << "Time spent: " << time_spent << " seconds" << std::endl;
  }
  return 0;
}
