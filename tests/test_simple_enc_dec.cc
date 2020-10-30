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

// Simple test: encode/decode with default config.
// Does not require Google Test.

#include <cstdio>

#include "imageio/image_dec.h"
#include "src/wp2/decode.h"
#include "src/wp2/encode.h"

int main(int argc, const char* argv[]) {
  (void)argc;
  (void)argv;

  WP2Status status;
  WP2::ArgbBuffer src;
  status = WP2::ReadImage("testdata/source0.ppm", &src);
  if (status != WP2_STATUS_OK) {
    fprintf(stderr, "ReadImage failed: %i\n", status);
    return 1;
  }

  WP2::MemoryWriter memory_writer;
  status = WP2::Encode(src, &memory_writer);
  if (status != WP2_STATUS_OK) {
    fprintf(stderr, "Encode failed: %i\n", status);
    return 1;
  }

  WP2::ArgbBuffer output;
  status = WP2::Decode(memory_writer.mem_, memory_writer.size_, &output);
  if (status != WP2_STATUS_OK) {
    fprintf(stderr, "Decode failed: %i\n", status);
    return 1;
  }

  printf("TEST OK.\n");
  return 0;
}
