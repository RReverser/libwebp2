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
//  Fuzzing of imageio (still and animated image decoding).
//
// Author: Yannis Guyon (yguyon@google.com)

#include <cstddef>
#include <cstdint>

#include "./fuzz_utils.h"
#include "imageio/anim_image_dec.h"
#include "src/wp2/base.h"

extern "C" int LLVMFuzzerTestOneInput(const uint8_t* data, size_t size) {
  // Truncate lengthy inputs.
  static constexpr uint32_t kMaxSize = 1024 * 1024;
  if (size > kMaxSize) size = kMaxSize;

  // Input data is an encoded animation (webp, gif etc.) that may be invalid.
  // It can also be a still image (jpeg, png etc.).
  WP2::ArgbBuffer buffer;
  WP2::ImageReader image_reader(data, size, &buffer, WP2::FileFormat::AUTO,
                                WP2::LogLevel::QUIET,
                                WP2::testing::kMaxNumPixels);
  bool is_last_frame = true;
  while (image_reader.ReadFrame(&is_last_frame) == WP2_STATUS_OK &&
         !is_last_frame) {
  }
  return 0;
}
