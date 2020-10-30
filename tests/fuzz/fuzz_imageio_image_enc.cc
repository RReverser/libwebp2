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
//  Fuzzing of imageio (still image encoding).
//
// Author: Yannis Guyon (yguyon@google.com)

#include <cstddef>
#include <cstdint>
#include <iostream>
#include <string>

#include "imageio/file_format.h"
#include "imageio/image_dec.h"
#include "imageio/image_enc.h"
#include "src/wp2/base.h"
#include "src/wp2/format_constants.h"
#include "tests/fuzz/fuzz_utils.h"
#include "tests/include/helpers.h"

extern "C" int LLVMFuzzerTestOneInput(const uint8_t* data, size_t size) {
  // Input is considered random data and converted to parameters.
  WP2::testing::FuzzedParameters params(data, size);

  // Create a randomly sized valid buffer with random pixels.
  WP2::ArgbBuffer buffer = params.ExtractArgbBuffer(
      /*max_num_pixels=*/256 * 256,  // Not too big to avoid timeouts.
      /*max_width_or_height=*/WP2::kMaxBufferDimension);
  if (buffer.IsEmpty()) return 0;
  FuzzerTemporaryFile dst(data, 0);
  const std::string dst_path = dst.filename();
  const WP2::FileFormat file_format = params.ExtractOutputFormat();

  // Write it in a randomly chosen valid format.
  size_t file_size;
  const WP2Status write_status = WP2::SaveImage(
      buffer, dst_path.c_str(), /*overwrite=*/true, file_format, &file_size);
  if (write_status != WP2_STATUS_OK) {
    std::cerr << "Unable to write image to " << dst_path << " in format "
              << (int)file_format << ": " << WP2GetStatusMessage(write_status)
              << std::endl;
    abort();
  }

  // PGM is written as 16 bits but read at 8 so it's disabled. BMP is unhandled.
  if (file_format != WP2::FileFormat::PGM &&
      file_format != WP2::FileFormat::BMP) {
    // Read the freshly written file to ensure it's valid.
    const WP2Status read_again_status =
        WP2::ReadImage(dst_path.c_str(), &buffer, /*file_size=*/nullptr,
                       file_format, WP2::LogLevel::QUIET);
    if (read_again_status != WP2_STATUS_OK) {
      std::cerr << "Unable to read freshly written image " << dst_path << ": "
                << WP2GetStatusMessage(read_again_status) << std::endl;
      abort();
    }
  }
  return 0;
}
