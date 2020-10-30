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
//  Fuzzing tools.
//
// Author: Yannis Guyon (yguyon@google.com)

#ifndef WP2_TESTS_FUZZ_FUZZ_UTILS_H_
#define WP2_TESTS_FUZZ_FUZZ_UTILS_H_

#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

#include "imageio/file_format.h"
#include "src/utils/utils.h"
#include "src/wp2/base.h"
#include "src/wp2/decode.h"
#include "src/wp2/encode.h"

// This header can be found in some GitHub oss-fuzz repositories.
#include "fuzzer_temp_file.h"

namespace WP2 {
namespace testing {

// Limit read image size to avoid out-of-memory and timeout issues.
// Hint:  kMaxNumPixels * kMaxBytesPerPixel (4) * kMaxNumFullCanvasAlloc (3)
//        should be below kMaxFuzzMemory (2 GB).
constexpr size_t kMaxNumPixels = 1u << 26;

// Used to extract test parameters from a bitstream of any size and content.
class FuzzedParameters {
 public:
  FuzzedParameters(const uint8_t data[], size_t data_size)
      : data_(data), data_size_(data_size) {}

  // Returns a value in range [0..max_value] or [min_value..max_value].
  // Only powers of two are uniformly distributed.
  uint32_t ExtractUInt32(uint32_t max_value);
  uint32_t ExtractUInt32(uint32_t min_value, uint32_t max_value);

  bool ExtractBool();

  // Returns the name, bytes and size of a valid image from testdata folder or
  // false if an error occurred.
  bool ExtractSource(const char** const file_name,
                     const uint8_t** const file_data, size_t* const file_size);

  // Returns a valid output format.
  FileFormat ExtractOutputFormat();

  // Returns a randomly sized buffer with random pixels. Returns an empty buffer
  // if there is not enough data or if the allocation failed.
  ArgbBuffer ExtractArgbBuffer(uint32_t max_num_pixels = 256 * 256,
                               uint32_t max_width_or_height = kImageDimMax,
                               WP2SampleFormat format = WP2_Argb_32);

  // Returns a (possibly invalid) wp2 bitstream, starting with the wp2 tag.
  bool TryExtractWP2Bitstream(Data* const encoded_data);

  // Extracts main encoding settings (quality, speed etc.) and sets 'config'.
  void ExtractSimpleConfig(uint32_t max_speed, EncoderConfig* const config);
  // Extracts other settings (partition_set etc.) and sets 'config'.
  void ExtractAdvancedConfig(bool fast, EncoderConfig* const config);

  // Extracts decoding settings.
  void ExtractConfig(DecoderConfig* const config);

 protected:
  const uint8_t* const data_;
  const size_t data_size_;
  size_t bit_pos_ = 0;
};

// Verifies that 'original' can be encoded with 'config' or aborts.
// Also checks that the decoded image is similar if 'compare_after_decoding'.
int TestImageEncConfig(const ArgbBuffer& original, const EncoderConfig& config,
                       bool compare_after_decoding, float expected_distortion);

int TestAnimEncConfig(const std::vector<ArgbBuffer>& original_frames,
                      const std::vector<uint32_t>& durations_ms,
                      uint32_t loop_count, const EncoderConfig& config,
                      bool compare_after_decoding, float expected_distortion);

// Verifies that 'encoded_data' does not crash when decoded with 'config'.
// Aborts if it's not decoded without error in case of 'expected_success'.
int TestImageDecConfig(DataView encoded_data, const DecoderConfig& config,
                       bool expected_success);

// Exports 'original' and 'decompressed' to PNG files at 'path' for debugging.
void SaveBeforeAfter(const ArgbBuffer& original, const ArgbBuffer& decompressed,
                     const std::string& path);

}  // namespace testing
}  // namespace WP2

#endif  // WP2_TESTS_FUZZ_FUZZ_UTILS_H_
