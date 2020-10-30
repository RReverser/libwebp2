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

#include "./fuzz_utils.h"

#include <algorithm>
#include <cassert>
#include <vector>

#include "src/wp2/format_constants.h"
#include "imageio/image_enc.h"

#include "examples/example_utils.h"
#include "imageio/imageio_util.h"
#include "tests/include/helpers.h"

namespace WP2 {
namespace testing {

//------------------------------------------------------------------------------

uint32_t FuzzedParameters::ExtractUInt32(uint32_t max_value) {
  uint32_t extracted_value = 0;
  int max_extracted_value = 0;

  while (bit_pos_ < 8 * data_size_ && max_extracted_value < max_value) {
    extracted_value = (extracted_value << 1);
    max_extracted_value = (max_extracted_value << 1) | 1u;

    const uint8_t mask = (uint8_t)1 << (bit_pos_ & 7);
    if (data_[bit_pos_ >> 3] & mask) extracted_value |= 1;

    ++bit_pos_;
  }

  // Also handles the case where (max_value + 1) would overflow.
  if (max_value < max_extracted_value) {
    return extracted_value % (max_value + 1);
  }
  return extracted_value;
}

uint32_t FuzzedParameters::ExtractUInt32(uint32_t min_value,
                                         uint32_t max_value) {
  assert(max_value >= min_value);
  return min_value + ExtractUInt32(max_value - min_value);
}

bool FuzzedParameters::ExtractBool() { return (ExtractUInt32(1) == 1); }

//------------------------------------------------------------------------------

bool FuzzedParameters::ExtractSource(const char** const file_name,
                                     const uint8_t** const file_data,
                                     size_t* const file_size) {
  std::vector<std::string> files;
  Data data;
  if (WP2GetFilesIn(GetTestDataPath(""), &files, /*recursive=*/false) &&
      files.size() > 0) {
    *file_name = files[ExtractUInt32((uint32_t)files.size() - 1)].c_str();
    if (IoUtilReadFile(*file_name, &data) == WP2_STATUS_OK) {
      *file_data = data.bytes;
      *file_size = data.size;
      return true;
    }
  }
  return false;
}

FileFormat FuzzedParameters::ExtractOutputFormat() {
  const FileFormat kOutputFormats[] = {
      FileFormat::PNG, FileFormat::TIFF, FileFormat::BMP,
      FileFormat::PAM, FileFormat::PGM,  FileFormat::PPM,
  };
  const uint32_t kNumOutputFormats =
      sizeof(kOutputFormats) / sizeof(kOutputFormats[0]);

  return kOutputFormats[ExtractUInt32(kNumOutputFormats - 1)];
}

//------------------------------------------------------------------------------

namespace {

// Returns all unique divisors of 'number', unsorted.
std::vector<uint32_t> GetDivisors(uint32_t number, uint32_t max_divisor) {
  std::vector<uint32_t> divisors;
  uint32_t divisor = (uint32_t)std::lround(
      std::ceil(number / (double)std::min(max_divisor, number)));
  for (; divisor * divisor < number; ++divisor) {
    if (number % divisor == 0) {
      if (divisor > max_divisor) return divisors;
      divisors.push_back(divisor);
      divisors.push_back(number / divisor);
    }
  }
  if (divisor * divisor == number && divisor <= max_divisor) {
    divisors.push_back(divisor);
  }
  return divisors;
}

}  // namespace

ArgbBuffer FuzzedParameters::ExtractArgbBuffer(uint32_t max_num_pixels,
                                               uint32_t max_width_or_height,
                                               WP2SampleFormat format) {
  assert(max_num_pixels <= kMaxBufferDimension);
  const size_t previous_bit_pos = bit_pos_;

  // Some bytes will be used for extracting the size and for the alignment.
  const size_t max_num_bits_for_size = (sizeof(uint32_t) + sizeof(uint8_t)) * 8;
  assert(bit_pos_ <= data_size_ * 8);

  // Check if there is enough data left.
  size_t num_available_bits = data_size_ * 8 - bit_pos_;
  if (num_available_bits < max_num_bits_for_size) return ArgbBuffer(format);
  num_available_bits -= max_num_bits_for_size;

  uint32_t num_pixels =
      std::min((uint32_t)num_available_bits / (8 * WP2FormatBpp(format)),
               max_num_pixels);
  if (num_pixels < 1) return ArgbBuffer(format);

  // Compute random valid width and height based on the number of pixels.
  std::vector<uint32_t> divisors =
      GetDivisors(num_pixels, /*max_divisor=*/max_width_or_height);
  if (divisors.empty()) {  // No valid divisor respecting 'max_width_or_height'.
    num_pixels = max_width_or_height;
    divisors = GetDivisors(num_pixels, /*max_divisor=*/max_width_or_height);
  }
  const uint32_t width = divisors[ExtractUInt32((uint32_t)divisors.size() - 1)];
  const uint32_t height = num_pixels / width;

  bit_pos_ = (bit_pos_ + 7u) & ~7u;  // Align.
  const uint32_t stride = width * WP2FormatBpp(format);
  assert(data_size_ - (bit_pos_ >> 3) >= height * stride);

  // Copy bytes and consider them as ARGB. Allocation may fail.
  ArgbBuffer buffer(format);
  if (buffer.Import(WP2_ARGB_32, width, height, data_ + (bit_pos_ >> 3),
                    stride) != WP2_STATUS_OK) {
    bit_pos_ = previous_bit_pos;  // Unread everything in case of failure.
    return ArgbBuffer(format);
  }
  bit_pos_ += (height * stride) * 8;
  return buffer;
}

//------------------------------------------------------------------------------

bool FuzzedParameters::TryExtractWP2Bitstream(Data* const encoded_data) {
  bit_pos_ = (bit_pos_ + 7u) & ~7u;  // Align.
  const size_t byte_pos = bit_pos_ / 8u;
  const size_t num_available_bytes = data_size_ - byte_pos;
  if (encoded_data->Resize(3u + num_available_bytes, /*keep_bytes=*/false) !=
      WP2_STATUS_OK) {
    return false;
  }
  encoded_data->bytes[0] = (uint8_t)((kSignature >> 0) & 0xffu);
  encoded_data->bytes[1] = (uint8_t)((kSignature >> 8) & 0xffu);
  encoded_data->bytes[2] = (uint8_t)((kSignature >> 16) & 0xffu);
  memcpy(encoded_data->bytes + 3, data_ + byte_pos, num_available_bytes);
  return true;
}

//------------------------------------------------------------------------------

void FuzzedParameters::ExtractSimpleConfig(uint32_t max_speed,
                                           EncoderConfig* const config) {
  config->quality = (float)ExtractUInt32(kMaxQuality);
  config->alpha_quality = (float)ExtractUInt32(kMaxQuality);
  config->speed = (int)ExtractUInt32(max_speed);
  config->decoding_orientation = (Orientation)ExtractUInt32(3);
  // config->create_preview = ExtractBool();  // Too slow.
  config->thread_level = (int)ExtractUInt32((1u << 31) - 1u);
}

void FuzzedParameters::ExtractAdvancedConfig(bool fast,
                                             EncoderConfig* const config) {
  config->transfer_function = (TransferFunction)ExtractUInt32(
      (uint32_t)WP2_TF_ARIB_STD_B67_HLG);  // Unimplemented.
  config->segment_id_mode = (EncoderConfig::SegmentIdMode)ExtractUInt32(2);
  for (uint32_t segment = 0; segment < kMaxNumSegments; ++segment) {
    config->segment_factors[segment] = (float)ExtractUInt32(kMaxLossyQuality);
  }
  config->tile_shape = (TileShape)ExtractUInt32((uint32_t)TILE_SHAPE_AUTO);

  constexpr PartitionMethod kFastEnoughPM[] = {
      MULTIPASS_PARTITIONING, ALL_4X4_PARTITIONING,   ALL_8X8_PARTITIONING,
      ALL_16X16_PARTITIONING, ALL_32X32_PARTITIONING, AUTO_PARTITIONING};
  config->partition_method = kFastEnoughPM[ExtractUInt32(
      sizeof(kFastEnoughPM) / sizeof(kFastEnoughPM[0]) - 1)];
  config->partition_set = (PartitionSet)ExtractUInt32(NUM_PARTITION_SETS - 1);
  config->partition_snapping = ExtractBool();

  config->csp_type = (Csp)ExtractUInt32(kNumCspTypes - 1);
  config->uv_mode = (EncoderConfig::UVMode)ExtractUInt32(
      (uint32_t)EncoderConfig::NumUVMode - 1);
  config->preprocessing = (int)ExtractUInt32(2);             // Unimplemented.
  config->preprocessing_strength = (int)ExtractUInt32(100);  // Unimplemented.
  config->use_delta_palette = ExtractBool();                 // Unimplemented.
  config->low_memory = ExtractBool();                        // Unimplemented.
  config->enable_alpha_filter = ExtractBool();

  if (fast) {
    config->sns = 0.f;
  } else {
    // config->pass = (int)ExtractUInt32(1, 10);  // Too slow.
    config->sns = (float)ExtractUInt32(100);
    config->error_diffusion = (int)ExtractUInt32(100);
    config->segments = (int)ExtractUInt32(1, kMaxNumSegments);
    config->use_random_matrix = ExtractBool();  // Experimental.
    config->store_grain = ExtractBool();        // Experimental.
    config->tune_perceptual = ExtractBool();    // Experimental.
  }
}

//------------------------------------------------------------------------------

void FuzzedParameters::ExtractConfig(DecoderConfig* const config) {
  config->thread_level = ExtractUInt32((1u << 31) - 1u);
  config->enable_deblocking_filter = ExtractBool();
  config->enable_directional_filter = ExtractBool();
  config->enable_restoration_filter = ExtractBool();
  config->grain_amplitude = (uint8_t)ExtractUInt32(100);
  config->incremental_mode = (DecoderConfig::IncrementalMode)ExtractUInt32(1);
}

//------------------------------------------------------------------------------

void SaveBeforeAfter(const ArgbBuffer& original, const ArgbBuffer& decompressed,
                     const std::string& path) {
  if (SaveImage(original, (path + "_original.png").c_str(),
                /*overwrite=*/true) != WP2_STATUS_OK ||
      SaveImage(decompressed, (path + "_decompressed.png").c_str(),
                /*overwrite=*/true) != WP2_STATUS_OK) {
    abort();
  }
}

//------------------------------------------------------------------------------

}  // namespace testing
}  // namespace WP2
