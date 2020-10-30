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
// Generic image-type guessing.
//
// Author: Skal (pascal.massimino@gmail.com)

#include "./image_dec.h"

#include <cstddef>
#include <cstdint>

#include "./anim_image_dec.h"

namespace WP2 {

static constexpr uint32_t kFormatNumBytes = 12;

static inline uint32_t GetBE16(const uint8_t buf[]) {
  return ((uint32_t)buf[0] << 8) | buf[1];
}

static inline uint32_t GetBE32(const uint8_t buf[]) {
  return ((uint32_t)buf[0] << 24) | (buf[1] << 16) | (buf[2] << 8) | buf[3];
}

FileFormat GuessImageFormat(const uint8_t* const data, size_t data_size) {
  if (data != nullptr && data_size >= kFormatNumBytes) {
    const uint32_t magic1 = GetBE32(data + 0);
    const uint32_t magic2 = GetBE16(data + 4);
    const uint32_t magic3 = GetBE32(data + 8);
    if (magic1 == 0x89504E47U) {
      return FileFormat::PNG;
    } else if (magic1 >= 0xFFD8FF00U && magic1 <= 0xFFD8FFFFU) {
      return FileFormat::JPEG;
    } else if (magic1 == 0x49492A00 || magic1 == 0x4D4D002A) {
      return FileFormat::TIFF;
    } else if (magic1 == 0x47494638 && (magic2 == 0x3761 || magic2 == 0x3961)) {
      return FileFormat::GIF;
    } else if (magic1 == 0x52494646 && magic3 == 0x57454250) {
      return FileFormat::WEBP;
    } else if ((magic1 >> 8) == 0xf4ff6f) {
      return FileFormat::WP2;
    } else if ((magic1 >> 16) == 0x424d) {
      return FileFormat::BMP;
    } else if (((magic1 >> 24) & 0xff) == 'P') {
      const int type = (magic1 >> 16) & 0xff;
      // we only support 'P5 -> P7' for now.
      if (type == '5') return FileFormat::PGM;
      if (type == '6') return FileFormat::PPM;
      if (type == '7') return FileFormat::PAM;
    } else if (magic1 == GetBE32((const uint8_t*)"YUV4")) {
      // Try to find the color space tag in available 'data'.
      const char* const header = reinterpret_cast<const char*>(data);
      const size_t header_size = std::min(data_size, (size_t)64);
      if (strnstr(header, header_size, " C")) {
        if (strnstr(header, header_size, " C420")) return FileFormat::Y4M_420;
        if (strnstr(header, header_size, " C444")) return FileFormat::Y4M_444;
      } else if (data_size < 64) {
        // No color tag yet. Assume 420, there are probably more bytes to come.
        return FileFormat::Y4M_420;
      }
    }
  }
  return FileFormat::UNSUPPORTED;
}

WP2Status ReadImage(const uint8_t* const data, size_t data_size,
                    ArgbBuffer* const buffer, FileFormat format,
                    LogLevel log_level) {
  ImageReader image_reader(data, data_size, buffer, format, log_level);
  return image_reader.ReadFrame();
}

WP2Status ReadImage(const char* const file_path, ArgbBuffer* const buffer,
                    size_t* const file_size, FileFormat format,
                    LogLevel log_level) {
  Data data;
  WP2_CHECK_STATUS(IoUtilReadFile(file_path, &data));
  if (file_size != nullptr) *file_size = data.size;
  return ReadImage(data.bytes, data.size, buffer, format, log_level);
}

}  // namespace WP2
