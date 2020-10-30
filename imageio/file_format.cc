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
//  File extension utils.
//
// Author: Yannis Guyon (yguyon@google.com)

#include "./file_format.h"

#include <cctype>
#include <cstring>

#include "src/utils/utils.h"

namespace WP2 {

const char* GetFileFormatStr(FileFormat file_format) {
  static constexpr const char* const kFileFormatStr[] = {
      "PNG", "JPEG", "WEBP", "WP2",     "TIFF",    "GIF",  "BMP",
      "PAM", "PGM",  "PPM",  "Y4M_420", "Y4M_444", "AUTO", "UNSUPPORTED"};
  STATIC_ASSERT_ARRAY_SIZE(kFileFormatStr, (int)FileFormat::UNSUPPORTED + 1);
  return kFileFormatStr[(int)file_format];
}

namespace {

constexpr struct FormatExtension {
  FileFormat format;
  const char* extension;
} kFormatExtension[] = {
    // Top row for the same extension (or format) represent the default format
    // (or extension).
    {FileFormat::PNG, "png"},
    {FileFormat::JPEG, "jpg"},
    {FileFormat::JPEG, "jpeg"},
    {FileFormat::BMP, "bmp"},
    {FileFormat::PAM, "pam"},
    {FileFormat::PGM, "pgm"},
    {FileFormat::PPM, "ppm"},
    {FileFormat::TIFF, "tiff"},
    {FileFormat::TIFF, "tif"},
    {FileFormat::GIF, "gif"},
    {FileFormat::WEBP, "webp"},
    {FileFormat::WP2, "wp2"},
    {FileFormat::Y4M_420, "y4m"},
    {FileFormat::Y4M_444, "y4m"},
};

bool IsExtension(const char* const str,
                 const char* const lower_case_extension) {
  for (size_t i = 0; std::tolower(str[i]) == lower_case_extension[i]; ++i) {
    if (lower_case_extension[i] == '\0') return true;
  }
  return false;
}

#ifdef _WIN32
constexpr char kPathSeparator = (char)'\\';
#else
constexpr char kPathSeparator = (char)'/';
#endif  // _WIN32

}  // namespace

FileFormat GetFormatFromExtension(const char* const file_path) {
  const char* const last_separator = std::strrchr(file_path, kPathSeparator);
  const char* const file_name =
      (last_separator != nullptr) ? (last_separator + 1) : file_path;
  const char* const last_dot = std::strrchr(file_name, '.');
  if (last_dot == nullptr) return FileFormat::UNSUPPORTED;
  for (const FormatExtension& f : kFormatExtension) {
    if (IsExtension(last_dot + 1, f.extension)) return f.format;
  }
  return FileFormat::UNSUPPORTED;
}

const char* GetExtensionFromFormat(FileFormat format) {
  for (const FormatExtension& f : kFormatExtension) {
    if (format == f.format) return f.extension;
  }
  if (format == FileFormat::AUTO) return "png";
  return nullptr;  // UNSUPPORTED
}

bool IsCustomColorSpace(FileFormat file_format) {
  return (file_format == FileFormat::Y4M_420 ||
          file_format == FileFormat::Y4M_444);
}

bool IsChromaSubsampled(FileFormat file_format) {
  return (file_format == FileFormat::Y4M_420);
}

}  // namespace WP2
