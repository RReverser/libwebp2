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
// Image savers

#include "./image_enc.h"

#include <cstring>

#include "./imageio_util.h"
#include "src/utils/utils.h"

namespace WP2 {

//------------------------------------------------------------------------------

WP2Status SaveImage(const ArgbBuffer& buffer, const char* const file_path,
                    bool overwrite, FileFormat format,
                    size_t* const file_size) {
  WP2_CHECK_OK(file_path != nullptr, WP2_STATUS_NULL_PARAMETER);
  WP2_CHECK_OK(!buffer.IsEmpty(), WP2_STATUS_INVALID_PARAMETER);

  FILE* fout = nullptr;
  const bool use_stdout = (strcmp(file_path, "-") == 0);
  if (use_stdout) {  // std::ftell() does not work with pipes.
    WP2_CHECK_OK(file_size == nullptr, WP2_STATUS_INVALID_PARAMETER);
  }

#ifdef HAVE_WINCODEC_H
  if (format != PNG) {
#else
  {
#endif
    const char* const mode = (overwrite ? "wb" : "wxb");
    fout =
        use_stdout ? IoUtilSetBinaryMode(stdout) : std::fopen(file_path, mode);
    WP2_CHECK_OK(fout != nullptr, (!use_stdout && FileExists(file_path))
                                      ? WP2_STATUS_BAD_WRITE
                                      : WP2_STATUS_INVALID_PARAMETER);
  }

  if (format == FileFormat::AUTO) {
    format = use_stdout ? FileFormat::PNG : GetFormatFromExtension(file_path);
  }

  WP2Status status;
  if (format == FileFormat::PNG) {
    status =
        WritePNG(buffer, fout, file_path, use_stdout, overwrite, file_size);
  } else if (format == FileFormat::TIFF) {
    status = WriteTIFF(buffer, fout);
  } else if (format == FileFormat::BMP) {
    status = WriteBMP(buffer, fout);
  } else if (format == FileFormat::PGM) {
    status = Write16bAsPGM(buffer, fout);
  } else if (format == FileFormat::PPM) {
    status = WritePPM(buffer, fout);
  } else if (format == FileFormat::PAM) {
    status = WritePAM(buffer, fout);
  } else {
    status = WP2_STATUS_UNSUPPORTED_FEATURE;
  }

  if (fout != nullptr && file_size != nullptr) *file_size = std::ftell(fout);
  if (fout != nullptr && fout != stdout) std::fclose(fout);
  return status;
}

//------------------------------------------------------------------------------

}  // namespace WP2
