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
//  Utility functions used by the image decoders.
//

#include "./imageio_util.h"

#if defined(_WIN32)
#include <fcntl.h>  // for _O_BINARY
#include <io.h>     // for _setmode()

#include <cstdio>
#endif
#include <algorithm>
#include <cstdint>
#include <cstring>

#include "src/dsp/dsp.h"
#include "src/dsp/math.h"
#include "src/utils/utils.h"

namespace WP2 {

// useful struct for RAII
class FileCloser {
 public:
  explicit FileCloser(FILE* const file) : file_(file) {}
  ~FileCloser() {
    if (file_ != nullptr && file_ != stdin && file_ != stdout) fclose(file_);
  }
 private:
  FILE* file_;
};

// -----------------------------------------------------------------------------
// File I/O

FILE* IoUtilSetBinaryMode(FILE* const file) {
#if defined(_WIN32)
  if (_setmode(_fileno(file), _O_BINARY) == -1) {
    fprintf(stderr, "Failed to reopen file in O_BINARY mode.\n");
    return nullptr;
  }
#endif
  return file;
}

WP2Status IoUtilReadFromStdin(Data* const data, size_t max_num_bytes) {
  static constexpr size_t kBlockSize = 16384;  // default initial size
  const size_t block_size =
      (max_num_bytes != 0) ? std::min(max_num_bytes, kBlockSize) : kBlockSize;

  WP2_CHECK_OK(data != nullptr, WP2_STATUS_NULL_PARAMETER);
  WP2_CHECK_OK(IoUtilSetBinaryMode(stdin) != nullptr,
               WP2_STATUS_INVALID_PARAMETER);

  size_t num_read_bytes = 0, max_size = 0;
  while (!feof(stdin)) {
    // We double the buffer size each time and read as much as possible.
    size_t extra_size = (max_size == 0) ? block_size : max_size;
    if (max_num_bytes != 0) {
      assert(max_size <= max_num_bytes);
      extra_size = std::min(max_num_bytes - max_size, extra_size);
    }

    WP2_CHECK_STATUS(data->Resize(max_size + extra_size,
                                  /*keep_bytes=*/(num_read_bytes > 0)));

    num_read_bytes += fread(data->bytes + num_read_bytes, 1, extra_size, stdin);
    max_size += extra_size;
    assert(num_read_bytes <= max_size);
    if (num_read_bytes < max_size) break;
    if (max_num_bytes != 0) {
      assert(num_read_bytes <= max_num_bytes);
      if (num_read_bytes >= max_num_bytes) break;
    }
  }
  WP2_CHECK_OK(ferror(stdin) == 0, WP2_STATUS_BAD_READ);
  // Trim unnecessary bytes.
  WP2_CHECK_STATUS(data->Resize(num_read_bytes, /*keep_bytes=*/true));
  return WP2_STATUS_OK;
}

static WP2Status IoUtilFileSize(FILE* const in, size_t* const size) {
  assert(in != nullptr && size != nullptr);
  fseek(in, 0, SEEK_END);
  const long fsize = ftell(in);  // NOLINT (long)
  fseek(in, 0, SEEK_SET);
  if (fsize < 0) return WP2_STATUS_BAD_READ;
  if ((long)(size_t)fsize != fsize) return WP2_STATUS_BAD_READ;  // NOLINT
  *size = (size_t)fsize;
  return WP2_STATUS_OK;
}

WP2Status IoUtilFileSize(const char* const file_path, size_t* const size) {
  *size = 0;
  WP2_CHECK_OK(file_path != nullptr && size != nullptr,
               WP2_STATUS_NULL_PARAMETER);
  FILE* const in = fopen(file_path, "rb");
  WP2_CHECK_OK(in != nullptr, WP2_STATUS_INVALID_PARAMETER);
  FileCloser closer(in);
  return IoUtilFileSize(in, size);
}

WP2Status IoUtilReadFile(const char* const file_path, Data* const data,
                         size_t max_num_bytes) {
  if ((file_path == nullptr) || !strcmp(file_path, "-")) {
    return IoUtilReadFromStdin(data, max_num_bytes);
  }

  WP2_CHECK_OK(data != nullptr, WP2_STATUS_NULL_PARAMETER);

  FILE* const in = fopen(file_path, "rb");
  WP2_CHECK_OK(in != nullptr, WP2_STATUS_INVALID_PARAMETER);
  FileCloser closer(in);
  size_t file_size;
  WP2_CHECK_STATUS(IoUtilFileSize(in, &file_size));
  if (max_num_bytes != 0) file_size = std::min(file_size, max_num_bytes);
  WP2_CHECK_STATUS(data->Resize(file_size, /*keep_bytes=*/false));

  static constexpr size_t kNumReadElems = 1;
  const bool success =
      (file_size == 0 ||
       fread(data->bytes, data->size, kNumReadElems, in) == kNumReadElems);
  if (!success) {
    data->Clear();
    return WP2_STATUS_BAD_READ;
  }
  return WP2_STATUS_OK;
}

WP2Status IoUtilWriteFile(const uint8_t* const data, size_t data_size,
                          const char* const file_path, bool overwrite,
                          size_t* const file_size) {
  const int to_stdout = (file_path == nullptr) || !strcmp(file_path, "-");
  if (to_stdout) {  // std::ftell() does not work with pipes.
    WP2_CHECK_OK(file_size == nullptr, WP2_STATUS_INVALID_PARAMETER);
  }

  WP2_CHECK_OK(data != nullptr, WP2_STATUS_NULL_PARAMETER);
  const char* const mode = (overwrite ? "wb" : "wxb");
  FILE* const out =
      to_stdout ? IoUtilSetBinaryMode(stdout) : fopen(file_path, mode);
  WP2_CHECK_OK(out != nullptr, (!to_stdout && FileExists(file_path))
                                   ? WP2_STATUS_BAD_WRITE
                                   : WP2_STATUS_INVALID_PARAMETER);
  FileCloser closer(out);

  static constexpr size_t kNumWrittenElems = 1;
  const bool success =
      (fwrite(data, data_size, kNumWrittenElems, out) == kNumWrittenElems);
  if (file_size != nullptr) *file_size = std::ftell(out);
  return success ? WP2_STATUS_OK : WP2_STATUS_BAD_WRITE;
}

bool FileExists(const char* const file_path) {
  FILE* const out = fopen(file_path, "ab");
  if (out == nullptr) return false;
  fclose(out);
  return true;
}

//------------------------------------------------------------------------------
// string variant

WP2Status IoUtilReadFile(const char* const file_path,
                         std::string* const data_string) {
  if ((file_path == nullptr) || !strcmp(file_path, "-")) {
    return IoUtilReadFromStdin(data_string);
  }

  WP2_CHECK_OK(data_string != nullptr, WP2_STATUS_NULL_PARAMETER);
  data_string->clear();

  FILE* const in = fopen(file_path, "rb");
  WP2_CHECK_OK(in != nullptr, WP2_STATUS_INVALID_PARAMETER);
  FileCloser closer(in);

  size_t data_size;
  WP2_CHECK_STATUS(IoUtilFileSize(in, &data_size));
  WP2_CHECK_ALLOC_OK(data_size <= WP2_MAX_ALLOCABLE_MEMORY);
  data_string->resize(data_size);
  WP2_CHECK_ALLOC_OK(data_string->size() >= data_size);

  void* const bytes = &(*data_string)[0];  // string::data() is const in C++11
  static constexpr size_t kNumReadElems = 1;
  const bool success =
      (fread(bytes, data_string->size(), kNumReadElems, in) == kNumReadElems);
  return success ? WP2_STATUS_OK : WP2_STATUS_BAD_READ;
}

WP2Status IoUtilReadFromStdin(std::string* const data_string) {
  static constexpr size_t kBlockSize = 16384;  // default initial size

  WP2_CHECK_OK(data_string != nullptr, WP2_STATUS_NULL_PARAMETER);
  WP2_CHECK_OK(IoUtilSetBinaryMode(stdin) != nullptr,
               WP2_STATUS_INVALID_PARAMETER);
  data_string->clear();

  while (!feof(stdin)) {
    const size_t current_size = data_string->size();
    const size_t extra_size = current_size ? kBlockSize : current_size;
    WP2_CHECK_ALLOC_OK(current_size + extra_size <= WP2_MAX_ALLOCABLE_MEMORY);
    data_string->resize(current_size + extra_size);
    WP2_CHECK_ALLOC_OK(data_string->size() >= current_size + extra_size);

    void* const bytes = &(*data_string)[current_size];
    data_string->resize(current_size + fread(bytes, 1, extra_size, stdin));
    if (data_string->size() < current_size + extra_size) break;
  }
  WP2_CHECK_OK(ferror(stdin) == 0, WP2_STATUS_BAD_READ);
  return WP2_STATUS_OK;
}

WP2Status IoUtilWriteFile(const std::string& data, const char* const file_path,
                          bool overwrite) {
  return IoUtilWriteFile((const uint8_t*)data.data(), data.size(), file_path,
                         overwrite);
}

//------------------------------------------------------------------------------

#if defined(WP2_HAVE_WEBP)
WP2Status ConvertToWebPPicture(const WP2::ArgbBuffer& in,
                               WebPPicture* const out, bool use_yuv) {
  WP2_CHECK_OK(WebPPictureInit(out), WP2_STATUS_VERSION_MISMATCH);
  out->use_argb = 1;
  out->width = in.width;
  out->height = in.height;
  WP2_CHECK_OK(WebPPictureAlloc(out), WP2_STATUS_OUT_OF_MEMORY);
  // TODO(skal): sharp-yuv conversion ?
  // Convert from Argb_32 to WebP's ARGB
  WP2ArgbConverterInit();
  for (int y = 0; y < out->height; ++y) {
    // BGRA order for conversion from uint8_t[4] to uint32_t[1] endianness.
    WP2ArgbConvertTo[WP2_BGRA_32]((const uint8_t*)in.GetRow(y), out->width,
                                  (uint8_t*)(out->argb + y * out->argb_stride));
  }
  if (use_yuv) {
    WP2_CHECK_OK(WebPPictureSharpARGBToYUVA(out), WP2_STATUS_OUT_OF_MEMORY);
  }
  return WP2_STATUS_OK;
}

WP2Status ConvertFromWebPPicture(const WebPPicture& in,
                                 WP2::ArgbBuffer* const out) {
  WP2_CHECK_OK(out != nullptr, WP2_STATUS_NULL_PARAMETER);
  // TODO(skal): potential endian problem
  return out->Import(WP2_bgrA_32, in.width, in.height, (const uint8_t*)in.argb,
                     in.argb_stride * sizeof(*in.argb));
}
#endif  // WP2_HAVE_WEBP

//------------------------------------------------------------------------------

}  // namespace WP2
