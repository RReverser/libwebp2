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

#include <cstdint>
#include <cstdio>

#ifdef WP2_HAVE_PNG
#include <png.h>
#include <setjmp.h>  // note: this must be included *after* png.h
#endif  // WP2_HAVE_PNG

#ifdef HAVE_WINCODEC_H
#ifdef __MINGW32__
#define INITGUID  // Without this GUIDs are declared extern and fail to link
#endif
#define CINTERFACE
#define COBJMACROS
#define _WIN32_IE 0x500  // Workaround bug in shlwapi.h when compiling C++
                         // code with COBJMACROS.
#include <ole2.h>  // CreateStreamOnHGlobal()
#include <shlwapi.h>
#include <windows.h>
#include <wincodec.h>
#endif  // HAVE_WINCODEC_H

#include "./image_enc.h"
#include "src/dsp/dsp.h"
#include "src/utils/utils.h"
#include "src/utils/vector.h"

namespace WP2 {

#ifdef HAVE_WINCODEC_H

namespace {

#define IFS(fn)                                                   \
  do {                                                            \
    if (SUCCEEDED(hr)) {                                          \
      hr = (fn);                                                  \
      if (FAILED(hr)) fprintf(stderr, #fn " failed %08lx\n", hr); \
    }                                                             \
  } while (0)

#ifdef __cplusplus
#define MAKE_REFGUID(x) (x)
#else
#define MAKE_REFGUID(x) &(x)
#endif

HRESULT CreateOutputStream(const char* file_path, bool write_to_mem,
                           bool overwrite, IStream** stream) {
  HRESULT hr = S_OK;
  if (write_to_mem) {
    // Output to a memory buffer. This is freed when 'stream' is released.
    IFS(CreateStreamOnHGlobal(nullptr, TRUE, stream));
  } else {
    const DWORD create_mode = (overwrite ? STGM_CREATE : STGM_FAILIFTHERE);
    IFS(SHCreateStreamOnFileA(file_path, STGM_WRITE | create_mode, stream));
  }
  if (FAILED(hr)) {
    fprintf(stderr, "Error opening output file %s (%08lx)\n", file_path, hr);
  }
  return hr;
}

HRESULT WriteUsingWIC(const char* file_path, bool use_stdout, bool overwrite,
                      REFGUID container_guid, uint8_t* rgb, int stride,
                      uint32_t width, uint32_t height, int has_alpha,
                      size_t* const file_size) {
  HRESULT hr = S_OK;
  IWICImagingFactory* factory = nullptr;
  IWICBitmapFrameEncode* frame = nullptr;
  IWICBitmapEncoder* encoder = nullptr;
  IStream* stream = nullptr;
  WICPixelFormatGUID pixel_format =
      has_alpha ? GUID_WICPixelFormat32bppBGRA : GUID_WICPixelFormat24bppBGR;

  if (file_path == nullptr || rgb == nullptr) return E_INVALIDARG;

  IFS(CoInitialize(nullptr));
  IFS(CoCreateInstance(
      MAKE_REFGUID(CLSID_WICImagingFactory), nullptr, CLSCTX_INPROC_SERVER,
      MAKE_REFGUID(IID_IWICImagingFactory), (LPVOID*)&factory));
  if (hr == REGDB_E_CLASSNOTREG) {
    fprintf(stderr,
            "Couldn't access Windows Imaging Component (are you running "
            "Windows XP SP3 or newer?). PNG support not available. "
            "Use -ppm or -pgm for available PPM and PGM formats.\n");
  }
  IFS(CreateOutputStream(file_path, use_stdout, overwrite, &stream));
  IFS(IWICImagingFactory_CreateEncoder(factory, container_guid, nullptr,
                                       &encoder));
  IFS(IWICBitmapEncoder_Initialize(encoder, stream, WICBitmapEncoderNoCache));
  IFS(IWICBitmapEncoder_CreateNewFrame(encoder, &frame, nullptr));
  IFS(IWICBitmapFrameEncode_Initialize(frame, nullptr));
  IFS(IWICBitmapFrameEncode_SetSize(frame, width, height));
  IFS(IWICBitmapFrameEncode_SetPixelFormat(frame, &pixel_format));
  IFS(IWICBitmapFrameEncode_WritePixels(frame, height, stride, height * stride,
                                        rgb));
  IFS(IWICBitmapFrameEncode_Commit(frame));
  IFS(IWICBitmapEncoder_Commit(encoder));

  if (SUCCEEDED(hr) && use_stdout) {
    HGLOBAL image;
    IFS(GetHGlobalFromStream(stream, &image));
    if (SUCCEEDED(hr)) {
      HANDLE std_output = GetStdHandle(STD_OUTPUT_HANDLE);
      DWORD mode;
      const BOOL update_mode = GetConsoleMode(std_output, &mode);
      const void* const image_mem = GlobalLock(image);
      DWORD bytes_written = 0;

      // Clear output processing if necessary, then output the image.
      if (update_mode) SetConsoleMode(std_output, 0);
      if (!WriteFile(std_output, image_mem, (DWORD)GlobalSize(image),
                     &bytes_written, nullptr) ||
          bytes_written != GlobalSize(image)) {
        hr = E_FAIL;
      } else if (file_size != nullptr) {
        *file_size = (size_t)bytes_written;
      }
      if (update_mode) SetConsoleMode(std_output, mode);
      GlobalUnlock(image);
    }
  }

  if (frame != nullptr) IUnknown_Release(frame);
  if (encoder != nullptr) IUnknown_Release(encoder);
  if (factory != nullptr) IUnknown_Release(factory);
  if (stream != nullptr) IUnknown_Release(stream);
  return hr;
}

}  // namespace

WP2Status WritePNG(const ArgbBuffer& buffer, FILE*, const char* const file_path,
                   bool use_stdout, bool overwrite, size_t* const file_size) {
  WP2_CHECK_OK(file_path != nullptr, WP2_STATUS_NULL_PARAMETER);
  WP2_CHECK_OK(buffer.format == WP2_Argb_32, WP2_STATUS_UNSUPPORTED_FEATURE);

  const uint32_t width = buffer.width;
  const uint32_t height = buffer.height;
  const uint8_t* const Argb = buffer.Argb;
  const uint32_t stride = buffer.stride;
  const int ok = SUCCEEDED(WriteUsingWIC(
      file_path, use_stdout, overwrite, MAKE_REFGUID(GUID_ContainerFormatPng),
      rgb, stride, width, height, 1, file_size));
  return ok ? WP2_STATUS_OK : WP2_STATUS_BAD_WRITE;
}

#elif defined(WP2_HAVE_PNG)  // !HAVE_WINCODEC_H

namespace {

void PNGAPI PNGErrorFunction(png_structp png, png_const_charp dummy) {
  (void)dummy;  // remove variable-unused warning
  longjmp(png_jmpbuf(png), 1);
}

}  // namespace

WP2Status WritePNG(const ArgbBuffer& buffer, FILE* fout,
                   const char*, bool, bool, size_t*) {
  WP2_CHECK_OK(fout != nullptr, WP2_STATUS_NULL_PARAMETER);
  WP2_CHECK_OK(buffer.format == WP2_Argb_32, WP2_STATUS_UNSUPPORTED_FEATURE);

  const uint32_t width = buffer.width;
  const uint32_t height = buffer.height;
  const bool has_alpha = buffer.HasTransparency();
  volatile png_bytep row = (png_bytep)buffer.GetRow(0);
  Vector_u8 ARGB;
  WP2_CHECK_ALLOC_OK(ARGB.resize(width * (has_alpha ? 4u : 3u)));
  const uint32_t stride = buffer.stride;
  volatile png_structp png;
  volatile png_infop info;
  png_uint_32 y;
  WP2Status status;

  png = png_create_write_struct(PNG_LIBPNG_VER_STRING, nullptr,
                                PNGErrorFunction, nullptr);
  if (png == nullptr) return WP2_STATUS_OUT_OF_MEMORY;

  info = png_create_info_struct(png);
  if (info == nullptr) {
    png_destroy_write_struct((png_structpp)&png, nullptr);
    return WP2_STATUS_OUT_OF_MEMORY;
  }
  if (setjmp(png_jmpbuf(png))) {
    png_destroy_write_struct((png_structpp)&png, (png_infopp)&info);
    status = WP2_STATUS_BAD_WRITE;
    goto End;
  }
  png_init_io(png, fout);
  png_set_IHDR(png, info, width, height, 8,
               has_alpha ? PNG_COLOR_TYPE_RGBA : PNG_COLOR_TYPE_RGB,
               PNG_INTERLACE_NONE, PNG_COMPRESSION_TYPE_DEFAULT,
               PNG_FILTER_TYPE_DEFAULT);
  png_write_info(png, info);

  WP2ArgbConverterInit();

  for (y = 0; y < height; ++y) {
    WP2ArgbConvertTo[has_alpha ? WP2_RGBA_32 : WP2_RGB_24](row, width,
                                                           ARGB.data());
    png_bytep ptr[1] = {ARGB.data()};
    png_write_rows(png, ptr, 1);
    row += stride;
  }
  png_write_end(png, info);
  status = WP2_STATUS_OK;

 End:
  png_destroy_write_struct((png_structpp)&png, (png_infopp)&info);
  return status;
}

#else  // !HAVE_WINCODEC_H && !WP2_HAVE_PNG

WP2Status WritePNG(const ArgbBuffer&, FILE* fout, const char*, bool, bool,
                   size_t*) {
  WP2_CHECK_OK(fout != nullptr, WP2_STATUS_NULL_PARAMETER);

  fprintf(stderr, "PNG support not compiled. Please install the libpng "
                  "development package before building.\n");
  fprintf(stderr, "You can run with -ppm flag to decode in PPM format.\n");
  return WP2_STATUS_BAD_WRITE;
}

#endif

}  // namespace WP2
