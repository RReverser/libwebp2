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
// PNG decode.
//
// Author: Skal (pascal.massimino@gmail.com)

#include <cstdio>

#ifdef HAVE_CONFIG_H
#include "src/wp2/config.h"
#endif

#include "./anim_image_dec.h"

#ifdef WP2_HAVE_PNG
#include <png.h>
#include <setjmp.h>  // note: this must be included *after* png.h
#include <cstdlib>
#include <cstring>

#include "./imageio_util.h"
#include "src/utils/utils.h"
#include "src/wp2/base.h"

namespace {

//------------------------------------------------------------------------------

void PNGAPI ErrorFunctionNoLog(png_structp png, png_const_charp error) {
  (void)error;
  longjmp(png_jmpbuf(png), 1);
}

void PNGAPI ErrorFunction(png_structp png, png_const_charp error) {
  if (error != NULL) {
    fprintf(stderr, "libpng error: %s\n", error);
  } else {
    fprintf(stderr, "unknown libpng error\n");
  }
  ErrorFunctionNoLog(png, error);
}

//------------------------------------------------------------------------------

// Converts the NULL terminated 'hexstring' which contains 2-byte character
// representations of hex values to raw data.
// 'hexstring' may contain values consisting of [A-F][a-f][0-9] in pairs,
// e.g., 7af2..., separated by any number of newlines.
// 'expected_length' is the anticipated processed size.
// On success the 'payload' is filled with its length equivalent to
// 'expected_length' and WP2_STATUS_OK is returned. An error is returned if the
// allocation failed, if the processed length is less than 'expected_length' or
// if any character aside from those above is encountered.
WP2Status HexStringToBytes(const char* const hexstring, size_t expected_length,
                           WP2::Data* const payload) {
  WP2_CHECK_STATUS(payload->Resize(expected_length, /*keep_bytes=*/false));

  const char* src = hexstring;
  size_t actual_length = 0;
  for (uint8_t* dst = payload->bytes;
       actual_length < expected_length && *src != '\0'; ++src) {
    if (*src == '\n') continue;
    const char val[3]{*src++, *src, '\0'};  // One hex byte string ("0F").
    char* end;
    *dst++ = (uint8_t)std::strtol(val, &end, 16);  // strtol NOLINT
    WP2_CHECK_OK(end == val + 2, WP2_STATUS_BITSTREAM_ERROR);
    ++actual_length;
  }

  WP2_CHECK_OK(actual_length == expected_length, WP2_STATUS_BITSTREAM_ERROR);
  return WP2_STATUS_OK;
}

WP2Status ProcessRawProfile(const char* const profile, size_t profile_len,
                            WP2::Data* const payload, WP2::LogLevel log_level) {
  WP2_CHECK_OK(profile != nullptr, WP2_STATUS_NULL_PARAMETER);
  WP2_CHECK_OK(profile_len > 0, WP2_STATUS_NOT_ENOUGH_DATA);

  // ImageMagick formats 'raw profiles' as
  // '\n<name>\n<length>(%8lu)\n<hex payload>\n'.
  const char* src = profile;
  if (*src != '\n') {
    if (log_level >= WP2::LogLevel::DEFAULT) {
      fprintf(stderr, "Malformed raw profile, expected '\\n' got '\\x%.2X'\n",
              *src);
    }
    return WP2_STATUS_BITSTREAM_ERROR;
  }
  ++src;
  // Skip the profile name and extract the length.
  while (*src != '\0' && *src++ != '\n') {
  }
  char* end;
  const size_t expected_length = (size_t)std::strtol(src, &end, 10);  // NOLINT
  if (*end != '\n') {
    if (log_level >= WP2::LogLevel::DEFAULT) {
      fprintf(stderr, "Malformed raw profile, expected '\\n' got '\\x%.2X'\n",
              *end);
    }
    return WP2_STATUS_BITSTREAM_ERROR;
  }
  ++end;

  // 'end' now points to the profile payload.
  WP2_CHECK_STATUS(HexStringToBytes(end, expected_length, payload));
  return WP2_STATUS_OK;
}

WP2Status ProcessProfile(const char* const profile, size_t profile_len,
                         WP2::Data* const payload) {
  WP2_CHECK_OK(profile != nullptr, WP2_STATUS_NULL_PARAMETER);
  WP2_CHECK_OK(profile_len > 0, WP2_STATUS_NOT_ENOUGH_DATA);

  return payload->CopyFrom((const uint8_t*)profile, profile_len);
}

//------------------------------------------------------------------------------

png_size_t GetTextLength(png_const_textp const text_chunk) {
#ifdef PNG_iTXt_SUPPORTED
  if (text_chunk->compression == PNG_ITXT_COMPRESSION_NONE ||
      text_chunk->compression == PNG_ITXT_COMPRESSION_zTXt) {
    return text_chunk->itxt_length;
  }
#endif
  return text_chunk->text_length;
}

// Returns the associated 'metadata' payload and if it 'is_raw_profile' based on
// the PNG text chunk 'key'. Returns null if unhandled.
WP2::Data* GetPayload(const char* const key, WP2::Metadata* const metadata,
                      bool* const is_raw_profile) {
  // www.sno.phy.queensu.ca/~phil/exiftool/TagNames/PNG.html#TextualData
  // See also: ExifTool on CPAN. Exiftool puts exif data in APP1 chunk, too.
  if (!std::strcmp(key, "Raw profile type exif") ||
      !std::strcmp(key, "Raw profile type APP1")) {
    *is_raw_profile = true;
    return &metadata->exif;
  }
  if (!std::strcmp(key, "Raw profile type xmp")) {
    *is_raw_profile = true;
    return &metadata->xmp;
  }
  // XMP Specification Part 3, Section 3 #PNG
  if (!std::strcmp(key, "XML:com.adobe.xmp")) {
    *is_raw_profile = false;
    return &metadata->xmp;
  }
  return nullptr;
}

// Looks for metadata in 'info' and extracts it to 'metadata'.
// Returns WP2_STATUS_OK on success.
WP2Status ExtractMetadataFromPNG(png_const_structrp png, png_infop const info,
                                 WP2::Metadata* const metadata,
                                 WP2::LogLevel log_level) {
  // Look for EXIF / XMP metadata.
  png_textp text_chunk = NULL;
  const int num_text_chunks = png_get_text(png, info, &text_chunk, NULL);
  for (int i = 0; i < num_text_chunks; ++i, ++text_chunk) {
    bool is_raw_profile = false;
    WP2::Data* const payload =
        GetPayload((const char*)text_chunk->key, metadata, &is_raw_profile);

    if (payload == nullptr) {
      if (log_level >= WP2::LogLevel::VERBOSE) {
        fprintf(stderr, "Ignoring unhandled '%s'\n", text_chunk->key);
      }
    } else if (!payload->IsEmpty()) {
      if (log_level >= WP2::LogLevel::VERBOSE) {
        fprintf(stderr, "Ignoring additional '%s'\n", text_chunk->key);
      }
    } else {
      const png_size_t text_length = GetTextLength(text_chunk);
      WP2Status status;
      if (is_raw_profile) {
        status = ProcessRawProfile(text_chunk->text, text_length, payload,
                                   log_level);
      } else {
        status = ProcessProfile(text_chunk->text, text_length, payload);
      }
      if (status != WP2_STATUS_OK) {
        if (log_level >= WP2::LogLevel::DEFAULT) {
          fprintf(stderr, "Failed to process: '%s'\n", text_chunk->key);
          fprintf(stderr, "Error extracting PNG metadata!\n");
        }
        return status;
      }
    }
  }

  if (metadata->iccp.IsEmpty()) {
    // Look for an ICC profile.
#if ((PNG_LIBPNG_VER_MAJOR << 8) | PNG_LIBPNG_VER_MINOR << 0) < \
    ((1 << 8) | (5 << 0))
    png_charp profile;
#else  // >= libpng 1.5.0
    png_bytep profile;
#endif
    png_uint_32 length;

    png_charp name;
    int comp_type;
    if (png_get_iCCP(png, info, &name, &comp_type, &profile, &length) ==
        PNG_INFO_iCCP) {
      const WP2Status status =
          ProcessProfile((const char*)profile, length, &metadata->iccp);
      if (status != WP2_STATUS_OK) {
        if (log_level >= WP2::LogLevel::DEFAULT) {
          fprintf(stderr, "Failed to process: '%s'\n", text_chunk->key);
          fprintf(stderr, "Error extracting PNG metadata!\n");
        }
        return status;
      }
    }
  }
  return WP2_STATUS_OK;
}

//------------------------------------------------------------------------------

struct PNGReadContext {
  const uint8_t* data;
  size_t data_size;
  png_size_t offset;
};

void ReadFunc(png_structp png_ptr, png_bytep data, png_size_t length) {
  PNGReadContext* const ctx = (PNGReadContext*)png_get_io_ptr(png_ptr);
  if (ctx->data_size - ctx->offset < length) {
    png_error(png_ptr, "ReadFunc: invalid read length (overflow)!");
  }
  std::memcpy(data, ctx->data + ctx->offset, length);
  ctx->offset += length;
}

}  // namespace

// -----------------------------------------------------------------------------

namespace WP2 {

class ImageReaderPNG : public ImageReader::Impl {
 public:
  ImageReaderPNG(const uint8_t* data, size_t data_size,
                 ArgbBuffer* const buffer, LogLevel log_level,
                 size_t max_num_pixels)
      : ImageReader::Impl(buffer, data, data_size, log_level, max_num_pixels) {}

  WP2Status ReadFrame(bool* const is_last,
                      uint32_t* const duration_ms) override {
    WP2_CHECK_STATUS(CheckData());

    // The following are volatile because they need to keep their value for
    // deletion in case of longjmp() error handling.
    png_infop volatile head_info = NULL;
    png_infop volatile end_info = NULL;
    tmp_rgb_ = nullptr;

    png_error_ptr error_function =
        (log_level_ >= LogLevel::DEFAULT) ? ErrorFunction : ErrorFunctionNoLog;
    png_structp const png = png_create_read_struct(PNG_LIBPNG_VER_STRING, NULL,
                                                   error_function, NULL);
    WP2_CHECK_ALLOC_OK(png != NULL);

    if (setjmp(png_jmpbuf(png))) {
      // On a libpng error, error_function() is called, longjmping here.
      png_destroy_read_struct((png_structpp)&png, (png_infopp)&head_info,
                              (png_infopp)&end_info);
      WP2Free(tmp_rgb_);
      tmp_rgb_ = nullptr;  // for good measure
      return WP2_STATUS_BITSTREAM_ERROR;
    }

    if ((head_info = png_create_info_struct(png)) == NULL ||
        (end_info = png_create_info_struct(png)) == NULL) {
      png_destroy_read_struct((png_structpp)&png, (png_infopp)&head_info,
                              (png_infopp)&end_info);
      return WP2_STATUS_OUT_OF_MEMORY;
    }
    PNGReadContext context = {data_, data_size_, 0};
    png_set_read_fn(png, &context, ReadFunc);
    png_read_info(png, head_info);

    WP2Status status = ReadCanvas(png, head_info);
    WP2Free(tmp_rgb_);
    tmp_rgb_ = nullptr;

    // Look for metadata at both the beginning and the end of the PNG file,
    // giving preference to the head.
    if (status == WP2_STATUS_OK) {
      status = ExtractMetadataFromPNG(png, head_info, &buffer_->metadata,
                                      log_level_);
    }
    if (status == WP2_STATUS_OK) {
      png_read_end(png, end_info);
      status =
          ExtractMetadataFromPNG(png, end_info, &buffer_->metadata, log_level_);
    }

    png_destroy_read_struct((png_structpp)&png, (png_infopp)&head_info,
                            (png_infopp)&end_info);
    if (status == WP2_STATUS_OK) {
      *is_last = true;
      *duration_ms = ImageReader::kInfiniteDuration;
    }
    return status;
  }

 protected:
  // Reads pixels into 'buffer_'. Allocates 'tmp_rgb_'.
  WP2Status ReadCanvas(png_structp const png, png_inforp const head_info) {
    int color_type, bit_depth, interlaced;
    png_uint_32 width, height;
    WP2_CHECK_OK(png_get_IHDR(png, head_info, &width, &height, &bit_depth,
                              &color_type, &interlaced, NULL, NULL),
                 WP2_STATUS_BITSTREAM_ERROR);
    WP2_CHECK_STATUS(CheckDimensions((uint32_t)width, (uint32_t)height));

    png_set_strip_16(png);
    png_set_packing(png);
    if (color_type == PNG_COLOR_TYPE_PALETTE) {
      png_set_palette_to_rgb(png);
    } else if (color_type == PNG_COLOR_TYPE_GRAY ||
               color_type == PNG_COLOR_TYPE_GRAY_ALPHA) {
      if (bit_depth < 8) png_set_expand_gray_1_2_4_to_8(png);
      png_set_gray_to_rgb(png);
    }

    bool has_alpha = !!(color_type & PNG_COLOR_MASK_ALPHA);
    if (png_get_valid(png, head_info, PNG_INFO_tRNS)) {
      // sometimes the tRNS chunk is misplaced and libpng behaves as
      // if there's no transparency. See pngrutil.c:1084
      png_bytep trans = nullptr;
      int num_trans = 0;
      if (!png_get_tRNS(png, head_info, &trans, &num_trans, nullptr)) {
        num_trans = 0;
      }
      has_alpha = (num_trans > 0);
      if (has_alpha) png_set_tRNS_to_alpha(png);
    }

    // Apply gamma correction if needed.
    double image_gamma = 1 / 2.2;
    if (png_get_gAMA(png, head_info, &image_gamma)) {
      const double screen_gamma = 2.2;
      png_set_gamma(png, screen_gamma, image_gamma);
    }

    const int num_passes = png_set_interlace_handling(png);
    png_read_update_info(png, head_info);

    const uint32_t depth = (has_alpha ? 4 : 3) * sizeof(*tmp_rgb_);
    size_t stride = 0;
    WP2_CHECK_STATUS(MultFitsIn(width, depth, &stride));
    // Make sure whole size fits in size_t.
    WP2_CHECK_STATUS(MultFitsIn<size_t>(stride, height));

    WP2_CHECK_STATUS(buffer_->Resize((uint32_t)width, (uint32_t)height));

    // For interlaced PNG, libpng needs to retain the full image in memory.
    tmp_rgb_ = (uint8_t*)WP2Malloc(num_passes > 1 ? height : 1, stride);
    WP2_CHECK_ALLOC_OK(tmp_rgb_ != nullptr);

    for (int p = 0; p < num_passes; ++p) {
      png_bytep row = tmp_rgb_;
      for (uint32_t y = 0; y < height; ++y) {
        png_read_row(png, row, NULL);
        if (p == num_passes - 1) {
          WP2_CHECK_STATUS(
              buffer_->ImportRow(has_alpha ? WP2_RGBA_32 : WP2_RGB_24, y, row));
        }
        if (num_passes > 1) row += stride;
      }
    }

    return WP2_STATUS_OK;
  }

 private:
  // Easier to store it here than as an argument of ReadCanvas().
  // Volatile because allocation happens between setjmp() and longjmp() and
  // needs to be freed afterwards.
  uint8_t* volatile tmp_rgb_;
};

void ImageReader::SetImplPNG(ArgbBuffer* const buffer, LogLevel log_level,
                             size_t max_num_pixels) {
  impl_.reset(new (WP2Allocable::nothrow) ImageReaderPNG(
      data_.bytes, data_.size, buffer, log_level, max_num_pixels));
  if (impl_ == nullptr) status_ = WP2_STATUS_OUT_OF_MEMORY;
}

}  // namespace WP2

#else  // !WP2_HAVE_PNG

namespace WP2 {

void ImageReader::SetImplPNG(ArgbBuffer* const buffer, LogLevel log_level,
                             size_t max_num_pixels) {
  (void)buffer;
  (void)max_num_pixels;
  if (log_level >= WP2::LogLevel::DEFAULT) {
    fprintf(stderr, "PNG support not compiled. Please install the libpng "
                    "development package before building.\n");
  }
  status_ = WP2_STATUS_VERSION_MISMATCH;
}

}  // namespace WP2

#endif  // WP2_HAVE_PNG

// -----------------------------------------------------------------------------
