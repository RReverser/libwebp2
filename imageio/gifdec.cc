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
// GIF decoding.
// Note: Legacy sequential functions are used, not DGifSlurp(). From giflib:
// "If you are handling large images on an extremely memory-limited machine,
//  you may need to use the functions for sequential read and write."
//
// Author: Yannis Guyon (yguyon@google.com)

#include <algorithm>
#include <cassert>
#include <cstdio>
#include <cstring>

#include "./anim_image_dec.h"

#ifdef HAVE_CONFIG_H
#include "wp2/config.h"
#endif

//------------------------------------------------------------------------------

#ifdef WP2_HAVE_GIF

#include <gif_lib.h>

#include "src/dsp/dsp.h"
#include "src/utils/utils.h"
#include "src/wp2/base.h"

// GIFLIB_MAJOR is only defined in libgif >= 4.2.0.
#if defined(GIFLIB_MAJOR) && defined(GIFLIB_MINOR)
#define LOCAL_GIF_VERSION ((GIFLIB_MAJOR << 8) | GIFLIB_MINOR)
#define LOCAL_GIF_PREREQ(maj, min) (LOCAL_GIF_VERSION >= (((maj) << 8) | (min)))
#else
#define LOCAL_GIF_VERSION 0
#define LOCAL_GIF_PREREQ(maj, min) 0
#endif

//------------------------------------------------------------------------------

namespace WP2 {
namespace {

// Returns the WP2Status corresponding to the D_GIF_ERR and logs it if not OK.
WP2Status ErrorToStatus(const GifFileType* const gif, int gif_error,
                        LogLevel log_level) {
#if LOCAL_GIF_PREREQ(5, 0)
  if (gif != NULL) gif_error = gif->Error;
#else
  (void)gif;
#endif
  if (gif_error == D_GIF_SUCCEEDED) return WP2_STATUS_OK;

  if (log_level >= LogLevel::DEFAULT) {
    // libgif 4.2.0 has retired PrintGifError() and added GifErrorString().
#if LOCAL_GIF_PREREQ(4, 2)
#if LOCAL_GIF_PREREQ(5, 0)
    // Static string actually, hence the const char* cast.
    const char* error_str = (const char*)GifErrorString(gif_error);
#else
    const char* error_str = (const char*)GifErrorString();
#endif
    if (error_str == NULL) error_str = "Unknown error";
    fprintf(stderr, "GIFLib Error %d: %s\n", gif_error, error_str);
#else
    fprintf(stderr, "GIFLib Error %d: ", gif_error);
    PrintGifError();
    fprintf(stderr, "\n");
#endif
  }

  switch (gif_error) {
    case D_GIF_ERR_OPEN_FAILED:
    case D_GIF_ERR_READ_FAILED:
    case D_GIF_ERR_CLOSE_FAILED:
    case D_GIF_ERR_NOT_READABLE:
      return WP2_STATUS_BAD_READ;
    case D_GIF_ERR_EOF_TOO_SOON:
      return WP2_STATUS_NOT_ENOUGH_DATA;
    case D_GIF_ERR_DATA_TOO_BIG:
      return WP2_STATUS_FILE_TOO_BIG;
    case D_GIF_ERR_NOT_ENOUGH_MEM:
      return WP2_STATUS_OUT_OF_MEMORY;
    default:
      return WP2_STATUS_BITSTREAM_ERROR;
  }
}

constexpr int kInvalidIndex = -1;
constexpr int kDefaultDurationMs = 100;
constexpr Argb32b kTransparentColor{0x00u, 0x00u, 0x00u, 0x00u};
constexpr Argb32b kWhiteColor{0xFFu, 0xFFu, 0xFFu, 0xFFu};

enum class DisposeMethod { NONE, BACKGROUND, RESTORE_PREVIOUS };

bool IsConversionNeeded(WP2SampleFormat format) {
  return (format != WP2_Argb_32 && format != WP2_ARGB_32);
}

//------------------------------------------------------------------------------

// Extracts background color (to dispose frames with).
Argb32b GetBackgroundColor(const ColorMapObject* const color_map,
                           int bgcolor_index, int transparent_index) {
  if (transparent_index != kInvalidIndex &&
      bgcolor_index == transparent_index) {
    return kTransparentColor;  // Special case.
  } else if (color_map == NULL || color_map->Colors == NULL ||
             bgcolor_index >= color_map->ColorCount) {
    // Invalid background color index. Assuming white background.
    return kWhiteColor;
  } else {
    const GifColorType c = color_map->Colors[bgcolor_index];
    return {0xFFu, c.Red, c.Green, c.Blue};
  }
}

// Converts 'indexed_color_pixels' and copies to 'pixels' in 'pixel_format'.
WP2Status ConvertIndexedColors(const GifFileType* const decoder,
                               const uint8_t* const indexed_color_pixels,
                               uint32_t* const temp_argb_pixels,
                               uint32_t num_pixels, int transparent_index,
                               uint8_t* const pixels,
                               WP2SampleFormat pixel_format) {
  const ColorMapObject* const color_map =
      (decoder->Image.ColorMap != NULL) ? decoder->Image.ColorMap
                                        : decoder->SColorMap;
  if (color_map == NULL) return WP2_STATUS_OK;  // Nothing to do.
  WP2_CHECK_OK(color_map->Colors != NULL && color_map->ColorCount > 0,
               WP2_STATUS_BITSTREAM_ERROR);

  assert(IsConversionNeeded(pixel_format) == (temp_argb_pixels != nullptr));
  if (IsConversionNeeded(pixel_format)) {
    WP2ArgbConverterInit();
    // For transparent pixels that reuse the previous value.
    WP2ArgbConvertFrom[pixel_format]((const uint8_t*)pixels, num_pixels,
                                     (uint8_t*)temp_argb_pixels);
  }
  uint32_t* const dst_pixels =
      (IsConversionNeeded(pixel_format) ? temp_argb_pixels : (uint32_t*)pixels);

  for (uint32_t i = 0; i < num_pixels; ++i) {
    if (indexed_color_pixels[i] == transparent_index) {
      // Do nothing. A tranparent pixel implies keeping previous color.
    } else if (indexed_color_pixels[i] < color_map->ColorCount) {
      const GifColorType color = color_map->Colors[indexed_color_pixels[i]];
      dst_pixels[i] = (uint32_t)0xffu | (color.Red << 8) |
                      (color.Green << 16) | (color.Blue << 24);
    } else {
      return WP2_STATUS_BITSTREAM_ERROR;
    }
  }

  if (IsConversionNeeded(pixel_format)) {
    WP2ArgbConvertTo[pixel_format]((const uint8_t*)temp_argb_pixels, num_pixels,
                                   pixels);
  }
  return WP2_STATUS_OK;
}

// Decodes pixels into 'frame'.
WP2Status DecodeFrame(GifFileType* const decoder, int transparent_index,
                      ArgbBuffer* const frame, uint8_t* const tmp_indexed_row,
                      uint32_t* const tmp_argb_row, LogLevel log_level) {
  if (decoder->Image.Interlace) {
    static constexpr uint32_t kNumPasses = 4;
    static constexpr uint32_t kOffsets[kNumPasses] = {0, 4, 2, 1};
    static constexpr uint32_t kJumps[kNumPasses] = {8, 8, 4, 2};
    for (uint32_t pass = 0; pass < kNumPasses; ++pass) {
      uint32_t y = kOffsets[pass];
      // It is not written in GIF specs if an interlaced image can be <5px tall.
      if (y < frame->height) {
        uint8_t* row = (uint8_t*)frame->GetRow(y);
        const uint32_t jump = kJumps[pass] * frame->stride;
        for (; y < frame->height; y += kJumps[pass], row += jump) {
          WP2_CHECK_STATUS(ErrorToStatus(
              decoder, DGifGetLine(decoder, tmp_indexed_row, frame->width),
              log_level));
          WP2_CHECK_STATUS(ConvertIndexedColors(
              decoder, tmp_indexed_row, tmp_argb_row, frame->width,
              transparent_index, row, frame->format));
        }
      }
    }
  } else {
    uint8_t* row = (uint8_t*)frame->GetRow(0);
    for (uint32_t y = 0; y < frame->height; ++y, row += frame->stride) {
      WP2_CHECK_STATUS(ErrorToStatus(
          decoder, DGifGetLine(decoder, tmp_indexed_row, frame->width),
          log_level));
      WP2_CHECK_STATUS(ConvertIndexedColors(
          decoder, tmp_indexed_row, tmp_argb_row, frame->width,
          transparent_index, row, frame->format));
    }
  }
  return WP2_STATUS_OK;
}

//------------------------------------------------------------------------------

// Decodes frame duration, dispose method, transparent color index.
WP2Status ReadGraphicsExtension(const GifByteType* const data,
                                uint32_t* const frame_duration_ms,
                                DisposeMethod* const frame_dispose_method_,
                                int* const transparent_index) {
  static constexpr int kDisposeShift = 2;
  static constexpr int kDisposeMask = 0x07;
  static constexpr int kTransparentMask = 0x01;

  const int flags = data[1];
  const int dispose_raw = (flags >> kDisposeShift) & kDisposeMask;
  const int duration_raw = data[2] | (data[3] << 8);  // In 10 ms units.
  WP2_CHECK_OK(data[0] == 4, WP2_STATUS_BITSTREAM_ERROR);
  *frame_duration_ms = (uint32_t)duration_raw * 10;
  if (dispose_raw == 2) {
    *frame_dispose_method_ = DisposeMethod::BACKGROUND;
  } else if (dispose_raw == 3) {
    *frame_dispose_method_ = DisposeMethod::RESTORE_PREVIOUS;
  } else {
    *frame_dispose_method_ = DisposeMethod::NONE;
  }
  *transparent_index = (flags & kTransparentMask) ? data[4] : kInvalidIndex;
  return WP2_STATUS_OK;
}

// Decodes loop count, metadata.
WP2Status ReadApplicationExtension(GifFileType* const decoder,
                                   const GifByteType* const data,
                                   uint32_t* const loop_count,
                                   LogLevel log_level) {
  if (data[0] != 11) return WP2_STATUS_OK;  // Chunk is too short.
  if (!memcmp(data + 1, "NETSCAPE2.0", 11) ||
      !memcmp(data + 1, "ANIMEXTS1.0", 11)) {
    GifByteType* extension_data = NULL;
    WP2_CHECK_STATUS(ErrorToStatus(
        decoder, DGifGetExtensionNext(decoder, &extension_data), log_level));
    WP2_CHECK_OK(extension_data != NULL, WP2_STATUS_BITSTREAM_ERROR);
    WP2_CHECK_OK(extension_data[0] >= 3 && extension_data[1] == 1,
                 WP2_STATUS_BITSTREAM_ERROR);
    assert(loop_count != nullptr);
    *loop_count = (uint32_t)extension_data[2] | (extension_data[3] << 8);
    if (*loop_count == 0) *loop_count = ImageReader::kInfiniteLoopCount;
  }
  return WP2_STATUS_OK;
}

}  // namespace

//------------------------------------------------------------------------------

// Will be used by giflib.
static int InputFunction(GifFileType* decoder, GifByteType* data_source,
                         int num_bytes_to_read);

class ImageReaderGIF : public ImageReader::Impl {
 public:
  ImageReaderGIF(const uint8_t* data, size_t data_size,
                 ArgbBuffer* const buffer, LogLevel log_level,
                 size_t max_num_pixels)
      : ImageReader::Impl(buffer, data, data_size, log_level, max_num_pixels),
        previous_canvas_((buffer_ != nullptr) ? buffer_->format : WP2_Argb_32),
        frame_(previous_canvas_.format) {}
  ~ImageReaderGIF() override {
    if (decoder_ != NULL) {
#if LOCAL_GIF_PREREQ(5, 1)
      int error = D_GIF_SUCCEEDED;
      DGifCloseFile(decoder_, &error);
      (void)ErrorToStatus(nullptr, error, log_level_);  // For error logging.
#else
      DGifCloseFile(decoder_);
#endif
    }
    WP2Free(tmp_indexed_row_);
    WP2Free(tmp_argb_row_);
  }

  // Returns the number of bytes successfully copied (at most 'num_bytes').
  size_t TryRead(size_t num_bytes, GifByteType* const into) {
    assert(num_read_bytes_ <= data_size_);
    num_bytes = std::min(num_bytes, data_size_ - num_read_bytes_);
    if (num_bytes > 0) {
      memcpy((void*)into, (const void*)(data_ + num_read_bytes_), num_bytes);
      num_read_bytes_ += num_bytes;
    }
    return num_bytes;
  }

 protected:
  // Allocates 'buffer_', 'previous_canvas_' and temporary arrays.
  WP2Status InitCanvas() {
    // Fix some broken GIF global headers that report 0 x 0 dimension.
    if (decoder_->SWidth == 0 || decoder_->SHeight == 0) {
      decoder_->Image.Left = 0;
      decoder_->Image.Top = 0;
      decoder_->SWidth = decoder_->Image.Width;
      decoder_->SHeight = decoder_->Image.Height;
      WP2_CHECK_OK(decoder_->SWidth > 0 && decoder_->SHeight > 0,
                   WP2_STATUS_BITSTREAM_ERROR);
    }
    WP2_CHECK_STATUS(CheckDimensions((uint32_t)decoder_->SWidth,
                                     (uint32_t)decoder_->SHeight));

    // Background color. Some decoders might always use transparent instead,
    // or might update it if the palette changes.
    background_color_ = GetBackgroundColor(
        decoder_->SColorMap, decoder_->SBackGroundColor, transparent_index_);

    // Allocate canvas and fill it with background color.
    WP2_CHECK_STATUS(buffer_->Resize((uint32_t)decoder_->SWidth,
                                     (uint32_t)decoder_->SHeight));
    buffer_->Fill(background_color_);
    WP2_CHECK_STATUS(previous_canvas_.CopyFrom(*buffer_));

    // Allocate temp rows (used for DGifGetLine()).
    tmp_indexed_row_ =
        (uint8_t*)WP2Malloc(buffer_->width, sizeof(*tmp_indexed_row_));
    WP2_CHECK_ALLOC_OK(tmp_indexed_row_ != nullptr);
    if (IsConversionNeeded(buffer_->format)) {
      tmp_argb_row_ =
          (uint32_t*)WP2Malloc(buffer_->width, sizeof(*tmp_argb_row_));
      WP2_CHECK_ALLOC_OK(tmp_argb_row_ != nullptr);
    }
    return WP2_STATUS_OK;
  }

  // Clears the canvas depending on 'previous_frame_dispose_method_'.
  WP2Status DisposeCanvas() {
    if (previous_frame_dispose_method_ == DisposeMethod::BACKGROUND) {
      frame_.Fill(background_color_);
    } else if (previous_frame_dispose_method_ ==
               DisposeMethod::RESTORE_PREVIOUS) {
      const uint8_t* const previous_pixels_inside_frame_rect_ =
          (const uint8_t*)previous_canvas_.GetRow(frame_rect_.y) +
          frame_rect_.x * WP2FormatBpp(previous_canvas_.format);
      WP2_CHECK_STATUS(frame_.Import(
          previous_canvas_.format, frame_rect_.width, frame_rect_.height,
          previous_pixels_inside_frame_rect_, previous_canvas_.stride));
    }
    WP2_CHECK_STATUS(previous_canvas_.CopyFrom(*buffer_));
    return WP2_STATUS_OK;
  }

 public:
  WP2Status ReadFrame(bool* const is_last,
                      uint32_t* const duration_ms) override {
    WP2_CHECK_STATUS(CheckData());

    GifRecordType record_type = TERMINATE_RECORD_TYPE;
    if (decoder_ == NULL) {
      // First frame.
      // Metadata must be decoded now but Giflib DGifGetRecordType() doesn't
      // work without DGifGetLine(), hence the use of ReadGIFMetadata().
      WP2_CHECK_STATUS(ReadGIFMetadata(data_, data_size_, &buffer_->metadata));

      int error = D_GIF_SUCCEEDED;
      decoder_ = DGifOpen((void*)this, InputFunction, &error);
      WP2_CHECK_STATUS(ErrorToStatus(decoder_, error, log_level_));
      WP2_CHECK_ALLOC_OK(decoder_ != NULL);
      WP2_CHECK_STATUS(ErrorToStatus(
          decoder_, DGifGetRecordType(decoder_, &record_type), log_level_));
    } else {
      // Not the first frame, we already know the next 'record_type'.
      record_type = IMAGE_DESC_RECORD_TYPE;
      ++frame_index_;
    }

    // In GIF format the duration is stored after the frame. The following
    // records are parsed till the next frame to know if it was the last one.
    bool just_decoded_frame = false;
    while (record_type != TERMINATE_RECORD_TYPE) {
      if (record_type == IMAGE_DESC_RECORD_TYPE) {
        if (just_decoded_frame) {
          *is_last = false;
          return WP2_STATUS_OK;
        }

        WP2_CHECK_STATUS(
            ErrorToStatus(decoder_, DGifGetImageDesc(decoder_), log_level_));

        if (frame_index_ == 0) {
          WP2_CHECK_STATUS(InitCanvas());
        } else {
          WP2_CHECK_STATUS(DisposeCanvas());
        }

        // Some broken GIF can have sub-rect with zero width/height.
        if (decoder_->Image.Width == 0 || decoder_->Image.Height == 0) {
          decoder_->Image.Width = decoder_->SWidth;
          decoder_->Image.Height = decoder_->SHeight;
        }

        frame_rect_ = {
            (uint32_t)decoder_->Image.Left,
            (uint32_t)decoder_->Image.Top,
            (uint32_t)decoder_->Image.Width,
            (uint32_t)decoder_->Image.Height
        };
        WP2_CHECK_STATUS(frame_.SetView(*buffer_, frame_rect_));
        WP2_CHECK_STATUS(DecodeFrame(decoder_, transparent_index_, &frame_,
                                     tmp_indexed_row_, tmp_argb_row_,
                                     log_level_));
        just_decoded_frame = true;

        // In GIF, graphic control extensions are optional, so we may not get
        // one before reading the next frame. To handle this case, we reset
        // frame properties to reasonable defaults.
        previous_frame_dispose_method_ = frame_dispose_method_;
        frame_dispose_method_ = DisposeMethod::NONE;
        *duration_ms = kDefaultDurationMs;
        transparent_index_ = kInvalidIndex;
      } else if (record_type == EXTENSION_RECORD_TYPE) {
        int extension;
        GifByteType* data = NULL;
        WP2_CHECK_STATUS(ErrorToStatus(
            decoder_, DGifGetExtension(decoder_, &extension, &data),
            log_level_));
        if (data != NULL) {
          if (extension == GRAPHICS_EXT_FUNC_CODE) {
            uint32_t frame_duration_ms = 0;
            WP2_CHECK_STATUS(ReadGraphicsExtension(data, &frame_duration_ms,
                                                   &frame_dispose_method_,
                                                   &transparent_index_));
            // Avoid instant frames from unconsistent GIF images.
            *duration_ms = (frame_duration_ms > 0) ? frame_duration_ms
                                                   : kDefaultDurationMs;
          } else if (extension == APPLICATION_EXT_FUNC_CODE) {
            WP2_CHECK_STATUS(ReadApplicationExtension(
                decoder_, data, &loop_count_, log_level_));
          }

          do {
            WP2_CHECK_STATUS(ErrorToStatus(
                decoder_, DGifGetExtensionNext(decoder_, &data), log_level_));
          } while (data != NULL);
        }
      } else {
        // Skipping over unknown record type.
      }

      WP2_CHECK_STATUS(ErrorToStatus(
          decoder_, DGifGetRecordType(decoder_, &record_type), log_level_));
    }

    WP2_CHECK_OK(just_decoded_frame, WP2_STATUS_BITSTREAM_ERROR);
    *is_last = true;
    return WP2_STATUS_OK;
  }

 protected:
  GifFileType* decoder_ = NULL;
  size_t num_read_bytes_ = 0;

  ArgbBuffer previous_canvas_;
  ArgbBuffer frame_;  // Is a subview of buffer_, which represents the canvas.
  Rectangle frame_rect_;
  uint32_t frame_index_ = 0;
  DisposeMethod frame_dispose_method_ = DisposeMethod::NONE;
  DisposeMethod previous_frame_dispose_method_ = DisposeMethod::NONE;
  int transparent_index_ = kInvalidIndex;
  Argb32b background_color_ = kTransparentColor;
  uint8_t* tmp_indexed_row_ = nullptr;
  uint32_t* tmp_argb_row_ = nullptr;
};

// -----------------------------------------------------------------------------

static int InputFunction(GifFileType* decoder, GifByteType* data_source,
                         int num_bytes_to_read) {
  ImageReaderGIF* const user_data = (ImageReaderGIF*)decoder->UserData;
  assert(user_data != nullptr);
  return (int)user_data->TryRead((size_t)num_bytes_to_read, data_source);
}

// -----------------------------------------------------------------------------

void ImageReader::SetImplGIF(ArgbBuffer* const buffer,
                             LogLevel log_level, size_t max_num_pixels) {
  impl_.reset(new (WP2Allocable::nothrow) ImageReaderGIF(
      data_.bytes, data_.size, buffer, log_level, max_num_pixels));
  if (impl_ == nullptr) status_ = WP2_STATUS_OUT_OF_MEMORY;
}

}  // namespace WP2

#else  // !WP2_HAVE_GIF

void WP2::ImageReader::SetImplGIF(WP2::ArgbBuffer* const buffer,
                                  WP2::LogLevel log_level,
                                  size_t max_num_pixels) {
  (void)buffer;
  (void)max_num_pixels;
  if (log_level >= WP2::LogLevel::DEFAULT) {
    fprintf(stderr, "GIF support not compiled. Please install the libgif "
                    "development package before building.\n");
  }
  status_ = WP2_STATUS_VERSION_MISMATCH;
}

#endif  // WP2_HAVE_GIF

// -----------------------------------------------------------------------------
