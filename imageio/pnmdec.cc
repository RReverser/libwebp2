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
// (limited) PNM decoder

#include <cassert>
#include <cstdio>
#include <cstring>

#include "./anim_image_dec.h"
#include "./imageio_util.h"
#include "src/utils/utils.h"

namespace {

enum PNMFlags {
  WIDTH_FLAG  = 1 << 0,
  HEIGHT_FLAG = 1 << 1,
  DEPTH_FLAG  = 1 << 2,
  MAXVAL_FLAG = 1 << 3,
  TUPLE_FLAG  = 1 << 4,
  ALL_NEEDED_FLAGS = WIDTH_FLAG | HEIGHT_FLAG | DEPTH_FLAG | MAXVAL_FLAG
};

// See https://en.wikipedia.org/wiki/Netpbm_format
enum class PNMType {
  UNHANDLED     = 0,
  GRAY_MAP      = 5,  // .pgm
  PIX_MAP       = 6,  // .ppm
  ARBITRARY_MAP = 7,  // .pam
};

struct PNMInfo {
  const uint8_t* data = nullptr;
  size_t data_size = 0;
  uint32_t width = 0;
  uint32_t height = 0;
  int bytes_per_px = 0;
  int depth = 0;         // 1, 2, 3, or 4
  int max_value = 0;
  PNMType type = PNMType::UNHANDLED;
  int seen_flags = 0;
};

// -----------------------------------------------------------------------------
// PNM decoding

constexpr size_t kMaxLineSize = 1024;
constexpr size_t kMinPNMHeaderSize = 3;

WP2Status ReadLine(const uint8_t* const data, size_t* const data_pos,
                   size_t data_size, char out[kMaxLineSize + 1],
                   size_t* const out_size) {
  size_t i = 0;
  do {
    for (i = 0; i < kMaxLineSize && *data_pos < data_size; ++i) {
      out[i] = data[(*data_pos)++];
      if (out[i] == '\n') break;
    }
    // Skip empty lines and comments.
  } while (*data_pos < data_size && (i == 0 || out[0] == '#'));
  out[i] = 0;  // safety sentinel
  *out_size = i;
  return (*out_size > 0 ? WP2_STATUS_OK : WP2_STATUS_NOT_ENOUGH_DATA);
}

WP2Status FlagError(const char flag[], WP2::LogLevel log_level) {
  if (log_level >= WP2::LogLevel::DEFAULT) {
    fprintf(stderr, "PAM header error: flags '%s' already seen.\n", flag);
  }
  return WP2_STATUS_BITSTREAM_ERROR;
}

// inspired from http://netpbm.sourceforge.net/doc/pam.html
WP2Status ReadPAMFields(PNMInfo* const info, size_t* const data_pos,
                        WP2::LogLevel log_level) {
  assert(info != nullptr);
  int expected_depth = -1;
  while (true) {
    char out[kMaxLineSize + 1];
    size_t out_size;
    WP2_CHECK_STATUS(
        ReadLine(info->data, data_pos, info->data_size, out, &out_size));
    uint32_t tmp;
    if (sscanf(out, "WIDTH %u", &tmp) == 1) {
      if (info->seen_flags & WIDTH_FLAG) return FlagError("WIDTH", log_level);
      info->seen_flags |= WIDTH_FLAG;
      info->width = tmp;
    } else if (sscanf(out, "HEIGHT %u", &tmp) == 1) {
      if (info->seen_flags & HEIGHT_FLAG) return FlagError("HEIGHT", log_level);
      info->seen_flags |= HEIGHT_FLAG;
      info->height = tmp;
    } else if (sscanf(out, "DEPTH %u", &tmp) == 1) {
      if (info->seen_flags & DEPTH_FLAG) return FlagError("DEPTH", log_level);
      info->seen_flags |= DEPTH_FLAG;
      info->depth = tmp;
    } else if (sscanf(out, "MAXVAL %u", &tmp) == 1) {
      if (info->seen_flags & MAXVAL_FLAG) return FlagError("MAXVAL", log_level);
      info->seen_flags |= MAXVAL_FLAG;
      info->max_value = tmp;
    } else if (!strcmp(out, "TUPLTYPE RGB_ALPHA")) {
      expected_depth = 4;
      info->seen_flags |= TUPLE_FLAG;
    } else if (!strcmp(out, "TUPLTYPE RGB")) {
      expected_depth = 3;
      info->seen_flags |= TUPLE_FLAG;
    } else if (!strcmp(out, "TUPLTYPE GRAYSCALE_ALPHA")) {
      expected_depth = 2;
      info->seen_flags |= TUPLE_FLAG;
    } else if (!strcmp(out, "TUPLTYPE GRAYSCALE")) {
      expected_depth = 1;
      info->seen_flags |= TUPLE_FLAG;
    } else if (!strcmp(out, "ENDHDR")) {
      break;
    } else {
      if (log_level >= WP2::LogLevel::DEFAULT) {
        static const char kEllipsis[] = " ...";
        if (out_size > 20) {
          sprintf(out + 20 - strlen(kEllipsis), kEllipsis);  // NOLINT
        }
        fprintf(stderr, "PAM header error: unrecognized entry [%s]\n", out);
      }
      return WP2_STATUS_BITSTREAM_ERROR;
    }
  }
  if (!(info->seen_flags & ALL_NEEDED_FLAGS)) {
    if (log_level >= WP2::LogLevel::DEFAULT) {
      fprintf(stderr, "PAM header error: missing tags%s%s%s%s\n",
              (info->seen_flags & WIDTH_FLAG) ? "" : " WIDTH",
              (info->seen_flags & HEIGHT_FLAG) ? "" : " HEIGHT",
              (info->seen_flags & DEPTH_FLAG) ? "" : " DEPTH",
              (info->seen_flags & MAXVAL_FLAG) ? "" : " MAXVAL");
    }
    return WP2_STATUS_BITSTREAM_ERROR;
  }
  if (expected_depth != -1 && info->depth != expected_depth) {
    if (log_level >= WP2::LogLevel::DEFAULT) {
      fprintf(stderr, "PAM header error: expected DEPTH %d but got DEPTH %d\n",
              expected_depth, info->depth);
    }
    return WP2_STATUS_BITSTREAM_ERROR;
  }
  return WP2_STATUS_OK;
}

WP2Status ReadHeader(PNMInfo* const info, size_t* const data_pos,
                     WP2::LogLevel log_level) {
  WP2_CHECK_OK(info != nullptr && info->data != nullptr && data_pos != nullptr,
               WP2_STATUS_NULL_PARAMETER);
  WP2_CHECK_OK(info->data_size >= kMinPNMHeaderSize,
               WP2_STATUS_NOT_ENOUGH_DATA);

  char out[kMaxLineSize + 1];
  size_t out_size;

  WP2_CHECK_STATUS(
      ReadLine(info->data, data_pos, info->data_size, out, &out_size));
  int type;
  WP2_CHECK_OK(sscanf(out, "P%d", &type) == 1, WP2_STATUS_BITSTREAM_ERROR);
  WP2_CHECK_OK(type >= 1 && type <= 7, WP2_STATUS_BITSTREAM_ERROR);
  if (type < 5) {
    if (log_level >= WP2::LogLevel::DEFAULT) {
      fprintf(stderr, "Unsupported P%d PNM format.\n", type);
    }
    return WP2_STATUS_UNSUPPORTED_FEATURE;
  }
  info->type = (PNMType)type;

  if (info->type == PNMType::GRAY_MAP || info->type == PNMType::PIX_MAP) {
    WP2_CHECK_STATUS(
        ReadLine(info->data, data_pos, info->data_size, out, &out_size));
    WP2_CHECK_OK(sscanf(out, "%d %d", &info->width, &info->height) == 2,
                 WP2_STATUS_BITSTREAM_ERROR);
    WP2_CHECK_STATUS(
        ReadLine(info->data, data_pos, info->data_size, out, &out_size));
    WP2_CHECK_OK(sscanf(out, "%d", &info->max_value) == 1,
                 WP2_STATUS_BITSTREAM_ERROR);

    // finish initializing missing fields
    info->depth = (info->type == PNMType::GRAY_MAP) ? 1 : 3;
  } else if (info->type == PNMType::ARBITRARY_MAP) {
    WP2_CHECK_STATUS(ReadPAMFields(info, data_pos, log_level));
  } else {
    assert(false);  // Already checked above.
  }

  // perform some basic numerical validation
  WP2_CHECK_OK(info->width > 0 && info->height > 0 &&
                   (info->depth >= 1 && info->depth <= 4) &&
                   info->max_value > 0 && info->max_value < 65536,
               WP2_STATUS_BITSTREAM_ERROR);
  info->bytes_per_px = info->depth * (info->max_value > 255 ? 2 : 1);
  return WP2_STATUS_OK;
}

}  // namespace

// -----------------------------------------------------------------------------

namespace WP2 {

class ImageReaderPNM : public ImageReader::Impl {
 public:
  ImageReaderPNM(const uint8_t* data, size_t data_size,
                 ArgbBuffer* const buffer, LogLevel log_level,
                 size_t max_num_pixels)
      : ImageReader::Impl(buffer, data, data_size, log_level, max_num_pixels) {}

  WP2Status ReadFrame(bool* const is_last,
                      uint32_t* const duration_ms) override {
    WP2_CHECK_STATUS(CheckData());  // Some basic validations.

    PNMInfo info;
    info.data = data_;
    info.data_size = data_size_;
    size_t data_pos = 0;
    WP2_CHECK_STATUS(ReadHeader(&info, &data_pos, log_level_));
    WP2_CHECK_STATUS(CheckDimensions(info.width, info.height));

    const uint64_t sample_size = (info.max_value > 255 ? 2 : 1);
    const uint64_t pixel_bytes =
        (uint64_t)info.width * info.height * info.bytes_per_px;
    if (data_size_ < data_pos + pixel_bytes) {
      if (log_level_ >= LogLevel::DEFAULT) {
        fprintf(stderr, "Truncated PNM file (P%d).\n", (int)info.type);
      }
      return WP2_STATUS_BITSTREAM_ERROR;
    }

    // ArgbBuffer::Resize() will also check that dimensions are acceptable.
    WP2_CHECK_STATUS(buffer_->Resize(info.width, info.height));

    Data rgb;  // Temp row.
    const uint32_t round = info.max_value >> 1;
    const WP2SampleFormat format =
         (info.depth == 1 || info.depth == 3) ? WP2_RGB_24 : WP2_RGBA_32;
    if (info.depth <= 2) {   // Convert grayscale(/A) -> RGB(A)
      const uint32_t channels = (info.depth == 1) ? 1 : 2;
      WP2_CHECK_STATUS(
          rgb.Resize(info.width * (2 + info.depth), /*keep_bytes=*/false));
      for (uint32_t j = 0; j < info.height; ++j) {
        const uint8_t* const in = data_ + data_pos;
        for (uint32_t k = 0, i = 0; i < info.width * channels; ++i) {
          uint32_t v = (sample_size == 2) ? 256u * in[2 * i + 0] + in[2 * i + 1]
                     : in[i];
          if (info.max_value != 255) v = (v * 255u + round) / info.max_value;
          if ((i % channels) == 0) {
            rgb.bytes[k + 0] = rgb.bytes[k + 1] = rgb.bytes[k + 2] = v;
            k += 3;
          } else {
            rgb.bytes[k] = v;   // alpha
            k += 1;
          }
        }
        WP2_CHECK_STATUS(buffer_->ImportRow(format, j, rgb.bytes));
        data_pos += info.width * sample_size * channels;
      }
    } else if (info.depth == 3 || info.depth == 4) {   // RGB or RGBA
      if (sample_size == 2 || info.max_value != 255) {
        WP2_CHECK_STATUS(
            rgb.Resize(info.width * info.depth, /*keep_bytes=*/false));
      }
      for (uint32_t j = 0; j < info.height; ++j) {
        const uint8_t* in = data_ + data_pos;
        data_pos += info.bytes_per_px * info.width;
        if (sample_size == 2 || info.max_value != 255) {
          for (uint32_t i = 0; i < info.width * info.depth; ++i) {
            const uint32_t v =
                 (sample_size == 2) ? 256u * in[2 * i + 0] + in[2 * i + 1] :
                 in[i];
            rgb.bytes[i] = (v * 255u + round) / info.max_value;
          }
          in = rgb.bytes;
        }
        WP2_CHECK_STATUS(buffer_->ImportRow(format, j, in));
      }
    } else {
      assert(false);  // Already checked in ReadHeader().
    }

    *is_last = true;
    *duration_ms = ImageReader::kInfiniteDuration;
    return WP2_STATUS_OK;
  }
};

void ImageReader::SetImplPNM(ArgbBuffer* const buffer, LogLevel log_level,
                             size_t max_num_pixels) {
  impl_.reset(new (WP2Allocable::nothrow) ImageReaderPNM(
      data_.bytes, data_.size, buffer, log_level, max_num_pixels));
  if (impl_ == nullptr) status_ = WP2_STATUS_OUT_OF_MEMORY;
}

}  // namespace WP2

// -----------------------------------------------------------------------------
