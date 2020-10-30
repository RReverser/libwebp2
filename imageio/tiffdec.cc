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
// TIFF decode.
//
// Author: Skal (pascal.massimino@gmail.com)

#include "./anim_image_dec.h"

#ifdef HAVE_CONFIG_H
#include "src/wp2/config.h"
#endif

#include <cstdio>
#include <cstring>

#ifdef WP2_HAVE_TIFF
#include <tiffio.h>

#include "./imageio_util.h"

namespace {

// -----------------------------------------------------------------------------

// Fills 'metadata' with extracted ICCP and/or XMP chunks.
WP2Status ExtractMetadataFromTIFF(TIFF* const tiff,
                                  WP2::Metadata* const metadata,
                                  WP2::LogLevel log_level) {
  void* tag_data = nullptr;
  uint32 tag_data_len = 0;
  if (TIFFGetField(tiff, TIFFTAG_ICCPROFILE, &tag_data_len, &tag_data)) {
    WP2_CHECK_STATUS(
        metadata->iccp.CopyFrom((const uint8_t*)tag_data, tag_data_len));
  }
  if (TIFFGetField(tiff, TIFFTAG_XMLPACKET, &tag_data_len, &tag_data)) {
    WP2_CHECK_STATUS(
        metadata->xmp.CopyFrom((const uint8_t*)tag_data, tag_data_len));
  }

  // TODO(jzern): To extract the raw EXIF directory some parsing of it would be
  // necessary to determine the overall size. In addition, value offsets in
  // individual directory entries may need to be updated as, depending on the
  // type, they are file based.
  // Exif 2.2 Section 4.6.2 Tag Structure
  // TIFF Revision 6.0 Part 1 Section 2 TIFF Structure #Image File Directory
  toff_t exif_ifd_offset;
  if (TIFFGetField(tiff, TIFFTAG_EXIFIFD, &exif_ifd_offset) &&
      log_level >= WP2::LogLevel::VERBOSE) {
    fprintf(stderr, "Warning: EXIF extraction from TIFF is unsupported.\n");
  }
  return WP2_STATUS_OK;
}

// -----------------------------------------------------------------------------

// Ad-hoc structure to supply read-from-memory functionalities.
struct MyData {
  const uint8_t* data;
  toff_t size;
  toff_t pos;
};

int MyClose(thandle_t opaque) {
  (void)opaque;
  return 0;
}

toff_t MySize(thandle_t opaque) {
  const MyData* const my_data = (MyData*)opaque;
  return my_data->size;
}

toff_t MySeek(thandle_t opaque, toff_t offset, int whence) {
  MyData* const my_data = (MyData*)opaque;
  offset += (whence == SEEK_CUR) ? my_data->pos
          : (whence == SEEK_SET) ? 0
          : my_data->size;
  if (offset > my_data->size) return (toff_t)-1;
  my_data->pos = offset;
  return offset;
}

int MyMapFile(thandle_t opaque, void** base, toff_t* size) {
  (void)opaque;
  (void)base;
  (void)size;
  return 0;
}
void MyUnmapFile(thandle_t opaque, void* base, toff_t size) {
  (void)opaque;
  (void)base;
  (void)size;
}

tsize_t MyRead(thandle_t opaque, void* dst, tsize_t size) {
  MyData* const my_data = (MyData*)opaque;
  if (my_data->pos + size > my_data->size) {
    size = my_data->size - my_data->pos;
  }
  if (size > 0) {
    memcpy(dst, my_data->data + my_data->pos, size);
    my_data->pos += size;
  }
  return size;
}

// -----------------------------------------------------------------------------

// Auto-freeing raw data struct.
struct TIFFData {
  ~TIFFData() { Free(); }
  WP2Status Alloc(tsize_t size, WP2::LogLevel log_level) {
    Free();
    bytes = (uint32*)_TIFFmalloc((tsize_t)size);
    if (bytes == NULL) {
      if (log_level >= WP2::LogLevel::DEFAULT) {
        fprintf(stderr, "Error allocating TIFF memory!\n");
      }
      return WP2_STATUS_OUT_OF_MEMORY;
    }
    return WP2_STATUS_OK;
  }
  void Free() {
    if (bytes != NULL) {
      _TIFFfree(bytes);
      bytes = NULL;
    }
  }
  uint32_t* bytes = NULL;
};

// -----------------------------------------------------------------------------

}  // namespace

namespace WP2 {

class ImageReaderTIFF : public ImageReader::Impl {
 public:
  ImageReaderTIFF(const uint8_t* data, size_t data_size,
                  ArgbBuffer* const buffer, LogLevel log_level,
                  size_t max_num_pixels)
      : ImageReader::Impl(buffer, data, data_size, log_level, max_num_pixels) {}

  WP2Status ReadFrame(bool* const is_last,
                      uint32_t* const duration_ms) override {
    WP2_CHECK_STATUS(CheckData());

    MyData my_data = {data_, (toff_t)data_size_, 0};
    TIFF* const tiff =
        TIFFClientOpen("Memory", "r", &my_data, MyRead, MyRead, MySeek, MyClose,
                       MySize, MyMapFile, MyUnmapFile);
    if (tiff == NULL) {
      if (log_level_ >= LogLevel::DEFAULT) {
        fprintf(stderr, "Error! Cannot parse TIFF file\n");
      }
      return WP2_STATUS_BITSTREAM_ERROR;
    }

    WP2Status status = DecodeTIFF(tiff, duration_ms);

    if (status == WP2_STATUS_OK) {
      status = ExtractMetadataFromTIFF(tiff, &buffer_->metadata, log_level_);
      if (status != WP2_STATUS_OK && log_level_ >= LogLevel::DEFAULT) {
        fprintf(stderr, "Error! Cannot extract TIFF metadata\n");
      }
    }

    if (status != WP2_STATUS_OK && log_level_ >= LogLevel::DEFAULT) {
      fprintf(stderr, "Status: %d (%s)\n", status, WP2GetStatusMessage(status));
    }

    TIFFClose(tiff);
    if (status == WP2_STATUS_OK) {
      *is_last = true;
      *duration_ms = ImageReader::kInfiniteDuration;
    }
    return status;
  }

 protected:
  // Reads pixels into 'buffer_'.
  WP2Status DecodeTIFF(TIFF* const tiff, uint32_t* const duration_ms) {
    const tdir_t dircount = TIFFNumberOfDirectories(tiff);
    if (dircount > 1 && log_level_ >= LogLevel::VERBOSE) {
      fprintf(stderr,
              "Warning: multi-directory TIFF files are not supported.\n"
              "Only the first will be used, %d will be ignored.\n",
              dircount - 1);
    }
    uint16_t samples_per_px = 0;
    if (!TIFFGetFieldDefaulted(tiff, TIFFTAG_SAMPLESPERPIXEL,
                               &samples_per_px)) {
      if (log_level_ >= LogLevel::DEFAULT) {
        fprintf(stderr,
                "Error! Cannot retrieve TIFF samples-per-pixel info.\n");
      }
      return WP2_STATUS_BITSTREAM_ERROR;
    }
    if (samples_per_px < 3 || samples_per_px > 4) {
      if (log_level_ >= LogLevel::DEFAULT) {
        fprintf(stderr, "Error! Unsupported TIFF samples-per-pixel %u\n",
                samples_per_px);
      }
      return WP2_STATUS_UNSUPPORTED_FEATURE;
    }

    uint32_t image_width, image_height;
    if (!(TIFFGetField(tiff, TIFFTAG_IMAGEWIDTH, &image_width) &&
          TIFFGetField(tiff, TIFFTAG_IMAGELENGTH, &image_height))) {
      if (log_level_ >= LogLevel::DEFAULT) {
        fprintf(stderr, "Error! Cannot retrieve TIFF image dimensions.\n");
      }
      return WP2_STATUS_BITSTREAM_ERROR;
    }
    WP2_CHECK_STATUS(CheckDimensions(image_width, image_height));

    uint32_t tile_width, tile_height;
    if (TIFFGetField(tiff, TIFFTAG_TILEWIDTH, &tile_width) &&
        TIFFGetField(tiff, TIFFTAG_TILELENGTH, &tile_height)) {
      WP2_CHECK_STATUS(CheckDimensions(tile_width, tile_height));
    }

    TIFFData raster;
    uint32_t stride = 0;
    WP2_CHECK_STATUS(MultFitsIn(image_width, sizeof(*raster.bytes), &stride));

    uint16_t extra_samples = 0;
    uint16_t* extra_samples_ptr = nullptr;
    if (samples_per_px > 3 &&
        !TIFFGetField(tiff, TIFFTAG_EXTRASAMPLES, &extra_samples,
                      &extra_samples_ptr)) {
      if (log_level_ >= LogLevel::DEFAULT) {
        fprintf(stderr, "Error! Cannot retrieve TIFF ExtraSamples info.\n");
      }
      return WP2_STATUS_BITSTREAM_ERROR;
    }

    // _Tiffmalloc uses a signed type for size.
    tsize_t alloc_size = 0;
    WP2_CHECK_STATUS(MultFitsIn(stride, image_height, &alloc_size));
    WP2_CHECK_STATUS(raster.Alloc((tsize_t)alloc_size, log_level_));

    if (!TIFFReadRGBAImageOriented(tiff, image_width, image_height,
                                   raster.bytes, ORIENTATION_TOPLEFT, 1)) {
      if (log_level_ >= LogLevel::DEFAULT) {
        fprintf(stderr, "Error! Cannot read TIFF RGBA samples.\n");
      }
      return WP2_STATUS_BITSTREAM_ERROR;
    }

    // TIFF data is ABGR. It is also always returned pre-multiplied.
#ifdef WORDS_BIGENDIAN
    TIFFSwabArrayOfLong(raster, width * height);
#endif
    const bool is_unassociated =
        (extra_samples > 0) && (extra_samples_ptr != nullptr) &&
        (extra_samples_ptr[0] == EXTRASAMPLE_UNASSALPHA);
    WP2_CHECK_STATUS(buffer_->Import(
        is_unassociated ? WP2_RGBA_32 : WP2_rgbA_32, image_width, image_height,
        (const uint8_t*)raster.bytes, stride));
    return WP2_STATUS_OK;
  }
};

void ImageReader::SetImplTIFF(ArgbBuffer* const buffer, LogLevel log_level,
                              size_t max_num_pixels) {
  impl_.reset(new (WP2Allocable::nothrow) ImageReaderTIFF(
      data_.bytes, data_.size, buffer, log_level, max_num_pixels));
  if (impl_ == nullptr) status_ = WP2_STATUS_OUT_OF_MEMORY;
}

}  // namespace WP2

#else  // !WP2_HAVE_TIFF

void WP2::ImageReader::SetImplTIFF(WP2::ArgbBuffer* const buffer,
                                   WP2::LogLevel log_level,
                                   size_t max_num_pixels) {
  (void)buffer;
  (void)max_num_pixels;
  if (log_level >= WP2::LogLevel::DEFAULT) {
    fprintf(stderr, "TIFF support not compiled. Please install the libtiff "
                    "development package before building.\n");
  }
  status_ = WP2_STATUS_VERSION_MISMATCH;
}

#endif  // WP2_HAVE_TIFF

// -----------------------------------------------------------------------------
