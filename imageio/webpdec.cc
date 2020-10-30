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
// WebP decode.
//
// Author: Skal (pascal.massimino@gmail.com)

#include <cstdio>

#include "./anim_image_dec.h"

#ifdef HAVE_CONFIG_H
#include "wp2/config.h"
#endif

#ifdef WP2_HAVE_WEBP
#include "./imageio_util.h"
#include "webp/decode.h"
#include "webp/demux.h"

namespace {

// -----------------------------------------------------------------------------
// Format conversion

WEBP_CSP_MODE FormatToWebPMode(WP2SampleFormat format) {
  switch (format) {
    case WP2_Argb_32:
      return MODE_Argb;
    case WP2_ARGB_32:
      return MODE_ARGB;
    case WP2_rgbA_32:
      return MODE_rgbA;
    case WP2_RGBA_32:
      return MODE_RGBA;
    case WP2_bgrA_32:
      return MODE_bgrA;
    case WP2_BGRA_32:
      return MODE_BGRA;
    default:
      // Warning: WP2_RGB_24 is not the same as MODE_RGB (565).
      return MODE_LAST;
  }
}

// -----------------------------------------------------------------------------
// Metadata

WP2Status Copy(const WebPChunkIterator& iter, WP2::Data* const dst) {
  return dst->CopyFrom((const uint8_t*)iter.chunk.bytes, iter.chunk.size);
}

WP2Status ExtractMetadata(const uint8_t* const data, size_t data_size,
                          WP2::Metadata* const metadata) {
  WebPData webp_data = {data, data_size};
  WebPDemuxer* const demux = WebPDemux(&webp_data);
  if (demux == NULL) return WP2_STATUS_OUT_OF_MEMORY;

  const uint32_t flags = WebPDemuxGetI(demux, WEBP_FF_FORMAT_FLAGS);
  WP2Status status = WP2_STATUS_OK;
  WebPChunkIterator chunk_iter;
  if ((flags & ICCP_FLAG) && WebPDemuxGetChunk(demux, "ICCP", 1, &chunk_iter)) {
    status = Copy(chunk_iter, &metadata->iccp);
    WebPDemuxReleaseChunkIterator(&chunk_iter);
    if (status != WP2_STATUS_OK) goto End;
  }
  if ((flags & EXIF_FLAG) && WebPDemuxGetChunk(demux, "EXIF", 1, &chunk_iter)) {
    status = Copy(chunk_iter, &metadata->exif);
    WebPDemuxReleaseChunkIterator(&chunk_iter);
    if (status != WP2_STATUS_OK) goto End;
  }
  if ((flags & XMP_FLAG) && WebPDemuxGetChunk(demux, "XMP ", 1, &chunk_iter)) {
    status = Copy(chunk_iter, &metadata->xmp);
    WebPDemuxReleaseChunkIterator(&chunk_iter);
    if (status != WP2_STATUS_OK) goto End;
  }
 End:
  WebPDemuxDelete(demux);
  return status;
}

//------------------------------------------------------------------------------
// WebP decoding

}  // namespace

namespace WP2 {

class ImageReaderWEBP : public ImageReader::Impl {
 public:
  ImageReaderWEBP(const uint8_t* data, size_t data_size,
                  ArgbBuffer* const buffer, LogLevel log_level,
                  size_t max_num_pixels)
      : ImageReader::Impl(buffer, data, data_size, log_level, max_num_pixels) {}
  ~ImageReaderWEBP() override {
    WebPAnimDecoderDelete((WebPAnimDecoder*)decoder_);
  }

  WP2Status ReadFrame(bool* const is_last,
                      uint32_t* const duration_ms) override {
    if (decoder_ == nullptr) {
      WP2Status status = WP2_STATUS_OK;
      WebPDecoderConfig config;
      WebPDecBuffer* const output_buffer = &config.output;
      WebPBitstreamFeatures* const bitstream = &config.input;

      WP2_CHECK_STATUS(CheckData());

      if (!WebPInitDecoderConfig(&config)) {
        return WP2_STATUS_VERSION_MISMATCH;
      }

      if (WebPGetFeatures(data_, data_size_, bitstream) != VP8_STATUS_OK) {
        status = WP2_STATUS_BITSTREAM_ERROR;
        goto End;
      }
      WP2_CHECK_STATUS(CheckDimensions((uint32_t)bitstream->width,
                                       (uint32_t)bitstream->height));

      status = buffer_->Resize((uint32_t)bitstream->width,
                               (uint32_t)bitstream->height);
      if (status != WP2_STATUS_OK) goto End;

      if (!bitstream->has_animation) {
        // Not an animation, regular decoding.
        output_buffer->colorspace = FormatToWebPMode(buffer_->format);
        if (output_buffer->colorspace == MODE_LAST) {
          status = WP2_STATUS_INVALID_COLORSPACE;
          goto End;
        }
        output_buffer->u.RGBA.rgba = (uint8_t*)buffer_->GetRow(0);
        output_buffer->u.RGBA.stride = buffer_->stride;
        output_buffer->u.RGBA.size = buffer_->stride * buffer_->height;
        output_buffer->is_external_memory = 1;

        if (WebPDecode(data_, data_size_, &config) != VP8_STATUS_OK) {
          status = WP2_STATUS_BITSTREAM_ERROR;
          goto End;
        }
        *is_last = true;
        *duration_ms = ImageReader::kInfiniteDuration;
      } else {
        // Animation, decode first frame.
        const WebPData webp_data = {data_, data_size_};
        decoder_ = WebPAnimDecoderNew(&webp_data, NULL);
        if (decoder_ == NULL) {
          status = WP2_STATUS_OUT_OF_MEMORY;
          goto End;
        }

        WebPAnimInfo info;
        if (WebPAnimDecoderGetInfo(decoder_, &info) == 0) {
          status = WP2_STATUS_BITSTREAM_ERROR;
          goto End;
        }
        loop_count_ =
            (info.loop_count == 0) ? WP2::ImageReader::kInfiniteLoopCount
                                   : info.loop_count;

        status = DecodeFrame(is_last, duration_ms);
      }

      // Exif/XMP/ICC extraction.
      if (status == WP2_STATUS_OK) {
        const WP2Status status_metadata =
            ExtractMetadata(data_, data_size_, &buffer_->metadata);
        if (status_metadata != WP2_STATUS_OK) status = status_metadata;
      }

     End:
      WebPFreeDecBuffer(output_buffer);
      if (status != WP2_STATUS_OK) {
        if (log_level_ >= LogLevel::DEFAULT) {
          fprintf(stderr, "Decoding failed.\nStatus: %d (%s)\n",
                  status, WP2GetStatusMessage(status));
        }
        buffer_->Deallocate();
      }
      return status;
    } else {
      // More frames.
      return DecodeFrame(is_last, duration_ms);
    }
  }

 protected:
  WP2Status DecodeFrame(bool* const is_last, uint32_t* const duration_ms) {
    uint8_t* buf = nullptr;
    int timestamp = last_timestamp_;

    // Decode frames until one with a duration > 0.
    while (WebPAnimDecoderHasMoreFrames(decoder_) &&
        timestamp <= last_timestamp_) {
      if (!WebPAnimDecoderGetNext(decoder_, &buf, &timestamp) ||
          timestamp < last_timestamp_) {
        return WP2_STATUS_BITSTREAM_ERROR;
      }
    }
    if (buf == nullptr) return WP2_STATUS_BITSTREAM_ERROR;  // No data = error.

    WP2_CHECK_STATUS(buffer_->Import(WP2_RGBA_32, buffer_->width,
                                     buffer_->height, buf, buffer_->width * 4));
    *is_last = !WebPAnimDecoderHasMoreFrames(decoder_);
    *duration_ms = (uint32_t)(timestamp - last_timestamp_);
    last_timestamp_ = timestamp;
    return WP2_STATUS_OK;
  }

  WebPAnimDecoder* decoder_ = nullptr;
  int last_timestamp_ = 0;
};

void ImageReader::SetImplWEBP(ArgbBuffer* const buffer,
                              LogLevel log_level, size_t max_num_pixels) {
  impl_.reset(new (WP2Allocable::nothrow) ImageReaderWEBP(
      data_.bytes, data_.size, buffer, log_level, max_num_pixels));
  if (impl_ == nullptr) status_ = WP2_STATUS_OUT_OF_MEMORY;
}

}  // namespace WP2

#else  // !WP2_HAVE_WEBP

void WP2::ImageReader::SetImplWEBP(WP2::ArgbBuffer* const buffer,
                                   WP2::LogLevel log_level,
                                   size_t max_num_pixels) {
  (void)buffer;
  (void)max_num_pixels;
  if (log_level >= WP2::LogLevel::DEFAULT) {
    fprintf(stderr, "WEBP support not compiled. Please install the libwebp "
                    "development package before building.\n");
  }
  status_ = WP2_STATUS_VERSION_MISMATCH;
}

#endif  // WP2_HAVE_WEBP

// -----------------------------------------------------------------------------
