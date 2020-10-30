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
// JPEG decode.
//
// Author: Skal (pascal.massimino@gmail.com)

#include <cstdio>

#ifdef HAVE_CONFIG_H
#include "src/wp2/config.h"
#endif

#include "./anim_image_dec.h"

#ifdef WP2_HAVE_JPEG
#include <jerror.h>
#include <jpeglib.h>
#include <csetjmp>
#include <cstdlib>
#include <cstring>

#include "./imageio_util.h"
#include "src/utils/utils.h"
#include "src/wp2/base.h"

namespace {

constexpr int kNumOutputComponents = 3;

// -----------------------------------------------------------------------------
// Metadata processing

// JPEG application-specific markers used for EXIF, XMP and ICCP.
// See https://en.wikipedia.org/wiki/JPEG
#ifndef JPEG_APP1
#define JPEG_APP1 (JPEG_APP0 + 1)
#endif
#ifndef JPEG_APP2
#define JPEG_APP2 (JPEG_APP0 + 2)
#endif

// Struct used for storing a JPEG metadata signature.
struct MetadataSignature {
  int code;              // JPEG application-specific marker.
  const uint8_t* bytes;  // Signature.
  size_t size;           // Signature length.
  size_t skip_size;      // Signature length plus some header bytes.
};

// Exif 2.2 Section 4.7.2 Interoperability Structure of APP1 ...
const MetadataSignature kEXIFSignature = {
    JPEG_APP1, (const uint8_t*)"Exif\0", 6, 6};

// XMP Specification Part 3 Section 3 Embedding XMP Metadata ... #JPEG
// TODO(jzern) Add support for 'ExtendedXMP'
const MetadataSignature kXMPSignature = {
    JPEG_APP1, (const uint8_t*)"http://ns.adobe.com/xap/1.0/", 29, 29};

// ICC.1:2010-12 (4.3.0.0) Annex B.4 Embedding ICC Profiles in JPEG files
const MetadataSignature kICCPSignature = {
    JPEG_APP2, (const uint8_t*)"ICC_PROFILE", 12, 14};

#undef JPEG_APP1
#undef JPEG_APP2

// Returns true if 'marker_code' and 'marker_data' match the 'signature'.
bool MarkerHasSignature(int marker_code, WP2::DataView marker_data,
                        const MetadataSignature& signature) {
  return (marker_code == signature.code &&
          marker_data.size >= signature.skip_size &&
          !std::memcmp(marker_data.bytes, signature.bytes, signature.size));
}

// Extracts chunk from 'marker_data' to 'metadata'.
WP2Status ExtractMetadataFromMarker(WP2::DataView marker_data,
                                    const MetadataSignature& signature,
                                    WP2::Data* const metadata,
                                    WP2::LogLevel log_level) {
  if (!metadata->IsEmpty()) {
    if (log_level >= WP2::LogLevel::VERBOSE) {
      fprintf(stderr, "Ignoring additional '%s' marker\n", signature.bytes);
    }
    return WP2_STATUS_OK;
  }
  assert(marker_data.size >= signature.size);  // Tested in MarkerHasSignature()
  return metadata->CopyFrom(marker_data.bytes + signature.size,
                            marker_data.size - signature.size);
}

// -----------------------------------------------------------------------------
// ICCP metadata processing

struct ICCPSegment {
  const uint8_t* data;
  size_t data_length;
  size_t seq;  // Sequence number [1, 255] for use in reassembly.
};

// Extracts ICC profile segment from 'marker_data' to 'iccp_segments'.
// If 'iccp_segment_expected_count' is 0, it will be set to the segment's
// extracted count. Otherwise it will be compared to it.
// 'iccp_max_seq' will be set to the segment's extracted sequence number if it's
// the highest seen so far.
// Returns WP2_STATUS_OK on success.
WP2Status StoreICCPFromMarker(WP2::DataView marker_data,
                              ICCPSegment iccp_segments[255],
                              size_t* const iccp_segment_expected_count,
                              size_t* const iccp_max_seq,
                              WP2::LogLevel log_level) {
  // ICC_PROFILE\0<seq><count>; 'seq' starts at 1.
  const size_t seq = marker_data.bytes[kICCPSignature.size];
  const size_t count = marker_data.bytes[kICCPSignature.size + 1];
  const size_t segment_size = marker_data.size - kICCPSignature.skip_size;

  if (segment_size == 0) {
    if (log_level >= WP2::LogLevel::DEFAULT) {
      fprintf(stderr, "[ICCP] size (%zu) cannot be 0!\n", segment_size);
    }
    return WP2_STATUS_BITSTREAM_ERROR;
  } else if (count == 0) {
    if (log_level >= WP2::LogLevel::DEFAULT) {
      fprintf(stderr, "[ICCP] count (%zu) cannot be 0!\n", count);
    }
    return WP2_STATUS_BITSTREAM_ERROR;
  } else if (seq == 0) {
    if (log_level >= WP2::LogLevel::DEFAULT) {
      fprintf(stderr, "[ICCP] invalid sequence (%zu)!\n", seq);
    }
    return WP2_STATUS_BITSTREAM_ERROR;
  }

  if (*iccp_segment_expected_count == 0) {
    *iccp_segment_expected_count = count;
  } else if (*iccp_segment_expected_count != count) {
    if (log_level >= WP2::LogLevel::DEFAULT) {
      fprintf(stderr, "[ICCP] Inconsistent segment count (%zu / %zu)!\n",
              *iccp_segment_expected_count, count);
    }
    return WP2_STATUS_BITSTREAM_ERROR;
  }

  ICCPSegment* const segment = iccp_segments + seq - 1;
  if (segment->data_length != 0) {
    if (log_level >= WP2::LogLevel::DEFAULT) {
      fprintf(stderr, "[ICCP] Duplicate segment number (%zu)!\n", seq);
    }
    return WP2_STATUS_BITSTREAM_ERROR;
  }

  segment->data = marker_data.bytes + kICCPSignature.skip_size;
  segment->data_length = segment_size;
  segment->seq = seq;
  if (seq > *iccp_max_seq) *iccp_max_seq = seq;
  return WP2_STATUS_OK;
}

// Returns WP2_STATUS_OK if everything's fine.
WP2Status VerifyICCP(size_t iccp_segment_count,
                     size_t iccp_segment_expected_count, size_t iccp_max_seq,
                     WP2::LogLevel log_level) {
  if (iccp_segment_expected_count != iccp_segment_count) {
    if (log_level >= WP2::LogLevel::DEFAULT) {
      fprintf(stderr,
              "[ICCP] Segment count: %zu does not match expected: %zu!\n",
              iccp_segment_count, iccp_segment_expected_count);
    }
    return WP2_STATUS_BITSTREAM_ERROR;
  }
  if (iccp_max_seq != iccp_segment_count) {
    if (log_level >= WP2::LogLevel::DEFAULT) {
      fprintf(stderr,
              "[ICCP] Discontinuous segments, expected: %zu actual: %zu!\n",
              iccp_segment_count, iccp_max_seq);
    }
    return WP2_STATUS_BITSTREAM_ERROR;
  }
  return WP2_STATUS_OK;
}

// Merges all ICC profile segments from 'iccp_segments' to 'metadata'.
// Returns WP2_STATUS_OK on success.
WP2Status MergeICCP(ICCPSegment iccp_segments[255],
                    size_t iccp_segment_count, WP2::Data* const metadata) {
  // The segments may appear out of order in the file, sort them based on
  // sequence number before assembling the payload.
  qsort(iccp_segments, iccp_segment_count, sizeof(*iccp_segments),
        [](const void* a, const void* b) {
          const size_t seq_a = ((const ICCPSegment*)a)->seq;
          const size_t seq_b = ((const ICCPSegment*)b)->seq;
          return (seq_a < seq_b) ? -1 : (seq_a > seq_b) ? 1 : 0;
        });

  size_t iccp_total_size = 0;
  for (size_t i = 0; i < iccp_segment_count; ++i) {
    iccp_total_size += iccp_segments[i].data_length;
  }
  WP2_CHECK_STATUS(metadata->Resize(iccp_total_size, /*keep_bytes=*/false));

  uint8_t* dst = metadata->bytes;
  for (size_t i = 0; i < iccp_segment_count; ++i) {
    memcpy(dst, iccp_segments[i].data, iccp_segments[i].data_length);
    dst += iccp_segments[i].data_length;
  }
  return WP2_STATUS_OK;
}

// -----------------------------------------------------------------------------
// Metadata processing

// Returns WP2_STATUS_OK on success.
WP2Status ExtractMetadataFromJPEG(j_decompress_ptr dinfo,
                                  WP2::Metadata* const metadata,
                                  WP2::LogLevel log_level) {
  metadata->Clear();

  // Treat ICC profiles separately as they may be segmented and out of order.
  size_t iccp_segment_count = 0;
  size_t iccp_segment_expected_count = 0;
  size_t iccp_max_seq = 0;
  ICCPSegment iccp_segments[255];
  memset(iccp_segments, 0, sizeof(iccp_segments));

  jpeg_saved_marker_ptr marker;
  for (marker = dinfo->marker_list; marker != NULL; marker = marker->next) {
    const int marker_code = marker->marker;
    const WP2::DataView marker_data = {marker->data, marker->data_length};

    if (MarkerHasSignature(marker_code, marker_data, kEXIFSignature)) {
      WP2_CHECK_STATUS(ExtractMetadataFromMarker(marker_data, kEXIFSignature,
                                                 &metadata->exif, log_level));
    } else if (MarkerHasSignature(marker_code, marker_data, kXMPSignature)) {
      WP2_CHECK_STATUS(ExtractMetadataFromMarker(marker_data, kXMPSignature,
                                                 &metadata->xmp, log_level));
    } else if (MarkerHasSignature(marker_code, marker_data, kICCPSignature)) {
      WP2_CHECK_STATUS(StoreICCPFromMarker(marker_data, iccp_segments,
                                           &iccp_segment_expected_count,
                                           &iccp_max_seq, log_level));
      ++iccp_segment_count;
    } else if (log_level >= WP2::LogLevel::VERBOSE) {
      fprintf(stderr, "Ignoring marker '%d'\n", marker->marker);
    }
  }

  if (iccp_segment_count > 0) {
    WP2_CHECK_STATUS(VerifyICCP(iccp_segment_count, iccp_segment_expected_count,
                                iccp_max_seq, log_level));
    WP2_CHECK_STATUS(
        MergeICCP(iccp_segments, iccp_segment_count, &metadata->iccp));
  }
  return WP2_STATUS_OK;
}

// -----------------------------------------------------------------------------
// Error handling

struct JPEGErrorContext {
  struct jpeg_error_mgr pub;
  jmp_buf setjmp_buffer;
};

void ErrorExitFunctionNoLog(j_common_ptr dinfo) {
  JPEGErrorContext* const myerr = (JPEGErrorContext*)dinfo->err;
  longjmp(myerr->setjmp_buffer, 1);
}

void ErrorExitFunction(j_common_ptr dinfo) {
  dinfo->err->output_message(dinfo);
  ErrorExitFunctionNoLog(dinfo);
}

// -----------------------------------------------------------------------------
// Context

struct JPEGReadContext {
  struct jpeg_source_mgr pub;
  const uint8_t* data;
  size_t data_size;
};

void ContextInit(j_decompress_ptr cinfo) {
  JPEGReadContext* const ctx = (JPEGReadContext*)cinfo->src;
  ctx->pub.next_input_byte = ctx->data;
  ctx->pub.bytes_in_buffer = ctx->data_size;
}

boolean ContextFill(j_decompress_ptr cinfo) {
  // We shouldn't get here.
  ERREXIT(cinfo, JERR_FILE_READ);
  return FALSE;
}

void ContextSkip(j_decompress_ptr cinfo, long jump_size) {  // NOLINT (long)
  JPEGReadContext* const ctx = (JPEGReadContext*)cinfo->src;
  size_t jump = (size_t)jump_size;
  if (jump > ctx->pub.bytes_in_buffer) {  // Don't overflow the buffer.
    jump = ctx->pub.bytes_in_buffer;
  }
  ctx->pub.bytes_in_buffer -= jump;
  ctx->pub.next_input_byte += jump;
}

void ContextTerm(j_decompress_ptr cinfo) {
  (void)cinfo;
}

void ContextSetup(volatile struct jpeg_decompress_struct* const cinfo,
                  JPEGReadContext* const ctx) {
  assert((void*)ctx == (void*)&ctx->pub);
  cinfo->src = &ctx->pub;
  ctx->pub.init_source = ContextInit;
  ctx->pub.fill_input_buffer = ContextFill;
  ctx->pub.skip_input_data = ContextSkip;
  ctx->pub.resync_to_restart = jpeg_resync_to_restart;
  ctx->pub.term_source = ContextTerm;
  ctx->pub.bytes_in_buffer = 0;
  ctx->pub.next_input_byte = NULL;
}

}  // namespace

// -----------------------------------------------------------------------------

namespace WP2 {

class ImageReaderJPEG : public ImageReader::Impl {
 public:
  ImageReaderJPEG(const uint8_t* data, size_t data_size,
                  ArgbBuffer* const buffer, LogLevel log_level,
                  size_t max_num_pixels)
      : ImageReader::Impl(buffer, data, data_size, log_level, max_num_pixels) {}

  WP2Status ReadFrame(bool* const is_last,
                      uint32_t* const duration_ms) override {
    WP2_CHECK_STATUS(CheckData());

    // The following are volatile because they need to keep their value for
    // deletion in case of longjmp() error handling.
    volatile jpeg_decompress_struct dinfo;
    memset((j_decompress_ptr)&dinfo, 0, sizeof(dinfo));  // for setjmp sanity
    uint8_t* volatile tmp_rgb = nullptr;

    JPEGReadContext ctx;
    memset(&ctx, 0, sizeof(ctx));
    ctx.data = data_;
    ctx.data_size = data_size_;

    JPEGErrorContext jerr;
    memset(&jerr, 0, sizeof(jerr));
    dinfo.err = jpeg_std_error(&jerr.pub);
    jerr.pub.error_exit = (log_level_ >= LogLevel::DEFAULT)
                              ? ErrorExitFunction
                              : ErrorExitFunctionNoLog;

    if (setjmp(jerr.setjmp_buffer)) {
      // On a libjpeg error, error_exit() is called, longjmping here.
      jpeg_destroy_decompress((j_decompress_ptr)&dinfo);
      WP2Free(tmp_rgb);
      return WP2_STATUS_BITSTREAM_ERROR;
    }

    jpeg_create_decompress((j_decompress_ptr)&dinfo);

    // Avoid out-of-memory issues.
    if (MultFitsIn(max_num_pixels_, kNumOutputComponents,
                   &dinfo.mem->max_memory_to_use) != WP2_STATUS_OK) {
      // If 'max_num_pixels_ * kNumOutputComponents' does not fit in the
      // 'max_memory_to_use' variable, leave it to its default value so memory
      // is not limited to a lower-than-specified bound.
    }
    ContextSetup(&dinfo, &ctx);

    jpeg_read_header((j_decompress_ptr)&dinfo, TRUE);

    dinfo.out_color_space = JCS_RGB;
    dinfo.do_fancy_upsampling = TRUE;

    WP2Status status = ReadCanvas(&dinfo, &tmp_rgb);
    WP2Free(tmp_rgb);
    tmp_rgb = nullptr;

    if (status == WP2_STATUS_OK) {
      status = ExtractMetadataFromJPEG((j_decompress_ptr)&dinfo,
                                       &buffer_->metadata, log_level_);
    }

    if (status == WP2_STATUS_OK) {
      jpeg_finish_decompress((j_decompress_ptr)&dinfo);
    }
    jpeg_destroy_decompress((j_decompress_ptr)&dinfo);
    if (status == WP2_STATUS_OK) {
      *is_last = true;
      *duration_ms = ImageReader::kInfiniteDuration;
    }
    return status;
  }

 protected:
  // Reads pixels into 'buffer_'.
  // Allocates '*tmp_rgb' which is volatile because allocation happens between
  // setjmp() and longjmp() and needs to be freed afterwards.
  WP2Status ReadCanvas(volatile jpeg_decompress_struct* const dinfo,
                       uint8_t* volatile* const tmp_rgb) {
    jpeg_start_decompress((j_decompress_ptr)dinfo);

    WP2_CHECK_OK(dinfo->output_components == kNumOutputComponents,
                 WP2_STATUS_UNSUPPORTED_FEATURE);

    const uint32_t width = dinfo->output_width;
    const uint32_t height = dinfo->output_height;
    WP2_CHECK_STATUS(CheckDimensions(width, height));

    const int depth = dinfo->output_components * sizeof(**tmp_rgb);
    size_t stride = 0;
    WP2_CHECK_STATUS(MultFitsIn(width, depth, &stride));
    // Make sure whole size fits in size_t.
    WP2_CHECK_STATUS(MultFitsIn<size_t>(stride, height));

    WP2_CHECK_STATUS(buffer_->Resize(width, height));

    *tmp_rgb = (uint8_t*)WP2Malloc(1, stride);
    WP2_CHECK_ALLOC_OK(*tmp_rgb != nullptr);
    JSAMPROW buffer[1];
    buffer[0] = (JSAMPLE*)*tmp_rgb;

    while (dinfo->output_scanline < dinfo->output_height) {
      const uint32_t row = dinfo->output_scanline;
      WP2_CHECK_OK(jpeg_read_scanlines((j_decompress_ptr)dinfo, buffer, 1) == 1,
                   WP2_STATUS_BITSTREAM_ERROR);
      WP2_CHECK_STATUS(buffer_->ImportRow(WP2_RGB_24, row, *tmp_rgb));
    }
    return WP2_STATUS_OK;
  }
};

void ImageReader::SetImplJPEG(ArgbBuffer* const buffer,
                              LogLevel log_level, size_t max_num_pixels) {
  impl_.reset(new (WP2Allocable::nothrow) ImageReaderJPEG(
      data_.bytes, data_.size, buffer, log_level, max_num_pixels));
  if (impl_ == nullptr) status_ = WP2_STATUS_OUT_OF_MEMORY;
}

}  // namespace WP2

#else  // !WP2_HAVE_JPEG

namespace WP2 {

void ImageReader::SetImplJPEG(ArgbBuffer* const buffer, LogLevel log_level,
                              size_t max_num_pixels) {
  (void)buffer;
  (void)max_num_pixels;
  if (log_level >= WP2::LogLevel::DEFAULT) {
    fprintf(stderr, "JPEG support not compiled. Please install the libjpeg "
                    "development package before building.\n");
  }
  status_ = WP2_STATUS_VERSION_MISMATCH;
}

}  // namespace WP2

#endif  // WP2_HAVE_JPEG

// -----------------------------------------------------------------------------
