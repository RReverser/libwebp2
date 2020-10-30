// Copyright 2020 Google LLC
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
//  Same as imageio but can handle formats in custom color space (y4m).
//
// Author: Yannis Guyon (yguyon@google.com)

#include "extras/ccsp_imageio.h"

#include "extras/extras.h"
#include "imageio/image_dec.h"
#include "imageio/image_enc.h"
#include "imageio/imageio_util.h"

namespace WP2 {

//------------------------------------------------------------------------------

CCSPImageReader::CCSPImageReader(const uint8_t* const data, size_t data_size,
                                 YUVPlane* const ccsp_buffer,
                                 CSPMtx* const ccsp_to_rgb,
                                 Metadata* const metadata,
                                 FileFormat data_format, LogLevel log_level,
                                 size_t max_num_pixels)
    : data_{data, data_size},
      ccsp_buffer_(ccsp_buffer),
      ccsp_to_rgb_(ccsp_to_rgb),
      metadata_(metadata) {
  SetCCSPImpl(data_format, log_level, max_num_pixels);
}

//------------------------------------------------------------------------------

// Proxy to copy the output of ImageReader from ArgbBuffer to CCSPBuffer.
class RGBToCCSPImpl : public ImageReader::Impl {
 public:
  RGBToCCSPImpl(const uint8_t* data, size_t data_size, FileFormat data_format,
                LogLevel log_level, size_t max_num_pixels,
                YUVPlane* const ccsp_buffer, CSPMtx* const ccsp_to_rgb,
                Metadata* const metadata)
      : Impl(nullptr, nullptr, 0, log_level, max_num_pixels),  // Unused.
        rgb_buffer_(WP2_Argb_32),
        image_reader_(data, data_size, &rgb_buffer_, data_format, log_level,
                      max_num_pixels),
        ccsp_buffer_(ccsp_buffer),
        ccsp_to_rgb_(ccsp_to_rgb),
        metadata_(metadata) {}

  WP2Status ReadFrame(bool* const is_last,
                      uint32_t* const duration_ms) override {
    const WP2Status status =
        image_reader_.ReadFrame(is_last, duration_ms);
    if (status == WP2_STATUS_OK) {
      CSPTransform csp_transform;  // For API needs.
      constexpr int16_t kV = 1 << CSPTransform::kMtxShift;
      constexpr int16_t kRGBToRGBMatrix[9] = {kV, 0, 0, 0, kV, 0, 0, 0, kV};
      constexpr int16_t kRgbAvg[3] = {0, 0, 0};
      if (!csp_transform.Init(kRGBToRGBMatrix, kRgbAvg)) assert(false);
      WP2_CHECK_STATUS(ccsp_buffer_->Import(
          rgb_buffer_, /*import_alpha=*/rgb_buffer_.HasTransparency(),
          csp_transform, /*resize_if_needed=*/true));
      // Samples are already in Argb 8 bits.
      *ccsp_to_rgb_ = CSPMtx({1, 0, 0, 0, 1, 0, 0, 0, 1}, 0);  // Identity.
      if (metadata_ != nullptr) {
        WP2_CHECK_STATUS(metadata_->CopyFrom(rgb_buffer_.metadata));
      }
    }
    return status;
  }

 private:
  ArgbBuffer rgb_buffer_;
  ImageReader image_reader_;
  YUVPlane* const ccsp_buffer_;
  CSPMtx* const ccsp_to_rgb_;
  Metadata* const metadata_;
};

void CCSPImageReader::SetCCSPImpl(FileFormat data_format, LogLevel log_level,
                                  size_t max_num_pixels) {
  if (data_format == FileFormat::AUTO) {
    data_format = GuessImageFormat(data_.bytes, data_.size);
  }

  if (data_format == FileFormat::Y4M_420 ||
      data_format == FileFormat::Y4M_444) {
    SetImplY4M(ccsp_buffer_, ccsp_to_rgb_, metadata_, log_level,
               max_num_pixels);
  } else if (data_format == FileFormat::WP2) {
    SetImplWP2(ccsp_buffer_, ccsp_to_rgb_, metadata_, log_level,
               max_num_pixels);
  } else {
    impl_.reset(new (WP2Allocable::nothrow) RGBToCCSPImpl(
        data_.bytes, data_.size, data_format, log_level, max_num_pixels,
        ccsp_buffer_, ccsp_to_rgb_, metadata_));
    if (impl_ == nullptr) status_ = WP2_STATUS_OUT_OF_MEMORY;
  }
}

//------------------------------------------------------------------------------

WP2Status CCSPImageReader::ReadFrame(bool* const is_last,
                                     uint32_t* const duration_ms) {
  // Fill values in case of error.
  if (is_last != nullptr) *is_last = true;
  if (duration_ms != nullptr) *duration_ms = ImageReader::kInfiniteDuration;

  WP2_CHECK_STATUS(status_);
  // Return an error if there is no more frame to decode. Otherwise it could be
  // misleading: whether the last frame should be displayed again or not.
  WP2_CHECK_OK(!is_done_, WP2_STATUS_INVALID_PARAMETER);
  bool is_last_frame = false;  // Initialized to something in case of error.
  uint32_t frame_duration_ms = 0;
  status_ = impl_->ReadFrame(&is_last_frame, &frame_duration_ms);
  if (status_ == WP2_STATUS_OK) {
    if (is_last != nullptr) *is_last = is_last_frame;
    if (duration_ms != nullptr) *duration_ms = frame_duration_ms;
    if (is_last_frame) is_done_ = true;
  }
  return status_;
}

uint32_t CCSPImageReader::GetLoopCount() const {
  return (impl_ != nullptr && status_ == WP2_STATUS_OK)
             ? impl_->GetLoopCount()
             : ImageReader::kInfiniteLoopCount;
}

//------------------------------------------------------------------------------

WP2Status ReadImage(const uint8_t* const data, size_t data_size,
                    YUVPlane* const ccsp_buffer, CSPMtx* const ccsp_to_rgb,
                    Metadata* const metadata, FileFormat format,
                    LogLevel log_level) {
  CCSPImageReader image_reader(data, data_size, ccsp_buffer, ccsp_to_rgb,
                               metadata, format, log_level);
  return image_reader.ReadFrame();
}

WP2Status ReadImage(const char* const file_path, YUVPlane* const ccsp_buffer,
                    CSPMtx* const ccsp_to_rgb, Metadata* const metadata,
                    size_t* const file_size, FileFormat format,
                    LogLevel log_level) {
  Data data;
  WP2_CHECK_STATUS(IoUtilReadFile(file_path, &data));
  if (file_size != nullptr) *file_size = data.size;
  return ReadImage(data.bytes, data.size, ccsp_buffer, ccsp_to_rgb, metadata,
                   format, log_level);
}

static constexpr uint32_t kUnknownBitDepth = 0;
static uint32_t GetY4MBitDepth(const char* const header, size_t header_size) {
  const char* tag = strnstr(header, header_size, " C420");
  if (tag == nullptr) tag = strnstr(header, header_size, " C444");
  if (tag == nullptr) return kUnknownBitDepth;
  tag += 5;
  header_size -= tag - header;
  if (tag[0] == ' ' || tag[0] == '\n' || tag[0] == '\r') return 8;
  if (strlcmp(tag, "jpeg", header_size) == 0) return 8;  // C***jpeg is 8 bits
  if (strlcmp(tag, "p8", header_size) == 0) return 8;    // unlikely
  if (strlcmp(tag, "p9", header_size) == 0) return 9;    // unlikely
  if (strlcmp(tag, "p10", header_size) == 0) return 10;
  if (strlcmp(tag, "p11", header_size) == 0) return 11;  // unlikely
  if (strlcmp(tag, "p12", header_size) == 0) return 12;
  return kUnknownBitDepth;
}

WP2Status ReadBitDepth(const uint8_t* const data, size_t data_size,
                       uint32_t* const bit_depth) {
  WP2_CHECK_OK(data != nullptr, WP2_STATUS_NULL_PARAMETER);
  WP2_CHECK_OK(data_size >= 12, WP2_STATUS_NOT_ENOUGH_DATA);
  WP2_CHECK_OK(bit_depth != nullptr, WP2_STATUS_NULL_PARAMETER);
  const FileFormat format = GuessImageFormat(data, data_size);
  WP2_CHECK_OK(format != FileFormat::UNSUPPORTED,
               WP2_STATUS_UNSUPPORTED_FEATURE);
  if (IsCustomColorSpace(format)) {
    *bit_depth = GetY4MBitDepth(reinterpret_cast<const char*>(data),
                                std::min(data_size, (size_t)64));
    WP2_CHECK_OK(*bit_depth != kUnknownBitDepth,
                 WP2_STATUS_UNSUPPORTED_FEATURE);
  } else {
    *bit_depth = 8;
  }
  return WP2_STATUS_OK;
}

WP2Status ReadBitDepth(const char* const file_path, uint32_t* const bit_depth) {
  Data data;
  WP2_CHECK_STATUS(IoUtilReadFile(file_path, &data));
  return ReadBitDepth(data.bytes, data.size, bit_depth);
}

//------------------------------------------------------------------------------

WP2Status SaveImage(const YUVPlane& buffer, const CSPMtx& ccsp_to_rgb,
                    uint32_t file_num_bits, const char* file_path,
                    bool overwrite, FileFormat format, const Metadata& metadata,
                    const SamplingTaps& filter, size_t* const file_size) {
  WP2_CHECK_OK(file_path != nullptr, WP2_STATUS_NULL_PARAMETER);
  WP2_CHECK_OK(!buffer.IsEmpty(), WP2_STATUS_INVALID_PARAMETER);
  const bool use_stdout = (strcmp(file_path, "-") == 0);
  if (use_stdout) {  // std::ftell() does not work with pipes.
    WP2_CHECK_OK(file_size == nullptr, WP2_STATUS_INVALID_PARAMETER);
  }

  if (format == FileFormat::AUTO) {
    format =
        use_stdout ? FileFormat::Y4M_444 : GetFormatFromExtension(file_path);
  }
  WP2_CHECK_OK(format != FileFormat::UNSUPPORTED,
               WP2_STATUS_UNSUPPORTED_FEATURE);

  WP2_CHECK_OK(ccsp_to_rgb.shift <= 16, WP2_STATUS_INVALID_PARAMETER);

  if (!IsCustomColorSpace(format)) {  // Output Argb samples.
    WP2_CHECK_OK(file_num_bits == 8, WP2_STATUS_INVALID_COLORSPACE);
    ArgbBuffer argb_buffer;
    WP2_CHECK_STATUS(buffer.Export(ccsp_to_rgb, /*resize_if_needed=*/true,
                                   &argb_buffer, &SamplingTaps::kUpSmooth));
    WP2_CHECK_STATUS(argb_buffer.metadata.CopyFrom(metadata));
    WP2_CHECK_STATUS(SaveImage(argb_buffer, file_path, overwrite, format));
    return WP2_STATUS_OK;
  }

  const char* const mode = (overwrite ? "wb" : "wxb");
  FILE* const fout =
      use_stdout ? IoUtilSetBinaryMode(stdout) : fopen(file_path, mode);
  WP2_CHECK_OK(fout != nullptr, (!use_stdout && FileExists(file_path))
                                    ? WP2_STATUS_BAD_WRITE
                                    : WP2_STATUS_INVALID_PARAMETER);

  WP2Status status;
  if (format == FileFormat::Y4M_420 || format == FileFormat::Y4M_444) {
    if (!metadata.IsEmpty()) {
      status = WP2_STATUS_INVALID_CONFIGURATION;
    } else {
      const SamplingTaps* const downsampling =
          (format == FileFormat::Y4M_420) ? &filter : nullptr;
      YUVPlane ycbcr;
      status =
          ToYCbCr(buffer, ccsp_to_rgb, file_num_bits, downsampling, &ycbcr);
      if (status == WP2_STATUS_OK) {
        status = WriteY4M(ycbcr, file_num_bits, fout);
        if (file_size != nullptr) *file_size = std::ftell(fout);
      }
    }
  } else {
    status = WP2_STATUS_UNSUPPORTED_FEATURE;
  }

  if (fout != nullptr && fout != stdout) fclose(fout);
  return status;
}

//------------------------------------------------------------------------------

WP2Status ToRGBA(const YUVPlane& ycbcr,
                 const SamplingTaps* const upsample, ArgbBuffer* const argb) {
  WP2_CHECK_OK(argb != nullptr, WP2_STATUS_NULL_PARAMETER);
  return ycbcr.Export(CSPMtx(kYCbCrToRGBMatrix, kYCbCrToRGBShift),
                      /*resize_if_needed=*/true, argb,
                      ycbcr.IsDownsampled() ? upsample : nullptr);
}

}  // namespace WP2
