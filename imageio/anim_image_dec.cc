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
//  ImageReader implementation.
//
// Author: Yannis Guyon (yguyon@google.com)

#include "./anim_image_dec.h"

#include "./image_dec.h"
#include "src/wp2/base.h"

namespace WP2 {

//------------------------------------------------------------------------------

void ImageReader::SetImpl(FileFormat data_format, LogLevel log_level,
                          size_t max_num_pixels, ArgbBuffer* const buffer) {
  if (data_format == FileFormat::AUTO) {
    data_format = GuessImageFormat(data_.bytes, data_.size);
  }

  if (data_format == FileFormat::PNG) {
    SetImplPNG(buffer, log_level, max_num_pixels);
  } else if (data_format == FileFormat::JPEG) {
    SetImplJPEG(buffer, log_level, max_num_pixels);
  } else if ((data_format == FileFormat::PAM) ||
             (data_format == FileFormat::PGM) ||
             (data_format == FileFormat::PPM)) {
    SetImplPNM(buffer, log_level, max_num_pixels);
  } else if (data_format == FileFormat::TIFF) {
    SetImplTIFF(buffer, log_level, max_num_pixels);
  } else if (data_format == FileFormat::GIF) {
    SetImplGIF(buffer, log_level, max_num_pixels);
  } else if (data_format == FileFormat::WEBP) {
    SetImplWEBP(buffer, log_level, max_num_pixels);
  } else if (data_format == FileFormat::WP2) {
    SetImplWP2(buffer, log_level, max_num_pixels);
  } else if (data_format == FileFormat::BMP) {
    SetImplBMP(buffer, log_level, max_num_pixels);
  } else {
    status_ = WP2_STATUS_UNSUPPORTED_FEATURE;
  }
  if (buffer != nullptr) buffer->Deallocate();
}

//------------------------------------------------------------------------------

ImageReader::ImageReader(const char* const file_path, ArgbBuffer* const buffer,
                         FileFormat data_format, LogLevel log_level,
                         size_t max_num_pixels)
    : data_{nullptr, 0} {
  const WP2Status file_status = IoUtilReadFile(file_path, &file_bytes_);
  if (file_status == WP2_STATUS_OK) {
    data_ = {file_bytes_.bytes, file_bytes_.size};
    SetImpl(data_format, log_level, max_num_pixels, buffer);
  } else {
    status_ = file_status;
  }
}

ImageReader::ImageReader(const uint8_t* const data, size_t data_size,
                         ArgbBuffer* const buffer, FileFormat data_format,
                         LogLevel log_level, size_t max_num_pixels)
    : data_{data, data_size} {
  SetImpl(data_format, log_level, max_num_pixels, buffer);
}

ImageReader::ImageReader(const std::string& data, ArgbBuffer* const buffer,
                         FileFormat data_format, LogLevel log_level,
                         size_t max_num_pixels)
    : data_{(const uint8_t*)data.data(), data.size()} {
  SetImpl(data_format, log_level, max_num_pixels, buffer);
}

//------------------------------------------------------------------------------

WP2Status ImageReader::ReadFrame(bool* const is_last,
                                 uint32_t* const duration_ms) {
  // Fill values in case of error.
  if (is_last != nullptr) *is_last = true;
  if (duration_ms != nullptr) *duration_ms = kInfiniteDuration;

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

uint32_t ImageReader::GetLoopCount() const {
  return (impl_ != nullptr && status_ == WP2_STATUS_OK) ? impl_->GetLoopCount()
                                                        : kInfiniteLoopCount;
}

WP2Status ImageReader::Impl::CheckDimensions(uint32_t width,
                                             uint32_t height) const {
  WP2_CHECK_OK(width > 0 && height > 0, WP2_STATUS_BAD_DIMENSION);
  WP2_CHECK_OK(width <= kMaxBufferDimension, WP2_STATUS_BAD_DIMENSION);
  WP2_CHECK_OK(height <= kMaxBufferDimension, WP2_STATUS_BAD_DIMENSION);
  WP2_CHECK_OK((uint64_t)width * height <= std::min((uint64_t)kMaxBufferArea,
                                                    (uint64_t)max_num_pixels_),
               WP2_STATUS_BAD_DIMENSION);
  return WP2_STATUS_OK;
}

WP2Status ImageReader::Impl::CheckData() const {
  WP2_CHECK_OK(data_ != nullptr, WP2_STATUS_NULL_PARAMETER);
  WP2_CHECK_OK(data_size_ > 0, WP2_STATUS_NOT_ENOUGH_DATA);
  WP2_CHECK_OK(buffer_ != nullptr, WP2_STATUS_NULL_PARAMETER);
  return WP2_STATUS_OK;
}

//------------------------------------------------------------------------------

}  // namespace WP2
