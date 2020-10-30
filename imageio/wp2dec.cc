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
// WP2 decode.

#include <cstdio>

#include "./anim_image_dec.h"
#include "./imageio_util.h"
#include "src/wp2/decode.h"

namespace WP2 {

//------------------------------------------------------------------------------
// WP2 decoding

class ImageReaderWP2 : public ImageReader::Impl {
 public:
  ImageReaderWP2(const uint8_t* data, size_t data_size,
                 ArgbBuffer* const buffer, LogLevel log_level,
                 size_t max_num_pixels)
      : ImageReader::Impl(buffer, data, data_size, log_level, max_num_pixels),
        decoder_(data, data_size, DecoderConfig::kDefault, buffer) {}

  WP2Status ReadFrame(bool* const is_last,
                      uint32_t* const duration_ms) override {
    if (frame_index_ == 0) {
      WP2_CHECK_STATUS(CheckData());

      // Check width/height before allocating anything and get loop count.
      BitstreamFeatures features;
      WP2_CHECK_STATUS(features.Read(data_, data_size_));
      WP2_CHECK_STATUS(CheckDimensions(features.width, features.height));
      loop_count_ = (features.loop_count == WP2::kInfiniteLoop)
                        ? ImageReader::kInfiniteLoopCount
                        : features.loop_count;
    }

    // Attempt to decode the next frame.
    uint32_t frame_duration_ms;
    WP2Status status;
    if (decoder_.ReadFrame(&frame_duration_ms)) {
      status = WP2_STATUS_OK;
      assert(frame_index_ == decoder_.GetCurrentFrameIndex());
      // Features are not null because ReadFrame() returned true.
      const BitstreamFeatures& features =
          *decoder_.TryGetDecodedFeatures();
      const FrameFeatures& frame_features =
          *decoder_.TryGetFrameDecodedFeatures(frame_index_);

      *is_last = frame_features.is_last;
      *duration_ms = features.is_animation ? frame_duration_ms
                                           : ImageReader::kInfiniteDuration;

      if (frame_features.is_last) {
        // Finish decoding: extract metadata if any and if enough bytes.
        if (decoder_.ReadFrame()) assert(false);  // No more frame.
        status = decoder_.GetStatus();
      } else if (frame_index_ == 0) {
        // Metadata must be extracted during the first ReadFrame() call.
        // ICCP was already decoded after header.
        const WP2Status metadata_status =
            ExtractMetadata(data_, data_size_, &buffer_->metadata.exif,
                            /*iccp=*/nullptr, &buffer_->metadata.xmp);
        if (metadata_status != WP2_STATUS_OK) status = metadata_status;
      }
    } else if (decoder_.GetStatus() == WP2_STATUS_OK) {
      // Decoding is done, there is no more frame, ReadFrame() should not have
      // been called.
      status = WP2_STATUS_INVALID_PARAMETER;
    } else {
      // Error or not enough data.
      status = decoder_.GetStatus();
    }

    if (status != WP2_STATUS_OK) {
      if (log_level_ >= LogLevel::DEFAULT) {
        fprintf(stderr, "Decoding failed.\nStatus: %d (%s)\n",
                status, WP2GetStatusMessage(status));
      }
    } else {
      ++frame_index_;
    }
    return status;
  }

 protected:
  ArrayDecoder decoder_;
  uint32_t frame_index_ = 0;
};

void ImageReader::SetImplWP2(ArgbBuffer* const buffer,
                             LogLevel log_level, size_t max_num_pixels) {
  impl_.reset(new (WP2Allocable::nothrow) ImageReaderWP2(
      data_.bytes, data_.size, buffer, log_level, max_num_pixels));
  if (impl_ == nullptr) status_ = WP2_STATUS_OUT_OF_MEMORY;
}

// -----------------------------------------------------------------------------

}  // namespace WP2
