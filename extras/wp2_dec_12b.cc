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
// WP2 decode into RGB 12 bits to keep as much precision as possible.

#include "extras/ccsp_imageio.h"
#include "extras/extras.h"
#include "src/utils/data_source.h"
#include "src/utils/utils.h"
#include "src/wp2/decode.h"

namespace WP2 {

//------------------------------------------------------------------------------

// Identity matrix but keep 12 bits of precision (unsigned) compared to 8.
constexpr CSPMtx kRGBToRGB12({16, 0, 0, 0, 16, 0, 0, 0, 16}, 0);

class ImageReaderWP212b : public ImageReader::Impl {
 public:
  ImageReaderWP212b(const uint8_t* data, size_t data_size,
                    YUVPlane* const ccsp_buffer, CSPMtx* const ccsp_to_rgb,
                    Metadata* const metadata, LogLevel log_level,
                    size_t max_num_pixels)
      : ImageReader::Impl(nullptr, data, data_size, log_level, max_num_pixels),
        ccsp_buffer_(ccsp_buffer),
        ccsp_to_rgb_(ccsp_to_rgb),
        metadata_(metadata) {}

  WP2Status ReadFrame(bool* const is_last,
                      uint32_t* const duration_ms) override {
    // Check width/height before allocating anything.
    BitstreamFeatures features;
    WP2_CHECK_STATUS(features.Read(data_, data_size_));
    // TODO(yguyon): Is it better to output only the first frame or an error?
    WP2_CHECK_OK(!features.is_animation, WP2_STATUS_UNSUPPORTED_FEATURE);
    *is_last = true;
    *duration_ms = ImageReader::kInfiniteDuration;
    WP2_CHECK_STATUS(CheckDimensions(features.width, features.height));
    loop_count_ = ImageReader::kInfiniteLoopCount;

    // Allocate a padded buffer to avoid doing it internally.
    // TODO(yguyon): Check if 'has_alpha' is needed beforehand (decode GParams).
    WP2_CHECK_STATUS(ccsp_buffer_->Resize(features.width, features.height,
                                          /*pad=*/kPredWidth,
                                          /*has_alpha=*/true));

    // Decode.
    WP2_CHECK_STATUS(Decode(data_, data_size_,
      kRGBToRGB12.mtx(), kRGBToRGB12.shift,
      ccsp_buffer_->Y.Row(0), ccsp_buffer_->Y.Step(), ccsp_buffer_->Y.Size(),
      ccsp_buffer_->U.Row(0), ccsp_buffer_->U.Step(), ccsp_buffer_->U.Size(),
      ccsp_buffer_->V.Row(0), ccsp_buffer_->V.Step(), ccsp_buffer_->V.Size(),
      ccsp_buffer_->HasAlpha() ? ccsp_buffer_->A.Row(0) : nullptr,
      ccsp_buffer_->A.Step(), ccsp_buffer_->A.Size(), metadata_));

    bool is_opaque = true;
    for (uint32_t y = 0; y < ccsp_buffer_->A.h_ && is_opaque; ++y) {
      for (uint32_t x = 0; x < ccsp_buffer_->A.w_ && is_opaque; ++x) {
        if (ccsp_buffer_->A.At(x, y) != kAlphaMax) is_opaque = false;
      }
    }
    if (is_opaque) ccsp_buffer_->A.Clear();
    // Remove padding without deallocating.
    WP2_CHECK_STATUS(ccsp_buffer_->SetView(
        *ccsp_buffer_, {0, 0, features.width, features.height}));

    *ccsp_to_rgb_ = CSPMtx({1, 0, 0, 0, 1, 0, 0, 0, 1}, 4);  // Go back to 8b.
    return WP2_STATUS_OK;
  }

 protected:
  // Pointers to user structs.
  YUVPlane* const ccsp_buffer_;
  CSPMtx* const ccsp_to_rgb_;
  Metadata* const metadata_;
};

void CCSPImageReader::SetImplWP2(YUVPlane* const buffer,
                                 CSPMtx* const ccsp_to_rgb,
                                 Metadata* const metadata, LogLevel log_level,
                                 size_t max_num_pixels) {
  impl_.reset(new (WP2Allocable::nothrow) ImageReaderWP212b(
      data_.bytes, data_.size, buffer, ccsp_to_rgb, metadata, log_level,
      max_num_pixels));
  if (impl_ == nullptr) status_ = WP2_STATUS_OUT_OF_MEMORY;
}

//------------------------------------------------------------------------------

}  // namespace WP2
