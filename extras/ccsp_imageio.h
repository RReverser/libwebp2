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

#ifndef WP2_EXTRAS_CCSP_IMAGEIO_H_
#define WP2_EXTRAS_CCSP_IMAGEIO_H_

#include <cstddef>
#include <cstdint>

#include "imageio/anim_image_dec.h"
#include "imageio/file_format.h"
#include "src/utils/csp.h"
#include "src/utils/plane.h"
#include "src/wp2/base.h"

namespace WP2 {

//------------------------------------------------------------------------------

// Same as ImageReader but can also read custom color spaces and/or higher bit
// depths.
class CCSPImageReader {
 public:
  static constexpr uint32_t kMinBitDepth = 8, kMaxBitDepth = 12;

  // 'data' must not change till the last ReadFrame() call.
  // The 'ccsp_to_rgb_matrix' and the 'ccsp_to_rgb_shift' are also returned and
  // can be applied in this order to the output 'ccsp_buffer' samples to obtain
  // RGB. If not null, 'metadata' is filled with EXIF, XMP and/or ICCP data, if
  // any. Clears 'ccsp_buffer'.
  CCSPImageReader(const uint8_t* const data, size_t data_size,
                  YUVPlane* const ccsp_buffer, CSPMtx* const ccsp_to_rgb,
                  Metadata* const metadata = nullptr,
                  FileFormat data_format = FileFormat::AUTO,
                  LogLevel log_level = LogLevel::DEFAULT,
                  size_t max_num_pixels = kMaxBufferArea);

  WP2Status ReadFrame(bool* const is_last = nullptr,
                      uint32_t* const duration_ms = nullptr);
  uint32_t GetLoopCount() const;

 protected:
  void SetImplY4M(YUVPlane* const ccsp_buffer, CSPMtx* const ccsp_to_rgb,
                  Metadata* const metadata, LogLevel log_level,
                  size_t max_num_pixels);
  void SetImplWP2(YUVPlane* const ccsp_buffer, CSPMtx* const ccsp_to_rgb,
                  Metadata* const metadata, LogLevel log_level,
                  size_t max_num_pixels);
  void SetCCSPImpl(FileFormat data_format, LogLevel log_level,
                   size_t max_num_pixels);

  std::unique_ptr<ImageReader::Impl> impl_;
  const DataView data_;  // Input. Points to user memory.
  WP2Status status_ = WP2_STATUS_OK;
  bool is_done_ = false;

  // Final output, points to the user data structs.
  YUVPlane* const ccsp_buffer_;
  CSPMtx* const ccsp_to_rgb_;
  Metadata* const metadata_;
};

//------------------------------------------------------------------------------

// Decodes 'data' into 'ccsp_buffer'. Image type (y4m, png etc.) is guessed from
// the bitstream if 'format' is AUTO.
// The 'ccsp_to_rgb_matrix' and the 'ccsp_to_rgb_shift' are also returned and
// can be applied in this order to the output 'buffer' samples to obtain RGB.
// Returns WP2_STATUS_OK upon success.
WP2Status ReadImage(const uint8_t* const data, size_t data_size,
                    YUVPlane* const ccsp_buffer, CSPMtx* const ccsp_to_rgb,
                    Metadata* const metadata = nullptr,
                    FileFormat format = FileFormat::AUTO,
                    LogLevel log_level = LogLevel::DEFAULT);

// Reads the file at 'file_path' and calls the function above.
// If 'file_size' is not null, the file's original size is returned.
WP2Status ReadImage(const char* const file_path, YUVPlane* const ccsp_buffer,
                    CSPMtx* const ccsp_to_rgb,
                    Metadata* const metadata = nullptr,
                    size_t* file_size = nullptr,
                    FileFormat format = FileFormat::AUTO,
                    LogLevel log_level = LogLevel::DEFAULT);

// Decodes enough 'data' to get the 'bit_depth' of the image.
WP2Status ReadBitDepth(const uint8_t* const data, size_t data_size,
                       uint32_t* const bit_depth);
WP2Status ReadBitDepth(const char* const file_path, uint32_t* const bit_depth);

//------------------------------------------------------------------------------

// The custom color space is specified by a matrix and a right shift applied in
// this order on 'buffer' samples to obtain Argb 8 bits.
// Returns WP2_STATUS_OK or an error. Returns written 'file_size' if not null.
WP2Status SaveImage(const YUVPlane& buffer, const CSPMtx& ccsp_to_rgb,
                    uint32_t file_num_bits, const char* file_path,
                    bool overwrite = false,
                    FileFormat format = FileFormat::AUTO,
                    const Metadata& metadata = Metadata(),
                    const SamplingTaps& filter = SamplingTaps::kDownSharp,
                    size_t* const file_size = nullptr);

// Converts samples to MPEG-2 YCbCr color space.
// Halves Cb and Cr resolution if 'downsample' is not null.
WP2Status ToYCbCr(const YUVPlane& ccsp, const CSPMtx& ccsp_to_rgb,
                  uint32_t ycbcr_num_bits, const SamplingTaps* const downsample,
                  YUVPlane* const ycbcr);
// This specialized version uses kRGBToYCbCrMatrix[] (hence: reduced range)
WP2Status ToYCbCr(const ArgbBuffer& rgb, uint32_t ycbcr_num_bits,
                  const SamplingTaps* const downsample, YUVPlane* const ycbcr);

// Expects 'buffer' in MPEG-2 YCbCr color space.
WP2Status WriteY4M(const YUVPlane& buffer, uint32_t num_bits, FILE* const fout);

//------------------------------------------------------------------------------

// Converts YCbCr content of 'yuv' into 'argb'. If 'ycbcr' is downsampled, the
// filter 'upsample' will be used for upsampling into 'argb'.
// This function is the inverse of ToYCbCr().
WP2Status ToRGBA(const YUVPlane& ycbcr, const SamplingTaps* const upsample,
                 ArgbBuffer* const argb);

//------------------------------------------------------------------------------

}  // namespace WP2

#endif  // WP2_EXTRAS_CCSP_IMAGEIO_H_
