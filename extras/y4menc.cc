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
// y4m (raw YCbCr bytes) encoding.

#include <cstdio>

#include "extras/ccsp_imageio.h"
#include "extras/extras.h"
#include "src/dsp/math.h"

namespace WP2 {

//------------------------------------------------------------------------------

WP2Status ToYCbCr(const YUVPlane& ccsp, const CSPMtx& ccsp_to_rgb,
                  uint32_t ycbcr_num_bits, const SamplingTaps* const downsample,
                  YUVPlane* const ycbcr) {
  WP2_CHECK_OK(!ccsp.IsEmpty(), WP2_STATUS_INVALID_PARAMETER);
  WP2_CHECK_OK(ccsp_to_rgb.shift <= 16, WP2_STATUS_INVALID_PARAMETER);

  WP2_CHECK_OK(ycbcr_num_bits >= 8 && ycbcr_num_bits <= 12,
               WP2_STATUS_INVALID_PARAMETER);

  // Combine CCSP->RGB and RGB->YCbCr into CCSP->YCbCr.
  int32_t ccsp_to_ycbcr_matrix[9] = {};
  Multiply(kRGBToYCbCrMatrix, ccsp_to_rgb.mtx(), ccsp_to_ycbcr_matrix);
  const uint32_t ccsp_to_ycbcr_shift =
      ccsp_to_rgb.shift + kRGBToYCbCrShift - (ycbcr_num_bits - 8);

  // TODO(yguyon): Bypass sampling and matrix multiplication if it is ~identity.
  YUVPlane source;
  if (ccsp.IsDownsampled()) {
    WP2_CHECK_STATUS(ycbcr->UpsampleFrom(
        ccsp, (downsample != nullptr && *downsample == SamplingTaps::kDownAvg)
                  ? SamplingTaps::kUpNearest
                  : SamplingTaps::kUpSmooth));
    WP2_CHECK_STATUS(source.SetView(*ycbcr));
  } else {
    WP2_CHECK_STATUS(ycbcr->Resize(ccsp.GetWidth(), ccsp.GetHeight()));
    WP2_CHECK_STATUS(source.SetView(ccsp));
  }

  for (uint32_t y = 0; y < source.GetHeight(); ++y) {
    for (uint32_t x = 0; x < source.GetWidth(); ++x) {
      // Internal precision is 64 bits as input samples are 16b and matrix 32b.
      Multiply<int16_t, int32_t, int16_t, int64_t>(
          source.Y.At(x, y), source.U.At(x, y), source.V.At(x, y),
          ccsp_to_ycbcr_matrix, ccsp_to_ycbcr_shift,
          &ycbcr->Y.At(x, y), &ycbcr->U.At(x, y), &ycbcr->V.At(x, y));
    }
  }

  if (downsample != nullptr) WP2_CHECK_STATUS(ycbcr->Downsample(*downsample));
  return WP2_STATUS_OK;
}

//------------------------------------------------------------------------------

WP2Status ToYCbCr(const ArgbBuffer& rgb, uint32_t ycbcr_num_bits,
                  const SamplingTaps* const downsample, YUVPlane* const ycbcr) {
  WP2_CHECK_OK(!rgb.IsEmpty(), WP2_STATUS_INVALID_PARAMETER);

  WP2_CHECK_OK(ycbcr_num_bits >= 8 && ycbcr_num_bits <= 12,
               WP2_STATUS_INVALID_PARAMETER);
  const uint32_t shift = kRGBToYCbCrShift - (ycbcr_num_bits - 8);

  WP2_CHECK_STATUS(ycbcr->Resize(rgb.width, rgb.height));
  for (uint32_t y = 0; y < rgb.height; ++y) {
    const uint8_t* const row = (const uint8_t*)rgb.GetRow(y);
    for (uint32_t x = 0; x < rgb.width; ++x) {
      const uint8_t* argb = &row[x * WP2FormatBpp(rgb.format)];
      Multiply(argb[1], argb[2], argb[3], kRGBToYCbCrMatrix, shift,
               &ycbcr->Y.At(x, y), &ycbcr->U.At(x, y), &ycbcr->V.At(x, y));
    }
  }

  if (downsample != nullptr) WP2_CHECK_STATUS(ycbcr->Downsample(*downsample));
  return WP2_STATUS_OK;
}

//------------------------------------------------------------------------------

WP2Status WriteY4M(const YUVPlane& buffer, uint32_t num_bits,
                   FILE* const fout) {
  WP2_CHECK_OK(fout != nullptr, WP2_STATUS_NULL_PARAMETER);
  WP2_CHECK_OK(!buffer.IsEmpty(), WP2_STATUS_INVALID_PARAMETER);
  WP2_CHECK_OK(num_bits >= CCSPImageReader::kMinBitDepth &&
                   num_bits <= CCSPImageReader::kMaxBitDepth,
               WP2_STATUS_INVALID_PARAMETER);

  const char* const color_space_tag = buffer.IsDownsampled() ? "C420" : "C444";
  constexpr const char* const kBitDepthTag[] = {"", "p9", "p10", "p11", "p12"};
  STATIC_ASSERT_ARRAY_SIZE(kBitDepthTag, CCSPImageReader::kMaxBitDepth -
                                             CCSPImageReader::kMinBitDepth + 1);
  const char* const bit_depth_tag =
      kBitDepthTag[num_bits - CCSPImageReader::kMinBitDepth];

  fprintf(fout, "YUV4MPEG2 W%u H%u F25:1 Ip A1:1 %s%s\nFRAME\n", buffer.Y.w_,
          buffer.Y.h_, color_space_tag, bit_depth_tag);

  // 'buffer' should be in YCbCr color space from unsigned alpha-premultiplied
  // RGB samples, with meaningful values of Y in [16:235] and Cb, Cr in
  // [-128+16:-128+240] (constants for a bit depth of 8).
  const uint32_t shift = num_bits - 8;
  const int32_t max = (1 << num_bits) - 1;
  const int32_t offset[3] = {16 << shift, 128 << shift, 128 << shift};

  for (Channel channel : {kYChannel, kUChannel, kVChannel}) {
    const Plane16& plane = buffer.GetChannel(channel);
    Vector_u8 row;
    WP2_CHECK_ALLOC_OK(row.resize(plane.w_ * (num_bits > 8 ? 2 : 1)));
    for (uint32_t y = 0; y < plane.h_; ++y) {
      const int16_t* const src = plane.Row(y);
      // Offset to unsigned before saving.
      if (num_bits <= 8) {
        for (uint32_t x = 0; x < plane.w_; ++x) {
          row[x] = (uint8_t)Clamp(src[x] + offset[channel], 0, max);
        }
      } else {
        for (uint32_t x = 0; x < plane.w_; ++x) {
          const uint16_t v = (uint16_t)Clamp(src[x] + offset[channel], 0, max);
          row[2 * x + 0] = (v >> 0) & 0xff;
          row[2 * x + 1] = (v >> 8) & 0xff;
        }
      }
      // write row-by-row, for speed
      WP2_CHECK_OK(fwrite(&row[0], row.size(), 1, fout) == 1,
                          WP2_STATUS_BAD_WRITE);
    }
  }

  return WP2_STATUS_OK;
}

//------------------------------------------------------------------------------

}  // namespace WP2
