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
// Image savers

#include <cstdint>
#include <cstdio>

#include "./image_enc.h"
#include "src/dsp/dsp.h"
#include "src/utils/utils.h"
#include "src/utils/vector.h"

namespace WP2 {

static WP2Status WritePPMPAM(const ArgbBuffer& buffer, int alpha,
                             FILE* const fout) {
  WP2_CHECK_OK(fout != nullptr, WP2_STATUS_NULL_PARAMETER);
  WP2_CHECK_OK(buffer.format == WP2_Argb_32, WP2_STATUS_UNSUPPORTED_FEATURE);

  const uint32_t width = buffer.width;
  const uint32_t height = buffer.height;
  const uint8_t* row = (const uint8_t*)buffer.GetRow(0);
  const uint32_t stride = buffer.stride;
  const size_t bytes_per_px = alpha ? 4 : 3;
  WP2SampleFormat fmt;
  Vector_u8 dest;
  uint32_t y;

  WP2_CHECK_OK(row != nullptr, WP2_STATUS_NULL_PARAMETER);
  WP2_CHECK_ALLOC_OK(dest.resize(width * bytes_per_px));
  WP2ArgbConverterInit();

  if (alpha) {
    fprintf(fout,
            "P7\nWIDTH %u\nHEIGHT %u\nDEPTH 4\nMAXVAL 255\n"
            "TUPLTYPE RGB_ALPHA\nENDHDR\n",
            width, height);
    fmt = WP2_RGBA_32;
  } else {
    fprintf(fout, "P6\n%u %u\n255\n", width, height);
    fmt = WP2_RGB_24;
  }
  for (y = 0; y < height; ++y) {
    WP2ArgbConvertTo[fmt](row, width, dest.data());
    WP2_CHECK_OK(fwrite(dest.data(), width, bytes_per_px, fout) == bytes_per_px,
                 WP2_STATUS_BAD_WRITE);
    row += stride;
  }
  return WP2_STATUS_OK;
}

WP2Status WritePPM(const ArgbBuffer& buffer, FILE* const fout) {
  return WritePPMPAM(buffer, 0, fout);
}

WP2Status WritePAM(const ArgbBuffer& buffer, FILE* const fout) {
  return WritePPMPAM(buffer, 1, fout);
}

//------------------------------------------------------------------------------
// Raw PGM

// Save 16b mode (RGBA4444, RGB565, ...) for debugging purpose.
WP2Status Write16bAsPGM(const ArgbBuffer& buffer, FILE* const fout) {
  WP2_CHECK_OK(fout != nullptr, WP2_STATUS_NULL_PARAMETER);
  WP2_CHECK_OK(buffer.format == WP2_Argb_32, WP2_STATUS_UNSUPPORTED_FEATURE);

  const uint32_t width = buffer.width;
  const uint32_t height = buffer.height;
  const uint8_t* rgba = (const uint8_t*)buffer.GetRow(0);
  const uint32_t stride = buffer.stride;
  const uint32_t bytes_per_px = 2;
  uint32_t y;

  WP2_CHECK_OK(rgba != nullptr, WP2_STATUS_NULL_PARAMETER);

  fprintf(fout, "P5\n%u %u\n255\n", width * bytes_per_px, height);
  for (y = 0; y < height; ++y) {
    WP2_CHECK_OK(fwrite(rgba, width, bytes_per_px, fout) == bytes_per_px,
                 WP2_STATUS_BAD_WRITE);
    rgba += stride;
  }
  return WP2_STATUS_OK;
}

}  // namespace WP2
