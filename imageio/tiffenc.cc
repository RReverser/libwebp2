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

#include "./image_enc.h"
#include "./imageio_util.h"
#include "src/dsp/dsp.h"
#include "src/utils/utils.h"
#include "src/utils/vector.h"

namespace WP2 {

namespace {

void PutLE16(uint8_t* const dst, uint32_t value) {
  dst[0] = (value >> 0) & 0xff;
  dst[1] = (value >> 8) & 0xff;
}

void PutLE32(uint8_t* const dst, uint32_t value) {
  PutLE16(dst + 0, (value >> 0) & 0xffff);
  PutLE16(dst + 2, (value >> 16) & 0xffff);
}

}  // namespace

#define NUM_IFD_ENTRIES 15
#define EXTRA_DATA_SIZE 16
// 10b for signature/header + n * 12b entries + 4b for IFD terminator:
#define EXTRA_DATA_OFFSET (10 + 12 * NUM_IFD_ENTRIES + 4)
#define TIFF_HEADER_SIZE (EXTRA_DATA_OFFSET + EXTRA_DATA_SIZE)

WP2Status WriteTIFF(const ArgbBuffer& buffer, FILE* fout) {
  WP2_CHECK_OK(fout != nullptr, WP2_STATUS_NULL_PARAMETER);
  WP2_CHECK_OK(!buffer.IsEmpty(), WP2_STATUS_INVALID_PARAMETER);
  WP2_CHECK_OK(buffer.format == WP2_Argb_32, WP2_STATUS_UNSUPPORTED_FEATURE);

  // TODO(skal): should be has_alpha = WP2ArgbBufferHasTransparency(buffer);
  // but tests are failing with that.
  const int has_alpha = 1;
  const uint32_t width = buffer.width;
  const uint32_t height = buffer.height;
  const uint8_t* Argb = (const uint8_t*)buffer.GetRow(0);
  const uint32_t stride = buffer.stride;
  const WP2SampleFormat fmt = has_alpha ? WP2_rgbA_32 : WP2_RGB_24;
  const uint8_t bytes_per_px = has_alpha ? 4 : 3;
  // For non-alpha case, we omit tag 0x152 (ExtraSamples).
  const uint8_t num_ifd_entries = has_alpha ? NUM_IFD_ENTRIES
                                            : NUM_IFD_ENTRIES - 1;
  uint8_t tiff_header[TIFF_HEADER_SIZE] = {
    0x49, 0x49, 0x2a, 0x00,   // little endian signature
    8, 0, 0, 0,               // offset to the unique IFD that follows
    // IFD (offset = 8). Entries must be written in increasing tag order.
    num_ifd_entries, 0,       // Number of entries in the IFD (12 bytes each).
    0x00, 0x01, 4, 0, 1, 0, 0, 0, 0, 0, 0, 0,    //  10: Width  (TBD)
    0x01, 0x01, 4, 0, 1, 0, 0, 0, 0, 0, 0, 0,    //  22: Height (TBD)
    0x02, 0x01, 3, 0, bytes_per_px, 0, 0, 0,     //  34: BitsPerSample: 8888
        EXTRA_DATA_OFFSET + 0, 0, 0, 0,
    0x03, 0x01, 3, 0, 1, 0, 0, 0, 1, 0, 0, 0,    //  46: Compression: none
    0x06, 0x01, 3, 0, 1, 0, 0, 0, 2, 0, 0, 0,    //  58: Photometric: RGB
    0x11, 0x01, 4, 0, 1, 0, 0, 0,                //  70: Strips offset:
        TIFF_HEADER_SIZE, 0, 0, 0,               //      data follows header
    0x12, 0x01, 3, 0, 1, 0, 0, 0, 1, 0, 0, 0,    //  82: Orientation: topleft
    0x15, 0x01, 3, 0, 1, 0, 0, 0,                //  94: SamplesPerPixels
        bytes_per_px, 0, 0, 0,
    0x16, 0x01, 4, 0, 1, 0, 0, 0, 0, 0, 0, 0,    // 106: Rows per strip (TBD)
    0x17, 0x01, 4, 0, 1, 0, 0, 0, 0, 0, 0, 0,    // 118: StripByteCount (TBD)
    0x1a, 0x01, 5, 0, 1, 0, 0, 0,                // 130: X-resolution
        EXTRA_DATA_OFFSET + 8, 0, 0, 0,
    0x1b, 0x01, 5, 0, 1, 0, 0, 0,                // 142: Y-resolution
        EXTRA_DATA_OFFSET + 8, 0, 0, 0,
    0x1c, 0x01, 3, 0, 1, 0, 0, 0, 1, 0, 0, 0,    // 154: PlanarConfiguration
    0x28, 0x01, 3, 0, 1, 0, 0, 0, 2, 0, 0, 0,    // 166: ResolutionUnit (inch)
    0x52, 0x01, 3, 0, 1, 0, 0, 0, 1, 0, 0, 0,    // 178: ExtraSamples: rgbA
    0, 0, 0, 0,                                  // 190: IFD terminator
    // EXTRA_DATA_OFFSET:
    8, 0, 8, 0, 8, 0, 8, 0,      // BitsPerSample
    72, 0, 0, 0, 1, 0, 0, 0      // 72 pixels/inch, for X/Y-resolution
  };

  size_t bytes_per_row = 0;
  WP2_CHECK_STATUS(MultFitsIn(width, bytes_per_px, &bytes_per_row));
  uint32_t strip_byte_count = 0;
  WP2_CHECK_STATUS(MultFitsIn(bytes_per_row, height, &strip_byte_count));

  // Fill placeholders in IFD:
  PutLE32(tiff_header + 10 + 8, width);
  PutLE32(tiff_header + 22 + 8, height);
  PutLE32(tiff_header + 106 + 8, height);
  PutLE32(tiff_header + 118 + 8, strip_byte_count);
  if (!has_alpha) PutLE32(tiff_header + 178, 0);  // IFD terminator

  Vector_u8 dest;
  WP2_CHECK_ALLOC_OK(dest.resize(bytes_per_row));
  WP2ArgbConverterInit();

  // write header
  WP2_CHECK_OK(fwrite(tiff_header, sizeof(tiff_header), 1, fout) == 1,
               WP2_STATUS_BAD_WRITE);
  // write pixel values
  for (uint32_t y = 0; y < height; ++y) {
    WP2ArgbConvertTo[fmt](Argb, width, dest.data());
    WP2_CHECK_OK(fwrite(dest.data(), bytes_per_px, width, fout) == width,
                 WP2_STATUS_BAD_WRITE);
    Argb += stride;
  }
  return WP2_STATUS_OK;
}

#undef TIFF_HEADER_SIZE
#undef EXTRA_DATA_OFFSET
#undef EXTRA_DATA_SIZE
#undef NUM_IFD_ENTRIES

}  // namespace WP2
