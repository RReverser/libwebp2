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
//  Simple GIF decoder to extract only metadata and skip everything else
//  as fast as possible.
//
// Author: Yannis Guyon (yguyon@google.com)

#include <cstring>

#include "./anim_image_dec.h"
#include "src/utils/data_source.h"
#include "src/utils/utils.h"

namespace WP2 {
namespace {

//------------------------------------------------------------------------------

// Reads 'expected_size' bytes from 'data_source' if 'expected_bytes' are
// matched. Returns false if nothing was read.
bool MatchBytes(DataSource* const data_source, const char* const expected_bytes,
                size_t expected_size) {
  const uint8_t* data = nullptr;
  if (!data_source->TryGetNext(expected_size, &data) ||
      std::memcmp((const void*)data, (const void*)expected_bytes,
                  expected_size) != 0) {
    return false;
  }
  data_source->MarkNumBytesAsRead(expected_size);
  return true;
}

bool MatchByte(DataSource* const data_source, uint8_t expected_byte) {
  return MatchBytes(data_source, (const char*)&expected_byte, 1);
}

// Reads next byte from 'data_source' to 'value'.
WP2Status ReadByte(DataSource* const data_source, uint8_t* const value) {
  const uint8_t* data = nullptr;
  WP2_CHECK_OK(data_source->TryReadNext(1, &data), WP2_STATUS_NOT_ENOUGH_DATA);
  if (value != nullptr) *value = *data;
  return WP2_STATUS_OK;
}

//------------------------------------------------------------------------------

// Discards the next 'num_bytes_to_discard' from 'data_source'.
WP2Status DiscardBytes(DataSource* const data_source,
                       size_t num_bytes_to_discard) {
  const uint8_t* data = nullptr;
  WP2_CHECK_OK(data_source->TryReadNext(num_bytes_to_discard, &data),
               WP2_STATUS_NOT_ENOUGH_DATA);
  return WP2_STATUS_OK;
}

WP2Status DiscardSubBlock(DataSource* const data_source,
                          bool* const last_sub_block) {
  uint8_t sub_block_size = 0;
  WP2_CHECK_STATUS(ReadByte(data_source, &sub_block_size));
  *last_sub_block = (sub_block_size == 0);
  WP2_CHECK_STATUS(DiscardBytes(data_source, sub_block_size));
  return WP2_STATUS_OK;
}

WP2Status DiscardSubBlocks(DataSource* const data_source) {
  bool last_sub_block = true;
  do {
    WP2_CHECK_STATUS(DiscardSubBlock(data_source, &last_sub_block));
  } while (!last_sub_block);
  return WP2_STATUS_OK;
}

WP2Status DiscardColorTable(DataSource* const data_source, uint32_t flags) {
  const bool color_table_present = !!(flags & 0x80u);
  if (color_table_present) {
    const uint32_t color_table_size = flags & 0x07u;
    WP2_CHECK_STATUS(DiscardBytes(data_source, 3 * (2u << color_table_size)));
  }
  return WP2_STATUS_OK;
}

//------------------------------------------------------------------------------

// Reads and check header.
WP2Status ReadHeader(DataSource* const data_source) {
  WP2_CHECK_OK(MatchBytes(data_source, "GIF", 3), WP2_STATUS_BITSTREAM_ERROR);
  WP2_CHECK_OK(
      MatchBytes(data_source, "87a", 3) || MatchBytes(data_source, "89a", 3),
      WP2_STATUS_BITSTREAM_ERROR);
  return WP2_STATUS_OK;
}

WP2Status DiscardLogicalScreenDescriptor(DataSource* const data_source) {
  // Skip 7-byte header except 5th byte, needed to skip color table.
  WP2_CHECK_STATUS(DiscardBytes(data_source, 4));
  uint8_t flags = 0;
  WP2_CHECK_STATUS(ReadByte(data_source, &flags));
  WP2_CHECK_STATUS(DiscardBytes(data_source, 2));
  WP2_CHECK_STATUS(DiscardColorTable(data_source, flags));
  return WP2_STATUS_OK;
}

//------------------------------------------------------------------------------

// Reads XMP or ICCP to 'metadata'.
// Usual case (including ICC profile): In each sub-block, the first byte
// specifies its size in bytes (0 to 255) and the rest of the bytes contain the
// data. Special case for XMP data: In each sub-block, the first byte is also
// part of the XMP payload. XMP in GIF also has a 257 byte padding data. See the
// XMP specification for details.
WP2Status ReadMetadata(DataSource* const data_source, bool is_xmp,
                       WP2::Data* const metadata) {
  metadata->Clear();
  while (!MatchByte(data_source, 0x00)) {
    const uint8_t* data = nullptr;
    WP2_CHECK_OK(data_source->TryGetNext(1, &data),
                 WP2_STATUS_NOT_ENOUGH_DATA);  // Not read yet.
    const uint8_t sub_block_size = *data;
    assert(sub_block_size > 0);  // Verified by MatchByte().
    WP2_CHECK_OK(data_source->TryReadNext(1u + sub_block_size, &data),
                 WP2_STATUS_NOT_ENOUGH_DATA);

    // If it is an XMP chunk, its size ('data[0]') must be included.
    WP2_CHECK_STATUS(metadata->Append(data + (is_xmp ? 0 : 1),
                                      sub_block_size + (is_xmp ? 1 : 0)));
  }
  if (is_xmp) {
    // XMP padding data is 0x01, 0xff, 0xfe ... 0x01, 0x00.
    static constexpr size_t xmp_padding_size = 257;
    if (metadata->size > xmp_padding_size) metadata->size -= xmp_padding_size;
  }
  return WP2_STATUS_OK;
}

// Reads metadata. Discards anything else.
WP2Status ReadExtension(DataSource* const data_source,
                        WP2::Metadata* const metadata) {
  if (MatchByte(data_source, 0xFF)) {
    if (MatchByte(data_source, 0x0B)) {
      if (MatchBytes(data_source, "XMP DataXMP", 11)) {
        if (metadata->xmp.IsEmpty()) {
          WP2_CHECK_STATUS(
              ReadMetadata(data_source, /*is_xmp=*/true, &metadata->xmp));
        } else {
          WP2_CHECK_STATUS(DiscardSubBlocks(data_source));
        }
      } else if (MatchBytes(data_source, "ICCRGBG1012", 11)) {
        if (metadata->iccp.IsEmpty()) {
          WP2_CHECK_STATUS(
              ReadMetadata(data_source, /*is_xmp=*/false, &metadata->iccp));
        } else {
          WP2_CHECK_STATUS(DiscardSubBlocks(data_source));
        }
      } else {
        WP2_CHECK_STATUS(DiscardBytes(data_source, 11));
        WP2_CHECK_STATUS(DiscardSubBlocks(data_source));
      }
    } else {
      WP2_CHECK_STATUS(DiscardSubBlocks(data_source));
    }
  } else {
    WP2_CHECK_STATUS(DiscardBytes(data_source, 1));
    WP2_CHECK_STATUS(DiscardSubBlocks(data_source));
  }
  return WP2_STATUS_OK;
}

//------------------------------------------------------------------------------

// Discards a frame.
WP2Status DiscardImage(DataSource* const data_source) {
  // Skip 9-byte header except 9th byte, needed to skip color table.
  WP2_CHECK_STATUS(DiscardBytes(data_source, 8));
  uint8_t flags = 0;
  WP2_CHECK_STATUS(ReadByte(data_source, &flags));
  WP2_CHECK_STATUS(DiscardColorTable(data_source, flags));
  // Skip pixels.
  WP2_CHECK_STATUS(DiscardBytes(data_source, 1));
  WP2_CHECK_STATUS(DiscardSubBlocks(data_source));
  return WP2_STATUS_OK;
}

// Reads or discards next segment from 'data_source'.
WP2Status ReadLogicalSegment(DataSource* const data_source,
                             WP2::Metadata* const metadata) {
  if (MatchByte(data_source, 0x21)) {
    WP2_CHECK_STATUS(ReadExtension(data_source, metadata));
  } else if (MatchByte(data_source, 0x2C)) {
    WP2_CHECK_STATUS(DiscardImage(data_source));
  } else {
    WP2_CHECK_STATUS(ReadByte(data_source, nullptr));  // for NOT_ENOUGH_DATA
    return WP2_STATUS_BITSTREAM_ERROR;
  }
  return WP2_STATUS_OK;
}

//------------------------------------------------------------------------------

// Internal version.
WP2Status ReadGIFMetadata(DataSource* const data_source,
                          Metadata* const metadata) {
  WP2_CHECK_OK(data_source != nullptr, WP2_STATUS_NULL_PARAMETER);
  WP2_CHECK_OK(metadata != nullptr, WP2_STATUS_NULL_PARAMETER);
  WP2_CHECK_STATUS(ReadHeader(data_source));
  WP2_CHECK_STATUS(DiscardLogicalScreenDescriptor(data_source));
  while (!MatchByte(data_source, 0x3B)) {
    WP2_CHECK_STATUS(ReadLogicalSegment(data_source, metadata));
  }
  return WP2_STATUS_OK;
}

}  // namespace

WP2Status ReadGIFMetadata(const uint8_t* const data, size_t data_size,
                          Metadata* const metadata) {
  ExternalDataSource data_source(data, data_size);
  return ReadGIFMetadata(&data_source, metadata);
}

//------------------------------------------------------------------------------

}  // namespace WP2
