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
//  extra helper functions
//
// Author: Skal (pascal.massimino@gmail.com)

#ifndef WP2_EXTRAS_EXTRAS_H_
#define WP2_EXTRAS_EXTRAS_H_

#include <string>
#include "src/wp2/base.h"

#define WP2_EXTRAS_ABI_VERSION 0x0001    // MAJOR(8b) + MINOR(8b)

// Returns the version number of the extras library, packed in hexadecimal using
// 8bits for each of major/minor/revision. E.g: v2.5.7 is 0x020507.
WP2_EXTERN int WP2GetExtrasVersion();

namespace WP2 {

//------------------------------------------------------------------------------
// Buffer rescaling

WP2_EXTERN WP2Status ArgbBufferRescale(ArgbBuffer* buffer,
                                       uint32_t width, uint32_t height);

//------------------------------------------------------------------------------
// Matrices that can be used as is with Encode() from encode.h.
// The custom color space input samples take at most 16 bits per channel.
// The output is RGB clamped on 8 bits of theoretical precision; the matrix
// elements and the shifts aim for that range, even if more bits are kept for
// increased internal precision.

// Identity matrix, stays on 8 bits RGB (unsigned).
static constexpr int16_t kRGB8ToRGBMatrix[9] = {1, 0, 0, 0, 1, 0, 0, 0, 1};
static constexpr uint32_t kRGB8ToRGBShift = 0;  // Scaled by 2^0 = 1, alreaby 8b

// Identity matrix for 10 bits of precision (unsigned).
static constexpr int16_t kRGB10ToRGBMatrix[9] = {1, 0, 0, 0, 1, 0, 0, 0, 1};
static constexpr uint32_t kRGB10ToRGBShift = 2;  // Divide by 4, from 10b to 8b

// Lossless luma/chroma YCoCg conversion on 10 bits (sign included).
static constexpr int16_t kYCoCgToRGBMatrix[9] = {1, 1, -1, 1, 0, 1, 1, -1, -1};
static constexpr uint32_t kYCoCgToRGBShift = 2;  // From 10b to 8b

// Lossy luma/chroma SDTV BT.601 Y'UV conversion.
//  1.0  0.00000  1.13983      with 12b fixed-point precision (x 1<<12)
//  1.0 -0.39465 -0.58060      so that Y'UV 11b + matrix 12b - shift 15b = 8b
//  1.0  2.03211  0.00000      of output precision (unsigned)
static constexpr int16_t kYpUVToRGBMatrix[9] = {
    4096, 0, 4669, 4096, -1616, -2378, 4096, 8324, 0,
};
static constexpr uint32_t kYpUVToRGBShift = 15;

// YCbCr->RGB matrix from http://www.mir.com/DMG/ycbcr.html, with reduced range
// (opposite of kRGBToYCbCrMatrix[]).
static constexpr int16_t kYCbCrToRGBMatrix[9] = {4769,  0,    6537, 4769, -1605,
                                                 -3330, 4769, 8263, 0};
static constexpr uint32_t kYCbCrToRGBShift = 12;

//------------------------------------------------------------------------------
// Matrices that can be used as is with Decode() from decode.h.
// The input is RGB with 8 bits of precision.
// The custom color space output samples take at most 11 bits per channel.

// Identity matrix, stays on 8 bits RGB (unsigned).
static constexpr int16_t kRGBToRGB8Matrix[9] = {1, 0, 0, 0, 1, 0, 0, 0, 1};
static constexpr uint32_t kRGBToRGB8Shift = 0;  // Scaled by 2^0 = 1

// Identity matrix but keep 10 bits of precision (unsigned).
static constexpr int16_t kRGBToRGB10Matrix[9] = {4, 0, 0, 0, 4, 0, 0, 0, 4};
static constexpr uint32_t kRGBToRGB10Shift = 0;  // Keep x4 factor.

// Lossless luma/chroma YCoCg conversion on 10 bits (sign included).
static constexpr int16_t kRGBToYCoCgMatrix[9] = {1, 2, 1, 2, 0, -2, -1, 2, -1};
static constexpr uint32_t kRGBToYCoCgShift = 0;  // Keep x4 range.

// Lossy luma/chroma SDTV BT.601 Y'UV conversion.
//   0.299    0.587    0.114      with 12 bits fixed-point precision (x 1<<12)
//  -0.14713 -0.28886  0.436      so that RGB 8b + matrix 12b - shift 9b = 11b
//   0.615   -0.51499 -0.10001    of output precision (sign included)
static constexpr int16_t kRGBToYpUVMatrix[9] = {1225, 2404, 467,   -603, -1183,
                                                1786, 2519, -2109, -410};
static constexpr uint32_t kRGBToYpUVShift = 9;

// RGB->YCbCr matrix taken from http://www.mir.com/DMG/ycbcr.html (x 1<<12),
// [0:255]->[16:235/240] range reduction included.
static constexpr int16_t kRGBToYCbCrMatrix[] = {1052, 2065, 401,   -607, -1192,
                                                1799, 1799, -1506, -293};
constexpr uint32_t kRGBToYCbCrShift = 12;

//------------------------------------------------------------------------------
// Textual I/O of PreviewData

class PreviewData;

// Converts PreviewData to and from base64 string.
std::string PreviewToBase64(const PreviewData& data);
WP2Status PreviewFromBase64(const std::string& in, PreviewData* preview);

// Converts PreviewData to and from textural description.
std::string PreviewToText(const PreviewData& data);
WP2Status PreviewFromText(const std::string& in, PreviewData* preview);

// print some infos to a string
std::string PrintPreview(const PreviewData& data, bool reduced = true);

//------------------------------------------------------------------------------

}    // namespace WP2

#endif  /* WP2_EXTRAS_EXTRAS_H_ */
