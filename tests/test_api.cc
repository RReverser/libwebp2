// Copyright 2020 Google LLC
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//      http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

// API tests.

#include <iostream>

#include "src/wp2/base.h"
#include "src/wp2/encode.h"

#include "include/helpers.h"

namespace WP2 {
namespace {

void FillWholeBuffer(ArgbBuffer* const buffer) {
  buffer->Fill({0xEFu, 0x48u, 0xACu, 0x02u});
}

//------------------------------------------------------------------------------

TEST(EncodeTest, Simple) {
  ArgbBuffer input;
  ASSERT_WP2_OK(input.Resize(256, 256));
  FillWholeBuffer(&input);

  MemoryWriter output;
  ASSERT_WP2_OK(Encode(input, &output));
}

//------------------------------------------------------------------------------

TEST(EncodeTest, BadDimensions) {
  ArgbBuffer input;
  ASSERT_EQ(input.Resize(1, kMaxBufferDimension + 1), WP2_STATUS_BAD_DIMENSION);
  ASSERT_EQ(input.Resize(kMaxBufferDimension + 1, 1), WP2_STATUS_BAD_DIMENSION);

  ASSERT_WP2_OK(input.Resize(256, 256));
  input.Fill(Rectangle(0, 0, 257, 256), Argb32b({0xFFu, 0xFFu, 0xFFu, 0xFFu}));
  input.Fill(Rectangle(255, 0, 2, 1), Argb32b({0xFFu, 0xFFu, 0xFFu, 0xFFu}));

  MemoryWriter output;
  ASSERT_WP2_OK(input.Resize(1, kImageDimMax + 1));
  FillWholeBuffer(&input);
  ASSERT_EQ(Encode(input, &output), WP2_STATUS_BAD_DIMENSION);

  ASSERT_WP2_OK(input.Resize(kImageDimMax + 1, 1));
  FillWholeBuffer(&input);
  ASSERT_EQ(Encode(input, &output), WP2_STATUS_BAD_DIMENSION);
}

//------------------------------------------------------------------------------

TEST(EncodeTest, BadFormat) {
  MemoryWriter output;
  for (uint32_t format_ind = 0; format_ind < (uint32_t)WP2_FORMAT_NUM;
       ++format_ind) {
    WP2SampleFormat format = (WP2SampleFormat)format_ind;
    // MemoryWriter is in WP2_Argb_32.
    if (format == WP2_Argb_32) continue;
    ArgbBuffer input(format);
    ASSERT_WP2_OK(input.Resize(1, 1));
    for (uint32_t i = 0; i < WP2FormatBpp(format); ++i) {
      ((uint8_t*)input.GetRow(0))[i] = (uint8_t)0xFF;
    }
    ASSERT_EQ(Encode(input, &output), WP2_STATUS_INVALID_COLORSPACE);
  }
}

//------------------------------------------------------------------------------
// For coverage of helpers.cc.

TEST(EncodeTest, Include) {
  EncoderConfig config;
  std::cout << config << std::endl;
}

//------------------------------------------------------------------------------

}  // namespace
}  // namespace WP2
