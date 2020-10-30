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
//   Block sizes and partitions
//
// Author: Vincent Rabaud (vincent.rabaud@google.com)

#include "src/common/lossy/block_size.h"

#include <algorithm>

#include "src/utils/utils.h"
#include "src/dsp/math.h"

namespace WP2 {

//------------------------------------------------------------------------------

// clang-format off
const uint32_t BlockWidth[] = {
  1, 2, 4,
  1, 2, 4, 8,
  1, 2, 4, 8,
     2, 4, 8
};

const uint32_t BlockHeight[] = {
  1, 1, 1,
  2, 2, 2, 2,
  4, 4, 4, 4,
     8, 8, 8
};
// clang-format on

BlockSize GetBlockSize(uint32_t width, uint32_t height) {
  constexpr uint32_t kNumWidths = WP2Log2Ceil_k(kMaxBlockSize) + 1;
  constexpr uint32_t kNumHeights = kNumWidths;
  constexpr BlockSize kBlockSizes[kNumWidths][kNumHeights] = {
      {BLK_4x4, BLK_4x8, BLK_4x16, BLK_4x16},
      {BLK_8x4, BLK_8x8, BLK_8x16, BLK_8x32},
      {BLK_16x4, BLK_16x8, BLK_16x16, BLK_16x32},
      {BLK_16x4, BLK_32x8, BLK_32x16, BLK_32x32}};
  return kBlockSizes[WP2Log2Floor(std::min(width, kMaxBlockSize))]
                    [WP2Log2Floor(std::min(height, kMaxBlockSize))];
}

BlockSize GetSnappedBlockSize(uint32_t x, uint32_t y,
                              uint32_t width, uint32_t height) {
  constexpr uint32_t kMaxSizePerPos[] = {8, 1, 2, 1, 4, 1, 2, 1};
  STATIC_ASSERT_ARRAY_SIZE(kMaxSizePerPos, kMaxBlockSize);
  return GetBlockSize(std::min(width, kMaxSizePerPos[x % kMaxBlockSize]),
                      std::min(height, kMaxSizePerPos[y % kMaxBlockSize]));
}

// clang-format off
const uint32_t kNumBlocks[] = {
  1, 2,  4,
  2, 4,  8, 16,
  4, 8, 16, 32,
    16, 32, 64
};
// clang-format on
uint32_t NumPix(BlockSize dim) {
  return kNumBlocks[dim] * kMinBlockSizePix * kMinBlockSizePix;
}

// clang-format off
const BlockSize kAllBlockSizes[] = {
  BLK_4x4,  BLK_8x4,  BLK_16x4,
  BLK_4x8,  BLK_8x8,  BLK_16x8,  BLK_32x8,
  BLK_4x16, BLK_8x16, BLK_16x16, BLK_32x16,
            BLK_8x32, BLK_16x32, BLK_32x32
};

// clang-format on

STATIC_ASSERT_ARRAY_SIZE(BlockWidth, BLK_LAST);
STATIC_ASSERT_ARRAY_SIZE(BlockHeight, BLK_LAST);
STATIC_ASSERT_ARRAY_SIZE(kNumBlocks, BLK_LAST);
STATIC_ASSERT_ARRAY_SIZE(kHalfDim, BLK_LAST);
STATIC_ASSERT_ARRAY_SIZE(kFullDim, BLK_LAST);
STATIC_ASSERT_ARRAY_SIZE(kAllBlockSizes, BLK_LAST);

// clang-format off
const char* const kDimNames[BLK_LAST + 1] = {
  "BLK_4x4  ", "BLK_8x4  ", "BLK_16x4 ",
  "BLK_4x8  ", "BLK_8x8  ", "BLK_16x8 ", "BLK_32x8 ",
  "BLK_4x16 ", "BLK_8x16 ", "BLK_16x16", "BLK_32x16",
               "BLK_8x32 ", "BLK_16x32", "BLK_32x32",
  "BLK_LAST "
};
// clang-format on

STATIC_ASSERT_ARRAY_SIZE(kDimNames, BLK_LAST + 1);

//------------------------------------------------------------------------------

BlockSize GetSplitSize(BlockSize size, bool split) {
  constexpr BlockSize kSplitTf[] = {
      /*BLK_4x4=*/BLK_4x4,     /*BLK_8x4=*/BLK_4x4,
      /*BLK_16x4=*/BLK_4x4,    /*BLK_4x8=*/BLK_4x4,
      /*BLK_8x8=*/BLK_4x4,     /*BLK_16x8=*/BLK_8x8,
      /*BLK_32x8=*/BLK_8x8,    /*BLK_4x16=*/BLK_4x4,
      /*BLK_8x16=*/BLK_8x8,    /*BLK_16x16=*/BLK_8x8,
      /*BLK_32x16=*/BLK_16x16, /*BLK_8x32=*/BLK_8x8,
      /*BLK_16x32=*/BLK_16x16, /*BLK_32x32=*/BLK_16x16};
  return split ? kSplitTf[size] : size;
}

uint32_t GetNumTransformsInBlock(BlockSize size, bool split) {
  constexpr uint32_t kNumSplits[] = {
      /*BLK_4x4=*/1,   /*BLK_8x4=*/2,
      /*BLK_16x4=*/4,  /*BLK_4x8=*/2,
      /*BLK_8x8=*/4,   /*BLK_16x8=*/2,
      /*BLK_32x8=*/4,  /*BLK_4x16=*/4,
      /*BLK_8x16=*/2,  /*BLK_16x16=*/4,
      /*BLK_32x16=*/2, /*BLK_8x32=*/4,
      /*BLK_16x32=*/2, /*BLK_32x32=*/4};
  return split ? kNumSplits[size] : 1;
}

//------------------------------------------------------------------------------

// clang-format off
const uint32_t kNumCoeffs[TRF_LAST] = {
  4,   8,  16,
  8,  16,  32,  64,
  16, 32,  64, 128,  256,
      64, 128, 256,  512,
          256, 512, 1024
};
const uint32_t TrfWidth[TRF_LAST] = {
  2, 4, 8,
  2, 4, 8, 16,
  2, 4, 8, 16, 32,
     4, 8, 16, 32,
        8, 16, 32
};
const uint32_t TrfHeight[TRF_LAST] = {
  2,  2,  2,
  4,  4,  4,  4,
  8,  8,  8,  8,  8,
     16, 16, 16, 16,
         32, 32, 32
};
const TrfSize kHalfDim[BLK_LAST] = {
  TRF_2x2,  TRF_4x2,  TRF_8x2,
  TRF_2x4,  TRF_4x4,  TRF_8x4,  TRF_16x4,
  TRF_2x8,  TRF_4x8,  TRF_8x8,  TRF_16x8,
            TRF_4x16, TRF_8x16, TRF_16x16
};
const TrfSize kFullDim[BLK_LAST] = {
  TRF_4x4,  TRF_8x4,  TRF_16x4,
  TRF_4x8,  TRF_8x8,  TRF_16x8,  TRF_32x8,
  TRF_4x16, TRF_8x16, TRF_16x16, TRF_32x16,
            TRF_8x32, TRF_16x32, TRF_32x32
};
const TrfSize kAllTrfSizes[TRF_LAST] = {
  TRF_2x2, TRF_4x2,  TRF_8x2,
  TRF_2x4, TRF_4x4,  TRF_8x4,  TRF_16x4,
  TRF_2x8, TRF_4x8,  TRF_8x8,  TRF_16x8,  TRF_32x8,
           TRF_4x16, TRF_8x16, TRF_16x16, TRF_32x16,
                     TRF_8x32, TRF_16x32, TRF_32x32,
};
const uint32_t TrfLog2[kMaxBlockSizePix + 1] {
  255,  0, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3,
  4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 5
};
// clang-format on

STATIC_ASSERT_ARRAY_SIZE(kAllTrfSizes, TRF_LAST)
STATIC_ASSERT_ARRAY_SIZE(kNumCoeffs, TRF_LAST);

//------------------------------------------------------------------------------

// clang-format off
const char* const kTDimNames[TRF_LAST + 1] = {
  "TRF_2x2  ", "TRF_4x2  ", "TRF_8x2  ",
  "TRF_2x4  ", "TRF_4x4  ", "TRF_8x4  ", "TRF_16x4 ",
  "TRF_2x8  ", "TRF_4x8  ", "TRF_8x8  ", "TRF_16x8 ", "TRF_32x8 ",
               "TRF_4x16 ", "TRF_8x16 ", "TRF_16x16", "TRF_32x16",
                            "TRF_8x32 ", "TRF_16x32", "TRF_32x32",
  "TRF_LAST "
};
// clang-format on

STATIC_ASSERT_ARRAY_SIZE(kTDimNames, TRF_LAST + 1);

}  // namespace WP2
