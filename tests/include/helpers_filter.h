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
// -----------------------------------------------------------------------------

// Helper functions for filter tests.

#ifndef WP2_TESTS_INCLUDE_HELPERS_FILTER_H_
#define WP2_TESTS_INCLUDE_HELPERS_FILTER_H_

#include <vector>

#include "src/common/lossy/block.h"
#include "src/dec/filters/block_map_filter.h"
#include "src/utils/plane.h"
#include "src/utils/vector.h"

namespace WP2 {
namespace testing {

//------------------------------------------------------------------------------

int16_t GetMinValueSigned(uint32_t num_precision_bits);
int16_t GetMaxValueSigned(uint32_t num_precision_bits);

// Fills Y, U and V with extreme but single values.
YUVPlane CreatePlainPixels(uint32_t width, uint32_t height,
                           int16_t min_value, int16_t max_value);

// Adds noise.
void Noise(int32_t min_value, int32_t max_value, size_t seed, int32_t strength,
           Plane16* const plane);
void Noise(size_t seed, int32_t strength, ArgbBuffer* const canvas);

// Pixel by pixel comparison.
bool AreEqual(const Plane16& a, const Plane16& b);
bool AreEqual(const YUVPlane& a, const YUVPlane& b);

// Returns true if all channels match.
bool AreEqual(const ArgbBuffer& a, const ArgbBuffer& b);

//------------------------------------------------------------------------------

// Creates as many segments as 'qualities' with YUV bounds in [min, max_value].
Vector<Segment> CreateSegments(uint32_t quality_hint,
                               int16_t min_value, int16_t max_value);

// Helper for 'blocks->RegisterBlock()' with a CodedBlock mockup.
void RegisterBlock(uint16_t x, uint16_t y, BlockSize dim, uint32_t segment_id,
                   double res_den, double min_bpp, bool with_lossy_alpha,
                   FilterBlockMap* const blocks);

// Creates randomly sized blocks with round-robin segment selection.
// 'res_den' is the proportion of non-zero residuals and 'min_bpp' is the
// minimum number of bits per pixel (as deduced from the bitstream).
void CreateBlocks(uint32_t tile_x, uint32_t tile_y,
                  const Vector<Segment>& segments, double res_den,
                  double min_bpp, uint32_t num_precision_bits,
                  YUVPlane* const pixels, FilterBlockMap* const blocks);

//------------------------------------------------------------------------------

}  // namespace testing
}  // namespace WP2

#endif  // WP2_TESTS_INCLUDE_HELPERS_FILTER_H_
