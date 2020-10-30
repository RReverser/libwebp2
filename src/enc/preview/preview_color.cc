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
//  Preview color collection
//
// Author: Skal (pascal.massimino@gmail.com)

#include "./preview_enc.h"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstring>
#include <iterator>
#include <numeric>

#include "src/common/color_precision.h"
#include "src/common/preview/preview.h"

namespace WP2 {

//------------------------------------------------------------------------------

namespace {

class ColorAccumulator {
 public:
  void Add(Argb32b color) { Add(color.a, color.r, color.g, color.b); }
  void Add(uint32_t a, uint32_t r, uint32_t g, uint32_t b) {
    sum_[0] += a;
    sum_[1] += r;
    sum_[2] += g;
    sum_[3] += b;
    ++count_;
  }

  bool IsEmpty() const { return (count_ == 0); }

  Argb32b GetAverage() const {
    assert(!IsEmpty());
    const double s = count_ ? (1. / count_) : 0.;
    return {    // round down! otherwise, 255.9 turns to 256.f (= 0u)
        (uint8_t)(sum_[0] * s), (uint8_t)(sum_[1] * s),
        (uint8_t)(sum_[2] * s), (uint8_t)(sum_[3] * s)
    };
  }
  bool operator<(const ColorAccumulator& other) const {
    if (count_ != other.count_) return (count_ > other.count_);
    for (int c = 0; c < 4; ++c) {
      if (sum_[c] != other.sum_[c]) return (sum_[c] > other.sum_[c]);
    }
    return ((uintptr_t)this < (uintptr_t)&other);  // last resort
  }

 private:
  uint32_t count_ = 0;
  uint64_t sum_[4] = {0, 0, 0, 0};
};

//------------------------------------------------------------------------------

// Weighted color difference.
float GetSquareDistance(AYCoCg19b c1, AYCoCg19b c2) {
  float square_distance = 0.3f * (c1.y  - c2.y)  * (c1.y  - c2.y) +
                          0.1f * (c1.co - c2.co) * (c1.co - c2.co) +
                          0.6f * (c1.cg - c2.cg) * (c1.cg - c2.cg);
  if (c1.a != c2.a) square_distance += 0.2f * 63 * 63;
  return square_distance;
}

// Pixel to pixel distance.
constexpr float GetSquareDistance(float x1, float y1, float x2, float y2) {
  return (x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2);
}

//------------------------------------------------------------------------------

// Weighted color average.
void MergeColors(AYCoCg19b c1, uint16_t c1_weight,
                 AYCoCg19b c2, uint16_t c2_weight,
                 AYCoCg19b* const into, uint16_t* const into_weight) {
  *into_weight = c1_weight + c2_weight;
  if (*into_weight == 0) {
    into->a  = (uint8_t)((c1.a  + c2.a  + 1) >> 1);
    into->y  = (uint8_t)((c1.y  + c2.y  + 1) >> 1);
    into->co = (uint8_t)((c1.co + c2.co + 1) >> 1);
    into->cg = (uint8_t)((c1.cg + c2.cg + 1) >> 1);
  } else {
    const uint32_t round = *into_weight >> 1;
    into->a =  (uint8_t)((c1.a  * c1_weight + c2.a  * c2_weight + round) /
                         *into_weight);
    into->y =  (uint8_t)((c1.y  * c1_weight + c2.y  * c2_weight + round) /
                         *into_weight);
    into->co = (uint8_t)((c1.co * c1_weight + c2.co * c2_weight + round) /
                         *into_weight);
    into->cg = (uint8_t)((c1.cg * c1_weight + c2.cg * c2_weight + round) /
                         *into_weight);
  }
}

}  // namespace

//------------------------------------------------------------------------------

// Works with premultiplied RGB values because alpha is averaged too.
Argb32b GetAverageColor(const ArgbBuffer& canvas, uint32_t radius,
                        uint32_t pixel_x, uint32_t pixel_y) {
  const uint32_t min_x =
      (uint32_t)std::max(0, (int32_t)pixel_x - (int32_t)radius);
  const uint32_t min_y =
      (uint32_t)std::max(0, (int32_t)pixel_y - (int32_t)radius);
  const uint32_t max_x = std::min(canvas.width - 1, pixel_x + radius);
  const uint32_t max_y = std::min(canvas.height - 1, pixel_y + radius);

  assert(canvas.format == WP2_Argb_32 || canvas.format == WP2_Argb_38);
  ColorAccumulator accumulator;
  const uint32_t channel_depth = WP2FormatBpc(canvas.format);
  assert(channel_depth == 1 || channel_depth == 2);
  if (channel_depth == 1) {
    for (uint32_t y = min_y; y <= max_y; ++y) {
      const uint8_t* const row = canvas.GetRow8(y);
      for (uint32_t x = min_x; x <= max_x; ++x) {
        accumulator.Add(row[x * 4 + 0], row[x * 4 + 1], row[x * 4 + 2],
                        row[x * 4 + 3]);
      }
    }
  } else {
    for (uint32_t y = min_y; y <= max_y; ++y) {
      const uint16_t* const row = canvas.GetRow16(y);
      for (uint32_t x = min_x; x <= max_x; ++x) {
        accumulator.Add(row[x * 4 + 0], row[x * 4 + 1], row[x * 4 + 2],
                        row[x * 4 + 3]);
      }
    }
  }
  return accumulator.GetAverage();
}

static int16_t GetAverageColor(const Plane16& plane) {
  const int64_t num_pixels = plane.w_ * plane.h_;
  int64_t accumulator = 0;
  assert(kMaxYuvBits + WP2Log2Ceil_k(num_pixels) < sizeof(accumulator) * 8 - 1);

  for (uint32_t y = 0; y < plane.h_; ++y) {
    const int16_t* const row = plane.Row(y);
    accumulator = std::accumulate(row, row + plane.w_, accumulator);
  }
  return DivRound(accumulator, num_pixels);
}

WP2Status PreviewData::ReduceColors(uint32_t max_num_colors) {
  if (max_num_colors >= palette_.size()) return WP2_STATUS_OK;  // nothing to do
  Vector_u16 counts;
  WP2_CHECK_STATUS(GetPaletteHistogram(vertices_, palette_, &counts));

  Vector<AYCoCg19b> old_palette;
  WP2_CHECK_ALLOC_OK(old_palette.copy_from(palette_));

  while (palette_.size() > max_num_colors) {
    const auto least_used_color_it =
        std::min_element(counts.begin(), counts.end());
    assert(least_used_color_it != counts.end());
    const auto least_used_color_index =
        std::distance(counts.begin(), least_used_color_it);

    const AYCoCg19b least_used_color = palette_[least_used_color_index];
    const uint16_t least_used_color_count = counts[least_used_color_index];
    palette_.erase(&palette_[least_used_color_index]);
    counts.erase(&counts[least_used_color_index]);

    const uint16_t closest_color_index =
        FindClosestColorIndex(palette_, least_used_color);
    MergeColors(palette_[closest_color_index], counts[closest_color_index],
                least_used_color, least_used_color_count,
                &palette_[closest_color_index], &counts[closest_color_index]);
  }
  // note: at this point, vertices_[].color_index are mostly invalid!
  return WP2_STATUS_OK;
}

//------------------------------------------------------------------------------

uint16_t PreviewData::GetBestIndex(const ArgbBuffer& canvas,
                                   uint32_t x, uint32_t y) const {
  const uint32_t round_x = (grid_width_ - 1) >> 1;
  const uint32_t round_y = (grid_height_ - 1) >> 1;
  x = (x * (canvas.width - 1) + round_x) / (grid_width_ - 1);
  y = (y * (canvas.height - 1) + round_y) / (grid_height_ - 1);
  const Argb32b color = GetAverageColor(canvas, /*radius=*/5, x, y);
  return FindClosestColorIndex(palette_, ToAYCoCg19b(color));
}

WP2Status MakePaletteValid(Vector<AYCoCg19b>* const palette) {
  assert(palette != nullptr);
  const AYCoCg19b kDefaultPalette[kPreviewMinNumColors] = {
    ToAYCoCg19b({0xff, 0x00, 0x00, 0x00}),
    ToAYCoCg19b({0xff, 0xff, 0xff, 0xff})
  };
  if (palette->empty()) {
    for (const auto& color : kDefaultPalette) {
      WP2_CHECK_ALLOC_OK(palette->push_back(color));
    }
  } else {
    const AYCoCg19b color = palette->back();
    for (uint32_t i = palette->size(); i < kPreviewMinNumColors; ++i) {
      WP2_CHECK_ALLOC_OK(palette->push_back(color));
    }
  }
  return WP2_STATUS_OK;
}

WP2Status PreviewData::CollectColors(const ArgbBuffer& canvas,
                                     uint32_t max_colors) {
  const uint32_t num_colors = vertices_.size() + 4u;
  const uint32_t w = grid_width_ - 1;
  const uint32_t h = grid_height_ - 1;
  Vector<ColorAccumulator> centroids;
  WP2_CHECK_ALLOC_OK(centroids.resize(num_colors));

  for (uint32_t y = 0; y <= h; ++y) {
    for (uint32_t x = 0; x <= w; ++x) {
      // Initial value take care of the corner cases where grid 2x2 or
      // there's no vertices_:
      uint32_t closest_vertex_index = (x == w);
      float square_distance_to_closest_vertex = 0.f;
      // Find the Voronoi cell.
      for (uint32_t i = 0; i < vertices_.size(); ++i) {
        const float square_distance =
            GetSquareDistance(x, y, vertices_[i].x, vertices_[i].y);
        if (i == 0 || square_distance < square_distance_to_closest_vertex) {
          square_distance_to_closest_vertex = square_distance;
          closest_vertex_index = i;
        }
      }
      // Accumulate.
      const uint32_t pixel_x = x * (canvas.width - 1) / w;
      const uint32_t pixel_y = y * (canvas.height - 1) / h;
      centroids[closest_vertex_index].Add(
          GetAverageColor(canvas, /*radius=*/5, pixel_x, pixel_y));
    }
  }

  // pick for palette the centroids with most counts
  palette_.reset();
  WP2_CHECK_ALLOC_OK(palette_.reserve(num_colors));
  // note, we only need Select(max_colors), actually:
  std::sort(centroids.begin(), centroids.end());
  for (const auto& c : centroids) {
    if (!c.IsEmpty()) {
      WP2_CHECK_ALLOC_OK(palette_.push_back(ToAYCoCg19b(c.GetAverage())));
      if (palette_.size() == vertices_.size()) break;
    }
  }

  // add the corners
  for (uint32_t y : {0u, canvas.height - 1}) {
    for (uint32_t x : {0u, canvas.width - 1}) {
      const Argb32b c = GetAverageColor(canvas, /*radius=*/15, x, y);
      WP2_CHECK_ALLOC_OK(palette_.push_back(ToAYCoCg19b(c)));
    }
  }

  // Make sure we have the minimum required number of colors
  WP2_CHECK_STATUS(MakePaletteValid(&palette_));
  // ...and reduce the palette if needed.
  WP2_CHECK_STATUS(ReduceColors(max_colors));

  // find closest matching color for all vertices
  for (VertexIndexedColor& vertex : vertices_) {
    vertex.color_index = GetBestIndex(canvas, vertex.x, vertex.y);
  }
  // And corners too:
  corners_[0] = GetBestIndex(canvas, 0, 0);
  corners_[1] = GetBestIndex(canvas, w, 0);
  corners_[2] = GetBestIndex(canvas, 0, h);
  corners_[3] = GetBestIndex(canvas, w, h);

  return WP2_STATUS_OK;
}

WP2Status PreviewData::CollectColors2(const ArgbBuffer& canvas,
                                      uint32_t max_colors) {
  const uint32_t num_colors = vertices_.size() + 4u;

  Vector<ColorAccumulator> centroids;
  WP2_CHECK_ALLOC_OK(centroids.resize(num_colors));

  Triangulation triangulation;
  WP2_CHECK_STATUS(triangulation.Create(*this));
  for (const auto& triangle : triangulation.GetTriangles()) {
    Argb32b colors[3];
    WP2_CHECK_STATUS(
        triangulation.GetBestColors(*this, canvas, triangle, colors));
    centroids[triangle.v0 + 4].Add(colors[0]);
    centroids[triangle.v1 + 4].Add(colors[1]);
    centroids[triangle.v2 + 4].Add(colors[2]);
  }

  // pick for palette the centroids with most counts
  palette_.reset();
  WP2_CHECK_ALLOC_OK(palette_.reserve(num_colors));
  // TODO(skal): use partial_sort(..max_colors...) if possible
  std::sort(centroids.begin(), centroids.end());
  for (const auto& c : centroids) {
    if (!c.IsEmpty()) {
      WP2_CHECK_ALLOC_OK(palette_.push_back(ToAYCoCg19b(c.GetAverage())));
      if (palette_.size() == vertices_.size()) break;
    }
  }

  // Make sure we have the minimum required number of colors
  WP2_CHECK_STATUS(MakePaletteValid(&palette_));
  // ...and reduce the palette if needed.
  WP2_CHECK_STATUS(ReduceColors(max_colors));

  // find closest matching color for all vertices
  for (VertexIndexedColor& vertex : vertices_) {
    vertex.color_index = GetBestIndex(canvas, vertex.x, vertex.y);
  }
  // And corners too:
  const uint32_t w = grid_width_ - 1;
  const uint32_t h = grid_height_ - 1;
  corners_[0] = GetBestIndex(canvas, 0, 0);
  corners_[1] = GetBestIndex(canvas, w, 0);
  corners_[2] = GetBestIndex(canvas, 0, h);
  corners_[3] = GetBestIndex(canvas, w, h);

  return WP2_STATUS_OK;
}

//------------------------------------------------------------------------------

WP2Status GetPaletteHistogram(const Vector<VertexIndexedColor>& vertices,
                              const Vector<AYCoCg19b>& palette,
                              Vector_u16* const counts) {
  WP2_CHECK_ALLOC_OK(counts->resize(palette.size()));
  std::fill(counts->begin(), counts->end(), 0);
  for (const VertexIndexedColor& v : vertices) ++counts->at(v.color_index);
  return WP2_STATUS_OK;
}

//------------------------------------------------------------------------------

uint16_t FindClosestColorIndex(const Vector<AYCoCg19b>& palette,
                               AYCoCg19b color) {
  float closest_color_square_distance = std::numeric_limits<float>::max();
  uint16_t closest_color_index = 0;
  for (uint16_t i = 0; i < palette.size(); ++i) {
    const float square_distance = GetSquareDistance(palette[i], color);
    if (i == 0 || square_distance < closest_color_square_distance) {
      closest_color_square_distance = square_distance;
      closest_color_index = i;
    }
  }
  return closest_color_index;
}

//------------------------------------------------------------------------------

RGB12b GetPreviewColor(const ArgbBuffer& canvas) {
  return ToRGB12b(GetAverageColor(canvas, kMaxBufferDimension, 0, 0));
}

RGB12b GetPreviewColor(const YUVPlane& canvas,
                       const CSPTransform& csp_transform) {
  assert(!canvas.IsEmpty());
  Ayuv38b color;
  color.a = canvas.A.IsEmpty() ? kAlphaMax
                               : (uint8_t)Clamp<int16_t>(
                                     GetAverageColor(canvas.A), 0, kAlphaMax);
  color.y = GetAverageColor(canvas.Y);
  color.u = GetAverageColor(canvas.U);
  color.v = GetAverageColor(canvas.V);
  // TODO(yguyon): Verify it makes sense with premultiplied alpha.
  return ToRGB12b(csp_transform.ToRGB(color));
}

RGB12b GetPreviewColor(const YUVPlane& canvas, const CSPMtx& ccsp_to_rgb) {
  assert(!canvas.IsEmpty());
  int16_t ccsp[3];
  ccsp[kYChannel] = GetAverageColor(canvas.Y);
  ccsp[kUChannel] = GetAverageColor(canvas.U);
  ccsp[kVChannel] = GetAverageColor(canvas.V);
  // TODO(yguyon): Verify it makes sense to average premultiplied samples.

  uint8_t argb[4] = {};
  argb[0] = canvas.A.IsEmpty() ? kAlphaMax
                               : (uint8_t)Clamp<int16_t>(
                                     GetAverageColor(canvas.A), 0, kAlphaMax);
  Multiply(ccsp[0], ccsp[1], ccsp[2], ccsp_to_rgb.mtx(), ccsp_to_rgb.shift,
           &argb[1], &argb[2], &argb[3]);
  return ToRGB12b(Argb32b{argb[0], argb[1], argb[2], argb[3]});
}

//------------------------------------------------------------------------------

}  // namespace WP2
