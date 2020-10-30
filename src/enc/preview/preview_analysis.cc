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
//  Preview vertices collection
//
// Author: Skal (pascal.massimino@gmail.com)

#include "./preview_enc.h"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstring>

#include "src/common/preview/preview.h"
#include "src/utils/random.h"

namespace WP2 {

//------------------------------------------------------------------------------

namespace {

// Converts pixels to opaque gray levels. Transparency gives black.
WP2Status ConvertToGrayLevels(const ArgbBuffer& src, Plane8u* const dst) {
  assert(src.width < (1u << 16) && src.height < (1u << 16));
  WP2_CHECK_STATUS(dst->Resize(src.width, src.height));

  for (uint32_t y = 0; y < src.height; ++y) {
    uint8_t* const dst_y = dst->Row(y);
    const uint8_t* src_pixel = (const uint8_t*)src.GetRow(y);
    for (uint32_t x = 0; x < src.width; ++x, src_pixel += 4) {
      const uint8_t r = src_pixel[1], g = src_pixel[2], b = src_pixel[3];
      dst_y[x] = (uint8_t)std::lround(0.2126 * r + 0.7152 * g + 0.0722 * b);
    }
  }
  return WP2_STATUS_OK;
}

// Resizes and fills 'dst' with edge intensity. Uses a Sobel kernel.
WP2Status DetectEdges(const ArgbBuffer& canvas, Plane16* const dst) {
  Plane8u gray;
  WP2_CHECK_STATUS(ConvertToGrayLevels(canvas, &gray));

  dst->Fill(0);
  for (uint32_t y = 1; y < canvas.height - 1; ++y) {
    const uint8_t* const src = gray.Row(y);
    const uint32_t dst_y = y * dst->h_ / canvas.height;
    auto* const sobel = dst->Row(dst_y);
    for (uint32_t x = 1; x < canvas.width - 1; ++x) {
      const int32_t neighbors_sum = (src +          1)[x] +
                                    (src -          1)[x] +
                                    (src + gray.step_)[x] +
                                    (src - gray.step_)[x];
      const int32_t diff = std::abs(neighbors_sum - 4 * (int32_t)src[x]);
      // rescale from canvas to grid coordinate
      const uint32_t dst_x = x * dst->w_ / canvas.width;
      sobel[dst_x] = std::min(sobel[dst_x] + diff, 32767);
    }
  }
  return WP2_STATUS_OK;
}

//------------------------------------------------------------------------------

constexpr float GetPotential(float square_distance) {
  return 1.f / square_distance;  // Alternative: pow(square_distance, -2.);
}

// Allocates and fills the 'potential_table' with values inversely proportional
// to the distance to the origin.
void CreatePotentialTable(Planef* const potential_table) {
  static constexpr float kMinimumDistance = 4.f;
  static constexpr float kMinimumDistancePotential =
      GetPotential(kMinimumDistance * kMinimumDistance);
  static constexpr float kMaximumDistance = 20.f;
  static constexpr float kMaximumDistancePotential = 0.f;

  for (uint32_t y = 0; y < potential_table->h_; ++y) {
    float* const dst = potential_table->Row(y);
    for (uint32_t x = 0; x < potential_table->w_; ++x) {
      const float square_distance = (float)y * y + x * x;
      float potential;
      if (x == 0 && y == 0) {
        potential = 1e20f;  // High potential on exact spot.
      } else if (square_distance < kMinimumDistance * kMinimumDistance) {
        potential = kMinimumDistancePotential;  // Clamp max.
      } else if (square_distance > kMaximumDistance * kMaximumDistance) {
        potential = kMaximumDistancePotential;  // Clamp min.
      } else {
        potential = GetPotential(square_distance);
      }
      dst[x] = potential;
    }
  }
}

// Increases each 'repulsive' cell depending on the value of the
// 'potential_table' corresponding to the distance to 'point'.
void RepulseFutureVerticesFrom(const VertexIndexedColor& point,
                               const Planef& potential_table,
                               Planef* const repulsive) {
  assert(point.x < potential_table.w_ && point.y < potential_table.h_);
  for (uint32_t y = 0; y < repulsive->h_; ++y) {
    const uint32_t distance_y = (y > point.y) ? y - point.y : point.y - y;
    float* const dst_row = repulsive->Row(y);
    const float* const src_row = potential_table.Row(distance_y);
    for (uint32_t x = 0; x < repulsive->w_; ++x) {
      const uint32_t distance_x = (x > point.x) ? x - point.x : point.x - x;
      dst_row[x] += src_row[distance_x];
    }
  }
}

//------------------------------------------------------------------------------

float GetDiff(const ArgbBuffer& reference, const ArgbBuffer& canvas,
              uint32_t pixel_x, uint32_t pixel_y,
              uint32_t grid_w, uint32_t grid_h,
              int radius = 0) {
  assert(canvas.width == reference.width);
  assert(canvas.height == reference.height);
  const uint32_t cx = DivRound(pixel_x * (canvas.width - 1),  grid_w - 1);
  const uint32_t cy = DivRound(pixel_y * (canvas.height - 1), grid_h - 1);
  float diff = 0.f;
  const int y_limit = canvas.height - 1;
  const int x_limit = canvas.width - 1;
  for (int dy = -radius; dy <= radius; ++dy) {
    const uint32_t y = (uint32_t)Clamp((int)cy + dy, 0, y_limit);
    for (int dx = -radius; dx <= radius; ++dx) {
      const uint32_t x = (uint32_t)Clamp((int)cx + dx, 0, x_limit);
      const uint8_t* const argb_reference = reference.GetRow8(y) + 4 * x;
      const uint8_t* const argb_canvas = canvas.GetRow8(y) + 4 * x;
      for (uint32_t c = 0; c < 4; ++c) {
        diff += std::abs((float)argb_reference[c] - argb_canvas[c]);
      }
    }
  }
  const uint32_t dim = 2 * radius + 1;
  return diff / (dim * dim);
}

}  // namespace

//------------------------------------------------------------------------------

WP2Status CollectVerticesFromEdgeSelectionAndRepulsion(
    const ArgbBuffer& canvas, uint32_t max_num_vertices,
    PreviewData* const preview) {
  assert(preview != nullptr);
  Vector<VertexIndexedColor>& vertices = preview->vertices_;
  vertices.clear();

  static constexpr float kEdgeScale = 0.05f;

  Plane16 edges;
  WP2_CHECK_STATUS(edges.Resize(preview->grid_width_, preview->grid_height_));
  WP2_CHECK_STATUS(DetectEdges(canvas, &edges));

  Planef repulsive;
  WP2_CHECK_STATUS(repulsive.Resize(edges.w_, edges.h_));
  repulsive.Fill(0.f);

  Planef potential_table;
  WP2_CHECK_STATUS(potential_table.Resize(edges.w_, edges.h_));
  CreatePotentialTable(&potential_table);

  WP2_CHECK_ALLOC_OK(vertices.resize(max_num_vertices));
  uint32_t num_vertices = 0;
  for (auto& vtx : vertices) {
    float max_attraction = std::numeric_limits<float>::lowest();
    bool found = false;
    for (uint16_t y = 0; y < edges.h_; ++y) {
      const auto* const edges_row = edges.Row(y);
      const float* const repulsive_row = repulsive.Row(y);
      for (uint16_t x = 0; x < edges.w_; ++x) {
        if (preview->IsCorner(x, y)) continue;
        const float attraction = kEdgeScale * edges_row[x] - repulsive_row[x];
        if (attraction > max_attraction) {
          max_attraction = attraction;
          vtx = { x, y, 0 };
          found = true;
        }
      }
    }
    if (!found) break;
    ++num_vertices;
    RepulseFutureVerticesFrom(vtx, potential_table, &repulsive);
  }
  WP2_CHECK_ALLOC_OK(vertices.resize(num_vertices));
  return WP2_STATUS_OK;
}

//------------------------------------------------------------------------------

WP2Status CollectVerticesFromColorDiffMaximization(
    const ArgbBuffer& canvas, uint32_t max_num_vertices,
    PreviewData* const preview) {
  assert(preview != nullptr);
  static constexpr uint32_t avg_radius = 3;

  auto& vertices = preview->vertices_;
  auto& palette = preview->palette_;

  WP2_CHECK_ALLOC_OK(vertices.reserve(max_num_vertices));
  WP2_CHECK_ALLOC_OK(palette.reserve(max_num_vertices));

  ArgbBuffer result;
  WP2_CHECK_STATUS(result.Resize(canvas.width, canvas.height));

  const uint32_t w = preview->grid_width_, h = preview->grid_height_;
  Vector_u8 used;  // Vector does not handle bool very well
  WP2_CHECK_ALLOC_OK(used.resize(w * h));
  std::fill(used.begin(), used.end(), 0u);

  palette.reset();
  WP2_CHECK_STATUS(MakePaletteValid(&palette));

  // Init corners with black color.
  for (auto& c : preview->corners_) c = 0;

  Triangulation triangulation;
  WP2_CHECK_STATUS(triangulation.Create(*preview));

  // Greedily add vertices one at a time.
  for (uint32_t i = 0; i < max_num_vertices; ++i) {
    uint32_t num_covered_pixels = 0;
    WP2_CHECK_STATUS(
        triangulation.FillTriangles(*preview, &result, &num_covered_pixels));
    assert(num_covered_pixels >= result.width * result.height);

    float best_diff = -1.f;
    uint32_t best_offset = 0;
    VertexIndexedColor best_vtx = { 0, 0, (uint16_t)palette.size() };
    for (uint16_t y = 0; y < h; ++y) {
      for (uint16_t x = 0; x < w; ++x) {
        if (preview->IsCorner(x, y)) continue;
        const uint32_t offset = x + y * w;
        if (used[offset]) continue;
        const float diff = GetDiff(canvas, result, x, y, w, h);
        assert(diff >= 0.f);
        if (diff > best_diff) {
          best_diff = diff;
          best_vtx.x = x;
          best_vtx.y = y;
          best_offset = offset;
        }
      }
    }
    if (best_diff < 0) break;
    used[best_offset] = 1u;
    const uint32_t cx = best_vtx.x * (canvas.width - 1) / (w - 1);
    const uint32_t cy = best_vtx.y * (canvas.height - 1) / (h - 1);
    const AYCoCg19b col =
        ToAYCoCg19b(GetAverageColor(canvas, avg_radius, cx, cy));
    WP2_CHECK_ALLOC_OK(palette.push_back(col));
    WP2_CHECK_ALLOC_OK(vertices.push_back(best_vtx));
    WP2_CHECK_STATUS(triangulation.Insert(vertices, i));
  }
  return WP2_STATUS_OK;
}

//------------------------------------------------------------------------------

WP2Status CollectVerticesFromRandomEdgeSelection(
    const ArgbBuffer& canvas, uint32_t max_num_vertices,
    PreviewData* const preview) {
  assert(preview != nullptr);
  Vector<VertexIndexedColor>& vertices = preview->vertices_;
  vertices.clear();

  Plane16 edges;
  WP2_CHECK_STATUS(edges.Resize(preview->grid_width_, preview->grid_height_));
  WP2_CHECK_STATUS(DetectEdges(canvas, &edges));

  // Remove isolated vertices.
  for (uint16_t y = 1; y < edges.h_ - 1; ++y) {
    auto* const row = edges.Row(y);
    for (uint16_t x = 1; x < edges.w_ - 1; ++x) {
      const int nb = (int)((row[x] > row[x + edges.step_]) +
                           (row[x] > row[x + 1]) +
                           (row[x] > row[x - 1]) +
                           (row[x] > row[x - edges.step_]));
      if (nb >= 3) row[x] = 0;
    }
  }

  struct Candidate {
    uint16_t x, y, score;
    bool operator<(const Candidate& other) const {
      return (score > other.score);
    }
  };

  const size_t num_best_candidates = 2 * max_num_vertices;
  Vector<Candidate> candidates_heap;
  WP2_CHECK_ALLOC_OK(candidates_heap.reserve(num_best_candidates));

  // Discard but the best 'num_kept_vertices' candidates.
  for (uint16_t y = 0; y < edges.h_; ++y) {
    auto* const row = edges.Row(y);
    for (uint16_t x = 0; x < edges.w_; ++x) {
      if (preview->IsCorner(x, y)) continue;
      if (row[x] > 0) {
        if (candidates_heap.size() == candidates_heap.capacity()) {
          std::pop_heap(candidates_heap.begin(), candidates_heap.end());
          candidates_heap.pop_back();
        }
        WP2_CHECK_ALLOC_OK(candidates_heap.push_back(
            {x, y, (uint16_t)row[x]}, /*resize_if_needed=*/false));
        std::push_heap(candidates_heap.begin(), candidates_heap.end());
      }
    }
  }
  max_num_vertices = std::min(candidates_heap.size(), (size_t)max_num_vertices);

  // Randomly pick at most 'max_num_vertices' from the 'num_best_candidates'.
  const size_t random_seed = edges.w_ + edges.h_;
  Shuffle(&candidates_heap[0], candidates_heap.end(), random_seed);

  WP2_CHECK_ALLOC_OK(vertices.resize(max_num_vertices));
  for (uint32_t i = 0; i < max_num_vertices; ++i) {
    const Candidate candidate = candidates_heap[i];
    vertices[i] = {candidate.x, candidate.y, 0};
  }
  return WP2_STATUS_OK;
}

//------------------------------------------------------------------------------

}  // namespace WP2
