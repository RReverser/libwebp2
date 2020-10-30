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
//  Preview triangulation
//
// Author: Skal (pascal.massimino@gmail.com)

#include "./preview.h"

namespace WP2 {

//------------------------------------------------------------------------------

bool Circle::Contains(int64_t x, int64_t y) const {
  const int64_t dx = (center_x - x * i_s);
  const int64_t dy = (center_y - y * i_s);
  const int64_t square_distance = dx * dx + dy * dy;
  // Note the '<=': this is important for replacing vertex!
  return (square_distance <= square_radius);
}

static bool TryGetCircumcircleOfTriangle(VertexIndexedColor v0,
                                         VertexIndexedColor v1,
                                         VertexIndexedColor v2,
                                         Circle* const circumcircle) {
  const int64_t dx10 = (int64_t)v1.x - v0.x;
  const int64_t dy10 = (int64_t)v1.y - v0.y;
  const int64_t dx20 = (int64_t)v2.x - v0.x;
  const int64_t dy20 = (int64_t)v2.y - v0.y;
  const int64_t inner_surface = 2 * (dx10 * dy20 - dx20 * dy10);
  if (inner_surface == 0) return false;

  circumcircle->i_s = inner_surface;
  const int64_t u = dx10 * (v0.x + v1.x) + dy10 * (v0.y + v1.y);
  const int64_t v = dx20 * (v0.x + v2.x) + dy20 * (v0.y + v2.y);
  circumcircle->center_x = +u * dy20 - v * dy10;
  circumcircle->center_y = -u * dx20 + v * dx10;
  const int64_t dx = (circumcircle->center_x - v0.x * circumcircle->i_s);
  const int64_t dy = (circumcircle->center_y - v0.y * circumcircle->i_s);
  circumcircle->square_radius = dx * dx + dy * dy;
  return true;
}

//------------------------------------------------------------------------------

WP2Status Triangulation::InitCorners(const PreviewData& data) {
  const uint16_t w = data.grid_width_ - 1;
  const uint16_t h = data.grid_height_ - 1;
  // should be consistent with GetVertex() :
  static constexpr vtx_t top_left = -1, top_right = -2;
  static constexpr vtx_t bot_left = -3, bot_right = -4;

  WP2_CHECK_ALLOC_OK(triangles_.resize(2));
  triangles_[0] = {top_left, top_right, bot_right, Circle() };
  triangles_[1] = {top_left, bot_right, bot_left, Circle() };
  corners_[0] = { 0, 0, data.corners_[0] };
  corners_[1] = { w, 0, data.corners_[1] };
  corners_[2] = { 0, h, data.corners_[2] };
  corners_[3] = { w, h, data.corners_[3] };
  WP2_CHECK_OK(
    TryGetCircumcircleOfTriangle(corners_[0], corners_[1], corners_[3],
                                 &triangles_[0].circumcircle) &&
    TryGetCircumcircleOfTriangle(corners_[0], corners_[3], corners_[2],
                                 &triangles_[1].circumcircle),
    WP2_STATUS_INVALID_PARAMETER);
  return WP2_STATUS_OK;
}

VertexIndexedColor Triangulation::GetVertex(
    const Vector<VertexIndexedColor>& vertices, vtx_t index) const {
  return (index < 0) ? corners_[-1 - index] : vertices[index];
}

//------------------------------------------------------------------------------

struct Edge { vtx_t from, to; };

static WP2Status InsertEdge(const Edge& new_edge, Vector<Edge>& edges) {
  for (uint32_t i = 0; i < edges.size(); ++i) {
    // There should not be any duplicate.
    assert(edges[i].from != new_edge.from || edges[i].to != new_edge.to);
    // However opposite edges cancel each other.
    if (edges[i].from == new_edge.to && edges[i].to == new_edge.from) {
      edges[i] = edges.back();
      edges.pop_back();
      return WP2_STATUS_OK;
    }
  }
  WP2_CHECK_ALLOC_OK(edges.push_back(new_edge));
  return WP2_STATUS_OK;
}

//------------------------------------------------------------------------------
// Bowyer–Watson algorithm
// (faster than Lawson's flipping one if you can save the circles per-triangle)

WP2Status Triangulation::Insert(const Vector<VertexIndexedColor>& vertices,
                                vtx_t new_vertex_index) {
  const VertexIndexedColor new_vertex = vertices[new_vertex_index];
  Vector<Edge> removed_edges;
  for (uint32_t t = 0; t < triangles_.size();) {
    const Triangle& triangle = triangles_[t];
    if (triangle.circumcircle.Contains(new_vertex.x, new_vertex.y)) {
      WP2_CHECK_STATUS(InsertEdge({triangle.v0, triangle.v1}, removed_edges));
      WP2_CHECK_STATUS(InsertEdge({triangle.v1, triangle.v2}, removed_edges));
      WP2_CHECK_STATUS(InsertEdge({triangle.v2, triangle.v0}, removed_edges));
      triangles_[t] = triangles_.back();
      triangles_.pop_back();
    } else {
      ++t;
    }
  }

  for (const Edge& edge : removed_edges) {
    Triangle triangle = { new_vertex_index, edge.from, edge.to, Circle() };
    if (TryGetCircumcircleOfTriangle(GetVertex(vertices, triangle.v0),
                                     GetVertex(vertices, triangle.v1),
                                     GetVertex(vertices, triangle.v2),
                                     &triangle.circumcircle)) {
      WP2_CHECK_ALLOC_OK(triangles_.push_back(triangle));
    }
  }
  return WP2_STATUS_OK;
}

//------------------------------------------------------------------------------

WP2Status Triangulation::Create(const PreviewData& data) {
  WP2_CHECK_STATUS(InitCorners(data));
  for (uint32_t i = 0; i < data.vertices_.size(); ++i) {
    WP2_CHECK_STATUS(Insert(data.vertices_, i));
  }
  return WP2_STATUS_OK;
}

//------------------------------------------------------------------------------

}  // namespace WP2
