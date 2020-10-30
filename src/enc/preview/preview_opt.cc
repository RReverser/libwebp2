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
//  Preview optimizer
//
// Author: Skal (pascal.massimino@gmail.com)

#include "./preview_enc.h"

#include <algorithm>
#include <cstdio>
#include <iterator>

#include "src/common/preview/preview.h"
#include "src/utils/random.h"
#include "src/wp2/base.h"
#include "src/wp2/encode.h"

namespace WP2 {

//------------------------------------------------------------------------------

void PreviewDataEnc::Print() const {
  fprintf(stderr,
      "%zu vertices, %zu colors, loss: %.1f, score: %.1f size: %u\n",
      vertices_.size(), palette_.size(),
      triangulation_loss_, score_, num_bytes_);
}

//------------------------------------------------------------------------------

bool PreviewData::CopyFrom(const PreviewData& other) {
  Reset();
  grid_width_  = other.grid_width_;
  grid_height_ = other.grid_height_;
  use_noise_   = other.use_noise_;
  std::copy(other.corners_, other.corners_ + 4, corners_);
  return vertices_.copy_from(other.vertices_) &&
         palette_.copy_from(other.palette_);
}

WP2Status PreviewDataEnc::ImportFrom(const PreviewData& data) {
  score_ = 0;
  triangulation_loss_ = 0;
  num_bytes_ = 0;
  WP2_CHECK_ALLOC_OK(CopyFrom(data));
  return WP2_STATUS_OK;
}

WP2Status PreviewDataEnc::ExportTo(PreviewData* const data) const {
  assert(data != nullptr);
  WP2_CHECK_ALLOC_OK(data->CopyFrom(*this));
  return WP2_STATUS_OK;
}

//------------------------------------------------------------------------------

WP2Status PreviewDataEnc::GetScore(float bigger_size_penalty) {
  // Encode() will sort vertices and palette, but it shouldn't change the result
  ANSEnc enc;
  WP2_CHECK_STATUS(Encode(&enc));
  WP2_CHECK_STATUS(enc.Assemble());
  num_bytes_ = enc.BufferSize();

  WP2_CHECK_STATUS(triangulation_.Create(*this));

  uint32_t num_covered_pixels = 0;
  WP2_CHECK_STATUS(triangulation_.GetLoss(*this, canvas_,
                                          &triangulation_loss_,
                                          &num_covered_pixels));
  // Make sure all canvas_' pixels were tested (fast good enough test).
  // TODO(yguyon): Make sure each pixel is tested exactly once.
  assert(num_covered_pixels >= canvas_.width * canvas_.height);
  // Scale loss to make it comparable to the size.
  triangulation_loss_ /= num_covered_pixels;

  score_ = triangulation_loss_ + bigger_size_penalty * num_bytes_;

  return WP2_STATUS_OK;
}

//------------------------------------------------------------------------------

namespace {

template <typename T>
constexpr T Move(T value, bool positively, int32_t distance, int32_t min_value,
                 int32_t max_value) {
  return (T)(positively ? (((int32_t)value + distance > max_value)
                               ? max_value : (int32_t)value + distance)
                        : (((int32_t)value - distance < max_value)
                               ? min_value : (int32_t)value - distance));
}

bool IsMutationHappening(UniformIntDistribution* const random,
                         uint8_t probability) {
  return (random->Get<uint8_t>(0u, 99u) < probability);
}

}  // namespace

bool PreviewDataEnc::MoveVertex(UniformIntDistribution* const rnd) {
  const uint32_t vtx_idx = rnd->Get<uint32_t>(0, vertices_.size() - 1u);
  VertexIndexedColor& vertex = vertices_[vtx_idx];
  const bool horizontally = rnd->FlipACoin();
  const bool positively = rnd->FlipACoin();
  static constexpr uint16_t distance = 1;
  uint32_t x = vertex.x, y = vertex.y;
  if (horizontally) {
    x = Move(x, positively, distance, 0, grid_width_ - 1);
  } else {
    y = Move(y, positively, distance, 0, grid_height_ - 1);
  }
  if (IsCorner(x, y)) return false;
  vertex.x = x;
  vertex.y = y;
  return true;
}

WP2Status PreviewDataEnc::AddVertex(UniformIntDistribution* const rnd) {
  VertexIndexedColor vertex = {
      rnd->Get<uint16_t>(0u, (uint16_t)grid_width_ - 1u),
      rnd->Get<uint16_t>(0u, (uint16_t)grid_height_ - 1u),
      0
  };
  if (!IsCorner(vertex.x, vertex.y)) {
    WP2_CHECK_ALLOC_OK(vertices_.push_back(vertex));
  }
  return WP2_STATUS_OK;
}

void PreviewDataEnc::RemoveVertex(UniformIntDistribution* const rnd) {
  const uint32_t vertex_index = rnd->Get<uint32_t>(0, vertices_.size() - 1u);
  vertices_[vertex_index] = vertices_.back();
  vertices_.pop_back();
}

void PreviewDataEnc::MoveColorIndex(UniformIntDistribution* const rnd) {
  const uint32_t col = rnd->Get<uint16_t>(0u, palette_.size() - 1);
  const int idx = rnd->Get<int>(-4, (int)vertices_.size() - 1);
  SetColorIndex(idx, col);
}

//------------------------------------------------------------------------------

void PreviewDataEnc::MoveColor(UniformIntDistribution* const rnd) {
  AYCoCg19b& color = palette_[rnd->Get<uint16_t>(0u, palette_.size() - 1u)];
  const uint16_t channel = rnd->Get<uint16_t>(0u, 2u);
  const bool positively = rnd->FlipACoin();
  const uint16_t distance = 1;
  const uint8_t min = 0, max = AYCoCg19b::kYCoCgMax;
  AYCoCg19b tmp = color;
  if (channel == 0) {
    tmp.y = Move(tmp.y, positively, distance, min, max);
  } else if (channel == 1) {
    tmp.co = Move(tmp.co, positively, distance, min, max);
  } else {
    tmp.cg = Move(tmp.cg, positively, distance, min, max);
  }
  // Check for dup if palette is becoming small
  if (palette_.size() == kPreviewMinNumColors) {
    for (const auto& c : palette_) {
      if (tmp.y == c.y && tmp.co == c.co &&
          tmp.cg == c.cg && tmp.a == c.a) {
        return;
      }
    }
  }
  color = tmp;  // good to go
}

WP2Status PreviewDataEnc::AddColor(UniformIntDistribution* const rnd) {
  const AYCoCg19b color = { rnd->Get<uint8_t>(0u, AYCoCg19b::kAlphaMax),
                            rnd->Get<uint8_t>(0u, AYCoCg19b::kYCoCgMax),
                            rnd->Get<uint8_t>(0u, AYCoCg19b::kYCoCgMax),
                            rnd->Get<uint8_t>(0u, AYCoCg19b::kYCoCgMax) };
  WP2_CHECK_ALLOC_OK(palette_.push_back(color));
  const uint32_t last_index = palette_.size() - 1;
  // Only reassign vertices that are closest to the new color.
  for (VertexIndexedColor& vertex : vertices_) {
    if (GetBestIndex(canvas_, vertex.x, vertex.y) == last_index) {
      vertex.color_index = last_index;
    }
  }
  // and corners too
  const uint32_t w = grid_width_ - 1;
  const uint32_t h = grid_height_ - 1;
  if (GetBestIndex(canvas_, 0, 0) == last_index) corners_[0] = last_index;
  if (GetBestIndex(canvas_, w, 0) == last_index) corners_[1] = last_index;
  if (GetBestIndex(canvas_, 0, h) == last_index) corners_[2] = last_index;
  if (GetBestIndex(canvas_, w, h) == last_index) corners_[3] = last_index;
  return WP2_STATUS_OK;
}

WP2Status PreviewDataEnc::RemoveColor(UniformIntDistribution* const rnd) {
  Vector_u16 palette_histogram;
  WP2_CHECK_STATUS(
            GetPaletteHistogram(vertices_, palette_, &palette_histogram));
  const uint32_t least_used_color_index = (uint32_t)std::abs(
      std::distance(palette_histogram.begin(),
                    std::min_element(palette_histogram.begin(),
                                     palette_histogram.end())));
  const AYCoCg19b color = palette_[least_used_color_index];
  palette_[least_used_color_index] = palette_.back();
  palette_.pop_back();
  // find the replacement color
  const uint16_t new_color_idx = FindClosestColorIndex(palette_, color);
  // Only reassign vertices that used this color.
  for (VertexIndexedColor& vertex : vertices_) {
    if (vertex.color_index == least_used_color_index) {
      vertex.color_index = new_color_idx;
    } else if (vertex.color_index == palette_.size()) {
      vertex.color_index = least_used_color_index;  // Was swapped.
    }
  }
  for (auto& c : corners_) {
    if (c == least_used_color_index) {
      c = new_color_idx;
    } else if (c == palette_.size()) {
      c = least_used_color_index;
    }
  }
  return WP2_STATUS_OK;
}

//------------------------------------------------------------------------------

WP2Status PreviewData::RemoveDuplicateVertices() {
  if (vertices_.empty()) return WP2_STATUS_OK;
  // Sort by 'y' then by 'x' for easier duplicate check.
  std::sort(vertices_.begin(), vertices_.end());

  uint32_t last = 0;
  for (uint32_t i = 1; i < vertices_.size(); ++i) {
    if (vertices_[i].x != vertices_[last].x ||
        vertices_[i].y != vertices_[last].y) {
      ++last;
      if (last < i) vertices_[last] = vertices_[i];
    }
  }
  if (++last < vertices_.size()) {
    WP2_CHECK_ALLOC_OK(vertices_.resize(last));
  }
  return WP2_STATUS_OK;
}

//------------------------------------------------------------------------------

// Assigns the color producing the smallest local loss to each vertex.
WP2Status PreviewDataEnc::OptimizeColorIndices() {
  WP2_CHECK_STATUS(triangulation_.Create(*this));
  Vector<Triangle> triangles_containing_vertex;
  WP2_CHECK_ALLOC_OK(triangles_containing_vertex.reserve(
      triangulation_.GetTriangles().size()));

  // TODO(skal): randomize order?
  for (int vtx_idx = -4; vtx_idx < (int)vertices_.size(); ++vtx_idx) {
    triangles_containing_vertex.clear();
    for (const Triangle& triangle : triangulation_.GetTriangles()) {
      if (triangle.v0 == vtx_idx || triangle.v1 == vtx_idx ||
          triangle.v2 == vtx_idx) {
        WP2_CHECK_ALLOC_OK(triangles_containing_vertex.push_back(
            triangle, /*resize_if_needed=*/false));
      }
    }

    float best_color_loss = 0.f;
    uint32_t best_color_index = 0;
    for (uint32_t color_index = 0; color_index < palette_.size();
         ++color_index) {
      float loss = 0.f;
      uint32_t num_covered_pixels = 0;
      SetColorIndex(vtx_idx, color_index);
      WP2_CHECK_STATUS(
        triangulation_.GetPartialLoss(*this, canvas_,
                                      triangles_containing_vertex,
                                      &loss, &num_covered_pixels));
      if (color_index == 0 || loss < best_color_loss) {
        best_color_loss = loss;
        best_color_index = color_index;
      }
    }
    SetColorIndex(vtx_idx, best_color_index);
  }
  return WP2_STATUS_OK;
}

//------------------------------------------------------------------------------

WP2Status PreviewData::Optimize(const ArgbBuffer& canvas,
                                const PreviewConfig& config, bool log) {
  WP2_CHECK_OK(config.IsValid(), WP2_STATUS_INVALID_CONFIGURATION);

  // Starting solution.
  PreviewDataEnc best(canvas);
  WP2_CHECK_STATUS(best.ImportFrom(*this));

  // Set optimization parameters.
  const uint32_t max_num_temp_colors =
      std::min((uint32_t)best.palette_.size() * 2u, kPreviewMaxNumColors);
  float bigger_size_penalty = config.importance_of_size_over_quality;

  const size_t random_seed = 3 * grid_width_ + 5 * grid_height_ +
      17 * best.vertices_.size() + 23 * best.palette_.size();
  UniformIntDistribution random(random_seed);

  // Get starting score.
  WP2_CHECK_STATUS(best.GetScore(bigger_size_penalty));
  if (log) {
    best.Print();
    fprintf(stderr, "Starting score: %f\n", best.Score());
  }
  float best_triangulation_score = best.Score();

  // Try some mutations, and keep the solution if it's better.
  uint32_t num_solutions_since_last_improvement = 0;
  for (uint32_t iteration = 0; iteration < config.num_iterations; ++iteration) {
    WP2_CHECK_STATUS(best.ExportTo(this));

    if (num_solutions_since_last_improvement >
        config.GetMaxNumConsecutiveSolutionsWithoutImprovement()) {
      break;
    }
    bool maybe_duplicated_vertices = false;
    for (uint32_t m = 0; m < config.num_mutations_per_iteration; ++m) {
      // Move a vertex.
      if (!best.vertices_.empty() &&
          IsMutationHappening(&random, config.proba_vertex_move)) {
        maybe_duplicated_vertices |= best.MoveVertex(&random);
      }

      // Add a vertex.
      if (IsMutationHappening(&random, config.proba_vertex_add)) {
        WP2_CHECK_STATUS(best.AddVertex(&random));
        maybe_duplicated_vertices = true;
      }

      // Remove a vertex.
      if (!best.vertices_.empty() &&
          IsMutationHappening(&random, config.proba_vertex_sub)) {
        best.RemoveVertex(&random);
      }

      // Re-assign a vertex color index.
      if (!best.vertices_.empty() &&
          IsMutationHappening(&random, config.proba_color_index_move)) {
        best.MoveColorIndex(&random);
      }

      // Change the channel Y, Co or Cg of a color.
      if (!best.palette_.empty() &&
          IsMutationHappening(&random, config.proba_color_move)) {
        best.MoveColor(&random);
      }

      // Add a color.
      if (best.palette_.size() < max_num_temp_colors &&
          IsMutationHappening(&random, config.proba_color_add)) {
        WP2_CHECK_STATUS(best.AddColor(&random));
      }

      // Remove the least used color.
      if (best.palette_.size() > kPreviewMinNumColors &&
          IsMutationHappening(&random, config.proba_color_sub)) {
        WP2_CHECK_STATUS(best.CollectColors2(canvas, best.palette_.size()));
        WP2_CHECK_STATUS(best.RemoveColor(&random));
      }
    }
    if (maybe_duplicated_vertices) {
      WP2_CHECK_STATUS(best.RemoveDuplicateVertices());
    }

    if (config.optimize_color_indices_every_n_iterations > 0 &&
        (iteration % config.optimize_color_indices_every_n_iterations) == 0) {
      WP2_CHECK_STATUS(best.OptimizeColorIndices());
    }

    // Compute current solution's score.
    WP2_CHECK_STATUS(best.GetScore(bigger_size_penalty));
    const float score = best.Score();
    const uint32_t num_bytes = best.Size();
    if (log) best.Print();

    // Pseudo simulated annealing:
    const float score_tolerance = config.score_tolerance *
                                  (config.num_iterations - iteration) /
                                  config.num_iterations;
    if (score < best_triangulation_score + score_tolerance) {
      // Take solution even if not strictly better but keep the lowest score.
      // TODO(yguyon): Return strictly best solution if nothing better at end?
      if (score < best_triangulation_score) {
        if (log) {
          fprintf(stderr, "New best solution, keeping it. "
                          "    score:%.1f -> %.1f\n",
                          best_triangulation_score, score);
        }
        best_triangulation_score = score;
        num_solutions_since_last_improvement = 0;
      } else {
        ++num_solutions_since_last_improvement;
        if (log) fprintf(stderr, "New (close enough) solution, keeping it.\n");
      }

      // Adapt 'bigger_size_penalty' to converge towards 'target_num_bytes'.
      if (config.target_num_bytes > 0 &&
          iteration > 70 * config.num_iterations / 100 &&
          bigger_size_penalty > 1.f / 2000.f) {
        const float mult = 0.9;
        if (num_bytes < config.target_num_bytes * 0.97) {
          bigger_size_penalty *= mult;
        } else if (num_bytes > config.target_num_bytes * 1.03) {
          bigger_size_penalty /= mult;
        }
      }
    } else {   // restore previous best
      WP2_CHECK_STATUS(best.ImportFrom(*this));
    }
  }
  return WP2_STATUS_OK;
}

//------------------------------------------------------------------------------

}  // namespace WP2
