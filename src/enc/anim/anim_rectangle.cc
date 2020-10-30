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
//  RectangleGroup implementation.
//
// Author: Yannis Guyon (yguyon@google.com)

#include <algorithm>
#include <cassert>
#include <utility>

#include "./anim_enc.h"

namespace WP2 {

//------------------------------------------------------------------------------

namespace {

Rectangle BoundingRectangle(const Rectangle& a, const Rectangle& b) {
  const uint32_t left_x = std::min(a.x, b.x);
  const uint32_t top_y = std::min(a.y, b.y);
  const uint32_t right_x = std::max(a.x + a.width, b.x + b.width);
  const uint32_t bottom_y = std::max(a.y + a.height, b.y + b.height);
  return {left_x, top_y, right_x - left_x, bottom_y - top_y};
}

uint64_t GetMergeAreaIncrease(const Rectangle& a, const Rectangle& b) {
  const uint64_t area_removed = a.GetArea() + b.GetArea();
  const uint64_t area_added = BoundingRectangle(a, b).GetArea();
  // Overlap is possible.
  return (area_added > area_removed) ? (area_added - area_removed) : 0u;
}

}  // namespace

//------------------------------------------------------------------------------

// 'max_num_rectangles_' is clamped to make sure there is no oveflow but still
// assert that it's the expected value.
RectangleGroup::RectangleGroup(uint32_t max_num_rectangles)
    : max_num_rectangles_(
          std::min(max_num_rectangles, kMaxStoredRectangles - 1)) {
  assert(max_num_rectangles < kMaxStoredRectangles);
  assert(max_num_rectangles_ >= 1);
}

void RectangleGroup::Clear() { num_rectangles_ = 0; }

//------------------------------------------------------------------------------

void RectangleGroup::MergeRectangles(uint32_t a_index, uint32_t b_index) {
  if (a_index != b_index) {
    assert(a_index < num_rectangles_ && b_index < num_rectangles_);
    if (a_index > b_index) std::swap(a_index, b_index);

    rectangles_[a_index] =
        BoundingRectangle(rectangles_[a_index], rectangles_[b_index]);

    --num_rectangles_;
    if (b_index != num_rectangles_) {
      rectangles_[b_index] = rectangles_[num_rectangles_];
    }
  }
}

bool RectangleGroup::TryToMergeARectangleWith(
    uint32_t rectangle_index, uint32_t* const other_rectangle_index) {
  // Browse from the end, it's likely that the last two rectangles are adjacent.
  uint32_t i = num_rectangles_;
  while (i-- > 0) {
    if (i != rectangle_index) {
      const uint64_t merge_area_increase =
          GetMergeAreaIncrease(rectangles_[i], rectangles_[rectangle_index]);
      if (merge_area_increase == 0) {
        MergeRectangles(i, rectangle_index);
        *other_rectangle_index = std::min(i, rectangle_index);
        return true;
      }
    }
  }
  return false;
}

//------------------------------------------------------------------------------

void RectangleGroup::AddPoint(uint32_t x, uint32_t y) {
  assert(x <= kMaxBufferDimension && y <= kMaxBufferDimension);

  // If it's contained in an existing rectangle, nothing to do.
  for (uint32_t i = 0; i < num_rectangles_; ++i) {
    if (rectangles_[i].Contains(x, y)) return;
  }

  // Otherwise try to expand an existing rectangle by the new point without
  // increasing the sum of the covered areas. Ripple changes if found.
  rectangles_[num_rectangles_] = {x, y, 1, 1};
  ++num_rectangles_;

  uint32_t rectangle_index = num_rectangles_ - 1;
  uint32_t other_rectangle_index;
  while (TryToMergeARectangleWith(rectangle_index, &other_rectangle_index)) {
    rectangle_index = other_rectangle_index;
  }

  if (num_rectangles_ > max_num_rectangles_) {
    // As last resort merge the pair minimizing the increased covered area.
    std::pair<uint32_t, uint32_t> best_merge = {0, 0};
    uint64_t best_merge_area_increase = 0;
    for (uint32_t i = 0; i < num_rectangles_; ++i) {
      for (uint32_t j = i + 1; j < num_rectangles_; ++j) {
        const uint64_t merge_area_increase =
            GetMergeAreaIncrease(rectangles_[i], rectangles_[j]);
        if (best_merge.first == best_merge.second ||
            merge_area_increase < best_merge_area_increase) {
          best_merge = {i, j};
          best_merge_area_increase = merge_area_increase;
        }
      }
    }
    MergeRectangles(best_merge.first, best_merge.second);

    // Ripple changes.
    rectangle_index = best_merge.first;
    while (TryToMergeARectangleWith(rectangle_index, &other_rectangle_index)) {
      rectangle_index = other_rectangle_index;
    }
  }
}

//------------------------------------------------------------------------------

void RectangleGroup::SetRectangle(const Rectangle& rectangle) {
  assert(rectangle.x <= kMaxBufferDimension &&
         rectangle.y <= kMaxBufferDimension &&
         rectangle.width <= kMaxBufferDimension &&
         rectangle.height <= kMaxBufferDimension);

  rectangles_[0] = rectangle;
  num_rectangles_ = 1;
}

//------------------------------------------------------------------------------

uint32_t RectangleGroup::GetNumRectangles() const { return num_rectangles_; }

const Rectangle& RectangleGroup::GetRectangle(uint32_t index) const {
  assert(index < num_rectangles_);
  return rectangles_[index];
}

Rectangle RectangleGroup::GetBoundingRectangle() const {
  if (num_rectangles_ == 0) return Rectangle();
  Rectangle bounding_rectangle = rectangles_[0];
  for (uint32_t i = 1; i < num_rectangles_; ++i) {
    bounding_rectangle = BoundingRectangle(rectangles_[i], bounding_rectangle);
  }
  return bounding_rectangle;
}

uint64_t RectangleGroup::GetCumulativeArea() const {
  uint64_t cumulative_area = 0;
  for (uint32_t i = 0; i < num_rectangles_; ++i) {
    cumulative_area += rectangles_[i].GetArea();
  }
  return cumulative_area;
}

//------------------------------------------------------------------------------

}  // namespace WP2
