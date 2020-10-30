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
//  Preview encoder
//
// Author: Skal (pascal.massimino@gmail.com)

#include "./preview_enc.h"

#include <algorithm>

#include "src/common/color_precision.h"
#include "src/common/preview/preview.h"
#include "src/dsp/math.h"
#include "src/utils/ans.h"
#include "src/utils/ans_utils.h"
#include "src/utils/utils.h"
#include "src/wp2/encode.h"

namespace WP2 {

namespace {

//------------------------------------------------------------------------------

constexpr uint32_t AbsUInt32(uint8_t a, uint8_t b) {
  return (uint32_t)((a >= b) ? (a - b) : (b - a));
}

void BlurPixel(const ArgbBuffer& src, uint32_t radius,
               uint32_t pixel_x, uint32_t pixel_y, uint32_t threshold,
               ArgbBuffer* const dst) {
  const Argb32b color = GetAverageColor(src, radius, pixel_x, pixel_y);

  assert(src.format == WP2_Argb_32 && src.format == dst->format);
  const uint8_t* const src_pixel =
      (const uint8_t*)src.GetRow(pixel_y) + pixel_x * 4;
  uint8_t* const dst_pixel = (uint8_t*)dst->GetRow(pixel_y) + pixel_x * 4;

  if (AbsUInt32(src_pixel[0], color.a) < threshold) dst_pixel[0] = color.a;
  if (AbsUInt32(src_pixel[1], color.r) < threshold) dst_pixel[1] = color.r;
  if (AbsUInt32(src_pixel[2], color.g) < threshold) dst_pixel[2] = color.g;
  if (AbsUInt32(src_pixel[3], color.b) < threshold) dst_pixel[3] = color.b;
}

WP2Status Blur(ArgbBuffer* const buffer, uint32_t radius, uint32_t threshold) {
  if (radius < 1) return WP2_STATUS_OK;

  ArgbBuffer original(buffer->format);
  WP2_CHECK_STATUS(original.CopyFrom(*buffer));

  for (uint32_t y = 0; y < buffer->height; ++y) {
    for (uint32_t x = 0; x < buffer->width; ++x) {
      BlurPixel(original, radius, x, y, threshold, buffer);
    }
  }
  return WP2_STATUS_OK;
}

//------------------------------------------------------------------------------
// 2x Horizontal averaging. Does not reallocate, stride stays the same.
bool HalfWidthDownscale(ArgbBuffer* const canvas, uint32_t min_w) {
  assert(canvas->format == WP2_Argb_32);
  assert(min_w > 1);
  const uint32_t half_width = std::max(1u, canvas->width / 2);
  if (half_width < min_w) return true;
  const uint32_t kChannels = 4;
  for (uint32_t y = 0; y < canvas->height; ++y) {
    const uint8_t* src = canvas->GetRow8(y);
    uint8_t* dst = canvas->GetRow8(y);
    for (uint32_t x = 0; x < half_width; ++x) {
      for (uint32_t c = 0; c < kChannels; ++c) {
        dst[c] = RightShiftRound(src[c] + src[c + kChannels], 1);
      }
      src += 2 * kChannels;
      dst += kChannels;
    }
  }
  canvas->width = half_width;
  return false;
}

// 2x Vertical averaging. Does not reallocate, stride stays the same.
bool HalfHeightDownscale(ArgbBuffer* const canvas, uint32_t min_h) {
  assert(canvas->format == WP2_Argb_32);
  assert(min_h > 1);
  const uint32_t half_height = std::max(1u, canvas->height / 2);
  if (half_height < min_h) return true;
  const uint32_t kChannels = 4;
  for (uint32_t y = 0; y < half_height; ++y) {
    const uint8_t* const src1 = canvas->GetRow8(2 * y + 0);
    const uint8_t* const src2 = canvas->GetRow8(2 * y + 1);
    uint8_t* const dst = canvas->GetRow8(y);
    for (uint32_t x = 0; x < kChannels * canvas->width; ++x) {
      dst[x] = RightShiftRound(src1[x] + src2[x], 1);
    }
  }
  canvas->height = half_height;
  return false;
}

//------------------------------------------------------------------------------


// Encode the 'value' with 'stats'.
void PutAValue(ANSEnc* enc, int32_t value, ValueStats* const stats,
               bool is_positive, WP2_OPT_LABEL) {
  assert(!is_positive || value >= 0);
  assert((uint32_t)std::abs(value) <= AYCoCg19b::kYCoCgMax);

  if (enc->PutABit((value == 0) ? 1 : 0, &stats->zero, label)) return;
  const bool is_negative = value < 0;
  const uint32_t abs_value = (uint32_t)std::abs(value) - 1;
  for (uint32_t i = 0, j = 1; j <= AYCoCg19b::kYCoCgMax; j <<= 1) {
    assert(i < AYCoCg19b::kYCoCgBitDepth);
    enc->PutABit((abs_value & j) ? 1 : 0, &stats->bits[i++], label);
  }
  if (!is_positive) enc->PutABit(is_negative ? 1 : 0, &stats->sign, label);
}

}  // namespace

//------------------------------------------------------------------------------
// I/O

void PreviewData::EncodeColorIndex(uint16_t idx, ANSEnc* const enc,
                                   ANSBinSymbol* const stats,
                                   uint16_t* const prediction) const {
  assert(idx < palette_.size());
  const bool color_match = (idx == *prediction);
  if (enc->PutABit(!color_match, stats, "color_index_match")) {
    const uint32_t coded_idx = (idx < *prediction) ? idx : idx - 1u;
    enc->PutRValue(coded_idx, palette_.size() - 1, "residual_color_index");
    *prediction = idx;
  }
}

WP2Status PreviewData::EncodePalette(ANSEnc* const enc) const {
  WP2_CHECK_OK(palette_.size() >= kPreviewMinNumColors &&
               palette_.size() <= kPreviewMaxNumColors,
               WP2_STATUS_INVALID_PARAMETER);
  ANSDebugPrefix prefix(enc, "palette");
  enc->PutRange(palette_.size(), kPreviewMinNumColors, kPreviewMaxNumColors,
                "num_colors");

  bool has_alpha = false;
  for (uint32_t i = 0; i < palette_.size() && !has_alpha; ++i) {
    if (!palette_[i].a) has_alpha = true;
  }
  enc->PutBit((uint32_t)has_alpha, kPreviewOpaqueProba, "has_alpha");

  AYCoCg19b pred = kPreviewStartingColorPrediction;
  ANSBinSymbol alpha(2, 2);  // TODO(yguyon): why 2, 2
  ValueStats stat_yco, stat_cg;
  for (AYCoCg19b color : palette_) {
    if (has_alpha) {
      enc->PutABit((uint32_t)(color.a ^ pred.a), &alpha, "residual_alpha");
    }
    const int32_t y = (int32_t)color.y - pred.y;
    const int32_t co = (int32_t)color.co - pred.co;
    const int32_t cg = (int32_t)color.cg - pred.cg;
    // The 'cg' value is assumed to be sorted in increasing order.
    // TODO(yguyon): why sort by 'cg'?
    PutAValue(enc, y, &stat_yco, /*is_positive=*/false, "residual_y");
    PutAValue(enc, co, &stat_yco, /*is_positive=*/false, "residual_co");
    PutAValue(enc, cg, &stat_cg, /*is_positive=*/true, "residual_cg");
    pred = color;
  }
  return WP2_STATUS_OK;
}

//------------------------------------------------------------------------------

WP2Status PreviewData::EncodeVertices(ANSEnc* const enc) const {
  uint32_t num_grid_points_left = grid_width_ * grid_height_ - 4u;

  const uint32_t min_num_vertices =
      std::max(std::max((uint32_t)palette_.size(), 4u) - 4u,
               kPreviewMinNumVertices);
  const uint32_t max_num_vertices =
      std::min(kPreviewMaxNumVertices, num_grid_points_left);
  WP2_CHECK_OK(vertices_.size() >= min_num_vertices &&
               vertices_.size() <= max_num_vertices,
               WP2_STATUS_INVALID_PARAMETER);

  ANSDebugPrefix prefix(enc, "vertices");
  enc->PutRange(vertices_.size(), min_num_vertices, max_num_vertices,
                "num_vertices");

  WP2_CHECK_OK(vertices_.size() <= num_grid_points_left,
               WP2_STATUS_INVALID_PARAMETER);

  // Encode corners (thus only color_index).
  ANSBinSymbol stat_idx(2, 2);  // TODO(yguyon): why 2, 2
  uint16_t predictor = 0;
  for (const auto c : corners_) EncodeColorIndex(c, enc, &stat_idx, &predictor);
  if (vertices_.empty()) return WP2_STATUS_OK;

  // Encode other points.
  uint32_t vertex_index = 0;
  for (uint32_t y = 0; y < grid_height_; ++y) {
    for (uint32_t x = 0; x < grid_width_; ++x) {
      if (IsCorner(x, y)) continue;
      const VertexIndexedColor& vertex = vertices_[vertex_index];
      const bool position_match = (vertex.y == y && vertex.x == x);
      const uint32_t num_vertices_left = vertices_.size() - vertex_index;
      const uint32_t proba_different_position =
          PROBA_MAX - (PROBA_MAX * num_vertices_left / num_grid_points_left);
      if (enc->PutBit((uint32_t)position_match, proba_different_position,
                      "position_match")) {
        EncodeColorIndex(vertex.color_index, enc, &stat_idx, &predictor);
        if (++vertex_index == vertices_.size()) return WP2_STATUS_OK;
      }
      --num_grid_points_left;
    }
  }
  return WP2_STATUS_INVALID_PARAMETER;  // Vertices outside the grid.
}

//------------------------------------------------------------------------------

void PreviewData::SetColorIndex(int vtx_idx, uint32_t index) {
  assert(vtx_idx >= -4 && vtx_idx < (int)vertices_.size());
  assert(index < palette_.size());
  if (vtx_idx >= 0) {
    vertices_[vtx_idx].color_index = index;
  } else {
    corners_[-1 - vtx_idx] = index;
  }
}

WP2Status PreviewData::SortPalette() {
  struct OrderedColor {
    AYCoCg19b color;
    uint32_t old_index;
    bool used_by_a_vertex;

    bool operator<(const OrderedColor& other) const {
      if (color.cg != other.color.cg) return (color.cg < other.color.cg);
      if (color.co != other.color.co) return (color.co < other.color.co);
      if (color.y != other.color.y) return (color.y < other.color.y);
      return (color.a < other.color.a);
    }
    bool operator==(const OrderedColor& other) const {
      return (color.a == other.color.a) &&
             (color.y == other.color.y) &&
             (color.co == other.color.co) &&
             (color.cg == other.color.cg);
    }
  };

  Vector<OrderedColor> ordered_palette;
  WP2_CHECK_ALLOC_OK(ordered_palette.resize(palette_.size()));
  for (uint32_t i = 0; i < ordered_palette.size(); ++i) {
    ordered_palette[i].color = palette_[i];
    ordered_palette[i].old_index = i;
    ordered_palette[i].used_by_a_vertex = false;
  }

  // Keep only used colors.
  for (const auto& vertex : vertices_) {
    WP2_CHECK_OK(vertex.color_index < ordered_palette.size(),
                 WP2_STATUS_INVALID_PARAMETER);
    ordered_palette[vertex.color_index].used_by_a_vertex = true;
  }
  for (const auto& v : corners_) ordered_palette[v].used_by_a_vertex = true;

  ordered_palette.erase(
      std::remove_if(
          ordered_palette.begin(), ordered_palette.end(),
          [](const OrderedColor& color) { return !color.used_by_a_vertex; }),
      ordered_palette.end());

  // Sort it for easier duplicate removal.
  std::sort(ordered_palette.begin(), ordered_palette.end());

  // Remember colors old indices.
  Vector_u32 old_index_to_new_index;
  WP2_CHECK_ALLOC_OK(old_index_to_new_index.resize(palette_.size()));

  palette_.clear();  // No need to reserve.
  const OrderedColor* previous_valid_color = nullptr;
  for (const OrderedColor& ordered_color : ordered_palette) {
    if (previous_valid_color != nullptr &&
        ordered_color == *previous_valid_color) {
      // Don't keep duplicates.
      old_index_to_new_index[ordered_color.old_index] =
          old_index_to_new_index[previous_valid_color->old_index];
    } else {
      old_index_to_new_index[ordered_color.old_index] = palette_.size();
      WP2_CHECK_ALLOC_OK(
          palette_.push_back(ordered_color.color, /*resize_if_needed=*/false));
      previous_valid_color = &ordered_color;
    }
  }

  // Map the vertices to the new colors.
  for (VertexIndexedColor& vertex : vertices_) {
    vertex.color_index = old_index_to_new_index[vertex.color_index];
  }
  for (auto& v : corners_) v = old_index_to_new_index[v];
  WP2_CHECK_STATUS(MakePaletteValid(&palette_));
  return WP2_STATUS_OK;
}

//------------------------------------------------------------------------------

WP2Status PreviewData::SortVertices() {
  // Sort by 'y' then by 'x' for easier duplicate check.
  std::sort(vertices_.begin(), vertices_.end());

  // Verify integrity.
  for (uint32_t i = 0; i < vertices_.size(); ++i) {
    const VertexIndexedColor& vertex = vertices_[i];
    WP2_CHECK_OK(!IsCorner(vertex.x, vertex.y), WP2_STATUS_INVALID_PARAMETER);
    // Is it inside the grid?
    WP2_CHECK_OK(!IsOutside(vertex.x, vertex.y), WP2_STATUS_INVALID_PARAMETER);
    // Is it different from the previous one?
    WP2_CHECK_OK(i == 0 ||
                 (vertex.x != vertices_[i - 1].x ||
                  vertex.y != vertices_[i - 1].y),
                 WP2_STATUS_INVALID_PARAMETER);
  }
  return WP2_STATUS_OK;
}

//------------------------------------------------------------------------------

WP2Status PreviewData::Encode(ANSEnc* const enc) const {
  WP2_CHECK_STATUS(const_cast<PreviewData*>(this)->SortVertices());
  WP2_CHECK_STATUS(const_cast<PreviewData*>(this)->SortPalette());
  {
    ANSDebugPrefix prefix(enc, "preview");
    assert(grid_width_ >= kPreviewMinGridSize);
    assert(grid_height_ >= kPreviewMinGridSize);
        enc->PutRange(grid_width_, kPreviewMinGridSize, kPreviewMaxGridSize,
                  "grid_dim");
    enc->PutRange(grid_height_, kPreviewMinGridSize, kPreviewMaxGridSize,
                  "grid_dim");
    enc->PutBit((uint32_t)use_noise_, kPreviewNoiseProba, "use_noise");
    WP2_CHECK_STATUS(EncodePalette(enc));
    WP2_CHECK_STATUS(EncodeVertices(enc));
  }
  WP2_CHECK_STATUS(enc->GetStatus());
  WP2_CHECK_STATUS(enc->Assemble());
  return WP2_STATUS_OK;
}

//------------------------------------------------------------------------------

namespace {

// Adjust buffer size and grid size.
WP2Status AdjustDimensions(const PreviewConfig& config,
                           ArgbBuffer* const buffer,
                           PreviewData* const preview) {
  uint32_t buffer_w = buffer->width, buffer_h = buffer->height;
  uint32_t grid_w = config.grid_size, grid_h = config.grid_size;

  if (buffer_w >= buffer_h) {
    grid_h = DivRound(grid_w * buffer_h, buffer_w);
  } else {
    grid_w = DivRound(grid_h * buffer_w, buffer_h);
  }
  grid_w = std::min(grid_w, kPreviewMaxGridSize);
  grid_h = std::min(grid_h, kPreviewMaxGridSize);
  grid_w = std::max(kPreviewMinGridSize, std::min(grid_w, buffer_w / 2));
  grid_h = std::max(kPreviewMinGridSize, std::min(grid_h, buffer_h / 2));

  preview->grid_width_ = grid_w;
  preview->grid_height_ = grid_h;

  // printf("Canvas: %d,%d,  grid:%d,%d ", buffer_w, buffer_h, grid_w, grid_h);

  buffer_w = std::lround(grid_w * config.grid_density);
  buffer_h = std::lround(grid_h * config.grid_density);
  // Downscale by 2x in alternate direction until reaching the dimension
  // just above the final one
  // TODO(skal): we could do a last pass of refinement with a simple
  // bilinear interpolation.
  bool done = false;
  while (!done) {
    done = HalfWidthDownscale(buffer, buffer_w) ||
           HalfHeightDownscale(buffer, buffer_h);
  }

  // printf("final Canvas: %d,%d\n", buffer->width, buffer->height);
  return WP2_STATUS_OK;
}

}  // namespace

WP2Status PreviewData::Generate(const ArgbBuffer& canvas,
                                const PreviewConfig& config) {
  WP2_CHECK_OK(grid_width_ >= kPreviewMinGridSize &&
               grid_width_ <= kPreviewMaxGridSize,
               WP2_STATUS_INVALID_PARAMETER);
  WP2_CHECK_OK(grid_height_ >= kPreviewMinGridSize &&
               grid_height_ <= kPreviewMaxGridSize,
               WP2_STATUS_INVALID_PARAMETER);

  const uint32_t max_num_vertices = std::min(config.num_vertices,
                                             grid_width_ * grid_height_ - 4);
  // Find best vertices.
  if (config.analysis_method ==
      PreviewConfig::AnalysisMethod::kEdgeSelectionAndRepulsion) {
    WP2_CHECK_STATUS(CollectVerticesFromEdgeSelectionAndRepulsion(
        canvas, max_num_vertices, this));
  } else if (config.analysis_method ==
             PreviewConfig::AnalysisMethod::kColorDiffMaximization) {
    WP2_CHECK_STATUS(CollectVerticesFromColorDiffMaximization(
        canvas, max_num_vertices, this));
  } else {
    WP2_CHECK_STATUS(CollectVerticesFromRandomEdgeSelection(
        canvas, max_num_vertices, this));
  }
  WP2_CHECK_STATUS(CollectColors2(canvas, config.num_colors));
  WP2_CHECK_STATUS(Optimize(canvas, config, /* log= */false));
  return WP2_STATUS_OK;
}

//------------------------------------------------------------------------------
// Entry point

WP2Status EncodePreview(const ArgbBuffer& buffer, const PreviewConfig& config,
                        Writer* const output) {
  WP2_CHECK_OK(buffer.format == WP2_Argb_32, WP2_STATUS_INVALID_COLORSPACE);
  WP2_CHECK_OK(config.IsValid(), WP2_STATUS_INVALID_CONFIGURATION);
  WP2_CHECK_OK(output != nullptr, WP2_STATUS_NULL_PARAMETER);

  // Prepare canvas.
  ArgbBuffer canvas(buffer.format);
  WP2_CHECK_STATUS(canvas.CopyFrom(buffer));

  PreviewData preview;
  WP2_CHECK_STATUS(AdjustDimensions(config, &canvas, &preview));

  WP2_CHECK_STATUS(Blur(&canvas, config.blur_radius, config.blur_threshold));

  // Generate and optimize.
  WP2_CHECK_STATUS(preview.Generate(canvas, config));

  // Create bitstream.
  ANSEnc enc;
  WP2_CHECK_STATUS(preview.Encode(&enc));
  WP2_CHECK_OK(output->Append(enc.Buffer(), enc.BufferSize()),
               WP2_STATUS_BAD_WRITE);
  return WP2_STATUS_OK;
}

//------------------------------------------------------------------------------

}  // namespace WP2
