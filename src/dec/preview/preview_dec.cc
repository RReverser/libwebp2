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
//  Preview decoder
//
// Author: Yannis Guyon (yguyon@google.com)

#include "./preview_dec.h"
#include "src/utils/data_source.h"

namespace WP2 {

//------------------------------------------------------------------------------

namespace {

int32_t ReadAValue(ANSDec* const dec, ValueStats* const stats, bool is_positive,
                   WP2_OPT_LABEL) {
  if (dec->ReadABit(&stats->zero, label)) return 0u;
  int32_t abs_value = 1;
  for (uint32_t i = 0, j = 1; j <= AYCoCg19b::kYCoCgMax; j <<= 1) {
    assert(i < AYCoCg19b::kYCoCgBitDepth);
    if (dec->ReadABit(&stats->bits[i++], label)) abs_value += j;
  }
  if (!is_positive && dec->ReadABit(&stats->sign, label)) return -abs_value;
  return abs_value;
}

}  // namespace

WP2Status PreviewData::DecodePalette(ANSDec* const dec) {
  ANSDebugPrefix prefix(dec, "palette");
  const uint32_t num_colors =
      dec->ReadRange(kPreviewMinNumColors, kPreviewMaxNumColors, "num_colors");

  const bool has_alpha = (bool)dec->ReadBit(kPreviewOpaqueProba, "has_alpha");

  WP2_CHECK_ALLOC_OK(palette_.resize(num_colors));
  AYCoCg19b pred = kPreviewStartingColorPrediction;
  ANSBinSymbol alpha(2, 2);
  ValueStats stat_yco, stat_cg;
  for (AYCoCg19b& color : palette_) {
    color.a =
        !has_alpha || (bool)(dec->ReadABit(&alpha, "residual_alpha") ^ pred.a);
    color.y =
        (uint8_t)(pred.y + ReadAValue(dec, &stat_yco, false, "residual_y"));
    color.co =
        (uint8_t)(pred.co + ReadAValue(dec, &stat_yco, false, "residual_co"));
    color.cg =
        (uint8_t)(pred.cg + ReadAValue(dec, &stat_cg, true, "residual_cg"));
    pred = color;
  }
  return dec->GetStatus();
}

//------------------------------------------------------------------------------

uint16_t PreviewData::DecodeColorIndex(ANSDec* const dec,
                                       ANSBinSymbol* const stats,
                                       uint16_t* const prediction) {
  if (dec->ReadABit(stats, "color_index_match")) {
    // TODO(yguyon): Map color index better (if only unused colors left..)
    const uint16_t coded_idx =
        dec->ReadRValue(palette_.size() - 1, "residual_color_index");
    *prediction = (coded_idx <= *prediction) ? coded_idx : coded_idx  + 1;
  }
  return *prediction;
}

WP2Status PreviewData::DecodeVertices(ANSDec* const dec) {
  uint32_t num_grid_points_left = grid_width_ * grid_height_ - 4u;
  const uint32_t min_num_vertices =
      std::max(std::max((uint32_t)palette_.size(), 4u) - 4u,
               kPreviewMinNumVertices);
  const uint32_t max_num_vertices =
      std::min(kPreviewMaxNumVertices, num_grid_points_left);

  ANSDebugPrefix prefix(dec, "vertices");
  const uint32_t num_vertices =
      dec->ReadRange(min_num_vertices, max_num_vertices, "num_vertices");
  WP2_CHECK_ALLOC_OK(vertices_.resize(num_vertices));

  // Decode corners.
  ANSBinSymbol stat_idx(2, 2);
  uint16_t predictor = 0;
  for (auto& c : corners_) c = DecodeColorIndex(dec, &stat_idx, &predictor);
  if (num_vertices == 0) return WP2_STATUS_OK;

  uint32_t vertex_index = 0;

  // Decode other points.
  for (uint16_t y = 0; y < grid_height_; ++y) {
    for (uint16_t x = 0; x < grid_width_; ++x) {
      if (IsCorner(x, y)) continue;
      const uint32_t num_vertices_left = vertices_.size() - vertex_index;
      const uint32_t proba_different_position =
          PROBA_MAX - (PROBA_MAX * num_vertices_left / num_grid_points_left);

      if (dec->ReadBit(proba_different_position, "position_match")) {
        VertexIndexedColor& vertex = vertices_[vertex_index];
        vertex.x = x;
        vertex.y = y;
        vertex.color_index = DecodeColorIndex(dec, &stat_idx, &predictor);
        if (++vertex_index == vertices_.size()) return dec->GetStatus();
      }
      --num_grid_points_left;
    }
  }
  // TODO(yguyon): Remove the possibility of BITSTREAM_ERROR
  return WP2_STATUS_BITSTREAM_ERROR;  // Vertices outside the grid.
}

//------------------------------------------------------------------------------

void PreviewData::Reset() {
  grid_width_ = 0;
  grid_height_ = 0;
  use_noise_ = false;
  vertices_.clear();
  palette_.clear();
  corners_[0] = corners_[1] = corners_[2] = corners_[3] = 0;
}

WP2Status PreviewData::Decode(ANSDec* const dec) {
  Reset();
  ANSDebugPrefix prefix(dec, "preview");
  grid_width_ = dec->ReadRange(kPreviewMinGridSize, kPreviewMaxGridSize,
                               "grid_dim");
  grid_height_ = dec->ReadRange(kPreviewMinGridSize, kPreviewMaxGridSize,
                               "grid_dim");
  use_noise_ = (bool)dec->ReadBit(kPreviewNoiseProba, "use_noise");
  WP2_CHECK_STATUS(DecodePalette(dec));
  WP2_CHECK_STATUS(DecodeVertices(dec));
  return dec->GetStatus();
}

WP2Status PreviewData::Decode(const uint8_t data[], uint32_t data_size) {
  ExternalDataSource source(data, data_size);
  ANSDec dec(&source);
  WP2_CHECK_STATUS(Decode(&dec));
  return WP2_STATUS_OK;
}

//------------------------------------------------------------------------------

WP2Status DecodePreview(const uint8_t* const data, uint32_t data_size,
                        ArgbBuffer* const output_buffer) {
  if (output_buffer != nullptr) {
    WP2_CHECK_OK(!output_buffer->IsEmpty(), WP2_STATUS_INVALID_PARAMETER);
    WP2_CHECK_OK(output_buffer->format == WP2_Argb_32,
                 WP2_STATUS_INVALID_COLORSPACE);

    PreviewData preview;
    WP2_CHECK_STATUS(preview.Decode(data, data_size));

    Triangulation triangulation;
    WP2_CHECK_STATUS(triangulation.Create(preview));

    uint32_t num_covered_pixels = 0;
    WP2_CHECK_STATUS(
      triangulation.FillTriangles(preview, output_buffer, &num_covered_pixels));
    // Make sure all 'canvas' pixels were drawn (fast good enough test).
    // TODO(yguyon): Make sure each pixel is drawn exactly once.
    assert(num_covered_pixels >= output_buffer->width * output_buffer->height);
    // TODO(skal): use_noise_
  }
  return WP2_STATUS_OK;
}

//------------------------------------------------------------------------------

}  // namespace WP2
