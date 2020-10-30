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
//  Preview rasterizer
//
// Author: Skal (pascal.massimino@gmail.com)

#include "./preview.h"

#include <cmath>

#include "src/dsp/math.h"
#include "src/utils/utils.h"
#include "src/utils/wiener.h"

namespace WP2 {

// Used for drawing.
struct VertexArgb {
  uint16_t x, y;
  Argb32b color;
};

// Maps grid coordinates to canvas pixels and draws triangles.
class Rasterizer {
 public:
  Rasterizer() : vertical_span_{0, 0} {}

  WP2Status Init(const ArgbBuffer& canvas, const PreviewData& preview) {
    img_w_ = canvas.width;
    img_h_ = canvas.height;
    if (img_w_ > (1u << 31) / img_h_) return WP2_STATUS_INVALID_PARAMETER;
    mult_w_ = DivCeil(img_w_ << kBits, preview.grid_width_ - 1u);
    mult_h_ = DivCeil(img_h_ << kBits, preview.grid_height_ - 1u);
    WP2_CHECK_ALLOC_OK(row_spans_.resize(img_h_));
    return WP2_STATUS_OK;
  }

  // Fills the triangle 'v0,v1,v2' in 'canvas' with interpolated colors.
  void FillTriangle(VertexArgb v0, VertexArgb v1, VertexArgb v2,
                    ArgbBuffer* const canvas, uint32_t* const num_pixels) {
    RecordTriangleSpan(&v0, &v1, &v2);
    DrawScans(v0, v1, v2, canvas, num_pixels);
  }

  // Same as FillTriangle() but returns the loss from 'reference' instead.
  float GetLossInTriangle(const ArgbBuffer& reference,
                          VertexArgb v0, VertexArgb v1, VertexArgb v2,
                          uint32_t* const num_pixels) {
    RecordTriangleSpan(&v0, &v1, &v2);
    return GetLossInScans(reference, v0, v1, v2, num_pixels);
  }

  void CollectColors(const ArgbBuffer& reference,
                     VertexArgb& v0, VertexArgb& v1, VertexArgb& v2);

 protected:
  // Increases vertical span and sets horizontal spans covered by triangle.
  void RecordTriangleSpan(VertexArgb* const v0,
                          VertexArgb* const v1,
                          VertexArgb* const v2);

  // Increases vertical span and sets horizontal spans covered by segment 'v0v1'
  void RecordEdgeSpan(VertexArgb v0, VertexArgb v1);

  // Paints all pixels in the area defined by the spans.
  void DrawScans(VertexArgb v0, VertexArgb v1, VertexArgb v2,
                 ArgbBuffer* const canvas,
                 uint32_t* const num_pixels) const;

  // Same as DrawScans() but returns the loss from 'canvas' instead.
  float GetLossInScans(const ArgbBuffer& reference,
                       VertexArgb v0, VertexArgb v1, VertexArgb v2,
                       uint32_t* const num_pixels) const;

  // Pixels are drawn in range [from, to)
  struct Span { uint16_t from, to; };

  // fixed-point precision when converting from grid-space to canvas-space
  const uint32_t kBits = 17u;
  const uint32_t kRound = (1u << kBits >> 1);
  uint32_t mult_w_, mult_h_;
  uint32_t img_w_, img_h_;    // width / height limit

  Span vertical_span_;
  Vector<Span> row_spans_;
};

// TODO(yguyon): Remove double from these functions, it's not portable enough
//------------------------------------------------------------------------------

void Rasterizer::RecordEdgeSpan(VertexArgb v0, VertexArgb v1) {
  uint32_t from_y = v0.y;
  uint32_t to_y = v1.y;
  if (from_y == to_y) return;  // Lines are drawn in range [from_y, to_y)

  const double slope = ((double)v1.x - v0.x) / ((double)v1.y - v0.y);
  const bool inverted = (from_y > to_y);
  double x = inverted ? v1.x : v0.x;
  if (inverted) std::swap(from_y, to_y);
  from_y = Clamp(from_y, 0u, img_h_ - 1);
  to_y = Clamp(to_y, 0u, img_h_);
  if (from_y < vertical_span_.from) vertical_span_.from = from_y;
  if (to_y > vertical_span_.to) vertical_span_.to = to_y;

  if (inverted) {
    for (uint32_t y = from_y; y < to_y; ++y) {
      row_spans_[y].from = Clamp((uint32_t)x, 0u, img_w_ - 1);
      x += slope;
    }
  } else {
    for (uint32_t y = from_y; y < to_y; ++y) {
      row_spans_[y].to = Clamp((uint32_t)x, 0u, img_w_);
      x += slope;
    }
  }
}

//------------------------------------------------------------------------------

static void CalcArea(VertexArgb v0, VertexArgb v1, VertexArgb v2,
                     double iM[4]) {
  const double dx1 = v1.x - v0.x, dy1 = v1.y - v0.y;
  const double dx2 = v2.x - v0.x, dy2 = v2.y - v0.y;
  const double area = dx1 * dy2 - dx2 * dy1;
  const double iArea = (area == 0.) ? 0. : 1. / area;
  iM[0] =  dy2 * iArea;
  iM[1] = -dy1 * iArea;
  iM[2] = -dx2 * iArea;
  iM[3] =  dx1 * iArea;
}

static void CalcGradient(const double iM[4],
                         float v0x, float v0y,
                         float f0, float f1, float f2,
                         double grad[3]) {
  const double df1 = f1 - f0;
  const double df2 = f2 - f0;
  grad[0] = iM[0] * df1 + iM[1] * df2;
  grad[1] = iM[2] * df1 + iM[3] * df2;
  grad[2] = f0 - grad[0] * v0x - grad[1] * v0y;
}

static uint8_t RoundToUInt8(double v) {
  return (uint8_t)((v < 0.) ? 0u : (v > 255.) ? 255u : std::lround(v));
}

static void GetGradients(VertexArgb v0, VertexArgb v1, VertexArgb v2,
                         double gradA[3],
                         double gradR[3], double gradG[3], double gradB[3]) {
  double iM[4];
  CalcArea(v0, v1, v2, iM);
  CalcGradient(iM, v0.x, v0.y, v0.color.a, v1.color.a, v2.color.a, gradA);
  CalcGradient(iM, v0.x, v0.y, v0.color.r, v1.color.r, v2.color.r, gradR);
  CalcGradient(iM, v0.x, v0.y, v0.color.g, v1.color.g, v2.color.g, gradG);
  CalcGradient(iM, v0.x, v0.y, v0.color.b, v1.color.b, v2.color.b, gradB);
}

void Rasterizer::DrawScans(VertexArgb v0, VertexArgb v1, VertexArgb v2,
                           ArgbBuffer* const canvas,
                           uint32_t* const num_pixels) const {
  assert(canvas->width == img_w_ && canvas->height == img_h_);
  double gradA[3], gradR[3], gradG[3], gradB[3];
  GetGradients(v0, v1, v2, gradA, gradR, gradG, gradB);
  double base_A = gradA[1] * vertical_span_.from + gradA[2];
  double base_R = gradR[1] * vertical_span_.from + gradR[2];
  double base_G = gradG[1] * vertical_span_.from + gradG[2];
  double base_B = gradB[1] * vertical_span_.from + gradB[2];

  for (uint32_t y = vertical_span_.from; y < vertical_span_.to; ++y) {
    uint8_t* const argb = canvas->GetRow8(y);
    const Span& horizontal_span = row_spans_[y];

    double cur_A = base_A + gradA[0] * horizontal_span.from;
    double cur_R = base_R + gradR[0] * horizontal_span.from;
    double cur_G = base_G + gradG[0] * horizontal_span.from;
    double cur_B = base_B + gradB[0] * horizontal_span.from;
    base_A += gradA[1];
    base_R += gradR[1];
    base_G += gradG[1];
    base_B += gradB[1];

    for (uint32_t x = horizontal_span.from; x < horizontal_span.to; ++x) {
      uint8_t* const argb_pixel = argb + x * 4;
      argb_pixel[0] = RoundToUInt8(cur_A);
      argb_pixel[1] = RoundToUInt8(cur_R);
      argb_pixel[2] = RoundToUInt8(cur_G);
      argb_pixel[3] = RoundToUInt8(cur_B);
      cur_A += gradA[0];
      cur_R += gradR[0];
      cur_G += gradG[0];
      cur_B += gradB[0];
      ++*num_pixels;
    }
#if 0   // uncomment to show the triangles:
    argb[4 * horizontal_span.from + 0] = 0xff;
    argb[4 * horizontal_span.from + 1] = 0xff;
    argb[4 * horizontal_span.from + 2] = 0xff;
    argb[4 * horizontal_span.from + 3] = 0xff;
#endif
  }
}

static uint8_t ColorFromGradient(const float g[3], int x, int y) {
  const int v = (int)(g[0] * x + g[1] * y + g[2]);
  return Clamp<uint8_t>(v, 0, 255);
}

static void GetColor(const float grad[4][3], int x, int y, Argb32b* const c) {
  c->a = ColorFromGradient(grad[0], x, y);
  c->r = ColorFromGradient(grad[1], x, y);
  c->g = ColorFromGradient(grad[2], x, y);
  c->b = ColorFromGradient(grad[3], x, y);
}

static void GetRawColor(const ArgbBuffer& buf, int x, int y, Argb32b* const c) {
  x = std::max(0, std::min(x, (int)buf.width - 1));
  y = std::max(0, std::min(y, (int)buf.height - 1));
  *c = ToArgb32b(buf.GetRow8(y) + 4 * x);
}

void Rasterizer::CollectColors(const ArgbBuffer& reference,
                               VertexArgb& v0, VertexArgb& v1, VertexArgb& v2) {
  RecordTriangleSpan(&v0, &v1, &v2);  // modifies x/y !

  WienerOptimizer<2, 4> opt;  // linearly predict four A,R,G,B values from x,y
  for (uint16_t y = vertical_span_.from; y < vertical_span_.to; ++y) {
    const uint8_t* const argb = reference.GetRow8(y);
    const Span& horizontal_span = row_spans_[y];
    for (uint16_t x = horizontal_span.from; x < horizontal_span.to; ++x) {
      const uint8_t* const src = argb + x * 4;
      const int16_t xy[2] = { (int16_t)x, (int16_t)y };
      const int16_t samples[4] = { src[0], src[1], src[2], src[3] };
      opt.AddSample(xy, samples);
    }
  }
  float mtx[4][2 + 1];
  if (opt.Optimize(mtx)) {
    GetColor(mtx, v0.x, v0.y, &v0.color);
    GetColor(mtx, v1.x, v1.y, &v1.color);
    GetColor(mtx, v2.x, v2.y, &v2.color);
  } else {
    GetRawColor(reference, v0.x, v0.y, &v0.color);
    GetRawColor(reference, v1.x, v1.y, &v1.color);
    GetRawColor(reference, v2.x, v2.y, &v2.color);
  }
}

//------------------------------------------------------------------------------

static uint32_t Square(double a, uint8_t b) {
  const int32_t v = (int32_t)RoundToUInt8(a) - b;
  return (uint32_t)(v * v);
}

float Rasterizer::GetLossInScans(const ArgbBuffer& reference,
                                 VertexArgb v0, VertexArgb v1, VertexArgb v2,
                                 uint32_t* const num_pixels) const {
  assert(reference.width == img_w_ && reference.height == img_h_);
  double gradA[3], gradR[3], gradG[3], gradB[3];
  GetGradients(v0, v1, v2, gradA, gradR, gradG, gradB);
  double base_A = gradA[1] * vertical_span_.from + gradA[2];
  double base_R = gradR[1] * vertical_span_.from + gradR[2];
  double base_G = gradG[1] * vertical_span_.from + gradG[2];
  double base_B = gradB[1] * vertical_span_.from + gradB[2];

  uint32_t a_loss = 0, r_loss = 0, g_loss = 0, b_loss = 0;
  for (uint32_t y = vertical_span_.from; y < vertical_span_.to; ++y) {
    const uint8_t* const argb_reference = reference.GetRow8(y);
    const Span& horizontal_span = row_spans_[y];

    double cur_A = base_A + gradA[0] * horizontal_span.from;
    double cur_R = base_R + gradR[0] * horizontal_span.from;
    double cur_G = base_G + gradG[0] * horizontal_span.from;
    double cur_B = base_B + gradB[0] * horizontal_span.from;
    base_A += gradA[1];
    base_R += gradR[1];
    base_G += gradG[1];
    base_B += gradB[1];

    for (uint32_t x = horizontal_span.from; x < horizontal_span.to; ++x) {
      const uint8_t* argb_pixel_reference = argb_reference + 4 * x;
      a_loss += Square(cur_A, argb_pixel_reference[0]);
      r_loss += Square(cur_R, argb_pixel_reference[1]);
      g_loss += Square(cur_G, argb_pixel_reference[2]);
      b_loss += Square(cur_B, argb_pixel_reference[3]);
      cur_A += gradA[0];
      cur_R += gradR[0];
      cur_G += gradG[0];
      cur_B += gradB[0];
      ++*num_pixels;
    }
  }
  return a_loss + .3f * r_loss + .6f * g_loss + .1f * b_loss;
}

//------------------------------------------------------------------------------

void Rasterizer::RecordTriangleSpan(VertexArgb* const v0, VertexArgb* const v1,
                                    VertexArgb* const v2) {
  assert(row_spans_.size() == img_h_);

  // Triangles will be drawn top-left based, meaning in ranges [x,X) [y,Y).
  // Vertices coordinates are rescaled to fill the entire img_w_ x img_h_ area.
  v0->x = (v0->x * mult_w_) >> kBits;
  v1->x = (v1->x * mult_w_) >> kBits;
  v2->x = (v2->x * mult_w_) >> kBits;
  v0->y = (v0->y * mult_h_) >> kBits;
  v1->y = (v1->y * mult_h_) >> kBits;
  v2->y = (v2->y * mult_h_) >> kBits;

  vertical_span_.from = 0x7fff;  // Reset which rows will be drawn.
  vertical_span_.to = 0;

  RecordEdgeSpan(*v0, *v1);
  RecordEdgeSpan(*v1, *v2);
  RecordEdgeSpan(*v2, *v0);
}

//------------------------------------------------------------------------------

WP2Status Triangulation::FillTriangles(const PreviewData& preview,
                                       ArgbBuffer* const canvas,
                                       uint32_t* const num_pixels) const {
  Rasterizer rasterizer;
  WP2_CHECK_STATUS(rasterizer.Init(*canvas, preview));

  const Vector<VertexIndexedColor>& vertices = preview.vertices_;
  const Vector<AYCoCg19b>& palette = preview.palette_;
  for (const Triangle& triangle : triangles_) {
    const VertexIndexedColor vtx0 = GetVertex(vertices, triangle.v0);
    const VertexIndexedColor vtx1 = GetVertex(vertices, triangle.v1);
    const VertexIndexedColor vtx2 = GetVertex(vertices, triangle.v2);
    VertexArgb v0{vtx0.x, vtx0.y, ToArgb32b(palette[vtx0.color_index])};
    VertexArgb v1{vtx1.x, vtx1.y, ToArgb32b(palette[vtx1.color_index])};
    VertexArgb v2{vtx2.x, vtx2.y, ToArgb32b(palette[vtx2.color_index])};
    rasterizer.FillTriangle(v0, v1, v2, canvas, num_pixels);
  }
  return WP2_STATUS_OK;
}

WP2Status Triangulation::GetPartialLoss(const PreviewData& preview,
                                        const ArgbBuffer& reference,
                                        const Vector<Triangle>& triangles,
                                        float* const loss,
                                        uint32_t* const num_pixels) const {
  Rasterizer rasterizer;
  WP2_CHECK_STATUS(rasterizer.Init(reference, preview));
  *loss = 0.f;

  const Vector<VertexIndexedColor>& vertices = preview.vertices_;
  const Vector<AYCoCg19b>& palette = preview.palette_;
  for (const Triangle& triangle : triangles) {
    const VertexIndexedColor vtx0 = GetVertex(vertices, triangle.v0);
    const VertexIndexedColor vtx1 = GetVertex(vertices, triangle.v1);
    const VertexIndexedColor vtx2 = GetVertex(vertices, triangle.v2);
    VertexArgb v0{vtx0.x, vtx0.y, ToArgb32b(palette[vtx0.color_index])};
    VertexArgb v1{vtx1.x, vtx1.y, ToArgb32b(palette[vtx1.color_index])};
    VertexArgb v2{vtx2.x, vtx2.y, ToArgb32b(palette[vtx2.color_index])};
    *loss += rasterizer.GetLossInTriangle(reference, v0, v1, v2, num_pixels);
  }
  return WP2_STATUS_OK;
}

WP2Status Triangulation::GetLoss(const PreviewData& preview,
                                 const ArgbBuffer& reference,
                                 float* const loss,
                                 uint32_t* const num_pixels) const {
  return GetPartialLoss(preview, reference, triangles_, loss, num_pixels);
}

//------------------------------------------------------------------------------

WP2Status Triangulation::GetBestColors(const PreviewData& preview,
                                       const ArgbBuffer& reference,
                                       const Triangle& triangle,
                                       Argb32b colors[3]) const {
  Rasterizer rasterizer;
  WP2_CHECK_STATUS(rasterizer.Init(reference, preview));
  const Vector<VertexIndexedColor>& vertices = preview.vertices_;
  const VertexIndexedColor vtx0 = GetVertex(vertices, triangle.v0);
  const VertexIndexedColor vtx1 = GetVertex(vertices, triangle.v1);
  const VertexIndexedColor vtx2 = GetVertex(vertices, triangle.v2);
  VertexArgb v0 = {vtx0.x, vtx0.y, {0, 0, 0, 0}};
  VertexArgb v1 = {vtx1.x, vtx1.y, {0, 0, 0, 0}};
  VertexArgb v2 = {vtx2.x, vtx2.y, {0, 0, 0, 0}};
  rasterizer.CollectColors(reference, v0, v1, v2);
  colors[0] = v0.color;
  colors[1] = v1.color;
  colors[2] = v2.color;
  return WP2_STATUS_OK;
}

//------------------------------------------------------------------------------

}  // namespace WP2
