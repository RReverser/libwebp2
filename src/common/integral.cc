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
//  'Integral' struct to quickly compute mean and variance over rectangles.
//
// Author: Skal (pascal.massimino@gmail.com)
//

#include "src/common/integral.h"

namespace WP2 {

void Integral::GetArea(const Plane<input_t>& p, uint32_t x0, uint32_t y0,
                       integral_t* const M, integral_t* const M2) const {
  *M  = 0;
  *M2 = 0;
  for (uint32_t y = 0; y < size_; ++y) {
    for (uint32_t x = 0; x < size_; ++x) {
      const integral_t v = p.AtClamped(x0 * size_ + x, y0 * size_ + y);
      *M  += v;
      *M2 += v * v;
    }
  }
}

void Integral::GetArea(const YUVPlane& yuv, float uv_mix,
                       uint32_t x0, uint32_t y0,
                       integral_t* const M, integral_t* const M2) const {
  if (uv_mix == 0.f) return GetArea(yuv.Y, x0, y0, M, M2);
  *M = 0;
  *M2 = 0;
  for (uint32_t y = 0; y < size_; ++y) {
    for (uint32_t x = 0; x < size_; ++x) {
      const integral_t v_y = yuv.Y.AtClamped(x0 * size_ + x, y0 * size_ + y);
      const integral_t v_uv = yuv.U.AtClamped(x0 * size_ + x, y0 * size_ + y)
                            + yuv.V.AtClamped(x0 * size_ + x, y0 * size_ + y);
      const integral_t v = (integral_t)(v_y * (1.f - uv_mix) + v_uv * uv_mix);
      *M += v;
      *M2 += v * v;
    }
  }
}

WP2Status Integral::Allocate(uint32_t w, uint32_t h, uint32_t size) {
  w_ = w;
  h_ = h;
  size_ = size;
  wstride_ = w + 1;
  WP2_CHECK_ALLOC_OK(w1_.resize(wstride_ * (h + 1)));
  WP2_CHECK_ALLOC_OK(w2_.resize(wstride_ * (h + 1)));

  for (uint32_t i = 0; i <= w_; ++i) {
    w1_[i] = 0;
    w2_[i] = 0;
  }
  return WP2_STATUS_OK;
}

void Integral::SetW(uint32_t i, uint32_t j, integral_t M, integral_t M2) {
  const size_t off = pos(i, j);
  w1_[off] = M + w1(i, j - 1) + w1(i - 1, j) - w1(i - 1, j - 1);
  w2_[off] = M2 + w2(i, j - 1) + w2(i - 1, j) - w2(i - 1, j - 1);
}

void Integral::AddValues(const YUVPlane& yuv, float uv_mix) {
  assert(w_ <= (yuv.Y.w_ + size_ - 1) / size_);
  assert(h_ <= (yuv.Y.h_ + size_ - 1) / size_);
  for (uint32_t j = 1; j <= h_; ++j) {
    w1_[pos(0, j)] = 0;
    w2_[pos(0, j)] = 0;
    for (uint32_t i = 1; i <= w_; ++i) {
      integral_t M, M2;
      GetArea(yuv, uv_mix, i - 1, j - 1, &M, &M2);
      SetW(i, j, M, M2);
    }
  }
}

void Integral::AddValues(const Plane16& plane) {
  assert(w_ <= (plane.w_ + size_ - 1) / size_);
  assert(h_ <= (plane.h_ + size_ - 1) / size_);
  for (uint32_t j = 1; j <= h_; ++j) {
    w1_[pos(0, j)] = 0;
    w2_[pos(0, j)] = 0;
    for (uint32_t i = 1; i <= w_; ++i) {
      integral_t M, M2;
      GetArea(plane, i - 1, j - 1, &M, &M2);
      SetW(i, j, M, M2);
    }
  }
}

}    // namespace WP2
