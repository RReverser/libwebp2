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
//   Random matrices
//
// Author: Skal (pascal.massimino@gmail.com)

#include <cassert>
#include <cmath>
#include <cstdio>
#include <random>

#include "src/common/lossy/rnd_mtx.h"
#include "src/common/lossy/block.h"

using std::mt19937;

#if defined WP2_COLLECT_MTX_SCORES
WP2::Stats<uint32_t> scores("rnd_mtx-score-per-block-type");
#endif

namespace WP2 {

//------------------------------------------------------------------------------

static inline uint32_t Uniform(mt19937* const gen, uint32_t range) {
  const uint32_t max_range = mt19937::max() - mt19937::min();
  const double scale = 1. * range / max_range;
  return (uint32_t)std::round(scale * ((*gen)() - mt19937::min()));
}

static bool IsHighFrequency(uint32_t x, uint32_t y) {
  return (x + y > 3);
}

// TODO(skal): super slow!!
uint32_t FillRandomMatrix(uint32_t id, uint32_t sparsity, uint32_t freq_cut_off,
                          uint32_t sx, uint32_t sy, int16_t mtx[]) {
  mt19937 gen(id + 1 /*seed*/);
  assert(mtx != nullptr && freq_cut_off > 0);
  uint32_t num_ones = 0;

  const double norm = 1. / (sx * sx + sy * sy);
  const double alpha = 0.2 * 1024. / freq_cut_off;
  for (uint32_t n = 0; n < sx * sy; ++n) {
    const uint32_t i = n % sx;
    const uint32_t j = n / sx;
    uint32_t v = 0;
    if (IsHighFrequency(i, j)) {
      v = Uniform(&gen, 0xffu);
      const float r = pow(norm * (i * i + j * j), 0.3);
      const float p = exp(-r * alpha);
      v *= p;
      v = (v > sparsity) ? (v & 1 ? -1 : 1) : 0;
    }
    mtx[n] = v;
    num_ones += (v != 0);
  }
  return num_ones;
}

static int16_t DecideValue(int16_t v0, int16_t v1) {
  if (v0 != 0 || v1 != 0) {
    return (v0 < 0) ? -1 : 1;
  }
  return 0;
}

uint32_t DownscaleMatrix(const int16_t* src, int16_t* dst,
                         uint32_t sx, uint32_t sy, bool vertically) {
  uint32_t num_ones = 0;
  for (uint32_t y = 0; y < sy; ++y) {
    if (!vertically) {
      for (uint32_t x = 0; x < sx; ++x, ++dst) {
        const int16_t v0 = src[x + sx * (2 * y + 0)];
        const int16_t v1 = src[x + sx * (2 * y + 1)];
        *dst = IsHighFrequency(x, y) ? DecideValue(v0, v1) : 0;
        num_ones += (*dst != 0);
      }
    } else {
      for (uint32_t x = 0; x < sx; ++x, ++dst) {
        const int16_t v0 = src[2 * x + 0 + (2 * sx) * y];
        const int16_t v1 = src[2 * x + 1 + (2 * sx) * y];
        *dst = IsHighFrequency(x, y) ? DecideValue(v0, v1) : 0;
        num_ones += (*dst != 0);
      }
    }
  }
  return num_ones;
}

void AddRandomMatrix(uint32_t size, const int16_t rnd[], int16_t coeffs[]) {
  for (uint32_t i = 0; i < size; ++i) {
    if (rnd[i]) {
      if (coeffs[i]) {
        const int16_t v = std::abs(coeffs[i]) + 1;
        coeffs[i] = (coeffs[i] < 0) ? -v : v;
      } else {
        coeffs[i] = rnd[i];
      }
    }
  }
}

void SubtractRandomMatrix(uint32_t size, const int16_t rnd[], int16_t coeffs[],
                          uint32_t* const num_coeffs) {
  *num_coeffs = 0;
  for (uint32_t i = 0; i < size; ++i) {
    const uint32_t c = std::abs(coeffs[i]);
    if (c == 0) {  /* nothing to do */
    } else if (c == 1) {
      if (rnd[i]) {
        coeffs[i] = 0;  // replace +1/-1 with +1/-1 from rnd[]
      }
    } else if (rnd[i]) {   // large coeffs
      const int16_t v = c - 1;
      coeffs[i] = (coeffs[i] < 0) ? -v : v;
    }
    if (coeffs[i] != 0) *num_coeffs = i + 1;
  }
}

uint32_t RandomMatrixScore(const int16_t rnd[],
                           const int16_t coeffs[], uint32_t sx, uint32_t sy) {
  int32_t count = 0;
  for (uint32_t i = 0; i < sx * sy; ++i) {
    if (rnd[i]) {
      const uint32_t v = std::abs(coeffs[i]);
      count += (v == 0) ? 0 : (v == 1) ? 2 : 0;
    }
  }
  return (count < 0) ? 0 : count + 5;
}

void PrintRandomMatrix(uint32_t sx, uint32_t sy, const int16_t m[],
                       bool full_coeffs) {
  for (uint32_t j = 0; j < sy; ++j) {
    for (uint32_t i = 0; i < sx; ++i) {
      const int16_t v = m[i + j * sx];
      if (full_coeffs) {
        printf("%3d", v);
      } else {
        assert(v >= -1 && v <= 1);
        printf("%c", "- +"[1 + v]);
      }
    }
    printf("\n");
  }
  printf("----- (%dx%d)\n", sx, sy);
}

//------------------------------------------------------------------------------

template<class T>
static WP2Status InitMtx(T* const mtx, uint32_t size,
                         uint32_t sparsity, uint32_t freq_cut_off) {
  WP2_CHECK_ALLOC_OK(mtx->resize(size));
  if (size >= kRndMtxBadIdx) return WP2_STATUS_INVALID_PARAMETER;
  for (uint32_t i = 0; i < size; ++i) {
    (*mtx)[i].FillRandom(i, sparsity, freq_cut_off);
  }
  return WP2_STATUS_OK;
}

WP2Status RndMtxSet::Init(uint32_t num_mtx, const QuantMtx& quant) {
  num_mtx_ = num_mtx;

  // TODO(skal): deduce or learn sparsity/cut-off from quant
  (void)quant;
  // 'sparsity' controls the probability of being non-zero
  constexpr uint32_t sparsity = 13;  // good values at default quality
  // cutoff controls the size of the zero'd hi-frequencies
  constexpr uint32_t cut_off = 40;

  // initial source of random matrix set
  WP2_CHECK_STATUS(InitMtx(&mtx32x32_, num_mtx_, sparsity, cut_off));
  // downscaled version (for speed)
  WP2_CHECK_STATUS(Downscale(mtx32x32_, &mtx16x32_));
  WP2_CHECK_STATUS(Downscale(mtx32x32_, &mtx32x16_));
  WP2_CHECK_STATUS(Downscale(mtx32x16_, &mtx16x16_));
  WP2_CHECK_STATUS(Downscale(mtx16x16_, &mtx8x16_));
  WP2_CHECK_STATUS(Downscale(mtx16x16_, &mtx16x8_));
  WP2_CHECK_STATUS(Downscale(mtx16x8_, &mtx8x8_));
  WP2_CHECK_STATUS(Downscale(mtx8x8_, &mtx4x8_));
  WP2_CHECK_STATUS(Downscale(mtx8x8_, &mtx8x4_));
  WP2_CHECK_STATUS(Downscale(mtx8x4_, &mtx4x4_));
  // TODO(skal): remove duplicates, esp. for small dimensions

  return WP2_STATUS_OK;
}

void RndMtxSet::Reset() {
  num_mtx_ = 0;
  mtx4x4_.reset();
  mtx4x8_.reset();
  mtx8x4_.reset();
  mtx8x8_.reset();
  mtx16x8_.reset();
  mtx8x16_.reset();
  mtx16x16_.reset();
  mtx32x16_.reset();
  mtx16x32_.reset();
  mtx32x32_.reset();
}

WP2Status RndMtxSet::CopyFrom(const RndMtxSet& src) {
  Reset();
  num_mtx_ = src.num_mtx_;
  if (num_mtx_ > 0) {
    WP2_CHECK_ALLOC_OK(mtx4x4_.copy_from(src.mtx4x4_));
    WP2_CHECK_ALLOC_OK(mtx4x8_.copy_from(src.mtx4x8_));
    WP2_CHECK_ALLOC_OK(mtx8x4_.copy_from(src.mtx8x4_));
    WP2_CHECK_ALLOC_OK(mtx8x8_.copy_from(src.mtx8x8_));
    WP2_CHECK_ALLOC_OK(mtx16x8_.copy_from(src.mtx16x8_));
    WP2_CHECK_ALLOC_OK(mtx8x16_.copy_from(src.mtx8x16_));
    WP2_CHECK_ALLOC_OK(mtx16x16_.copy_from(src.mtx16x16_));
    WP2_CHECK_ALLOC_OK(mtx32x16_.copy_from(src.mtx32x16_));
    WP2_CHECK_ALLOC_OK(mtx16x32_.copy_from(src.mtx16x32_));
    WP2_CHECK_ALLOC_OK(mtx32x32_.copy_from(src.mtx32x32_));
  }
  return WP2_STATUS_OK;
}

template<class T>
bool CodedBlock::CheckUseMtx(BlockSize dim, const T& mtx) {
  if (blk_.dim() == dim && !GetCodingParams(kYChannel)->split_tf) {
    // TODO(skal): later: mtx_[1] = ...
    mtx_[0] = FindAndSubtractBestMatch(mtx, coeffs_[kYChannel][0],
                                       &num_coeffs_[kYChannel][0]);
    return (mtx_[0] != kRndMtxBadIdx);
  }
  return false;
}

void RndMtxSet::DecideUseRndMtx(CodedBlock* const cb) const {
  assert(cb->mtx_set_ == this);
  cb->use_mtx_ = cb->CheckUseMtx(BLK_4x4, mtx4x4_) ||
                 cb->CheckUseMtx(BLK_4x8, mtx4x8_) ||
                 cb->CheckUseMtx(BLK_8x4, mtx8x4_) ||
                 cb->CheckUseMtx(BLK_8x8, mtx8x8_) ||
                 cb->CheckUseMtx(BLK_16x8, mtx16x8_) ||
                 cb->CheckUseMtx(BLK_8x16, mtx8x16_) ||
                 cb->CheckUseMtx(BLK_16x16, mtx16x16_) ||
                 cb->CheckUseMtx(BLK_32x16, mtx32x16_) ||
                 cb->CheckUseMtx(BLK_16x32, mtx16x32_) ||
                 cb->CheckUseMtx(BLK_32x32, mtx32x32_);
}

bool RndMtxSet::IsOk(uint32_t id, BlockSize dim) const {
  if (dim == BLK_4x4 && id < mtx4x4_.size()) return true;
  if (dim == BLK_4x8 && id < mtx4x8_.size()) return true;
  if (dim == BLK_8x4 && id < mtx8x4_.size()) return true;
  if (dim == BLK_8x8 && id < mtx8x8_.size()) return true;
  if (dim == BLK_16x8 && id < mtx16x8_.size()) return true;
  if (dim == BLK_8x16 && id < mtx8x16_.size()) return true;
  if (dim == BLK_16x16 && id < mtx16x16_.size()) return true;
  if (dim == BLK_32x16 && id < mtx32x16_.size()) return true;
  if (dim == BLK_16x32 && id < mtx16x32_.size()) return true;
  if (dim == BLK_32x32 && id < mtx32x32_.size()) return true;
  return false;
}

void RndMtxSet::AddMatrix(uint32_t id, BlockSize dim, int16_t dst[]) const {
  if (dim == BLK_4x4) mtx4x4_[id].AddTo(dst);
  if (dim == BLK_4x8) mtx4x8_[id].AddTo(dst);
  if (dim == BLK_8x4) mtx8x4_[id].AddTo(dst);
  if (dim == BLK_8x8) mtx8x8_[id].AddTo(dst);
  if (dim == BLK_16x8) mtx16x8_[id].AddTo(dst);
  if (dim == BLK_8x16) mtx8x16_[id].AddTo(dst);
  if (dim == BLK_16x16) mtx16x16_[id].AddTo(dst);
  if (dim == BLK_32x16) mtx32x16_[id].AddTo(dst);
  if (dim == BLK_16x32) mtx16x32_[id].AddTo(dst);
  if (dim == BLK_32x32) mtx32x32_[id].AddTo(dst);
}

}    // namespace WP2
