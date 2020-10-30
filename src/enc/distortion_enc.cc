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
// Calculation of distortion
//
// Author: Skal (pascal.massimino@gmail.com)

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdlib>

#include "src/dsp/dsp.h"
#include "src/utils/csp.h"
#include "src/utils/plane.h"
#include "src/utils/utils.h"
#include "src/wp2/base.h"

namespace {

// Max value returned in case of exact similarity.
const double kNoDistortion_dB = 99.;
// If there is the slightest difference, clamp to a value strictly lower than
// 99 dB because it is used as "pristine lossless level".
const double kMinDistortion_dB = kNoDistortion_dB - 0.1;

// 'v' contains the total sum of the squared differences.
// 'size' is the number of pixels and 'max' is the maximum difference value.
double GetPSNR(double v, double size, double max) {
  assert(max > 0.);
  assert(v >= 0.);
  if (v == 0.) return kNoDistortion_dB;
  assert(size > 0.);
  v /= size * max * max;
  if (v <= 0.) return kMinDistortion_dB;  // For very small distortion.
  return std::min(-4.3429448 * std::log(v), kMinDistortion_dB);
}

// 'v' is closer to 0 for a higher distortion. 'size' is the number of pixels.
double GetLogSSIM(double v, double size) {
  assert(v >= 0.);
  assert(v <= size);
  if (v >= size) return kNoDistortion_dB;
  assert(size > 0.);
  v = 1. - v / size;
  if (v <= 0.) return kMinDistortion_dB;  // For very small distortion.
  return std::min(-10.0 * log10(v), kMinDistortion_dB);
}

// -----------------------------------------------------------------------------

template<uint32_t channel_step>
void WP2SSIMGetClipped(const uint8_t* const src1, size_t step1,
                       const uint8_t* const src2, size_t step2,
                       uint32_t xo, uint32_t yo, uint32_t W, uint32_t H,
                       WP2DistoStats* const stats);

template<>
void WP2SSIMGetClipped<4>(const uint8_t* const src1, size_t step1,
                          const uint8_t* const src2, size_t step2,
                          uint32_t xo, uint32_t yo, uint32_t W, uint32_t H,
                          WP2DistoStats* const stats) {
  return WP2SSIMGetClipped4x8u(src1, step1, src2, step2, xo, yo, W, H, stats);
}

template<>
void WP2SSIMGetClipped<1>(const uint8_t* const src1, size_t step1,
                          const uint8_t* const src2, size_t step2,
                          uint32_t xo, uint32_t yo, uint32_t W, uint32_t H,
                          WP2DistoStats* const stats) {
  return WP2SSIMGetClipped8u(src1, step1, src2, step2, xo, yo, W, H, stats);
}

template<uint32_t channel_step>
void WP2SSIMGetClipped(const int16_t* const src1, size_t step1,
                       const int16_t* const src2, size_t step2,
                       uint32_t xo, uint32_t yo, uint32_t W, uint32_t H,
                       WP2DistoStats* const stats) {
  static_assert(channel_step == 1, "invalid channel step");
  return WP2SSIMGetClipped16s(src1, step1, src2, step2, xo, yo, W, H,
                              stats);
}

// -----------------------------------------------------------------------------

template<uint32_t channel_step>
void WP2SSIMGet(const uint8_t* const src1, size_t step1,
                const uint8_t* const src2, size_t step2,
                WP2DistoStats* const stats);

template<>
void WP2SSIMGet<4>(const uint8_t* const src1, size_t step1,
                   const uint8_t* const src2, size_t step2,
                   WP2DistoStats* const stats) {
  return WP2SSIMGet4x8u(src1, step1, src2, step2, stats);
}

template<>
void WP2SSIMGet<1>(const uint8_t* const src1, size_t step1,
                   const uint8_t* const src2, size_t step2,
                   WP2DistoStats* const stats) {
  return WP2SSIMGet8u(src1, step1, src2, step2, stats);
}

template<uint32_t channel_step>
void WP2SSIMGet(const int16_t* const src1, size_t step1,
                const int16_t* const src2, size_t step2,
                WP2DistoStats* const stats) {
  static_assert(channel_step == 1, "invalid channel step");
  return WP2SSIMGet16s(src1, step1, src2, step2, stats);
}

// -----------------------------------------------------------------------------

typedef double (*WP2SSIMGetFromStatsFunc)(uint32_t bit_depth,
                                          const WP2DistoStats& stats);

template <class T, uint32_t channel_step,
          WP2SSIMGetFromStatsFunc get_from_stats = WP2SSIMCalculation>
double AccumulateSSIM(uint32_t w, uint32_t h, uint32_t bit_depth,
                      const T* const src, uint32_t src_step,
                      const T* const ref, uint32_t ref_step) {
  const uint32_t pad = kWP2SSIMMargin - kWP2SSIMKernel;
  const uint32_t w0 = std::min(w, kWP2SSIMKernel);
  const uint32_t w1 = std::max(pad, w) - pad;
  const uint32_t h0 = std::min(h, kWP2SSIMKernel);
  const uint32_t h1 = std::max(kWP2SSIMKernel, h) - kWP2SSIMKernel;
  uint32_t x, y;
  double sum = 0.;
  for (y = 0; y < h0; ++y) {
    for (x = 0; x < w; ++x) {
      WP2DistoStats stats;
      WP2SSIMGetClipped<channel_step>(src, src_step, ref, ref_step,
                                      x, y, w, h, &stats);
      sum += get_from_stats(bit_depth, stats);
    }
  }
  for (; y < h1; ++y) {
    for (x = 0; x < w0; ++x) {
      WP2DistoStats stats;
      WP2SSIMGetClipped<channel_step>(src, src_step, ref, ref_step,
                                      x, y, w, h, &stats);
      sum += get_from_stats(bit_depth, stats);
    }
    const size_t y_src_off = (y - kWP2SSIMKernel) * src_step;
    const size_t y_ref_off = (y - kWP2SSIMKernel) * ref_step;
    for (; x < w1; ++x) {
      const size_t off1 = channel_step * (x - kWP2SSIMKernel) + y_src_off;
      const size_t off2 = channel_step * (x - kWP2SSIMKernel) + y_ref_off;
      WP2DistoStats stats;
      WP2SSIMGet<channel_step>(src + off1, src_step, ref + off2, ref_step,
                               &stats);
      sum += get_from_stats(bit_depth, stats);
    }
    for (; x < w; ++x) {
      WP2DistoStats stats;
      WP2SSIMGetClipped<channel_step>(src, src_step, ref, ref_step,
                                      x, y, w, h, &stats);
      sum += get_from_stats(bit_depth, stats);
    }
  }
  for (; y < h; ++y) {
    for (x = 0; x < w; ++x) {
      WP2DistoStats stats;
      WP2SSIMGetClipped<channel_step>(src, src_step, ref, ref_step,
                                      x, y, w, h, &stats);
      sum += get_from_stats(bit_depth, stats);
    }
  }
  return sum;
}

template<WP2SSIMGetFromStatsFunc get_from_stats>
double AccumulateSSIM(const WP2::ArgbBuffer& a, const WP2::ArgbBuffer& b,
                      uint32_t channel) {
  return AccumulateSSIM<uint8_t, /*CHANNEL_STRIDE=*/4, get_from_stats>(
      a.width, a.height, /*bit_depth=*/8,
      (const uint8_t*)a.GetRow(0) + channel, a.stride,
      (const uint8_t*)b.GetRow(0) + channel, b.stride);
}

double AccumulateSSIM(const WP2::ArgbBuffer& a, const WP2::ArgbBuffer& b,
                      uint32_t channel) {
  return AccumulateSSIM<WP2SSIMCalculation>(a, b, channel);
}

double AccumulateCsSSIM(const WP2::ArgbBuffer& a, const WP2::ArgbBuffer& b,
                        uint32_t channel) {
  return AccumulateSSIM<WP2CsSSIMCalculation>(a, b, channel);
}

}  // namespace

//------------------------------------------------------------------------------
// MSSSIM tools

static constexpr uint32_t kNumMSSSIMPasses = 5;
static constexpr float kMSSSIMExponents[kNumMSSSIMPasses] = {
  // Simoncelli & al: https://ece.uwaterloo.ca/~z70wang/publications/msssim.pdf
  //    original: 0.0448f, 0.2856f, 0.3001f, 0.2363f, 0.13333f
  // Modified exponents that weight more the 1rst layer, and are more regular:
  0.35f, 0.23f, 0.19f, 0.15f, 0.08f
};

//------------------------------------------------------------------------------
// local-min distortion
//
// For every pixel in the *reference* picture, we search for the local best
// match in the compressed image. This is not a symmetrical measure.

namespace {

// Search radius. Shouldn't be too large.
constexpr uint32_t kRadius = 2u;

double AccumulateLSIM(const uint8_t* const src, uint32_t src_stride,
                      const uint8_t* const ref, uint32_t ref_stride,
                      uint32_t x, uint32_t y, uint32_t w, uint32_t h) {
  const uint32_t y_0 = std::max(y, kRadius) - kRadius;
  const uint32_t y_1 = std::min(h - 1, y + kRadius);
  const uint32_t x_0 = std::max(x, kRadius) - kRadius;
  const uint32_t x_1 = std::min(w - 1, x + kRadius);
  double best_sse = 255. * 255.;
  const double value = (double)ref[y * ref_stride + 4 * x];
  for (uint32_t j = y_0; j <= y_1; ++j) {
    const uint8_t* const s = src + j * src_stride;
    for (uint32_t i = x_0; i <= x_1; ++i) {
      const double diff = s[4 * i] - value;
      const double sse = diff * diff;
      if (sse < best_sse) best_sse = sse;
    }
  }
  return best_sse;
}

double AccumulateLSIM(const WP2::ArgbBuffer& src_buffer,
                      const WP2::ArgbBuffer& ref_buffer,
                      uint32_t channel) {
  const uint32_t w = src_buffer.width, h = ref_buffer.height;
  const uint8_t* const src = (const uint8_t*)src_buffer.GetRow(0) + channel;
  const size_t src_stride = src_buffer.stride;
  const uint8_t* const ref = (const uint8_t*)ref_buffer.GetRow(0) + channel;
  const size_t ref_stride = ref_buffer.stride;
  double total_sse = 0.;
  for (uint32_t y = 0; y < h; ++y) {
    for (uint32_t x = 0; x < w; ++x) {
      total_sse += AccumulateLSIM(src, src_stride, ref, ref_stride, x, y, w, h);
    }
  }
  return total_sse;
}

//------------------------------------------------------------------------------
// PSNR-HVS-M

// Implementation based on libaom's, originally written by Gregory Maxwell.
/*
 [1] Nikolay Ponomarenko, Flavia Silvestri, Karen Egiazarian, Marco Carli,
   Jaakko Astola, Vladimir Lukin, "On between-coefficient contrast masking
   of DCT basis functions", CD-ROM Proceedings of the Third
   International Workshop on Video Processing and Quality Metrics for Consumer
   Electronics VPQM-07, Scottsdale, Arizona, USA, 25-26 January, 2007, 4 p.
*/

constexpr int32_t Square(int32_t v) { return v * v; }

// The "CSF" masks, normalized to a maximum of 1.
constexpr float kMaskY[64] = {
  0.707107, 1.000000, 0.910460, 0.647842, 0.437644, 0.296179, 0.203577, .142570,
  1.000000, 0.848508, 0.894231, 0.736766, 0.537328, 0.379415, 0.267584, .190557,
  0.910460, 0.894231, 0.586549, 0.476847, 0.382396, 0.292942, 0.219082, .162654,
  0.647842, 0.736766, 0.476847, 0.337452, 0.264452, 0.210943, 0.166115, .129150,
  0.437644, 0.537328, 0.382396, 0.264452, 0.196055, 0.154089, 0.123575, .099098,
  0.296179, 0.379415, 0.292942, 0.210943, 0.154089, 0.118036, 0.093888, .076012,
  0.203577, 0.267584, 0.219082, 0.166115, 0.123575, 0.093888, 0.073737, .059452,
  0.142570, 0.190557, 0.162654, 0.129150, 0.099098, 0.076012, 0.059452, .047632,
};
constexpr float kMaskU[64] = {
  0.776648, 1.000000, 0.480685, 0.467268, 0.426770, 0.364938, 0.303670, .249968,
  1.000000, 0.644234, 0.493198, 0.561579, 0.540895, 0.477208, 0.404920, .337658,
  0.480685, 0.493198, 0.397731, 0.417047, 0.419163, 0.390151, 0.345352, .297155,
  0.467268, 0.561579, 0.417047, 0.350023, 0.325845, 0.305370, 0.278533, .247362,
  0.426770, 0.540895, 0.419163, 0.325845, 0.274940, 0.246065, 0.223518, .201486,
  0.364938, 0.477208, 0.390151, 0.305370, 0.246065, 0.209154, 0.184641, .165418,
  0.303670, 0.404920, 0.345352, 0.278533, 0.223518, 0.184641, 0.158178, .139126,
  0.249968, 0.337658, 0.297155, 0.247362, 0.201486, 0.165418, 0.139126, .120098
};
constexpr float kMaskV[64] = {
  0.776648, 1.000000, 0.480685, 0.422929, 0.386274, 0.330309, 0.274855, .226248,
  1.000000, 0.644234, 0.446398, 0.508290, 0.489569, 0.431926, 0.366497, .305618,
  0.480685, 0.446398, 0.359990, 0.377473, 0.379388, 0.353129, 0.312582, .268958,
  0.422929, 0.508290, 0.377473, 0.316810, 0.294926, 0.276394, 0.252103, .223890,
  0.386274, 0.489569, 0.379388, 0.294926, 0.248851, 0.222716, 0.202308, .182367,
  0.330309, 0.431926, 0.353129, 0.276394, 0.222716, 0.189308, 0.167120, .149721,
  0.274855, 0.366497, 0.312582, 0.252103, 0.202308, 0.167120, 0.143168, .125925,
  0.226248, 0.305618, 0.268958, 0.223890, 0.182367, 0.149721, 0.125925, .108702
};

const float kCos[8][8] = {
  {1.000000e+00,  1.000000e+00,  1.000000e+00,  1.000000e+00,
      1.000000e+00, 1.000000e+00, 1.000000e+00, 1.000000e+00, },
  {9.807853e-01,  8.314696e-01,  5.555702e-01,  1.950903e-01,
     -1.950903e-01, -5.555702e-01, -8.314696e-01, -9.807853e-01, },
  {9.238795e-01,  3.826834e-01, -3.826834e-01, -9.238795e-01,
     -9.238795e-01, -3.826834e-01, 3.826834e-01, 9.238795e-01, },
  {8.314696e-01, -1.950903e-01, -9.807853e-01, -5.555702e-01,
      5.555702e-01, 9.807853e-01, 1.950903e-01, -8.314696e-01, },
  {7.071068e-01, -7.071068e-01, -7.071068e-01,  7.071068e-01,
      7.071068e-01, -7.071068e-01, -7.071068e-01, 7.071068e-01, },
  {5.555702e-01, -9.807853e-01,  1.950903e-01,  8.314696e-01,
     -8.314696e-01, -1.950903e-01, 9.807853e-01, -5.555702e-01, },
  {3.826834e-01, -9.238795e-01,  9.238795e-01, -3.826834e-01,
     -3.826834e-01, 9.238795e-01, -9.238795e-01, 3.826834e-01, },
  {1.950903e-01, -5.555702e-01,  8.314696e-01, -9.807853e-01,
      9.807853e-01, -8.314696e-01, 5.555702e-01, -1.950903e-01, },
};

void Dct8(int32_t dct[64]) {
#if 0
  for (uint32_t j = 0; j < 8; ++j) {
    printf("  {");
    for (uint32_t i = 0; i < 8; ++i) {
      printf("%.6e, ", cos(((2.0 * i + 1) * j * M_PI) / 16.));
    }
    printf("},\n");
  }
  exit(1);
#endif
  int32_t tmp[64];
  for (uint32_t k = 0, v = 0; v < 8; ++v) {
    for (uint32_t u = 0; u < 8; ++u, ++k) {
      tmp[k] = WP2SlowDct8x8(dct, kCos[u], kCos[v]);
    }
  }
  memcpy(dct, tmp, sizeof(tmp));
}

float CalcPSNRHVS(const WP2::Plane16& a, const WP2::Plane16& b,
                  const float mask[64]) {
  assert(a.w_ == b.w_);
  assert(a.h_ == b.h_);

  WP2TransformInit();

  const uint32_t step = 7;
  float total = 0.;
  uint32_t count = 0;
  for (uint32_t y = 0; y + 8 <= a.h_; y += step) {
    for (uint32_t x = 0; x + 8 <= a.w_; x += step) {
      int32_t dct[2][64];
      float means[2] = { 0, 0 };  // mean for each sources
      // and same for all 2x2 sub-blocks
      float sub_means[2][4] = { {0}, {0} };
      // Import and compute means.
      for (uint32_t j = 0, k = 0; j < 8; ++j) {
        const int16_t* const src_a = &a.At(x, y + j);
        const int16_t* const src_b = &b.At(x, y + j);
        for (uint32_t i = 0; i < 8; ++i, ++k) {
          const uint32_t sub_blk = ((j & 4) >> 1) + (i >> 2);
          dct[0][k] = src_a[i];
          dct[1][k] = src_b[i];
          for (uint32_t s : {0, 1}) {
            means[s] += dct[s][k];
            sub_means[s][sub_blk] += dct[s][k];
          }
        }
      }
      for (uint32_t s : {0, 1}) {
        means[s] /= 64.;
        for (int32_t sub = 0; sub <= 3; ++sub) sub_means[s][sub] /= 16.;
      }

      // Compute variances.
      float vars[2] = { 0, 0 };
      float sub_vars[2][4] = { {0}, {0} };
      for (uint32_t j = 0, k = 0; j < 8; ++j) {
        for (uint32_t i = 0; i < 8; ++i, ++k) {
          const uint32_t sub_blk = ((j & 4) >> 1) + (i >> 2);
          for (uint32_t s : {0, 1}) {
            vars[s] += Square(dct[s][k] - means[s]);
            sub_vars[s][sub_blk] += Square(dct[s][k] - sub_means[s][sub_blk]);
          }
        }
      }

      float masks[2];
      for (uint32_t s : {0, 1}) {
        Dct8(dct[s]);
        masks[s] = 0.;
        if (vars[s] > 0.) {
          // equation (2): Em(D) = Ew(D) * delta(D) / 16
          vars[s] = (sub_vars[s][0] + sub_vars[s][1] +
                     sub_vars[s][2] + sub_vars[s][3]) * 0.25 / vars[s];
          // equation (1), C_00 = 0
          uint64_t sum = 0;
          for (uint32_t k = 1; k < 64; ++k) {
            sum += (uint64_t)dct[s][k] * dct[s][k] * mask[k];
          }
          masks[s] = std::sqrt(sum * vars[s]) / 256. / 32.;   // E-norm
        }
      }
      const float final_mask = std::max(masks[0], masks[1]);
      for (uint32_t k = 0; k < 64; ++k) {
        float err = std::fabs(dct[0][k] - dct[1][k]);
        if (k > 0) err = std::max(err * mask[k] - final_mask, 0.f);
        total += err * err;
        ++count;
      }
    }
  }
  if (count > 0) total /= count;
  return total;
}

}    // namespace

//------------------------------------------------------------------------------

namespace WP2 {

WP2Status ArgbBuffer::GetDistortion(const ArgbBuffer& ref,
                                    const Rectangle& window,
                                    MetricType metric_type,
                                    float result[5]) const {
  ArgbBuffer view, ref_view;
  WP2_CHECK_STATUS(view.SetView(*this, window));
  WP2_CHECK_STATUS(ref_view.SetView(ref, window));
  return view.GetDistortion(ref_view, metric_type, result);
}

WP2Status ArgbBuffer::GetDistortion(const ArgbBuffer& ref,
                                    MetricType metric_type,
                                    float result[5]) const {
  return GetDistortion(ref, metric_type, /*include_alpha_in_all=*/true, result);
}

WP2Status ArgbBuffer::GetDistortion(const ArgbBuffer& ref,
                                    MetricType metric_type,
                                    bool include_alpha_in_all,
                                    float result[5]) const {
  WP2_CHECK_OK(!IsEmpty() && !ref.IsEmpty() && (result != nullptr),
               WP2_STATUS_NULL_PARAMETER);
  WP2_CHECK_OK((width == ref.width) && (height == ref.height),
               WP2_STATUS_BAD_DIMENSION);
  result[0] = result[1] = result[2] = result[3] = result[4] = 0.f;

  const double pix_size = width * height;
  const double max = 255.;  // Maximum value per pixel per channel.
  double total_distortion = 0.;
  const bool has_transparency = HasTransparency();
  const bool include_alpha = has_transparency && include_alpha_in_all;
  const uint32_t pix_channel_size = pix_size * (include_alpha ? 4 : 3);
  if (metric_type == PSNR) {
    WP2PSNRInit();
    uint64_t result_64b[4] = {0, 0, 0, 0};
    for (uint32_t y = 0; y < height; ++y) {
      if (WP2FormatBpc(ref.format) == 1) {
        const uint8_t* const ptr1 = GetRow8(y);
        const uint8_t* const ptr2 = ref.GetRow8(y);
        WP2SumSquaredError4x8u(ptr1, ptr2, width, result_64b);
      } else {
        assert(WP2FormatBpc(ref.format) == 2);
        const uint16_t* const ptr1 = GetRow16(y);
        const uint16_t* const ptr2 = ref.GetRow16(y);
        WP2SumSquaredError2x16u(ptr1, ptr2, width, result_64b);
      }
    }
    for (uint32_t c = 0; c < 4; ++c) {
      if (c > 0 || include_alpha) total_distortion += (float)result_64b[c];
      result[c] = (float)GetPSNR(result_64b[c], pix_size, max);
    }
    result[4] = (float)GetPSNR(total_distortion, pix_channel_size, max);
  } else if (metric_type == SSIM) {
    assert(WP2FormatBpc(ref.format) == 1);  // not supported yet
    WP2SSIMInit();
    for (uint32_t c = 0; c < 4; ++c) {
      const double disto = AccumulateSSIM(*this, ref, c);
      if (c > 0 || include_alpha) total_distortion += disto;
      result[c] = (float)GetLogSSIM(disto, pix_size);
    }
    result[4] = (float)GetLogSSIM(total_distortion, pix_channel_size);
  } else if (metric_type == MSSSIM) {
    assert(WP2FormatBpc(ref.format) == 1);  // not supported yet
    WP2SSIMInit();
    ArgbBuffer buf1, buf2;  // working copies, to be down-scaled
    WP2_CHECK_STATUS(buf1.CopyFrom(ref));
    WP2_CHECK_STATUS(buf2.CopyFrom(*this));
    double acc[4] = {1., 1., 1., 1.};
    for (uint32_t pass = 0; pass < kNumMSSSIMPasses - 1; ++pass) {
      const double norm = 1. / (buf1.width * buf1.height);
      for (uint32_t c = 0; c < 4; ++c) {
        const double sum = norm * AccumulateCsSSIM(buf1, buf2, c);
        acc[c] *= pow(sum, kMSSSIMExponents[pass]);
      }
      buf1.SimpleHalfDownsample();
      buf2.SimpleHalfDownsample();
    }
    double total = 0.;
    for (uint32_t c = 0; c < 4; ++c) {
      const double norm = 1. / (buf1.width * buf1.height);
      const double sum = norm * AccumulateSSIM(buf1, buf2, c);
      acc[c] *= pow(sum, kMSSSIMExponents[kNumMSSSIMPasses - 1]);
      if (c > 0 || include_alpha) total += acc[c];
      result[c] = GetLogSSIM(acc[c], 1.);
    }
    result[4] = GetLogSSIM(total, include_alpha ? 4. : 3.);
  } else if (metric_type == LSIM) {
    for (uint32_t c = 0; c < 4; ++c) {
      const double disto = AccumulateLSIM(*this, ref, c);
      if (c > 0 || include_alpha) total_distortion += disto;
      result[c] = (float)GetPSNR(disto, pix_size, max);
    }
    result[4] = (float)GetPSNR(total_distortion, pix_channel_size, max);
  } else {
    // We are working in YUV space.
    CSPTransform transf;
    transf.InitYCbCr();
    YUVPlane yuv0, yuv1;
    WP2_CHECK_STATUS(yuv0.Import(*this, has_transparency, transf,
                                 /*resize_if_needed=*/true));
    WP2_CHECK_STATUS(
        yuv1.Import(ref, has_transparency, transf, /*resize_if_needed=*/true));
    if (metric_type == PSNRHVS) {
      // note: for alpha, we use the Y-mask.
      result[0] = has_transparency ? CalcPSNRHVS(yuv0.A, yuv1.A, kMaskY) : 0;
      result[1] = CalcPSNRHVS(yuv0.Y, yuv1.Y, kMaskY);
      result[2] = CalcPSNRHVS(yuv0.U, yuv1.U, kMaskU);
      result[3] = CalcPSNRHVS(yuv0.V, yuv1.V, kMaskV);
      for (uint32_t c = 0; c < 4; ++c) {
        if (c > 0 || include_alpha) total_distortion += result[c];
        result[c] = (float)GetPSNR(result[c], pix_size, max);
      }
      result[4] = (float)GetPSNR(total_distortion, pix_channel_size, max);
    } else if (metric_type == PSNR_YUV) {
      WP2_CHECK_STATUS(yuv0.GetDistortion(
          yuv1, transf.GetYUVPrecisionBits() + 1, PSNR, result));
    } else if (metric_type == SSIM_YUV) {
      WP2_CHECK_STATUS(yuv0.GetDistortion(
          yuv1, transf.GetYUVPrecisionBits() + 1, SSIM, result));
    } else {
      return WP2_STATUS_INVALID_PARAMETER;
    }
  }
  return WP2_STATUS_OK;
}

WP2Status ArgbBuffer::GetDistortionBlackOrWhiteBackground(
    const ArgbBuffer& ref, MetricType metric_type, float* result) const {
  // Since we use premultiplied values, no compositing is the same as
  // compositing on a black background.
  WP2_CHECK_STATUS(
      GetDistortion(ref, metric_type, /*include_alpha_in_all=*/false, result));
  if (!ref.HasTransparency()) {
    return WP2_STATUS_OK;
  }

  // Composite on a white background.
  // TODO(mayrla): we could composite directly inside GetDistortion() to use
  // less memory. Compositing over pure white could also be simplified.
  ArgbBuffer ref_white, this_white;
  WP2_CHECK_STATUS(ref_white.CopyFrom(ref));
  WP2_CHECK_STATUS(this_white.CopyFrom(*this));
  WP2_CHECK_STATUS(ref_white.CompositeOver(Argb32b{255, 255, 255, 255}));
  WP2_CHECK_STATUS(this_white.CompositeOver(Argb32b{255, 255, 255, 255}));
  float result_white_background[5];
  WP2_CHECK_STATUS(this_white.GetDistortion(ref_white, metric_type,
                                            /*include_alpha_in_all=*/false,
                                            result_white_background));

  // Result is the worst of the two.
  if (result_white_background[4] < result[4]) {
    // We intentionally omit copying result_white_background[0] which will be
    // 99 anyway as the image composited to a white background has no alpha.
    result[1] = result_white_background[1];
    result[2] = result_white_background[2];
    result[3] = result_white_background[3];
    result[4] = result_white_background[4];
  }
  return WP2_STATUS_OK;
}

WP2Status ArgbBuffer::GetDistortionBlackOrWhiteBackground(
    const ArgbBuffer& ref, const Rectangle& window, MetricType metric_type,
    float* result) const {
  ArgbBuffer view, ref_view;
  WP2_CHECK_STATUS(view.SetView(*this, window));
  WP2_CHECK_STATUS(ref_view.SetView(ref, window));
  return view.GetDistortionBlackOrWhiteBackground(ref_view, metric_type,
                                                  result);
}

WP2Status ArgbBuffer::GetDistortion(const ArgbBuffer& ref,
                                    uint32_t x, uint32_t y,
                                    MetricType metric_type,
                                    float* result) const {
  WP2_CHECK_OK(!IsEmpty() && !ref.IsEmpty(), WP2_STATUS_NULL_PARAMETER);
  WP2_CHECK_OK((width == ref.width) && (height == ref.height),
               WP2_STATUS_BAD_DIMENSION);
  WP2_CHECK_OK(x < width && y < height, WP2_STATUS_INVALID_PARAMETER);

  double total = 0.;
  if (metric_type == PSNR) {
    const uint8_t* const ptr1 = (const uint8_t*)GetRow(y);
    const uint8_t* const ptr2 = (const uint8_t*)ref.GetRow(y);
    for (uint32_t c = 0; c < 4; ++c) {
      const int diff = ptr1[4 * x + c] - ptr2[4 * x + c];
      total += diff * diff;
    }
    total = GetPSNR(total, 4., /*max=*/255.);
  } else if (metric_type == MSSSIM) {
    WP2SSIMInit();
    // This is extremely costly: we extract a 16x16 buffer around (x,y),
    // and call the regular MSSSIM function on it.
    const uint32_t xx = std::max(0, (int)x - 8);
    const uint32_t yy = std::max(0, (int)y - 8);
    const uint32_t w = std::min(xx + 16, width) - xx;
    const uint32_t h = std::min(yy + 16, height) - yy;
    ArgbBuffer tmp1, tmp2;
    WP2_CHECK_STATUS(tmp1.SetView(ref, {xx, yy, w, h}));
    WP2_CHECK_STATUS(tmp2.SetView(*this, {xx, yy, w, h}));
    float results[5];
    WP2_CHECK_STATUS(tmp1.GetDistortion(tmp2, MSSSIM, results));
    total = results[4];
  } else {
    WP2SSIMInit();
    // we need absolute top-left pointers
    const uint8_t* const ptr1 = (const uint8_t*)GetRow(0);
    const uint8_t* const ptr2 = (const uint8_t*)ref.GetRow(0);
    if (metric_type == SSIM) {
      for (uint32_t c = 0; c < 4; ++c) {
        WP2DistoStats stats;
        WP2SSIMGetClipped<4>(ptr1 + c, stride,
                             ptr2 + c, ref.stride, x, y, width, height, &stats);
        total += WP2SSIMCalculation(/*bit_depth=*/8, stats);
      }
      total = GetLogSSIM(total, 4.);
    } else if (metric_type == LSIM) {
      for (uint32_t c = 0; c < 4; ++c) {
        total = AccumulateLSIM(ptr1 + c, stride,
                               ptr2 + c, ref.stride, x, y, width, height);
      }
      total = GetPSNR(total, 4., /*max=*/255.);
    } else {
      return WP2_STATUS_UNSUPPORTED_FEATURE;  // TODO(skal) ?
    }
  }

  if (result != nullptr) *result = (float)total;
  return WP2_STATUS_OK;
}

//------------------------------------------------------------------------------

template <>
WP2Status Plane16::GetSSIM(const Plane16& src, uint32_t bit_depth,
                           double* const score) const {
  WP2_CHECK_OK(w_ == src.w_ && h_ == src.h_, WP2_STATUS_BAD_DIMENSION);
  WP2SSIMInit();
  *score = AccumulateSSIM<int16_t, /*channel_step=*/1>(
      w_, h_, bit_depth, Row(0), Step(), src.Row(0), src.Step());
  return WP2_STATUS_OK;
}

template <>
WP2Status Plane<uint8_t>::GetSSIM(const Plane<uint8_t>& src, uint32_t bit_depth,
                                  double* const score) const {
  WP2_CHECK_OK(w_ == src.w_ && h_ == src.h_, WP2_STATUS_BAD_DIMENSION);
  WP2SSIMInit();
  *score = AccumulateSSIM<uint8_t, /*channel_step=*/1>(
      w_, h_, bit_depth, Row(0), Step(), src.Row(0), src.Step());
  return WP2_STATUS_OK;
}

WP2Status YUVPlane::GetDistortion(const YUVPlane& ref, uint32_t bit_depth,
                                  MetricType metric_type,
                                  float result[5]) const {
  WP2_CHECK_OK(!IsEmpty(), WP2_STATUS_NULL_PARAMETER);
  WP2_CHECK_OK(!ref.IsEmpty(), WP2_STATUS_NULL_PARAMETER);
  WP2_CHECK_OK(result != nullptr, WP2_STATUS_NULL_PARAMETER);
  WP2_CHECK_OK(A.IsEmpty() == ref.A.IsEmpty(), WP2_STATUS_BAD_DIMENSION);
  WP2_CHECK_OK(GetWidth() == ref.GetWidth() && GetHeight() == ref.GetHeight() &&
                   IsDownsampled() == ref.IsDownsampled(),
               WP2_STATUS_BAD_DIMENSION);
  result[0] = result[1] = result[2] = result[3] = result[4] = 0.f;

  double total_distortion = 0.;
  const uint32_t num_pixels = GetWidth() * GetHeight();
  const uint32_t num_channels = ref.HasAlpha() ? 4 : 3;
  const double uv_multiplier = (double)num_pixels / (U.w_ * U.h_);
  if (metric_type == PSNR || metric_type == PSNR_YUV) {
    // Expected output is AYUV.
    uint64_t result_64b[4] = {0, 0, 0, 0};
    for (Channel channel : {kYChannel, kUChannel, kVChannel, kAChannel}) {
      const Plane16& p = GetChannel(channel);
      const Plane16& ref_p = ref.GetChannel(channel);
      if (channel == kAChannel && A.IsEmpty()) continue;
      WP2_CHECK_STATUS(p.GetSSE(ref_p, &result_64b[channel]));
    }

    // Hack to still have a "valid" general distortion with a bit depth and/or
    // subsampling mix: bring all squared errors to 8 bits and full sampling.
    constexpr double kMax = (1u << 8) - 1u;
    const double range_scale = std::exp2(8. - bit_depth);
    const double yuv_scale = range_scale * range_scale;  // Squared error.
    const double uv_scale = yuv_scale * uv_multiplier;
    const double result_double[] = {
        result_64b[0] * yuv_scale, result_64b[1] * uv_scale,
        result_64b[2] * uv_scale, (double)result_64b[3]};

    for (Channel channel : {kYChannel, kUChannel, kVChannel, kAChannel}) {
      if (channel != kAChannel || !A.IsEmpty()) {
        total_distortion += (float)result_double[channel];
      }
      // Expected output is AYUV, but Channel enum is YUVA.
      result[(channel + 1) % 4] =
          (float)GetPSNR(result_double[channel], num_pixels, kMax);
    }
    result[4] =
        (float)GetPSNR(total_distortion, num_pixels * num_channels, kMax);
  } else if (metric_type == SSIM || metric_type == SSIM_YUV) {
    for (Channel channel : {kYChannel, kUChannel, kVChannel, kAChannel}) {
      if (channel == kAChannel && A.IsEmpty()) {
        result[0] = (float)kNoDistortion_dB;
      } else {
        const Plane16& p = GetChannel(channel);
        const Plane16& ref_p = ref.GetChannel(channel);
        double disto;
        WP2_CHECK_STATUS(p.GetSSIM(ref_p,
            (channel == kAChannel) ? kAlphaBits : bit_depth, &disto));
        if (channel == kUChannel || channel == kVChannel) {
          disto *= uv_multiplier;  // Simulate 'num_pixels' for 420.
        }
        total_distortion += disto;
        // Expected output is AYUV, but Channel enum is YUVA.
        result[(channel + 1) % 4] = (float)GetLogSSIM(disto, num_pixels);
      }
    }
    result[4] = (float)GetLogSSIM(total_distortion, num_pixels * num_channels);
  } else {
    return WP2_STATUS_UNSUPPORTED_FEATURE;
  }
  return WP2_STATUS_OK;
}

//------------------------------------------------------------------------------

}    // namespace WP2
