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
// Colorspace utilities
//
// Author: Skal (pascal.massimino@gmail.com)

#ifndef WP2_UTILS_CSP_H_
#define WP2_UTILS_CSP_H_

#include <array>
#include <cstdint>

#include "src/wp2/base.h"
#include "src/wp2/format_constants.h"

namespace WP2 {

//------------------------------------------------------------------------------

class ANSEnc;
class ANSDec;

// Contains a matrix and its fixed-point precision (the number of bits to
// right-shift the samples with after multiplication).
struct CSPMtx {
  CSPMtx() = default;
  CSPMtx(const CSPMtx&) = default;
  constexpr CSPMtx(const std::array<int16_t, 9>& m, uint32_t s)
      : matrix(m), shift(s) {}
  constexpr CSPMtx(const int16_t m[9], uint32_t s)
      : matrix{m[0], m[1], m[2], m[3], m[4], m[5], m[6], m[7], m[8]},
        shift(s) {}

  const int16_t* mtx() const { return matrix.data(); }
  bool operator==(const CSPMtx& other) const {
    return matrix == other.matrix && shift == other.shift;
  }
  bool operator!=(const CSPMtx& other) const { return !operator==(other); }
  CSPMtx& operator=(const CSPMtx&) = default;

  std::array<int16_t, 9> matrix;
  uint32_t shift;
};

//------------------------------------------------------------------------------

// Custom YUV <-> RGB transform matrix
// YUV samples are internally represented using yuv_bits_ (including sign).
// The maximum precision is kMaxYuvBits (+ 1 bit for the sign).
// The inverse matrix (YUV to RGB) is im_[].
// The direct matrix m_[] (RGB to YUV) is always derived from im_[], and not
// computed 'as is' nor stored. It uplifts from kRgbBits input to yuv_bits_
// internal representation. m_ and im_ elements use kMtxBits in memory and
// are stored as fixed point with a precision of kMtxShift.

class CSPTransform {
 public:
  // kMtxBits is high for accuracy (> kMtxShift) but low enough to keep elements
  // within 16 bits so that the 16b x 16b = 32b mult is easy to optimize.
  static constexpr uint32_t kMtxBits = 15;   // Sign not included.
  static constexpr uint32_t kMtxShift = 12;  // Applied after matrix mul
                                             // (fixed-point precision).

  CSPTransform() { InitYCoCg(); }
  Csp GetType() const { return type_; }

  // YUV <-> RGB conversion with several type overloads.
  void ToYUV(int16_t r, int16_t g, int16_t b,
             int16_t* const y, int16_t* const u, int16_t* const v) const;
  void ToRGB(int16_t y, int16_t u, int16_t v,
             int16_t* const r, int16_t* const g, int16_t* const b) const;

  void ToYUV(uint8_t r, uint8_t g, uint8_t b,
             int16_t* const y, int16_t* const u, int16_t* const v) const;
  void ToRGB(int16_t y, int16_t u, int16_t v,
             uint8_t* const r, uint8_t* const g, uint8_t* const b) const;

  void ToYUV(const int16_t rgb[3], int16_t yuv[3]) const;
  void ToRGB(const int16_t yuv[3], int16_t rgb[3]) const;

  void ToYUV(const uint8_t rgb[3], int16_t yuv[3]) const;
  void ToRGB(const int16_t yuv[3], uint8_t rgb[3]) const;

  Ayuv38b ToYUV(const Argb32b& color) const;
  Argb32b ToRGB(const Ayuv38b& color) const;

  // Converts the channels 'c0, c1, c2' from a custom color space to RGB using
  // 'ccsp_to_rgb_matrix'. Copies alpha if not null and clamp RGB to [0:alpha].
  static WP2Status CustomToArgb(uint32_t width, uint32_t height,
                                const int16_t* c0_buffer, uint32_t c0_step,
                                const int16_t* c1_buffer, uint32_t c1_step,
                                const int16_t* c2_buffer, uint32_t c2_step,
                                const int16_t* a_buffer, uint32_t a_step,
                                const CSPMtx& ccsp_to_rgb,
                                ArgbBuffer* const argb);

  // Converts the channels 'c0, c1, c2' from a custom color space to
  // 'CSPTransform::type_' color space using 'ccsp_to_rgb_matrix'.
  // Input 'c0/1/2_buffer' can be the same as output 'y/u/v_buffer' (inplace).
  WP2Status CustomToYUV(uint32_t width, uint32_t height,
                        const int16_t* c0_buffer, uint32_t c0_step,
                        const int16_t* c1_buffer, uint32_t c1_step,
                        const int16_t* c2_buffer, uint32_t c2_step,
                        const CSPMtx& ccsp_to_rgb, int16_t* y_buffer,
                        uint32_t y_step, int16_t* u_buffer, uint32_t u_step,
                        int16_t* v_buffer, uint32_t v_step) const;

  // Converts the channels 'y, u, v' from the 'CSPTransform::type_' color space
  // to a custom color space defined by 'rgb_to_ccsp_matrix'.
  // Input 'y/u/v_buffer' can be the same as output 'c0/1/2_buffer' (inplace).
  WP2Status YUVToCustom(uint32_t width, uint32_t height,
                        const int16_t* y_buffer, uint32_t y_step,
                        const int16_t* u_buffer, uint32_t u_step,
                        const int16_t* v_buffer, uint32_t v_step,
                        const CSPMtx& rgb_to_ccsp, int16_t* c0_buffer,
                        uint32_t c0_step, int16_t* c1_buffer, uint32_t c1_step,
                        int16_t* c2_buffer, uint32_t c2_step) const;

  // Returns false if the yuv_to_rgb_matrix[] is not invertible.
  WP2_NO_DISCARD bool Init(const int16_t yuv_to_rgb_matrix[9],
                           const int16_t rgb_avg[3] = kDefaultRgbAvg);

  void InitYCoCg();  // lossless
  void InitYCbCr();  // regular YCbCr
  void InitYIQ();    // YIQ
  // Scans input 'image' to deduce optimal matrix.
  // Returns WP2_STATUS_INVALID_COLORSPACE if it failed to find a correct one.
  WP2Status Optimize(const ArgbBuffer& image);

  // Calls one of the above depending on 'csp_type'.
  // 'argb_buffer' cannot be empty if it is kCustom.
  WP2Status Init(Csp csp_type, const ArgbBuffer& argb_buffer = ArgbBuffer());

  // Encodes 'm' then decodes it as 'dm'. In case the former is not perfectly
  // shaped as the expected eigenmatrix, 'dm' might have suffered some loss
  // during encoding. The maximum 'error' between the two is returned.
  static WP2Status MakeEigenMatrixEncodable(const int16_t m[9], int16_t dm[9],
                                            int32_t* const error);

  // Writes this instance to the bitstream.
  WP2Status Write(ANSEnc* const enc) const;
  WP2Status Read(ANSDec* const dec);

  // Returns the number of bits used to store YUV samples, without the sign.
  uint32_t GetYUVPrecisionBits() const { return yuv_bits_ - 1; }
  // Returns min/max YUV values that this transform produces.
  int16_t GetYUVMin() const { return min_yuv_value_; }
  int16_t GetYUVMax() const { return max_yuv_value_; }
  uint32_t GetYUVRange() const { return 1 + max_yuv_value_ - min_yuv_value_; }

  // debug
  void Print() const;
  bool CheckRoundTrip(const int16_t yuv[3], const int16_t rgb[3]) const;
  const int16_t* GetRgbToYuvMatrix() const { return m_; }
  const int16_t* GetYuvToRgbMatrix() const { return im_; }
  const int16_t* GetRgbAverage() const { return avg_; }
  // In-place apply the transform to 'image' area defined by {x,y,w,h}.
  // Output is yuva packed as 8b-unsigned (not s16).
  // 'num_precision_bits' includes the sign.
  void Apply(ArgbBuffer* const image, uint32_t num_precision_bits,
             uint32_t x, uint32_t y, uint32_t w, uint32_t h) const;

 private:
  // Input RGB and output RGB are in range [0..255].
  static constexpr int16_t kDefaultRgbAvg[3] = {
      (1u << (kRgbBits - 1)), (1u << (kRgbBits - 1)), (1u << (kRgbBits - 1))};

  // Sets yuv_bits_ and min_yuv_value_/max_yuv_value_.
  void SetYUVBounds(int16_t yuv_min, int16_t yuv_max);
  bool ComputeRGBToYUVMatrix();  // Computes m_[] from im_[].
  // deduce best YUV bounds from m_ and avg_
  void StoreCustomYUVBounds();

  Csp type_;

  uint32_t yuv_bits_;  // Depends on type. Includes sign.
  int16_t min_yuv_value_;
  int16_t max_yuv_value_;

  int16_t m_[9];    // RGB to YUV matrix, scaled by 1<<kMtxShift.
  int16_t im_[9];   // YUV to RGB matrix, scaled by 1<<kMtxShift.
  int16_t avg_[3];  // Input RGB average. Used as translation.
};

//------------------------------------------------------------------------------

}   // namespace WP2

#endif   // WP2_UTILS_CSP_H_
