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
//  AOM / AV1 wrapper
//
// Author: Skal (pascal.massimino@gmail.com)

#ifndef WP2_EXTRAS_AOM_UTILS_H_
#define WP2_EXTRAS_AOM_UTILS_H_

#include <string>
#include <vector>

#include "src/wp2/base.h"

namespace WP2 {

struct ParamsAV1 {
  float quality = 75.f;         // in [0..100]
  float alpha_quality = 100.f;  // in [0..100]
  uint32_t size = 0;            // target size if not zero (predates 'quality')
  uint32_t speed = 5;           // in [0..9]
  uint32_t threads = 1;
  uint32_t pass = 2;
  bool use_yuv444 = true;
  float filter_strength = 50.f;  // in [0..100]
  MetricType tuning = SSIM;      // either PSNR or SSIM
  uint32_t bittrace = 0;         // 0=off, 1=bits, 2=bytes
  bool compact_trace = false;    // compact traces into fewer categories
  bool draw_blocks = false;      // for visual debug
  bool draw_transforms = false;
  const char* draw_mi_data = nullptr;
  bool draw_mi_number = false;  // draw value of mi data as red square (9 bits)
};

struct QuantAV1 {
  // Dequantization values per segment, one for DC, one for AC.
  int16_t y_dequant[8][2];
  int16_t u_dequant[8][2];
  int16_t v_dequant[8][2];
};

WP2Status CompressAV1(const ArgbBuffer& ref, const ParamsAV1& params,
                      ArgbBuffer* const decoded, std::string* const out,
                      double timing[2] = nullptr,
                      std::vector<WP2::Rectangle>* const blocks = nullptr,
                      std::vector<WP2::Rectangle>* const transforms = nullptr,
                      QuantAV1* const quant = nullptr);

WP2Status CompressAVIF(const ArgbBuffer& ref, const ParamsAV1& params,
                       ArgbBuffer* const decoded, std::string* const out,
                       double timing[2] = nullptr);

}  // namespace WP2

#endif  // WP2_EXTRAS_AOM_UTILS_H_
