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
//   AOM residual decoding, branched from libgav1, with as little changes as
//   possible.
//
// Author: Vincent Rabaud (vrabaud@google.com)

#ifndef WP2_DEC_RESIDUALS_DEC_AOM_H_
#define WP2_DEC_RESIDUALS_DEC_AOM_H_

#include <array>

#include "src/common/lossy/aom/array_2d.h"
#include "src/common/lossy/block.h"
#include "src/common/lossy/residuals.h"
#include "src/common/lossy/residuals_aom.h"
#include "src/common/lossy/transforms.h"
#include "src/dec/symbols_dec.h"
#include "src/utils/ans.h"
#include "src/utils/utils.h"

namespace WP2 {
namespace libgav1 {

using Allocable = WP2Allocable;

//------------------------------------------------------------------------------
// From libgav1/src/tile/tile.h

// Class specialized in reading transform residuals.
class AOMResidualReader : public ResidualIO, public AOMResidualIO {
 public:
  void Init() { ResidualIO::Init(/*use_aom_coeffs=*/true); }

  // Reads residual coefficients.
  WP2Status ReadCoeffs(Channel channel, uint32_t x_pix, uint32_t y_pix,
                       WP2::BlockSize dim, bool is420, TransformPair tx_type,
                       bool first_is_dc, int max_num_coeffs, ANSDec* const dec,
                       SymbolReader* const sr, AOMContext* const aom_context,
                       int16_t* const coeffs, uint32_t* const num_coeffs);

 private:
  // Returns the number of non-zero coefficients that were read or -1 on error.
  int ReadTransformCoefficients(Plane plane, int start_x, int start_y,
                                TransformSize tx_size, TransformPair tx_type,
                                bool first_is_dc, uint32_t max_num_coeffs,
                                int16_t* const coeffs,
                                AOMContext* const aom_context);

  // This function specializes the parsing of DC coefficient by removing some of
  // the branches when i == 0 (since scan[0] is always 0 and scan[i] is always
  // non-zero for all other possible values of i). |dc_category| is an output
  // parameter that is populated when |is_dc_coefficient| is true.
  // |coefficient_level| is an output parameter which accumulates the
  // coefficient level.
  template <bool is_dc_coefficient>
  bool ReadSignAndApplyDequantization(
      int adjusted_tx_width_log2, uint16_t dc_sign_cluster,
      int16_t* const coeffs, int8_t* dc_category,
      int* coefficient_level);  // Part of 5.11.39.
  int ReadCoeffBaseRange(int clamped_tx_size_context, int cdf_context,
                         int plane_type);  // Part of 5.11.39.

  ANSDec* dec_;
  SymbolReader* sr_;
};

}  // namespace libgav1
}  // namespace WP2

#endif /* WP2_DEC_RESIDUALS_DEC_AOM_H_ */
