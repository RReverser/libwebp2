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
//   AOM residual encoding, reverse-engineered from libgav1 decoder.
//
// Author: Vincent Rabaud (vrabaud@google.com)

#ifndef WP2_ENC_RESIDUALS_ENC_AOM_H_
#define WP2_ENC_RESIDUALS_ENC_AOM_H_

#include "src/common/lossy/block.h"
#include "src/common/lossy/residuals_aom.h"
#include "src/common/lossy/transforms.h"
#include "src/enc/symbols_enc.h"
#include "src/utils/ans.h"

namespace WP2 {
namespace libgav1 {

// Class writing AOM residuals.
// Each member function matches one in AOMResidualReader.
class AOMResidualWriter : public ResidualIO, public AOMResidualIO {
 public:
  void Init() { ResidualIO::Init(/*use_aom_coeffs=*/true); }

  // Writes coefficients of block of size 'dim' at 'x_pix, y_pix' .
  // 'tx_type' is written only if 'implicit_tx' is 'kUnknownTf'.
  // If 'do_update', the 'aom_context' is modified.
  static WP2Status WriteCoeffs(uint32_t x_pix, uint32_t y_pix,
                               WP2::BlockSize dim, TransformPair tx_type,
                               bool first_is_dc, bool is420,
                               const int16_t* const coeffs, Channel channel,
                               bool do_update, SymbolManager* const sm,
                               ANSEncBase* const enc,
                               AOMContext* const aom_context);

 private:
  static void WriteTransformCoefficients(
      const int16_t* const res, const AOMContext& aom_context, bool is420,
      Plane plane, int start_x, int start_y, TransformSize tx_size,
      TransformPair tx_type, bool first_is_dc, SymbolManager* const sm,
      ANSEncBase* const enc, int* const coefficient_level,
      int8_t* const dc_category);
  template <bool is_dc_coefficient>
  static void WriteSignAndApplyDequantization(
      int16_t res, uint32_t dc_sign_cluster, int* coefficient_level,
      SymbolManager* const sm, ANSEncBase* const enc);
  static uint32_t WriteCoeffBaseRange(uint32_t level,
                                      int clamped_tx_size_context,
                                      int cdf_context, int plane_type,
                                      SymbolManager* const sm,
                                      ANSEncBase* const enc);
};

}  // namespace libgav1
}  // namespace WP2

#endif /* WP2_ENC_RESIDUALS_ENC_AOM_H_ */
