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
// Coding tools configuration
//
// Author: Skal (pascal.massimino@gmail.com)

#include "src/common/constants.h"
#include "src/wp2/encode.h"

#ifdef HAVE_CONFIG_H
#include "src/wp2/config.h"
#endif

namespace WP2 {

//------------------------------------------------------------------------------

bool EncoderConfig::IsValid() const {
  bool not_ok = false;
  // Tested in declaration order.
  not_ok |= (quality < 0 || quality > kMaxQuality);
  not_ok |= (target_psnr < 0);
  not_ok |= (alpha_quality < 0 || alpha_quality > (int)kMaxQuality);
  not_ok |= (speed < 0 || speed > 9);
  not_ok |= (preview_size > 0 && preview == nullptr);
  not_ok |= (transfer_function > (int)WP2_TF_ARIB_STD_B67_HLG);
  not_ok |= (pass < 1 || pass > 10);
  not_ok |= (sns < 0 || sns > 100);
  not_ok |= (error_diffusion < 0 || error_diffusion > 100);
  not_ok |= (segments < 1 || segments > (int)kMaxNumSegments);
  not_ok |= (segment_id_mode < 0 || segment_id_mode > SEGMENT_ID_IMPLICIT);
  for (const float s : segment_factors) {
    not_ok |= (s < 0.f || s > kMaxLossyQuality);
  }
  not_ok |= (tile_shape < 0 || tile_shape > TILE_SHAPE_AUTO);
  not_ok |= (partition_method < 0 || partition_method >= NUM_PARTITION_METHODS);
  not_ok |= (partition_set < 0 || partition_set >= NUM_PARTITION_SETS);
  not_ok |= ((int)uv_mode >= NumUVMode || (int)uv_mode < UVModeAdapt);
  not_ok |= (preprocessing < 0 || preprocessing > 7);
  not_ok |= (preprocessing_strength < 0 || preprocessing_strength > 100);
  not_ok |= (thread_level < 0 || thread_level > 2147483647);
  not_ok |= (use_neural_compression != 0 && use_neural_compression != 1);
  return !not_ok;
}

//------------------------------------------------------------------------------

}  // namespace WP2
