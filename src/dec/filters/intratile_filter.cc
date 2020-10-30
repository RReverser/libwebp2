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
// Intratile filter.
//
// Author: Yannis Guyon (yguyon@google.com)

#include "src/dec/filters/intratile_filter.h"

#include "src/utils/utils.h"
#include "src/wp2/decode.h"

namespace WP2 {

//------------------------------------------------------------------------------

bool IsIntratileFilterEnabled(const DecoderConfig& config,
                              const BitstreamFeatures& features,
                              const GlobalParams& gparams) {
  return (IsDeblockingFilterEnabled(config, features) ||
          IsDirectionalFilterEnabled(config, features) ||
          IsRestorationFilterEnabled(config) ||
          IsGrainFilterEnabled(config, gparams.segments_) ||
          IsAlphaFilterEnabled(config, gparams));
}

// The number of pixel rows above a block that may be used for prediction. These
// pixels should not be filtered before that.
static constexpr uint32_t kMaxNumPredRows = 1u;

//------------------------------------------------------------------------------

IntratileFilter::IntratileFilter(const DecoderConfig& config,
                                 const BitstreamFeatures& features,
                                 const GlobalParams& gparams,
                                 const FilterBlockMap& blocks,
                                 const RstrFltParams& params)
    : config_(config),
      blocks_(blocks),
      enabled_(IsIntratileFilterEnabled(config, features, gparams)),
      deblocking_filter_(config, features, blocks),
      directional_filter_(config, features, blocks),
      restoration_filter_(config, blocks, params),
      alpha_filter_(config, gparams, blocks),
      grain_filter_(config, gparams.segments_, blocks) {
  (void)config_;
}

WP2Status IntratileFilter::Allocate() {
  if (!enabled_) return WP2_STATUS_OK;

  WP2_CHECK_STATUS(deblocking_filter_.Allocate());
  WP2_CHECK_STATUS(directional_filter_.Allocate());
  WP2_CHECK_STATUS(restoration_filter_.Allocate());
  WP2_CHECK_STATUS(alpha_filter_.Allocate());
  WP2_CHECK_STATUS(grain_filter_.Allocate());
  return WP2_STATUS_OK;
}

//------------------------------------------------------------------------------

uint32_t IntratileFilter::Apply(uint32_t num_decoded_yuv_rows) {
  if (!enabled_) return num_decoded_yuv_rows;

  uint32_t num_filtered_rows = num_decoded_yuv_rows;
  if (num_filtered_rows < blocks_.tile_rect_.height) {
    // Do not deblock pixels that might still be used for prediction.
    num_filtered_rows = SafeSub(num_filtered_rows, kMaxNumPredRows);
  }
  num_filtered_rows = deblocking_filter_.Deblock(num_filtered_rows);
  num_filtered_rows = directional_filter_.Smooth(num_filtered_rows);
  num_filtered_rows = restoration_filter_.Enhances(num_filtered_rows);
  num_filtered_rows = alpha_filter_.Smooth(num_filtered_rows);

  // TODO(skal): The GrainFilter happens now (YUV) but the IntertileFilter might
  //             blur it later (RGB). Either remove it or find a solution.
  num_filtered_rows = grain_filter_.Apply(num_filtered_rows);

  return num_filtered_rows;
}

//------------------------------------------------------------------------------

}  // namespace WP2
