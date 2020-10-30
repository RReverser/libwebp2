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
//  Preview config
//
// Author: Yannis Guyon (yguyon)

#include "src/enc/preview/preview_enc.h"

#include <cassert>
#include <cmath>

namespace WP2 {

//------------------------------------------------------------------------------

static float Lerp(uint32_t min_value, uint32_t max_value, float value) {
  return min_value * (1.f - value) + max_value * value;
}
static uint32_t ILerp(uint32_t min_value, uint32_t max_value, float value) {
  return (uint32_t)std::lround(Lerp(min_value, max_value, value));
}

PreviewConfig::PreviewConfig(float quality, int speed) {
  const float q = quality / 100.0f;                  // [0,1]
  const float s = speed / 9.0f;                      // [0,1]
  const float qs = ((1 + q) * (1 + s) - 1.f) / 3.f;  // [0,1]
  assert(q >= 0.f && q <= 1.f && s >= 0.f && s <= 1.f);

  num_colors = ILerp(6, 10, std::sqrt(qs));
  grid_density = Lerp(1.0, 3.0, s);
  num_iterations = ILerp(500u, 3000u, s);

  assert(IsValid());
}

//------------------------------------------------------------------------------

bool PreviewConfig::IsValid() const {
  if (num_vertices < kPreviewMinNumVertices) return false;
  if (num_vertices > kPreviewMaxNumVertices) return false;
  if (num_colors < kPreviewMinNumColors) return false;
  if (num_colors > kPreviewMaxNumColors) return false;
  if (grid_size < kPreviewMinGridSize) return false;
  if (grid_size > kPreviewMaxGridSize) return false;

  if (grid_density < kPreviewMinGridDensity) return false;
  if (grid_density > kPreviewMaxGridDensity) return false;
  if (blur_radius > kPreviewMaxBlurRadius) return false;
  if (blur_threshold < kPreviewMinBlurThreshold) return false;
  if (blur_threshold > kPreviewMaxBlurThreshold) return false;
  if (num_iterations > kPreviewMaxNumIterations) return false;

  // TODO(yguyon): Verify all settings
  return true;
}

//------------------------------------------------------------------------------

}  // namespace WP2
