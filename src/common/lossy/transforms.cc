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
//  Transforms
//
#include "src/common/lossy/transforms.h"

namespace WP2 {

const WP2TransformType kTfX[]{
    kDct, kAdst, kDct,  kAdst, kIdentity, kDct,      kIdentity};
const WP2TransformType kTfY[]{
    kDct, kAdst, kAdst, kDct,  kDct,      kIdentity, kIdentity};
STATIC_ASSERT_ARRAY_SIZE(kTfX, kNumTransformPairs);
STATIC_ASSERT_ARRAY_SIZE(kTfY, kNumTransformPairs);

TransformClass GetTransformClass(TransformPair transform) {
  if (kTfX[transform] == kIdentity && kTfY[transform] == kIdentity) {
    return TransformClass::kTwoD;
  } else if (kTfX[transform] == kIdentity) {
    return TransformClass::kVertical;
  } else if (kTfY[transform] == kIdentity) {
    return TransformClass::kHorizontal;
  } else {
    return TransformClass::kTwoD;
  }
}

}  // namespace WP2
