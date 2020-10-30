// Copyright 2020 Google LLC
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//      http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.
// -----------------------------------------------------------------------------
//
// Test library version compared to header.

#include "extras/extras.h"
#include "include/helpers.h"
#include "src/wp2/base.h"

namespace WP2 {
namespace {

TEST(VersionTest, Base) {
  ASSERT_TRUE(WP2CheckVersion());
}

TEST(VersionTest, Extras) {
  ASSERT_FALSE(
      WP2_ABI_IS_INCOMPATIBLE(WP2GetExtrasVersion(), WP2_EXTRAS_ABI_VERSION));
}

}  // namespace
}  // namespace WP2
