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
// SuspendableDataSource implementation.
//
// Author: Yannis Guyon (yguyon@google.com)

#include "src/utils/data_source_context.h"

#if defined(WP2_USE_CONTEXT_SWITCH) && (WP2_USE_CONTEXT_SWITCH > 0)

namespace WP2 {

void SuspendableDataSource::SetContext(LocalContext* const context) {
  context_ = context;
}

bool SuspendableDataSource::HasEnoughDataToResume() const {
  return (GetNumNextBytes() >= num_requested_bytes_);
}

void SuspendableDataSource::Reset() {
  DataSource::Reset();
  context_ = nullptr;
  num_requested_bytes_ = 0;
}

bool SuspendableDataSource::Fetch(size_t num_requested_bytes) {
  num_requested_bytes_ = num_requested_bytes;
  // It should be verified before resuming but in case there is still not enough
  // data to decode, the while loop will keep yielding.
  while (!HasEnoughDataToResume()) {
    if ((context_ == nullptr) || !context_->Yield()) return false;
  }
  return true;
}

}  // namespace WP2

#endif  // WP2_USE_CONTEXT_SWITCH
