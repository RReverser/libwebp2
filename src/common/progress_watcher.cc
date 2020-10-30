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
// ProgressWatcher implementation.
//
// Author: Yannis Guyon (yguyon@google.com)

#include "src/common/progress_watcher.h"

#include <cmath>

#include "src/common/constants.h"
#include "src/utils/utils.h"

namespace WP2 {

WP2Status ProgressWatcher::Set(float progress) {
  WP2_CHECK_STATUS(thread_lock_.Acquire());
  if (!aborted_) {
    assert(progress >= 0.f && progress <= 1.f);
    progress_ = progress;
    if (hook_ != nullptr) aborted_ = !hook_->OnUpdate(progress_);
  }
  thread_lock_.Release();
  return aborted_ ? WP2_STATUS_USER_ABORT : WP2_STATUS_OK;
}

WP2Status ProgressWatcher::Rewind() { return Set(0.f); }

WP2Status ProgressWatcher::Finish() {
  assert(progress_ < 1.f);
  assert(std::abs(progress_ + kProgressDecEnd - 1.f) < 0.001f);
  return Set(1.f);
}

WP2Status ProgressWatcher::AdvanceBy(float step) {
  WP2_CHECK_STATUS(thread_lock_.Acquire());
  if (!aborted_) {
    progress_ += scale_ * step;
    assert(progress_ >= 0.f && progress_ <= 1.f);
    if (hook_ != nullptr) aborted_ = !hook_->OnUpdate(progress_);
  }
  const bool aborted = aborted_;  // Prevent threaded data race.
  thread_lock_.Release();
  return aborted ? WP2_STATUS_USER_ABORT : WP2_STATUS_OK;
}

}  // namespace WP2
