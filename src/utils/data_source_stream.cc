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
// StreamDataSource implementation.
//
// Author: Yannis Guyon (yguyon@google.com)

#include "src/utils/data_source_stream.h"

#include <algorithm>

#include "src/dsp/math.h"
#include "src/utils/utils.h"

namespace WP2 {

WP2Status StreamDataSource::AppendAsExternal(const uint8_t* bytes,
                                             size_t size) {
  if (num_new_bytes_to_discard_ > 0) {
    const size_t num_discarded_new_bytes =
        std::min(size, num_new_bytes_to_discard_);
    bytes += num_discarded_new_bytes;
    size -= num_discarded_new_bytes;
    num_new_bytes_to_discard_ -= num_discarded_new_bytes;
  }
  if (size > 0) {
    assert(bytes != nullptr);
    WP2_CHECK_STATUS(AppendExternalToInternal());  // Make room if needed.
    external_bytes_ = bytes;
    num_external_bytes_ = size;

    // Only one contiguous array is returned by TryGetNext(), so multiple
    // buffers are merged. Thus StreamDataSource contains either external or
    // internal data but never both at the same time.
    if (!internal_data_.empty()) {
      WP2_CHECK_STATUS(AppendExternalToInternal());
    } else {
      available_bytes_ = external_bytes_;
      num_available_bytes_ = num_external_bytes_;
    }
  }
  return WP2_STATUS_OK;
}

WP2Status StreamDataSource::AppendExternalToInternal() {
  if (num_external_bytes_ > 0) {
    assert(external_bytes_ != nullptr);
    const uint8_t* const external_bytes = external_bytes_;
    const size_t num_external_bytes = num_external_bytes_;
    // Remove external reference even if error. AppendExternalToInternal() might
    // get called again and external_bytes_ might not be valid anymore.
    external_bytes_ = available_bytes_ = nullptr;
    num_external_bytes_ = num_available_bytes_ = 0;

    const size_t internal_size = internal_data_.size();
    WP2_CHECK_ALLOC_OK(
        internal_data_.resize(internal_size + num_external_bytes));
    memcpy(internal_data_.data() + internal_size, external_bytes,
           num_external_bytes);
    available_bytes_ = internal_data_.data();
    num_available_bytes_ = internal_data_.size();
  }
  return WP2_STATUS_OK;
}

void StreamDataSource::Reset() {
  DataSource::Reset();
  external_bytes_ = nullptr;
  num_external_bytes_ = 0;
  internal_data_.clear();
  num_new_bytes_to_discard_ = 0;
}

bool StreamDataSource::Fetch(size_t num_requested_bytes) { return false; }

void StreamDataSource::OnDiscard(size_t num_bytes) {
  if (num_external_bytes_ > 0) {
    if (num_bytes >= num_external_bytes_) {
      num_new_bytes_to_discard_ =
          SafeAdd(num_new_bytes_to_discard_, num_bytes - num_external_bytes_);
      external_bytes_ = nullptr;
      num_external_bytes_ = 0;
    } else {
      assert(external_bytes_ != nullptr);
      external_bytes_ += num_bytes;
      num_external_bytes_ -= num_bytes;
    }
    available_bytes_ = external_bytes_;
    num_available_bytes_ = num_external_bytes_;
  } else {
    if (num_bytes >= internal_data_.size()) {
      num_new_bytes_to_discard_ =
          SafeAdd(num_new_bytes_to_discard_, num_bytes - internal_data_.size());
      internal_data_.clear();
    } else {
      const size_t new_size = internal_data_.size() - num_bytes;
      memmove(internal_data_.data(), internal_data_.data() + num_bytes,
              new_size);
      const bool success = internal_data_.resize(new_size);
      (void)success;
      assert(success);
    }
    available_bytes_ = internal_data_.data();
    num_available_bytes_ = internal_data_.size();
  }
}

}  // namespace WP2
