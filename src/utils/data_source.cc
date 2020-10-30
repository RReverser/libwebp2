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
// DataSource implementation.
//
// Author: Yannis Guyon (yguyon@google.com)

#include "src/utils/data_source.h"

#include <algorithm>
#include <cassert>

namespace WP2 {

//------------------------------------------------------------------------------

bool DataSource::TryGetNext(size_t num_bytes, const uint8_t** data_pointer) {
  if (num_bytes > GetNumNextBytes()) {
    if (!Fetch(num_bytes)) return false;
    // Fetch must return false if there isn't enough data.
    assert(num_bytes <= GetNumNextBytes());
  }
  assert(data_pointer != nullptr);
  *data_pointer = available_bytes_ + num_read_bytes_;
  return true;
}

size_t DataSource::GetNumNextBytes() const {
  return num_available_bytes_ - std::min(num_available_bytes_, num_read_bytes_);
}

size_t DataSource::GetNumReadBytes() const {
  return num_read_bytes_;
}

size_t DataSource::GetNumDiscardedBytes() const {
  return num_discarded_bytes_;
}

void DataSource::MarkNumBytesAsRead(size_t num_bytes) {
  num_read_bytes_ += num_bytes;
}

void DataSource::UnmarkAllReadBytes() { num_read_bytes_ = 0; }

void DataSource::Discard(size_t num_bytes) {
  if (num_bytes >= num_available_bytes_) {
    available_bytes_ = nullptr;
    num_available_bytes_ = 0;
  } else {
    available_bytes_ += num_bytes;
    num_available_bytes_ -= num_bytes;
  }
  num_read_bytes_ -= std::min(num_bytes, num_read_bytes_);
  OnDiscard(num_bytes);
  num_discarded_bytes_ += num_bytes;
}

//------------------------------------------------------------------------------

const uint8_t* DataSource::DataHandle::GetBytes() const {
  if (data_source_ != nullptr && data_source_->available_bytes_ != nullptr &&
      offset_ >= data_source_->num_discarded_bytes_) {
    const size_t offset_in_available_bytes =
        offset_ - data_source_->num_discarded_bytes_;
    if (data_source_->num_available_bytes_ >=
        offset_in_available_bytes + size_) {
      return data_source_->available_bytes_ + offset_in_available_bytes;
    }
  }
  return nullptr;
}

bool DataSource::TryGetNext(size_t num_bytes, DataHandle* const data_handle) {
  const uint8_t* data_pointer;
  if (TryGetNext(num_bytes, &data_pointer)) {
    const size_t offset =
        (data_pointer - available_bytes_) + num_discarded_bytes_;
    *data_handle = DataHandle(this, offset, num_bytes);
    return true;
  }
  return false;
}

bool DataSource::TryReadNext(size_t num_bytes, DataHandle* const data_handle) {
  if (TryGetNext(num_bytes, data_handle)) {
    MarkNumBytesAsRead(num_bytes);
    return true;
  }
  return false;
}

void DataSource::Reset() {
  num_read_bytes_ = 0;
  num_discarded_bytes_ = 0;
  available_bytes_ = nullptr;
  num_available_bytes_ = 0;
}

//------------------------------------------------------------------------------

void ExternalDataSource::Update(const uint8_t* data, size_t data_size) {
  if (num_discarded_bytes_ >= data_size) {
    available_bytes_ = nullptr;
    num_available_bytes_ = 0;
  } else {
    available_bytes_ = data + num_discarded_bytes_;
    num_available_bytes_ = data_size - num_discarded_bytes_;
  }
}

// With this simple version of DataSource there is no way of getting more data
// through Fetch(). There is nothing to discard OnDisposal() either.
bool ExternalDataSource::Fetch(size_t) { return false; }
void ExternalDataSource::OnDiscard(size_t) {}

//------------------------------------------------------------------------------

}  // namespace WP2
