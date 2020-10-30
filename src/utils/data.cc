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
// Data and Metadata
//
// Author: Skal (pascal.massimino@gmail.com)

#include <cstring>  // for memcpy()

#include "src/utils/utils.h"
#include "src/wp2/base.h"

namespace WP2 {

//------------------------------------------------------------------------------

Data::Data(Data&& other) noexcept { TrivialMoveCtor(this, &other); }

void Data::Clear() {
  WP2Free(bytes);
  bytes = nullptr;
  size = 0;
}

WP2Status Data::Resize(size_t data_size, bool keep_bytes) {
  if (size == data_size) return WP2_STATUS_OK;
  if (data_size == 0) {
    Clear();
  } else if (keep_bytes) {
    uint8_t* const new_bytes = (uint8_t*)WP2Realloc(bytes, 1, data_size);
    WP2_CHECK_ALLOC_OK(new_bytes != nullptr);
    bytes = new_bytes;
    size = data_size;
  } else {
    Clear();
    bytes = (uint8_t*)WP2Malloc(data_size, sizeof(*bytes));
    WP2_CHECK_ALLOC_OK(bytes != nullptr);
    size = data_size;
  }
  return WP2_STATUS_OK;
}

WP2Status Data::CopyFrom(const uint8_t* data, size_t data_size) {
  if (data_size == 0) {
    Clear();
  } else {
    WP2_CHECK_OK(data != nullptr, WP2_STATUS_INVALID_PARAMETER);
    WP2_CHECK_STATUS(Resize(data_size, /*keep_bytes=*/false));
    std::memcpy(bytes, data, size);
  }
  return WP2_STATUS_OK;
}

WP2Status Data::Append(const uint8_t* data, size_t data_size) {
  if (data_size > 0) {
    WP2_CHECK_OK(data != nullptr, WP2_STATUS_INVALID_PARAMETER);
    WP2_CHECK_STATUS(Resize(size + data_size, /*keep_bytes=*/true));
    std::memcpy(bytes + size - data_size, data, data_size);
  }
  return WP2_STATUS_OK;
}

void swap(Data& a, Data& b) {
  std::swap(a.bytes, b.bytes);
  std::swap(a.size, b.size);
}

//------------------------------------------------------------------------------

bool Metadata::IsEmpty() const {
  if (iccp.size > 0 || exif.size > 0 || xmp.size > 0) {
    return false;
  }
  return true;
}

void Metadata::Clear() {
  exif.Clear();
  iccp.Clear();
  xmp.Clear();
}

WP2Status Metadata::CopyFrom(const Metadata& src) {
  WP2_CHECK_STATUS(exif.CopyFrom(src.exif.bytes, src.exif.size));
  WP2_CHECK_STATUS(iccp.CopyFrom(src.iccp.bytes, src.iccp.size));
  WP2_CHECK_STATUS(xmp.CopyFrom(src.xmp.bytes, src.xmp.size));
  return WP2_STATUS_OK;
}

void swap(Metadata& a, Metadata& b) {
  swap(a.exif, b.exif);
  swap(a.iccp, b.iccp);
  swap(a.xmp, b.xmp);
}

//------------------------------------------------------------------------------

}  // namespace WP2
