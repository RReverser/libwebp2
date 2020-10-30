// Copyright 2019 Google LLC
// Copyright 2019 The libgav1 Authors
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

#ifndef WP2_COMMON_LOSSY_AOM_ARRAY_2D_H_
#define WP2_COMMON_LOSSY_AOM_ARRAY_2D_H_

#include <cstddef>
#include <cstring>
#include <memory>
#include <new>
#include <type_traits>

#include "src/wp2/base.h"

namespace WP2 {
namespace libgav1 {

// Exposes a 1D allocated memory buffer as a 2D array.
template <typename T>
class Array2DView {
 public:
  Array2DView() = default;
  Array2DView(int rows, int columns, T* const data) {
    Reset(rows, columns, data);
  }

  // Copyable and Movable.
  Array2DView(const Array2DView& rhs) = default;
  Array2DView& operator=(const Array2DView& rhs) = default;

  void Reset(int rows, int columns, T* const data) {
    rows_ = rows;
    columns_ = columns;
    data_ = data;
  }

  int rows() const { return rows_; }
  int columns() const { return columns_; }

  T* operator[](int row) { return const_cast<T*>(GetRow(row)); }

  const T* operator[](int row) const { return GetRow(row); }

 private:
  const T* GetRow(int row) const {
    const ptrdiff_t offset = static_cast<ptrdiff_t>(row) * columns_;
    return (row < rows_) ? data_ + offset : nullptr;
  }

  int rows_ = 0;
  int columns_ = 0;
  T* data_ = nullptr;
};

// Allocates and owns the contiguous memory and exposes an Array2DView of
// dimension |rows| x |columns|.
template <typename T>
class Array2D {
 public:
  Array2D() = default;

  WP2_NO_DISCARD bool CopyFrom(const Array2D<T>& other) {
    if (!Reset(other.rows(), other.columns(), /*zero_initialize=*/false)) {
      return false;
    }
    size_ = other.size_;
    if (std::is_trivially_copyable<T>::value &&
        std::is_trivially_destructible<T>::value) {
      std::memcpy((void*)data_.get(), (const void*)other.data_.get(),
                  size_ * sizeof(T));
    } else {
      std::copy(other.data_.get(), other.data_.get() + other.size_,
                data_.get());
    }
    return true;
  }

  WP2_NO_DISCARD bool Reset(int rows, int columns,
                            bool zero_initialize = true) {
    const size_t size = rows * columns;
    // If T is not a trivial type, we should always reallocate the data_
    // buffer, so that the destructors of any existing objects are invoked.
    if (!std::is_trivial<T>::value || size_ < size) {
      // Note: This invokes the global operator new if T is a non-class type,
      // such as integer or enum types, or a class type that is not derived
      // from libgav1::Allocable, such as std::unique_ptr. If we enforce a
      // maximum allocation size or keep track of our own heap memory
      // consumption, we will need to handle the allocations here that use the
      // global operator new.
      if (zero_initialize) {
        data_.reset(new (std::nothrow) T[size]());
      } else {
        data_.reset(new (std::nothrow) T[size]);
      }
      if (data_ == nullptr) {
        size_ = 0;
        return false;
      }
      size_ = size;
    } else if (zero_initialize) {
      // Cast the data_ pointer to void* to avoid the GCC -Wclass-memaccess
      // warning. The memset is safe because T is a trivial type.
      void* dest = data_.get();
      memset(dest, 0, sizeof(T) * size);
    }
    data_view_.Reset(rows, columns, data_.get());
    return true;
  }

  int rows() const { return data_view_.rows(); }
  int columns() const { return data_view_.columns(); }
  T* data() { return data_.get(); }

  T* operator[](int row) { return data_view_[row]; }

  const T* operator[](int row) const { return data_view_[row]; }

 private:
  std::unique_ptr<T[]> data_ = nullptr;
  size_t size_ = 0;
  Array2DView<T> data_view_;
};

}  // namespace libgav1
}  // namespace WP2

#endif  // WP2_COMMON_LOSSY_AOM_ARRAY_2D_H_
