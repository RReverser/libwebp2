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
// Wiener-optimal filtering
//
// Author: Skal (pascal.massimino@gmail.com)

#include "src/utils/wiener.h"

#include <cmath>
#include <utility>

namespace WP2 {

const float kEpsilon = 1e-10f;  // tolerance for non-invertibility

bool LinearSolve(size_t n, double A[], double y[], float x[]) {
  // Partial pivoting (idempotent).
  for (size_t i = n - 1; i > 0; i--) {
    if (A[(i - 1) * n] < A[i * n]) {
      for (size_t j = 0; j < n; ++j) {
        std::swap(A[i * n + j], A[(i - 1) * n + j]);
      }
      std::swap(y[i], y[i - 1]);
    }
  }
  // Forward elimination
  for (size_t k = 0; k < n - 1; k++) {
    if (fabs(A[k * n + k]) < kEpsilon) {
      return false;
    }
    const double norm = 1. / A[k * n + k];
    for (size_t i = k + 1; i < n; ++i) {
      const double c = norm * A[i * n + k];
      for (size_t j = 0; j < n; ++j) {
        A[i * n + j] -= c * A[k * n + j];
      }
      y[i] -= c * y[k];
    }
  }
  // Backward substitution
  for (size_t i = n; i-- > 0; ) {
    if (fabs(A[i * n + i]) < kEpsilon) return false;
    double c = 0.;
    for (size_t j = i + 1; j < n; ++j) {
      c += A[i * n + j] * x[j];
    }
    x[i] = (y[i] - c) / A[i * n + i];
  }
  return true;
}

}    // namespace WP2
