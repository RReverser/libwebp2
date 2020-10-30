// Copyright 2018 Google LLC
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
// Binary that tests some asserts with ANS
// Author: vrabaud@google.com (Vincent Rabaud)

#include <cassert>

#include "src/utils/ans.h"
#include "src/utils/ans_utils.h"
#include "src/utils/vector.h"

int main(int argc, const char* argv[]) {
  if (argc != 2) return 1;
  const int c = atoi(argv[1]);    // NOLINT (atoi)
  WP2::ANSEnc enc;
  // All cases <0 should execute normally.
  // All cases >0 should fail because of an assert.
  // The 0 case assert(0). It is used to check if the executable is built with
  // assert support or not.
  switch (c) {
    case -4:
      enc.PutUValue(0, 0, "zero");
      break;
    case -3: {
      // Compiles because ctors are implicitly-defined.
      struct Foo {
        Foo() = default;
        // Foo(const Foo&) = delete;  // Does not compile.
        // Foo(Foo&&) = delete;       // Does not compile.
      };
      WP2::Vector<Foo> a;
      if (!a.resize(1)) assert(false);         // Uses default and move ctors.
      if (!a.push_back(Foo())) assert(false);  // Uses copy ctor.
      break;
    }
    case -2:
      // Proba should be in [0,PROBA_MAX].
      enc.PutBit(0, PROBA_MAX, "bit");
      enc.PutBit(1, 0, "bit");
      break;
    case -1:
      // Range has to fit in [1..(1 << WP2::kANSMaxRangeBits)]
      enc.PutRValue(1, (1 << WP2::kANSMaxRangeBits), "rvalue");
      break;
    case 0:
      assert(0);
      break;
    case 1:
      // Range has to fit in [1..(1 << WP2::kANSMaxRangeBits)]
      enc.PutRValue(1, (1 << WP2::kANSMaxRangeBits) + 1, "rvalue");
      break;
    case 2:
      // Value has to be < range.
      enc.PutRValue(10, 1, "rvalue");
      break;
    case 3:
      // Value has to fit in U bits.
      enc.PutUValue(10, 2, "uvalue");
      break;
    case 4:
      // Value has to fit in IO_BITS.
      enc.PutUValue(1 << IO_BITS, IO_BITS + 1, "uvalue");
      break;
    case 5:
      // Dict has to be valid.
      enc.PutASymbol(0, (WP2::ANSAdaptiveSymbol*)nullptr, "id");
      break;
    case 6:
      // Proba should fit in [0,PROBA_MAX].
      enc.PutBit(1, PROBA_MAX + 1, "bit");
      break;
    case 7:
      enc.PutUValue(1, 0, "zero");
      break;
    case 8:
      enc.PutUValue(1, WP2::kANSMaxUniformBits + 1, "zero");
      break;
    default:
      return 1;
  }
  // Test for asserts in the Emit functions.
  const WP2Status status = enc.Assemble();
  (void)status;
  assert(status == WP2_STATUS_OK);

  return 0;
}
