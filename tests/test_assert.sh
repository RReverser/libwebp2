#!/bin/sh
# Copyright 2020 Google LLC
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#      http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

#
# tests that asserts trigger properly

set -ex

if [ "$#" -ge 1 ]; then
  # eval so that the passed in directory can contain variables.
  BINARY_DIR=$(eval echo ${1})
else
  # assume "tests" is the current directory
  BINARY_DIR=.
fi

ASSERT_BIN=${BINARY_DIR}/test_assert_bin

if [ ! -f "${ASSERT_BIN}" ]; then
  echo "${ASSERT_BIN} not found."
  exit 1
fi

# Check if we have assert support.
if ${ASSERT_BIN} 0 >/dev/null 2>/dev/null; then
  echo "No assert support."
  exit 0
fi

for i in $(seq -4 -1; seq 1 8)
do
  printf "Case %s ... " "${i}"
  # Check that the binary does not crash.
  if ${ASSERT_BIN} "${i}" >/dev/null 2>/dev/null; then
    echo "PASSED."
    # If it does not crash, i has to be < 0.
    if [ "${i}" -gt 0 ]; then
      exit 1
    fi
  else
    # If it does crash, i has to be > 0.
    echo "CRASHED."
    if [ "${i}" -lt 0 ]; then
      exit 1
    fi
  fi
done

echo "TEST OK"
exit 0
