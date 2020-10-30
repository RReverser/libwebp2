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


# tests for encoding / decoding extra formats

set -ex

if [ "$#" -ge 1 ]; then
  # eval so that the passed in directory can contain variables.
  BINARY_DIR=$(eval echo ${1})
else
  # assume "tests" is the current directory
  BINARY_DIR=..
fi
if [ "$#" -ge 2 ]; then
  TESTDATA_DIR=$(eval echo ${2})
else
  TESTDATA_DIR=./testdata
fi

EWP2=${BINARY_DIR}/ewp2
TMP_Y4M=/tmp/test_enc_ewp2.y4m
TMP_WP2=/tmp/test_enc_ewp2.wp2

# Encode/decode wp2 from/to 4:4:4 with 12 bits
${EWP2} "${TESTDATA_DIR}/ccsp/source3_64x48_C444p12.y4m" \
        -o "${TMP_WP2}" -q 80 -speed 5
${EWP2} "${TMP_WP2}" -o "${TMP_Y4M}" -depth 12 \
        -psnr "${TESTDATA_DIR}/ccsp/source3_64x48_C444p12.y4m" 34.7

# Encode/decode wp2 from/to 4:2:0 with 8 bits
${EWP2} "${TESTDATA_DIR}/ccsp/source3_64x48_C420p8.y4m" \
        -o "${TMP_WP2}" -q 100 -speed 4 -sampling nearest
${EWP2} "${TMP_WP2}" -o "${TMP_Y4M}" -depth 8 -420 -sampling nearest \
        -psnr "${TESTDATA_DIR}/ccsp/source3_64x48_C420p8.y4m" 55.0
${EWP2} "${TMP_WP2}" -o "${TMP_Y4M}" -depth 8 -420 -sampling nearest \
        -ssim "${TESTDATA_DIR}/ccsp/source3_64x48_C420p8.y4m" 35.0

echo "TEST OK"
exit 0
