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
TMP_Y4M=/tmp/test_io_ewp2.y4m
TMP_PNG=/tmp/test_io_ewp2.png
TMP_WEBP=/tmp/test_io_ewp2.webp

# Convert y4m to y4m losslessly
${EWP2} "${TESTDATA_DIR}/ccsp/source3_C444p8.y4m" \
        -o "${TMP_Y4M}" -depth 8 -444 \
        -psnr "${TESTDATA_DIR}/ccsp/source3_C444p8.y4m" 99.0

${EWP2} "${TESTDATA_DIR}/ccsp/source3_C420p12.y4m" \
        -o "${TMP_Y4M}" -depth 12 -420 -sampling nearest \
        -psnr "${TESTDATA_DIR}/ccsp/source3_C420p12.y4m" 99.0

# No output
${EWP2} "${TMP_Y4M}" -psnr "${TMP_Y4M}" 99.0 -size
${EWP2} "${TMP_Y4M}" -psnr "${TMP_Y4M}" 99.0 -v
${EWP2} "${TMP_Y4M}" -psnr "${TMP_Y4M}" 99.0 -size -v -rm A

# Some loss during up/downsampling
${EWP2} "${TESTDATA_DIR}/ccsp/source3_C420p12.y4m" \
        -o "${TMP_Y4M}" -depth 12 -420 -sampling smooth \
        -psnr "${TESTDATA_DIR}/ccsp/source3_C420p12.y4m" 50.0

# Encode/decode WebP from JPG
${EWP2} "${TESTDATA_DIR}/source3.jpg" -o "${TMP_WEBP}" -q 100 -speed 2
${EWP2} "${TMP_WEBP}" -o "${TMP_PNG}" -ssim "${TESTDATA_DIR}/source3.jpg" 99.0

# Encode/decode WebP from/to PNG (some premultiplied-alpha loss)
${EWP2} "${TESTDATA_DIR}/source1.png" -o "${TMP_WEBP}" -q 100 -speed 2
${EWP2} "${TMP_WEBP}" -o "${TMP_PNG}" -ssim "${TESTDATA_DIR}/source1.png" 38.0
${EWP2} "${TESTDATA_DIR}/source1.png" -o "${TMP_WEBP}" -q 0 -speed 0
${EWP2} "${TMP_WEBP}" -o "${TMP_PNG}" -ssim "${TESTDATA_DIR}/source1.png" 8.0

echo "TEST OK"
exit 0
