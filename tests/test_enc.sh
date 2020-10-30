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


# tests for encoding / decoding

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

CWP2=${BINARY_DIR}/cwp2
DWP2=${BINARY_DIR}/dwp2
GET_DISTO=${BINARY_DIR}/get_disto

TMP_WP2=/tmp/test_enc.wp2
TMP_PNG=/tmp/test_enc.png

echo "lossy test"
${CWP2} -o ${TMP_WP2} "${TESTDATA_DIR}/source1.png" -q 75 -speed 2 -summary \
        -mt 2147483647 -orientation 2 -crop 7 8 497 304 -segment_mode implicit
${DWP2} -mt 2147483647 ${TMP_WP2} -o ${TMP_PNG}
${GET_DISTO} ${TMP_PNG} "${TESTDATA_DIR}/source1.png" -lsim \
             -crop_source 7 8 497 304 -scale -ref ${TMP_WP2}

echo "EXIF + XMP test"
${CWP2} -o ${TMP_WP2} "${TESTDATA_DIR}/test_exif_xmp.webp" -psnrhvs -q 10 \
        -summary -z 4 -mt -crop 0 0 128 128 -segment_mode explicit
${DWP2} ${TMP_WP2} -o ${TMP_PNG}
${GET_DISTO} ${TMP_WP2} "${TESTDATA_DIR}/test_exif_xmp.webp" \
             -ssim -crop_source 0 0 0 0 -crop 0 0 128 128

echo "lossless get_disto test"
${CWP2} -o ${TMP_WP2} "${TESTDATA_DIR}/source3.jpg" -q 100 -short -speed 0
${GET_DISTO} ${TMP_WP2} "${TESTDATA_DIR}/source3.jpg" -exact -scale ${TMP_WP2} \
             -gray -psnrhvs

echo "TEST OK"
exit 0
