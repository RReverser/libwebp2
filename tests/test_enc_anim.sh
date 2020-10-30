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
# ------------------------------------------------------------------------------

# tests for encoding / decoding animations

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
DST_DIR=/tmp/test_enc_anim
rm -fr "${DST_DIR}"
mkdir -p "${DST_DIR}"

# Input images are expected to be of equal sizes.
echo "-f test"
${CWP2} \
  -f "${TESTDATA_DIR}/source0.pgm" 50 \
     "${TESTDATA_DIR}/source1.png" 75 \
     "${TESTDATA_DIR}/source2.tiff" 100 \
  -o "${DST_DIR}/anim.wp2" -speed 0

echo "-frames test"
${DWP2} \
  "${DST_DIR}/anim.wp2" \
  -o "${DST_DIR}/" -frames -force -v -bt \
  -no_dblk_filter -no_drct_filter -no_rstr_filter
${GET_DISTO} "${TESTDATA_DIR}/source0.pgm" \
  "${DST_DIR}/frame00000_50ms.png" -min 32.0
${GET_DISTO} "${TESTDATA_DIR}/source1.png" \
  "${DST_DIR}/frame00001_75ms.png" -min 32.0
${GET_DISTO} "${TESTDATA_DIR}/source2.tiff" \
  "${DST_DIR}/frame00002_100ms.png" -min 32.0

echo "parse file names for frame order and duration"
${CWP2} "${DST_DIR}/" -r -frames -o "${DST_DIR}/anim.wp2" -speed 0

echo "from animated wp2 to animated wp2"
${CWP2} "${DST_DIR}/anim.wp2" -o "${DST_DIR}/anim_encoded_twice.wp2" -speed 0

echo "one-frame animation"
${CWP2} -f "${TESTDATA_DIR}/source3_222x167.jpg" 50 -speed 0

echo "decode still image as an animation"
${CWP2} "${TESTDATA_DIR}/source3_222x167.jpg" -o "${DST_DIR}/still.wp2" -speed 0
${DWP2} "${DST_DIR}/still.wp2" -o "${DST_DIR}/" -frames -force -BT -v \
        -dblk_filter -drct_filter -rstr_filter

echo "test file paths that begin with a dash"
cp -f "${TESTDATA_DIR}/source0.pgm" "${DST_DIR}/-source0.pgm"
cp -f "${TESTDATA_DIR}/source1.png" "${DST_DIR}/-source1.png"
cd "${DST_DIR}/"
${CWP2} -f -- -source0.pgm 500 -- -source1.png 750 -o -anim.wp2 -speed 0

# Clean up
rm -fr "${DST_DIR}"

echo "TEST OK"
exit 0
