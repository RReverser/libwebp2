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

# tests for encoding / decoding directories

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
TMP_TESTDATA_DIR=/tmp/test_enc_dir_testdata
DST_DIR=/tmp/test_enc_dir
rm -fr "${TMP_TESTDATA_DIR}"
mkdir -p ${TMP_TESTDATA_DIR}
rm -fr "${DST_DIR}"
mkdir -p ${DST_DIR}

# Copy some files so that it doesn't take ages compressing the whole folder.
cp -f "${TESTDATA_DIR}/source0.pgm" "${TMP_TESTDATA_DIR}/"
cp -f "${TESTDATA_DIR}/source2.tiff" "${TMP_TESTDATA_DIR}/"
cp -f "${TESTDATA_DIR}/source1_1x1.png" "${TMP_TESTDATA_DIR}/"
cp -f "${TESTDATA_DIR}/source1_64x48.png" "${TMP_TESTDATA_DIR}/"

echo "input folder test"
${CWP2} ${TMP_TESTDATA_DIR} -r -o ${DST_DIR} -force -speed 0
${DWP2} ${DST_DIR} -r -o ${DST_DIR} -force

echo "multithreaded decoding"
${DWP2} ${DST_DIR} -r -o ${DST_DIR} -force -mt
${DWP2} ${DST_DIR} -r -o ${DST_DIR} -force -mt 1
${DWP2} ${DST_DIR} -r -o ${DST_DIR} -mt 4 -force

echo "several input files test"
${CWP2} "${TESTDATA_DIR}/source1.png" "${TESTDATA_DIR}/source3.jpg" \
        -o ${DST_DIR} -force -speed 0
${DWP2} "${DST_DIR}/source1.wp2" "${DST_DIR}/source3.wp2" \
        -o ${DST_DIR} -pam -force

echo "inplace test"
${CWP2} "${DST_DIR}/source1.pam" "${DST_DIR}/source3.pam" \
        -inplace -force -speed 0
${DWP2} "${DST_DIR}/source1.wp2" "${DST_DIR}/source3.wp2" \
        -inplace -force

set +e  # Negative tests.

echo "Negative tests"
${DWP2} "${DST_DIR}/source1.png" "${DST_DIR}/source3.png" && die
${CWP2} "${TESTDATA_DIR}/source1.png" "${TESTDATA_DIR}/source3.jpg" \
        -o "${DST_DIR}/not_a_directory.wp2" && die
${DWP2} "${DST_DIR}/source1.wp2" "${DST_DIR}/source3.wp2" \
        -o "${DST_DIR}/not_a_directory.png" && die
${DWP2} "${DST_DIR}/source1.wp2" -mt 4294967295 && die
${DWP2} "${DST_DIR}/source3.wp2" -mt -1 && die

echo "test -d option with multiple files"
${CWP2} "${TESTDATA_DIR}/source1_64x48.png" "${TESTDATA_DIR}/source1.png" \
        -o "${TMP_TESTDATA_DIR}" -d "${DST_DIR}/tmp.webp" && die
${CWP2} "${TESTDATA_DIR}/source1_64x48.png"  -d "${DST_DIR}/" && die

echo "test good/bad files mix"
printf "%b" '\xf4\xff\x6f\xff\xff\xff\xff\xff\xff\xff\xff\xff' \
    > "${DST_DIR}/bad.wp2"
${DWP2} "${DST_DIR}/source1.wp2" "${DST_DIR}/bad.wp2" \
        -o ${DST_DIR} -pam -force -mt && die
${DWP2} "${DST_DIR}/bad.wp2" "${DST_DIR}/source3.wp2" \
        -o ${DST_DIR} -pam -mt -force && die
${DWP2} ${DST_DIR} -r -o ${DST_DIR} -force -mt 1 && die

# Clean up
rm -fr "${TMP_TESTDATA_DIR}"

echo "TEST OK"
exit 0
