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


# tests for image formats I/O

set -ex

quiet="yes"  # change to 'no' if needed

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

TMP_IMG=/tmp/test_io
TMP_WP2="${TMP_IMG}.wp2"

formats_to_test="png ppm tiff pam bmp"
sources1="source0.ppm alpha_ramp.tiff source4.webp"

do_cmd() {
  cmd=$*
  if [ "x${quiet}" = "xno" ]; then
    echo "running: ${cmd}"
    ${cmd}
  else
    ${cmd} > /dev/null
  fi
}

check_roundtrip() {
  echo "=== $1 ==="
  min_psnr=$2
  min_msssim=$3
  do_cmd ${CWP2} "${1}" $4 -o "${TMP_WP2}" -psnr -short -speed 0
  for fmt in ${formats_to_test}; do
    do_cmd ${DWP2} -quiet -v "${TMP_WP2}" -${fmt} -o "${TMP_IMG}.${fmt}"
    if [ "x${fmt}" != "xppm" ]; then
      ${GET_DISTO} "${TMP_IMG}.${fmt}" "${1}" -min ${min_psnr}
    else  # ppm doesn't handle alpha, so psnr is non-sense
      ${GET_DISTO} "${TMP_IMG}.${fmt}" "${1}" -min 2
    fi
    ${GET_DISTO} "${TMP_WP2}" "${1}" -min ${min_psnr}
    # test more exotic metric and options :
    ${GET_DISTO} "${TMP_WP2}" "${1}" -msssim -half -min ${min_msssim}
  done
}

# Minimal quality checks
for file in ${sources1}; do
  check_roundtrip "${TESTDATA_DIR}/${file}" 33.1 19.25
done

check_roundtrip "${TESTDATA_DIR}/source3_222x167.jpg" 28.8 18.8 "-pass 3"

check_roundtrip "${TESTDATA_DIR}/source1_64x48.png" 29.3 19.8 "-alpha_q 80"

${DWP2} -i "${TMP_WP2}"

# MSAN check on specific crash vectors found by the fuzzer.
${CWP2} ${TESTDATA_DIR}/specific/testcase-6544524677218304 || echo "OK"
${CWP2} ${TESTDATA_DIR}/specific/testcase-6544524677218304 -q 100 || echo "OK"

echo "ALL OK"
exit 0
