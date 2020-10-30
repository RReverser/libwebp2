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
TESTDATA_DIR=${TESTDATA_DIR}/specific

CWP2=${BINARY_DIR}/cwp2
GET_DISTO=${BINARY_DIR}/get_disto

TMP_WP2=/tmp/test_io_pnm.wp2

# check exotic PNM files
do_check_psnr() {
  local file=$1
  shift
  local min_psnr=$1
  shift
  ${GET_DISTO} -min "${min_psnr}" "${TMP_WP2}" "${TESTDATA_DIR}/${file}" "$@"
}

for v in 37 255 1001 65535; do
  min_psnr=41  # low quality only for max_value = 37
  if [ "x${v}" != "x37" ]; then
    min_psnr=99
  fi
  opt="-q 100 -quiet"

  ${CWP2} "${TESTDATA_DIR}/bug449_${v}.pgm" ${opt} -o "${TMP_WP2}"
  do_check_psnr bug449_ref_gray.png ${min_psnr}

  ${CWP2} "${TESTDATA_DIR}/bug449_${v}_g.pam" ${opt} -o "${TMP_WP2}"
  do_check_psnr bug449_ref_gray_alpha.png ${min_psnr}

  ${CWP2} "${TESTDATA_DIR}/bug449_${v}.ppm" ${opt} -o "${TMP_WP2}"
  do_check_psnr bug449_ref_noalpha.png ${min_psnr}

  ${CWP2} "${TESTDATA_DIR}/bug449_${v}.pam" ${opt} -o "${TMP_WP2}"
  do_check_psnr bug449_ref.png ${min_psnr}
done

echo "ALL OK"
exit 0
