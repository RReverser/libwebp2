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
#
# tests for command lines

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
VWP2=${BINARY_DIR}/vwp2
SRC_FILE1="${TESTDATA_DIR}/source1.png"
SRC_FILE2="${TESTDATA_DIR}/source3.jpg"
SRC_FILE3="${TESTDATA_DIR}/source1_64x48.png"
TMP_FILE=/tmp/cmd_test.wp2
TMP_FILE2=/tmp/cmd_test2.png
TMP_FILE3=/tmp/cmd_test3.png

# simple coverage of command line arguments

${CWP2}
${CWP2} -version
${CWP2} -h
${CWP2} ${SRC_FILE3}
cat ${SRC_FILE3} | ${CWP2} -- -

${DWP2} -version
${DWP2} -h

${VWP2}
${VWP2} -version
${VWP2} -h
${VWP2} -q 30 -nomt -noasm
${VWP2} -quants 30,64,12,64,13,432,12,125,32 -nomt -noasm
${VWP2} -info -d -alt /tmp -ps 3 -speed 3 -segments 1 -pass 2 \
        -nometadata -bt 2 -q 50 -set_block 1 2 6 9 -h

${GET_DISTO} -version
${GET_DISTO} -h

# test -d option
${CWP2} -o ${TMP_FILE} ${SRC_FILE3} -speed 1 -d ${TMP_FILE2} -bt -create_preview
${DWP2} ${TMP_FILE} -o ${TMP_FILE3} -histos -count_blocks -short
${GET_DISTO} ${TMP_FILE2} ${TMP_FILE3} -exact

# 'preview file'
echo "asdfvkladfsvkal" > ${TMP_FILE}
${CWP2} ${TMP_FILE3} -preview ${TMP_FILE}

# lossy vdebug test
${CWP2} -o ${TMP_FILE} ${SRC_FILE1} -crop 30 31 1 20 -speed 0 -diffusion 50 \
        -quants "60,43, 23, 23,54 ,23" -BT -histos -tile_shape 3
for vdebug in `seq -1 33`; do
  ${DWP2} ${TMP_FILE} -mt -vdebug ${vdebug} &> /dev/null
done

# lossless vdebug test
${CWP2} -o ${TMP_FILE} ${SRC_FILE1} -crop 30 31 1 20 -speed 0 -q 99 -sns 30 -v
for vdebug in `seq -1 20`; do
  ${DWP2} ${TMP_FILE} -vdebug ${vdebug} &> /dev/null
done

# grain test
${CWP2} -o ${TMP_FILE} ${SRC_FILE1} -crop 30 31 200 200 -speed 0 \
        -grain -noasm -segment_mode auto -perceptual
${GET_DISTO} ${TMP_FILE} ${SRC_FILE1} -scale -crop_source 30 31 200 200 \
             -prefix "this is a test" -suffix "good"

for grain in 0 37 78 100; do
  ${DWP2} ${TMP_FILE} -grain ${grain} &> /dev/null
done
${DWP2} ${TMP_FILE} -grain 10 -noasm &> /dev/null

# metrics identity test
for metric in psnr ssim msssim lsim psnrhvs; do
  ${GET_DISTO} ${SRC_FILE2} ${SRC_FILE2} -exact -${metric}
  ${GET_DISTO} ${SRC_FILE3} ${SRC_FILE3} -exact -${metric}
  ${CWP2} -${metric} ${SRC_FILE1} -crop 0 0 1 1 -target_size 2000
done

# misc flag test for 'get_disto'
${GET_DISTO} ${SRC_FILE3} ${SRC_FILE3} -exact -gray -scale -o ${TMP_FILE}.webp
${GET_DISTO}
${GET_DISTO} ${SRC_FILE1}

# negative tests
set +e
${GET_DISTO} ${SRC_FILE1} ${SRC_FILE1} -o
${GET_DISTO} ${SRC_FILE1} ${SRC_FILE1} -crop 4 2 1
${GET_DISTO} ${SRC_FILE1} ${SRC_FILE1} -bad_option

${CWP2} -tile_shape 10
${CWP2} -csp
${CWP2} -uv_mode
${CWP2} -pm
${CWP2} -ps
${CWP2} -transfer
${CWP2} -target_size
${CWP2} -target_psnr
${CWP2} -loop
${CWP2} -segments
${CWP2} -pre
${CWP2} -hint
${CWP2} -preds
${CWP2} -bad_option
${CWP2} -preview /tmp
${CWP2} -preview
${CWP2} ${TMP_FILE2} ${TMP_FILE3} -preview ${TMP_FILE}  # multiple input
${CWP2} -
${CWP2} -grain    # missing input file name

${DWP2}
${DWP2} -bad_option
${DWP2} -v

echo "TEST OK"
exit 0
