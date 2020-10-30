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
# Tests for SWIG examples.
# For this google3 version, binaries don't have .py extension and are not run
# through python3 command.

set -ex

if [ "$#" -ge 1 ]; then
  # 'eval' so that the passed in directory can contain variables.
  BINARY_DIR=$(eval echo ${1})
else
  BINARY_DIR=../swig/examples
fi
if [ "$#" -ge 2 ]; then
  GET_DISTO_DIR=$(eval echo ${2})
else
  GET_DISTO_DIR=../swig/examples
fi
if [ "$#" -ge 3 ]; then
  TESTDATA_DIR=$(eval echo ${3})
else
  TESTDATA_DIR=./testdata
fi

CWP2="python3 ${BINARY_DIR}/pycwp2.py"
DWP2="python3 ${BINARY_DIR}/pydwp2.py"
GET_DISTO="${GET_DISTO_DIR}/get_disto"
SRC_FILE0="${TESTDATA_DIR}/source0.ppm"
SRC_FILE1="${TESTDATA_DIR}/source1.png"
SRC_FILE2="${TESTDATA_DIR}/animation1.webp"
TMP_FILE="/tmp/test_swig_examples.wp2"
DST_FILE="/tmp/test_swig_examples.png"
FRAMES_FOLDER="/tmp/test_swig_examples/"
mkdir -p ${FRAMES_FOLDER}

${CWP2} --version
${DWP2} --version

# Image
${CWP2} ${SRC_FILE0} -o ${TMP_FILE} -s 0
${DWP2} ${TMP_FILE} -o ${DST_FILE}
${GET_DISTO} ${SRC_FILE0} ${DST_FILE}

# Animation
${CWP2} --frames ${SRC_FILE0} 1000 ${SRC_FILE1} 1000 -o ${TMP_FILE} -s 0
${DWP2} ${TMP_FILE} -o ${DST_FILE} -f ${FRAMES_FOLDER}
${GET_DISTO} ${SRC_FILE0} "${FRAMES_FOLDER}/frame0_1000ms.png"
${GET_DISTO} ${SRC_FILE1} "${FRAMES_FOLDER}/frame1_1000ms.png"

${CWP2} ${SRC_FILE2} -o ${TMP_FILE} -s 0
${DWP2} ${TMP_FILE} -o ${DST_FILE} -f ${FRAMES_FOLDER}
