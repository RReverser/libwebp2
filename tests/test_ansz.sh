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
# tests for ansz

set -ex

if [ "$#" -ge 1 ]; then
  # eval so that the passed in directory can contain variables.
  BINARY_DIR=$(eval echo ${1})
else
  # assume "tests" is the current directory
  BINARY_DIR=.
fi
if [ "$#" -ge 2 ]; then
  TESTDATA_DIR=$(eval echo ${2})
else
  TESTDATA_DIR=./testdata
fi
if [ "$#" -ge 3 ]; then
  ANSZCC_DIR=$(eval echo ${3})
else
  ANSZCC_DIR=.
fi

ANSZ=${BINARY_DIR}/ansz
TMP_FILE1=/tmp/ansz_test.ansz
TMP_FILE2=/tmp/ansz_test2

for infile in ${ANSZCC_DIR}/ansz.cc ${TESTDATA_DIR}/*.*; do
  ${ANSZ} "${infile}" -o ${TMP_FILE1} -check
  ${ANSZ} -d ${TMP_FILE1} -o ${TMP_FILE2}
  diff ${TMP_FILE2} "${infile}"
  ${ANSZ} "${infile}" -o ${TMP_FILE1} -check -o0
  ${ANSZ} -d ${TMP_FILE1} -o ${TMP_FILE2}
  diff ${TMP_FILE2} "${infile}"
  echo "${infile} OK"
done

rm -f ${TMP_FILE1} ${TMP_FILE2}

echo "TEST OK"
exit 0
