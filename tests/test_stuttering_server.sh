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
# Test for stuttering_http_server.py.

set -ex

if [ "$#" -ge 1 ]; then
  # 'eval' so that the passed in directory can contain variables.
  BINARY_DIR=$(eval echo ${1})
else
  BINARY_DIR=./tools
fi
if [ "$#" -ge 2 ]; then
  TESTDATA_DIR=$(eval echo ${2})
else
  TESTDATA_DIR=./testdata
fi

SRC_FOLDER="${TESTDATA_DIR}"
SRC_FILE="source1_64x48.png"
DST_FOLDER="/tmp/test_stuttering_server/"
mkdir -p ${DST_FOLDER}

STUTTERING_SERVER="python3 ${BINARY_DIR}/stuttering_http_server.py"
PORT=8080

# Option '--test' makes the binary exit after the first transfer.
${STUTTERING_SERVER} --ip "localhost" --port ${PORT} --folder "${SRC_FOLDER}" \
                     --test & SERVER_PID="$!"
# Make sure the server is killed when this test exits.
trap "kill -0 ${SERVER_PID} && kill -s KILL ${SERVER_PID}" \
     0 HUP INT QUIT ABRT KILL TERM
sleep 10  # Wait for server to come alive.
wget -O "${DST_FOLDER}/${SRC_FILE}" \
     "http://localhost:${PORT}/${SRC_FILE}?bytes=512&ms=100"
# Verify dest file is the same as src file.
cmp "${SRC_FOLDER}/${SRC_FILE}" "${DST_FOLDER}/${SRC_FILE}"

# Wait for the http server to exit by itself or timeout and fail the test.
tail --pid=${SERVER_PID} -f /dev/null
