#!/bin/bash
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

# Fail on any error.
set -e
# Display commands being run.
set -x

SECONDS=0

rm -rf build
mkdir build
cd build

# CMake release
cmake -D CMAKE_BUILD_TYPE=RelWithDebInfo -D WP2_ENABLE_TESTS=ON \
-D WP2_BUILD_EXAMPLES=ON -D WP2_BUILD_EXTRAS=ON -D WP2_BUILD_SWIG_PY=ON ..
make -j 8
ctest -j 8

rm -rf ./*

# CMake debug
cmake -D CMAKE_BUILD_TYPE=Debug -D WP2_ENABLE_TESTS=ON \
-D WP2_BUILD_EXAMPLES=ON -D WP2_BUILD_EXTRAS=ON ..
make -j 8
# ctest -j 8  # Tests take too much time to run in Debug.

rm -rf ./*

duration=$SECONDS
echo "Build and tests took ${duration} seconds."
