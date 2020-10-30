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

# Pull and run Docker image
WP2_ROOT=$(pwd)/git/wp2
docker run -v "${WP2_ROOT}":/wp2 \
 gcr.io/google.com/webm-project/wp2_docker_image \
 bash -c "cd /wp2 && ./kokoro/ubuntu/build.sh"
