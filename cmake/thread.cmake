#
# Copyright 2020 Google LLC.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#      https://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#

# Check for thread libraries.
include(CheckCSourceCompiles)

find_package(Threads)

if(Threads_FOUND)
  if(CMAKE_USE_PTHREADS_INIT)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -pthread")
  endif()
  foreach(PTHREAD_TEST HAVE_PTHREAD_PRIO_INHERIT PTHREAD_CREATE_UNDETACHED)
    check_c_source_compiles(
      "
        #include <pthread.h>
        int main (void) {
          int attr = ${PTHREAD_TEST};
          return attr;
        }
      "
      ${PTHREAD_TEST})
  endforeach()
  list(APPEND WP2_DEP_LIBRARIES ${CMAKE_THREAD_LIBS_INIT})
endif()

set(WP2_USE_THREAD ${Threads_FOUND})
