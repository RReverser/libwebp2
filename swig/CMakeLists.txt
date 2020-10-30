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

find_package(SWIG REQUIRED)
include(UseSWIG)

find_package(PythonLibs 3 REQUIRED)

set_property(SOURCE libwebp2.i PROPERTY CPLUSPLUS ON)
swig_add_library(
  libwebp2
  TYPE
  SHARED
  LANGUAGE
  python
  OUTPUT_DIR
  ${CMAKE_CURRENT_BINARY_DIR}
  OUTFILE_DIR
  ${CMAKE_CURRENT_BINARY_DIR}
  SOURCES
  ${CMAKE_CURRENT_SOURCE_DIR}/libwebp2.i)
target_include_directories(${SWIG_MODULE_libwebp2_REAL_NAME}
                           PRIVATE ${PROJECT_SOURCE_DIR} ${PYTHON_INCLUDE_DIRS})
target_link_libraries(${SWIG_MODULE_libwebp2_REAL_NAME}
                      PRIVATE imageio webp2 ${PYTHON_LIBRARIES})

# There are warnings relative to 'swig_this' and 'SWIG_init'. Removing them. See
# https://github.com/swig/swig/issues/1278
get_property(
  FLAGS
  TARGET ${SWIG_MODULE_libwebp2_REAL_NAME}
  PROPERTY COMPILE_FLAGS)
set(FLAGS "${FLAGS} -Wno-shadow -Wno-missing-declarations")
set_property(TARGET ${SWIG_MODULE_libwebp2_REAL_NAME} PROPERTY COMPILE_FLAGS
                                                             ${FLAGS})
