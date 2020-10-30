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

add_test(
  NAME test_io
  COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/test_io.sh ${CMAKE_BINARY_DIR}
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
add_test(
  NAME test_io_pnm
  COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/test_io_pnm.sh ${CMAKE_BINARY_DIR}
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
add_test(
  NAME test_cmd
  COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/test_cmd.sh ${CMAKE_BINARY_DIR}
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
add_test(
  NAME test_enc
  COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/test_enc.sh ${CMAKE_BINARY_DIR}
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
add_test(
  NAME test_enc_anim
  COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/test_enc_anim.sh ${CMAKE_BINARY_DIR}
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
add_test(
  NAME test_enc_dir
  COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/test_enc_dir.sh ${CMAKE_BINARY_DIR}
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})

if(WP2_BUILD_EXTRAS)
  add_test(
    NAME test_io_ewp2
    COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/test_io_ewp2.sh ${CMAKE_BINARY_DIR}
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
  add_test(
    NAME test_enc_ewp2
    COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/test_enc_ewp2.sh ${CMAKE_BINARY_DIR}
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
endif()

add_executable(ansz ansz.cc)
target_link_libraries(ansz imageio webp2)
add_test(
  NAME test_ansz
  COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/test_ansz.sh ${CMAKE_CURRENT_BINARY_DIR}
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})

add_executable(test_assert_bin test_assert.cc)
target_link_libraries(test_assert_bin webp2 imageio)
add_test(NAME test_assert COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/test_assert.sh
                                  ${CMAKE_CURRENT_BINARY_DIR})

# Test with no Google Test dependency.
add_executable(test_simple_enc_dec test_simple_enc_dec.cc)
target_link_libraries(test_simple_enc_dec imageio webp2)
add_test(NAME test_simple_enc_dec
         COMMAND ${CMAKE_CURRENT_BINARY_DIR}/test_simple_enc_dec
         WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})

# ##############################################################################
# Dependencies

# Google Test default directories.
set(GTEST_SOURCE_DIR
    "${PROJECT_SOURCE_DIR}/../googletest"
    CACHE STRING "Path to Google Test source directory.")
set(GTEST_BUILD_DIR
    "${CMAKE_BINARY_DIR}/googletest"
    CACHE STRING
          "Path to directory where Google Test will be configured and built.")

# Try the Ubuntu path if we do not have a googletest directory.
if(NOT EXISTS "${GTEST_SOURCE_DIR}/CMakeLists.txt")
  if(EXISTS "/usr/src/gtest/CMakeLists.txt")
    set(GTEST_SOURCE_DIR "/usr/src/googletest")
  endif()
endif()

# Confirm Google Test can be built.
include(GoogleTest)
if(EXISTS "${GTEST_SOURCE_DIR}/CMakeLists.txt")
  message(STATUS "Using ${GTEST_SOURCE_DIR}/CMakeLists.txt.")
  add_subdirectory(${GTEST_SOURCE_DIR} ${GTEST_BUILD_DIR} EXCLUDE_FROM_ALL)
  # There will be a compilation warning: gtest.cc:2973:13: warning: format
  # string is not a string literal So we remove "-Werror" for that target. Also
  # a "no previous declaration" warning is shown if compiling from libgtest-dev:
  # we use -Wno-missing-declarations.
  get_property(
    FLAGS
    TARGET gmock
    PROPERTY COMPILE_FLAGS)
  set(FLAGS "${FLAGS} -Wno-format-nonliteral -Wno-missing-declarations")
  set_property(TARGET gtest PROPERTY COMPILE_FLAGS ${FLAGS})
  set_property(TARGET gmock PROPERTY COMPILE_FLAGS ${FLAGS})
  set_property(TARGET gmock_main PROPERTY COMPILE_FLAGS ${FLAGS})
  set(GTEST_FOUND TRUE)
  set(GTEST_INCLUDE_DIRS "${GTEST_SOURCE_DIR}/googletest/include")
  set(GMOCK_INCLUDE_DIRS "${GTEST_SOURCE_DIR}/googlemock/include")
  set(GTEST_BOTH_LIBRARIES gtest gtest_main)
  foreach(TARGET_NAME gtest gtest_main)
    get_property(
      tmp
      TARGET ${TARGET_NAME}
      PROPERTY INTERFACE_SYSTEM_INCLUDE_DIRECTORIES)
    set_property(TARGET ${TARGET_NAME} PROPERTY INTERFACE_INCLUDE_DIRECTORIES
                                                ${tmp})
    # Reset INTERFACE_SYSTEM_INCLUDE_DIRECTORIES to allow warnings.
    set_property(TARGET ${TARGET_NAME}
                 PROPERTY INTERFACE_SYSTEM_INCLUDE_DIRECTORIES "")
  endforeach()
else()
  # If find_package works, GTEST_FOUND, GTEST_INCLUDE_DIRS, GTEST_BOTH_LIBRARIES
  # and GTEST_LIBRARIES will be set.
  find_package(GTest QUIET)
  if(GTEST_FOUND)
    message(STATUS "Found installed CMake Google Test package.")
  endif()
endif()

# ##############################################################################
# GTests

if(GTEST_FOUND AND NOT WP2_BUILD_WP2_JS)
  # Helpers.
  add_library(
    tests_helpers STATIC ${CMAKE_CURRENT_SOURCE_DIR}/include/helpers.cc
                         ${CMAKE_CURRENT_SOURCE_DIR}/include/helpers.h)
  target_link_libraries(
    tests_helpers
    PUBLIC webp2
    INTERFACE ${GTEST_BOTH_LIBRARIES}
    PRIVATE example_utils imageio)
  target_include_directories(tests_helpers PUBLIC ${GTEST_INCLUDE_DIRS}
                                                  ${GMOCK_INCLUDE_DIRS})

  add_library(
    tests_helpers_filter STATIC
    ${CMAKE_CURRENT_SOURCE_DIR}/include/helpers_filter.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/include/helpers_filter.h)
  target_link_libraries(tests_helpers_filter ${GTEST_BOTH_LIBRARIES} webp2
                        tests_helpers)

  add_library(
    tests_helpers_incr STATIC
    ${CMAKE_CURRENT_SOURCE_DIR}/include/helpers_incr.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/include/helpers_incr.h)
  target_link_libraries(tests_helpers_incr webp2 tests_helpers)

  macro(add_wp2_test TEST_FILE NUM_SHARDS)
    string(REPLACE "/" "_" TEST_NAME ${TEST_FILE})
    if(${NUM_SHARDS} GREATER 1)
      # gtest_discover_tests() could also be used but it logs too many tests.
      math(EXPR TEST_MAX_SHARD "${NUM_SHARDS} - 1")
      foreach(SHARD RANGE 0 ${TEST_MAX_SHARD})
        set(TEST_NAME_SHARDED ${TEST_NAME}_${SHARD})
        add_executable(${TEST_NAME_SHARDED} ${TEST_FILE}.cc)
        target_include_directories(${TEST_NAME_SHARDED}
                                   PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/..)
        target_link_libraries(${TEST_NAME_SHARDED} PRIVATE tests_helpers
                                                           ${ARGN})
        add_test(
          NAME ${TEST_NAME_SHARDED}
          COMMAND ${CMAKE_CURRENT_BINARY_DIR}/${TEST_NAME_SHARDED}
          WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})

        # GoogleTest needs two environment variables for sharding:
        set_property(
          TEST ${TEST_NAME_SHARDED}
          APPEND
          PROPERTY ENVIRONMENT "GTEST_TOTAL_SHARDS=${NUM_SHARDS};"
                   "GTEST_SHARD_INDEX=${SHARD}")
      endforeach(SHARD)
    else()
      add_executable(${TEST_NAME} ${TEST_FILE}.cc)
      target_include_directories(${TEST_NAME}
                                 PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/..)
      target_link_libraries(${TEST_NAME} PRIVATE tests_helpers ${ARGN})
      add_test(
        NAME ${TEST_NAME}
        COMMAND ${CMAKE_CURRENT_BINARY_DIR}/${TEST_NAME}
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
    endif()
  endmacro()

  add_wp2_test(test_alloc 1)
  add_wp2_test(test_alpha_filter 1 tests_helpers_filter)
  add_wp2_test(test_analysis 1)
  add_wp2_test(test_anim 4)
  add_wp2_test(test_anim_incr 4 tests_helpers_incr)
  add_wp2_test(test_anim_frame 4)
  add_wp2_test(test_anim_rectangle 1)
  add_wp2_test(test_ans 1 example_utils)
  add_wp2_test(test_aom_residuals 4)
  add_wp2_test(test_api 1)
  add_wp2_test(test_block_size 1)
  add_wp2_test(test_buffer 4)
  add_wp2_test(test_coded_block 1)
  add_wp2_test(test_color_precision 1)
  add_wp2_test(test_context 1)
  add_wp2_test(test_csp 1)
  add_wp2_test(test_data_source_stream 1)
  add_wp2_test(test_data_source 1)
  add_wp2_test(test_deblocking_filter 4 tests_helpers_filter)
  add_wp2_test(test_dec 4)
  add_wp2_test(test_decoder_api 1)
  add_wp2_test(test_diffusion 1)
  add_wp2_test(test_directional_filter 4 tests_helpers_filter)
  add_wp2_test(test_dsp 1)
  add_wp2_test(test_dsp_speed 1)
  add_wp2_test(test_example_utils 1 example_utils)
  add_wp2_test(test_formats 1)
  add_wp2_test(test_front_mgr 1)
  add_wp2_test(test_global_params 1)
  add_wp2_test(test_header 1)
  add_wp2_test(test_imageio_anim 1)
  add_wp2_test(test_imageio_util 1)
  add_wp2_test(test_incr_dec 4 tests_helpers_incr)
  add_wp2_test(test_intertile_filter 1 tests_helpers_filter)
  add_wp2_test(test_math 1)
  add_wp2_test(test_optim 1)
  add_wp2_test(test_orientation 1)
  add_wp2_test(test_partitions 4)
  add_wp2_test(test_partition_score_func 1)
  add_wp2_test(test_pessimization 1)
  add_wp2_test(test_plane 1 tests_helpers_filter)
  add_wp2_test(test_predictor 1)
  add_wp2_test(test_preview_api 1)
  add_wp2_test(test_progress 4)
  add_wp2_test(test_score 1)
  add_wp2_test(test_quant_mtx 1)
  add_wp2_test(test_residuals 1)
  add_wp2_test(test_restoration_filter 8 tests_helpers_filter)
  add_wp2_test(test_segment 1)
  add_wp2_test(test_simple_encode 1)
  add_wp2_test(test_slow_enc 1)
  add_wp2_test(test_speeds 4)
  add_wp2_test(test_symbols 1)
  add_wp2_test(test_syntax_writer 1)
  add_wp2_test(test_thread 1)
  add_wp2_test(test_tile_shapes 1)
  add_wp2_test(test_utils 1)
  add_wp2_test(test_uv_mode 4)
  add_wp2_test(test_vector 1)

  add_wp2_test(lossless/test_color_cache 1)
  add_wp2_test(lossless/test_group4 1)
  add_wp2_test(lossless/test_histogram 1)
  add_wp2_test(lossless/test_palette 1)

  if(WP2_BUILD_EXTRAS)
    add_wp2_test(test_custom_csp 4 extras)
    add_wp2_test(test_custom_csp_dec 4 extras)
    add_wp2_test(test_custom_csp_enc 4 extras)
    add_wp2_test(test_custom_csp_enc_dec 4 extras)
    add_wp2_test(test_distortion 4 extras tests_helpers_filter)
    add_wp2_test(test_extras 1 extras)
    add_wp2_test(test_imageio 1 extras)
    add_wp2_test(test_imageio_conversion 1 extras)
    add_wp2_test(test_transf 1 extras)
    add_wp2_test(test_version 1 extras)
  endif()
else()
  if(NOT WP2_BUILD_WP2_JS)
    message(WARNING "Didn't find Google Test, disabled some tests.")
  endif()
endif()

# ##############################################################################
# SWIG

if(WP2_BUILD_SWIG_PY)
  # To run the tests, 'execute' permission must be set.
  file(
    COPY test_stuttering_server.sh test_swig_examples.sh test_swig_python.py
    DESTINATION ${CMAKE_CURRENT_BINARY_DIR}
    FILE_PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ
                     GROUP_EXECUTE WORLD_EXECUTE)

  # test_stuttering_server doesn't work in a Docker environment (network issue).
  if(WP2_ENABLE_TEST_HTTP_SERVER)
    add_test(
      NAME test_stuttering_server
      COMMAND ${CMAKE_CURRENT_BINARY_DIR}/test_stuttering_server.sh
              ${CMAKE_SOURCE_DIR}/tests/tools
              ${CMAKE_CURRENT_SOURCE_DIR}/testdata
      WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
  endif()

  add_test(
    NAME test_swig_examples
    COMMAND ${CMAKE_CURRENT_BINARY_DIR}/test_swig_examples.sh
            ${CMAKE_SOURCE_DIR}/swig/examples ${CMAKE_BINARY_DIR}
            ${CMAKE_CURRENT_SOURCE_DIR}/testdata
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})

  add_test(
    NAME test_swig_python
    COMMAND ${CMAKE_CURRENT_BINARY_DIR}/test_swig_python.py
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})

  # Add the path to libwebp2 module in PYTHONPATH environment variable.
  set_tests_properties(
    test_swig_examples test_swig_python
    PROPERTIES ENVIRONMENT "PYTHONPATH=\$PYTHONPATH:${CMAKE_BINARY_DIR}/swig")
endif()

