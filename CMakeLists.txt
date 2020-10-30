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

cmake_minimum_required(VERSION 3.1.3)
project(
  wp2
  LANGUAGES C CXX
  VERSION 0.1.0)
set(CMAKE_C_STANDARD 99)
set(CMAKE_CXX_STANDARD 11)

# Force some policies to adopt the new behavior.
if(POLICY CMP0072)
  # Use the new OpenGL.
  cmake_policy(SET CMP0072 NEW)
endif(POLICY CMP0072)
if(POLICY CMP0048)
  # Have project manage the version.
  cmake_policy(SET CMP0048 NEW)
endif(POLICY CMP0048)

# Options for coder / decoder executables.
option(WP2_ENABLE_SIMD "Enable any SIMD optimization." ON)
set(WP2_ENABLE_CONTEXT_SWITCH
    AUTO
    CACHE STRING "Use context switch optimization for incremental decoding.")
set_property(CACHE WP2_ENABLE_CONTEXT_SWITCH PROPERTY STRINGS AUTO ON OFF)
option(WP2_BUILD_TESTS "Build test targets." ON)
option(WP2_ENABLE_TESTS "Enable test targets." ON)
option(WP2_ENABLE_TEST_HTTP_SERVER "Enable stuttering_http_server test." OFF)
set(AOM_SRC_PATH CACHE STRING "Path to AOM sources.")
set(AOM_BUILD_PATH CACHE STRING "Path to AOM build directory.")
set(AVIF_SRC_PATH CACHE STRING "Path to AVIF sources.")
set(AVIF_BUILD_PATH CACHE STRING "Path to AVIF build directory.")
set(SJPEG_INSTALL_PATH CACHE STRING
                             "Path where sjpeg is installed (for rd-curve).")
# Debug options.
option(
  WP2_HAVE_AOM_DBG
  "Specifies AOM was built with -DCONFIG_ACCOUNTING=1 -DCONFIG_INSPECTION=1."
  OFF)
option(WP2_ENABLE_ANS_CHECK "Enable enc/dec ANS comparison." OFF)
option(WP2_ENABLE_ANS_DEEP_CHECK "Enable more verbose ANS and pixel comparison."
       OFF)
option(WP2_ENABLE_BITTRACE "Enable tracing." OFF)
option(WP2_ENABLE_TRACE "Enable tracing." OFF)
# Other build options.
option(WP2_BUILD_EXAMPLES "Build the cwp2 / dwp2 / vwp2 command line tools." ON)
option(WP2_BUILD_EXTRAS "Build the extras/ library." ON)
option(WP2_BUILD_WITH_CLANG_TIDY "Check the style using clang tidy: slower."
       OFF)
option(WP2_BUILD_WP2_JS "Emscripten build of wp2.js." OFF)
option(WP2_BUILD_SWIG_PY "Build SWIG python interface." OFF)

set(WP2_DEP_LIBRARIES)
set(WP2_DEP_INCLUDE_DIRS)
set(WP2_DEFINITIONS)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE
      "Release"
      CACHE STRING "Build type: Release, Debug or RelWithDebInfo" FORCE)
endif()
if(WP2_ENABLE_ANS_CHECK OR WP2_ENABLE_ANS_DEEP_CHECK)
  list(APPEND WP2_DEFINITIONS WP2_ENC_DEC_MATCH)
  if(WP2_ENABLE_ANS_DEEP_CHECK)
    list(APPEND WP2_DEFINITIONS WP2_ENC_DEC_DEEP_MATCH)
  endif()
endif()
if(WP2_ENABLE_BITTRACE)
  list(APPEND WP2_DEFINITIONS WP2_BITTRACE)
endif()
if(WP2_ENABLE_TRACE)
  list(APPEND WP2_DEFINITIONS WP2_TRACE)
endif()
if(NOT WP2_ENABLE_CONTEXT_SWITCH STREQUAL "AUTO")
  if(WP2_ENABLE_CONTEXT_SWITCH)
    list(APPEND WP2_DEFINITIONS WP2_USE_CONTEXT_SWITCH=1)
  else()
    list(APPEND WP2_DEFINITIONS WP2_USE_CONTEXT_SWITCH=0)
  endif()
endif()

# disable some configurations for WP2_JS
if(WP2_BUILD_WP2_JS)
  set(WP2_ENABLE_SIMD OFF)
  set(WP2_ENABLE_TESTS OFF)
endif()

if(NOT WP2_BUILD_TESTS)
  set(WP2_ENABLE_TESTS OFF)
endif()

# ##############################################################################
# Android only.
if(ANDROID)
  include_directories(${ANDROID_NDK}/sources/android/cpufeatures)
  add_library(cpufeatures STATIC
              ${ANDROID_NDK}/sources/android/cpufeatures/cpu-features.c)
  target_link_libraries(cpufeatures dl)
  set(WP2_DEP_LIBRARIES ${WP2_DEP_LIBRARIES} cpufeatures)
  set(WP2_DEP_INCLUDE_DIRS ${WP2_DEP_INCLUDE_DIRS}
                           ${ANDROID_NDK}/sources/android/cpufeatures)
endif()

# Check for SIMD extensions.
include(${CMAKE_CURRENT_LIST_DIR}/cmake/cpu.cmake)

# Check for thread support
include(${CMAKE_CURRENT_LIST_DIR}/cmake/thread.cmake)
if(WP2_USE_THREAD)
  list(APPEND WP2_DEFINITIONS WP2_USE_THREAD)
endif()

# Set compiler warnings.
include(${CMAKE_CURRENT_LIST_DIR}/cmake/compiler.cmake)

# Strip symbols when compiling in release to reduce binary size.
if(NOT WIN32 AND
   NOT APPLE AND
   (CMAKE_CXX_COMPILER_ID MATCHES "Clang" OR
    CMAKE_CXX_COMPILER_ID MATCHES "GNU"))
  set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} -s")
  set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -s")
endif()

# ##############################################################################
# WP2 various sub-libraries

# wp2_dsp library.
add_library(
  wp2_dsp OBJECT
  ${CMAKE_CURRENT_SOURCE_DIR}/src/dsp/alpha.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/src/dsp/ans_dsp.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/src/dsp/argb_converter.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/src/dsp/cpu.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/src/dsp/cdef_dsp.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/src/dsp/csp_dsp.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/src/dsp/deblocking_filter_dsp.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/src/dsp/dec_dsp.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/src/dsp/enc_dsp.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/src/dsp/grain_c.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/src/dsp/intra_c.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/src/dsp/math.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/src/dsp/math.h
  ${CMAKE_CURRENT_SOURCE_DIR}/src/dsp/ssim_dsp.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/src/dsp/quant_c.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/src/dsp/quant_sse.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/src/dsp/score_dsp.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/src/dsp/transf_c.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/src/dsp/transf_sse.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/src/dsp/wiener_dsp.cc
  # lossless
  ${CMAKE_CURRENT_SOURCE_DIR}/src/dsp/lossless/lossless_common.h
  ${CMAKE_CURRENT_SOURCE_DIR}/src/dsp/lossless/lossless_dsp_enc.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/src/dsp/lossless/lossless.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/src/dsp/lossless/lossless.h)
target_include_directories(wp2_dsp PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
target_compile_definitions(wp2_dsp PUBLIC ${WP2_DEFINITIONS})

# ANS library.
add_library(
  wp2_ans OBJECT
  ${CMAKE_CURRENT_SOURCE_DIR}/src/utils/ans.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/src/utils/ans.h
  ${CMAKE_CURRENT_SOURCE_DIR}/src/utils/ans_cost_table.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/src/utils/ans_utils.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/src/utils/ans_utils.h)
target_include_directories(wp2_ans PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
target_compile_definitions(wp2_ans PUBLIC ${WP2_DEFINITIONS})

# Preview library. TODO(yguyon): Split it as enc/dec
add_library(
  wp2_preview OBJECT
  ${CMAKE_CURRENT_SOURCE_DIR}/src/common/preview/preview.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/src/common/preview/preview.h
  ${CMAKE_CURRENT_SOURCE_DIR}/src/common/preview/preview_rasterizer.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/src/dec/preview/preview_dec.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/src/dec/preview/preview_dec.h
  ${CMAKE_CURRENT_SOURCE_DIR}/src/enc/preview/preview_analysis.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/src/enc/preview/preview_color.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/src/enc/preview/preview_config.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/src/enc/preview/preview_enc.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/src/enc/preview/preview_enc.h
  ${CMAKE_CURRENT_SOURCE_DIR}/src/enc/preview/preview_opt.cc)
target_include_directories(wp2_preview PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
target_compile_definitions(wp2_preview PUBLIC ${WP2_DEFINITIONS})

# utils library.
add_library(
  wp2_utils OBJECT
  ${CMAKE_CURRENT_SOURCE_DIR}/src/common/color_precision.h
  ${CMAKE_CURRENT_SOURCE_DIR}/src/common/constants.h
  ${CMAKE_CURRENT_SOURCE_DIR}/src/common/integral.h
  ${CMAKE_CURRENT_SOURCE_DIR}/src/common/integral.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/src/common/global_params.h
  ${CMAKE_CURRENT_SOURCE_DIR}/src/common/global_params.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/src/common/progress_watcher.h
  ${CMAKE_CURRENT_SOURCE_DIR}/src/common/progress_watcher.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/src/utils/argb_buffer.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/src/utils/context_switch.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/src/utils/context_switch.h
  ${CMAKE_CURRENT_SOURCE_DIR}/src/utils/csp.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/src/utils/csp.h
  ${CMAKE_CURRENT_SOURCE_DIR}/src/utils/data.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/src/utils/data_source.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/src/utils/data_source.h
  ${CMAKE_CURRENT_SOURCE_DIR}/src/utils/data_source_context.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/src/utils/data_source_context.h
  ${CMAKE_CURRENT_SOURCE_DIR}/src/utils/data_source_stream.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/src/utils/data_source_stream.h
  ${CMAKE_CURRENT_SOURCE_DIR}/src/utils/front_mgr.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/src/utils/front_mgr.h
  ${CMAKE_CURRENT_SOURCE_DIR}/src/utils/memory.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/src/utils/orientation.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/src/utils/orientation.h
  ${CMAKE_CURRENT_SOURCE_DIR}/src/utils/plane.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/src/utils/plane.h
  ${CMAKE_CURRENT_SOURCE_DIR}/src/utils/quantizer.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/src/utils/quantizer.h
  ${CMAKE_CURRENT_SOURCE_DIR}/src/utils/random.h
  ${CMAKE_CURRENT_SOURCE_DIR}/src/utils/stats.h
  ${CMAKE_CURRENT_SOURCE_DIR}/src/utils/status.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/src/utils/thread_utils.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/src/utils/thread_utils.h
  ${CMAKE_CURRENT_SOURCE_DIR}/src/utils/utils.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/src/utils/utils.h
  ${CMAKE_CURRENT_SOURCE_DIR}/src/utils/vector.h
  ${CMAKE_CURRENT_SOURCE_DIR}/src/utils/wiener.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/src/utils/wiener.h)
target_include_directories(wp2_utils PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
target_compile_definitions(wp2_utils PUBLIC ${WP2_DEFINITIONS})

# 'swapcontext' is deprecated on MacOS with gcc.
if(APPLE)
  set_source_files_properties(
    ${CMAKE_CURRENT_SOURCE_DIR}/src/utils/context_switch.cc PROPERTIES
    COMPILE_FLAGS -Wno-deprecated-declarations)
endif()

# wp2_enc_dec library.
add_library(
  wp2_enc_dec OBJECT
  # common lossy/lossless
  ${CMAKE_CURRENT_SOURCE_DIR}/src/dec/symbols_dec.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/src/dec/symbols_dec.h
  ${CMAKE_CURRENT_SOURCE_DIR}/src/enc/symbols_enc.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/src/enc/symbols_enc.h
  ${CMAKE_CURRENT_SOURCE_DIR}/src/dec/incr/decoder_api.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/src/dec/incr/decoder_context.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/src/dec/incr/decoder_context.h
  ${CMAKE_CURRENT_SOURCE_DIR}/src/dec/incr/decoder_info.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/src/dec/incr/decoder_info.h
  ${CMAKE_CURRENT_SOURCE_DIR}/src/dec/incr/decoder_skip.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/src/dec/incr/decoder_skip.h
  ${CMAKE_CURRENT_SOURCE_DIR}/src/dec/incr/decoder_stages.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/src/dec/incr/decoder_state.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/src/dec/incr/decoder_state.h
  ${CMAKE_CURRENT_SOURCE_DIR}/src/dec/incr/decoder_tiles.cc
  # common encode/decoder
  ${CMAKE_CURRENT_SOURCE_DIR}/src/common/filters/rstr_flt_params.h
  ${CMAKE_CURRENT_SOURCE_DIR}/src/common/filters/rstr_flt_params.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/src/common/lossy/aom/array_2d.h
  ${CMAKE_CURRENT_SOURCE_DIR}/src/common/lossy/aom/cdfs.inc
  ${CMAKE_CURRENT_SOURCE_DIR}/src/common/lossy/block.h
  ${CMAKE_CURRENT_SOURCE_DIR}/src/common/lossy/block.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/src/common/lossy/block_size.h
  ${CMAKE_CURRENT_SOURCE_DIR}/src/common/lossy/block_size.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/src/common/lossy/block_size_io.h
  ${CMAKE_CURRENT_SOURCE_DIR}/src/common/lossy/block_size_io.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/src/common/lossy/context.h
  ${CMAKE_CURRENT_SOURCE_DIR}/src/common/lossy/context.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/src/common/lossy/predictor.h
  ${CMAKE_CURRENT_SOURCE_DIR}/src/common/lossy/predictor.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/src/common/lossy/quant_mtx.h
  ${CMAKE_CURRENT_SOURCE_DIR}/src/common/lossy/quant_mtx.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/src/common/lossy/residuals.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/src/common/lossy/residuals.h
  ${CMAKE_CURRENT_SOURCE_DIR}/src/common/lossy/residuals_aom.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/src/common/lossy/residuals_aom.h
  ${CMAKE_CURRENT_SOURCE_DIR}/src/common/lossy/segment.h
  ${CMAKE_CURRENT_SOURCE_DIR}/src/common/lossy/segment.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/src/common/lossy/transforms.h
  ${CMAKE_CURRENT_SOURCE_DIR}/src/common/lossy/transforms.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/src/common/lossy/rnd_mtx.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/src/common/symbols.h
  ${CMAKE_CURRENT_SOURCE_DIR}/src/common/symbols.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/src/common/vdebug.cc
  # lossy
  ${CMAKE_CURRENT_SOURCE_DIR}/src/dec/alpha_dec.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/src/dec/anim_dec.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/src/dec/bypass_dec.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/src/dec/features_dec.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/src/dec/filters/alpha_filter.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/src/dec/filters/alpha_filter.h
  ${CMAKE_CURRENT_SOURCE_DIR}/src/dec/filters/block_map_filter.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/src/dec/filters/block_map_filter.h
  ${CMAKE_CURRENT_SOURCE_DIR}/src/dec/filters/deblocking_filter.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/src/dec/filters/deblocking_filter.h
  ${CMAKE_CURRENT_SOURCE_DIR}/src/dec/filters/directional_filter.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/src/dec/filters/directional_filter.h
  ${CMAKE_CURRENT_SOURCE_DIR}/src/dec/filters/grain_filter.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/src/dec/filters/grain_filter.h
  ${CMAKE_CURRENT_SOURCE_DIR}/src/dec/filters/intertile_filter.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/src/dec/filters/intertile_filter.h
  ${CMAKE_CURRENT_SOURCE_DIR}/src/dec/filters/intratile_filter.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/src/dec/filters/intratile_filter.h
  ${CMAKE_CURRENT_SOURCE_DIR}/src/dec/filters/restoration_filter.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/src/dec/filters/restoration_filter.h
  ${CMAKE_CURRENT_SOURCE_DIR}/src/dec/lossless_dec.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/src/dec/lossy_dec.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/src/dec/main_dec.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/src/dec/residuals_dec.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/src/dec/residuals_dec_aom.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/src/dec/residuals_dec_aom.h
  ${CMAKE_CURRENT_SOURCE_DIR}/src/dec/syntax_dec.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/src/dec/tile_dec.h
  ${CMAKE_CURRENT_SOURCE_DIR}/src/dec/tile_dec.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/src/enc/alpha_enc.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/src/enc/analysis_enc.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/src/enc/anim/anim_enc.h
  ${CMAKE_CURRENT_SOURCE_DIR}/src/enc/anim/anim_enc.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/src/enc/anim/anim_frame.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/src/enc/anim/anim_rectangle.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/src/enc/bypass_enc.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/src/enc/config_enc.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/src/enc/distortion_enc.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/src/enc/lossless_enc.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/src/enc/lossy_enc.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/src/enc/main_enc.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/src/enc/nearlossless_enc.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/src/enc/partition_enc.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/src/enc/partition_score_func.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/src/enc/partitioner.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/src/enc/predictor_enc.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/src/enc/residuals_enc.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/src/enc/residuals_enc_aom.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/src/enc/residuals_enc_aom.h
  ${CMAKE_CURRENT_SOURCE_DIR}/src/enc/screen_content/screen_enc.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/src/enc/screen_content/screen_enc.h
  ${CMAKE_CURRENT_SOURCE_DIR}/src/enc/segment_enc.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/src/enc/syntax_enc.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/src/enc/writer_enc.cc
  # lossless
  ${CMAKE_CURRENT_SOURCE_DIR}/src/common/lossless/color_cache.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/src/common/lossless/color_cache.h
  ${CMAKE_CURRENT_SOURCE_DIR}/src/dec/lossless/group4_dec.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/src/dec/lossless/losslessi_dec.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/src/dec/lossless/losslessi_dec.h
  ${CMAKE_CURRENT_SOURCE_DIR}/src/enc/lossless/analysisl_enc.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/src/enc/lossless/backward_references_enc.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/src/enc/lossless/backward_references_cost_enc.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/src/enc/lossless/backward_references_enc.h
  ${CMAKE_CURRENT_SOURCE_DIR}/src/enc/lossless/group4_enc.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/src/enc/lossless/histogram_enc.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/src/enc/lossless/histogram_enc.h
  ${CMAKE_CURRENT_SOURCE_DIR}/src/enc/lossless/palette.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/src/enc/lossless/palette.h
  ${CMAKE_CURRENT_SOURCE_DIR}/src/enc/lossless/predictorl_enc.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/src/enc/lossless/losslessi_enc.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/src/enc/lossless/losslessi_enc.h)
target_include_directories(wp2_enc_dec PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
target_compile_definitions(wp2_enc_dec PUBLIC ${WP2_DEFINITIONS})

# ##############################################################################
# the WP2 library eventually.

add_library(
  webp2
  $<TARGET_OBJECTS:wp2_ans> $<TARGET_OBJECTS:wp2_preview>
  $<TARGET_OBJECTS:wp2_dsp> $<TARGET_OBJECTS:wp2_utils>
  $<TARGET_OBJECTS:wp2_enc_dec>)
if(WP2_DEP_LIBRARIES)
  target_link_libraries(webp2 ${WP2_DEP_LIBRARIES})
endif()
if(WP2_USE_THREAD)
  target_link_libraries(webp2 Threads::Threads)
endif()
target_compile_definitions(webp2 PUBLIC ${WP2_DEFINITIONS})

# Make sure the OBJECT libraries are built with position independent code (it is
# not ON by default).
set_target_properties(wp2_dsp wp2_ans wp2_preview wp2_utils wp2_enc_dec webp2
                      PROPERTIES POSITION_INDEPENDENT_CODE ON)

# Set the version numbers.
set_target_properties(webp2 PROPERTIES VERSION ${PROJECT_VERSION} SOVERSION
                                     ${PROJECT_VERSION})
target_include_directories(
  webp2 INTERFACE $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
                $<INSTALL_INTERFACE:include>)

# Install the different headers and libraries.
install(
  FILES ${CMAKE_CURRENT_SOURCE_DIR}/src/wp2/base.h
        ${CMAKE_CURRENT_SOURCE_DIR}/src/wp2/debug.h
        ${CMAKE_CURRENT_SOURCE_DIR}/src/wp2/decode.h
        ${CMAKE_CURRENT_SOURCE_DIR}/src/wp2/encode.h
        ${CMAKE_CURRENT_SOURCE_DIR}/src/wp2/format_constants.h
  DESTINATION include/wp2)

# ##############################################################################

# build the imageio library
add_library(
  imageio
  ${CMAKE_CURRENT_SOURCE_DIR}/imageio/anim_image_dec.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/imageio/anim_image_dec.h
  ${CMAKE_CURRENT_SOURCE_DIR}/imageio/bmpenc.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/imageio/bmpdec.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/imageio/file_format.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/imageio/file_format.h
  ${CMAKE_CURRENT_SOURCE_DIR}/imageio/gifdec.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/imageio/gifdec_metadata.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/imageio/image_dec.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/imageio/image_dec.h
  ${CMAKE_CURRENT_SOURCE_DIR}/imageio/image_enc.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/imageio/image_enc.h
  ${CMAKE_CURRENT_SOURCE_DIR}/imageio/imageio_util.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/imageio/imageio_util.h
  ${CMAKE_CURRENT_SOURCE_DIR}/imageio/jpegdec.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/imageio/pngdec.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/imageio/pngenc.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/imageio/pnmdec.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/imageio/pnmenc.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/imageio/tiffdec.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/imageio/tiffenc.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/imageio/webpdec.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/imageio/webpenc.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/imageio/wp2dec.cc)
target_link_libraries(imageio webp2)

# Find the standard image libraries.
foreach(I_LIB PNG JPEG TIFF GIF WebP)
  find_package(${I_LIB} QUIET)
  set(WP2_HAVE_${I_LIB} ${${I_LIB}_FOUND})
  if(${I_LIB}_FOUND)
    if(${I_LIB} STREQUAL "JPEG")
      target_link_libraries(imageio ${JPEG_LIBRARIES})
      target_include_directories(imageio PRIVATE ${JPEG_INCLUDE_DIRS})
    elseif(I_LIB STREQUAL "GIF")
      if(WP2_HAVE_GIF)
        # GIF find_package only locates the header and library, it doesn't fail
        # compile tests when detecting the version, but falls back to 3 (as of
        # at least cmake 3.7.2). Make sure the library links to avoid incorrect
        # detection when cross compiling.
        cmake_push_check_state()
        set(CMAKE_REQUIRED_LIBRARIES ${GIF_LIBRARIES})
        set(CMAKE_REQUIRED_INCLUDES ${GIF_INCLUDE_DIR})
        check_c_source_compiles(
          "
          #include <gif_lib.h>
          int main(void) {
            (void)DGifOpenFileHandle;
            return 0;
          }
          "
          GIF_COMPILES)
        cmake_pop_check_state()
        if(GIF_COMPILES)
          target_link_libraries(imageio ${GIF_LIBRARIES})
          target_include_directories(imageio PRIVATE ${GIF_INCLUDE_DIR})
        else()
          unset(WP2_HAVE_GIF)
        endif()
      endif()
    elseif(I_LIB STREQUAL "WebP")
      if(TARGET WebP::webp)
        message(STATUS "Found new CMake WebP")
        target_link_libraries(imageio WebP::webp WebP::webpdemux)
      else()
        message(STATUS "Found old CMake WebP")
        target_link_libraries(imageio ${WEBP_LIBRARIES})
        target_include_directories(imageio PRIVATE ${WEBP_INCLUDE_DIRS})
      endif()
      set(WP2_HAVE_WEBP ${WebP_FOUND})
    else()
      target_link_libraries(imageio ${I_LIB}::${I_LIB})
    endif()
  elseif(I_LIB STREQUAL "WebP")
    message(STATUS "Looking for non-CMake WebP")
    # Find WebP in case it was not found with CMake.
    find_path(WEBP_INCLUDE_DIR NAMES webp/decode.h)
    find_path(WEBP_DEMUX_INCLUDE_DIR NAMES webp/demux.h)
    find_library(WEBP_LIBRARY NAMES webp)
    find_library(WEBP_DEMUX_LIBRARY NAMES webpdemux)
    if(WEBP_INCLUDE_DIR
       AND WEBP_DEMUX_INCLUDE_DIR
       AND WEBP_LIBRARY
       AND WEBP_DEMUX_LIBRARY)
      message(STATUS "Found non-cmake WebP: ${WEBP_LIBRARY}")
      target_include_directories(
        imageio
        PRIVATE ${WEBP_INCLUDE_DIR} ${WEBP_DEMUX_INCLUDE_DIR}
        PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>)
      target_link_libraries(imageio ${WEBP_LIBRARY} ${WEBP_DEMUX_LIBRARY})
      set(WebP_FOUND TRUE)
    endif()
    set(WP2_HAVE_WEBP ${WebP_FOUND})
  endif()
endforeach()

target_include_directories(
  imageio PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>)

if(WP2_HAVE_JPEG)
  target_compile_definitions(imageio PUBLIC WP2_HAVE_JPEG)
endif()
if(WP2_HAVE_PNG)
  target_compile_definitions(imageio PUBLIC WP2_HAVE_PNG)
endif()
if(WP2_HAVE_WEBP)
  target_compile_definitions(imageio PUBLIC WP2_HAVE_WEBP)
endif()
if(WP2_HAVE_TIFF)
  target_compile_definitions(imageio PUBLIC WP2_HAVE_TIFF)
endif()
if(WP2_HAVE_GIF)
  target_compile_definitions(imageio PUBLIC WP2_HAVE_GIF)
endif()
if(WP2_USE_THREAD)
  target_link_libraries(imageio Threads::Threads)
endif()
set_target_properties(imageio PROPERTIES POSITION_INDEPENDENT_CODE ON)

# ##############################################################################

find_program(iwyu_path NAMES include-what-you-use iwyu)
# Uncomment below, install iwpu, and compile with clang to have warnings for
# extra includes.
if(iwyu_path)
  set_property(TARGET wp2_ans PROPERTY CXX_INCLUDE_WHAT_YOU_USE ${iwyu_path})
  set_property(TARGET wp2_dsp PROPERTY CXX_INCLUDE_WHAT_YOU_USE ${iwyu_path})
  set_property(TARGET wp2_enc_dec PROPERTY CXX_INCLUDE_WHAT_YOU_USE
                                           ${iwyu_path})
  set_property(TARGET wp2_preview PROPERTY CXX_INCLUDE_WHAT_YOU_USE
                                           ${iwyu_path})
  set_property(TARGET wp2_utils PROPERTY CXX_INCLUDE_WHAT_YOU_USE ${iwyu_path})
  set_property(TARGET webp2 PROPERTY CXX_INCLUDE_WHAT_YOU_USE ${iwyu_path})
  set_property(TARGET imageio PROPERTY CXX_INCLUDE_WHAT_YOU_USE ${iwyu_path})
endif()

# ##############################################################################

# Build the executables if asked for.
if(WP2_BUILD_EXAMPLES)
  # utils
  add_library(
    example_utils STATIC
    ${CMAKE_CURRENT_SOURCE_DIR}/examples/example_utils.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/examples/example_utils.h
    ${CMAKE_CURRENT_SOURCE_DIR}/examples/stopwatch.h)
  target_link_libraries(example_utils imageio webp2)
  if(WP2_DEFINITIONS)
    target_compile_definitions(example_utils PUBLIC ${WP2_DEFINITIONS})
  endif()

  # cwp2
  add_executable(cwp2 ${CMAKE_CURRENT_SOURCE_DIR}/examples/cwp2.cc)
  allownonliteral(cwp2)
  target_link_libraries(cwp2 webp2 imageio example_utils)

  # dwp2
  add_executable(dwp2 ${CMAKE_CURRENT_SOURCE_DIR}/examples/dwp2.cc)
  allownonliteral(dwp2)
  target_link_libraries(dwp2 webp2 imageio example_utils)

  # vwp2 (only if OpenGL is present)

  # Find the OpenGL/GLUT libraries.
  if(NOT WP2_BUILD_WP2_JS)
    find_package(OpenGL)
    find_package(GLUT)
  endif()

  if(OPENGL_FOUND
     AND GLUT_FOUND
     AND NOT WP2_BUILD_WP2_JS)
    add_executable(vwp2 ${CMAKE_CURRENT_SOURCE_DIR}/examples/vwp2.cc)
    target_include_directories(vwp2 PRIVATE ${OPENGL_INCLUDE_DIR}
                                            ${GLUT_INCLUDE_DIR})
    target_link_libraries(vwp2 ${OPENGL_LIBRARIES} ${GLUT_LIBRARIES} webp2
                          imageio example_utils)
    if(AOM_SRC_PATH)
      target_link_libraries(vwp2 extras)
      target_compile_definitions(vwp2 PRIVATE WP2_HAVE_AOM)
      if(WP2_HAVE_AOM_DBG)
        target_compile_definitions(vwp2 PRIVATE WP2_HAVE_AOM_DBG)
      endif()
    endif()
    if(WP2_HAVE_WEBP)
      target_compile_definitions(vwp2 PRIVATE WP2_HAVE_WEBP)
    endif()
    # check whether we need to include GLUT/glut.h or GL/glut.h
    include(CheckIncludeFileCXX)
    check_include_file_cxx(GLUT/glut.h HAVE_GLUT_GLUT_H)
    check_include_file_cxx(GL/glut.h HAVE_GL_GLUT_H)
    if(HAVE_GLUT_GLUT_H)
      list(APPEND WP2_DEFINITIONS HAVE_GLUT_GLUT_H)
    elseif(HAVE_GL_GLUT_H)
      list(APPEND WP2_DEFINITIONS HAVE_GL_GLUT_H)
    endif()
    target_compile_definitions(
      vwp2 PRIVATE ${OPENGL_DEFINITIONS} WP2_HAVE_OPENGL ${GLUT_DEFINITIONS}
      ${WP2_DEFINITIONS})

    # Force to link against pthread.
    set(THREADS_PREFER_PTHREAD_FLAG ON)
    find_package(Threads)
    include(CheckCXXSourceCompiles)
    include(CMakePushCheckState)
    cmake_push_check_state()
    set(CMAKE_REQUIRED_FLAGS "-Wl,--no-as-needed")
    check_cxx_source_compiles("int main(void){return 0;}" FLAG_NO_AS_NEEDED)
    cmake_pop_check_state()
    if(FLAG_NO_AS_NEEDED)
      target_link_libraries(vwp2 "-Wl,--no-as-needed")
    endif()
    if(APPLE)
      target_compile_options(vwp2 PUBLIC "-Wno-deprecated-declarations")
    endif()
  else() # Build empty vwp2 for test_cmd.
    add_executable(vwp2 ${CMAKE_CURRENT_SOURCE_DIR}/examples/vwp2.cc)
    target_link_libraries(vwp2 webp2 imageio example_utils)
  endif()
  allownonliteral(vwp2)

  # get_disto
  add_executable(get_disto ${CMAKE_CURRENT_SOURCE_DIR}/examples/get_disto.cc)
  target_link_libraries(get_disto webp2 imageio example_utils)
  allownonliteral(get_disto)

  # rd_curve
  add_executable(rd_curve ${CMAKE_CURRENT_SOURCE_DIR}/examples/rd_curve.cc)
  target_link_libraries(rd_curve webp2 imageio example_utils extras)
  if(WP2_HAVE_WEBP)
    target_compile_definitions(rd_curve PRIVATE WP2_HAVE_WEBP)
  endif()
  if(SJPEG_INSTALL_PATH)
    target_link_libraries(rd_curve ${SJPEG_INSTALL_PATH}/lib/libsjpeg.a)
    target_include_directories(rd_curve PRIVATE ${SJPEG_INSTALL_PATH}/include/)
    target_compile_definitions(rd_curve PRIVATE WP2_HAVE_SJPEG)
  endif()
  allownonliteral(rd_curve)

  install(TARGETS cwp2 dwp2 RUNTIME DESTINATION bin)

  # ############################################################################
  # Tools (that could be moved/removed later)

  # libaom needs to be compiled with CONFIG_INSPECTION now:
  # ~~~
  #   cd aom/build
  #   cmake ../ -DCONFIG_ACCOUNTING=1 -DCONFIG_INSPECTION=1
  # ~~~

  if((AOM_SRC_PATH AND AOM_BUILD_PATH) OR AVIF_BUILD_PATH)
    add_executable(av1enc extras/av1enc.cc)
    target_link_libraries(av1enc webp2 extras example_utils imageio)
  endif()

  add_executable(mk_preview extras/mk_preview.cc)
  target_link_libraries(mk_preview webp2 extras example_utils imageio)

  add_executable(ewp2 extras/ewp2.cc)
  target_link_libraries(ewp2 webp2 extras example_utils imageio)
endif()

# ##############################################################################

if(WP2_BUILD_EXTRAS)
  add_library(
    extras
    ${CMAKE_CURRENT_SOURCE_DIR}/extras/aom_utils.h
    ${CMAKE_CURRENT_SOURCE_DIR}/extras/aom_utils.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/extras/ccsp_imageio.h
    ${CMAKE_CURRENT_SOURCE_DIR}/extras/ccsp_imageio.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/extras/extras.h
    ${CMAKE_CURRENT_SOURCE_DIR}/extras/extras.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/extras/wp2_dec_12b.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/extras/y4mdec.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/extras/y4menc.cc)
  set_target_properties(extras PROPERTIES POSITION_INDEPENDENT_CODE ON)
  target_include_directories(extras PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
  target_link_libraries(extras example_utils webp2)
  # AOM code needs WebP code
  if(AOM_BUILD_PATH
     AND AOM_SRC_PATH
     AND WP2_HAVE_WEBP)
    target_link_directories(extras PUBLIC ${AOM_BUILD_PATH})
    target_link_libraries(extras libaom.a)
    target_include_directories(extras PRIVATE ${AOM_SRC_PATH}/)
    target_compile_definitions(extras PRIVATE WP2_HAVE_AOM WP2_HAVE_WEBP)
    if(WP2_HAVE_AOM_DBG)
      if(NOT AOM_SRC_PATH OR NOT AOM_BUILD_PATH)
        message(FATAL_ERROR "AOM_SRC_PATH and AOM_BUILD_PATH need to be set.")
      endif()
      target_include_directories(extras SYSTEM PRIVATE ${AOM_SRC_PATH}
                                                       ${AOM_BUILD_PATH})
      target_compile_definitions(extras PUBLIC WP2_HAVE_AOM_DBG)
    endif()
  endif()
  # Deal with AVIF.
  if (AVIF_BUILD_PATH)
    target_include_directories(extras SYSTEM PRIVATE ${AVIF_SRC_PATH}/include)
    target_link_directories(extras PUBLIC ${AVIF_BUILD_PATH})
    target_link_libraries(extras avif m Threads::Threads)
    target_compile_definitions(extras PUBLIC WP2_HAVE_AVIF)
  endif()
endif()

# ##############################################################################

if(WP2_BUILD_WP2_JS)
  # JavaScript version
  add_executable(wp2_js ${CMAKE_CURRENT_SOURCE_DIR}/extras/wp2_to_sdl.cc)
  target_link_libraries(wp2_js webp2 SDL)
  target_include_directories(wp2_js PRIVATE ${CMAKE_CURRENT_BINARY_DIR}
                                            ${CMAKE_CURRENT_SOURCE_DIR}/src)
  set(WP2_HAVE_SDL 1)
  set_target_properties(
    wp2_js
    PROPERTIES LINK_FLAGS
               "-s EXPORTED_FUNCTIONS='[\"_WP2ToSDL\"]' -s INVOKE_RUN=0 \
       -s EXTRA_EXPORTED_RUNTIME_METHODS='[\"cwrap\"]'")
  set_target_properties(wp2_js PROPERTIES OUTPUT_NAME webp2)
  target_compile_definitions(wp2_js PUBLIC EMSCRIPTEN WP2_HAVE_SDL)

  # WASM version
  add_executable(wp2_wasm ${CMAKE_CURRENT_SOURCE_DIR}/extras/wp2_to_sdl.cc)
  target_link_libraries(wp2_wasm webp2 SDL)
  target_include_directories(wp2_wasm PRIVATE ${CMAKE_CURRENT_BINARY_DIR}
                                              ${CMAKE_CURRENT_SOURCE_DIR}/src)
  set(WP2_HAVE_SDL 1)
  set_target_properties(
    wp2_wasm
    PROPERTIES
      LINK_FLAGS
      "-s WASM=1 -s ALLOW_MEMORY_GROWTH=1 \
       -s 'BINARYEN_METHOD=\"native-wasm\"' \
       -s EXPORTED_FUNCTIONS='[\"_WP2ToSDL\"]' -s INVOKE_RUN=0 \
       -s EXTRA_EXPORTED_RUNTIME_METHODS='[\"cwrap\"]'")
  target_compile_definitions(wp2_wasm PUBLIC EMSCRIPTEN WP2_HAVE_SDL)

  target_compile_definitions(wp2_dsp PUBLIC EMSCRIPTEN)
endif()

# ##############################################################################

if(WP2_BUILD_SWIG_PY)
  add_subdirectory(swig)
endif()

# ##############################################################################

if(WP2_BUILD_TESTS)
  if(WP2_ENABLE_TESTS)
    enable_testing()
  endif()
  # Build tests even if they are not run.
  add_subdirectory(tests)
endif()

# ##############################################################################

# Install the exported targets.
install(
  TARGETS webp2
  EXPORT ${PROJECT_NAME}_targets
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib)

if(ANDROID)
  install(
    TARGETS cpufeatures
    EXPORT ${PROJECT_NAME}_targets
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib)
endif()

set(ConfigPackageLocation lib/wp2/cmake/)

# Install the exports.
install(
  EXPORT wp2_targets
  NAMESPACE ${PROJECT_NAME}::
  FILE ${PROJECT_NAME}Targets.cmake
  DESTINATION ${ConfigPackageLocation})

# Create the CMake version file.
include(CMakePackageConfigHelpers)
write_basic_package_version_file(
  "${CMAKE_CURRENT_BINARY_DIR}/wp2ConfigVersion.cmake"
  VERSION ${PROJECT_VERSION}
  COMPATIBILITY AnyNewerVersion)

# Create the Config file.
include(CMakePackageConfigHelpers)
set(ConfigPackageLocation lib/${PROJECT_NAME}/cmake/)
configure_package_config_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/cmake/${PROJECT_NAME}Config.cmake.in
  ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}Config.cmake INSTALL_DESTINATION
  ${ConfigPackageLocation})

# Install the generated CMake files.
install(FILES "${CMAKE_CURRENT_BINARY_DIR}/wp2ConfigVersion.cmake"
              "${CMAKE_CURRENT_BINARY_DIR}/wp2Config.cmake"
        DESTINATION ${ConfigPackageLocation})

# ##############################################################################
# Man pages.

install(FILES ${CMAKE_CURRENT_SOURCE_DIR}/man/cwp2.1
              ${CMAKE_CURRENT_SOURCE_DIR}/man/dwp2.1 DESTINATION share/man/man1)
