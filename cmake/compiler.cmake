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

# Find Clang Tidy.
if(WP2_BUILD_WITH_CLANG_TIDY)
  find_program(
    CLANG_TIDY_EXE
    NAMES "clang-tidy"
    DOC "Path to clang-tidy executable")
  if(NOT CLANG_TIDY_EXE)
    message(STATUS "clang-tidy not found.")
  else()
    message(STATUS "clang-tidy found: ${CLANG_TIDY_EXE}")
    set(CMAKE_CXX_CLANG_TIDY
        "clang-tidy"
        "-config={Checks: 'readability-braces-around-statements, \
clang-analyzer-deadcode.DeadStores, \
-clang-analyzer-core.CallAndMessage', \
CheckOptions: [{key: readability-braces-around-statements.ShortStatementLines, \
value: 1}]}")
  endif()
endif()

# Deal with the compiler flags.
if(CMAKE_CXX_COMPILER_ID MATCHES "Clang" OR CMAKE_CXX_COMPILER_ID MATCHES "GNU")
  set(EXTRA_FLAGS "-Wall -Wmissing-declarations -Wshadow -Wformat=2")
  # Disable RTTI (useful for Skia).
  set(EXTRA_FLAGS "${EXTRA_FLAGS} -Werror -Wsign-compare -fno-rtti")
  if(CMAKE_CXX_COMPILER_ID MATCHES "GNU")
    if(CMAKE_CXX_COMPILER_VERSION VERSION_GREATER "5.3")
      set(EXTRA_FLAGS "${EXTRA_FLAGS} -Wsubobject-linkage")
    endif()
  endif()
  if(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    set(EXTRA_FLAGS
        "${EXTRA_FLAGS} -Wno-tautological-constant-out-of-range-compare")
  endif()
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${EXTRA_FLAGS}")
endif()
if(CMAKE_C_COMPILER_ID MATCHES "Clang" OR CMAKE_CXX_COMPILER_ID MATCHES "Clang")
  set(EXTRA_C_FLAGS "${EXTRA_FLAGS} -Wextra -Wold-style-definition")
  set(EXTRA_C_FLAGS "${EXTRA_C_FLAGS} -Wmissing-prototypes")
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${EXTRA_C_FLAGS}")
endif()

# prevents the warning: "warning: format string is not a string literal".
function(AllowNonLiteral varName)
  get_property(
    FLAGS
    TARGET ${varName}
    PROPERTY COMPILE_FLAGS)
  set(FLAGS "${FLAGS} -Wno-format-nonliteral")
  set_property(TARGET ${varName} PROPERTY COMPILE_FLAGS ${FLAGS})
endfunction()
