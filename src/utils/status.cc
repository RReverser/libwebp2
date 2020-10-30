// Copyright 2019 Google LLC
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//     https://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.
// -----------------------------------------------------------------------------
//
// Misc. common utility functions
//
// Author: Skal (pascal.massimino@gmail.com)

#include "src/utils/utils.h"

static const char* kStatusMessage[WP2_STATUS_LAST + 1] = {
  "WP2_STATUS_OK",
  "WP2_STATUS_VERSION_MISMATCH",
  "WP2_STATUS_OUT_OF_MEMORY",
  "WP2_STATUS_INVALID_PARAMETER",
  "WP2_STATUS_NULL_PARAMETER",
  "WP2_STATUS_BAD_DIMENSION",
  "WP2_STATUS_USER_ABORT",
  "WP2_STATUS_UNSUPPORTED_FEATURE",

  "WP2_STATUS_BITSTREAM_ERROR",
  "WP2_STATUS_NOT_ENOUGH_DATA",
  "WP2_STATUS_BAD_READ",
  "WP2_STATUS_NEURAL_DECODE_FAILURE",

  "WP2_STATUS_BITSTREAM_OUT_OF_MEMORY",
  "WP2_STATUS_INVALID_CONFIGURATION",
  "WP2_STATUS_BAD_WRITE",
  "WP2_STATUS_FILE_TOO_BIG",
  "WP2_STATUS_INVALID_COLORSPACE",
  "WP2_STATUS_NEURAL_ENCODE_FAILURE",

  "WP2_STATUS_LAST"
};

static const char* kStatusText[WP2_STATUS_LAST + 1] = {
  "OK. no error.",
  "mismatch between user and system library version",
  "memory error allocating objects",
  "a parameter value is invalid",
  "a pointer parameter is NULL",
  "picture has invalid width/height",
  "abort request by user",
  "unsupported feature",

  "bitstream has syntactic error",
  "premature end-of-file during decoding",
  "error while reading bytes",
  "neural decode failure",

  "memory error while flushing bits",
  "encoding configuration is invalid",
  "error while flushing bytes",
  "file is bigger than 4G",
  "encoder called with improper colorspace",
  "neural encode failure",

  "last",
};

const char* WP2GetStatusMessage(WP2Status status) {
  return kStatusMessage[status];
}

const char* WP2GetStatusText(WP2Status status) {
  return kStatusText[status];
}

//------------------------------------------------------------------------------
