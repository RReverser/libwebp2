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

#include <algorithm>
#include <cstdarg>
#include <cstdio>
#include <cstring>
#include <memory>

//------------------------------------------------------------------------------

int WP2GetVersion() { return WP2_VERSION; }
int WP2GetABIVersion() { return WP2_ABI_VERSION; }

void WP2ErrorLog(const char* const file, int line, WP2Status s,
                 const char* const extra_message) {
#if defined(WP2_TRACE) || defined(WP2_ERROR_TRACE)
  assert(wp2_error_tracer != nullptr);
  wp2_error_tracer->Log(file, line, s, extra_message);
#else
  (void)file;
  (void)line;
  (void)s;
  (void)extra_message;
#endif  // defined(WP2_TRACE) || defined(WP2_ERROR_TRACE)
}

#if defined(WP2_TRACE) || defined(WP2_ERROR_TRACE)
void WP2ErrorTracer::Log(const char* const file, int line, WP2Status s,
                         const char* const extra_message) {
  // TODO(skal): mutex?
  fprintf(stderr, "Status check error");
  if (file != nullptr) {
    fprintf(stderr, " at %s:%d", file, line);
  }
  fprintf(stderr, ". Got: %s.\n", WP2GetStatusMessage(s));
  if (extra_message != nullptr) {
    fprintf(stderr, "%s.\n", extra_message);
  }
}

static WP2ErrorTracer default_error_tracer;
WP2ErrorTracer* wp2_error_tracer = &default_error_tracer;
#endif  // defined(WP2_TRACE) || defined(WP2_ERROR_TRACE)

//------------------------------------------------------------------------------

[[gnu::format (printf, 1, 0)]]
void WP2Print(const char* const fmt, ...) {
  va_list arg;
  va_start(arg, fmt);
  vprintf(fmt, arg);
  va_end(arg);
}

#if defined(WP2_BITTRACE) || defined(WP2_TRACE) || defined(WP2_ENC_DEC_MATCH)
[[gnu::format(printf, 2, 3)]]
void WP2SPrint(std::string* const str, const char* const fmt, ...) {
  va_list list;
  va_start(list, fmt);
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wformat-nonliteral"
  const size_t size = vsnprintf(nullptr, 0, fmt, list);
  va_end(list);
  std::unique_ptr<char[]> tmp(new char[1 + size]);
  va_start(list, fmt);
  vsprintf(tmp.get(), fmt, list);
  va_end(list);
  *str = std::string(tmp.get(), tmp.get() + size);
#pragma GCC diagnostic pop
}
#endif

// provide a default printing mechanism in case none is supplied
#if defined(WP2_TRACE) && !defined(WP2_HAVE_WP2TRACE)
// we use gnu::format to avoid error message "format not a string literal..."
[[gnu::format (printf, 1, 0)]]
void WP2Trace(const char* format, const std::string prefix, ...) {
  fprintf(stderr, "%s", prefix.c_str());
  va_list list;
  va_start(list, prefix);
  vfprintf(stderr, format, list);
  fprintf(stderr, "\n");
  va_end(list);
}
#endif   // WP2_TRACE && !WP2_HAVE_WP2TRACE

//------------------------------------------------------------------------------

namespace WP2 {

const char* strnstr(const char* str, size_t count, const char* const target) {
  if (str == nullptr || target == nullptr) return nullptr;
  const size_t target_len = std::strlen(target);
  if (target_len == 0) return str;  // Empty 'target' matches everything.
  for (; count >= target_len && *str != '\0'; ++str, --count) {
    if (std::strncmp(str, target, target_len) == 0) return str;
  }
  return nullptr;
}

int strlcmp(const char* lhs, const char* rhs, size_t count) {
  return std::strncmp(lhs, rhs, std::min(count, std::strlen(rhs)));
}

}  // namespace WP2

//------------------------------------------------------------------------------
