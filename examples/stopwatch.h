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
//  Helper functions to measure elapsed time.
//
// Author: Mikolaj Zalewski (mikolajz@google.com)

#ifndef WP2_EXAMPLES_STOPWATCH_H_
#define WP2_EXAMPLES_STOPWATCH_H_

#include "src/wp2/base.h"

#if defined _WIN32 && !defined __GNUC__
#include <windows.h>

static inline double GetStopwatchTime() {
  LARGE_INTEGER watch;
  LARGE_INTEGER freq;
  if (!QueryPerformanceCounter(&watch))
    return 0.0;
  if (!QueryPerformanceFrequency(&freq))
    return 0.0;
  if (freq.QuadPart == 0)
    return 0.0;
  return watch.QuadPart / (double)(freq.QuadPart);
}

#else    /* !_WIN32 */
#include <sys/time.h>

static inline double GetStopwatchTime() {
  struct timeval watch;
  gettimeofday(&watch, NULL);
  const double sec = (double)(watch.tv_sec);
  const double usec = (double)(watch.tv_usec);
  return sec + usec / 1000000.0;
}

#endif   /* _WIN32 */

#endif  /* WP2_EXAMPLES_STOPWATCH_H_ */
