// Copyright 2020 Google LLC
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//      http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

// Worker and ThreadLock test.

#include <tuple>

#include "examples/stopwatch.h"
#include "include/helpers.h"
#include "src/dsp/dsp.h"
#include "src/utils/thread_utils.h"
#include "src/utils/utils.h"

namespace WP2 {
namespace {

const double kRunDuration = 0.5;

class TestWorker : public WorkerBase {
 public:
  TestWorker(int* const concurrency_value, ThreadLock* const lock,
             int expected_value, bool sanitize)
      : concurrency_value_(concurrency_value),
        lock_(lock),
        expected_value_(expected_value),
        found_unexpected_value_(false),
        sanitize_(sanitize) {}

  WP2Status ExecuteSanitizeThread() {
    if (lock_ != nullptr) WP2_CHECK_STATUS(lock_->Acquire());
    *concurrency_value_ = expected_value_;
    found_unexpected_value_ = false;
    const double start = GetStopwatchTime();
    while ((GetStopwatchTime() - start) < kRunDuration) {
      found_unexpected_value_ |= (*concurrency_value_ != expected_value_);
      *concurrency_value_ = expected_value_;
    }
    if (lock_ != nullptr) lock_->Release();
    return WP2_STATUS_OK;
  }

  WP2_TSAN_IGNORE_FUNCTION
  WP2Status ExecuteNoSanitizeThread() {
    if (lock_ != nullptr) WP2_CHECK_STATUS(lock_->Acquire());
    *concurrency_value_ = expected_value_;
    found_unexpected_value_ = false;
    const double start = GetStopwatchTime();
    while ((GetStopwatchTime() - start) < kRunDuration) {
      found_unexpected_value_ |= (*concurrency_value_ != expected_value_);
      *concurrency_value_ = expected_value_;
    }
    if (lock_ != nullptr) lock_->Release();
    return WP2_STATUS_OK;
  }

  WP2Status Execute() override {
    if (sanitize_) return ExecuteSanitizeThread();
    return ExecuteNoSanitizeThread();
  }

 public:
  int* const concurrency_value_;
  ThreadLock* const lock_;
  const int expected_value_;
  bool found_unexpected_value_;
  const bool sanitize_;
};

class ThreadTest : public ::testing::TestWithParam<std::tuple<bool, bool>> {};

TEST_P(ThreadTest, MultiThread) {
  const bool threaded = std::get<0>(GetParam());
  const bool use_lock = std::get<1>(GetParam());

  int concurrency_value = 0;
  ThreadLock lock;
  const bool expect_the_unexpected = (threaded && !use_lock);
  const bool sanitize = !expect_the_unexpected;
  TestWorker first_worker(&concurrency_value, use_lock ? &lock : nullptr,
                          /*expected_value=*/1, sanitize);
  TestWorker second_worker(&concurrency_value, use_lock ? &lock : nullptr,
                           /*expected_value=*/-1, sanitize);

  EXPECT_WP2_OK(first_worker.Start(threaded));
  EXPECT_WP2_OK(second_worker.Start(threaded));
  EXPECT_WP2_OK(first_worker.End());
  EXPECT_WP2_OK(second_worker.End());

#ifdef WP2_USE_THREAD
  EXPECT_EQ(first_worker.found_unexpected_value_, expect_the_unexpected);
  EXPECT_EQ(second_worker.found_unexpected_value_, expect_the_unexpected);
#else
  EXPECT_FALSE(first_worker.found_unexpected_value_);
  EXPECT_FALSE(second_worker.found_unexpected_value_);
#endif  // WP2_USE_THREAD
}

INSTANTIATE_TEST_SUITE_P(ThreadTestInstantiation, ThreadTest,
                         ::testing::Combine(::testing::Bool() /* threaded */,
                                            ::testing::Bool() /* use_lock */));

}  // namespace
}  // namespace WP2
