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
// LocalContext and MainContext implementation.
//
// Author: Yannis Guyon (yguyon@google.com)

#include "./context_switch.h"

#if defined(WP2_USE_CONTEXT_SWITCH) && (WP2_USE_CONTEXT_SWITCH > 0)

#include <cstdlib>

#if defined(WP2_TRACE)
#include "src/utils/utils.h"
#define LOG_CONTEXT_ERROR(message)                                      \
  do {                                                                  \
    WP2ErrorLog(__FILE__, __LINE__, WP2_STATUS_OUT_OF_MEMORY, message); \
  } while (0)
#else
#define LOG_CONTEXT_ERROR(message)  do {} while (0)
#endif  // WP2_TRACE

namespace WP2 {

//------------------------------------------------------------------------------

void LocalContext::Reset() {
  running_ = false;
  opened_ = false;
  closed_ = true;
  should_close_ = false;
  inter_context_data_ = nullptr;
#ifdef _WIN32
  local_fiber_ = nullptr;
  main_fiber_ = nullptr;
#else
  local_context_.uc_link = nullptr;
  local_context_.uc_stack.ss_sp = nullptr;
  local_context_.uc_stack.ss_size = 0;
  main_context_.uc_link = nullptr;
  main_context_.uc_stack.ss_sp = nullptr;
  main_context_.uc_stack.ss_size = 0;
#endif
}

bool LocalContext::Yield() {
  if (should_close_ || error_) return false;

  if (!running_ || !opened_ || closed_) {
    LOG_CONTEXT_ERROR("Context state is invalid.");
    error_ = true;
    return false;
  }

#ifdef _WIN32
  SwitchToFiber(main_fiber_);
#else
  if (swapcontext(&local_context_, &main_context_) != 0) {
    LOG_CONTEXT_ERROR("Can not swap context.");
    error_ = true;
    return false;
  }
#endif

  if (!running_ || !opened_ || closed_) {
    LOG_CONTEXT_ERROR("Context state is invalid.");
    error_ = true;
  }
  return !error_ && !should_close_;
}

void LocalContext::Close() {
  if (!running_ || !opened_ || closed_) {
    LOG_CONTEXT_ERROR("Context state is invalid.");
    error_ = true;
  }

  closed_ = true;

#ifdef _WIN32
  SwitchToFiber(main_fiber_);
#else
  if (swapcontext(&local_context_, &main_context_) != 0) {
    LOG_CONTEXT_ERROR("Can not Close() from context.");
    return;  // Let's see if just returning works.
  }
#endif

  LOG_CONTEXT_ERROR("Yield() called after Close().");
  error_ = true;
}

//------------------------------------------------------------------------------

bool MainContext::CreateLocalContext(
    void (*suspendable_function)(LocalContext*), void* inter_context_data) {
  if (!context_.closed_) {
    LOG_CONTEXT_ERROR("Another context is already created.");
    return false;
  } else if (suspendable_function == nullptr) {
    LOG_CONTEXT_ERROR("Parameter is invalid.");
    return false;
  }

  context_.Reset();
  context_.closed_ = false;
  context_.error_ = false;
  context_.inter_context_data_ = inter_context_data;

#ifdef _WIN32
  context_.main_fiber_ = ConvertThreadToFiber(nullptr);
  if (context_.main_fiber_ == nullptr) {
    LOG_CONTEXT_ERROR("Can not allocate calling fiber.");
    context_.Reset();
    context_.error_ = true;
    return false;
  }
  context_.local_fiber_ = CreateFiber(
      (SIZE_T)LocalContext::kCallstackSize,
      (LPFIBER_START_ROUTINE)suspendable_function, (LPVOID)&context_);
  if (context_.local_fiber_ == nullptr) {
    LOG_CONTEXT_ERROR("Can not allocate called fiber.");
    context_.Reset();
    context_.error_ = true;
    return false;
  }
#else
  if (getcontext(&context_.local_context_) != 0) {
    LOG_CONTEXT_ERROR("Can not allocate called context.");
    context_.Reset();
    context_.error_ = true;
    return false;
  }
  context_.local_context_.uc_link = nullptr;
  context_.local_context_.uc_stack.ss_sp = context_.stack_buffer_;
  context_.local_context_.uc_stack.ss_size =
      sizeof(context_.stack_buffer_) / sizeof(context_.stack_buffer_[0]);
  context_.local_context_.uc_stack.ss_flags = 0;
  makecontext(&context_.local_context_,
              reinterpret_cast<void (*)()>(suspendable_function), 1, &context_);
#endif
  return true;
}

bool MainContext::Resume() {
  if (context_.closed_) {
    // This can be used in a loop, no need to print an error.
    CloseLocalContext();  // Clean up now, in case it's not done later on.
    return false;
  } else if (context_.should_close_) {
    LOG_CONTEXT_ERROR("Resume() called after Close().");
    return false;
  } else if (context_.running_) {
    LOG_CONTEXT_ERROR("Resume() called but already running.");
    return false;
  }
  if (context_.error_) return false;

  context_.running_ = true;

#ifdef _WIN32
  context_.opened_ = true;
  SwitchToFiber(context_.local_fiber_);
#else
  bool first_time_distant_context_entered = !context_.opened_;
  context_.opened_ = true;
  if (swapcontext(&context_.main_context_, &context_.local_context_) != 0) {
    context_.running_ = false;
    if (first_time_distant_context_entered) context_.opened_ = false;
    LOG_CONTEXT_ERROR("Can not swap context.");
    context_.error_ = true;
    return false;
  }
#endif

  context_.running_ = false;
  if (context_.closed_) {
#ifdef _WIN32
    DeleteFiber(context_.local_fiber_);
#endif
    context_.Reset();
  }
  return true;
}

void MainContext::CloseLocalContext() {
  if (!context_.opened_ || context_.closed_) {
#ifdef _WIN32
    if (context_.local_fiber_ != nullptr) DeleteFiber(context_.local_fiber_);
#endif
    context_.Reset();
    return;
  }
  context_.should_close_ = true;  // Skip all Yield() and aim for Close().
  context_.running_ = true;

#ifdef _WIN32
  SwitchToFiber(context_.local_fiber_);
  DeleteFiber(context_.local_fiber_);
#else
  if (swapcontext(&context_.main_context_, &context_.local_context_) != 0) {
    LOG_CONTEXT_ERROR("Can not swap context.");
    context_.running_ = false;
    context_.error_ = true;
  }
#endif
  context_.Reset();
}

//------------------------------------------------------------------------------

}  // namespace WP2

#endif  // WP2_USE_CONTEXT_SWITCH
