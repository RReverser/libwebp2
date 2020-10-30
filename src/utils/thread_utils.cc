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
// Multi-threaded worker
//
// Author: Skal (pascal.massimino@gmail.com)
//         James Zern (jzern@google.com)

#include "src/utils/thread_utils.h"

#if defined(WP2_USE_THREAD)
#if defined(_WIN32)
#include <windows.h>
#include <process.h>
#else
#include <pthread.h>
#endif
#endif

#include <cassert>

#include "src/utils/utils.h"

namespace WP2 {

#ifdef WP2_USE_THREAD

#if defined(_WIN32)

typedef HANDLE pthread_t;
typedef CRITICAL_SECTION pthread_mutex_t;
typedef struct {
  uint32_t stack_size;
} pthread_attr_t;

#if _WIN32_WINNT >= 0x0600  // Windows Vista / Server 2008 or greater
#define USE_WINDOWS_CONDITION_VARIABLE
typedef CONDITION_VARIABLE pthread_cond_t;
#else
typedef struct {
  HANDLE waiting_sem_;
  HANDLE received_sem_;
  HANDLE signal_event_;
} pthread_cond_t;
#endif  // _WIN32_WINNT >= 0x600

#ifndef WINAPI_FAMILY_PARTITION
#define WINAPI_PARTITION_DESKTOP 1
#define WINAPI_FAMILY_PARTITION(x) x
#endif

#if !WINAPI_FAMILY_PARTITION(WINAPI_PARTITION_DESKTOP)
#define USE_CREATE_THREAD
#endif

#endif  // _WIN32

typedef struct {
  pthread_mutex_t mutex_;
  pthread_cond_t  condition_;
  pthread_t       thread_;
} WorkerImpl;

#if defined(_WIN32)

//------------------------------------------------------------------------------
// simplistic pthread emulation layer

// _beginthreadex requires __stdcall
#define THREADFN unsigned int __stdcall
#define THREAD_RETURN(val) (unsigned int)((DWORD_PTR)val)

#if _WIN32_WINNT >= 0x0501  // Windows XP or greater
#define WaitForSingleObject(obj, timeout) \
  WaitForSingleObjectEx(obj, timeout, FALSE /*bAlertable*/)
#endif

static int pthread_attr_init(pthread_attr_t* attr) {
  if (attr != NULL) attr->stack_size = 0;
  return 0;
}
static int pthread_attr_destroy(pthread_attr_t* attr) { return 0; }
static int pthread_attr_setstacksize(pthread_attr_t* attr, size_t stack_size) {
  if (attr != NULL) attr->stack_size = stack_size;
  return 1;
}

static int pthread_create(pthread_t* const thread, const pthread_attr_t* attr,
                          unsigned int (__stdcall *start)(void*), void* arg) {
  const size_t stack_size = (attr != NULL) ? attr->stack_size : 0u;
#ifdef USE_CREATE_THREAD
  *thread = CreateThread(NULL,        /* lpThreadAttributes */
                         stack_size,  /* dwStackSize */
                         start,
                         arg,
                         0,           /* dwCreationFlags */
                         NULL);       /* lpThreadId */
#else
  *thread = (pthread_t)_beginthreadex(NULL,        /* void *security */
                                      stack_size,  /* unsigned stack_size */
                                      start,
                                      arg,
                                      0,           /* unsigned initflag */
                                      NULL);       /* unsigned *thrdaddr */
#endif
  if (*thread == NULL) return 1;
  SetThreadPriority(*thread, THREAD_PRIORITY_ABOVE_NORMAL);
  return 0;
}

static int pthread_join(pthread_t thread, void** value_ptr) {
  (void)value_ptr;
  return (WaitForSingleObject(thread, INFINITE) != WAIT_OBJECT_0 ||
          CloseHandle(thread) == 0);
}

// Mutex
static int pthread_mutex_init(pthread_mutex_t* const mutex, void* mutexattr) {
  (void)mutexattr;
#if _WIN32_WINNT >= 0x0600  // Windows Vista / Server 2008 or greater
  InitializeCriticalSectionEx(mutex, 0 /*dwSpinCount*/, 0 /*Flags*/);
#else
  InitializeCriticalSection(mutex);
#endif
  return 0;
}

static int pthread_mutex_lock(pthread_mutex_t* const mutex) {
  EnterCriticalSection(mutex);
  return 0;
}

static int pthread_mutex_unlock(pthread_mutex_t* const mutex) {
  LeaveCriticalSection(mutex);
  return 0;
}

static int pthread_mutex_destroy(pthread_mutex_t* const mutex) {
  DeleteCriticalSection(mutex);
  return 0;
}

// Condition
static int pthread_cond_destroy(pthread_cond_t* const condition) {
  int ok = 1;
#ifdef USE_WINDOWS_CONDITION_VARIABLE
  (void)condition;
#else
  ok &= (CloseHandle(condition->waiting_sem_) != 0);
  ok &= (CloseHandle(condition->received_sem_) != 0);
  ok &= (CloseHandle(condition->signal_event_) != 0);
#endif
  return !ok;
}

static int pthread_cond_init(pthread_cond_t* const condition, void* cond_attr) {
  (void)cond_attr;
#ifdef USE_WINDOWS_CONDITION_VARIABLE
  InitializeConditionVariable(condition);
#else
  condition->waiting_sem_ = CreateSemaphore(NULL, 0, 1, NULL);
  condition->received_sem_ = CreateSemaphore(NULL, 0, 1, NULL);
  condition->signal_event_ = CreateEvent(NULL, FALSE, FALSE, NULL);
  if (condition->waiting_sem_ == NULL ||
      condition->received_sem_ == NULL ||
      condition->signal_event_ == NULL) {
    pthread_cond_destroy(condition);
    return 1;
  }
#endif
  return 0;
}

static int pthread_cond_signal(pthread_cond_t* const condition) {
  int ok = 1;
#ifdef USE_WINDOWS_CONDITION_VARIABLE
  WakeConditionVariable(condition);
#else
  if (WaitForSingleObject(condition->waiting_sem_, 0) == WAIT_OBJECT_0) {
    // a thread is waiting in pthread_cond_wait: allow it to be notified
    ok = SetEvent(condition->signal_event_);
    // wait until the event is consumed so the signaler cannot consume
    // the event via its own pthread_cond_wait.
    ok &= (WaitForSingleObject(condition->received_sem_, INFINITE) !=
           WAIT_OBJECT_0);
  }
#endif
  return !ok;
}

static int pthread_cond_wait(pthread_cond_t* const condition,
                             pthread_mutex_t* const mutex) {
  int ok;
#ifdef USE_WINDOWS_CONDITION_VARIABLE
  ok = SleepConditionVariableCS(condition, mutex, INFINITE);
#else
  // note that there is a consumer available so the signal isn't dropped in
  // pthread_cond_signal
  if (!ReleaseSemaphore(condition->waiting_sem_, 1, NULL)) return 1;
  // now unlock the mutex so pthread_cond_signal may be issued
  pthread_mutex_unlock(mutex);
  ok = (WaitForSingleObject(condition->signal_event_, INFINITE) ==
        WAIT_OBJECT_0);
  ok &= ReleaseSemaphore(condition->received_sem_, 1, NULL);
  pthread_mutex_lock(mutex);
#endif
  return !ok;
}

#else  // !_WIN32
# define THREADFN void*
# define THREAD_RETURN(val) val
#endif  // _WIN32

//------------------------------------------------------------------------------

static THREADFN ThreadLoop(void* ptr) {
  Worker* const worker = (Worker*)ptr;
  assert(worker != nullptr);
  WorkerImpl* const impl = (WorkerImpl*)worker->GetImpl();
  assert(impl != nullptr);
  bool done = false;
  while (!done) {
    pthread_mutex_lock(&impl->mutex_);
    while (worker->state_ == Worker::OK) {   // wait in idling mode
      pthread_cond_wait(&impl->condition_, &impl->mutex_);
    }
    if (worker->state_ == Worker::WORK) {
      worker->status_ = worker->Execute();
      worker->state_ = Worker::OK;
    } else if (worker->state_ == Worker::NOT_OK) {   // finish the worker
      done = true;
    }
    // signal to the main thread that we're done (for Sync())
    pthread_mutex_unlock(&impl->mutex_);
    pthread_cond_signal(&impl->condition_);
  }
  return THREAD_RETURN(NULL);    // Thread is finished
}

// main thread state control
static void ChangeState(WorkerBase* const worker,
                        decltype(WorkerBase::state_) new_state) {
  // No-op when attempting to change state on a thread that didn't come up.
  // Checking state_ without acquiring the lock first would result in a data
  // race.
  assert(worker != nullptr);
  WorkerImpl* const impl = (WorkerImpl*)worker->GetImpl();
  if (impl == nullptr) return;

  pthread_mutex_lock(&impl->mutex_);
  if (worker->state_ >= Worker::OK) {
    // wait for the worker to finish
    while (worker->state_ != Worker::OK) {
      pthread_cond_wait(&impl->condition_, &impl->mutex_);
    }
    // assign new status and release the working thread if needed
    if (new_state != Worker::OK) {
      worker->state_ = new_state;
      // Note the associated mutex does not need to be held when signaling the
      // condition. Unlocking the mutex first may improve performance in some
      // implementations, avoiding the case where the waiting thread can't
      // reacquire the mutex when woken.
      pthread_mutex_unlock(&impl->mutex_);
      pthread_cond_signal(&impl->condition_);
      return;
    }
  }
  pthread_mutex_unlock(&impl->mutex_);
}

#endif  // WP2_USE_THREAD

//------------------------------------------------------------------------------

struct ThreadLock::Impl : public WP2Allocable {
#ifdef WP2_USE_THREAD
  pthread_mutex_t mutex;
#endif
};

ThreadLock::ThreadLock() {
#ifdef WP2_USE_THREAD
  // If new Impl() fails, Acquire() will return WP2_STATUS_OUT_OF_MEMORY.
  impl_ = new (WP2Allocable::nothrow) Impl();
  if ((impl_ != nullptr) && (pthread_mutex_init(&impl_->mutex, nullptr) != 0)) {
    delete impl_;
    impl_ = nullptr;
  }
#else
  impl_ = nullptr;
#endif
}

ThreadLock::~ThreadLock() {
#ifdef WP2_USE_THREAD
  if (impl_ != nullptr) {
    pthread_mutex_destroy(&impl_->mutex);
    delete impl_;
  }
#endif
}

WP2Status ThreadLock::Acquire() {
#ifdef WP2_USE_THREAD
  WP2_CHECK_ALLOC_OK(impl_ != nullptr);
  WP2_CHECK_OK(pthread_mutex_lock(&impl_->mutex) == 0,
               WP2_STATUS_OUT_OF_MEMORY);
#endif
  return WP2_STATUS_OK;
}

void ThreadLock::Release() {
#ifdef WP2_USE_THREAD
  if (impl_ != nullptr) pthread_mutex_unlock(&impl_->mutex);
#endif
}

//------------------------------------------------------------------------------

WP2Status WorkerBase::Start(bool do_mt, uint32_t stack_size) {
  const auto& I = GetWorkerInterface();
  WP2_CHECK_OK(I.Reset(this, do_mt, stack_size), WP2_STATUS_OUT_OF_MEMORY);
  I.Launch(this);
  return WP2_STATUS_OK;
}

WP2Status WorkerBase::End() {
  const auto& I = GetWorkerInterface();
  const WP2Status status = I.Sync(this);
  I.End(this);
  WP2_CHECK_STATUS(status);
  return WP2_STATUS_OK;
}

//------------------------------------------------------------------------------

static class : public WorkerInterface {
  WP2Status Sync(WorkerBase* const worker) const override {
    if (worker == nullptr) return WP2_STATUS_NULL_PARAMETER;
#ifdef WP2_USE_THREAD
    ChangeState(worker, Worker::OK);
#endif
    assert(worker->state_ <= Worker::OK);
    return worker->status_;
  }

  bool Reset(WorkerBase* const worker,
             bool do_mt, uint32_t stack_size) const override {
    if (worker == nullptr) return false;
    bool ok = true;
    worker->status_ = WP2_STATUS_OK;
    if (worker->state_ < Worker::OK) {
      worker->state_ = Worker::OK;
#ifdef WP2_USE_THREAD
      if (do_mt) {
        WorkerImpl* const impl = (WorkerImpl*)WP2Calloc(1, sizeof(WorkerImpl));
        worker->SetImpl(impl);
        if (impl == nullptr) return false;
        if (pthread_mutex_init(&impl->mutex_, NULL)) goto Error;
        if (pthread_cond_init(&impl->condition_, NULL)) {
          pthread_mutex_destroy(&impl->mutex_);
          goto Error;
        }
        pthread_attr_t attr;
        ok = (pthread_attr_init(&attr) == 0);
        if (ok && stack_size > 0) {
          ok = (pthread_attr_setstacksize(&attr, stack_size) == 0);
        }
        if (ok) {
          pthread_mutex_lock(&impl->mutex_);
          ok = (pthread_create(&impl->thread_, &attr, ThreadLoop, worker) == 0);
          if (ok) worker->state_ = Worker::OK;
          pthread_mutex_unlock(&impl->mutex_);
        }
        if (pthread_attr_destroy(&attr) != 0) ok = false;
        if (!ok) {
          pthread_mutex_destroy(&impl->mutex_);
          pthread_cond_destroy(&impl->condition_);
   Error:
          WP2Free(impl);
          worker->SetImpl(nullptr);
          worker->state_ = Worker::NOT_OK;
          return false;
        }
      }
#endif
    } else if (worker->state_ == Worker::WORK) {
      ok = (Sync(worker) == WP2_STATUS_OK);
    }
    assert(!ok || worker->state_ == Worker::OK);
    return ok;
  }

  void Launch(WorkerBase* const worker) const override {
    if (worker == nullptr) return;   // !?
    if (worker->state_ == Worker::OK) {
#ifdef WP2_USE_THREAD
      if (worker->GetImpl() != nullptr) {
        ChangeState(worker, Worker::WORK);
        return;
      }
      // fall through
#endif
      worker->status_ = worker->Execute();   // direct call
    }
  }

  void End(WorkerBase* const worker) const override {
    if (worker == nullptr) return;   // !?
#ifdef WP2_USE_THREAD
    WorkerImpl* const impl = (WorkerImpl*)worker->GetImpl();
    if (impl != nullptr) {
      ChangeState(worker, Worker::NOT_OK);
      pthread_join(impl->thread_, NULL);
      pthread_mutex_destroy(&impl->mutex_);
      pthread_cond_destroy(&impl->condition_);
      WP2Free(impl);
      worker->SetImpl(nullptr);
      assert(worker->state_ == Worker::NOT_OK);
      return;
    }
#endif
    worker->state_ = Worker::NOT_OK;
    assert(worker->GetImpl() == nullptr);
  }
} kDefaultInterface;

//------------------------------------------------------------------------------

static const WorkerInterface* g_worker_interface = &kDefaultInterface;

void SetWorkerInterface(const WorkerInterface& winterface) {
  g_worker_interface = &winterface;
}

const WorkerInterface& GetWorkerInterface() {
  return *g_worker_interface;
}

//------------------------------------------------------------------------------

}    // namespace WP2
