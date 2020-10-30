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
// -----------------------------------------------------------------------------

// Context test.

#include "src/utils/context_switch.h"

#if defined(WP2_USE_CONTEXT_SWITCH) && (WP2_USE_CONTEXT_SWITCH > 0)

#include "gtest/gtest.h"

namespace WP2 {
namespace {

//------------------------------------------------------------------------------

struct InterContextData {
 public:
  int data = 0;
};

void IncrementData(LocalContext* const context,
                   InterContextData* const inter_context_data) {
  inter_context_data->data += 1;
  if (!context->Yield()) return;
  inter_context_data->data += 2;
  if (!context->Yield()) return;
  inter_context_data->data += 4;
  if (!context->Yield()) return;
  inter_context_data->data += 8;
  if (!context->Yield()) return;
  inter_context_data->data += 16;
}

void RunContext(LocalContext* const context) {
  IncrementData(context, (InterContextData*)context->GetInterContextData());
  context->Close();
}

void RunBadContext(LocalContext* const context) {
  context->Close();
  context->Yield();  // Nothing should be called after Close().
  context->Close();
}

//------------------------------------------------------------------------------

TEST(ContextTest, Simple) {
  InterContextData inter_context_data;
  MainContext main_context;

  ASSERT_TRUE(
      main_context.CreateLocalContext(&RunContext, &inter_context_data));

  ASSERT_TRUE(main_context.Resume());
  EXPECT_EQ(inter_context_data.data, 1);
  ASSERT_TRUE(main_context.Resume());
  EXPECT_EQ(inter_context_data.data, 1 + 2);
  ASSERT_TRUE(main_context.Resume());
  EXPECT_EQ(inter_context_data.data, 1 + 2 + 4);
  ASSERT_TRUE(main_context.Resume());
  main_context.CloseLocalContext();

  EXPECT_EQ(inter_context_data.data, 1 + 2 + 4 + 8);
}

//------------------------------------------------------------------------------

TEST(ContextTest, Loop) {
  InterContextData inter_context_data;
  MainContext main_context;
  ASSERT_TRUE(main_context.Ok());
  ASSERT_TRUE(
      main_context.CreateLocalContext(&RunContext, &inter_context_data));

  while (main_context.Resume()) {
    ASSERT_TRUE(main_context.Ok());
  }
  main_context.CloseLocalContext();

  ASSERT_TRUE(main_context.Ok());
  EXPECT_EQ(inter_context_data.data, 1 + 2 + 4 + 8 + 16);
}

//------------------------------------------------------------------------------

TEST(ContextTest, ExpectErrors) {
  MainContext main_context;

  ASSERT_FALSE(main_context.CreateLocalContext(nullptr, nullptr));
  ASSERT_FALSE(main_context.Resume());

  {
    InterContextData inter_context_data;
    ASSERT_TRUE(
        main_context.CreateLocalContext(&RunContext, &inter_context_data));
    main_context.CloseLocalContext();
    ASSERT_FALSE(main_context.Resume());
  }

  {
    InterContextData inter_context_data;
    ASSERT_TRUE(
        main_context.CreateLocalContext(&RunContext, &inter_context_data));
    ASSERT_TRUE(main_context.Resume());
    main_context.CloseLocalContext();
    ASSERT_FALSE(main_context.Resume());
  }

  {
    InterContextData inter_context_data;
    ASSERT_TRUE(
        main_context.CreateLocalContext(&RunContext, &inter_context_data));
    while (inter_context_data.data != (1 + 2 + 4 + 8 + 16)) {
      ASSERT_TRUE(main_context.Resume());
    }
    ASSERT_FALSE(main_context.Resume());
  }
}

//------------------------------------------------------------------------------

TEST(ContextTest, BadFunction) {
  {
    MainContext main_context;
    ASSERT_TRUE(main_context.CreateLocalContext(&RunBadContext, nullptr));
    ASSERT_TRUE(main_context.Resume());
    ASSERT_FALSE(main_context.Resume());
    main_context.CloseLocalContext();
  }
}

//------------------------------------------------------------------------------

}  // namespace
}  // namespace WP2

#endif  // WP2_USE_CONTEXT_SWITCH
