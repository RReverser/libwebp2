// Copyright 2020 Google LLC
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
// implementation of the 1990 CCITT Group4 spec encoder.
// https://en.wikipedia.org/wiki/Group_4_compression
//
// Author: Vincent Rabaud (vrabaud@google.com)

#include "src/common/lossless/color_cache.h"
#include "src/enc/lossless/losslessi_enc.h"
#include "src/enc/symbols_enc.h"

namespace WP2L {

static uint16_t GetColor(const uint16_t* const argb, int x,
                         uint16_t fallback = 0) {
  if (x < 0) return fallback;
  // Get the color from the green channel.
  return argb[x * 4 + 2];
}

static void WriteColor(uint16_t current_color, uint16_t expected_new_color,
                       uint16_t actual_new_color, MoveToFrontCache* const mtf,
                       WP2::ANSEncBase* const enc,
                       WP2::SymbolManager* const sm) {
  assert(actual_new_color != current_color);
  assert(expected_new_color != current_color);
  const bool color_change = (actual_new_color != expected_new_color);
  sm->Process(/*cluster=*/0, kSymbolG4ColorChange, color_change, "color_change",
              enc);

  const uint32_t current_color_idx = mtf->GetIndex(current_color);
  const uint32_t expected_new_color_idx = mtf->GetIndex(expected_new_color);
  const uint32_t actual_new_color_idx = mtf->GetIndex(actual_new_color);

  mtf->MoveToFront(actual_new_color);

  if (!color_change) return;

  uint16_t coded_color = actual_new_color_idx;
  if (actual_new_color_idx > current_color_idx) --coded_color;
  if (actual_new_color_idx > expected_new_color_idx) --coded_color;
  // TODO(maryla): optionally use a color cache? or a 'move-to-front' queue for
  //  most used colors.
  sm->Process(/*cluster=*/0, kSymbolG4NewColor, coded_color, "color_change",
              enc);
}

// Implementation taken from https://www.itu.int/rec/T-REC-T.6-198811-I/en.
// See the above for diagrams showing the meaning of a0/a1/b1/b2.
static void Impl(const uint16_t* const argb, uint32_t width, uint32_t height,
                 MoveToFrontCache* const mtf, WP2::ANSEncBase* const enc,
                 WP2::SymbolManager* const sm, EncodeInfo* const encode_info) {
  WP2::ANSDebugPrefix prefix(enc, "group4");

  const uint32_t num_colors = mtf->num_colors();
  const uint16_t first_color = GetColor(argb, /*x=*/0);

  for (uint32_t y = 0; y < height; ++y) {
    const uint16_t* prev = (y > 0) ? &argb[(y - 1) * width * 4] : nullptr;
    const uint16_t* cur = &argb[y * width * 4];

    mtf->MoveToFront(first_color);

    // Start from the virtual pixel before the first column.
    int a0_prev = 0;
    for (int a0 = -1; a0 < (int)width - 1;) {
      // All rows are assumed to start with first_color. We could also use the
      // first color of the row above, but it's not actually better.
      const uint16_t a0_color = GetColor(cur, a0, first_color);

      // Find a1: next pixel changing color on the current line.
      int a1 = a0 + 1;
      while (a1 < (int)width && GetColor(cur, a1) == a0_color) ++a1;
      const uint16_t a1_color = a1 < (int)width ? GetColor(cur, a1) : 0;

      // Find b1: next changing pixel of color different from a0_color
      // on the previous line.
      int b1 = a0 + 1;
      if (y == 0) {
        b1 = width;
      } else {
        while (b1 < (int)width &&
               !(GetColor(prev, b1) != GetColor(prev, b1 - 1, first_color) &&
                 GetColor(prev, b1) != a0_color)) {
          ++b1;
        }
      }
      const uint64_t b1_color =
          (b1 < (int)width) ? GetColor(prev, b1, !a0_color) : !a0_color;
      // Find b2: next changing pixel after b1 on the previous line.
      int b2 = b1 + 1;
      if (y == 0) {
        b2 = width;
      } else {
        while (b2 < (int)width && GetColor(prev, b2) == GetColor(prev, b1)) {
          ++b2;
        }
      }
      const uint16_t b2_color =
          (b2 < (int)width) ? GetColor(prev, b2, a0_color) : a0_color;

      // Recent colors on the row above are likely to be reused.
      mtf->MoveToFront(b1_color);
      mtf->MoveToFront(b2_color);

      const uint32_t a0b1_dist = b1 - a0;
      const uint32_t a0a1_max = b2 - a0;
      const bool horizontal_mode_allowed =
          (a0b1_dist > kGroup4Window + 1) ||
          (a0b1_dist + kGroup4Window < a0a1_max);
      const bool pass_mode_allowed = (b2 + 1) <= (int)width;
      const uint32_t max_mode = horizontal_mode_allowed + pass_mode_allowed;

      if (b2 < a1) {
        assert(pass_mode_allowed);
        // Pass mode.
        sm->Process(/*cluster=*/0, kSymbolG4Type,
                    std::min((uint32_t)Group4Mode::kPass, max_mode), max_mode,
                    "mode", enc);
        a0 = b2;
      } else {
        int dist = a1 - b1;
        if ((uint32_t)std::abs(dist) <= kGroup4Window) {
          // Vertical mode.
          if (max_mode > 0) {
            sm->Process(/*cluster=*/0, kSymbolG4Type,
                        (uint32_t)Group4Mode::kVertical, max_mode, "mode", enc);
          }
          // Restricting the range as below is actually worse.
          const int min_dist = -(int)kGroup4Window;
          const int max_dist = (int)std::min(b1 + kGroup4Window, width) - b1;

          sm->Process(/*cluster=*/0, kSymbolG4VerticalDist, dist - min_dist,
                      max_dist - min_dist, "vdist", enc);

          if (a1 < (int)width && num_colors > 2) {
            WriteColor(a0_color, /*expected_new_color=*/b1_color,
                       /*actual_new_color=*/a1_color, mtf, enc, sm);
          }

          a0 = a1;
        } else {
          // Horizontal mode.
          assert(horizontal_mode_allowed);
          sm->Process(/*cluster=*/0, kSymbolG4Type,
                      std::min((uint32_t)Group4Mode::kHorizontal, max_mode),
                      max_mode, "mode", enc);
          const char* const str[2] = {"hdist_bg_minus1", "hdist_fg_minus1"};
          // Note1: some pixel locations are impossible since they would be
          // coded in vertical mode, but removing those values by shifting all
          // values above does not give better results.
          // Note2: we use two clusters: one for the first_color, one for all
          // the other colors. Having more clusters (e.g. 3 or 4 clusters for
          // images with 3 or 4 colors) does not seem to help.
          // TODO(maryla): around 14% of images do benefit from having only one
          // cluster.
          sm->Process(
              /*cluster=*/(uint32_t)(a0_color != first_color),
              kSymbolG4HorizontalDist, a1 - a0 - 1, a0a1_max - 1,
              str[a0_color != first_color], enc);
          const int a1a2_max = (int)width - a1;
          int a2 = a1 + 1;
          if (a1a2_max > 0) {
            while (a2 < (int)width && GetColor(cur, a2) == a1_color) {
              ++a2;
            }
            sm->Process(
                /*cluster=*/(uint32_t)(a0_color == first_color),
                kSymbolG4HorizontalDist, a2 - a1 - 1, a1a2_max - 1,
                str[a0_color == first_color], enc);

            if (num_colors > 2) {
              WriteColor(a0_color, /*expected_new_color=*/b1_color,
                         /*actual_new_color=*/a1_color, mtf, enc, sm);

              if (a2 < (int)width) {
                const uint16_t expected_a2_color =
                    (b2_color != a1_color) ? b2_color : a0_color;
                const uint16_t a2_color = GetColor(cur, a2);
                WriteColor(a1_color, expected_a2_color,
                           /*actual_new_color=*/a2_color, mtf, enc, sm);
              }
            }
          }
          a0 = a2;
        }
      }
      // Update the costs.
      if (encode_info != nullptr && !encode_info->bits_per_pixel.empty()) {
        for (uint32_t x = a0_prev; x < std::min((uint32_t)a0, width); ++x) {
          // Set an arbitrary cost of 5.
          encode_info->bits_per_pixel[y * width + x] = 5. / (a0 - a0_prev);
        }
      }
      a0_prev = a0;
    }
    if (encode_info != nullptr && !encode_info->line_tokens.empty()) {
      encode_info->line_tokens[y] = enc->NumTokens();
    }
    // Update the costs.
    if (encode_info != nullptr && !encode_info->bits_per_pixel.empty() &&
        a0_prev < (int)width) {
      // All leftover pixels are deduced.
      std::fill(&encode_info->bits_per_pixel[0] + y * width + a0_prev,
                &encode_info->bits_per_pixel[0] + y * width + width, 0);
    }
  }
}

WP2Status Group4Encode(const uint16_t* const argb, uint32_t width,
                       uint32_t height, uint32_t num_colors,
                       bool use_move_to_front, WP2::ANSEncBase* const enc,
                       EncodeInfo* const encode_info) {
  WP2::SymbolsInfo info;
  InitGroup4(width, num_colors, &info);

  MoveToFrontCache mtf;
  WP2_CHECK_STATUS(mtf.Init(use_move_to_front, num_colors));

  // Get symbol statistics.
  WP2::SymbolRecorder sr;
  WP2_CHECK_STATUS(sr.Allocate(info, /*num_records=*/0));
  WP2::ANSEncNoop noop;
  Impl(argb, width, height, &mtf, &noop, &sr, /*encode_info=*/nullptr);

  // Write the headers.
  WP2::SymbolWriter sw;
  WP2_CHECK_STATUS(sw.Init(info));
  WP2_CHECK_STATUS(sw.Allocate());

  enc->AddDebugPrefix("GlobalHeader");
  enc->AddDebugPrefix("transforms");
  WP2::ANSDictionaries dicts;
  const uint32_t max_nnz = width * height;

  {
    WP2::ANSDebugPrefix prefix(enc, "group4");
    for (uint32_t s = 0; s < (uint32_t)kSymbolG4Num; ++s) {
      if (num_colors <= 2 &&
          (s == kSymbolG4ColorChange || s == kSymbolG4NewColor)) {
        continue;
      }
      WP2_CHECK_STATUS(
          sw.WriteHeader(max_nnz, s, sr, kSymbolGroup4Names[s], enc, &dicts));
    }
  }
  enc->PopDebugPrefix();
  enc->PopDebugPrefix();

  // Write the data.
  mtf.Reset();
  Impl(argb, width, height, &mtf, enc, &sw, encode_info);

  return WP2_STATUS_OK;
}

//------------------------------------------------------------------------------

}  // namespace WP2L
