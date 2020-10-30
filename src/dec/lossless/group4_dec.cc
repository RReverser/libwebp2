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
#include "src/common/symbols.h"
#include "src/dec/lossless/losslessi_dec.h"
#include "src/dec/wp2_dec_i.h"
#include "src/utils/plane.h"

namespace WP2L {


static void SetColor(int16_t color, uint32_t x, int16_t* const dst) {
  // Write to the image to the green channel (channel 2).
  dst[4 * x + 0] = dst[4 * x + 1] = dst[4 * x + 3] = 0;
  dst[4 * x + 2] = color;
}

static uint16_t GetColor(const int16_t* const argb, int x,
                         int16_t fallback = 0) {
  if (x < 0) return fallback;
  // Get the color from the green channel.
  return argb[x * 4 + 2];
}

static uint16_t ChangeColor(uint32_t num_colors, uint16_t current_color,
                            uint16_t expected_new_color,
                            MoveToFrontCache* const mtf,
                            WP2::SymbolReader* const sr) {
  assert(current_color != expected_new_color);
  if (num_colors == 2) return !current_color;
  const bool color_change = sr->Read(kSymbolG4ColorChange, "color_change");
  if (!color_change) {
    mtf->MoveToFront(expected_new_color);
    return expected_new_color;
  }
  uint32_t new_color_idx = sr->Read(kSymbolG4NewColor, "color_change");
  const uint32_t current_color_idx = mtf->GetIndex(current_color);
  const uint32_t expected_new_color_idx = mtf->GetIndex(expected_new_color);
  if (new_color_idx >= std::min(expected_new_color_idx, current_color_idx)) {
    ++new_color_idx;
  }
  if (new_color_idx >= std::max(expected_new_color_idx, current_color_idx)) {
    ++new_color_idx;
  }
  const uint32_t new_color = mtf->GetColor(new_color_idx);
  mtf->MoveToFront(new_color);
  return new_color;
}

// Implementation taken from https://www.itu.int/rec/T-REC-T.6-198811-I/en.
// See the above for diagrams showing the meaning of a0/a1/b1/b2.
WP2Status Decoder::DecodeGroup4(uint32_t width, uint32_t last_row,
                                WP2::Planef* const bits_per_pixel,
                                int16_t** src_out) {
  assert(num_colors_ >= 2);
  WP2::ANSDebugPrefix prefix(dec_, "group4");

  WP2::SymbolReader* const sr = &hdr_.sr_;
  int16_t* const data = pixels_.data();
  int16_t* src = data + 4 * last_pixel_;

  const int16_t first_color = group4_first_color_;

  const bool get_cost = (bits_per_pixel != nullptr || config_.info != nullptr);
  for (uint32_t y = last_pixel_ / width; y < last_row; ++y) {
    const int16_t* prev = (y > 0) ? &data[(y - 1) * width * 4] : nullptr;

    // Define the virtual pixel before the first column.
    int16_t current_color = first_color;
    mtf_.MoveToFront(current_color);

    int a0 = -1;
    uint32_t x = 0;
    while (x < width) {
      double cost = 0.;
      const uint32_t start_x = x;
      double* cost_ptr = get_cost ? &cost : nullptr;

      const int16_t a0_color = current_color;

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
          (b1 < (int)width) ? GetColor(prev, b1) : !a0_color;

      // Find b2: next changing pixel after b1 on the previous line.
      uint32_t b2 = b1 + 1;
      if (y == 0) {
        b2 = width;
      } else {
        while (b2 < width && GetColor(prev, b2) == GetColor(prev, b1)) ++b2;
      }
      const uint16_t b2_color =
          (b2 < width) ? GetColor(prev, b2, a0_color) : a0_color;

      mtf_.MoveToFront(b1_color);
      mtf_.MoveToFront(b2_color);

      const uint32_t a0b1_dist = b1 - a0;
      const uint32_t a0a1_max = b2 - a0;
      // Horizontal mode allowed if all the possible pixel locations aren't
      // covered by the vertical mode.
      const bool horizontal_mode_allowed =
          (a0b1_dist > kGroup4Window + 1) ||
          (a0b1_dist + kGroup4Window < a0a1_max);
      // Pass mode allowed if it's not at the right edge of the tile.
      const bool pass_mode_allowed = (b2 + 1) <= width;
      const uint32_t max_mode = horizontal_mode_allowed + pass_mode_allowed;

      Group4Mode mode = Group4Mode::kVertical;
      if (max_mode == 2) {
        mode = (Group4Mode)sr->Read(kSymbolG4Type, "mode", cost_ptr);
      } else if (max_mode == 1) {
        int32_t mode_int;
        WP2_CHECK_STATUS(sr->TryRead(/*cluster=*/0, kSymbolG4Type, max_mode,
                                     "mode", &mode_int, cost_ptr));
        mode = (mode_int == 0)           ? Group4Mode::kVertical
               : horizontal_mode_allowed ? Group4Mode::kHorizontal
                                         : Group4Mode::kPass;
      }

      switch (mode) {
        case Group4Mode::kPass:
          // b2 can be equal to width so make sure we also do not get past the
          // end of the row.
          while (x < width && x <= b2) SetColor(a0_color, x++, src);
          a0 = b2;
          break;
        case Group4Mode::kHorizontal: {
          const char* const str[2] = {"hdist_bg_minus1", "hdist_fg_minus1"};
          int32_t a0a1_minus1;
          WP2_CHECK_STATUS(sr->TryRead(
              a0_color != first_color, kSymbolG4HorizontalDist, a0a1_max - 1,
              str[a0_color != first_color], &a0a1_minus1, cost_ptr));
          const uint32_t a1 = a0 + a0a1_minus1 + 1;
          while (x < a1) SetColor(a0_color, x++, src);
          if (x < width) {
            const int a1a2_max = (int)width - a1;
            int32_t a1a2_minus1;
            WP2_CHECK_STATUS(sr->TryRead(
                a0_color == first_color, kSymbolG4HorizontalDist, a1a2_max - 1,
                str[a0_color == first_color], &a1a2_minus1, cost_ptr));
            const uint32_t a2 = a1 + a1a2_minus1 + 1;
            assert(a2 <= width);
            const int16_t a1_color =
                ChangeColor(num_colors_, a0_color, b1_color, &mtf_, sr);
            current_color = a1_color;
            while (x < a2) SetColor(a1_color, x++, src);
            if (x < width) {
              const uint16_t expected_a2_color =
                  (b2_color != a1_color) ? b2_color : a0_color;
              const int16_t a2_color = ChangeColor(
                  num_colors_, a1_color, expected_a2_color, &mtf_, sr);
              current_color = a2_color;
              SetColor(a2_color, x++, src);
            }
            a0 = a2;
          }
          break;
        }
        case Group4Mode::kVertical: {
          // min_dist could also be optimized upon (e.g. when at the beginning
          // of the line) but that would change the statistics as some values
          // would be offset. E.g., if all values are the same, they would be
          // shifted by a non-constant min_dist which would result in different
          // values.
          const int min_dist = -(int)kGroup4Window;
          const int max_dist = (int)std::min(b1 + kGroup4Window, width) - b1;
          int32_t a1b1;
          WP2_CHECK_STATUS(sr->TryRead(/*cluster=*/0, kSymbolG4VerticalDist,
                                       max_dist - min_dist, "vdist", &a1b1,
                                       cost_ptr));

          const uint32_t a1 = b1 + a1b1 + min_dist;
          // We need to keep x < width as we do not optimize upon min_dist.
          while (x < std::min(width, a1)) SetColor(a0_color, x++, src);
          if (x < width) {
            const int16_t a1_color =
                ChangeColor(num_colors_, a0_color, b1_color, &mtf_, sr);
            current_color = a1_color;
            SetColor(a1_color, x++, src);
          }
          a0 = a1;
          break;
        }
      }

      if (get_cost) {
        RegisterSymbolForVDebug((int)mode, /*is_group4=*/true,
                                y * width + start_x, x - start_x, cost,
                                config_.info);
        if (bits_per_pixel != nullptr) {
          const float bits = cost / (x - start_x);
          bits_per_pixel->Fill({start_x, y, x - start_x, 1}, bits);
        }
      }
    }
    src += 4 * width;
    WP2_CHECK_STATUS(ProcessRows(y, config_.info));
  }

  if (src_out != nullptr) *src_out = src;

  return WP2_STATUS_OK;
}

}  // namespace WP2L
