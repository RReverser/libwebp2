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
//  Screen content encoding (e.g. screenshots)
//
// Author: Maryla (maryla@google.com)

#include "src/enc/screen_content/screen_enc.h"

#include "src/common/constants.h"
#include "src/common/lossless/color_cache.h"
#include "src/enc/anim/anim_enc.h"
#include "src/enc/preview/preview_enc.h"
#include "src/enc/wp2_enc_i.h"
#include "src/utils/utils.h"
#include "src/wp2/encode.h"

namespace WP2 {
namespace {

// Split the image into a grid of kScreenBlockSize x kScreenBlockSize blocks,
// and decide between lossy and lossless for each block.
// TODO(maryla): we could try using the real paritioning so the blocks follow
// image edges.
static constexpr uint32_t kScreenBlockSize = 32;

uint32_t GetDisto(const ArgbBuffer& left, const ArgbBuffer& right,
                  const Rectangle& window) {
  uint64_t disto_per_channel[4] = {0, 0, 0, 0};
  for (uint32_t y = 0; y < window.height; ++y) {
    const uint8_t* const ptr1 = (const uint8_t*)left.GetRow(window.y + y);
    const uint8_t* const ptr2 = (const uint8_t*)right.GetRow(window.y + y);
    WP2SumSquaredError4x8u(&ptr1[window.x * 4], &ptr2[window.x * 4],
                           window.width, disto_per_channel);
  }
  return disto_per_channel[0] + disto_per_channel[1] + disto_per_channel[2] +
         disto_per_channel[3];
}

bool HasOnlyOneColor(const ArgbBuffer& input, uint32_t i, uint32_t j) {
  const uint32_t first_pix =
      ((uint32_t*)input.GetRow(j * kScreenBlockSize))[i * kScreenBlockSize];
  for (uint32_t y = 0; y < kScreenBlockSize; ++y) {
    if (j * kScreenBlockSize + y >= input.height) break;
    const uint32_t* const row =
        (uint32_t*)input.GetRow(j * kScreenBlockSize + y);
    for (uint32_t x = 0; x < kScreenBlockSize; ++x) {
      if (i * kScreenBlockSize + x >= input.width) break;
      const uint32_t pix = row[i * kScreenBlockSize + x];
      if (pix != first_pix) return false;
    }
  }
  return true;
}

// Computes RD score per pixel.
WP2Status GetScores(const EncoderConfig& config, const ArgbBuffer& input,
                    float lambda, Planef* const scores) {
  MemoryWriter memory_writer;
  WP2_CHECK_STATUS(Encode(input, &memory_writer, config));

  ArgbBuffer output;
  DecoderConfig dec_config;
  DecoderInfo info;
  dec_config.info = &info;
  Planef bits_per_pixel;
  info.bits_per_pixel = &bits_per_pixel;
  WP2_CHECK_STATUS(
      Decode(memory_writer.mem_, memory_writer.size_, &output, dec_config));

  const uint32_t block_height = DivCeil(input.height, kScreenBlockSize);
  const uint32_t block_width = DivCeil(input.width, kScreenBlockSize);
  WP2_CHECK_STATUS(scores->Resize(block_width, block_height));
  scores->Fill(0.f);
  for (uint32_t y = 0; y < block_height; ++y) {
    for (uint32_t x = 0; x < block_width; ++x) {
      const Rectangle window{
          x * kScreenBlockSize, y * kScreenBlockSize,
          std::min(kScreenBlockSize, input.width - x * kScreenBlockSize),
          std::min(kScreenBlockSize, input.height - y * kScreenBlockSize)};

      for (uint32_t j = 0; j < window.height; ++j) {
        for (uint32_t i = 0; i < window.width; i++) {
          scores->At(x, y) += bits_per_pixel.At(window.x + i, window.y + j);
        }
      }

      const uint32_t disto = GetDisto(output, input, window);
      const float score = scores->At(x, y) * lambda + disto;
      scores->At(x, y) = score;
    }
  }
  return WP2_STATUS_OK;
}

WP2Status DecideLossy(const ArgbBuffer& input, const Planef& lossy_scores,
                      const Planef& lossless_scores,
                      Plane<bool>* const is_lossy) {
  const uint32_t block_width = lossy_scores.w_;
  const uint32_t block_height = lossy_scores.h_;
  WP2_CHECK_STATUS(is_lossy->Resize(block_width, block_height));
  // Choose between lossy and lossless based on score.
  for (uint32_t y = 0; y < block_height; ++y) {
    for (uint32_t x = 0; x < block_width; ++x) {
      if (HasOnlyOneColor(input, x, y)) {
        // For monochrome areas, use the same mode as the neighboring blocks.
        is_lossy->At(x, y) = (x > 0) ? is_lossy->At(x - 1, y)
                                     : (y > 0) ? is_lossy->At(x, y - 1) : false;
      } else {
        is_lossy->At(x, y) = (lossy_scores.At(x, y) < lossless_scores.At(x, y));
      }
    }
  }

  // Second step: remove small holes. If at least 3 direct neighbours share the
  // same mode (lossy or lossless), the current block gets set to that mode.
  bool something_changed;
  uint16_t iter = 0;
  const uint16_t kMaxIter = 5;
  do {
    something_changed = false;
    for (uint32_t y = 0; y < block_height; ++y) {
      for (uint32_t x = 0; x < block_width; ++x) {
        uint16_t num_total = 0;
        uint16_t num_lossy = 0;
        if (x > 0) {
          ++num_total;
          num_lossy += is_lossy->At(x - 1, y);
        }
        if (x + 1 < block_width) {
          ++num_total;
          num_lossy += is_lossy->At(x + 1, y);
        }
        if (y > 0) {
          ++num_total;
          num_lossy += is_lossy->At(x, y - 1);
        }
        if (y + 1 < block_height) {
          ++num_total;
          num_lossy += is_lossy->At(x, y + 1);
        }
        const uint16_t num_lossless = num_total - num_lossy;
        if (!is_lossy->At(x, y) &&
            ((num_total >= 3 && num_lossless <= 1) || (num_lossless == 0))) {
          is_lossy->At(x, y) = true;
          something_changed = true;
        } else if (is_lossy->At(x, y) &&
                   ((num_total >= 3 && num_lossy <= 1) || (num_lossy == 0))) {
          is_lossy->At(x, y) = false;
          something_changed = true;
        }
      }
    }
  } while (something_changed && ++iter < kMaxIter);
  return WP2_STATUS_OK;
}

Rectangle ScaleRect(const Rectangle& rect, uint32_t multiplier) {
  return {rect.x * multiplier, rect.y * multiplier, rect.width * multiplier,
          rect.height * multiplier};
}

}  // namespace

WP2Status EncodeScreenContent(const ArgbBuffer& input, Writer* const output,
                              const EncoderConfig& config) {
#if !defined(WP2_BITTRACE)
  fprintf(stderr, "WARNING! WP2_BITTRACE required for screen content\n");
  return Encode(input, output, config);
#endif

  if (config.quality > kMaxLossyQuality || config.quality <= 50) {
    return Encode(input, output, config);
  }
  const uint32_t block_width = DivCeil(input.width, kScreenBlockSize);
  const uint32_t block_height = DivCeil(input.height, kScreenBlockSize);
  // Empirical value. TODO(maryla): tune.
  const float lambda = std::max(85.f, -42.f * config.quality + 3400.f);

  const uint32_t kNearLosslessQuality = 98;

  EncoderConfig lossless_config = config;
  lossless_config.quality = kNearLosslessQuality;
  Planef lossless_scores;
  WP2_CHECK_STATUS(GetScores(lossless_config, input, lambda, &lossless_scores));
  EncoderConfig lossy_config = config;
  // We don't actually need to keep lossy and lossless_score, we could subtract
  // one from the other directly. But this is simpler/clearer.
  Planef lossy_scores;
  WP2_CHECK_STATUS(GetScores(lossy_config, input, lambda, &lossy_scores));

  Plane<bool> is_lossy;
  WP2_CHECK_STATUS(
      DecideLossy(input, lossy_scores, lossless_scores, &is_lossy));

  RectangleGroup lossy_rect_group(1);
  RectangleGroup lossless_rect_group(1);
  for (uint32_t y = 0; y < block_height; ++y) {
    for (uint32_t x = 0; x < block_width; ++x) {
      if (is_lossy.At(x, y)) {
        lossy_rect_group.AddPoint(x, y);
      } else {
        lossless_rect_group.AddPoint(x, y);
      }
    }
  }

  if (lossy_rect_group.GetNumRectangles() == 0 ||
      lossless_rect_group.GetNumRectangles() == 0) {
    if (lossless_rect_group.GetNumRectangles() > 0) {
      return Encode(input, output, lossless_config);
    }
    return Encode(input, output, lossy_config);
  }

  Subframe lossless_subframe;
  Subframe lossy_subframe;
  WP2_CHECK_STATUS(lossy_subframe.rgb_pixels.CopyFrom(input));
  WP2_CHECK_STATUS(lossless_subframe.rgb_pixels.CopyFrom(input));
  lossy_subframe.window =
      ScaleRect(lossy_rect_group.GetRectangle(0), kScreenBlockSize)
          .ClipWith(lossy_subframe.GetFrameRect());
  lossless_subframe.window =
      ScaleRect(lossless_rect_group.GetRectangle(0), kScreenBlockSize)
          .ClipWith(lossless_subframe.GetFrameRect());

  lossless_subframe.duration_ms = 0;
  lossless_subframe.dispose = true;
  lossy_subframe.duration_ms = 1;
  lossy_subframe.blend = true;

  // Remove lossy areas from the lossless frame and vice versa.
  for (uint32_t y = 0; y < block_height; ++y) {
    for (uint32_t x = 0; x < block_width; ++x) {
      const bool lossy = is_lossy.At(x, y);
      const uint32_t w =
          std::min(input.width - x * kScreenBlockSize, kScreenBlockSize);
      for (uint32_t j = 0; j < kScreenBlockSize; ++j) {
        if (y * kScreenBlockSize + j >= input.height) break;
        uint8_t* lossless_row = (uint8_t*)(lossless_subframe.rgb_pixels.GetRow(
            y * kScreenBlockSize + j));
        uint8_t* lossy_row = (uint8_t*)(lossy_subframe.rgb_pixels.GetRow(
            y * kScreenBlockSize + j));
        if (lossy) {
          std::fill(&lossless_row[(x * kScreenBlockSize) * 4],
                    &lossless_row[(x * kScreenBlockSize + w) * 4], 0);
        } else {
          std::fill(&lossy_row[(x * kScreenBlockSize) * 4],
                    &lossy_row[(x * kScreenBlockSize + w) * 4], 0);
        }
      }
    }
  }

  const bool has_alpha = input.HasTransparency();
  const bool has_icc = (input.metadata.iccp.size > 0);
  const bool has_trailing_data =
      (input.metadata.xmp.size > 0) || (input.metadata.exif.size > 0);
  WP2_CHECK_STATUS(EncodeHeader(config, input.width, input.height, has_alpha,
                                /*is_anim=*/true, /*loop_count=*/1,
                                kTransparentArgb38b, GetPreviewColor(input),
                                has_icc, has_trailing_data, output));

  SubframeEncoder frame_encode(GetQualityHint(config.quality), has_alpha);
  WP2_CHECK_STATUS(frame_encode.EncodeSubframe(lossless_config,
                                               lossless_subframe,
                                               /*is_last=*/false, output));

  WP2_CHECK_STATUS(frame_encode.EncodeSubframe(lossy_config, lossy_subframe,
                                               /*is_last=*/true, output));

  return WP2_STATUS_OK;
}

}  // namespace WP2
