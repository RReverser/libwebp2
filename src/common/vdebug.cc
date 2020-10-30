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
//   Everything related to visual debug
//
// Author: Skal (pascal.massimino@gmail.com)

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <numeric>
#include <string>
#include <unordered_set>

#include "src/common/color_precision.h"
#include "src/common/lossy/block.h"
#include "src/common/lossy/block_size.h"
#include "src/common/lossy/predictor.h"
#include "src/dec/filters/alpha_filter.h"
#include "src/dec/filters/deblocking_filter.h"
#include "src/dec/filters/directional_filter.h"
#include "src/dec/filters/grain_filter.h"
#include "src/dec/filters/intertile_filter.h"
#include "src/dec/filters/restoration_filter.h"
#include "src/dec/tile_dec.h"
#include "src/dec/wp2_dec_i.h"
#include "src/dsp/dsp.h"
#include "src/dsp/lossless/lossless_common.h"
#include "src/dsp/math.h"
#include "src/enc/partitioner.h"
#include "src/enc/wp2_enc_i.h"
#include "src/utils/utils.h"
#include "src/wp2/base.h"
#include "src/wp2/decode.h"
#include "src/wp2/format_constants.h"

namespace WP2 {

//------------------------------------------------------------------------------

// Declared in segment.h.
const Argb32b kSegmentColors[kMaxNumSegments]{
    {0xff, 0x30, 0xf0, 0xd0},  // Teal
    {0xff, 0x3d, 0x70, 0xf0},  // Blue
    {0xff, 0xb3, 0x30, 0xf0},  // Purple
    {0xff, 0xe4, 0x30, 0xf0},  // Pink
    {0xff, 0xf0, 0x30, 0x50},  // Red
    {0xff, 0xf0, 0x8c, 0x4d},  // Orange
    {0xff, 0xf0, 0xd0, 0x30},  // Yellow
    {0xff, 0x20, 0xf0, 0x60},  // Green
};

//------------------------------------------------------------------------------

// Sets pixel at (x, y) with 'color' if it fits in 'debug_output' boundaries.
static void MarkPixelForVDebug(int32_t x, int32_t y, Argb32b color,
                               ArgbBuffer* const debug_output) {
  assert(debug_output->format == WP2_Argb_32);
  if (x >= 0 && x < (int32_t)debug_output->width &&  // Within bounds.
      y >= 0 && y < (int32_t)debug_output->height) {
    ToUInt8(color, (uint8_t*)debug_output->GetRow(y) + x * 4);
  }
}

// Sets pixels with 'color' (or respectively 'bg_color' if not fully
// transparent), in the rectangle (w, h) starting at (x, y), for non-space
// characters in 'mask' (or respectively spaces).
static void DrawPixelArray(uint32_t x, uint32_t y, uint32_t w, uint32_t h,
                           Argb32b color, Argb32b bg_color,
                           const std::vector<std::string>& mask, uint32_t scale,
                           ArgbBuffer* const debug_output) {
  w += x;
  h += y;
  for (uint32_t l = 0; l < mask.size(); ++l) {
    for (uint32_t ls = 0; ls < scale; ++ls) {
      const uint32_t final_y = y + l * scale + ls;
      if (final_y >= h) break;

      for (uint32_t c = 0; mask[l][c] != '\0'; ++c) {
        for (uint32_t cs = 0; cs < scale; ++cs) {
          const uint32_t final_x = x + c * scale + cs;
          if (final_x >= w) break;

          uint8_t* const pixel =
              (uint8_t*)debug_output->GetRow(final_y) + final_x * 4;
          if (mask[l][c] != ' ') {
            ToUInt8(color, pixel);
          } else if (bg_color.a > 0) {
            ToUInt8(bg_color, pixel);
          }
        }
      }
    }
  }
}

//------------------------------------------------------------------------------

static void DrawUndefinedPattern(ArgbBuffer* const debug_output) {
  // Easily recognizable striped background.
  for (uint32_t y = 0; y < debug_output->height; ++y) {
    for (uint32_t x = 0; x < debug_output->width; ++x) {
      ToUInt8(((x % 10) == (y % 10)) ? Argb32b{0xff, 0x22, 0x22, 0x22}
                                     : Argb32b{0xff, 0x00, 0x00, 0x00},
              (uint8_t*)debug_output->GetRow(y) + x * 4);
    }
  }

  // Display "undefined" if the 'debug_output' content is not overwritten.
  DrawPixelArray(debug_output->width / 3, debug_output->height / 3,
                 debug_output->width - debug_output->width / 3,
                 debug_output->height - debug_output->height / 3,
                 Argb32b{0xff, 0x55, 0x55, 0x55}, Argb32b{0, 0, 0, 0},
                 {
                     "            #        # #             #",
                     "#  # # #   ##  ##   #    # #   ##   ##",
                     "#  # ## # # # #### ### # ## # #### # #",
                     "# ## #  # # # #     #  # #  # #    # #",
                     " # # #  #  ##  ##   #  # #  #  ##   ##",
                 },
                 /*scale=*/Clamp(debug_output->width / 128u, 1u, 4u),
                 debug_output);
}

// Called once 'features' are known to initialize 'DecoderConfig::info'.
WP2Status SetupDecoderInfo(const BitstreamFeatures& features,
                           const DecoderConfig& config) {
  if (VDMatch(config, "")) {  // Matches any visual debug.
    ArgbBuffer* const debug_output = &config.info->debug_output;
    if (VDMatch(config, "original")) {
      // Keep the content of the buffer but check the size.
      WP2_CHECK_OK(!debug_output->IsEmpty() &&
                       debug_output->width == features.raw_width &&
                       debug_output->height == features.raw_height,
                   WP2_STATUS_BAD_DIMENSION);
      // TODO(yguyon): Rotate 'debug_output' if 'orientation' is not kOriginal
    } else {
      WP2_CHECK_STATUS(
          debug_output->Resize(features.raw_width, features.raw_height));
      DrawUndefinedPattern(debug_output);
    }
  }
  if (config.info != nullptr) {  // Clear output values.
#if defined(WP2_BITTRACE)
    config.info->bit_traces.clear();
    config.info->blocks.clear();
    // Register the 'header_size' into WP2_BITTRACE.
    LabelStats* const stats =
        &config.info->bit_traces["BitstreamFeatures (not ANS)"];
    stats->bits += features.header_size * 8.;
    stats->num_occurrences += 1;
    // TODO(yguyon): Add WP2BitCounts to HeaderDec to have finer labels rather
    //               than storing the whole BitstreamFeatures as a Symbol.
    stats->type = LabelStats::Type::Symbol;
    stats->histo[0] += 1;
#endif  // WP2_BITTRACE
    config.info->header_size = features.header_size;
    config.info->selection_info.clear();
  }
  return WP2_STATUS_OK;
}

void RegisterChunkSize(const DecoderConfig& config, uint32_t chunk_size,
                       uint32_t chunk_size_size) {
  (void)config, (void)chunk_size, (void)chunk_size_size;
#if defined(WP2_BITTRACE)
  if (config.info != nullptr) {
    LabelStats* const stats = &config.info->bit_traces["ChunkSize (not ANS)"];
    stats->bits += chunk_size_size * 8.;
    stats->num_occurrences += 1;
    stats->type = LabelStats::Type::RValue;
    stats->histo[chunk_size] += 1;
  }
#endif  // WP2_BITTRACE
}

template<class V>
static void CopyPredictorNames(const V& predictors,
                               std::vector<std::string>* const names) {
  names->clear();
  names->reserve(predictors.size());
  for (Predictor* predictor : predictors) {
    names->push_back(predictor->GetName());
  }
}

void FillDecoderInfo(const GlobalParams& gparams, const DecoderConfig& config) {
  if (config.info != nullptr) {  // Clear output values.
    config.info->num_segments = gparams.segments_.size();
    CopyPredictorNames(gparams.y_preds_, &config.info->y_predictors);
    CopyPredictorNames(gparams.uv_preds_, &config.info->uv_predictors);
    CopyPredictorNames(gparams.a_preds_, &config.info->a_predictors);
    config.info->explicit_segment_ids = gparams.explicit_segment_ids_;
  }
}

WP2Status SetupEncoderInfo(uint32_t width, uint32_t height,
                           const EncoderConfig& config) {
  if (VDMatch(config, "")) {  // Matches any visual debug.
    ArgbBuffer* const debug_output = &config.info->debug_output;
    WP2_CHECK_STATUS(debug_output->Resize(width, height));
    DrawUndefinedPattern(debug_output);
  }
  if (config.info != nullptr) {
#if defined(WP2_BITTRACE)
    config.info->blocks.clear();
#endif
    config.info->selection_info.clear();
  }
  return WP2_STATUS_OK;
}

//------------------------------------------------------------------------------

// Returns true if the 'str' contains the 'token'. Tokens are separated by '/'.
static bool VDMatch(const char str[], const char token[]) {
  if (str != nullptr && token != nullptr) {
    const size_t token_length = std::strlen(token);
    if (token_length == 0) return true;

    for (const char* pos = str; (pos = std::strstr(pos, token)) != nullptr;
         pos += token_length) {
      if ((pos == str || pos[-1] == '/') &&  // Beginning of a token.
          (pos[token_length] == '\0' || pos[token_length] == '/')) {  // End.
        return true;  // A complete token is matched.
      }
    }
  }
  return false;
}

static Channel VDChannel(const char visual_debug[]) {
  if (VDMatch(visual_debug, "y")) return kYChannel;
  if (VDMatch(visual_debug, "u")) return kUChannel;
  if (VDMatch(visual_debug, "v")) return kVChannel;
  assert(VDMatch(visual_debug, "a"));
  return kAChannel;
}

// Returns true if an 'index' exists as a token in 'str' and outputs it.
static bool VDIndex(const char str[], uint32_t* const index) {
  if (str != nullptr) {
    // For each 'token' in 'str'.
    for (const char* token = str; *token != '\0'; ++token) {
      if (*token == '/') continue;  // Empty token.
      bool keep_token = true;
      uint32_t token_value = 0;
      static constexpr uint32_t max_value =
          std::numeric_limits<decltype(token_value)>::max() / 10;

      // For each character in 'token'.
      for (const char* pos = token; *pos != '\0'; ++pos) {
        if (*pos == '/') break;  // Token ended.
        if (*pos < '0' || *pos > '9') keep_token = false;
        if (keep_token && token_value > max_value) keep_token = false;
        if (keep_token) token_value = token_value * 10 + (*pos - '0');
      }
      if (keep_token) {
        *index = token_value;
        return true;
      }
    }
  }
  return false;
}

// Returns true if the point 'config.info->selection' is in 'rect'.
template <typename TConfig>
static bool VDSelected(const Rectangle& rect, const TConfig& config) {
  return (config.info != nullptr && config.info->selection.width > 0 &&
          config.info->selection.height > 0 &&
          config.info->selection.x >= rect.x &&
          config.info->selection.x < rect.x + rect.width &&
          config.info->selection.y >= rect.y &&
          config.info->selection.y < rect.y + rect.height);
}

// Declared in wp2_dec_i.h.
bool VDMatch(const DecoderConfig& config, const char token[]) {
  if (config.info != nullptr) return VDMatch(config.info->visual_debug, token);
  return false;
}
Channel VDChannel(const DecoderConfig& config) {
  return VDChannel(config.info->visual_debug);
}
bool VDSelected(uint32_t tile_x, uint32_t tile_y, const Rectangle& rect,
                const DecoderConfig& config) {
  return VDSelected({tile_x + rect.x, tile_y + rect.y, rect.width, rect.height},
                    config);
}

// Declared in wp2_enc_i.h.
bool VDMatch(const EncoderConfig& config, const char token[]) {
  if (config.info != nullptr) return VDMatch(config.info->visual_debug, token);
  return false;
}
Channel VDChannel(const EncoderConfig& config) {
  return VDChannel(config.info->visual_debug);
}
bool VDSelected(uint32_t tile_x, uint32_t tile_y, const Rectangle& rect,
                const EncoderConfig& config) {
  return VDSelected({tile_x + rect.x, tile_y + rect.y, rect.width, rect.height},
                    config);
}
float VDGetParam(const EncoderConfig& config, float default_value) {
  float value;
  return (config.info != nullptr && config.info->visual_debug != nullptr &&
          std::sscanf(config.info->visual_debug, "p=%f", &value) == 1)
             ? value
             : default_value;
}

static bool IsHorizontal(const DecoderConfig& config) {
  if (VDMatch(config, "horizontal")) return true;
  assert(VDMatch(config, "vertical"));
  return false;
}

//------------------------------------------------------------------------------

// Appends a formatted literal to a std::string. Avoids <sstream>.
#define WP2SAppend(str_ptr, ...)                                         \
  do {                                                                   \
    const size_t size = std::snprintf(nullptr, 0, __VA_ARGS__) + 1;      \
    (str_ptr)->resize((str_ptr)->size() + size);                         \
    std::snprintf((char*)&(str_ptr)->at((str_ptr)->size() - size), size, \
                  __VA_ARGS__);                                          \
    (str_ptr)->pop_back(); /* Remove ending '\0' */                      \
  } while (false)

// Simple tool for gathering and displaying min, max and average values.
class MinMaxAvg {
 public:
  void Add(double value) {
    if (num_values_ == 0) min_ = max_ = sum_ = value;
    min_ = std::min(min_, value), max_ = std::max(max_, value), sum_ += value;
    ++num_values_;
  }

  uint32_t GetNumValues() const { return num_values_; }
  double GetMin() const { return min_; }
  double GetMax() const { return max_; }
  double GetAvg() const { return (num_values_ > 0) ? sum_ / num_values_ : 0; }

  std::string ToString() const {
    std::string str;
    WP2SAppend(&str, "min: %.3f, avg: %.3f, avg: %.3f (%u)", min_, GetAvg(),
               max_, num_values_);
    return str;
  }

 private:
  uint32_t num_values_ = 0;
  double sum_ = 0, min_ = 0, max_ = 0;
};

//------------------------------------------------------------------------------

// Maps a value in [0; num_values) to the [-512;512] range.
static int16_t ToPixel10(uint32_t value, uint32_t num_values) {
  if (num_values <= 1) return 0;
  return (int16_t)(DivRound(value * 1024u, num_values - 1u) - 512);
}

// Maps a value in [0; num_values) to the [0;255] range.
static int16_t ToPixel8(uint32_t value, uint32_t num_values) {
  if (num_values <= 1) return 0;
  return (uint8_t)DivRound(value * 255u, num_values - 1u);
}

void CodedBlock::StoreTransform(const DecoderConfig& config,
                                uint32_t tile_x, uint32_t tile_y,
                                ArgbBuffer* const debug_output) const {
  const CodingParams& params = GetCodingParams(kYChannel);
  uint32_t value, max_value;
  if (VDMatch(config, "xy")) {
    value = params.tf;
    max_value = kNumTransformPairs - 1;
  } else {
    value = (VDMatch(config, "x") ? params.tf_x() : params.tf_y());
    max_value = kNumTransforms - 1;
  }
  value = max_value - value;  // It makes more sense to flip values on screen.

  const bool all_zero = !HasCoeffs(kYChannel);
  const Rectangle r_px{x_pix(), y_pix(),
                       std::min(w_pix(), debug_output->width - x_pix()),
                       std::min(h_pix(), debug_output->height - y_pix())};
  const uint8_t r = (uint8_t)(value * 255u / max_value);
  const uint8_t g = all_zero ? 0x00 : 0xCC;
  const uint8_t b = r;
  debug_output->Fill(r_px, Argb32b{0xff, r, g, b});

  if (VDSelected(tile_x, tile_y, r_px, config)) {
    debug_output->DrawRect(r_px, Argb32b{0xff, 0xff, 0x00, 0x00});
    std::string* const str_ptr = &config.info->selection_info;
    WP2SAppend(str_ptr, "Block at %u, %u px (w %u, h %u)\n", x_pix() + tile_x,
               y_pix() + tile_y, w_pix(), h_pix());
    if (all_zero) WP2SAppend(str_ptr, "All zero\n");
    WP2SAppend(str_ptr, "Transform x: %s\n", WP2TransformNames[params.tf_x()]);
    WP2SAppend(str_ptr, "Transform y: %s\n", WP2TransformNames[params.tf_y()]);
  }
}

// If the block is currently selected, outputs the pixel values from 'pixels'
// (in local block coordinates) to the debug selection_info.
static void OutputPixelValuesWithContext(const DecoderConfig& config,
                                         uint32_t tile_x, uint32_t tile_y,
                                         const CodedBlock& cb,
                                         const Plane16& pixels) {
  if (VDSelected(tile_x, tile_y, cb.blk().rect_pix(), config)) {
    // Print pixels with the same formatting as "prediction/modes".
    std::string* const str_ptr = &config.info->selection_info;
    WP2SAppend(str_ptr, "Block at %u, %u px (w %u, h %u)\n",
               tile_x + cb.x_pix(), tile_y + cb.y_pix(), cb.w_pix(),
               cb.h_pix());
    WP2SAppend(str_ptr, "\n\n\n\n");  // See DisplayPredictionMode().
    const int16_t* const context = cb.GetContext(
        VDChannel(config), /*fill_in=*/true, /*extend_right=*/false);
    const int16_t* const context_str = cb.GetContext(
        VDChannel(config), /*fill_in=*/true, /*extend_right=*/true);
    *str_ptr +=
        GetContextAndBlockPixelsStr(context, context_str, cb.w_pix(),
                                    cb.h_pix(), pixels.Row(0), pixels.Step());
  }
}

void CodedBlock::StoreResiduals(const DecoderConfig& config, uint32_t tile_x,
                                uint32_t tile_y, const QuantMtx& quant,
                                Channel channel,
                                Plane16* const dst_plane) const {
  const CodedBlock::CodingParams& params = GetCodingParams(channel);
  const BlockSize split_size = GetSplitSize(dim(), params.split_tf);
  const uint32_t split_w = BlockWidthPix(split_size);
  const uint32_t split_h = BlockHeightPix(split_size);
  const TrfSize tf_size = tdim(channel);
  const bool reduced_transform =
      (channel == kUChannel || channel == kVChannel) && is420_;
  uint32_t tf_i = 0;
  std::string* const str_ptr = &config.info->selection_info;
  for (uint32_t split_y = 0; split_y < h_pix(); split_y += split_h) {
    for (uint32_t split_x = 0; split_x < w_pix(); split_x += split_w) {
      int32_t res[kMaxBlockSizePix2];
      // TODO(skal): rnd_mtx
      quant.Dequantize(coeffs_[channel][tf_i], num_coeffs_[channel][tf_i],
                       tf_size, res);
      WP2InvTransform2D(res, params.tf_x(), params.tf_y(), split_w, split_h,
                        res, reduced_transform);

      const bool selected = VDSelected(
          tile_x, tile_y,
          {x_pix() + split_x, y_pix() + split_y, split_w, split_h}, config);
      for (uint32_t y = 0; y < split_h; ++y) {
        for (uint32_t x = 0; x < split_w; ++x) {
          dst_plane->At(x_pix() + split_x + x, y_pix() + split_y + y) =
              res[x + split_w * y];
          if (selected) WP2SAppend(str_ptr, "%4d ", res[x + split_w * y]);
        }
        if (selected) WP2SAppend(str_ptr, "\n");
      }
      ++tf_i;
    }
  }
}

void CodedBlock::StoreOriginalResiduals(
    const EncoderConfig& config, uint32_t tile_pos_x, uint32_t tile_pos_y,
    int32_t original_res[kMaxNumTransformsPerBlock][kMaxBlockSizePix2],
    Plane16* const dst_plane) const {
  const BlockSize split_size =
      GetSplitSize(dim(), GetCodingParams(VDChannel(config)).split_tf);
  const uint32_t split_w = BlockWidthPix(split_size);
  const uint32_t split_h = BlockHeightPix(split_size);
  uint32_t tf_i = 0;
  std::string* const str_ptr = &config.info->selection_info;
  for (uint32_t split_y = 0; split_y < h_pix(); split_y += split_h) {
    for (uint32_t split_x = 0; split_x < w_pix(); split_x += split_w) {
      const bool selected = VDSelected(
          tile_pos_x, tile_pos_y,
          {x_pix() + split_x, y_pix() + split_y, split_w, split_h}, config);
      for (uint32_t k = 0, y = 0; y < split_h; ++y) {
        for (uint32_t x = 0; x < split_w; ++x, ++k) {
          const int32_t value = original_res[tf_i][k];
          dst_plane->At(x_pix() + split_x + x, y_pix() + split_y + y) = value;
          if (selected) WP2SAppend(str_ptr, "%4d ", value);
        }
        if (selected) WP2SAppend(str_ptr, "\n");
      }
      ++tf_i;
    }
  }
}

static Argb32b CtxtToColor(int16_t context) {
  const uint8_t v = Clamp(DivRound((context + 512) * 255, 1023), 0x00, 0xff);
  return Argb32b{0xff, v, v, v};
}

// Draws the context used for the prediction of 'cb'.
static void DisplayPredContext(const CodedBlock& cb, Channel channel,
                               bool draw_red_border,
                               ArgbBuffer* const debug_output) {
  const uint32_t x = cb.x_pix(), y = cb.y_pix(), w = cb.w_pix(), h = cb.h_pix();
  const Rectangle rect = {x, y, std::min(w, debug_output->width - x),
                                std::min(h, debug_output->height - y)};

  if (draw_red_border) {  // Encircle the left, top, right context lines.
    const uint32_t bx = SafeSub(x, 2u), by = SafeSub(y, 2u);
    const Argb32b red = {0xff, 0xff, 0, 0};
    debug_output->DrawRect({bx, by, x - bx + 1, y + rect.height - by}, red);
    debug_output->DrawRect({bx, by, x + rect.width - bx, y - by + 1}, red);
    debug_output->DrawRect({x + w - 1, by, 3, y + rect.height - by}, red);
    debug_output->DrawRect(
        {x + w - 1, by, 2 + ContextWithTRSize(w, h) - (h + 1 + w), 3}, red);
  }

  const int16_t* const context =
      cb.GetContext(channel, /*fill_in=*/true, /*extend_right=*/false);
  if (x >= 1) {
    for (uint32_t j = 0; j < rect.height; ++j) {  // Left
      debug_output->Fill({x - 1, y + j, 1, 1}, CtxtToColor(context[h - j - 1]));
    }
    if (y >= 1) {  // Top-left
      debug_output->Fill({x - 1, y - 1, 1, 1}, CtxtToColor(context[h]));
    }
  }
  if (y >= 1) {
    for (uint32_t i = 0; i < rect.width; ++i) {  // Top
      debug_output->Fill({x + i, y - 1, 1, 1}, CtxtToColor(context[h + 1 + i]));
    }
    if (x + w < debug_output->width) {  // Top-right
      debug_output->Fill({x + w, y - 1, 1, 1}, CtxtToColor(context[h + 1 + w]));
    }
  }
  if (x + w < debug_output->width) {
    for (uint32_t j = 0; j < rect.height; ++j) {  // Right
      debug_output->Fill({x + w, y + j, 1, 1},
                         CtxtToColor(context[h + 1 + w + 1 + j]));
    }
  }

  // Extended right
  const int16_t* const tr_context =
      cb.GetContext(channel, /*fill_in=*/true, /*extend_right=*/true);
  const uint32_t tr_size = ContextWithTRSize(w, h);
  for (uint32_t i = h + 1 + w + 1, tr_x = x + w + 1; i < tr_size; ++i, ++tr_x) {
    debug_output->Fill({tr_x, y - 1, 1, 1}, CtxtToColor(tr_context[i]));
  }
}

void CodedBlock::StorePredictionScore(const EncoderConfig& config,
                                      const Rectangle& tile_rect,
                                      Channel channel, const Predictor& pred,
                                      TransformPair tf, float dist,
                                      float lambda, float res_rate,
                                      float pred_rate, float tf_rate,
                                      float score, bool is_best) const {
  if (!VDMatch(config, "prediction-scores")) return;
  if (VDChannel(config) != channel) return;

  ArgbBuffer debug_output;  // Tile view.
  WP2_ASSERT_STATUS(debug_output.SetView(config.info->debug_output, tile_rect));

  const Rectangle rect = {  // Block rect.
      x_pix(), y_pix(), std::min(w_pix(), debug_output.width  - x_pix()),
                        std::min(h_pix(), debug_output.height - y_pix())};
  const bool selected = VDSelected(tile_rect.x, tile_rect.y, rect, config);

  debug_output.Fill(rect, Argb32b{0xff, 0x00, 0x22, 0x00});
  if (selected) {
    debug_output.DrawRect(rect, Argb32b{0xff, 0xff, 0x33, 0x33});

    // Check the formula did not change.
    assert(std::abs(dist + lambda * (res_rate + pred_rate + tf_rate) - score) <
           0.001f);

    // Print at least all better predictors. Skip worse ones if too many.
    std::string* const str_ptr = &config.info->selection_info;
    if (is_best || std::count(str_ptr->begin(), str_ptr->end(), '\n') < 30) {
      WP2SAppend(str_ptr, "Score %9.1f = dist %9.1f", score, dist);
      WP2SAppend(str_ptr,
                 " + lambda %9.1f * (res %7.4f + pred %7.4f + tf %7.4f)",
                 lambda, res_rate, pred_rate, tf_rate);
      WP2SAppend(str_ptr, " %s (%s - %s/%s x%u)\n", is_best ? "best" : "    ",
                 pred.GetName().c_str(), WP2TransformNames[kTfX[tf]],
                 WP2TransformNames[kTfY[tf]], GetNumTransforms(channel));
    }
  }
  DisplayPredContext(*this, channel, /*draw_red_border=*/selected,
                     &debug_output);
}

// Prints the selected block prediction to 'config.info->selection_info' and
// fills either the 'raw_prediction' (values) or the 'debug_output' (mode).
static void DisplayPredictionMode(const DecoderConfig& config,
                                  const Rectangle& tile_rect,
                                  const CodedBlock& cb, uint32_t num_preds,
                                  const Predictor& p, Channel channel,
                                  const Rectangle& block_rect, bool selected,
                                  Plane16* const raw_prediction,
                                  ArgbBuffer* const debug_output) {
  if (raw_prediction != nullptr) {
    Plane16 dst_view;
    WP2_ASSERT_STATUS(dst_view.SetView(*raw_prediction, cb.blk().rect_pix()));
    const Predictor& pred = *cb.GetCodingParams(channel).pred;
    pred.Predict(cb, channel, dst_view.Row(0), dst_view.Step());
  } else {
    assert(debug_output != nullptr);
    const uint8_t c = ToPixel8(p.mode(), num_preds);
    debug_output->Fill(block_rect, Argb32b{255, c, c, c});
    p.Draw(block_rect, debug_output);
    if (selected) {
      debug_output->DrawRect(block_rect, Argb32b{0xff, 0xff, 0x33, 0x33});
    }
  }

  if (selected) {
    std::string* const str_ptr = &config.info->selection_info;
    WP2SAppend(str_ptr, "Prediction block at %u, %u px (w %u, h %u)\n\n",
               block_rect.x + tile_rect.x, block_rect.y + tile_rect.y,
               block_rect.width, block_rect.height);

    WP2SAppend(str_ptr, " %u: %s\n", p.mode(), p.GetName().c_str());

    if (!VDMatch(config, "short")) {
      const std::string prediction = p.GetPredStr(cb, channel);
      WP2SAppend(str_ptr, "\nPrediction:\n%s", prediction.c_str());
    }
  }
}

void CodedBlock::StorePredictionModes(const DecoderConfig& config,
                                      const Rectangle& tile_rect,
                                      Channel channel, const Predictors& preds,
                                      Plane16* const raw_prediction,
                                      ArgbBuffer* const debug_output) const {
  const Rectangle r_px{x_pix(), y_pix(),
                       std::min(w_pix(), tile_rect.width - x_pix()),
                       std::min(h_pix(), tile_rect.height - y_pix())};

  // Whole block
  const bool selected = VDSelected(tile_rect.x, tile_rect.y, r_px, config);

  const CodingParams& params = GetCodingParams(channel);
  DisplayPredictionMode(config, tile_rect, *this, preds.GetMaxMode(),
                        *params.pred, channel, r_px, selected, raw_prediction,
                        debug_output);

  if (selected && debug_output != nullptr) {
    DisplayPredContext(*this, channel, /*draw_red_border=*/true, debug_output);
  }
}

void CodedBlock::AppendOriginalPixels(const DecoderConfig& config,
                                      uint32_t tile_x, uint32_t tile_y,
                                      const CSPTransform& csp_tranform,
                                      ArgbBuffer* const debug_output) const {
  if (VDSelected(tile_x, tile_y, {x_pix(), y_pix(), w_pix(), h_pix()},
                 config)) {
    const Channel channel = VDChannel(config);
    YUVPlane original_yuv;  // 'debug_output' should be the original RGB pixels.
    WP2_ASSERT_STATUS(
        original_yuv.Import(*debug_output, /*import_alpha=*/true, csp_tranform,
                            /*resize_if_needed=*/true, /*pad=*/kPredWidth));
    // Use a dummy CodedBlock to get the 'context' out of 'original_yuv'.
    CodedBlock cb;
    cb.SetDimDefault(blk(), /*full_left_ctx=*/true);
    if (cb.x_pix() + cb.w_pix() < original_yuv.Y.w_) {
      cb.right_occupancy_ = cb.y() + cb.h();  // Full right context.
    }
    ContextCache context_cache;
    cb.SetContextInput(original_yuv, &context_cache);

    YUVPlane view;
    WP2_ASSERT_STATUS(view.SetView(original_yuv, blk().rect_pix()));
    OutputPixelValuesWithContext(config, tile_x, tile_y, cb,
                                 view.GetChannel(channel));
  }
}

void CodedBlock::AppendCompressedPixels(const DecoderConfig& config,
                                        uint32_t tile_x, uint32_t tile_y,
                                        ArgbBuffer* const debug_output) const {
  if (VDMatch(config, "r") || VDMatch(config, "g") || VDMatch(config, "b")) {
    return;
  }
  if (VDSelected(tile_x, tile_y, blk().rect_pix(), config)) {
    OutputPixelValuesWithContext(config, tile_x, tile_y, *this,
                                 out_.GetChannel(VDChannel(config)));
    WP2SAppend(
        &config.info->selection_info,
        "Compressed values are before applying post processing filters\n");
  }
}

void SyntaxReader::StoreBitCost(const DecoderConfig& config,
                                uint32_t tile_x, uint32_t tile_y,
                                const Block& block,
                                Plane16* const dst_plane) const {
  (void)config, (void)tile_x, (void)tile_y, (void)block, (void)dst_plane;
#if defined(WP2_BITTRACE)
  const Rectangle& rect = block.rect_pix();
  const std::string prefix =
      VDMatch(config, "overall") ? "" : kCoeffsStr[VDChannel(config)];

  double num_bits = 0.;
  for (const auto& bt : dec_->GetBitTracesCustom()) {
    if (bt.first.size() >= prefix.size() &&
        std::strncmp(bt.first.c_str(), prefix.c_str(), prefix.size()) == 0) {
      num_bits += bt.second.bits;
    }
  }

  // Max value is 6bpp (which is a lot!). Use 8b fixed point precision.
  const int16_t value = ToPixel10(
      (uint32_t)std::lround(num_bits * 256. / rect.GetArea()), 6u << 8);
  dst_plane->Fill(rect, value);

  if (VDSelected(tile_x, tile_y, rect, config)) {
    dst_plane->DrawRect(rect, 127);
    std::string* const str_ptr = &config.info->selection_info;
    WP2SAppend(str_ptr, "\nBlock at %u, %u px\n\n", rect.x, rect.y);
    WP2SAppend(str_ptr, "%f bits over %u pixels\n", num_bits,
               (uint32_t)rect.GetArea());
    WP2SAppend(str_ptr, " = %f bpp\n", num_bits / rect.GetArea());
  }
#endif  // defined(WP2_BITTRACE)
}

void CodedBlock::StoreCoeffMethod(const DecoderConfig& config,
                                  Plane16* const dst_plane) const {
  const Channel channel = VDChannel(config);
  const BlockSize split_size =
      GetSplitSize(dim(), GetCodingParams(channel).split_tf);
  const uint32_t split_w = BlockWidthPix(split_size);
  const uint32_t split_h = BlockHeightPix(split_size);
  uint32_t tf_i = 0;
  for (uint32_t split_y = 0; split_y < h_pix(); split_y += split_h) {
    for (uint32_t split_x = 0; split_x < w_pix(); split_x += split_w) {
      dst_plane->Fill(
          {x_pix() + split_x, y_pix() + split_y, split_w - 1, split_h - 1},
          ToPixel10((uint32_t)method_[channel][tf_i],
                    (uint32_t)EncodingMethod::kNumMethod));
      ++tf_i;
    }
  }
}

// Similar to CflPredictor::LinearRegression, but uses the actual block
// instead of the block context.
void CodedBlock::CflLinearRegression(Channel channel, int16_t yuv_min,
                                     int16_t yuv_max, float* const a,
                                     float* const b,
                                     std::string* const debug_str) const {
  int16_t a_int, b_int;
  CflPredictor cfl_predictor(yuv_min, yuv_max);
  cfl_predictor.BlockLinearRegression(channel, *this, &a_int, &b_int);
  *a = (float)a_int / (1 << CflPredictor::kAPrecShift);
  *b = (float)b_int / (1 << CflPredictor::kBPrecShift);

  if (debug_str != nullptr) {
    WP2SAppend(debug_str, "Chroma = %.2f * luma + %.2f\n", *a, *b);
  }
}

void CodedBlock::StoreBestCflSlope(Channel channel, int16_t yuv_min,
                                   int16_t yuv_max, Plane16* const dst_plane,
                                   std::string* const debug_str) const {
  float a, b;
  CflLinearRegression(channel, yuv_min, yuv_max, &a, &b, debug_str);
  dst_plane->Fill(blk_.rect_pix(), (int16_t)std::lround(a * 100.f));
}

void CodedBlock::StoreCflSlope(Channel channel, int16_t yuv_min,
                               int16_t yuv_max, Plane16* const dst_plane,
                               std::string* const debug_str) const {
  int16_t a, b;
  CflPredictor cfl_predictor(yuv_min, yuv_max);
  cfl_predictor.ContextLinearRegression(channel, *this, &a, &b);
  dst_plane->Fill(
      blk_.rect_pix(),
      (int16_t)std::lround(a * 100.f / (1 << CflPredictor::kAPrecShift)));

  if (debug_str != nullptr) {
    WP2SAppend(debug_str, "Chroma = %.2f * luma + %.2f\n",
               (float)a / (1 << CflPredictor::kAPrecShift),
               (float)b / (1 << CflPredictor::kBPrecShift));
  }
}

void CodedBlock::StoreCflIntercept(Channel channel, int16_t yuv_min,
                                   int16_t yuv_max, Plane16* const dst_plane,
                                   std::string* const debug_str) const {
  int16_t a, b;
  CflPredictor cfl_predictor(yuv_min, yuv_max);
  cfl_predictor.ContextLinearRegression(channel, *this, &a, &b);
  dst_plane->Fill(
      blk_.rect_pix(),
      (int16_t)std::lround(b / (1 << (CflPredictor::kBPrecShift + 1))));

  if (debug_str != nullptr) {
    WP2SAppend(debug_str, "Chroma = %.2f * luma + %.2f\n",
               (float)a / (1 << CflPredictor::kAPrecShift),
               (float)b / (1 << CflPredictor::kBPrecShift));
  }
}

void CodedBlock::StoreBestCflIntercept(Channel channel, int16_t yuv_min,
                                       int16_t yuv_max,
                                       Plane16* const dst_plane,
                                       std::string* const debug_str) const {
  float a, b;
  CflLinearRegression(channel, yuv_min, yuv_max, &a, &b, debug_str);
  dst_plane->Fill(blk_.rect_pix(), (int16_t)std::lround(b / 2.f));
}

void CodedBlock::StoreBestCflPrediction(Channel channel, int16_t yuv_min,
                                        int16_t yuv_max,
                                        Plane16* const dst_plane,
                                        std::string* const debug_str) const {
  const Plane16& luma = out_.GetChannel(kYChannel);

  float a, b;
  CflLinearRegression(channel, yuv_min, yuv_max, &a, &b, debug_str);

  for (uint32_t y = 0; y < h_pix(); ++y) {
    for (uint32_t x = 0; x < w_pix(); ++x) {
      const int32_t v = std::lround(a * luma.At(x, y) + b);
      dst_plane->At(x_pix() + x, y_pix() + y) = (int16_t)v;
    }
  }
}

void CodedBlock::StoreBestCflResiduals(Channel channel, int16_t yuv_min,
                                       int16_t yuv_max,
                                       Plane16* const dst_plane,
                                       std::string* const debug_str) const {
  const Plane16& luma = out_.GetChannel(kYChannel);
  const Plane16& chroma = in_.GetChannel(channel);

  float a, b;
  CflLinearRegression(channel, yuv_min, yuv_max, &a, &b, debug_str);

  for (uint32_t y = 0; y < h_pix(); ++y) {
    for (uint32_t x = 0; x < w_pix(); ++x) {
      const int32_t v = std::lround(a * luma.At(x, y) + b);
      dst_plane->At(x_pix() + x, y_pix() + y) = (int16_t)(chroma.At(x, y) - v);
    }
  }
}

//------------------------------------------------------------------------------

void CodedBlock::StoreErrorDiffusion(const EncoderConfig& config,
                                     uint32_t tile_x, uint32_t tile_y,
                                     Plane16* const dst_plane) const {
  if (!VDMatch(config, "error-diffusion")) return;
  const Channel channel = VDChannel(config);

  if (VDSelected(tile_x, tile_y, blk().rect_pix(), config)) {
    WP2SAppend(&config.info->selection_info,
               "Block at %4u, %4u (%2u x %2u px)\n", tile_x + x_pix(),
               tile_y + y_pix(), w_pix(), h_pix());
    WP2SAppend(&config.info->selection_info,
               "Propagated error = %d | New error = %d\n", dc_error_[channel],
               dc_error_next_[channel]);
  }

  const int16_t error = VDMatch(config, "propagated-error")
                            ? dc_error_[channel]
                            : dc_error_next_[channel];
  dst_plane->Fill(blk_.rect_pix(), error / 2);
}

//------------------------------------------------------------------------------

void Block::Draw(YUVPlane* const yuv) const {
  const uint32_t x0 = x_pix();
  const uint32_t y0 = y_pix();
  for (uint32_t i = 0; i < w_pix(); ++i) {
    yuv->Y.At(x0 + i, y0) = 0x1ff;
    yuv->U.At(x0 + i, y0) = 0x100;
    yuv->V.At(x0 + i, y0) = 0x100;
  }
  for (uint32_t j = 0; j < h_pix(); ++j) {
    yuv->Y.At(x0, y0 + j) = 0x1ff;
    yuv->U.At(x0, y0 + j) = 0x100;
    yuv->V.At(x0, y0 + j) = 0x100;
  }
}

//------------------------------------------------------------------------------

WP2Status CodedBlock::Draw(const DecoderConfig& config, uint32_t tile_x,
                           uint32_t tile_y, const GlobalParams& gparams,
                           ArgbBuffer* const debug_output) const {
  const uint32_t x = x_pix(), bw = w_pix();
  const uint32_t y = y_pix(), bh = h_pix();
  const uint32_t w = std::min(bw, debug_output->width - x);  // Border
  const uint32_t h = std::min(bh, debug_output->height - y);
  const uint32_t mid_w = std::min(bw - 2, debug_output->width - (x + 1));
  const uint32_t mid_h = std::min(bh - 2, debug_output->height - (y + 1));

  if (VDMatch(config, "partition")) {
    if (VDMatch(config, "split-tf")) {
      const BlockSize split_size =
          GetSplitSize(dim(), GetCodingParams(kYChannel).split_tf);
      const uint32_t split_w = BlockWidthPix(split_size);
      const uint32_t split_h = BlockHeightPix(split_size);
      for (uint32_t split_y = 0; split_y < h_pix(); split_y += split_h) {
        for (uint32_t split_x = 0; split_x < w_pix(); split_x += split_w) {
          debug_output->Fill({x + split_x, y + split_y, 1, split_h},
                             Argb32b{0xff, 0x00, 0x30, 0xf0});
          debug_output->Fill({x + split_x, y + split_y, split_w, 1},
                             Argb32b{0xff, 0x00, 0x40, 0xf0});
        }
      }
    }
    debug_output->Fill({x, y, 1, h}, Argb32b{0xff, 0xff, 0xff, 0x30});
    debug_output->Fill({x, y, w, 1}, Argb32b{0xff, 0xff, 0xff, 0x40});
  } else if (VDMatch(config, "segment-ids") ||
             VDMatch(config, "quantization")) {
    const bool is_alpha =
        VDMatch(config, "quantization") && (VDChannel(config) == kAChannel);
    if (is_alpha && !HasLossyAlpha()) {
      debug_output->Fill({x, y, w, h}, Argb32b{255, 0, 0, 0});
      return WP2_STATUS_OK;
    }
    const Segment& segment =
        is_alpha ? gparams.a_segments_[0] : gparams.segments_[id_];

    Argb32b color;
    if (VDMatch(config, "segment-ids")) {
      color = kSegmentColors[id_];
      if (mid_w > 0 && mid_h > 0) {
        debug_output->Fill({x + 1, y + 1, mid_w, mid_h}, color);
      }
    } else {
      // showing quantization
      uint16_t quant[kNumQuantZones];
      segment.GetQuantSteps(quant, VDChannel(config));
      // Quant expressed in [0 .. 1]
      const float normalized_quant = (float)quant[0] / (float)kMaxQuantStep;

      // Gamma it not really necessary but just to show I studied colors.
      constexpr float kGamma = 2.2f;
      const float red = std::pow(normalized_quant, 1.f / kGamma);
      const float green = std::pow((1.f - normalized_quant), 1.f / kGamma);
      // Max out saturation.
      const float mult = std::min(1.f / red, 1.f / green) * 255.f;

      color = Argb32b{255, (uint8_t)std::min((uint32_t)(red * mult), 255u),
                      (uint8_t)std::min((uint32_t)(green * mult), 255u), 0};
      debug_output->Fill({x, y, w, h}, color);
    }

    if (VDSelected(tile_x, tile_y, {x, y, bw, bh}, config)) {
      std::string* const str_ptr = &config.info->selection_info;
      WP2SAppend(str_ptr, "\nSegment id: %u\n", id_);
      if (segment.use_quality_factor_) {
        WP2SAppend(str_ptr, " quant_factor: %d\n", segment.quality_factor_);
      }
      if (is_alpha) {
        WP2SAppend(str_ptr, " A quant_steps:   %d %d %d %d\n",
                   segment.quant_steps_[kAChannel][0],
                   segment.quant_steps_[kAChannel][1],
                   segment.quant_steps_[kAChannel][2],
                   segment.quant_steps_[kAChannel][3]);
      } else {
        WP2SAppend(str_ptr, " Y quant_steps:   %d %d %d %d\n",
                   segment.quant_steps_[kYChannel][0],
                   segment.quant_steps_[kYChannel][1],
                   segment.quant_steps_[kYChannel][2],
                   segment.quant_steps_[kYChannel][3]);
        WP2SAppend(str_ptr, " U quant_steps: %d %d %d %d\n",
                   segment.quant_steps_[kUChannel][0],
                   segment.quant_steps_[kUChannel][1],
                   segment.quant_steps_[kUChannel][2],
                   segment.quant_steps_[kUChannel][3]);
        WP2SAppend(str_ptr, " V quant_steps: %d %d %d %d\n",
                   segment.quant_steps_[kVChannel][0],
                   segment.quant_steps_[kVChannel][1],
                   segment.quant_steps_[kVChannel][2],
                   segment.quant_steps_[kVChannel][3]);
      }

      debug_output->Fill({x, y, w, h}, Argb32b{255, 0, 0, 0});
      if (mid_w > 0 && mid_h > 0) {
        debug_output->Fill({x + 1, y + 1, mid_w, mid_h}, color);
      }
    }
  } else if (VDMatch(config, "is420")) {
    debug_output->DrawRect(
        {x, y, w, h}, is420_ ? Argb32b{0xff, 0x11, 0xff, 0x11}
                             : Argb32b{0xff, 0xdd, 0xdd, 0xdd});
    if (is420_ && mid_w > 0 && mid_h > 0) {
      DrawPixelArray(x + 1, y + 1, mid_w, mid_h,
                     Argb32b{0xff, 0xaa, 0xff, 0xaa},
                     Argb32b{0xff, 0x00, 0x00, 0x00},
                     {
                         "             ",
                         "   #  #   #  ",
                         "  ## # # # # ",
                         " # #   # # # ",
                         " ###  #  # # ",
                         "   # ###  #  ",
                         "             ",
                     },
                     1, debug_output);
    }
  }
  return WP2_STATUS_OK;
}

//------------------------------------------------------------------------------

void CodedBlock::ToBlockInfo(BlockInfo* const blk) const {
  const CodingParams& y_params = GetCodingParams(kYChannel);

  blk->rect = {x_pix(), y_pix(), w_pix(), h_pix()};
  blk->segment_id = id_;
  blk->is420 = is420_;
  blk->tf_x = (uint8_t)y_params.tf_x();
  blk->tf_y = (uint8_t)y_params.tf_y();
  blk->y_context_is_constant = y_context_is_constant_;

  blk->y_pred = y_params.pred->mode();
  blk->has_lossy_alpha = HasLossyAlpha();
  if (HasLossyAlpha()) {
    const CodingParams& a_params = GetCodingParams(kAChannel);
    blk->a_pred = a_params.pred->mode();
  }
  blk->uv_pred = GetCodingParams(kUChannel).pred->mode();
  std::memcpy(blk->pred_scores, pred_scores_, sizeof(pred_scores_));

  std::memset(blk->coeffs, 0, sizeof(blk->coeffs));
  std::fill(blk->split_tf, blk->split_tf + 4, false);
  std::fill(blk->encoding_method[0], blk->encoding_method[0] + 4 * 4, -1);
  for (Channel c : {kYChannel, kUChannel, kVChannel, kAChannel}) {
    if (c == kAChannel && !HasLossyAlpha()) continue;
    blk->split_tf[c] = GetCodingParams(c).split_tf;
    for (uint32_t i = 0; i < GetNumTransforms(c); ++i) {
      blk->encoding_method[c][i] = (int8_t)method_[c][i];
      std::copy(coeffs_[c][i], coeffs_[c][i] + NumCoeffsPerTransform(c),
                blk->coeffs[c][i]);
    }
  }
  blk->bits = 0;  // unknown for now
}

//------------------------------------------------------------------------------

void CodedBlock::Store420Scores(const EncoderConfig& config, uint32_t pos_x,
                                uint32_t pos_y, float lambda_u, float lambda_v,
                                bool reduced, uint32_t disto, float rate_u,
                                float rate_v) {
  if (!VDMatch(config, "is420-scores")) return;
  if (VDSelected(pos_x, pos_y, blk().rect_pix(), config)) {
    WP2SAppend(&config.info->selection_info,
               "%s score: %f = disto %d + lambda_u %f * rate_u %f + lambda_v "
               "%f * rate_v %f\n",
               reduced ? "420" : "444",
               disto + lambda_u * rate_u + lambda_v * rate_v, disto,
               lambda_u, rate_u, lambda_v, rate_v);
  }
}

WP2Status CodedBlock::Store420Decision(const EncoderConfig& config,
                                       uint32_t pos_x, uint32_t pos_y,
                                       Debug420Decision decision) const {
  if (!VDMatch(config, "is420-scores")) return WP2_STATUS_OK;

  if (VDSelected(pos_x, pos_y, blk().rect_pix(), config)) {
    WP2SAppend(&config.info->selection_info,
               decision == Debug420Decision::k444EarlyExit
                   ? "Decision: 444 early exit (DC only)\n"
                   : decision == Debug420Decision::k444 ? "Decision: 444\n"
                                                        : "Decision: 420\n");
  }

  ArgbBuffer& debug_output = config.info->debug_output;
  Argb32b colors[3] = {Argb32b{0xff, 0xff, 0x11, 0xff},
                       Argb32b{0xff, 0xff, 0x11, 0x11},
                       Argb32b{0xff, 0x11, 0xff, 0x11}};
  const Rectangle rect = {x_pix(), y_pix(),
                          std::min(w_pix(), debug_output.width - x_pix()),
                          std::min(h_pix(), debug_output.height - y_pix())};
  debug_output.Fill(rect, colors[(int)decision]);
  debug_output.DrawRect(rect, Argb32b{0xff, 0x00, 0x00, 0x00});
  return WP2_STATUS_OK;
}

WP2Status CodedBlock::StoreLambdaMult(const EncoderConfig& config,
                                      uint32_t pos_x, uint32_t pos_y) const {
  if (!VDMatch(config, "lambda-mult")) return WP2_STATUS_OK;
  const bool selected = VDSelected(pos_x, pos_y, blk().rect_pix(), config);

  if (selected) {
    WP2SAppend(&config.info->selection_info, "lambda-mult: %.3f", lambda_mult_);
  }

  ArgbBuffer& debug_output = config.info->debug_output;
  pos_x += x_pix();
  pos_y += y_pix();
  const Rectangle rect = {pos_x, pos_y,
                          std::min(w_pix(), debug_output.width - pos_x),
                          std::min(h_pix(), debug_output.height - pos_y)};
  const uint8_t gray = (uint8_t)Clamp((lambda_mult_ - 1.0f) * 500.f + 128.f,
                                      0.f, 255.f);
  debug_output.Fill(rect, Argb32b{0xff, gray, gray, gray});
  if (selected) {
    debug_output.DrawRect(rect, Argb32b{0xff, 0xc0, 0x30, 0x00});
  }
  return WP2_STATUS_OK;
}

//------------------------------------------------------------------------------

void Tile::Draw(const DecoderConfig& config) {
  Argb32b color{0xff, 0xbb, 0xbb, 0x66};

  if (VDSelected(rect, config)) {
    color.r = color.g = color.b = 0xff;  // Highlight the tile's border.
    std::string* const str_ptr = &config.info->selection_info;
    WP2SAppend(str_ptr, "\nTile at %u, %u px (%u x %u)\n", rect.x, rect.y,
               rect.width, rect.height);
    assert(chunk_size_is_known);
    WP2SAppend(str_ptr, "Tile size: %zu bytes\n", chunk_size);
  }

  config.info->debug_output.DrawRect(rect, color);
}

//------------------------------------------------------------------------------

// Copies pixels in 'area'.
static void CopyArea(const YUVPlane& from, const Rectangle& area,
                     YUVPlane* const to) {
  for (Channel c : {kYChannel, kUChannel, kVChannel, kAChannel}) {
    if (c == kAChannel && !from.HasAlpha()) continue;

    const Plane16& from_plane = from.GetChannel(c);
    Plane16* const to_plane = &to->GetChannel(c);

    for (uint32_t y = 0; y < area.height; ++y) {
      for (uint32_t x = 0; x < area.width; ++x) {
        const uint32_t px = area.x + x, py = area.y + y;
        to_plane->At(px, py) = from_plane.At(px, py);
      }
    }
  }
}

// Saves the unfiltered 'canvas' to compare with the filtered ones.
static void SavePixelsForVDebugDiff(const YUVPlane& canvas, uint32_t num_rows,
                                    YUVPlane* const unfiltered_pixels,
                                    uint32_t* const num_debug_rows) {
  if (unfiltered_pixels->Y.w_ != canvas.Y.w_ ||
      unfiltered_pixels->Y.h_ != canvas.Y.h_) {
    WP2_ASSERT_STATUS(unfiltered_pixels->Resize(canvas.Y.w_, canvas.Y.h_,
                                                /*pad=*/1, canvas.HasAlpha()));
  }

  // Copy the unfiltered available rows.
  assert(*num_debug_rows <= num_rows && num_rows <= canvas.Y.h_);
  CopyArea(canvas,
           {0, *num_debug_rows, canvas.Y.w_, num_rows - *num_debug_rows},
           unfiltered_pixels);
  *num_debug_rows = num_rows;
}

// Set pixels of 'debug_output' in 'rect' as the difference between the pixels
// of 'unfiltered_pixels' and 'filtered_pixels', scaled by 'multiplier'.
static void ApplyVDebugDiff(const YUVPlane& unfiltered_pixels,
                            const YUVPlane& filtered_pixels,
                            int32_t yuv_max_value, int32_t multiplier,
                            const Rectangle& rect,
                            ArgbBuffer* const debug_output) {
  assert(unfiltered_pixels.GetWidth() == filtered_pixels.GetWidth());
  assert(unfiltered_pixels.GetHeight() == filtered_pixels.GetHeight());
  assert(unfiltered_pixels.GetWidth() >= rect.width);
  assert(unfiltered_pixels.GetHeight() >= rect.height);
  YUVPlane diff_pixels;
  WP2_ASSERT_STATUS(diff_pixels.Resize(rect.width, rect.height));

  for (Channel c : {kYChannel, kUChannel, kVChannel}) {
    const Plane16& unfiltered = unfiltered_pixels.GetChannel(c);
    const Plane16& filtered = filtered_pixels.GetChannel(c);
    Plane16* const diff = &diff_pixels.GetChannel(c);
    for (uint32_t y = 0; y < rect.height; ++y) {
      for (uint32_t x = 0; x < rect.width; ++x) {
        diff->At(x, y) =
            Clamp((filtered.At(x, y) - unfiltered.At(x, y)) * multiplier,
                  -yuv_max_value, yuv_max_value);
      }
    }
  }
  // TODO(yguyon): Use the actual CSPTransform
  ArgbBuffer debug_output_view(debug_output->format);
  WP2_ASSERT_STATUS(debug_output_view.SetView(*debug_output, rect));
  WP2_ASSERT_STATUS(diff_pixels.Export(
      CSPTransform(), /*resize_if_needed=*/false, &debug_output_view));
}
static void ApplyVDebugAlphaDiff(const YUVPlane& unfiltered_pixels,
                                 const YUVPlane& filtered_pixels,
                                 const Rectangle& rect,
                                 ArgbBuffer* const debug_output) {
  if (!unfiltered_pixels.HasAlpha()) return;
  assert(unfiltered_pixels.GetWidth() == filtered_pixels.GetWidth());
  assert(unfiltered_pixels.GetHeight() == filtered_pixels.GetHeight());
  assert(unfiltered_pixels.GetWidth() >= rect.width);
  assert(unfiltered_pixels.GetHeight() >= rect.height);

  Channel c = kAChannel;
  const Plane16& unfiltered = unfiltered_pixels.GetChannel(c);
  const Plane16& filtered = filtered_pixels.GetChannel(c);
  const int32_t kMultiplier = 4;  // Multiplier to make things more visible.
  for (uint32_t y = 0; y < rect.height; ++y) {
    uint8_t* row = (uint8_t*)debug_output->GetRow(rect.y + y);
    for (uint32_t x = 0; x < rect.width; ++x) {
      const int32_t diff = filtered.At(x, y) - unfiltered.At(x, y);
      const int32_t v = diff * kMultiplier * 255 / (int32_t)kAlphaMax + 128;
      const uint32_t abs_x = rect.x + x;
      row[4 * abs_x + 0] = 255;
      row[4 * abs_x + 1] = row[4 * abs_x + 2] = row[4 * abs_x + 3] =
          Clamp<int32_t>(v, 0, 255);
    }
  }
}

static void ApplyVDebugDiff(const ArgbBuffer& unfiltered_pixels,
                            const ArgbBuffer& filtered_pixels,
                            int32_t multiplier,
                            ArgbBuffer* const debug_output) {
  assert(unfiltered_pixels.format == filtered_pixels.format);
  assert(unfiltered_pixels.width == filtered_pixels.width);
  assert(unfiltered_pixels.height == filtered_pixels.height);
  assert(unfiltered_pixels.format == debug_output->format);
  assert(unfiltered_pixels.width == debug_output->width);
  assert(unfiltered_pixels.height == debug_output->height);
  const uint32_t num_channels = WP2FormatBpp(debug_output->format);
  const uint32_t alpha_channel = WP2AlphaChannelIndex(debug_output->format);

  for (uint32_t y = 0; y < unfiltered_pixels.height; ++y) {
    for (uint32_t x = 0; x < unfiltered_pixels.width; ++x) {
      for (uint32_t c = 0; c < num_channels; ++c) {
        const int32_t unfiltered =
            ((uint8_t*)unfiltered_pixels.GetRow(y))[x * num_channels + c];
        const int32_t filtered =
            ((uint8_t*)filtered_pixels.GetRow(y))[x * num_channels + c];

        int32_t diff = (filtered - unfiltered) * multiplier;
        if (c == alpha_channel) {
          diff = 255 - std::abs(diff);  // Opaque by default.
        } else {
          diff = 128 + diff;  // Grey by default.
        }

        uint8_t* const pixel =
            (uint8_t*)debug_output->GetRow(y) + x * num_channels + c;
        *pixel = (uint8_t)Clamp(diff, 0, 255);
      }
    }
  }
}

//------------------------------------------------------------------------------

void FilterBlockMap::ApplyVDebug(const DecoderConfig& config,
                                 ArgbBuffer* const debug_output) {
  const int vd = VDMatch(config, "bpp") ? -1 : VDChannel(config);

  std::unordered_set<uint32_t> done;  // NOLINT (absl::flat_hash_set)
  for (uint32_t by = 0; by < tile_rect_.height / kMinBlockSizePix; ++by) {
    for (uint32_t bx = 0; bx < tile_rect_.width / kMinBlockSizePix; ++bx) {
      const FilterBlockMap::BlockFeatures& block_features =
          features_[by * max_num_blocks_x_ + bx];
      if (!done.insert(block_features.index).second) continue;

      const uint32_t x = bx * kMinBlockSizePix, y = by * kMinBlockSizePix;
      const uint32_t w =
          std::min(BlockWidthPix(block_features.size), tile_rect_.width - x);
      const uint32_t h =
          std::min(BlockHeightPix(block_features.size), tile_rect_.height - y);

      const uint8_t color =
          (vd == -1) ? block_features.min_bpp : block_features.res_den[vd];
      debug_output->Fill({x, y, w, h}, Argb32b{255, color, 0, 0});
      debug_output->DrawRect({x, y, w, h}, Argb32b{255, 128, 128, 128});

      if (VDSelected(tile_rect_.x, tile_rect_.y, {x, y, w, h}, config)) {
        debug_output->DrawRect({x, y, w, h}, Argb32b{255, 255, 255, 255});

        std::string* const str_ptr = &config.info->selection_info;
        WP2SAppend(str_ptr, "\nBlock at %u, %u px\n\n", x, y);

        WP2SAppend(str_ptr, "Min bits-per-pixel: %f\n",
                   block_features.min_bpp / 128.);
        for (Channel c : {kYChannel, kUChannel, kVChannel, kAChannel}) {
          WP2SAppend(str_ptr, "Residual density (%s): %f\n", kChannelStr[c],
                     block_features.res_den[c] / 255.);
        }
      }
    }
  }
}

//------------------------------------------------------------------------------

void DeblockingFilter::SavePixelsForVDebug(uint32_t num_rows) {
  if (VDMatch(config_, "diff")) {
    SavePixelsForVDebugDiff(*blocks_.pixels_, num_rows, &unfiltered_pixels_,
                            &num_debug_rows_);
  } else if (VDMatch(config_, "strength")) {
    if (num_rows > 0 && num_debug_rows_ == 0) {
      // Registered pixels won't fill the whole canvas so fill it with black.
      config_.info->debug_output.Fill(
          blocks_.tile_rect_, Argb32b{0xff, 0x00, 0x00, 0x00});
    }
    num_debug_rows_ = num_rows;
  }
}

void DeblockingFilter::RegisterPixelsForVDebug(
    const DecoderConfig& config, uint32_t q0_x, uint32_t q0_y,
    uint32_t max_half, uint32_t filtered_half, bool horizontal_filtering,
    Channel channel, uint32_t strength, uint32_t sharpness,
    bool before_deblocking, bool deblocked, uint32_t num_bits,
    const int16_t* const q0, uint32_t step) {
  if (VDMatch(config, "strength") && VDChannel(config) == channel &&
      IsHorizontal(config) == horizontal_filtering) {
    ArgbBuffer* const debug_output = &config.info->debug_output;

    if (!before_deblocking) {
      // Strength heatmap on R, inverted sharpness heatmap on G.
      const uint32_t r = DivRound(strength * 255, kDblkMaxStrength);
      const uint32_t g =
          DivRound((kDblkMaxSharpness - sharpness) * 255, kDblkMaxSharpness);

      // Display gradients in filtered direction.
      for (uint32_t i = filtered_half; i >= 1; --i) {
        debug_output->Fill(
            horizontal_filtering ? Rectangle{q0_x - i, q0_y, i * 2, 1}
                                 : Rectangle{q0_x, q0_y - i, 1, i * 2},
            deblocked ? Argb32b{255, (uint8_t)(r / i), (uint8_t)(g / i), 0}
                      : Argb32b{255, (uint8_t)(200 / i), (uint8_t)(200 / i),
                                (uint8_t)(255 / i)});
      }
    }

    const Rectangle rect =
        horizontal_filtering
            ? Rectangle{q0_x - filtered_half, q0_y, 2 * filtered_half, 1}
            : Rectangle{q0_x, q0_y - filtered_half, 1, 2 * filtered_half};
    if (VDSelected(rect, config)) {
      debug_output->DrawRect(rect, Argb32b{255, 255, 255, 255});

      std::string* const str_ptr = &config.info->selection_info;
      if (before_deblocking) {
        WP2SAppend(str_ptr, "\nEdge near %u, %u px (%s)\n\n", q0_x, q0_y,
                   kChannelStr[channel]);

        WP2SAppend(str_ptr, "  Half-length: %u\n", filtered_half);
        WP2SAppend(str_ptr, "  Strength: %u\n", strength);
        WP2SAppend(str_ptr, "  Sharpness: %u\n", sharpness);
        WP2SAppend(str_ptr, "  Deblock threshold: %d\n\n",
                   DeblockThresholdFromSharpness(sharpness, num_bits));

        for (int i = -(int32_t)max_half; i < (int32_t)max_half; ++i) {
          if (i == 0) WP2SAppend(str_ptr, " |");
          WP2SAppend(str_ptr, " %4d", q0[i * (int32_t)step]);
        }
        WP2SAppend(str_ptr, "\n");
      } else {  // !before_deblocking
        for (int i = -(int32_t)max_half; i < (int32_t)max_half; ++i) {
          if (i == 0) WP2SAppend(str_ptr, " |");
          WP2SAppend(str_ptr, " %4d", q0[i * (int32_t)step]);
        }
        WP2SAppend(str_ptr, "\n\n");
        WP2SAppend(str_ptr, "  Deblocked: %s\n", deblocked ? "yes" : "no");
      }
    }
  }
}

void DeblockingFilter::ApplyVDebug(uint32_t num_rows) {
  if (VDMatch(config_, "diff") && num_rows == blocks_.tile_rect_.height) {
    if (VDMatch(config_, "a")) {
      ApplyVDebugAlphaDiff(unfiltered_pixels_, *blocks_.pixels_,
                           blocks_.tile_rect_, &config_.info->debug_output);
    } else {
      const int32_t yuv_max_value =
          (1 << (blocks_.num_precision_bits_ - 1)) - 1;
      ApplyVDebugDiff(unfiltered_pixels_, *blocks_.pixels_, yuv_max_value, 16,
                      blocks_.tile_rect_, &config_.info->debug_output);
    }
  }
}

//------------------------------------------------------------------------------

void GrainFilter::SavePixelsForVDebug(uint32_t num_rows) {
  if (VDMatch(config_, "diff")) {
    SavePixelsForVDebugDiff(*blocks_.pixels_, num_rows, &unfiltered_pixels_,
                            &num_debug_rows_);
  }
}

void GrainFilter::ApplyVDebug(uint32_t num_rows) {
  if (VDMatch(config_, "diff") && num_rows == blocks_.tile_rect_.height) {
    const int32_t yuv_max_value = (1 << (blocks_.num_precision_bits_ - 1)) - 1;
    ApplyVDebugDiff(unfiltered_pixels_, *blocks_.pixels_, yuv_max_value, 2,
                    blocks_.tile_rect_, &config_.info->debug_output);
  }
}

//------------------------------------------------------------------------------

void DirectionalFilter::SavePixelsForVDebug(uint32_t num_rows) {
  if (VDMatch(config_, "diff")) {
    SavePixelsForVDebugDiff(*blocks_.pixels_, num_rows, &unfiltered_pixels_,
                            &num_debug_rows_);
  }
}

void DirectionalFilter::RegisterPixelsForVDebug(uint32_t from_x, uint32_t to_x,
                                                uint32_t from_y, uint32_t to_y,
                                                uint32_t primary_strength,
                                                uint32_t secondary_strength,
                                                uint32_t direction,
                                                uint32_t variance) {
  ArgbBuffer* const debug_output = &config_.info->debug_output;
  const Rectangle rect{blocks_.tile_rect_.x + from_x,
                       blocks_.tile_rect_.y + from_y,
                       to_x - from_x + 1,
                       to_y - from_y + 1};

  if (VDMatch(config_, "strength")) {
    const uint8_t primary_color =
        (uint8_t)DivRound(primary_strength * 255, kMaxPriStr);
    const uint8_t secondary_color =
        (uint8_t)DivRound(secondary_strength * 255, kMaxSecStr);

    debug_output->Fill(
        rect, Argb32b{255, primary_color, secondary_color, 0});
  } else if (VDMatch(config_, "variance")) {
    const uint8_t color =
        (uint8_t)Clamp(DivRound(WP2Log2Floor(variance >> 6) * 255, 16), 0, 255);

    debug_output->Fill(rect, Argb32b{255, 0, 0, color});
  } else if (VDMatch(config_, "direction")) {
    assert(direction < kDrctFltNumDirs);
    // Light checkerboard to distinguish areas where the direction is computed.
    debug_output->Fill(rect, (((from_x + from_y) / kDrctFltSize) & 1u)
                                     ? Argb32b{255, 40, 40, 40}
                                     : Argb32b{255, 0, 0, 0});

    // Draw the oriented tap pattern in the middle of the block.
    for (int32_t dist = -(int32_t)kDrctFltTapDist;
         dist <= (int32_t)kDrctFltTapDist; ++dist) {
      const int32_t cx =
          blocks_.tile_rect_.x + from_x + (int32_t)kDrctFltSize / 2;
      const int32_t cy =
          blocks_.tile_rect_.y + from_y + (int32_t)kDrctFltSize / 2;
      MarkPixelForVDebug(cx, cy, {255, 255, 255, 255}, debug_output);
      if (dist != 0) {
        const int32_t s = (dist < 0) ? -1 : 1, d = std::abs(dist) - 1;
        // Primary taps.
        MarkPixelForVDebug(cx + s * kDrctFltTapPos[direction][d][0],
                           cy + s * kDrctFltTapPos[direction][d][1],
                           {255, 255, 255, 127}, debug_output);
        // Secondary taps.
        MarkPixelForVDebug(
            cx + s * kDrctFltTapPos[(direction + 2u) % kDrctFltNumDirs][d][0],
            cy + s * kDrctFltTapPos[(direction + 2u) % kDrctFltNumDirs][d][1],
            {255, 127, 127, 255}, debug_output);
        MarkPixelForVDebug(
            cx + s * kDrctFltTapPos[(direction + 6u) % kDrctFltNumDirs][d][0],
            cy + s * kDrctFltTapPos[(direction + 6u) % kDrctFltNumDirs][d][1],
            {255, 127, 127, 255}, debug_output);
      }
    }
  }
}

void DirectionalFilter::ApplyVDebug(uint32_t num_rows) {
  if (VDMatch(config_, "diff") && num_rows == blocks_.tile_rect_.height) {
    const int32_t yuv_max_value = (1 << (blocks_.num_precision_bits_ - 1)) - 1;
    ApplyVDebugDiff(unfiltered_pixels_, *blocks_.pixels_, yuv_max_value, 16,
                    blocks_.tile_rect_, &config_.info->debug_output);
  }
}

//------------------------------------------------------------------------------

void RestorationFilter::SavePixelsForVDebug(uint32_t num_rows) {
  if (VDMatch(config_, "diff")) {
    SavePixelsForVDebugDiff(*blocks_.pixels_, num_rows, &unfiltered_pixels_,
                            &num_debug_rows_);
  }
}

void RestorationFilter::RegisterPixelsForVDebug(uint32_t from_x, uint32_t to_x,
                                                uint32_t from_y,
                                                uint32_t to_y) {
  if (VDMatch(config_, "strength")) {
    ArgbBuffer* const debug_output = &config_.info->debug_output;
    const Rectangle rect{blocks_.tile_rect_.x + from_x,
                         blocks_.tile_rect_.y + from_y,
                         to_x - from_x + 1,
                         to_y - from_y + 1};
    debug_output->Fill(rect, Argb32b{255, 128, 128, 128});

    if (VDSelected(rect, config_)) {
      MinMaxAvg min_max_avg_strength;
      for (uint32_t y = from_y; y <= to_y; ++y) {
        for (uint32_t x = from_x; x <= to_x; ++x) {
          const uint8_t strength =
              filter_strength_map_[(y - from_y) * blocks_.tile_rect_.width +
                                   (x - from_x)];
          const uint8_t color = (uint8_t)DivRound(strength * 255, 63);
          MarkPixelForVDebug(blocks_.tile_rect_.x + x, blocks_.tile_rect_.y + y,
                             Argb32b{255, color, 0, 0}, debug_output);
          min_max_avg_strength.Add(strength / 63.);
        }
      }

      // Highlight the border.
      debug_output->Fill(rect, Argb32b{255, 255, 255, 255});
      std::string* const str_ptr = &config_.info->selection_info;
      WP2SAppend(str_ptr, "\nStrength %s\n\n",
                 min_max_avg_strength.ToString().c_str());

      const uint32_t area = from_y / kWieFltHeight;
      for (Channel c : {kYChannel, kUChannel, kVChannel}) {
        WP2SAppend(str_ptr, "Channel %s\n", kChannelStr[c]);
        for (uint32_t h_or_v : {0, 1}) {
          int32_t tap_weights[kWieFltNumTaps];
          std::copy(params_.half_tap_weights[c][area][h_or_v],
                    params_.half_tap_weights[c][area][h_or_v] + kWieFltTapDist,
                    tap_weights);
          WienerHalfToFullWgts(tap_weights, tap_weights);
          WP2SAppend(str_ptr, "  %s",
                     (h_or_v == 0) ? "Horizontal:" : "Vertical:  ");
          for (int32_t weight : tap_weights) {
            WP2SAppend(str_ptr, " %6.3f",
                       (double)weight / (1u << (kWieFltNumBitsTapWgts - 1)));
          }
          *str_ptr += '\n';
        }
      }
    }
  }
}

void RestorationFilter::ApplyVDebug(uint32_t num_rows) {
  if (VDMatch(config_, "diff") && num_rows == blocks_.tile_rect_.height) {
    const int32_t yuv_max_value = (1 << (blocks_.num_precision_bits_ - 1)) - 1;
    ApplyVDebugDiff(unfiltered_pixels_, *blocks_.pixels_, yuv_max_value, 16,
                    blocks_.tile_rect_, &config_.info->debug_output);
  }
}

//------------------------------------------------------------------------------

void IntertileFilter::SavePixelsForVDebug(uint32_t num_rows) {
  if (VDMatch(*config_, "diff")) {
    SavePixelsForVDebugDiff(canvas_, num_rows, &unfiltered_pixels_,
                            &num_debug_rows_);
  } else if (VDMatch(*config_, "strength")) {
    if (num_rows > 0 && num_debug_rows_ == 0) {
      // Registered pixels won't fill the whole canvas so fill it with black.
      config_->info->debug_output.Fill(Argb32b{0xff, 0x00, 0x00, 0x00});
    }
    num_debug_rows_ = num_rows;
  }
}

void IntertileFilter::ApplyVDebug(uint32_t num_rows) {
  if (VDMatch(*config_, "diff") && (num_rows == canvas_.Y.h_)) {
    const FilterBlockMap& blocks = tiles_layout_->GetTileAt(0, 0).block_map;
    const int32_t yuv_max_value = (1 << (blocks.num_precision_bits_ - 1)) - 1;
    ApplyVDebugDiff(unfiltered_pixels_, canvas_, yuv_max_value, 16,
                    {0, 0, canvas_.Y.w_, canvas_.Y.h_},
                    &config_->info->debug_output);
  }
}

//------------------------------------------------------------------------------

void AlphaFilter::SavePixelsForVDebug(uint32_t num_rows) {
  if (VDMatch(config_, "diff")) {
    SavePixelsForVDebugDiff(*blocks_.pixels_, num_rows, &unfiltered_pixels_,
                            &num_debug_rows_);
  }
}

void AlphaFilter::ApplyVDebug(uint32_t num_rows) {
  if (VDMatch(config_, "diff") && num_rows == blocks_.tile_rect_.height) {
    ApplyVDebugAlphaDiff(unfiltered_pixels_, *blocks_.pixels_,
                         blocks_.tile_rect_, &config_.info->debug_output);
  }
}

//------------------------------------------------------------------------------

// Sets all R, G and B channels to the value of a single 'channel'.
static void ToGray(const ArgbBuffer& from, uint32_t channel,
                   ArgbBuffer* const to) {
  assert(from.width == to->width && from.height == to->height);
  assert(from.format == WP2_Argb_32 && to->format == WP2_Argb_32);
  uint32_t kAlpha = 0;
  for (uint32_t y = 0; y < to->height; ++y) {
    const uint8_t* const from_row = (const uint8_t*)from.GetRow(y);
    uint8_t* const to_row = (uint8_t*)to->GetRow(y);
    for (size_t x = 0; x < to->width; ++x) {
      const uint8_t* const from_pixel = from_row + x * 4;
      uint8_t* const to_pixel = to_row + x * 4;
      uint32_t gray = from_pixel[channel];
      if (channel != kAlpha) {  // Unpremultiply by alpha.
        uint8_t alpha = from_pixel[kAlpha];
        assert(gray <= alpha);
        gray = (alpha > 0) ? ((gray * 255 + alpha / 2) / alpha) : 0;
      }
      to_pixel[0] = 255;
      to_pixel[1] = to_pixel[2] = to_pixel[3] = (uint8_t)gray;
    }
  }
}

// Handles VisualDebug for "compressed" and "original". Declared in wp2_dec_i.h.
void ApplyVDebugBeforeAfter(const DecoderConfig& config,
                            const CSPTransform& csp_tranform,
                            const Tile& tile,
                            ArgbBuffer* const debug_output) {
  if (VDMatch(config, "original/diff")) {
    ArgbBuffer argb;
    WP2_ASSERT_STATUS(argb.Resize(debug_output->width, debug_output->height));
    YUVPlane non_padded;
    WP2_ASSERT_STATUS(non_padded.SetView(
        tile.yuv_output, {0, 0, tile.rect.width, tile.rect.height}));
    WP2_ASSERT_STATUS(
        non_padded.Export(csp_tranform, /*resize_if_needed=*/false, &argb));
    ApplyVDebugDiff(*debug_output, argb, 2, debug_output);
  } else if (VDMatch(config, "y") || VDMatch(config, "u") ||
             VDMatch(config, "v")) {
    const uint32_t yuv_bits = csp_tranform.GetYUVPrecisionBits() + 1;
    if (VDMatch(config, "original")) {
      csp_tranform.Apply(debug_output, yuv_bits, 0, 0, debug_output->width,
                         debug_output->height);
      ToGray(*debug_output, VDChannel(config) + 1, debug_output);
    } else if (VDMatch(config, "compressed")) {
      const Plane16& plane = tile.yuv_output.GetChannel(VDChannel(config));
      WP2_ASSERT_STATUS(plane.ToGray(debug_output, yuv_bits));
    }
  } else {
    const bool a = VDMatch(config, "a"), r = VDMatch(config, "r"),
               g = VDMatch(config, "g"), b = VDMatch(config, "b");
    if (a || r || g || b) {
      const uint32_t channel = a ? 0 : r ? 1 : g ? 2 : 3;
      if (VDMatch(config, "original")) {
        ToGray(*debug_output, channel, debug_output);
      } else if (VDMatch(config, "compressed")) {
        ArgbBuffer argb;
        WP2_ASSERT_STATUS(
            argb.Resize(debug_output->width, debug_output->height));
        YUVPlane non_padded_yuv;
        WP2_ASSERT_STATUS(non_padded_yuv.SetView(
            tile.yuv_output, {0, 0, tile.rect.width, tile.rect.height}));
        WP2_ASSERT_STATUS(non_padded_yuv.Export(
            csp_tranform, /*resize_if_needed=*/false, &argb));
        ToGray(argb, channel, debug_output);
      }
    }
  }

  if (VDMatch(config, "histogram") && VDSelected(tile.rect, config)) {
    // A gray image of either Y,U,V,R,G or B should already be in
    // 'debug_output'. The histogram will be based on that.
    ArgbBuffer selected_area(debug_output->format);
    Rectangle rect = config.info->selection.ClipWith(tile.rect);
    assert(rect.GetArea() != 0);
    rect.x -= tile.rect.x;
    rect.y -= tile.rect.y;
    WP2_ASSERT_STATUS(selected_area.SetView(*debug_output, rect));

    std::array<uint32_t, 32> histogram{0};
    uint8_t min_value, max_value;
    const uint32_t num_buckets =
        GetHistogram(selected_area, /*channel=*/1, histogram.size(),
                     histogram.data(), &min_value, &max_value);

    debug_output->DrawRect(rect, Argb32b{255, 255, 128, 128});
    std::string* const str_ptr = &config.info->selection_info;
    WP2SAppend(str_ptr, "\nArea %u x %u at %u, %u\n", rect.width, rect.height,
               rect.x + tile.rect.x, rect.y + tile.rect.y);
    WP2SAppend(str_ptr, "\nHistogram [%u, %u]\n", min_value, max_value);
    for (uint32_t bucket = 0; bucket < num_buckets; ++bucket) {
      const float percent = 100.f * histogram[bucket] / rect.GetArea();
      WP2SAppend(str_ptr, "  Bucket %2u: %5.2f%% (%8u)  ", bucket, percent,
                 histogram[bucket]);
      *str_ptr += std::string(std::lround(percent), '|');
      WP2SAppend(str_ptr, "\n");
    }

    std::array<uint32_t, histogram.size()> clusters{0};
    std::array<uint32_t, kMaxNumSegments> centers{0};
    const uint32_t num_clusters = ClusterHistogram(
        histogram.data(), num_buckets, /*max_clusters=*/kMaxNumSegments,
        clusters.data(), centers.data());
    WP2SAppend(str_ptr, "\nNum clusters: %u\n", num_clusters);
  }
}

//------------------------------------------------------------------------------

// Fills the background of the current tile with luma and display some info.
WP2Status PartitionScoreFunc::ClearVDebug() const {
  if (VDMatch(*config_, "encoder/partition")) {
    ArgbBuffer debug_output;
    WP2_CHECK_STATUS(
        debug_output.SetView(config_->info->debug_output, tile_rect_));
    if (VDMatch(*config_, "method")) {
      const PartitionMethod pm =
          FinalPartitionMethod(*config_, tile_rect_.width, tile_rect_.height);
      if (VDSelected(tile_rect_, *config_)) {
        debug_output.Fill(Argb32b{255, 255, 255, 255});
        std::string* const str_ptr = &config_->info->selection_info;
        WP2SAppend(str_ptr, "\nTile at %4u, %4u (%2u x %2u px): %s\n",
                   tile_rect_.x, tile_rect_.y, tile_rect_.width,
                   tile_rect_.height, kPartitionMethodString[pm]);
      } else {
        WP2_CHECK_STATUS(src_->Y.ToGray(&debug_output, kMaxYuvBits));
      }
      debug_output.DrawRect(
          {0, 0, debug_output.width, debug_output.height},
          Argb32b{255, (uint8_t)(pm * 255 / NUM_PARTITION_METHODS),
                  (uint8_t)(100 + pm * 100 / NUM_PARTITION_METHODS),
                  (uint8_t)(50 + pm * 50 / NUM_PARTITION_METHODS)});
    } else if (VDMatch(*config_, "score")) {
      WP2_CHECK_STATUS(src_->Y.ToGray(&debug_output, kMaxYuvBits));
    }
  }
  return WP2_STATUS_OK;
}

// Should be called after ClearVDebug() and for each block score computed.
bool PartitionScoreFunc::RegisterScoreForVDebug(const Block& block, float score,
                                                bool ending_new_line) const {
  if (!VDMatch(*config_, "encoder/partition/score") ||
      !VDSelected(
          tile_rect_.x, tile_rect_.y,
          {block.x_pix(), block.y_pix(), kMinBlockSizePix, kMinBlockSizePix},
          *config_)) {
    return false;
  }
  // Highlight the block with a color dependent on its size.
  const uint8_t color_w = (uint8_t)(191 + block.w() * 64 / kMaxBlockSize);
  const uint8_t color_h = (uint8_t)(191 + block.h() * 64 / kMaxBlockSize);
  config_->info->debug_output.Fill(
      {tile_rect_.x + block.x_pix() + 1, tile_rect_.y + block.y_pix() + 1,
       block.w_pix() - 2, block.h_pix() - 2},
      Argb32b{255, color_w, color_h, 0});

  std::string* const str_ptr = &config_->info->selection_info;
  WP2SAppend(str_ptr, "Block at %4u, %4u (%2u x %2u px), score: %9.3f",
             tile_rect_.x + block.x_pix(), tile_rect_.y + block.y_pix(),
             block.w_pix(), block.h_pix(), score);
  if (ending_new_line) WP2SAppend(str_ptr, "\n");
  return true;
}

// Should be called after ClearVDebug() and for each block actually encoded.
void Partitioner::RegisterOrderForVDebug(uint32_t pass, uint32_t pass_index,
                                         const Block& block,
                                         uint32_t block_index,
                                         uint32_t max_num_blocks) const {
  if (VDMatch(*config_, "encoder/partition")) {
    ArgbBuffer* const debug_output = &config_->info->debug_output;
    Rectangle rect;
    rect.x = tile_rect_.x + block.x_pix();
    rect.y = tile_rect_.y + block.y_pix();
    rect.width = std::min(block.w_pix(), debug_output->width - rect.x);
    rect.height = std::min(block.h_pix(), debug_output->height - rect.y);

    if (VDMatch(*config_, "pass") || VDMatch(*config_, "order")) {
      Argb32b color;
      if (VDMatch(*config_, "pass")) {
        color.a = 255;
        color.r = color.g = color.b =
            (uint8_t)(255 - Clamp(pass_index * 12u, 0u, 255u));
      } else {
        color.a = 255;
        color.r = color.g = color.b =
            (uint8_t)(55u + DivRound(200u * block_index, max_num_blocks));
      }
      debug_output->Fill(rect, color);
      uint32_t selected_pass_index;
      if (VDMatch(*config_, "pass")) {
        if (VDIndex(config_->info->visual_debug, &selected_pass_index) &&
            pass_index == selected_pass_index) {
          debug_output->DrawRect(rect, Argb32b{255, 0, 0, 255});
        } else if (VDMatch(*config_, "all")) {
          const Argb32b kPassColor[] = {{255, 141, 202, 175},
                                        {255, 140, 73, 192},
                                        {255, 141, 208, 85},
                                        {255, 205, 104, 150},
                                        {255, 133, 148, 196}};
          STATIC_ASSERT_ARRAY_SIZE(kPassColor,
                                   (int)MultiScoreFunc::Pass::Any + 1);
          debug_output->DrawRect(rect, kPassColor[pass]);
        }
      }

      if (VDSelected(rect, *config_)) {
        // Highlight the border.
        debug_output->DrawRect(rect, Argb32b{255, 255, 0, 0});

        std::string* const str_ptr = &config_->info->selection_info;
        WP2SAppend(
            str_ptr, "Block at %u, %u (%u x %u px), pass: %u, order: %u\n",
            rect.x, rect.y, rect.width, rect.height, pass_index, block_index);
      }
    } else if (VDMatch(*config_, "score")) {
      // To have a reference, borders of actual encoded blocks are drawn.
      // Selected blocks are colored in RegisterScoreForVDebug().
      debug_output->DrawRect(rect, Argb32b{255, 64, 128, 64});
    }
  }
}

bool BlockScoreFunc::RegisterScoreForVDebug(
    const Block blocks[], uint32_t num_blocks, const float rate[4],
    const float disto[4], float total_rate, float total_disto, float score) {
  Rectangle hull = blocks[0].rect();  // Print surrounding block if possible.
  for (uint32_t i = 1; i < num_blocks; ++i) {
    hull = hull.MergeWith(blocks[i].rect());
  }
  const Block hull_block(hull.x, hull.y, GetBlockSize(hull.width, hull.height));
  if (!PartitionScoreFunc::RegisterScoreForVDebug(hull_block, score,
                                                  /*ending_new_line=*/false)) {
    return false;
  }

  std::string* const str_ptr = &config_->info->selection_info;
  WP2SAppend(str_ptr, "  rate %8.2f disto %8.2f", total_rate, total_disto);
  if (num_blocks > 1) {  // Print sub-blocks stats.
    WP2SAppend(str_ptr, "  %u blks: ", num_blocks);
    for (uint32_t i = 0; i < num_blocks; ++i) {
      WP2SAppend(str_ptr, " %ux%u(r%.2f,d%.2f)", blocks[i].w_pix(),
                 blocks[i].h_pix(), rate[i], disto[i]);
    }
  }
  WP2SAppend(str_ptr, "\n");
  return true;
}

bool AreaScoreFunc::RegisterScoreForVDebug(
    BlockSize grid_size, const Vector<CodedBlock>& area_blocks, float score,
    float disto, float rate) const {
  // Only display something when the current 'area_' is selected.
  if (!VDMatch(*config_, "encoder/partition/score") ||
      !VDSelected(tile_rect_.x, tile_rect_.y, area_, *config_)) {
    return false;
  }

  std::string* const str_ptr = &config_->info->selection_info;
  if (grid_size == BLK_LAST) {  // Meaning default partitioning (happens once).
    WP2SAppend(str_ptr, "\nArea at %4u, %4u (%2u x %2u px)\n",
               tile_rect_.x + area_.x, tile_rect_.y + area_.y,
               area_.width, area_.height);
    WP2SAppend(str_ptr, "\nscore: %9.3f = disto %7.3f + rate %7.3f\n", score,
               disto, rate);
    WP2SAppend(str_ptr, "  (default)\n");
    ArgbBuffer debug_output;
    WP2_ASSERT_STATUS(
        debug_output.SetView(config_->info->debug_output, tile_rect_));
    WP2_ASSERT_STATUS(buffer_.Export(
        gparams_->transf_, /*resize_if_needed=*/false, &debug_output));
    debug_output.DrawRect(area_, Argb32b{255, 255, 0, 0});
    for (const CodedBlock& cb : area_blocks) {
      debug_output.DrawRect(cb.blk().rect_pix(), Argb32b{255, 128, 255, 0});
    }
  } else {  // Grid partitioning (happens for each of several block sizes).
    WP2SAppend(str_ptr, "\nscore: %9.3f = disto %7.3f + rate %7.3f\n", score,
               disto, rate);
    WP2SAppend(str_ptr, "  (grid of %2u x %2u px)\n", BlockWidthPix(grid_size),
               BlockHeightPix(grid_size));
  }
  return true;
}

bool SubAreaScoreFunc::RegisterScoreForVDebug(
    const Block& block, const Vector<CodedBlock>& area_remaining_blocks,
    float score, float disto, float rate) const {
  if (PartitionScoreFunc::RegisterScoreForVDebug(block, score)) {
    std::string* const str_ptr = &config_->info->selection_info;
    WP2SAppend(str_ptr, "  disto %7.3f + rate %7.3f%s\n", disto, rate,
               (default_block_.dim() == BLK_LAST) ? " (default)" : "");
    return true;
  }
  return false;
}

bool TileScoreFunc::RegisterScoreForVDebug(const char label[],
                                           const Block& block,
                                           float score) const {
  if (VDMatch(*config_, "encoder/partition/score")) {
    std::string* const str_ptr = &config_->info->selection_info;
    WP2SAppend(str_ptr, "\n%s: score %.6f, %4zu blocks, %4zu-byte tile", label,
               score, blocks_.size(),
               enc_tiles_layout_.tiles.front().enc.BufferSize());
    if (block.dim() != BLK_LAST) {
      WP2SAppend(str_ptr, ", on block at %3u,%3u (%2ux%2u)", block.x_pix(),
                 block.y_pix(), block.w_pix(), block.h_pix());
    }
    WP2SAppend(str_ptr, ", tile distortion %.4f dB", distortion_[4]);
    return true;
  }
  return false;
}

//------------------------------------------------------------------------------

void ExhaustivePartitioner::RegisterScoreForVDebug(
    float best_partition_score, size_t best_partition_size,
    size_t num_iterations) const {
  if (VDMatch(*config_, "encoder/partition/score")) {
    std::string* const str_ptr = &config_->info->selection_info;
    WP2SAppend(str_ptr, "\n\nNum iterations: %zu\n", num_iterations);
    WP2SAppend(str_ptr, "Best partition score: %8.3f size: %4zu\n",
               best_partition_score, best_partition_size);
  }
}

//------------------------------------------------------------------------------

void MultiPassPartitioner::RegisterPassForVDebug(
    MultiScoreFunc::Pass pass, BlockSize block_size,
    uint32_t num_chosen_blocks) const {
  if (!VDMatch(*config_, "encoder/partition/pass")) return;
  if (!VDSelected(tile_rect_, *config_)) return;

  const char* const kPassStr[] = {"FlatLumaAlpha", "NarrowStdDev",
                                  "GoodQuantDCT", "Direction", "Any"};
  STATIC_ASSERT_ARRAY_SIZE(kPassStr, (int)MultiScoreFunc::Pass::Any + 1);

  std::string* const str_ptr = &config_->info->selection_info;
  WP2SAppend(str_ptr, "Pass: %2u - %s - %2ux%2u - %4u blocks\n", pass_index_,
             kPassStr[(int)pass], BlockWidthPix(block_size),
             BlockHeightPix(block_size), num_chosen_blocks);
}

//------------------------------------------------------------------------------

// Displays the normalized value of the selected VDebug menu.
void MultiScoreFunc::DrawValueThresholdColor(
    uint32_t cell_width, uint32_t cell_height, const Block& block,
    ArgbBuffer* const debug_output) const {
  float value = 0.f;
  if (VDMatch(*config_, "luma-alpha-gradient")) {
    value = GetLumaAlphaGradient(block) / GetLumaAlphaGradientThreshold(block);
  } else if (VDMatch(*config_, "narrow-std-dev")) {
    value = GetStdDevRange(block) / GetStdDevRangeThreshold(block);
  } else if (VDMatch(*config_, "quant-dct")) {
    value = GetQuantDCT(block) / GetQuantDCTThreshold(block);
  } else if (VDMatch(*config_, "direction")) {
    value = GetDirection(block) / GetDirectionThreshold(block);
  } else {
    assert(false);
  }
  // Display as blue if passing threshold, otherwise as red. Dark blue means
  // close to the threshold, dark red means far from it.
  const float color[4] = {255.f, (value <= 1.f) ? 0.f : 55.f + 200.f / value,
                          0.f,
                          (value <= 1.f) ? 55.f + 200.f * (1.f - value) : 0.f};

  for (uint32_t sub_y = 0; sub_y < cell_height; ++sub_y) {
    if (block.y_pix() + sub_y >= debug_output->height) break;
    uint8_t* const row = (uint8_t*)debug_output->GetRow(block.y_pix() + sub_y);
    for (uint32_t sub_x = 0; sub_x < cell_width; ++sub_x) {
      if (block.x_pix() + sub_x >= debug_output->width) break;
      uint8_t* const pixel = &row[(block.x_pix() + sub_x) * 4];
      for (uint32_t c : {1, 2, 3}) {
        pixel[c] = (uint8_t)std::lround((pixel[c] + color[c]) * 0.5f);
      }
    }
  }
}

// Displays the lossy luma.
WP2Status MultiScoreFunc::DrawLossyLuma(const Block& block,
                                        const Rectangle& block_rect,
                                        bool draw_side,
                                        ArgbBuffer* const debug_output) const {
  Plane16 p;
  WP2_CHECK_STATUS(p.Resize(block.w_pix(), block.h_pix()));
  int32_t coeffs[kMaxBlockSizePix2];
  QuantizeCoeffs(kYChannel, block, coeffs);
  std::copy(coeffs, coeffs + block.w_pix() * block.h_pix(), p.Row(0));
  ArgbBuffer block_debug_output;
  WP2_CHECK_STATUS(block_debug_output.SetView(*debug_output, block_rect));
  WP2_CHECK_STATUS(p.ToGray(&block_debug_output, kMaxYuvBits + 1u));

  if (draw_side) {
    // If there is enough space, also display the original luma.
    if (block.x_pix() >= block.w_pix()) {
      GetCoeffs(kYChannel, block, coeffs);
      std::copy(coeffs, coeffs + block.w_pix() * block.h_pix(), p.Row(0));
      WP2_CHECK_STATUS(block_debug_output.SetView(
          *debug_output, {block.x_pix() - block.w_pix(), block.y_pix(),
                          block.w_pix(), block.h_pix()}));
      WP2_CHECK_STATUS(p.ToGray(&block_debug_output, kMaxYuvBits + 1u));
    }

    // If there is enough space, also display the quantized luma for the
    // sub-blocks, for comparison during the merging pass.
    if (VDMatch(*config_, "merge-flat") &&
        block.x_pix() + 2 * block.w_pix() <= debug_output->width &&
        block.y_pix() + block.h_pix() <= debug_output->height) {
      QuantizeCoeffs(kYChannel, block, BLK_8x8, coeffs);
      std::copy(coeffs, coeffs + block.w_pix() * block.h_pix(), p.Row(0));
      WP2_CHECK_STATUS(block_debug_output.SetView(
          *debug_output, {block.x_pix() + block.w_pix(), block.y_pix(),
                          block.w_pix(), block.h_pix()}));
      WP2_CHECK_STATUS(p.ToGray(&block_debug_output, kMaxYuvBits + 1u));
    }
  }
  return WP2_STATUS_OK;
}

// Displays the selected block border and info.
WP2Status MultiScoreFunc::DrawSelection(const Block& block,
                                        const Rectangle& block_rect,
                                        ArgbBuffer* const debug_output) const {
  Rectangle border = {SafeSub(block_rect.x, 1u), SafeSub(block_rect.y, 1u), 0,
                      0};
  border.width =
      std::min(block_rect.x + block_rect.width + 1u, debug_output->width) -
      border.x;
  border.height =
      std::min(block_rect.y + block_rect.height + 1u, debug_output->height) -
      border.y;
  debug_output->DrawRect(border, Argb32b{255, 255, 255, 0});

  std::string* const str_ptr = &config_->info->selection_info;
  WP2SAppend(str_ptr, "\nBlock at %4u, %4u (%2u x %2u px)\n",
             tile_rect_.x + block.x_pix(), tile_rect_.y + block.y_pix(),
             block.w_pix(), block.h_pix());
  const uint8_t segment_id = AssignSegmentId(
      *config_, *gparams_, {tile_rect_.x, tile_rect_.y, src_->Y.w_, src_->Y.h_},
      block);
  WP2SAppend(str_ptr, "Segment id: %u\n", segment_id);
  WP2SAppend(str_ptr, "YA gradient: %5.3f = %5.2f / %5.2f\n",
             GetLumaAlphaGradient(block) / GetLumaAlphaGradientThreshold(block),
             GetLumaAlphaGradient(block), GetLumaAlphaGradientThreshold(block));
  WP2SAppend(str_ptr, "Std dev:     %5.3f = %5.2f / %5.2f\n",
             GetStdDevRange(block) / GetStdDevRangeThreshold(block),
             GetStdDevRange(block), GetStdDevRangeThreshold(block));
  WP2SAppend(str_ptr, "Quant DCT:   %5.3f = %5.2f / %5.2f\n",
             GetQuantDCT(block) / GetQuantDCTThreshold(block),
             GetQuantDCT(block), GetQuantDCTThreshold(block));
  WP2SAppend(str_ptr, "Direction:   %5.3f = %6.1f / %6.1f\n",
             GetDirection(block) / GetDirectionThreshold(block),
             GetDirection(block), GetDirectionThreshold(block));
  return WP2_STATUS_OK;
}

// Displays the orientations computed by MultiScoreFunc::ComputeDirection().
static void DrawDirection(const Rectangle& rect, uint32_t direction,
                          Argb32b color, ArgbBuffer* const debug_output) {
  // 4x4 mask per direction
  constexpr uint8_t kPatterns[kDrctFltNumDirs][4][4] = {
      {{0, 0, 0, 1}, {0, 0, 1, 0}, {0, 1, 0, 0}, {1, 0, 0, 0}},
      {{0, 0, 0, 0}, {0, 0, 1, 1}, {1, 1, 0, 0}, {0, 0, 0, 0}},
      {{0, 0, 0, 0}, {1, 1, 1, 1}, {1, 1, 1, 1}, {0, 0, 0, 0}},
      {{0, 0, 0, 0}, {1, 1, 0, 0}, {0, 0, 1, 1}, {0, 0, 0, 0}},
      {{1, 0, 0, 0}, {0, 1, 0, 0}, {0, 0, 1, 0}, {0, 0, 0, 1}},
      {{0, 1, 0, 0}, {0, 1, 0, 0}, {0, 0, 1, 0}, {0, 0, 1, 0}},
      {{0, 1, 1, 0}, {0, 1, 1, 0}, {0, 1, 1, 0}, {0, 1, 1, 0}},
      {{0, 0, 1, 0}, {0, 0, 1, 0}, {0, 1, 0, 0}, {0, 1, 0, 0}}};
  assert(direction < 8 && rect.width <= 4 && rect.height <= 4);

  for (uint32_t y = 0; y < rect.height; ++y) {
    for (uint32_t x = 0; x < rect.width; ++x) {
      if (kPatterns[direction][y][x] == 1) {
        ToUInt8(color, (uint8_t*)debug_output->GetRow(rect.y + y) +
                           (rect.x + x) * WP2FormatBpp(debug_output->format));
      }
    }
  }
}

WP2Status MultiScoreFunc::DrawVDebug() const {
  if (!VDMatch(*config_, "encoder/partition/multi")) return WP2_STATUS_OK;
  ArgbBuffer debug_output;
  WP2_CHECK_STATUS(
      debug_output.SetView(config_->info->debug_output, tile_rect_));

  // Image processing output.
  if (VDMatch(*config_, "spread-5x5")) {
    WP2_CHECK_STATUS(spread_.Y.ToGray(&debug_output, kMaxYuvBits + 1u));
    if (VDSelected(0, 0, tile_rect_, *config_)) {
      const uint32_t x = config_->info->selection.x - tile_rect_.x;
      const uint32_t y = config_->info->selection.y - tile_rect_.y;
      debug_output.DrawRect({x, y, 1, 1}, Argb32b{255, 255, 0, 0});
      std::string* const str_ptr = &config_->info->selection_info;
      WP2SAppend(str_ptr, "\nspread-5x5 at %4u, %4u: %4d\n",
                 config_->info->selection.x, config_->info->selection.y,
                 spread_.Y.At(x, y));
    }
    return WP2_STATUS_OK;
  }

  // Display the general orientation of all 4x4 blocks in the grid.
  if (VDMatch(*config_, "direction-4x4")) {
    // Display the luma source as background.
    WP2_CHECK_STATUS(src_->Y.ToGray(&debug_output, kMaxYuvBits + 1u));
    for (uint32_t y = 0; y + 4 <= debug_output.height; y += 4) {
      for (uint32_t x = 0; x + 4 <= debug_output.width; x += 4) {
        const Rectangle block_rect = {x, y,
                                      std::min(4u, debug_output.width - x),
                                      std::min(4u, debug_output.height - y)};
        const uint32_t i = x / 4 + y / 4 * num_block_cols_;
        if (direction_certainty_[i] > 0) {  // Do not draw random directions.
          const uint8_t c =
              Clamp(direction_certainty_[i] * 255u / 3u, 0u, 255u);
          const Argb32b color = {255, c, (uint8_t)(c / 2), (uint8_t)(255 - c)};
          DrawDirection(block_rect, direction_[i], color, &debug_output);
        }
        if (VDSelected(tile_rect_.x, tile_rect_.y, block_rect, *config_)) {
          debug_output.DrawRect(block_rect, Argb32b{255, 0, 255, 0});
          const uint32_t angle_deg =
              (10 - (int32_t)direction_[i]) * 180 / kDrctFltNumDirs;
          std::string* const str_ptr = &config_->info->selection_info;
          WP2SAppend(str_ptr, "Block at %4u, %4u: \n", x, y);
          WP2SAppend(str_ptr, "  angle %u degrees\n", angle_deg);
          WP2SAppend(str_ptr, "  certainty %u/3\n", direction_certainty_[i]);
        }
      }
    }
    return WP2_STATUS_OK;
  }

  // Original luma.
  WP2_CHECK_STATUS(src_->Y.ToGray(&debug_output, kMaxYuvBits + 1u));

  // clang-format off
  Block block(0, 0, VDMatch(*config_, "32x32") ? BLK_32x32 :
                    VDMatch(*config_, "16x16") ? BLK_16x16 :
                    VDMatch(*config_, "8x8")   ? BLK_8x8   : BLK_4x4);
  // clang-format on
  const bool snapped =
      (VDMatch(*config_, "raw-quant-dct") || config_->partition_snapping);
  const uint32_t cell_width = snapped ? block.w_pix() : kMinBlockSizePix;
  const uint32_t cell_height = snapped ? block.h_pix() : kMinBlockSizePix;
  Block selected_block = Block();

  // Display all blocks.
  for (uint32_t y = 0; y + block.h_pix() <= debug_output.height;
       y += cell_height) {
    for (uint32_t x = 0; x + block.w_pix() <= debug_output.width;
         x += cell_width) {
      block.SetXY(x / kMinBlockSizePix, y / kMinBlockSizePix);
      const uint32_t max_width = debug_output.width - block.x_pix();
      const uint32_t max_height = debug_output.height - block.y_pix();
      const Rectangle cell_rect = {block.x_pix(), block.y_pix(),
                                   std::min(cell_width, max_width),
                                   std::min(cell_height, max_height)};
      const Rectangle block_rect = {block.x_pix(), block.y_pix(),
                                    std::min(block.w_pix(), max_width),
                                    std::min(block.h_pix(), max_height)};
      if (VDMatch(*config_, "raw-quant-dct")) {
        WP2_CHECK_STATUS(DrawLossyLuma(block, block_rect, /*draw_side=*/false,
                                       &debug_output));
      } else {
        DrawValueThresholdColor(cell_width, cell_height, block, &debug_output);
      }
      if (VDSelected(
              tile_rect_.x, tile_rect_.y,
              VDMatch(*config_, "raw-quant-dct") ? block_rect : cell_rect,
              *config_)) {
        selected_block = block;
      }
    }
  }

  // Display selected block on top.
  if (selected_block.dim() != BLK_LAST) {
    const uint32_t max_width = debug_output.width - selected_block.x_pix();
    const uint32_t max_height = debug_output.height - selected_block.y_pix();
    const Rectangle block_rect = {selected_block.x_pix(),
                                  selected_block.y_pix(),
                                  std::min(selected_block.w_pix(), max_width),
                                  std::min(selected_block.h_pix(), max_height)};

    if (VDMatch(*config_, "quant-dct") || VDMatch(*config_, "merge-flat")) {
      WP2_CHECK_STATUS(DrawLossyLuma(selected_block, block_rect,
                                     /*draw_side=*/true, &debug_output));
    }
    WP2_CHECK_STATUS(DrawSelection(selected_block, block_rect, &debug_output));
  }
  return WP2_STATUS_OK;
}

//------------------------------------------------------------------------------

void SegmentationGridVDebug(const EncoderConfig& config, const Rectangle& rect,
                            uint32_t tile_pos_x, uint32_t tile_pos_y,
                            float value) {
  const uint8_t gray = (uint8_t)std::lround(value * 10);
  ArgbBuffer* const debug_output = &config.info->debug_output;
  config.info->debug_output.Fill(
      {tile_pos_x + rect.x, tile_pos_y + rect.y,
       std::min(rect.width, debug_output->width - tile_pos_x - rect.x),
       std::min(rect.height, debug_output->height - tile_pos_y - rect.y)},
      Argb32b{255, gray, gray, gray});

  if (VDSelected(tile_pos_x, tile_pos_y, rect, config)) {
    config.info->debug_output.Fill(
        {tile_pos_x + rect.x, tile_pos_y + rect.y, rect.width, rect.height},
        Argb32b{255, 128, 128, 0});

    std::string* const str_ptr = &config.info->selection_info;
    WP2SAppend(str_ptr, "\nRect at %4u, %4u (%2u x %2u px), value: %9.3f\n",
               tile_pos_x + rect.x, tile_pos_y + rect.y, rect.width,
               rect.height, value);
  }
}

void SegmentationBlockVDebug(const EncoderConfig& config, const CodedBlock& cb,
                             uint32_t tile_pos_x, uint32_t tile_pos_y,
                             float score, const Segment& segment) {
  SegmentationGridVDebug(config,
                         {cb.x_pix(), cb.y_pix(), cb.w_pix(), cb.h_pix()},
                         tile_pos_x, tile_pos_y, score);

  config.info->debug_output.Fill(
      {tile_pos_x + cb.x_pix() + cb.w_pix() - 1,
       tile_pos_y + cb.y_pix() + cb.h_pix() - 1, 1, 1},
      kSegmentColors[cb.id_]);

  const Rectangle rect{cb.x_pix(), cb.y_pix(), cb.w_pix(), cb.h_pix()};

  if (VDSelected(tile_pos_x, tile_pos_y, rect, config)) {
    std::string* const str_ptr = &config.info->selection_info;
    WP2SAppend(str_ptr, "Segment %d risk %f\n", cb.id_, segment.risk_);
  }
}

//------------------------------------------------------------------------------

void SyntaxWriter::PutRawPixels(const Block& block, const YUVPlane& pixels,
                                ANSEnc* const enc) {
  enc->PutRValue(block.x(), SizeBlocks(tile_rect_.width), "DEEP_MATCH_X");
  enc->PutRValue(block.y(), SizeBlocks(tile_rect_.height), "DEEP_MATCH_Y");
  enc->PutRValue(block.w(), kMaxBlockSize + 1, "DEEP_MATCH_W");
  enc->PutRValue(block.h(), kMaxBlockSize + 1, "DEEP_MATCH_H");
  for (Channel channel : {kYChannel, kUChannel, kVChannel, kAChannel}) {
    if (pixels.GetChannel(channel).IsEmpty()) continue;
    for (uint32_t y = 0; y < block.h_pix(); ++y) {
      for (uint32_t x = 0; x < block.w_pix(); ++x) {
        enc->PutSUValue(pixels.GetChannel(channel).At(x, y),
                        gparams_->transf_.GetYUVPrecisionBits() + 1,
                        kChannelStr[channel]);
      }
    }
  }
}

void SyntaxReader::ReadAndCompareRawPixels(const Block& block,
                                           const YUVPlane& pixels,
                                           ANSDec* const dec) {
  const uint32_t block_x = dec->ReadRValue(SizeBlocks(width_), "DEEP_MATCH_X");
  const uint32_t block_y = dec->ReadRValue(SizeBlocks(height_), "DEEP_MATCH_Y");
  const uint32_t block_w = dec->ReadRValue(kMaxBlockSize + 1, "DEEP_MATCH_W");
  const uint32_t block_h = dec->ReadRValue(kMaxBlockSize + 1, "DEEP_MATCH_H");
  if (block_x != block.x() || block_y != block.y()) assert(false);
  if (block_w != block.w() || block_h != block.h()) assert(false);
  for (Channel channel : {kYChannel, kUChannel, kVChannel, kAChannel}) {
    if (pixels.GetChannel(channel).IsEmpty()) continue;
    for (uint32_t y = 0; y < block.h_pix(); ++y) {
      for (uint32_t x = 0; x < block.w_pix(); ++x) {
        const int32_t px = dec->ReadSUValue(
            gparams_->transf_.GetYUVPrecisionBits() + 1, kChannelStr[channel]);
        if (px != pixels.GetChannel(channel).At(x, y)) assert(false);
      }
    }
  }
}

//------------------------------------------------------------------------------

}  // namespace WP2

namespace WP2L {

static bool VDSelectedLossless(uint32_t pos, uint32_t length,
                               const WP2::DecoderInfo& decoder_info,
                               const WP2::Rectangle& tile_rect) {
  const WP2::Rectangle& selection = decoder_info.selection;
  if (selection.width == 0 || selection.height == 0) return false;

  const uint32_t x = selection.x;
  const uint32_t y = selection.y;
  if (!tile_rect.Contains(x, y)) return false;
  const uint32_t selected_pos =
      (y - tile_rect.y) * tile_rect.width + (x - tile_rect.x);
  return (selected_pos >= pos && selected_pos < pos + length);
}

static void Fill(uint32_t pos, uint32_t length, const WP2::Rectangle& tile_rect,
                 WP2::Argb32b color, WP2::ArgbBuffer* out) {
  const uint32_t tile_width = tile_rect.width;
  const uint32_t start_x = tile_rect.x + pos % tile_width;
  const uint32_t start_y = tile_rect.y + pos / tile_width;
  const uint32_t end_y = tile_rect.y + (pos + length) / tile_width;
  if (start_y == end_y) {
    out->Fill({start_x, start_y, length, 1}, color);
    return;
  }
  out->Fill({start_x, start_y, tile_width - pos % tile_width, 1}, color);
  if (end_y - start_y > 1) {
    out->Fill({tile_rect.x, start_y + 1, tile_width, end_y - start_y - 1},
              color);
  }
  out->Fill({tile_rect.x, end_y, (pos + length) % tile_width, 1}, color);
}

static const char* const kSymbolTypeNames[] = {"Literal", "Copy", "Cache"};
STATIC_ASSERT_ARRAY_SIZE(kSymbolTypeNames, kSymbolTypeNum);
static const char* const kGroup4ModeNames[] = {
    "Group4 Vertical", "Group4 Horizontal", "Group4 Pass"};
static const WP2::Argb32b kSymbolColors[2][3] = {
    {{255, 255, 0, 0}, {255, 0, 255, 0}, {255, 0, 0, 255}},       // SymbolType
    {{255, 255, 255, 0}, {255, 0, 255, 255}, {255, 255, 0, 255}}  // Group4
};

void Decoder::RegisterSymbolForVDebug(int symbol_type, bool is_group4,
                                      uint32_t pos, uint32_t length,
                                      double cost,
                                      WP2::DecoderInfo* const decoder_info) {
  assert(length > 0);
  if (decoder_info == nullptr) return;
  if (!WP2::VDMatch(decoder_info->visual_debug, "bits-per-pixel/overall") &&
      !WP2::VDMatch(decoder_info->visual_debug, "lossless/symbols")) {
    return;
  }

  if (WP2::VDMatch(decoder_info->visual_debug, "lossless/symbols")) {
    // Only show once per tile.
    if (pos == 0 && VDSelected(0, 0, tile_->rect, config_)) {
      std::string* const str_ptr = &decoder_info->selection_info;
      WP2SAppend(str_ptr, "Tile at (%d, %d)\n", tile_->rect.x, tile_->rect.y);
      if (use_group4_) {
        WP2SAppend(str_ptr, "Group 4\n");
        WP2SAppend(str_ptr, "Move-to-front cache: %s\n",
                   mtf_.enabled() ? "enabled" : "disabled");
      } else {
        uint32_t cache_bits = hdr_.cache_config_.cache_bits;
        if (cache_bits > 0) {
          WP2SAppend(str_ptr, "Color cache: %d bits\n", cache_bits);
        } else {
          WP2SAppend(str_ptr, "No color cache\n");
        }
      }
    }

    Fill(pos, length, tile_->rect, kSymbolColors[is_group4][symbol_type],
         &decoder_info->debug_output);
  }

  const uint32_t tile_width = tile_->rect.width;
  if (VDSelectedLossless(pos, length, *decoder_info, tile_->rect)) {
    if (WP2::VDMatch(decoder_info->visual_debug, "lossless/symbols")) {
      Fill(pos, length, tile_->rect, {255, 255, 255, 255},
           &decoder_info->debug_output);
    }
    std::string* const str_ptr = &decoder_info->selection_info;
    WP2SAppend(str_ptr, "Pixel (%d, %d) Symbol at (%u, %u)\n",
               decoder_info->selection.x, decoder_info->selection.y,
               tile_->rect.x + pos % tile_width,
               tile_->rect.y + pos / tile_width);
    WP2SAppend(str_ptr, "Type %s length %d cost %f bpp %f\n",
               is_group4 ? kGroup4ModeNames[symbol_type]
                         : kSymbolTypeNames[symbol_type],
               length, cost, cost / length);
  }
}

static const char* const kTransformNames[] = {
    "Predictor", "CrossColor", "SubstractGreen", "ColorIndexing", "Group4"};
STATIC_ASSERT_ARRAY_SIZE(kTransformNames, NUM_TRANSFORMS);

void Decoder::RegisterUnprocessedRowForVDebug(
    uint32_t num_rows, uint32_t num_bits, const int16_t* const rows,
    WP2::DecoderInfo* const decoder_info) {
  if (decoder_info == nullptr) return;
  if (!WP2::VDMatch(decoder_info->visual_debug, "lossless/transformed")) {
    return;
  }

  // Only show once per tile.
  if (last_row_ == 0 && VDSelected(0, 0, tile_->rect, config_)) {
    std::string* const str_ptr = &decoder_info->selection_info;
    WP2SAppend(str_ptr, "Tile at (%d, %d): %d transform(s)\n", tile_->rect.x,
               tile_->rect.y, next_transform_);
    for (uint32_t i = 0; i < next_transform_; ++i) {
      const Transform& transform = transforms_[i];
      WP2SAppend(str_ptr, "Transform #%d: %s", i,
                 kTransformNames[transform.type_]);
      if (transform.type_ == CROSS_COLOR_TRANSFORM ||
          transform.type_ == PREDICTOR_TRANSFORM) {
        WP2SAppend(str_ptr, " Bits: %d", transform.bits_);
      } else if (transform.type_ == COLOR_INDEXING_TRANSFORM ||
                 transform.type_ == GROUP4) {
        WP2SAppend(str_ptr, " Num colors: %d", num_colors_);
      }
      WP2SAppend(str_ptr, "\n");
    }
  }

  const bool opaque = WP2::VDMatch(decoder_info->visual_debug, "force-opaque");
  for (uint32_t y = 0; y < num_rows; ++y) {
    const uint32_t abs_y = tile_->rect.y + last_row_ + y;
    uint8_t* const row = (uint8_t*)decoder_info->debug_output.GetRow(abs_y);
    for (uint32_t x = 0; x < tile_->rect.width; ++x) {
      const uint32_t abs_x = tile_->rect.x + x;
      const uint32_t offset = (y * tile_->rect.width + x) * 4;
      if (decoder_info->selection.Contains(abs_x, abs_y)) {
        std::string* const str_ptr = &decoder_info->selection_info;
        WP2SAppend(str_ptr, "Pixel (%d, %d) argb (%u, %u, %u, %u)\n",
                   decoder_info->selection.x, decoder_info->selection.y,
                   rows[offset + 0], rows[offset + 1], rows[offset + 2],
                   rows[offset + 3]);
      }
      for (uint32_t c = 0; c < 4; ++c) {
        row[abs_x * 4 + c] =
            WP2::RightShiftRound(rows[offset + c], num_bits - 8);
      }
      if (opaque) row[abs_x * 4 + 0] = 255;
    }
  }

  if (VDSelected(0, 0, tile_->rect, config_)) {
    decoder_info->debug_output.DrawRect(tile_->rect,
                                        WP2::Argb32b{255, 255, 0, 0});
  }
}

void Decoder::HeaderVDebug() {
  WP2::DecoderInfo* const decoder_info = config_.info;
  if (decoder_info == nullptr) return;

  if (!WP2::VDMatch(decoder_info->visual_debug, "lossless/clusters")) {
    return;
  }

  const uint32_t bits = hdr_.histogram_subsample_bits_;
  std::string* const str_ptr = &decoder_info->selection_info;
  const uint32_t num_clusters = hdr_.sr_.symbols_info().NumClusters(0);
  const WP2::Rectangle& tile_rect = tile_->rect;
  if (VDSelected(0, 0, tile_rect, config_)) {
    WP2SAppend(str_ptr, "Tile at (%d, %d): %d cluster(s) %d histo_bits\n",
               tile_rect.x, tile_rect.y, num_clusters, bits);
  }

  WP2::ArgbBuffer* const debug_output = &decoder_info->debug_output;
  if (num_clusters == 1) {
    debug_output->Fill(tile_rect, WP2::Argb32b{255, 0, 0, 0});
  } else {
    const uint32_t num_blocks_x = WP2LSubSampleSize(tile_rect.width, bits);
    const uint32_t num_blocks_y = WP2LSubSampleSize(tile_rect.height, bits);
    for (uint32_t y = 0; y < num_blocks_y; ++y) {
      for (uint32_t x = 0; x < num_blocks_x; ++x) {
        const WP2::Rectangle rect{
            tile_rect.x + (x << bits), tile_rect.y + (y << bits),
            std::min(1u << bits, tile_rect.width - (x << bits)),
            std::min(1u << bits, tile_rect.height - (y << bits))};
        const uint32_t cluster =
            hdr_.histogram_image_[4 * (num_blocks_x * y + x) + 2];
        const uint8_t gray = WP2::ToPixel8(cluster, num_clusters);
        debug_output->Fill(rect, WP2::Argb32b{255, gray, gray, gray});
        if (VDSelected(0, 0, rect, config_)) {
          WP2SAppend(str_ptr, "Rect (%d, %d) %d x %d : cluster %d\n", rect.x,
                     rect.y, rect.width, rect.height, cluster);
          debug_output->DrawRect(rect, WP2::Argb32b{255, 255, 0, 0});
        }
      }
    }
  }

  if (VDSelected(0, 0, tile_rect, config_)) {
    debug_output->DrawRect(tile_rect, WP2::Argb32b{255, 0, 255, 0});
  }
}

}  // namespace WP2L
