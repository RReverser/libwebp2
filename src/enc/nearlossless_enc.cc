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
//  Near-lossless pre-processing
//
// Author: Skal (pascal.massimino@gmail.com)
//

#include "src/dsp/dsp.h"
#include "src/enc/analysis.h"
#include "src/utils/vector.h"

namespace WP2 {

// returns the number of bits to remove (q=99 -> 1bit, q=96 -> 4bit)
static uint32_t NumBitsToRemove(float quality) {
  return (uint32_t)std::min(4.f, std::max(1.f, 100.f - quality));
}

static void BuildMap(uint8_t map[256], uint32_t bits) {
  const uint32_t mask = (1u << bits) - 1;
  for (uint32_t i = 0; i < 256; ++i) {
    const uint32_t biased = i + (mask >> 1) + ((i >> bits) & 1);
    map[i] = (biased > 0xff) ? 0xff : (biased & ~mask);
  }
}

static uint32_t SmoothValue(uint32_t v, const uint8_t map[]) {
  return (map[(v >> 24) & 0xff] << 24) | (map[(v >> 16) & 0xff] << 16) |
         (map[(v >>  8) & 0xff] <<  8) | (map[(v >>  0) & 0xff] <<  0);
}

static bool IsFar(uint32_t a, uint32_t b, int limit) {
  for (uint32_t k = 0; k < 32; k += 8) {
    const int delta =
        (int)((a >> k) & 0xff) - (int)((b >> k) & 0xff);
    if (delta >= limit || delta <= -limit) {
      return true;
    }
  }
  return false;
}

static uint32_t Smooth(const uint32_t* const prev,
                       const uint32_t* const cur,
                       const uint32_t* const next,
                       uint32_t x, uint32_t limit,
                       const uint8_t map[]) {
  uint32_t v = cur[x];
  if (IsFar(v, cur[x - 1], limit) || IsFar(v, cur[x + 1], limit) ||
      IsFar(v, prev[x], limit) || IsFar(v, next[x], limit)) {
    v = SmoothValue(v, map);
  }
  return v;
}

WP2Status PreprocessNearLossless(const ArgbBuffer& in_buffer,
                                 const EncoderConfig& config, bool is_alpha,
                                 ArgbBuffer* const out_buffer) {
  assert(out_buffer != nullptr);
  assert(in_buffer.format == WP2_Argb_32);
  assert(out_buffer->format == WP2_Argb_32);
  const uint32_t w = in_buffer.width, h = in_buffer.height;
  WP2_CHECK_STATUS(out_buffer->Resize(w, h));
  Vector_u32 tmp;
  WP2_CHECK_ALLOC_OK(tmp.resize(2 * w));   // temporary rotating buffer

  const ArgbBuffer* src = &in_buffer;  // source
  const float quality = is_alpha ? config.alpha_quality : config.quality;
  for (uint32_t bits = NumBitsToRemove(quality); bits > 0; --bits) {
    const uint32_t limit = std::max(1u, ((1u << bits) - 1u) / 2u);
    uint8_t map[256];
    BuildMap(map, bits);  // could be cached per 'bits' value [1..4]
    uint32_t* prev = &tmp[0 * w];
    uint32_t* cur  = &tmp[1 * w];
    const uint32_t* in = (const uint32_t*)src->GetRow(0);
    const uint32_t* in_p = in;
    for (uint32_t y = 0; y < h; ++y) {
      const uint32_t* const in_n =
          (const uint32_t*)src->GetRow(y + 1 < h ? y + 1 : h - 1);
      cur[0] = in[0];
      cur[w - 1] = in[w - 1];
      for (uint32_t x = 1; x + 1 < w; x += 1) {
        cur[x] = Smooth(in_p, in, in_n, x, limit, map);
      }
      if (y > 0) memcpy(out_buffer->GetRow(y - 1), prev, w * sizeof(*prev));
      std::swap(cur, prev);
      in_p = in;
      in = in_n;
    }
    memcpy(out_buffer->GetRow(h - 1), prev, w * sizeof(*prev));
    src = out_buffer;  // use output buffer as source for next passes
  }
  return WP2_STATUS_OK;
}

}  // namespace WP2
