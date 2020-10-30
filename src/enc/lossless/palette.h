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
// Color Cache for WebP Lossless
//
// Author: Vincent Rabaud (vrabaud@google.com)

#ifndef WP2_ENC_LOSSLESS_PALETTE_H_
#define WP2_ENC_LOSSLESS_PALETTE_H_

#include "src/utils/ans.h"
#include "src/utils/vector.h"
#include "src/wp2/base.h"
#include "src/wp2/format_constants.h"

namespace WP2L {

// Forward declare.
class Encoder;
class InternalParams;

class Palette {
 public:
  Palette() : num_colors_(0) {}
  void Init(WP2SampleFormat format, bool has_alpha) {
    format_ = format;
    has_alpha_ = has_alpha;
  }
  // Populates the palette with colors from the given 'pic'.
  // If the unique color count is more than kMaxPaletteSize (or if there is only
  // one color), the resulting palette will be disabled (Size() returns 0).
  WP2Status AnalyzeAndCreate(const WP2::ArgbBuffer& pic);
  // Remaps argb values in 'pic' to packed palettes entries in dst[].
  // We assume that all values in 'pic' have a corresponding entry in the
  // palette.
  // Should be called after calling FindBestMethod() since the order of colors
  // can be changed on write.
  WP2Status Apply(const WP2::ArgbBuffer& pic, uint16_t* dst) const;
  // Returns the number of colors in the palette. Returns 0 if AnalyzeAndCreate
  // hasn't been called or if it was called on an image with more than
  // kMaxPaletteSize unique colors or less than two colors.
  size_t Size() const { return num_colors_; }
  // Try to find the best writing method 'method_'. Might reorder the colors_.
  // The storage cost cost_ will be calculated.
  // The 'speed' parameter impacts the depth of analyzis to get a better cost.
  WP2Status FindBestMethod(Encoder* const encoder, int speed);
  // Returns the bit cost of the current palette, set by FindBestMethod().
  float GetCost() const { return cost_; }
  // Writes the palette to the provided bitwriter 'enc' according to method_.
  WP2Status Write(WP2::ANSEncBase* const enc, WP2::ANSDictionaries* const dicts,
                  Encoder* const encoder) const;

  // for debugging
  const WP2::Vector_s16& GetColors() const { return colors_; }

 private:
  bool IsGrayScale() const;  // returns true if palette is grayscale (r==b==g)
  // try coding a palette 'colors' and if it has a better cost_, set it as the
  // new colors_.
  WP2Status TryPalette(WP2::Vector_s16& colors,
                       Encoder* const encoder, InternalParams& params,
                       bool also_as_image = true);

 private:
  WP2SampleFormat format_;
  bool has_alpha_ = false;
  bool is_grayscale_ = false;
  WP2::Vector_s16 colors_;
  uint32_t num_colors_ = 0;     // 0 = palette is disabled
  uint32_t max_palette_size_ = 0;   // depends on picture
  bool full_optimize_ = false;  // depending on size, do full size-optimization

  // The following is updated after calling AnalyzeAndCreate().
  enum PaletteStorageMethod {
    kMethodVerbatim = 0,   // store as-is
    kMethodImage,          // Store as an image.
    // Store as a list of values in range [0:previous value] or [previous
    // value:channel_max] (whichever is best according to the sign).
    kMethodList
  } method_ = kMethodVerbatim;
  float cost_;
};

}  // namespace WP2L

#endif  // WP2_ENC_LOSSLESS_PALETTE_H_
