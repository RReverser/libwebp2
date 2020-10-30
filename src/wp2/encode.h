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
//   WP2 encoder: main interface
//
// Author: Skal (pascal.massimino@gmail.com)

#ifndef WP2_WP2_ENCODE_H_
#define WP2_WP2_ENCODE_H_

#include <cassert>
#include <string>

#include "./base.h"
#include "./debug.h"
#include "./format_constants.h"

namespace WP2 {

//------------------------------------------------------------------------------
// Coding parameters

struct EncoderConfig {
 public:
  static const EncoderConfig kDefault;

  EncoderConfig() { assert(IsValid()); }

  // Returns true if all parameters are within their valid ranges.
  bool IsValid() const;

 public:
  // Compression parameters:
  float quality = 75.0f;  // Range: [0 = smallest file .. 100 = lossless]
                          // Quality in (95-100) range will use near-lossless.
                          // Quality 100 is strictly lossless.
  size_t target_size = 0;  // If non-zero, set the desired target size in bytes.
                           // Takes precedence over the 'quality' parameter.
  float target_psnr = 0.f;  // If non-zero, specifies the minimal distortion to
                            // try to achieve. Precedence over 'target_size'.

  float alpha_quality = 100.f;  // Range: [0 = smallest size .. 100 = lossless]
  int speed = 5;  // Quality/speed trade-off. Range: [0=fast .. 9=slower-better]

  // Side parameters:
  // Set whether the image will be rotated during decoding.
  Orientation decoding_orientation = Orientation::kOriginal;
  // Add a heavily compressed preview to be decoded and displayed before final
  // pixels (small size overhead up to kMaxPreviewSize).
  // Takes precedence over 'preview_size / preview[]' params.
  bool create_preview = false;
  // If preview_size>0, incorporate 'preview' data as preview in the bitstream.
  size_t preview_size = 0;
  const uint8_t* preview = nullptr;

  TransferFunction transfer_function = WP2_TF_ITU_R_BT2020_10BIT;

  // Parameters related to lossy compression only:
  int pass = 1;  // Number of entropy-analysis passes. Range: [1..10]
  // Spatial noise shaping strength in [0(=off), 100]
  // Affects how we spread noise between 'risky' areas (where noise is easily
  // visible) and easier areas (where it's less visible). A high SNS
  // value leads to skewing noise more towards areas where it should be less
  // visible. In general this improves SSIM but worsens PSNR.
  float sns = 50.f;
  int error_diffusion = 0;  // error diffusion strength [0=off, 100=max]

  int segments = 4;  // Max number of segments. Range: [1..kMaxNumSegments]
  // Selector for explicit or implicit segment-id.
  typedef enum {
    SEGMENT_ID_AUTO,     // use ID_EXPLICIT above a quality threshold
    SEGMENT_ID_EXPLICIT,
    SEGMENT_ID_IMPLICIT
  } SegmentIdMode;
  SegmentIdMode segment_id_mode = SEGMENT_ID_AUTO;
  // If non-zero, forces a quantization factor for the corresponding segment.
  float segment_factors[kMaxNumSegments] = { 0.f };

  // Size of tiles (width/height) in pixels. Each tile is compressed
  // independently, possibly in parallel.
  TileShape tile_shape = TILE_SHAPE_AUTO;

  // Algorithm for dividing the image into blocks.
  PartitionMethod partition_method = AUTO_PARTITIONING;
  // The set of allowed block sizes for image partitioning.
  PartitionSet partition_set = SOME_RECTS;
  // If true, use binary space partitioning instead of floating partition.
  bool partition_snapping = true;

  Csp csp_type = Csp::kYCoCg;  // Colorspace.

  typedef enum {
    UVModeAdapt = 0,  // Mix of 420 and 444 (per block)
    UVMode420,        // All blocks 420
    UVMode444,        // All blocks 444
    UVModeAuto,       // Choose any of the above automatically
    NumUVMode         // End-of-list marker
  } UVMode;
  UVMode uv_mode = UVModeAuto;  // Default sub-sampling mode for U/V planes.

  int preprocessing = 0;           // Preprocessing filter.
  int preprocessing_strength = 0;  // Range: [0 .. 100]

  bool use_random_matrix = false;  // Experimental.
  bool store_grain = false;        // Experimental: store grain info

  bool tune_perceptual = false;    // Experimental: tune for SSIM

  // Parameters related to lossless compression only:
  bool use_delta_palette = false;  // Reserved for future lossless feature.

  // 0 to disable multi-threading. Any other value represents the maximum number
  // of simultaneous extra threads (can be as large as 2^31 - 1).
  // No impact on encoded bytes.
  int thread_level = 0;

  // Memory usage reduction (but CPU use increase). No impact on encoded bytes.
  bool low_memory = false;

  // Neural compression:
  int use_neural_compression = 0;       // Neural network compression.
  const char* graphdef_path = nullptr;  // Directory holding encoder / decoder
                                        // graphdefs structure:
                                        // base/qq/[en|de]coder.pbbin

  // Enable alpha post processing filter (edge preserving blur).
  bool enable_alpha_filter = true;

  // TODO(maryla): remove from public API. Should match kDefaultQuantMultiplier.
  int u_quant_multiplier = 7;  // 0.875
  int v_quant_multiplier = 7;

  EncoderInfo* info = nullptr;   // If not null, report internal stats and info.
};

//------------------------------------------------------------------------------
// Picture hint

// Enumerate some predefined picture types for Encode(), depending on the
// source picture type.
typedef enum PictureHint {
  HINT_NONE = 0,     // nothing particular (default)
  HINT_PICTURE,      // digital picture, like portrait, inner shot
  HINT_PHOTO,        // outdoor photograph, with natural lighting
  HINT_DRAWING,      // hand or line drawing, with high-contrast details
  HINT_ICON,         // small-sized colorful images
  HINT_TEXT          // text-like
} PictureHint;

//------------------------------------------------------------------------------
// Abstract interface. Inherit to write the encoding output somewhere.

class Writer {
 public:
  virtual ~Writer() = default;

  // Will be called when 'data_size' bytes of encoded 'data' are ready.
  // Can be used to save the output of the encoding.
  // Should return false if an error happened.
  WP2_NO_DISCARD virtual bool Append(const void* data, size_t data_size) = 0;

  // TODO(yguyon):
  // Will be regularly called during the encoding process. Can be used to show
  // the progress to the user and add the possibility to stop the encoding.
  // 'percent' is in range [0..100].
  // Should return false if an error happened or if the encoding has to stop.
  WP2_NO_DISCARD virtual bool UpdateProgress(int percent) {
    (void)percent;
    return true;
  }
};

// Writes the encoding output to memory.
class MemoryWriter : public Writer {
 public:
  ~MemoryWriter() override { Reset(); }
  void Reset();  // Releases memory.
  WP2_NO_DISCARD bool Append(const void* data, size_t data_size) override;

  uint8_t* mem_ = nullptr;  // Final buffer (of capacity 'max_size' >= 'size').
  size_t size_ = 0;         // Final size.
  size_t max_size_ = 0;     // Total capacity.
};

// Writes the encoding output to a string.
class StringWriter : public Writer {
 public:
  explicit StringWriter(std::string* const str) : str_(str) {}
  WP2_NO_DISCARD bool Append(const void* data, size_t data_size) override;

 private:
  std::string* const str_;
};

// Counts the bytes in the encoding output.
class Counter : public Writer {
 public:
  WP2_NO_DISCARD bool Append(const void* data, size_t data_size) override;

  size_t total_size_ = 0;
};

//------------------------------------------------------------------------------
// Main encoding function

// Encodes an image. Returns WP2_STATUS_OK if no error occurred.
WP2_NO_DISCARD WP2Status
Encode(const ArgbBuffer& input, Writer* output,
       const EncoderConfig& config = EncoderConfig::kDefault,
       PictureHint picture_hint = HINT_NONE);

//------------------------------------------------------------------------------
// Specialized call

// Encodes an image of 'width' by 'height' pixels of channels 'c0, c1, c2' in a
// custom color space (converted from alpha-premultiplied RGB samples), 10b max.
// The '_step' parameters specify the distances in samples between scanlines and
// must be multiples of two. The 'ccsp_to_rgb_matrix' and the following
// 'ccsp_to_rgb_shift' are assumed to output RGB with 8 unsigned bits of
// precision. Alpha must have a precision of 8 unsigned bits unless 'a_buffer'
// is null. Returns WP2_STATUS_OK or an error.
WP2_NO_DISCARD WP2Status Encode(
    uint32_t width, uint32_t height,
    const int16_t* c0_buffer, uint32_t c0_step,
    const int16_t* c1_buffer, uint32_t c1_step,
    const int16_t* c2_buffer, uint32_t c2_step,
    const int16_t* a_buffer, uint32_t a_step,
    const int16_t ccsp_to_rgb_matrix[9], uint32_t ccsp_to_rgb_shift,
    Writer* output, const EncoderConfig& config = EncoderConfig::kDefault,
    const Metadata& metadata = Metadata());

//------------------------------------------------------------------------------
// AnimationEncoder

class AnimationEncoder {
 public:
  AnimationEncoder();
  AnimationEncoder(const AnimationEncoder&) = delete;
  AnimationEncoder(AnimationEncoder&&) noexcept;
  ~AnimationEncoder();

  // Copies data from 'pixels' to an internal buffer. 'pixels.metadata' is
  // ignored. Frames must have a 'duration_ms' of at least 1 millisecond.
  // 'force_dispose' may increase the quality of the current frame at the cost
  // of a bigger file. Returns WP2_STATUS_OK if no error occurred.
  WP2_NO_DISCARD WP2Status AddFrame(const ArgbBuffer& pixels,
                                    uint32_t duration_ms,
                                    bool force_dispose = false,
                                    PictureHint picture_hint = HINT_NONE);

  // Same as above but for a custom color space. See Encode().
  // All input frames must have the same color space within an animation.
  WP2_NO_DISCARD WP2Status AddFrame(
      uint32_t width, uint32_t height,
      const int16_t* c0_buffer, uint32_t c0_step,
      const int16_t* c1_buffer, uint32_t c1_step,
      const int16_t* c2_buffer, uint32_t c2_step,
      const int16_t* a_buffer, uint32_t a_step,
      const int16_t ccsp_to_rgb_matrix[9], uint32_t ccsp_to_rgb_shift,
      uint32_t duration_ms, bool force_dispose = false);

  // Encodes all frames. Call it again to compress the same frames with a
  // different 'config' and/or to another 'output'.
  // Returns WP2_STATUS_OK if no error occurred.
  WP2_NO_DISCARD WP2Status
  Encode(Writer* output, const EncoderConfig& config = EncoderConfig::kDefault,
         uint8_t loop_count = kInfiniteLoop,
         const Metadata& metadata = Metadata()) const;

  struct State;

 private:
  State* state_;  // Pointer to implementation.
};

}  // namespace WP2

#endif /* WP2_WP2_ENCODE_H_ */
