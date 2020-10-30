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
//  Common functions
//
// Author: Skal (pascal.massimino@gmail.com)

#ifndef WP2_WP2_BASE_H_
#define WP2_WP2_BASE_H_

#define WP2_VERSION     0x000001    // MAJOR(8b) + MINOR(8b) + REVISION(8b)
#define WP2_ABI_VERSION   0x0001    // MAJOR(8b) + MINOR(8b)

//------------------------------------------------------------------------------

#include <cassert>
#include <cinttypes>
#include <cstddef>  // for size_t

#ifndef WP2_EXTERN
// This explicitly marks library functions and allows for changing the
// signature for e.g., Windows DLL builds.
# if defined(__GNUC__)
#  define WP2_EXTERN extern __attribute__ ((visibility ("default")))
# else
#  define WP2_EXTERN extern
# endif  /* __GNUC__ >= 4 */
#endif  /* WP2_EXTERN */

#if defined(__GNUC__)
// With C++17, we could use the [[nodiscard]] attribute.
#define WP2_NO_DISCARD __attribute__((warn_unused_result))
#else
#define WP2_NO_DISCARD
#endif

// Macro to check ABI compatibility (same major revision number)
#define WP2_ABI_IS_INCOMPATIBLE(a, b) (((a) >> 8) != ((b) >> 8))

// Return the library's version number, packed in hexadecimal using 8bits for
// each of major/minor/revision. E.g: v2.5.7 is 0x020507.
WP2_EXTERN int WP2GetVersion();

// Return the library's ABI version number, similar to WP2GetVersion()
WP2_EXTERN int WP2GetABIVersion();

// Return true if the link and compile version are compatible.
// Should be checked first before using the library.
static inline bool WP2CheckVersion() {
  return !WP2_ABI_IS_INCOMPATIBLE(WP2GetABIVersion(), WP2_ABI_VERSION);
}

//------------------------------------------------------------------------------
// Status enumeration

typedef enum
#if __cplusplus >= 201703L
    [[nodiscard]]  // Triggers a warning on calls not using the returned value
                   // of functions returning a WP2Status by value.
#endif
{ WP2_STATUS_OK = 0,
  WP2_STATUS_VERSION_MISMATCH,
  WP2_STATUS_OUT_OF_MEMORY,            // memory error allocating objects
  WP2_STATUS_INVALID_PARAMETER,        // a parameter value is invalid
  WP2_STATUS_NULL_PARAMETER,           // a pointer parameter is NULL
  WP2_STATUS_BAD_DIMENSION,            // picture has invalid width/height
  WP2_STATUS_USER_ABORT,               // abort request by user
  WP2_STATUS_UNSUPPORTED_FEATURE,      // unsupported feature
  // the following are specific to decoding:
  WP2_STATUS_BITSTREAM_ERROR,          // bitstream has syntactic error
  WP2_STATUS_NOT_ENOUGH_DATA,          // premature EOF during decoding
  WP2_STATUS_BAD_READ,                 // error while reading bytes
  WP2_STATUS_NEURAL_DECODE_FAILURE,    // neural decoder failed
  // the following are specific to encoding:
  WP2_STATUS_BITSTREAM_OUT_OF_MEMORY,  // memory error while flushing bits
  WP2_STATUS_INVALID_CONFIGURATION,    // encoding configuration is invalid
  WP2_STATUS_BAD_WRITE,                // error while flushing bytes
  WP2_STATUS_FILE_TOO_BIG,             // file is bigger than 4G
  WP2_STATUS_INVALID_COLORSPACE,       // encoder called with bad colorspace
  WP2_STATUS_NEURAL_ENCODE_FAILURE,    // neural encoder failed

  WP2_STATUS_LAST                      // list terminator. always last.
} WP2Status;

// returns a printable string for the name of the status
WP2_EXTERN const char* WP2GetStatusMessage(WP2Status status);
// returns a short description of the status
WP2_EXTERN const char* WP2GetStatusText(WP2Status status);

//------------------------------------------------------------------------------
//  WP2SampleFormat

// The format of pixels passed for import into ArgbBuffer.
// These enum names describe the samples order are they are laid out
// in memory. For instance WP2_ARGB_32 will refer to 8b samples
// ordered as: A,R,G,B,A,R,G,B, ...
// The 'X' means that the alpha values are to be ignored and assumed to be set
// to 0xff everywhere (the image has no transparency). This assumption is not
// checked nor enforced.
// Non-capital names (e.g.:MODE_Argb) relates to pre-multiplied RGB channels.
typedef enum {
  // 32b/pixel formats
  WP2_Argb_32 = 0,
  WP2_ARGB_32,
  WP2_XRGB_32,

  WP2_rgbA_32,
  WP2_RGBA_32,
  WP2_RGBX_32,

  WP2_bgrA_32,
  WP2_BGRA_32,
  WP2_BGRX_32,

  // 24b/pixel formats
  WP2_RGB_24,
  WP2_BGR_24,

  // HDR format: 8 bits for A, 10 per RGB. Internally each is stored as uint16.
  WP2_Argb_38,

  // TODO(skal): RGB565_16? rgbA4444_16? RGBA4444_16?
  WP2_FORMAT_NUM
} WP2SampleFormat;

// Returns the number of bytes per sample for the given format.
WP2_EXTERN uint32_t WP2FormatBpp(WP2SampleFormat format);

// Returns the number of bytes used per channel internally.
WP2_EXTERN uint32_t WP2FormatBpc(WP2SampleFormat format);

// Returns the index of the alpha channel, or the number of channels if none.
WP2_EXTERN uint32_t WP2AlphaChannelIndex(WP2SampleFormat format);

namespace WP2 {

//------------------------------------------------------------------------------
// distortion type

typedef enum {
  PSNR,
  PSNRHVS,
  PSNR_YUV,  // PSNR_YUV is PSNR in YCbCr space for now (could be any other YUV)
  SSIM,
  SSIM_YUV,
  MSSSIM,
  LSIM,
  NUM_METRIC_TYPES
} MetricType;

//------------------------------------------------------------------------------
// Color structs

// The red, green and blue channels are premultiplied by the alpha component.
struct Argb32b {
  uint8_t a;        // Precision:  8 bits
  uint8_t r, g, b;  // Precision:  8 bits
};

struct Argb38b {
  uint8_t a;         // Precision:  8 bits
  uint16_t r, g, b;  // Precision:  10 bits
};

struct RGB12b {
  uint8_t r, g, b;  // Precision:  4 bits
};

struct Ayuv38b {
  uint8_t a;        // Precision:  8 bits
  int16_t y, u, v;  // Precision:  10 bits
};

// Colorspaces.
enum class Csp { kYCoCg, kYCbCr, kCustom, kYIQ };
constexpr uint32_t kNumCspTypes = 4;

//------------------------------------------------------------------------------
// Data and Metadata container

// Struct owning a chunk of memory which will be freed upon destruction.
struct Data {
  Data() = default;
  Data(Data&& other) noexcept;
  Data(const Data&) = delete;
  ~Data() { Clear(); }

  bool IsEmpty() const { return (size == 0); }

  // Deallocates the memory, and re-initializes the structure.
  void Clear();
  // Increases or decreases allocated memory to 'data_size' bytes.
  // If 'keep_bytes' is true, the 'bytes' will contain the same data as before.
  WP2_NO_DISCARD WP2Status Resize(size_t data_size, bool keep_bytes);
  // Resizes to 'data_size' then copies all 'data'.
  WP2_NO_DISCARD WP2Status CopyFrom(const uint8_t* data, size_t data_size);
  // Resizes to 'size + data_size' then copies all 'data'.
  WP2_NO_DISCARD WP2Status Append(const uint8_t* data, size_t data_size);

  uint8_t* bytes = nullptr;
  size_t size = 0;
};

void swap(Data& a, Data& b);

//------------------------------------------------------------------------------

struct Metadata {
  Metadata() = default;
  Metadata(const Metadata&) = delete;
  Metadata(Metadata&&) = default;

  Data exif;
  Data iccp;
  Data xmp;

  // Returns true is there is no metadata stored.
  bool IsEmpty() const;
  // Deallocate the memory, and re-initialize the structure.
  void Clear();
  // Copy exif, iccp and xmp. Returns WP2_STATUS_OK or an error.
  WP2_NO_DISCARD WP2Status CopyFrom(const Metadata& src);
};

void swap(Metadata& a, Metadata& b);

//------------------------------------------------------------------------------
// Rectangle struct. Rectangle(0, 0, 1, 1) represents the top-left most pixel.

struct Rectangle {
  uint32_t x = 0;
  uint32_t y = 0;
  uint32_t width = 0;
  uint32_t height = 0;

  Rectangle() = default;
  Rectangle(uint32_t rect_x, uint32_t rect_y,
            uint32_t rect_width, uint32_t rect_height)
      : x(rect_x), y(rect_y), width(rect_width), height(rect_height) {}
  bool operator==(const Rectangle& other) const {
    return x == other.x && y == other.y &&
           width == other.width && height == other.height;
  }

  inline uint64_t GetArea() const { return (uint64_t)width * height; }
  inline bool Contains(uint32_t px, uint32_t py) const {
    return (px >= x && px < x + width) && (py >= y && py < y + height);
  }
  inline bool Contains(const WP2::Rectangle& r) const {
    if (r.GetArea() == 0) return Contains(r.x, r.y);
    return Contains(r.x, r.y) &&
           Contains(r.x + r.width - 1, r.y + r.height - 1);
  }
  inline bool Intersects(const WP2::Rectangle& r) const {
    return ((x >= r.x && x < r.x + r.width) || (r.x >= x && r.x < x + width)) &&
           ((y >= r.y && y < r.y + r.height) || (r.y >= y && r.y < y + height));
  }
  Rectangle ClipWith(const Rectangle& r) const;
  Rectangle MergeWith(const Rectangle& r) const;
};

enum class Orientation { kOriginal = 0, k90Clockwise, k180, k270Clockwise };

//------------------------------------------------------------------------------
// Main Input/Output buffer structure

struct ArgbBuffer {
 public:
  explicit ArgbBuffer(WP2SampleFormat format_in = WP2_Argb_32) noexcept;
  ArgbBuffer(const ArgbBuffer&) = delete;  // No invisible copy, for safety.
  ArgbBuffer(ArgbBuffer&&) noexcept;
  ~ArgbBuffer() { Deallocate(); }

  // Returns true if this buffer is empty (has no pixel).
  bool IsEmpty() const { return (pixels == nullptr); }
  // Deallocates this buffer's memory (including metadata) and resets it empty.
  void Deallocate();

  // Returns an opaque pointer to the first channel element of row 'y'.
  const void* GetRow(uint32_t y) const;
  void* GetRow(uint32_t y);
  const uint8_t* GetRow8(uint32_t y) const;
  uint8_t* GetRow8(uint32_t y);
  const uint16_t* GetRow16(uint32_t y) const;
  uint16_t* GetRow16(uint32_t y);

  // Sets the format but only if the buffer is empty.
  WP2Status SetFormat(WP2SampleFormat format_in);

  // Returns true if this buffer points to a non-owned memory (View).
  bool IsView() const { return is_external_memory; }
  // Returns true if this buffer points to slow memory.
  bool IsSlow() const { return is_slow_memory; }

  // Resizes this buffer to dimensions 'new_width' x 'new_height'.
  // If it's a view, success depends on whether it's already exactly well-sized.
  // Otherwise it will try to allocate the exact size, if not already the case.
  // The content of this buffer either stays the same or is uninitialized.
  // Previous metadata is discarded and deallocated.
  // Returns WP2_STATUS_OK upon success (no memory error or invalid arguments).
  WP2_NO_DISCARD WP2Status Resize(uint32_t new_width, uint32_t new_height);

  // Calls Resize() with dimensions of 'src' and copies the content of 'src'
  // into this buffer. Only WP2_Argb_32 format is accepted. Previous metadata is
  // discarded and deallocated, new metadata is not copied.
  // Returns WP2_STATUS_OK upon success.
  WP2_NO_DISCARD WP2Status CopyFrom(const ArgbBuffer& src);

  // Same as CopyFrom() but allows this buffer to be in any format. 'src' should
  // be in WP2_Argb_32 format. Previous metadata is discarded and deallocated,
  // new metadata is not copied. Returns WP2_STATUS_OK upon success.
  WP2_NO_DISCARD WP2Status ConvertFrom(const ArgbBuffer& src);

  // Calls Resize() and copies 'samples' in 'input_format' into this buffer,
  // converting if necessary (only towards WP2_Argb_32).
  // Returns WP2_STATUS_OK upon success.
  WP2_NO_DISCARD WP2Status Import(WP2SampleFormat input_format,
                                  uint32_t new_width, uint32_t new_height,
                                  const uint8_t* samples,
                                  uint32_t samples_stride);

  // Copies a number of pixels equal to this buffer's width from 'samples' in
  // 'input_format' into this buffer at 'row_index', converting if necessary
  // (only towards WP2_Argb_32). Buffer should be already allocated.
  // Returns WP2_STATUS_OK upon success.
  WP2_NO_DISCARD WP2Status ImportRow(WP2SampleFormat input_format,
                                     uint32_t row_index,
                                     const uint8_t* samples);

  // Fills the area or just a rectangle border defined by 'window' with 'color'.
  // The rectangle will be clipped to the dimension of the buffer.
  void Fill(const Rectangle& window, Argb32b color);
  void Fill(const Rectangle& window, Argb38b color);
  void Fill(Argb32b color);  // fills the whole buffer
  void DrawRect(const Rectangle& window, Argb32b color);
  // 'color' must contain WP2FormatBpp(format) bytes.
  void Fill(const Rectangle& window, const uint8_t color[]);
  void Fill(const Rectangle& window, const uint16_t color[]);

  // Composites this buffer on the given background color.
  WP2_NO_DISCARD WP2Status CompositeOver(Argb32b background);
  // Composites this buffer on the given background buffer. Buffer size must
  // match.
  WP2_NO_DISCARD WP2Status CompositeOver(const ArgbBuffer& background);
  // Composites the given buffer on top of this buffer. Buffer size must match.
  // If 'window' is not empty, only this area is composited.
  WP2_NO_DISCARD WP2Status CompositeUnder(const ArgbBuffer& foreground,
                                          const Rectangle& window = {});

  // Deallocates then initiliazes this buffer as a view pointing to 'src'.
  // Previous metadata is discarded and deallocated, new metadata is not copied
  // nor pointed to. The view may modify the pixels of 'src'.
  // Returns WP2_STATUS_OK upon success.
  WP2_NO_DISCARD WP2Status SetView(const ArgbBuffer& src);
  // Same as above but the view is a 'window' of 'src' buffer.
  WP2_NO_DISCARD WP2Status SetView(const ArgbBuffer& src,
                                   const Rectangle& window);

  // Deallocates then initiliazes this buffer as a view pointing to the external
  // 'samples'. If this buffer 'is_slow', fewer read/write operations will
  // happen during encoding/decoding, at the expense of processing time and/or
  // internal memory usage. Returns WP2_STATUS_OK upon success.
  WP2_NO_DISCARD WP2Status SetExternal(uint32_t new_width, uint32_t new_height,
                                       uint8_t* samples, uint32_t new_stride,
                                       bool is_slow = false);
  WP2_NO_DISCARD WP2Status SetExternal(uint32_t new_width, uint32_t new_height,
                                       uint16_t* samples, uint32_t new_stride,
                                       bool is_slow = false);

  // Scans the buffer for the presence of non fully opaque alpha values.
  // Returns true in such case. Otherwise returns false (indicating that the
  // alpha plane can be ignored altogether e.g.).
  bool HasTransparency() const;

  // Computes PSNR, SSIM, LSIM or PSNRHVS distortion metric between two images.
  // Results are in dB, stored in result[] in the A/R/G/B/All order.
  // For PSNRHVS, the result[] content is A/Y/U/V/All.
  // Warning: this function is rather CPU-intensive.
  // Returns WP2_STATUS_OK upon success (no memory allocation error, etc.)
  WP2_NO_DISCARD WP2Status GetDistortion(const ArgbBuffer& ref,
                                         MetricType metric_type,
                                         float result[5]) const;
  // Same as above but computes distortion metric on a sub window.
  WP2_NO_DISCARD WP2Status GetDistortion(const ArgbBuffer& ref,
                                         const Rectangle& window,
                                         MetricType metric_type,
                                         float result[5]) const;
  // Same as above, but computes the distortion at a given x,y point (slow!).
  // metric_type should not be PSNRHVS.
  WP2_NO_DISCARD WP2Status GetDistortion(const ArgbBuffer& ref, uint32_t x,
                                         uint32_t y, MetricType metric_type,
                                         float* result) const;

  // Computes the distortion when the image is composited on a white or black
  // background and returns the worst of the two. For images with no
  // transparency, distortion is only computed once since the background doesn't
  // matter.
  WP2_NO_DISCARD WP2Status GetDistortionBlackOrWhiteBackground(
      const ArgbBuffer& ref, MetricType metric_type, float result[5]) const;
  // Same as above but computes distortion metric on a sub window.
  WP2_NO_DISCARD WP2Status GetDistortionBlackOrWhiteBackground(
      const ArgbBuffer& ref, const Rectangle& window, MetricType metric_type,
      float result[5]) const;

  // Returns WP2_STATUS_OK upon success (same format expected for both).
  WP2_NO_DISCARD WP2Status Swap(ArgbBuffer* other);

  // Simple in-place 2x down-sampler. Cheap! (box-average)
  void SimpleHalfDownsample();

 public:
  // If 'include_alpha_in_all', result[4] will take alpha into account (if any).
  // Otherwise alpha can still have an influence since we use premultiplied
  // RGB, except if the image is all black.
  WP2_NO_DISCARD WP2Status GetDistortion(const ArgbBuffer& ref,
                                         MetricType metric_type,
                                         bool include_alpha_in_all,
                                         float result[5]) const;

  // Should always be WP2_Argb_32 for encoding.
  // TODO(vrabaud) make it private and use a getter.
  const WP2SampleFormat format;

  uint32_t width = 0;   // Number of pixels on the horizontal axis.
  uint32_t height = 0;  // Number of pixels on the vertical axis.
  uint32_t stride = 0;  // Difference in bytes from a scanline to the next one.

  Metadata metadata;    // ICC / XMP / EXIF metadata.

 private:
  void FillImpl(const Rectangle& window, const void* color);
  WP2Status SetExternalImpl(uint32_t new_width, uint32_t new_height,
                            void* samples, uint32_t new_stride, bool is_slow);

  void* pixels = nullptr;           // Pointer to samples.
  bool is_external_memory = false;  // If false, 'pixels' is 'private_memory'.
  void* private_memory = nullptr;   // Internally allocated memory.
  bool is_slow_memory = false;  // If true, the memory is considered 'slow' and
                                // the number of reads/writes will be reduced.
};

//------------------------------------------------------------------------------

// Can be inherited to get notified of progress and/or to abort encoding or
// decoding.
class ProgressHook {
 public:
  virtual ~ProgressHook() = default;
  // Will be regularly called with 'progress' advancing from 0.f (decoding
  // didn't start) to 1.f (decoding is done). For animations, 'progress' will
  // be reset to 0.f at each frame beginning; otherwise it never decreases.
  // Return true to continue, false to stop decoding.
  virtual bool OnUpdate(float progress) = 0;
};

//------------------------------------------------------------------------------

}   // namespace WP2

#endif  /* WP2_WP2_BASE_H_ */
