// Copyright 2018 Google LLC
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
//  Libwp2 SWIG interface definition.
//
// Author: Yannis Guyon (yguyon@google.com)

%module libwebp2

%{
#include "imageio/anim_image_dec.h"
#include "imageio/image_enc.h"
#include "imageio/image_dec.h"
#include "imageio/imageio_util.h"
#include "src/wp2/base.h"
#include "src/wp2/decode.h"
#include "src/wp2/encode.h"
%}

//------------------------------------------------------------------------------
// Typemaps

%include "stdint.i"
%include "typemaps.i"

// Maps uint8_t array and target language bytes structure.
#ifdef SWIGPYTHON
// Python bytes to C array.
%typemap(in) (const uint8_t* input_data, size_t input_data_size) {
  $1 = (uint8_t*)PyBytes_AsString($input);
  $2 = PyBytes_Size($input);
}

// C array to Python bytes.
%typemap(in, numinputs=0)
    (const uint8_t** output_data, size_t* output_data_size)
    (uint8_t* temp_output_data, size_t temp_output_data_size) {
  temp_output_data = nullptr;
  $1 = &temp_output_data;
  temp_output_data_size = 0;
  $2 = &temp_output_data_size;
}
%typemap(argout) (const uint8_t** output_data, size_t* output_data_size) {
  Py_XDECREF($result);
  $result = PyBytes_FromStringAndSize((const char*)*$1,
                                      (*$1 == nullptr) ? 0 : *$2);
}
#endif  // SWIGPYTHON

// Raises an exception when a function returns a bad WP2Status.
%typemap(out) WP2Status %{
  if ($1 != WP2_STATUS_OK) {
    SWIG_exception_fail(SWIG_RuntimeError,
        ("$name returned error " + std::to_string($1) + ": " +
         WP2GetStatusText($1)).c_str());
  }
  resultobj = SWIG_From_int(static_cast<int>(result));
%}

//------------------------------------------------------------------------------
// Base

extern int WP2GetVersion();
extern int WP2GetABIVersion();
extern bool WP2CheckVersion();

typedef enum {
  WP2_STATUS_OK = 0,
  WP2_STATUS_VERSION_MISMATCH,
  WP2_STATUS_OUT_OF_MEMORY,
  WP2_STATUS_INVALID_PARAMETER,
  WP2_STATUS_NULL_PARAMETER,
  WP2_STATUS_BAD_DIMENSION,
  WP2_STATUS_USER_ABORT,
  WP2_STATUS_UNSUPPORTED_FEATURE,
  WP2_STATUS_BITSTREAM_ERROR,
  WP2_STATUS_NOT_ENOUGH_DATA,
  WP2_STATUS_BAD_READ,
  WP2_STATUS_NEURAL_DECODE_FAILURE,
  WP2_STATUS_BITSTREAM_OUT_OF_MEMORY,
  WP2_STATUS_INVALID_CONFIGURATION,
  WP2_STATUS_BAD_WRITE,
  WP2_STATUS_FILE_TOO_BIG,
  WP2_STATUS_INVALID_COLORSPACE,
  WP2_STATUS_NEURAL_ENCODE_FAILURE
} WP2Status;

extern const char* WP2GetStatusText(WP2Status status);

typedef enum {
  WP2_Argb_32, WP2_ARGB_32, WP2_XRGB_32,
  WP2_rgbA_32, WP2_RGBA_32, WP2_RGBX_32,
  WP2_bgrA_32, WP2_BGRA_32, WP2_BGRX_32,
  WP2_RGB_24, WP2_BGR_24
} WP2SampleFormat;

namespace WP2 {

struct Rectangle {
  uint32_t x, y, width, height;
};

//------------------------------------------------------------------------------
// Reading

struct ArgbBuffer {
  const WP2SampleFormat format;
  uint32_t width;
  uint32_t height;
  uint32_t stride;

  ArgbBuffer(WP2SampleFormat format_in = WP2_Argb_32);
  bool IsEmpty() const;
  void* GetRow(uint32_t y);
  WP2Status CopyFrom(const ArgbBuffer& src);
  WP2Status ConvertFrom(const ArgbBuffer& src);
};

// Gets a pointer to data.
%extend ArgbBuffer {
void GetBytes(const uint8_t** output_data, size_t* const output_data_size) {
  if ((output_data != nullptr) && (output_data_size != nullptr)) {
    *output_data = (uint8_t*)$self->GetRow(0);
    *output_data_size = $self->height * $self->stride;
  }
}
};

WP2Status ReadImage(const char* file_path, ArgbBuffer* buffer,
                    size_t* file_size = nullptr);

%apply bool* OUTPUT { bool* const is_last };
%apply uint32_t* OUTPUT { uint32_t* const duration_ms };
class ImageReader {
 public:
  ImageReader(const char* const file_path, ArgbBuffer* const buffer);

  static constexpr uint32_t kInfiniteDuration = 0;
  WP2Status ReadFrame(bool* const is_last, uint32_t* const duration_ms);
};

//------------------------------------------------------------------------------
// Writing

WP2Status IoUtilWriteFile(const uint8_t* data, size_t data_size,
                          const char* const file_name, bool overwrite);
WP2Status SaveImage(const ArgbBuffer& buffer, const char* const file_path,
                    bool overwrite = false);

//------------------------------------------------------------------------------
// Encoding

class MemoryWriter {
 public:
  uint8_t* mem_;
  size_t size_;
};

// Gets a pointer to data.
%extend MemoryWriter {
void GetBytes(const uint8_t** output_data, size_t* const output_data_size) {
  if ((output_data != nullptr) && (output_data_size != nullptr)) {
    *output_data = $self->mem_;
    *output_data_size = $self->size_;
  }
}
};

struct EncoderConfig {
  float quality;
  int speed;
  int thread_level;

  bool IsValid() const;
};

WP2Status Encode(const ArgbBuffer& input, MemoryWriter* output,
                 const EncoderConfig& config);

class AnimationEncoder {
 public:
  WP2Status AddFrame(const ArgbBuffer& Argb, uint32_t duration_ms,
                     bool force_dispose = false);
  WP2Status Encode(MemoryWriter* memory_writer,
                   const EncoderConfig& config) const;
};

//------------------------------------------------------------------------------
// Decoding

struct BitstreamFeatures {
  uint32_t width;
  uint32_t height;
  int is_animation;
  int has_preview;
  int quality_hint;
  int has_icc;
  int has_trailing_data;
  uint8_t loop_count;
  size_t header_size;
};

struct FrameFeatures {
  uint32_t duration_ms;
  Rectangle window;
  bool is_last;
  uint32_t last_dispose_frame_index;
};

WP2Status Decode(const uint8_t* input_data, size_t input_data_size,
                 ArgbBuffer* output_buffer);

class Decoder {
 public:
  bool ReadFrame();
  WP2Status GetStatus() const;
  bool Failed() const;
  const BitstreamFeatures* TryGetDecodedFeatures() const;
  const FrameFeatures* TryGetFrameDecodedFeatures(uint32_t frame_index) const;
  uint32_t GetCurrentFrameIndex() const;
  uint32_t GetNumFrameDecodedFeatures() const;
  uint32_t GetNumAvailableFrames() const;
  Rectangle GetDecodedArea() const;
  const ArgbBuffer& GetPixels() const;
 protected:
  Decoder();
};

// Simple language-independant interface to access frame duration.
%extend Decoder {
  int GetFrameDurationMs() const {
    const WP2::FrameFeatures* current_frame_features =
        $self->TryGetFrameDecodedFeatures($self->GetCurrentFrameIndex());
    if (current_frame_features != nullptr) {
      return (int)current_frame_features->duration_ms;
    }
    return -1;
  }
};

class ArrayDecoder : public Decoder {
 public:
  ArrayDecoder();
  void SetInput(const uint8_t* input_data, size_t input_data_size);
};

class StreamDecoder : public Decoder {
 public:
  StreamDecoder();
  void AppendInput(const uint8_t* input_data, size_t input_data_size);
};

//------------------------------------------------------------------------------

}
