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

// Common includes for tests.
// Helper for accessing test data files.

#ifndef WP2_TESTS_INCLUDE_HELPERS_H_
#define WP2_TESTS_INCLUDE_HELPERS_H_

#include <cmath>
#include <iostream>
#include <random>
#include <string>
#include <vector>

#include "src/common/lossy/block_size.h"
#include "src/utils/random.h"
#include "src/utils/vector.h"
#include "src/utils/plane.h"
#include "src/dsp/dsp.h"
#include "src/wp2/decode.h"
#include "src/wp2/encode.h"

#include "gmock/gmock.h"
#include "gtest/gtest.h"

#if !defined(__has_feature)
#define __has_feature(x) 0
#endif

#if __has_feature(memory_sanitizer)
#define WP2_HAVE_MSAN 1
#include <sanitizer/msan_interface.h>
#endif

// Helper to assert that a function returning a WP2Status returned OK.
#define ASSERT_WP2_OK(expression) ASSERT_EQ(expression, WP2_STATUS_OK)
#define EXPECT_WP2_OK(expression) EXPECT_EQ(expression, WP2_STATUS_OK)

// Undocumented, just for testing.
// WARNING WARNING! This function is NOT thread safe, and should be used
// in a single-thread testing environment!
extern void WP2DspReset();

namespace WP2 {
namespace testing {

//------------------------------------------------------------------------------

// the following is only initialized once:
static const WP2CPUInfo GetCPUInfo = WP2GetCPUInfo;

static bool GetCPUInfoNoSSE(WP2CPUFeature feature) {
  if (feature == kSSE4_1 || feature == kSSE4_2 || feature == kSSE ||
      feature == kAVX) {
    return false;
  }
  return GetCPUInfo(feature);
}

static bool GetCPUInfoNoAVX(WP2CPUFeature feature) {
  if (feature == kAVX) return false;
  return GetCPUInfo(feature);
}

static bool GetCPUInfoForceSlowSSSE3(WP2CPUFeature feature) {
  if (feature == kSlowSSSE3 && GetCPUInfo(kSSE3)) {
    return true;  // we have SSE3 -> force SlowSSSE3
  }
  return GetCPUInfo(feature);
}

static bool GetCPUInfoOnlyC(WP2CPUFeature feature) { return false; }

// To be used with '::testing::ValuesIn(kWP2CPUInfos)'
// Warning! No multi-thread safe! don't use in *sharded* tests!
const WP2CPUInfo kWP2CpuInfos[] = {
  GetCPUInfoOnlyC, GetCPUInfoForceSlowSSSE3,
  GetCPUInfoNoSSE, GetCPUInfoNoAVX,
  GetCPUInfo
};

// Same, but as struct rather than typedef
static struct WP2CPUInfoStruct {
  const WP2CPUInfo cpu_info;
  const char* const name;
} const kWP2CpuInfoStructs[] = {
  { GetCPUInfoOnlyC, "OnlyC" },
  { GetCPUInfoForceSlowSSSE3, "ForceSlowSSSE3" },
  { GetCPUInfoNoSSE, "NoSSE" },
  { GetCPUInfoNoAVX, "NoAVX" },
  { GetCPUInfo, "Regular GetCPUInfo" },
};

//------------------------------------------------------------------------------

std::string GetTestDataPath(const std::string& file_name);
void GetTestDataPaths(const std::vector<const char*>& file_names,
                      std::vector<std::string>* const file_paths);

std::string GetTempDataPath(const std::string& file_path);

//------------------------------------------------------------------------------

bool HasSameData(const Data& a, const Data& b);

//------------------------------------------------------------------------------

WP2Status ReadImages(const std::vector<std::string>& file_paths,
                     std::vector<ArgbBuffer>* const frames);

WP2Status ReadAnimation(const uint8_t* const data, size_t size,
                        std::vector<ArgbBuffer>* const frames,
                        std::vector<uint32_t>* const frame_durations = nullptr,
                        uint32_t* const loop_count = nullptr);

WP2Status ReadAnimation(const std::string& file_path,
                        std::vector<ArgbBuffer>* const frames,
                        std::vector<uint32_t>* const frame_durations = nullptr,
                        uint32_t* const loop_count = nullptr);

//------------------------------------------------------------------------------

WP2Status CompressImage(
    const std::string& file_name, Writer* const encoded_data,
    ArgbBuffer* original_image = nullptr,
    float quality = EncoderConfig::kDefault.quality,
    int speed = EncoderConfig::kDefault.speed,
    int thread_level = EncoderConfig::kDefault.thread_level,
    uint32_t num_downsamplings = 0);

WP2Status CompressAnimation(
    const std::vector<const char*>& frames_file_names,
    const std::vector<uint32_t>& durations_ms, Writer* const encoded_data,
    std::vector<ArgbBuffer>* frames = nullptr,
    float quality = EncoderConfig::kDefault.quality,
    int speed = EncoderConfig::kDefault.speed,
    int thread_level = EncoderConfig::kDefault.thread_level,
    uint32_t num_downsamplings = 0);

//------------------------------------------------------------------------------

// Returns a lower-bound expected PSNR.
float GetExpectedDistortion(
    float quality = EncoderConfig::kDefault.quality,
    float alpha_quality = EncoderConfig::kDefault.alpha_quality);

float GetExpectedDistortion(const EncoderConfig& encoder_config);

// Returns true if 'src' and 'dec' are similar of at least 'expected_distortion'
::testing::AssertionResult Compare(
    const ArgbBuffer& src, const ArgbBuffer& dec, const std::string& file_name,
    float expected_distortion = 99.f, MetricType metric = PSNR,
    Orientation decoding_orientation = Orientation::kOriginal);

::testing::AssertionResult Compare(const YUVPlane& src, const YUVPlane& dec,
                                   uint32_t bit_depth,
                                   const std::string& file_name,
                                   float expected_distortion = 99.f,
                                   MetricType metric = PSNR);

//------------------------------------------------------------------------------

// Returns WP2_STATUS_INVALID_PARAMETER if the comparison fails.
WP2Status EncodeDecodeCompare(
    const std::string& file_name,
    const EncoderConfig& encoder_config = EncoderConfig::kDefault,
    const DecoderConfig& decoder_config = DecoderConfig::kDefault);

// Encodes the image with 'file_name' cropped by 'window' with the
// 'encoder_config'. Returns the distortion of decoding it with config 'b' minus
// the one with 'a'.
WP2Status GetDistortionDiff(const std::string& file_name,
                            const Rectangle& window,
                            const EncoderConfig& encoder_config,
                            const DecoderConfig& decoder_config_a,
                            const DecoderConfig& decoder_config_b,
                            float disto_diff[5], MetricType metric = PSNR);
// Encodes the image with 'file_name'cropped by 'window'  with encoder config
// 'a' and 'b''. Returns the distortion of encoding (then decoding) the image
// 'file_name' with config 'b' minus the one with 'a'.
WP2Status GetDistortionDiff(const std::string& file_name,
                            const Rectangle& window,
                            const EncoderConfig& encoder_config_a,
                            const EncoderConfig& encoder_config_b,
                            const DecoderConfig& decoder_config,
                            float disto_diff[5], MetricType metric = PSNR);

// Simple one-liner for debugging; avoids the inclusion of image_enc.h.
void DumpImage(const ArgbBuffer& image, const std::string& file_path);

//------------------------------------------------------------------------------

// Creates a random block partition of a tile of size w x h.
void CreatePartition(uint32_t w, uint32_t h, bool snapped,
                     UniformIntDistribution* const gen,
                     VectorNoCtor<Block>* const blocks);

//------------------------------------------------------------------------------
// struct for speeding up pseudo-random filling of buffers

template<size_t N, typename T> struct PrecalculatedRandom {
 public:
  PrecalculatedRandom(uint32_t min, uint32_t max) {
    for (auto& b : buffer) b = min + rng.GetSigned(max - min);
  }
  void Fill(T* dst, uint32_t w, uint32_t h, uint32_t step) {
    assert(w < N);
    for (uint32_t j = 0; j < h; ++j) {
      const uint32_t offset = GetUnsigned(N - w);
      for (uint32_t i = 0; i < w; ++i) dst[i] = buffer[offset + i];
      dst += step;
    }
  }
  void Fill(T* dst, uint32_t w) { Fill(dst, w, 1, /*step=*/0); }

  // this version is not using the precalc'd buffer[], so is slower
  void FreshFill(uint32_t min, uint32_t max,
                 T* dst, uint32_t w, uint32_t h, uint32_t step) {
    for (uint32_t y = 0; y < h; ++y) {
      for (uint32_t x = 0; x < w; ++x) {
        dst[x] = min + rng.GetUnsigned(max - min);
      }
      dst += step;
    }
  }

  uint32_t GetUnsigned(uint32_t max) { return rng.GetUnsigned(max); }

 private:
  T buffer[N];
  PseudoRNG rng;
};

//------------------------------------------------------------------------------

static inline uint32_t SillyCRC(uint32_t crc, uint32_t v) {
  v += 1;
  return ((crc >> 3) + (17 + v)) ^ ((crc << 5) * (23 + 7 * v));
}

}  // namespace testing
}  // namespace WP2

// Displays an error message instead of a number, even in ASSERT_EQ() output.
inline std::ostream& operator<<(std::ostream& out, const WP2Status& status) {
  out << WP2GetStatusMessage(status);
  return out;
}

// Displays EncoderConfig as code that is ready-to-paste, for easy reproduction.
std::ostream& operator<<(std::ostream& out, const WP2::EncoderConfig& config);
std::ostream& operator<<(std::ostream& out, const WP2::DecoderConfig& config);

//------------------------------------------------------------------------------

#endif  // WP2_TESTS_INCLUDE_HELPERS_H_
