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

// Test premultiplied, unpremultiplied, reordered and/or opaque decoding output.
// There is no loss during format conversion because the source ('original') is
// premultiplied before encoding and the 'premul <-> unpremul <-> premul'
// operation is stable.

#include <string>
#include <tuple>
#include <vector>

#include "imageio/image_dec.h"
#include "include/helpers.h"
#include "src/dsp/dsp.h"
#include "src/wp2/decode.h"
#include "src/wp2/encode.h"

namespace WP2 {
namespace {

//------------------------------------------------------------------------------

// Copies data from 'src' to 'dst' and convert it if formats differ.
void Convert(const ArgbBuffer& src, ArgbBuffer* const dst) {
  if (src.format == WP2_Argb_32) {
    ASSERT_WP2_OK(dst->ConvertFrom(src));
  } else {
    ASSERT_EQ(dst->format, WP2_Argb_32);
    ASSERT_WP2_OK(dst->Resize(src.width, src.height));

    WP2ArgbConverterInit();
    for (uint32_t y = 0; y < src.height; ++y) {
      const uint8_t* const src_row = (const uint8_t*)src.GetRow(y);
      uint8_t* const dst_row = (uint8_t*)dst->GetRow(y);
      WP2ArgbConvertFrom[src.format](src_row, src.width, dst_row);
    }
  }
}

// Reference implementation, without the tricks used in argb_converter.cc
void Premultiply(uint8_t* argb, uint32_t width) {
  for (uint32_t i = 0; i < width; ++i, argb += 4) {
    const uint32_t alpha = argb[0];
    if (alpha < 255u) {
      argb[1] = (uint8_t)((argb[1] * alpha + 127u) / 255u);
      argb[2] = (uint8_t)((argb[2] * alpha + 127u) / 255u);
      argb[3] = (uint8_t)((argb[3] * alpha + 127u) / 255u);
    }
  }
}

// Copies alpha from 'src' to 'dst' (must both be WP2_Argb_32).
void CopyAlpha(const ArgbBuffer& src, ArgbBuffer* const dst, bool premultiply) {
  ASSERT_EQ(src.format, WP2_Argb_32);
  ASSERT_EQ(dst->format, WP2_Argb_32);
  ASSERT_EQ(src.width, dst->width);
  ASSERT_EQ(src.height, dst->height);

  for (uint32_t y = 0; y < src.height; ++y) {
    const uint8_t* const src_row = (const uint8_t*)src.GetRow(y);
    uint8_t* const dst_row = (uint8_t*)dst->GetRow(y);
    for (uint32_t x = 0; x < src.width; ++x) dst_row[x * 4] = src_row[x * 4];
    if (premultiply) Premultiply(dst_row, src.width);
  }
}

// Go back to Argb32 for comparison and copy alpha if it was lost.
ArgbBuffer ConvertOutputBackToOriginalFormat(const ArgbBuffer& original,
                                             const ArgbBuffer& output) {
  ArgbBuffer output_in_original_format(original.format);
  Convert(output, &output_in_original_format);
  if (output.format == WP2_RGB_24 || output.format == WP2_BGR_24) {
    CopyAlpha(original, &output_in_original_format, /*premultiply=*/true);
  } else if (output.format == WP2_XRGB_32 || output.format == WP2_RGBX_32 ||
             output.format == WP2_BGRX_32) {
    CopyAlpha(original, &output_in_original_format, /*premultiply=*/false);
  }
  return output_in_original_format;
}

//------------------------------------------------------------------------------
// Basic format conversion only tests

TEST(FormatTest, ConvertSimple) {
  const std::string& file_name = "1x1";
  ArgbBuffer original(WP2_Argb_32);
  ASSERT_WP2_OK(original.Resize(1, 1));
  uint8_t* const a = (uint8_t*)original.GetRow(0) + 0;
  uint8_t* const r = (uint8_t*)original.GetRow(0) + 1;
  uint8_t* const g = (uint8_t*)original.GetRow(0) + 2;
  uint8_t* const b = (uint8_t*)original.GetRow(0) + 3;
  *a = 0xAA;
  *r = 0x55;
  *g = 0xAA;
  *b = 0x00;

  for (WP2SampleFormat format = WP2_Argb_32; format < WP2_FORMAT_NUM;
       format = (WP2SampleFormat)(format + 1)) {
    if (WP2FormatBpc(format) != WP2FormatBpc(WP2_Argb_32)) continue;
    ArgbBuffer converted(format);
    Convert(original, &converted);

    ArgbBuffer reconverted =
        ConvertOutputBackToOriginalFormat(original, converted);
    ASSERT_TRUE(testing::Compare(original, reconverted, file_name));
  }
}

TEST(FormatTest, ConvertImage) {
  const std::string& file_name = "source1_64x48.png";
  ArgbBuffer original(WP2_Argb_32);
  ASSERT_WP2_OK(
      ReadImage(testing::GetTestDataPath(file_name).c_str(), &original));

  for (WP2SampleFormat format = WP2_Argb_32; format < WP2_FORMAT_NUM;
       format = (WP2SampleFormat)(format + 1)) {
    if (WP2FormatBpc(format) != WP2FormatBpc(WP2_Argb_32)) continue;
    ArgbBuffer converted(format);
    Convert(original, &converted);

    ArgbBuffer reconverted =
        ConvertOutputBackToOriginalFormat(original, converted);
    ASSERT_TRUE(testing::Compare(original, reconverted, file_name));
  }
}

//------------------------------------------------------------------------------
// Encode and decode tests

typedef std::tuple<std::vector<const char*>, WP2SampleFormat, float> Param;
class FormatTest : public ::testing::TestWithParam<Param> {};

//------------------------------------------------------------------------------

TEST_P(FormatTest, EncodeDecodeImage) {
  const char* src_file_name = std::get<0>(GetParam()).front();
  const WP2SampleFormat output_format = std::get<1>(GetParam());
  const float quality = std::get<2>(GetParam());

  ArgbBuffer src;
  MemoryWriter data;
  ASSERT_WP2_OK(
      testing::CompressImage(src_file_name, &data, &src, quality, /*speed=*/2));

  ArgbBuffer output(output_format);
  ASSERT_WP2_OK(Decode(data.mem_, data.size_, &output));

  ArgbBuffer output_original_format =
      ConvertOutputBackToOriginalFormat(src, output);
  ASSERT_TRUE(testing::Compare(src, output_original_format, src_file_name,
                               testing::GetExpectedDistortion(quality)));
}

//------------------------------------------------------------------------------

TEST_P(FormatTest, EncodeDecodeAnimation) {
  const std::vector<const char*>& file_names = std::get<0>(GetParam());
  const std::vector<uint32_t> durations_ms(file_names.size(), 100u);
  const WP2SampleFormat output_format = std::get<1>(GetParam());
  const float quality = std::get<2>(GetParam());
  const uint32_t speed = 0, thread_level = 0;

  std::vector<ArgbBuffer> frames;
  MemoryWriter encoded_data;
  ASSERT_WP2_OK(testing::CompressAnimation(file_names, durations_ms,
                                           &encoded_data, &frames, quality,
                                           speed, thread_level));

  ArgbBuffer decoded_frame(output_format);
  ArrayDecoder decoder(encoded_data.mem_, encoded_data.size_,
                       DecoderConfig::kDefault, &decoded_frame);

  for (size_t i = 0; i < frames.size(); ++i) {
    ASSERT_TRUE(decoder.ReadFrame());

    ArgbBuffer decoded_frame_original_format =
        ConvertOutputBackToOriginalFormat(frames[i], decoded_frame);
    EXPECT_TRUE(testing::Compare(frames[i], decoded_frame_original_format,
                                 file_names[i],
                                 testing::GetExpectedDistortion(quality)));
  }
}

//------------------------------------------------------------------------------

constexpr float kMinExpectedSimilarityForPreview = 10.f;
constexpr float kMinExpectedSimilarityForPreviewColor = 5.f;

TEST_P(FormatTest, EncodeDecodePreview) {
  const char* file_name = std::get<0>(GetParam()).front();
  const WP2SampleFormat output_format = std::get<1>(GetParam());
  const float quality = std::get<2>(GetParam());

  ArgbBuffer original_image;
  ASSERT_WP2_OK(
      ReadImage(testing::GetTestDataPath(file_name).c_str(), &original_image));

  EncoderConfig config;
  config.quality = quality;
  config.speed = 2;
  config.create_preview = true;
  MemoryWriter memory_writer;
  EXPECT_WP2_OK(Encode(original_image, &memory_writer, config));

  ArgbBuffer decompressed_preview(output_format);
  ASSERT_WP2_OK(ExtractPreview(memory_writer.mem_, memory_writer.size_,
                               &decompressed_preview));

  ArgbBuffer decompressed_preview_original_format =
      ConvertOutputBackToOriginalFormat(original_image, decompressed_preview);
  float disto[5];
  ASSERT_WP2_OK(
      decompressed_preview_original_format.GetDistortion(original_image,
                                                         PSNR, disto));
  EXPECT_GE(disto[4], config.create_preview
                          ? kMinExpectedSimilarityForPreview
                          : kMinExpectedSimilarityForPreviewColor);
}

//------------------------------------------------------------------------------
// Input images are expected to be of equal sizes inside a Param.

INSTANTIATE_TEST_SUITE_P(
    FormatTestInstantiation, FormatTest,
    ::testing::Values(
        Param({"source1_64x48.png"}, WP2_ARGB_32, /*quality=*/100.0f),
        Param({"source1_64x48.png"}, WP2_rgbA_32, /*quality=*/40.0f)));

INSTANTIATE_TEST_SUITE_P(
    ExhaustiveFormatTestInstantiation, FormatTest,
    ::testing::Combine(::testing::Values(std::vector<const char*>{
                           "source1_1x1.png"}),
                       ::testing::Range(WP2_Argb_32, WP2_Argb_38),
                       /*quality=*/::testing::Values(10.f, 100.f)));

// This one takes a while to run so it is disabled.
// Can still be run with flag --test_arg=--gunit_also_run_disabled_tests
INSTANTIATE_TEST_SUITE_P(
    DISABLED_HeavyFormatTestInstantiation, FormatTest,
    ::testing::Combine(
        ::testing::Values(
            std::vector<const char*>{"source1_64x48.png"},
            std::vector<const char*>{"alpha_ramp.png", "alpha_ramp.lossy.webp",
                                     "alpha_ramp.webp"},
            std::vector<const char*>{"source0.pgm", "source1.png",
                                     "source2.tiff", "source3.jpg"}),
        ::testing::Range(WP2_Argb_32, WP2_Argb_38),
        ::testing::Values(0.f, 70.f, 96.f)));

//------------------------------------------------------------------------------

static struct FmtTestCase {
  WP2SampleFormat fmt;
  const char* name;
  size_t seed;
  bool has_alpha;
} kFmtTestCases[] = {
  { WP2_Argb_32, "WP2_Argb_32", 412412, true },
  { WP2_ARGB_32, "WP2_ARGB_32", 655473, true },
  { WP2_XRGB_32, "WP2_XRGB_32", 653134, false },
  { WP2_rgbA_32, "WP2_rgbA_32", 812463, true },
  { WP2_RGBA_32, "WP2_RGBA_32", 547271, true },
  { WP2_RGBX_32, "WP2_RGBX_32", 872353, false },
  { WP2_bgrA_32, "WP2_bgrA_32", 274353, true },
  { WP2_BGRA_32, "WP2_BGRA_32", 897534, true },
  { WP2_BGRX_32, "WP2_BGRX_32", 414125, false },
};

typedef std::tuple<FmtTestCase, testing::WP2CPUInfoStruct> TestParam;
static void PrintTo(const TestParam& p, std::ostream* os) {
  *os << "{" << std::get<0>(p).name << ", " << std::get<1>(p).name << "}";
}

class SpeedTest : public ::testing::TestWithParam<TestParam> {
  void SetUp() override {
    WP2DspReset();
    WP2GetCPUInfo = std::get<1>(GetParam()).cpu_info;
    WP2ArgbConverterInit();
    WP2PSNRInit();
    test_ = std::get<0>(GetParam());
    bpp_ = WP2FormatBpp(test_.fmt);
    convert_to_ = WP2ArgbConvertTo[test_.fmt];
    convert_from_ = WP2ArgbConvertFrom[test_.fmt];
  }

 protected:
#if defined(NDEBUG)
  static constexpr uint32_t kNumTest = 1000u;
#else
  static constexpr uint32_t kNumTest = 100u;
#endif
  static constexpr uint32_t kMaxWidth = 1024u;
  static constexpr uint32_t kMaxHeight = 64u;
  static constexpr int32_t kMaxSize = kMaxWidth * kMaxHeight;
  std::vector<uint8_t> argb_in_, tmp_, argb_out_;
  FmtTestCase test_;
  uint32_t bpp_;
  bool no_alpha_;
  WP2ArgbConverterF convert_to_, convert_from_;

  void Init(UniformIntDistribution* const random) {
    argb_in_.resize(4 * kMaxSize, 0);
    for (uint32_t i = 0; i < kMaxSize; ++i) {
      const uint8_t A = test_.has_alpha ? random->Get(0u, 255u) : 0xffu;
      argb_in_[4 * i + 0] = A;
      argb_in_[4 * i + 1] = (A > 0) ? random->Get<uint32_t>(0u, A) : 0u;
      argb_in_[4 * i + 2] = (A > 0) ? random->Get<uint32_t>(0u, A) : 0u;
      argb_in_[4 * i + 3] = (A > 0) ? random->Get<uint32_t>(0u, A) : 0u;
    }
    argb_out_ = argb_in_;
    tmp_.resize(bpp_ * kMaxWidth);
  }

  bool DoTest() {
    UniformIntDistribution random(test_.seed);
    Init(&random);
    for (uint32_t test = 0; test < kNumTest; ++test) {
      const uint32_t width = 1 + random.Get(0u, kMaxWidth - 1u);
      const size_t stride = 4 * width;
      for (uint32_t y = 0; y < kMaxHeight; ++y) {
        convert_to_(&argb_in_[y * stride], width, &tmp_[0]);
        convert_from_(&tmp_[0], width, &argb_out_[y * stride]);
      }
    }
    if (WP2SumSquaredError8u(&argb_in_[0], &argb_out_[0], kMaxSize * 4)) {
      for (int i = 0; i < 12; ++i) {
        printf("[%d %d %d]", argb_in_[i], tmp_[i], argb_out_[i]);
      }
      printf("\n");
      return false;
    }
    return true;
  }
};

TEST_P(SpeedTest, TestBackForth) {
  EXPECT_TRUE(DoTest());
}

INSTANTIATE_TEST_SUITE_P(SpeedTestInstantiation, SpeedTest,
    ::testing::Combine(::testing::ValuesIn(kFmtTestCases),
                       ::testing::ValuesIn(testing::kWP2CpuInfoStructs)));

//------------------------------------------------------------------------------

}  // namespace
}  // namespace WP2
