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

// Measure wp2 codec format pessimization level by computing the ratio of random
// and mutated bitstreams that decompress fine compared to those that are
// invalid.

#include <cstdio>
#include <limits>
#include <tuple>

#include "imageio/image_dec.h"
#include "include/helpers.h"
#include "src/common/constants.h"
#include "src/dec/wp2_dec_i.h"
#include "src/enc/wp2_enc_i.h"
#include "src/utils/random.h"
#include "src/wp2/decode.h"
#include "src/wp2/encode.h"

// Define WP2_ERROR_TRACE to gather stats about where the random bitstreams fail
// to decode. Do not use multithreading, even only at tile or test level.
#if defined(WP2_TRACE) || defined(WP2_ERROR_TRACE)
#include <unordered_map>

#include "examples/example_utils.h"

class ErrorTracer : public WP2ErrorTracer {
 public:
  ErrorTracer() { wp2_error_tracer = this; }
  ~ErrorTracer() {
    wp2_error_tracer = nullptr;
    WriteToFile("/tmp/wp2_error_trace.txt");
  }

  void Log(const char* const file, int line, WP2Status, const char*) override {
    // Lighten the callstack.
    const std::string call = WP2GetFileName(file) + ":" + std::to_string(line);
    if (!callstack_.empty()) callstack_ += ", ";
    callstack_ += call;  // Reversed call order but more convenient text output.
  }
  void Flush() {
    ++callstack_count_[callstack_];
    callstack_.clear();
  }
  void WriteToFile(const char* const file_path) {
    FILE* const file = fopen(file_path, "w");
    if (file != nullptr) {
      for (const auto& it : callstack_count_) {
        fprintf(file, "%u, %s\n", it.second, it.first.c_str());
      }
      fclose(file);
    }
  }

 private:
  std::string callstack_;
  std::unordered_map<std::string, uint32_t> callstack_count_;
};
#else
class ErrorTracer {
 public:
  void Flush() {}
};
#endif  // defined(WP2_TRACE) || defined(WP2_ERROR_TRACE)

namespace WP2 {
namespace {

//------------------------------------------------------------------------------

// Inspired from GlobalAnalysis() but analysis-driven values are replaced by
// random ones.
WP2Status GenerateRandomGlobalParams(const EncoderConfig& config,
                                     bool image_has_alpha,
                                     UniformIntDistribution* const random,
                                     GlobalParams* const gparams) {
  // TODO(yguyon): Also test fully random GlobalParams.
  gparams->Reset();
  gparams->type_ = DecideGlobalParamsType(config);

  gparams->has_alpha_ = image_has_alpha && random->FlipACoin();
  if (gparams->type_ != GlobalParams::GP_LOSSLESS) {
    gparams->partition_set_ = config.partition_set;
    gparams->partition_snapping_ = config.partition_snapping;
    gparams->explicit_segment_ids_ =
        (config.segment_id_mode == WP2::EncoderConfig::SEGMENT_ID_AUTO)
            ? (config.quality > 63.f)
            : (config.segment_id_mode ==
               WP2::EncoderConfig::SEGMENT_ID_EXPLICIT);
    gparams->use_rnd_mtx_ = config.use_random_matrix;
    WP2_CHECK_STATUS(gparams->transf_.Init(config.csp_type));

    const uint32_t num_segments = random->Get(
        1u, std::min((uint32_t)config.segments,
                     GetMaxNumSegments(gparams->explicit_segment_ids_,
                                       GetQualityHint(config.quality),
                                       config.partition_set)));
    WP2_CHECK_ALLOC_OK(gparams->segments_.resize(num_segments));
    for (uint32_t i = 0; i < num_segments; ++i) {
      gparams->segments_[i].risk_ = random->Get(0, 127) / 127.f;
      gparams->segments_[i].risk_class_ = i;
      gparams->segments_[i].SetQuantizationFactor(
          gparams->transf_, GetQualityHint(config.quality),
          kDefaultQuantMultiplier, kDefaultQuantMultiplier,
          random->Get(1, 128) / 1.28f, random->Get(1, 128) / 1.28f * 40.f,
          ALL_RECTS);
      gparams->segments_[i].use_grain_ = false;  // Nah.
    }
    gparams->maybe_use_lossy_alpha_ =
        gparams->has_alpha_ && (config.speed > 0) &&
        (config.alpha_quality <= kMaxLossyQuality) && random->FlipACoin();
    gparams->enable_alpha_filter_ =
        gparams->maybe_use_lossy_alpha_ && config.enable_alpha_filter;
    WP2_CHECK_ALLOC_OK(gparams->a_segments_.resize(1));
    const float quant_factor =
        99.f * (kMaxLossyQuality - config.alpha_quality) / kMaxLossyQuality;
    const float lambda = std::max(85.f, -42.f * config.alpha_quality + 3400.f);
    gparams->a_segments_.front().SetQuantizationFactor(
        gparams->transf_, GetQualityHint(config.quality),
        /*u_quant_multiplier=*/0, /*v_quant_multiplier=*/0, quant_factor,
        lambda, ALL_RECTS);

    WP2_CHECK_STATUS(gparams->InitFixedPredictors());
  }
  return WP2_STATUS_OK;
}

//------------------------------------------------------------------------------

// Generates a random wp2 bitstream from 'seed' and 'config' to 'output'.
// It contains the header (defaults values), a semi-random GlobalParams and a
// randomly-sized tile chunk with random bytes.
WP2Status RandomEncode(size_t seed, const EncoderConfig& config,
                       Writer* const output) {
  EncoderConfig final_config = config;
  final_config.tile_shape = TILE_SHAPE_SQUARE_512;
  WP2_CHECK_OK(final_config.IsValid(), WP2_STATUS_INVALID_CONFIGURATION);
  WP2_CHECK_OK(!final_config.create_preview, WP2_STATUS_INVALID_CONFIGURATION);

  UniformIntDistribution random(seed);
  const uint32_t width = random.Get(1, 512);
  const uint32_t height = random.Get(1, 512);
  const bool has_alpha = random.FlipACoin();

  // Append proxy header data, only config and dimensions are relevant.
  WP2_CHECK_STATUS(EncodeHeader(final_config, width, height, has_alpha,
                                /*is_anim=*/false, /*loop_count=*/0,
                                /*background_color=*/kTransparentArgb38b,
                                /*preview_color=*/RGB12b(), /*has_icc=*/false,
                                /*has_trailing_data=*/false, output));

  // Append randomly generated global parameters.
  GlobalParams gparams;
  WP2_CHECK_STATUS(
      GenerateRandomGlobalParams(config, has_alpha, &random, &gparams));
  WP2_CHECK_STATUS(EncodeGLBL(config, gparams, GetQualityHint(config.quality),
                              has_alpha, output));

  // Append a randomly sized random bitstream representing a single tile.
  // Testing several tiles at once would not bring useful information as they
  // are isolated. Clamp to ~4 bytes per pixel.
  const uint32_t num_bytes = random.Get(10u, width * height * 4u);
  {
    uint8_t buf[kMaxVarIntLength];
    const size_t buf_size = WriteVarInt(num_bytes - 1, kChunkSizeMax, buf);
    WP2_CHECK_ALLOC_OK(output->Append(buf, buf_size));
  }
  Vector_u8 bytes;
  WP2_CHECK_ALLOC_OK(bytes.resize(num_bytes));
  for (uint8_t& byte : bytes) byte = random.Get<uint8_t>(0, 255);
  WP2_CHECK_ALLOC_OK(output->Append(bytes.data(), bytes.size()));
  return WP2_STATUS_OK;
}

//------------------------------------------------------------------------------

// For command line output.
struct BitstreamStats {
  size_t length = 0;           // Of the bitstream, header included.
  BitstreamFeatures features;  // Extracted header.
  bool is_valid = false;       // Can be decoded without an error.

  static void PrintTitleRow() {
    printf("Length   Width  Height  Bpp     Valid\n");
  }
  void PrintRow() const {
    printf("%7zu, %5u, %6u, %6.2f, %s\n", length, features.width,
           features.height, length / (float)(features.width * features.height),
           is_valid ? "Yes" : "No");
  }
};

// From a wp2 'bitstream', checks if it 'is_valid' and can extract 'stats'.
WP2Status GetBistreamStats(const uint8_t bitstream[], size_t bitstream_length,
                           bool* const is_valid,
                           BitstreamStats* const stats = nullptr) {
  ArgbBuffer buffer;
  const WP2Status status = Decode(bitstream, bitstream_length, &buffer);
  if (status != WP2_STATUS_BITSTREAM_ERROR) WP2_CHECK_STATUS(status);
  *is_valid = (status == WP2_STATUS_OK);

  if (stats != nullptr) {
    stats->length = bitstream_length;
    WP2_CHECK_STATUS(stats->features.Read(bitstream, bitstream_length));
    stats->is_valid = *is_valid;
  }
  return WP2_STATUS_OK;
}

//------------------------------------------------------------------------------

// Creates 'num_tries' random bitstreams and counts the number of valid ones.
WP2Status GetPessimizationProportion(
    size_t seed_generator_seed, const EncoderConfig& config, uint32_t num_tries,
    uint32_t* const num_valid,
    std::vector<BitstreamStats>* const all_stats = nullptr) {
  ErrorTracer error_tracer;
  UniformIntDistribution seed_generator(seed_generator_seed);
  if (num_valid != nullptr) *num_valid = 0;
  for (uint32_t i = 0; i < num_tries; ++i) {
    const size_t seed = seed_generator.Get(0u, 1234567890u);
    bool is_valid;
    if (all_stats != nullptr) all_stats->emplace_back();
    BitstreamStats* const stats =
        (all_stats != nullptr) ? &all_stats->back() : nullptr;
    MemoryWriter writer;
    WP2_CHECK_STATUS(RandomEncode(seed, config, &writer));
    WP2_CHECK_STATUS(
        GetBistreamStats(writer.mem_, writer.size_, &is_valid, stats));
    error_tracer.Flush();
    if (num_valid != nullptr && is_valid) ++*num_valid;
  }
  return WP2_STATUS_OK;
}

//------------------------------------------------------------------------------

// Tests that a valid bitstream, well, is.
TEST(PessimizationTest, NotRandom) {
  constexpr const char file_name[] = "source1_64x48.png";
  ArgbBuffer src;
  ASSERT_WP2_OK(ReadImage(testing::GetTestDataPath(file_name).c_str(), &src));
  MemoryWriter writer;
  ASSERT_WP2_OK(Encode(src, &writer));
  bool is_valid;
  ASSERT_WP2_OK(GetBistreamStats(writer.mem_, writer.size_, &is_valid));
  EXPECT_TRUE(is_valid);
}

// Try different configurations.
class PessimizationTest : public ::testing::TestWithParam<
                              std::tuple<size_t, float, uint32_t, bool>> {};

TEST_P(PessimizationTest, RandomBitstream) {
  const size_t seed = std::get<0>(GetParam());
  EncoderConfig config = EncoderConfig::kDefault;
  config.quality = config.alpha_quality = std::get<1>(GetParam());
  const uint32_t num_tries = std::get<2>(GetParam());
  const bool print_stats = std::get<3>(GetParam());

  uint32_t num_valid;
  std::vector<BitstreamStats> all_stats;
  ASSERT_WP2_OK(GetPessimizationProportion(seed, config, num_tries, &num_valid,
                                           print_stats ? &all_stats : nullptr));

  const float percentage_of_valid_bitstreams = 100.f * num_valid / num_tries;
  EXPECT_GE(percentage_of_valid_bitstreams, 0.f);
  // TODO(yguyon): 0% seems buggy?

  if (print_stats) {
    BitstreamStats::PrintTitleRow();
    for (const BitstreamStats& stats : all_stats) stats.PrintRow();
    printf("Valid bitstreams: %u/%u (%.3f %%)\n", num_valid, num_tries,
           percentage_of_valid_bitstreams);
  }
}

INSTANTIATE_TEST_SUITE_P(
    PessimizationTestInstantiation, PessimizationTest,
    ::testing::Combine(::testing::Values(42),           // seed
                       ::testing::Values(75.f, 100.f),  // quality
                       ::testing::Values(50),           // num tries
                       ::testing::Values(false)         // print stats
                       ));

//------------------------------------------------------------------------------

// For a given wp2 'bitstream', returns the position of the first tile byte
// (corresponding to the first byte of the variable integer of the tile chunk
// size).
WP2Status GetTileBytesPosition(const uint8_t bitstream[],
                               size_t bitstream_length,
                               size_t* const tile_bytes_position,
                               uint32_t* const tile_width,
                               uint32_t* const tile_height) {
  ExternalDataSource data_source(bitstream, bitstream_length);
  BitstreamFeatures features;
  WP2_CHECK_STATUS(DecodeHeader(&data_source, &features));
  WP2_CHECK_STATUS(SkipPreview(&data_source, features));
  WP2_CHECK_STATUS(SkipICC(&data_source, features));
  AnimationFrame frame;
  uint32_t frame_index = 0;
  WP2_CHECK_STATUS(DecodeANMF(&data_source, features, frame_index, &frame));
  GlobalParams gparams;
  WP2_CHECK_STATUS(
      DecodeGLBL(&data_source, DecoderConfig::kDefault, features, &gparams));
  *tile_bytes_position =
      data_source.GetNumDiscardedBytes() + data_source.GetNumReadBytes();
  *tile_width = features.tile_width;
  *tile_height = features.tile_height;
  return WP2_STATUS_OK;
}

// Mutates one byte randomly located in the 'bytes'.
void RandomByteChange(uint8_t bytes[], size_t num_bytes,
                      UniformIntDistribution* const random) {
  assert(num_bytes > 0);
  const uint32_t i = random->Get(0u, (uint32_t)num_bytes - 1u);
  bytes[i] = random->Get<uint8_t>(0, 255);
}

// Randomly changes something in the 'tile_bytes' while preserving the integrity
// of the chunk size.
WP2Status MutateBitstream(uint32_t tile_width, uint32_t tile_height,
                          uint8_t tile_bytes[], size_t tile_num_bytes,
                          UniformIntDistribution* const random) {
  ExternalDataSource data_source(tile_bytes, tile_num_bytes);
  const uint64_t max_size = kMaxNumBytesPerPixel * tile_width * tile_height + 1;
  size_t tile_chunk_size;
  WP2_CHECK_OK(ReadVarInt(&data_source, max_size, &tile_chunk_size),
               WP2_STATUS_NOT_ENOUGH_DATA);
  const size_t var_int_num_bytes = data_source.GetNumReadBytes();
  WP2_CHECK_OK(var_int_num_bytes + tile_chunk_size == tile_num_bytes,
               WP2_STATUS_NOT_ENOUGH_DATA);

  RandomByteChange(tile_bytes + var_int_num_bytes,
                   tile_num_bytes - var_int_num_bytes, random);
  // TODO(yguyon): Add more operations
  return WP2_STATUS_OK;
}

//------------------------------------------------------------------------------

// Mutate a valid bistream rather than generating completely random ones.
class ImageMutationTest
    : public ::testing::TestWithParam<
          std::tuple<const char*, size_t, float, uint32_t, bool>> {};

TEST_P(ImageMutationTest, FromValidImage) {
  const char* const file_name = std::get<0>(GetParam());
  const size_t seed = std::get<1>(GetParam());
  EncoderConfig config = EncoderConfig::kDefault;
  config.quality = config.alpha_quality = std::get<2>(GetParam());
  if (config.tile_shape == TILE_SHAPE_AUTO) {
    config.tile_shape = TILE_SHAPE_SQUARE_512;
  }
  const uint32_t max_tile_size = kMaxTileSize;
  const uint32_t num_tries = std::get<3>(GetParam());
  const bool print_stats = std::get<4>(GetParam());

  MemoryWriter writer;
  {
    // Read and encode an image, only keep the bitstream.
    ArgbBuffer src;
    ASSERT_WP2_OK(ReadImage(testing::GetTestDataPath(file_name).c_str(), &src));
    ASSERT_WP2_OK(src.SetView(src, {0, 0, std::min(src.width, max_tile_size),
                                    std::min(src.height, max_tile_size)}));
    ASSERT_WP2_OK(Encode(src, &writer, config));
  }
  bool is_valid = false;
  ASSERT_WP2_OK(GetBistreamStats(writer.mem_, writer.size_, &is_valid));
  ASSERT_TRUE(is_valid);  // No mutation happened yet.

  // Do not mutate anything prior to the tile data.
  size_t tile_bytes_position = 0;
  uint32_t tile_width = 0;
  uint32_t tile_height = 0;
  ASSERT_WP2_OK(GetTileBytesPosition(writer.mem_, writer.size_,
                                     &tile_bytes_position, &tile_width,
                                     &tile_height));

  // Mutate the valid bitstream 'num_tries' times (start from the original each
  // time).
  UniformIntDistribution random(seed);
  uint32_t num_valid = 0;
  Data data;
  ErrorTracer error_tracer;
  for (uint32_t i = 0; i < num_tries; ++i) {
    ASSERT_WP2_OK(data.CopyFrom(writer.mem_, writer.size_));
    ASSERT_WP2_OK(MutateBitstream(tile_width, tile_height,
                                  data.bytes + tile_bytes_position,
                                  data.size - tile_bytes_position, &random));
    ASSERT_WP2_OK(GetBistreamStats(data.bytes, data.size, &is_valid));
    error_tracer.Flush();
    if (is_valid) ++num_valid;
  }

  const float percentage_of_valid_bitstreams = 100.f * num_valid / num_tries;
  EXPECT_GE(percentage_of_valid_bitstreams, 0.f);
  // TODO(yguyon): 2% seems super low.

  if (print_stats) {
    printf("Valid mutations: %u/%u (%.3f %%), tile bytes: %zu/%zu\n", num_valid,
           num_tries, percentage_of_valid_bitstreams,
           writer.size_ - tile_bytes_position, writer.size_);
  }
}

INSTANTIATE_TEST_SUITE_P(
    ImageMutationTestInstantiationLossy, ImageMutationTest,
    ::testing::Combine(::testing::Values("source1_64x48.png",
                                         "source1_32x32.png"),
                       ::testing::Values(42),    // seed
                       ::testing::Values(75.f),  // quality
                       ::testing::Values(50),    // num tries
                       ::testing::Values(false)  // print stats
                       ));

INSTANTIATE_TEST_SUITE_P(
    ImageMutationTestInstantiationLossless, ImageMutationTest,
    ::testing::Combine(::testing::Values("source1_64x48.png",
                                         "source1_32x32.png"),
                       ::testing::Values(12),     // seed
                       ::testing::Values(100.f),  // quality
                       ::testing::Values(50),     // num tries
                       ::testing::Values(false)   // print stats
                       ));

//------------------------------------------------------------------------------

}  // namespace
}  // namespace WP2
