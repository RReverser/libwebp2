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

#include "./helpers.h"

#include <cassert>

#include "examples/example_utils.h"
#include "imageio/anim_image_dec.h"
#include "imageio/image_dec.h"
#include "imageio/image_enc.h"
#include "src/utils/front_mgr.h"
#include "src/utils/orientation.h"

namespace WP2 {
namespace testing {

//------------------------------------------------------------------------------

std::string GetTestDataPath(const std::string& file_name) {
  return "testdata/" + file_name;
}

void GetTestDataPaths(const std::vector<const char*>& file_names,
                      std::vector<std::string>* const file_paths) {
  file_paths->resize(file_names.size());
  for (size_t i = 0; i < file_names.size(); ++i) {
    file_paths->at(i) = GetTestDataPath(file_names[i]);
  }
}

std::string GetTempDataPath(const std::string& file_path) {
  return "/tmp/" + WP2GetFileName(file_path);
}

//------------------------------------------------------------------------------

bool HasSameData(const Data& a, const Data& b) {
  if (a.size != b.size) return false;
  if (a.size == 0) return true;
  if (a.bytes == nullptr || b.bytes == nullptr) return false;
  return (std::memcmp((void*)a.bytes, (void*)b.bytes, a.size) == 0);
}

//------------------------------------------------------------------------------

WP2Status ReadImages(const std::vector<std::string>& file_paths,
                     std::vector<ArgbBuffer>* const frames) {
  frames->resize(file_paths.size());
  for (size_t i = 0; i < frames->size(); ++i) {
    WP2_CHECK_STATUS(ReadImage(file_paths[i].c_str(), &frames->at(i)));
  }
  return WP2_STATUS_OK;
}

WP2Status ReadAnimation(const uint8_t* const data, size_t size,
                        std::vector<ArgbBuffer>* const frames,
                        std::vector<uint32_t>* const frame_durations,
                        uint32_t* const loop_count) {
  ArgbBuffer buffer;
  ImageReader image_reader(data, size, &buffer);
  frames->clear();
  if (frame_durations != nullptr) frame_durations->clear();

  bool is_last_frame = false;
  do {
    uint32_t duration_ms = 0;
    WP2_CHECK_STATUS(image_reader.ReadFrame(&is_last_frame, &duration_ms));
    assert(!buffer.IsEmpty() && !buffer.IsView());
    frames->emplace_back(buffer.format);
    WP2_CHECK_STATUS(frames->back().CopyFrom(buffer));
    if (frame_durations != nullptr) frame_durations->push_back(duration_ms);
  } while (!is_last_frame);

  if (loop_count != nullptr) *loop_count = image_reader.GetLoopCount();
  return WP2_STATUS_OK;
}

WP2Status ReadAnimation(const std::string& file_path,
                        std::vector<ArgbBuffer>* const frames,
                        std::vector<uint32_t>* const frame_durations,
                        uint32_t* const loop_count) {
  Data data;
  WP2_CHECK_STATUS(IoUtilReadFile(file_path.c_str(), &data));
  WP2_CHECK_STATUS(ReadAnimation(data.bytes, data.size, frames, frame_durations,
                                 loop_count));
  return WP2_STATUS_OK;
}

//------------------------------------------------------------------------------

WP2Status CompressImage(const std::string& file_name,
                        Writer* const encoded_data, ArgbBuffer* original_image,
                        float quality, int speed, int thread_level,
                        uint32_t num_downsamplings) {
  ArgbBuffer internal_original_image;
  if (original_image == nullptr) original_image = &internal_original_image;
  WP2_CHECK_STATUS(
      ReadImage(GetTestDataPath(file_name).c_str(), original_image));
  for (uint32_t i = 0; i < num_downsamplings; ++i) {
    original_image->SimpleHalfDownsample();
  }

  EncoderConfig config;
  config.quality = quality;
  config.speed = speed;
  config.thread_level = thread_level;
  return Encode(*original_image, encoded_data, config);
}

WP2Status CompressAnimation(const std::vector<const char*>& frames_file_names,
                            const std::vector<uint32_t>& durations_ms,
                            Writer* const encoded_data,
                            std::vector<ArgbBuffer>* frames, float quality,
                            int speed, int thread_level,
                            uint32_t num_downsamplings) {
  assert(frames_file_names.size() == durations_ms.size());
  std::vector<std::string> file_paths;
  GetTestDataPaths(frames_file_names, &file_paths);
  std::vector<ArgbBuffer> internal_frames;
  if (frames == nullptr) frames = &internal_frames;
  WP2_CHECK_STATUS(ReadImages(file_paths, frames));
  for (uint32_t i = 0; i < num_downsamplings; ++i) {
    for (ArgbBuffer& frame : *frames) frame.SimpleHalfDownsample();
  }

  AnimationEncoder animation_encoder;
  for (size_t i = 0; i < frames->size(); ++i) {
    WP2_CHECK_STATUS(
        animation_encoder.AddFrame(frames->at(i), durations_ms[i]));
  }

  EncoderConfig config;
  config.quality = quality;
  config.speed = speed;
  config.thread_level = thread_level;
  return animation_encoder.Encode(encoded_data, config);
}

//------------------------------------------------------------------------------

float GetExpectedDistortion(float quality, float alpha_quality) {
  if (quality <= kMaxLossyQuality || alpha_quality <= kMaxLossyQuality) {
    return (15.f + quality * 0.06f + alpha_quality * 0.02f);  // Lossy
  } else if (quality < kMaxQuality || alpha_quality < kMaxQuality) {
    quality -= kMaxLossyQuality;
    alpha_quality -= kMaxLossyQuality;
    return (30.f + 1.f * quality + 1.f * alpha_quality);  // Near-lossless
  } else {
    return 99.f;  // Lossless
  }
}

float GetExpectedDistortion(const EncoderConfig& encoder_config) {
  return GetExpectedDistortion(encoder_config.quality,
                               encoder_config.alpha_quality);
}

::testing::AssertionResult Compare(const ArgbBuffer& src, const ArgbBuffer& dec,
                                   const std::string& file_name,
                                   float expected_distortion, MetricType metric,
                                   Orientation decoding_orientation) {
  if (src.width != RotateWidth(decoding_orientation, dec.width, dec.height) ||
      src.height != RotateHeight(decoding_orientation, dec.width, dec.height)) {
    return ::testing::AssertionFailure()
           << "Original image dimensions " << src.width << "x" << src.height
           << " are different than the decompressed ones " << dec.width << "x"
           << dec.height << " with file " << file_name;
  }

  WP2::ArgbBuffer dec_oriented(dec.format);
  if (decoding_orientation != WP2::Orientation::kOriginal) {
    if (WP2::RotateBuffer(WP2::GetInverseOrientation(decoding_orientation), dec,
                          &dec_oriented) != WP2_STATUS_OK) {
      return ::testing::AssertionFailure() << "Failed to rotate buffer";
    }
  } else {
    if (dec_oriented.SetView(dec) != WP2_STATUS_OK) {
      return ::testing::AssertionFailure() << "Failed to set view";
    }
  }

  float disto[5];
  if (dec_oriented.GetDistortion(src, metric, disto) != WP2_STATUS_OK) {
    return ::testing::AssertionFailure()
           << "Can not get distortion with file " << file_name;
  }
  if (disto[4] < expected_distortion) {
    auto assertion_failure =
        ::testing::AssertionFailure()
        << "Distortion is too high: " << disto[0] << ", " << disto[1] << ", "
        << disto[2] << ", " << disto[3] << ", " << disto[4] << " with file "
        << file_name << " (expected at least " << expected_distortion << ")";
    for (uint32_t y = 0; y < src.height; ++y) {
      const uint8_t* a = (const uint8_t*)src.GetRow(y);
      const uint8_t* b = (const uint8_t*)dec_oriented.GetRow(y);
      for (uint32_t x = 0; x < src.width; ++x) {
        const int32_t aa = a[0], ar = a[1], ag = a[2], ab = a[3];
        const int32_t ba = b[0], br = b[1], bg = b[2], bb = b[3];
        if (aa != ba || ar != br || ag != bg || ab != bb) {
          assertion_failure << std::endl
                            << "First different pixel at " << x << ", " << y
                            << ":   " << aa << ", " << ar << ", " << ag << ", "
                            << ab << " vs " << ba << ", " << br << ", " << bg
                            << ", " << bb;
          return assertion_failure;  // Don't spam.
        }
        a += 4;
        b += 4;
      }
    }
    return assertion_failure;
  }
  return ::testing::AssertionSuccess();
}

::testing::AssertionResult Compare(const YUVPlane& src, const YUVPlane& dec,
                                   uint32_t bit_depth,
                                   const std::string& file_name,
                                   float expected_distortion,
                                   MetricType metric) {
  // some easy checks
  if ((src.Y.w_ != dec.Y.w_ || src.Y.h_ != dec.Y.h_) ||
      (src.U.w_ != dec.U.w_ || src.U.h_ != dec.U.h_) ||
      (src.V.w_ != dec.V.w_ || src.V.h_ != dec.V.h_) ||
      src.A.IsEmpty() != dec.A.IsEmpty()) {
    return ::testing::AssertionFailure() << "Different image dimensions";
  }

  float disto[5];
  if (dec.GetDistortion(src, bit_depth, metric, disto) != WP2_STATUS_OK) {
    return ::testing::AssertionFailure()
           << "Can not get distortion with file " << file_name;
  }
  if (disto[4] < expected_distortion) {
    auto assertion_failure =
        ::testing::AssertionFailure()
        << "Distortion is too high: " << disto[0] << ", " << disto[1] << ", "
        << disto[2] << ", " << disto[3] << ", " << disto[4] << " with file "
        << file_name << " (expected at least " << expected_distortion << ")";
    for (Channel channel : {kYChannel, kUChannel, kVChannel, kAChannel}) {
      const Plane16& src_plane = src.GetChannel(channel);
      const Plane16& dec_plane = dec.GetChannel(channel);
      for (uint32_t y = 0; y < src_plane.h_; ++y) {
        for (uint32_t x = 0; x < src_plane.w_; ++x) {
          if (src_plane.At(x, y) != dec_plane.At(x, y)) {
            assertion_failure << std::endl
                              << "First different pixel in "
                              << kChannelStr[channel] << " at " << x << ", "
                              << y << ":   " << src_plane.At(x, y) << " vs "
                              << dec_plane.At(x, y);
            y = src_plane.h_;
            break;  // Don't spam.
          }
        }
      }
    }
    return assertion_failure;
  }
  return ::testing::AssertionSuccess();
}

//------------------------------------------------------------------------------

WP2Status EncodeDecodeCompare(const std::string& file_name,
                              const EncoderConfig& encoder_config,
                              const DecoderConfig& decoder_config) {
  ArgbBuffer original;
  WP2_CHECK_STATUS(
      ReadImage(testing::GetTestDataPath(file_name).c_str(), &original));

  MemoryWriter memory_writer;
  WP2_CHECK_STATUS(Encode(original, &memory_writer, encoder_config));

  ArgbBuffer decompressed;
  WP2_CHECK_STATUS(Decode(memory_writer.mem_, memory_writer.size_,
                          &decompressed, decoder_config));

  WP2_CHECK_OK(testing::Compare(original, decompressed, file_name,
                                testing::GetExpectedDistortion(encoder_config)),
               WP2_STATUS_INVALID_PARAMETER);
  return WP2_STATUS_OK;
}

WP2Status GetDistortionDiff(const std::string& file_name,
                            const Rectangle& window,
                            const EncoderConfig& encoder_config,
                            const DecoderConfig& decoder_config_a,
                            const DecoderConfig& decoder_config_b,
                            float disto_diff[5], MetricType metric) {
  ArgbBuffer original, cropped;
  WP2_CHECK_STATUS(ReadImage(GetTestDataPath(file_name).c_str(), &original));
  WP2_CHECK_STATUS(cropped.SetView(original, window));
  MemoryWriter writer;
  WP2_CHECK_STATUS(Encode(cropped, &writer, encoder_config));
  ArgbBuffer a, b;
  WP2_CHECK_STATUS(Decode(writer.mem_, writer.size_, &a, decoder_config_a));
  WP2_CHECK_STATUS(Decode(writer.mem_, writer.size_, &b, decoder_config_b));
  float disto_a[5], disto_b[5];
  WP2_CHECK_STATUS(a.GetDistortion(cropped, metric, disto_a));
  WP2_CHECK_STATUS(b.GetDistortion(cropped, metric, disto_b));
  for (int i = 0; i < 5; ++i) disto_diff[i] = disto_b[i] - disto_a[i];
  return WP2_STATUS_OK;
}

WP2Status GetDistortionDiff(const std::string& file_name,
                            const Rectangle& window,
                            const EncoderConfig& encoder_config_a,
                            const EncoderConfig& encoder_config_b,
                            const DecoderConfig& decoder_config,
                            float disto_diff[5], MetricType metric) {
  ArgbBuffer original, cropped;
  WP2_CHECK_STATUS(ReadImage(GetTestDataPath(file_name).c_str(), &original));
  WP2_CHECK_STATUS(cropped.SetView(original, window));
  MemoryWriter writer;
  ArgbBuffer a, b;
  WP2_CHECK_STATUS(Encode(cropped, &writer, encoder_config_a));
  WP2_CHECK_STATUS(Decode(writer.mem_, writer.size_, &a, decoder_config));
  writer.Reset();
  WP2_CHECK_STATUS(Encode(cropped, &writer, encoder_config_b));
  WP2_CHECK_STATUS(Decode(writer.mem_, writer.size_, &b, decoder_config));
  float disto_a[5], disto_b[5];
  WP2_CHECK_STATUS(a.GetDistortion(cropped, metric, disto_a));
  WP2_CHECK_STATUS(b.GetDistortion(cropped, metric, disto_b));
  for (int i = 0; i < 5; ++i) disto_diff[i] = disto_b[i] - disto_a[i];
  return WP2_STATUS_OK;
}

void DumpImage(const ArgbBuffer& image, const std::string& file_path) {
  WP2_ASSERT_STATUS(SaveImage(image, file_path.c_str(), /*overwrite=*/true));
}

//------------------------------------------------------------------------------

// Uniformly draws a random BlockSize whose width/height fit in bw/bh (in
// kMinBlockSize units).
static WP2::BlockSize UniformBlockSizeDrawInBounds(
    uint32_t bw, uint32_t bh, UniformIntDistribution* const gen) {
  BlockSize dims[BLK_LAST + 1];
  uint32_t size = 0;
  for (uint32_t dim = 0; (BlockSize)dim != BLK_LAST; ++dim) {
    if (WP2::BlockWidth[dim] <= bw && WP2::BlockHeight[dim] <= bh) {
      dims[size++] = (BlockSize)dim;
    }
  }
  return dims[gen->Get(0u, size - 1)];
}

//------------------------------------------------------------------------------

void CreatePartition(uint32_t w, uint32_t h, bool snapped,
                     UniformIntDistribution* const gen,
                     VectorNoCtor<Block>* const blocks) {
  const uint32_t bwidth = SizeBlocks(w), bheight = SizeBlocks(h);
  WP2::FrontMgrLexico m;
  WP2_ASSERT_STATUS(m.Init(ALL_RECTS, snapped, w, h));
  while (!m.Done()) {
    // Find the beginning of the next block in lexicographic order.
    uint32_t bx = 0;
    uint32_t by = bheight;
    for (uint32_t x = 0; x < bwidth; ++x) {
      const uint32_t y = m.GetOccupancy(x);
      if (y < by) {
        bx = x;
        by = y;
      }
    }
    // Find max width/height.
    uint32_t max_bw = 1, max_bh = std::min(kMaxBlockSize, bheight - by);
    while (bx + max_bw < bwidth && max_bw < kMaxBlockSize &&
           m.GetOccupancy(bx + max_bw) == m.GetOccupancy(bx)) {
      ++max_bw;
    }

    if (snapped) {
      const BlockSize snapped_block_size =
          GetSnappedBlockSize(bx, by, max_bw, max_bh);
      max_bw = BlockWidth[snapped_block_size];
      max_bh = BlockHeight[snapped_block_size];
    }
    const BlockSize dim = UniformBlockSizeDrawInBounds(max_bw, max_bh, gen);
    const uint32_t bh = BlockHeight[dim];

    const Block block(bx, by, (BlockSize)dim);

    if (!blocks->push_back(block)) assert(false);
    m.Use(block);
    if (!m.IsDone(bx, by)) assert(false);
    if (m.IsDone(bx, by + bh)) assert(false);
  }
}

//------------------------------------------------------------------------------

}  // namespace testing
}  // namespace WP2

#define PRINT_SETTING1(var_name)                                             \
  do {                                                                       \
    out << std::endl << "config." #var_name " = " << config.var_name << ";"; \
    if (config.var_name != default_config.var_name) {                        \
      out << "  // Default is " << default_config.var_name;                  \
    }                                                                        \
  } while (0);

#define PRINT_SETTING2(var_name, var_type)                                  \
  do {                                                                      \
    out << std::endl                                                        \
        << "config." #var_name " = (" #var_type ")" << (int)config.var_name \
        << ";";                                                             \
    if (config.var_name != default_config.var_name) {                       \
      out << "  // Default is " << (int)default_config.var_name;            \
    }                                                                       \
  } while (0);

#define PRINT_SETTINGP(var_name, var_type)                         \
  do {                                                             \
    out << std::endl << "config." #var_name " = ";                 \
    if (config.var_name == nullptr) {                              \
      out << "nullptr;";                                           \
    } else {                                                       \
      out << "[\" #var_type \" instance];  // Default is nullptr"; \
    }                                                              \
  } while (0);

std::ostream& operator<<(std::ostream& out, const WP2::EncoderConfig& config) {
  const WP2::EncoderConfig default_config = WP2::EncoderConfig::kDefault;
  out << std::endl << "EncoderConfig config = EncoderConfig::kDefault;";
  PRINT_SETTING1(quality);
  PRINT_SETTING1(target_size);
  PRINT_SETTING1(target_psnr);
  PRINT_SETTING1(alpha_quality);
  PRINT_SETTING1(speed);
  PRINT_SETTING2(decoding_orientation, Orientation);
  PRINT_SETTING2(create_preview, bool);
  PRINT_SETTING2(transfer_function, TransferFunction);
  PRINT_SETTING1(pass);
  PRINT_SETTING2(sns, float);
  PRINT_SETTING1(error_diffusion);
  PRINT_SETTING1(segments);
  PRINT_SETTING2(partition_method, PartitionMethod);
  PRINT_SETTING2(partition_set, PartitionSet);
  PRINT_SETTING2(partition_snapping, bool);
  PRINT_SETTING2(csp_type, Csp);
  PRINT_SETTING2(uv_mode, EncoderConfig::UVMode);
  PRINT_SETTING1(preprocessing);
  PRINT_SETTING1(preprocessing_strength);
  PRINT_SETTING2(use_random_matrix, bool);
  PRINT_SETTING2(store_grain, bool);
  PRINT_SETTING2(use_delta_palette, bool);
  PRINT_SETTING1(thread_level);
  PRINT_SETTING2(low_memory, bool);
  PRINT_SETTING1(use_neural_compression);
  PRINT_SETTINGP(info, EncoderInfo);
  PRINT_SETTING2(enable_alpha_filter, bool);
  out << std::endl;
  return out;
}

std::ostream& operator<<(std::ostream& out, const WP2::DecoderConfig& config) {
  const WP2::DecoderConfig default_config = WP2::DecoderConfig::kDefault;
  out << std::endl << "DecoderConfig config = DecoderConfig::kDefault;";
  PRINT_SETTING1(thread_level);
  PRINT_SETTING2(enable_deblocking_filter, bool);
  PRINT_SETTING2(enable_directional_filter, bool);
  PRINT_SETTING2(enable_restoration_filter, bool);
  PRINT_SETTING1(grain_amplitude);
  PRINT_SETTING2(incremental_mode, DecoderConfig::IncrementalMode);
  PRINT_SETTINGP(progress_hook, ProgressHook);
  PRINT_SETTINGP(info, DecoderInfo);
  out << std::endl;
  return out;
}

#undef PRINT_SETTING1
#undef PRINT_SETTING2
#undef PRINT_SETTINGP
