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

#include "include/helpers.h"
#include "src/wp2/base.h"
#include "src/wp2/decode.h"
#include "src/wp2/encode.h"

namespace WP2 {
namespace {

//------------------------------------------------------------------------------
// Helper functions to generate sample bitstreams that can be quickly decoded.

Data GenerateBitstream() {
  MemoryWriter writer;
  if (testing::CompressImage("test_exif_xmp.webp", &writer, nullptr,
                             EncoderConfig::kDefault.quality,
                             EncoderConfig::kDefault.speed,
                             EncoderConfig::kDefault.thread_level,
                             /*num_downsamplings=*/2) != WP2_STATUS_OK) {
    abort();
  }
  Data data;  // Transfer ownership.
  using std::swap;
  swap(data.bytes, writer.mem_);
  swap(data.size, writer.size_);
  writer.max_size_ = 0;
  return data;
}

Data GenerateAnimatedBitstream() {
  MemoryWriter writer;
  if (testing::CompressAnimation(
          {"source0.pgm", "source1.png"}, {50, 50}, &writer, nullptr,
          EncoderConfig::kDefault.quality, EncoderConfig::kDefault.speed,
          EncoderConfig::kDefault.thread_level,
          /*num_downsamplings=*/2) != WP2_STATUS_OK) {
    abort();
  }
  Data data;  // Transfer ownership.
  using std::swap;
  swap(data.bytes, writer.mem_);
  swap(data.size, writer.size_);
  writer.max_size_ = 0;
  return data;
}

std::vector<Data> GenerateAnimatedBitstreamChunks() {
  const Data bitstream = GenerateAnimatedBitstream();
  std::vector<Data> chunks;
  const size_t chunk_size = bitstream.size / 10;
  for (size_t i = 0; i < bitstream.size; i += chunk_size) {
    chunks.emplace_back(Data());
    if (chunks.back().CopyFrom(bitstream.bytes + i,
                               std::min(chunk_size, bitstream.size - i)) !=
        WP2_STATUS_OK) {
      abort();
    }
  }
  return chunks;
}

//------------------------------------------------------------------------------
// Still images
// Note: These are examples using the Decoder class. The Decode() function
//       should be enough for these simple use cases (see decode.h).

// Decode the still image or the first frame contained in a full bitstream.
TEST(DecoderAPITest, StillImage) {
  const Data bitstream = GenerateBitstream();

  ArrayDecoder decoder(bitstream.bytes, bitstream.size);
  if (decoder.ReadFrame()) {
    // Success! 'decoder.GetPixels()' can be accessed here.
  } else {
    // Failure! See 'decoder.GetStatus()'.
  }
}

// Decode the still image or the first frame contained in a full bitstream
// with specific settings.
TEST(DecoderAPITest, StillImageConfig) {
  const Data bitstream = GenerateBitstream();

  DecoderConfig config;
  config.thread_level = 11;  // Decode using 12 threads (including this one).
  ArrayDecoder decoder(bitstream.bytes, bitstream.size, config);
  if (decoder.ReadFrame()) { /* Success! */ }
}

// Decode the still image or the first frame contained in a full bitstream
// into a user-allocated pixel array.
TEST(DecoderAPITest, StillImageToUserMem) {
  const Data bitstream = GenerateBitstream();

  // Begin by extracting the width and height of the image.
  BitstreamFeatures features;
  ASSERT_WP2_OK(features.Read(bitstream.bytes, bitstream.size));

  // Pixels are stored as A,r,g,b,A,r,g,b,A,r,g,b... so there are 'width*4'
  // bytes per row.
  const uint32_t num_bytes_per_row = features.width * 4;
  std::vector<uint8_t> pixels(num_bytes_per_row * features.height);
  ArgbBuffer buffer;
  ASSERT_WP2_OK(buffer.SetExternal(features.width, features.height,
                                   pixels.data(), num_bytes_per_row));

  ArrayDecoder decoder(bitstream.bytes, bitstream.size,
                       DecoderConfig::kDefault, &buffer);
  if (decoder.ReadFrame()) { /* Success! */ }
}

//------------------------------------------------------------------------------
// Metadata

// Decode the still image or the first frame contained in a full bitstream
// and the metadata located at the end of the bitstream (XMP, EXIF).
// Note: This is an example using the Decoder class. The Decode() and/or
//       ExtractMetadata() functions should be enough for this simple use case.
TEST(DecoderAPITest, StillImageWithMetadata) {
  const Data bitstream = GenerateBitstream();

  ArrayDecoder decoder(bitstream.bytes, bitstream.size);
  ASSERT_TRUE(decoder.ReadFrame());          // Decode one frame.
  decoder.SkipNumNextFrames(kMaxNumFrames);  // Ignore everything till XMP/EXIF.
  decoder.ReadFrame();                       // Decode XMP/EXIF, if any.
  if (decoder.GetStatus() == WP2_STATUS_OK) {
    // Success! Access ICC, XMP, EXIF chunks in 'decoder.GetPixels().metadata'.
    ASSERT_FALSE(decoder.GetPixels().metadata.IsEmpty());
  }
}

//------------------------------------------------------------------------------
// Animations

// Decode all frames contained in a full bitstream and display them.
TEST(DecoderAPITest, AnimationToScreen) {
  const Data bitstream = GenerateAnimatedBitstream();

  ArrayDecoder decoder(bitstream.bytes, bitstream.size);
  uint32_t duration_ms;
  while (decoder.ReadFrame(&duration_ms)) {
    // Display 'decoder.GetPixels()' for 'duration_ms'.
  }
  if (decoder.GetStatus() != WP2_STATUS_OK) {
    // An error happened or data is missing.
  }
}

// Decode all frames contained in a full bitstream into an array of frames.
TEST(DecoderAPITest, AnimationToVector) {
  const Data bitstream = GenerateAnimatedBitstream();

  ArrayDecoder decoder(bitstream.bytes, bitstream.size);
  std::vector<ArgbBuffer> frames;
  while (decoder.ReadFrame()) {
    frames.emplace_back(ArgbBuffer());
    ASSERT_WP2_OK(frames.back().CopyFrom(decoder.GetPixels()));
  }
  if (decoder.GetStatus() != WP2_STATUS_OK) { /* Error */ }
}

// Count the number of frames and their features before decoding any pixel.
TEST(DecoderAPITest, AnimationFeatures) {
  const Data bitstream = GenerateAnimatedBitstream();

  ArrayDecoder decoder(bitstream.bytes, bitstream.size);
  decoder.SkipNumNextFrames(kMaxNumFrames);  // Disable pixel decoding.
  decoder.ReadFrame();                       // Decode features.

  if (decoder.GetStatus() == WP2_STATUS_OK) {
    // Number of frames is given by:        decoder.GetNumFrameDecodedFeatures()
    // Frame features can be accessed with: decoder.TryGetFrameDecodedFeatures()

    // Rewind and start decoding pixels.
    decoder.Rewind();
    while (decoder.ReadFrame()) { /* ... */ }
  }
}

//------------------------------------------------------------------------------
// Streamed animations

constexpr size_t kStep = 100;

// Decode all frames contained in a growing bitstream.
TEST(DecoderAPITest, AnimationArray) {
  const Data bitstream = GenerateAnimatedBitstream();
  size_t num_bytes = 0;

  ArrayDecoder decoder;
  while (num_bytes < bitstream.size) {
    num_bytes = std::min(num_bytes + kStep, bitstream.size);  // Simulate stream

    decoder.SetInput(bitstream.bytes, num_bytes);  // Growing input array
    while (decoder.ReadFrame()) { /* Frame is ready */ }
    if (decoder.Failed()) break;
  }
  if (decoder.GetStatus() != WP2_STATUS_OK) { /* Error */ }
}

// Decode all frames contained in a stream available chunk by chunk.
TEST(DecoderAPITest, AnimationStream) {
  std::vector<Data> chunks = GenerateAnimatedBitstreamChunks();

  StreamDecoder decoder;
  for (Data& chunk : chunks) {
    decoder.AppendInput(chunk.bytes, chunk.size);
    // 'chunk' is destroyed now but it is ok, AppendInput() copied it.
    chunk.Clear();

    while (decoder.ReadFrame()) { /* Frame is ready */ }
    if (decoder.Failed()) break;
  }
  if (decoder.GetStatus() != WP2_STATUS_OK) { /* Error */ }
}

// Decode all frames contained in a stream available chunk by chunk and avoid
// copying chunks if possible.
TEST(DecoderAPITest, AnimationStreamPersistent) {
  std::vector<Data> chunks = GenerateAnimatedBitstreamChunks();

  StreamDecoder decoder;
  for (Data& chunk : chunks) {
    // 'chunk' is not yet destroyed, signal that fact to AppendInput() to avoid
    // or shorten a copy if possible.
    decoder.AppendInput(chunk.bytes, chunk.size,
                        /*data_is_persistent=*/true);

    while (decoder.ReadFrame()) { /* Frame is ready */ }
    if (decoder.Failed()) break;

    // AppendInput() expects 'new_chunk' to be persistent till the next call to
    // AppendInput(), so call it again.
    decoder.AppendInput(nullptr, 0);
    // 'chunk' is destroyed now but it is ok, AppendInput() copied the part of
    // it needed to continue the decoding.
    chunk.Clear();
  }
  if (decoder.GetStatus() != WP2_STATUS_OK) { /* Error */ }
}

//------------------------------------------------------------------------------

// TODO(yguyon): Example for decoding an animation 'loop_count' times
// TODO(yguyon): Example for CustomDecoder

//------------------------------------------------------------------------------

}  // namespace
}  // namespace WP2
