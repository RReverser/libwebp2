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

// Incremental decoding test function.

#include "./helpers_incr.h"

#include <cstdlib>

#include "./helpers.h"
#include "src/utils/orientation.h"

namespace WP2 {

//------------------------------------------------------------------------------

namespace {

class ExternalCustomDecoder : public CustomDecoder {
 public:
  ExternalCustomDecoder(const DecoderConfig& config, ArgbBuffer* output_buffer)
      : CustomDecoder(config, output_buffer) {}

  void SetData(const uint8_t* data, size_t size) {
    data_ = data;
    size_ = size;
  }

 private:
  void Discard(size_t num_bytes) override { num_discarded_bytes_ += num_bytes; }

  void Reset() override {
    data_ = nullptr;
    size_ = 0;
    num_discarded_bytes_ = 0;
  }

 protected:
  const uint8_t* data_ = nullptr;
  size_t size_ = 0;
  size_t num_discarded_bytes_ = 0;
};

// Invalidates the buffer between CustomDecoder::Read() calls.
class UnstableCustomDecoder : public ExternalCustomDecoder {
 public:
  UnstableCustomDecoder(const DecoderConfig& config, ArgbBuffer* output_buffer)
      : ExternalCustomDecoder(config, output_buffer) {}

  void MakePreviousDataDisappear() {
    if (disappearing_input_.IsEmpty()) return;
    memset(disappearing_input_.bytes, 0, disappearing_input_.size);
    disappearing_input_.Clear();
  }

 private:
  void Fetch(size_t num_requested_bytes, const uint8_t** available_bytes,
             size_t* num_available_bytes) override {
    (void)num_requested_bytes;
    if (num_discarded_bytes_ <= size_) {
      WP2_ASSERT_STATUS(disappearing_input_.CopyFrom(
          data_ + num_discarded_bytes_, size_ - num_discarded_bytes_));
    } else {
      // All bytes were discarded (or more if chunks with corrupted sizes were
      // skipped).
      disappearing_input_.Clear();
    }
    *available_bytes = disappearing_input_.bytes;
    *num_available_bytes = disappearing_input_.size;
  }

  Data disappearing_input_;
};

// Changes the buffer between DataSource::Fetch() calls. Since it only returns
// at most 'num_requested_bytes', the buffer should change between most of the
// calls to TryGetNext() too.
class SwappingCustomDecoder : public ExternalCustomDecoder {
 public:
  SwappingCustomDecoder(const DecoderConfig& config, ArgbBuffer* output_buffer)
      : ExternalCustomDecoder(config, output_buffer) {}

 private:
  void Fetch(size_t num_requested_bytes, const uint8_t** available_bytes,
             size_t* num_available_bytes) override {
    swap(current_input_, previous_input_);
    if (!previous_input_.IsEmpty()) {
      // Make sure previous input is garbage but is not deleted yet.
      memset(previous_input_.bytes, 0, previous_input_.size);
    }

    if (num_discarded_bytes_ <= size_) {
      // Input from two Fetch() calls ago is freed and replaced.
      WP2_ASSERT_STATUS(current_input_.CopyFrom(
          data_ + num_discarded_bytes_,
          std::min(size_ - num_discarded_bytes_, num_requested_bytes)));
    } else {
      // All bytes were discarded (or more if chunks with corrupted sizes were
      // skipped).
      current_input_.Clear();
    }
    *available_bytes = current_input_.bytes;
    *num_available_bytes = current_input_.size;
  }

  Data current_input_;
  Data previous_input_;
};

void VerifyFrameFeatures(const FrameFeatures* const frame_features,
                         uint32_t frame_index, const ArgbBuffer& output) {
  ASSERT_NE(frame_features, nullptr);
  ASSERT_LE(frame_features->last_dispose_frame_index, frame_index);
  if (frame_index == 0 ||
      frame_features->last_dispose_frame_index == frame_index) {
    // Disposed frame or not an animation.
    ASSERT_EQ(frame_features->window.x, 0u);
    ASSERT_EQ(frame_features->window.y, 0u);
    if (!output.IsEmpty()) {
      ASSERT_EQ(frame_features->window.width, output.width);
      ASSERT_EQ(frame_features->window.height, output.height);
    }
  } else {
    ASSERT_GT(frame_features->window.GetArea(), 0u);
    if (!output.IsEmpty()) {
      ASSERT_LE(frame_features->window.x + frame_features->window.width,
                output.width);
      ASSERT_LE(frame_features->window.y + frame_features->window.height,
                output.height);
    }
  }
}

// Copies pixels 'from' 'rect' 'to' 'rect'.
WP2Status ImportRect(const ArgbBuffer& from, const Rectangle& rect,
                     ArgbBuffer* const to) {
  WP2_CHECK_OK(to != nullptr, WP2_STATUS_NULL_PARAMETER);
  if (rect.GetArea() == 0) return WP2_STATUS_OK;
  ArgbBuffer from_view(from.format), to_view(to->format);
  WP2_CHECK_STATUS(from_view.SetView(from, rect));
  WP2_CHECK_STATUS(to_view.SetView(*to, rect));
  for (uint32_t y = 0; y < rect.height; ++y) {
    WP2_CHECK_STATUS(
        to_view.ImportRow(from_view.format, y, from_view.GetRow8(y)));
  }
  return WP2_STATUS_OK;
}

}  // namespace

//------------------------------------------------------------------------------

namespace testing {

// Like assert() but is not silent if NDEBUG is defined (-c opt).
#define CHECK_OR_DIE(cond)     \
  do {                         \
    if (!(cond)) std::abort(); \
  } while (0)

WP2Status DecodeIncremental(const DecoderConfig& config, DataView input,
                            const IncrementalDecodingTestSetup& setup,
                            std::vector<ArgbBuffer>* const decoded_frames,
                            std::vector<uint32_t>* const decoded_durations_ms) {
  const size_t incr_size_step = std::max((size_t)1, setup.incr_size_step);

  // Setup the chosen Decoder instance.
  ArgbBuffer output;
  ArrayDecoder array_idec(config, &output);
  StreamDecoder stream_idec(config, &output);
  UnstableCustomDecoder unstable_custom_idec(config, &output);
  SwappingCustomDecoder swapping_custom_idec(config, &output);
  Decoder* idec;
  if (setup.decoder_type == DecoderType::kArray) {
    idec = (Decoder*)&array_idec;
  } else if (setup.decoder_type == DecoderType::kStream) {
    idec = (Decoder*)&stream_idec;
  } else if (setup.decoder_type == DecoderType::kUnstableCustom) {
    idec = (Decoder*)&unstable_custom_idec;
  } else {
    idec = (Decoder*)&swapping_custom_idec;
  }

  // State of the decoding.
  uint32_t frame_index = 0;
  Rectangle decoded_area;
  Rectangle previous_decoded_area(0, 0, 0, 0);
  ArgbBuffer incremental_output;
  bool rewinded = false;
  bool skipped_frames = false;
  uint32_t expected_total_num_frames = 0;

  // Start, continue, rewind decoding as long as there is something left to do.
  size_t available_input_size = 0;
  uint32_t current_action_index = 0;
  while (available_input_size < input.size ||
         idec->GetStatus() == WP2_STATUS_NOT_ENOUGH_DATA ||
         current_action_index < setup.actions.size()) {
    WP2_CHECK_OK(!idec->Failed(), idec->GetStatus());
    size_t previous_input_size = available_input_size;
    available_input_size += incr_size_step;
    if (available_input_size > input.size) available_input_size = input.size;

    // Rewind or skip frames as planned.
    while (current_action_index < setup.actions.size() &&
           (setup.actions[current_action_index].bistream_position <=
                available_input_size ||
            previous_input_size >= input.size)) {
      const DecoderAction& action = setup.actions[current_action_index];
      if (action.type == DecoderAction::Type::kRewind ||
          action.type == DecoderAction::Type::kRewindKeepBytes) {
        const bool output_buffer_changed = (action.value > 0);
        if (output_buffer_changed) output.Deallocate();
        idec->Rewind(output_buffer_changed);
        previous_input_size = 0;
        if (action.type != DecoderAction::Type::kRewindKeepBytes) {
          available_input_size = 0;
        }
        frame_index = 0;
        previous_decoded_area = Rectangle(0, 0, 0, 0);
        incremental_output.Deallocate();
        if (decoded_frames != nullptr) decoded_frames->clear();
        if (decoded_durations_ms != nullptr) decoded_durations_ms->clear();
        rewinded = true;
        skipped_frames = false;
      } else if (action.type == DecoderAction::Type::kSkip) {
        const uint32_t num_frames_to_skip = action.value;
        idec->SkipNumNextFrames(num_frames_to_skip);
        if (num_frames_to_skip > 0) {
          frame_index = std::min(frame_index + action.value, kMaxNumFrames);
          previous_decoded_area = Rectangle(0, 0, 0, 0);
          incremental_output.Deallocate();
          skipped_frames = true;
        }
      }
      ++current_action_index;
    }

    // Update the input.
    if (setup.decoder_type == DecoderType::kArray) {
      array_idec.SetInput(input.bytes, available_input_size);
    } else if (setup.decoder_type == DecoderType::kStream) {
      stream_idec.AppendInput(input.bytes + previous_input_size,
                              available_input_size - previous_input_size,
                              /*data_is_persistent=*/true);
    } else if (setup.decoder_type == DecoderType::kUnstableCustom) {
      unstable_custom_idec.SetData(input.bytes, available_input_size);
    } else {
      swapping_custom_idec.SetData(input.bytes, available_input_size);
    }

    // Decode all available bytes until a non-skipped frame is ready.
    uint32_t duration_ms = 0;
    const bool decoded_a_full_frame = idec->ReadFrame(&duration_ms);
    WP2_CHECK_OK(!idec->Failed(), idec->GetStatus());

    if (setup.decoder_type == DecoderType::kUnstableCustom) {
      // Exercise with an invalid buffer until Fetch() is called (during the
      // next ReadFrame() call).
      unstable_custom_idec.MakePreviousDataDisappear();
    }

    // Verify basic values.
    if (idec->TryGetDecodedFeatures() != nullptr &&
        idec->TryGetDecodedFeatures()->is_animation) {
      const uint32_t num_known_frames = idec->GetNumFrameDecodedFeatures();
      if (!skipped_frames || decoded_a_full_frame) {
        CHECK_OR_DIE(idec->GetCurrentFrameIndex() == frame_index);
        if (!rewinded) {
          CHECK_OR_DIE(frame_index <= num_known_frames);
          if (decoded_a_full_frame) {
            CHECK_OR_DIE(frame_index < num_known_frames);
          }
        }
      }
      for (uint32_t i = 0; i < num_known_frames; ++i) {
        VerifyFrameFeatures(idec->TryGetFrameDecodedFeatures(i), i, output);
      }
      if (idec->GetStatus() == WP2_STATUS_OK) {
        for (uint32_t i = 0; i < num_known_frames; ++i) {
          CHECK_OR_DIE(idec->TryGetFrameDecodedFeatures(i)->is_last ==
                       (i == num_known_frames - 1u));
        }
        CHECK_OR_DIE(idec->GetNumFrameDecodedFeatures() ==
                     idec->GetNumAvailableFrames());
        if (expected_total_num_frames > 0u) {
          CHECK_OR_DIE(idec->GetNumAvailableFrames() ==
                       expected_total_num_frames);
        }
      } else {
        CHECK_OR_DIE(idec->GetNumFrameDecodedFeatures() >=
                     idec->GetNumAvailableFrames());
        if (!rewinded && available_input_size < input.size &&
            available_input_size + kANMFHeaderSize >= input.size &&
            !idec->TryGetDecodedFeatures()->has_trailing_data) {
          // Close to the end, no metadata afterwards, expecting one more frame.
          if (expected_total_num_frames == 0u) {
            expected_total_num_frames = idec->GetNumAvailableFrames() + 1u;
          } else {
            CHECK_OR_DIE(expected_total_num_frames ==
                         idec->GetNumAvailableFrames() + 1u);
          }
        }
      }
    } else {
      CHECK_OR_DIE(idec->GetCurrentFrameIndex() == frame_index);
      if (!rewinded) {
        CHECK_OR_DIE(idec->GetNumFrameDecodedFeatures() <= 1u);
      }
    }
    decoded_area = idec->GetDecodedArea();
    const BitstreamFeatures* const features = idec->TryGetDecodedFeatures();
    if (features != nullptr) {
      if (!output.IsEmpty()) {
        CHECK_OR_DIE(features->width == output.width);
        CHECK_OR_DIE(features->height == output.height);
      }
      const Rectangle unoriented =
          RotateRectangle(GetInverseOrientation(features->orientation),
                          features->width, features->height, decoded_area);
      CHECK_OR_DIE(unoriented.x == 0u);
      CHECK_OR_DIE(unoriented.y == 0u);
    }
    CHECK_OR_DIE(decoded_area.width >= previous_decoded_area.width);
    CHECK_OR_DIE(decoded_area.width <= output.width);
    CHECK_OR_DIE(decoded_area.height >= previous_decoded_area.height);
    CHECK_OR_DIE(decoded_area.height <= output.height);

    // Save incremental output to make sure data doesn't change later.
    if (decoded_area.GetArea() > 0) {
      CHECK_OR_DIE(features != nullptr);
      if (incremental_output.IsEmpty()) {
        WP2_CHECK_STATUS(
            incremental_output.Resize(features->width, features->height));
      }

      if (previous_decoded_area.GetArea() > 0) {
        // Check the pixels from last loop before overwriting them.
        float disto[5];
        WP2_CHECK_STATUS(incremental_output.GetDistortion(
            output, previous_decoded_area, PSNR, disto));
        CHECK_OR_DIE(disto[4] == 99.f);
      }
      WP2_CHECK_STATUS(ImportRect(output, decoded_area, &incremental_output));
    }

    if (idec->TryGetDecodedFeatures() != nullptr) {
      if (decoded_a_full_frame) {
        CHECK_OR_DIE(decoded_area.width == output.width);
        CHECK_OR_DIE(decoded_area.height == output.height);
        // Compare incremental data and the current canvas.
        CHECK_OR_DIE(testing::Compare(incremental_output, output,
                                      "incrementally decoded"));
        if (decoded_frames != nullptr) {
          decoded_frames->push_back(ArgbBuffer());
          WP2_CHECK_STATUS(decoded_frames->back().CopyFrom(output));
        }
        if (idec->TryGetDecodedFeatures()->is_animation) {
          if (decoded_durations_ms != nullptr) {
            decoded_durations_ms->push_back(duration_ms);
          }
          CHECK_OR_DIE(
              idec->TryGetFrameDecodedFeatures(frame_index)->duration_ms ==
              duration_ms);
        }
        if (!idec->TryGetFrameDecodedFeatures(frame_index)->is_last) {
          ++frame_index;
        }
        previous_decoded_area = Rectangle(0, 0, 0, 0);
      } else {
        if (!idec->TryGetDecodedFeatures()->has_trailing_data && !rewinded &&
            frame_index == 0 && available_input_size < input.size) {
          CHECK_OR_DIE(idec->GetNumAvailableFrames() == 0u);
        }
        previous_decoded_area = decoded_area;
      }
    }

    if (current_action_index >= setup.actions.size() &&
        previous_input_size >= input.size && !decoded_a_full_frame) {
      // Nothing will happen further so avoid an infinite loop.
      WP2_CHECK_STATUS(idec->GetStatus());
    }
  }
  return WP2_STATUS_OK;
}

#undef CHECK_OR_DIE

}  // namespace testing

//------------------------------------------------------------------------------

}  // namespace WP2
