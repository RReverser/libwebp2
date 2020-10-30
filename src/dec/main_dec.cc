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
// main entry for the decoder
//
// Author: Skal (pascal.massimino@gmail.com)

#include <algorithm>
#include <cstddef>

#include "src/common/color_precision.h"
#include "src/common/progress_watcher.h"
#include "src/dec/filters/intertile_filter.h"
#include "src/dec/wp2_dec_i.h"
#include "src/utils/data_source.h"
#include "src/utils/orientation.h"
#include "src/utils/utils.h"
#include "src/wp2/base.h"
#include "src/wp2/decode.h"

namespace WP2 {

const DecoderConfig DecoderConfig::kDefault;

//------------------------------------------------------------------------------

namespace {

// Returns true if the 'buffer_size' and 'buffer_stride' can hold an image of
// 'width' by 'height' with 'num_bytes_per_pixel'.
bool IsBufferLengthValid(uint32_t width, uint32_t height,
                         uint32_t num_bytes_per_pixel, uint32_t buffer_stride,
                         size_t buffer_size) {
  assert(width > 0 && height > 0);
  const size_t min_stride = width * num_bytes_per_pixel;
  const size_t min_buffer_size =
      (size_t)(height - 1u) * buffer_stride + min_stride;
  return (buffer_stride <= buffer_size && buffer_stride >= min_stride &&
          buffer_size >= min_buffer_size);
}

// Read input data into each Tile::private_input.
WP2Status ReadTiles(const DecoderConfig& config, DataSource* const data_source,
                    TilesLayout* const tiles_layout) {
  for (Tile& tile : tiles_layout->tiles) {
    const size_t size_before = data_source->GetNumReadBytes();
    WP2_CHECK_STATUS(DecodeTileChunkSize(*tiles_layout->gparams, tile.rect,
                                         data_source, &tile.chunk_size));
    tile.chunk_size_is_known = true;
    RegisterChunkSize(config, tile.chunk_size,
                      data_source->GetNumReadBytes() - size_before);
    WP2_CHECK_OK(data_source->TryReadNext(tile.chunk_size, &tile.data_handle),
                 WP2_STATUS_NOT_ENOUGH_DATA);
  }
  for (Tile& tile : tiles_layout->tiles) {
    const uint8_t* const tile_data = tile.data_handle.GetBytes();
    // Maybe a DataSource::TryReadNext() invalidated some previously read data.
    WP2_CHECK_OK(tile_data != nullptr, WP2_STATUS_BITSTREAM_OUT_OF_MEMORY);
    tile.private_input = ExternalDataSource(tile_data, tile.chunk_size);
    tile.input = &tile.private_input;
  }
  return WP2_STATUS_OK;
}

// Decode the data in each Tile::input.
WP2Status DecodeTiles(const DecoderConfig& config,
                      DataSource* const data_source,
                      TilesLayout* const tiles_layout,
                      Vector<TileDecoder>* const workers) {
  WP2_CHECK_STATUS(ReadTiles(config, data_source, tiles_layout));
  WP2Status status = WP2_STATUS_OK;

  // Decode tiles in parallel or sequentially.
  for (TileDecoder& worker : *workers) {
    // The last worker runs in the main thread, after starting others.
    const bool threaded = (&worker != &workers->back());
    status = worker.Start(threaded);
    if (status != WP2_STATUS_OK) break;
  }

  // Loop over all threads and close them properly. No early return.
  for (TileDecoder& worker : *workers) {
    const WP2Status worker_status = worker.End();
    if (status == WP2_STATUS_OK) status = worker_status;
  }
  return status;
}

//------------------------------------------------------------------------------

WP2Status InternalDecode(const DecoderConfig& config,
                         DataSource* const data_source,
                         ArgbBuffer* const rgb_output,
                         YUVPlane* const yuv_output, Metadata* const metadata,
                         CSPTransform* const csp_transform = nullptr) {
  WP2DecDspInit();
  BitstreamFeatures features;
  ProgressWatcher progress(config.progress_hook);
  WP2_CHECK_STATUS(progress.Rewind());

  // Decode the header (width, height, etc.).
  WP2_CHECK_STATUS(DecodeHeader(data_source, &features));
  WP2_CHECK_STATUS(progress.AdvanceBy(kProgressDecHeader));
  WP2_CHECK_STATUS(SetupDecoderInfo(features, config));

  // Preview is only accessible through ExtractPreview().
  WP2_CHECK_STATUS(SkipPreview(data_source, features));
  if (features.has_preview) {
    WP2_CHECK_STATUS(progress.AdvanceBy(kProgressDecPreview));
  }

  // Exactly one output format.
  const bool rgb_output_wanted = (rgb_output != nullptr);
  const bool yuv_output_wanted = (yuv_output != nullptr);
  WP2_CHECK_OK(rgb_output_wanted != yuv_output_wanted,
               WP2_STATUS_INVALID_PARAMETER);

  // Placeholders for missing alpha, unoriented planes, padding etc.
  YUVPlane padded_yuv, non_padded_yuv;
  Plane16 alpha_plane;

  // Reset metadata and decode ICC if any.
  if (metadata != nullptr) metadata->Clear();
  if (rgb_output_wanted) {
    // ArgbBuffer::Resize() clears even metadata so do it before anything else.
    WP2_CHECK_STATUS(
        rgb_output->Resize(features.raw_width, features.raw_height));
  }
  if (metadata != nullptr) {
    WP2_CHECK_STATUS(DecodeICC(data_source, features, &metadata->iccp));
  } else {
    WP2_CHECK_STATUS(SkipICC(data_source, features));
  }
  if (features.has_icc) WP2_CHECK_STATUS(progress.AdvanceBy(kProgressDecICC));

  // Decode the first frame if it is an animation, the whole image otherwise.
  AnimationFrame frame;
  uint32_t frame_index = 0;
  bool is_preframe;
  do {
    WP2_CHECK_STATUS(DecodeANMF(data_source, features, frame_index, &frame));
    is_preframe = (frame.duration_ms == 0);
    // We cannot have more than kMaxNumPreframes preframes.
    WP2_CHECK_OK(!(is_preframe && frame_index >= kMaxNumPreframes),
                 WP2_STATUS_BITSTREAM_ERROR);
    // Preframes are not supported when outputting yuv.
    WP2_CHECK_OK(!(is_preframe && yuv_output_wanted),
                 WP2_STATUS_UNSUPPORTED_FEATURE);
    // The first frame must dispose.
    WP2_CHECK_OK(frame_index > 0 || frame.dispose, WP2_STATUS_BITSTREAM_ERROR);
    if (features.is_animation && frame_index == 0) {
      WP2_CHECK_STATUS(progress.AdvanceBy(kProgressDecANMF));
    }

    // From now on, progress will advance as pixels are available in the buffer.
    // If this is a pre-frame, we assume there's one more frame, but there
    // might actually be even more.
    const float num_pixels = (float)features.width * features.height *
                             (frame_index + 1 + (is_preframe ? 1 : 0));
    const float progress_range_for_pixels =
        1.f - progress.GetProgress() - kProgressDecEnd -
        (features.has_trailing_data ? kProgressDecMetadata : 0.f);
    progress.SetAdvancementScale(progress_range_for_pixels / num_pixels);
    WP2_CHECK_STATUS(progress.AdvanceBy(num_pixels - frame.window.GetArea()));

    TilesLayout tiles_layout;

    // Global parameters for this frame (has an alpha component, etc.)
    GlobalParams gparams;
    tiles_layout.gparams = &gparams;
    WP2_CHECK_STATUS(DecodeGLBL(data_source, config, features, &gparams));
    if (csp_transform != nullptr) *csp_transform = gparams.transf_;
    FillDecoderInfo(gparams, config);
    const uint32_t padded_frame_width = Pad(frame.window.width, kPredWidth);
    const uint32_t padded_frame_height = Pad(frame.window.height, kPredWidth);

    // Setup a YUV buffer if needed.
    if (yuv_output_wanted) {
      // The right and bottom tiles might need padding for easier block
      // decoding. Depending on the 'frame.window', dimensions and orientation
      // of the image, a temporary 'padded_yuv' buffer may be allocated.
      const uint32_t min_width =
          std::max(frame.window.x + padded_frame_width, features.raw_width);
      const uint32_t min_height =
          std::max(frame.window.y + padded_frame_height, features.raw_height);
      if (!yuv_output->IsEmpty() && (yuv_output->GetWidth() >= min_width &&
                                     yuv_output->GetHeight() >= min_height)) {
        WP2_CHECK_STATUS(
            padded_yuv.SetView(*yuv_output, {0, 0, min_width, min_height}));
        if (!gparams.has_alpha_) {
          // No alpha plane to decode.
          padded_yuv.A.Clear();
        } else if (padded_yuv.A.IsEmpty()) {
          // Create an alpha plane that will be decoded but not output.
          // YUVPlane should not contain a mix of owned memory and views.
          WP2_CHECK_STATUS(alpha_plane.Resize(min_width, min_height));
          WP2_CHECK_STATUS(
              padded_yuv.A.SetView(alpha_plane, {0, 0, min_width, min_height}));
        }
      } else {
        WP2_CHECK_STATUS(padded_yuv.Resize(min_width, min_height, /*pad=*/1,
                                           gparams.has_alpha_));
      }

      WP2_CHECK_STATUS(non_padded_yuv.SetView(
          padded_yuv, {0, 0, features.raw_width, features.raw_height}));
    }

    const bool rgb_output_needed =
        rgb_output_wanted || (gparams.type_ == GlobalParams::GP_LOSSLESS ||
                              gparams.type_ == GlobalParams::GP_BOTH);
    const bool yuv_output_needed =
        yuv_output_wanted || (gparams.type_ == GlobalParams::GP_LOSSY ||
                              gparams.type_ == GlobalParams::GP_BOTH);

    // Create views or buffers of the dimensions of the 'frame.window'.
    ArgbBuffer rgb_frame(rgb_output_wanted ? rgb_output->format : WP2_Argb_32);

    if (rgb_output_needed) {
      if (rgb_output_wanted && !frame.blend) {
        WP2_CHECK_STATUS(rgb_frame.SetView(*rgb_output, frame.window));
      } else {
        WP2_CHECK_STATUS(
            rgb_frame.Resize(frame.window.width, frame.window.height));
      }
    }

    YUVPlane yuv_frame;
    if (yuv_output_needed) {
      if (yuv_output_wanted) {
        WP2_CHECK_STATUS(yuv_frame.SetView(
            padded_yuv, {frame.window.x, frame.window.y,
                         padded_frame_width, padded_frame_height}));
      } else {
        WP2_CHECK_STATUS(
            yuv_frame.Resize(frame.window.width, frame.window.height,
                             /*pad=*/kPredWidth, gparams.has_alpha_));
      }
    }

    WP2_CHECK_STATUS(GetTilesLayout(frame.window.width, frame.window.height,
                                    features.tile_width, features.tile_height,
                                    &rgb_frame, &yuv_frame, &tiles_layout));
    // All tiles should be immediately decodable.
    tiles_layout.num_assignable_tiles = (uint32_t)tiles_layout.tiles.size();

    {  // Actual tile decoding.
      // TODO(maryla): support multiple frames (add the bpps?). This will
      // only contain information for the last frame.
      if (config.info != nullptr && config.info->bits_per_pixel != nullptr) {
        WP2_CHECK_STATUS(config.info->bits_per_pixel->Resize(
            frame.window.width, frame.window.height));
      }
      Vector<TileDecoder> workers;
      WP2_CHECK_STATUS(
          SetupWorkers(features, config, &progress, &tiles_layout, &workers));
      WP2_CHECK_STATUS(
          DecodeTiles(config, data_source, &tiles_layout, &workers));
    }

    // In case of the YUV decoding stepping on the background borders because
    // of padding, FillBorders() is done afterwards.
    if (frame.dispose) {
      if (rgb_output_wanted) {
        WP2_CHECK_STATUS(FillBorders(features, frame, rgb_output));
      } else {
        WP2_CHECK_STATUS(
            FillBorders(features, frame, gparams.transf_, &non_padded_yuv));
      }
    }

    // Filter the edges of the lossy adjacent tiles.
    if (!yuv_frame.IsEmpty()) {
      IntertileFilter intertile_filter;
      intertile_filter.Init(config, features, gparams, tiles_layout,
                            &yuv_frame);
      intertile_filter.Deblock(frame.window.height);
      assert(intertile_filter.GetNumFilteredRows() == frame.window.height);
    }

    // Per-tile YUV <-> RGB conversion, if any.
    // TODO(maryla): in case there are preframes, conversion could be done
    // after we finished decoding the final frame if the yuv padding does not
    // stomp over pixels.
    for (Tile& tile : tiles_layout.tiles) {
      if (rgb_output_wanted) {
        if (tile.output_is_yuv) {
          assert(!tile.rgb_output.IsEmpty());
          assert(!tile.yuv_output.IsEmpty());
          YUVPlane yuv;  // Non-padded.
          WP2_CHECK_STATUS(yuv.SetView(
              tile.yuv_output, {0, 0, tile.rect.width, tile.rect.height}));
          WP2_CHECK_STATUS(yuv.Export(
              gparams.transf_, /*resize_if_needed=*/false, &tile.rgb_output));
        }
      } else {  // yuv_output_wanted
        if (!tile.output_is_yuv) {
          assert(!tile.rgb_output.IsEmpty());
          assert(!tile.yuv_output.IsEmpty());
          WP2_CHECK_STATUS(tile.yuv_output.Import(
              tile.rgb_output, gparams.has_alpha_, gparams.transf_,
              /*resize_if_needed=*/false, /*pad=*/kPredWidth));
        }
      }
    }
    if (yuv_output_wanted && !gparams.has_alpha_ && !yuv_output->A.IsEmpty()) {
      yuv_output->A.Fill(kAlphaMax);
    }
    ++frame_index;

    if (frame.blend) {
      if (rgb_output_wanted) {
        if (frame.dispose) {
          rgb_output->Fill(ToArgb32b(features.background_color));
        }
        ArgbBuffer output_view;
        WP2_CHECK_STATUS(output_view.SetView(*rgb_output, frame.window));
        WP2_CHECK_STATUS(output_view.CompositeUnder(rgb_frame));
      } else {
        // TODO(maryla): handle blending for yuv output.
        return WP2_STATUS_UNSUPPORTED_FEATURE;
      }
    }
  } while (is_preframe);

  // Skip all remaining frames (DecodeTiles() is not called).
  while (!frame.is_last) {
    WP2_CHECK_STATUS(DecodeANMF(data_source, features, frame_index, &frame));
    WP2_CHECK_STATUS(SkipTiles(data_source, features, frame.window));
    ++frame_index;
  }
  progress.SetAdvancementScale(1.);

  // Decode or skip trailing metadata if any (EXIF, XMP).
  if (metadata != nullptr) {
    WP2_CHECK_STATUS(
        DecodeMetadata(data_source, features, &metadata->exif, &metadata->xmp));
  } else {
    WP2_CHECK_STATUS(SkipMetadata(data_source, features));
  }
  if (features.has_trailing_data) {
    WP2_CHECK_STATUS(progress.AdvanceBy(kProgressDecMetadata));
  }

  // TODO(yguyon): rotate directly during decoding
  if (rgb_output_wanted) {
    WP2_CHECK_STATUS(RotateBuffer(features.orientation, rgb_output));
  } else {  // yuv_output_wanted
    // If 'padded_yuv' is a view, it means 'yuv_output' already contains the
    // right output and must only be rotated. Otherwise 'padded_yuv' owns the
    // pixels and the view 'non_padded_yuv' should be rotated to 'yuv_output'.
    if (padded_yuv.IsView()) {
      WP2_CHECK_STATUS(RotateBuffer(features.orientation, yuv_output));
    } else {
      WP2_CHECK_STATUS(
          RotateBuffer(features.orientation, non_padded_yuv, yuv_output));
    }
  }
  WP2_CHECK_STATUS(progress.Finish());
  return WP2_STATUS_OK;
}

}  // namespace

//------------------------------------------------------------------------------

WP2Status Decode(const uint8_t* data, size_t data_size, ArgbBuffer* output,
                 const DecoderConfig& config) {
  WP2_CHECK_OK(output != nullptr, WP2_STATUS_NULL_PARAMETER);
  const WP2SampleFormat format_default =
      (WP2FormatBpc(output->format) == 1) ? WP2_Argb_32 : WP2_Argb_38;
  ArgbBuffer tmp_output(format_default);
  ArgbBuffer* real_output = nullptr;

  if ((output->format != format_default) || output->IsSlow()) {
    // work on tmp buffer
    real_output = output;
    output = &tmp_output;  // tmp working buffer
  }
  ExternalDataSource data_source(data, data_size);
  WP2Status status = InternalDecode(config, &data_source, output,
                                    /*yuv_output=*/nullptr, &output->metadata);
  if (status != WP2_STATUS_OK) {
    output->Deallocate();
    // Because it's not incremental decoding, missing data means an error.
    if (status == WP2_STATUS_NOT_ENOUGH_DATA) return WP2_STATUS_BITSTREAM_ERROR;
    return status;
  }
  if (real_output != nullptr) {
    WP2_CHECK_STATUS(real_output->ConvertFrom(*output));
    swap(real_output->metadata, output->metadata);
  }
  return WP2_STATUS_OK;
}

WP2Status Decode(const std::string& data, ArgbBuffer* output,
                 const DecoderConfig& config) {
  return Decode((const uint8_t*)data.data(), data.size(), output, config);
}

WP2Status DecodeArgb_32(const uint8_t* input_data, size_t input_data_size,
                        uint8_t* output_buffer, uint32_t output_stride,
                        size_t output_buffer_size,
                        const DecoderConfig& config) {
  if (config.progress_hook != nullptr) {
    WP2_CHECK_OK(config.progress_hook->OnUpdate(0.f), WP2_STATUS_USER_ABORT);
    // OnUpdate() will be called twice with 'progress' at 0; not a big deal.
  }
  ArgbBuffer output;
  BitstreamFeatures features;
  ExternalDataSource data_source(input_data, input_data_size);
  const uint32_t bpp = WP2FormatBpp(output.format);
  WP2_CHECK_STATUS(DecodeHeader(&data_source, &features));
  WP2_CHECK_OK(IsBufferLengthValid(features.width, features.height,
                                   bpp, output_stride, output_buffer_size),
               WP2_STATUS_BAD_DIMENSION);
  WP2_CHECK_STATUS(output.SetExternal(features.width, features.height,
                                      output_buffer, output_stride));
  return Decode(input_data, input_data_size, &output, config);
}

WP2Status Decode(const uint8_t* input_data, size_t input_data_size,
                 const int16_t rgb_to_ccsp_matrix[9],
                 uint32_t rgb_to_ccsp_shift, int16_t* c0_buffer,
                 uint32_t c0_step, size_t c0_buffer_size, int16_t* c1_buffer,
                 uint32_t c1_step, size_t c1_buffer_size, int16_t* c2_buffer,
                 uint32_t c2_step, size_t c2_buffer_size, int16_t* a_buffer,
                 uint32_t a_step, size_t a_buffer_size, Metadata* metadata,
                 const DecoderConfig& config) {
  WP2_CHECK_OK(input_data != nullptr, WP2_STATUS_NULL_PARAMETER);
  WP2_CHECK_OK(rgb_to_ccsp_matrix != nullptr, WP2_STATUS_NULL_PARAMETER);
  const CSPMtx rgb_to_ccsp(rgb_to_ccsp_matrix, rgb_to_ccsp_shift);

  if (config.progress_hook != nullptr) {
    WP2_CHECK_OK(config.progress_hook->OnUpdate(0.f), WP2_STATUS_USER_ABORT);
  }
  WP2_CHECK_OK(rgb_to_ccsp.shift <= 16, WP2_STATUS_INVALID_PARAMETER);

  YUVPlane output;
  BitstreamFeatures features;
  ExternalDataSource data_source(input_data, input_data_size);
  WP2_CHECK_STATUS(DecodeHeader(&data_source, &features));

  WP2_CHECK_OK(c0_buffer != nullptr, WP2_STATUS_NULL_PARAMETER);
  WP2_CHECK_OK(c1_buffer != nullptr, WP2_STATUS_NULL_PARAMETER);
  WP2_CHECK_OK(c2_buffer != nullptr, WP2_STATUS_NULL_PARAMETER);


  // Try to fit the worst case into the output buffer to avoid temp buffers.
  uint32_t view_width = std::max(Pad(features.width, kPredWidth),
                                 Pad(features.raw_width, kPredWidth));
  uint32_t view_height = std::max(Pad(features.height, kPredWidth),
                                  Pad(features.raw_height, kPredWidth));

  const uint32_t bpp = sizeof(c0_buffer[0]);  // Number of bytes per pixel.
  if (!IsBufferLengthValid(view_width, view_height, bpp, c0_step * bpp,
                           c0_buffer_size) ||
      !IsBufferLengthValid(view_width, view_height, bpp, c1_step * bpp,
                           c1_buffer_size) ||
      !IsBufferLengthValid(view_width, view_height, bpp, c2_step * bpp,
                           c2_buffer_size) ||
      (a_buffer != nullptr &&
       !IsBufferLengthValid(view_width, view_height, bpp,
                            a_step * bpp, a_buffer_size))) {
    // Not possible so settle for the minimum valid output buffer size.
    view_width = features.width;
    view_height = features.height;
    WP2_CHECK_OK(IsBufferLengthValid(view_width, view_height,
                                     bpp, c0_step * bpp, c0_buffer_size),
                 WP2_STATUS_BAD_DIMENSION);
    WP2_CHECK_OK(IsBufferLengthValid(view_width, view_height,
                                     bpp, c1_step * bpp, c1_buffer_size),
                 WP2_STATUS_BAD_DIMENSION);
    WP2_CHECK_OK(IsBufferLengthValid(view_width, view_height,
                                     bpp, c2_step * bpp, c2_buffer_size),
                 WP2_STATUS_BAD_DIMENSION);
    if (a_buffer != nullptr) {
      WP2_CHECK_OK(IsBufferLengthValid(view_width, view_height,
                                       bpp, a_step * bpp, a_buffer_size),
                   WP2_STATUS_BAD_DIMENSION);
    }
  }

  WP2_CHECK_STATUS(
      output.Y.SetView(c0_buffer, view_width, view_height, c0_step));
  WP2_CHECK_STATUS(
      output.U.SetView(c1_buffer, view_width, view_height, c1_step));
  WP2_CHECK_STATUS(
      output.V.SetView(c2_buffer, view_width, view_height, c2_step));
  if (a_buffer != nullptr) {
    WP2_CHECK_STATUS(
        output.A.SetView(a_buffer, view_width, view_height, a_step));
  }

  data_source = ExternalDataSource(input_data, input_data_size);
  CSPTransform csp_transform;
  WP2_CHECK_STATUS(InternalDecode(config, &data_source, /*rgb_output=*/nullptr,
                                  &output, metadata, &csp_transform));

  WP2_CHECK_STATUS(csp_transform.YUVToCustom(
      features.width, features.height, output.Y.Row(0), output.Y.Step(),
      output.U.Row(0), output.U.Step(), output.V.Row(0), output.V.Step(),
      rgb_to_ccsp, output.Y.Row(0), output.Y.Step(), output.U.Row(0),
      output.U.Step(), output.V.Row(0), output.V.Step()));
  return WP2_STATUS_OK;
}

//------------------------------------------------------------------------------

WP2Status ExtractPreview(const uint8_t* const data, size_t data_size,
                         ArgbBuffer* const output_buffer) {
  WP2_CHECK_OK(output_buffer != nullptr, WP2_STATUS_NULL_PARAMETER);

  ExternalDataSource data_source(data, data_size);
  BitstreamFeatures features;

  WP2_CHECK_STATUS(DecodeHeader(&data_source, &features));
  if (features.has_preview) {
    // TODO(yguyon): rotate directly during decoding
    if (output_buffer->IsEmpty()) {
      WP2_CHECK_STATUS(output_buffer->Resize(features.width, features.height));
    }

    if (output_buffer->format != WP2_Argb_32 || output_buffer->IsSlow() ||
        features.orientation == Orientation::k90Clockwise ||
        features.orientation == Orientation::k270Clockwise) {
      // Create another buffer to leave intact the layout of 'output_buffer'
      // that would be modified by a rotation of 90 or 270 degrees.
      // This is also needed if output is not in Argb32 format or is slow.
      ArgbBuffer raw_output_buffer(WP2_Argb_32);
      WP2_CHECK_STATUS(raw_output_buffer.Resize(
          RotateWidth(GetInverseOrientation(features.orientation),
                      output_buffer->width, output_buffer->height),
          RotateHeight(GetInverseOrientation(features.orientation),
                       output_buffer->width, output_buffer->height)));

      WP2_CHECK_STATUS(
          DecodePreview(&data_source, features, &raw_output_buffer));

      WP2_CHECK_STATUS(RotateSubBuffer(
          features.orientation, raw_output_buffer,
          {0, 0, raw_output_buffer.width, raw_output_buffer.height},
          output_buffer));
    } else {
      // The layout of 'output_buffer' will not be modified, only its content.
      WP2_CHECK_STATUS(DecodePreview(&data_source, features, output_buffer));
      WP2_CHECK_STATUS(RotateBuffer(features.orientation, output_buffer));
    }
  } else {
    if (output_buffer->IsEmpty()) {
      WP2_CHECK_STATUS(output_buffer->Resize(features.width, features.height));
    }
    std::array<uint8_t, 4> bytes;
    ToUInt8(ToArgb32b(features.preview_color), bytes.data());
    if (output_buffer->format != WP2_Argb_32) {
      const std::array<uint8_t, 4> bytes_argb32b = bytes;
      WP2ArgbConverterInit();
      WP2ArgbConvertTo[output_buffer->format](bytes_argb32b.data(), 1,
                                              bytes.data());
    }
    output_buffer->Fill({0, 0, output_buffer->width, output_buffer->height},
                        bytes.data());
  }
  return WP2_STATUS_OK;
}

//------------------------------------------------------------------------------

WP2Status ExtractMetadata(const uint8_t* data, size_t data_size,
                          Data* const exif, Data* const iccp, Data* const xmp) {
  if (exif != nullptr) exif->Clear();
  if (iccp != nullptr) iccp->Clear();
  if (xmp != nullptr) xmp->Clear();
  ExternalDataSource data_source(data, data_size);
  BitstreamFeatures features;

  WP2_CHECK_STATUS(DecodeHeader(&data_source, &features));
  WP2_CHECK_STATUS(SkipPreview(&data_source, features));
  WP2_CHECK_STATUS(DecodeICC(&data_source, features, iccp));

  if (features.has_trailing_data && (exif != nullptr || xmp != nullptr)) {
    // Skip all frames because there is trailing data.
    AnimationFrame frame;
    uint32_t frame_index = 0;
    do {
      WP2_CHECK_STATUS(DecodeANMF(&data_source, features, frame_index, &frame));
      WP2_CHECK_STATUS(SkipTiles(&data_source, features, frame.window));
      ++frame_index;
    } while (!frame.is_last);

    WP2_CHECK_STATUS(DecodeMetadata(&data_source, features, exif, xmp));
  }
  return WP2_STATUS_OK;
}

//------------------------------------------------------------------------------

}  // namespace WP2
