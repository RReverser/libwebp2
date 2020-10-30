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
// main entry for the encoder
//
// Author: Skal (pascal.massimino@gmail.com)

#include <cassert>
#include <cstring>

#include "src/common/color_precision.h"
#include "src/common/constants.h"
#include "src/dec/tile_dec.h"
#include "src/enc/partitioner.h"
#include "src/enc/preview/preview_enc.h"
#include "src/enc/wp2_enc_i.h"
#include "src/utils/ans.h"
#include "src/utils/ans_utils.h"
#include "src/utils/utils.h"
#include "src/wp2/base.h"
#include "src/wp2/encode.h"
#include "src/wp2/format_constants.h"

namespace WP2 {

const EncoderConfig EncoderConfig::kDefault;

//------------------------------------------------------------------------------

WP2Status TileEncoder::AssignNextTile() {
  tile_ = nullptr;

  // Get first unassigned tile, and increment the counter if valid.
  WP2_CHECK_STATUS(tiles_layout_->assignment_lock.Acquire());
  const uint32_t tile_index = tiles_layout_->first_unassigned_tile_index;
  if (tile_index < tiles_layout_->tiles.size()) {
    tile_ = &tiles_layout_->tiles[tile_index];
    ++tiles_layout_->first_unassigned_tile_index;
  }
  tiles_layout_->assignment_lock.Release();
  return WP2_STATUS_OK;
}

WP2Status TileEncoder::Execute() {
  while (tile_ != nullptr) {
    ANSEnc& enc = tile_->enc;
    assert(enc.Buffer() == nullptr);
    if (tiles_layout_->gparams->type_ == GlobalParams::GP_BOTH) {
      ANSDebugPrefix prefix(&enc, "GlobalHeader");
      enc.PutBool(use_lossless_, "use_lossless");
    }
    if (use_lossless_) {
      WP2_CHECK_STATUS(LosslessEncode(&enc));
    } else {
      const Rectangle padded_tile_rect = {tile_->rect.x, tile_->rect.y,
                                          Pad(tile_->rect.width, kPredWidth),
                                          Pad(tile_->rect.height, kPredWidth)};
      VectorNoCtor<Block> partition;
      WP2_CHECK_STATUS(AddForcedBlocks(*config_, padded_tile_rect, &partition));
      WP2_CHECK_STATUS(LossyEncode(partition, &enc));
    }

    WP2_CHECK_STATUS(enc.Assemble());

    WP2_CHECK_STATUS(AssignNextTile());
  }
  return WP2_STATUS_OK;
}

//------------------------------------------------------------------------------

size_t WriteVarInt(uint64_t size, uint64_t size_upper,
                   uint8_t dst[kMaxVarIntLength]) {
  uint64_t n = 0;
  if (size >= kUpperVarInt) return 0;
  assert(size < size_upper);
  // TODO(yguyon): No need to write anything when 'size_upper' is 1
  //               or remove possibility of 'size' being 0 (pessimization)
  while (size_upper >= 256 && size >= 128) {
    dst[n++] = (size & 0x7f) | 0x80;
    size >>= 7;
    size_upper >>= 7;
  }
  dst[n++] = size;
  return n;
}

namespace {

WP2Status WriteChunk(DataView data, Writer* const output) {
  if (data.size == 0) return WP2_STATUS_OK;
  uint8_t buf[kMaxVarIntLength];
  const size_t buf_size = WriteVarInt(data.size - 1, kChunkSizeMax, buf);
  WP2_CHECK_OK(output->Append(buf, buf_size), WP2_STATUS_BAD_WRITE);
  WP2_CHECK_OK(output->Append(data.bytes, data.size), WP2_STATUS_BAD_WRITE);
  return WP2_STATUS_OK;
}

bool WriteTag(const uint32_t tag, Writer* const output) {
  const uint8_t data[3] = {(uint8_t)(tag >> 0), (uint8_t)(tag >> 8),
                           (uint8_t)(tag >> 16)};
  return output->Append(data, sizeof(data));
}

//------------------------------------------------------------------------------

// Fills the 'tiles_layout'. Dimensions are in pixels, not padded.
WP2Status GetEncTilesLayout(uint32_t frame_width, uint32_t frame_height,
                            uint32_t tile_width, uint32_t tile_height,
                            const ArgbBuffer& frame_rgb,
                            const YUVPlane& frame_yuv,
                            const GlobalParams& gparams,
                            EncTilesLayout* const tiles_layout) {
  tiles_layout->tile_width = tile_width;
  tiles_layout->tile_height = tile_height;
  const uint32_t num_tiles =
      GetNumTiles(frame_width, frame_height, tile_width, tile_height,
                  &tiles_layout->num_tiles_x, &tiles_layout->num_tiles_y);
  WP2_CHECK_ALLOC_OK(tiles_layout->tiles.resize(num_tiles));

  uint32_t tile_y = 0, tile_x = 0;
  for (EncTile& tile : tiles_layout->tiles) {
    tile.rect = {tile_x * tile_width, tile_y * tile_height,
                 std::min(frame_width - tile_x * tile_width, tile_width),
                 std::min(frame_height - tile_y * tile_height, tile_height)};

    // Set views on RGB and YUV input buffers, one or both may be used.
    assert(!frame_rgb.IsEmpty() || !frame_yuv.IsEmpty());
    if (!frame_rgb.IsEmpty()) {
      assert(frame_rgb.width == frame_width &&
             frame_rgb.height == frame_height);
      WP2_CHECK_STATUS(tile.rgb_input.SetFormat(frame_rgb.format));
      WP2_CHECK_STATUS(tile.rgb_input.SetView(frame_rgb, tile.rect));
    }
    if (!frame_yuv.IsEmpty()) {
      assert(frame_yuv.GetWidth() == Pad(frame_width, kPredWidth) &&
             frame_yuv.GetHeight() == Pad(frame_height, kPredWidth));
      WP2_CHECK_STATUS(tile.yuv_input.SetView(
          frame_yuv,
          {tile.rect.x, tile.rect.y,
           std::min(frame_yuv.GetWidth() - tile.rect.x, tile_width),
           std::min(frame_yuv.GetHeight() - tile.rect.y, tile_height)}));
    }

    if (++tile_x == tiles_layout->num_tiles_x) {
      ++tile_y;
      tile_x = 0;
    }
  }
  tiles_layout->first_unassigned_tile_index = 0;
  tiles_layout->gparams = &gparams;
  return WP2_STATUS_OK;
}

WP2Status SetupWorkers(uint32_t width, uint32_t height,
                       const ArgbBuffer& rgb_buffer, const YUVPlane& yuv_buffer,
                       const EncoderConfig& config, const GlobalParams& gparams,
                       EncTilesLayout* const tiles_layout,
                       Vector<TileEncoder>* const workers) {
  // The number of workers is limited by the number of threads and tiles.
  const size_t num_workers = (uint32_t)std::min(
      (size_t)config.thread_level + 1, tiles_layout->tiles.size());
  WP2_CHECK_ALLOC_OK(workers->resize(num_workers));

  for (TileEncoder& worker : *workers) {
    worker.status_ = WP2_STATUS_OK;
    worker.config_ = &config;
    worker.use_lossless_ =
        (!config.use_neural_compression && config.quality > kMaxLossyQuality);
    worker.tiles_layout_ = tiles_layout;
    WP2_CHECK_STATUS(worker.AssignNextTile());
  }
  return WP2_STATUS_OK;
}

WP2Status CodeTiles(EncTilesLayout* const tiles_layout,
                    Vector<TileEncoder>* const workers, Writer* const output) {
  const bool do_mt = (workers->size() > 1);
  WP2Status status = WP2_STATUS_OK;
  for (TileEncoder& tile : *workers) {
    status = tile.Start(do_mt);
    if (status != WP2_STATUS_OK) break;
  }
  for (TileEncoder& tile : *workers) {
    const WP2Status thread_status = tile.End();
    if (thread_status != WP2_STATUS_OK) {
      status = thread_status;
      continue;
    }
    if (status == WP2_STATUS_OK) status = tile.status_;
  }
  WP2_CHECK_STATUS(status);
  for (EncTile& tile : tiles_layout->tiles) {
    const uint32_t max_num_bytes =
        GetTileMaxNumBytes(*tiles_layout->gparams, tile.rect);
    {
      uint8_t buffer[kMaxVarIntLength];
      const size_t size =
          WriteVarInt(std::min(tile.enc.BufferSize(), (size_t)max_num_bytes),
                      max_num_bytes + 1, buffer);
      WP2_CHECK_ALLOC_OK(output->Append(buffer, size));
    }

    if (tile.enc.BufferSize() >= max_num_bytes) {
      // Discard ANS encoding because it takes at least as many bytes as raw
      // pixels.
      WP2_CHECK_STATUS(BypassTileEnc(*tiles_layout->gparams, &tile, output));
    } else {
      // TODO(yguyon): Check if 'tile.enc.BufferSize()' can even be 0.
      WP2_CHECK_ALLOC_OK(
          output->Append(tile.enc.Buffer(), tile.enc.BufferSize()));
    }
    tile.enc.WipeOut();
  }
  return status;
}

//------------------------------------------------------------------------------

// Get background flag from 'background_color'.
uint8_t BackgroundColorToBackgroundBits(Argb38b background_color) {
  if (background_color.a == 0x00u) return kBackgroundTransparent;
  if (background_color.a == 0xFFu) {
    if (background_color.r == 0x03FFu && background_color.g == 0x03FFu &&
        background_color.b == 0x03FFu) {
      return kBackgroundWhite;
    }
    if (background_color.r == 0x0000u && background_color.g == 0x0000u &&
        background_color.b == 0x0000u) {
      return kBackgroundBlack;
    }
  }
  return kBackgroundCustom;
}

bool HasPreview(const EncoderConfig& config) {
  return (config.create_preview || config.preview_size > 0);
}

}  // anonymous namespace

//------------------------------------------------------------------------------

// Encode header to 'output'.
WP2Status EncodeHeader(const EncoderConfig& config,
                       uint32_t width, uint32_t height,
                       bool has_alpha, bool is_anim, uint32_t loop_count,
                       Argb38b background_color, RGB12b preview_color,
                       bool has_icc, bool has_trailing_data,
                       Writer* const output) {
  WP2_CHECK_OK(CheckPremultiplied(background_color),
               WP2_STATUS_INVALID_PARAMETER);
  const bool use_bt2020 =
      (config.transfer_function == WP2_TF_ITU_R_BT2020_10BIT);
  const uint8_t background_bits =
      BackgroundColorToBackgroundBits(background_color);

  {
    uint8_t header[kHeaderMaxSize];
    HeaderEnc henc(header, sizeof(header));

    // signature
    henc.PutBits(kSignature, 24, "signature");
    // infos
    henc.PutBits(width - 1, kImageDimNumBits, "width_m1");
    henc.PutBits(height - 1, kImageDimNumBits, "height_m1");
    henc.PutBits((uint32_t)config.decoding_orientation, 2, "orientation");
    henc.PutBits(has_alpha ? 1 : 0, 1, "has_alpha");
    henc.PutBits(is_anim ? 1 : 0, 1, "is_animation");

    henc.PutBits(ToUInt32(preview_color), 12, "preview_color");
    henc.PutBits(GetQualityHint(config.quality), kQualityHintNumBits,
                 "quality_hint");

    henc.PutBits(HasPreview(config) ? 1 : 0, 1, "has_preview");
    henc.PutBits(has_icc ? 1 : 0, 1, "has_icc");
    henc.PutBits(has_trailing_data ? 1 : 0, 1, "has_trailing_data");
    henc.PutBits(use_bt2020 ? 1 : 0, 1, "default_transfer_function");

    static_assert(TileShape::TILE_SHAPE_AUTO <= (1 << kTileShapeBits),
                  "invalid TILE_SHAPE_AUTO value");
    henc.PutBits(FinalTileShape(config), kTileShapeBits, "tile_shape");
    henc.PutBits(0x2u, 2, "filler");  // reserved

    if (is_anim) {
      henc.PutBits(loop_count, kLoopCountNumBits, "loop_count");
      henc.PutBits(background_bits, kBackgroundNumBits, "background");
      if (background_bits == kBackgroundCustom) {
        henc.PutBits(background_color.a, 8, "background");
        henc.PutBits(background_color.r, 10, "background");
        henc.PutBits(background_color.g, 10, "background");
        henc.PutBits(background_color.b, 10, "background");
        henc.PutBits(0x0, 2, "filler");
      }
    }
    if (!use_bt2020) {
      henc.PutBits((uint32_t)config.transfer_function - 1, 4,
                   "transfer_function");
      henc.PutBits(0x0, 4, "filler");
    }
    henc.Align();
    WP2_CHECK_OK(henc.Ok(), WP2_STATUS_BAD_WRITE);

    assert(henc.Used() >= kHeaderMinSize && henc.Used() <= kHeaderMaxSize);
    WP2_CHECK_OK(output->Append(header, henc.Used()), WP2_STATUS_BAD_WRITE);
  }
  return WP2_STATUS_OK;
}

WP2Status EncodePreview(const ArgbBuffer& buffer, const EncoderConfig& config,
                        Writer* const output) {
  // TODO(skal): add size limit on preview chunk length?
  if (config.create_preview) {
    MemoryWriter preview;
    WP2_CHECK_STATUS(EncodePreview(
        buffer, PreviewConfig(config.quality, config.speed), &preview));
    WP2_CHECK_STATUS(WriteChunk({preview.mem_, preview.size_}, output));
  } else if (config.preview_size > 0) {
    WP2_CHECK_OK(config.preview != nullptr, WP2_STATUS_INVALID_PARAMETER);
    WP2_CHECK_STATUS(
        WriteChunk({config.preview, (size_t)config.preview_size}, output));
  }
  return WP2_STATUS_OK;
}

WP2Status EncodeICC(DataView iccp, Writer* const output) {
  if (iccp.size > 0) {
    WP2_CHECK_STATUS(WriteChunk(iccp, output));
  }
  return WP2_STATUS_OK;
}

WP2Status EncodeGLBL(const EncoderConfig& config, const GlobalParams& gparams,
                     uint32_t quality_hint, bool image_has_alpha,
                     Writer* const output) {
  ANSEnc enc;
  WP2_CHECK_STATUS(gparams.Write(quality_hint, image_has_alpha, &enc));
  WP2_CHECK_STATUS(enc.Assemble());
  WP2_CHECK_STATUS(WriteChunk({enc.Buffer(), enc.BufferSize()}, output));
  return WP2_STATUS_OK;
}

WP2Status EncodeTiles(uint32_t width, uint32_t height,
                      const ArgbBuffer& rgb_buffer, const YUVPlane& yuv_buffer,
                      const CSPTransform& transf, const EncoderConfig& config,
                      uint32_t quality_hint, bool image_has_alpha,
                      Writer* const output) {
  WP2EncDspInit();
  GlobalParams gparams;
  FeatureMap features;
  gparams.features_ = &features;
  ArgbBuffer pre_processed_buffer;
  const ArgbBuffer* buffer = &rgb_buffer;

  // near-lossless pre-processing
  constexpr uint32_t kMinNearLosslessDimension = 16;    // minimum dimension
  if (yuv_buffer.IsEmpty() &&
      rgb_buffer.width >= kMinNearLosslessDimension &&
      rgb_buffer.height >= kMinNearLosslessDimension &&
      config.quality > kMaxLossyQuality && config.quality < kMaxQuality) {
    WP2_CHECK_STATUS(PreprocessNearLossless(
        rgb_buffer, config, /*is_alpha=*/false, &pre_processed_buffer));
    buffer = &pre_processed_buffer;
  }

  WP2_CHECK_STATUS(
      GlobalAnalysis(*buffer, yuv_buffer, transf, config, &gparams));

  // TODO(skal): analysis and num_mtx_ decision...
  WP2_CHECK_STATUS(gparams.InitRndMtxSet());
  WP2_CHECK_STATUS(gparams.Init());

  WP2_CHECK_STATUS(
      EncodeGLBL(config, gparams, quality_hint, image_has_alpha, output));

  const uint32_t tile_width = TileWidth(FinalTileShape(config), width);
  const uint32_t tile_height =
      TileHeight(FinalTileShape(config), /*image_width=*/width);

  EncTilesLayout tiles_layout;
  WP2_CHECK_STATUS(GetEncTilesLayout(width, height, tile_width, tile_height,
                                     *buffer, yuv_buffer, gparams,
                                     &tiles_layout));

  Vector<TileEncoder> workers;
  WP2_CHECK_STATUS(SetupWorkers(width, height, *buffer, yuv_buffer, config,
                                gparams, &tiles_layout, &workers));
  WP2_CHECK_STATUS(CodeTiles(&tiles_layout, &workers, output));
  return WP2_STATUS_OK;
}

WP2Status EncodeMetadata(const Metadata& metadata, Writer* const output) {
  // put a terminating signature to count extra chunks and
  // finish off with trailing metadata
  const bool has_xmp = (metadata.xmp.size > 0);
  const bool has_exif = (metadata.exif.size > 0);
  const uint32_t content_bits = (has_xmp ? 1u : 0u) | (has_exif ? 2u : 0u);
  if (content_bits > 0) {
    WP2_CHECK_OK(WriteTag(kTagMask | (content_bits << 16), output),
                 WP2_STATUS_BAD_WRITE);
    WP2_CHECK_STATUS(
        WriteChunk({metadata.xmp.bytes, metadata.xmp.size}, output));
    WP2_CHECK_STATUS(
        WriteChunk({metadata.exif.bytes, metadata.exif.size}, output));
  }
  return WP2_STATUS_OK;
}

TileShape FinalTileShape(const EncoderConfig& config) {
  if (config.tile_shape == TILE_SHAPE_AUTO) {
    return config.quality > kMaxLossyQuality ? TILE_SHAPE_SQUARE_256
                                             : TILE_SHAPE_SQUARE_512;
  }
  return config.tile_shape;
}

PartitionMethod FinalPartitionMethod(const EncoderConfig& config,
                                     uint32_t tile_width,
                                     uint32_t tile_height) {
  if (config.partition_method == AUTO_PARTITIONING) {
    if (tile_width <= kMinBlockSizePix && tile_height <= kMinBlockSizePix) {
      return ALL_4X4_PARTITIONING;  // Skip any setup.
    }
    if (tile_width < 32 && tile_height < 32) {
      // Small enough to use slow methods.
      if (config.speed == 0) return ALL_16X16_PARTITIONING;
      if (config.speed <= 2) return MULTIPASS_PARTITIONING;
      if (!config.partition_snapping) return MULTIPASS_PARTITIONING;
      if (config.speed <= 5) return AREA_ENCODE_PARTITIONING;
      if (config.speed <= 7) return TILE_ENCODE_PARTITIONING;
      if (tile_width <= 16 && tile_height <= 16) return EXHAUSTIVE_PARTITIONING;
      return TILE_ENCODE_PARTITIONING;
    }
    if (config.speed == 0) return ALL_16X16_PARTITIONING;
    if (config.speed <= 6) return MULTIPASS_PARTITIONING;
    if (!config.partition_snapping) return MULTIPASS_PARTITIONING;
    return AREA_ENCODE_PARTITIONING;
  }
  return config.partition_method;
}

//------------------------------------------------------------------------------
// Still image encoding function

WP2Status Encode(const ArgbBuffer& input, Writer* output,
                 const EncoderConfig& config, PictureHint picture_hint) {
  WP2_CHECK_OK(output != nullptr, WP2_STATUS_NULL_PARAMETER);

  WP2_CHECK_OK(config.IsValid(), WP2_STATUS_INVALID_CONFIGURATION);

  WP2_CHECK_OK(input.format == WP2_Argb_32 || (input.format == WP2_Argb_38 &&
                                               config.quality == kMaxQuality),
               WP2_STATUS_INVALID_COLORSPACE);
  WP2_CHECK_OK((input.width > 0) && (input.height > 0) &&
                   (input.width <= kImageDimMax) &&
                   (input.height <= kImageDimMax),
               WP2_STATUS_BAD_DIMENSION);

  (void)picture_hint;
  WP2_CHECK_OK(output != nullptr, WP2_STATUS_NULL_PARAMETER);
  WP2_CHECK_STATUS(SetupEncoderInfo(input.width, input.height, config));

  const RGB12b preview_color = GetPreviewColor(input);
  const bool has_alpha = input.HasTransparency();
  const bool has_icc = (input.metadata.iccp.size > 0);
  const bool has_trailing_data =
      (input.metadata.xmp.size > 0) || (input.metadata.exif.size > 0);
  WP2_CHECK_STATUS(EncodeHeader(config, input.width, input.height,
                                has_alpha, /*is_anim=*/false, /*loop_count=*/0,
                                kTransparentArgb38b, preview_color, has_icc,
                                has_trailing_data, output));

  WP2_CHECK_STATUS(EncodePreview(input, config, output));
  WP2_CHECK_STATUS(
      EncodeICC({input.metadata.iccp.bytes, input.metadata.iccp.size}, output));

  const GlobalParams::Type type = DecideGlobalParamsType(config);
  const bool yuv_is_needed = (type != GlobalParams::GP_LOSSLESS);

  YUVPlane yuv_input;
  CSPTransform csp_transform;
  if (yuv_is_needed) {
    WP2_CHECK_STATUS(csp_transform.Init(config.csp_type, input));
    WP2_CHECK_STATUS(yuv_input.Import(input, has_alpha, csp_transform,
                                      /*resize_if_needed=*/true,
                                      /*pad=*/kPredWidth));
  }

  WP2_CHECK_STATUS(
      EncodeTiles(input.width, input.height, input, yuv_input, csp_transform,
                  config, GetQualityHint(config.quality), has_alpha, output));
  WP2_CHECK_STATUS(EncodeMetadata(input.metadata, output));
  return WP2_STATUS_OK;
}

//------------------------------------------------------------------------------

WP2Status Encode(uint32_t width, uint32_t height,
                 const int16_t* c0_buffer, uint32_t c0_step,
                 const int16_t* c1_buffer, uint32_t c1_step,
                 const int16_t* c2_buffer, uint32_t c2_step,
                 const int16_t* a_buffer, uint32_t a_step,
                 const int16_t ccsp_to_rgb_matrix[9],
                 uint32_t ccsp_to_rgb_shift, Writer* output,
                 const EncoderConfig& config, const Metadata& metadata) {
  WP2_CHECK_OK(ccsp_to_rgb_matrix != nullptr, WP2_STATUS_NULL_PARAMETER);
  WP2_CHECK_OK(output != nullptr, WP2_STATUS_NULL_PARAMETER);
  const CSPMtx ccsp_to_rgb(ccsp_to_rgb_matrix, ccsp_to_rgb_shift);

  WP2_CHECK_OK(config.IsValid(), WP2_STATUS_INVALID_CONFIGURATION);

  WP2_CHECK_OK((width > 0 && height > 0) &&
                   (width <= kImageDimMax && height <= kImageDimMax),
               WP2_STATUS_BAD_DIMENSION);
  const bool has_alpha = (a_buffer != nullptr);

  WP2_CHECK_OK(c0_buffer != nullptr, WP2_STATUS_NULL_PARAMETER);
  WP2_CHECK_OK(c1_buffer != nullptr, WP2_STATUS_NULL_PARAMETER);
  WP2_CHECK_OK(c2_buffer != nullptr, WP2_STATUS_NULL_PARAMETER);

  WP2_CHECK_OK(c0_step >= width && c1_step >= width && c2_step >= width,
               WP2_STATUS_BAD_DIMENSION);
  if (has_alpha) {
    WP2_CHECK_OK(a_step >= width, WP2_STATUS_BAD_DIMENSION);
  }

  WP2_CHECK_OK(ccsp_to_rgb.shift <= 16, WP2_STATUS_INVALID_PARAMETER);
  WP2_CHECK_STATUS(SetupEncoderInfo(width, height, config));

  const GlobalParams::Type type = DecideGlobalParamsType(config);
  const bool rgb_is_needed =
      (type != GlobalParams::GP_LOSSY ||   // for lossless
       config.csp_type == Csp::kCustom ||  // for CSPTransform::Optimize()
       config.create_preview);             // for EncodePreview()
  // TODO(yguyon): Some of these cases could also be done directly in YUV space
  //               instead of needing RGB conversion.
  const bool yuv_is_needed =
      (type != GlobalParams::GP_LOSSLESS);  // lossy, higher precision than RGB

  ArgbBuffer rgb_input(WP2_Argb_32);
  YUVPlane yuv_input;
  CSPTransform csp_transform;
  if (rgb_is_needed) {
    WP2_CHECK_STATUS(rgb_input.Resize(width, height));
    WP2_CHECK_STATUS(CSPTransform::CustomToArgb(
        width, height, c0_buffer, c0_step, c1_buffer, c1_step, c2_buffer,
        c2_step, a_buffer, a_step, ccsp_to_rgb, &rgb_input));
  }
  if (yuv_is_needed) {
    WP2_CHECK_STATUS(csp_transform.Init(config.csp_type, rgb_input));
    WP2_CHECK_STATUS(
        yuv_input.Resize(width, height, /*pad=*/kPredWidth, has_alpha));

    WP2_CHECK_STATUS(csp_transform.CustomToYUV(
        width, height, c0_buffer, c0_step, c1_buffer, c1_step, c2_buffer,
        c2_step, ccsp_to_rgb, yuv_input.Y.Row(0), yuv_input.Y.Step(),
        yuv_input.U.Row(0), yuv_input.U.Step(), yuv_input.V.Row(0),
        yuv_input.V.Step()));
    if (has_alpha) {
      Plane16 non_padded_alpha;
      WP2_CHECK_STATUS(
          non_padded_alpha.SetView(yuv_input.A, {0, 0, width, height}));
      non_padded_alpha.From(a_buffer, a_step);
    }

    WP2_CHECK_STATUS(yuv_input.FillPad(width, height));
  }

  const RGB12b preview_color =
      yuv_input.IsEmpty() ? GetPreviewColor(rgb_input)
                          : GetPreviewColor(yuv_input, csp_transform);
  const bool has_icc = (metadata.iccp.size > 0);
  const bool has_trailing_data =
      (metadata.xmp.size > 0 || metadata.exif.size > 0);
  WP2_CHECK_STATUS(EncodeHeader(
      config, width, height, has_alpha, /*is_anim=*/false, /*loop_count=*/0,
      kTransparentArgb38b, preview_color, has_icc, has_trailing_data, output));

  WP2_CHECK_STATUS(EncodePreview(rgb_input, config, output));
  WP2_CHECK_STATUS(
      EncodeICC({metadata.iccp.bytes, metadata.iccp.size}, output));

  WP2_CHECK_STATUS(
      EncodeTiles(width, height, rgb_input, yuv_input, csp_transform, config,
                  GetQualityHint(config.quality), has_alpha, output));
  WP2_CHECK_STATUS(EncodeMetadata(metadata, output));
  return WP2_STATUS_OK;
}

//------------------------------------------------------------------------------

}  // namespace WP2
