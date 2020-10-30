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
// WP2 lossy decoding of the alpha plane.
//
// Author: Maryla (maryla@google.com)

#include "src/dec/wp2_dec_i.h"
#include "src/utils/ans_utils.h"

namespace WP2 {

AlphaReader::AlphaReader(ANSDec* const dec, const Tile& tile)
    : dec_(dec), lossless_decoder_() {
  alpha_tile_.chunk_size_is_known = tile.chunk_size_is_known;
  alpha_tile_.chunk_size = tile.chunk_size;
  alpha_tile_.input = tile.input;
  alpha_tile_.rect = tile.rect;
  alpha_tile_.num_decoded_rows = 0;
}

WP2Status AlphaReader::Allocate() {
  return alpha_tile_.rgb_output.Resize(alpha_tile_.rect.width,
                                       alpha_tile_.rect.height);
}

WP2Status AlphaReader::ReadHeader(const GlobalParams& gparams) {
  ANSDebugPrefix prefix(dec_, "Alpha");
  alpha_mode_ = gparams.maybe_use_lossy_alpha_
                    ? (AlphaMode)dec_->ReadRValue(kAlphaModeNum, "alpha_mode")
                    : kAlphaModeLossless;
  gparams_ = &gparams;
  if (alpha_mode_ == kAlphaModeLossy) {
    WP2_CHECK_STATUS(
        mode_predictor_.Init(alpha_tile_.rect.width, alpha_tile_.rect.height));
  }
  return WP2_STATUS_OK;
}

void AlphaReader::GetBlockHeader(SymbolReader* const sr, CodedBlock* const cb) {
  if (alpha_mode_ == kAlphaModeLossy) {
    cb->alpha_mode_ = (BlockAlphaMode)mode_predictor_.Read(*cb, sr);
  } else {
    cb->alpha_mode_ = kBlockAlphaLossless;
  }
}

WP2Status AlphaReader::GetBlock(CodedBlock* const cb) {
  ANSDebugPrefix prefix(dec_, "Alpha");
  if (alpha_mode_ == kAlphaModeLossless) {
    WP2_CHECK_STATUS(ReadLossless(cb));
  }
  return WP2_STATUS_OK;
}

void AlphaReader::Reconstruct(CodedBlock* const cb) const {
  switch (cb->alpha_mode_) {
    case kBlockAlphaLossless:
      ReconstructLossless(cb, &cb->out_.A);
      break;
    case kBlockAlphaLossy:
      BlockCoeffs32 res;
      cb->Dequantize(gparams_->a_segments_[0].quant_a_, kAChannel, &res);
      cb->Reconstruct(kAChannel, /*reduced_transform=*/false, &res);
      break;
    case kBlockAlphaFullTransp:
      cb->out_.A.Fill(0);
      break;
    case kBlockAlphaFullOpaque:
      cb->out_.A.Fill(kAlphaMax);
      break;
  }
}

void AlphaReader::ReconstructLossless(CodedBlock* const cb,
                                      Plane16* const plane) const {
  // Copy from buffer to block.
  const uint32_t block_width = cb->w_pix();
  const uint32_t block_height = cb->h_pix();
  for (uint32_t y = 0, y_offset = cb->y_pix(); y < block_height;
       ++y, ++y_offset) {
    for (uint32_t x = 0, x_offset = cb->x_pix(); x < block_width;
         ++x, ++x_offset) {
      if (x_offset >= alpha_tile_.rect.width) {
        plane->At(x, y) = plane->At(x - 1, y);
      } else if (y_offset >= alpha_tile_.rect.height) {
        plane->At(x, y) = plane->At(x, y - 1);
      } else {
        const uint8_t* const alpha_row =
            (uint8_t*)alpha_tile_.rgb_output.GetRow(y_offset);
        plane->At(x, y) = alpha_row[4 * x_offset + 1];
      }
    }
  }
}

WP2Status AlphaReader::ReadLossless(CodedBlock* const cb) {
  const uint32_t max_line =
      std::min(cb->y_pix() + cb->h_pix(), alpha_tile_.rect.height);

  if (next_line_ < max_line) {
    // Read as many new lines as we need for this block.
    const uint32_t lines_to_read = max_line - next_line_;
    const size_t num_used_bytes_before = dec_->GetNumUsedBytes();

    if (first_block_) {
      // Initialize the lossless decoder and read the header.
      lossless_decoder_.Init(DecoderConfig(), *gparams_, /*progress=*/nullptr,
                             dec_, &alpha_tile_);
      WP2_CHECK_STATUS(lossless_decoder_.DecodeHeader());

      first_block_ = false;
    }

    WP2_CHECK_STATUS(lossless_decoder_.DecodeLines(lines_to_read));

    const size_t min_num_used_bytes = ANSDec::GetMinNumUsedBytesDiff(
        num_used_bytes_before, dec_->GetNumUsedBytes());
    cb->alpha_lossless_bytes_ = min_num_used_bytes;

    next_line_ = max_line;
  } else {
    cb->alpha_lossless_bytes_ = 0;
  }

  return WP2_STATUS_OK;
}

}  // namespace WP2
