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
// ContextTileDecoder implementation.
//
// Author: Yannis Guyon (yguyon@google.com)

#include "src/dec/incr/decoder_context.h"

#if defined(WP2_USE_CONTEXT_SWITCH) && (WP2_USE_CONTEXT_SWITCH > 0)

namespace WP2 {

//------------------------------------------------------------------------------

void ContextTileDecoder::SuspendingDecode(LocalContext* const context) {
  InterContextData* const inter_context_data =
      (InterContextData*)context->GetInterContextData();
  inter_context_data->input.SetContext(context);  // Enable suspension.

  inter_context_data->status = inter_context_data->worker.Start(false);

  const WP2Status worker_status = inter_context_data->worker.End();
  if (inter_context_data->status == WP2_STATUS_OK) {
    inter_context_data->status = worker_status;  // Return first error.
  }

  inter_context_data->input.SetContext(nullptr);  // Disable suspension.
  context->Close();
}

//------------------------------------------------------------------------------

void ContextTileDecoder::DecodePartialTile() {
  const uint8_t* partial_tile_data = partial_tile_->data_handle.GetBytes();
  const size_t partial_tile_data_size = partial_tile_->data_handle.GetSize();
  inter_context_data_.input.Update(partial_tile_data, partial_tile_data_size);
  if (inter_context_data_.input.HasEnoughDataToResume() &&
      !main_context_.Resume()) {
    // Not being able to Resume() is probably related to memory.
    inter_context_data_.status = WP2_STATUS_OUT_OF_MEMORY;
  }
}

//------------------------------------------------------------------------------

WP2Status ContextTileDecoder::Init(const BitstreamFeatures& features,
                                   const DecoderConfig& config,
                                   TilesLayout* const tiles_layout,
                                   Tile* const partial_tile) {
  auto& worker = inter_context_data_.worker;   // shortcut
  worker.features_ = &features;
  worker.config_ = &config;
  worker.tiles_layout_ = tiles_layout;
  worker.gparams_ = tiles_layout->gparams;
  worker.tile_ = partial_tile;
  worker.self_assigning_ = false;   // Isolated worker.
  inter_context_data_.input = SuspendableDataSource();  // Clear.
  partial_tile_ = partial_tile;
  partial_tile_->input = &inter_context_data_.input;
  return WP2_STATUS_OK;
}

WP2Status ContextTileDecoder::Start() {
  if (main_context_.CreateLocalContext(&SuspendingDecode,
                                       &inter_context_data_)) {
    DecodePartialTile();
  } else {
    // Not being able to CreateLocalContext() is probably related to memory.
    inter_context_data_.status = WP2_STATUS_OUT_OF_MEMORY;
  }
  return inter_context_data_.status;
}

WP2Status ContextTileDecoder::Continue() {
  DecodePartialTile();

  if ((inter_context_data_.status == WP2_STATUS_OK) &&
      (!main_context_.IsLocalContextClosed())) {
    return WP2_STATUS_NOT_ENOUGH_DATA;  // Tile is incomplete.
  }
  return inter_context_data_.status;
}

void ContextTileDecoder::Clear() {
  main_context_.CloseLocalContext();
  inter_context_data_.worker.tiles_layout_ = nullptr;
  inter_context_data_.worker.tile_ = nullptr;
  if (partial_tile_ != nullptr) partial_tile_->input = nullptr;
  partial_tile_ = nullptr;
  inter_context_data_.input = SuspendableDataSource();  // Clear.
  inter_context_data_.status = WP2_STATUS_OK;
}

Tile* ContextTileDecoder::GetPartialTile() const {
  return partial_tile_;
}

}  // namespace WP2

#endif  // WP2_USE_CONTEXT_SWITCH
