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
// WP2 decoder: internal header.
//
// Author: Skal (pascal.massimino@gmail.com)

#ifndef WP2_DEC_WP2_DEC_I_H_
#define WP2_DEC_WP2_DEC_I_H_

#include "src/common/global_params.h"
#include "src/common/lossy/block.h"
#include "src/common/lossy/block_size_io.h"
#include "src/common/lossy/context.h"
#include "src/common/lossy/residuals.h"
#include "src/common/lossy/rnd_mtx.h"
#include "src/common/progress_watcher.h"
#include "src/dec/lossless/losslessi_dec.h"
#include "src/dec/symbols_dec.h"
#include "src/dec/residuals_dec_aom.h"
#include "src/dec/tile_dec.h"
#include "src/dsp/dsp.h"
#include "src/utils/ans.h"
#include "src/utils/csp.h"
#include "src/utils/data_source.h"
#include "src/utils/front_mgr.h"
#include "src/wp2/decode.h"
#include "src/wp2/format_constants.h"

namespace WP2 {

//------------------------------------------------------------------------------
// Decoding functions.
// "dec" is the input (ANS decoder), with "features" and "config".
// "picture" is the output (buffer view of the tile to decode).
// "num_decoded_rows" is for incremental decoding; kept up-to-date.
//   A row must be left intact once it is included in num_decoded_rows.

WP2Status LossyDecode(const BitstreamFeatures& features,
                      const DecoderConfig& config,
                      TilesLayout* const tiles_layout,
                      ProgressWatcher* const progress, ANSDec* const dec,
                      Tile* const tile);

WP2Status LosslessDecode(const BitstreamFeatures& features,
                         const DecoderConfig& config,
                         const GlobalParams& gparams,
                         ProgressWatcher* const progress, ANSDec* const dec,
                         Tile* const tile);

WP2Status NeuralDecode(ArgbBuffer* const picture, ANSDec* const dec,
                       const BitstreamFeatures& features);

//------------------------------------------------------------------------------

// These internal methods are to avoid polluting public decode.h
void ResetDecoderInfo(DecoderInfo* const info);
WP2Status MergeDecoderInfo(DecoderInfo* const from, DecoderInfo* const to);

//------------------------------------------------------------------------------
// Internal decode.

// Output of an ANMF chunk header decoding.
struct AnimationFrame {
  bool dispose = true;       // Fill the whole canvas with the background color
                             // before decoding this frame.
                             // Leave the canvas as is if false.
  bool blend = false;        // Blend frame into current canvas, based on alpha.
                             // Overwrite pixels (including alpha) if false.
  uint32_t duration_ms = 0;  // Milliseconds during which this frame is visible.
  Rectangle window;          // The canvas area modified by this frame.
  bool is_last = false;      // If true, there's no frame after this one.
};

// Initialize and register debug data.
WP2Status SetupDecoderInfo(const BitstreamFeatures& features,
                           const DecoderConfig& config);
void RegisterChunkSize(const DecoderConfig& config, uint32_t chunk_size,
                       uint32_t chunk_size_size);

// Decode steps are split into functions for incremental decoding. These
// functions only call DataSource::MarkNumBytesAsRead(). This way the calling
// function can decide whether to retry later (UnmarkAllReadBytes()) or not
// (Discard(GetNumNextBytes())).

// Fills DecoderInfo in 'config' based on provided GlobalParams.
void FillDecoderInfo(const GlobalParams& gparams, const DecoderConfig& config);

WP2Status DecodeHeader(DataSource* const data_source,
                       BitstreamFeatures* const features);

WP2Status DecodePreview(DataSource* const data_source,
                        const BitstreamFeatures& features,
                        ArgbBuffer* const output_buffer);
WP2Status SkipPreview(DataSource* const data_source,
                      const BitstreamFeatures& features);

WP2Status DecodeICC(DataSource* const data_source,
                    const BitstreamFeatures& features, Data* const iccp);

WP2Status DecodeANMF(DataSource* const data_source,
                     const BitstreamFeatures& features, uint32_t frame_index,
                     AnimationFrame* const frame);

// if gparams = nullptr, just skip over the data without parsing
WP2Status DecodeGLBL(DataSource* const data_source,
                     const DecoderConfig& config,
                     const BitstreamFeatures& features,
                     GlobalParams* const gparams = nullptr);

// Fills the area outside 'frame.window' with the 'features.background_color'.
WP2Status FillBorders(const BitstreamFeatures& features,
                      const AnimationFrame& frame, ArgbBuffer* const output);
WP2Status FillBorders(const BitstreamFeatures& features,
                      const AnimationFrame& frame,
                      const CSPTransform& csp_transform,
                      YUVPlane* const output);

WP2Status DecodeTileChunkSize(const GlobalParams& params,
                              const Rectangle& tile_rect,
                              DataSource* const data_source,
                              size_t* const tile_chunk_size);

WP2Status DecodeMetadata(DataSource* const data_source,
                         const BitstreamFeatures& features, Data* const exif,
                         Data* const xmp);

// Returns true if variable is successfully read, false otherwise.
bool ReadVarInt(DataSource* const data_source, uint64_t buf_size_upper,
                size_t* const size);

// Extracts a chunk. 'chunk_data' is not a copy and must be used before reading
// from the 'data_source' again. If 'chunk_data' is null, only the chunk size is
// read and the chunk itself is skipped (marked as read).
WP2Status ReadChunk(DataSource* const data_source, DataView* const chunk_data);
// Reads the chunk size and set that many following bytes as read.
WP2Status SkipChunk(DataSource* const data_source);

// These two Skip() functions do not attempt to read anything past the chunk
// size; if data length check through reading is necessary but data copy should
// be avoided, call the Decode() functions instead with null Data pointers.
WP2Status SkipICC(DataSource* const data_source,
                  const BitstreamFeatures& features);
WP2Status SkipMetadata(DataSource* const data_source,
                       const BitstreamFeatures& features);

//------------------------------------------------------------------------------

// Class specialized in reading transform residuals.
class ResidualReader : public ResidualIO {
 public:
  ResidualReader() = default;
  // Updates the symbol reader with what will be needed to read residuals.
  WP2Status ReadHeader(SymbolReader* const sr, uint32_t num_coeffs_max_y,
                       uint32_t num_coeffs_max_uv, uint32_t num_transforms,
                       bool has_alpha, bool has_lossy_alpha);

  // Reads residual coefficients.
  WP2Status ReadCoeffs(Channel channel, ANSDec* const dec,
                       SymbolReader* const sr, CodedBlock* const cb,
                       libgav1::AOMContext* const aom_context,
                       BlockInfo* const info);

 private:
  WP2Status ReadHeaderForResidualSymbols(uint32_t num_coeffs_max,
                                         Channel channel,
                                         SymbolReader* const sr);
  WP2Status ReadCoeffsMethod01(Channel channel, uint32_t tf_i,
                               ANSDec* const dec, SymbolReader* const sr,
                               CodedBlock* const cb, BlockInfo* const info);

  uint32_t num_channels_;
  libgav1::AOMResidualReader aom_reader_;
};

//------------------------------------------------------------------------------
// AlphaReader: Class for reading the alpha plane.

class AlphaReader : public WP2Allocable {
 public:
  // Image width and height are in pixels.
  explicit AlphaReader(ANSDec* const dec, const Tile& tile);
  WP2Status Allocate();

  WP2Status ReadHeader(const GlobalParams& gparams);

  // Reads alpha information for the given block from the stream.
  void GetBlockHeader(SymbolReader* const sr, CodedBlock* const cb);
  // Reads alpha information for the given block from the stream.
  WP2Status GetBlock(CodedBlock* const cb);
  // Reconstructs the alpha plane for the given block.
  void Reconstruct(CodedBlock* const cb) const;
  // Copies the lossless alpha data for this block to the given plane (at 0, 0).
  void ReconstructLossless(CodedBlock* const cb, Plane16* const plane) const;

  AlphaMode GetAlphaMode() const { return alpha_mode_; }

 private:
  WP2Status ReadLossless(CodedBlock* const cb);

  const GlobalParams* gparams_ = nullptr;
  // ANS reader
  ANSDec* const dec_;
  // Alpha buffer. Alpha values are stored as a grayscale image.
  // TODO(maryla): this is pretty memory inefficient: in reality we don't need
  // four channels, and we don't need to keep the whole alpha image in memory
  // at once (we never need more than the max height of a block).
  Tile alpha_tile_;
  // Whether alpha is encoded with pure lossy, pure lossless, or a mix.
  AlphaMode alpha_mode_;
  // Next line to read in pixels.
  uint32_t next_line_ = 0;
  // Whether the next processed block is the first block of the image.
  bool first_block_ = true;

  // Block alpha mode predictor when using lossy.
  AlphaModePredictor mode_predictor_;

  // *** Members variables for lossless decoding. ***
  WP2L::Decoder lossless_decoder_;
};

//------------------------------------------------------------------------------
// SyntaxReader: Class for parsing syntactic elements.

class SyntaxReader {
 public:
  SyntaxReader(ANSDec* const dec, const Rectangle& rect);

  WP2Status ReadHeader(const BitstreamFeatures& features,
                       const GlobalParams& gparams, const Tile& tile);

  // Reads 'cb' size, segment id, coeffs etc. and sets 'out' as the
  // reconstructed pixels buffer destination.
  WP2Status GetBlock(CodedBlock* const cb, YUVPlane* const out,
                     BlockInfo* const info);
  BlockSize GetBlockSize(const FrontMgrNxNBase& mgr);

  PartitionSet GetPartitionSet() const { return gparams_->partition_set_; }
  bool GetPartitionSnap() const { return gparams_->partition_snapping_; }

  const AlphaReader& alpha_reader() const { return *alpha_reader_; }

  void ApplyDecoderConfig(const DecoderConfig& config);

 private:
  void ReadSplitTransform(Channel channel, CodedBlock* const cb);
  void ReadHasCoeffs(Channel channel, CodedBlock* const cb);
  void ReadTransform(Channel channel, CodedBlock* const cb);

  // Fill the is420_ data in the CodedBlock.
  void ReadIs420(CodedBlock* const cb);

  WP2Status LoadDictionaries();

  WP2Status ReadPredModes(CodedBlock* const cb);

  const GlobalParams* gparams_ = nullptr;

  ANSDec* const dec_;

  uint32_t num_blocks_, num_transforms_;
  ChromaSubsampling chroma_subsampling_;

  SegmentIdPredictor segment_ids_;

  ResidualReader residual_reader_;

  std::unique_ptr<AlphaReader> alpha_reader_;

  SymbolsInfo symbols_info_;
  SymbolReader sr_;
  BlockContext context_;
  uint32_t width_;
  uint32_t height_;
  bool use_aom_coeffs_;

 public:
  // Visual debug: stores bit cost for a block
  void StoreBitCost(const WP2::DecoderConfig& config, uint32_t tile_x,
                    uint32_t tile_y, const Block& block,
                    Plane16* const dst_plane) const;
  // Decodes the 'block' boundaries and the 'pixels'. Check that they match.
  void ReadAndCompareRawPixels(const Block& block, const YUVPlane& pixels,
                               ANSDec* const dec);
};

//------------------------------------------------------------------------------

// Returns true if 'DecoderInfo::visual_debug' in the 'config' exists and
// contains the 'token'. Tokens are separated by '/'.
bool VDMatch(const DecoderConfig& config, const char token[]);
// Returns the channel corresponding to the current debug view. Asserts that
// the debug view contains the name of a channel ("y", "u", "v" or "a") in its
// path.
Channel VDChannel(const DecoderConfig& config);
// Returns true if the given rectangle is selected in visual debug.
bool VDSelected(uint32_t tile_x, uint32_t tile_y, const Rectangle& rect,
                const DecoderConfig& config);

// Handles VisualDebug for "compressed" and "original".
void ApplyVDebugBeforeAfter(const DecoderConfig& config,
                            const CSPTransform& csp_tranform,
                            const Tile& tile,
                            ArgbBuffer* const debug_output);

//------------------------------------------------------------------------------

}  // namespace WP2

#endif  // WP2_DEC_WP2_DEC_I_H_
