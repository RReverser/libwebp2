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
//   WP2 encoder: internal header.
//
// Author: Skal (pascal.massimino@gmail.com)

#ifndef WP2_WP2_ENC_I_H_
#define WP2_WP2_ENC_I_H_

#include <memory>

#include "src/common/global_params.h"
#include "src/common/lossy/block.h"
#include "src/common/lossy/block_size_io.h"
#include "src/common/lossy/context.h"
#include "src/common/lossy/predictor.h"
#include "src/common/lossy/residuals.h"
#include "src/common/lossy/rnd_mtx.h"
#include "src/dsp/dsp.h"
#include "src/enc/analysis.h"
#include "src/enc/lossless/losslessi_enc.h"
#include "src/enc/residuals_enc_aom.h"
#include "src/enc/symbols_enc.h"
#include "src/utils/ans.h"
#include "src/utils/front_mgr.h"
#include "src/utils/plane.h"
#include "src/utils/thread_utils.h"
#include "src/utils/wiener.h"
#include "src/wp2/encode.h"

namespace WP2 {

// Encodes 'size' (< 'size_upper') into a variable length integer.
// Returns the number of bytes used.
size_t WriteVarInt(uint64_t size, uint64_t size_upper,
                   uint8_t dst[kMaxVarIntLength]);

// Resets 'EncoderConfig::info' if any. Not thread-safe.
WP2Status SetupEncoderInfo(uint32_t width, uint32_t height,
                           const EncoderConfig& config);

// Writes WP2 header to 'output' based on 'config' and other arugments.
WP2Status EncodeHeader(const EncoderConfig& config,
                       uint32_t width, uint32_t height,
                       bool has_alpha, bool is_anim, uint32_t loop_count,
                       Argb38b background_color, RGB12b preview_color,
                       bool has_icc, bool has_trailing_data,
                       Writer* const output);

// Creates and compresses a preview from 'buffer' based on 'config' and writes
// it to 'output'.
WP2Status EncodePreview(const ArgbBuffer& buffer, const EncoderConfig& config,
                        Writer* const output);

// Writes 'iccp' to 'output'.
WP2Status EncodeICC(DataView iccp, Writer* const output);

// Writes global parameters 'gparams' to 'output'.
WP2Status EncodeGLBL(const EncoderConfig& config, const GlobalParams& gparams,
                     uint32_t quality_hint, bool image_has_alpha,
                     Writer* const output);

// Compresses an image into tiles based on 'config' and writes them to 'output'.
// The input is either 'rgb_buffer' or 'yuv_buffer' depending on what is needed.
WP2Status EncodeTiles(uint32_t width, uint32_t height,
                      const ArgbBuffer& rgb_buffer, const YUVPlane& yuv_buffer,
                      const CSPTransform& transf, const EncoderConfig& config,
                      uint32_t quality_hint, bool image_has_alpha,
                      Writer* const output);

// Writes 'metadata' to 'output'.
WP2Status EncodeMetadata(const Metadata& metadata, Writer* const output);

// Returns the TileShape, converting from TILE_SHAPE_AUTO to a concrete
// one if necessary.
TileShape FinalTileShape(const EncoderConfig& config);

// Returns the partition method used at encoding.
// 'tile_width' and 'tile_height' are in pixels, padded or not.
PartitionMethod FinalPartitionMethod(const EncoderConfig& config,
                                     uint32_t tile_width, uint32_t tile_height);

//------------------------------------------------------------------------------

// Choose encoding settings based on the given 'config'.
ChromaSubsampling DecideChromaSubsampling(const EncoderConfig& config,
                                          bool more_than_one_block);
bool DecideAOMCoeffs(const EncoderConfig& config, const Rectangle& tile_rect);
WP2Status DecideTransforms(const EncoderConfig& config,
                           Vector<TransformPair>* const transforms,
                           Vector<TransformPair>* const transforms_subset);

//------------------------------------------------------------------------------

struct EncTile {
  // Tile position and dimensions within the frame in pixels, not padded.
  Rectangle rect;

  // Views of this tile.
  ArgbBuffer rgb_input;  // Not padded.
  YUVPlane yuv_input;    // Padded.

  // Output. Assemble() to get the encoded bytes.
  ANSEnc enc;
};

struct EncTilesLayout {
  uint32_t num_tiles_x = 0;
  uint32_t num_tiles_y = 0;
  uint32_t tile_width = 0;
  uint32_t tile_height = 0;
  Vector<EncTile> tiles;
  uint32_t first_unassigned_tile_index = 0;  // Next tile to assign a worker to.
  ThreadLock assignment_lock;                // In case of concurrency.

  // Frame's global parameters (externally owned).
  const GlobalParams* gparams = nullptr;
};

class TileEncoder : public WorkerBase {
 public:
  // Assigns the next 'tile_' to this worker if possible.
  WP2Status AssignNextTile();
  // As long as 'tile_' is not null, decode it and AssignNextTile().
  WP2Status Execute() override;

  // Compresses the tile into 'enc'.
  WP2Status LossyEncode(const VectorNoCtor<Block>& forced_partition,
                        ANSEnc* const enc);
  WP2Status LosslessEncode(ANSEnc* const enc);

  const EncoderConfig* config_ = nullptr;
  bool use_lossless_ = false;
  EncTilesLayout* tiles_layout_ = nullptr;  // Used to get next assignable tile.
  EncTile* tile_ = nullptr;                 // Currently assigned tile.

  // Debug output.
  EncoderInfo* info_ = nullptr;
};

// Encodes pixels as row-ordered raw samples.
WP2Status BypassTileEnc(const GlobalParams& gparams, EncTile* const tile,
                        Writer* const output);

WP2Status NeuralEncode(const ArgbBuffer& buffer, const EncoderConfig& config,
                       ANSEnc* const enc);

//------------------------------------------------------------------------------

// Class specialized in writing transform residuals to a stream.
class ResidualWriter : public ResidualIO {
 public:
  // Initializes the underlying memory.
  WP2Status Init(bool use_aom_coeffs, bool has_alpha);

  // Default copy is enough to create a deep clone of the full state.
  void CopyFrom(const ResidualWriter& other) { operator=(other); }

  // Finds the best encoding method for the given 'num_coeffs'.
  // If 'cost' is not nullptr, it is set to the cost of storing the residuals
  // with this method (coefficients only, not including dictionaries).
  static WP2Status FindBestEncodingMethod(
      TrfSize dim, const int16_t* const coeffs, uint32_t num_coeffs,
      bool first_is_dc, Channel channel, uint32_t num_channels,
      SymbolCounter* const counter, EncodingMethod* const encoding_method,
      float* const cost = nullptr);

  // Records residuals from 'cb'.
  WP2Status RecordCoeffs(const CodedBlock& cb, Channel channel,
                         SymbolRecorder* const recorder,
                         libgav1::AOMContext* const aom_context);

  // Writes the data (e.g. dictionaries) needed to interpret the symbols that
  // will later be written one by one.
  WP2Status WriteHeader(uint32_t num_coeffs_max_y, uint32_t num_coeffs_max_uv,
                        uint32_t num_transforms, bool has_lossy_alpha,
                        const SymbolRecorder& recorder,
                        ANSDictionaries* const dicts, ANSEncBase* const enc,
                        SymbolWriter* const sw);

  // Writes residuals from 'cb'.
  WP2Status WriteCoeffs(const CodedBlock& cb, Channel channel,
                          ANSEncBase* enc, SymbolManager* const sm,
                          libgav1::AOMContext* const aom_context);

  // Returns a (very) crude estimation of the number of bits needed to code
  // the residuals.
  // The logic for storing the residuals is the same as WriteCoeffs
  // except symbol statistics and adaptive symbols do not use global statistics
  // but fixed ones.
  static WP2Status GetPseudoRate(Channel channel, uint32_t num_channels,
                                 TrfSize dim, const int16_t* const coeffs,
                                 uint32_t num_coeffs, bool first_is_dc,
                                 SymbolCounter* const counter,
                                 float* const pseudo_rate);
  // Returns a slightly less crude estimation of the number of bits needed to
  // code the residuals. Much slower than the pseudo rate above.
  static WP2Status GetRate(Channel channel, uint32_t num_channels, TrfSize dim,
                           const int16_t* const coeffs, uint32_t num_coeffs,
                           bool first_is_dc, SymbolCounter* const counter,
                           float* const rate,
                           EncodingMethod* const encoding_method = nullptr);
  static WP2Status GetRateAOM(const CodedBlock& cb, Channel channel,
                              const libgav1::AOMContext& aom_context,
                              SymbolCounter* const counter, float* const rate);

 private:
  // Calls WriteHeader on SymbolWriter for symbols that do depend on a
  // specific residual method.
  WP2Status WriteHeaderForResidualSymbols(Channel channel,
                                          uint32_t num_coeffs_max,
                                          const SymbolRecorder& recorder,
                                          ANSEncBase* const enc,
                                          SymbolWriter* const sw,
                                          ANSDictionaries* const dicts_in);

  // Stores the DC separately.
  // 'range' is the maximum absolute value of 'v'.
  // 'can_be_zero' specifies whether v can be 0.
  static void StoreDC(Channel channel, uint32_t num_channels,
                      ANSEncBase* const enc, SymbolManager* const sm, int16_t v,
                      bool can_be_zero);
  // Stores residuals.
  // The residuals in 'coeffs' are stored in 'enc' and 'sm'.
  static void StoreCoeffs(const int16_t* const coeffs, uint32_t num_coeffs,
                          bool first_is_dc, TrfSize dim, Channel channel,
                          uint32_t num_channels, EncodingMethod method,
                          ANSEncBase* enc, SymbolManager* const sm,
                          bool is_pseudo_rate);
  uint32_t num_channels_;
  libgav1::AOMResidualWriter aom_writer_;
};

// Class for writing the alpha plane.
class AlphaWriter {
 public:
  // Setup. 'tile_rect' is in pixels, not padded.
  WP2Status Init(const EncoderConfig& config, const GlobalParams& gparams,
                 const BlockContext& context, const YUVPlane& yuv,
                 const Rectangle& tile_rect);

  // Deep copy. 'yuv, context, dicts' must be the clones of 'other' instances.
  WP2Status CopyFrom(const AlphaWriter& other, const BlockContext& context);

  // Decides cb->alpha_mode_. If lossy is used,
  // coefficients will be quantized and an encoding method decided.
  WP2Status DecideAlpha(CodedBlock* const cb,
                        const ResidualWriter& residual_writer,
                        Counters* const counters);

  WP2Status ResetRecord();
  WP2Status Record(const CodedBlock& cb);
  void WriteBlockBeforeCoeffs(const CodedBlock& cb, SymbolManager* const sm,
                              ANSEncBase* const enc);
  WP2Status Write(const CodedBlock& cb, ANSEnc* const enc);

  // Freezes and writes the dictionaries.
  WP2Status WriteHeader(uint32_t num_coeffs_max, ANSEncBase* const enc);

  AlphaMode GetAlphaMode() const { return alpha_mode_; }

 private:
  WP2Status WriteLossless(const CodedBlock& cb, ANSEnc* const enc);

  // Decides whether to use lossy residuals or plain black/white for a given
  // block. If using residuals, encoding params (predictors/transform) and
  // residuals are computed.
  WP2Status ProcessLossy(CodedBlock* const cb,
                         const ResidualWriter& residual_writer,
                         const Vector_f& bits_per_pixel,
                         Counters* const counters);

  BlockAlphaMode GetBlockAlphaMode(const CodedBlock& cb,
                                   const Plane16& alpha) const;

  const GlobalParams* gparams_ = nullptr;
  Rectangle tile_rect_;  // In pixels, not padded.
  EncoderConfig config_;
  AlphaMode alpha_mode_;  // Global alpha mode.

  // *** Members variables for lossy encoding. ***
  // Mode predictor for individual blocks.
  AlphaModePredictor mode_predictor_;

  // *** Members variables for lossless encoding. ***
  // TODO(maryla): ArgbBuffer uses 4 channels even though we only need one...
  ArgbBuffer alpha_;
  // Separate encoder for lossless.
  ANSEnc lossless_enc_;
  WP2L::EncodeInfo lossless_encode_info_;
  // Next line to write.
  uint32_t next_line_ = 0;
  const BlockContext* context_;
};

// Class for holding ANS encoder and dictionaries.
class SyntaxWriter : public WP2Allocable {
 public:
  // 'has_alpha' is the same as buffer.HasTransparency. It is just here to
  // save CPU. 'tile_rect' is in pixels, not padded.
  WP2Status Init(ANSDictionaries* const dicts, const EncoderConfig& config,
                 const GlobalParams& gparams, const YUVPlane& yuv,
                 ChromaSubsampling chroma_subsampling,
                 const Rectangle& tile_rect, uint32_t num_blocks,
                 bool use_aom_coeffs);
  // Deep copy. 'dicts' must be copied beforehand.
  WP2Status CopyFrom(const SyntaxWriter& other, ANSDictionaries* const dicts);
  // Initializes some data based on results from the previous pass. To be called
  // at the beginning of every pass.
  WP2Status InitPass();
  // Records the block's content, except for size and alpha.
  WP2Status Record(const CodedBlock& cb);
  // Records the block's size.
  WP2Status RecordSize(const FrontMgrNxNBase& mgr, BlockSize dim);
  // Writes the headers (frozen dictionaries, etc.).
  WP2Status WriteHeader(ANSEncBase* const enc);
  // These write the block header (segment id, preds, etc.).
  void WriteBlockBeforeCoeffs(const CodedBlock& cb, bool update_ymodes,
                              SymbolManager* const sm, ANSEncBase* const enc);
  // Writes the predictors used by the given block for the given channel.
  static void WriteYAPredictors(const CodedBlock& cb, Channel channel,
                                const YModePredictor* const ymode_predictor,
                                bool update_ymodes, SymbolManager* const sm,
                                ANSEncBase* const enc);
  static void WriteUVPredictors(const CodedBlock& cb, Channel channel,
                                SymbolManager* const sm, ANSEncBase* const enc);
  // Writes whether the block is split into smaller square transforms.
  static void WriteSplitTransform(const CodedBlock& cb, Channel channel,
                                  SymbolManager* const sm,
                                  ANSEncBase* const enc);
  // Writes whether at least one coeff is not zero.
  static void WriteHasCoeffs(const CodedBlock& cb, Channel channel,
                             SymbolManager* const sm, ANSEncBase* const enc);
  // Writes the transform if it is not implicit nor all zero.
  static void WriteTransform(const CodedBlock& cb, Channel channel,
                             SymbolManager* const sm, ANSEncBase* const enc);

  WP2Status WriteBlocks(const Vector<CodedBlock>& cblocks,
                        FrontMgrNxNBase* const mgr, ANSEnc* const enc);

  SymbolWriter* symbol_writer() { return &symbol_writer_; }
  const SymbolRecorder& symbol_recorder() const { return symbol_recorder_; }
  Counters* counters() const { return &counters_; }
  const BlockContext& context() const { return context_; }

  // Finds the best encoding method for the given coeffs
  WP2Status FindBestEncodingMethods(CodedBlock* const cb);

  // Decide alpha mode and encoding params if needed
  WP2Status DecideAlpha(CodedBlock* const cb);
  WP2Status RecordAlpha(const CodedBlock& cb);

  // Fills 'segment_ids_' with default values taken from 'gparams_'.
  WP2Status SetInitialSegmentIds();

  ChromaSubsampling chroma_subsampling() const { return chroma_subsampling_; }

  // Final writing of the block, except for its size.
  WP2Status WriteBlock(const CodedBlock& cb, uint32_t block_index,
                       ANSEnc* const enc);

 private:
  // Resets all records, called by InitPass.
  WP2Status ResetRecord();

  void RecordBlockHeader(const CodedBlock& cb);

  const EncoderConfig* config_;
  ANSDictionaries* dicts_;
  AlphaWriter alpha_writer_;

  // Actual number of blocks and transforms if known or upper bound.
  uint32_t num_blocks_, num_transforms_;
  // Number of times Record() has been called.
  uint32_t recorded_blocks_;
  // Pass number.
  uint32_t pass_number_ = 0u;

  Rectangle tile_rect_;  // In pixels, not padded.

  ChromaSubsampling chroma_subsampling_;

  const GlobalParams* gparams_;
  bool use_aom_coeffs_;

  SegmentIdPredictor segment_ids_;
  BlockContext context_;

  SymbolWriter symbol_writer_;
  SymbolRecorder symbol_recorder_;
  mutable Counters counters_;  // temporary data made to be shared
  ResidualWriter residual_writer_;

  // Debugging.
  // Set to true to print debug information for checking accuracy of estimated
  // rates (bit costs).
  static const bool kDebugPrintRate = false;
  // Estimated residual rate per block, per channel.
  Vector<std::array<float, 3>> residual_rate_;
  // Encodes the 'block' boundaries and the 'pixels'.
  void PutRawPixels(const Block& block, const YUVPlane& pixels,
                    ANSEnc* const enc);
};

//------------------------------------------------------------------------------

// Returns true if 'EncoderInfo::visual_debug' in the 'config' exists and
// contains the 'token'. Tokens are separated by '/'.
bool VDMatch(const EncoderConfig& config, const char token[]);
// Returns the channel corresponding to the current debug view. Asserts that
// the debug view contains the name of a channel ("y", "u", "v" or "a") in its
// path.
Channel VDChannel(const EncoderConfig& config);
// Returns true if the given rectangle is selected in visual debug.
bool VDSelected(uint32_t tile_x, uint32_t tile_y, const Rectangle& rect,
                const EncoderConfig& config);
// Returns the parameter if 'config->info.visual_debug' matches "p="+parameter
// or 'default_value' otherwise.
float VDGetParam(const EncoderConfig& config, float default_value);

//------------------------------------------------------------------------------

}  // namespace WP2

#endif  /* WP2_WP2_ENC_I_H_ */
