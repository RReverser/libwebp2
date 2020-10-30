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
//   All about blocks
//
// Author: Skal (pascal.massimino@gmail.com)

#ifndef WP2_COMMON_LOSSY_BLOCK_H_
#define WP2_COMMON_LOSSY_BLOCK_H_

#include "src/common/constants.h"
#include "src/common/global_params.h"
#include "src/common/lossy/block_size.h"
#include "src/common/lossy/context.h"
#include "src/common/lossy/predictor.h"
#include "src/common/lossy/quant_mtx.h"
#include "src/common/lossy/residuals.h"
#include "src/common/lossy/residuals_aom.h"
#include "src/common/lossy/transforms.h"
#include "src/dsp/dsp.h"
#include "src/utils/plane.h"
#include "src/utils/vector.h"
#include "src/wp2/decode.h"
#include "src/wp2/encode.h"
#include "src/wp2/format_constants.h"

namespace WP2 {

class FrontMgrBase;
class SymbolReader;
class SymbolCounter;
class SymbolRecorder;
class SymbolWriter;
class SymbolManager;
class Segment;
class RndMtxSet;

//------------------------------------------------------------------------------

// Class containing several SymbolCounters for re-use (to avoid too many
// mallocs).
class Counters : public WP2Allocable {
 public:
  WP2Status Init(const SymbolRecorder& recorder);
  WP2Status CopyFrom(const Counters& other, const SymbolRecorder& recorder);
  const SymbolRecorder* recorder() const { return recorder_; }
  // Returns a pointer to the counter with symbols related to transforms.
  SymbolCounter* transform() { return transform_.get(); }
  // Returns a pointer to the counter with symbols related to residuals.
  SymbolCounter* residuals() { return residuals_.get(); }
  // Returns a pointer to the counter with symbols related to AOM residual.
  SymbolCounter* residuals_aom() { return residuals_aom_.get(); }
  // Returns a pointer to the counter with symbols related to predictors.
  SymbolCounter* predictor() { return predictor_.get(); }

 private:
  std::unique_ptr<SymbolCounter> transform_;
  std::unique_ptr<SymbolCounter> residuals_;
  std::unique_ptr<SymbolCounter> residuals_aom_;
  std::unique_ptr<SymbolCounter> predictor_;
  const SymbolRecorder* recorder_ = nullptr;
};

//------------------------------------------------------------------------------
// Context passed around and updated once we know the right parameters for each
// block. Only contains AOM context for now.

class BlockContext {
 public:
  // 'use_aom' defines the usage of AOM residuals. 'width' and 'height' are in
  // pixels.
  WP2Status Init(bool use_aom, uint32_t width, uint32_t height) {
    use_aom_ = use_aom;
    if (use_aom) WP2_CHECK_STATUS(aom_.Init(width, height));
    WP2_CHECK_STATUS(modes_.InitMap(width, height));
    return WP2_STATUS_OK;
  }
  WP2Status CopyFrom(const BlockContext& other) {
    use_aom_ = other.use_aom_;
    if (use_aom_) WP2_CHECK_STATUS(aom_.CopyFrom(other.aom_));
    WP2_CHECK_STATUS(modes_.CopyFrom(other.modes_));
    return WP2_STATUS_OK;
  }
  void Reset() {
    if (use_aom_) aom_.Reset();
    modes_.Reset();
  }
  // Returns a pointer to the AOM context, or nullptr if AOM residuals are not
  // used.
  const libgav1::AOMContext* aom() const { return use_aom_ ? &aom_ : nullptr; }
  libgav1::AOMContext* aom() { return use_aom_ ? &aom_ : nullptr; }
  bool use_aom() const { return use_aom_; }
  const YModePredictor& ymodes() const { return modes_; }
  YModePredictor& ymodes() { return modes_; }

 private:
  bool use_aom_ = false;
  libgav1::AOMContext aom_;
  YModePredictor modes_;
};

//------------------------------------------------------------------------------
// Context around a block (full context or 4x4 sub context).

// Class that caches contexts for speed-up.
class ContextCache {
 public:
  // Gets the large context for a given block.
  // TODO(vrabaud) fill_in == false and extend_right == true seems unused.
  const int16_t* GetLarge(const CodedBlock& cb, Channel channel, bool fill_in,
                          bool extend_right);
  // Sets contexts for all parameters as uncomputed.
  void Reset();

 private:
  // Large context per channel, per fill_in, per extend right.
  int16_t large_[4][2][2][kMaxContextSize];
  bool is_computed_[4][2][2] = {{{false, false}, {false, false}},
                                {{false, false}, {false, false}},
                                {{false, false}, {false, false}},
                                {{false, false}, {false, false}}};
};

//------------------------------------------------------------------------------
// Struct to reduce memory footprint. A block can have at most 4 transforms and
// each transform may contain kMaxBlockSizePix2 coeffs but there are no more
// than kMaxBlockSizePix2 per block too. Allocate only what is necessary while
// keeping the convenient index access.

static constexpr uint32_t kMaxNumTransformsPerBlock = 4;
// Splits 1 and 2 are swapped to store half in [0] and half in [1].
static constexpr uint32_t kBlockCoeffOffset[] = {
    0, kMaxBlockSizePix2 / 2, kMaxBlockSizePix2 / 4, kMaxBlockSizePix2 * 3 / 4};
STATIC_ASSERT_ARRAY_SIZE(kBlockCoeffOffset, kMaxNumTransformsPerBlock);

template <typename T>
struct BlockCoeffs {
  T* operator[](uint32_t index) { return &memory[kBlockCoeffOffset[index]]; }
  const T* operator[](uint32_t index) const {
    return &memory[kBlockCoeffOffset[index]];
  }

 private:
  T memory[kMaxBlockSizePix2];  // Internal storage.
};
typedef BlockCoeffs<int16_t> BlockCoeffs16;
typedef BlockCoeffs<int32_t> BlockCoeffs32;

//------------------------------------------------------------------------------
// CodedBlock

class CodedBlock {
 public:
  CodedBlock() = default;
  CodedBlock(const CodedBlock&) = delete;
  CodedBlock(CodedBlock&&) = default;
  virtual ~CodedBlock() = default;

  // Sets range of legal YUV values (inclusive, used for reconstruction) to
  // ['yuv_min', 'yuv_max'].
  void SetRange(int16_t yuv_min, int16_t yuv_max);

  // Set the views on the source and destination buffers based on 'blk_'.
  // Position and size must be set prior to calling these.
  void SetSrcInput(const YUVPlane& in);  // Used by GetDisto().
  // Used by GetContext(). If context_cache == nullptr, the previous instance
  // will be used. It must be != nullptr the first time it is used.
  void SetContextInput(const YUVPlane& in,
                       ContextCache* const context_cache = nullptr);
  void SetReconstructedOutput(YUVPlane* const out);  // Used by Optimize*() etc.

  // Sets the U/V predictor for the whole block.
  void SetUVPredictor(const Predictor* const pred);
  // Sets the predictor for a given 'channel' of the whole block.
  void SetPredictor(Channel channel, const Predictor* const pred);
  // To be called when LumaContextIsUniform().
  void SetLumaUniformPredictor(const Predictors& preds);

  // Assigns the best of the 'transforms' (based on full pred+tf+quant+recons
  // score) to the CodingParams of 'channel'. The final predictor must be
  // selected before calling FindBestTf().
  WP2Status FindBestTf(const EncoderConfig& config, const Rectangle& tile_rect,
                       const QuantMtx& quant, const BlockContext& context,
                       Channel channel, uint32_t num_channels, bool reduced,
                       const Vector<TransformPair>& transforms,
                       Counters* const counters);
  // Assigns the best of the 'preds' (based on pred score only) to the
  // CodingParams of 'channel'.
  void FindBestPred(const Predictors& preds, Channel channel, bool reduced);
  // Assigns the best pair of 'preds/transforms' (based on full
  // pred+tf+quant+recons score) to the CodingParams of 'channel'.
  WP2Status FindBestPredTf(const EncoderConfig& config,
                           const Rectangle& tile_rect, const Predictors& preds,
                           const Segment& segment, const BlockContext& context,
                           Channel channel, uint32_t num_channels, bool reduced,
                           const Vector<TransformPair>& transforms,
                           Counters* const counters);
  // Selects whether to 'split_tf' or not based on residuals rate.
  WP2Status FindBestSplitTf(const EncoderConfig& config,
                            const Rectangle& tile_rect, const QuantMtx& quant,
                            const BlockContext& context, Channel channel,
                            uint32_t num_channels, bool reduced,
                            Counters* const counters);
  // Calls FindBestPred/Tf() depending on 'config.speed' and 'context' then
  // Quantize(). All 'transforms' will be tested but only 'transforms_subset'
  // will be tried for each of the 'preds'.
  WP2Status OptimizeModesLuma(const EncoderConfig& config,
                              const Rectangle& tile_rect, bool has_alpha,
                              const Predictors& preds, const Segment& segment,
                              const BlockContext& context,
                              const Vector<TransformPair>& transforms,
                              const Vector<TransformPair>& transforms_subset,
                              Counters* const counters);

  // Decides whether to use yuv420 based on coeffs_[kUChannel or kVChannel].
  WP2Status DecideChromaSubsampling(const EncoderConfig& config,
                                    uint32_t tile_pos_x, uint32_t tile_pos_y,
                                    bool has_alpha, const QuantMtx& quant_u,
                                    const QuantMtx& quant_v,
                                    Counters* const counters);

  // Returns true if the context of the block (neighboring pixels used for
  // prediction) all have the same value for the given channel.
  bool ContextIsConstant(Channel channel) const;

  // Calls FindBestUVModes() and DecideChromaSubsampling() then Quantize().
  WP2Status OptimizeModesChroma(
    const EncoderConfig& config, const Rectangle& tile_rect, bool has_alpha,
    const FrontMgrBase& mgr, const UVPredictors& preds, const Segment& segment,
    const BlockContext& context, ChromaSubsampling chroma_subsampling,
    DCDiffusionMap* const dc_error_u, DCDiffusionMap* const dc_error_v,
    Counters* const counters);

  // Computes the prediction to 'out'.
  void PredictBlock(Channel channel, int16_t* out, uint32_t step) const;

  // Subtracts the 'prediction' from 'in_' and stores the residuals in 'res'.
  void GetResiduals(Channel channel, const Plane16& prediction,
                    BlockCoeffs32* const res) const;

  // Applies the TransformPair to 'res', quantizes the 'coeffs_' and
  // reconstructs the pixels to 'out_' using the 'prediction'.
  void TransformAndReconstruct(const QuantMtx& quant, Channel channel,
                               bool reduced_transform,
                               const Plane16& prediction,
                               BlockCoeffs32* const res);

  // Quantizes and reconstructs channel 'channel'. If 'reduced_transform' is
  // true, 2x downsampling of residuals is performed.
  void Quantize(const QuantMtx& quant, Channel channel, bool reduced_transform);

  // Dequantizes coeffs_[] into res[]. dim_, use_mtx_ and is420_ must be set.
  void Dequantize(const QuantMtx& quant, Channel channel,
                  BlockCoeffs32* const res);

  // Reconstructs 'channel' from dequantized residuals 'res' and 'prediction'
  // to 'out_'. 'prediction' can also point to 'out_' or be empty.
  // Upon return, res[] contains the inverse-transformed coeffs.
  void Reconstruct(Channel channel, bool reduced_transform,
                   BlockCoeffs32* const res,
                   const Plane16& prediction = Plane16()) const;

  // Visual debug: stores transform at block position
  void StoreTransform(const WP2::DecoderConfig& config, uint32_t tile_x,
                      uint32_t tile_y, ArgbBuffer* const debug_output) const;

  // Visual debug: stores residuals at block position
  void StoreResiduals(const DecoderConfig& config, uint32_t tile_x,
                      uint32_t tile_y, const QuantMtx& quant, Channel channel,
                      Plane16* const dst_plane) const;
  // Visual debug: stores original residuals (before quantization) at block
  // position.
  void StoreOriginalResiduals(
      const EncoderConfig& config, uint32_t tile_pos_x, uint32_t tile_pos_y,
      int32_t original_res[kMaxNumTransformsPerBlock][kMaxBlockSizePix2],
      Plane16* const dst_plane) const;
  // Visual debug: stores prediction scores
  void StorePredictionScore(const EncoderConfig& config,
                            const Rectangle& tile_rect, Channel channel,
                            const Predictor& pred, TransformPair tf, float dist,
                            float lambda, float res_rate, float pred_rate,
                            float tf_rate, float score, bool is_best) const;
  // Visual debug: stores prediction modes
  void StorePredictionModes(const WP2::DecoderConfig& config,
                            const Rectangle& tile_rect, Channel channel,
                            const Predictors& preds,
                            Plane16* const raw_prediction,
                            ArgbBuffer* const debug_output) const;
  // Visual debug: append the values of the original pixels as text
  void AppendOriginalPixels(const DecoderConfig& config, uint32_t tile_x,
                            uint32_t tile_y, const CSPTransform& csp_tranform,
                            ArgbBuffer* const debug_output) const;
  // Visual debug: append the values of the original pixels as text
  void AppendCompressedPixels(const DecoderConfig& config, uint32_t tile_x,
                              uint32_t tile_y,
                              ArgbBuffer* const debug_output) const;
  // Visual debug: stores coeff methods
  void StoreCoeffMethod(const WP2::DecoderConfig& config,
                        Plane16* const dst_plane) const;
  // Visual debug: stores the best possible prediction of the given channel,
  // assuming we could transmit the best slope/bias instead of deducing them
  // from context.
  void StoreBestCflPrediction(Channel channel, int16_t yuv_min, int16_t yuv_max,
                              Plane16* const dst_plane,
                              std::string* const debug_str) const;
  // Visual debug: stores raw (unquantized) residuals after ideal chroma from
  // luma prediction as above.
  void StoreBestCflResiduals(Channel channel, int16_t yuv_min, int16_t yuv_max,
                             Plane16* const dst_plane,
                             std::string* const debug_str) const;
  // Visual debug: stores the best slope for Chroma from luma ('a' in 'chroma =
  // a * luma + b'), after doing a linear regression between the output luma and
  // input chroma.
  void StoreBestCflSlope(Channel channel, int16_t yuv_min, int16_t yuv_max,
                         Plane16* const dst_plane,
                         std::string* const debug_str) const;
  // Visual debug: stores the best intercept for Chroma from luma ('b' in
  // 'chroma = a * luma + b'), after doing a linear regression between the
  // output luma and input chroma.
  void StoreBestCflIntercept(Channel channel, int16_t yuv_min, int16_t yuv_max,
                             Plane16* const dst_plane,
                             std::string* const debug_str) const;
  // Visual debug: stores the slope for Chroma from luma ('a' in 'chroma =
  // a * luma + b'), after doing a linear regression on the context.
  void StoreCflSlope(Channel channel, int16_t yuv_min, int16_t yuv_max,
                     Plane16* const dst_plane,
                     std::string* const debug_str) const;
  // Visual debug: stores the intercept for Chroma from luma ('b' in 'chroma =
  // a * luma + b'), after doing a linear regression on the context.
  void StoreCflIntercept(Channel channel, int16_t yuv_min, int16_t yuv_max,
                         Plane16* const dst_plane,
                         std::string* const debug_str) const;
  // Visual debug: stores information around the chroma subsampling decision.
  void Store420Scores(const EncoderConfig& config, uint32_t pos_x,
                      uint32_t pos_y, float lambda_u, float lambda_v,
                      bool reduced, uint32_t disto,
                      float rate_u, float rate_v);
  enum class Debug420Decision { k444EarlyExit, k444, k420 };
  WP2Status Store420Decision(const EncoderConfig& config, uint32_t pos_x,
                             uint32_t pos_y, Debug420Decision decision) const;

  // Visual debug: store lambda multiplier map.
  WP2Status StoreLambdaMult(const EncoderConfig& config,
                            uint32_t pos_x, uint32_t pos_y) const;

  // Visual debug: stores propagated error or new error for U/V error diffusion.
  void StoreErrorDiffusion(const EncoderConfig& config, uint32_t tile_x,
                           uint32_t tile_y, Plane16* const dst_plane) const;

  // VisualDebug
  WP2Status Draw(const DecoderConfig& config, uint32_t tile_x, uint32_t tile_y,
                 const GlobalParams& gparams,
                 ArgbBuffer* const debug_output) const;

  // Exports debug info during decoding
  void ToBlockInfo(struct BlockInfo* const blk) const;

  static constexpr int16_t kMissing = std::numeric_limits<int16_t>::max();

  // Retrieves boundary samples for the whole block.
  // If 'extend_right' is set to true, the context is extended horizontally to
  // the right instead of along the right boundary.
  // If 'fill_in' is true, missing context pixels are filled in
  // with neighboring values. Otherwise, they are set to kMissing.
  virtual const int16_t* GetContext(Channel channel, bool fill_in,
                                    bool extend_right) const;
  void SetContextCache(ContextCache* const context_cache);
  void ResetContextCache() const;

  // Returns left and right occupancy, in kMinBlockSizePix units.
  // If top_context_extent is not nullptr, fills it with the possible extent of
  // the top context.
  void GetOccupancy(int8_t* const left, int8_t* const right,
                    int8_t* const top_context_extent) const;

  // Coding parameters for a given channel.
  struct CodingParams {
    // Predictor.
    const Predictor* pred = nullptr;
    // Some predictors (e.g. angle predictors) have sub-modes.
    uint8_t pred_sub_mode = 0;
    // Frequency domain transform.
    TransformPair tf = kDctDct;
    // Whether to split the block into smaller transforms or not.
    bool split_tf = false;

    // Returns the frequency domain transform for the horizontal axis.
    WP2TransformType tf_x() const { return kTfX[tf]; }
    // Returns the frequency domain transform for the vertical axis.
    WP2TransformType tf_y() const { return kTfY[tf]; }
  };

  const CodingParams& GetCodingParams(Channel channel) const;
  CodingParams* GetCodingParams(Channel channel);

  // Returns the transform if it is predetermined (e.g. by the predictor) and
  // does not need to be signaled. Otherwise returns 'kUnknownTf'.
  TransformPair GetImplicitTf(Channel channel) const;

  // Returns the sum of square error between 'p_in_' and 'p_out_' for the pixels
  // inside 'tile_rect' (px, not padded).
  uint32_t GetDisto(Channel channel, const Rectangle& tile_rect) const;

  // Returns the number of sub-transforms.
  uint32_t GetNumTransforms(Channel channel) const;

  // Returns the number of coefficients per split transform, taking is420 into
  // account.
  uint32_t NumCoeffsPerTransform(Channel channel) const;

  // Returns true if the transform implies that the first coeff is DC.
  bool IsFirstCoeffDC(Channel channel) const;

  // Returns true if there is at least one non zero coeff in the block.
  bool HasCoeffs(Channel channel) const;

  // Accessors and setters (padded).
  uint32_t x() const { return blk_.x(); }
  uint32_t x_pix() const { return blk_.x_pix(); }
  uint32_t y() const { return blk_.y(); }
  uint32_t y_pix() const { return blk_.y_pix(); }
  uint32_t w() const { return blk_.w(); }
  uint32_t w_pix() const { return blk_.w_pix(); }
  uint32_t h() const { return blk_.h(); }
  uint32_t h_pix() const { return blk_.h_pix(); }

  // Non-padded dimensions deduced from 'tile_rect' (px, not padded).
  uint32_t visible_w_pix(const Rectangle& tile_rect) const;
  uint32_t visible_h_pix(const Rectangle& tile_rect) const;

  BlockSize dim() const { return blk_.dim(); }
  // Returns the size of the transform but taking into account the potential
  // size reduction due to is420_ with U/V.
  TrfSize tdim(Channel channel) const {
    return GetTransform(
        GetSplitSize(dim(), GetCodingParams(channel).split_tf),
        (channel == kUChannel || channel == kVChannel) && is420_);
  }
  Block blk() const { return blk_; }
  // Sets position/dimension based on 'block', occupancy based on 'mgr'.
  void SetDim(const Block& block, const FrontMgrBase& mgr);
  // When no mgr is provided, left block will be considered available and the
  // right block not available.
  void SetDimDefault(const Block& block, bool full_left_ctx = false);
  void SetXY(const FrontMgrBase& mgr, uint32_t x, uint32_t y);
  // Returns the plane we use to compute the context.
  const Plane16& GetContextPlane(Channel channel) const {
    return predict_from_.GetChannel(channel);
  }

  // Pixel transfer. 'yuv' is the full plane, not the block's sub-view.
  // No bound-check or clipping is performed.
  void ExtractFrom(const YUVPlane& yuv, Channel channel) const;

  // Returns an estimation of the number of bits necessary to encode the
  // predictor modes.
  WP2Status PredictorRate(const YModePredictor* const ymode_predictor,
                          Channel channel, SymbolCounter* const counter,
                          float* const rate) const;

  // Returns an estimation of the number of bits necessary to encode the
  // transform.
  WP2Status TransformRate(Channel channel, uint32_t num_transforms,
                          SymbolCounter* const counter,
                          float* const rate) const;
  WP2Status SplitTransformRate(Channel channel, SymbolCounter* const counter,
                               float* const rate) const;

  // Estimated cost in bits of storing residuals for the given channel.
  WP2Status ResidualRate(const BlockContext& context, Channel channel,
                         uint32_t num_channels, SymbolCounter* const counter,
                         SymbolCounter* const aom_counter, float* rate) const;
  // Gets the different distortion/rates of a block as well as the final score.
  WP2Status GetRates(const Rectangle& tile_rect, const BlockContext& context,
                     const QuantMtx& quant, Channel channel,
                     uint32_t num_channels, float tf_rate,
                     Counters* const counters, uint32_t* const dist,
                     float* const res_rate, float* const pred_rate,
                     float* const score) const;

  // Returns true if this block has alpha residuals.
  bool HasLossyAlpha() const;

  friend ContextCache;

 private:
  void CflLinearRegression(Channel channel, int16_t yuv_min, int16_t yuv_max,
                           float* const a, float* const b,
                           std::string* const debug_str) const;

  void GetContext(Channel channel, bool fill_in, bool extend_right,
                  int16_t context[kMaxContextSize]) const;

 public:
  uint32_t left_occupancy_;   // occupancy along the left edge
  uint32_t right_occupancy_;  // occupancy along the right edge
  // Length of the top context beyond the top-right corner of the block.
  uint32_t top_context_extent_;
  uint8_t id_;      // segment id
  bool is420_;      // reduced transform

  uint8_t mtx_[2];        // index of rnd_mtx to use
  bool use_mtx_ = false;  // true if rnd_mtx should be used
  const RndMtxSet* mtx_set_ = nullptr;
  template<class T> bool CheckUseMtx(BlockSize dim, const T& mtx);

  // Input and output samples must be views of buffers external to CodedBlock
  // pointing to the rectangle covered by 'blk_'. Except for 'predict_from_',
  // they will not be accessed outside of these bounds.
  YUVPlane in_;            // source samples (input), within blk_
  mutable YUVPlane out_;   // coded samples (output), within blk_

  // Quantized residuals for Y/U/V/A, aka coefficients.
  BlockCoeffs16 coeffs_[4];
  // Size of the coefficients for Y/U/V/A per split transform.
  // If size == 0, all coefficients are 0. If size != 0, the coefficient at
  // size - 1 is the last non-zero one.
  uint32_t num_coeffs_[4][kMaxNumTransformsPerBlock];
  BlockAlphaMode alpha_mode_ = kBlockAlphaLossless;
  // Approximate number of bytes taken up by lossless alpha in this block.
  size_t alpha_lossless_bytes_ = 0;
  bool y_context_is_constant_ = false;  // only used for encoding
  float pred_scores_[4] = {0};  // per channel, only available during encoding

  // Encoding method for the Y, U, V and A residuals per split transform.
  EncodingMethod method_[4][kMaxNumTransformsPerBlock];

  // Diffusion errors, indexed by Y/U/V channels.
  int16_t dc_error_[3] = {0};
  int16_t dc_error_next_[3] = {0};

  // Min/max YUV values, used to clip output when reconstructing.
  int16_t yuv_min_;
  int16_t yuv_max_;

  // Chroma from luma parameters for U, V, A.
  int16_t cfl_[3] = {0};

  // Per-block lambda modulation. Larger multiplier will favor rate over
  // distortion for the block: should be assigned to highly textured blocks.
  float lambda_mult_ = 1.f;

#if defined(WP2_BITTRACE)  // extra store for debugging
  // Original uncompressed residuals in spatial domain.
  int32_t (*original_res_)[4][kMaxNumTransformsPerBlock][kMaxBlockSizePix2] =
      nullptr;
#endif

 private:
  ContextCache* context_cache_ = nullptr;
  YUVPlane predict_from_;  // context (input), within 1px left/top/right of blk_

  CodingParams params_[3];  // Y/UV/A
  Block blk_;      // position and dimension
};

}    // namespace WP2

#endif  /* WP2_COMMON_LOSSY_BLOCK_H_ */
