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
// Predictors for blocks.
//
// Authors: Skal (pascal.massimino@gmail.com)

#ifndef WP2_COMMON_LOSSY_PREDICTOR_H_
#define WP2_COMMON_LOSSY_PREDICTOR_H_

#include <memory>
#include <string>

#include "src/common/constants.h"
#include "src/common/lossy/transforms.h"
#include "src/dsp/dsp.h"
#include "src/utils/ans.h"
#include "src/utils/plane.h"
#include "src/utils/utils.h"
#include "src/utils/vector.h"
#include "src/utils/wiener.h"

namespace WP2 {

//------------------------------------------------------------------------------

// Returns a string containing the values of the pixels inside the 'block' of
// 'w' by 'h' pixels and its surrounding 'context'.
std::string GetContextAndBlockPixelsStr(const int16_t context[],
                                        const int16_t context_tr[], uint32_t w,
                                        uint32_t h, const int16_t block[],
                                        uint32_t step);

//------------------------------------------------------------------------------
// Predictors

typedef WienerOptimizer<kContextSize, kPredWidth * kPredHeight> PredOptimizer;

class CodedBlock;
class SymbolManager;
class Counters;
class SymbolReader;

// A predictor predicts a kPredWidth x kPredHeight block's values using
// heuristics and local context information.
class Predictor : public WP2Allocable {
 public:
  // Context of the prediction: each predictor can belong to a context (mainly
  // horizontal, mainly vertical ...). AV1 defines the mapping between the
  // predictor and a mode in Intra_Mode_Context.
  enum class ModeContext {
    kDC,
    kVertical,
    kHorizontal,
    kTopRight,
    kLeft,  // for left angle predictors that are not horizontal.
    kNum
  };

  virtual ~Predictor() = default;

  // Predicts values for the given block, on channel 'channel', for the
  // whole block.
  // Stores the result to 'output', incremented by 'step' for each row.
  virtual void Predict(const CodedBlock& cb, Channel channel, int16_t output[],
                       uint32_t step) const = 0;

  // Returns true if it should be processed again everytime luma changes.
  virtual bool DependsOnLuma() const { return false; }

  // Returns true if the transform is set by this predictor.
  virtual bool OverridesTransform(Channel channel,
                                  TransformPair* const tf) const {
    if (channel == kUChannel || channel == kVChannel) {
      // When overriding this function for UV, use Mode_To_Txfm from the AV1
      // specs.
      if (tf != nullptr) *tf = kDctDct;
      return true;
    }
    return false;
  }

  // The three following methods are for predictors that signal some
  // parameters in the bitstream.

  // Computes parameters needed for this transform and stores them in the block.
  // This should be called on the encoder side before using the predictor.
  virtual void ComputeParams(CodedBlock* const cb, Channel channel) const {
    (void)cb, (void)channel;
  }
  virtual void WriteParams(const CodedBlock& cb, Channel channel,
                           SymbolManager* const sm,
                           ANSEncBase* const enc) const {
    (void)cb, (void)channel, (void)sm, (void)enc;
  }
  virtual void ReadParams(CodedBlock* const cb, Channel channel,
                          SymbolReader* const sm, ANSDec* const dec) const {
    (void)cb, (void)channel, (void)sm, (void)dec;
  }
  // Setter/Getter for the mode of a predictor. For now, each predictor is its
  // own mode, and we set it to the index in the vector containing all
  // predictors.
  void SetMode(uint32_t mode) { mode_ = mode; }
  virtual uint32_t mode() const { return mode_; }
  virtual ModeContext mode_context() const { return ModeContext::kNum; }

  // Debugging methods.
  virtual std::string GetName() const = 0;
  virtual std::string GetInfoStr() const { return ""; }
  virtual std::string GetFakePredStr() const = 0;
  virtual std::string GetPredStr(const CodedBlock& cb,
                                 Channel channel) const = 0;
  virtual void Draw(const WP2::Rectangle& rect,
                    ArgbBuffer* const debug_output) const = 0;

  static constexpr uint32_t kInvalidMode = std::numeric_limits<uint32_t>::max();

 private:
  uint32_t mode_ = kInvalidMode;
};

// This vector will *own* the Predictor pointers and destruct them
// when going out of scope.
struct PredictorVector : public Vector<Predictor*> {
  PredictorVector() = default;
  PredictorVector(PredictorVector&&) = default;
  ~PredictorVector() { reset(); }
  void reset();
};

//------------------------------------------------------------------------------
// Predictors.

// Class containing the different predictors.
struct EncoderConfig;
class BlockContext;
class Segment;
class SymbolRecorder;
class Predictors {
 public:
  Predictors() = default;
  Predictors(Predictors&&) = default;
  void reset() {
    preds_.reset();
    preds_no_angle_.reset();
  }
  bool empty() const { return preds_.empty(); }
  uint32_t size() const { return preds_.size(); }
  const Predictor* operator[](size_t i) const { return preds_[i]; }
  PredictorVector::const_iterator begin() const { return preds_.begin(); }
  PredictorVector::const_iterator end() const { return preds_.end(); }
  // Finds the best predictor for a block 'cb'. This predictor has to give a
  // score better than the current one 'best_score'.
  // This is speed dependent: it looks at all non-angle predictors, and for
  // angle predictors, it looks at the main angles but only refine some of them
  // depending on the speed.
  // For Y/A channel, best_score needs to be given and the function will try to
  // find a better score. If found, best_factor will be replaced and
  // best_predictor set.
  WP2Status FindBest(const EncoderConfig& config, const Rectangle& tile_rect,
                     uint32_t num_channels, const BlockContext& context,
                     const Segment& segment, bool reduced,
                     const Vector<TransformPair>& transforms,
                     Counters* const counters, CodedBlock* const cb,
                     float* const best_score,
                     const Predictor** const best_predictor,
                     TransformPair* const best_tf) const;
  // Get the index of the firt predictor with a given mode.
  const Predictor* GetFirstWithMode(uint32_t mode) const;
  const Predictor* GetWithMode(uint32_t mode, uint32_t sub_mode) const;
  // Get the maximum mode contained in a predictor vector.
  uint32_t GetMaxMode() const;

 protected:
  enum class Pred {
    kDcAll,
    kDcLeft,
    kDcTop,
    kSmooth2D,
    kSmoothVertical,
    kSmoothHorizontal,
    kTrueMotion,
    kGradient,
    kAngle45,
    kAngle67,
    kAngle90,
    kAngle113,
    kAngle135,
    kAngle157,
    kAngle180,
    kAngle203,
    kCfl,
    kSignalingCfl,
    kZero,
    kNum,
  };

  // Fills the different angle and non-angle predictors.
  WP2Status FillImpl(const Pred* const mapping, uint32_t num_modes,
                     int16_t min_value, int16_t max_value,
                     const CSPTransform* const transform);

  PredictorVector preds_;
  Vector<const Predictor*> preds_no_angle_;
  std::array<const Predictor*, kAnglePredNum> preds_main_angle_;
  Channel channel_;

 private:
  // At each index, contains a pointer to the first Predictor with
  // mode == index.
  std::array<const Predictor*, kYBasePredNum + kAnglePredNum> first_with_mode_;
};

class YAPredictors : public Predictors {
 public:
  YAPredictors() = default;
};

class YPredictors : public YAPredictors {
 public:
  YPredictors();
  YPredictors(YPredictors&&) = default;
  // Fills the object with the non-angle and angle predictors.
  WP2Status Fill(int16_t min_value, int16_t max_value);
};

class APredictors : public YAPredictors {
 public:
  APredictors();
  APredictors(APredictors&&) = default;
  // Fills the object with the non-angle and angle predictors.
  WP2Status Fill(int16_t min_value, int16_t max_value,
                 const CSPTransform* const transform);
};

class UVPredictors : public Predictors {
 public:
  UVPredictors();
  UVPredictors(UVPredictors&&) = default;
  // Fills the object with the non-angle and angle predictors.
  WP2Status Fill(int16_t min_value, int16_t max_value);
};

//------------------------------------------------------------------------------

// Context-based predictor base class.
class ContextPredictor : public Predictor {
 public:
  void Predict(const CodedBlock& cb, Channel channel, int16_t output[],
               uint32_t step) const override;

  std::string GetFakePredStr() const override;
  std::string GetPredStr(const CodedBlock& cb, Channel channel) const override;

 protected:
  virtual void Predict(const int16_t context[], uint32_t w, uint32_t h,
                       int16_t output[], uint32_t step) const = 0;
  std::string GetPredStr(const int16_t context[], const int16_t context_tr[],
                         bool extend_right, uint32_t width,
                         uint32_t height) const;
  bool extend_context_ = false;  // whether we extend the top context right
};

//------------------------------------------------------------------------------
// CFL predictors.

// Chroma-form-luma predictor. Predicts chroma as being:
// chroma = a * luma + b
// where the a and b parameters are entirely deduced from context pixels.
// Also used for alpha.
class CflPredictor : public Predictor {
 public:
  CflPredictor(int16_t min_value, int16_t max_value);

  void Predict(const CodedBlock& cb, Channel channel, int16_t output[],
               uint32_t step) const override;

  bool DependsOnLuma() const override { return true; }

  std::string GetName() const override;
  std::string GetFakePredStr() const override;
  std::string GetPredStr(const CodedBlock& cb, Channel channel) const override;
  void Draw(const WP2::Rectangle& rect,
            ArgbBuffer* const debug_output) const override {}

  // Does a simple linear regression based on the context for chroma and luma,
  // assuming chroma = a * luma + b
  void ContextLinearRegression(Channel channel, const CodedBlock& cb,
                               int16_t* a, int16_t* b) const;
  // Does a simple linear regression based on the reconstructed luma and input
  // chroma, assuming chroma = a * luma + b
  void BlockLinearRegression(Channel channel, const CodedBlock& cb, int16_t* a,
                             int16_t* b) const;

  static constexpr uint32_t kAPrecShift = 8;  // In bits.
  static constexpr uint32_t kBPrecShift = 2;

 protected:
  virtual void DoPredict(const CodedBlock& cb, Channel channel,
                         uint32_t x_start, uint32_t y_start, uint32_t w,
                         uint32_t h, int16_t* output, uint32_t step) const;

  const int16_t min_value_;
  const int16_t max_value_;
};

// Like CflPredictor, but signals in the bitstream the difference between
// the predicted slope ('a' in 'chroma = a * chroma + b') and the actual one.
// More exact but obviously more costly.
class SignalingCflPredictor : public CflPredictor {
 public:
  SignalingCflPredictor(int16_t min_value, int16_t max_value);

  std::string GetName() const override;
  std::string GetPredStr(const CodedBlock& cb, Channel channel) const override;

  void ComputeParams(CodedBlock* const cb, Channel channel) const override;
  void WriteParams(const CodedBlock& cb, Channel channel,
                   SymbolManager* const sm,
                   ANSEncBase* const enc) const override;
  void ReadParams(CodedBlock* const cb, Channel channel, SymbolReader* const sr,
                  ANSDec* const dec) const override;

  bool DependsOnLuma() const override { return true; }

  constexpr static uint32_t kAResBits = 10;

 protected:
  constexpr static uint32_t kAResPrecision = 6;

  void DoPredict(const CodedBlock& cb, Channel channel, uint32_t unused_x,
                 uint32_t unused_y, uint32_t w, uint32_t h, int16_t* output,
                 uint32_t step) const override;
  void GetParams(const CodedBlock& cb, Channel channel, uint32_t w, uint32_t h,
                 int32_t* a, int32_t* b) const;

 private:
  int16_t GetParam(const CodedBlock& cb, Channel channel) const;
};

//------------------------------------------------------------------------------

// Adds predictors for the alpha plane.
WP2Status InitAlphaPredictors(const CSPTransform& transform,
                              APredictors* const preds);
// Total number of angle predictors.
static constexpr uint32_t kDirectionalPredNumYA =
    kAnglePredNum * (2 * kDirectionalMaxAngleDeltaYA + 1);
static constexpr uint32_t kDirectionalPredNumUV =
    kAnglePredNum * (2 * kDirectionalMaxAngleDeltaUV + 1);

// Number of prediction modes (8 'base' predictors + directional).
static constexpr uint32_t kYPredModeNum = kYBasePredNum + kAnglePredNum;
static constexpr uint32_t kAPredModeNum = kABasePredNum + kAnglePredNum;
static constexpr uint32_t kUVPredModeNum = kUVBasePredNum + kAnglePredNum;
// Number of predictors (number of modes and their sub-modes).
static constexpr uint32_t kYPredNum = kYBasePredNum + kDirectionalPredNumYA;
static constexpr uint32_t kAPredNum = kABasePredNum + kDirectionalPredNumYA;
static constexpr uint32_t kUVPredNum = kUVBasePredNum + kDirectionalPredNumUV;

}  // namespace WP2

#endif  // WP2_COMMON_LOSSY_PREDICTOR_H_
