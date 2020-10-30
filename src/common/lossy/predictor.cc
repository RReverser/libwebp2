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
//  Predictors
//
// Author: Skal (pascal.massimino@gmail.com)
//

#include "src/common/lossy/predictor.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <limits>
#include <numeric>

#include "src/common/constants.h"
#include "src/common/lossy/block.h"
#include "src/dec/symbols_dec.h"
#include "src/dsp/dsp.h"
#include "src/dsp/math.h"
#include "src/enc/symbols_enc.h"
#include "src/utils/utils.h"
#include "src/utils/vector.h"

#define LOG 0

namespace WP2 {

//------------------------------------------------------------------------------
// PredictorVector

void PredictorVector::reset() {
  for (auto& p : *this) delete p;  // we own the predictors
  Vector<Predictor*>::reset();
}

//------------------------------------------------------------------------------
// For debugging only.

// Appends a formatted literal to a std::string.
#define WP2SAppend(str_ptr, ...)                                         \
  do {                                                                   \
    const size_t size = std::snprintf(nullptr, 0, __VA_ARGS__) + 1;      \
    (str_ptr)->resize((str_ptr)->size() + size);                         \
    std::snprintf((char*)&(str_ptr)->at((str_ptr)->size() - size), size, \
                  __VA_ARGS__);                                          \
    (str_ptr)->pop_back(); /* Remove ending '\0' */                      \
  } while (false)

std::string GetContextAndBlockPixelsStr(const int16_t context[],
                                        const int16_t context_tr[], uint32_t w,
                                        uint32_t h, const int16_t block[],
                                        uint32_t step) {
  std::string str;
  WP2SAppend(&str, "%4d | ", context[h]);
  for (uint32_t i = 0; i < w; ++i) WP2SAppend(&str, "%4d ", context[h + 1 + i]);
  WP2SAppend(&str, "| %4d", context[h + 1 + w]);
  uint32_t num_extra_dashes = 0;
  const uint32_t tr_size = ContextWithTRSize(w, h);
  uint32_t tr_i = h + 1 + w + 1;
  if (context_tr != nullptr) {
    str += " ";
    // Display the first extended right context values.
    for (uint32_t x = 0; x < 7 && tr_i < tr_size;
         ++x, ++tr_i, ++num_extra_dashes) {
      WP2SAppend(&str, "%4d ", context_tr[tr_i]);
    }
  }
  str += "\n-----+-";
  for (uint32_t i = 0; i < w; ++i) str += "-----";
  str += "+-----";
  for (uint32_t i = 0; i < num_extra_dashes; ++i) str += "-----";
  str += "\n";
  for (uint32_t y = 0; y < h || (context_tr != nullptr && tr_i < tr_size);
       ++y) {
    if (y < h) {
      WP2SAppend(&str, "%4d | ", context[h - 1 - y]);
      for (uint32_t x = 0; x < w; ++x) {
        WP2SAppend(&str, "%4d ", block[x + y * step]);
      }
      WP2SAppend(&str, "| %4d", context[h + 1 + w + 1 + y]);
    } else {
      WP2SAppend(&str, "%*s", 7 + w * 5 + 6, " ");
    }
    // Display remaining extended right context in small batches on the right.
    if (context_tr != nullptr) {
      if (tr_i < tr_size) {
        WP2SAppend(&str, "               ");
        for (uint32_t x = 0; x < 4 && tr_i < tr_size; ++x, ++tr_i) {
          WP2SAppend(&str, " %4d", context_tr[tr_i]);
        }
      }
    }
    WP2SAppend(&str, "\n");
  }
  return str;
}

const int16_t kFakeContext[kContextSize] = {
    2,  4,  6,  8,   // left
    10,              // top left
    12, 14, 16, 18,  // top
    20,              // top-right
    22, 24, 26, 28   // right
};

//------------------------------------------------------------------------------

// Dummy predictor that always predicts zero.
class ZeroPredictor : public Predictor {
 public:
  void Predict(const CodedBlock& cb, Channel, int16_t output[],
               uint32_t step) const override {
    for (uint32_t y = 0; y < cb.h_pix(); ++y) {
      std::fill(output + y * step, output + y * step + cb.w_pix(), 0);
    }
  }

  std::string GetName() const override { return "zero predictor"; }
  std::string GetFakePredStr() const override { return "zero predictor"; }
  std::string GetPredStr(const CodedBlock& cb, Channel channel) const override {
    return "zero predictor";
  }
  void Draw(const WP2::Rectangle& rect,
            ArgbBuffer* const debug_output) const override {}
};

//------------------------------------------------------------------------------

void ContextPredictor::Predict(const CodedBlock& cb, Channel channel,
                               int16_t output[], uint32_t step) const {
  const int16_t* const context = cb.GetContext(
      channel, /*fill_in=*/true, /*extend_right=*/extend_context_);
  Predict(context, cb.w_pix(), cb.h_pix(), output, step);
}

std::string ContextPredictor::GetFakePredStr() const {
  return GetPredStr(kFakeContext, kFakeContext, /*extend_right=*/false,
                    kPredWidth, kPredHeight);
}

std::string ContextPredictor::GetPredStr(const CodedBlock& cb,
                                         Channel channel) const {
  const int16_t* const context =
      cb.GetContext(channel, /*fill_in=*/true, /*extend_right=*/false);
  // Always display the top right context.
  const int16_t* const context_tr =
      cb.GetContext(channel, /*fill_in=*/true, /*extend_right=*/true);

  return GetPredStr(context, context_tr, extend_context_, cb.w_pix(),
                    cb.h_pix());
}

std::string ContextPredictor::GetPredStr(const int16_t context[],
                                         const int16_t context_tr[],
                                         bool extend_right, uint32_t width,
                                         uint32_t height) const {
  int16_t output[kMaxBlockSizePix2];
  Predict(extend_right ? context_tr : context, width, height,
          output, /*step=*/width);
  return GetContextAndBlockPixelsStr(context, context_tr, width, height,
                                     output, /*step=*/width);
}

//------------------------------------------------------------------------------

// The DC predictor sets all pixels to the average of top and/or left context.
// TODO(vrabaud) remove Large in the name once the other DCPredictor class from
//               context is removed.
class LargeDCPredictor : public ContextPredictor {
 public:
  enum Type { kAll, kLeft, kTop };
  explicit LargeDCPredictor(Type type) : type_(type) {}
  ModeContext mode_context() const override { return ModeContext::kDC; }

  std::string GetName() const override {
    const char* const kTypeStr[]{"all", "left", "top"};
    std::string name = "DC predictor (avg of ";
    name += kTypeStr[type_];
    name += " context)";
    return name;
  }
  void Draw(const WP2::Rectangle& rect,
            ArgbBuffer* const debug_output) const override {}

 protected:
  void Predict(const int16_t context[], uint32_t w, uint32_t h,
               int16_t output[], uint32_t step) const override {
    // DC predictors will never predict values outside of the input range so we
    // don't need to pass in accurate min/max bounds.
    const int16_t min = std::numeric_limits<int16_t>::min();
    const int16_t max = std::numeric_limits<int16_t>::max();
    BasePredictors[BPRED_DC + type_](context, w, h, min, max, output, step);
  }

  const Type type_;
};

//------------------------------------------------------------------------------

// TrueMotion predictor.
class TMPredictor : public ContextPredictor {
 public:
  TMPredictor(int16_t min_value, int16_t max_value)
      : min_value_(min_value), max_value_(max_value) {}
  ModeContext mode_context() const override {
    // Just like 135 degrees.
    return ModeContext::kLeft;
  }

  std::string GetName() const override { return "TM predictor"; }
  void Draw(const WP2::Rectangle& rect,
            ArgbBuffer* const debug_output) const override {}

 protected:
  void Predict(const int16_t context[], uint32_t w, uint32_t h,
               int16_t output[], uint32_t step) const override {
    BasePredictors[BPRED_TM](context, w, h, min_value_, max_value_, output,
                             step);
  }

 private:
  const int16_t min_value_;
  const int16_t max_value_;
};

//------------------------------------------------------------------------------

// Smooth predictor that interpolates pixels horizontally, vertically or a mix
// of both.
class SmoothPredictor : public ContextPredictor {
 public:
  enum class SmoothType {
    k2DSmooth,
    kVerticalSmooth,
    kHorizontalSmooth,
    kNumSmoothType
  };
  explicit SmoothPredictor(SmoothType type) : type_(type) {}
  ModeContext mode_context() const override {
    if (type_ == SmoothType::k2DSmooth) {
      return ModeContext::kDC;
    } else if (type_ == SmoothType::kVerticalSmooth) {
      return ModeContext::kVertical;
    } else {
      return ModeContext::kHorizontal;
    }
  }

  std::string GetName() const override {
    const char* const kTypeStr[]{"2D", "vertical", "horizontal"};
    std::string name = "smooth predictor (";
    name += kTypeStr[(int)type_];
    name += ")";
    return name;
  }
  void Draw(const WP2::Rectangle& rect,
            ArgbBuffer* const debug_output) const override {}

 protected:
  void Predict(const int16_t context[], uint32_t w, uint32_t h,
               int16_t output[], uint32_t step) const override {
    // Smooth predictors will never predict values outside of the input range so
    // we don't need to pass in accurate min/max bounds.
    const int16_t min = std::numeric_limits<int16_t>::min();
    const int16_t max = std::numeric_limits<int16_t>::max();
    const BasePredictor p = (type_ == SmoothType::k2DSmooth)
                                ? BPRED_SMOOTH
                                : (type_ == SmoothType::kVerticalSmooth)
                                      ? BPRED_SMOOTH_V
                                      : BPRED_SMOOTH_H;
    BasePredictors[p](context, w, h, min, max, output, step);
  }

  const SmoothType type_;
};

//------------------------------------------------------------------------------

static uint32_t GetMaxAngleDelta(Channel channel) {
  return (channel == kUChannel || channel == kVChannel)
             ? kDirectionalMaxAngleDeltaUV
             : kDirectionalMaxAngleDeltaYA;
}

// Angle predictor.
class AnglePredictor : public ContextPredictor {
 public:
  // Use AV1 convention to name angle predictors.
  enum class Type {
    // Angles must be in the same order as in Predictors::Pred
    D45_PRED,
    D67_PRED,
    V_PRED,  // 90 degrees
    D113_PRED,
    D135_PRED,
    D157_PRED,
    H_PRED,  // 180 degrees
    D203_PRED,
    Num
  };
  AnglePredictor(Type type, Channel channel, int16_t angle_delta,
                 int16_t min_value, int16_t max_value)
      : min_value_(min_value),
        max_value_(max_value),
        type_(type),
        angle_delta_(angle_delta),
        angle_idx_(AngleIdx(type, channel, angle_delta)) {
    extend_context_ = (angle_idx_ < kAngle_90);
  }

  ModeContext mode_context() const override {
    if (type_ == Type::D45_PRED || type_ == Type::D67_PRED) {
      return ModeContext::kTopRight;
    } else if (type_ == Type::V_PRED) {
      return ModeContext::kVertical;
    } else if (type_ == Type::D113_PRED || type_ == Type::D135_PRED ||
               type_ == Type::D157_PRED || type_ == Type::D203_PRED) {
      return ModeContext::kLeft;
    } else if (type_ == Type::H_PRED) {
      return ModeContext::kHorizontal;
    } else {
      assert(false);
    }
    return ModeContext::kDC;
  }

  // Returns the angle index for precalculated angles.
  // Index goes from 0 (35.35 degrees) to 55 (212.14 degrees).
  static uint8_t AngleIdx(Type type, Channel channel, int16_t angle_step) {
    const int16_t max_delta = kDirectionalMaxAngleDeltaYA;
    if (channel == kUChannel || channel == kVChannel) {
      // Note that the U/V prediction angles are not exactly even spaced.
      // But they map to already existing angles for Y/A.
      // TODO(skal): should the YA / UV angles be made independent?
      angle_step *= 2;
    }
    const int idx = (int)type * (2 * max_delta + 1) + max_delta + angle_step;
    assert(idx >= 0 && idx < (int)kNumDirectionalAngles);
    return (uint8_t)idx;
  }

  float GetAngleDegrees() const {
    return 45.f + 22.5f * ((int)angle_idx_ - 3) / 7.;
  }

  void WriteParams(const CodedBlock& cb, Channel channel,
                   SymbolManager* const sm,
                   ANSEncBase* const enc) const override {
    if (channel == kVChannel) return;
    (void)cb;
    (void)sm;
    // By storing the angle step, we are then able, at decoding time, to figure
    // out the sub-mode of an angle predictor.
    const uint32_t max_delta = GetMaxAngleDelta(channel);
    enc->PutRValue(max_delta + angle_delta_, 2 * max_delta + 1, "angle_step");
  }

  void ReadParams(CodedBlock* const cb, Channel channel, SymbolReader* const sm,
                  ANSDec* const dec) const override {
    if (channel == kVChannel) return;
    (void)sm;
    CodedBlock::CodingParams* const params = cb->GetCodingParams(channel);
    // TODO(vrabaud) explore using a symbol for the step, like in AV1.
    const uint32_t max_delta = GetMaxAngleDelta(channel);
    const uint32_t angle_step =
        dec->ReadRValue(2 * max_delta + 1, "angle_step");
    params->pred_sub_mode = angle_step;
  }

  std::string GetName() const override {
    const int angle_deg = std::round(GetAngleDegrees());
    return "angle predictor (" + std::to_string(angle_deg) +
           " degrees, idx=" + std::to_string((int)angle_idx_) + ")";
  }

  void Draw(const WP2::Rectangle& rect,
            ArgbBuffer* const debug_output) const override {
    const float angle = (M_PI / 180.f) * GetAngleDegrees();
    const float cos = std::cos(angle);
    const float sin = std::sin(angle);
    const uint32_t pixel_depth = WP2FormatBpp(debug_output->format);
    const uint32_t l =
        std::sqrt(rect.width * rect.width + rect.height * rect.height) / 2;
    for (float i = 0; i < l; i += 0.5) {
      const int y = std::round(rect.y + (rect.height - 1) / 2.f - i * sin);
      const int x = std::round(rect.x + (rect.width - 1) / 2.f + i * cos);
      if (!rect.Contains(x, y)) break;
      const uint8_t color[4] = {255, 255, 0, 0};
      uint8_t* const dst_pixel =
          (uint8_t*)debug_output->GetRow(y) + x * pixel_depth;
      std::memcpy(dst_pixel, color, pixel_depth);
    }
  }

 protected:
  void Predict(const int16_t context[], uint32_t w, uint32_t h,
               int16_t output[], uint32_t step) const override {
    const uint32_t context_size =
      (angle_idx_ < kAngle_90) ? ContextWithTRSize(w, h) : ContextSize(w, h);
    for (uint32_t i = 0; i < context_size; ++i) {
      assert(context[i] >= min_value_ && context[i] <= max_value_);
    }
    // SimpleAnglePredictor is faster and currently gives better results than
    // BaseAnglePredictor.
    const uint32_t w_idx = TrfLog2[w];
    const uint32_t h_idx = TrfLog2[h];
    SimpleAnglePredictor(angle_idx_, context, w_idx, h_idx, output, step);
  }
  const int16_t min_value_;
  const int16_t max_value_;

 private:
  Type type_;
  const int16_t angle_delta_;
  uint8_t angle_idx_;
};

//------------------------------------------------------------------------------

// Fuse predictor.
class FusePredictor : public ContextPredictor {
 public:
  FusePredictor(float strength, int16_t min_value, int16_t max_value)
      : strength_(strength), min_value_(min_value), max_value_(max_value) {
    PrecomputeLargeWeightTable(strength_, table_);
  }
  ModeContext mode_context() const override {
    // TODO(vrabaud) find the right one.
    return ModeContext::kDC;
  }

  std::string GetName() const override {
    std::string name = "fuse predictor (strength: ";
    name += std::to_string(strength_);
    name += ")";
    return name;
  }

  void Draw(const WP2::Rectangle& rect,
            ArgbBuffer* const debug_output) const override {}

 protected:
  void Predict(const int16_t context[], uint32_t w, uint32_t h,
               int16_t output[], uint32_t step) const override {
    BaseFusePredictor(table_, context, w, h, output, step, min_value_,
                      max_value_);
  }

  const float strength_;
  LargeWeightTable table_;
  const int16_t min_value_;
  const int16_t max_value_;
};

//------------------------------------------------------------------------------

// Extrapolates a gradient from the bottom-left, top-left and top-right corners.
// This is sometimes more powerful than other predictors because it can guess
// values outside the ones in the context.
class GradientPredictor : public ContextPredictor {
 public:
  GradientPredictor(int16_t min_value, int16_t max_value)
      : min_value_(min_value), max_value_(max_value) {}
  ModeContext mode_context() const override {
    // TODO(vrabaud) find the right one.
    // Just like 135 degrees.
    return ModeContext::kLeft;
  }

  std::string GetName() const override { return "gradient predictor"; }
  void Draw(const WP2::Rectangle& rect,
            ArgbBuffer* const debug_output) const override {}
  bool OverridesTransform(Channel channel,
                          TransformPair* const tf) const override {
    if (tf != nullptr) *tf = kAdstAdst;
    return true;
  }

 protected:
  void Predict(const int16_t c[], uint32_t w, uint32_t h, int16_t output[],
               uint32_t step) const override {
    // Average the three known (or extrapolated) corners.
    const int32_t bottom_left = DivRound(c[0] + c[1] + c[2], 3);
    const int32_t top_left =
        DivRound(c[h - 2] + c[h - 1] + c[h] + c[h + 1] + c[h + 2], 5);
    const int32_t top_right =
        DivRound(c[h + w - 1] + c[h + w] + c[h + w + 1], 3);
    // Extrapolate the middle value and from it, the bottom-right corner.
    const int32_t middle = DivRound(bottom_left + top_right, 2);
    const int32_t bottom_right =
        Clamp(middle - top_left + middle, min_value_, max_value_);

    // Create a gradient by bidimensional interpolation.
    const int32_t max_x = (int32_t)w - 1, max_y = (int32_t)h - 1;
    for (int32_t y = 0; y <= max_y; ++y) {
      const int32_t left =
          DivRound(top_left * (max_y - y) + bottom_left * y, max_y);
      const int32_t right =
          DivRound(top_right * (max_y - y) + bottom_right * y, max_y);
      for (int32_t x = 0; x <= max_x; ++x) {
        output[x] = DivRound(left * (max_x - x) + right * x, max_x);
        assert(output[x] >= min_value_ && output[x] <= max_value_);
      }
      output += step;
    }
  }

  const int32_t min_value_;
  const int32_t max_value_;
};

//------------------------------------------------------------------------------
// Alpha predictor.

void AdjustAlphaPrediction(const CodedBlock& cb, CSPTransform transform,
                           int16_t output[], uint32_t step);

// Predictor for the alpha plane, that wraps a normal predictor, and adjusts
// its prediction to take into account the fact RGB is premultiplied, therefore
// alpha >= max(r, g, b)
// This predictor assumes that YUV is encoded/reconstructed before alpha.
// The supplied Predictor is not owned, so we don't delete it.
class AlphaPredictor : public Predictor {
 public:
  // Takes ownership of 'pred'.
  AlphaPredictor(Predictor* const pred, const CSPTransform& transform)
      : pred_(pred), transform_(transform) {}

  ~AlphaPredictor() override { delete pred_; }

  void Predict(const CodedBlock& cb, Channel channel, int16_t output[],
               uint32_t step) const override {
    assert(channel == kAChannel);  // Should only be used for alpha.
    pred_->Predict(cb, channel, output, step);
    AdjustAlphaPrediction(cb, transform_, output, step);
  }

  std::string GetName() const override {
    return "alpha predictor based on " + pred_->GetName();
  }
  std::string GetFakePredStr() const override {
    printf("alpha predictor based on:\n");
    return pred_->GetFakePredStr();
  }
  std::string GetPredStr(const CodedBlock& cb, Channel channel) const override {
    std::string str = pred_->GetPredStr(cb, channel);
    str += "==> adjusted by alpha predictor:\n";
    int16_t output[kMaxBlockSizePix2];
    Predict(cb, channel, output, cb.w_pix());
    const int16_t* const context = cb.GetContext(
      channel, /*fill_in=*/true, /*extend_right=*/false);
    str +=
        GetContextAndBlockPixelsStr(context, /*context_tr=*/nullptr, cb.w_pix(),
                                    cb.h_pix(), output, cb.w_pix());
    return str;
  }
  void Draw(const WP2::Rectangle& rect,
            ArgbBuffer* const debug_output) const override {
    pred_->Draw(rect, debug_output);
  }

  uint32_t mode() const override { return pred_->mode(); }
  void ComputeParams(CodedBlock* const cb, Channel channel) const override {
    pred_->ComputeParams(cb, channel);
  }
  void WriteParams(const CodedBlock& cb, Channel channel,
                   SymbolManager* const sm,
                   ANSEncBase* const enc) const override {
    pred_->WriteParams(cb, channel, sm, enc);
  }
  void ReadParams(CodedBlock* const cb, Channel channel, SymbolReader* const sm,
                  ANSDec* const dec) const override {
    pred_->ReadParams(cb, channel, sm, dec);
  }
  bool DependsOnLuma() const override { return true; }

 protected:
  Predictor* const pred_;
  const CSPTransform& transform_;
};

//------------------------------------------------------------------------------

const Predictor* Predictors::GetFirstWithMode(uint32_t mode) const {
  return first_with_mode_[mode];
}

const Predictor* Predictors::GetWithMode(uint32_t mode,
                                         uint32_t sub_mode) const {
  uint32_t ind = 0;
  for (const auto& p : preds_) {
    assert(p->mode() != Predictor::kInvalidMode);
    if (p->mode() == mode) break;
    ++ind;
  }
  if (ind < preds_.size()) {
    return preds_[ind + sub_mode];
  } else {
    assert(false);
    return nullptr;
  }
}

// Get the maximum mode contained in a predictor vector.
uint32_t Predictors::GetMaxMode() const {
  return preds_no_angle_.size() + kAnglePredNum - 1;
}

WP2Status Predictors::FillImpl(const Pred* const mapping, uint32_t num_modes,
                               int16_t min_value, int16_t max_value,
                               const CSPTransform* const transform) {
  const int max_delta = GetMaxAngleDelta(channel_);
  WP2_CHECK_ALLOC_OK(preds_.reserve(
      num_modes + 2 * max_delta * (uint32_t)AnglePredictor::Type::Num));
  assert(num_modes >= (uint32_t)AnglePredictor::Type::Num);
  WP2_CHECK_ALLOC_OK(
      preds_no_angle_.reserve(num_modes - (uint32_t)AnglePredictor::Type::Num));

  for (uint8_t i = 0; i < num_modes; ++i) {
    Pred p = mapping[i];
    WP2::Predictor* pred = nullptr;

    bool is_angle = false;
    switch (p) {
      case Pred::kDcAll:
        pred = new (WP2Allocable::nothrow)
            LargeDCPredictor(LargeDCPredictor::kAll);
        break;
      case Pred::kDcLeft:
        pred = new (WP2Allocable::nothrow)
            LargeDCPredictor(LargeDCPredictor::kLeft);
        break;
      case Pred::kDcTop:
        pred = new (WP2Allocable::nothrow)
            LargeDCPredictor(LargeDCPredictor::kTop);
        break;
      case Pred::kSmooth2D:
        pred = new (WP2Allocable::nothrow)
            SmoothPredictor(SmoothPredictor::SmoothType::k2DSmooth);
        break;
      case Pred::kSmoothVertical:
        pred = new (WP2Allocable::nothrow)
            SmoothPredictor(SmoothPredictor::SmoothType::kVerticalSmooth);
        break;
      case Pred::kSmoothHorizontal:
        pred = new (WP2Allocable::nothrow)
            SmoothPredictor(SmoothPredictor::SmoothType::kHorizontalSmooth);
        break;
      case Pred::kTrueMotion:
        pred = new (WP2Allocable::nothrow) TMPredictor(min_value, max_value);
        break;
      case Pred::kGradient:
        pred =
            new (WP2Allocable::nothrow) GradientPredictor(min_value, max_value);
        break;
      case Pred::kAngle45:
      case Pred::kAngle67:
      case Pred::kAngle90:
      case Pred::kAngle113:
      case Pred::kAngle135:
      case Pred::kAngle157:
      case Pred::kAngle180:
      case Pred::kAngle203:
        is_angle = true;
        for (int step = -max_delta; step <= max_delta; ++step) {
          is_angle = true;
          // We assume they're in the same order in ContextPredictor and in
          // AnglePredictor::Type
          const auto type =
              (AnglePredictor::Type)((int)p - (int)Pred::kAngle45);
          pred = new (WP2Allocable::nothrow)
              AnglePredictor(type, channel_, step, min_value, max_value);
          WP2_CHECK_ALLOC_OK(pred != nullptr);
          pred->SetMode(i);
          if (channel_ == kAChannel) {
            AlphaPredictor* const alpha_pred =
                new (WP2Allocable::nothrow) AlphaPredictor(pred, *transform);
            if (alpha_pred == nullptr) {
              delete pred;
              return WP2_STATUS_OUT_OF_MEMORY;
            }
            pred = alpha_pred;
          }
          if (!preds_.push_back(pred)) {
            delete pred;
            return WP2_STATUS_OUT_OF_MEMORY;
          }
          // Store the main angle.
          if (step == 0) preds_main_angle_[(uint32_t)type] = pred;
        }
        break;
      case Pred::kCfl:
        pred = new (WP2Allocable::nothrow) CflPredictor(min_value, max_value);
        break;
      case Pred::kSignalingCfl:
        pred = new (WP2Allocable::nothrow)
            SignalingCflPredictor(min_value, max_value);
        break;
      case Pred::kZero:
        pred = new (WP2Allocable::nothrow) ZeroPredictor();
        break;
      default:
        assert(false);
        break;
    }

    if (is_angle) continue;

    WP2_CHECK_ALLOC_OK(pred != nullptr);
    pred->SetMode(i);

    if (channel_ == kAChannel) {
      WP2::Predictor* const alpha_pred =
          new (WP2Allocable::nothrow) AlphaPredictor(pred, *transform);
      if (alpha_pred == nullptr) {
        delete pred;
        return WP2_STATUS_OUT_OF_MEMORY;
      }
      pred = alpha_pred;
    }

    if (!preds_.push_back(pred)) {
      delete pred;
      return WP2_STATUS_OUT_OF_MEMORY;
    }
    WP2_CHECK_ALLOC_OK(preds_no_angle_.push_back(pred));
  }
  // Pre-compute the first predictors with a given mode.
  assert(GetMaxMode() < first_with_mode_.size());
  for (uint32_t mode = 0; mode <= GetMaxMode(); ++mode) {
    first_with_mode_[mode] = GetWithMode(mode, 0);
  }
  return WP2_STATUS_OK;
}

// Struct computing some predictor scores.
class PredictorScorer {
 public:
  PredictorScorer(const EncoderConfig& config, const Rectangle& tile_rect,
                  const BlockContext& context, const Segment& segment,
                  Channel channel, uint32_t num_channels, bool reduced,
                  Counters* const counters, float best_score)
      : config_(config),
        tile_rect_(tile_rect),
        context_(context),
        segment_(segment),
        channel_(channel),
        num_channels_(num_channels),
        reduced_(reduced),
        counters_(counters),
        score_(0.f),
        dist_(0),
        best_score_(best_score),  // Needed for VDebug 'is_best'.
        best_dist_(std::numeric_limits<uint32_t>::max()),
        best_predictor_(nullptr),
        best_tf_(kUnknownTf) {}

  // Computes the score and dist of the given predictor 'p' applied to 'cb'.
  WP2Status UpdateUVPredScore(const Predictor& p, CodedBlock* const cb) {
    uint32_t dist[4];
    float scores[4], res_rate[4], pred_rate[4];
    for (Channel channel : {kUChannel, kVChannel}) {
      const QuantMtx& quant = segment_.GetQuant(channel);
      cb->Quantize(quant, channel, reduced_);
      WP2_CHECK_STATUS(cb->GetRates(
          tile_rect_, context_, quant, channel, num_channels_,
          /*tf_rate=*/0.f, counters_, &dist[channel], &res_rate[channel],
          &pred_rate[channel], &scores[channel]));
    }

    score_ = (scores[kUChannel] + scores[kVChannel]) / 2.f;
    dist_ = (dist[kUChannel] + dist[kVChannel]) / 2;

    // Store scores for debug.
    for (Channel channel : {kUChannel, kVChannel}) {
      const QuantMtx& quant = segment_.GetQuant(channel);
      cb->StorePredictionScore(
          config_, tile_rect_, channel, p, cb->GetCodingParams(channel)->tf,
          dist[channel], cb->lambda_mult_ * quant.lambda2, res_rate[channel],
          pred_rate[channel], /*tf_rate=*/0.f, scores[channel],
          /*is_best=*/(score_ < best_score_));
    }
    return WP2_STATUS_OK;
  }

  // Computes the score and dist of the given predictor 'p' and transform 'tf'
  // applied to 'cb'.
  WP2Status UpdateYAPredScore(const Predictor& p, const Plane16& prediction,
                              TransformPair tf, uint32_t num_transforms,
                              BlockCoeffs32* const res, CodedBlock* const cb) {
    assert(channel_ == kYChannel || channel_ == kAChannel);
    CodedBlock::CodingParams* const params = cb->GetCodingParams(channel_);
    params->tf = tf;
    assert(cb->GetImplicitTf(channel_) == kUnknownTf ||
           cb->GetImplicitTf(channel_) == tf);

    const QuantMtx& quant = segment_.GetQuant(channel_);
    float res_rate, pred_rate, tf_rate;

    cb->TransformAndReconstruct(quant, channel_, reduced_, prediction, res);

    WP2_CHECK_STATUS(cb->TransformRate(channel_, num_transforms,
                                       counters_->transform(), &tf_rate));
    WP2_CHECK_STATUS(cb->GetRates(tile_rect_, context_, quant, channel_,
                                  num_channels_, tf_rate, counters_, &dist_,
                                  &res_rate, &pred_rate, &score_));

    // Store score for debug.
    cb->StorePredictionScore(config_, tile_rect_, channel_, p, params->tf,
                             dist_, cb->lambda_mult_ * quant.lambda1, res_rate,
                             pred_rate, tf_rate, score_,
                             /*is_best=*/(score_ < best_score_));
    return WP2_STATUS_OK;
  }

  // Computes the score of the predictior 'p' and stores its score, dist etc. if
  // it is the best.
  WP2Status UpdatePredScore(const Predictor* const p, const Plane16& prediction,
                            TransformPair tf, uint32_t num_tf,
                            BlockCoeffs32* const res, CodedBlock* const cb) {
    if (channel_ == kUChannel) {
      WP2_CHECK_STATUS(UpdateUVPredScore(*p, cb));
    } else {
      WP2_CHECK_STATUS(UpdateYAPredScore(*p, prediction, tf, num_tf, res, cb));
    }
    if (score_ < best_score_) {
      best_score_ = score_;
      best_dist_ = dist_;
      best_predictor_ = p;
      best_tf_ = tf;
    }
    return WP2_STATUS_OK;
  }

  float GetLastScore() const { return score_; }
  uint32_t GetLastDist() const { return dist_; }
  float GetBestScore() const { return best_score_; }
  uint32_t GetBestDist() const { return best_dist_; }
  const Predictor* GetBestPredictor() const { return best_predictor_; }
  TransformPair GetBestTf() const { return best_tf_; }

 private:
  const EncoderConfig& config_;
  const Rectangle tile_rect_;
  const BlockContext& context_;
  const Segment& segment_;
  const Channel channel_;
  const uint32_t num_channels_;
  const bool reduced_;
  Counters* counters_;

  float score_;
  uint32_t dist_;
  float best_score_;
  uint32_t best_dist_;
  const Predictor* best_predictor_;
  TransformPair best_tf_;
};

//------------------------------------------------------------------------------

WP2Status Predictors::FindBest(
    const EncoderConfig& config, const Rectangle& tile_rect,
    uint32_t num_channels, const BlockContext& context, const Segment& segment,
    bool reduced, const Vector<TransformPair>& transforms,
    Counters* const counters, CodedBlock* const cb, float* const best_score,
    const Predictor** const best_predictor,
    TransformPair* const best_tf) const {
  *best_predictor = nullptr;
  // This function should only be called once for both U and V.
  assert(channel_ != kVChannel);
  const bool is_uv = (channel_ == kUChannel);
  const uint32_t num_tf = std::max((uint32_t)transforms.size(), 1u);
  PredictorScorer scorer(
      config, tile_rect, context, segment, channel_, num_channels, reduced,
      counters, is_uv ? std::numeric_limits<float>::max() : *best_score);

  int16_t prediction[kMaxBlockSizePix2];
  Plane16 pred_view;
  WP2_CHECK_STATUS(pred_view.SetView(prediction, kMaxBlockSizePix,
                                     kMaxBlockSizePix, kMaxBlockSizePix));
  BlockCoeffs32 res;
  TransformPair overridden_tf;

  // Test the non-angle predictors.
  for (const Predictor* const p : preds_no_angle_) {
    cb->SetPredictor(channel_, p);
    if (!is_uv) {
      cb->PredictBlock(channel_, pred_view.Row(0), pred_view.Step());
      cb->GetResiduals(channel_, pred_view, &res);
    }
    if (p->OverridesTransform(channel_, &overridden_tf)) {
      WP2_CHECK_STATUS(scorer.UpdatePredScore(p, pred_view, overridden_tf,
                                              /*num_tf=*/1, &res, cb));
    } else {
      for (TransformPair tf : transforms) {
        WP2_CHECK_STATUS(
            scorer.UpdatePredScore(p, pred_view, tf, num_tf, &res, cb));
      }
    }
  }
  const float best_score_no_angle = scorer.GetBestScore();
  const uint32_t best_dist_no_angle = scorer.GetBestDist();

  // Test the main angle predictors.
  struct AngleScore {
    float score;
    uint32_t dist;
    const Predictor* p;
    TransformPair tf;
  };
  VectorNoCtor<AngleScore> angle_scores;
  WP2_CHECK_ALLOC_OK(
      angle_scores.reserve((uint32_t)AnglePredictor::Type::Num * num_tf));

  for (const Predictor* const p : preds_main_angle_) {
    cb->SetPredictor(channel_, p);
    if (!is_uv) {
      cb->PredictBlock(channel_, pred_view.Row(0), pred_view.Step());
      cb->GetResiduals(channel_, pred_view, &res);
    }
    if (p->OverridesTransform(channel_, &overridden_tf)) {
      WP2_CHECK_STATUS(scorer.UpdatePredScore(p, pred_view, overridden_tf,
                                              /*num_tf=*/1, &res, cb));
      angle_scores.push_back_no_resize(
          {scorer.GetLastScore(), scorer.GetLastDist(), p, overridden_tf});
    } else {
      for (TransformPair tf : transforms) {
        WP2_CHECK_STATUS(
            scorer.UpdatePredScore(p, pred_view, tf, num_tf, &res, cb));
        angle_scores.push_back_no_resize(
            {scorer.GetLastScore(), scorer.GetLastDist(), p, tf});
      }
    }
  }

  std::sort(angle_scores.begin(), angle_scores.end(),
            [](const AngleScore& a, const AngleScore& b) {
              return (a.score < b.score);
            });

  // Choose to refine all the ones within a percentage of the lowest score. This
  // helps refine predictors that give similar scores.
  const float min_allowed_score =
      angle_scores[0].score +
      (angle_scores.back().score - angle_scores[0].score) * config.speed / 9.f +
      1;
  assert(config.speed != 9.f || angle_scores.back().score <= min_allowed_score);
  for (const AngleScore& s : angle_scores) {
    // We refine every predictor with a better score, or distortion and make
    // sure we refine all the ones within a percentage of the lowest score.
    if (s.score <= best_score_no_angle || s.dist <= best_dist_no_angle ||
        s.score <= min_allowed_score) {
      // Refine the angle predictor.
      for (const Predictor* const p : preds_) {
        if (p->mode() == s.p->mode() && p != s.p) {
          cb->SetPredictor(channel_, p);
          if (!is_uv) {
            cb->PredictBlock(channel_, pred_view.Row(0), pred_view.Step());
            cb->GetResiduals(channel_, pred_view, &res);
          }
          const bool ovrd_tf = p->OverridesTransform(channel_, &overridden_tf);
          WP2_CHECK_STATUS(scorer.UpdatePredScore(
              p, pred_view, s.tf, ovrd_tf ? 1 : num_tf, &res, cb));
        }
      }
    }
  }

  if (scorer.GetBestScore() < *best_score) {
    assert(scorer.GetBestPredictor() != nullptr);
    *best_score = scorer.GetBestScore();
    *best_predictor = scorer.GetBestPredictor();
    *best_tf = scorer.GetBestTf();
  }
  return WP2_STATUS_OK;
}

//------------------------------------------------------------------------------

YPredictors::YPredictors() { channel_ = kYChannel; }

WP2Status YPredictors::Fill(int16_t min_value, int16_t max_value) {
  constexpr Pred kYPreds[] = {
      Pred::kDcAll,      Pred::kDcLeft,         Pred::kDcTop,
      Pred::kSmooth2D,   Pred::kSmoothVertical, Pred::kSmoothHorizontal,
      Pred::kTrueMotion, Pred::kGradient,       Pred::kAngle45,
      Pred::kAngle67,    Pred::kAngle90,        Pred::kAngle113,
      Pred::kAngle135,   Pred::kAngle157,       Pred::kAngle180,
      Pred::kAngle203};

  STATIC_ASSERT_ARRAY_SIZE(kYPreds, kYPredModeNum);
  WP2_CHECK_STATUS(FillImpl(kYPreds, kYPredModeNum, min_value, max_value,
                            /*transform=*/nullptr));
  return WP2_STATUS_OK;
}

//------------------------------------------------------------------------------

APredictors::APredictors() { channel_ = kAChannel; }

WP2Status APredictors::Fill(int16_t min_value, int16_t max_value,
                            const CSPTransform* const transform) {
  constexpr Pred kAPreds[] = {Pred::kDcAll,          Pred::kSmooth2D,
                              Pred::kSmoothVertical, Pred::kSmoothHorizontal,
                              Pred::kTrueMotion,     Pred::kAngle45,
                              Pred::kAngle67,        Pred::kAngle90,
                              Pred::kAngle113,       Pred::kAngle135,
                              Pred::kAngle157,       Pred::kAngle180,
                              Pred::kAngle203,       Pred::kCfl,
                              Pred::kSignalingCfl,   Pred::kZero};

  assert(transform != nullptr);
  STATIC_ASSERT_ARRAY_SIZE(kAPreds, kAPredModeNum);
  WP2_CHECK_STATUS(
      FillImpl(kAPreds, kAPredModeNum, min_value, max_value, transform));
  return WP2_STATUS_OK;
}

//------------------------------------------------------------------------------

UVPredictors::UVPredictors() { channel_ = kUChannel; }

WP2Status UVPredictors::Fill(int16_t min_value, int16_t max_value) {
  constexpr Pred kUVPreds[] = {
      Pred::kCfl,      Pred::kDcAll,    Pred::kSignalingCfl, Pred::kSmooth2D,
      Pred::kAngle45,  Pred::kAngle67,  Pred::kAngle90,      Pred::kAngle113,
      Pred::kAngle135, Pred::kAngle157, Pred::kAngle180,     Pred::kAngle203};

  STATIC_ASSERT_ARRAY_SIZE(kUVPreds, kUVPredModeNum);
  WP2_CHECK_STATUS(FillImpl(kUVPreds, kUVPredModeNum, min_value, max_value,
                            /*transform=*/nullptr));
  return WP2_STATUS_OK;
}

//------------------------------------------------------------------------------

CflPredictor::CflPredictor(int16_t min_value, int16_t max_value)
    : min_value_(min_value), max_value_(max_value) {}

static void LinearRegression(const int16_t* const luma,
                             const int16_t* const chroma, uint32_t size,
                             int16_t* const a, int16_t* const b) {
  int32_t num_values = 0;
  int32_t l_sum = 0;
  int32_t uv_sum = 0;
  int32_t l_uv_sum = 0;
  int32_t l_l_sum = 0;
  for (uint32_t i = 0; i < size; ++i) {
    if (luma[i] == CodedBlock::kMissing) continue;
    ++num_values;
    const int32_t l = luma[i], uv = chroma[i];
    l_sum += l;
    uv_sum += uv;
    l_uv_sum += l * uv;
    l_l_sum += l * l;
  }
  if (num_values == 0) {
    *a = *b = 0;
    return;
  }
  const int64_t num = (int64_t)l_uv_sum * num_values - (int64_t)l_sum * uv_sum;
  const int64_t den = (int64_t)l_l_sum * num_values - (int64_t)l_sum * l_sum;
  if (num == 0 || den == 0) {
    *a = 0;
    *b = (int16_t)DivRound(LeftShift(uv_sum, CflPredictor::kBPrecShift),
                           num_values);    // no need to s-clamp to 16bit
  } else {
     const int32_t tmp_a = DivRound(LeftShift(num, CflPredictor::kAPrecShift),
                                    den);
     *a = (int16_t)ClampToSigned(tmp_a, 16);
     int64_t tmp_b = den * uv_sum - num * l_sum;
     tmp_b = DivRound(LeftShift(tmp_b, CflPredictor::kBPrecShift),
                                num_values * den);
     *b = (int16_t)ClampToSigned(tmp_b, 16);
  }
}

void CflPredictor::ContextLinearRegression(Channel channel,
                                           const CodedBlock& cb,
                                           int16_t* const a,
                                           int16_t* const b) const {
  const int16_t* const luma_context =
      cb.GetContext(kYChannel, /*fill_in=*/false, /*extend_right=*/false);
  const int16_t* const chroma_context =
      cb.GetContext(channel, /*fill_in=*/false, /*extend_right=*/false);
  const uint32_t context_size = ContextSize(cb.w_pix(), cb.h_pix());
  LinearRegression(luma_context, chroma_context, context_size, a, b);
}

void CflPredictor::BlockLinearRegression(Channel channel, const CodedBlock& cb,
                                         int16_t* const a,
                                         int16_t* const b) const {
  const Plane16& luma_plane = cb.out_.GetChannel(kYChannel);
  const Plane16& chroma_plane = cb.in_.GetChannel(channel);
  int16_t luma[kMaxBlockSizePix2], chroma[kMaxBlockSizePix2];

  for (uint32_t k = 0, y = 0; y < cb.h_pix(); ++y) {
    const int16_t* const y_src = luma_plane.Row(y);
    const int16_t* const c_src = chroma_plane.Row(y);
    for (uint32_t x = 0; x < cb.w_pix(); ++x, ++k) {
      luma[k] = y_src[x];
      chroma[k] = c_src[x];
    }
  }
  LinearRegression(luma, chroma, cb.h_pix() * cb.w_pix(), a, b);
}

void CflPredictor::DoPredict(const CodedBlock& cb, Channel channel,
                             uint32_t x_start, uint32_t y_start, uint32_t w,
                             uint32_t h, int16_t* output, uint32_t step) const {
  const Plane16& luma = cb.GetContextPlane(kYChannel);
  int16_t a, b;
  ContextLinearRegression(channel, cb, &a, &b);
  const int32_t shifted_b =
      ChangePrecision<int32_t>(b, kBPrecShift, kAPrecShift) +
      (1 << kAPrecShift >> 1);   // we pre-incorporate the rounding constant

  for (uint32_t y = 0; y < h; ++y) {
    const int16_t* const src = &luma.At(x_start, y_start + y);
    for (uint32_t x = 0; x < w; ++x) {
      const int32_t v = ((int32_t)a * src[x] + shifted_b) >> kAPrecShift;
      output[x] = (int16_t)Clamp<int32_t>(v, min_value_, max_value_);
    }
    output += step;
  }
}

void CflPredictor::Predict(const CodedBlock& cb, Channel channel,
                           int16_t output[], uint32_t step) const {
  assert(channel != kYChannel);
  DoPredict(cb, channel, 0, 0, cb.w_pix(), cb.h_pix(), output, step);
}

std::string CflPredictor::GetName() const {
  return "chroma-from-luma predictor";
}

std::string CflPredictor::GetFakePredStr() const {
  return "(fake prediction not supported)\n";
}

std::string CflPredictor::GetPredStr(const CodedBlock& cb,
                                     Channel channel) const {
  std::string str;
  int16_t a, b;
  ContextLinearRegression(channel, cb, &a, &b);
  const float fa = (float)a / (1 << kAPrecShift);
  const float fb = (float)b / (1 << kBPrecShift);
  const Plane16& luma = cb.GetContextPlane(kYChannel);
  WP2SAppend(&str, "Chroma = %.2f * Luma + %.2ff\n", fa, fb);
  int16_t output[kMaxBlockSizePix2];
  Predict(cb, channel, output, kMaxBlockSizePix);
  // Display subset of values to fit screen.
  WP2SAppend(&str, "       LUMA          |            chroma\n");
  for (uint32_t j = 0; j < kPredHeight; ++j) {
    for (uint32_t i = 0; i < kPredWidth; ++i) {
      WP2SAppend(&str, "%4d ", luma.At(i, j));
    }
    WP2SAppend(&str, " | ");
    for (uint32_t i = 0; i < kPredWidth; ++i) {
      WP2SAppend(&str, "% 6.1f ", fa * luma.At(i, j) + fb);
    }
    WP2SAppend(&str, "\n");
  }
  if (cb.w() > 1 || cb.h() > 1) WP2SAppend(&str, "(cropped)\n");
  return str;
}

SignalingCflPredictor::SignalingCflPredictor(int16_t min_value,
                                             int16_t max_value)
    : CflPredictor(min_value, max_value) {}

void SignalingCflPredictor::GetParams(const CodedBlock& cb, Channel channel,
                                      uint32_t w, uint32_t h, int32_t* a,
                                      int32_t* b) const {
  const Plane16& luma = cb.GetContextPlane(kYChannel);
  int16_t predicted_a, predicted_b;
  ContextLinearRegression(channel, cb, &predicted_a, &predicted_b);
  *a = predicted_a + GetParam(cb, channel);

  int64_t luma_sum = 0;
  for (uint32_t y = 0; y < h; ++y) {
    const int16_t* const y_src = luma.Row(y);
    for (uint32_t x = 0; x < w; ++x) luma_sum += y_src[x];
  }
  const int32_t average_luma = DivRound<int64_t>(luma_sum, w * h);
  *b = ChangePrecision<int32_t>(predicted_b, kBPrecShift, kAPrecShift) -
       GetParam(cb, channel) * average_luma;
}

void SignalingCflPredictor::DoPredict(const CodedBlock& cb, Channel channel,
                                      uint32_t, uint32_t, uint32_t w,
                                      uint32_t h, int16_t* output,
                                      uint32_t step) const {
  const Plane16& luma = cb.GetContextPlane(kYChannel);
  int32_t a, b;
  GetParams(cb, channel, w, h, &a, &b);
  b += (1 << kAPrecShift >> 1);   // pre-incorporate the rounding constant

  for (uint32_t y = 0; y < h; ++y) {
    const int16_t* const y_src = luma.Row(y);
    for (uint32_t x = 0; x < w; ++x) {
      const int32_t v = (a * y_src[x] + b) >> kAPrecShift;
      output[x] = (int16_t)Clamp<int32_t>(v, min_value_, max_value_);
    }
    output += step;
  }
}

void SignalingCflPredictor::ComputeParams(CodedBlock* const cb,
                                          Channel channel) const {
  int16_t predicted_a, predicted_b;
  ContextLinearRegression(channel, *cb, &predicted_a, &predicted_b);

  int16_t best_a, best_b;
  BlockLinearRegression(channel, *cb, &best_a, &best_b);
  const int16_t a_res =
      ChangePrecision(ClampToSigned(best_a - predicted_a, kAResBits),
                      kAPrecShift, kAResPrecision);
  cb->cfl_[channel - 1] = a_res;
}

void SignalingCflPredictor::WriteParams(const CodedBlock& cb, Channel channel,
                                        SymbolManager* const sm,
                                        ANSEncBase* const enc) const {
  const int16_t a_res =
      ChangePrecision(GetParam(cb, channel), kAPrecShift, kAResPrecision);
  const uint32_t a_unsigned = std::abs(a_res);
  const uint32_t sign = a_res < 0 ? 1 : 0;
  sm->Process(kSymbolCflSlope, a_unsigned, "cfl_slope", enc);
  enc->PutUValue(sign, 1, "cfl_slope_sign");
}

void SignalingCflPredictor::ReadParams(CodedBlock* const cb, Channel channel,
                                       SymbolReader* const sr,
                                       ANSDec* const dec) const {
  const uint32_t a_unsigned =
      sr->Read(/*cluster=*/0, kSymbolCflSlope, "cfl_slope");
  cb->cfl_[channel - 1] = a_unsigned;
  const uint32_t sign = dec->ReadUValue(1, "cfl_slope_sign");
  if (sign) {
    cb->cfl_[channel - 1] *= -1;
  }
}

std::string SignalingCflPredictor::GetName() const {
  return "signaling chroma-from-luma";
}

std::string SignalingCflPredictor::GetPredStr(const CodedBlock& cb,
                                              Channel channel) const {
  const Plane16& luma = cb.GetContextPlane(kYChannel);
  std::string str;
  int16_t predicted_a, predicted_b;
  ContextLinearRegression(channel, cb, &predicted_a, &predicted_b);
  int32_t a, b;
  GetParams(cb, channel, cb.w(), cb.h_pix(), &a, &b);

  const float fa = (float)a / (1 << kAPrecShift);
  const float fb = (float)b / (1 << kAPrecShift);
  WP2SAppend(&str, "Predicted a: %.2f with correction: %.2f\n",
             (float)predicted_a / (1 << kAPrecShift), fa);
  WP2SAppend(&str, "Chroma = %.2f * Luma + %.2ff\n", fa, fb);
  int16_t output[kMaxBlockSizePix2];
  Predict(cb, channel, output, kMaxBlockSizePix);
  // Display subset of values to fit screen.
  WP2SAppend(&str, "       LUMA          |            chroma\n");
  for (uint32_t j = 0; j < kPredHeight; ++j) {
    for (uint32_t i = 0; i < kPredWidth; ++i) {
      WP2SAppend(&str, "%4d ", luma.At(i, j));
    }
    WP2SAppend(&str, " | ");
    for (uint32_t i = 0; i < kPredWidth; ++i) {
      WP2SAppend(&str, "% 6.1f ", fa * luma.At(i, j) + fb);
    }
    WP2SAppend(&str, "\n");
  }
  if (cb.w() > 1 || cb.h() > 1) WP2SAppend(&str, "(cropped)\n");
  return str;
}

int16_t SignalingCflPredictor::GetParam(const CodedBlock& cb,
                                        Channel channel) const {
  return ChangePrecision(cb.cfl_[channel - 1], kAResPrecision, kAPrecShift);
}

//------------------------------------------------------------------------------

void AdjustAlphaPrediction(const CodedBlock& cb, CSPTransform transform,
                           int16_t output[], uint32_t step) {
  for (uint16_t j = 0; j < cb.h_pix(); ++j) {
    const int16_t* const y_src = cb.GetContextPlane(kYChannel).Row(j);
    const int16_t* const u_src = cb.GetContextPlane(kUChannel).Row(j);
    const int16_t* const v_src = cb.GetContextPlane(kVChannel).Row( j);
    for (uint16_t i = 0; i < cb.w_pix(); ++i) {
      int16_t r, g, b;
      transform.ToRGB(y_src[i], u_src[i], v_src[i], &r, &g, &b);

      // Alpha is larger than the largest R, G or B value.
      // TODO(maryla): since this is based on reconstructed RGB, there might be
      // some noise. Add a safety margin at lower quality?
      const int16_t min_alpha_value = std::max({r, g, b, (int16_t)0});
      output[j * step + i] = std::max(output[j * step + i], min_alpha_value);
    }
  }
}

// TODO(vrabaud) Remove
WP2Status InitAlphaPredictors(const CSPTransform& transform,
                              APredictors* const preds) {
  const uint32_t kMinValue = 0;
  const uint32_t kMaxValue = kAlphaMax;

  preds->reset();
  WP2_CHECK_STATUS(preds->Fill(kMinValue, kMaxValue, &transform));
  assert(preds->size() == kAPredNum);

  return WP2_STATUS_OK;
}

//------------------------------------------------------------------------------

}  // namespace WP2
