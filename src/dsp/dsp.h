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
//   Speed-critical functions.
//
// Author: Skal (pascal.massimino@gmail.com)

#ifndef WP2_DSP_DSP_H_
#define WP2_DSP_DSP_H_

#ifdef HAVE_CONFIG_H
#include "src/wp2/config.h"
#endif

#include "src/wp2/base.h"
#include "src/wp2/format_constants.h"

//------------------------------------------------------------------------------
// CPU detection

#if defined(__GNUC__)
# define LOCAL_GCC_VERSION ((__GNUC__ << 8) | __GNUC_MINOR__)
# define LOCAL_GCC_PREREQ(maj, min) \
    (LOCAL_GCC_VERSION >= (((maj) << 8) | (min)))
#else
# define LOCAL_GCC_VERSION 0
# define LOCAL_GCC_PREREQ(maj, min) 0
#endif

// We are inheriting constructors, so we need the following
// minimal versions of the compilers.
#if defined(__clang__)
// clang 3.3, cf https://clang.llvm.org/cxx_status.html
static_assert(LOCAL_GCC_PREREQ(3, 3), "clang needs to be at version >= 3.3");
#elif defined(__GNUC__)
// gcc 4.7, cf https://gcc.gnu.org/projects/cxx-status.html
static_assert(LOCAL_GCC_PREREQ(4, 8), "gcc needs to be at version >= 4.8");
#elif defined(_MSC_VER)
// MSVC 2015, cf https://msdn.microsoft.com/en-us/library/hh567368.aspx
static_assert(_MSC_VER >= 1900, "Visual Studio needs to be at version >= 2015");
#endif

#ifndef __has_builtin
# define __has_builtin(x) 0
#endif

// for now, none of the optimizations below are available in emscripten
#if !defined(EMSCRIPTEN)

#if defined(_MSC_VER) && (defined(_M_X64) || defined(_M_IX86))
#define WP2_MSC_SSE   // Visual C++ SSE4.2 targets
#endif

// WP2_HAVE_* are used to indicate the presence of the instruction set in dsp
// files without intrinsics, allowing the corresponding Init() to be called.
// Files containing intrinsics will need to be built targeting the instruction
// set so should succeed on one of the earlier tests.
#if defined(__SSE4_2__) || defined(WP2_MSC_SSE) || defined(WP2_HAVE_SSE)
#define WP2_USE_SSE    // incorporate all instructions up to SSE4.2 (inclusive)
#include <smmintrin.h>
#endif

#if defined(__AVX2__) || defined(WP2_HAVE_AVX2)
#define WP2_USE_AVX2
#include <immintrin.h>
#endif

#if defined(__ANDROID__) && defined(__ARM_ARCH_7A__)
#define WP2_ANDROID_NEON  // Android targets that might support NEON
#endif

// The intrinsics currently cause compiler errors with arm-nacl-gcc and the
// inline assembly would need to be modified for use with Native Client.
#if (defined(__ARM_NEON__) || defined(WP2_ANDROID_NEON) || \
     defined(__aarch64__) || defined(WP2_HAVE_NEON)) && \
    !defined(__native_client__)
#define WP2_USE_NEON
#endif

#if defined(_MSC_VER) && defined(_M_ARM)
#define WP2_USE_NEON
#define WP2_USE_INTRINSICS
#endif

#if defined(WP2_USE_NEON)
#include <arm_neon.h>
#endif

#if defined(__mips__) && !defined(__mips64) && \
    defined(__mips_isa_rev) && (__mips_isa_rev >= 1) && (__mips_isa_rev < 6)
#define WP2_USE_MIPS32
#if (__mips_isa_rev >= 2)
#define WP2_USE_MIPS32_R2
#if defined(__mips_dspr2) || (__mips_dsp_rev >= 2)
#define WP2_USE_MIPS_DSP_R2
#endif
#endif
#endif

#if defined(__mips_msa) && defined(__mips_isa_rev) && (__mips_isa_rev >= 5)
#define WP2_USE_MSA
#endif

#endif  /* EMSCRIPTEN */

// This macro prevents thread_sanitizer from reporting known concurrent writes.
#define WP2_TSAN_IGNORE_FUNCTION
#if defined(__has_feature)
#if __has_feature(thread_sanitizer)
#undef WP2_TSAN_IGNORE_FUNCTION
#define WP2_TSAN_IGNORE_FUNCTION __attribute__((no_sanitize_thread))
#endif
#endif

#define WP2_UBSAN_IGNORE_UNDEF
#define WP2_UBSAN_IGNORE_UNSIGNED_OVERFLOW
#if defined(__clang__) && defined(__has_attribute)
#if __has_attribute(no_sanitize)
// This macro prevents the undefined behavior sanitizer from reporting
// failures. This is only meant to silence unaligned loads on platforms that
// are known to support them.
#undef WP2_UBSAN_IGNORE_UNDEF
#define WP2_UBSAN_IGNORE_UNDEF \
  __attribute__((no_sanitize("undefined")))

// This macro prevents the undefined behavior sanitizer from reporting
// failures related to unsigned integer overflows. This is only meant to
// silence cases where this well defined behavior is expected.
#undef WP2_UBSAN_IGNORE_UNSIGNED_OVERFLOW
#define WP2_UBSAN_IGNORE_UNSIGNED_OVERFLOW \
  __attribute__((no_sanitize("unsigned-integer-overflow")))
#endif
#endif

typedef enum {
  kSSE,        // everything up to 4.2
  kAVX,
  kAVX2,
  kNEON,
  kMIPS32,
  kMIPSdspR2,
  kMSA,
  // some particular fine-grained CPU types
  kSSE2, kSSE3, kSSE4_1, kSSE4_2,
  kSlowSSSE3,  // special feature for slow SSSE3 architectures
} WP2CPUFeature;

// returns true if the CPU supports the feature.
typedef bool (*WP2CPUInfo)(WP2CPUFeature feature);
WP2_EXTERN WP2CPUInfo WP2GetCPUInfo;
#define WP2_CPUINFO_IS_DEFINED

//------------------------------------------------------------------------------
// Init stub generator

// Defines an init function stub to ensure each module exposes a symbol,
// avoiding a compiler warning.
#define WP2_DSP_INIT_STUB(func) \
  extern void func(); \
  WP2_TSAN_IGNORE_FUNCTION void func() {}

//------------------------------------------------------------------------------
// Encoding

static constexpr uint32_t WP2QBits = 16;  // fixed-point precision for quant
// Quantize: coeffs[] = (res[] * iq[] + bias[]) >> WP2Bits (with sign handling).
// iq[], bias[], res[] and coeffs[] are assumed to be pointing to at
// least 8 elements (uint32_t or int16_t) of memory underneath.
extern void (*WP2Quantize)(const uint32_t iq[], const uint32_t bias[],
                           const int32_t res[], int16_t coeffs[],
                           uint32_t len);
// Performs unclamped out[] = in[] * dequant[] for i in [0..len), and pad the
// rest of out[] with 0s. It's assumed that in[] is padded with 0s for i >= len.
extern void (*WP2Dequantize)(const int16_t in[], const int16_t dequants[],
                             int32_t out[], uint32_t len , uint32_t max_len);

void WP2QuantizeInit();

// must be called before using any of the above
void WP2EncDspInit();

//------------------------------------------------------------------------------
// Decoding

// must be called before anything using the above
void WP2DecDspInit();

//------------------------------------------------------------------------------
// Common utils

typedef void (*WP2ArgbConverterF)(const uint8_t* src,
                                  uint32_t width, uint8_t* dst);
// functions to convert from any format to internal Argb
extern WP2ArgbConverterF WP2ArgbConvertFrom[WP2_FORMAT_NUM];
// functions to convert from internal Argb to any format
extern WP2ArgbConverterF WP2ArgbConvertTo[WP2_FORMAT_NUM];

void WP2ArgbConverterInit();


// Applies the matrix mtx[] to (y,u,v)[] + offset[]. Alias-safe.
extern void (*WP2YUVToCustom)(
  const int16_t* y, const int16_t* u, const int16_t* v,
  const int16_t offset[3], const int16_t mtx[9],
  int16_t* dst0, int16_t* dst1, int16_t* dst2, uint32_t width);

typedef void (*WP2YUVToArgbFunc)(
  const int16_t* y, const int16_t* u, const int16_t* v, const int16_t* a,
  const int16_t avg[3], const int16_t mtx[9], uint8_t* argb, uint32_t width);
// converts to RGB + 0xff-alpha.
// mtx[] is assumed to be 12-bit fixed-point precision matrix.
extern WP2YUVToArgbFunc WP2YUVToXRGB;  // the 'a' pointer is ignored here
extern WP2YUVToArgbFunc WP2YUVToArgb;  // Convert to pre-multiplied Argb
extern WP2YUVToArgbFunc WP2YUVToARGB;  // Convert to un-multiplied ARGB

void WP2CSPConverterInit();

//------------------------------------------------------------------------------
// ANS

namespace WP2 {

// Updates the discrete cumulative distribution function.
// The number of elements in 'cdf_base', 'cdf_var' and 'cumul' must be a
// multiple of 8 at least as big as 'n'.
typedef void (*ANSUpdateCDFFunc)(uint32_t n, const uint16_t cdf_base[],
                                 const uint16_t cdf_var[], uint32_t mult,
                                 uint16_t cumul[]);
extern ANSUpdateCDFFunc ANSUpdateCDF;

void ANSInit();

}  // namespace WP2

//------------------------------------------------------------------------------
// Transforms

enum WP2TransformType { kDct = 0, kAdst, kHadamard, kIdentity, kNumTransforms };
const char* const WP2TransformNames[] = {"DCT", "ADST", "Hadamard", "Identity"};

typedef void (*WP2TransposeF)(const int32_t* in, uint32_t w, uint32_t h,
                              int32_t* out);

typedef void (*WP2TransformF)(const int32_t* input, int32_t* output);

// indexed on size=2,4,8,16,32
extern WP2TransformF WP2FwdDct[5], WP2FwdAdst[5], WP2FwdHadamard[5];
extern WP2TransformF WP2InvDct[5], WP2InvAdst[5], WP2InvHadamard[5];

extern WP2TransposeF WP2Transpose;

void WP2TransformInit();

// These functions work inplace when 'src' and 'dst' point to the same address.
// Entry call for 1D transforms.
void WP2InvTransform(const int32_t src[], WP2TransformType type, uint32_t size,
                     int32_t dst[]);
void WP2Transform(const int32_t src[], WP2TransformType type, uint32_t size,
                  int32_t dst[]);
// If 'reduced' is true, transforms at half the resolution.
void WP2Transform2D(const int32_t src[], WP2TransformType tf_x,
                    WP2TransformType tf_y, uint32_t size_x, uint32_t size_y,
                    int32_t dst[], bool reduced = false);
void WP2InvTransform2D(const int32_t src[], WP2TransformType tf_x,
                       WP2TransformType tf_y, uint32_t size_x, uint32_t size_y,
                       int32_t dst[], bool reduced = false);
// Converts full transform coefficients array to half-transform ones.
// coeffs[] is replaced by its top-left quadrant values divided by 4.
void WP2ReduceCoeffs(const int32_t src[], uint32_t size_x, uint32_t size_y,
                     int32_t dst[]);

// computes (int32)(256.f * sum_i_j{ in[i + 8 * j] * cos_x[i] * cos_y[j] })
// Warning! Do not expect bit-exact same result amongst several implementation,
// because of rounding order of floating-point computations. The SSE version
// is different (off-by-one) from the C-version in ~1% of cases, e.g.
extern int32_t (*WP2SlowDct8x8)(const int32_t in[64],
                                const float cos_x[8], const float cos_y[8]);

//------------------------------------------------------------------------------
// Deblocking filter

namespace WP2 {

static constexpr uint32_t kDblkMaxStrength = 63;
static constexpr uint32_t kDblkMaxSharpness = 7;
// max number of filtered pixels on each side of the edge
static constexpr uint32_t kDblkMaxHalf = 8;

void DblkFilterInit();

// Applies the deblocking filter on a line of 2 * 'half_length' across the
// edge. 'q0' is the top/left-most pixel on the bottom/right side of the edge
// and is 'step' units away from its immediate neighbors.
// The higher 'filter_strength' and the lower 'filter_sharpness' are, the
// higher the variance threshold (under which the filter is applied) is.
// The 'filter_strength' also influences the number of actually filtered
// pixels and how much they are modified (for narrow filter only).
// 'num_precision_bits' includes the sign. 'min' and 'max' are color boundaries.
// Returns true if the line was deblocked, false if pixels were above deblocking
// threshold (no deblocking).
typedef bool (*DeblockLineF)(uint32_t filter_strength, int32_t threshold,
                             uint32_t half_length, bool is_chroma, int32_t min,
                             int32_t max, int32_t step, int16_t* q0);
extern DeblockLineF DeblockLine;
// Same as above but doesn't actually deblock the line, just returns true if
// it would be deblocked when calling DeblockLine.
extern DeblockLineF WouldDeblockLine;

// Returns the number of pixels in [1:half] that are flat on each side of 'q0'.
// For half = 3,4, the samples q0[-4..3] must all be available for read.
// For 4 < half <= 8, the samples q0[-8..7] must be readable.
extern uint32_t (*MeasureFlatLength)(uint32_t filter_sharpness, uint32_t half,
                                     const int16_t* q0, int32_t step);

// Returns the sample threshold deduced from sharpness and precision bits (sign
// included), to be passed to DeblockLineF functions.
int32_t DeblockThresholdFromSharpness(uint32_t filter_sharpness,
                                      uint32_t num_precision_bits);

}  // namespace WP2

//------------------------------------------------------------------------------
// Directional filter

namespace WP2 {

// There are 8 directions, 0 being lower-left/upper-right and 2 being
// horizontal (with top-left origin).
static constexpr uint32_t kDrctFltNumDirs = 8;

// Filtering is applied on blocks of 8x8 pixels.
static constexpr uint32_t kDrctFltSize = 8;

// Kernel size in pixels (around the filtered pixel so 2+1+2 in total).
static constexpr uint32_t kDrctFltTapDist = 2;

// Tap position (x, y) for a given direction [0:7] and a given distance [0:1].
static constexpr int8_t kDrctFltTapPos[kDrctFltNumDirs][kDrctFltTapDist][2] = {
    {{1, -1}, {2, -2}}, {{1, 0}, {2, -1}}, {{1, 0}, {2, 0}}, {{1, 0}, {2, 1}},
    {{1, 1}, {2, 2}},   {{0, 1}, {1, 2}},  {{0, 1}, {0, 2}}, {{0, 1}, {-1, 2}}};

void DrctFilterInit();

// CDEF direction function signature. Section 7.15.2.
// |src| is a pointer to the source block. Pixel size is set by |bitdepth|
// with |step| given in int16_t units. |direction| and |variance| are output
// parameters and must not be nullptr.
typedef void (*CdefDirectionFunc)(const int16_t* src, int32_t step,
                                  uint32_t bitdepth, uint32_t* const direction,
                                  uint32_t* const variance);
extern CdefDirectionFunc CdefDirection4x4;
extern CdefDirectionFunc CdefDirection8x8;

// Fills the extended input buffer with unfiltered pixels or unknown values.
typedef void (*CdefPadFunc)(const int16_t* const src, int32_t src_step,
                            const int16_t* const left, int32_t left_step,
                            const int16_t* const top, int32_t top_step,
                            int32_t width, int32_t height, int32_t n_left,
                            int32_t n_right, int32_t n_top, int32_t n_bottom,
                            int16_t* tmp, int32_t tmp_step);
extern CdefPadFunc CdefPad;

// CDEF filtering function signature. Section 7.15.3.
// |src| is a pointer to the padded input block.
// |block_width|, |block_height| are the width/height in px of the input block.
// |primary_strength|, |secondary_strength|, and |damping| are parameters.
// |direction| is the filtering direction. |dst| is the output buffer.
typedef void (*CdefFilteringFunc)(const int16_t* src, int32_t src_step,
                                  uint32_t bitdepth, int block_width,
                                  int block_height, int primary_strength,
                                  int secondary_strength, int damping,
                                  int direction, int16_t* const dst,
                                  int32_t dst_step);
extern CdefFilteringFunc CdefFiltering;

}  // namespace WP2

//------------------------------------------------------------------------------
// Restoration filter

namespace WP2 {

// Wiener kernel size in pixels (around the filtered pixel so 3+1+3 in total).
// TODO(yguyon): Use a dynamic tap_dist (2 for chroma etc.)
static constexpr uint32_t kWieFltTapDist = 3;
static constexpr uint32_t kWieFltNumTaps = 2 * kWieFltTapDist + 1;
static constexpr uint32_t kWieFltNumBitsTapWgts = 8;  // Including the sign.
static constexpr uint32_t kWieFltNumBitsOverflow = 2;

// Including the sign.
static constexpr uint32_t kWieNumBitsPerTapWgt[kWieFltTapDist]{4, 5, 6};

// Dimensions of the internal buffer for the Wiener filter.
static constexpr uint32_t kWieFltWidth = kMaxTileSize;  // AV1 is 256*1.5 (why?)
static constexpr uint32_t kWieFltHeight = 64;           // Like AV1.

static constexpr uint32_t kWieFltBufWidth = kWieFltWidth + 2 * kWieFltTapDist;
static constexpr uint32_t kWieFltBufHeight = kWieFltHeight + 2 * kWieFltTapDist;
static constexpr uint32_t kWieFltBufSize = kWieFltBufHeight * kWieFltBufWidth;

void WienerFilterInit();

// Fills 'full_tap_weights' from 'half_tap_weights' by symmetry and knowing it
// must give a unit vector.
void WienerHalfToFullWgts(const int32_t half_tap_weights[kWieFltTapDist],
                          int32_t full_tap_weights[kWieFltNumTaps]);

// Applies the Wiener filter on  'dst' of dimensions 'width x height', which
// must fit within the internal buffer (along with tap margins). Available
// surrounding margins of size 'n_*' pixels should be at most 'kWieFltTapDist'.
// 'left' can be nullptr or contain 'n_left x height' pixels. Same for 'right'.
// 'top' can be nullptr or contain 'n_top x (n_left + height + n_right)' pixels.
// Same for 'bottom'. The 'strength_map' must have the same dimension as 'dst'.
// 'num_precision_bits' includes the sign.
typedef WP2Status (*WienerFilterF)(
    uint32_t width, uint32_t height, const int16_t* const left,
    size_t left_step, const int16_t* const right, size_t right_step,
    const int16_t* const top, size_t top_step, const int16_t* const bottom,
    size_t bottom_step, uint32_t n_left, uint32_t n_right, uint32_t n_top,
    uint32_t n_bottom, const int32_t tap_weights_h[kWieFltNumTaps],
    const int32_t tap_weights_v[kWieFltNumTaps],
    const uint8_t* const strength_map, size_t strength_step,
    uint32_t num_precision_bits, size_t dst_step, int16_t* const dst);

extern WienerFilterF WienerFilter;

}  // namespace WP2

//------------------------------------------------------------------------------
// Grain

namespace WP2 {

void GrainFilterInit();

class PseudoRNG;
typedef void (*AddGrainF)(int16_t* samples, size_t step, PseudoRNG* const rng,
                          uint32_t amp, uint32_t cut_off);
typedef void (*GenerateGrainF)(PseudoRNG* const rng,
                               uint32_t amp, uint32_t cut_off, int32_t dst[]);

extern AddGrainF AddGrain4x4;
extern GenerateGrainF GenerateGrain4x4;

}  // namespace WP2

//------------------------------------------------------------------------------
// Predictions

namespace WP2 {

// 'step' is in *dst units, not bytes.
typedef void (*LPredictF)(const int16_t* ctx, uint32_t bw, uint32_t bh,
                          int16_t min_value, int16_t max_value, int16_t* dst,
                          size_t step);
typedef void (*AnglePredictorF)(uint8_t angle, const int16_t* ctx, uint32_t bw,
                                uint32_t bh, int16_t* dst, size_t step);

// generic predictors (dc, smooth, ...)
extern LPredictF BasePredictors[BPRED_LAST];

// Computes a weight based on the distance from the context point (x, y) to the
// direction axis (defined by the cos/sin of its angle, axes X going right and Y
// up, VP9 convention) going through the block point (x0, y0).
// Helper function used by angle predictors.
float ProjWeight(float x0, float y0, float c, float s, float x, float y);

// Generic directional predictors. 'angle_deg' is in degrees.
extern void BaseAnglePredictor(float angle_deg, const int16_t ctx[],
                               uint32_t bw, uint32_t bh, int16_t* dst,
                               size_t step, int16_t min_value,
                               int16_t max_value);
// Simple version of the generic directional predictors: interpolations are
// closer to how AV1 does it.
extern void SimpleAnglePredictor(uint8_t angle_idx, const int16_t ctx[],
                                 uint32_t log2_bw, uint32_t log2_bh,
                                 int16_t* dst, size_t step);

extern void (*AnglePredInterpolate)(const int16_t* src, int32_t frac,
                                    int16_t* dst, uint32_t len);

// Number of angle subdivision within a 22.5 degrees span.
constexpr uint32_t kDirectionalMaxDelta = (2 * kDirectionalMaxAngleDeltaYA + 1);
// number of indexed angled: 8 based angles: {45.0, 67.5, ... 202.5} + 7 deltas
constexpr uint32_t kNumDirectionalAngles = 8 * kDirectionalMaxDelta;
// angle_idx goes for 0 (=35.35.. degrees) to 55 (=212.14 degrees).
// Some notable values:
constexpr uint32_t kAngle_45 = 3;
constexpr uint32_t kAngle_90 = 17;
constexpr uint32_t kAngle_135 = 31;
constexpr uint32_t kAngle_180 = 45;
static_assert(kAngle_90 - kAngle_45 == 2 * kDirectionalMaxDelta,
              "bad kAngle 90-45");
static_assert(kAngle_135 - kAngle_90 == 2 * kDirectionalMaxDelta,
              "bad kAngle 135-90");
static_assert(kAngle_180 - kAngle_135 == 2 * kDirectionalMaxDelta,
              "bad kAngle 180-135");

// Precalculation of 'fuse' weight. 'strength' in [0..1].
static constexpr uint32_t LargeWeightTableDim = 16u;
typedef float LargeWeightTable[LargeWeightTableDim *  LargeWeightTableDim];
extern void PrecomputeLargeWeightTable(float strength, LargeWeightTable table);

// "fuse" weights computation. 'w,h,x,y' in pixels.
extern uint32_t ComputeFuseWeights(uint32_t w, uint32_t h, int32_t x, int32_t y,
                                   const LargeWeightTable table,
                                   float weights[kMaxContextSize]);

// Generic "fuse" prediction. 'bw,bh' in pixels.
extern void BaseFusePredictor(const LargeWeightTable table, const int16_t ctx[],
                              uint32_t bw, uint32_t bh, int16_t* dst,
                              size_t step, int16_t min_value,
                              int16_t max_value);

// Generic "Paeth" prediction.
extern void BasePaethPredictor(const int16_t* ctx, uint32_t bw, uint32_t bh,
                               int16_t min_value, int16_t max_value,
                               int16_t* dst, size_t step);

// Subtract prediction 'pred' from 'src'. 'len' is assumed to be multiple of 4.
// The 'src - pred' difference value is assumed to fit in int16_t.
extern void (*SubtractRow)(const int16_t src[], const int16_t pred[],
                           int32_t dst[], uint32_t len);

// Add prediction 'res' to 'src', clamp and store in 'dst'.
// 'len' is assumed to be multiple of 4.
// 'res[]' and 'res[] + src[]' is assumed to fit in int16_t.
extern void (*AddRow)(const int16_t src[], const int32_t res[],
                      int32_t min, int32_t max, int16_t dst[], uint32_t len);

// Equivalent version for blocks instead of row (same restrictions)
typedef void (*SubtractBlockFunc)(const int16_t src[], uint32_t src_step,
                                  const int16_t pred[], uint32_t pred_step,
                                  int32_t dst[], uint32_t dst_step,
                                  uint32_t height);
typedef void (*AddBlockFunc)(const int16_t src[], uint32_t src_step,
                             const int32_t res[], uint32_t res_step,
                             int32_t min, int32_t max,
                             int16_t dst[], uint32_t dst_step,
                             uint32_t height);

// the array is indexed by log2(width / 4): 0=4, 1=8, 2=16, 3=32
extern SubtractBlockFunc SubtractBlock[4];
extern AddBlockFunc AddBlock[4];

// Will initialize the above functions.
void PredictionInit();

}  // namespace WP2

//------------------------------------------------------------------------------
// Scoring

namespace WP2 {

// returns the min/max value of a 4x4 blocks located at 'src'
extern void (*GetBlockMinMax)(const int16_t* src, uint32_t step,
                              uint32_t num_blocks,
                              int16_t min[], int16_t max[]);

// returns the Min/Max values over a 5x5 block. It is assumed that the whole
// 8x5 block memory is readable.
extern void (*GetBlockMinMax_5x5)(const int16_t* src, uint32_t step,
                                  int16_t* const min, int16_t* const max);

extern void GetBlockMinMaxGeneric(const int16_t* src, uint32_t step,
                                  uint32_t w, uint32_t h,
                                  int16_t* const min, int16_t* const max);

extern void ScoreDspInit();

}  // namespace WP2

//------------------------------------------------------------------------------
// SSIM / PSNR

// total size of the kernel: 2 * kWP2SSIMKernel + 1
static constexpr uint32_t kWP2SSIMKernel = 3;
static constexpr uint32_t kWP2SSIMWeightSum = 16 * 16;   // sum{kWeight}^2

// some implementation (like SSE2) will process more than the kernel's size of
// data in the non-clipped fast version. This value specifies how many
// RGBA samples must be readable when calling WP2SSIMGetFunc(), all
// implementations considered.
static constexpr uint32_t kWP2SSIMMargin = 8;
static_assert(kWP2SSIMMargin >= 2 * kWP2SSIMKernel + 1,
              "invalid kWP2SSIMMargin value");

// struct for accumulating statistical moments
struct WP2DistoStats {
  uint32_t w = 0;             // sum(w_i) : sum of weights
  int32_t xm = 0, ym = 0;     // sum(w_i * x_i), sum(w_i * y_i)
  uint64_t xxm = 0, yym = 0;  // sum(w_i * x_i * x_i), sum(w_i * y_i * y_i)
  int64_t xym = 0;            // sum(w_i * x_i * y_i)
};

// Compute the final SSIM value
extern double WP2SSIMCalculation(uint32_t bit_depth,
                                 const WP2DistoStats& stats);
// This version only computes the (contrast * saturation) part of the SSIM
extern double WP2CsSSIMCalculation(uint32_t bit_depth,
                                   const WP2DistoStats& stats);

// These functions assume an xstep of '4'.
// src1 and src2 are absolute pointers to top-left corner.
typedef void (*WP2SSIMGetClipped8uFunc)(const uint8_t* src1, size_t step1,
                                        const uint8_t* src2, size_t step2,
                                        // center position:
                                        uint32_t xo, uint32_t yo,
                                        // plane dimension:
                                        uint32_t W, uint32_t H,
                                        WP2DistoStats* const stats);
typedef void (*WP2SSIMGetClipped16sFunc)(const int16_t* src1, size_t step1,
                                         const int16_t* src2, size_t step2,
                                         // center position:
                                         uint32_t xo, uint32_t yo,
                                         // plane dimension:
                                         uint32_t W, uint32_t H,
                                         WP2DistoStats* const stats);

// This version is called with the guarantee that you can load 8 bytes and
// 8 rows at offset src1 and src2 (when xstep is '4')
typedef void (*WP2SSIMGet8uFunc)(const uint8_t* src1, size_t step1,
                                 const uint8_t* src2, size_t step2,
                                 WP2DistoStats* const stats);
typedef void (*WP2SSIMGet16sFunc)(const int16_t* src1, size_t step1,
                                  const int16_t* src2, size_t step2,
                                  WP2DistoStats* const stats);

extern WP2SSIMGet8uFunc WP2SSIMGet4x8u;                // unclipped / unchecked
extern WP2SSIMGet8uFunc WP2SSIMGet8u;                  // unclipped / unchecked
extern WP2SSIMGet16sFunc WP2SSIMGet16s;                // unclipped / unchecked
extern WP2SSIMGetClipped8uFunc WP2SSIMGetClipped4x8u;  // with clipping
extern WP2SSIMGetClipped8uFunc WP2SSIMGetClipped8u;    // with clipping
extern WP2SSIMGetClipped16sFunc WP2SSIMGetClipped16s;  // with clipping

void WP2SSIMInit();

//------------------------------------------------------------------------------
// alpha-related functions

// returns true if src[0..len) contains 'value'
extern bool (*WP2HasValue8b)(const uint8_t* src, size_t len, uint8_t value);
extern bool (*WP2HasValue16b)(const int16_t* src, size_t len, int16_t value);

// returns true if src[0..len) contains other values than 'value'
extern bool (*WP2HasOtherValue8b)(
    const uint8_t* src, size_t len, uint8_t value);
extern bool (*WP2HasOtherValue16b)(
    const int16_t* src, size_t len, int16_t value);

// return true if src[4 * i + 0] is different from 'value'
extern bool (*WP2HasOtherValue8b32b)(
    const uint8_t* src, size_t len, uint8_t value);

void WP2AlphaInit();

//------------------------------------------------------------------------------

// Return the sum of squared error for a continuous row of samples.
extern uint64_t (*WP2SumSquaredError8u)(const uint8_t* src1,
                                        const uint8_t* src2,
                                        uint32_t len);
// the int16_t in src1[] and src2[] are assumed in range [-1024, 1023]
extern uint64_t (*WP2SumSquaredError16s)(const int16_t* src1,
                                         const int16_t* src2,
                                         uint32_t len);
// Accumulate 32b samples channels by channels. result[] is not reset to 0.
extern void (*WP2SumSquaredError4x8u)(const uint8_t* src1,
                                      const uint8_t* src2,
                                      uint32_t len, uint64_t result[4]);
// Accumulate 64b samples channels by channels. result[] is not reset to 0.
extern void (*WP2SumSquaredError2x16u)(const uint16_t* src1,
                                       const uint16_t* src2, uint32_t len,
                                       uint64_t result[4]);
// returns the SSE for a w x h block.
// 16b input is assumed to be in range [-1024,1023].
extern uint64_t (*WP2SumSquaredErrorBlock)(const int16_t* src1, uint32_t step1,
                                           const int16_t* src2, uint32_t step2,
                                           uint32_t w, uint32_t h);

void WP2PSNRInit();

//------------------------------------------------------------------------------

#endif  /* WP2_DSP_DSP_H_ */
