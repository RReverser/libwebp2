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
// AOM / AV1 wrappers
//
// Author: Skal (pascal.massimino@gmail.com)

#include <cstdint>
#include <cstring>
#include <cmath>

#include "extras/aom_utils.h"
#include "extras/ccsp_imageio.h"
#include "src/utils/utils.h"
#include "src/wp2/decode.h"
#include "src/wp2/encode.h"
#include "imageio/imageio_util.h"
#include "examples/example_utils.h"

#ifdef WP2_HAVE_AOM
#include "aom/aom_decoder.h"
#include "aom/aom_encoder.h"
#include "aom/aomcx.h"
#include "aom/aomdx.h"
#if defined(WP2_HAVE_AOM_DBG)
#include "stats/aomstats.h"
// to benefit from these, one should build libaom using the config:
//  cmake ../ -DCONFIG_ACCOUNTING=1 -DCONFIG_INSPECTION=1
#include "av1/decoder/accounting.h"
#include "av1/decoder/inspection.h"
#include "av1/common/common_data.h"
#endif
#include "examples/stopwatch.h"
#if !defined(WP2_HAVE_WEBP)
#error "WEBP support is required by WP2_HAVE_AOM"
#else
#include "webp/decode.h"
#include "webp/encode.h"
#endif
#endif  // WP2_HAVE_AOM

#ifdef WP2_HAVE_AVIF
#include "avif/avif.h"
#endif

//------------------------------------------------------------------------------

#if defined(WP2_HAVE_AOM)
#if (AOM_ENCODER_ABI_VERSION == (8 + 1 + AOM_CODEC_ABI_VERSION))
#define USE_OLD_AOM 1
#endif
static bool encode_frame(aom_codec_ctx_t* const codec, aom_image_t* const img,
                         int pts, aom_enc_frame_flags_t flags,
                         std::string* const out, bool report_error = true) {
  const aom_codec_err_t err =
      aom_codec_encode(codec, img, pts, /*duration=*/1, flags);
  if (err != AOM_CODEC_OK) {
    if (report_error) {
      fprintf(stderr, "aom_codec_encode() failed ['%s']\n",
              aom_codec_err_to_string(err));
    }
    return false;
  }

  aom_codec_iter_t iter = NULL;
  const aom_codec_cx_pkt_t* pkt = NULL;
  bool got_pkts = false;
  while ((pkt = aom_codec_get_cx_data(codec, &iter)) != NULL) {
    got_pkts = true;
    if (pkt->kind == AOM_CODEC_CX_FRAME_PKT) {
      out->append((const char*)pkt->data.frame.buf, pkt->data.frame.sz);
    }
  }
  return got_pkts;
}

#endif  // WP2_HAVE_AOM

#if defined(WP2_HAVE_AOM_DBG)

namespace {

void inspect_cb(void* pbi, void* data) {
  insp_frame_data* const frame_data = (insp_frame_data*)data;
  ifd_inspect(frame_data, pbi, /*skip_non_transform=*/0);
}

// Given a 'mi_data_name', returns the value of the variable with that name in
// the instance 'mi' of the 'insp_mi_data' struct.
int16_t GetMiData(const insp_mi_data& mi, const std::string& mi_data_name) {
  if (mi_data_name == "mode") return mi.mode;
  if (mi_data_name == "uv_mode") return mi.uv_mode;
#if defined(USE_OLD_AOM)
  if (mi_data_name == "sb_type") return mi.sb_type;
#else
  if (mi_data_name == "sb_type") return mi.bsize;
#endif
  if (mi_data_name == "skip") return mi.skip;
  if (mi_data_name == "segment_id") return mi.segment_id;
  if (mi_data_name == "dual_filter_type") return mi.dual_filter_type;
  if (mi_data_name == "filter[0]") return mi.filter[0];
  if (mi_data_name == "filter[1]") return mi.filter[1];
  if (mi_data_name == "tx_type") return mi.tx_type;
  if (mi_data_name == "tx_size") return mi.tx_size;
  if (mi_data_name == "cdef_level") return mi.cdef_level;
  if (mi_data_name == "cdef_strength") return mi.cdef_strength;
  if (mi_data_name == "cfl_alpha_idx") return mi.cfl_alpha_idx;
  if (mi_data_name == "cfl_alpha_sign") return mi.cfl_alpha_sign;
  if (mi_data_name == "current_qindex") return mi.current_qindex;
  if (mi_data_name == "compound_type") return mi.compound_type;
  if (mi_data_name == "motion_mode") return mi.motion_mode;
  if (mi_data_name == "intrabc") return mi.intrabc;
  if (mi_data_name == "palette") return mi.palette;
  if (mi_data_name == "uv_palette") return mi.uv_palette;
  assert(false);
  return 0;
}

// Draws a 3x3 red square representing the 'value' as 9 bits.
void DrawNumber(uint32_t value, int x, int y, WP2::ArgbBuffer* const out) {
  for (uint32_t j = 0; j < 3; ++j) {
    for (uint32_t i = 0; i < 3; ++i) {
      if (x + i < out->width && y + j < out->height) {
        ((uint32_t*)out->GetRow(y + j))[(x + i)] =
            (value & 1) ? 0xffff8888u : 0xffff0000u;
      }
      value >>= 1;
    }
  }
}

// Fills rectangles with uniform colors representing the values of
// 'mi_data_name' in 'data' blocks.
void DrawMiData(const insp_frame_data& data, const std::string& mi_data_name,
                bool draw_mi_number, WP2::ArgbBuffer* const out) {
  const uint32_t num_cols = data.mi_cols;
  const uint32_t num_rows = data.mi_rows;
  assert(out != nullptr);

  const auto minmax = std::minmax_element(
      data.mi_grid, data.mi_grid + num_cols * num_rows,
      [&](const insp_mi_data& lhs, const insp_mi_data& rhs) {
        return (GetMiData(lhs, mi_data_name) < GetMiData(rhs, mi_data_name));
      });
  const int16_t min = GetMiData(*minmax.first, mi_data_name);
  const int16_t max = GetMiData(*minmax.second, mi_data_name);

  // To remember already painted areas.
  std::vector<int32_t> values(num_cols * num_rows, max + 1);

  for (uint32_t r = 0; r < num_rows && r * 4 < out->height; ++r) {
    for (uint32_t c = 0; c < num_cols && c * 4 < out->width; ++c) {
      const insp_mi_data& mi = data.mi_grid[r * num_cols + c];
      const int16_t mi_data = GetMiData(mi, mi_data_name);

      const uint8_t intensity = (mi_data - min) * 0xff / std::max(1, max - min);
      const uint32_t color = 0xff000000u + 0x010101u * intensity;

      if (values[r * num_cols + c] == max + 1) {  // Not yet drawn.
#if defined(USE_OLD_AOM)
        const BLOCK_SIZE bsize = (BLOCK_SIZE)mi.sb_type;
#else
        const BLOCK_SIZE bsize = (BLOCK_SIZE)mi.bsize;
#endif
        const uint32_t w = mi_size_wide[bsize], h = mi_size_high[bsize];
        const uint32_t min_x = c * 4, min_y = r * 4;
        out->Fill({min_x, min_y, w * 4, h * 4},
                  WP2::Argb32b{0xff, intensity, intensity, intensity});

        if (draw_mi_number) {
          DrawNumber((uint32_t)(mi_data - min), min_x + 1, min_y + 1, out);
        }

        for (uint32_t sub_r = r; sub_r < std::min(num_rows, r + h); ++sub_r) {
          for (uint32_t sub_c = c; sub_c < std::min(num_cols, c + w); ++sub_c) {
            values[sub_r * num_cols + sub_c] = mi_data;
          }
        }
      } else if (mi_data != values[r * num_cols + c]) {
        // Already drawn but non-uniform values were found within a block.
        const uint32_t min_x = c * 4, min_y = r * 4;
        const uint32_t max_x = std::min((c + 1) * 4, out->width);
        const uint32_t max_y = std::min((r + 1) * 4, out->height);

        for (uint32_t y = min_y; y < max_y; ++y) {
          for (uint32_t x = min_x; x < max_x; ++x) {
            if ((x + y) & 1) continue;  // Semi-transparent checkerboard.
            ((uint32_t*)out->GetRow(y))[x] = color;
          }
        }
      }
    }
  }
}

void DrawQuant(const insp_frame_data& data, int channel,
               WP2::ArgbBuffer* const out) {
  const uint32_t num_cols = data.mi_cols;
  const uint32_t num_rows = data.mi_rows;
  assert(out != nullptr);

  for (uint32_t r = 0; r < num_rows && r * 4 < out->height; ++r) {
    for (uint32_t c = 0; c < num_cols && c * 4 < out->width; ++c) {
      const insp_mi_data& mi = data.mi_grid[r * num_cols + c];
      const int16_t segment = GetMiData(mi, "segment_id");
      int16_t quant = (channel == 0)
                          ? data.y_dequant[segment][0]
                          : (channel == 1) ? data.u_dequant[segment][0]
                                           : data.v_dequant[segment][0];
      // Raw quant values are multiplied by 4. We shift them back to
      // make them comparable with those used by wp2.
      quant /= 4;

      // Cap quant at roughly the same as kMaxQuantStep.
      constexpr int16_t kCap = 8192;
      const float normalized_quant = std::min(quant, kCap) / (float)kCap;
      // Same conversion as in vdebug.cc
      constexpr float kGamma = 2.2f;
      const float red = std::pow(normalized_quant, 1.f / kGamma);
      const float green = std::pow((1.f - normalized_quant), 1.f / kGamma);
      // Max out saturation.
      const float mult = std::min(1.f / red, 1.f / green);
      WP2::Argb32b color;
      color.r = std::min((int32_t)(red * mult * 0xff), 0xff);
      color.g = std::min((int32_t)(green * mult * 0xff), 0xff);
      color.b = 0;
      color.a = 0xff;
      out->Fill({c * 4, r * 4, 4, 4}, color);
    }
  }
}

std::vector<WP2::Rectangle> ExtractBoxes(const insp_frame_data& data,
                                         uint32_t width, uint32_t height,
                                         bool extract_transform_boxes) {
  std::vector<WP2::Rectangle> boxes;
  const uint32_t w = (uint32_t)data.mi_cols;
  const uint32_t h = (uint32_t)data.mi_rows;
  std::vector<bool> done(w * h, false);
  for (uint32_t y = 0; y < h && y * 4u < height; ++y) {
    for (uint32_t x = 0; x < w && x * 4u < width; ++x) {
      if (done[y * w + x]) continue;

      const insp_mi_data& mi = data.mi_grid[y * w + x];
#if defined(USE_OLD_AOM)
      const BLOCK_SIZE bsize = (BLOCK_SIZE)mi.sb_type;
#else
      const BLOCK_SIZE bsize = (BLOCK_SIZE)mi.bsize;
#endif
      const uint32_t bw = mi_size_wide[bsize];
      const uint32_t bh = mi_size_high[bsize];

      if (extract_transform_boxes) {
        // For lossy intra luma only. See av1_loopfilter.c, get_transform_size()
        const TX_SIZE tx_size = (TX_SIZE)mi.tx_size;
        const uint32_t tw_px = tx_size_wide[tx_size];
        const uint32_t th_px = tx_size_high[tx_size];
        const uint32_t tw = tw_px / 4, th = th_px / 4;

        for (uint32_t ty = 0; ty < bh && (y + ty) * 4u < height; ty += th) {
          for (uint32_t tx = 0; tx < bw && (x + tx) * 4u < width; tx += tw) {
            const uint32_t tx_px = (x + tx) * 4u, ty_px = (y + ty) * 4u;
            boxes.emplace_back(tx_px, ty_px,  // Transform position/dimensions
                               std::min(tw_px, width - tx_px),
                               std::min(th_px, height - ty_px));
          }
        }
      } else {
        const uint32_t x_px = x * 4u, y_px = y * 4u;
        const uint32_t bw_px = bw * 4u, bh_px = bh * 4u;
        boxes.emplace_back(x_px, y_px,  // Block position/dimensions
                           std::min(bw_px, width - x_px),
                           std::min(bh_px, height - y_px));
      }

      for (uint32_t yd = y; yd < std::min(h, y + bh); ++yd) {
        for (uint32_t xd = x; xd < std::min(w, x + bw); ++xd) {
          done[yd * w + xd] = true;
        }
      }
    }
  }
  return boxes;
}

void DrawBoxes(const std::vector<WP2::Rectangle>& boxes,
               WP2::Argb32b color, WP2::ArgbBuffer* const out) {
  assert(out != nullptr);
  for (const auto& rect : boxes) out->DrawRect(rect, color);
}

// category id, similar to those reported by 'dwebp'
static const char* const kCtgCoeffs        = "coeffs";
static const char* const kCtgPaletteCoeffs = "palette-coeffs";
static const char* const kCtgPredModes     = "pred-modes";
static const char* const kCtgSegmentId     = "segment-id";
static const char* const kCtgHeader        = "header";
static const char* const kCtgPaletteHeader = "palette-header";
static const char* const kCtgBlockSize     = "block-size";
static std::map<std::string, const char*> kCategories = {
  // coeffs
  {"read_coeffs_reverse_2d",          kCtgCoeffs},
  {"av1_read_coeffs_txb",             kCtgCoeffs},
  {"read_skip_txfm",                  kCtgCoeffs},
  {"read_coeffs_reverse",             kCtgCoeffs},
  {"read_golomb",                     kCtgCoeffs},
  // pred-modes
  {"read_intra_mode",                 kCtgPredModes},
  {"read_intra_mode_uv",              kCtgPredModes},
  {"read_angle_delta",                kCtgPredModes},
  {"read_filter_intra_mode_info",     kCtgPredModes},
  // segment-id
  {"read_segment_id",                 kCtgSegmentId},
  // block-size
  {"av1_read_tx_type",                kCtgBlockSize},
  {"read_selected_tx_size",           kCtgBlockSize},
  {"read_partition",                  kCtgHeader},
  {"read_skip",                       kCtgHeader},
  // everything about filters
  {"read_wiener_filter",              kCtgHeader},
  {"loop_restoration_read_sb_coeffs", kCtgHeader},
  {"read_cdef",                       kCtgHeader},
  {"cfl:signs",                       kCtgHeader},
  {"cfl:alpha_u",                     kCtgHeader},
  {"cfl:alpha_v",                     kCtgHeader},
  {"read_sgrproj_filter",             kCtgHeader},
  // palette
  {"decode_color_map_tokens",         kCtgPaletteCoeffs},
  {"av1_read_uniform",                kCtgPaletteHeader},
  {"read_palette_colors_y",           kCtgPaletteHeader},
  {"read_palette_colors_uv",          kCtgPaletteHeader},
  {"read_palette_mode_info",          kCtgPaletteHeader},
};

}  // namespace

#endif  // WP2_HAVE_AOM_DBG

namespace WP2 {

// small class to wrap plain-C destruction call
template<class T> class Cleaner {
 public:
  Cleaner(T* const obj, void (*dtor)(T*)) : obj_(obj), dtor_(dtor) {}
  ~Cleaner() { dtor_(obj_); }
  T* const obj_;
  void (*dtor_)(T*);
};

#if defined(WP2_HAVE_AOM)

//------------------------------------------------------------------------------

// From twopass_encoder.c / aomenc.c
static WP2Status wp2_get_frame_stats(aom_codec_ctx_t* ctx,
                                     const aom_image_t* img,
                                     aom_codec_pts_t pts, unsigned int duration,
                                     aom_enc_frame_flags_t flags,
                                     aom_fixed_buf_t* stats,
                                     bool* const got_stats) {
  *got_stats = false;
  aom_codec_iter_t iter = NULL;
  const aom_codec_cx_pkt_t* pkt = NULL;
  const aom_codec_err_t res = aom_codec_encode(ctx, img, pts, duration, flags);
  if (res != AOM_CODEC_OK) {
    fprintf(stderr, "Failed to get frame stats.\n");
    return WP2_STATUS_BAD_WRITE;
  }

  while ((pkt = aom_codec_get_cx_data(ctx, &iter)) != NULL) {
    *got_stats = true;

    if (pkt->kind == AOM_CODEC_STATS_PKT) {
      const uint8_t* const pkt_buf =
          (const uint8_t*)pkt->data.twopass_stats.buf;
      const size_t pkt_size = pkt->data.twopass_stats.sz;
      stats->buf = realloc(stats->buf, stats->sz + pkt_size);
      memcpy((uint8_t*)stats->buf + stats->sz, pkt_buf, pkt_size);
      stats->sz += pkt_size;
    }
  }

  return WP2_STATUS_OK;
}

//------------------------------------------------------------------------------

// Converts 'ref' image into AV1 struct 'raw'.
static WP2Status CreateAV1Img(const ArgbBuffer& ref, aom_image_t* const raw,
                              bool use_yuv444 = false) {
  YUVPlane yuv;
  WP2_CHECK_STATUS(ToYCbCr(ref, /*ycbcr_num_bits=*/8,
                           use_yuv444 ? nullptr : &SamplingTaps::kDownSharp,
                           &yuv));

  const int w = ref.width;
  const int h = ref.height;
  if (!aom_img_alloc(raw, use_yuv444 ? AOM_IMG_FMT_I444 : AOM_IMG_FMT_I420,
                     w, h, /*align=*/16)) {
    fprintf(stderr, "aom_img_alloc() failed!");
    return WP2_STATUS_OUT_OF_MEMORY;
  }
  // Match y4m format with luma in [16:235] range.
  yuv.Y.To(raw->planes[AOM_PLANE_Y], raw->stride[AOM_PLANE_Y], 16);
  yuv.U.To(raw->planes[AOM_PLANE_U], raw->stride[AOM_PLANE_U], 128);
  yuv.V.To(raw->planes[AOM_PLANE_V], raw->stride[AOM_PLANE_V], 128);
  if (ref.HasTransparency()) {
    fprintf(stderr, "# WARNING! Alpha channel is discarded with AV1!!\n");
  }

  return WP2_STATUS_OK;
}

//------------------------------------------------------------------------------

// Initializes the 'codec'. Uses 'params', 'pass' and 'stats' to configure it.
static WP2Status CreateCodec(int w, int h, const ParamsAV1& params,
                             aom_enc_pass pass, aom_fixed_buf_t stats,
                             aom_codec_ctx_t* const codec) {
  const aom_codec_iface_t* const itf = aom_codec_av1_cx();
  if (itf == NULL) {
    fprintf(stderr, "aom_codec_av1_cx() failed!\n");
    return WP2_STATUS_INVALID_PARAMETER;
  }

  aom_codec_enc_cfg_t cfg;
  if (aom_codec_enc_config_default(itf, &cfg, 0) != AOM_CODEC_OK) {
    fprintf(stderr, "aom_codec_enc_config_default() failed!\n");
    return WP2_STATUS_INVALID_PARAMETER;
  }
  cfg.g_w = w;
  cfg.g_h = h;
  cfg.g_timebase.num = 1;
  cfg.g_timebase.den = 25;
  cfg.g_bit_depth = AOM_BITS_8;
  cfg.g_usage =
      (params.speed == 0) ? AOM_USAGE_REALTIME : AOM_USAGE_GOOD_QUALITY;
  cfg.g_profile = (params.use_yuv444 ? 1 : 0);
  cfg.g_threads = params.threads;
  cfg.rc_end_usage = AOM_Q;
  cfg.g_pass = pass;
  cfg.rc_twopass_stats_in = stats;

  if (aom_codec_enc_init(codec, itf, &cfg, 0) != AOM_CODEC_OK) {
    fprintf(stderr, "aom_codec_enc_init() failed!\n");
    return WP2_STATUS_OUT_OF_MEMORY;
  }

  if (params.quality > 95) {
    if (aom_codec_control(codec, AV1E_SET_LOSSLESS, 1) != AOM_CODEC_OK) {
      fprintf(stderr, "aom_codec_control(AV1E_SET_LOSSLESS) failed\n");
      return WP2_STATUS_INVALID_PARAMETER;
    }
  } else {
    const int level = std::lround(63 - params.quality * 63 / 95);
    if (aom_codec_control(codec, AOME_SET_CQ_LEVEL, level) != AOM_CODEC_OK) {
      fprintf(stderr, "aom_codec_control(AOME_SET_CQ_LEVEL) failed\n");
      return WP2_STATUS_INVALID_PARAMETER;
    }
  }
  // The only metrics AV1 can be tuned for are PSNR, SSIM or VMAF.
  WP2_CHECK_OK(params.tuning == PSNR || params.tuning == SSIM,
               WP2_STATUS_INVALID_CONFIGURATION);
  const int tuning = (params.tuning == PSNR) ? 0 : 1;
  if (aom_codec_control(codec, AOME_SET_TUNING, tuning) != AOM_CODEC_OK) {
    fprintf(stderr, "aom_codec_control(AOME_SET_TUNING) failed\n");
    return WP2_STATUS_INVALID_PARAMETER;
  }
  // AV1's speed param goes opposite to WP2's
  const int speed = std::max(8 - (int)params.speed * 8 / 9, 0);
  if (aom_codec_control(codec, AOME_SET_CPUUSED, speed) != AOM_CODEC_OK) {
    fprintf(stderr,
            "aom_codec_control(AOME_SET_CPUUSED) failed for speed setting\n");
    return WP2_STATUS_INVALID_PARAMETER;
  }
  return WP2_STATUS_OK;
}

//------------------------------------------------------------------------------

static constexpr int kFlags = AOM_EFLAG_FORCE_KF;

static WP2Status FirstPassAV1(const ArgbBuffer& ref, const ParamsAV1& params,
                              double* const timing,
                              aom_fixed_buf_t* const stats) {
  aom_image_t raw;
  WP2_CHECK_STATUS(CreateAV1Img(ref, &raw, params.use_yuv444));
  const auto C1 = Cleaner<aom_image>(&raw, aom_img_free);

  aom_codec_ctx_t codec;
  // TODO(yguyon): Check possible leaks of itf/cfg/codec
  WP2_CHECK_STATUS(CreateCodec((int)ref.width, (int)ref.height, params,
                               AOM_RC_FIRST_PASS, /*stats=*/{NULL, 0}, &codec));

  const double start_time = GetStopwatchTime();
  bool got_stats;
  // Calculate frame statistics.
  WP2_CHECK_STATUS(wp2_get_frame_stats(&codec, &raw, /*pts=*/0, /*duration=*/1,
                                       kFlags, stats, &got_stats));
  // Flush encoder.
  do {
    WP2_CHECK_STATUS(wp2_get_frame_stats(
        &codec, NULL, /*pts=*/0, /*duration=*/1, kFlags, stats, &got_stats));
  } while (got_stats);
  if (timing != nullptr) *timing = GetStopwatchTime() - start_time;
  if (aom_codec_destroy(&codec) != AOM_CODEC_OK) {
    fprintf(stderr, "aom_codec_destroy() failed!\n");
    return WP2_STATUS_BITSTREAM_ERROR;
  }
  return WP2_STATUS_OK;
}

//------------------------------------------------------------------------------

// Encodes 'ref' into 'out', decodes it into 'decoded'.
static WP2Status CompressAV1(const ArgbBuffer& ref, const ParamsAV1& params,
                             aom_fixed_buf_t stats, ArgbBuffer* const decoded,
                             std::string* const out, double timing[],
                             std::vector<WP2::Rectangle>* const blocks,
                             std::vector<WP2::Rectangle>* const transforms,
                             QuantAV1* const quant) {
  (void)blocks, (void)transforms, (void)quant;
  WP2_CHECK_OK(decoded != nullptr, WP2_STATUS_NULL_PARAMETER);
  WP2_CHECK_OK(out != nullptr, WP2_STATUS_NULL_PARAMETER);
  out->clear();

  // encode
  {
    aom_image_t raw;
    WP2_CHECK_STATUS(CreateAV1Img(ref, &raw, params.use_yuv444));
    auto C1 = Cleaner<aom_image>(&raw, aom_img_free);

    aom_codec_ctx_t codec;
    // TODO(yguyon): Check possible leaks of itf/cfg/codec
    WP2_CHECK_STATUS(
        CreateCodec((int)ref.width, (int)ref.height, params,
                    (stats.buf != NULL) ? AOM_RC_LAST_PASS : AOM_RC_ONE_PASS,
                    stats, &codec));
    // encode
    const double start_time = GetStopwatchTime();
    if (encode_frame(&codec, &raw, 0 /* frame count */, kFlags, out, false)) {
      fprintf(stderr, "!?! first calls to encode_frame() should 'fail'!\n");
      return WP2_STATUS_BITSTREAM_ERROR;
    }
    while (encode_frame(&codec, NULL, -1, 0, out)) { /* flush encoder */ }
    if (timing != nullptr) timing[0] = GetStopwatchTime() - start_time;
    if (aom_codec_destroy(&codec) != AOM_CODEC_OK) {
      fprintf(stderr, "aom_codec_destroy() failed!\n");
      return WP2_STATUS_BITSTREAM_ERROR;
    }
  }

  // decode back
  WP2_CHECK_STATUS(decoded->Resize(ref.width, ref.height));
  {
    const aom_codec_iface_t* const itf = aom_codec_av1_dx();
    if (itf == NULL) {
      fprintf(stderr, "aom_codec_av1_dx() failed!\n");
      return WP2_STATUS_INVALID_PARAMETER;
    }

    aom_codec_ctx_t codec;
    if (aom_codec_dec_init(&codec, itf, NULL, 0) != AOM_CODEC_OK) {
      fprintf(stderr, "aom_codec_dec_init() failed\n");
      return WP2_STATUS_OUT_OF_MEMORY;
    }

#if defined(WP2_HAVE_AOM_DBG)
    insp_frame_data frame_data;
    ifd_init(&frame_data, ref.width, ref.height);
    WP2_CHECK_OK(frame_data.mi_grid != NULL, WP2_STATUS_OUT_OF_MEMORY);
    auto C4 = Cleaner<insp_frame_data>(&frame_data, ifd_clear);

    struct aom_inspect_init inspect = { inspect_cb, &frame_data };
    const bool use_inspect =
        (params.draw_mi_data != nullptr || params.draw_blocks ||
         params.draw_transforms || blocks != nullptr || transforms != nullptr ||
         quant != nullptr) &&
        (aom_codec_control(&codec, AV1_SET_INSPECTION_CALLBACK, &inspect) ==
         AOM_CODEC_OK);
#endif    // WP2_HAVE_AOM_DBG

    const double start_time = GetStopwatchTime();
    aom_codec_dec_cfg cfg;
    memset(&cfg, 0, sizeof(cfg));
    cfg.allow_lowbitdepth = 1;
    if (aom_codec_decode(&codec,
                         (const uint8_t*)out->data(), out->size(), &cfg)
          != AOM_CODEC_OK) {
      fprintf(stderr, "aom_codec_decode() failed.\n");
      return WP2_STATUS_OUT_OF_MEMORY;
    }
    if (timing != nullptr) timing[1] = GetStopwatchTime() - start_time;
    aom_image_t* img = NULL;
    aom_codec_iter_t iter = NULL;
    bool ok = false;
    while ((img = aom_codec_get_frame(&codec, &iter)) != NULL) {
      YUVPlane yuv;
      WP2_CHECK_STATUS(yuv.Resize(ref.width, ref.height, /*pad=*/1,
                                  /*has_alpha=*/false, !params.use_yuv444));
      yuv.Y.From(img->planes[AOM_PLANE_Y], img->stride[AOM_PLANE_Y], -16);
      yuv.U.From(img->planes[AOM_PLANE_U], img->stride[AOM_PLANE_U], -128);
      yuv.V.From(img->planes[AOM_PLANE_V], img->stride[AOM_PLANE_V], -128);
      WP2_CHECK_STATUS(
          ToRGBA(yuv, params.use_yuv444 ? nullptr : &SamplingTaps::kUpSmooth,
                 decoded));
      ok = true;
      break;
    }
    if (!ok) {
      fprintf(stderr, "Round-trip error!\n");
      return WP2_STATUS_OUT_OF_MEMORY;
    }
#if defined(WP2_HAVE_AOM_DBG)
#if defined(WP2_BITTRACE)
    if (params.bittrace > 0) {
      Accounting* acc = NULL;
      if (aom_codec_control(&codec, AV1_GET_ACCOUNTING, &acc) == AOM_CODEC_OK) {
        std::map<std::string, WP2BitCounts::mapped_type> accum;
        const uint32_t num_syms = acc->syms.num_syms;
        for (uint32_t i = 0; i < num_syms; ++i) {
          const AccountingSymbol& sym = acc->syms.syms[i];
          const double bits = sym.bits / 8.;
          const char* str = acc->syms.dictionary.strs[sym.id];
          const uint32_t num_samples = sym.samples;
          if (params.compact_trace) {
            if (kCategories[str] == nullptr) {
              fprintf(stderr, "Unmapped category! Missing line:\n");
              fprintf(stderr, "  {\"%s\",    kCtgXXX},\n", str);
            } else {
              str = kCategories[str];
            }
          }
          accum[str].bits += bits;
          accum[str].num_occurrences += num_samples;
        }
        std::vector<WP2TraceType> values;
        for (const auto& p : accum) values.push_back(p);
        const bool use_bytes = (params.bittrace == 2);
        PrintBitTraceCluster(values, true, use_bytes, 0);
      } else {
        fprintf(stderr, "Error! AV1 accounting (bit traces) not compiled in\n");
      }
    }
#endif  // WP2_BITTRACE
    if (use_inspect) {
      if (params.draw_mi_data != nullptr) {
        if (std::strcmp(params.draw_mi_data, "quant") == 0) {
          DrawQuant(frame_data, /*channel=*/0, decoded);
        } else {
          DrawMiData(frame_data, params.draw_mi_data, params.draw_mi_number,
                     decoded);
        }
      }
      // Print transforms before blocks because transforms are smaller.
      if (params.draw_transforms || transforms != nullptr) {
        std::vector<WP2::Rectangle> t =
            ExtractBoxes(frame_data, ref.width, ref.height,
                         /*extract_transform_boxes=*/true);
        if (params.draw_transforms) {
          DrawBoxes(t, {0xff, 0x08, 0x80, 0xf0}, decoded);
        }
        if (transforms != nullptr) std::swap(*transforms, t);
      }
      if (params.draw_blocks || blocks != nullptr) {
        std::vector<WP2::Rectangle> b =
            ExtractBoxes(frame_data, ref.width, ref.height,
                         /*extract_transform_boxes=*/false);
        if (params.draw_blocks) DrawBoxes(b, {0xff, 0x08, 0xf0, 0x80}, decoded);
        if (blocks != nullptr) std::swap(*blocks, b);
      }
      if (quant != nullptr) {
        for (uint32_t i = 0; i < 8; ++i) {
          for (uint32_t j = 0; j < 2; ++j) {
            // Raw quant values are multiplied by 4. We shift them back to
            // make them comparable with those used by wp2.
            quant->y_dequant[i][j] = frame_data.y_dequant[i][j] >> 2;
            quant->u_dequant[i][j] = frame_data.u_dequant[i][j] >> 2;
            quant->v_dequant[i][j] = frame_data.v_dequant[i][j] >> 2;
          }
        }
      }
    }
#endif  // WP2_HAVE_AOM_DBG

    if (aom_codec_destroy(&codec) != AOM_CODEC_OK) {
      fprintf(stderr, "aom_codec_destroy() failed!\n");
      return WP2_STATUS_OUT_OF_MEMORY;
    }
  }

  return WP2_STATUS_OK;
}

//------------------------------------------------------------------------------

WP2Status CompressAV1(const ArgbBuffer& ref, const ParamsAV1& params,
                      ArgbBuffer* const decoded, std::string* const out,
                      double timing[],
                      std::vector<WP2::Rectangle>* const blocks,
                      std::vector<WP2::Rectangle>* const transforms,
                      QuantAV1* const quant) {
  if (params.pass < 2) {
    return CompressAV1(ref, params, /*stats=*/{NULL, 0}, decoded, out, timing,
                       blocks, transforms, quant);
  } else {
    aom_fixed_buf_t stats = {NULL, 0};
    double first_pass_timing;
    WP2_CHECK_STATUS(FirstPassAV1(ref, params, &first_pass_timing, &stats));
    WP2_CHECK_STATUS(CompressAV1(ref, params, stats, decoded, out, timing,
                                 blocks, transforms, quant));
    if (timing != nullptr) timing[0] += first_pass_timing;
    free(stats.buf);
    return WP2_STATUS_OK;
  }
}

#else

WP2Status CompressAV1(const ArgbBuffer& ref, const ParamsAV1& params,
                      ArgbBuffer* const decoded, std::string* const out,
                      double timing[],
                      std::vector<WP2::Rectangle>* const blocks,
                      std::vector<WP2::Rectangle>* const transforms,
                      QuantAV1* const quant) {
  (void)ref;
  (void)params;
  (void)decoded;
  (void)out;
  (void)timing;
  (void)blocks;
  (void)transforms;
  (void)quant;
  fprintf(stderr, "CompressAV1: AOM support not enabled at compile time.\n");
  return WP2_STATUS_UNSUPPORTED_FEATURE;
}

#endif  // WP2_HAVE_AOM

#ifdef WP2_HAVE_AVIF

static float ToAvifQuantizer(float quality) {
  // Note AVIF_QUANTIZER_LOSSLESS doesn't actually seem lossless.
  if (quality > 95) return AVIF_QUANTIZER_LOSSLESS;
  return AVIF_QUANTIZER_WORST_QUALITY +
         quality *
             (AVIF_QUANTIZER_BEST_QUALITY - AVIF_QUANTIZER_WORST_QUALITY) / 100;
}

WP2Status CompressAVIF(const ArgbBuffer& ref, const ParamsAV1& params,
                       ArgbBuffer* const decoded, std::string* const out,
                       double timing[2]) {
  WP2_CHECK_OK(decoded != nullptr, WP2_STATUS_NULL_PARAMETER);
  WP2_CHECK_OK(out != nullptr, WP2_STATUS_NULL_PARAMETER);
  out->clear();

  YUVPlane yuv;
  WP2_CHECK_STATUS(
      ToYCbCr(ref, 8, params.use_yuv444 ? nullptr : &SamplingTaps::kDownSharp,
              &yuv));

  const uint32_t w = ref.width;
  const uint32_t h = ref.height;
  constexpr int kDepth = 8;
  const avifPixelFormat format =
    params.use_yuv444 ? AVIF_PIXEL_FORMAT_YUV444 :  AVIF_PIXEL_FORMAT_YUV420;
  avifImage* const image = avifImageCreate(w, h, kDepth, format);
  const auto C1 = Cleaner<avifImage>(image, avifImageDestroy);

  const bool has_alpha = ref.HasTransparency();
  if (has_alpha) {
    avifImageAllocatePlanes(image, AVIF_PLANES_YUV | AVIF_PLANES_A);
  } else {
    avifImageAllocatePlanes(image, AVIF_PLANES_YUV);
  }

  // import into yuv
  yuv.Y.To(image->yuvPlanes[AVIF_CHAN_Y], image->yuvRowBytes[AVIF_CHAN_Y]);
  yuv.U.To(image->yuvPlanes[AVIF_CHAN_U], image->yuvRowBytes[AVIF_CHAN_U], 128);
  yuv.V.To(image->yuvPlanes[AVIF_CHAN_V], image->yuvRowBytes[AVIF_CHAN_V], 128);
  if (has_alpha) {
    WP2_CHECK_STATUS(yuv.A.Resize(w, h));
    for (uint32_t y = 0; y < ref.height; ++y) {
      auto* const row_yuv = yuv.A.Row(y);
      const uint8_t* const row = (const uint8_t*)ref.GetRow(y);
      for (uint32_t x = 0; x < ref.width; ++x) {
        const uint8_t* argb = &row[x * WP2FormatBpp(ref.format)];
        row_yuv[x] = argb[0];
      }
    }
    yuv.A.To(image->alphaPlane, image->alphaRowBytes);
  }

  // Encode.
  avifRWData avif_enc = AVIF_DATA_EMPTY;
  const auto C3 = Cleaner<avifRWData>(&avif_enc, avifRWDataFree);
  avifEncoder* const encoder = avifEncoderCreate();
  const auto C4 = Cleaner<avifEncoder>(encoder, avifEncoderDestroy);
  encoder->maxThreads = 1;
  // Scale to go from AVIF_QUANTIZER_BEST_QUALITY to
  // AVIF_QUANTIZER_WORST_QUALITY.
  encoder->minQuantizer = ToAvifQuantizer(params.quality);
  encoder->maxQuantizer = ToAvifQuantizer(params.quality);
  encoder->minQuantizerAlpha = ToAvifQuantizer(params.alpha_quality);
  encoder->maxQuantizerAlpha = ToAvifQuantizer(params.alpha_quality);
  encoder->speed = (int)AVIF_SPEED_FASTEST +
                   ((int)AVIF_SPEED_SLOWEST - (int)AVIF_SPEED_FASTEST) *
                       (int)params.speed / 9;
  double start_time = GetStopwatchTime();
  const avifResult encodeResult = avifEncoderWrite(encoder, image, &avif_enc);
  if (timing != nullptr) timing[0] = GetStopwatchTime() - start_time;
  if (encodeResult == AVIF_RESULT_OK) {
    // output contains a valid .avif file's contents
    *out = std::string(reinterpret_cast<char*>(avif_enc.data), avif_enc.size);
  } else {
    fprintf(stderr, "ERROR: Failed to encode: %s\n",
            avifResultToString(encodeResult));
  }

  // Decode back.
  avifImage* const avif_dec = avifImageCreateEmpty();
  avifDecoder* const decoder = avifDecoderCreate();
  start_time = GetStopwatchTime();
  const avifResult decodeResult =
      avifDecoderRead(decoder, avif_dec, (avifROData*)&avif_enc);
  if (timing != nullptr) timing[1] = GetStopwatchTime() - start_time;
  avifDecoderDestroy(decoder);
  if (decodeResult == AVIF_RESULT_OK) {
    // re-use 'yuv' plane
    yuv.Y.From(avif_dec->yuvPlanes[AVIF_CHAN_Y],
               avif_dec->yuvRowBytes[AVIF_CHAN_Y]);
    yuv.U.From(avif_dec->yuvPlanes[AVIF_CHAN_U],
               avif_dec->yuvRowBytes[AVIF_CHAN_U], -128);
    yuv.V.From(avif_dec->yuvPlanes[AVIF_CHAN_V],
               avif_dec->yuvRowBytes[AVIF_CHAN_V], -128);
    if (has_alpha) {
      yuv.A.From(avif_dec->alphaPlane, avif_dec->alphaRowBytes);
    }
    WP2_CHECK_STATUS(ToRGBA(yuv, &SamplingTaps::kUpSmooth, decoded));
  } else {
    fprintf(stderr, "ERROR: Failed to encode: %s\n",
            avifResultToString(decodeResult));
  }

  return WP2_STATUS_OK;
}

#else

WP2Status CompressAVIF(const ArgbBuffer& ref, const ParamsAV1& params,
                       ArgbBuffer* const decoded, std::string* const out,
                       double timing[2]) {
  (void)ref;
  (void)params;
  (void)decoded;
  (void)out;
  (void)timing;
  fprintf(stderr, "CompressAVIF: AVIF support not enabled at compile time.\n");
  return WP2_STATUS_UNSUPPORTED_FEATURE;
}

#endif  // WP2_HAVE_AVIF

}  // namespace WP2
