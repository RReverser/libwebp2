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
// WP2 visualization tool
// usage: vwp2 input.png [-q quality][...]
//
// Author: Skal (pascal.massimino@gmail.com)

#if defined(__unix__) || defined(__CYGWIN__)
#define _POSIX_C_SOURCE 200112L  // for setenv
#endif

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <numeric>
#include <string>
#include <vector>

#ifdef HAVE_CONFIG_H
#include "src/wp2/config.h"
#endif

#include "examples/example_utils.h"
#include "imageio/image_dec.h"
#include "imageio/image_enc.h"
#include "imageio/imageio_util.h"
#include "src/common/color_precision.h"
#include "src/common/constants.h"
#include "src/common/lossy/block_size.h"
#include "src/dec/wp2_dec_i.h"
#include "src/dsp/dsp.h"
#include "src/dsp/math.h"
#include "src/enc/wp2_enc_i.h"
#include "src/utils/plane.h"
#include "src/utils/random.h"
#include "src/utils/utils.h"
#include "src/wp2/decode.h"
#include "src/wp2/encode.h"

#if defined(WP2_HAVE_AOM)
#include "extras/aom_utils.h"
#endif

#if defined(WP2_HAVE_OPENGL)
#define GL_SILENCE_DEPRECATION  // to avoid deprecation warnings on MacOS
#if defined(HAVE_GLUT_GLUT_H)
#include <OpenGL/gl.h>
#include <GLUT/glut.h>
#else
#include <GL/gl.h>
#include <GL/glut.h>
#endif

#if defined(WP2_HAVE_GLUT) && defined(FREEGLUT)
#include <GL/freeglut.h>
#endif

#if defined(_MSC_VER) && _MSC_VER < 1900
#define snprintf _snprintf
#endif

namespace {

//------------------------------------------------------------------------------
// Visual Debug: set DecoderInfo::visual_debug to get extra visual data into
// DecoderInfo::debug_output during decoding.

struct VDToken {
  VDToken(const char* const token_name,  // NOLINT - not explicit
          std::vector<VDToken>&& next_tokens = std::vector<VDToken>())
      : name(token_name),
        next(std::move(next_tokens)),
        longest_next_name_size(
            next.empty() ? 0
                         : (int32_t)std::max_element(
                               next.begin(), next.end(),
                               [](const VDToken& lhs, const VDToken& rhs) {
                                 return lhs.name.size() < rhs.name.size();
                               })->name.size()) {}

  const std::string name;
  const std::vector<VDToken> next;       // Direct children.
  const int32_t longest_next_name_size;  // Longest name size amongst 'next'.
};

std::vector<VDToken> GetVDTokensCfl() {
  return {"best-prediction", "best-residuals", "best-slope", "best-intercept"};
}

std::vector<VDToken> GetVDTokensChannelCommon() {
  return {{"compressed"},
          {"prediction", {{"raw"}, {"modes", {"long", "short"}}}},
          {"residuals"},
          {"encoder", {"original-residuals", "prediction-scores"}}};
}

std::vector<VDToken> GetVDTokensAlpha() {
  std::vector<VDToken> r = GetVDTokensChannelCommon();
  r.emplace_back<VDToken>({"is_lossy"});
  r.emplace_back<VDToken>({"lossy"});
  r.emplace_back<VDToken>({"lossless"});
  return r;
}

std::vector<VDToken> GetVDTokensUV() {
  std::vector<VDToken> r = GetVDTokensChannelCommon();
  r.emplace_back("coeff-method");
  r.emplace_back<VDToken>(
      {"chroma-from-luma", {"prediction", "slope", "intercept"}});
  return r;
}

std::vector<VDToken> GetVDTokensY() {
  std::vector<VDToken> r = GetVDTokensChannelCommon();
  r.emplace_back("coeff-method");
  return r;
}

// Tree of possible 'visual_debug' paths.
const VDToken kVisualDebug  // NOLINT - non-trivially destructible global var
    {"",
     {{"decompressed"},  // No vdebug.
      {"blocks",
       {{"partition", {"split-tf", "blocks-only"}},
        "segment-ids",
        {"quantization", {"y", "u", "v", "a"}},
        "is420",
        {"encoder", {"is420-scores", "lambda-mult"}}}},
      {"transform", {"xy", "x", "y"}},
      {"y", GetVDTokensY()},
      {"u", GetVDTokensUV()},
      {"v", GetVDTokensUV()},
      {"a", GetVDTokensAlpha()},
      {"filters",
       {{"filter-block-map", {{"bpp"}, {"res", {"y", "u", "v", "a"}}}},
        {"deblocking-filter",
         {{"diff", {"yuv", "a"}},
          {"strength",
           {{"horizontal", {"y", "u", "v", "a"}},
            {"vertical", {"y", "u", "v", "a"}}}}}},
        {"directional-filter", {"diff", "strength", "direction", "variance"}},
        {"restoration-filter", {"diff", "strength"}},
        {"intertile-filter",
         {"diff",
          {"strength",
           {{"horizontal", {"y", "u", "v", "a"}},
            {"vertical", {"y", "u", "v", "a"}}}}}},
        {"grain-filter", {"diff"}},
        {"alpha-filter", {"diff"}}}},
      {"encoder",
       {{"partition",
         {"method",
          {"pass",
           {"all", "0",  "1",  "2",  "3",  "4",  "5",  "6",  "7",  "8", "9",
            "10",  "11", "12", "13", "14", "15", "16", "17", "18", "19"}},
          "order",
          "score",
          {"multi",
           {{"raw-quant-dct", {"32x32", "16x16", "8x8", "4x4"}},
            {"luma-alpha-gradient", {"32x32", "16x16", "8x8", "4x4"}},
            {"narrow-std-dev", {"32x32", "16x16", "8x8", "4x4"}},
            {"quant-dct", {"32x32", "16x16", "8x8", "4x4"}},
            {"direction", {"32x32", "16x16", "8x8", "4x4"}},
            {"analysis", {"spread-5x5", "direction-4x4"}}}}}},
        {"segmentation", {"variance", "score", "id", "blocks"}},
        {"error-diffusion",
         {{"u", {"propagated-error", "new-error"}},
          {"v", {"propagated-error", "new-error"}}}},
        {"chroma-from-luma",
         {{"u", GetVDTokensCfl()},
          {"v", GetVDTokensCfl()},
          {"a", GetVDTokensCfl()}}}}},
      {"compressed", {"a", "r", "g", "b"}},
      // Requires the reference image to be copied to DecoderInfo::debug_output.
      {"original",
       {"diff",
        "a",
        "r",
        "g",
        "b",
        "y",
        "u",
        "v",
        {"histogram", {"y", "u", "v"}}}},
      {"lossless",
       {{"symbols", "clusters", {"transformed", {"regular", "force-opaque"}}}}},
      // Requires WP2_BITTRACE to be defined during the compilation.
      {"bits-per-pixel", {{"overall"}, {"coeffs", {"y", "u", "v"}}}},
      {"error-map",
       {"PSNR",
        "SSIM",
        "MSSSIM",
        "LSIM",
        {"diff-with-alt", {"PSNR", "SSIM", "MSSSIM", "LSIM"}}}}}};

// Verifies that 'str' starts with 'token' and that the ending or the next
// character (if any) is the separator '/'.
bool StrStartsWithToken(const std::string& str, const std::string& token) {
  return (std::strncmp(str.c_str(), token.c_str(), token.size()) == 0 &&
          (str.size() <= token.size() || str[token.size()] == '/' ||
           (!token.empty() && token.back() == '/')));
}

// Recursively retrieves the paths to all leaves.
void GetAllLeavesPaths(const VDToken& node, const std::string& current_path,
                       std::vector<std::string>* const paths) {
  for (const VDToken& next : node.next) {
    if (next.next.empty()) {
      paths->push_back(current_path + next.name);
    } else {
      GetAllLeavesPaths(next, current_path + next.name + '/', paths);
    }
  }
}

// Returns the path to the first leaf matching 'mask' or to an adjacent one
// depending on 'leaf_offset'.
std::string GetLeaf(const std::string& mask, int32_t leaf_offset = 0) {
  // Printing all paths to an array then parsing it is not great but simple.
  std::vector<std::string> leaves_paths;
  GetAllLeavesPaths(kVisualDebug, "", &leaves_paths);

  for (uint32_t i = 0; i < leaves_paths.size(); ++i) {
    if (StrStartsWithToken(leaves_paths[i], mask)) {
      return leaves_paths[((int32_t)(i + leaves_paths.size()) + leaf_offset) %
                          leaves_paths.size()];
    }
  }
  return "";
}

//------------------------------------------------------------------------------

constexpr uint32_t kNumSparks = 1000;

// Unfortunate global variables. Gathered into a struct for comfort.
struct Params {
  Params() : current_file(~0u) {
    enc_config.thread_level = 64;
    enc_config.info = &einfo;
    dec_config.thread_level = 64;
    dec_config.info = &dinfo;
#if defined(WP2_BITTRACE)
    einfo.store_blocks = true;  // Preemptively to avoid (slow) reencoding.
#endif
  }
  WP2_NO_DISCARD bool SetCurrentFile(size_t file_number);
  WP2_NO_DISCARD bool SetBitstream(const char* const file_name);
  WP2_NO_DISCARD bool SetAltFile(const char* const file_name);
  WP2_NO_DISCARD bool SetAltImage();
  WP2_NO_DISCARD bool DumpCurrentCanvas(const char* const file_path);
  void ShiftAlt();   // move alt1 to alt2
  WP2_NO_DISCARD bool EncodeImage();
  WP2_NO_DISCARD bool DecodeOutput();
  WP2_NO_DISCARD bool EncodeAndDecode();
  WP2_NO_DISCARD bool DecodeAndSwap(WP2::ArgbBuffer* const buffer);
  WP2_NO_DISCARD bool GetOriginal();
  WP2_NO_DISCARD bool GetCompressed();
  WP2_NO_DISCARD bool ComputeYUVDistortion();

  const WP2::ArgbBuffer& GetBuffer() const;
  WP2_NO_DISCARD bool GetCurrentCanvas(WP2::ArgbBuffer* const buffer);

  WP2_NO_DISCARD bool CompressWebP(bool match_size);
  float webp_distortion = 0.f;
  uint32_t webp_size = 0;
  WP2_NO_DISCARD bool CompressAV1(bool copy_partition, float av1_quality,
                                  size_t* const av1_file_size = nullptr);
  WP2_NO_DISCARD bool CompressAV1ToMatch(bool copy_partition,
                                         size_t target_file_size);

  WP2Status status = WP2_STATUS_OK;

  enum Show {
    kDebug,
    kOriginal,
    kCompressed,
    kPreview,
    kPreviewColor,
    kAlt1, kAlt2,
    kWebP,
    kInfo,
    kHelp,
    kMenu
  };
  Show show = kDebug;
  int mouse_x = 20, mouse_y = 20;  // In pixels, starting from top left.
  int mouse_x_down = 20, mouse_y_down = 20;  // Position at last GLUT_DOWN.
  bool moved_since_last_down = false;
  std::string message;
  int message_end_time;
  bool show_interface = true;

  float sparks_x[kNumSparks] = {0};
  float sparks_y[kNumSparks] = {0};
  uint32_t sparks_index = 0;
  bool display_sparks = false;

  enum {
    kDisplayNone,
    kDisplayHeader,
    kDisplayYCoeffs,
    kDisplayUCoeffs,
    kDisplayVCoeffs,
    kDisplayACoeffs,
    kDisplayPredModes,
    kDisplayGrid,
    kDisplayNum
  } display_block = kDisplayNone;
  void DisplayBlockInfo(std::vector<std::string>* const msg);

  const char* partition_file_path = "/tmp/partition.txt";
  const char* dump_wp2_path = "/tmp/dump.wp2";
  const char* dump_png_path = "/tmp/dump.png";
  // If not empty, is the current block being drawn with the mouse.
  WP2::Rectangle forcing_block;
  // Same as EncoderInfo::force_partition but not yet encoded.
  std::vector<WP2::Rectangle> force_partition;
  std::vector<WP2::EncoderInfo::ForcedSegment> force_segment;
  std::vector<WP2::EncoderInfo::ForcedPredictor> force_predictor;
  void ClearForcedElements();
  bool IsDebugViewVisible() const;
  bool IsPartitionVisible() const;
  bool IsSegmentVisible() const;
  bool IsPredictorVisible() const;
  void DisplayPartition() const;
  void DisplayForcedSegments() const;
  void DisplayForcedPredictors() const;
  uint32_t getLineThickness() const;
  bool ReadPartition(bool read_only = true);
  // Checks if either the encoding or decoding visual debug config match the
  // given token.
  bool VDMatchEncDec(const char* token) const;

  bool IsHamburgerMenuVisible() const;

  void ApplyZoomLevel();

  enum Background {
    kCheckerboard,
    kWhite,
    kBlack,
    kPink,
    kBackgroundNum,
  } background = kCheckerboard;

  // Debug environment value. Can be retrieved with 'atoi(getenv("WP2DBG"))'.
  uint32_t wp2dbg_value = 0;

  size_t current_file;
  std::vector<std::string> files;
  std::string input;   // input file
  std::string output;  // currently encoded file

  float riskiness = 0.f;

  uint32_t enc_duration_ms;  // encoding time in ms
  uint32_t dec_duration_ms;  // decoding time in ms

  float quality = 75.;
  float alpha_quality = 100.;
  WP2::EncoderConfig enc_config;  // encoding parameters
  WP2::DecoderConfig dec_config;
  std::string menu_selection = "decompressed";  // hovered menu selection
  std::string visual_debug = "decompressed";    // displayed debug output

  bool view_only = false;
  WP2::ArgbBuffer in;      // original samples in RGB
  WP2::ArgbBuffer in_yuv;  // original alpha, Y, U or V samples
  // image shown on spacebar (can be a view on 'in' or 'in_yuv')
  WP2::ArgbBuffer original;
  WP2::ArgbBuffer out;  // recompressed samples
  WP2::ArgbBuffer out_yuv;  // recompressed samples
  WP2::ArgbBuffer preview;
  WP2::Argb32b preview_color;
  WP2::ArgbBuffer webp;  // result of WebP compression
  // alternate comparison pictures & messages
  WP2::ArgbBuffer alt1, alt2;
  std::string alt1_msg, alt2_msg;

  WP2::EncoderInfo einfo;
  WP2::DecoderInfo dinfo;
  std::string original_selection_info;
  WP2::MetricType distortion_metric = WP2::PSNR;
  float distortion[5] = {0.f}, yuv_distortion[3] = {0.f};
  int bit_trace = 0;
  uint32_t bit_trace_level = 0;
  uint32_t visual_bit_trace_level = 3;

  uint32_t width, height;
  uint32_t viewport_width, viewport_height;
  int zoom_level = 0;  // Multiply width, height by 2^zoom_level.

  bool crop = false;
  WP2::Rectangle crop_area;
  bool keep_metadata = false;
} kParams;  // NOLINT (non trivially destructible global var)

void SetMenuSelection(const std::string& selection) {
  kParams.menu_selection = selection;
  assert(!kParams.menu_selection.empty());
  kParams.visual_debug = kParams.menu_selection;
  kParams.dinfo.visual_debug = kParams.visual_debug.c_str();
  kParams.einfo.visual_debug = kParams.visual_debug.c_str();

  if (!WP2::VDMatch(kParams.enc_config, "encoder")) {
    kParams.einfo.visual_debug = nullptr;
  }
  if (WP2::VDMatch(kParams.dec_config, "decompressed") ||
      WP2::VDMatch(kParams.dec_config, "encoder")) {
    // "decompressed" is in fact the regular output.
    // "encoder" is reserved to EncoderConfig.
    kParams.dinfo.visual_debug = nullptr;
  }
}

//------------------------------------------------------------------------------
// Drawing tools

// Displays 'text' starting at pixel (x, y) in OpenGL coordinates [-1:1].
void PrintString(float x, float y, const std::string& text,
                 bool small = false) {
  glRasterPos2f(x, y);
  void* const font = small ? GLUT_BITMAP_8_BY_13 : GLUT_BITMAP_9_BY_15;
  for (char c : text) glutBitmapCharacter(font, c);
}

// Converts from pixel sizes to relative sizes as used by OpenGL.
float ToRelativeX(int32_t pixels) {
  return pixels * 2.f / kParams.viewport_width;
}
float ToRelativeY(int32_t pixels) {
  return pixels * 2.f / kParams.viewport_height;
}
// Converts from pixel positions to absolute positions as used by OpenGL.
float ToAbsoluteX(int32_t pixels) { return ToRelativeX(pixels) - 1; }
float ToAbsoluteY(int32_t pixels) { return 1 - ToRelativeY(pixels); }

// Converts from viewport pixel position to image pixel coordinates.
uint32_t ToImageX(int viewport_x) {
  return WP2::DivRound(
      (uint32_t)WP2::Clamp<int>(viewport_x, 0, kParams.viewport_width - 1) *
          (kParams.width - 1),
      kParams.viewport_width - 1);
}
uint32_t ToImageY(int viewport_y) {
  return WP2::DivRound(
      (uint32_t)WP2::Clamp<int>(viewport_y, 0, kParams.viewport_height - 1) *
          (kParams.height - 1),
      kParams.viewport_height - 1);
}

// Floors to block-aligned coordinate.
uint32_t AlignWithBlocks(uint32_t coord) {
  return coord - coord % WP2::kMinBlockSizePix;
}

void DisplaySparks() {
  if (!kParams.display_sparks) return;

  static WP2::UniformIntDistribution random(/*seed=*/0u);

  const float mouse_x =
      ToRelativeX(kParams.mouse_x + random.Get<int32_t>(-5, 5)) - 1.f;
  const float mouse_y = -(ToRelativeY(kParams.mouse_y) - 1.f);
  kParams.sparks_x[kParams.sparks_index] = mouse_x;
  kParams.sparks_y[kParams.sparks_index] = mouse_y;
  kParams.sparks_index = (kParams.sparks_index + 1) % kNumSparks;

  glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
  for (uint32_t i = 0; i < kNumSparks; ++i) {
    if (kParams.sparks_x[i] != 0 && kParams.sparks_y[i] != 0) {
      kParams.sparks_y[i] -= ToRelativeY(2);

      glColor4f(random.Get<int32_t>(100, 255) / 255.f,
                random.Get<int32_t>(100, 255) / 255.f,
                random.Get<int32_t>(100, 255) / 255.f, 1.f);
      glRectf(kParams.sparks_x[i] - ToRelativeX(1),
              kParams.sparks_y[i] - ToRelativeY(1),
              kParams.sparks_x[i] + ToRelativeX(1),
              kParams.sparks_y[i] + ToRelativeY(1));
      if (kParams.sparks_y[i] < -1) {
        kParams.sparks_x[i] = kParams.sparks_y[i] = 0;
      }
    }
  }
}

//------------------------------------------------------------------------------
// Hamburger menu

constexpr int32_t kMenuCharWidth = 9;
constexpr int32_t kMenuItemHeight = 20;
constexpr int32_t kMenuBorderWidth = 2;

// Renders a submenu at (x, y) with children of 'prefix'. Highlights the node
// matching the 'current_selection' or sets it to the hovered entry.
// Returns true if the mouse is inside the submenu.
bool DisplaySubMenu(int32_t x, int32_t y, const VDToken& prefix,
                    std::string* const current_selection) {
  const int32_t longest_name_size = prefix.longest_next_name_size;
  const int32_t width = (longest_name_size + 3) * kMenuCharWidth +
                                  2 * kMenuBorderWidth;
  const int32_t height = prefix.next.size() * kMenuItemHeight;
  const int32_t left = x, right = x + width - 1;  // Inclusive.

  const int32_t mx = kParams.mouse_x, my = kParams.mouse_y;
  const bool mouse_is_nearby = (mx >= left && mx <= right + 50 &&
                                my + 150 >= y && my <= y + height + 150);
  const bool mouse_in_column = (mouse_is_nearby && mx <= right);
  const bool mouse_in_submenu =
      (mouse_in_column && my >= y && my <= y + height - 1);
  bool was_hovered = false;  // True if an element is hovered by the mouse.

  // Background.
  const float kMenuAlpha = .8;
  glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
  glColor4f(0.f, 0.f, 0.f, kMenuAlpha);
  glRectf(ToAbsoluteX(x), ToAbsoluteY(y), ToAbsoluteX(x + width),
          ToAbsoluteY(y + height));

  // For each branch in the current VisualDebug tree node.
  int32_t top = y;
  for (const VDToken& token : prefix.next) {
    const int32_t bottom = top + kMenuItemHeight - 1;  // Inclusive.

    // Set 'current_selection' if hovered by the mouse.
    const bool is_hovered = (mouse_in_submenu && my >= top && my <= bottom);
    if (is_hovered) {
      assert(!was_hovered);  // Check that only one item is hovered.
      *current_selection = token.name;
      if (!token.next.empty()) *current_selection += "/";
      was_hovered = true;
    }

    // Do not set as selected if another entry in submenu is hovered.
    bool is_selected =
        (is_hovered || (!mouse_in_column &&
                        StrStartsWithToken(*current_selection, token.name)));

    std::string displayed_name = token.name;
    if (!token.next.empty()) {
      displayed_name.append(longest_name_size - token.name.size(), ' ');
      displayed_name += " /";
    }
    if (is_hovered) {
      glColor4f(1.0f, 1.0f, 1.0f, kMenuAlpha);
    } else if (is_selected) {
      glColor4f(0.2f, 0.8f, 0.2f, kMenuAlpha);
    } else {
      glColor4f(0.7f, 0.7f, 0.7f, kMenuAlpha);
    }
    PrintString(ToAbsoluteX(left + kMenuBorderWidth), ToAbsoluteY(bottom - 5),
                displayed_name);

    // Display submenu if any and selected.
    if (is_selected && !token.next.empty()) {
      std::string suffix =
          current_selection->substr(/*pos*/ token.name.size() + 1);
      *current_selection =
          current_selection->substr(/*pos*/ 0, /*n*/ token.name.size() + 1);
      was_hovered |= DisplaySubMenu(right + 1, top, token, &suffix);
      *current_selection += suffix;
    }

    top += kMenuItemHeight;
  }

  if (mouse_in_submenu) assert(was_hovered);
  if (mouse_is_nearby && !was_hovered) *current_selection = "";
  return was_hovered || mouse_is_nearby;
}

// Renders the menu and selected or hovered submenus.
void DisplayMenu() {
  const bool is_hamburger_hovered =
      (kParams.mouse_x >= 0 && kParams.mouse_x < kMenuItemHeight &&
       kParams.mouse_y >= 0 && kParams.mouse_y < kMenuItemHeight);

  // Show the hamburger menu.
  if (kParams.IsHamburgerMenuVisible()) {
    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
    glColor4f(0.f, 0.f, 0.f, .8f);
    glRectf(-1, 1, -1 + ToRelativeX(kMenuItemHeight),
            1 - ToRelativeY(kMenuItemHeight));
    glColor4f(1.f, 1.f, 1.f, .8f);

    for (uint32_t i = 0; i < 3; ++i) {
      glRectf(-1 + ToRelativeX(2), 1 - ToRelativeY(i * 5 + 3),
              -1 + ToRelativeX(kMenuItemHeight - 2),
              1 - ToRelativeY(i * 5 + 5));
    }

    // Show menu items if hamburger menu is hovered.
    if (is_hamburger_hovered) kParams.show = Params::kMenu;
  } else if (kParams.show == Params::kMenu) {
    kParams.show = Params::kDebug;
    kParams.menu_selection = kParams.visual_debug;
  }

  // Show menu items.
  if (kParams.show == Params::kMenu) {
    const bool menu_hovered = DisplaySubMenu(0, kMenuItemHeight, kVisualDebug,
                                             &kParams.menu_selection);
    if (!menu_hovered && !is_hamburger_hovered) {
      kParams.show = Params::kDebug;
      kParams.menu_selection = kParams.visual_debug;
    }
  }
}

//------------------------------------------------------------------------------
// Messages

constexpr uint32_t kFontWidthPix = 9;

// Split 'str' by the line break delimiter into 'msg' entries.
void SplitIntoMessages(const std::string& str,
                       std::vector<std::string>* const msg) {
  std::string::size_type from = 0, to;
  while (from < str.size()) {
    if ((to = str.find('\n', from)) == std::string::npos) to = str.size();
    msg->emplace_back(str.substr(from, to - from));
    from = to + 1;
  }
  // Remove any empty line at the back.
  while (!msg->empty() && msg->back().empty()) msg->pop_back();
}

void PrintMessages(const std::vector<std::string>& msg, bool print_lower,
                   const float text_color[4], const float bg_color[4],
                   bool small = false, bool outline = false) {
  if (msg.empty()) return;

  const float line_height = ToRelativeY(19);
  const float line_start =
      1.f - ToRelativeY(print_lower ? 25 : 5) - line_height;
  const float left_start = -1.f + ToRelativeX(10);

  // Background unless fully transparent
  if (bg_color[3] > 0) {
    const float y_top = line_start + line_height;
    const float x_left = left_start - ToRelativeX(5);
    float HY = msg.size() * line_height;
    float HX = 0;
    for (const std::string& str : msg) HX = std::max(HX, (float)str.size());
    HX *= ToRelativeX(kFontWidthPix);
    HX += ToRelativeX(10);
    HY += ToRelativeY(10);
    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
    glColor4f(bg_color[0], bg_color[1], bg_color[2], bg_color[3]);
    glRectf(x_left, y_top, x_left + HX, y_top - HY);
  }

  if (outline) {
    glColor4f(bg_color[0], bg_color[1], bg_color[2], text_color[3]);
    for (int dy = -1; dy <= 1; dy += 2) {
      for (int dx = -1; dx <= 1; dx += 2) {
        for (size_t i = 0; i < msg.size(); ++i) {
          const float position = line_start - i * line_height;
          PrintString(left_start + ToRelativeX(dx), position + ToRelativeY(dy),
                      msg[i], small);
        }
      }
    }
  }

  // Text
  glColor4f(text_color[0], text_color[1], text_color[2], text_color[3]);
  for (size_t i = 0; i < msg.size(); ++i) {
    const float position = line_start - i * line_height;
    PrintString(left_start, position, msg[i], small);
  }
}

#if defined(WP2_BITTRACE)
constexpr uint32_t kBitTraceMaxLevels = 10;
constexpr float kInfoAlpha = 0.8f;

// Displays a sort of bar graph at the bottom of the screen showing the relative
// size of different syntax elements.
// 'traces' contains the size in bits taken by each (group of) symbols.
void DisplayBitTraces(const std::map<std::string, float>& traces,
                      double total_size) {
  const float line_height = ToRelativeY(35);

  const float mouse_x = ToRelativeX(kParams.mouse_x) - 1;
  const float mouse_y = -(ToRelativeY(kParams.mouse_y) - 1);
  glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);

  // First pass draws the bars, second pass draws the labels, third pass draws
  // hovered label.
  for (uint32_t pass = 0; pass < 3; ++pass) {
    std::map<std::string, float> left_pos;
    uint32_t index[kBitTraceMaxLevels] = {0};

    for (const auto& v : traces) {
      const float relative_size = v.second / total_size;
      const float width = relative_size * 2;
      const int level = std::count(v.first.begin(), v.first.end(), '/');
      std::string prefix = "";
      if (level > 0) {
        prefix = v.first.substr(0, v.first.find_last_of("/"));
      }

      const float bottom = -1 + line_height * level;
      const float top = bottom + line_height;
      const float left = left_pos[prefix] - 1;
      const float right = left + width;

      if (pass == 0) {
        // First pass : draw the bars.
        const GLfloat color = index[level] % 2;
        glColor4f(color, color, color, kInfoAlpha);
        glRectf(left, top, right, bottom);
        // Add a black line at the top.
        glColor4f(0.f, 0.f, 0.f, kInfoAlpha);
        glRectf(left, top, right, top - ToRelativeY(1));
      } else if (pass == 1) {
        // Second pass: draw the labels.
        glColor4f(1.f, 0.f, 0.f, kInfoAlpha);
        const float top_pos = bottom + ToRelativeY(10);

        std::string title = v.first;
        if (level > 0) {
          title = title.substr(title.find_last_of("/") + 1);
        }
        const uint32_t num_chars = std::max(
            0.f,
            std::floor((width - ToRelativeX(4)) / ToRelativeX(kFontWidthPix)));
        if (num_chars > 1) {
          if (num_chars < title.length()) {
            title = title.substr(0, num_chars - 1);
            title += ".";  // To show that it's truncated.
          }
          PrintString(left_pos[prefix] - 1 + ToRelativeX(2), top_pos, title);
        }
      } else if (pass == 2 && mouse_x > left && mouse_x < right &&
                 mouse_y < top && mouse_y > bottom) {
        // Third pass: draw hover label.
        std::string title = v.first;
        if (level > 0) {
          title = title.substr(title.find_last_of("/") + 1);
        }
        title += SPrintf(" (%.f B, %.2f%%)", title.c_str(), v.second / 8.f,
                         relative_size * 100);
        const float title_width = title.size() * ToRelativeX(kFontWidthPix);
        const float title_x =
            (mouse_x + title_width > 1) ? 1 - title_width : mouse_x;
        // Draw a black outline.
        glColor4f(0.f, 0.f, 0.f, 1.f);
        static constexpr int kBorderWidth = 2;
        for (int dy = -kBorderWidth; dy <= kBorderWidth; ++dy) {
          for (int dx = -kBorderWidth; dx <= kBorderWidth; ++dx) {
            PrintString(title_x + ToRelativeX(dx), mouse_y + ToRelativeY(dy),
                        title);
          }
        }
        glColor4f(1.f, 1.f, 0.f, 1.f);
        PrintString(title_x, mouse_y, title);
      }
      ++index[level];
      left_pos[v.first] = left_pos[prefix];
      left_pos[prefix] += width;
    }
  }
}
#endif

//------------------------------------------------------------------------------

std::vector<std::string> GetHelp() {
  std::vector<std::string> m;
  // Adding a new shortcut? These are still available: j, p, u, y
  m.emplace_back("Keyboard shortcuts:");
  m.emplace_back(
      "  'i' ............... overlay file information ('I' for YUV disto)");
  m.emplace_back("  'v'/'V' ........... cycle through visual debugging");
  m.emplace_back("  'e'/'E' ........... shortcut for error map");
  m.emplace_back(
      "  space ............. show the original uncompressed picture");
  m.emplace_back("  tab ............... show the compressed picture");
  m.emplace_back("  alt+1 / alt+2 ..... show alternate picture(s)");
  m.emplace_back(
      "  'a' ............... saves current canvas as alternate picture");
  m.emplace_back(std::string("  'A' ............... saves current canvas to ") +
                 kParams.dump_png_path);
  m.emplace_back(std::string("  alt+'a' ........... loads ") +
                 kParams.dump_png_path + std::string(" as alternate picture"));
  m.emplace_back(
      "  up/down ........... change the compression factor by +/- 1 units");
  m.emplace_back(
      "  left/right ........ change the compression factor by +/- 10 units");
  m.emplace_back("  ctrl + arrows ..... change alpha compression factor");
  m.emplace_back("  's' ............... change encoding speed parameter");
  m.emplace_back("  't' ............... toggle segment id mode");
  m.emplace_back(
      "  '1' ............... cycle colorspace type "
      "(YCoCg,YCbCr,Custom,YIQ)");
  m.emplace_back("  '2' ............... cycle UV-Mode (Adapt,420,Sharp,444)");
  m.emplace_back("  '3' ............... change tile size");
  m.emplace_back("  '4' ............... change the partition method");
  m.emplace_back("  '5' ............... change the partition set");
  m.emplace_back("  '6' ............... toggle block snapping");
  m.emplace_back("  '7' ............... increase number of segments");
  m.emplace_back("  '8' ............... increase SNS value");
  m.emplace_back("  '9' ............... increase error diffusion strength");
  m.emplace_back("  '0' ............... toggle perceptual tuning");
  m.emplace_back("  'm' ............... toggle use of random matrix");
  m.emplace_back("  'x' ............... toggle the env variable WP2DBG");
  m.emplace_back("  'f'/'F'/alt+'f/F' . toggle filters");
  m.emplace_back("  'g'/'G' ........... decrease/reset grain amplitude");
  m.emplace_back("  'b' ............... display block side-info");
  m.emplace_back("  'k' ............... toggle metadata");
  m.emplace_back("  'o' ............... change background color (for alpha)");
  m.emplace_back("  '\\' ............... show the preview if available");
  m.emplace_back("  '|' ............... show the preview color if available");
  m.emplace_back("  'z'/'Z' ........... zoom in/out");
  m.emplace_back("  'n'/'N' ........... change distortion metric");
#if defined(WP2_HAVE_WEBP)
  m.emplace_back("  'w'/'W' ........... show WebP equivalent in size / disto");
#endif
#if defined(WP2_HAVE_AOM)
  m.emplace_back("  'r' ............... show AV1 equivalent in quality factor");
  m.emplace_back("  alt+'r' ........... show AV1 with same file size (slow)");
#endif
#if defined(WP2_HAVE_AOM_DBG)
  m.emplace_back("  'R' ............... copy AV1 partition (luma transforms)");
#endif
  m.emplace_back(std::string("  'l'/'L' ........... "
                             "load or dump the current bitstream as ") +
                 kParams.dump_wp2_path);
  m.emplace_back("  'h' ............... show this help message");
  m.emplace_back("  'H' ............... toggle interface display");
  m.emplace_back("  'q' / 'Q' / ESC ... quit");

  m.emplace_back("");
  m.emplace_back("Block layout edition:");
  m.emplace_back("  Editor is shown by bringing the \"blocks/partitioning\"");
  m.emplace_back("  menu ('v' key) or the \"display block\" mode ('b' key).");
  m.emplace_back("  Force blocks position and size by drawing rectangles with");
  m.emplace_back("  the right mouse button. Remove one with right or middle ");
  m.emplace_back("  click. Press 'c' to convert rectangles (dashed outline)");
  m.emplace_back("  into actual forced encoded blocks (magenta outline).");
  m.emplace_back("  Press 'd' to dump or 'D' to load the blocks in file at");
  m.emplace_back(std::string("  ") + kParams.partition_file_path);
  m.emplace_back("Segment edition: In segment-ids view, right click to cycle");
  m.emplace_back("  through segments, then press 'c' to apply. ");
  m.emplace_back("  Middle click to remove.");
  m.emplace_back("Predictor edition: In prediction view, right click to cycle");
  m.emplace_back("  through predictors, then press 'c' to apply. Middle click");
  m.emplace_back("  to remove. A red square means predictor was ignored ");
  m.emplace_back("  because context is constant.");
  return m;
}

//------------------------------------------------------------------------------
// Info

const char* const kCSPString[] = {"YCoCg", "YCbCr", "Custom", "YIQ"};
const char* const kUVModeString[] = {"UVAdapt", "UV420", "UV444", "UVAuto"};
const char* const kTileShapeString[] = {"128", "256", "512", "Wide", "Auto"};
const char* const kPartitionSetString[] = {
    "Small squares",  "Small rects", "All rects", "Thick rects",
    "Medium squares", "All squares", "Some rects"};
const char* const kOnOff[] = {"Off", "On"};
const char* const kSegmentModes[] = {"Auto", "Explicit", "Implicit"};

STATIC_ASSERT_ARRAY_SIZE(kCSPString, WP2::kNumCspTypes);
STATIC_ASSERT_ARRAY_SIZE(kUVModeString, WP2::EncoderConfig::NumUVMode);
STATIC_ASSERT_ARRAY_SIZE(kPartitionSetString, WP2::NUM_PARTITION_SETS);

void PrintInfo() {
  DisplaySparks();
  std::vector<std::string> msg;
  kParams.DisplayBlockInfo(&msg);
  kParams.DisplayPartition();
  kParams.DisplayForcedSegments();
  kParams.DisplayForcedPredictors();
  DisplayMenu();

  if (kParams.show == Params::kMenu ||
      (kParams.forcing_block.width > 0 && kParams.forcing_block.height > 0)) {
    // Do not print messages if something else is displayed instead.
    return;
  }

  bool print_lower = kParams.IsHamburgerMenuVisible();
  float text_color[4] = {.0f, .0f, .0f, .9f};
  float bg_color[4] = {.8f, .8f, .9f, .9f};
  bool small = false;

  if (kParams.display_block != Params::kDisplayNone) {
    bg_color[0] = .8f;
    bg_color[1] = .7f;
    bg_color[2] = .7f;
    bg_color[3] = .8f;
    small = (kParams.display_block > Params::kDisplayHeader);
  } else if (kParams.show == Params::kHelp) {
    const std::vector<std::string> help = GetHelp();
    msg.insert(msg.end(), help.begin(), help.end());
  } else if (kParams.show == Params::kOriginal) {
    msg.emplace_back("- O R I G I N A L -");
    if (!kParams.original_selection_info.empty()) {
      SplitIntoMessages(kParams.original_selection_info, &msg);
      print_lower = true;  // Match standard kDebug text formatting.
      small = true;
    }
  } else if (kParams.show == Params::kCompressed) {
    msg.emplace_back("- Compressed -");
  } else if (kParams.show == Params::kPreview) {
    msg.emplace_back("- Preview -");
  } else if (kParams.show == Params::kPreviewColor) {
    msg.emplace_back("- Preview color -");
  } else if (kParams.show == Params::kWebP) {
    msg.emplace_back(SPrintf("- WebP [size=%u %s=%.2f] -", kParams.webp_size,
                             kWP2MetricNames[kParams.distortion_metric],
                             kParams.webp_distortion));
    msg.emplace_back(SPrintf(
        "[ref WP2: size=%u %s=%.2f]", kParams.output.size(),
        kWP2MetricNames[kParams.distortion_metric], kParams.distortion[4]));
  } else if (kParams.show == Params::kInfo) {
    char tmp[80];

    bg_color[0] = .8f;
    bg_color[1] = .9f;
    bg_color[2] = .9f;
    bg_color[3] = .8f;
    small = true;

    msg.emplace_back(kParams.files[kParams.current_file]);

    snprintf(tmp, sizeof(tmp), "Dimension: %d x %d",
             kParams.width, kParams.height);
    msg.emplace_back(tmp);
    if (kParams.viewport_width != kParams.width ||
        kParams.viewport_height != kParams.height) {
      snprintf(tmp, sizeof(tmp), " (displayed as %d x %d)",
               kParams.viewport_width, kParams.viewport_height);
      msg.back() += tmp;
    }

    snprintf(tmp, sizeof(tmp),
             "Quality: %3.1f [speed=%d      sns=%5.1f diffusion=%d",
             kParams.quality, kParams.enc_config.speed, kParams.enc_config.sns,
             kParams.enc_config.error_diffusion);
    msg.emplace_back(tmp);
    snprintf(tmp, sizeof(tmp),
             "              segments=%u/%d(%c) snap=%s]",
             kParams.dinfo.num_segments, kParams.enc_config.segments,
             kSegmentModes[(int)kParams.enc_config.segment_id_mode][0],
             kOnOff[kParams.enc_config.partition_snapping]);
    msg.emplace_back(tmp);

    if (kParams.in.HasTransparency()) {
      snprintf(tmp, sizeof(tmp), "Alpha quality: %.1f (%s=%.2f)",
               kParams.alpha_quality,
               kWP2MetricNames[kParams.distortion_metric],
               kParams.distortion[0]);
      msg.emplace_back(tmp);
    }

    snprintf(tmp, sizeof(tmp), "Size: %zu [%.2f bpp] (enc %u ms + dec %u ms)",
             kParams.output.size(),
             8.f * kParams.output.size() / (kParams.width * kParams.height),
             kParams.enc_duration_ms, kParams.dec_duration_ms);
    msg.emplace_back(tmp);
    snprintf(tmp, sizeof(tmp), " (%.1f%% of original)",
             100. * kParams.output.size() / kParams.input.size());
    msg.back() += tmp;

    if (kParams.in.metadata.IsEmpty()) {
      assert(kParams.out.metadata.IsEmpty());
    } else {
      const WP2::Metadata& mi = kParams.in.metadata;
      const WP2::Metadata& mo = kParams.out.metadata;
      snprintf(tmp, sizeof(tmp),  // 6 spaces to align with "Size: ".
               "      %s %zu B of metadata (ICCP %zu, EXIF %zu, XMP %zu)",
               kParams.keep_metadata ? "included" : "ignored",
               (mi.iccp.size + mi.exif.size + mi.xmp.size), mi.iccp.size,
               mi.exif.size, mi.xmp.size);
      msg.emplace_back(tmp);
      if (!kParams.keep_metadata) {
        assert(kParams.out.metadata.IsEmpty());
      } else if (mi.iccp.size != mo.iccp.size || mi.exif.size != mo.exif.size ||
                 mi.xmp.size != mo.xmp.size) {
        snprintf(tmp, sizeof(tmp), " decoded as ICCP %zu, EXIF %zu, XMP %zu",
                 mo.iccp.size, mo.exif.size, mo.xmp.size);
        msg.back() += tmp;
      }
    }

#if defined(WP2_ENC_DEC_MATCH)
    snprintf(tmp, sizeof(tmp),  // 6 spaces to align with "Size: ".
             "      Warning: size is inflated by WP2_ENC_DEC_MATCH");
    msg.emplace_back(tmp);
#endif

    snprintf(tmp, sizeof(tmp), "Filters: dblk=%s, drct=%s, rstr=%s, alpha=%s",
             kOnOff[kParams.dec_config.enable_deblocking_filter],
             kOnOff[kParams.dec_config.enable_directional_filter],
             kOnOff[kParams.dec_config.enable_restoration_filter],
             kOnOff[kParams.enc_config.enable_alpha_filter]);
    msg.emplace_back(tmp);

    snprintf(tmp, sizeof(tmp), "%s=%.2f  (a %.1f, r %.1f, g %.1f, b %.1f)",
             kWP2MetricNames[kParams.distortion_metric], kParams.distortion[4],
             kParams.distortion[0], kParams.distortion[1],
             kParams.distortion[2], kParams.distortion[3]);
    msg.emplace_back(tmp);
    if (kParams.yuv_distortion[0] > 0.f) {  // Maybe it is unavailable.
      snprintf(tmp, sizeof(tmp), "  (y %.1f, u %.1f, v %.1f)",
               kParams.yuv_distortion[0], kParams.yuv_distortion[1],
               kParams.yuv_distortion[2]);
      msg.back() += tmp;
    }

    snprintf(tmp, sizeof(tmp), "CSP: %s   UV: %s  Partition: %s, %s  Tile: %s",
             kCSPString[(int)kParams.enc_config.csp_type],
             kUVModeString[(int)kParams.enc_config.uv_mode],
             WP2::kPartitionMethodString[kParams.enc_config.partition_method],
             kPartitionSetString[kParams.enc_config.partition_set],
             kTileShapeString[kParams.enc_config.tile_shape]);
    msg.emplace_back(tmp);
    snprintf(tmp, sizeof(tmp), "Perceptual:%s  RndMtx: %s  Grain: %s Amp: %d",
             kOnOff[kParams.enc_config.tune_perceptual],
             kOnOff[kParams.enc_config.use_random_matrix],
             kOnOff[kParams.enc_config.store_grain],
             kParams.dec_config.grain_amplitude);
    msg.emplace_back(tmp);

#if defined(WP2_BITTRACE)
    double total_size = 0;
    for (const auto& p : kParams.dinfo.bit_traces) {
      total_size += p.second.bits;
    }

    std::map<std::string, float> traces;
    for (const auto& p : kParams.dinfo.bit_traces) {
      uint32_t level = 0;
      size_t i = 0;
      while (level < kParams.visual_bit_trace_level && i != std::string::npos) {
        i = p.first.find("/", (i == 0) ? 0 : i + 1);
        traces[p.first.substr(0, i)] += p.second.bits;
        ++level;
      }
    }
    DisplayBitTraces(traces, total_size);
#endif  // WP2_BITTRACE
  } else if (kParams.show == Params::kAlt1) {
    msg.emplace_back(kParams.alt1_msg);
  } else if (kParams.show == Params::kAlt2) {
    msg.emplace_back(kParams.alt2_msg);
  } else {
    bg_color[0] = .9f;
    bg_color[1] = .9f;
    bg_color[2] = .9f;
    bg_color[3] = .8f;
    small = true;

    if (WP2::VDMatch(kParams.enc_config, "")) {
      msg.emplace_back(kParams.einfo.visual_debug);
      SplitIntoMessages(kParams.einfo.selection_info, &msg);
    } else if (WP2::VDMatch(kParams.dec_config, "")) {
      msg.emplace_back(kParams.dinfo.visual_debug);
      SplitIntoMessages(kParams.dinfo.selection_info, &msg);
    }
  }

  if (!kParams.message.empty()) {
    msg.push_back(kParams.message);
    if (glutGet(GLUT_ELAPSED_TIME) > kParams.message_end_time) {
      kParams.message.clear();
    }
  }

  PrintMessages(msg, print_lower, text_color, bg_color, small);
}

// Adds a temporary 'message' on screen.
void SetMessage(const std::string& message, bool even_if_show_info = false,
                int duration_ms = 2000) {
  if (kParams.show != Params::kInfo || even_if_show_info) {
    kParams.message = message;
    kParams.message_end_time = glutGet(GLUT_ELAPSED_TIME) + duration_ms;
  }
}

// Same as above with different default arguments.
void AddInfo(const std::string& message) {
  SetMessage(message, /*even_if_show_info=*/true, /*duration_ms=*/5000);
}

void SetError(const std::string& message) {
  SetMessage(message, /*even_if_show_info=*/true, /*duration_ms=*/60000);

  // Reset outputs so that nothing misleading is displayed.
  const uint32_t width = kParams.in.width, height = kParams.in.height;
  for (WP2::ArgbBuffer* const buffer :
       {&kParams.out, &kParams.dinfo.debug_output}) {
    if (buffer->Resize(width, height) != WP2_STATUS_OK) assert(false);
    buffer->Fill({0, 0, 0, 0});
  }
#if defined(WP2_BITTRACE)
  kParams.einfo.blocks.clear();
  kParams.dinfo.bit_traces.clear();
  kParams.dinfo.blocks.clear();
#endif  // WP2_BITTRACE
  kParams.dinfo.header_size = 0;
  kParams.dinfo.selection_info.clear();
}

//------------------------------------------------------------------------------
// Timer callbacks

void DisplayLoop(int delay_ms) {
  glutTimerFunc((unsigned int)delay_ms, DisplayLoop, delay_ms);
  glutPostRedisplay();
}

//------------------------------------------------------------------------------
// error maps

WP2Status ComputeErrorMap(WP2::MetricType metric,
                          const WP2::ArgbBuffer& original,
                          const WP2::ArgbBuffer& decoded,
                          WP2::ArgbBuffer* const error_map,
                          WP2::DecoderInfo* const info = nullptr) {
  WP2_CHECK_OK(
      decoded.width == original.width && decoded.height == original.height,
      WP2_STATUS_BAD_DIMENSION);
  WP2_CHECK_STATUS(error_map->Resize(original.width, original.height));
  const double norm =
      (metric == WP2::PSNR) ? 5. : (metric == WP2::LSIM) ? 5. : 8.;
  for (uint32_t y = 0; y < error_map->height; ++y) {
    uint8_t* const dst = (uint8_t*)error_map->GetRow(y);
    for (uint32_t x = 0; x < error_map->width; ++x) {
      float d;
      WP2_CHECK_STATUS(decoded.GetDistortion(original, x, y, metric, &d));
      const uint8_t v =
          (uint8_t)std::lround(WP2::Clamp(255. - d * norm, 0., 255.));
      dst[4 * x + 0] = 255;  // TODO(skal): handle alpha properly
      dst[4 * x + 1] = dst[4 * x + 2] = dst[4 * x + 3] = v;

      if (info != nullptr && info->selection.x == x && info->selection.y == y) {
        info->selection_info =
            SPrintf("%s at %u, %u: %f\n", kWP2MetricNames[metric], x, y, d);
      }
    }
  }
  return WP2_STATUS_OK;
}

// Sets 'diff' to be a shadowed version of 'original' but with a greener or
// redder tint depending on 'error_a' or 'error_b' being the lowest.
void ComputeDiffPixel(const uint8_t original[4], uint8_t error_a,
                      uint8_t error_b, uint8_t diff[4]) {
  diff[0] = WP2::kAlphaMax;
  for (uint32_t i : {1, 2, 3}) {
    // Darkened/greyed out version of the original sample.
    diff[i] =
        WP2::RightShiftRound((63u << 1) + ((uint32_t)original[i] << 6), 7);
  }
  if (error_a > error_b) {
    diff[1] =
        (uint8_t)WP2::Clamp(diff[1] + (error_a - error_b) * 10u, 0u, 255u);
  } else if (error_b > error_a) {
    diff[2] =
        (uint8_t)WP2::Clamp(diff[2] + (error_b - error_a) * 10u, 0u, 255u);
  }
}

// Fills 'diff_map' with 'original' pixels that are greener if 'decoded' is
// closer than 'alt', redder otherwise.
WP2Status ComputeDiffMap(WP2::MetricType metric,
                         const WP2::ArgbBuffer& original,
                         const WP2::ArgbBuffer& decoded,
                         const WP2::ArgbBuffer& alt,
                         WP2::ArgbBuffer* const diff_map,
                         WP2::DecoderInfo* const info = nullptr) {
  if (info != nullptr) {
    float disto_dec[5] = {}, disto_alt[5] = {};
    WP2_CHECK_STATUS(decoded.GetDistortion(original, metric, disto_dec));
    WP2_CHECK_STATUS(alt.GetDistortion(original, metric, disto_alt));
    info->selection_info =
        SPrintf("%s of decoded: %5.2f (A %4.1f, R %4.1f, G %4.1f, B %4.1f)\n",
                kWP2MetricNames[metric], disto_dec[4], disto_dec[0],
                disto_dec[1], disto_dec[2], disto_dec[3]);
    info->selection_info +=
        SPrintf("%s of alt:     %5.2f (A %4.1f, R %4.1f, G %4.1f, B %4.1f)\n",
                kWP2MetricNames[metric], disto_alt[4], disto_alt[0],
                disto_alt[1], disto_alt[2], disto_alt[3]);
  }

  WP2::ArgbBuffer error_map_decoded, error_map_alt;
  WP2_CHECK_STATUS(
      ComputeErrorMap(metric, original, decoded, &error_map_decoded));
  WP2_CHECK_STATUS(ComputeErrorMap(metric, original, alt, &error_map_alt));
  WP2_CHECK_STATUS(diff_map->Resize(original.width, original.height));
  const uint32_t bpp = WP2FormatBpp(original.format);
  for (uint32_t y = 0; y < diff_map->height; ++y) {
    for (uint32_t x = 0; x < diff_map->width; ++x) {
      const uint32_t i = x * bpp + 1;
      const uint8_t err_dec = ((const uint8_t*)error_map_decoded.GetRow(y))[i];
      const uint8_t err_alt = ((const uint8_t*)error_map_alt.GetRow(y))[i];
      ComputeDiffPixel(&((const uint8_t*)original.GetRow(y))[x * bpp], err_dec,
                       err_alt, &((uint8_t*)diff_map->GetRow(y))[x * bpp]);

      if (info != nullptr && info->selection.x == x && info->selection.y == y) {
        info->selection_info +=
            SPrintf("Normalized diff at %3u,%3u: decoded %u alt %u\n", x, y,
                    err_dec, err_alt);
      }
    }
  }

  return WP2_STATUS_OK;
}

//------------------------------------------------------------------------------
// compression

bool Params::EncodeImage() {
  enc_config.quality = quality;
  enc_config.alpha_quality = alpha_quality;
  if (!enc_config.IsValid()) {
    SetError("Error: Invalid configuration");
    fprintf(stderr, "Invalid configuration\n");
    return false;
  }

  const double start = GetStopwatchTime();
  if (view_only) {
    output = input;
  } else {
    output.clear();
    einfo.force_partition = force_partition;
    einfo.force_segment = force_segment;
    einfo.force_predictor = force_predictor;

    WP2::StringWriter writer(&output);
    WP2::Metadata metadata;
    using std::swap;
    if (!keep_metadata) swap(metadata, in.metadata);
    const WP2Status enc_status = Encode(in, &writer, enc_config);
    if (!keep_metadata) swap(metadata, in.metadata);
    if (enc_status != WP2_STATUS_OK) {
      SetError("Error: Compression failed");
      fprintf(stderr, "Compression failed: %s\n", WP2GetStatusText(enc_status));
      return false;
    }
  }
  enc_duration_ms = (uint32_t)std::lround(1000. * (GetStopwatchTime() - start));
  return true;
}

bool Params::EncodeAndDecode() {
  if (!EncodeImage()) return false;
  return DecodeOutput();
}

bool Params::DecodeOutput() {
  const double start = GetStopwatchTime();

  if (view_only &&
      WP2::GuessImageFormat((const uint8_t*)output.data(), output.size()) !=
          WP2::FileFormat::WP2) {
    status = out.CopyFrom(in);
  } else {
    if (WP2::VDMatch(dec_config, "original")) {
      // The reference image is needed in 'debug_output' if the 'visual_debug'
      // is set to "original".
      status = dinfo.debug_output.CopyFrom(in);
      if (status != WP2_STATUS_OK) {
        SetError("Error: Out of memory");
        return false;
      }
    }

    status = Decode(output, &out, dec_config);
    if (bit_trace > 0) {
#if defined(WP2_BITTRACE)
      PrintBitTraces(dinfo, output.size(), /*sort_values=*/true,
                     /*use_bytes=*/bit_trace == 2, /*show_histograms=*/false,
                     /*short_version=*/false, bit_trace_level);
#endif
    }
  }
  dec_duration_ms = (uint32_t)std::lround(1000. * (GetStopwatchTime() - start));

  if (status != WP2_STATUS_OK) {
    SetError("Error: Decompression failed");
    fprintf(stderr, "Decompression failed\n");
    return false;
  }

  if (out.GetDistortionBlackOrWhiteBackground(in, distortion_metric,
                                              distortion) != WP2_STATUS_OK) {
    fprintf(stderr, "Error while computing the distortion.\n");
  }
  if (WP2::VDMatch(dec_config, "error-map")) {
    // clang-format off
    const WP2::MetricType metric =
            WP2::VDMatch(dec_config,   "PSNR") ?   WP2::PSNR :
            WP2::VDMatch(dec_config,   "LSIM") ?   WP2::LSIM :
            WP2::VDMatch(dec_config, "MSSSIM") ? WP2::MSSSIM :
                                                  WP2::SSIM;
    // clang-format on
    if (WP2::VDMatch(dec_config, "diff-with-alt")) {
      if (kParams.alt1.IsEmpty()) {
        SetError("alt-1 is empty");
      } else if (ComputeDiffMap(metric, kParams.in, kParams.out, kParams.alt1,
                                &kParams.dec_config.info->debug_output,
                                kParams.dec_config.info) != WP2_STATUS_OK) {
        SetError("Error: Diff map computation failed");
        fprintf(stderr, "Diff map computation failed\n");
        return false;
      }
    } else if (ComputeErrorMap(metric, kParams.in, kParams.out,
                               &kParams.dec_config.info->debug_output,
                               kParams.dec_config.info) != WP2_STATUS_OK) {
      SetError("Error: Error map computation failed");
      fprintf(stderr, "Error map computation failed\n");
      return false;
    }
  }
  if (!force_segment.empty() && !dinfo.explicit_segment_ids) {
    force_segment.clear();
  }

  if (!ComputeYUVDistortion()) return false;
  return true;
}

bool Params::ComputeYUVDistortion() {
  if (!view_only && show == kInfo) {  // Only displayed during 'kInfo'.
    const char* const kOrigVD[] = {"original/y", "original/u", "original/v"};
    const char* const kDecVD[] = {"y/compressed", "u/compressed",
                                  "v/compressed"};
    WP2::ArgbBuffer original_yuv, decompressed_yuv;
    for (WP2::Channel c : {WP2::kYChannel, WP2::kUChannel, WP2::kVChannel}) {
      WP2::DecoderInfo decoder_info;
      dec_config.info = &decoder_info;
      decoder_info.visual_debug = kOrigVD[c];
      if (!DecodeAndSwap(&original_yuv)) return false;
      decoder_info.visual_debug = kDecVD[c];
      if (!DecodeAndSwap(&decompressed_yuv)) return false;
      dec_config.info = &dinfo;
      float disto[5];
      WP2_ASSERT_STATUS(decompressed_yuv.GetDistortion(
          original_yuv, distortion_metric, disto));
      yuv_distortion[c] = disto[1];
    }
  } else {
    yuv_distortion[0] = yuv_distortion[1] = yuv_distortion[2] = 0.f;
  }
  return true;
}

bool Params::VDMatchEncDec(const char* token) const {
  return WP2::VDMatch(dec_config, token) || WP2::VDMatch(enc_config, token);
}

bool Params::DecodeAndSwap(WP2::ArgbBuffer* const buffer) {
  status = buffer->Swap(&dec_config.info->debug_output);  // Avoid realloc
  assert(status == WP2_STATUS_OK);
  status = dec_config.info->debug_output.CopyFrom(in);
  assert(status == WP2_STATUS_OK);
  status = Decode(output, &out, dec_config);
  if (status != WP2_STATUS_OK) {
    fprintf(stderr, "Decompression failed?!\n");
    return false;
  }
  status = buffer->Swap(&dec_config.info->debug_output);
  assert(status == WP2_STATUS_OK);
  return true;
}

bool Params::GetOriginal() {
  if (VDMatch(dec_config, "residuals")) {
    WP2::EncoderInfo encoder_info;
    WP2::Channel c = VDChannel(dec_config);
    const std::string kChannelName[] = {"y", "u", "v", "a"};
    const std::string vdebug = kChannelName[c] + "/encoder/original-residuals";
    encoder_info.visual_debug = vdebug.c_str();
    encoder_info.selection = dinfo.selection;
    enc_config.info = &encoder_info;

    // Avoid realloc
    if (in_yuv.Swap(&enc_config.info->debug_output) != WP2_STATUS_OK) {
      return false;
    }
    if (!EncodeImage()) return false;
    if (in_yuv.Swap(&enc_config.info->debug_output) != WP2_STATUS_OK) {
      return false;
    }

    original_selection_info = encoder_info.selection_info;
    enc_config.info = &einfo;
    if (original.SetView(in_yuv) != WP2_STATUS_OK) return false;
    return true;
  }

  WP2::DecoderInfo decoder_info;
  // clang-format off
  decoder_info.visual_debug = VDMatchEncDec("a") ? "original/a" :
                              VDMatchEncDec("y") ? "original/y" :
                              VDMatchEncDec("u") ? "original/u" :
                              VDMatchEncDec("v") ? "original/v" :
                              VDMatchEncDec("r") ? "original/r" :
                              VDMatchEncDec("g") ? "original/g" :
                              VDMatchEncDec("b") ? "original/b" : nullptr;
  // clang-format on
  decoder_info.selection =
      (VDMatchEncDec("prediction/modes/long") ||
       VDMatchEncDec("prediction/raw") || VDMatchEncDec("compressed"))
          ? dinfo.selection
          : WP2::Rectangle();
  if (decoder_info.visual_debug != nullptr) {
    dec_config.info = &decoder_info;
    if (!DecodeAndSwap(&in_yuv)) return false;
    original_selection_info = decoder_info.selection_info;
    dec_config.info = &dinfo;
    if (original.SetView(in_yuv) != WP2_STATUS_OK) return false;
  } else {
    original_selection_info.clear();
    if (original.SetView(in) != WP2_STATUS_OK) return false;
  }
  return true;
}

bool Params::GetCompressed() {
  WP2::DecoderInfo decoder_info;
  // clang-format off
  decoder_info.visual_debug = VDMatchEncDec("a") ? "a/compressed" :
                              VDMatchEncDec("y") ? "y/compressed" :
                              VDMatchEncDec("u") ? "u/compressed" :
                              VDMatchEncDec("v") ? "v/compressed" : nullptr;
  // clang-format on
  if (decoder_info.visual_debug != nullptr) {
    dec_config.info = &decoder_info;
    if (!DecodeAndSwap(&out_yuv)) return false;
    dec_config.info = &dinfo;
  } else {
    out_yuv.Deallocate();
  }
  return true;
}

bool Params::CompressWebP(bool match_size) {
#if defined(WP2_HAVE_WEBP)
  WebPConfig webp_config;
  if (!WebPConfigPreset(&webp_config, WEBP_PRESET_DEFAULT, quality)) {
    return false;
  }
  if (enc_config.thread_level > 0) ++webp_config.thread_level;
  webp_config.method = (enc_config.speed * 6 + 4) / 9;
  webp_config.segments = std::min(enc_config.segments, 4);
  webp_config.sns_strength = (int)std::lround(enc_config.sns);
  webp_config.pass = 6;
  if (match_size) {
    webp_config.target_size = (int)output.size();
  } else {
    if (distortion_metric != WP2::PSNR) return false;
    webp_config.target_PSNR = distortion[4];
  }
  if (!WebPValidateConfig(&webp_config)) return false;

  WP2::MemoryWriter writer;
  status = WP2::CompressWebP(in, webp_config, &writer);
  if (status != WP2_STATUS_OK) return false;
  webp_size = writer.size_;

  status =
      WP2::ReadImage(writer.mem_, writer.size_, &webp, WP2::FileFormat::WEBP);
  if (status != WP2_STATUS_OK) return false;

  float disto[5];
  status =
      webp.GetDistortionBlackOrWhiteBackground(in, distortion_metric, disto);
  if (status != WP2_STATUS_OK) return false;
  webp_distortion = disto[4];

  show = kWebP;
  return true;
#else
  return false;
#endif  // WP2_HAVE_WEBP
}

bool Params::CompressAV1(bool copy_partition, float av1_quality,
                         size_t* const av1_file_size) {
#if !defined(WP2_HAVE_AOM_DBG)
  if (copy_partition) return false;
#endif
#if defined(WP2_HAVE_AOM)
  WP2::ParamsAV1 cfg;
  cfg.quality = av1_quality;
  // The settings below seem to give the best results (yes, even threads).
  cfg.speed = 9;
  cfg.threads = 1;
  cfg.pass = 2;
  cfg.use_yuv444 =
      (kParams.enc_config.uv_mode != WP2::EncoderConfig::UVMode420);
  cfg.draw_blocks = (kParams.display_block != Params::kDisplayNone);
  cfg.draw_transforms = (kParams.display_block == Params::kDisplayYCoeffs);
  if (copy_partition) force_partition.clear();

  printf("Compressing AV1 at quality %.1f...\n", cfg.quality);
  std::string tmp;
  double timing[2];
  ShiftAlt();
  if (WP2::CompressAV1(
          in, cfg, &alt1, &tmp, timing, /*blocks=*/nullptr,
          /*transforms=*/(copy_partition ? &force_partition : nullptr)) !=
      WP2_STATUS_OK) {
    return false;
  }
  printf("Done compressing AV1 in just %.2f seconds (%zu bytes)!\n", timing[0],
         tmp.size());

  if (av1_file_size != nullptr) *av1_file_size = tmp.size();
  alt1_msg = SPrintf("AV1 - size = %zu", tmp.size());

  if (!cfg.draw_blocks && !cfg.draw_transforms) {
    // Only print distortion if the output image is unaltered by debug overlays.
    float disto[5];
    if (alt1.GetDistortionBlackOrWhiteBackground(in, distortion_metric,
                                                 disto) != WP2_STATUS_OK) {
      return false;
    }
    alt1_msg +=
        SPrintf("  %s = %.1f", kWP2MetricNames[distortion_metric], disto[4]);
  }

  if (copy_partition) {
    WP2ConvertPartition(in.width, in.height, enc_config.partition_set,
                        /*ignore_invalid=*/true, &force_partition);
  }

  show = kAlt1;
  return true;
#else
  (void)copy_partition;
  return false;
#endif  // WP2_HAVE_AOM
}

bool Params::CompressAV1ToMatch(bool copy_partition, size_t target_file_size) {
  // Begin by testing 0, it's probably a good choice.
  int lo_q = -100, hi_q = 100, mid_q;
  do {
    mid_q = (lo_q + hi_q) / 2;
    size_t file_size;
    if (!CompressAV1(copy_partition, mid_q, &file_size)) {
      AddInfo("AV1 not available");
      return false;
    }
    if (std::abs((float)file_size - target_file_size) <
        0.001f * target_file_size) {
      break;
    } else if (file_size < target_file_size) {
      lo_q = mid_q + 1;
    } else {
      hi_q = mid_q;
    }
  } while (lo_q >= 0 && lo_q < hi_q);
  return true;
}

//------------------------------------------------------------------------------

void Params::ApplyZoomLevel() {
  if (zoom_level >= 0) {
    viewport_width = width << zoom_level;
    viewport_height = height << zoom_level;
  } else {
    viewport_width = width >> -zoom_level;
    viewport_height = height >> -zoom_level;
  }
}

bool Params::SetCurrentFile(size_t file) {
  if (file >= files.size()) file = files.size() - 1;
  if (file != current_file) {
    ClearForcedElements();
    current_file = file;
    const std::string& file_name = files[file];
    status = WP2::IoUtilReadFile(file_name.c_str(), &input);
    if (status != WP2_STATUS_OK) {
      fprintf(stderr, "Could not read the input file %s\n", file_name.c_str());
      return false;
    }
    const uint8_t* data = (const uint8_t*)input.data();
    const size_t data_size = input.size();

    const WP2::FileFormat format = WP2::GuessImageFormat(data, data_size);
    if (format == WP2::FileFormat::UNSUPPORTED) {
      fprintf(stderr, "Unsupported input format\n");
      return false;
    }
    if (WP2::ReadImage(data, data_size, &in, format) != WP2_STATUS_OK) {
      fprintf(stderr, "Could not decode the input file %s\n",
              file_name.c_str());
      return false;
    }

    if (crop && in.SetView(in, crop_area) != WP2_STATUS_OK) {
      fprintf(stderr, "Error! Cropping operation failed. Skipping.");
    }

    preview_color = {0x00u, 0x00u, 0x00u, 0x00u};
    preview.Deallocate();
    if (format == WP2::FileFormat::WP2) {
      WP2::BitstreamFeatures wp2_features;
      if (wp2_features.Read(data, data_size) != WP2_STATUS_OK) {
        fprintf(stderr, "Inconsistent state: could not read WP2 features.\n");
        return false;
      }
      if (wp2_features.has_preview &&
          ExtractPreview(data, data_size, &preview) != WP2_STATUS_OK) {
        fprintf(stderr, "Could not decode preview of %s\n",
                file_name.c_str());
        preview.Deallocate();
      }
      preview_color = WP2::ToArgb32b(wp2_features.preview_color);
    }
    width = in.width;
    height = in.height;

    if (!EncodeAndDecode()) return false;
  }
  return (status == WP2_STATUS_OK);
}

void Params::ShiftAlt() {
  const WP2Status s = alt1.Swap(&alt2);
  (void)s;
  assert(s == WP2_STATUS_OK);
  std::swap(alt1_msg, alt2_msg);
  alt1_msg.clear();
  alt1.Deallocate();
}

bool Params::SetBitstream(const char* const file_name) {
  std::string data;
  status = WP2::IoUtilReadFile(file_name, &data);
  if (status != WP2_STATUS_OK) {
    fprintf(stderr, "Could not read the file %s\n", file_name);
    return false;
  }

  WP2::ArgbBuffer tmp;
  status = WP2::Decode(data, &tmp);
  if (status != WP2_STATUS_OK) {
    fprintf(stderr, "Could not decode the file %s\n", file_name);
    return false;
  }

  ClearForcedElements();
  alt1_msg.clear();
  alt1.Deallocate();
  alt2_msg.clear();
  alt2.Deallocate();

  view_only = true;
  using std::swap;
  swap(data, input);
  data.clear();
  status = tmp.Swap(&in);
  return (status == WP2_STATUS_OK);
}

bool Params::SetAltFile(const char* const file_name) {
  ShiftAlt();
  status = WP2::ReadImage(file_name, &alt1);
  if (status != WP2_STATUS_OK) {
    fprintf(stderr, "Could not decode the alternate file %s\n", file_name);
    return false;
  }
  const uint32_t w = alt1.width, h = alt1.height;
  if (w != width || h != height) {
    alt1.Deallocate();
    fprintf(stderr, "Alternate picture has incompatible dimensions "
                    " (%dx%d vs expected %dx%d)\n",
            w, h, width, height);
    return false;
  }
  alt1_msg = SPrintf("- Alt Pic (%s) -", file_name);
  return true;
}

// Copies the current canvas to the alternative buffer.
bool Params::SetAltImage() {
  if (show == Params::kAlt1 || show == Params::kAlt2) return false;
  WP2::ArgbBuffer tmp;
  if (GetCurrentCanvas(&tmp)) {
    ShiftAlt();
    status = tmp.IsView() ? alt1.CopyFrom(tmp) : alt1.Swap(&tmp);
    if (status == WP2_STATUS_OK) alt1_msg = "- Saved canvas -";
  }
  if (status != WP2_STATUS_OK) {
    fprintf(stderr, "Could not copy the alternate image\n");
    return false;
  }
  return true;
}

bool Params::DumpCurrentCanvas(const char* const file_path) {
  WP2::ArgbBuffer tmp;
  if (GetCurrentCanvas(&tmp)) {
    status = WP2::SaveImage(tmp, file_path, /*overwrite=*/true);
  }

  if (status != WP2_STATUS_OK) {
    fprintf(stderr, "Could not dump the current canvas\n");
    return false;
  }
  return true;
}

// Returns the currently displayed buffer.
const WP2::ArgbBuffer& Params::GetBuffer() const {
  if (show == Params::kOriginal) return original;
  if (show == Params::kPreview) return preview;
  if (show == Params::kWebP) return webp;
  if (show == Params::kAlt1 && !alt1.IsEmpty()) return alt1;
  if (show == Params::kAlt2 && !alt2.IsEmpty()) return alt2;
  if (show == Params::kCompressed) return out_yuv.IsEmpty() ? out : out_yuv;
  if (WP2::VDMatch(dec_config, "")) return dec_config.info->debug_output;
  if (WP2::VDMatch(enc_config, "")) return enc_config.info->debug_output;
  return out;
}

// Copies the current canvas to 'buffer'.
bool Params::GetCurrentCanvas(WP2::ArgbBuffer* const buffer) {
  if (show == Params::kPreviewColor) {
    ShiftAlt();
    status = buffer->Resize(out.width, out.height);
    if (status == WP2_STATUS_OK) buffer->Fill(preview_color);
  } else {
    status = buffer->SetView(GetBuffer());
  }
  return (status == WP2_STATUS_OK);
}

//------------------------------------------------------------------------------
// Blocks

// Same as glRectf() but with 'rect' in pixels, top-left being (0, 0).
void glRect(const WP2::Rectangle& rect, uint32_t thickness = 1) {
  const float x = 2.f * rect.x / kParams.width - 1.f;
  // Compensate OpenGL bottom-left origin by shifting by one pixel vertically.
  const float y = 1.f - 2.f * (rect.y + rect.height) / kParams.height -
       1.f / kParams.viewport_height;
  const float w = 2.f * rect.width / kParams.width;
  const float h = 2.f * rect.height / kParams.height;

  // Thickness
  const float px_w = 2.f / kParams.viewport_width;
  const float px_h = 2.f / kParams.viewport_height;
  for (uint32_t j = 0; j < thickness; ++j) {
    for (uint32_t i = 0; i < thickness; ++i) {
      glRectf(x + i * px_w, y - j * px_h, x + w + i * px_w, y + h - j * px_h);
    }
  }
}

// Helper class for the creation of a nice looking paragraph.
class ParagraphMaker {
 public:
  // Initializes with a potential prefix that will define an indentation for the
  // rest of the paragraph.
  void Init(const std::string& str = "") {
    strs_.clear();
    strs_.push_back(str);
    indentation_ = std::string(str.length(), ' ');
  }
  // Appends a string to the last string of the list.
  void Append(const std::string& str) {
    assert(!strs_.empty());
    strs_.back() += str;
  }
  // Appends a new element to the list.
  void AppendToList(const std::string& str) {
    assert(!strs_.empty());
    strs_.back() += str + ", ";
    if (strs_.back().size() >= 50u) strs_.push_back(indentation_);
  }
  // Push back all the elements of the list to the msg.
  void PushBackTo(std::vector<std::string>* const msg) {
    if (strs_.empty()) return;
    // Remove the last element if empty.
    if (strs_.back() == indentation_) strs_.pop_back();
    if (strs_.empty()) return;
    // Remove the last ', '.
    const uint32_t index = strs_.back().size() - 2;
    if (strs_.back().substr(index) == ", ") {
      strs_.back() = strs_.back().substr(0, index);
    }
    msg->insert(msg->end(), strs_.begin(), strs_.end());
  }

 private:
  std::vector<std::string> strs_;
  std::string indentation_;
};

void Params::DisplayBlockInfo(std::vector<std::string>* const msg) {
  if (display_block == kDisplayNone) return;
#if defined(WP2_BITTRACE)
  assert(dinfo.store_blocks && einfo.store_blocks &&
         einfo.blocks.size() == dinfo.blocks.size());
  // Sort encoded and decoded blocks for easy comparison.
  for (std::vector<WP2::BlockInfo>* blocks : {&dinfo.blocks, &einfo.blocks}) {
    std::sort(blocks->begin(), blocks->end(),
              [](const WP2::BlockInfo& l, const WP2::BlockInfo& r) {
                return (l.rect.x < r.rect.x) ||
                       (l.rect.x == r.rect.x && l.rect.y < r.rect.y);
              });
  }

  const char* const kEncNames[] = {"X", "Mtd0", "Mtd1", "DC  ", "Zero"};
  const int selection_x = (int)dinfo.selection.x;
  const int selection_y = (int)dinfo.selection.y;
  for (uint32_t index = 0; index < dinfo.blocks.size(); ++index) {
    const WP2::BlockInfo& b = dinfo.blocks[index];
    const WP2::BlockInfo& b_enc = einfo.blocks[index];
    assert(b.rect.x == b_enc.rect.x && b.rect.width == b_enc.rect.width &&
           b.rect.y == b_enc.rect.y && b.rect.height == b_enc.rect.height);
    assert(b.segment_id == b_enc.segment_id);
    assert(b.y_context_is_constant == b_enc.y_context_is_constant);

    const int bx = b.rect.x, by = b.rect.y;
    const int bw = b.rect.width, bh = b.rect.height;

    const bool block_selected = (bx <= selection_x && selection_x < bx + bw &&
                                 by <= selection_y && selection_y < by + bh);
    if (display_block == kDisplayGrid || !block_selected) {
      glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
      glColor4f(0.0f, 0.2f, 0.1f, 0.2f);
      // +1 width and height for a thinner border betweens blocks.
      glRect({b.rect.x, b.rect.y, b.rect.width, b.rect.height});
      continue;
    }

    bool forced = false;
    for (const WP2::Rectangle& rect : einfo.force_partition) {
      if (rect.x == b.rect.x && rect.y == b.rect.y) {
        assert(rect.width == b.rect.width && rect.height == b.rect.height);
        forced = true;
        break;
      }
    }
    msg->push_back(SPrintf("size: %2d x %2d  pos: %3d x %3d  %s",
                           bw, bh, bx, by, forced ? "(forced)" : ""));
    float disto[5];
    const WP2::Rectangle rect{(uint32_t)bx, (uint32_t)by, (uint32_t)bw,
                              (uint32_t)bh};
    WP2_ASSERT_STATUS(kParams.out.GetDistortionBlackOrWhiteBackground(
        kParams.in, rect, kParams.distortion_metric, disto));
    msg->push_back(
        SPrintf("bits:%5.1f (%.2f bpp) %s:%5.1f", b.bits, b.bits / (bw * bh),
                kWP2MetricNames[kParams.distortion_metric], disto[4]));

    // Isolate the bit costs about residuals (Y, U, V, A, remaining stuff).
    static constexpr uint32_t kAll = 4;
    std::map<const std::string, WP2::LabelStats> bit_traces[5];
    double bit_sums[5] = {0};

    for (const auto& bt : b.bit_traces) {
      std::string label = bt.first;
      uint32_t channel;
      for (channel = 0; channel < kAll; ++channel) {
        const std::string prefix = WP2::kCoeffsStr[channel];

        if (label.size() >= prefix.size() &&
            std::strncmp(label.c_str(), prefix.c_str(), prefix.size()) == 0) {
          // 'label' starts with 'prefix', remove 'prefix' from it.
          assert(label.size() > prefix.size() + 1 &&
                 label[prefix.size()] == '/');
          label = label.substr(prefix.size() + 1);
          break;
        }
      }

      // Merge by first token.
      const size_t slash_pos = label.find('/');
      if (slash_pos != std::string::npos) label = label.substr(0, slash_pos);

      bit_traces[channel][label].bits += bt.second.bits;
      bit_traces[channel][label].num_occurrences += bt.second.num_occurrences;
      bit_sums[channel] += bt.second.bits;
    }

    ParagraphMaker pm;
    if (display_block == kDisplayHeader) {
      pm.Init("Costs: ");
      for (uint32_t c = 0; c < (b.has_lossy_alpha ? kAll : 3); ++c) {
        pm.AppendToList(SPrintf("%s: %.1f", WP2::kCoeffsStr[c], bit_sums[c]));
      }
      for (const auto& p : bit_traces[kAll]) {
        pm.AppendToList(SPrintf("%s: %.1f", p.first.c_str(), p.second.bits));
      }
      pm.PushBackTo(msg);
      msg->push_back(SPrintf("segment id: %d (%s %s)", b.segment_id,
                             b.is420 ? "is420" : "",
                             b.has_lossy_alpha ? "has_alpha_res" : ""));
      msg->push_back(SPrintf("transform: %s %s",
                             WP2TransformNames[b.tf_x],
                             WP2TransformNames[b.tf_y]));
      for (uint32_t tf_i = 0; tf_i < 4; ++tf_i) {
        if (b.encoding_method[WP2::kYChannel][tf_i] == -1) continue;
        msg->push_back(
            SPrintf("Transform %u. Mthd: Y=%s U=%s V=%s (A=%s)", tf_i,
                    kEncNames[b.encoding_method[WP2::kYChannel][tf_i] + 1],
                    kEncNames[b.encoding_method[WP2::kUChannel][tf_i] + 1],
                    kEncNames[b.encoding_method[WP2::kVChannel][tf_i] + 1],
                    kEncNames[b.encoding_method[WP2::kAChannel][tf_i] + 1]));
      }
    } else if (display_block <= kDisplayACoeffs) {
      if (display_block == kDisplayACoeffs && !b.has_lossy_alpha) {
        msg->push_back("- no alpha -");
      } else {
        const uint32_t channel = (uint32_t)(display_block - kDisplayYCoeffs);
        msg->push_back(SPrintf("%s channel.", WP2::kChannelStr[channel]));
        for (uint32_t tf_i = 0; tf_i < 4; ++tf_i) {
          if (!b.residual_info[channel][tf_i].empty()) {
            pm.Init();
            // Display info about the residuals (index of last element ...).
            pm.Append(SPrintf("Transform %u. Residuals: ", tf_i));
            for (const std::string& str : b.residual_info[channel][tf_i]) {
              pm.AppendToList(str);
            }
            pm.PushBackTo(msg);
          }
        }
        // Display stats about residuals.
        {
          typedef decltype(b.bit_traces)::value_type StatType;
          std::vector<StatType> stats(bit_traces[channel].begin(),
                                      bit_traces[channel].end());
          // Sort those stats from lowest to highest usage.
          std::vector<uint32_t> indices(stats.size());
          std::iota(indices.begin(), indices.end(), 0);
          std::sort(indices.begin(), indices.end(),
                    [stats](uint32_t i1, uint32_t i2) {
                      return stats[i1].second.bits < stats[i2].second.bits;
                    });
          pm.Init("Costs: ");
          pm.Append(SPrintf("all: %.1f, ", bit_sums[channel]));
          // Display stats as bit cost (/num of occurences = individual bit
          // cost).
          for (uint32_t i : indices) {
            const auto& p = stats[i];
            pm.AppendToList(
                SPrintf("%s: %.1f (/%d=%.1f)", p.first.c_str(), p.second.bits,
                        p.second.num_occurrences,
                        p.second.bits / p.second.num_occurrences));
          }
          pm.PushBackTo(msg);
        }
        for (uint32_t tf_i = 0; tf_i < 4; ++tf_i) {
          if (b.encoding_method[channel][tf_i] == -1) continue;
          msg->push_back(
              SPrintf("Enc %u. Mthd: %s", tf_i,
                      kEncNames[b.encoding_method[channel][tf_i]] + 1));
          const uint32_t scale =
              (b.is420 && (channel == 1 || channel == 2)) ? 2 : 1;
          const WP2::BlockSize split_size =
              WP2::GetSplitSize(WP2::GetBlockSize(bw / WP2::kMinBlockSizePix,
                                                  bh / WP2::kMinBlockSizePix),
                                b.split_tf[channel]);
          const uint32_t w = WP2::BlockWidthPix(split_size) / scale;
          const uint32_t h = WP2::BlockHeightPix(split_size) / scale;
          const bool show_original = (kParams.show == Params::kOriginal);
          if (show_original || b.encoding_method[channel][tf_i] !=
                                   (int8_t)WP2::EncodingMethod::kAllZero) {
            int32_t res[WP2::kMaxBlockSizePix2];
            const uint32_t num_coeffs = w * h;
            std::copy(&b_enc.original_res[channel][tf_i][0],
                      &b_enc.original_res[channel][tf_i][0] + num_coeffs,
                      &res[0]);

            const bool reduced_transform =
                (channel == WP2::kUChannel || channel == WP2::kVChannel) &&
                b_enc.is420;
            WP2Transform2D(res, (WP2TransformType)b_enc.tf_x,
                           (WP2TransformType)b_enc.tf_y, w, h, res,
                           reduced_transform);

            const int32_t norm = (int32_t)num_coeffs;
            for (uint32_t j = 0; j < h; ++j) {
              std::string line;
              for (uint32_t i = 0; i < w; ++i) {
                const uint32_t idx = i + w * j;
                if (show_original) {
                  const int32_t v = res[idx];
                  if (std::abs(v / norm) > 0) {
                    line += SPrintf("%3d ", v / norm);
                  } else {
                    const char symbol = (std::abs(v) >= (norm * 0.50)) ? 'X'
                                      : (std::abs(v) >= (norm * 0.25)) ? 'x'
                                                                       : '.';
                    line += SPrintf("  %c ", symbol);
                  }
                } else {
                  const int32_t v = b.coeffs[channel][tf_i][idx];
                  line += v ? SPrintf("%3d ", v) : "  . ";
                }
              }
              msg->push_back(line);
            }
            if (show_original) {
              msg->push_back(SPrintf("original coeffs in units of %u", norm));
            }
          } else {
            msg->push_back(" - no coeffs -");
          }
        }
      }
    } else if (display_block == kDisplayPredModes) {
      msg->push_back(SPrintf("Chroma prediction mode: %d", b.uv_pred));
      msg->push_back(SPrintf("Luma prediction: %d", b.y_pred));
      if (b.has_lossy_alpha) {
        msg->push_back(SPrintf("Alpha prediction: %d", b.a_pred));
      } else {
        msg->push_back("-no alpha-");
      }

      msg->push_back(b.y_context_is_constant ? "Y context is constant"
                                             : "Y context is not constant");
      const int optimize_modes =
          (b.y_context_is_constant || enc_config.speed == 0) ? 0 :
                                     (enc_config.speed == 1) ? 1 : 2;
      msg->push_back(SPrintf("OptimizeModes%d() Y "
                             "score: %.1f (%.1f per pixel)",
                             optimize_modes,
                             b_enc.pred_scores[WP2::kYChannel],
                             b_enc.pred_scores[WP2::kYChannel] / (bw * bh)));
      msg->push_back(SPrintf("FindBestUVModes() score: %.1f (%.1f per pixel)",
                             b_enc.pred_scores[WP2::kUChannel],
                             b_enc.pred_scores[WP2::kUChannel] / (bw * bh)));
    }
    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
    glColor4f(1.f, 0.f, 0.f, 0.5f);  // red
    glRect(b.rect);
  }
#endif  // WP2_BITTRACE
}

void Params::ClearForcedElements() {
  force_partition.clear();
  force_segment.clear();
  force_predictor.clear();
  forcing_block = WP2::Rectangle();
}

bool Params::IsDebugViewVisible() const {
  return (show == kDebug || show == kInfo || show == kHelp || show == kMenu);
}

bool Params::IsPartitionVisible() const {
  return (display_block != kDisplayNone ||
          (IsDebugViewVisible() && WP2::VDMatch(dec_config, "partition"))) &&
         !IsSegmentVisible() && !IsPredictorVisible();
}

bool Params::IsSegmentVisible() const {
  return IsDebugViewVisible() && WP2::VDMatch(dec_config, "segment-ids");
}

bool Params::IsPredictorVisible() const {
  return IsDebugViewVisible() && WP2::VDMatch(dec_config, "prediction");
}

void Params::DisplayPartition() const {
  if (!IsPartitionVisible()) return;

  // Adapt to zoom.
  const uint32_t thick = getLineThickness();
  glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);

  // Encoded forced blocks.
  glColor4f(1.0f, 0.2f, 1.0f, 1.0f);
  for (const WP2::Rectangle& rect : einfo.force_partition) glRect(rect, thick);

  glLineStipple(2 * thick, 0xAAAA);
  glEnable(GL_LINE_STIPPLE);

  // All blocks, including those that haven't been applied yet.
  glColor4f(0.6f, 1.0f, 0.6f, 1.0f);
  for (const WP2::Rectangle& rect : force_partition) glRect(rect, thick);

  // If a block is currently being drawn with the mouse.
  if (forcing_block.width > 0 && forcing_block.height > 0) {
    glColor4f(1.0f, 1.0f, 1.0f, 1.0f);
    glRect(forcing_block, thick);
  }

  glDisable(GL_LINE_STIPPLE);
}

void Params::DisplayForcedSegments() const {
  if (!IsSegmentVisible()) return;

#if defined(WP2_BITTRACE)
  glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);

  const uint32_t thick = getLineThickness();
  glLineStipple(2 * thick, 0xAAAA);

  for (const auto& segment : force_segment) {
    for (const WP2::BlockInfo& block : einfo.blocks) {
      if (block.rect.Contains(segment.x, segment.y)) {
        const WP2::Rectangle rect{block.rect.x, block.rect.y,
                                  block.rect.width - 1, block.rect.height - 1};
        glColor4f(0.f, 0.f, 0.f, 1.f);
        glRect(rect, thick * 2);
        const WP2::Argb32b color = WP2::kSegmentColors[segment.segment_id];
        glColor4f(color.r / 255.f, color.g / 255.f, color.b / 255.f, 1.f);
        glRect(rect, thick);
        break;
      }
    }
  }
#endif
}

#ifdef WP2_BITTRACE
static const std::vector<std::string>& GetPredictors(
    const WP2::BlockInfo& block, WP2::Channel channel) {
  const WP2::DecoderInfo& info = *kParams.dec_config.info;
  return (channel == WP2::kYChannel)
             ? info.y_predictors
             : (channel == WP2::kAChannel) ? info.a_predictors
                                           : info.uv_predictors;
}
#endif

void Params::DisplayForcedPredictors() const {
  if (!IsPredictorVisible()) return;

#if defined(WP2_BITTRACE)
  WP2::Channel channel = WP2::VDChannel(kParams.dec_config);
  if (channel == WP2::kVChannel) channel = WP2::kUChannel;

  glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);

  for (const auto& predictor : force_predictor) {
    if (predictor.channel != channel) continue;
    for (const WP2::BlockInfo& block : einfo.blocks) {
      if (block.rect.Contains(predictor.x, predictor.y)) {
        if (channel == WP2::kYChannel && block.y_context_is_constant) {
          glColor4f(1.f, 0.f, 0.f, 1.f);
        } else {
          glColor4f(0.f, 1.f, 0.f, 1.f);
        }
        glRect(block.rect, 4.f);

        const std::vector<std::string>& predictors =
            GetPredictors(block, channel);
        const uint32_t predictor_id =
            std::min(predictor.predictor_id, (uint8_t)(predictors.size() - 1));

        const float color = (float)predictor_id / (predictors.size() - 1);
        glColor4f(color, color, color, 1.f);
        glRect(block.rect, 2.f);

        const std::string& title = predictors[predictor_id];
        glColor4f(1.f, 0.f, 1.f, 1.f);
        PrintString((predictor.x) * 2.f / kParams.width - 1,
                    1 - (predictor.y) * 2.f / kParams.height, title);

        break;
      }
    }
  }
#endif
}

uint32_t Params::getLineThickness() const {
  return (display_block != kDisplayNone)
             ? 1u
             : std::max(1u, std::min(viewport_width / width,
                                     viewport_height / height));
}

bool Params::IsHamburgerMenuVisible() const {
  return (show == kDebug || show == kMenu ||
          (show == kOriginal && display_block != kDisplayNone)) &&
         (forcing_block.width == 0 || forcing_block.height == 0);
}

//------------------------------------------------------------------------------

bool Params::ReadPartition(bool read_only) {
  force_partition.clear();
  if (!WP2ReadPartition(partition_file_path, read_only, &force_partition)) {
    return false;
  }
  WP2ConvertPartition(in.width, in.height, enc_config.partition_set,
                      /*ignore_invalid=*/false, &force_partition);
  return true;
}

//------------------------------------------------------------------------------
// Handlers

// T must be cast-able to and from 'int'
template <class T>
void Modify(T& value, int max, bool increment, int increment_step = 1) {
  int v = (int)value + (increment ? increment_step : -increment_step);
  v = (v + max) % max;  // deals with negative values
  value = (T)v;
}
// For ranges
template <class T>
void ModifyR(T& value, T max, bool increment, int increment_step = 1) {
  if (glutGetModifiers() == GLUT_ACTIVE_ALT) {
    // Toggle on/off (GLUT doesn't handle shift+alt+key: can't use 'increment')
    value = (T)(((int)value > (int)max / 2) ? 0 : (int)max - 1);
  } else {
    Modify(value, max, increment, increment_step);
  }
}

// Reshapes the window based on current image size and zoom level.
void ReshapeWindow() {
  kParams.ApplyZoomLevel();
  glutReshapeWindow(kParams.viewport_width, kParams.viewport_height);
  // Update viewport_width/height based on actual window dimensions. These can
  // be different from the ones we requested because of the size of the screen
  // or because the window is in fullscreen mode.
  kParams.viewport_width = glutGet(GLUT_WINDOW_WIDTH);
  kParams.viewport_height = glutGet(GLUT_WINDOW_HEIGHT);
  glutPostRedisplay();
}

void HandleKey(unsigned char key, int pos_x, int pos_y) {
  kParams.mouse_x = pos_x;
  kParams.mouse_y = pos_y;

  // Ignore key strokes if the menu is opened.
  if (kParams.show == Params::kMenu) return;

  const bool is_ctrl = (glutGetModifiers() == GLUT_ACTIVE_CTRL);
  const bool is_alt = (glutGetModifiers() == GLUT_ACTIVE_ALT);
  (void)is_ctrl;

  // first, the keys with special modifiers
  if (is_alt && key >= '0' && key <= '9') {
    if (key == '1' && !kParams.alt1.IsEmpty()) {
      kParams.show = Params::kAlt1;
    } else if (key == '2' && !kParams.alt2.IsEmpty()) {
      kParams.show = Params::kAlt2;
    }
    glutPostRedisplay();
  } else if (key == 'q' || key == 'Q' || key == 27 /* Esc */) {
    if (key == 27 /* Esc */ && (kParams.einfo.selection.width > 0 ||
                                kParams.einfo.selection.height > 0)) {
      kParams.einfo.selection = kParams.dinfo.selection = {0, 0, 0, 0};
      if (WP2::VDMatch(kParams.enc_config, "encoder")) {
        if (!kParams.EncodeAndDecode()) return;
      } else {
        if (!kParams.DecodeOutput()) return;
      }
      AddInfo("Selection cleared. Press escape again to quit.");
    } else {
#if defined(WP2_HAVE_GLUT) && defined(FREEGLUT)
      glutLeaveMainLoop();
#else
      exit(0);
#endif
    }
  } else if (key == ' ') {
    kParams.show = Params::kOriginal;
    if (!kParams.GetOriginal()) AddInfo("An error occurred");
    glutPostRedisplay();
  } else if (key == 'W' || key == 'w') {
    if (!kParams.CompressWebP(key == 'w')) AddInfo("WebP not available");
    glutPostRedisplay();
  } else if (key == 'r' || key == 'R') {
    const bool binary_search = is_alt;
    const bool copy_partition = (key == 'R');
    if (binary_search
            ? kParams.CompressAV1ToMatch(copy_partition, kParams.output.size())
            : kParams.CompressAV1(copy_partition, kParams.quality)) {
      AddInfo(copy_partition ? "Copied partition (press 'c' to use it)"
                             : "Stored AV1 as alt-1 image");
    } else {
      AddInfo("AV1 not available");
    }
    glutPostRedisplay();
  } else if (key == '\\' && !kParams.preview.IsEmpty()) {
    kParams.show = Params::kPreview;
    glutPostRedisplay();
  } else if (key == '|' && kParams.preview_color.a == 0xFFu) {
    kParams.show = Params::kPreviewColor;
    glutPostRedisplay();
  } else if (key == 13) {  // return
    if (!kParams.alt1.IsEmpty()) {
      kParams.show = Params::kAlt1;
      glutPostRedisplay();
    }
  } else if (key == '\t') {
    kParams.show = Params::kCompressed;
    if (!kParams.GetCompressed()) AddInfo("Unable to display compressed image");
    glutPostRedisplay();
  } else if (key == 'a' || key == 'A') {
    if (is_alt) {
      if (kParams.SetAltFile(kParams.dump_png_path)) {
        AddInfo(SPrintf("Loaded alt-1 from '%s'", kParams.dump_png_path));
      } else {
        AddInfo(SPrintf("Could not load '%s'", kParams.dump_png_path));
      }
    } else if (key == 'A') {
      if (kParams.DumpCurrentCanvas(kParams.dump_png_path)) {
        AddInfo(SPrintf("Dumped canvas to '%s'", kParams.dump_png_path));
      } else {
        AddInfo(SPrintf("Could not dump to '%s'", kParams.dump_png_path));
      }
    } else {
      if (kParams.SetAltImage()) {
        AddInfo("Saved canvas to alt-1 image");
      } else {
        AddInfo("Unable to save canvas");
      }
      glutPostRedisplay();
    }
  } else if (key == 'b' || key == 'B') {
    const bool need_update = (kParams.dinfo.store_blocks == 0);
    Modify(kParams.display_block, Params::kDisplayNum, key == 'b');
    kParams.dinfo.store_blocks =
        (kParams.display_block != Params::kDisplayNone);
    if (need_update) {
      if (!kParams.DecodeOutput()) return;
    }
    glutPostRedisplay();
  } else if (key == 'E' || key == 'e') {
    std::string selection =
        GetLeaf(kParams.menu_selection, /*leaf_offset=*/(key == 'e') ? 1 : -1);
    if (selection.find("error-map") != 0) {  // NOLINT
      selection = "error-map/PSNR";
    }
    SetMenuSelection(selection);
    if (!kParams.DecodeOutput()) return;
  } else if (key == 's' || key == 'S') {
    ModifyR(kParams.enc_config.speed, 10, key == 's');
    if (!kParams.EncodeAndDecode()) return;
    SetMessage(SPrintf("Speed: %d", kParams.enc_config.speed));
  } else if (key == 'o' || key == 'O') {
    ModifyR(kParams.background, Params::kBackgroundNum, key == 'o');
    const char* const kBackgroundName[] = {"Checkerboard", "Clearly White",
                                           "Just Black", "Not Pink"};
    STATIC_ASSERT_ARRAY_SIZE(kBackgroundName, Params::kBackgroundNum);
    SetMessage(SPrintf("Background: %s", kBackgroundName[kParams.background]));
  } else if (key == 'l') {
    if (kParams.SetBitstream(kParams.dump_wp2_path)) {
      AddInfo(SPrintf("Loaded file from '%s'", kParams.dump_wp2_path));
      if (!kParams.EncodeAndDecode()) return;  // 'view_only' so just copy.
    } else {
      AddInfo(SPrintf("Couldn't load file from '%s'", kParams.dump_wp2_path));
    }
  } else if (key == 'L') {
    if (kParams.output.empty() ||
        WP2::IoUtilWriteFile(kParams.output, kParams.dump_wp2_path,
                             /*overwrite=*/true) != WP2_STATUS_OK) {
      AddInfo(SPrintf("Couldn't save dump file '%s'", kParams.dump_wp2_path));
    } else {
      AddInfo(SPrintf("Saved bits as '%s'", kParams.dump_wp2_path));
    }
  } else if (key == 'M' || key == 'm') {
    kParams.enc_config.use_random_matrix =
        !kParams.enc_config.use_random_matrix;
    if (!kParams.EncodeAndDecode()) return;
    SetMessage(
        SPrintf("Use random matrix: %d", kParams.enc_config.use_random_matrix));
  } else if (key == 'x' || key == 'X') {
    Modify(kParams.wp2dbg_value, 2u, key == 'x');
    setenv("WP2DBG", (kParams.wp2dbg_value == 0) ? "0" : "1", 1);
    if (!kParams.EncodeAndDecode()) return;
    AddInfo(SPrintf("export WP2DBG=%d", kParams.wp2dbg_value));
  } else if (key == 'f' || key == 'F') {
    if (glutGetModifiers() == GLUT_ACTIVE_ALT) {
      kParams.dec_config.enable_restoration_filter =
          !kParams.dec_config.enable_restoration_filter;
    } else if (glutGetModifiers() == (GLUT_ACTIVE_ALT | GLUT_ACTIVE_SHIFT)) {
      kParams.enc_config.enable_alpha_filter =
          !kParams.enc_config.enable_alpha_filter;
    } else if (key == 'F') {
      kParams.dec_config.enable_directional_filter =
          !kParams.dec_config.enable_directional_filter;
    } else {
      kParams.dec_config.enable_deblocking_filter =
          !kParams.dec_config.enable_deblocking_filter;
    }
    if (!kParams.EncodeAndDecode()) return;
    SetMessage(
        SPrintf("Filters: dblk=%s, drct=%s, rstr=%s, alpha=%s",
                kOnOff[kParams.dec_config.enable_deblocking_filter],
                kOnOff[kParams.dec_config.enable_directional_filter],
                kOnOff[kParams.dec_config.enable_restoration_filter],
                kOnOff[kParams.enc_config.enable_alpha_filter]));
  } else if (key == 'g' || key == 'G') {
    if (key == 'G') {
      kParams.dec_config.grain_amplitude = 0;
    } else {
      ModifyR(kParams.dec_config.grain_amplitude, (uint8_t)101, false, 20);
    }
    kParams.enc_config.store_grain = (kParams.dec_config.grain_amplitude > 0);
    if (!kParams.EncodeAndDecode()) return;
    SetMessage(
        SPrintf("Grain amplitude: %d", kParams.dec_config.grain_amplitude));
  } else if (key == 'c' || key == 'C') {
    if (!kParams.EncodeAndDecode()) return;
    AddInfo(SPrintf("Forcing %zu blocks %zu segments %zu predictors",
                    kParams.einfo.force_partition.size(),
                    kParams.einfo.force_segment.size(),
                    kParams.einfo.force_predictor.size()));
  } else if (key == 'd') {
    EXIT_IF_FALSE(
        WP2WritePartition(kParams.force_partition, kParams.partition_file_path),
        "Error: Unable to write partition file");
    AddInfo(SPrintf("Dumped %zu blocks", kParams.force_partition.size()));
  } else if (key == 'D') {
    if (kParams.ReadPartition(/*read_only=*/true)) {
      AddInfo(SPrintf("Loaded %zu blocks", kParams.force_partition.size()));
    } else {
      AddInfo("Error: Unable to parse partition file");
    }
  } else if (key == '+' || key == '_') {
    const int file = (int)kParams.current_file + (key == '+' ? 1 : -1);
    if (file >= 0 && file < (int)kParams.files.size()) {
      if (kParams.SetCurrentFile(file)) {
        ReshapeWindow();
      } else {
        AddInfo(SPrintf("Unable to read file %d", file));
      }
    } else {
      AddInfo(SPrintf("File %d is not referenced", file));
    }
  } else if (key == 'i' || key == 'I') {
    if (kParams.show != Params::kInfo) {
      kParams.show = Params::kInfo;
      kParams.show_interface = true;
      kParams.display_block = Params::kDisplayNone;
      if (key == 'I') {
        if (!kParams.ComputeYUVDistortion()) return;
      } else {
        kParams.yuv_distortion[0] = 0.f;
        kParams.yuv_distortion[1] = 0.f;
        kParams.yuv_distortion[2] = 0.f;
      }
    } else {
      kParams.show = Params::kDebug;
    }
    glutPostRedisplay();
  } else if (key == 't' || key == 'T') {
    Modify(kParams.enc_config.segment_id_mode, 3u, key == 't');
    if (!kParams.EncodeAndDecode()) return;
    SetMessage(SPrintf("use_segment_id_mode: %s",
                       kSegmentModes[kParams.enc_config.segment_id_mode]));
  } else if (key == 'n' || key == 'N') {
    if (kParams.show != Params::kInfo) {
      kParams.show = Params::kInfo;
      kParams.show_interface = true;
      kParams.display_block = Params::kDisplayNone;
    }
    Modify(kParams.distortion_metric, WP2::NUM_METRIC_TYPES, key == 'n');
    if (!kParams.DecodeOutput()) return;
    glutPostRedisplay();
  } else if (key == 'h') {
    if (kParams.show != Params::kHelp) {
      kParams.show = Params::kHelp;
      kParams.show_interface = true;
      kParams.display_block = Params::kDisplayNone;
    } else {
      kParams.show = Params::kDebug;
    }
    glutPostRedisplay();
  } else if (key == 'H') {
    kParams.show_interface = !kParams.show_interface;
    glutPostRedisplay();
  } else if (key == '1' || key == '!') {
    Modify(kParams.enc_config.csp_type, WP2::kNumCspTypes, key == '1');
    if (!kParams.EncodeAndDecode()) return;
    SetMessage(SPrintf("CSP type: %s",
                       kCSPString[(uint32_t)kParams.enc_config.csp_type]));
  } else if (key == '2' || key == '@') {
    Modify(kParams.enc_config.uv_mode, WP2::EncoderConfig::NumUVMode,
           key == '2');
    if (!kParams.EncodeAndDecode()) return;
    SetMessage(
        SPrintf("UV Mode: %s", kUVModeString[kParams.enc_config.uv_mode]));
  } else if (key == '3' || key == '#') {
    Modify(kParams.enc_config.tile_shape, WP2::NUM_TILE_SHAPES,
           key == '3');
    if (!kParams.EncodeAndDecode()) return;
    SetMessage(SPrintf("Tile shape: %s",
                       kTileShapeString[kParams.enc_config.tile_shape]));
  } else if (key == '4' || key == '$') {
    Modify(kParams.enc_config.partition_method, WP2::NUM_PARTITION_METHODS,
           key == '4');
    if (!kParams.EncodeAndDecode()) return;
    SetMessage(SPrintf(
        "Partition method: %s",
        WP2::kPartitionMethodString[kParams.enc_config.partition_method]));
  } else if (key == '5' || key == '%') {
    Modify(kParams.enc_config.partition_set, WP2::NUM_PARTITION_SETS,
           key == '5');
    if (!kParams.EncodeAndDecode()) return;
    SetMessage(SPrintf("Partition set: %s",
                       kPartitionSetString[kParams.enc_config.partition_set]));
  } else if (key == '6' || key == '^') {
    kParams.enc_config.partition_snapping =
        !kParams.enc_config.partition_snapping;
    if (!kParams.EncodeAndDecode()) return;
    SetMessage(SPrintf("Partition snapping: %s",
                       kOnOff[kParams.enc_config.partition_snapping]));
  } else if (key == '7' || key == '&') {
    kParams.enc_config.segments -= 1;  // 'segments' value is 1-based
    ModifyR(kParams.enc_config.segments, (int)WP2::kMaxNumSegments, key == '7');
    kParams.enc_config.segments += 1;
    kParams.force_segment.clear();
    if (!kParams.EncodeAndDecode()) return;
    SetMessage(SPrintf("Num segments: %d", kParams.enc_config.segments));
  } else if (key == '8' || key == '*') {
    kParams.enc_config.sns += (key == '8') ? 7.f : -6.f;
    kParams.enc_config.sns = WP2::Clamp(kParams.enc_config.sns, 0.f, 100.f);
    if (!kParams.EncodeAndDecode()) return;
    SetMessage(SPrintf("SNS: %.1f", kParams.enc_config.sns));
  } else if (key == '9' || key == '(') {
    kParams.enc_config.error_diffusion += (key == '9') ? 7 : -6;
    kParams.enc_config.error_diffusion =
        WP2::Clamp(kParams.enc_config.error_diffusion, 0, 100);
    if (!kParams.EncodeAndDecode()) return;
    SetMessage(
        SPrintf("Error diffusion: %d", kParams.enc_config.error_diffusion));
  } else if (key == '0' || key == ')') {
    kParams.enc_config.tune_perceptual = !kParams.enc_config.tune_perceptual;
    if (!kParams.EncodeAndDecode()) return;
    SetMessage(SPrintf("Perceptual: %s",
                       kOnOff[kParams.enc_config.tune_perceptual]));
  } else if (key == 'k' || key == 'K') {
    kParams.keep_metadata = !kParams.keep_metadata;
    if (kParams.in.metadata.IsEmpty()) {
      SetMessage("Metadata: none");
    } else {
      if (!kParams.EncodeAndDecode()) return;
      SetMessage(kParams.keep_metadata ? "Metadata: on" : "Metadata: off");
    }
  } else if (key == 'v' || key == 'V') {
    SetMenuSelection(
        GetLeaf(kParams.menu_selection, /*leaf_offset=*/(key == 'v') ? 1 : -1));
    if (WP2::VDMatch(kParams.enc_config, "encoder")) {
      if (!kParams.EncodeAndDecode()) return;
    } else {
      if (!kParams.DecodeOutput()) return;
    }
  } else if (key == 'z' || key == 'Z') {
    kParams.zoom_level += (key == 'z') ? +1 : -1;  // Zoom in or out.
    ReshapeWindow();
  } else if (key == '?') {
    kParams.display_sparks = !kParams.display_sparks;
  }
}

void HandleKeyUp(unsigned char key, int pos_x, int pos_y) {
  (void)pos_x;
  (void)pos_y;
  const bool is_alt = (glutGetModifiers() == GLUT_ACTIVE_ALT);
  if (key == ' ' || key == '\\' || key == '|' || key == 13 || key == '\t') {
    kParams.show = Params::kDebug;
    kParams.out_yuv.Deallocate();
    glutPostRedisplay();
  } else if (key == 'w' || key == 'W') {
    kParams.show = Params::kDebug;
    kParams.webp.Deallocate();
    glutPostRedisplay();
  } else if (key == 'R' || key == 'r') {
    kParams.show = Params::kDebug;
    glutPostRedisplay();
  } else if (is_alt && (key == '1' || key == '2')) {
    kParams.show = Params::kDebug;
    glutPostRedisplay();
  }
}

uint32_t GetSqDist(uint32_t a, uint32_t b) {
  return (a < b) ? ((b - a) * (b - a)) : ((a - b) * (a - b));
}
uint32_t GetDist(const WP2::Rectangle& rect, WP2::BlockSize size) {
  return GetSqDist(rect.width, WP2::BlockWidthPix(size)) +
         GetSqDist(rect.height, WP2::BlockHeightPix(size));
}

void HandleMouseMove(int x, int y) {
  kParams.mouse_x = x;
  kParams.mouse_y = y;
}

void HandleMouseMoveButtonPressed(int x, int y) {
  kParams.mouse_x = x;
  kParams.mouse_y = y;
  if (x != kParams.mouse_x_down || y != kParams.mouse_y_down) {
    kParams.moved_since_last_down = true;
  }
  const uint32_t img_x = ToImageX(x), img_y = ToImageY(y);
  const uint32_t img_x_down = ToImageX(kParams.mouse_x_down);
  const uint32_t img_y_down = ToImageY(kParams.mouse_y_down);

  WP2::Rectangle* const rect = &kParams.forcing_block;
  if (rect->width > 0 && rect->height > 0) {
    // Maximum block position and size (both relative to top-left).
    const uint32_t max_w =
        WP2::SizeBlocks(kParams.width) * WP2::kMinBlockSizePix;
    const uint32_t max_h =
        WP2::SizeBlocks(kParams.height) * WP2::kMinBlockSizePix;
    const uint32_t max_x = WP2::SafeSub(max_w, WP2::kMinBlockSizePix);
    const uint32_t max_y = WP2::SafeSub(max_h, WP2::kMinBlockSizePix);

    rect->x = std::min(AlignWithBlocks(std::min(img_x_down, img_x)), max_x);
    rect->y = std::min(AlignWithBlocks(std::min(img_y_down, img_y)), max_y);
    rect->width = std::max(img_x_down, img_x) - rect->x + 1;
    rect->height = std::max(img_y_down, img_y) - rect->y + 1;

    // Find the valid block matching the drawn rectangle the most.
    WP2::BlockSize closest_size = WP2::BLK_4x4;
    const WP2::PartitionSet ps = kParams.enc_config.partition_set;
    for (uint32_t i = 0; i < WP2::GetNumBlockSizes(ps); ++i) {
      const WP2::BlockSize size = WP2::GetBlockSizes(ps)[i];
      if (rect->x + WP2::BlockWidthPix(size) <= max_w &&
          rect->y + WP2::BlockHeightPix(size) <= max_h &&
          rect->x + rect->width >= WP2::BlockWidthPix(size) &&
          rect->y + rect->height >= WP2::BlockHeightPix(size) &&
          GetDist(*rect, size) < GetDist(*rect, closest_size)) {
        closest_size = size;
      }
    }
    rect->width = WP2::BlockWidthPix(closest_size);
    rect->height = WP2::BlockHeightPix(closest_size);

    // Replace the rectangle so that the first click (down) is a corner.
    if (img_x < img_x_down) {
      rect->x = WP2::SafeSub(AlignWithBlocks(img_x_down), rect->width);
    }
    if (img_y < img_y_down) {
      rect->y = WP2::SafeSub(AlignWithBlocks(img_y_down), rect->height);
    }
  }
}

void HandleSegmentClick(int button, int state, uint32_t img_x_down,
                        uint32_t img_y_down, uint32_t img_x, uint32_t img_y) {
  (void)button;
  (void)state;
  (void)img_x_down;
  (void)img_y_down;
  (void)img_x;
  (void)img_y;

#if defined(WP2_BITTRACE)
  if (!kParams.IsSegmentVisible() || button == GLUT_LEFT_BUTTON ||
      state != GLUT_DOWN) {
    return;
  }

  if (!kParams.dinfo.explicit_segment_ids) {
    AddInfo("Cannot force segments at low quality (segment ids are implicit).");
    return;
  }

  uint32_t num_segments = kParams.dinfo.num_segments;
  // No point in forcing the segment if there's only one...
  if (num_segments == 1) return;

  uint32_t block_x = 0, block_y = 0;
  for (const WP2::BlockInfo block : kParams.einfo.blocks) {
    if (block.rect.Contains(img_x_down, img_y_down)) {
      block_x = block.rect.x;
      block_y = block.rect.y;
      break;
    }
  }

  bool found = false;
  for (uint32_t i = 0; i < kParams.force_segment.size(); ++i) {
    WP2::EncoderInfo::ForcedSegment* const forced_segment =
        &kParams.force_segment[i];
    if (forced_segment->x == block_x && forced_segment->y == block_y) {
      found = true;
      forced_segment->segment_id = (forced_segment->segment_id + 1);
      if (button == GLUT_MIDDLE_BUTTON ||
          forced_segment->segment_id == num_segments) {
        // Remove.
        kParams.force_segment[i] = kParams.force_segment.back();
        kParams.force_segment.pop_back();
      }
      break;
    }
  }
  if (button == GLUT_RIGHT_BUTTON && !found) {
    WP2::EncoderInfo::ForcedSegment forced = {block_x, block_y, 0};
    kParams.force_segment.emplace_back(forced);
  }
#endif
}

void HandlePredictorClick(int button, int state, uint32_t img_x_down,
                          uint32_t img_y_down, uint32_t img_x, uint32_t img_y) {
  (void)button;
  (void)state;
  (void)img_x_down;
  (void)img_y_down;
  (void)img_x;
  (void)img_y;

#if defined(WP2_BITTRACE)
  if (!kParams.IsPredictorVisible() || button == GLUT_LEFT_BUTTON ||
      state != GLUT_DOWN)
    return;

  WP2::Channel channel = WP2::VDChannel(kParams.dec_config);
  if (channel == WP2::kVChannel) channel = WP2::kUChannel;

  const WP2::BlockInfo* b = nullptr;
  for (const WP2::BlockInfo& block : kParams.einfo.blocks) {
    if (block.rect.Contains(img_x_down, img_y_down)) {
      b = &block;
      break;
    }
  }
  assert(b != nullptr);

  const uint32_t num_predictors = GetPredictors(*b, channel).size();

  bool found = false;
  for (uint32_t i = 0; i < kParams.force_predictor.size(); ++i) {
    auto* const forced_predictor = &kParams.force_predictor[i];
    if (forced_predictor->x == b->rect.x && forced_predictor->y == b->rect.y &&
        forced_predictor->channel == channel) {
      found = true;
      ++forced_predictor->predictor_id;
      if (button == GLUT_MIDDLE_BUTTON ||
          forced_predictor->predictor_id >= num_predictors) {
        // Remove.
        kParams.force_predictor[i] = kParams.force_predictor.back();
        kParams.force_predictor.pop_back();
      }
      break;
    }
  }
  if (button == GLUT_RIGHT_BUTTON && !found) {
    WP2::EncoderInfo::ForcedPredictor forced = {channel, b->rect.x, b->rect.y,
                                                0};
    kParams.force_predictor.emplace_back(forced);
  }
#endif
}

void HandlePartitionClick(int button, int state, uint32_t img_x_down,
                          uint32_t img_y_down, uint32_t img_x, uint32_t img_y) {
  if (button != GLUT_RIGHT_BUTTON && button != GLUT_MIDDLE_BUTTON) return;

  if (button == GLUT_RIGHT_BUTTON && state == GLUT_DOWN) {
    if (kParams.IsPartitionVisible()) {
      // Start drawing a block.
      kParams.forcing_block.x = AlignWithBlocks(img_x_down);
      kParams.forcing_block.y = AlignWithBlocks(img_y_down);
      kParams.forcing_block.width = WP2::kMinBlockSizePix;
      kParams.forcing_block.height = WP2::kMinBlockSizePix;
    } else {
      kParams.forcing_block = {0, 0, 0, 0};
    }
    return;
  }

  const bool delete_only = button == GLUT_MIDDLE_BUTTON;
  if ((kParams.forcing_block.width > 0 && kParams.forcing_block.height > 0 &&
       state == GLUT_UP) ||
      delete_only) {
    bool removed_block = false;
    if (!kParams.moved_since_last_down) {
      // Remove the block under the mouse.
      for (uint32_t i = 0; i < kParams.force_partition.size(); ++i) {
        if (kParams.force_partition[i].Contains(img_x, img_y)) {
          kParams.force_partition[i] = kParams.force_partition.back();
          kParams.force_partition.pop_back();
          removed_block = true;
          break;
        }
      }
    }
    if (!removed_block && !delete_only) {
      // Remove all blocks intersecting the new one.
      for (uint32_t i = 0; i < kParams.force_partition.size();) {
        if (kParams.force_partition[i].Intersects(kParams.forcing_block)) {
          kParams.force_partition[i] = kParams.force_partition.back();
          kParams.force_partition.pop_back();
        } else {
          ++i;
        }
      }
      kParams.force_partition.emplace_back(kParams.forcing_block);
    }
  }
  kParams.forcing_block = {0, 0, 0, 0};
}

void HandleMenuClick(int button, int state, uint32_t img_x_down,
                     uint32_t img_y_down, uint32_t img_x, uint32_t img_y) {
  if (button != GLUT_LEFT_BUTTON) return;

  if (state == GLUT_DOWN) {
    kParams.einfo.selection = kParams.dinfo.selection = {img_x, img_y, 1, 1};
    kParams.forcing_block = {0, 0, 0, 0};
  } else if (state == GLUT_UP) {
    if (kParams.show == Params::kMenu) {
      kParams.show = Params::kDebug;
      if (!kParams.menu_selection.empty()) {
        // Select the first available leaf (DFS) if not already the case.
        SetMenuSelection(GetLeaf(kParams.menu_selection));
      }
      kParams.einfo.selection = kParams.dinfo.selection = {0, 0, 0, 0};
      kParams.forcing_block = {0, 0, 0, 0};
      if (WP2::VDMatch(kParams.enc_config, "encoder")) {
        if (!kParams.EncodeAndDecode()) return;
      } else {
        if (!kParams.DecodeOutput()) return;
      }
    } else {
      const uint32_t min_x = std::min(img_x_down, img_x);
      const uint32_t min_y = std::min(img_y_down, img_y);
      kParams.einfo.selection = kParams.dinfo.selection = {
          min_x, min_y, std::max(img_x_down, img_x) - min_x + 1,
          std::max(img_y_down, img_y) - min_y + 1};
      kParams.forcing_block = {0, 0, 0, 0};

      // Refresh vdebugs that might handle mouse 'selection'.
      if (WP2::VDMatch(kParams.enc_config, "encoder")) {
        if (!kParams.EncodeAndDecode()) return;
      } else if (WP2::VDMatch(kParams.dec_config, "")) {
        if (!kParams.DecodeOutput()) return;
      }
    }
  }
}

void HandleMouseClick(int button, int state, int x, int y) {
  kParams.mouse_x = x;
  kParams.mouse_y = y;
  if (state == GLUT_DOWN) {
    kParams.mouse_x_down = x;
    kParams.mouse_y_down = y;
    kParams.moved_since_last_down = false;
  }
  const uint32_t img_x = ToImageX(x), img_y = ToImageY(y);
  const uint32_t img_x_down = ToImageX(kParams.mouse_x_down);
  const uint32_t img_y_down = ToImageY(kParams.mouse_y_down);

  HandleMenuClick(button, state, img_x_down, img_y_down, img_x, img_y);
  HandleSegmentClick(button, state, img_x_down, img_y_down, img_x, img_y);
  HandlePredictorClick(button, state, img_x_down, img_y_down, img_x, img_y);
  HandlePartitionClick(button, state, img_x_down, img_y_down, img_x, img_y);
}

void SetQuality(float incr) {
  if (glutGetModifiers() & GLUT_ACTIVE_ALT) incr *= 100.;  // Min or max
  const bool set_alpha =
      kParams.in.HasTransparency() && glutGetModifiers() & GLUT_ACTIVE_CTRL;
  float* const quality = set_alpha ? &kParams.alpha_quality : &kParams.quality;
  const float new_quality = WP2::Clamp(*quality + incr, 0.f, 100.f);
  if (*quality == new_quality) return;
  *quality = new_quality;
  if (!kParams.EncodeAndDecode()) return;
  SetMessage(
      SPrintf("%s: %.1f", (set_alpha ? "Alpha quality" : "Quality"), *quality));
}

void HandleSpecialKeys(int key, int pos_x, int pos_y) {
  (void)pos_x;
  (void)pos_y;
  switch (key) {
    default: return;
    case GLUT_KEY_UP: SetQuality(1); break;
    case GLUT_KEY_DOWN: SetQuality(-1); break;
    case GLUT_KEY_RIGHT: SetQuality(10); break;
    case GLUT_KEY_LEFT: SetQuality(-10); break;
  }
}

void HandleReshape(int width, int height) {
  glViewport(0, 0, width, height);
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();
  kParams.viewport_width = (uint32_t)width;
  kParams.viewport_height = (uint32_t)height;
}

//------------------------------------------------------------------------------
// Display

void DisplayCheckerboard(uint32_t width, uint32_t height) {
  static constexpr uint32_t kTileSize = 8;

  glMatrixMode(GL_PROJECTION);
  glPushMatrix();
  glLoadIdentity();
  glOrtho(0, width, height, 0, -1, 1);  // Vertical flip to start at top-left.

  for (uint32_t top = 0; top < height; top += kTileSize) {
    for (uint32_t left = 0; left < width; left += kTileSize) {
      const GLubyte color = !((top + left) & kTileSize) ? 255 : 200;
      glColor3ub(color, color, color);
      // We do not care about coordinates outside viewport, don't test it.
      glRecti(left, top, left + kTileSize, top + kTileSize);
    }
  }

  glPopMatrix();
}

void DisplayBackground(uint32_t width, uint32_t height) {
  if (kParams.background == Params::kCheckerboard) {
    DisplayCheckerboard(width, height);
    return;
  }

  switch (kParams.background) {
    case Params::kWhite:
      glColor3f(1.f, 1.f, 1.f);
      break;
    case Params::kBlack:
      glColor3f(0.f, 0.f, 0.f);
      break;
    case Params::kPink:
      glColor3f(1.f, 0.f, 1.f);
      break;
    default:
      assert(false);
  }
  glRectf(-1.f, -1.f, 1.f, 1.f);
}

void DisplayBuffer(const WP2::ArgbBuffer& buffer) {
  if (buffer.IsEmpty()) return;
  if (buffer.HasTransparency()) {
    DisplayBackground(buffer.width, buffer.height);
  }

  // convert from src->Argb Argb format to GL's BGRA format.
  std::vector<uint8_t> rgba(4u * buffer.width * buffer.height);
  for (uint32_t y = 0; y < buffer.height; ++y) {
    const uint8_t* const src = (const uint8_t*)buffer.GetRow(y);
    uint8_t* const dst = &rgba[y * 4u * buffer.width];
    WP2ArgbConvertTo[WP2_RGBA_32](src, buffer.width, dst);
  }
  glDrawPixels(buffer.width, buffer.height, GL_RGBA, GL_UNSIGNED_BYTE,
               reinterpret_cast<const GLvoid*>(&rgba[0]));
}

void DisplayColor(WP2::Argb32b color) {
  assert(color.a == 255u);
  glMatrixMode(GL_PROJECTION);
  glPushMatrix();
  glLoadIdentity();
  glOrtho(0, 1, 1, 0, -1, 1);  // Vertical flip to start at top-left.
  glColor3ub(color.r, color.g, color.b);
  glRecti(0, 0, 1, 1);
  glPopMatrix();
}

void HandleDisplay() {
  if (kParams.out.IsEmpty()) return;
  glPushMatrix();
  glPixelZoom((GLfloat)(+1. / kParams.width * kParams.viewport_width),
              (GLfloat)(-1. / kParams.height * kParams.viewport_height));
  glRasterPos2f(-1.f, 1.f);
  glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
  glPixelStorei(GL_UNPACK_ROW_LENGTH, kParams.width);

  if (kParams.show == Params::kPreviewColor) {
    DisplayColor(kParams.preview_color);
  } else {
    DisplayBuffer(kParams.GetBuffer());
  }

  if (kParams.show_interface) PrintInfo();
  glPopMatrix();
  glutSwapBuffers();
}

void StartDisplay(const char* window_title) {
  const uint32_t width = kParams.width;
  const uint32_t height = kParams.height;
  const uint32_t swidth = (uint32_t)glutGet(GLUT_SCREEN_WIDTH);
  const uint32_t sheight = (uint32_t)glutGet(GLUT_SCREEN_HEIGHT);

  (void)swidth;
  (void)sheight;
  glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGBA);
  glutInitWindowSize(width, height);
  glutCreateWindow(window_title);
  glutDisplayFunc(HandleDisplay);
  glutReshapeFunc(HandleReshape);
  glutIdleFunc(nullptr);
  glutKeyboardFunc(HandleKey);
  glutKeyboardUpFunc(HandleKeyUp);
  glutIgnoreKeyRepeat(GL_TRUE);
  glutPassiveMotionFunc(HandleMouseMove);
  glutMotionFunc(HandleMouseMoveButtonPressed);
  glutMouseFunc(HandleMouseClick);
  glutSpecialFunc(HandleSpecialKeys);
  glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
  glEnable(GL_BLEND);
  glClearColor(0., 0., 0., 1.);
  glClear(GL_COLOR_BUFFER_BIT);
  WP2ArgbConverterInit();
  DisplayLoop(40);
}

//------------------------------------------------------------------------------
// Main

void Help() {
  ProgramOptions opt;
  opt.Add("Visualizer for WebP2 (re-)compression, using OpenGL");
  opt.Add("");
  opt.Add("Usage:");
  opt.Add("  vwp2 in_file [options]");
  opt.Add("");
  opt.Add("Options are:");
  opt.Add("-q <float>", SPrintf("quality factor in (0..100) range, "
                                "from strongest compression to lossless "
                                "(%1.1f by default)",
                                WP2::EncoderConfig::kDefault.quality));
  opt.Add("-alpha_q <float>",
          SPrintf("alpha quality factor in (0..100) range, "
                  "from strongest compression to lossless "
                  "(%1.1f by default)",
                  WP2::EncoderConfig::kDefault.alpha_quality));
  opt.Add("-speed <int>", SPrintf("compression effort in (0..9) range, "
                                  "from fastest to slowest/best results "
                                  "(%d by default)",
                                  WP2::EncoderConfig::kDefault.speed));
  opt.Add("-pm <int>",
          SPrintf("partition method (0..%d)", WP2::NUM_PARTITION_METHODS - 1));
  opt.Add("-ps <int>",
          SPrintf("partition set (0..%d)", WP2::NUM_PARTITION_SETS - 1));
  opt.Add("-snap", "force quadtree-like partitioning");
  opt.Add("-sns <float>", "segmentation strength (0=off..100)");
  opt.Add("-uv_mode <int>", "UV subsampling mode (0=Adapt,1=420,2=444,3=Auto)");
  opt.Add("-csp <int>", "color space");
  opt.Add("-segments <int>", "number of segments ([1..8])");
  opt.Add("-segment_mode [mode]", "one of auto, explicit, implicit");
  opt.Add("-pass <int>", "number of entropy-analysis passes (1..10)");
  opt.Add("-[no_]dblk_filter", "enable/disable deblocking filter.");
  opt.Add("-[no_]drct_filter", "enable/disable directional filter.");
  opt.Add("-crop <x> <y> <w> <h>", "crop picture with the given rectangle");
  opt.Add("-set_block <x> <y> <w> <h>", "force one block (in px)");
  opt.Add("-set_partition <file>", "force blocks partition (x y w h per line)");
  opt.Add("-diffusion <int>", "error diffusion strength (0=off..100)");
  opt.Add("-set_dump_path <file>", "path to the dump file (for 'l'/'L' keys)");
  opt.Add("-quants <floats,...>", "force per-segment quantization factors");
  opt.Add("-nomt", "disable multi-threading");
  opt.Add("-[no]metadata", "include [or ignore] source metadata");
  opt.Add("-d", "decode-only: just show the .wp2 file");
  opt.Add("-alt <file>", "alternative image to compare with");
  opt.Add("-vbt <int>", "number of visual bittrace levels");
  opt.Add("-version", "print version number and exit");
  opt.Add("-info", "print info overlay");
  opt.Add("-title <text>", "set a custom title for the window");
  opt.Add("");
  // System options.
  opt.AddSystemOptionSection();
  opt.Print();
  // Keyboard options.
  const std::vector<std::string> help = GetHelp();
  for (const std::string& str : help) printf("%s\n", str.c_str());
}

}  // namespace

int main(int argc, char* argv[]) {
  EXIT_IF_FALSE(WP2CheckVersion(), "Error! Library version mismatch!");

  bool read_initial_partition = false;
  const char* alt_name = nullptr;
  const char* window_title = "WebP2 viewer";
  for (int c = 1; c < argc; ++c) {
    bool parse_error = false;
    if (!strcmp(argv[c], "-h") || !strcmp(argv[c], "-help")) {
      Help();
      return 0;
    } else if (!strcmp(argv[c], "-info")) {
      kParams.show = Params::kInfo;
    } else if (!strcmp(argv[c], "-d")) {
      kParams.view_only = true;
    } else if (!strcmp(argv[c], "-alt") && c + 1 < argc) {
      alt_name = argv[++c];
    } else if (!strcmp(argv[c], "-pm") && c + 1 < argc) {
      kParams.enc_config.partition_method =
          (WP2::PartitionMethod)ExUtilGetInt(argv[++c], &parse_error);
    } else if (!strcmp(argv[c], "-ps") && c + 1 < argc) {
      kParams.enc_config.partition_set =
          (WP2::PartitionSet)ExUtilGetInt(argv[++c], &parse_error);
    } else if (!strcmp(argv[c], "-snap")) {
      kParams.enc_config.partition_snapping = true;
    } else if (!strcmp(argv[c], "-sns") && c + 1 < argc) {
      kParams.enc_config.sns = ExUtilGetFloat(argv[++c], &parse_error);
    } else if (!strcmp(argv[c], "-uv_mode") && c + 1 < argc) {
      kParams.enc_config.uv_mode =
          (WP2::EncoderConfig::UVMode)ExUtilGetInt(argv[++c], &parse_error);
    } else if (!strcmp(argv[c], "-csp") && c + 1 < argc) {
      kParams.enc_config.csp_type =
          (WP2::Csp)ExUtilGetInt(argv[++c], &parse_error);
    } else if (!strcmp(argv[c], "-speed") && c + 1 < argc) {
      kParams.enc_config.speed = ExUtilGetInt(argv[++c], &parse_error);
    } else if (!strcmp(argv[c], "-segments") && c + 1 < argc) {
      kParams.enc_config.segments = ExUtilGetInt(argv[++c], &parse_error);
    } else if (!strcmp(argv[c], "-segment_mode") && c + 1 < argc) {
      ++c;
      if (!strcmp(argv[c], "auto")) {
        kParams.enc_config.segment_id_mode =
            WP2::EncoderConfig::SEGMENT_ID_AUTO;
      } else if (!strcmp(argv[c], "explicit")) {
        kParams.enc_config.segment_id_mode =
            WP2::EncoderConfig::SEGMENT_ID_EXPLICIT;
      } else if (!strcmp(argv[c], "implicit")) {
        kParams.enc_config.segment_id_mode =
            WP2::EncoderConfig::SEGMENT_ID_IMPLICIT;
      } else {
        EXIT_IF_FALSE(false, "unsupported segment mode '%s'", argv[c]);
      }
    } else if (!strcmp(argv[c], "-pass") && c + 1 < argc) {
      kParams.enc_config.pass = ExUtilGetInt(argv[++c], &parse_error);
    } else if (!strcmp(argv[c], "-nomt")) {
      kParams.enc_config.thread_level = 0;
      kParams.dec_config.thread_level = 0;
    } else if (!strcmp(argv[c], "-metadata") ||
               !strcmp(argv[c], "-nometadata")) {
      kParams.keep_metadata = !strcmp(argv[c], "-metadata");
    } else if (!strcmp(argv[c], "-bt") || !strcmp(argv[c], "-BT")) {
      kParams.bit_trace = !strcmp(argv[c], "-bt") ? 1 : 2;
      if (c + 1 < argc) {
        if (isdigit(argv[c + 1][0])) {
          kParams.bit_trace_level = ExUtilGetUInt(argv[++c], &parse_error);
        } else {
          kParams.bit_trace_level = 0;
        }
      }
    } else if (!strcmp(argv[c], "-vbt") && c + 1 < argc) {
      kParams.visual_bit_trace_level = ExUtilGetUInt(argv[++c], &parse_error);
#if defined(WP2_BITTRACE)
      if (kParams.visual_bit_trace_level > kBitTraceMaxLevels) {
        printf("-vbt must be at most '%d'\n", kBitTraceMaxLevels);
        parse_error = true;
      }
#else
      printf("WP2_BITTRACE required for -vbt\n");
      parse_error = true;
#endif  // WP2_BITTRACE
    } else if (!strcmp(argv[c], "--") && c + 1 < argc) {
      kParams.files.emplace_back(argv[++c]);
    } else if (!strcmp(argv[c], "-q") && c + 1 < argc) {
      kParams.quality = ExUtilGetFloat(argv[++c], &parse_error);
    } else if (!strcmp(argv[c], "-quants") && c + 1 < argc) {
      ExUtilGetFloats(argv[++c],
                      kParams.enc_config.segment_factors, WP2::kMaxNumSegments,
                      &parse_error);
    } else if (!strcmp(argv[c], "-alpha_q") && c + 1 < argc) {
      kParams.alpha_quality = ExUtilGetFloat(argv[++c], &parse_error);
    } else if (!strcmp(argv[c], "-dblk_filter") ||
               !strcmp(argv[c], "-no_dblk_filter")) {
      kParams.dec_config.enable_deblocking_filter =
          !strcmp(argv[c], "-dblk_filter");
    } else if (!strcmp(argv[c], "-drct_filter") ||
               !strcmp(argv[c], "-no_drct_filter")) {
      kParams.dec_config.enable_directional_filter =
          !strcmp(argv[c], "-drct_filter");
    } else if (!strcmp(argv[c], "-rstr_filter") ||
               !strcmp(argv[c], "-no_rstr_filter")) {
      kParams.dec_config.enable_restoration_filter =
          !strcmp(argv[c], "-rstr_filter");
    } else if (!strcmp(argv[c], "-title") && c + 1 < argc) {
      window_title = argv[++c];
    } else if (!strcmp(argv[c], "-crop") && c + 4 < argc) {
      kParams.crop = true;
      kParams.crop_area.x = ExUtilGetUInt(argv[++c], &parse_error);
      kParams.crop_area.y = ExUtilGetUInt(argv[++c], &parse_error);
      kParams.crop_area.width = ExUtilGetUInt(argv[++c], &parse_error);
      kParams.crop_area.height = ExUtilGetUInt(argv[++c], &parse_error);
    } else if (!strcmp(argv[c], "-set_block") && c + 4 < argc) {
      kParams.force_partition.emplace_back(
          ExUtilGetUInt(argv[c + 1], &parse_error),
          ExUtilGetUInt(argv[c + 2], &parse_error),
          ExUtilGetUInt(argv[c + 3], &parse_error),
          ExUtilGetUInt(argv[c + 4], &parse_error));
      c += 4;
    } else if (!strcmp(argv[c], "-set_partition") && c + 1 < argc) {
      kParams.partition_file_path = argv[++c];
      read_initial_partition = true;
    } else if (!strcmp(argv[c], "-diffusion") && c + 1 < argc) {
      kParams.enc_config.error_diffusion =
          ExUtilGetInt(argv[++c], &parse_error);
    } else if (!strcmp(argv[c], "-set_dump_path") && c + 1 < argc) {
      const std::string extension = WP2GetFileExtension(argv[++c]);
      if (extension == "wp2") {
        kParams.dump_wp2_path = argv[c];
      } else if (extension == "png") {
        kParams.dump_png_path = argv[c];
      } else {
        printf("Unhandled extension for -set_dump_path '%s'\n", argv[c]);
        parse_error = true;
      }
    } else if (argv[c][0] == '-') {
      bool must_stop;
      if (ProgramOptions::ParseSystemOptions(argv[c], &must_stop)) {
        if (must_stop) return 0;
      } else {
        printf("Unknown option '%s'\n", argv[c]);
        parse_error = true;
      }
    } else {
      kParams.files.emplace_back(argv[c]);
    }

    if (parse_error) {
      Help();
      return 1;
    }
  }

  if (kParams.files.empty()) {
    printf("missing input file(s)!!\n");
    Help();
    return 0;
  }

  if (!kParams.SetCurrentFile(0)) return 1;

  if (read_initial_partition) {
    EXIT_IF_FALSE(kParams.ReadPartition(/*read_only=*/false),
                  "Error: Unable to parse partition file");
  }

  if (alt_name != nullptr) {
    EXIT_IF_FALSE(kParams.SetAltFile(alt_name),
                  "Error: Unable to read alt file");
  }

#if defined(__unix__) || defined(__CYGWIN__)
  // Work around GLUT compositor bug.
  // https://bugs.launchpad.net/ubuntu/+source/freeglut/+bug/369891
  setenv("XLIB_SKIP_ARGB_VISUALS", "1", 1);
#endif

  // Start display (and timer)
  glutInit(&argc, argv);
#if defined(WP2_HAVE_GLUT) && defined(FREEGLUT)
  glutSetOption(GLUT_ACTION_ON_WINDOW_CLOSE, GLUT_ACTION_CONTINUE_EXECUTION);
#endif
  StartDisplay(window_title);

  glutMainLoop();

  // Should only be reached when using FREEGLUT:
  return 0;
}

#else  // !WP2_HAVE_OPENGL

int main(int argc, const char *argv[]) {
  fprintf(stderr, "OpenGL support not enabled in %s.\n", argv[0]);
  (void)argc;
  return 0;
}

#endif

//------------------------------------------------------------------------------
