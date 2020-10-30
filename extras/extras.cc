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
//  extra helper functions
//
// Author: Skal (pascal.massimino@gmail.com)

#include <cstdint>
#include <cstring>  // for memcpy()
#include <sstream>

#include "./extras.h"
#include "src/utils/utils.h"
#include "src/common/preview/preview.h"

//------------------------------------------------------------------------------

int WP2GetExtrasVersion() { return WP2_EXTRAS_ABI_VERSION; }

namespace WP2 {

//------------------------------------------------------------------------------

WP2Status ArgbBufferRescale(ArgbBuffer* buffer,
                            uint32_t width, uint32_t height) {
  ArgbBuffer rescaled_buffer;
  WP2_CHECK_STATUS(rescaled_buffer.Resize(width, height));

  // Simple point-sampling for now.
  uint8_t* dst = (uint8_t*)rescaled_buffer.GetRow(0);
  const uint8_t* const src = (const uint8_t*)buffer->GetRow(0);
  for (uint32_t j = 0; j < height; ++j) {
    const uint32_t J = j * buffer->height / height;
    for (uint32_t i = 0; i < width; ++i) {
      const uint32_t I = i * buffer->width / width;
      memcpy(&dst[4 * i], &src[J * buffer->stride + 4 * I], sizeof(uint32_t));
    }
    dst += rescaled_buffer.stride;
  }
  WP2_CHECK_STATUS(buffer->CopyFrom(rescaled_buffer));
  return WP2_STATUS_OK;
}

//------------------------------------------------------------------------------

namespace {

constexpr uint32_t kMaxLen = kPreviewMaxGridSize + 1;

constexpr uint8_t kBase64Chars[] =
    "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/";

// Inverse map: for (uint8_t i = 0; i < 64; ++i) kMap[kBase64Chars[i]] = i;
constexpr uint8_t kMap[] = {
  0x3e, 0xff, 0xff, 0xff, 0x3f, 0x34, 0x35, 0x36, 0x37, 0x38, 0x39, 0x3a, 0x3b,
  0x3c, 0x3d, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0x00, 0x01, 0x02, 0x03,
  0x04, 0x05, 0x06, 0x07, 0x08, 0x09, 0x0a, 0x0b, 0x0c, 0x0d, 0x0e, 0x0f, 0x10,
  0x11, 0x12, 0x13, 0x14, 0x15, 0x16, 0x17, 0x18, 0x19, 0xff, 0xff, 0xff, 0xff,
  0xff, 0xff, 0x1a, 0x1b, 0x1c, 0x1d, 0x1e, 0x1f, 0x20, 0x21, 0x22, 0x23, 0x24,
  0x25, 0x26, 0x27, 0x28, 0x29, 0x2a, 0x2b, 0x2c, 0x2d, 0x2e, 0x2f, 0x30, 0x31,
  0x32, 0x33
};
static_assert(sizeof(kMap) == 'z' - '+' + 1, "missing mapping elements.");

// Converts 3 bytes to 4 printable ASCII characters.
void Base64_3To4(const uint8_t in3[3], char out4[4 + 1]) {
  out4[0] = kBase64Chars[((in3[0] & 0xfc) >> 2)];
  out4[1] = kBase64Chars[((in3[0] & 0x03) << 4) + ((in3[1] & 0xf0) >> 4)];
  out4[2] = kBase64Chars[((in3[1] & 0x0f) << 2) + ((in3[2] & 0xc0) >> 6)];
  out4[3] = kBase64Chars[ (in3[2] & 0x3f)];
  out4[4] = '\0';
}

constexpr bool VerifyMapping1(size_t i = 0) {
  return (i >= sizeof(kBase64Chars) - 1) ||
         (kMap[(int)kBase64Chars[i] - '+'] == i && VerifyMapping1(i + 1));
}
static_assert(VerifyMapping1(), "Incoherent direct mapping!");

constexpr bool VerifyMapping2(size_t i = 0) {
  return (i >= sizeof(kMap)) ||
         ((kMap[i] == 0xff || kBase64Chars[kMap[i]] == '+' + i) &&
           VerifyMapping2(i + 1));
}
static_assert(VerifyMapping2(), "Incoherent reverse mapping!");

}  // namespace

std::string PreviewToBase64(const PreviewData& data) {
  ANSEnc enc;
  if (data.Encode(&enc) != WP2_STATUS_OK) return "";

  std::string out;
  const uint32_t size = enc.BufferSize();
  const uint8_t* const buf = enc.Buffer();
  // Base64 converts 3 bytes to 4 ASCII characters, hence the 4/3 size expansion
  out.reserve(DivCeil(size * 4, 3));
  uint8_t in3[3];
  char out4[4 + 1];
  int i = 0;
  for (uint32_t p = 0; p < size; ++p) {
    in3[i++] = buf[p];
    if (i == 3) {
      Base64_3To4(in3, out4);
      out += out4;
      i = 0;
    }
  }
  if (i > 0) {  // i = 1 or 2
    for (int j = i; j < 3; ++j) in3[j] = '\0';
    Base64_3To4(in3, out4);
    out += out4[0];
    out += out4[1];
    out += (i == 2) ? out4[2] : '=';
    out += '=';
  }
  return out;
}

WP2Status PreviewFromBase64(const std::string& in, PreviewData* preview) {
  WP2_CHECK_OK(preview != nullptr, WP2_STATUS_INVALID_PARAMETER);
  preview->Reset();

  Vector_u8 out;
  WP2_CHECK_ALLOC_OK(out.resize(in.size()));

  uint32_t size = 0;
  uint32_t val = 0;
  int nb_bits = -8;
  for (uint8_t c : in) {
    if (c < '+' || c > 'z') break;
    c -= '+';
    if (kMap[c] == 0xff) break;
    val = (val << 6) + kMap[c];
    nb_bits += 6;
    if (nb_bits >= 0) {
      WP2_CHECK_OK(size < out.size(), WP2_STATUS_BAD_READ);
      out[size++] = (uint8_t)(val >> nb_bits);
      nb_bits -= 8;
    }
  }
  return preview->Decode(out.data(), size);
}

std::string PreviewToText(const PreviewData& data) {
  char tmp[kMaxLen];
  std::string out;
  snprintf(tmp, kMaxLen, "%d %d %d %d\n",
           data.grid_width_, data.grid_height_,
           (int)data.palette_.size(), data.use_noise_ ? 1 : 0);
  out += tmp;
  for (auto& c : data.palette_) {
    snprintf(tmp, kMaxLen, "%d %d %d %d\n", c.y, c.co, c.cg, c.a);
    out += tmp;
  }

  for (uint32_t k = 0, y = 0; y < data.grid_height_; ++y) {
    for (uint32_t x = 0; x < data.grid_width_; ++x) {
      char c = '.';
      if (data.IsCorner(x, y)) {
        c = 'a' + data.corners_[(x != 0) + 2 * (y != 0)];
      } else if (k < data.vertices_.size() &&
                 data.vertices_[k].y == y && data.vertices_[k].x == x) {
        c = 'a' + data.vertices_[k].color_index;
        ++k;
      }
      out += c;
    }
    out += "\n";
  }
  return out;
}

WP2Status PreviewFromText(const std::string& in, PreviewData* preview) {
  WP2_CHECK_OK(preview != nullptr, WP2_STATUS_INVALID_PARAMETER);
  preview->Reset();

  std::istringstream s(in);
  uint32_t nb_colors, use_noise;
  s >> preview->grid_width_ >> preview->grid_height_ >> nb_colors >> use_noise;
  preview->use_noise_ = !!use_noise;
  WP2_CHECK_OK(preview->grid_width_ >= kPreviewMinGridSize &&
               preview->grid_width_ <= kPreviewMaxGridSize,
               WP2_STATUS_BAD_READ);
  WP2_CHECK_OK(preview->grid_height_ >= kPreviewMinGridSize &&
               preview->grid_height_ <= kPreviewMaxGridSize,
               WP2_STATUS_BAD_READ);
  WP2_CHECK_OK(nb_colors >= kPreviewMinNumColors &&
               nb_colors <= kPreviewMaxNumColors,
               WP2_STATUS_BAD_READ);
  WP2_CHECK_OK(!s.eof(), WP2_STATUS_BAD_READ);

  bool has_alpha = false;
  WP2_CHECK_ALLOC_OK(preview->palette_.resize(nb_colors));
  for (auto& c : preview->palette_) {
    WP2_CHECK_OK(!s.eof(), WP2_STATUS_BAD_READ);
    int y, co, cg, a;
    s >> y >> co >> cg >> a;
    c.y  = Clamp(y,  0, 63);
    c.co = Clamp(co, 0, 63);
    c.cg = Clamp(cg, 0, 63);
    c.a = (a != 0) ? 1 : 0;
    has_alpha |= (c.a == 0);
  }
  const uint32_t min_num_vertices =
      std::max(kPreviewMinNumVertices, std::max(nb_colors, 4u) - 4u);
  const uint32_t max_num_vertices =
      std::min(kPreviewMaxNumVertices,
               preview->grid_width_ * preview->grid_height_ - 4);
  WP2_CHECK_ALLOC_OK(preview->vertices_.reserve(max_num_vertices));
  for (uint16_t y = 0; y < preview->grid_height_; ++y) {
    for (uint16_t x = 0; x < preview->grid_width_; ++x) {
      if (!s.eof()) {  // otherwise: assume the rest of the grid is empty.
        char c;
        s >> c;
        const int idx = c - 'a';
        if (preview->IsCorner(x, y)) {
          WP2_CHECK_OK(idx >= 0 && idx < (int)nb_colors, WP2_STATUS_BAD_READ);
          preview->corners_[(x != 0) + 2 * (y != 0)] = idx;
        } else if (c != '.') {
          WP2_CHECK_OK(idx >= 0 && idx < (int)nb_colors, WP2_STATUS_BAD_READ);
          const VertexIndexedColor v = { x, y, (uint16_t)idx };
          WP2_CHECK_ALLOC_OK(preview->vertices_.push_back(v));
        }
      }
    }
  }
  WP2_CHECK_OK(preview->vertices_.size() >= min_num_vertices &&
               preview->vertices_.size() <= max_num_vertices,
               WP2_STATUS_BAD_READ);
  return WP2_STATUS_OK;
}

std::string PrintPreview(const PreviewData& data, bool reduced) {
  std::string out;
  char tmp[kMaxLen];
  snprintf(tmp, kMaxLen,
           "Grid: %d x %d, %zu vertices, %zu colors, use_noise:%d  ",
           data.grid_width_, data.grid_height_,
           data.vertices_.size(),
           data.palette_.size(), data.use_noise_);
  out += tmp;
  snprintf(tmp, kMaxLen, "Corners: [%d %d %d %d]\n",
           data.corners_[0], data.corners_[1],
           data.corners_[2], data.corners_[3]);
  out += tmp;
  if (!reduced) {
    for (uint32_t i = 0; i < data.palette_.size(); ++i) {
      const auto& col = data.palette_[i];
      snprintf(tmp, kMaxLen, "color #%d: y=%d co=%d cg=%d a=%d\n",
               i, col.y, col.co, col.cg, col.a);
      out += tmp;
    }
    out += "Vertices:\n";
    for (uint32_t i = 0; i < data.vertices_.size(); ++i) {
      const auto& v = data.vertices_[i];
      snprintf(tmp, kMaxLen, "#%d: pos=(%d,%d)   idx=%d\n",
               i, v.x, v.y, v.color_index);
      out += tmp;
    }
  }
  return out;
}

//------------------------------------------------------------------------------

}    // namespace WP2
