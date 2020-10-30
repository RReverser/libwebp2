// Copyright 2020 Google LLC
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
// Y4M decoding adapted from
//   chromium.googlesource.com/webm/libvpx/+/refs/heads/master/
//
// Author: Yannis Guyon (yguyon@google.com)

#include <cassert>
#include <cstdio>
#include <cstring>

#include "extras/ccsp_imageio.h"
#include "src/utils/data_source.h"
#include "src/utils/utils.h"

namespace WP2 {
namespace {

// -----------------------------------------------------------------------------
// Some constants and structs from vpx_image.h

#define VPX_IMG_FMT_PLANAR 0x100       /**< Image is a planar format. */
#define VPX_IMG_FMT_UV_FLIP 0x200      /**< V plane precedes U in memory. */
#define VPX_IMG_FMT_HIGHBITDEPTH 0x800 /**< Image uses 16bit framebuffer. */

typedef enum vpx_img_fmt {
  VPX_IMG_FMT_NONE,
  VPX_IMG_FMT_YV12 =
      VPX_IMG_FMT_PLANAR | VPX_IMG_FMT_UV_FLIP | 1, /**< planar YVU */
  VPX_IMG_FMT_I420 = VPX_IMG_FMT_PLANAR | 2,
  VPX_IMG_FMT_I422 = VPX_IMG_FMT_PLANAR | 5,
  VPX_IMG_FMT_I444 = VPX_IMG_FMT_PLANAR | 6,
  VPX_IMG_FMT_I440 = VPX_IMG_FMT_PLANAR | 7,
  VPX_IMG_FMT_I42016 = VPX_IMG_FMT_I420 | VPX_IMG_FMT_HIGHBITDEPTH,
  VPX_IMG_FMT_I42216 = VPX_IMG_FMT_I422 | VPX_IMG_FMT_HIGHBITDEPTH,
  VPX_IMG_FMT_I44416 = VPX_IMG_FMT_I444 | VPX_IMG_FMT_HIGHBITDEPTH,
  VPX_IMG_FMT_I44016 = VPX_IMG_FMT_I440 | VPX_IMG_FMT_HIGHBITDEPTH
} vpx_img_fmt_t; /**< alias for enum vpx_img_fmt */

/*!\brief List of supported color spaces */
typedef enum vpx_color_space {
  VPX_CS_UNKNOWN = 0,   /**< Unknown */
  VPX_CS_BT_601 = 1,    /**< BT.601 */
  VPX_CS_BT_709 = 2,    /**< BT.709 */
  VPX_CS_SMPTE_170 = 3, /**< SMPTE.170 */
  VPX_CS_SMPTE_240 = 4, /**< SMPTE.240 */
  VPX_CS_BT_2020 = 5,   /**< BT.2020 */
  VPX_CS_RESERVED = 6,  /**< Reserved */
  VPX_CS_SRGB = 7       /**< sRGB */
} vpx_color_space_t;    /**< alias for enum vpx_color_space */

/*!\brief List of supported color range */
typedef enum vpx_color_range {
  VPX_CR_STUDIO_RANGE = 0, /**< Y [16..235], UV [16..240] */
  VPX_CR_FULL_RANGE = 1    /**< YUV/RGB [0..255] */
} vpx_color_range_t;       /**< alias for enum vpx_color_range */

/**\brief Image Descriptor */
typedef struct vpx_image {
  vpx_img_fmt_t fmt;       /**< Image Format */
  vpx_color_space_t cs;    /**< Color Space */
  vpx_color_range_t range; /**< Color Range */

  /* Image storage dimensions */
  unsigned int w;         /**< Stored image width */
  unsigned int h;         /**< Stored image height */
  unsigned int bit_depth; /**< Stored image bit-depth */

  /* Image display dimensions */
  unsigned int d_w; /**< Displayed image width */
  unsigned int d_h; /**< Displayed image height */

  /* Image intended rendering dimensions */
  unsigned int r_w; /**< Intended rendering image width */
  unsigned int r_h; /**< Intended rendering image height */

  /* Chroma subsampling info */
  unsigned int x_chroma_shift; /**< subsampling order, X */
  unsigned int y_chroma_shift; /**< subsampling order, Y */

/* Image data pointers. */
#define VPX_PLANE_PACKED 0  /**< To be used for all packed formats */
#define VPX_PLANE_Y 0       /**< Y (Luminance) plane */
#define VPX_PLANE_U 1       /**< U (Chroma) plane */
#define VPX_PLANE_V 2       /**< V (Chroma) plane */
#define VPX_PLANE_ALPHA 3   /**< A (Transparency) plane */
  unsigned char *planes[4]; /**< pointer to the top left pixel for each plane */
  int stride[4];            /**< stride between rows for each plane */

  int bps; /**< bits per sample (for packed formats) */

  /*!\brief The following member may be set by the application to associate
   * data with this image.
   */
  void *user_priv;

  /* The following members should be treated as private. */
  unsigned char *img_data; /**< private */
  int img_data_owner;      /**< private */
  int self_allocd;         /**< private */

  void *fb_priv; /**< Frame buffer data associated with the image. */
} vpx_image_t;   /**< alias for struct vpx_image */

// -----------------------------------------------------------------------------
// Copied as is from y4minput.h

typedef struct y4m_input y4m_input;

/*The function used to perform chroma conversion.*/
typedef void (*y4m_convert_func)(y4m_input *_y4m, unsigned char *_dst,
                                 unsigned char *_src);

struct y4m_input {
  int pic_w;
  int pic_h;
  int fps_n;
  int fps_d;
  int par_n;
  int par_d;
  char interlace;
  int src_c_dec_h;
  int src_c_dec_v;
  int dst_c_dec_h;
  int dst_c_dec_v;
  char chroma_type[16];
  /*The size of each converted frame buffer.*/
  size_t dst_buf_sz;
  /*The amount to read directly into the converted frame buffer.*/
  size_t dst_buf_read_sz;
  /*The size of the auxilliary buffer.*/
  size_t aux_buf_sz;
  /*The amount to read into the auxilliary buffer.*/
  size_t aux_buf_read_sz;
  y4m_convert_func convert;
  unsigned char *dst_buf;
  unsigned char *aux_buf;
  enum vpx_img_fmt vpx_fmt;
  int bps;
  unsigned int bit_depth;
};

// -----------------------------------------------------------------------------
// Copied from y4minput.c
// - Replaced FILE by DataSource and adapted file_read().
// - Testing the LogLevel ('log' arg) before each fprintf().
// - Removed "static" keywords because this an anonymous namespace already.
// - Added "// NOLINT" to strcpy() to remove "use abseil" warning.
// - Testing width/height in y4m_input_open() to avoid integer overflows.

// Reads 'size' bytes from 'file' into 'buf'. Returns true on success.
int file_read(void *buf, size_t size, DataSource *file) {
  const uint8_t* data = nullptr;
  if (file->TryReadNext(size, &data)) {
    // TODO(yguyon): Copying data for API convenience. Adapting file_read()
    //               signature might avoid that.
    std::memcpy(buf, data, size);
    return true;
  }
  return false;
}

int y4m_parse_tags(y4m_input *_y4m, char *_tags) {
  int got_w;
  int got_h;
  int got_fps;
  int got_interlace;
  int got_par;
  int got_chroma;
  char *p;
  char *q;
  got_w = got_h = got_fps = got_interlace = got_par = got_chroma = 0;
  for (p = _tags;; p = q) {
    /*Skip any leading spaces.*/
    while (*p == ' ') p++;
    /*If that's all we have, stop.*/
    if (p[0] == '\0') break;
    /*Find the end of this tag.*/
    for (q = p + 1; *q != '\0' && *q != ' '; q++) {
    }
    /*Process the tag.*/
    switch (p[0]) {
      case 'W': {
        if (sscanf(p + 1, "%d", &_y4m->pic_w) != 1) return -1;
        got_w = 1;
        break;
      }
      case 'H': {
        if (sscanf(p + 1, "%d", &_y4m->pic_h) != 1) return -1;
        got_h = 1;
        break;
      }
      case 'F': {
        if (sscanf(p + 1, "%d:%d", &_y4m->fps_n, &_y4m->fps_d) != 2) {
          return -1;
        }
        got_fps = 1;
        break;
      }
      case 'I': {
        _y4m->interlace = p[1];
        got_interlace = 1;
        break;
      }
      case 'A': {
        if (sscanf(p + 1, "%d:%d", &_y4m->par_n, &_y4m->par_d) != 2) {
          return -1;
        }
        got_par = 1;
        break;
      }
      case 'C': {
        if (q - p > 16) return -1;
        memcpy(_y4m->chroma_type, p + 1, q - p - 1);
        _y4m->chroma_type[q - p - 1] = '\0';
        got_chroma = 1;
        break;
      }
        /*Ignore unknown tags.*/
    }
  }
  if (!got_w || !got_h || !got_fps) return -1;
  if (!got_interlace) _y4m->interlace = '?';
  if (!got_par) _y4m->par_n = _y4m->par_d = 0;
  /*Chroma-type is not specified in older files, e.g., those generated by
     mplayer.*/
  if (!got_chroma) strcpy(_y4m->chroma_type, "420");  // NOLINT
  return 0;
}
/*All anti-aliasing filters in the following conversion functions are based on
   one of two window functions:
  The 6-tap Lanczos window (for down-sampling and shifts):
   sinc(\pi*t)*sinc(\pi*t/3), |t|<3  (sinc(t)==sin(t)/t)
   0,                         |t|>=3
  The 4-tap Mitchell window (for up-sampling):
   7|t|^3-12|t|^2+16/3,             |t|<1
   -(7/3)|x|^3+12|x|^2-20|x|+32/3,  |t|<2
   0,                               |t|>=2
  The number of taps is intentionally kept small to reduce computational
   overhead and limit ringing.
  The taps from these filters are scaled so that their sum is 1, and the
  result is scaled by 128 and rounded to integers to create a filter whose
   intermediate values fit inside 16 bits.
  Coefficients are rounded in such a way as to ensure their sum is still 128,
   which is usually equivalent to normal rounding.
  Conversions which require both horizontal and vertical filtering could
   have these steps pipelined, for less memory consumption and better cache
   performance, but we do them separately for simplicity.*/
#define OC_MINI(_a, _b) ((_a) > (_b) ? (_b) : (_a))
#define OC_MAXI(_a, _b) ((_a) < (_b) ? (_b) : (_a))
#define OC_CLAMPI(_a, _b, _c) (OC_MAXI(_a, OC_MINI(_b, _c)))
/*420jpeg chroma samples are sited like:
  Y-------Y-------Y-------Y-------
  |       |       |       |
  |   BR  |       |   BR  |
  |       |       |       |
  Y-------Y-------Y-------Y-------
  |       |       |       |
  |       |       |       |
  |       |       |       |
  Y-------Y-------Y-------Y-------
  |       |       |       |
  |   BR  |       |   BR  |
  |       |       |       |
  Y-------Y-------Y-------Y-------
  |       |       |       |
  |       |       |       |
  |       |       |       |
  420mpeg2 chroma samples are sited like:
  Y-------Y-------Y-------Y-------
  |       |       |       |
  BR      |       BR      |
  |       |       |       |
  Y-------Y-------Y-------Y-------
  |       |       |       |
  |       |       |       |
  |       |       |       |
  Y-------Y-------Y-------Y-------
  |       |       |       |
  BR      |       BR      |
  |       |       |       |
  Y-------Y-------Y-------Y-------
  |       |       |       |
  |       |       |       |
  |       |       |       |
  We use a resampling filter to shift the site locations one quarter pixel (at
   the chroma plane's resolution) to the right.
  The 4:2:2 modes look exactly the same, except there are twice as many chroma
   lines, and they are vertically co-sited with the luma samples in both the
   mpeg2 and jpeg cases (thus requiring no vertical resampling).*/
void y4m_42xmpeg2_42xjpeg_helper(unsigned char *_dst,
                                        const unsigned char *_src, int _c_w,
                                        int _c_h) {
  int y;
  int x;
  for (y = 0; y < _c_h; y++) {
    /*Filter: [4 -17 114 35 -9 1]/128, derived from a 6-tap Lanczos
       window.*/
    for (x = 0; x < OC_MINI(_c_w, 2); x++) {
      _dst[x] = (unsigned char)OC_CLAMPI(
          0,
          (4 * _src[0] - 17 * _src[OC_MAXI(x - 1, 0)] + 114 * _src[x] +
           35 * _src[OC_MINI(x + 1, _c_w - 1)] -
           9 * _src[OC_MINI(x + 2, _c_w - 1)] + _src[OC_MINI(x + 3, _c_w - 1)] +
           64) >>
              7,
          255);
    }
    for (; x < _c_w - 3; x++) {
      _dst[x] = (unsigned char)OC_CLAMPI(
          0,
          (4 * _src[x - 2] - 17 * _src[x - 1] + 114 * _src[x] +
           35 * _src[x + 1] - 9 * _src[x + 2] + _src[x + 3] + 64) >>
              7,
          255);
    }
    for (; x < _c_w; x++) {
      _dst[x] = (unsigned char)OC_CLAMPI(
          0,
          (4 * _src[x - 2] - 17 * _src[x - 1] + 114 * _src[x] +
           35 * _src[OC_MINI(x + 1, _c_w - 1)] -
           9 * _src[OC_MINI(x + 2, _c_w - 1)] + _src[_c_w - 1] + 64) >>
              7,
          255);
    }
    _dst += _c_w;
    _src += _c_w;
  }
}
/*Handles both 422 and 420mpeg2 to 422jpeg and 420jpeg, respectively.*/
void y4m_convert_42xmpeg2_42xjpeg(y4m_input *_y4m, unsigned char *_dst,
                                         unsigned char *_aux) {
  int c_w;
  int c_h;
  int c_sz;
  int pli;
  /*Skip past the luma data.*/
  _dst += _y4m->pic_w * _y4m->pic_h;
  /*Compute the size of each chroma plane.*/
  c_w = (_y4m->pic_w + _y4m->dst_c_dec_h - 1) / _y4m->dst_c_dec_h;
  c_h = (_y4m->pic_h + _y4m->dst_c_dec_v - 1) / _y4m->dst_c_dec_v;
  c_sz = c_w * c_h;
  for (pli = 1; pli < 3; pli++) {
    y4m_42xmpeg2_42xjpeg_helper(_dst, _aux, c_w, c_h);
    _dst += c_sz;
    _aux += c_sz;
  }
}
/*This format is only used for interlaced content, but is included for
   completeness.
  420jpeg chroma samples are sited like:
  Y-------Y-------Y-------Y-------
  |       |       |       |
  |   BR  |       |   BR  |
  |       |       |       |
  Y-------Y-------Y-------Y-------
  |       |       |       |
  |       |       |       |
  |       |       |       |
  Y-------Y-------Y-------Y-------
  |       |       |       |
  |   BR  |       |   BR  |
  |       |       |       |
  Y-------Y-------Y-------Y-------
  |       |       |       |
  |       |       |       |
  |       |       |       |
  420paldv chroma samples are sited like:
  YR------Y-------YR------Y-------
  |       |       |       |
  |       |       |       |
  |       |       |       |
  YB------Y-------YB------Y-------
  |       |       |       |
  |       |       |       |
  |       |       |       |
  YR------Y-------YR------Y-------
  |       |       |       |
  |       |       |       |
  |       |       |       |
  YB------Y-------YB------Y-------
  |       |       |       |
  |       |       |       |
  |       |       |       |
  We use a resampling filter to shift the site locations one quarter pixel (at
   the chroma plane's resolution) to the right.
  Then we use another filter to move the C_r location down one quarter pixel,
   and the C_b location up one quarter pixel.*/
void y4m_convert_42xpaldv_42xjpeg(y4m_input *_y4m, unsigned char *_dst,
                                         unsigned char *_aux) {
  unsigned char *tmp;
  int c_w;
  int c_h;
  int c_sz;
  int pli;
  int y;
  int x;
  /*Skip past the luma data.*/
  _dst += _y4m->pic_w * _y4m->pic_h;
  /*Compute the size of each chroma plane.*/
  c_w = (_y4m->pic_w + 1) / 2;
  c_h = (_y4m->pic_h + _y4m->dst_c_dec_h - 1) / _y4m->dst_c_dec_h;
  c_sz = c_w * c_h;
  tmp = _aux + 2 * c_sz;
  for (pli = 1; pli < 3; pli++) {
    /*First do the horizontal re-sampling.
      This is the same as the mpeg2 case, except that after the horizontal
       case, we need to apply a second vertical filter.*/
    y4m_42xmpeg2_42xjpeg_helper(tmp, _aux, c_w, c_h);
    _aux += c_sz;
    switch (pli) {
      case 1: {
        /*Slide C_b up a quarter-pel.
          This is the same filter used above, but in the other order.*/
        for (x = 0; x < c_w; x++) {
          for (y = 0; y < OC_MINI(c_h, 3); y++) {
            _dst[y * c_w] = (unsigned char)OC_CLAMPI(
                0,
                (tmp[0] - 9 * tmp[OC_MAXI(y - 2, 0) * c_w] +
                 35 * tmp[OC_MAXI(y - 1, 0) * c_w] + 114 * tmp[y * c_w] -
                 17 * tmp[OC_MINI(y + 1, c_h - 1) * c_w] +
                 4 * tmp[OC_MINI(y + 2, c_h - 1) * c_w] + 64) >>
                    7,
                255);
          }
          for (; y < c_h - 2; y++) {
            _dst[y * c_w] = (unsigned char)OC_CLAMPI(
                0,
                (tmp[(y - 3) * c_w] - 9 * tmp[(y - 2) * c_w] +
                 35 * tmp[(y - 1) * c_w] + 114 * tmp[y * c_w] -
                 17 * tmp[(y + 1) * c_w] + 4 * tmp[(y + 2) * c_w] + 64) >>
                    7,
                255);
          }
          for (; y < c_h; y++) {
            _dst[y * c_w] = (unsigned char)OC_CLAMPI(
                0,
                (tmp[(y - 3) * c_w] - 9 * tmp[(y - 2) * c_w] +
                 35 * tmp[(y - 1) * c_w] + 114 * tmp[y * c_w] -
                 17 * tmp[OC_MINI(y + 1, c_h - 1) * c_w] +
                 4 * tmp[(c_h - 1) * c_w] + 64) >>
                    7,
                255);
          }
          _dst++;
          tmp++;
        }
        _dst += c_sz - c_w;
        tmp -= c_w;
        break;
      }
      case 2: {
        /*Slide C_r down a quarter-pel.
          This is the same as the horizontal filter.*/
        for (x = 0; x < c_w; x++) {
          for (y = 0; y < OC_MINI(c_h, 2); y++) {
            _dst[y * c_w] = (unsigned char)OC_CLAMPI(
                0,
                (4 * tmp[0] - 17 * tmp[OC_MAXI(y - 1, 0) * c_w] +
                 114 * tmp[y * c_w] + 35 * tmp[OC_MINI(y + 1, c_h - 1) * c_w] -
                 9 * tmp[OC_MINI(y + 2, c_h - 1) * c_w] +
                 tmp[OC_MINI(y + 3, c_h - 1) * c_w] + 64) >>
                    7,
                255);
          }
          for (; y < c_h - 3; y++) {
            _dst[y * c_w] = (unsigned char)OC_CLAMPI(
                0,
                (4 * tmp[(y - 2) * c_w] - 17 * tmp[(y - 1) * c_w] +
                 114 * tmp[y * c_w] + 35 * tmp[(y + 1) * c_w] -
                 9 * tmp[(y + 2) * c_w] + tmp[(y + 3) * c_w] + 64) >>
                    7,
                255);
          }
          for (; y < c_h; y++) {
            _dst[y * c_w] = (unsigned char)OC_CLAMPI(
                0,
                (4 * tmp[(y - 2) * c_w] - 17 * tmp[(y - 1) * c_w] +
                 114 * tmp[y * c_w] + 35 * tmp[OC_MINI(y + 1, c_h - 1) * c_w] -
                 9 * tmp[OC_MINI(y + 2, c_h - 1) * c_w] + tmp[(c_h - 1) * c_w] +
                 64) >>
                    7,
                255);
          }
          _dst++;
          tmp++;
        }
        break;
      }
    }
    /*For actual interlaced material, this would have to be done separately on
       each field, and the shift amounts would be different.
      C_r moves down 1/8, C_b up 3/8 in the top field, and C_r moves down 3/8,
       C_b up 1/8 in the bottom field.
      The corresponding filters would be:
       Down 1/8 (reverse order for up): [3 -11 125 15 -4 0]/128
       Down 3/8 (reverse order for up): [4 -19 98 56 -13 2]/128*/
  }
}
/*Perform vertical filtering to reduce a single plane from 4:2:2 to 4:2:0.
  This is used as a helper by several converation routines.*/
void y4m_422jpeg_420jpeg_helper(unsigned char *_dst,
                                       const unsigned char *_src, int _c_w,
                                       int _c_h) {
  int y;
  int x;
  /*Filter: [3 -17 78 78 -17 3]/128, derived from a 6-tap Lanczos window.*/
  for (x = 0; x < _c_w; x++) {
    for (y = 0; y < OC_MINI(_c_h, 2); y += 2) {
      _dst[(y >> 1) * _c_w] =
          OC_CLAMPI(0,
                    (64 * _src[0] + 78 * _src[OC_MINI(1, _c_h - 1) * _c_w] -
                     17 * _src[OC_MINI(2, _c_h - 1) * _c_w] +
                     3 * _src[OC_MINI(3, _c_h - 1) * _c_w] + 64) >>
                        7,
                    255);
    }
    for (; y < _c_h - 3; y += 2) {
      _dst[(y >> 1) * _c_w] =
          OC_CLAMPI(0,
                    (3 * (_src[(y - 2) * _c_w] + _src[(y + 3) * _c_w]) -
                     17 * (_src[(y - 1) * _c_w] + _src[(y + 2) * _c_w]) +
                     78 * (_src[y * _c_w] + _src[(y + 1) * _c_w]) + 64) >>
                        7,
                    255);
    }
    for (; y < _c_h; y += 2) {
      _dst[(y >> 1) * _c_w] = OC_CLAMPI(
          0,
          (3 * (_src[(y - 2) * _c_w] + _src[(_c_h - 1) * _c_w]) -
           17 * (_src[(y - 1) * _c_w] + _src[OC_MINI(y + 2, _c_h - 1) * _c_w]) +
           78 * (_src[y * _c_w] + _src[OC_MINI(y + 1, _c_h - 1) * _c_w]) +
           64) >>
              7,
          255);
    }
    _src++;
    _dst++;
  }
}
/*420jpeg chroma samples are sited like:
  Y-------Y-------Y-------Y-------
  |       |       |       |
  |   BR  |       |   BR  |
  |       |       |       |
  Y-------Y-------Y-------Y-------
  |       |       |       |
  |       |       |       |
  |       |       |       |
  Y-------Y-------Y-------Y-------
  |       |       |       |
  |   BR  |       |   BR  |
  |       |       |       |
  Y-------Y-------Y-------Y-------
  |       |       |       |
  |       |       |       |
  |       |       |       |
  422jpeg chroma samples are sited like:
  Y---BR--Y-------Y---BR--Y-------
  |       |       |       |
  |       |       |       |
  |       |       |       |
  Y---BR--Y-------Y---BR--Y-------
  |       |       |       |
  |       |       |       |
  |       |       |       |
  Y---BR--Y-------Y---BR--Y-------
  |       |       |       |
  |       |       |       |
  |       |       |       |
  Y---BR--Y-------Y---BR--Y-------
  |       |       |       |
  |       |       |       |
  |       |       |       |
  We use a resampling filter to decimate the chroma planes by two in the
   vertical direction.*/
void y4m_convert_422jpeg_420jpeg(y4m_input *_y4m, unsigned char *_dst,
                                        unsigned char *_aux) {
  int c_w;
  int c_h;
  int c_sz;
  int dst_c_w;
  int dst_c_h;
  int dst_c_sz;
  int pli;
  /*Skip past the luma data.*/
  _dst += _y4m->pic_w * _y4m->pic_h;
  /*Compute the size of each chroma plane.*/
  c_w = (_y4m->pic_w + _y4m->src_c_dec_h - 1) / _y4m->src_c_dec_h;
  c_h = _y4m->pic_h;
  dst_c_w = (_y4m->pic_w + _y4m->dst_c_dec_h - 1) / _y4m->dst_c_dec_h;
  dst_c_h = (_y4m->pic_h + _y4m->dst_c_dec_v - 1) / _y4m->dst_c_dec_v;
  c_sz = c_w * c_h;
  dst_c_sz = dst_c_w * dst_c_h;
  for (pli = 1; pli < 3; pli++) {
    y4m_422jpeg_420jpeg_helper(_dst, _aux, c_w, c_h);
    _aux += c_sz;
    _dst += dst_c_sz;
  }
}
/*420jpeg chroma samples are sited like:
  Y-------Y-------Y-------Y-------
  |       |       |       |
  |   BR  |       |   BR  |
  |       |       |       |
  Y-------Y-------Y-------Y-------
  |       |       |       |
  |       |       |       |
  |       |       |       |
  Y-------Y-------Y-------Y-------
  |       |       |       |
  |   BR  |       |   BR  |
  |       |       |       |
  Y-------Y-------Y-------Y-------
  |       |       |       |
  |       |       |       |
  |       |       |       |
  422 chroma samples are sited like:
  YBR-----Y-------YBR-----Y-------
  |       |       |       |
  |       |       |       |
  |       |       |       |
  YBR-----Y-------YBR-----Y-------
  |       |       |       |
  |       |       |       |
  |       |       |       |
  YBR-----Y-------YBR-----Y-------
  |       |       |       |
  |       |       |       |
  |       |       |       |
  YBR-----Y-------YBR-----Y-------
  |       |       |       |
  |       |       |       |
  |       |       |       |
  We use a resampling filter to shift the original site locations one quarter
   pixel (at the original chroma resolution) to the right.
  Then we use a second resampling filter to decimate the chroma planes by two
   in the vertical direction.*/
void y4m_convert_422_420jpeg(y4m_input *_y4m, unsigned char *_dst,
                                    unsigned char *_aux) {
  unsigned char *tmp;
  int c_w;
  int c_h;
  int c_sz;
  int dst_c_h;
  int dst_c_sz;
  int pli;
  /*Skip past the luma data.*/
  _dst += _y4m->pic_w * _y4m->pic_h;
  /*Compute the size of each chroma plane.*/
  c_w = (_y4m->pic_w + _y4m->src_c_dec_h - 1) / _y4m->src_c_dec_h;
  c_h = _y4m->pic_h;
  dst_c_h = (_y4m->pic_h + _y4m->dst_c_dec_v - 1) / _y4m->dst_c_dec_v;
  c_sz = c_w * c_h;
  dst_c_sz = c_w * dst_c_h;
  tmp = _aux + 2 * c_sz;
  for (pli = 1; pli < 3; pli++) {
    /*In reality, the horizontal and vertical steps could be pipelined, for
       less memory consumption and better cache performance, but we do them
       separately for simplicity.*/
    /*First do horizontal filtering (convert to 422jpeg)*/
    y4m_42xmpeg2_42xjpeg_helper(tmp, _aux, c_w, c_h);
    /*Now do the vertical filtering.*/
    y4m_422jpeg_420jpeg_helper(_dst, tmp, c_w, c_h);
    _aux += c_sz;
    _dst += dst_c_sz;
  }
}
/*420jpeg chroma samples are sited like:
  Y-------Y-------Y-------Y-------
  |       |       |       |
  |   BR  |       |   BR  |
  |       |       |       |
  Y-------Y-------Y-------Y-------
  |       |       |       |
  |       |       |       |
  |       |       |       |
  Y-------Y-------Y-------Y-------
  |       |       |       |
  |   BR  |       |   BR  |
  |       |       |       |
  Y-------Y-------Y-------Y-------
  |       |       |       |
  |       |       |       |
  |       |       |       |
  411 chroma samples are sited like:
  YBR-----Y-------Y-------Y-------
  |       |       |       |
  |       |       |       |
  |       |       |       |
  YBR-----Y-------Y-------Y-------
  |       |       |       |
  |       |       |       |
  |       |       |       |
  YBR-----Y-------Y-------Y-------
  |       |       |       |
  |       |       |       |
  |       |       |       |
  YBR-----Y-------Y-------Y-------
  |       |       |       |
  |       |       |       |
  |       |       |       |
  We use a filter to resample at site locations one eighth pixel (at the source
   chroma plane's horizontal resolution) and five eighths of a pixel to the
   right.
  Then we use another filter to decimate the planes by 2 in the vertical
   direction.*/
void y4m_convert_411_420jpeg(y4m_input *_y4m, unsigned char *_dst,
                                    unsigned char *_aux) {
  unsigned char *tmp;
  int c_w;
  int c_h;
  int c_sz;
  int dst_c_w;
  int dst_c_h;
  int dst_c_sz;
  int tmp_sz;
  int pli;
  int y;
  int x;
  /*Skip past the luma data.*/
  _dst += _y4m->pic_w * _y4m->pic_h;
  /*Compute the size of each chroma plane.*/
  c_w = (_y4m->pic_w + _y4m->src_c_dec_h - 1) / _y4m->src_c_dec_h;
  c_h = _y4m->pic_h;
  dst_c_w = (_y4m->pic_w + _y4m->dst_c_dec_h - 1) / _y4m->dst_c_dec_h;
  dst_c_h = (_y4m->pic_h + _y4m->dst_c_dec_v - 1) / _y4m->dst_c_dec_v;
  c_sz = c_w * c_h;
  dst_c_sz = dst_c_w * dst_c_h;
  tmp_sz = dst_c_w * c_h;
  tmp = _aux + 2 * c_sz;
  for (pli = 1; pli < 3; pli++) {
    /*In reality, the horizontal and vertical steps could be pipelined, for
       less memory consumption and better cache performance, but we do them
       separately for simplicity.*/
    /*First do horizontal filtering (convert to 422jpeg)*/
    for (y = 0; y < c_h; y++) {
      /*Filters: [1 110 18 -1]/128 and [-3 50 86 -5]/128, both derived from a
         4-tap Mitchell window.*/
      for (x = 0; x < OC_MINI(c_w, 1); x++) {
        tmp[x << 1] = (unsigned char)OC_CLAMPI(
            0,
            (111 * _aux[0] + 18 * _aux[OC_MINI(1, c_w - 1)] -
             _aux[OC_MINI(2, c_w - 1)] + 64) >>
                7,
            255);
        tmp[x << 1 | 1] = (unsigned char)OC_CLAMPI(
            0,
            (47 * _aux[0] + 86 * _aux[OC_MINI(1, c_w - 1)] -
             5 * _aux[OC_MINI(2, c_w - 1)] + 64) >>
                7,
            255);
      }
      for (; x < c_w - 2; x++) {
        tmp[x << 1] =
            (unsigned char)OC_CLAMPI(0,
                                     (_aux[x - 1] + 110 * _aux[x] +
                                      18 * _aux[x + 1] - _aux[x + 2] + 64) >>
                                         7,
                                     255);
        tmp[x << 1 | 1] = (unsigned char)OC_CLAMPI(
            0,
            (-3 * _aux[x - 1] + 50 * _aux[x] + 86 * _aux[x + 1] -
             5 * _aux[x + 2] + 64) >>
                7,
            255);
      }
      for (; x < c_w; x++) {
        tmp[x << 1] = (unsigned char)OC_CLAMPI(
            0,
            (_aux[x - 1] + 110 * _aux[x] + 18 * _aux[OC_MINI(x + 1, c_w - 1)] -
             _aux[c_w - 1] + 64) >>
                7,
            255);
        if ((x << 1 | 1) < dst_c_w) {
          tmp[x << 1 | 1] = (unsigned char)OC_CLAMPI(
              0,
              (-3 * _aux[x - 1] + 50 * _aux[x] +
               86 * _aux[OC_MINI(x + 1, c_w - 1)] - 5 * _aux[c_w - 1] + 64) >>
                  7,
              255);
        }
      }
      tmp += dst_c_w;
      _aux += c_w;
    }
    tmp -= tmp_sz;
    /*Now do the vertical filtering.*/
    y4m_422jpeg_420jpeg_helper(_dst, tmp, dst_c_w, c_h);
    _dst += dst_c_sz;
  }
}
/*Convert 444 to 420jpeg.*/
void y4m_convert_444_420jpeg(y4m_input *_y4m, unsigned char *_dst,
                                    unsigned char *_aux) {
  unsigned char *tmp;
  int c_w;
  int c_h;
  int c_sz;
  int dst_c_w;
  int dst_c_h;
  int dst_c_sz;
  int tmp_sz;
  int pli;
  int y;
  int x;
  /*Skip past the luma data.*/
  _dst += _y4m->pic_w * _y4m->pic_h;
  /*Compute the size of each chroma plane.*/
  c_w = (_y4m->pic_w + _y4m->src_c_dec_h - 1) / _y4m->src_c_dec_h;
  c_h = _y4m->pic_h;
  dst_c_w = (_y4m->pic_w + _y4m->dst_c_dec_h - 1) / _y4m->dst_c_dec_h;
  dst_c_h = (_y4m->pic_h + _y4m->dst_c_dec_v - 1) / _y4m->dst_c_dec_v;
  c_sz = c_w * c_h;
  dst_c_sz = dst_c_w * dst_c_h;
  tmp_sz = dst_c_w * c_h;
  tmp = _aux + 2 * c_sz;
  for (pli = 1; pli < 3; pli++) {
    /*Filter: [3 -17 78 78 -17 3]/128, derived from a 6-tap Lanczos window.*/
    for (y = 0; y < c_h; y++) {
      for (x = 0; x < OC_MINI(c_w, 2); x += 2) {
        tmp[x >> 1] = OC_CLAMPI(0,
                                (64 * _aux[0] + 78 * _aux[OC_MINI(1, c_w - 1)] -
                                 17 * _aux[OC_MINI(2, c_w - 1)] +
                                 3 * _aux[OC_MINI(3, c_w - 1)] + 64) >>
                                    7,
                                255);
      }
      for (; x < c_w - 3; x += 2) {
        tmp[x >> 1] = OC_CLAMPI(0,
                                (3 * (_aux[x - 2] + _aux[x + 3]) -
                                 17 * (_aux[x - 1] + _aux[x + 2]) +
                                 78 * (_aux[x] + _aux[x + 1]) + 64) >>
                                    7,
                                255);
      }
      for (; x < c_w; x += 2) {
        tmp[x >> 1] =
            OC_CLAMPI(0,
                      (3 * (_aux[x - 2] + _aux[c_w - 1]) -
                       17 * (_aux[x - 1] + _aux[OC_MINI(x + 2, c_w - 1)]) +
                       78 * (_aux[x] + _aux[OC_MINI(x + 1, c_w - 1)]) + 64) >>
                          7,
                      255);
      }
      tmp += dst_c_w;
      _aux += c_w;
    }
    tmp -= tmp_sz;
    /*Now do the vertical filtering.*/
    y4m_422jpeg_420jpeg_helper(_dst, tmp, dst_c_w, c_h);
    _dst += dst_c_sz;
  }
}
/*The image is padded with empty chroma components at 4:2:0.*/
void y4m_convert_mono_420jpeg(y4m_input *_y4m, unsigned char *_dst,
                                     unsigned char *_aux) {
  int c_sz;
  (void)_aux;
  _dst += _y4m->pic_w * _y4m->pic_h;
  c_sz = ((_y4m->pic_w + _y4m->dst_c_dec_h - 1) / _y4m->dst_c_dec_h) *
         ((_y4m->pic_h + _y4m->dst_c_dec_v - 1) / _y4m->dst_c_dec_v);
  memset(_dst, 128, c_sz * 2);
}
/*No conversion function needed.*/
void y4m_convert_null(y4m_input *_y4m, unsigned char *_dst,
                             unsigned char *_aux) {
  (void)_y4m;
  (void)_dst;
  (void)_aux;
}
int y4m_input_open(y4m_input *_y4m, DataSource *_fin, char *_skip, int _nskip,
                   int only_420, bool log) {
  char buffer[80] = { 0 };
  int ret;
  int i;
  /*Read until newline, or 80 cols, whichever happens first.*/
  for (i = 0; i < 79; i++) {
    if (_nskip > 0) {
      buffer[i] = *_skip++;
      _nskip--;
    } else {
      if (!file_read(buffer + i, 1, _fin)) return -1;
    }
    if (buffer[i] == '\n') break;
  }
  /*We skipped too much header data.*/
  if (_nskip > 0) return -1;
  if (i == 79) {
    if (log) fprintf(stderr, "Error parsing header; not a YUV2MPEG2 file?\n");
    return -1;
  }
  buffer[i] = '\0';
  if (memcmp(buffer, "YUV4MPEG", 8)) {
    if (log) fprintf(stderr, "Incomplete magic for YUV4MPEG file.\n");
    return -1;
  }
  if (buffer[8] != '2' && log) {
    fprintf(stderr, "Incorrect YUV input file version; YUV4MPEG2 required.\n");
  }
  ret = y4m_parse_tags(_y4m, buffer + 5);
  if (_y4m->pic_w < 1 || _y4m->pic_h < 1) {
    if (log) fprintf(stderr, "Invalid images dimensions.\n");
    return -1;
  }
  if ((size_t)_y4m->pic_w > kMaxBufferDimension ||
      (size_t)_y4m->pic_h > kMaxBufferDimension ||
      (size_t)_y4m->pic_w * (size_t)_y4m->pic_h > kMaxBufferDimension) {
    if (log) fprintf(stderr, "Image dimensions are too big.\n");
    return -1;
  }
  if (ret < 0) {
    if (log) fprintf(stderr, "Error parsing YUV4MPEG2 header.\n");
    return ret;
  }
  if (_y4m->interlace == '?') {
    if (log) {
      fprintf(stderr,
              "Warning: Input video interlacing format unknown; "
              "assuming progressive scan.\n");
    }
  } else if (_y4m->interlace != 'p') {
    if (log) {
      fprintf(stderr,
              "Input video is interlaced; "
              "Only progressive scan handled.\n");
    }
    return -1;
  }
  _y4m->vpx_fmt = VPX_IMG_FMT_I420;
  _y4m->bps = 12;
  _y4m->bit_depth = 8;
  if (strcmp(_y4m->chroma_type, "420") == 0 ||
      strcmp(_y4m->chroma_type, "420jpeg") == 0) {
    _y4m->src_c_dec_h = _y4m->dst_c_dec_h = _y4m->src_c_dec_v =
        _y4m->dst_c_dec_v = 2;
    _y4m->dst_buf_read_sz =
        _y4m->pic_w * _y4m->pic_h +
        2 * ((_y4m->pic_w + 1) / 2) * ((_y4m->pic_h + 1) / 2);
    /* Natively supported: no conversion required. */
    _y4m->aux_buf_sz = _y4m->aux_buf_read_sz = 0;
    _y4m->convert = y4m_convert_null;
  } else if (strcmp(_y4m->chroma_type, "420p10") == 0) {
    _y4m->src_c_dec_h = 2;
    _y4m->dst_c_dec_h = 2;
    _y4m->src_c_dec_v = 2;
    _y4m->dst_c_dec_v = 2;
    _y4m->dst_buf_read_sz =
        2 * (_y4m->pic_w * _y4m->pic_h +
             2 * ((_y4m->pic_w + 1) / 2) * ((_y4m->pic_h + 1) / 2));
    /* Natively supported: no conversion required. */
    _y4m->aux_buf_sz = _y4m->aux_buf_read_sz = 0;
    _y4m->convert = y4m_convert_null;
    _y4m->bit_depth = 10;
    _y4m->bps = 15;
    _y4m->vpx_fmt = VPX_IMG_FMT_I42016;
    if (only_420) {
      if (log) {
        fprintf(stderr, "Unsupported conversion from 420p10 to 420jpeg\n");
      }
      return -1;
    }
  } else if (strcmp(_y4m->chroma_type, "420p12") == 0) {
    _y4m->src_c_dec_h = 2;
    _y4m->dst_c_dec_h = 2;
    _y4m->src_c_dec_v = 2;
    _y4m->dst_c_dec_v = 2;
    _y4m->dst_buf_read_sz =
        2 * (_y4m->pic_w * _y4m->pic_h +
             2 * ((_y4m->pic_w + 1) / 2) * ((_y4m->pic_h + 1) / 2));
    /* Natively supported: no conversion required. */
    _y4m->aux_buf_sz = _y4m->aux_buf_read_sz = 0;
    _y4m->convert = y4m_convert_null;
    _y4m->bit_depth = 12;
    _y4m->bps = 18;
    _y4m->vpx_fmt = VPX_IMG_FMT_I42016;
    if (only_420) {
      if (log) {
        fprintf(stderr, "Unsupported conversion from 420p12 to 420jpeg\n");
      }
      return -1;
    }
  } else if (strcmp(_y4m->chroma_type, "420mpeg2") == 0) {
    _y4m->src_c_dec_h = _y4m->dst_c_dec_h = _y4m->src_c_dec_v =
        _y4m->dst_c_dec_v = 2;
    _y4m->dst_buf_read_sz = _y4m->pic_w * _y4m->pic_h;
    /*Chroma filter required: read into the aux buf first.*/
    _y4m->aux_buf_sz = _y4m->aux_buf_read_sz =
        2 * ((_y4m->pic_w + 1) / 2) * ((_y4m->pic_h + 1) / 2);
    _y4m->convert = y4m_convert_42xmpeg2_42xjpeg;
  } else if (strcmp(_y4m->chroma_type, "420paldv") == 0) {
    _y4m->src_c_dec_h = _y4m->dst_c_dec_h = _y4m->src_c_dec_v =
        _y4m->dst_c_dec_v = 2;
    _y4m->dst_buf_read_sz = _y4m->pic_w * _y4m->pic_h;
    /*Chroma filter required: read into the aux buf first.
      We need to make two filter passes, so we need some extra space in the
       aux buffer.*/
    _y4m->aux_buf_sz = 3 * ((_y4m->pic_w + 1) / 2) * ((_y4m->pic_h + 1) / 2);
    _y4m->aux_buf_read_sz =
        2 * ((_y4m->pic_w + 1) / 2) * ((_y4m->pic_h + 1) / 2);
    _y4m->convert = y4m_convert_42xpaldv_42xjpeg;
  } else if (strcmp(_y4m->chroma_type, "422jpeg") == 0) {
    _y4m->src_c_dec_h = _y4m->dst_c_dec_h = 2;
    _y4m->src_c_dec_v = 1;
    _y4m->dst_c_dec_v = 2;
    _y4m->dst_buf_read_sz = _y4m->pic_w * _y4m->pic_h;
    /*Chroma filter required: read into the aux buf first.*/
    _y4m->aux_buf_sz = _y4m->aux_buf_read_sz =
        2 * ((_y4m->pic_w + 1) / 2) * _y4m->pic_h;
    _y4m->convert = y4m_convert_422jpeg_420jpeg;
  } else if (strcmp(_y4m->chroma_type, "422") == 0) {
    _y4m->src_c_dec_h = 2;
    _y4m->src_c_dec_v = 1;
    if (only_420) {
      _y4m->dst_c_dec_h = 2;
      _y4m->dst_c_dec_v = 2;
      _y4m->dst_buf_read_sz = _y4m->pic_w * _y4m->pic_h;
      /*Chroma filter required: read into the aux buf first.
        We need to make two filter passes, so we need some extra space in the
         aux buffer.*/
      _y4m->aux_buf_read_sz = 2 * ((_y4m->pic_w + 1) / 2) * _y4m->pic_h;
      _y4m->aux_buf_sz =
          _y4m->aux_buf_read_sz + ((_y4m->pic_w + 1) / 2) * _y4m->pic_h;
      _y4m->convert = y4m_convert_422_420jpeg;
    } else {
      _y4m->vpx_fmt = VPX_IMG_FMT_I422;
      _y4m->bps = 16;
      _y4m->dst_c_dec_h = _y4m->src_c_dec_h;
      _y4m->dst_c_dec_v = _y4m->src_c_dec_v;
      _y4m->dst_buf_read_sz =
          _y4m->pic_w * _y4m->pic_h + 2 * ((_y4m->pic_w + 1) / 2) * _y4m->pic_h;
      /*Natively supported: no conversion required.*/
      _y4m->aux_buf_sz = _y4m->aux_buf_read_sz = 0;
      _y4m->convert = y4m_convert_null;
    }
  } else if (strcmp(_y4m->chroma_type, "422p10") == 0) {
    _y4m->src_c_dec_h = 2;
    _y4m->src_c_dec_v = 1;
    _y4m->vpx_fmt = VPX_IMG_FMT_I42216;
    _y4m->bps = 20;
    _y4m->bit_depth = 10;
    _y4m->dst_c_dec_h = _y4m->src_c_dec_h;
    _y4m->dst_c_dec_v = _y4m->src_c_dec_v;
    _y4m->dst_buf_read_sz = 2 * (_y4m->pic_w * _y4m->pic_h +
                                 2 * ((_y4m->pic_w + 1) / 2) * _y4m->pic_h);
    _y4m->aux_buf_sz = _y4m->aux_buf_read_sz = 0;
    _y4m->convert = y4m_convert_null;
    if (only_420) {
      if (log) {
        fprintf(stderr, "Unsupported conversion from 422p10 to 420jpeg\n");
      }
      return -1;
    }
  } else if (strcmp(_y4m->chroma_type, "422p12") == 0) {
    _y4m->src_c_dec_h = 2;
    _y4m->src_c_dec_v = 1;
    _y4m->vpx_fmt = VPX_IMG_FMT_I42216;
    _y4m->bps = 24;
    _y4m->bit_depth = 12;
    _y4m->dst_c_dec_h = _y4m->src_c_dec_h;
    _y4m->dst_c_dec_v = _y4m->src_c_dec_v;
    _y4m->dst_buf_read_sz = 2 * (_y4m->pic_w * _y4m->pic_h +
                                 2 * ((_y4m->pic_w + 1) / 2) * _y4m->pic_h);
    _y4m->aux_buf_sz = _y4m->aux_buf_read_sz = 0;
    _y4m->convert = y4m_convert_null;
    if (only_420) {
      if (log) {
        fprintf(stderr, "Unsupported conversion from 422p12 to 420jpeg\n");
      }
      return -1;
    }
  } else if (strcmp(_y4m->chroma_type, "411") == 0) {
    _y4m->src_c_dec_h = 4;
    _y4m->dst_c_dec_h = 2;
    _y4m->src_c_dec_v = 1;
    _y4m->dst_c_dec_v = 2;
    _y4m->dst_buf_read_sz = _y4m->pic_w * _y4m->pic_h;
    /*Chroma filter required: read into the aux buf first.
      We need to make two filter passes, so we need some extra space in the
       aux buffer.*/
    _y4m->aux_buf_read_sz = 2 * ((_y4m->pic_w + 3) / 4) * _y4m->pic_h;
    _y4m->aux_buf_sz =
        _y4m->aux_buf_read_sz + ((_y4m->pic_w + 1) / 2) * _y4m->pic_h;
    _y4m->convert = y4m_convert_411_420jpeg;
    if (log) fprintf(stderr, "Unsupported conversion from yuv 411\n");
    return -1;
  } else if (strcmp(_y4m->chroma_type, "444") == 0) {
    _y4m->src_c_dec_h = 1;
    _y4m->src_c_dec_v = 1;
    if (only_420) {
      _y4m->dst_c_dec_h = 2;
      _y4m->dst_c_dec_v = 2;
      _y4m->dst_buf_read_sz = _y4m->pic_w * _y4m->pic_h;
      /*Chroma filter required: read into the aux buf first.
        We need to make two filter passes, so we need some extra space in the
         aux buffer.*/
      _y4m->aux_buf_read_sz = 2 * _y4m->pic_w * _y4m->pic_h;
      _y4m->aux_buf_sz =
          _y4m->aux_buf_read_sz + ((_y4m->pic_w + 1) / 2) * _y4m->pic_h;
      _y4m->convert = y4m_convert_444_420jpeg;
    } else {
      _y4m->vpx_fmt = VPX_IMG_FMT_I444;
      _y4m->bps = 24;
      _y4m->dst_c_dec_h = _y4m->src_c_dec_h;
      _y4m->dst_c_dec_v = _y4m->src_c_dec_v;
      _y4m->dst_buf_read_sz = 3 * _y4m->pic_w * _y4m->pic_h;
      /*Natively supported: no conversion required.*/
      _y4m->aux_buf_sz = _y4m->aux_buf_read_sz = 0;
      _y4m->convert = y4m_convert_null;
    }
  } else if (strcmp(_y4m->chroma_type, "444p10") == 0) {
    _y4m->src_c_dec_h = 1;
    _y4m->src_c_dec_v = 1;
    _y4m->vpx_fmt = VPX_IMG_FMT_I44416;
    _y4m->bps = 30;
    _y4m->bit_depth = 10;
    _y4m->dst_c_dec_h = _y4m->src_c_dec_h;
    _y4m->dst_c_dec_v = _y4m->src_c_dec_v;
    _y4m->dst_buf_read_sz = 2 * 3 * _y4m->pic_w * _y4m->pic_h;
    _y4m->aux_buf_sz = _y4m->aux_buf_read_sz = 0;
    _y4m->convert = y4m_convert_null;
    if (only_420) {
      if (log) {
        fprintf(stderr, "Unsupported conversion from 444p10 to 420jpeg\n");
      }
      return -1;
    }
  } else if (strcmp(_y4m->chroma_type, "444p12") == 0) {
    _y4m->src_c_dec_h = 1;
    _y4m->src_c_dec_v = 1;
    _y4m->vpx_fmt = VPX_IMG_FMT_I44416;
    _y4m->bps = 36;
    _y4m->bit_depth = 12;
    _y4m->dst_c_dec_h = _y4m->src_c_dec_h;
    _y4m->dst_c_dec_v = _y4m->src_c_dec_v;
    _y4m->dst_buf_read_sz = 2 * 3 * _y4m->pic_w * _y4m->pic_h;
    _y4m->aux_buf_sz = _y4m->aux_buf_read_sz = 0;
    _y4m->convert = y4m_convert_null;
    if (only_420) {
      if (log) {
        fprintf(stderr, "Unsupported conversion from 444p12 to 420jpeg\n");
      }
      return -1;
    }
  } else if (strcmp(_y4m->chroma_type, "mono") == 0) {
    _y4m->src_c_dec_h = _y4m->src_c_dec_v = 0;
    _y4m->dst_c_dec_h = _y4m->dst_c_dec_v = 2;
    _y4m->dst_buf_read_sz = _y4m->pic_w * _y4m->pic_h;
    /*No extra space required, but we need to clear the chroma planes.*/
    _y4m->aux_buf_sz = _y4m->aux_buf_read_sz = 0;
    _y4m->convert = y4m_convert_mono_420jpeg;
  } else {
    if (log) {
      fprintf(stderr, "Unknown chroma sampling type: %s\n", _y4m->chroma_type);
    }
    return -1;
  }
  /*The size of the final frame buffers is always computed from the
     destination chroma decimation type.*/
  _y4m->dst_buf_sz =
      _y4m->pic_w * _y4m->pic_h +
      2 * ((_y4m->pic_w + _y4m->dst_c_dec_h - 1) / _y4m->dst_c_dec_h) *
          ((_y4m->pic_h + _y4m->dst_c_dec_v - 1) / _y4m->dst_c_dec_v);
  if (_y4m->bit_depth == 8)
    _y4m->dst_buf = (unsigned char *)malloc(_y4m->dst_buf_sz);
  else
    _y4m->dst_buf = (unsigned char *)malloc(2 * _y4m->dst_buf_sz);
  if (_y4m->aux_buf_sz > 0)
    _y4m->aux_buf = (unsigned char *)malloc(_y4m->aux_buf_sz);
  return 0;
}
void y4m_input_close(y4m_input *_y4m) {
  free(_y4m->dst_buf);
  free(_y4m->aux_buf);
}
int y4m_input_fetch_frame(y4m_input *_y4m, DataSource *_fin,
                          vpx_image_t *_img, bool log) {
  char frame[6];
  int pic_sz;
  int c_w;
  int c_h;
  int c_sz;
  int bytes_per_sample = _y4m->bit_depth > 8 ? 2 : 1;
  /*Read and skip the frame header.*/
  if (!file_read(frame, 6, _fin)) return 0;
  if (memcmp(frame, "FRAME", 5)) {
    if (log) fprintf(stderr, "Loss of framing in Y4M input data\n");
    return -1;
  }
  if (frame[5] != '\n') {
    char c;
    int j;
    for (j = 0; j < 79 && file_read(&c, 1, _fin) && c != '\n'; j++) {
    }
    if (j == 79) {
      if (log) fprintf(stderr, "Error parsing Y4M frame header\n");
      return -1;
    }
  }
  /*Read the frame data that needs no conversion.*/
  if (!file_read(_y4m->dst_buf, _y4m->dst_buf_read_sz, _fin)) {
    if (log) fprintf(stderr, "Error reading Y4M frame data.\n");
    return -1;
  }
  /*Read the frame data that does need conversion.*/
  if (!file_read(_y4m->aux_buf, _y4m->aux_buf_read_sz, _fin)) {
    if (log) fprintf(stderr, "Error reading Y4M frame data.\n");
    return -1;
  }
  /*Now convert the just read frame.*/
  (*_y4m->convert)(_y4m, _y4m->dst_buf, _y4m->aux_buf);
  /*Fill in the frame buffer pointers.
    We don't use vpx_img_wrap() because it forces padding for odd picture
     sizes, which would require a separate fread call for every row.*/
  memset(_img, 0, sizeof(*_img));
  /*Y4M has the planes in Y'CbCr order, which libvpx calls Y, U, and V.*/
  _img->fmt = _y4m->vpx_fmt;
  _img->w = _img->d_w = _y4m->pic_w;
  _img->h = _img->d_h = _y4m->pic_h;
  _img->x_chroma_shift = _y4m->dst_c_dec_h >> 1;
  _img->y_chroma_shift = _y4m->dst_c_dec_v >> 1;
  _img->bps = _y4m->bps;
  /*Set up the buffer pointers.*/
  pic_sz = _y4m->pic_w * _y4m->pic_h * bytes_per_sample;
  c_w = (_y4m->pic_w + _y4m->dst_c_dec_h - 1) / _y4m->dst_c_dec_h;
  c_w *= bytes_per_sample;
  c_h = (_y4m->pic_h + _y4m->dst_c_dec_v - 1) / _y4m->dst_c_dec_v;
  c_sz = c_w * c_h;
  _img->stride[VPX_PLANE_Y] = _img->stride[VPX_PLANE_ALPHA] =
      _y4m->pic_w * bytes_per_sample;
  _img->stride[VPX_PLANE_U] = _img->stride[VPX_PLANE_V] = c_w;
  _img->planes[VPX_PLANE_Y] = _y4m->dst_buf;
  _img->planes[VPX_PLANE_U] = _y4m->dst_buf + pic_sz;
  _img->planes[VPX_PLANE_V] = _y4m->dst_buf + pic_sz + c_sz;
  _img->planes[VPX_PLANE_ALPHA] = _y4m->dst_buf + pic_sz + 2 * c_sz;
  return 1;
}

}  // namespace

// -----------------------------------------------------------------------------
// The following is not copied from chromium.googlesource.com/webm/libvpx, this
// is the WP2 imageio implementation.

class ImageReaderY4M : public ImageReader::Impl {
 public:
  ImageReaderY4M(const uint8_t *data, size_t data_size,
                 YUVPlane *const ccsp_buffer, CSPMtx *const ccsp_to_rgb,
                 Metadata *const metadata, LogLevel log_level,
                 size_t max_num_pixels)
      : ImageReader::Impl(nullptr, data, data_size, log_level, max_num_pixels),
        ccsp_buffer_(ccsp_buffer),
        ccsp_to_rgb_(ccsp_to_rgb),
        metadata_(metadata),
        data_source_(data, data_size),
        info_({}),
        image_({}),
        opened_(false) {}
  ~ImageReaderY4M() override {
    if (opened_) y4m_input_close(&info_);
  }

  WP2Status ReadFrame(bool* const is_last,
                      uint32_t* const duration_ms) override {
    WP2_CHECK_OK(ccsp_buffer_ != nullptr, WP2_STATUS_NULL_PARAMETER);
    WP2_CHECK_OK(ccsp_to_rgb_ != nullptr, WP2_STATUS_NULL_PARAMETER);
    if (!opened_) {
      // y4m_input_open() does not seem to allocate while returning an error, so
      // always calling y4m_input_close() is not currently necessary but in case
      // the behavior changes, it is safer to always set 'opened_' to true.
      opened_ = true;
      WP2_CHECK_OK(y4m_input_open(&info_, &data_source_, /*_skip=*/nullptr,
                                  /*_nskip=*/0, /*only_420=*/false,
                                  log_level_ >= LogLevel::DEFAULT) >= 0,
                   WP2_STATUS_BITSTREAM_ERROR);
      WP2_CHECK_STATUS(CheckDimensions(info_.pic_w, info_.pic_h));
    }

    WP2_CHECK_OK(y4m_input_fetch_frame(&info_, &data_source_, &image_,
                                       log_level_ >= LogLevel::DEFAULT) >= 0,
                 WP2_STATUS_BITSTREAM_ERROR);

    // Most fields of 'vpx_image_t' are left to their default value, like
    // VPX_CS_UNKNOWN which should be YCbCr VPX_CS_BT_601.
    // VPX_CR_STUDIO_RANGE seems to be the correct default.
    // See https://en.wikipedia.org/wiki/YCbCr
    WP2_CHECK_OK(image_.cs == VPX_CS_UNKNOWN, WP2_STATUS_UNSUPPORTED_FEATURE);

    WP2_CHECK_OK(info_.bit_depth >= CCSPImageReader::kMinBitDepth &&
                     info_.bit_depth <= CCSPImageReader::kMaxBitDepth,
                 WP2_STATUS_UNSUPPORTED_FEATURE);
    const int w_shift = image_.x_chroma_shift, h_shift = image_.y_chroma_shift;
    if (image_.fmt == VPX_IMG_FMT_I444 || image_.fmt == VPX_IMG_FMT_I44416) {
      WP2_CHECK_OK(w_shift == 0 && h_shift == 0,
                   WP2_STATUS_UNSUPPORTED_FEATURE);
    } else {
      WP2_CHECK_OK(
          image_.fmt == VPX_IMG_FMT_I420 || image_.fmt == VPX_IMG_FMT_I42016,
          WP2_STATUS_UNSUPPORTED_FEATURE);
      WP2_CHECK_OK(w_shift == 1 && h_shift == 1,
                   WP2_STATUS_UNSUPPORTED_FEATURE);
    }

    // TODO(yguyon): Check if there is a way to know if the file was encoded
    //               with yuv444p12be or yuv444p12le (big/little endian)
    WP2_CHECK_OK((image_.w == (unsigned int)info_.pic_w &&
                  image_.h == (unsigned int)info_.pic_h),
                 WP2_STATUS_UNSUPPORTED_FEATURE);

    // From YCbCr to RGB 8 bits.
    const uint32_t shift = info_.bit_depth - CCSPImageReader::kMinBitDepth;
    // See vpx_color_range_t.
    int32_t offset[] = {0, LeftShift(-128, shift), LeftShift(-128, shift)};

    // YCbCr->RGB matrix taken from http://www.mir.com/DMG/ycbcr.html (x 1<<12)
    if (image_.range == VPX_CR_STUDIO_RANGE) {
      offset[kYChannel] = LeftShift(-16, shift);
      // Multiply Y by 255/(235-16) and UV by 255/(240-16).
      *ccsp_to_rgb_ = CSPMtx({4769, 0, 6537, 4769, -1605, -3330, 4769, 8263, 0},
                             12 + shift);
    } else {
      assert(image_.range == VPX_CR_FULL_RANGE);

      *ccsp_to_rgb_ = CSPMtx({4096, 0, 5743, 4096, -1410, -2925, 4096, 7258, 0},
                             12 + shift);
    }

    // TODO(yguyon): Handle alpha if y4m does
    constexpr uint32_t kToVPX[] = {VPX_PLANE_Y, VPX_PLANE_U, VPX_PLANE_V};
    for (Channel channel : {kYChannel, kUChannel, kVChannel}) {
      Plane16* const plane = &ccsp_buffer_->GetChannel(channel);
      if (channel == kUChannel || channel == kVChannel) {
        // Subsampled dimensions, rounded up.
        WP2_CHECK_STATUS(
            plane->Resize((image_.w + ((1u << w_shift) - 1u)) >> w_shift,
                          (image_.h + ((1u << h_shift) - 1u)) >> h_shift));
      } else {
        WP2_CHECK_STATUS(plane->Resize(image_.w, image_.h));
      }

      const uint8_t* src_row = image_.planes[kToVPX[channel]];
      int16_t* dst_row = plane->Row(0);
      for (uint32_t y = 0; y < plane->h_; ++y) {
        const int16_t* const src_row16b =
            reinterpret_cast<const int16_t *>(src_row);
        for (uint32_t x = 0; x < plane->w_; ++x) {
          if (info_.bit_depth == 8) {
            dst_row[x] = (int16_t)(src_row[x] + offset[channel]);
          } else {
            dst_row[x] = (int16_t)(src_row16b[x] + offset[channel]);
          }
        }
        src_row += image_.stride[kToVPX[channel]];
        dst_row += plane->Step();
      }
    }

    // TODO(yguyon): Handle metadata if y4m does
    // TODO(yguyon): Handle several frames/correct durations
    *is_last = true;
    *duration_ms = ImageReader::kInfiniteDuration;
    return WP2_STATUS_OK;
  }

 protected:
  // Pointers to user structs.
  YUVPlane* const ccsp_buffer_;
  CSPMtx* const ccsp_to_rgb_;
  Metadata* const metadata_;

  ExternalDataSource data_source_;
  y4m_input info_;
  vpx_image_t image_;
  bool opened_;
};

void CCSPImageReader::SetImplY4M(YUVPlane *const ccsp_buffer,
                                 CSPMtx *const ccsp_to_rgb,
                                 Metadata *const metadata, LogLevel log_level,
                                 size_t max_num_pixels) {
  impl_.reset(new (WP2Allocable::nothrow) ImageReaderY4M(
      data_.bytes, data_.size, ccsp_buffer, ccsp_to_rgb, metadata, log_level,
      max_num_pixels));
  if (impl_ == nullptr) status_ = WP2_STATUS_OUT_OF_MEMORY;
}

}  // namespace WP2

//------------------------------------------------------------------------------
