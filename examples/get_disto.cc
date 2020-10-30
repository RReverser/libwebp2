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
// Simple tool to load two webp/png/jpg/tiff files and compute PSNR/SSIM.
// This is mostly a wrapper around ArgbBuffer::GetDistortion().
//
// Author: Skal (pascal.massimino@gmail.com)

#include <cassert>
#include <cstdio>
#include <cstring>

#include "examples/example_utils.h"
#include "imageio/image_dec.h"
#include "imageio/image_enc.h"
#include "imageio/imageio_util.h"
#include "src/utils/plane.h"
#include "src/wp2/encode.h"

using WP2::ArgbBuffer;
using WP2::YUVPlane;

namespace {

void RescalePlane(ArgbBuffer* const pic1, int channel, uint32_t max) {
  const uint32_t factor = (max > 0) ? (255u << 16) / max : 0;
  for (uint32_t y = 0; y < pic1->height; ++y) {
    uint8_t* const ptr = (uint8_t*)pic1->GetRow(y) + channel;
    for (uint32_t x = 0; x < pic1->width; ++x) {
      ptr[4 * x] = (ptr[4 * x] * factor + (1 << 15)) >> 16;
    }
  }
}

// Return the max absolute difference (and replace pic1 with the max map)
uint32_t DiffScaleChannel(ArgbBuffer* const pic1, const ArgbBuffer* const pic2,
                          int channel, bool do_scaling) {
  uint32_t max = 0;
  for (uint32_t y = 0; y < pic1->height; ++y) {
    uint8_t* const ptr1 = (uint8_t*)pic1->GetRow(y) + channel;
    const uint8_t* const ptr2 = (const uint8_t*)pic2->GetRow(y) + channel;
    for (uint32_t x = 0; x < pic1->width; ++x) {
      const uint32_t diff = abs(ptr1[4 * x] - ptr2[4 * x]);
      if (diff > max) max = diff;
      ptr1[4 * x] = diff;
    }
  }

  if (do_scaling) RescalePlane(pic1, channel, max);
  return max;
}

// Convert an argb picture to luminance.
void ConvertToGray(ArgbBuffer* const pic) {
  assert(pic != nullptr);
  for (uint32_t y = 0; y < pic->height; ++y) {
    uint8_t* const row = (uint8_t*)pic->GetRow(y);
    for (uint32_t x = 0; x < pic->width; ++x) {
      const uint32_t r = row[4 * x + 1];
      const uint32_t g = row[4 * x + 2];
      const uint32_t b = row[4 * x + 3];
      // We use BT.709 for converting to luminance.
      const uint8_t Y = (uint8_t)(0.2126 * r + 0.7152 * g + 0.0722 * b + .5);
      row[4 * x + 1] = row[4 * x + 2] = row[4 * x + 3] = Y;
    }
  }
}

void Help() {
  ProgramOptions opt;
  opt.Add("Usage:");
  opt.Add("  get_disto [options] compressed_file orig_file");
  opt.Add("");
  opt.Add("Options:");
  opt.AddMetricOptions();
  opt.Add("-exact", "expect a perfect match");
  opt.Add("-min <float>", "check overall distortion is larger than value");
  opt.Add("-o <file>", "save the diff map as a WP2 lossless file");
  opt.Add("-scale", "scale the difference map to fit [0-255] range");
  opt.Add("-gray", "convert difference map to gray in the end");
  opt.Add("-half", "2x downscale before comparison");
  opt.Add("-crop <x> <y> <w> <h>",
          "Crop picture with given rectangle. Alternatively, one can use "
          "-crop_source to only crop the reference file (second argument).");
  opt.Add("-gray", "use grayscale for difference map (-scale)");
  opt.Add("-ref <string>", "specify reference file for bpp calculation");
  opt.Add("-prefix <string>", "add a prefix label");
  opt.Add("-suffix <string>", "add a suffix label");
  opt.Add("");
  opt.AddSystemOptionSection();
  opt.Add(
      " Also handles PNG, PPM, JPG, WEBP and TIFF files, in addition to WP2.");
  opt.Add("");
  opt.Print();
}

}    // namespace

int main(int argc, const char *argv[]) {
  EXIT_IF_FALSE(WP2CheckVersion(), "Error! Library version mismatch!");

  WP2::MetricType type = WP2::PSNR;
  bool scale = false;
  bool use_gray = false;
  bool exact = false;
  bool downscale = false;
  float min_disto = 0.;
  const char* prefix = nullptr;
  const char* suffix = nullptr;
  const char* path1 = nullptr;
  const char* path2 = nullptr;
  const char* ref_file = nullptr;
  const char* output = nullptr;
  enum { NO_CROP, CROP_BOTH, CROP_SOURCE } crop = NO_CROP;
  WP2::Rectangle crop_area;

  for (int c = 1; c < argc; ++c) {
    bool parse_error = false;
    const char* const arg = argv[c];
    if (!strcmp(arg, "-exact")) {
      exact = true;
    } else if (!strcmp(arg, "-scale")) {
      scale = true;
    } else if (!strcmp(arg, "-gray")) {
      use_gray = true;
    } else if (!strcmp(arg, "-min") && c + 1 < argc) {
      min_disto = ExUtilGetFloat(argv[++c], &parse_error);
    } else if (!strcmp(arg, "-half")) {
      downscale = true;
    } else if (!strcmp(arg, "-ref") && c + 1 < argc) {
      ref_file = argv[++c];
    } else if (!strcmp(arg, "-prefix") && c + 1 < argc) {
      prefix = argv[++c];
    } else if (!strcmp(arg, "-suffix") && c + 1 < argc) {
      suffix = argv[++c];
    } else if ((!strcmp(arg, "-crop") || !strcmp(arg, "-crop_source")) &&
               c < argc - 4) {
      crop = !strcmp(arg, "-crop") ? CROP_BOTH : CROP_SOURCE;
      crop_area.x = ExUtilGetUInt(argv[++c], &parse_error);
      crop_area.y = ExUtilGetUInt(argv[++c], &parse_error);
      crop_area.width = ExUtilGetUInt(argv[++c], &parse_error);
      crop_area.height = ExUtilGetUInt(argv[++c], &parse_error);
    } else if (!strcmp(arg, "-h")) {
      Help();
      return 0;
    } else if (!strcmp(arg, "-o")) {
      EXIT_IF_FALSE(c + 1 < argc, "missing file name after %s option.", arg);
      output = argv[++c];
    } else if (argv[c][0] == '-') {
      bool must_stop;
      if (ProgramOptions::ParseSystemOptions(argv[c], &must_stop)) {
        if (must_stop) return 0;
      } else if (ProgramOptions::ParseMetricOptions(argv[c], &type)) {
      } else {
        printf("Unknown option '%s'\n", argv[c]);
        parse_error = true;
      }
    } else if (path1 == nullptr) {
      path1 = arg;
    } else {
      path2 = arg;
    }
    if (parse_error) {
      Help();
      return 1;
    }
  }

  if (path1 == nullptr || path2 == nullptr) {
    fprintf(stderr, "Error: missing arguments.\n");
    Help();
    return 0;
  }

  size_t file_size1, file_size2;
  ArgbBuffer pic1_rgb, pic2_rgb;
  EXIT_IF_ERROR(WP2::ReadImage(path1, &pic1_rgb, &file_size1),
                "Can't decode input file '%s'.", path1);
  EXIT_IF_ERROR(WP2::ReadImage(path2, &pic2_rgb, &file_size2),
                "Can't decode input file '%s'.", path2);

  // save original dimensions (before crop) for bpp calculation.
  const uint32_t width = pic1_rgb.width, height = pic1_rgb.height;
  if (crop != NO_CROP) {
    if (crop == CROP_BOTH) {
      EXIT_IF_ERROR(pic1_rgb.SetView(pic1_rgb, crop_area),
                    "Cropping operation failed! Wrong parameters?");
    }
    EXIT_IF_ERROR(pic2_rgb.SetView(pic2_rgb, crop_area),
                  "Cropping operation failed! Wrong parameters?");
  }
  if (downscale) {
    pic1_rgb.SimpleHalfDownsample();
    pic2_rgb.SimpleHalfDownsample();
  }

  size_t ref_size = file_size1;
  if (ref_file != nullptr) {
    EXIT_IF_ERROR(WP2::IoUtilFileSize(ref_file, &ref_size),
                  "Could not read size of reference file");
  }
  const float bpp = 8.f * ref_size / (width * height);

  float disto[5];
  EXIT_IF_ERROR(
      pic1_rgb.GetDistortionBlackOrWhiteBackground(pic2_rgb, type, disto),
      "Error while computing the distortion.");

  if (exact) min_disto = 99.f;
  if (min_disto > 0.) {   // simply check min disto and exit
    const bool ok = (disto[4] >= min_disto);
    fprintf(stderr, "[%s] Check: %f > %f ? %s\n",
            path1, disto[4], min_disto, ok ? "OK" : "FAIL");
    return ok ? 0 : 1;
  }
  if (prefix != nullptr) printf("%s ", prefix);
  printf("%u %.2f    %.2f %.2f %.2f %.2f [ %.2f bpp %s]",
         (unsigned int)ref_size,
         disto[4], disto[0], disto[1], disto[2], disto[3],
         bpp, (crop != NO_CROP) ? "(cropped)" : "");
  if (suffix != nullptr) printf(" %s", suffix);
  printf("\n");

  if (output != nullptr) {
    fprintf(stderr, "max differences per channel: ");
    for (int n = 1; n < 4; ++n) {    // skip the alpha channel
      const int range = DiffScaleChannel(&pic1_rgb, &pic2_rgb, n, scale);
      if (range < 0) fprintf(stderr, "\nError computing diff map\n");
      fprintf(stderr, "[%d]", range);
    }
    fprintf(stderr, "\n");
    if (use_gray) ConvertToGray(&pic1_rgb);
    EXIT_IF_ERROR(WP2::SaveImage(pic1_rgb, output, /*overwrite=*/true),
                  "Error during saving [%s].", output);
  }
  return 0;
}
