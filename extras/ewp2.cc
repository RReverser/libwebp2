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
//  To and from Y4M conversion command line tool.
//
// Author: Yannis Guyon (yguyon@google.com)

#include <array>
#include <cstdio>
#include <cstring>
#include <string>

#include "examples/example_utils.h"
#include "extras/ccsp_imageio.h"
#include "extras/extras.h"
#include "imageio/image_dec.h"
#include "imageio/image_enc.h"
#include "src/wp2/decode.h"
#include "src/wp2/encode.h"

using FileFormat = WP2::FileFormat;

//------------------------------------------------------------------------------

namespace {

void Help() {
  ProgramOptions opt;
  opt.Add("Usage:");
  opt.Add("  ewp2 in_file [-o out_file] [options]");
  opt.Add("");
  opt.Add("Options valid if 'out_file' has the .wp2 or .webp extension:");
  opt.Add("-q <float>",
          SPrintf("quality factor (0:small..100:big), default=%1.1f",
                  WP2::EncoderConfig::kDefault.quality));
  opt.Add("-speed <int>",
          SPrintf("compression effort (0:fast..9:slower/better), default=%d",
                  WP2::EncoderConfig::kDefault.speed));
  opt.Add("-csp <int>", "color space (wp2 only)");
  opt.Add("");
  opt.Add(
      "Options valid if 'out_file' has the .y4m extension "
      "(4:4:4 8b by default):");
  opt.Add("-depth <int>",
          SPrintf("output bit depth per plane per sample (%u..%u)",
                  WP2::CCSPImageReader::kMinBitDepth,
                  WP2::CCSPImageReader::kMaxBitDepth));
  opt.Add("-444", "4:4:4 output samples (full chroma resolution)");
  opt.Add("-420", "4:2:0 output samples (half chroma resolution)");
  opt.Add("");
  opt.Add("Options valid if 'in_file' and/or 'out_file' are 4:2:0:");
  opt.Add("-sampling <filter>",
          "algorithm to use for up- and/or chroma downsampling: "
          "smooth (default), nearest");
  opt.Add("");
  opt.Add("Options valid if 'out_file' is specified:");
  opt.Add("-[no]metadata", "include [or exclude] input metadata in output");
  opt.Add("-[psnr|ssim] <file path> <float>",
          "error returned if output (or input if no given output) has a "
          "distortion worse than specified compared to a reference file");
  opt.Add("");
  opt.Add("Other options:");
  opt.Add("-v", "verbose");
  opt.Add("-size", "print size and bits-per-pixel");
  opt.Add("-rm [A|R|G|B|Y|U|V]",
          "set specified input channel(s) values to opaque/black/grey");
  opt.AddSystemOptionSection();
  opt.Print();
}

enum ChannelId { kA, kR, kG, kB, kY, kU, kV, kNumChans };

void PrintSize(const std::string& file_path, size_t size, uint32_t num_pixels,
               bool verbose, bool print_disto) {
  const double bpp = size * 8. / num_pixels;
  if (verbose) {
    printf("%s %zu bytes (%.3f bpp)\n", file_path.c_str(), size, bpp);
  } else {
    printf("%s %zu %.3f bpp%s", WP2GetFileName(file_path).c_str(), size, bpp,
           print_disto ? "" : "\n");
  }
}

}  // namespace

//------------------------------------------------------------------------------

int main(int argc, char* argv[]) {
  EXIT_IF_FALSE(WP2CheckVersion(), "Error! Library version mismatch!");
  if (argc < 2) {
    fprintf(stderr, "Error: Missing arguments.\n");
    Help();
    return 0;
  }

  std::string in_path;
  std::string out_path;
  bool keep_metadata = true;

  bool expecting_lossy_enc = false;
  bool expecting_wp2_enc = false;
  WP2::EncoderConfig wp2_enc_config;

  bool expecting_ccsp_enc = false;
  uint32_t ccsp_enc_bit_depth = 8;  // 4:4:4 8b by default
  uint32_t ccsp_enc_420 = false;
  const WP2::SamplingTaps* downsampling_filter = &WP2::SamplingTaps::kDownSharp;
  const WP2::SamplingTaps* upsampling_filter = &WP2::SamplingTaps::kUpSmooth;

  std::string ref_path;
  WP2::MetricType ref_metric = WP2::NUM_METRIC_TYPES;
  float ref_distortion = 99.f;

  bool verbose = false, print_size = false;
  bool rm_chan[kNumChans] = {};

  for (int c = 1; c < argc; ++c) {
    bool parse_error = false;
    if (!strcmp(argv[c], "-h") || !strcmp(argv[c], "-help")) {
      Help();
      return 0;
    } else if (!strcmp(argv[c], "-o") && c + 1 < argc) {
      EXIT_IF_FALSE(out_path.empty(), "Error! Output was already specified.");
      out_path = argv[++c];
    } else if (!strcmp(argv[c], "-q") && c + 1 < argc) {
      expecting_lossy_enc = true;
      wp2_enc_config.quality = ExUtilGetFloat(argv[++c], &parse_error);
    } else if (!strcmp(argv[c], "-speed") && c + 1 < argc) {
      expecting_lossy_enc = true;
      wp2_enc_config.speed = ExUtilGetInt(argv[++c], &parse_error);
    } else if (!strcmp(argv[c], "-csp") && c + 1 < argc) {
      expecting_lossy_enc = true;
      expecting_wp2_enc = true;
      wp2_enc_config.csp_type = (WP2::Csp)ExUtilGetInt(argv[++c], &parse_error);
    } else if (!strcmp(argv[c], "-depth") && c + 1 < argc) {
      expecting_ccsp_enc = true;
      ccsp_enc_bit_depth = ExUtilGetInt(argv[++c], &parse_error);
      EXIT_IF_FALSE(
          ccsp_enc_bit_depth >= WP2::CCSPImageReader::kMinBitDepth &&
              ccsp_enc_bit_depth <= WP2::CCSPImageReader::kMaxBitDepth,
          "Error! Bit depth is outside range.");
    } else if (!strcmp(argv[c], "-444") || !strcmp(argv[c], "-420")) {
      expecting_ccsp_enc = true;
      ccsp_enc_420 = !strcmp(argv[c], "-420");
    } else if (!strcmp(argv[c], "-sampling") && c + 1 < argc) {
      ++c;
      if (!strcmp(argv[c], "smooth")) {
        downsampling_filter = &WP2::SamplingTaps::kDownSharp;
        upsampling_filter = &WP2::SamplingTaps::kUpSmooth;
      } else if (!strcmp(argv[c], "nearest")) {
        downsampling_filter = &WP2::SamplingTaps::kDownAvg;
        upsampling_filter = &WP2::SamplingTaps::kUpNearest;
      } else {
        EXIT_IF_FALSE(false, "Error! Unsupported sampling '%s'", argv[c]);
      }
    } else if ((!strcmp(argv[c], "-psnr") || !strcmp(argv[c], "-ssim")) &&
               c + 2 < argc) {
      ref_metric = !strcmp(argv[c], "-psnr") ? WP2::PSNR : WP2::SSIM;
      ref_path = argv[++c];
      ref_distortion = ExUtilGetFloat(argv[++c], &parse_error);
    } else if (!strcmp(argv[c], "-v")) {
      verbose = true;
    } else if (!strcmp(argv[c], "-size")) {
      print_size = true;
    } else if (!strcmp(argv[c], "-rm") && c + 1 < argc) {
      constexpr char kChanChars[] = {'A', 'R', 'G', 'B', 'Y', 'U', 'V'};
      STATIC_ASSERT_ARRAY_SIZE(kChanChars, kNumChans);
      for (char character : std::string(argv[++c])) {
        const char* const pos = std::strchr(kChanChars, character);
        EXIT_IF_FALSE(pos != nullptr, "Error! Unsupported -rm '%c'", character);
        rm_chan[pos - kChanChars] = true;
      }
    } else if (!strcmp(argv[c], "-metadata")) {
      keep_metadata = true;
    } else if (!strcmp(argv[c], "-nometadata")) {
      keep_metadata = false;
    } else if (!strcmp(argv[c], "--") && c + 1 < argc) {
      EXIT_IF_FALSE(in_path.empty(), "Error! Input was already specified.");
      in_path = argv[++c];
    } else if (argv[c][0] == '-') {
      bool must_stop;
      EXIT_IF_FALSE(ProgramOptions::ParseSystemOptions(argv[c], &must_stop),
                    "Error! Unknown option '%s'", argv[c]);
      if (must_stop) return 0;
    } else {
      EXIT_IF_FALSE(in_path.empty(), "Error! Input was already specified.");
      in_path = argv[c];
    }
    if (parse_error) return 1;
  }

  const bool print_disto = !ref_path.empty();

  EXIT_IF_FALSE(!in_path.empty(), "Error! No input specified.");
  EXIT_IF_FALSE(!expecting_lossy_enc || !expecting_ccsp_enc,
                "Error! Output format mixup.");
  EXIT_IF_FALSE(
      (!expecting_lossy_enc && !expecting_ccsp_enc) || !out_path.empty(),
      "Error! No output specified.");
  EXIT_IF_FALSE(wp2_enc_config.IsValid(), "Error! Invalid configuration.");

  WP2::YUVPlane ccsp_buffer;
  WP2::CSPMtx ccsp_to_rgb = {};
  WP2::Metadata metadata;

  WP2::Data in_data;
  EXIT_IF_ERROR(WP2::IoUtilReadFile(in_path.c_str(), &in_data),
                "Error! Failed to read file '%s'.", in_path.c_str());
  const FileFormat in_format =
      WP2::GuessImageFormat(in_data.bytes, in_data.size);

  FileFormat out_format = WP2::GetFormatFromExtension(out_path.c_str());
  if (out_format == FileFormat::Y4M_420 || out_format == FileFormat::Y4M_444) {
    out_format = (ccsp_enc_420 ? FileFormat::Y4M_420 : FileFormat::Y4M_444);
  } else {
    EXIT_IF_FALSE(!ccsp_enc_420,
                  "Error! Output format %s does not support 4:2:0.",
                  WP2::GetFileFormatStr(out_format));
  }

  if (print_size) {
    EXIT_IF_FALSE(out_path != "-", "Error! Printing stdout size is unhandled.");
  }
  size_t output_file_size = 0;
  size_t* const file_size_ptr = print_size ? &output_file_size : nullptr;

  if (verbose) printf("Decoding '%s'\n", in_path.c_str());

  if (in_format == FileFormat::WP2 && (out_format == FileFormat::Y4M_444 ||
                                       out_format == FileFormat::Y4M_420)) {
    // Bypass conversion to and from RGB by decoding directly into YCbCr.
    WP2::BitstreamFeatures features;
    EXIT_IF_ERROR(features.Read(in_data.bytes, in_data.size),
                  "Error! Failed to decode the wp2 image header.");
    const uint32_t rgb_to_ycbcr_shift =
        WP2::kRGBToYCbCrShift - (ccsp_enc_bit_depth - 8);

    // Allocate a padded buffer to avoid doing it internally.
    // TODO(yguyon): Check if 'has_alpha' is needed beforehand (decode GParams).
    EXIT_IF_ERROR(
        ccsp_buffer.Resize(features.width, features.height,
                           /*pad=*/WP2::kPredWidth, /*has_alpha=*/true),
        "Error! Failed to allocate buffer.");

    EXIT_IF_ERROR(
        WP2::Decode(
            in_data.bytes, in_data.size, WP2::kRGBToYCbCrMatrix,
            rgb_to_ycbcr_shift,
            ccsp_buffer.Y.Row(0), ccsp_buffer.Y.Step(), ccsp_buffer.Y.Size(),
            ccsp_buffer.U.Row(0), ccsp_buffer.U.Step(), ccsp_buffer.U.Size(),
            ccsp_buffer.V.Row(0), ccsp_buffer.V.Step(), ccsp_buffer.V.Size(),
            ccsp_buffer.A.Row(0), ccsp_buffer.A.Step(), ccsp_buffer.A.Size(),
            keep_metadata ? &metadata : nullptr),
        "Error! Failed to decode the wp2 image.");

    EXIT_IF_FALSE(
        !rm_chan[kR] && !rm_chan[kG] && !rm_chan[kB],
        "Error! Input wp2 is decoded straight to YUV so -rm RGB is invalid.");
    if (rm_chan[kA]) ccsp_buffer.A.Clear();
    if (rm_chan[kY]) ccsp_buffer.Y.Fill(0);
    if (rm_chan[kU]) ccsp_buffer.U.Fill(0);
    if (rm_chan[kV]) ccsp_buffer.V.Fill(0);

    // Remove padding.
    WP2_ASSERT_STATUS(ccsp_buffer.SetView(
        ccsp_buffer, {0, 0, features.width, features.height}));

    if (WP2::IsChromaSubsampled(out_format)) {
      EXIT_IF_ERROR(ccsp_buffer.Downsample(*downsampling_filter),
                    "Error! Failed to downsample to 4:2:0.");
    }

    if (!out_path.empty()) {
      if (verbose) printf("Saving to '%s'\n", out_path.c_str());
      FILE* const fout = (out_path == "-") ? WP2::IoUtilSetBinaryMode(stdout)
                                           : std::fopen(out_path.c_str(), "wb");
      EXIT_IF_FALSE(fout != nullptr, "Error! Failed to write file '%s'.",
                    out_path.c_str());
      const WP2Status status =
          WP2::WriteY4M(ccsp_buffer, ccsp_enc_bit_depth, fout);
      if (file_size_ptr != nullptr) *file_size_ptr = std::ftell(fout);
      if (fout != stdout) fclose(fout);  // Close file even if error.
      EXIT_IF_ERROR(status, "Error! Failed to write image '%s'.",
                    out_path.c_str());
      if (!metadata.IsEmpty()) fprintf(stderr, "Warning: Ignored metadata.\n");
    }
  } else {  // General case: any format -> RGB -> any format.
    EXIT_IF_ERROR(
        WP2::ReadImage(in_data.bytes, in_data.size, &ccsp_buffer, &ccsp_to_rgb,
                       keep_metadata ? &metadata : nullptr),
        "Error! Failed to read input '%s'.", in_path.c_str());

    if (WP2::IsCustomColorSpace(in_format)) {
      EXIT_IF_FALSE(!rm_chan[kR] && !rm_chan[kG] && !rm_chan[kB],
                    "Error! Input is considered YUV so -rm RGB is invalid.");
    } else {
      EXIT_IF_FALSE(!rm_chan[kY] && !rm_chan[kU] && !rm_chan[kV],
                    "Error! Input is considered RGB so -rm YUV is invalid.");
    }
    if (rm_chan[kA]) ccsp_buffer.A.Clear();
    if (rm_chan[kY] || rm_chan[kR]) ccsp_buffer.Y.Fill(0);
    if (rm_chan[kU] || rm_chan[kG]) ccsp_buffer.U.Fill(0);
    if (rm_chan[kV] || rm_chan[kB]) ccsp_buffer.V.Fill(0);

    if (!out_path.empty()) {
      if (!WP2::IsChromaSubsampled(out_format) && ccsp_buffer.IsDownsampled()) {
        EXIT_IF_ERROR(ccsp_buffer.Upsample(*upsampling_filter),
                      "Error! Failed to upsample to 4:4:4.");
      }

      if (expecting_ccsp_enc) {
        EXIT_IF_FALSE(
            WP2::IsCustomColorSpace(out_format),
            "Error! Expected custom color space output instead of %s format.",
            WP2::GetFileFormatStr(out_format));
      } else if (expecting_lossy_enc) {
        EXIT_IF_FALSE(
            out_format == FileFormat::WP2 ||
                (!expecting_wp2_enc && out_format == FileFormat::WEBP),
            "Error! Expected wp2%s output instead of %s format.",
            expecting_wp2_enc ? "" : " or WebP",
            WP2::GetFileFormatStr(out_format));
      }

      if (verbose) printf("Saving to '%s'\n", out_path.c_str());
      if (out_format == FileFormat::WP2) {
        WP2::MemoryWriter writer;
        EXIT_IF_ERROR(
            WP2::Encode(ccsp_buffer.GetWidth(), ccsp_buffer.GetHeight(),
                        ccsp_buffer.Y.Row(0), ccsp_buffer.Y.Step(),
                        ccsp_buffer.U.Row(0), ccsp_buffer.U.Step(),
                        ccsp_buffer.V.Row(0), ccsp_buffer.V.Step(),
                        ccsp_buffer.HasAlpha() ? ccsp_buffer.A.Row(0) : nullptr,
                        ccsp_buffer.A.Step(), ccsp_to_rgb.mtx(),
                        ccsp_to_rgb.shift, &writer, wp2_enc_config, metadata),
            "Error! Failed to encode as wp2.");

        EXIT_IF_ERROR(
            WP2::IoUtilWriteFile(writer.mem_, writer.size_, out_path.c_str(),
                                 /*overwrite=*/true, file_size_ptr),
            "Error! Failed to write image '%s'.", out_path.c_str());
      } else if (out_format == FileFormat::WEBP) {
        WP2::ArgbBuffer argb_buffer;
        EXIT_IF_ERROR(ccsp_buffer.Export(ccsp_to_rgb, /*resize_if_needed=*/true,
                                         &argb_buffer, upsampling_filter),
                      "Error! Failed to convert to Argb.");
        if (!metadata.IsEmpty()) fprintf(stderr, "Warning: Ignored metadata\n");

        WP2::MemoryWriter writer;
        EXIT_IF_ERROR(WP2::CompressWebP(argb_buffer, wp2_enc_config.quality,
                                        wp2_enc_config.speed, &writer),
                      "Error! Failed to encode as WebP.");

        EXIT_IF_ERROR(
            WP2::IoUtilWriteFile(writer.mem_, writer.size_, out_path.c_str(),
                                 /*overwrite=*/true, file_size_ptr),
            "Error! Failed to write image '%s'.", out_path.c_str());
      } else {
        EXIT_IF_ERROR(
            WP2::SaveImage(ccsp_buffer, ccsp_to_rgb, ccsp_enc_bit_depth,
                           out_path.c_str(), /*overwrite=*/true, out_format,
                           metadata, *downsampling_filter, file_size_ptr),
            "Error! Failed to save image '%s'.", out_path.c_str());
      }
    }
  }

  if (out_path.empty() && verbose) printf("No specified output.\n");

  if (print_size) {
    const uint32_t num_px = ccsp_buffer.GetWidth() * ccsp_buffer.GetHeight();
    if (out_path.empty()) {
      PrintSize(in_path, in_data.size, num_px, verbose, print_disto);
    } else {
      PrintSize(out_path, output_file_size, num_px, verbose, print_disto);
    }
  }

  if (print_disto) {
    const std::string& cmp_path = out_path.empty() ? in_path : out_path;

    WP2::YUVPlane cmp, ref;
    WP2::CSPMtx cmp_matrix = {}, ref_matrix = {};

    if (cmp_path == in_path) {  // Compare 'ref' with input.
      using std::swap;
      swap(ccsp_buffer, cmp);
      cmp_matrix = ccsp_to_rgb;
    } else {  // Compare 'ref' with output.
      EXIT_IF_ERROR(
          WP2::ReadImage(cmp_path.c_str(), &cmp, &cmp_matrix),
          "Error! Failed to read '%s'.", cmp_path.c_str());
    }
    EXIT_IF_ERROR(
        WP2::ReadImage(ref_path.c_str(), &ref, &ref_matrix),
        "Error! Failed to read reference '%s'.", ref_path.c_str());

    EXIT_IF_FALSE(
        cmp_matrix == ref_matrix,
        "Error! Can only compute distortion for the same color space.");
    uint32_t cmp_bit_depth, ref_bit_depth;
    EXIT_IF_ERROR(WP2::ReadBitDepth(cmp_path.c_str(), &cmp_bit_depth),
                  "Error! Failed to read bit depth of '%s'.", cmp_path.c_str());
    EXIT_IF_ERROR(WP2::ReadBitDepth(ref_path.c_str(), &ref_bit_depth),
                  "Error! Failed to read bit depth of '%s'.", ref_path.c_str());
    // Otherwise the matrix or the shift should have been different.
    assert(cmp_bit_depth == ref_bit_depth);

    if (cmp.HasAlpha() != ref.HasAlpha()) {
      // Set image with missing alpha plane to opaque.
      if (!cmp.HasAlpha()) {
        EXIT_IF_ERROR(cmp.A.Resize(cmp.GetWidth(), cmp.GetHeight()),
                      "Error! Allocation failure.");
        cmp.A.Fill(WP2::kAlphaMax);
      } else {
        EXIT_IF_ERROR(ref.A.Resize(ref.GetWidth(), ref.GetHeight()),
                      "Error! Allocation failure.");
        ref.A.Fill(WP2::kAlphaMax);
      }
    }

    float distortion[5];
    EXIT_IF_ERROR(cmp.GetDistortion(ref, ref_bit_depth, ref_metric, distortion),
                  "Error! Failed to compute distortion.");

    if (verbose) {
      printf("Distortion %f (alpha %f, %f, %f, %f) between '%s' and '%s'\n",
             distortion[4], distortion[0], distortion[1], distortion[2],
             distortion[3], cmp_path.c_str(), ref_path.c_str());
    } else {
      printf("%s distortion %f   %f %f %f %f\n",
             print_size ? "   " : WP2GetFileName(cmp_path).c_str(),
             distortion[4], distortion[0], distortion[1], distortion[2],
             distortion[3]);
    }
    EXIT_IF_FALSE(
        distortion[4] >= ref_distortion,
        "Error! Distortion %f from %s to %s is worse than the required %f.",
        distortion[4], WP2GetFileName(ref_path).c_str(),
        WP2GetFileName(cmp_path).c_str(), ref_distortion);
  }

  return 0;
}

//------------------------------------------------------------------------------
