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
//  simple command line tools to create/manipulate/optimize preview bits
//
// Author: Skal (pascal.massimino@gmail.com)

#include <cstdio>
#include <cstring>
#include <string>
#include <vector>

#include "extras/extras.h"
#include "examples/example_utils.h"
#include "imageio/file_format.h"
#include "imageio/imageio_util.h"
#include "imageio/image_dec.h"
#include "imageio/image_enc.h"
#include "src/common/preview/preview.h"   // TODO(skal): make better include
#include "src/enc/preview/preview_enc.h"
#include "src/dec/preview/preview_dec.h"
#include "src/wp2/encode.h"

#define CHECK_OK(COND, STR...) do {   \
  if (!(COND)) {                      \
    fprintf(stderr, STR);             \
    return 1;                         \
  }                                   \
} while (false)
#define CHECK_STATUS(CALL, STR...) CHECK_OK(((CALL) == WP2_STATUS_OK), STR)

static int verbose = 1;
static bool check = false;
#define VERBOSE(STR...) if (verbose > 0) printf(STR)  // NOLINT (string literal)

namespace {
bool operator==(const WP2::AYCoCg19b& a, const WP2::AYCoCg19b& b) {
  return (a.y == b.y) && (a.co == b.co) && (a.cg == b.cg) && (a.a == b.a);
}
bool operator==(const WP2::VertexIndexedColor& a,
                const WP2::VertexIndexedColor& b) {
  return (a.x == b.x) && (a.y == b.y) && (a.color_index == b.color_index);
}
template<typename T> bool SameArray(const T a[], const T b[], size_t size) {
  for (size_t i = 0; i < size; ++i) {
    if (!(a[i] == b[i])) return false;
  }
  return true;
}

bool operator==(const WP2::PreviewData& A, const WP2::PreviewData& B) {
  return (A.grid_width_ == B.grid_width_) &&
         (A.grid_height_ == B.grid_height_) &&
         (A.palette_.size() == B.palette_.size()) &&
         (A.vertices_.size() == B.vertices_.size()) &&
         SameArray(A.corners_, B.corners_, 4) &&
         SameArray(A.palette_.data(), B.palette_.data(), A.palette_.size()) &&
         SameArray(A.vertices_.data(), B.vertices_.data(), A.vertices_.size());
}

}  // namespace

static void Help() {
  ProgramOptions opt;
  opt.Add("Usage:");
  opt.Add("   preview -convert [-dry] text_files... ");
  opt.Add("   preview image [bits|text] [options... ]"
          "[-o output] [-d output]");
  opt.Add("");
  opt.Add("PreviewConfiguration Options:");
  opt.Add(" -m [edge,rand,diff]", "analysis method");
  opt.Add(" -n <int>", "starting number of vertices");
  opt.Add(" -nc <int>", "starting number of colors");
  opt.Add(" -density <float>", "grid density");
  opt.Add(" -g <int>", "initial grid size");
  opt.Add(" -b <int>", "blur radius");
  opt.Add("");

  opt.Add("Optimization Options:");
  opt.Add(" -s <int>", "target size in bytes");
  opt.Add(" -i <int>",  "number of optim iterations");
  opt.Add(" -l <int>",  "size vs quality trade-off. "
          "0=favor quality, 100=favor size");
  opt.Add(" -p none", "reset all proba to 0");
  opt.Add(" -p <move> <int>",  "change the proba for the given <move>:");
  opt.Add("    . vmove: vertex move");
  opt.Add("    . vadd:  vertex insertion");
  opt.Add("    . vsub:  vertex removal");
  opt.Add("    . cmove: color move");
  opt.Add("    . cadd:  color insertion");
  opt.Add("    . csub:  color removal");
  opt.Add("    . cidx:  color index move");
  opt.Add("");

  opt.Add("Other Options:");
  opt.Add(" -txt", "use Text for preview data");
  opt.Add(" -info", "print full preview info");
  opt.Add(" -short", "print short preview info");
  opt.Add(" -w <int>", "output width for thumbnail");
  opt.Add(" -d <file_prefix>", "prefix for reconstructed thumbnail");
  opt.Add(" -o <file_prefix>", "prefix for output data");
  opt.Add(" -check", "perform round-trip check");
  opt.Add(" -noforce", "don't force file over-write");
  opt.Add(" -quiet", "reduced verbosity");
  opt.Add(" -verbose", "extra verbosity");

  opt.Print();
}

//------------------------------------------------------------------------------

static int ConvertTextFiles(int argc, const char* argv[]) {
  // pure txt->bits conversion mode
  uint32_t total_size = 0;
  uint32_t nb_files = 0;
  bool dry_run = false;
  for (int c = 2; c < argc; ++c) {
    if (!strcmp(argv[c], "-dry")) {
      dry_run = true;
      continue;
    } else if (!strcmp(argv[c], "-h")) {
      Help();
      return 0;
    }
    std::string in;
    CHECK_STATUS(WP2::IoUtilReadFile(argv[c], &in),
                 "bad input [%s]\n", argv[c]);
    WP2::PreviewData data;
    CHECK_STATUS(PreviewFromText(in, &data), "Error in text-to-preview\n");
    WP2::ANSEnc enc;
    CHECK_STATUS(data.Encode(&enc),
                 "Error encoding to binary! SHOULDN'T HAPPEN!!\n");
    if (!dry_run) {
      const std::string out_name = WP2RemoveFileExtension(argv[c]) + ".bits";
      CHECK_STATUS(WP2::IoUtilWriteFile(enc.Buffer(), enc.BufferSize(),
                                        out_name.c_str()),
                   "Couldn't save compressed bits [%s]!\n", out_name.c_str());
    }
    if (check) {
      WP2::PreviewData verif;
      CHECK_STATUS(verif.Decode(enc.Buffer(), enc.BufferSize()),
                   "Round-trip Decode() failed.\n");
      CHECK_STATUS(verif == data, "Round-trip verification failed!\n");
    }
    total_size += enc.BufferSize();
    ++nb_files;
  }
  VERBOSE("%u files, size: %u bytes\n", nb_files, total_size);
  return 0;
}

//------------------------------------------------------------------------------

static void ResetProbas(WP2::PreviewConfig* const config) {
  config->proba_vertex_move = 0;
  config->proba_vertex_add = 0;
  config->proba_vertex_sub = 0;
  config->proba_color_move = 0;
  config->proba_color_add = 0;
  config->proba_color_sub = 0;
  config->proba_color_index_move = 0;
}

int main(int argc, const char* argv[]) {
  if (argc < 2) {
    fprintf(stderr, "missing arguments!\n");
    Help();
    return 1;
  }
  if (argc > 1 && !strcmp(argv[1], "-convert")) {
    return ConvertTextFiles(argc, argv);
  }

  WP2::PreviewConfig config(75.f, 5);
  const char* img_file = nullptr;
  const char* bin_file = nullptr;
  const char* out_file = nullptr;
  const char* preview_file = nullptr;
  bool overwrite = true;
  uint32_t out_w = 256;   // width of decoded preview
  enum { TEXT, INFO, SHORT } out_format = SHORT;
  for (int c = 1; c < argc; ++c) {
    bool parse_error = false;
    if (!strcmp(argv[c], "-h")) {
      Help();
      return 0;
    } else if (!strcmp(argv[c], "-check")) {
      check = true;
    } else if (!strcmp(argv[c], "-noforce")) {
      overwrite = false;
    } else if (!strcmp(argv[c], "-quiet")) {
      verbose = 0;
      out_format = SHORT;
    } else if (!strcmp(argv[c], "-verbose")) {
      verbose = 2;
    } else if (!strcmp(argv[c], "-txt")) {
      out_format = TEXT;
    } else if (!strcmp(argv[c], "-info")) {
      out_format = INFO;
    } else if (!strcmp(argv[c], "-short")) {
      out_format = SHORT;
    } else if (!strcmp(argv[c], "-w") && c + 1 < argc) {
      out_w = std::max(ExUtilGetUInt(argv[++c], &parse_error), 32u);
    } else if (!strcmp(argv[c], "-o") && c + 1 < argc) {
      out_file = argv[++c];
    } else if (!strcmp(argv[c], "-d") && c + 1 < argc) {
      preview_file = argv[++c];
    } else if (!strcmp(argv[c], "-s") && c + 1 < argc) {
      config.target_num_bytes = ExUtilGetUInt(argv[++c], &parse_error);
    } else if (!strcmp(argv[c], "-l") && c + 1 < argc) {
      config.importance_of_size_over_quality =
          4.4f * ExUtilGetUInt(argv[++c], &parse_error) / 100.f;
    } else if (!strcmp(argv[c], "-p") && c + 1 < argc) {
      if (!strcmp(argv[c + 1], "none")) {
        ResetProbas(&config);
        c += 1;
      } else if (c + 2 < argc) {
        auto* const dst =
          !strcmp(argv[c + 1], "vmove") ? &config.proba_vertex_move :
          !strcmp(argv[c + 1], "vadd") ? &config.proba_vertex_add :
          !strcmp(argv[c + 1], "vsub") ? &config.proba_vertex_sub :
          !strcmp(argv[c + 1], "cmove") ? &config.proba_color_move :
          !strcmp(argv[c + 1], "cadd") ? &config.proba_color_add :
          !strcmp(argv[c + 1], "csub") ? &config.proba_color_sub :
          !strcmp(argv[c + 1], "cidx") ? &config.proba_color_index_move :
          nullptr;
        CHECK_OK(dst != nullptr, "Invalid -p move '%s'\n", argv[c + 1]);
        *dst = ExUtilGetUInt(argv[c + 2], &parse_error);
        c += 2;
      }
    } else if (!strcmp(argv[c], "-m") && c + 1 < argc) {
      ++c;
      if (!strcmp(argv[c], "edge")) {
        config.analysis_method =
          WP2::PreviewConfig::AnalysisMethod::kEdgeSelectionAndRepulsion;
      } else if (!strcmp(argv[c], "diff")) {
        config.analysis_method =
          WP2::PreviewConfig::AnalysisMethod::kColorDiffMaximization;
      } else if (!strcmp(argv[c], "rand")) {
        config.analysis_method =
          WP2::PreviewConfig::AnalysisMethod::kRandomEdgeSelection;
      } else {
        CHECK_OK(false, "Unknown analysis method '%s'\n", argv[c]);
      }
    } else if (!strcmp(argv[c], "-n") && c + 1 < argc) {
      config.num_vertices = ExUtilGetUInt(argv[++c], &parse_error);
    } else if (!strcmp(argv[c], "-nc") && c + 1 < argc) {
      config.num_colors = ExUtilGetUInt(argv[++c], &parse_error);
    } else if (!strcmp(argv[c], "-density") && c + 1 < argc) {
      config.grid_density = ExUtilGetUInt(argv[++c], &parse_error);
    } else if (!strcmp(argv[c], "-g") && c + 1 < argc) {
      config.grid_size = ExUtilGetUInt(argv[++c], &parse_error);
    } else if (!strcmp(argv[c], "-b") && c + 1 < argc) {
      config.blur_radius = ExUtilGetUInt(argv[++c], &parse_error);
    } else if (!strcmp(argv[c], "-i") && c + 1 < argc) {
      config.num_iterations = ExUtilGetUInt(argv[++c], &parse_error);
    } else if (c < argc && argv[c][0] != '-') {
      if (img_file == nullptr) {
        img_file = argv[c];
      } else {
        bin_file = argv[c];
      }
    } else {
      parse_error = true;
    }
    if (parse_error) {
      Help();
      return 1;
    }
  }
  CHECK_OK(img_file != nullptr, "Missing input image");
  CHECK_OK(config.IsValid(), "PreviewConfig is invalid.\n");

  WP2::ArgbBuffer original;
  CHECK_STATUS(WP2::ReadImage(img_file, &original),
               "Can't read image [%s]\n", img_file);

  WP2::PreviewData data;
  if (bin_file == nullptr) {
    WP2::MemoryWriter writer;
    CHECK_STATUS(WP2::EncodePreview(original, config, &writer),
                "Error in EncodePreview()\n");
    CHECK_STATUS(data.Decode(writer.mem_, writer.size_),
                 "Round-trip error. SHOULDN'T HAPPEN!!'\n");
  } else {
    // Read the starting state
    const bool as_txt = (WP2GetFileExtension(bin_file) == "txt");
    std::string bin;
    CHECK_STATUS(WP2::IoUtilReadFile(bin_file, &bin),
                 "bad bin file [%s]\n", bin_file);
    if (as_txt) {
      CHECK_STATUS(PreviewFromText(bin, &data), "Error in text-to-preview\n");
    } else {
      CHECK_STATUS(data.Decode((const uint8_t*)bin.data(), bin.size()),
                  "Can't decode binary data [%s]", bin_file);
    }
    // optimize further
    CHECK_STATUS(data.Optimize(original, config, /*log=*/true),
                 "Optimize() call failed\n");
  }

  // save binary/text result
  WP2::ANSEnc enc;
  CHECK_STATUS(data.Encode(&enc),
               "Error encoding to binary! SHOULDN'T HAPPEN!!\n");
  const uint32_t preview_size = enc.BufferSize();
  const uint8_t* const preview_data = enc.Buffer();
  VERBOSE("Creating preview for %s. Size:%u\n", img_file, preview_size);

  if (out_file != nullptr) {
    const bool as_txt = (WP2GetFileExtension(out_file) == "txt");
    if (as_txt) {
      const std::string s = PreviewToText(data);
      CHECK_STATUS(
          WP2::IoUtilWriteFile(
              (const uint8_t*)s.data(), s.size(), out_file, overwrite),
          "Couldn't write preview as text [%s]!\n", out_file);
    } else {   // raw output
      CHECK_STATUS(
          WP2::IoUtilWriteFile(
              preview_data, preview_size, out_file, overwrite),
          "Couldn't write preview as binary [%s]!\n", out_file);
    }
    VERBOSE("Saved [%s] as %s\n", out_file, as_txt ? "text" : "binary");
  }

  // save reconstructed preview
  if (preview_file != nullptr) {
    const uint32_t out_h = std::max(32u,
        WP2::DivRound(out_w * original.height, original.width));
    WP2::ArgbBuffer thumb;
    CHECK_STATUS(thumb.Resize(out_w, out_h), "Memory error\n");
    CHECK_STATUS(WP2::DecodePreview(preview_data, preview_size, &thumb),
                 "Error in DecodePreview()\n");
    CHECK_STATUS(WP2::SaveImage(thumb, preview_file, overwrite),
                 "Couldn't save decoded preview as [%s]!\n",
                 preview_file);
    VERBOSE("Saved [%s]\n", preview_file);
  }

  // print final info
  std::string infos = "\n";
  if (out_format == INFO) {
    infos = PrintPreview(data, /* reduced= */false);
  } else if (out_format == SHORT) {
    infos = PrintPreview(data, /* reduced= */true);
  } else if (out_format == TEXT) {
    infos = PreviewToText(data);
  }
  VERBOSE("%s", infos.c_str());

  return 0;
}
