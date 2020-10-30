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
//  Command-line tool for decoding a WP2 image.
//
// Author: Skal (pascal.massimino@gmail.com)

#include <cstdio>
#include <cstring>
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
#include "src/utils/thread_utils.h"

namespace {

constexpr uint32_t kDefaultThreadLevel = 20;  // When -mt is set.

//------------------------------------------------------------------------------

void Help() {
  ProgramOptions opt;
  opt.Add("Decodes the WP2 image file.");
  opt.Add("");
  opt.Add("Usage:");
  opt.Add("  dwp2 in_file [options] [-o out_file]");
  opt.Add("");
  opt.Add("Use the following options to force the conversion into:");
  opt.Add("-pam", "save the raw RGBA samples as a color PAM");
  opt.Add("-ppm", "save the raw RGB samples as a color PPM");
  opt.Add("-bmp", "save as uncompressed BMP format");
  opt.Add("-tiff", "save as uncompressed TIFF format");
  opt.Add("-png", "save the raw RGBA samples as PNG (default)");
  opt.Add("");
  opt.Add("Other options are:");
  opt.Add("-mt [int]", SPrintf("enable multi-threading, with %u extra threads "
                               "if unspecified (0 to disable, default is %u)",
                               kDefaultThreadLevel,
                               WP2::DecoderConfig::kDefault.thread_level));
  opt.Add("-grain <int>", "grain amplitude (if available): [0=off .. 100]");
  opt.Add(
      "-[no_]dblk_filter",
      SPrintf("enable or disable deblocking filter (%s by default)",
              WP2::DecoderConfig::kDefault.enable_deblocking_filter ? "on"
                                                                    : "off"));
  opt.Add(
      "-[no_]drct_filter",
      SPrintf("enable or disable directional filter (%s by default)",
              WP2::DecoderConfig::kDefault.enable_directional_filter ? "on"
                                                                     : "off"));
  opt.Add(
      "-[no_]rstr_filter",
      SPrintf("enable or disable restoration filter (%s by default)",
              WP2::DecoderConfig::kDefault.enable_restoration_filter ? "on"
                                                                     : "off"));
  opt.Add("-r", "recurse into input directories");
  opt.Add("-inplace",
          "output files in the same directories as input files (replaces -o)");
  opt.Add("-frames",
          "output frames named \"frame[index]_[duration]ms\" to the directory "
          "specified with -o");
  opt.Add("-force",
          "overwrite destination if it exists (default=off unless out_file is "
          "explicit)");
  opt.Add("-i", "just print bitstream info and exit");
  opt.Add("-v", "verbose (e.g. print encoding/decoding times)");
  opt.Add("-quiet", "quiet mode, don't print anything");
#ifdef WP2_BITTRACE
  opt.Add("-bt / -BT", "set bit trace level to 1 / 2");
  opt.Add("-histos", "display histograms per symbol when -bt is set");
#endif
  opt.Add("");
  opt.AddSystemOptionSection();
  opt.Print();
}

//------------------------------------------------------------------------------

struct DecodingSettings {
  // File writing.
  std::string out_path;
  bool inplace_output = false;
  bool output_frames_to_directory = false;
  bool output_to_directory = false;
  bool allow_overwrite = false;
  WP2::FileFormat out_format = WP2::FileFormat::AUTO;

  // Standard output settings.
  bool verbose = false;
  bool quiet = false;
  bool print_info = false;
  bool short_output = false;

  WP2::DecoderConfig config;
  const char* visual_debug = nullptr;

  // For WP2_BITTRACE.
  int bit_trace = 0;  // 0=none, 1=trace in bits, 2=trace in bytes
  uint32_t bit_trace_level = 0;
  bool show_histograms = false;
  bool count_blocks = false;
};

//------------------------------------------------------------------------------

// Get images to decode from input arguments.
std::vector<std::string> GetFilesToDecode(const std::vector<std::string>& paths,
                                          bool recursive_input, bool verbose,
                                          bool quiet) {
  std::vector<std::string> files_to_decode;
  for (const std::string& path : paths) {
    if (WP2IsDirectory(path)) {
      if (recursive_input) {
        std::vector<std::string> files_in_dir;
        EXIT_IF_FALSE(WP2GetFilesIn(path, &files_in_dir, /*recursive=*/true),
                      "Error! Opening directory %s failed.", path.c_str());

        // Only keep files with a wp2 extension.
        for (const std::string& file : files_in_dir) {
          WP2::Data data;
          if (WP2::IoUtilReadFile(file.c_str(), &data, /*max_num_bytes=*/12) !=
              WP2_STATUS_OK) {
            if (verbose) {
              printf("Could not read file %s: skipping.\n", file.c_str());
            }
          } else if (WP2::GuessImageFormat(data.bytes, data.size) !=
                     WP2::FileFormat::WP2) {
            if (verbose) {
              printf("File %s has no WP2 header: skipping.\n", file.c_str());
            }
          } else {
            files_to_decode.push_back(file);
          }
        }
      } else if (!quiet) {
        fprintf(stderr,
                "Ignoring directory %s because -r was not specified.\n",
                path.c_str());
      }
    } else {
      // No extension checking for explicit input files.
      files_to_decode.push_back(path);
    }
  }
  return files_to_decode;
}

uint8_t GetNumDigits(size_t number) {
  uint8_t num_digits = 1;
  for (number /= 10; number > 0; ++num_digits) number /= 10;
  return num_digits;
}

std::string GetFrameFilename(size_t index, uint32_t duration_ms) {
  std::string index_str = std::to_string(index);
  // Add leading zeros.
  const size_t max_num_digits = GetNumDigits(WP2::kMaxNumFrames - 1);
  index_str = std::string(max_num_digits - index_str.size(), '0') + index_str;
  return "frame" + index_str + "_" + std::to_string(duration_ms) + "ms";
}

//------------------------------------------------------------------------------

std::string PrintSummary(const WP2::BitstreamFeatures& bitstream,
                         size_t input_data_size, bool short_output) {
  const WP2::Argb32b prev_col = WP2::ToArgb32b(bitstream.preview_color);
  const uint32_t degrees = (uint32_t)bitstream.orientation * 90;

  std::string str;
  str += SPrintf("Dimension:     %u x %u\n", bitstream.width, bitstream.height);
  if (!short_output) {
    str += SPrintf("animation:     %s", bitstream.is_animation ? "yes" : "no");
    if (bitstream.is_animation) {
      str += SPrintf(" (loop count:%d)", bitstream.loop_count);
    }
    str += "\n";
    str += SPrintf("quality hint:  %u\n", bitstream.quality_hint);
    str += SPrintf("preview:       %s\n", bitstream.has_preview ? "yes" : "no");
    str += SPrintf("preview color: %u, %u, %u\n", prev_col.r, prev_col.g,
                   prev_col.b);
    str += SPrintf("rotated by:    %u degrees clockwise\n", degrees);
    str += SPrintf("transfer:      %d\n", (int)bitstream.transfer_function);
    str += SPrintf("ICC:           %s\n", bitstream.has_icc ? "yes" : "no");
    str += SPrintf("XMP/EXIF:      %s\n",
                   bitstream.has_trailing_data ? "yes" : "no");
    str += SPrintf("file size:     %zu\n", input_data_size);
  }
  return str;
}

//------------------------------------------------------------------------------

class ImageDecoder : public WP2::WorkerBase {
 public:
  ImageDecoder(const DecodingSettings& settings,
               const std::vector<std::string>& image_pool,
               size_t* const next_image_index, WP2::ThreadLock* const mutex)
      : settings_(settings),
        image_pool_(image_pool),
        next_image_index_(next_image_index),
        mutex_(mutex) {}

 protected:
  WP2Status SaveImage(const std::string& in_path, const WP2::ArgbBuffer& buffer,
                      std::string* const standard_out,
                      std::string* const standard_err) {
    if (!settings_.quiet) {
      *standard_out +=
          SPrintf("Decoded %s. Dimensions: %d x %d. Transparency: %s.\n",
                  in_path.c_str(), buffer.width, buffer.height,
                  buffer.HasTransparency() ? "yes" : "no");
    }
    if (!settings_.out_path.empty()) {
      std::string final_out_path;
      if (settings_.output_to_directory) {
        const char* const extension =
            GetExtensionFromFormat(settings_.out_format);
        WP2_CHECK_OK(extension != nullptr, WP2_STATUS_UNSUPPORTED_FEATURE);
        if (settings_.inplace_output) {
          final_out_path = WP2RemoveFileExtension(in_path) + '.' + extension;
        } else {
          final_out_path =
              WP2InputFileToOutputDir(in_path, settings_.out_path, extension);
        }
      } else {
        final_out_path = settings_.out_path;
      }

      const double start = GetStopwatchTime();
      const WP2Status status =
          WP2::SaveImage(buffer, final_out_path.c_str(),
                         settings_.allow_overwrite, settings_.out_format);
      if (!settings_.quiet && status != WP2_STATUS_OK) {
        *standard_err +=
            SPrintf("Could not save to '%s'\n", final_out_path.c_str());
        if (status == WP2_STATUS_BAD_WRITE && !settings_.allow_overwrite) {
          *standard_err +=
              "(hint: use -force flag to overwrite an existing file)\n";
        }
        return status;
      }

      if (!settings_.quiet) {
        *standard_out += SPrintf("Saved to '%s'\n", final_out_path.c_str());
      }
      if (settings_.verbose) {
        const double write_time = GetStopwatchTime() - start;
        *standard_out += SPrintf("Time to write output: %.3fs\n", write_time);
      }
    }
    return WP2_STATUS_OK;
  }

  WP2Status DecodeAndSaveImage(const std::string& in_path,
                               std::string* const standard_out,
                               std::string* const standard_err) {
    WP2::ArgbBuffer output_buffer(WP2_Argb_32);

    WP2::Data data;
    WP2_CHECK_STATUS(WP2::IoUtilReadFile(in_path.c_str(), &data));

    WP2::BitstreamFeatures bitstream;
    WP2_CHECK_STATUS(bitstream.Read(data.bytes, data.size));

    if (settings_.print_info) {
      *standard_out +=
          PrintSummary(bitstream, data.size, settings_.short_output);
      return WP2_STATUS_OK;
    }

    WP2::DecoderConfig config = settings_.config;
    WP2::DecoderInfo info;
    if (settings_.verbose || settings_.bit_trace || settings_.count_blocks ||
        settings_.visual_debug != nullptr) {
      info.visual_debug = settings_.visual_debug;
      config.info = &info;  // report extra statistics
      if (settings_.count_blocks) info.store_blocks = true;
    }

    const WP2::ArgbBuffer& buffer_to_save =
        (settings_.visual_debug != nullptr) ? info.debug_output : output_buffer;

    if (settings_.output_frames_to_directory) {
      WP2::ArrayDecoder decoder(data.bytes, data.size, config, &output_buffer);
      size_t num_frames = 0;
      double start = GetStopwatchTime();
      uint32_t duration_ms;
      while (decoder.ReadFrame(&duration_ms)) {
        if (settings_.verbose) {
          const double decode_time = GetStopwatchTime() - start;
          *standard_out +=
              SPrintf("Time to decode picture: %.3fs\n", decode_time);
          start = GetStopwatchTime();
        }

        const std::string frame_filename =
            GetFrameFilename(num_frames, duration_ms);
        WP2_CHECK_STATUS(SaveImage(frame_filename, buffer_to_save, standard_out,
                                   standard_err));
        ++num_frames;
      }

      WP2_CHECK_STATUS(decoder.GetStatus());

      // Output a frame even if it is a still image.
      if (!bitstream.is_animation) {
        assert(num_frames == 1);
        const std::string frame_filename =
            GetFrameFilename(/*index=*/0, WP2::kMaxFrameDurationMs);
        WP2_CHECK_STATUS(SaveImage(frame_filename, buffer_to_save, standard_out,
                                   standard_err));
      }
    } else {
      const double start = GetStopwatchTime();
      const WP2Status status =
          WP2::Decode(data.bytes, data.size, &output_buffer, config);
      if (settings_.verbose) {
        const double decode_time = GetStopwatchTime() - start;
        *standard_out +=
            SPrintf("Time to decode picture: %.3fs\n", decode_time);
      }
      WP2_CHECK_STATUS(status);

      WP2_CHECK_STATUS(
          SaveImage(in_path, buffer_to_save, standard_out, standard_err));
    }

    if (!settings_.quiet && (settings_.bit_trace || settings_.count_blocks)) {
      if (image_pool_.size() > 1) {
        *standard_err +=
            "Bit traces and block count are printed only for single images.\n";
      } else {
#if !defined(WP2_BITTRACE)
        *standard_err +=
            "Bit traces and block count are not available without WP2_BITTRACE "
            "compile flag.\n";
#else
        if (settings_.bit_trace) {
          PrintBitTraces(info, data.size, true, settings_.bit_trace == 2,
                         settings_.show_histograms, settings_.short_output,
                         settings_.bit_trace_level);
        }
        if (settings_.count_blocks) {
          printf("num-blocks : %zu\n", info.blocks.size());
        }
#endif
      }
    }
    return WP2_STATUS_OK;
  }

  WP2Status Execute() override {
    while (true) {
      // Get the next image to decode from the shared pool, if any.
      if (mutex_ != nullptr) WP2_CHECK_STATUS(mutex_->Acquire());
      if (*next_image_index_ >= image_pool_.size()) {
        // Nothing left to do.
        if (mutex_ != nullptr) mutex_->Release();
        break;
      }
      const std::string& in_path = image_pool_[*next_image_index_];
      ++*next_image_index_;
      if (mutex_ != nullptr) mutex_->Release();

      std::string standard_out, standard_err;
      const WP2Status status =
          DecodeAndSaveImage(in_path, &standard_out, &standard_err);

      if (!settings_.quiet) {
        // Acquire mutex to be sure to print in order without mixing lines.
        if (mutex_ != nullptr) WP2_CHECK_STATUS(mutex_->Acquire());
        printf("%s", standard_out.c_str());
        fprintf(stderr, "%s", standard_err.c_str());
        if (mutex_ != nullptr) mutex_->Release();
      }
      WP2_CHECK_STATUS(status);
    }
    return WP2_STATUS_OK;
  }

 protected:
  // Shared settings, const so no need to thread-lock it.
  const DecodingSettings& settings_;

  // Shared pool of images to decode.
  const std::vector<std::string>& image_pool_;
  size_t* const next_image_index_ = nullptr;
  WP2::ThreadLock* const mutex_ = nullptr;  // Acquire/release if not null.
};

//------------------------------------------------------------------------------

}  // namespace

int main(int argc, char* argv[]) {
  std::vector<std::string> in_paths;
  bool recursive_input = false;
  DecodingSettings settings = {};
  uint32_t thread_level = 0;
  int c;

  EXIT_IF_FALSE(WP2CheckVersion(), "Error! Library version mismatch!");

  for (c = 1; c < argc; ++c) {
    bool parse_error = false;
    if (!strcmp(argv[c], "-h") || !strcmp(argv[c], "-help")) {
      Help();
      return 0;
    } else if (!strcmp(argv[c], "-o") && c + 1 < argc) {
      EXIT_IF_FALSE(settings.out_path.empty() && !settings.inplace_output,
                    "Error! Output was already specified.");
      settings.out_path = argv[++c];
    } else if (!strcmp(argv[c], "-png")) {
      settings.out_format = WP2::FileFormat::PNG;
    } else if (!strcmp(argv[c], "-pam")) {
      settings.out_format = WP2::FileFormat::PAM;
    } else if (!strcmp(argv[c], "-ppm")) {
      settings.out_format = WP2::FileFormat::PPM;
    } else if (!strcmp(argv[c], "-bmp")) {
      settings.out_format = WP2::FileFormat::BMP;
    } else if (!strcmp(argv[c], "-tiff")) {
      settings.out_format = WP2::FileFormat::TIFF;
    } else if (!strcmp(argv[c], "-mt")) {
      if (c + 1 < argc && ExUtilTryGetUInt(argv[c + 1], &thread_level)) {
        ++c;
      } else {
        thread_level = kDefaultThreadLevel;
      }
    } else if (!strcmp(argv[c], "-grain") && c + 1 < argc) {
      settings.config.grain_amplitude = ExUtilGetUInt(argv[++c], &parse_error);
    } else if (!strcmp(argv[c], "-dblk_filter") ||
               !strcmp(argv[c], "-no_dblk_filter")) {
      settings.config.enable_deblocking_filter =
          !strcmp(argv[c], "-dblk_filter");
    } else if (!strcmp(argv[c], "-drct_filter") ||
               !strcmp(argv[c], "-no_drct_filter")) {
      settings.config.enable_directional_filter =
          !strcmp(argv[c], "-drct_filter");
    } else if (!strcmp(argv[c], "-rstr_filter") ||
               !strcmp(argv[c], "-no_rstr_filter")) {
      settings.config.enable_restoration_filter =
          !strcmp(argv[c], "-rstr_filter");
    } else if (!strcmp(argv[c], "-r")) {
      recursive_input = true;
    } else if (!strcmp(argv[c], "-inplace")) {
      EXIT_IF_FALSE(!settings.output_frames_to_directory,
                    "Error! -inplace is not compatible with -frames.");
      EXIT_IF_FALSE(settings.out_path.empty(),
                    "Error! Output was already specified.");
      settings.inplace_output = true;
    } else if (!strcmp(argv[c], "-frames")) {
      EXIT_IF_FALSE(!settings.inplace_output,
                    "Error! -frames is not compatible with -inplace.");
      settings.output_frames_to_directory = true;
    } else if (!strcmp(argv[c], "-force")) {
      settings.allow_overwrite = true;
    } else if (!strcmp(argv[c], "-i")) {
      settings.print_info = true;
    } else if (!strcmp(argv[c], "-vdebug") && c + 1 < argc) {
      settings.visual_debug = argv[++c];
    } else if (!strcmp(argv[c], "-bt") || !strcmp(argv[c], "-BT")) {
      settings.bit_trace = !strcmp(argv[c], "-bt") ? 1 : 2;
      if (c + 1 < argc) {
        if (isdigit(argv[c + 1][0])) {
          settings.bit_trace_level = ExUtilGetUInt(argv[++c], &parse_error);
        } else {
          settings.bit_trace_level = 0;
        }
      }
    } else if (!strcmp(argv[c], "-histos")) {
      settings.show_histograms = true;
    } else if (!strcmp(argv[c], "-count_blocks")) {
      settings.count_blocks = true;
    } else if (!strcmp(argv[c], "-quiet")) {
      settings.quiet = true;
    } else if (!strcmp(argv[c], "-short")) {
      settings.short_output = true;
    } else if (!strcmp(argv[c], "-v")) {
      settings.verbose = true;
    } else if (!strcmp(argv[c], "--")) {
      if (c + 1 < argc) in_paths.emplace_back(argv[++c]);
      break;
    } else if (argv[c][0] == '-') {
      bool must_stop;
      if (ProgramOptions::ParseSystemOptions(argv[c], &must_stop)) {
        if (must_stop) return 0;
      } else {
        printf("Unknown option '%s'\n", argv[c]);
        Help();
        return 1;
      }
    } else {
      in_paths.emplace_back(argv[c]);
    }

    if (parse_error) {
      Help();
      return 1;
    }
  }

  const std::vector<std::string> in_files = GetFilesToDecode(
      in_paths, recursive_input, settings.verbose, settings.quiet);
  if (in_files.empty()) {
    fprintf(stderr, "missing input file!!\n");
    Help();
    return 1;
  }

  if (settings.quiet) settings.verbose = false;

  settings.output_to_directory =
      (settings.inplace_output || WP2IsDirectory(settings.out_path));
  // Destination is explicitly specified, allow overwriting it.
  if (!settings.output_to_directory) settings.allow_overwrite = true;

  if (settings.output_frames_to_directory) {
    EXIT_IF_FALSE(in_files.size() == 1,
                  "Error! -frames requires a single input file.");
    EXIT_IF_FALSE(settings.output_to_directory,
                  "Error! -frames requires a single output directory.");
  } else {
    EXIT_IF_FALSE(
        in_files.size() == 1 || settings.output_to_directory ||
            settings.out_path.empty(),
        "Error! Several input images require the output to be an existing "
        "directory.");
  }

  // Use a pool of images to decode.
  size_t next_image_index = 0;
  WP2::ThreadLock mutex;

  // Create a pool of workers that will get jobs from 'in_files'.
  const size_t num_workers =
      std::min(1 + (size_t)thread_level, in_files.size());
  const bool use_mt_workers = num_workers > 1;

  // Enable multithreading at image level or at tile level but not both.
  // TODO(yguyon): This could be optimized depending on the number of tiles etc.
  settings.config.thread_level = use_mt_workers ? 0 : thread_level;

  std::vector<ImageDecoder> workers;
  workers.reserve(num_workers);

  WP2Status status = WP2_STATUS_OK;
  for (size_t i = 0; i < num_workers; ++i) {
    workers.emplace_back(settings, in_files, &next_image_index,
                         use_mt_workers ? &mutex : nullptr);
    status = workers.back().Start(use_mt_workers);
    if (status != WP2_STATUS_OK) break;
  }

  for (ImageDecoder& worker : workers) {
    const WP2Status end_status = worker.End();
    if (status == WP2_STATUS_OK) status = end_status;
    // Still wait for others to finish in case of error.
  }
  EXIT_IF_ERROR(status, "Error: %s", WP2GetStatusText(status));

  if (!settings.quiet && settings.out_path.empty()) {
    printf("Nothing written; use -o flag to save the result.\n");
  }
  return 0;
}

//------------------------------------------------------------------------------
