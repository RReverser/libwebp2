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
//  Helper functions used by examples.

#include "./example_utils.h"

#ifdef _WIN32
#include <windows.h>
#else
#include <dirent.h>
#endif  // _WIN32

#include <cstdio>
#include <string>
#include <vector>

//------------------------------------------------------------------------------

// Prints the distortion stored in 'values'.
static void PrintDistortion(WP2::MetricType metric, float values[5],
                            size_t data_size, bool short_output) {
  if (!short_output) {
    fprintf(stderr, "%s: ", kWP2MetricNames[metric]);
    if (metric == WP2::PSNR_YUV || metric == WP2::SSIM_YUV ||
        metric == WP2::PSNRHVS) {
      fprintf(stderr, "A: %2.2f  Y: %2.2f  U: %2.2f  V: %2.2f   ",
              values[0], values[1], values[2], values[3]);
    } else {
      fprintf(stderr, "A: %2.2f  R: %2.2f  G: %2.2f  B: %2.2f   ",
              values[0], values[1], values[2], values[3]);
    }
    fprintf(stderr, "Total size: %7d %s: %2.2f\n",
            (int)data_size, kWP2MetricNames[metric], values[4]);
  } else {
    fprintf(stderr, "%7d %.4f\n", (int)data_size, values[4]);
  }
}

void WP2PrintDistortion(const std::vector<WP2::ArgbBuffer>& refs,
                        const uint8_t* const data, size_t data_size,
                        WP2::MetricType metric, bool short_output) {
  WP2::ArrayDecoder decoder(data, data_size);
  size_t num_frames = 0;
  float avg_values[5] = {0};
  while (decoder.ReadFrame()) {
    EXIT_IF_FALSE(refs.size() > num_frames,
                  "Round trip error!! This shouldn't happen.");
    float values[5];
    EXIT_IF_ERROR(decoder.GetPixels().GetDistortionBlackOrWhiteBackground(
                      refs[num_frames], metric, values),
                  "Error while computing distortion!");
    if (!short_output) {
      fprintf(stderr, "Frame %3zu distortion ", num_frames);
      PrintDistortion(metric, values, data_size, short_output);
    }
    for (int i = 0; i < 5; ++i) avg_values[i] += values[i];
    ++num_frames;
  }
  EXIT_IF_ERROR(decoder.GetStatus(),
                "Round trip error!! This shouldn't happen.");

  if (num_frames > 0) {
    for (int i = 0; i < 5; ++i) avg_values[i] /= num_frames;
    if (!short_output) {
      fprintf(stderr, "Total avg distortion over %zu frames ", num_frames);
    }
    PrintDistortion(metric, avg_values, data_size, short_output);
  } else {
    float values[5];
    EXIT_IF_FALSE(refs.size() == 1,
                  "Round trip error!! This shouldn't happen.");
    EXIT_IF_ERROR(decoder.GetPixels().GetDistortionBlackOrWhiteBackground(
                      refs.front(), metric, values),
                  "Error while computing distortion!");
    PrintDistortion(metric, values, data_size, short_output);
  }
}

void WP2PrintDistortion(const WP2::ArgbBuffer& ref, const uint8_t* const data,
                        size_t data_size, WP2::MetricType metric,
                        bool short_output) {
  std::vector<WP2::ArgbBuffer> refs(1);
  if (refs.front().SetView(ref) != WP2_STATUS_OK) assert(false);
  WP2PrintDistortion(refs, data, data_size, metric, short_output);
}

//------------------------------------------------------------------------------

void WP2CollectBitTraces(const uint8_t* const data, size_t data_size,
                         int bit_trace, uint32_t bit_trace_level,
                         bool short_output, bool show_histograms) {
  WP2::DecoderConfig config;
  WP2::DecoderInfo info;
  config.info = &info;

  WP2::ArrayDecoder decoder(data, data_size, config);
  while (decoder.ReadFrame()) continue;  // Decode everything.
  EXIT_IF_ERROR(decoder.GetStatus(),
                "Round trip error!! This shouldn't happen.");

#if defined(WP2_BITTRACE)
  PrintBitTraces(info, data_size, /*sort_value=*/true,
                 /*use_bytes=*/bit_trace == 2, show_histograms,
                 /*short_version=*/short_output, bit_trace_level);
#else
  (void)bit_trace;
  (void)bit_trace_level;
  (void)short_output;
  fprintf(stderr,
          "Bit traces are not available without "
          "WP2_BITTRACE compile flag.\n");
#endif  // WP2_BITTRACE
}

//------------------------------------------------------------------------------

#ifdef _WIN32
static constexpr char kPathSeparator = (char)'\\';
#else
static constexpr char kPathSeparator = (char)'/';
#endif  // _WIN32

std::string WP2JoinPath(const std::string& prefix, const std::string& suffix) {
  std::string joined_path = prefix;
  if (!joined_path.empty() && joined_path.back() != kPathSeparator) {
    joined_path += kPathSeparator;
  }
  joined_path += suffix;
  return joined_path;
}

bool WP2IsDirectory(const std::string& path) {
  if (path.empty()) return false;
#ifdef _WIN32
  auto attributes = GetFileAttributes(path.c_str());
  if (attributes == INVALID_FILE_ATTRIBUTES) return false;
  return (attributes == FILE_ATTRIBUTE_DIRECTORY);
#else
  DIR* dir = opendir(path.c_str());
  if (dir == nullptr) return false;
  (void)closedir(dir);
  return true;
#endif  // _WIN32
}

bool WP2GetFilesIn(const std::string& dir_path,
                   std::vector<std::string>* const files, bool recursive) {
#ifdef _WIN32
  WIN32_FIND_DATA ffd;
  HANDLE handle = FindFirstFile(WP2JoinPath(dir_path, "*").c_str(), &ffd);
  if (handle == INVALID_HANDLE_VALUE) return false;

  do {
    const std::string file_path = WP2JoinPath(dir_path, ffd.cFileName);
    if (ffd.dwFileAttributes & FILE_ATTRIBUTE_DIRECTORY) {
      if (recursive) WP2GetFilesIn(file_path, files, recursive);
    } else {
      files->push_back(file_path);
    }
  } while (FindNextFile(handle, &ffd) != 0);

  if (GetLastError() != ERROR_NO_MORE_FILES) {
    FindClose(handle);
    return false;
  }
  FindClose(handle);
  return true;
#else
  DIR* dir = opendir(dir_path.c_str());
  if (dir == nullptr) return false;
  struct dirent* ent;
  while ((ent = readdir(dir)) != nullptr) {
    if ((strcmp(ent->d_name, ".") == 0) || (strcmp(ent->d_name, "..") == 0)) {
      continue;
    }
    const std::string ent_path = WP2JoinPath(dir_path, ent->d_name);
    if (recursive) {
      if (!WP2GetFilesIn(ent_path, files, recursive)) {
        files->push_back(ent_path);
      }
    } else if (!WP2IsDirectory(ent_path)) {
      files->push_back(ent_path);
    }
  }
  closedir(dir);
  return true;
#endif  // _WIN32
}

std::string WP2GetFileName(const std::string& path) {
  const size_t last_separator_pos = path.find_last_of(kPathSeparator);
  if (last_separator_pos == std::string::npos) return path;
  return path.substr(last_separator_pos + 1);
}

std::string WP2RemoveFileExtension(const std::string& path) {
  const size_t last_separator_pos = path.find_last_of(kPathSeparator);
  const size_t last_dot_pos = path.find_last_of('.');
  if ((last_dot_pos == std::string::npos) ||
      ((last_separator_pos != std::string::npos) &&
       (last_separator_pos > last_dot_pos))) {
    return path;
  }
  return path.substr(0, last_dot_pos);
}

std::string WP2GetFileExtension(const std::string& path) {
  const size_t last_separator_pos = path.find_last_of(kPathSeparator);
  const size_t last_dot_pos = path.find_last_of('.');
  if ((last_dot_pos == std::string::npos) ||
      (last_dot_pos + 1 == path.length()) ||
      ((last_separator_pos != std::string::npos) &&
       (last_separator_pos > last_dot_pos))) {
    return "";
  }
  return path.substr(last_dot_pos + 1);
}

std::string WP2InputFileToOutputDir(const std::string& input_file_path,
                                    const std::string& output_dir_path,
                                    const std::string& extension) {
  assert(!extension.empty());
  return WP2JoinPath(output_dir_path,
                     WP2RemoveFileExtension(WP2GetFileName(input_file_path)) +
                         '.' + extension);
}

//------------------------------------------------------------------------------

// Handy class just for collecting all params / work variable in a single place.
struct PartitionParams {
 public:
  const WP2::BlockSize* valid_blocks;
  uint32_t max_width, max_height;
  bool skip_invalid;
  uint32_t split_count;
  std::vector<WP2::Rectangle> partition;

 public:
  // search the valid dimension that is strictly smaller than 'd'.
  uint32_t SearchValidDim(uint32_t d, bool search_width) const {
    uint32_t new_d = WP2::kMinBlockSizePix;
    for (const WP2::BlockSize* b = valid_blocks; *b != WP2::BLK_LAST; ++b) {
      const uint32_t test_d = search_width ? WP2::BlockWidthPix(*b)
                                           : WP2::BlockHeightPix(*b);
      if (test_d > new_d && test_d < d) new_d = test_d;
    }
    return new_d;
  }
  void AddBox(const WP2::Rectangle& r) {
    uint32_t w = std::min(r.width,  max_width - r.x);
    uint32_t h = std::min(r.height, max_height - r.y);
    if (w < 4u || h < 4u) return;
    w = (w + 3u) & ~3u;
    h = (h + 3u) & ~3u;

    for (const WP2::BlockSize* b = valid_blocks; *b != WP2::BLK_LAST; ++b) {
      if (w == WP2::BlockWidthPix(*b) && h == WP2::BlockHeightPix(*b)) {
        // we found a valid block to push.
        partition.emplace_back(r.x, r.y, w, h);
        return;
      }
    }
    if (skip_invalid) return;
    // We have invalid block to recursively split along the largest dimension.
    split_count += 1;
    if (w >= h) {
      const uint32_t new_w = SearchValidDim(w, /*search_width=*/true);
      AddBox({r.x +     0, r.y,     new_w, h});
      AddBox({r.x + new_w, r.y, w - new_w, h});
    } else {
      const uint32_t new_h = SearchValidDim(h, /*search_width=*/false);
      AddBox({r.x, r.y +     0, w,     new_h});
      AddBox({r.x, r.y + new_h, w, h - new_h});
    }
  }
};

bool WP2ReadPartition(const char file_path[], bool read_only,
                      std::vector<WP2::Rectangle>* const partition) {
  std::FILE* const file = std::fopen(file_path, "r");
  if (file == nullptr) {
    // If not 'read_only', it's considered output, not an error.
    return !read_only;
  }
  char line[128];
  while (std::fgets(line, sizeof(line), file) != nullptr) {
    if (line[0] == '\0' || line[0] == '\n') continue;
    WP2::Rectangle r;
    if (std::sscanf(line, "%u %u %u %u",  // NOLINT (sscanf() not recommended)
                    &r.x, &r.y, &r.width, &r.height) != 4) {
      std::fclose(file);
      return false;  // Parsing failed.
    }
    partition->emplace_back(r);
  }
  std::fclose(file);
  return true;
}

void WP2ConvertPartition(uint32_t pic_width, uint32_t pic_height,
                         WP2::PartitionSet partition_set, bool ignore_invalid,
                         std::vector<WP2::Rectangle>* const partition) {
  PartitionParams p = {WP2::GetBlockSizes(partition_set),
                       pic_width,
                       pic_height,
                       ignore_invalid,
                       /*split_count=*/0,
                       /*partition=*/{}};
  for (WP2::Rectangle r : *partition) p.AddBox(r);
  if (p.split_count > 0) {
    fprintf(stderr, "Number of splits: %u\n", p.split_count);
  }
  std::swap(*partition, p.partition);
}

bool WP2WritePartition(const std::vector<WP2::Rectangle>& partition,
                       const char file_path[]) {
  std::FILE* const file = std::fopen(file_path, "w");
  if (file == nullptr) return false;
  for (const WP2::Rectangle& rect : partition) {
    fprintf(file, "%u %u %u %u\n", rect.x, rect.y, rect.width, rect.height);
  }
  std::fclose(file);
  return true;
}
//------------------------------------------------------------------------------
// Program options.

void ProgramOptions::Add(const char arg[], const std::string& desc) {
  content_.emplace_back(arg, desc);
}

void ProgramOptions::Add(const std::string& desc) {
  content_.emplace_back(nullptr, desc);
}

void ProgramOptions::AddSystemOptionSection() {
  Add("System options:");
  // TODO(vrabaud) add("-mt", "use multi-threading");
#ifndef WP2_DLL
  Add("-noasm", "disable all assembly optimizations");
#endif
  Add("-version", "print version number and exit");
  Add("-h / -help", "this help");
  Add("");
}

void ProgramOptions::AddMetricOptions() {
  Add("-ssim / -ssim_yuv", "print SSIM distortion in RGB / YCoCg");
  Add("-msssim", "print MSSSIM distortion in RGB");
  Add("-psnr / -psnr_yuv", "print PSNR distortion in RGB / YCoCg");
  Add("-lsim", "print LSIM distortion in RGB");
  Add("-psnrhvs", "print PSNRHVS distortion");
}

// Writes a string by wrapping it at 80 columns.
static void SafePrint(const char s_in[], uint32_t margin, uint32_t max = 80u) {
  margin = std::min(margin, max - 1);
  assert(margin < 80 && max <= 80);
  // print truncated word-delimited string
  while (margin + strlen(s_in) >= max) {
    char tmp[80 + 1] = { 0 };
    snprintf(tmp, max - margin, "%s", s_in);
    char* end = strrchr(tmp, ' ');
    if (end == NULL) end = tmp + max - margin;   // word is too long. Ellipsis?
    *end = 0;
    s_in += end + 1 - tmp;
    printf("%s\n%*s", tmp, margin, "");
  }
  printf("%s\n", s_in);
}

void ProgramOptions::Print() const {
  // Figure out the longest argument.
  size_t size_max = 0u;
  for (const auto& c : content_) {
    if (c.first != nullptr) size_max = std::max(size_max, strlen(c.first));
  }
  // left-justify and print
  for (uint32_t i = 0; i < content_.size(); ++i) {
    if (content_[i].first == nullptr) {   // Print normal text.
      printf("%s\n", content_[i].second.c_str());
    } else {                              // Print the arguments.
      const size_t padding = strlen(content_[i].first);
      printf("  %s %.*s.. ", content_[i].first, (int)(size_max - padding),
             "...............................................................");
      SafePrint(content_[i].second.c_str(), size_max + 6);
    }
  }
}

bool ProgramOptions::ParseSystemOptions(const char* const arg,
                                        bool* const must_stop) {
  *must_stop = false;
  if (!strcmp(arg, "-version")) {
    const auto version = (uint32_t)WP2GetVersion();
    printf("WebP2 version: %d.%d.%d\n", (version >> 16u) & 0xffu,
           (version >> 8u) & 0xffu, (version >> 0u) & 0xffu);
    *must_stop = true;
    return true;
  } else if (!strcmp(arg, "-noasm")) {
    ExUtilDisableSIMD();
    return true;
  }
  return false;
}

bool ProgramOptions::ParseMetricOptions(const char* const arg,
                                        WP2::MetricType* type) {
  if (!strcmp(arg, "-psnr")) {
    *type = WP2::PSNR;
  } else if (!strcmp(arg, "-psnrhvs")) {
    *type = WP2::PSNRHVS;
  } else if (!strcmp(arg, "-psnr_yuv")) {
    *type = WP2::PSNR_YUV;
  } else if (!strcmp(arg, "-ssim")) {
    *type = WP2::SSIM;
  } else if (!strcmp(arg, "-ssim_yuv")) {
    *type = WP2::SSIM_YUV;
  } else if (!strcmp(arg, "-msssim")) {
    *type = WP2::MSSSIM;
  } else if (!strcmp(arg, "-lsim")) {
    *type = WP2::LSIM;
  } else {
    return false;
  }
  return true;
}

//------------------------------------------------------------------------------

// Parse numerical values (and update *error in case of error).
int ExUtilGetInt(const char* const v, bool* const error, int base) {
  int n;
  if (!ExUtilTryGetInt(v, &n, base) && error != nullptr && !*error) {
    *error = true;
    fprintf(stderr, "Error! '%s' is not an integer.\n",
            (v != nullptr) ? v : "(null)");
  }
  return n;
}
uint32_t ExUtilGetUInt(const char* const v, bool* const error, int base) {
  return (uint32_t)ExUtilGetInt(v, error, base);
}

bool ExUtilTryGetInt(const char* const v, int* const n, int base) {
  char* end = nullptr;
  *n = (v != nullptr) ? (int)strtol(v, &end, base) : 0;  // NOLINT
  return (end != v);
}
bool ExUtilTryGetUInt(const char* const v, uint32_t* const n, int base) {
  char* end = nullptr;
  *n = (v != nullptr) ? (uint32_t)strtoul(v, &end, base) : 0;  // NOLINT
  return (end != v);
}

float ExUtilGetFloat(const char* const v, bool* const error) {
  char* end = nullptr;
  const float f = (v != nullptr) ? (float)strtod(v, &end) : 0.f;
  if (end == v && error != nullptr && !*error) {
    *error = true;
    fprintf(stderr, "Error! '%s' is not a floating point number.\n",
            (v != nullptr) ? v : "(null)");
  }
  return f;
}

uint32_t ExUtilGetFloats(const char* v, float values[], uint32_t num_values_max,
                         bool* const error, char separator) {
  uint32_t n = 0;
  bool parse_error = (v == nullptr);
  if (!parse_error) {
    while (n < num_values_max && strlen(v) > 0) {
      const float value = ExUtilGetFloat(v, &parse_error);
      if (parse_error) break;
      values[n++] = value;
      v = strchr(v, separator);
      if (v == NULL) break;
      v += 1;  // skip the separator
    }
  }
  if (error != nullptr) *error = *error || parse_error;
  return parse_error ? 0 : n;
}

//------------------------------------------------------------------------------
