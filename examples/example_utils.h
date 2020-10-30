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
//  Helper functions used by examples
//
// Author: Skal (pascal.massimino@gmail.com)

#ifndef WP2_EXAMPLES_EXAMPLE_UTILS_H_
#define WP2_EXAMPLES_EXAMPLE_UTILS_H_

#include <algorithm>
#include <cmath>
#include <cstdarg>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <map>
#include <string>
#include <vector>

#include "examples/stopwatch.h"
#include "src/common/lossy/block_size.h"
#include "src/common/lossy/block_size_io.h"
#include "src/utils/ans.h"
#include "src/utils/ans_utils.h"
#include "src/wp2/base.h"
#include "src/wp2/decode.h"

//------------------------------------------------------------------------------

// Parse numerical values (and update *error in case of error).
extern int ExUtilGetInt(const char* const v, bool* const error, int base = 0);
extern uint32_t ExUtilGetUInt(const char* const v, bool* const error,
                              int base = 0);

extern bool ExUtilTryGetInt(const char* const v, int* const n, int base = 0);
extern bool ExUtilTryGetUInt(const char* const v, uint32_t* const n,
                             int base = 0);

extern float ExUtilGetFloat(const char* v, bool* const error);

// Parse several float values into values[] (at most 'num_values_max') from a
// list. Return true in case of error ('parse_error').
extern uint32_t ExUtilGetFloats(const char* v,
                                float values[], uint32_t num_values_max,
                                bool* const error, char separator = ',');

//------------------------------------------------------------------------------

// ExUtilDisableSIMD() will disable any CPU acceleration
#if !defined(WP2_DLL)
#if !defined(WP2_CPUINFO_IS_DEFINED)
#include "src/dsp/dsp.h"  // for WP2GetCPUInfo definition
#endif  // WP2_CPUINFO_IS_DEFINED
static inline void ExUtilDisableSIMD() { WP2GetCPUInfo = nullptr; }
#else
static inline void ExUtilDisableSIMD() {}
#endif  // WP2_DLL

//------------------------------------------------------------------------------
// Distortion

static const char* const kWP2MetricNames[] = {
  "PSNR", "PSNRHVS", "PSNR_YUV", "SSIM", "SSIM_YUV", "MSSSIM", "LSIM"
};
STATIC_ASSERT_ARRAY_SIZE(kWP2MetricNames, WP2::NUM_METRIC_TYPES);

// Computes and prints the distortion of the compressed image stored in 'data'
// compared to 'refs' (single image or several frames if it's an animation).
void WP2PrintDistortion(const std::vector<WP2::ArgbBuffer>& refs,
                        const uint8_t* const data, size_t data_size,
                        WP2::MetricType metric, bool short_output);

// One-liner for still images because ArgbBuffer doesn't provide a copy ctor.
void WP2PrintDistortion(const WP2::ArgbBuffer& ref, const uint8_t* const data,
                        size_t data_size, WP2::MetricType metric,
                        bool short_output);

//------------------------------------------------------------------------------
// Bit trace

#if defined(WP2_BITTRACE)
// Report bit traces
typedef std::pair<std::string, WP2BitCounts::mapped_type> WP2TraceType;
static inline void PrintBitTraceCluster(const std::vector<WP2TraceType>& values,
                                        bool sort_values, bool use_bytes,
                                        bool show_histograms,
                                        uint32_t level = 0,
                                        uint32_t level_max = 0) {
  if (level > level_max) return;
  // Cluster the data. Keep track of each cluster and its summary.
  std::map<std::string, std::vector<WP2TraceType>> clusters;
  std::vector<WP2TraceType> cluster_totals;
  double data_size = 0.;
  for (const auto& v : values) {
    const size_t i = v.first.find("/");
    if (i == std::string::npos) {
      cluster_totals.push_back(v);
      data_size += v.second.bits;
    } else {
      // Store any non-root element in the matching cluster.
      clusters[v.first.substr(0, i)].emplace_back(v.first.substr(i + 1),
                                                  v.second);
    }
  }

  // Process the biggest clusters first.
  sort(cluster_totals.begin(), cluster_totals.end(),
       [=](const WP2TraceType& a, const WP2TraceType& b) {
         return sort_values ? (a.second.bits > b.second.bits)
                            : (a.first > b.first);
       });
  const char * const byte_str = use_bytes ? "bytes" : "bits";
  const float byte_factor = use_bytes ? 8.f : 1.f;
  for (const auto& v : cluster_totals) {
    // Print cluster_type  _alignment_: cluster_size bit/bytes.
    const int len = fprintf(stderr, "%*s%s%*s %5.2f %% : %9.1f %s ",
                            2 * level, "",
                            v.first.c_str(),
                            std::max(30 - (int)v.first.size(), 0), "",
                            100. * v.second.bits / data_size,
                            v.second.bits / byte_factor,
                            byte_str);
    // Print _alignment_(% total).
    fprintf(stderr, "%*savg: %6.2f %s/call (%6d calls)",
            std::max(60 - len, 0), "",
            v.second.bits / byte_factor /
                (v.second.num_occurrences == 0 ? 1 : v.second.num_occurrences),
            byte_str,
            v.second.num_occurrences);
    // Display the histogram
    if (show_histograms && !v.second.histo.empty()) {
      // Build the histogram.
      std::vector<uint32_t> counts(v.second.histo.rbegin()->first + 1, 0u);
      for (const auto& p : v.second.histo) counts[p.first] = p.second;
      // Display the type and entropy.
      const float entropy = WP2::ANSCountsQuantizedCost(
                                counts.data(), counts.data(), counts.size()) /
                            byte_factor;
      constexpr const char* types[] = {"Bit",    "ABit",    "RValue",
                                       "UValue", "ASymbol", "Symbol"};
      fprintf(stderr, "%*s %s, entr: %.1f %s [",
              5, "",
              types[(int)v.second.type],
              entropy,
              byte_str);
      for (uint32_t i = 0; i < counts.size(); ++i) {
        if (counts[i] == 0) {
          // Count how many zeros we have.
          uint32_t j = i + 1;
          while (j < counts.size() && counts[j] == 0) ++j;
          // Shrink the number of written 0s if 8 or more are found.
          if (j - i > 7) {
            fprintf(stderr, "%d zeros", j - i);
            i = j - 1;
          } else {
            fprintf(stderr, "0");
          }
        } else {
          fprintf(stderr, "%d", counts[i]);
        }
        fprintf(stderr, "%s", (i == counts.size() - 1) ? "" : ", ");
      }
      fprintf(stderr, "]");
    }
    fprintf(stderr, "\n");
    PrintBitTraceCluster(clusters[v.first], sort_values, use_bytes,
                         show_histograms, level + 1, level_max);
  }
}

static inline void PrintBitTraces(const WP2::DecoderInfo& info,
                                  size_t data_size, bool sort_values = true,
                                  bool use_bytes = false,
                                  bool show_histograms = false,
                                  bool short_version = false,
                                  uint32_t bit_trace_level = 0) {
  fprintf(stderr, "bit traces:\n");
  std::map<std::string, WP2BitCounts::mapped_type> acc;
  for (const auto& p : info.bit_traces) {
    int level = bit_trace_level;
    if (p.first.find(WP2::ANSDec::kPrefixStr) == 0) {
      const std::string str = p.first.substr(strlen(WP2::ANSDec::kPrefixStr));
      // The number of calls for a cluster is roughly the number of times a
      // prefix was added.
      acc[str].num_occurrences += p.second.num_occurrences;
      continue;
    }
    size_t i = 0;
    while (level >= 0 && i != std::string::npos) {
      i = p.first.find("/", (i == 0) ? 0 : i + 1);
      acc[p.first.substr(0, i)].bits += p.second.bits;
      // For the leaf, we know the number of calls.
      if (i == std::string::npos) acc[p.first.substr(0, i)] = p.second;
      --level;
    }
  }
  if (!short_version) {
    std::vector<WP2TraceType> values;
    for (const auto& p : acc) values.push_back(p);
    // Print the first level cluster, each cluster printing the sub-ones
    // recursively.
    PrintBitTraceCluster(values, sort_values, use_bytes, show_histograms,
                         /*level=*/0, bit_trace_level);
    double total = 0.;
    for (const auto& p : info.bit_traces) total += p.second.bits;
    fprintf(stderr, "---- ANS bytes estimated: %u, real: %u ---------------\n",
            (uint32_t)std::ceil(total / 8.), (uint32_t)data_size);
  }
  // For the short summary, only focus on the first level.
  for (const auto& v : acc) {
    if (v.first.find("/") != std::string::npos) continue;
    fprintf(stderr, "%s : %.1f %s\n", v.first.c_str(),
            v.second.bits / (use_bytes ? 8. : 1.),
            use_bytes ? "bytes" : "bits");
  }
}
#endif  // WP2_BITTRACE

void WP2CollectBitTraces(const uint8_t* const data, size_t data_size,
                         int bit_trace, uint32_t bit_trace_level,
                         bool short_output, bool show_histograms);

//------------------------------------------------------------------------------

// Print to an std::string like sprintf() would
[[gnu::format (printf, 1, 0)]]
static inline std::string SPrintf(const char* format, ...) {
  va_list args;
  va_start(args, format);
  std::string result;
  const size_t max_len = 512;
  char buffer[max_len];
  const size_t len = (size_t)vsnprintf(buffer, max_len, format, args);
  va_end(args);
  result.append(buffer, len);
  return result;
}

// verify the status of a 'call' and exit upon error with an error-message
#define EXIT_IF_ERROR(call, error_format, ...)                              \
  do {                                                                      \
    const WP2Status local_status = (call);                                  \
    if (local_status != WP2_STATUS_OK) {                                    \
      fprintf(stderr, (error_format), ##__VA_ARGS__);                       \
      fprintf(stderr, "\nStatus: %s\n", WP2GetStatusMessage(local_status)); \
      exit(1);                                                              \
    }                                                                       \
  } while (0)

#define EXIT_IF_FALSE(call, error_format, ...)        \
  do {                                                \
    if (!(call)) {                                    \
      fprintf(stderr, (error_format), ##__VA_ARGS__); \
      fprintf(stderr, "\n");                          \
      exit(1);                                        \
    }                                                 \
  } while (0)

//------------------------------------------------------------------------------
// Tools for working with directories.

// Join 'prefix' and 'suffix' with the OS' path separator, if not already there.
std::string WP2JoinPath(const std::string& prefix, const std::string& suffix);

// Get everything after the last path separator.
std::string WP2GetFileName(const std::string& path);

// Remove everything after and the last dot, if there is no separator after it.
std::string WP2RemoveFileExtension(const std::string& path);

// Returns everything after the last dot, if there is no separator after it.
std::string WP2GetFileExtension(const std::string& path);

// Returns true if 'path' is a directory that can be opened.
bool WP2IsDirectory(const std::string& path);

// Append the paths of all files in 'dir_path' to 'files'.
bool WP2GetFilesIn(const std::string& dir_path,
                   std::vector<std::string>* const files, bool recursive);

// Combine input file name and output directory path to form an output file path
// with the given 'extension'.
std::string WP2InputFileToOutputDir(const std::string& input_file_path,
                                    const std::string& output_dir_path,
                                    const std::string& extension);

//------------------------------------------------------------------------------
// Reading / Writing partition files

// Read the 'partition' from the file located at 'file_path'. Returns true in
// case of success (or if the file can not be opened and 'read_only' is true).
bool WP2ReadPartition(const char file_path[], bool read_only,
                      std::vector<WP2::Rectangle>* const partition);

// Adapts the 'partition' by either ignoring blocks that are not part of the
// 'partition_set' or by splitting them until they do.
void WP2ConvertPartition(uint32_t pic_width, uint32_t pic_height,
                         WP2::PartitionSet partition_set, bool ignore_invalid,
                         std::vector<WP2::Rectangle>* const partition);

// Writes the 'partition' to the file located at 'file_path'.
bool WP2WritePartition(const std::vector<WP2::Rectangle>& partition,
                       const char file_path[]);

//------------------------------------------------------------------------------
// Program options.

// Class helping the formatting of program options:
//  - end of lines are automatically added
//  - text is wrapped at 80 columns
//  - options are displayed as: "arg ...... arg description"
class ProgramOptions {
 public:
  // Adds generic text to the options.
  void Add(const std::string& desc);
  // Adds an arguments and its description.
  void Add(const char arg[], const std::string& desc);
  // Adds system options -noasm, -version, -h.
  void AddSystemOptionSection();
  // Adds metric options:  -ssim, -msssim, -psnr, -lsim, -psnrhvs.
  void AddMetricOptions();
  // Prints the compiled options with the formatting described above.
  void Print() const;
  // Parses a "main" argument for a system options. Returns true if one was
  // found.
  WP2_NO_DISCARD
  static bool ParseSystemOptions(const char* const arg, bool* const must_stop);
  // Parses the "main" argument for a metric option.
  // Leaves type untouched if nothing was found. Returns true if one was found.
  WP2_NO_DISCARD
  static bool ParseMetricOptions(const char* const arg, WP2::MetricType* type);

 private:
  // Contains pairs of argument and description. If it is just a
  // description, argument/first is "".
  std::vector<std::pair<const char*, const std::string>> content_;
};

#endif /* WP2_EXAMPLES_EXAMPLE_UTILS_H_ */
