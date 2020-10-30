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
//  simple command line to compute rate-distortion curves.
//
// Author: Skal (pascal.massimino@gmail.com)

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <sstream>
#include <vector>

#include "src/utils/utils.h"   // this deep link needs to be first
#include "examples/example_utils.h"
#include "examples/stopwatch.h"
#include "imageio/image_dec.h"
#include "imageio/image_enc.h"
#include "src/utils/thread_utils.h"
#include "src/utils/vector.h"
#include "src/wp2/base.h"
#include "src/wp2/encode.h"
#include "src/wp2/decode.h"

#ifdef WP2_HAVE_WEBP
#include "imageio/imageio_util.h"
#endif  // WP2_HAVE_WEBP

#include "extras/aom_utils.h"

#ifdef WP2_HAVE_SJPEG
#include "sjpeg.h"
#endif

namespace {

typedef enum { WP2_STATS = 0, WEBP_STATS, AV1_STATS,
               JPEG_STATS, NEURAL_STATS } StatType;

StatType kAllTypes[] = {WP2_STATS, WEBP_STATS, AV1_STATS, JPEG_STATS,
                        NEURAL_STATS};
static const char* const kTypeNames[] = {"WP2", "WebP", "AV1", "JPEG",
                                         "Neural"};

static const int kDefaultSteps = 8;  // default number of steps

static const float kDefaultQMin =  0.f;
static const float kDefaultQMax = 95.f;  // don't do lossless by default

enum class OutputFormat { kText, kGnuplot, kHtml, kLumascope };

//------------------------------------------------------------------------------

struct TestParams {
 public:
  std::string in_file = "";
  float qmin = kDefaultQMin;
  float qmax = kDefaultQMax;
  uint32_t qsteps = kDefaultSteps;
  bool quiet = false;
  WP2::EncoderConfig config;
  std::vector<WP2::Rectangle> force_partition;
  float q = 0;                       // in [0..100]
  bool use_mt = true;
  uint32_t stack_size = 50 * 4096u;   // must be a multiple of system's pagesize
  bool disable_filters = false;
  bool use_sharp_yuv = false;
  WP2::MetricType metric = WP2::PSNR;
  const char* neural_path = nullptr;
  // label for the gnuplot curve
  const char* label = kTypeNames[WP2_STATS];
};

//------------------------------------------------------------------------------

static void PrintHelpAndExit(int return_code) {
  const WP2::EncoderConfig kDefaultConfig = WP2::EncoderConfig::kDefault;

  ProgramOptions opt;
  opt.Add("Usage:");
  opt.Add(" rd_curve [options] input_file");
  opt.Add("");
  opt.Add("Options:");
  opt.Add("-speed", SPrintf("quality/speed trade-off "
                            "(0 = fast .. 9 = slow, default = %d)",
                            kDefaultConfig.speed));
  opt.Add("-tile_shape <int>",
          SPrintf("tile shape (0=128, 1=256, 2=512, 3=wide, default = %d)",
                  kDefaultConfig.tile_shape));
  opt.Add("-qmin", SPrintf("starting quality (%.1f)", kDefaultQMin));
  opt.Add("-qmax", SPrintf("ending quality (%.1f)", kDefaultQMax));
  opt.Add("-qsteps", SPrintf("number of points per curve (%d)", kDefaultSteps));
  opt.Add("-html", "emit HTML report");
  opt.Add("-gnuplot",
          "emit Gnuplot report (incompatible with -html and -lumascope)");
  opt.Add("-lumascope <string>",
          "emit json for go/lumascope with the passed in path prefix for image "
          "paths (incompatible with -gnuplot and -html)");
  opt.Add("-log file label", "add an extra log file to gnuplot curve");
  opt.Add("-save",
          "save output images as out.xxx.xxx.png "
          "(implied with -html and -lumascope)");
  opt.Add("-save_folder <string>", "folder for saving images (\"./\")");
  opt.AddMetricOptions();
  opt.Add("-webp", "compute WebP rd-curve (if compiled)");
  opt.Add("-av1", "compute AV1 rd-curve (if compiled)");
  opt.Add("-jpeg", "compute JPEG rd-curve (if compiled)");
  opt.Add("-nowp2", "skip WP2 rd-curve");
  opt.Add("-csp <int>", "color space to use");
  opt.Add("-sharp_yuv", "use slower / sharper YUV conversion");
  opt.Add("-pm <int>", SPrintf("partition method (0..%d, default = %d)",
                               WP2::NUM_PARTITION_METHODS - 1,
                               kDefaultConfig.partition_method));
  opt.Add("-ps <int>",
          SPrintf("partition set (0..%d, default = %d)",
                  WP2::NUM_PARTITION_SETS - 1, kDefaultConfig.partition_set));
  opt.Add("-uv_mode <int>", "UV-mode to use (0=Adapt,1=420,2=444,3=Auto)");
  opt.Add("-segments <int>",
          SPrintf("number of segments to use (%u)", kDefaultConfig.segments));
  opt.Add("-segment_mode [mode]", "one of auto, explicit, implicit");
  opt.Add("-pass <int>", SPrintf("number of passes (%u)", kDefaultConfig.pass));
  opt.Add("-nofilter", "disable all filters during decoding");
  opt.Add("-sns <int>", "SNS strength to use (0..100)");
  opt.Add("-[no_]perceptual", "turn perceptual tuning on/off");
  opt.Add("-snap", "snap blocks to power-of-two position");
  opt.Add("-set_partition <file>", "force blocks partition from file");
  opt.Add("-diffusion <int>", "error diffusion strength (0=off..100)");
  opt.Add("-crop <x> <y> <w> <h>", "crop picture with the given rectangle");
  opt.Add("-neural <path>", "use neural with the directory graphs");
  opt.Add("-label <string>", "string to use for gnuplot label");
  opt.Add("-gray", "convert source to gray-scale");
  opt.Add("-nomt", "don't use multi-threading");
  opt.Add("-ssize <int>", "Stack size for threads");
  opt.Add("-include <comma-list>",
          "only include given bit trace prefixes when computing size for WP2");
  opt.Add("-exclude <comma-list>",
          "exclude given bit trace prefixes when computing size for WP2 (takes "
          "precedence over -include)");
  opt.Add("");
  opt.AddSystemOptionSection();
  opt.Print();
  exit(return_code);
}

}  // namespace

//------------------------------------------------------------------------------

class Task : public WP2::Worker {
 public:
  Task() noexcept {};
  void Init(const TestParams& params,
            const WP2::ArgbBuffer& ref,  // needs to out-live the worker object
            const std::string prefix, StatType type,
            const std::vector<std::string>& include = {},
            const std::vector<std::string>& exclude = {}) {
    type_ = type;
    params_ = params;
    ref_ = &ref;
    prefix_ = prefix;
    include_ = include;
    exclude_ = exclude;
  }

  // Main task: encode, measure distortion, decode.
  WP2Status Execute() override {
    WP2::ArgbBuffer decoded;
    switch (type_) {
      default:
      case WP2_STATS:
        CompressWP2(&decoded);
        break;
      case WEBP_STATS:
        CompressWebP(&decoded);
        break;
      case AV1_STATS:
        CompressAV1(&decoded);
        break;
      case JPEG_STATS:
        CompressJPEG(&decoded);
        break;
      case NEURAL_STATS:
        CompressNeural(&decoded);
        break;
    }
    // finish up
    bpp_ = 8. * size_ / (ref_->width * ref_->height);

    // compute distortions
    float values[5];
    EXIT_IF_ERROR(
        decoded.GetDistortionBlackOrWhiteBackground(*ref_, WP2::PSNR, values),
        "GetDistortionBlackOrWhiteBackground() failed.");
    psnr_ = values[4];
    if (params_.metric != WP2::PSNR) {
      EXIT_IF_ERROR(decoded.GetDistortionBlackOrWhiteBackground(
                        *ref_, params_.metric, values),
                    "GetDistortionBlackOrWhiteBackground() failed.");
      disto_ = values[4];
    }
    if (!prefix_.empty()) {
      const std::string file = prefix_ + SPrintf("%d.png", type_);
      EXIT_IF_ERROR(WP2::SaveImage(decoded, file.c_str(), /*overwrite=*/true),
                    "Cannot save output [%s]", file.c_str());
    }
    return WP2_STATUS_OK;
  }

 public:
  StatType type_;
  TestParams params_;
  const WP2::ArgbBuffer* ref_;
  std::string prefix_;

  std::vector<std::string> include_;
  std::vector<std::string> exclude_;

  // Task stats
  size_t size_ = 0;  // Size in bytes.
  float bpp_ = 0;
  double psnr_ = 0., disto_ = 0.;
  float enc_time_ = 0., dec_time_ = 0.;

 private:
  // used by the two methods below
  void DoCompressWP2(WP2::ArgbBuffer* const decoded);
  void CompressWP2(WP2::ArgbBuffer* const decoded);
  void CompressNeural(WP2::ArgbBuffer* const decoded);

  void CompressWebP(WP2::ArgbBuffer* const decoded);
  void CompressAV1(WP2::ArgbBuffer* const decoded);
  void CompressJPEG(WP2::ArgbBuffer* const decoded);
};

//------------------------------------------------------------------------------

void Task::DoCompressWP2(WP2::ArgbBuffer* const decoded) {
  // Create working copy
  WP2::ArgbBuffer src;
  EXIT_IF_ERROR(src.CopyFrom(*ref_), "Can't copy buffer");

  params_.config.thread_level = (int)params_.use_mt;
  WP2::EncoderInfo einfo;
  einfo.force_partition = params_.force_partition;
  params_.config.info = &einfo;
  EXIT_IF_FALSE(params_.config.IsValid(),
                "Error! Invalid configuration. Some parameters are erroneous.");

  // Set up the memory writer
  WP2::MemoryWriter memory_writer;

  double start_time = GetStopwatchTime();
  EXIT_IF_ERROR(WP2::Encode(src, &memory_writer, params_.config),
                "Error! Cannot encode picture as WP2");
  enc_time_ = GetStopwatchTime() - start_time;
  size_ = memory_writer.size_;
  src.Deallocate();  // let's be nice to others
  params_.config.info = nullptr;

  // decode back
  WP2::DecoderConfig dec_config;
  dec_config.thread_level = params_.use_mt ? 1 : 0;
  if (params_.disable_filters) {
    dec_config.enable_deblocking_filter = false;
    dec_config.enable_directional_filter = false;
    dec_config.enable_restoration_filter = false;
  }
  WP2::DecoderInfo dinfo;
  dec_config.info = &dinfo;
  start_time = GetStopwatchTime();
  EXIT_IF_ERROR(
      WP2::Decode(memory_writer.mem_, memory_writer.size_, decoded, dec_config),
      "FATAL: Decoding failed for quality %f!", params_.config.quality);
  dec_time_ = GetStopwatchTime() - start_time;

#if defined(WP2_BITTRACE)
  if (!include_.empty() || !exclude_.empty()) {
    float total_size = 0;
    for (const auto& pair : dinfo.bit_traces) {
      const std::string& label = pair.first;
      const float size = pair.second.bits;

      for (const std::string& prefix : include_) {
        if (label.find(prefix) == 0) {
          total_size += size;
          break;
        }
      }
      for (const std::string& prefix : exclude_) {
        if (label.find(prefix) == 0) {
          total_size -= size;
          break;
        }
      }
    }

    if (!include_.empty()) {
      size_ = 0;
    }
    size_ += total_size / 8;
  }
#endif  // WP2_BITTRACE
}

void Task::CompressWP2(WP2::ArgbBuffer* const decoded) {
  params_.config.quality = params_.q;
  params_.config.alpha_quality = params_.q;
  DoCompressWP2(decoded);
}

void Task::CompressNeural(WP2::ArgbBuffer* const decoded) {
  assert(params_.neural_path != nullptr);
  params_.config.graphdef_path = params_.neural_path;
  params_.config.use_neural_compression = 1;
  // neural-net compression only accepts integer quality in [0..9] range
  params_.config.quality = (int)(params_.q * 9 / 100);
  DoCompressWP2(decoded);
}

//------------------------------------------------------------------------------

void Task::CompressWebP(WP2::ArgbBuffer* const decoded) {
#if defined(WP2_HAVE_WEBP)
  WebPConfig config;
  EXIT_IF_FALSE(WebPConfigPreset(&config, WEBP_PRESET_DEFAULT, params_.q),
                "Can't initialize WebP. Wrong library version?");
  if (params_.use_mt) ++config.thread_level;
  config.method = (params_.config.speed * 6 + 4) / 9;
  config.segments = params_.config.segments;
  config.sns_strength = params_.config.sns;
  config.pass = params_.config.pass;
  config.use_sharp_yuv = params_.use_sharp_yuv;

  // encode
  WP2::MemoryWriter writer;
  double start_time = GetStopwatchTime();
  EXIT_IF_ERROR(WP2::CompressWebP(*ref_, config, &writer),
                "Error! Could not encode WebP.");
  enc_time_ = GetStopwatchTime() - start_time;
  size_ = writer.size_;

  // decode back
  start_time = GetStopwatchTime();
  EXIT_IF_ERROR(WP2::ReadImage(writer.mem_, writer.size_, decoded,
                               WP2::FileFormat::WEBP),
                "Can't decode back the WebP output!");
  dec_time_ = GetStopwatchTime() - start_time;
#else
  (void)decoded;
#endif  // WP2_HAVE_WEBP
}

//------------------------------------------------------------------------------

void Task::CompressAV1(WP2::ArgbBuffer* const decoded) {
  WP2::ParamsAV1 p;
  p.quality = params_.q;
  p.use_yuv444 = (params_.config.uv_mode == WP2::EncoderConfig::UVMode444);
  if (params_.config.thread_level > 1) p.threads = 4;
  p.speed = params_.config.speed;
  p.pass = (params_.config.pass > 1) ? 2 : 1;

  std::string out;
  double timing[2];
  EXIT_IF_FALSE(WP2::CompressAV1(*ref_, p, decoded, &out, timing)
                    == WP2_STATUS_OK,
                "CompressAV1() failed!");
  enc_time_ = timing[0];
  dec_time_ = timing[1];
  size_ = out.size();
}

//------------------------------------------------------------------------------

void Task::CompressJPEG(WP2::ArgbBuffer* const decoded) {
#if defined(WP2_HAVE_SJPEG)
  sjpeg::EncoderParam param;
  param.SetQuality(params_.q);
  param.use_trellis = (params_.config.speed >= 3);
  param.adaptive_quantization = (params_.config.speed >= 2);
  param.Huffman_compress = (params_.config.speed >= 1);
  param.yuv_mode =
      (params_.config.uv_mode == WP2::EncoderConfig::UVMode444)
          ? SJPEG_YUV_444
          : (params_.config.uv_mode == WP2::EncoderConfig::UVMode420)
                ? SJPEG_YUV_420
                : params_.use_sharp_yuv ? SJPEG_YUV_SHARP : SJPEG_YUV_AUTO;

  const uint32_t w = ref_->width;
  const uint32_t h = ref_->height;
  const size_t stride = 3 * w;
  // Extract RGB
  WP2::Vector_u8 rgb;
  EXIT_IF_FALSE(rgb.resize(stride * h), "Malloc failure for rgb[]");
  for (uint32_t y = 0; y < h; ++y) {
    const uint8_t* const row = (const uint8_t*)ref_->GetRow(y);
    for (uint32_t x = 0; x < w; ++x) {
      rgb[y * stride + 3 * x + 0] = row[4 * x + 1];
      rgb[y * stride + 3 * x + 1] = row[4 * x + 2];
      rgb[y * stride + 3 * x + 2] = row[4 * x + 3];
    }
  }

  // encode
  std::string out;
  double start_time = GetStopwatchTime();
  EXIT_IF_FALSE(sjpeg::Encode(rgb.data(), w, h, stride, param, &out),
                "JPEG compression failed!");
  enc_time_ = GetStopwatchTime() - start_time;
  size_ = out.size();

  // decode back
  start_time = GetStopwatchTime();
  EXIT_IF_ERROR(WP2::ReadImage((const uint8_t*)out.data(), out.size(), decoded,
                               WP2::FileFormat::JPEG),
                "Can't decode back the JPEG output!");
  dec_time_ = GetStopwatchTime() - start_time;
#else
  (void)decoded;
#endif  // WP2_HAVE_SJPEG
}

//------------------------------------------------------------------------------

// struct for plotting extra log files
struct LogFile {
  const char* const name;
  const char* const label;
};

class OutputWriter {
 public:
  OutputWriter(const TestParams& params, const std::vector<StatType>& types)
      : params_(params), types_(types) {}
  virtual ~OutputWriter() {}

  virtual void Start() {}
  virtual void StartStep(uint32_t step, const Task* const task) {}
  virtual void TaskStats(StatType type, const Task* const task, bool first) {}
  virtual void EndStep(uint32_t step) { report_ += "\n"; }
  virtual void End() {}

  const std::string& GetOutput() const { return report_; }

 protected:
  bool DoType(StatType type) const {
    for (StatType t : types_) {
      if (t == type) return true;
    }
    return false;
  }

  const TestParams params_;
  const std::vector<StatType> types_;
  std::string report_;
};

class TextOutputWriter : public OutputWriter {
 public:
  TextOutputWriter(const TestParams& params, const std::vector<StatType>& types)
      : OutputWriter(params, types) {}

  void Start() override {
    report_ +=
        "# Q       {size (bytes), bpp, psnr (dB), SSIM*, "
        "enc-time (ms), dec-time (ms)}\n";
    report_ += "#     \t";
    for (auto type : kAllTypes) {
      if (DoType(type)) {
        report_ += SPrintf("|            %-10s             ", kTypeNames[type]);
      }
    }
    report_ += "\n";
  }

  void StartStep(uint32_t step, const Task* const task) override {
    report_ += SPrintf("%.1f  \t", task->params_.q);
  }

  void TaskStats(StatType type, const Task* const task, bool first) override {
    report_ +=
        SPrintf("   %6d %.3f %2.2lf %2.2lf %.2f %.2f", task->size_, task->bpp_,
                task->psnr_, task->disto_, task->enc_time_, task->dec_time_);
  }

  void End() override {
    report_ += SPrintf("# params: qmin:%.1f qmax:%.1f speed:%d\n\n",
                       params_.qmin, params_.qmax, params_.config.speed);
  }
};

class HtmlOutputWriter : public OutputWriter {
 public:
  HtmlOutputWriter(const TestParams& params, const std::vector<StatType>& types)
      : OutputWriter(params, types) {}

  void Start() override {
    report_ += "<html>\n<head><script>\n";
    const int max_type = types_.size();
    report_ += SPrintf("var max_type = %d;\n", max_type);
    report_ += "var types = [";
    for (auto type : kAllTypes) {
      if (DoType(type)) {
        report_ += SPrintf(" %d, ", type);
      }
    }
    report_ += " ];\n";
    report_ += "function hide(q) {\n";
    report_ += "  for (var t = 0; t < max_type; ++t) {\n";
    report_ += "    var id = 'img-' + q + '-' + types[t];\n";
    report_ += "    document.getElementById(id).innerHTML = \"\";\n";
    report_ += "  }\n";
    report_ += "}\n";
    report_ += "function go(q) {\n";
    report_ += "  for (var t = 0; t < max_type; ++t) {\n";
    report_ += "    var id = 'img-' + q + '-' + types[t];\n";
    report_ +=
        "    var html = \"<img src='out.\" + q + \".\""
        " + types[t] + \".png' ";
    report_ += "onclick='hide(\\\"\" + q + \"\\\");'>\";\n";
    report_ += "    document.getElementById(id).innerHTML = html;\n";
    report_ += "  }\n";
    report_ += "}\n";
    report_ += "</script></head>\n";
    report_ += SPrintf("<title>RD-curve for file %s</title>\n",
                       params_.in_file.c_str());
    report_ += SPrintf("<body>\n<center>RD-curve for file %s</center>\n",
                       params_.in_file.c_str());
    report_ += "<table border='1px'>\n";
    report_ += "<tr><td>Q</td>";
    for (auto type : kAllTypes) {
      if (DoType(type)) {
        report_ += SPrintf("<td colspan='6'>%s</td>", kTypeNames[type]);
      }
    }
    report_ += "</tr>";
    report_ += "\n";
  }

  void StartStep(uint32_t step, const Task* const task) override {
    report_ += "<tr><td>";
    report_ += SPrintf("%.1f  \t", task->params_.q);
  }

  void TaskStats(StatType type, const Task* const task, bool first) override {
    report_ += SPrintf(
        "<td>%6d</td><td>%.3f</td><td>%2.2lf</td><td>%2.2lf</td><td>%.2f</td>"
        "<td>%.2f</td>",
        task->size_, task->bpp_, task->psnr_, task->disto_, task->enc_time_,
        task->dec_time_);
  }

  void EndStep(uint32_t step) override {
    report_ += "</tr>\n";
    report_ += SPrintf(
        "</tr>\n<tr><td><input type='button' value='show' "
        "onclick='go(\"%.3d\");'></td>\n",
        step);
    for (auto type : kAllTypes) {
      if (DoType(type)) {
        report_ += SPrintf(
            "<td colspan='6'>"
            "<div id='img-%.3d-%d'></div></td>\n",
            step, type);
      }
    }
    report_ += "</tr>\n";
  }

  void End() override {
    report_ += "</table>\n";
    report_ += SPrintf("# params: qmin:%.1f qmax:%.1f speed:%d\n\n",
                       params_.qmin, params_.qmax, params_.config.speed);
    report_ += "</p></body></html>\n";
  }
};

class GnuplotOutputWriter : public OutputWriter {
 public:
  GnuplotOutputWriter(const TestParams& params,
                      const std::vector<StatType>& types,
                      const std::vector<LogFile> extra_log = {})
      : OutputWriter(params, types), extra_log_(extra_log) {}

  void Start() override {
    report_ += "#!/bin/sh\n";
    report_ += "cat << EOF > /tmp/cmd.tmp\n";
    report_ += "set style data linesp\n";
    report_ += SPrintf("set title 'RD-curve [%s]'\n", params_.in_file.c_str());
    report_ += "set xlabel 'bpp'\n";
    report_ += SPrintf("set ylabel '%s'\n", kWP2MetricNames[params_.metric]);
    report_ += "set terminal png size 1024,800\n";
    report_ += "set output '/tmp/rd-curve.png'\n";
    report_ += "plot ";
    const int skip = (params_.metric != WP2::PSNR) ? 2 : 1;
    int column = 3;

    std::string curves = "";
    for (auto type : kAllTypes) {
      if (DoType(type)) {
        const char* const label =
            (type == WP2_STATS) ? params_.label : kTypeNames[type];
        if (!curves.empty()) {
          curves += ", ";
        }
        curves += SPrintf("'/tmp/cmd.data' u %d:%d t '%s'", column,
                          column + skip, label);
        for (const auto& log : extra_log_) {
          curves += SPrintf(", '%s' u %d:%d t '%s / %s'", log.name, column,
                          column + skip, log.label, kTypeNames[type]);
        }
        column += 6;
      }
    }
    report_ += curves;

    report_ += "\nEOF\n";
    report_ += "cat << EOF > /tmp/cmd.data\n";
    report_ +=
        "# Q       {size (bytes), bpp, psnr (dB), SSIM*, "
        "enc-time (ms), dec-time (ms)}\n";
    report_ += "#     \t";
    for (auto type : kAllTypes) {
      if (DoType(type)) {
        report_ += SPrintf("|            %-10s             ", kTypeNames[type]);
      }
    }
    report_ += "\n";
  }

  void StartStep(uint32_t step, const Task* const task) override {
    report_ += SPrintf("%.1f  \t", task->params_.q);
  }

  void TaskStats(StatType type, const Task* const task, bool first) override {
    report_ +=
        SPrintf("   %6d %.3f %2.2lf %2.2lf %.2f %.2f", task->size_, task->bpp_,
                task->psnr_, task->disto_, task->enc_time_, task->dec_time_);
  }

  void End() override {
    report_ += SPrintf("# params: qmin:%.1f qmax:%.1f speed:%d\n\n",
                       params_.qmin, params_.qmax, params_.config.speed);
    report_ += "EOF\n";
    report_ += "gnuplot /tmp/cmd.tmp\n";
    if (!params_.quiet) report_ += "feh /tmp/rd-curve.png\n";
  }

 private:
  const std::vector<LogFile> extra_log_;
};

class LumascopeOutputWriter : public OutputWriter {
 public:
  LumascopeOutputWriter(const TestParams& params,
                        const std::vector<StatType>& types,
                        const std::string& lumascope_path)
      : OutputWriter(params, types), lumascope_path_(lumascope_path) {}

  void Start() override {
    report_ += "{\n";
    report_ += "\"version\": 2,\n";
    report_ += SPrintf("\"title\": \"RD-curve for file %s\",\n",
                       params_.in_file.c_str());
    report_ += "\"sets\": [{\n";
    report_ += "\"urls\": [\n";

    for (uint32_t step = 0; step <= params_.qsteps; ++step) {
      if (step > 0) {
        report_ += ",\n";
      }
      report_ += "[\n";

      std::string urls = "";
      for (auto type : kAllTypes) {
        if (DoType(type)) {
          if (!urls.empty()) {
            urls += ",\n";
          }
          urls += SPrintf("\"%s/out.%.3d.%d.png\"", lumascope_path_.c_str(),
                          step, type);
        }
      }
      report_ += urls;
      report_ += "\n]";
    }

    report_ += "\n],\n";  // End of URLs.
    report_ += "\"topTexts\": [\n";
  }

  void StartStep(uint32_t step, const Task* const task) override {
    if (step == 0) {
      report_ += "null,\n";   // in version 2, it's the default label
    } else {
      report_ += "],\n";
    }
    report_ += "[";
  }

  void TaskStats(StatType type, const Task* const task, bool first) override {
    if (!first) {
      report_ += ",\n";
    }
    const std::string names[] = {"WP2", "WEBP", "AV1", "JPEG", "NEURAL"};
    report_ += SPrintf("\"%s ; Size: %d ; bpp: %.3f ; PSNR : %2.2lf ;",
                       names[type].c_str(), task->size_,
                       task->bpp_, task->psnr_);
    if (params_.metric != WP2::PSNR) {
      report_ += SPrintf(" %s: %2.2lf ;", kWP2MetricNames[params_.metric],
                         task->disto_);
    }
    report_ += SPrintf(" enc: %.2f ; dec: %.2f\"",
                       task->enc_time_, task->dec_time_);
  }

  void End() override { report_ += "]\n]\n}]}\n"; }

 private:
  const std::string lumascope_path_;
};

static std::vector<std::string> Split(const std::string& s, char delim) {
  std::vector<std::string> result;
  std::stringstream ss(s);
  std::string item;
  while (getline(ss, item, delim)) {
    result.push_back(item);
  }
  return result;
}

//------------------------------------------------------------------------------

static int RdCurveMain(int argc, char* argv[], std::string* output_str) {
  EXIT_IF_FALSE(WP2CheckVersion(), "Error! Library version mismatch!");

  TestParams params;
  OutputFormat output_format = OutputFormat::kText;
  bool do_webp = false;
  bool do_jpeg = false;
  bool do_wp2 = true;
  bool do_av1 = false;
  bool do_neural = false;
  bool crop = false;
  bool convert_to_gray = false;
  WP2::Rectangle crop_area;
  bool save_output = false;
  const char* save_folder = ".";
  std::vector<std::string> include;
  std::vector<std::string> exclude;
  std::string lumascope_path;
  std::vector<LogFile> extra_log;
  const char* partition_file_path = nullptr;

  if (argc == 1) PrintHelpAndExit(0);

  // some default values
  for (int c = 1; c < argc; ++c) {
    bool parse_error = false;
    if (!strcmp(argv[c], "-h") || !strcmp(argv[c], "-help")) {
      PrintHelpAndExit(0);
    } else if (!strcmp(argv[c], "-speed") && c + 1 < argc) {
      params.config.speed = ExUtilGetInt(argv[++c], &parse_error);
    } else if (!strcmp(argv[c], "-tile_shape") && c + 1 < argc) {
      params.config.tile_shape =
          (WP2::TileShape)ExUtilGetInt(argv[++c], &parse_error);
    } else if (!strcmp(argv[c], "-csp") && c + 1 < argc) {
      params.config.csp_type = (WP2::Csp)ExUtilGetInt(argv[++c], &parse_error);
    } else if (!strcmp(argv[c], "-pm") && c + 1 < argc) {
      params.config.partition_method =
          (WP2::PartitionMethod)ExUtilGetInt(argv[++c], &parse_error);
    } else if (!strcmp(argv[c], "-ps") && c + 1 < argc) {
      params.config.partition_set =
          (WP2::PartitionSet)ExUtilGetInt(argv[++c], &parse_error);
    } else if (!strcmp(argv[c], "-snap")) {
      params.config.partition_snapping = true;
    } else if (!strcmp(argv[c], "-uv_mode") && c + 1 < argc) {
      params.config.uv_mode =
          (WP2::EncoderConfig::UVMode)ExUtilGetInt(argv[++c], &parse_error);
    } else if (!strcmp(argv[c], "-neural") && c + 1 < argc) {
      params.neural_path = argv[++c];
      do_neural = true;
    } else if (!strcmp(argv[c], "-label") && c + 1 < argc) {
      params.label = argv[++c];
    } else if (!strcmp(argv[c], "-crop") && c + 4 < argc) {
      crop = true;
      crop_area.x = ExUtilGetUInt(argv[++c], &parse_error);
      crop_area.y = ExUtilGetUInt(argv[++c], &parse_error);
      crop_area.width = ExUtilGetUInt(argv[++c], &parse_error);
      crop_area.height = ExUtilGetUInt(argv[++c], &parse_error);
    } else if (!strcmp(argv[c], "-nomt")) {
      params.use_mt = false;
    } else if (!strcmp(argv[c], "-ssize") && c + 1 < argc) {
      params.stack_size = ExUtilGetInt(argv[++c], &parse_error);
    } else if (!strcmp(argv[c], "-sharp_yuv")) {
      params.use_sharp_yuv = true;
    } else if (!strcmp(argv[c], "-gray")) {
      convert_to_gray = true;
    } else if (!strcmp(argv[c], "-perceptual")) {
      params.config.tune_perceptual = true;
    } else if (!strcmp(argv[c], "-no_perceptual")) {
      params.config.tune_perceptual = false;
    } else if (!strcmp(argv[c], "-sns") && c + 1 < argc) {
      params.config.sns = ExUtilGetInt(argv[++c], &parse_error);
    } else if (!strcmp(argv[c], "-segment_mode") && c + 1 < argc) {
      ++c;
      if (!strcmp(argv[c], "auto")) {
        params.config.segment_id_mode = WP2::EncoderConfig::SEGMENT_ID_AUTO;
      } else if (!strcmp(argv[c], "explicit")) {
        params.config.segment_id_mode = WP2::EncoderConfig::SEGMENT_ID_EXPLICIT;
      } else if (!strcmp(argv[c], "implicit")) {
        params.config.segment_id_mode = WP2::EncoderConfig::SEGMENT_ID_IMPLICIT;
      } else {
        fprintf(stderr, "unsupported segment mode '%s'\n", argv[c]);
        parse_error = true;
      }
    } else if (!strcmp(argv[c], "-diffusion") && c + 1 < argc) {
      params.config.error_diffusion = ExUtilGetInt(argv[++c], &parse_error);
    } else if (!strcmp(argv[c], "-segments") && c + 1 < argc) {
      params.config.segments = ExUtilGetUInt(argv[++c], &parse_error);
    } else if (!strcmp(argv[c], "-nofilter")) {
      params.disable_filters = true;
    } else if (!strcmp(argv[c], "-pass") && c + 1 < argc) {
      params.config.pass = ExUtilGetUInt(argv[++c], &parse_error);
    } else if (!strcmp(argv[c], "-qmin") && c + 1 < argc) {
      params.qmin = ExUtilGetFloat(argv[++c], &parse_error);
    } else if (!strcmp(argv[c], "-qmax") && c + 1 < argc) {
      params.qmax = ExUtilGetFloat(argv[++c], &parse_error);
    } else if (!strcmp(argv[c], "-qsteps") && c + 1 < argc) {
      params.qsteps = ExUtilGetUInt(argv[++c], &parse_error);
    } else if (!strcmp(argv[c], "-set_partition") && c + 1 < argc) {
      partition_file_path = argv[++c];
    } else if (!strcmp(argv[c], "-html")) {
      output_format = OutputFormat::kHtml;
      save_output = true;
    } else if (!strcmp(argv[c], "-gnuplot")) {
      output_format = OutputFormat::kGnuplot;
    } else if (!strcmp(argv[c], "-log") && c + 2 < argc) {
      extra_log.push_back({argv[c + 1], argv[c + 2]});
      c += 2;
    } else if (!strcmp(argv[c], "-lumascope") && c + 1 < argc) {
      output_format = OutputFormat::kLumascope;
      lumascope_path = argv[++c];
      save_output = true;
    } else if (!strcmp(argv[c], "-webp")) {
      do_webp = true;
    } else if (!strcmp(argv[c], "-jpeg")) {
      do_jpeg = true;
    } else if (!strcmp(argv[c], "-av1")) {
      do_av1 = true;
    } else if (!strcmp(argv[c], "-nowp2")) {
      do_wp2 = false;
    } else if (!strcmp(argv[c], "-save")) {
      save_output = true;
    } else if (!strcmp(argv[c], "-quiet")) {
      params.quiet = true;
    } else if (!strcmp(argv[c], "-save_folder") && c + 1 < argc) {
      save_folder = argv[++c];
    } else if (!strcmp(argv[c], "-exclude") && c + 1 < argc) {
      exclude = Split(argv[++c], ',');
#if !defined(WP2_BITTRACE)
      fprintf(stderr,
              "Error: -exclude requires the WP2_BITTRACE compile flag!\n");
      return 1;
#endif
    } else if (!strcmp(argv[c], "-include") && c + 1 < argc) {
      include = Split(argv[++c], ',');
#if !defined(WP2_BITTRACE)
      fprintf(stderr,
              "Error: -include requires the WP2_BITTRACE compile flag!\n");
      return 1;
#endif
    } else if (argv[c][0] == '-') {
      bool must_stop;
      if (ProgramOptions::ParseSystemOptions(argv[c], &must_stop)) {
        if (must_stop) return 0;
      } else if (ProgramOptions::ParseMetricOptions(argv[c], &params.metric)) {
      } else {
        fprintf(stderr, "Error! Unknown option '%s'\n", argv[c]);
        PrintHelpAndExit(-1);
      }
    } else {
      params.in_file = argv[c];
    }

    if (parse_error) PrintHelpAndExit(-1);
  }

  std::vector<StatType> types;
  if (do_wp2) types.push_back(WP2_STATS);
  if (do_webp) types.push_back(WEBP_STATS);
  if (do_av1) types.push_back(AV1_STATS);
  if (do_jpeg) types.push_back(JPEG_STATS);
  if (do_neural) types.push_back(NEURAL_STATS);

  if (params.in_file.empty()) {
    fprintf(stderr, "No input file specified!\n");
    PrintHelpAndExit(1);
  }

  // Read the input
  WP2::ArgbBuffer ref;
  EXIT_IF_ERROR(WP2::ReadImage(params.in_file.c_str(), &ref),
                "Error! Cannot read input picture file '%s'",
                params.in_file.c_str());

  ref.metadata.Clear();  // remove metadata for the test
  if (crop) {
    EXIT_IF_ERROR(ref.SetView(ref, crop_area), "Cropping operation failed.");
  }
  if (convert_to_gray) {
    for (uint32_t y = 0; y < ref.height; ++y) {
      uint8_t* const row = (uint8_t*)ref.GetRow(y);
      for (uint32_t x = 0; x < ref.width; ++x) {
        const uint32_t r = row[4 * x + 1];
        const uint32_t g = row[4 * x + 2];
        const uint32_t b = row[4 * x + 3];
        const uint8_t Y = (uint8_t)(0.2126 * r + 0.7152 * g + 0.0722 * b + .5);
        row[4 * x + 1] = row[4 * x + 2] = row[4 * x + 3] = Y;
      }
    }
  }

  if (partition_file_path != nullptr) {
    EXIT_IF_FALSE(WP2ReadPartition(partition_file_path, /*read_only=*/false,
                                   &params.force_partition),
                  "Error: Unable to parse partition file");
    WP2ConvertPartition(ref.width, ref.height, params.config.partition_set,
                        /*ignore_invalid=*/true, &params.force_partition);
  }

  std::unique_ptr<OutputWriter> output;
  switch (output_format) {
    case OutputFormat::kText:
      output.reset(new TextOutputWriter(params, types));
      break;
    case OutputFormat::kGnuplot:
      output.reset(new GnuplotOutputWriter(params, types, extra_log));
      break;
    case OutputFormat::kHtml:
      output.reset(new HtmlOutputWriter(params, types));
      break;
    case OutputFormat::kLumascope:
      output.reset(new LumascopeOutputWriter(params, types, lumascope_path));
      break;
  }

  output->Start();

  // set up multi-threading
  const auto& itf = WP2::GetWorkerInterface();

  const uint32_t num_tasks =
      (do_wp2 + do_webp + do_jpeg + do_av1 + do_neural) * (params.qsteps + 1);
  WP2::Vector<Task> tasks;
  EXIT_IF_FALSE(tasks.resize(num_tasks), "Task allocation failed.");

  Task* task = &tasks[0];
  for (uint32_t step = 0; step <= params.qsteps; ++step) {
    params.q = params.qmin + step * (params.qmax - params.qmin) / params.qsteps;
    std::string prefix;
    if (save_output) prefix = SPrintf("%s/out.%.3d.", save_folder, step);
    if (do_wp2) {
      task->Init(params, ref, prefix, WP2_STATS, include, exclude);
      EXIT_IF_FALSE(itf.Reset(task, params.use_mt, params.stack_size),
                    "Bad WP2 thread initialization!.");
      ++task;
    }
    if (do_webp) {
      task->Init(params, ref, prefix, WEBP_STATS);
      EXIT_IF_FALSE(itf.Reset(task, params.use_mt,
                              std::max(params.stack_size, 131072u)),
                    "Bad WebP thread initialization!.");
      ++task;
    }
    if (do_av1) {
      task->Init(params, ref, prefix, AV1_STATS);
      EXIT_IF_FALSE(itf.Reset(task, params.use_mt,
                              std::max(params.stack_size, 524288u)),
                    "Bad AV1 thread initialization!.");
      ++task;
    }
    if (do_jpeg) {
      task->Init(params, ref, prefix, JPEG_STATS);
      EXIT_IF_FALSE(itf.Reset(task, params.use_mt, params.stack_size),
                    "Bad JPEG thread initialization!.");
      ++task;
    }
    if (do_neural) {
      task->Init(params, ref, prefix, NEURAL_STATS);
      EXIT_IF_FALSE(itf.Reset(task, params.use_mt),
                    "Bad NEURAL thread initialization!.");
      ++task;
    }
  }

  // Main loop.
  for (auto& t : tasks) itf.Launch(&t);   // launch (or execute) all tasks
  for (auto& t : tasks) {
    EXIT_IF_ERROR(itf.Sync(&t), "Bad thread state!");
    itf.End(&t);
  }

  // Collect results.
  if (!tasks.empty()) {
    task = &tasks[0];
    for (uint32_t step = 0; step <= params.qsteps; ++step) {
      output->StartStep(step, task);

      bool first = true;
      for (auto type : kAllTypes) {
        if ((type == WP2_STATS && do_wp2) || (type == WEBP_STATS && do_webp) ||
            (type == AV1_STATS && do_av1) || (type == JPEG_STATS && do_jpeg) ||
            (type == NEURAL_STATS && do_neural)) {
          output->TaskStats(type, task, first);
          ++task;
          first = false;
        }
      }

      output->EndStep(step);
    }
  }
  output->End();

  // Done.
  *output_str = output->GetOutput();
  return 0;
}

#ifndef WP2_RD_CURVE_NO_MAIN
int main(int argc, char* argv[]) {
  std::string output_str;
  const int ret = RdCurveMain(argc, argv, &output_str);
  printf("%s", output_str.c_str());
  return ret;
}
#endif

//------------------------------------------------------------------------------
