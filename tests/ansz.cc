// Copyright 2018 Google LLC
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//      http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.
// -----------------------------------------------------------------------------
//
//  compressor/decompressor(s) based on ANS
//
//  ./ansz input [-o output.ansz] [-check] [-o0]
//  ./ansz output.ansz -o verif
//  diff input verif && echo "OK!"
//
// Author: skal@google.com (Pascal Massimino)

#include <cstdio>
#include <cstring>
#include <string>
#include <vector>

#include "imageio/imageio_util.h"
#include "src/dsp/dsp.h"
#include "src/utils/ans.h"
#include "src/utils/ans_utils.h"

namespace WP2 {
namespace {

enum CompressionMethod { ORDER_0, ORDER_N };

////////////////////////////////////////////////////////////////////////////////
// simple varint coding

void PutSize(size_t size, ANSEnc* const enc) {
  while (size >= 128) {
    enc->PutUValue((size & 127) | 128, 8, "uvalue");
    size >>= 7;
  }
  enc->PutUValue(size, 8, "uvalue");
}

size_t ReadSize(ANSDec* const dec) {
  size_t size = 0;
  for (int n = 0; n < 32; n += 7) {
    const size_t v = dec->ReadUValue(8, "uvalue");
    size |= (v & 0x7f) << n;
    if (v < 128) break;
  }
  return size;
}

////////////////////////////////////////////////////////////////////////////////
// Sliding context order-N compression

// Class for storing / updating the context tables.
struct Ctx {
  Ctx() {
    // TODO(skal): add some options to play with adaptation methods and speed.
    for (auto& c : ctx1) {
      c.InitFromUniform(16);
      c.SetAdaptationSpeed(ANSAdaptiveSymbol::Method::kAOM);
    }
    for (auto& c : ctx2) {
      c.InitFromUniform(16);
      c.SetAdaptationSpeed(ANSAdaptiveSymbol::Method::kAOM);
    }
    ctx = 0x00;
  }
  ANSAdaptiveSymbol* Ctx1() { return &ctx1[ctx]; }
  ANSAdaptiveSymbol* Ctx2() { return &ctx2[ctx]; }
  // we keep as context #CtxBits bits previously seen
  void Update(uint8_t nib) {
    assert(nib < 16);
    ctx = ((ctx << 4) | nib) & (kCtxSize - 1);
  }

 private:
  static constexpr uint32_t kCtxBits = 8;   // number of context bits
  static constexpr uint32_t kCtxSize = 1u << kCtxBits;
  // we use two sets of contexts for hi and lo nibbles
  ANSAdaptiveSymbol ctx1[kCtxSize], ctx2[kCtxSize];
  uint32_t ctx;
};

////////////////////////////////////////////////////////////////////////////////

static bool CompressN(const std::string& in, ANSEnc* const enc) {
  std::unique_ptr<Ctx> ctx(new Ctx);
  for (const uint8_t s : in) {
    const uint32_t nib1 = s >> 4, nib2 = s & 0xf;
    enc->PutASymbol(nib1, ctx->Ctx1(), "nib1");
    ctx->Update(nib1);
    enc->PutASymbol(nib2, ctx->Ctx2(), "nib2");
    ctx->Update(nib2);
  }
  return true;
}

static bool DecompressN(ANSDec* const dec, std::string* const out) {
  std::unique_ptr<Ctx> ctx(new Ctx);
  for (auto& s : *out) {
    const uint8_t nib1 = dec->ReadASymbol(ctx->Ctx1(), "nib1");
    ctx->Update(nib1);
    const uint8_t nib2 = dec->ReadASymbol(ctx->Ctx2(), "nib2");
    ctx->Update(nib2);
    s = (nib1 << 4) | nib2;
  }
  return (dec->GetStatus() == WP2_STATUS_OK);
}

////////////////////////////////////////////////////////////////////////////////
// Order-0 compression

const uint32_t kMaxFreq = kANSMaxRange;

static bool Compress0(const std::string& in, ANSEnc* const enc) {
  ANSDictionaries dicts;
  if (dicts.Add(256u) != WP2_STATUS_OK) return false;
  auto* const dict = dicts.back();
  Vector_u32 counts;
  if (!counts.resize(256)) return false;
  for (auto& c : counts) c = 0;
  for (const uint32_t s : in) {
    dict->RecordSymbol(s & 0xff);
    ++counts[s & 0xff];
  }

  // Quantize to fit in IO_BITS bits.
  if (ANSCountsQuantize(false, kMaxFreq, counts.size(), &counts[0], nullptr)) {
    if (dict->SetQuantizedCounts(counts) != WP2_STATUS_OK) return false;
  }
  std::vector<OptimizeArrayStorageStat> stats(counts.size());
  StoreVector(counts.data(), counts.size(), kMaxFreq, &stats[0], enc);
  if (dict->ToCodingTable() != WP2_STATUS_OK) {
    fprintf(stderr, "Error during ToCodingTable()\n");
    return false;
  }
  for (const auto& s : in) enc->PutSymbol(s & 0xff, *dict, "id");
  return true;
}

static bool Decompress0(ANSDec* const dec, std::string* const out) {
  VectorNoCtor<ANSSymbolInfo> codes;
  {
    Vector_u32 counts;
    if (!counts.resize(256)) return false;
    ReadVector(dec, kMaxFreq, counts);
    if (ANSCountsToSpreadTable(
        counts.data(), 256, ANS_LOG_TAB_SIZE, codes) != WP2_STATUS_OK) {
      return false;
    }
  }
  for (auto& s : *out) {
    s = dec->ReadSymbol(codes.data(), ANS_LOG_TAB_SIZE, "id");
  }
  return (dec->GetStatus() == WP2_STATUS_OK);
}

////////////////////////////////////////////////////////////////////////////////
// Generic calls

static bool Compress(const std::string& in, std::string* const out,
                     CompressionMethod method) {
  ANSEnc enc;
  PutSize(in.size(), &enc);
  enc.PutUValue(method == ORDER_0 ? 0 : 1, 1, "method");
  const bool ok = (method == 0) ? Compress0(in, &enc) : CompressN(in, &enc);
  if (!ok) return false;
  if (enc.Assemble() != WP2_STATUS_OK) return false;
  *out = std::string((const char*)enc.Buffer(), enc.BufferSize());
  return true;
}

static bool Decompress(const std::string& in, std::string* const out) {
  ExternalDataSource data_source((const uint8_t*)in.data(), in.size());
  ANSDec dec(&data_source);
  out->resize(ReadSize(&dec));
  const uint32_t method = dec.ReadUValue(1, "method");
  return (method == 0) ? Decompress0(&dec, out) : DecompressN(&dec, out);
}

}  // namespace
}  // namespace WP2

////////////////////////////////////////////////////////////////////////////////

int main(int argc, const char* argv[]) {
  WP2::ANSInit();
  bool decomp = false;
  bool check = false;
  WP2::CompressionMethod method = WP2::ORDER_N;   // default:
  const char *in_file = NULL, *out_file = NULL;
  for (int c = 1; c < argc; ++c) {
    if (!strcmp(argv[c], "-d")) {
      decomp = true;
    } else if (!strcmp(argv[c], "-check")) {
      check = true;
    } else if (!strcmp(argv[c], "-o0")) {
      method = WP2::ORDER_0;
    } else if (!strcmp(argv[c], "-o") && c + 1 < argc) {
      out_file = argv[++c];
    } else if (!strcmp(argv[c], "-h")) {
      printf("Usage %s in_file [-d] [-check] [-o0] [-o out_file]\n",
             argv[0]);
      return 0;
    } else {
      in_file = argv[c];
    }
  }
  if (in_file == NULL) {
    fprintf(stderr, "Missing filename!\n");
    return -1;
  }
  std::string in, out;
  if (WP2::IoUtilReadFile(in_file, &in) != WP2_STATUS_OK) {
    fprintf(stderr, "Error opening file: %s\n", in_file);
    return 1;
  }
  bool ok;
  if (!decomp) {
    ok = WP2::Compress(in, &out, method);
    if (ok && check) {
      std::string tmp;
      ok = WP2::Decompress(out, &tmp);
      if (!ok || tmp != in) {
        fprintf(stderr, "Validation failed!\n");
        return 2;
      } else {
        fprintf(stderr, "Check OK.\n");
      }
    }
  }
  if (decomp) ok = WP2::Decompress(in, &out);

  if (!ok) {
    fprintf(stderr, "Processing error (%s)\n",
            decomp ? "decompressing" : "compressing");
    return 2;
  }
  if (out_file != NULL) {
    if (WP2::IoUtilWriteFile(out, out_file, /*overwrite=*/true) !=
        WP2_STATUS_OK) {
      fprintf(stderr, "Error saving output file (%s, %d bytes)\n", out_file,
              (int)out.size());
      return 3;
    }
  }
  fprintf(stderr, "in: %d out:%d    (%.2f%%) %s\n",
    (int)in.size(), (int)out.size(), 100.f * out.size() / in.size(),
    in.size() < out.size() ? "***" : "");
  return 0;
}
