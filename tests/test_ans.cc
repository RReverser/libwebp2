// Copyright 2020 Google LLC
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

// test coding elements in ANS

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <random>
#include <vector>

#include "examples/example_utils.h"
#include "include/helpers.h"
#include "src/dsp/dsp.h"
#include "src/utils/ans.h"
#include "src/utils/ans_utils.h"
#include "src/utils/quantizer.h"
#include "src/utils/random.h"

namespace WP2 {
namespace {

constexpr uint32_t kBaseSeed = 463634;
constexpr uint32_t kTestFactor = 10;  // The higher the longer this test takes.
constexpr bool kVerbose = false;      // Additional debug info to stdout.

//------------------------------------------------------------------------------

// Draws a random bit and have its probability be 0 or PROBA_MAX N% of the time.
inline void DrawRandomBit(uint32_t N, UniformIntDistribution* const random,
                          uint32_t* const p, uint32_t* const b) {
  if ((N > 0) && (random->Get(0u, 100u) <= N)) {
    // Force to reach the bounds of probabilities.
    *p = PROBA_MAX * random->Get(0u, 1u);
  } else {
    *p = random->Get(0u, PROBA_MAX);
  }
  if (*p == 0) {
    *b = 1;
  } else if (*p == PROBA_MAX) {
    *b = 0;
  } else {
    *b = random->Get(0u, 1u);
  }
}

struct Token {
  enum Type { Bit, Bool, ASymbol, Symbol, Range, MinMax, Uniform, Signed, ABit,
             None };
  explicit Token(Type token_type = None, uint32_t token_first = 0,
                 uint32_t token_second = 0, uint32_t token_third = 0)
      : type(token_type),
        first(token_first),
        second(token_second),
        third(token_third) {}
  Type type;
  int32_t first = 0;
  int32_t second = 0;
  int32_t third = 0;
};

// Generate some random signal and append it to the ANS encoder.
// It creates 'n_elements' of type chosen randomly from the ones in
// 'used_types'.
// 'tokens' also gets updated.
void GenerateRandomEnc(size_t n_elements,
                       const std::vector<Token::Type>& used_types,
                       uint32_t max_asymbol, uint32_t n_dict,
                       uint32_t max_symbol, uint32_t max_range,
                       ANSBinSymbol* const bin_symbol,
                       ANSAdaptiveSymbol* const asymbol,
                       ANSDictionaries* const dicts,
                       UniformIntDistribution* const random,
                       std::vector<Token>* tokens, ANSEncBase* const enc) {
  const bool has_symbols = (std::find(used_types.begin(), used_types.end(),
                                      Token::Symbol) != used_types.end());
  uint32_t n_dict_ini = 0;
  if (has_symbols) {
    assert(n_dict > 0);
    n_dict_ini = dicts->size();
    for (uint32_t i = 0; i < n_dict; ++i) {
      ASSERT_WP2_OK(dicts->Add(max_symbol));
    }
  }
  const size_t size_ini = tokens->size();
  tokens->resize(size_ini + n_elements);
  // Pre-compute the signal.
  for (size_t i = size_ini; i < tokens->size(); ++i) {
    const uint32_t type_idx = random->Get(0u, (uint32_t)used_types.size() - 1u);
    const Token::Type t = used_types[type_idx];
    switch (t) {
      case Token::Bit: {
        uint32_t p, b;
        DrawRandomBit(10, random, &p, &b);
        (*tokens)[i] = Token(t, b, p);
        break;
      }
      case Token::ABit: {
        uint32_t b;
        const uint32_t thresh = random->Get(0u, 10u);
        b = (random->Get(0u, 10u) > thresh) ? 1 : 0;
        (*tokens)[i] = Token(t, b);
        break;
      }
      case Token::Bool: {
        const bool b = random->Get(0u, 1u);
        (*tokens)[i] = Token(t, b);
        break;
      }
      case Token::ASymbol: {
        const uint32_t s = random->Get(0u, max_asymbol - 1u);
        (*tokens)[i] = Token(t, s);
        break;
      }
      case Token::Symbol: {
        const uint32_t d = n_dict_ini + random->Get(0u, n_dict - 1u);
        const uint32_t s = random->Get(0u, max_symbol - 1u);
        (*dicts)[d]->RecordSymbol(s);
        (*tokens)[i] = Token(t, d, s);
        break;
      }
      case Token::Uniform:
      case Token::Signed: {
        const uint32_t U = random->Get(0u, kANSMaxUniformBits);
        const uint32_t v = random->Get(0u, (1u << U) - 1u);
        (*tokens)[i] = Token(t, v, U);
        break;
      }
      case Token::Range: {
        const uint32_t R =
            random->Get(1u, max_range ? max_range : kANSMaxRange);
        const uint32_t v = random->Get(0u, R - 1u);
        (*tokens)[i] = Token(t, v, R);
        break;
      }
      case Token::MinMax: {
        const uint32_t max =
            random->Get(0u, max_range ? max_range - 1u : kANSMaxRange);
        const uint32_t min = random->Get(0u, max);
        const uint32_t v = random->Get(min, max);
        (*tokens)[i] = Token(t, v, min, max);
        break;
      }
      default:
        assert(false);
    }
  }
  if (has_symbols) {
    EXPECT_WP2_OK(dicts->ToCodingTable());
  }
  for (size_t i = size_ini; i < tokens->size(); ++i) {
    const Token& tok = (*tokens)[i];
    switch (tok.type) {
      case Token::Bit: {
        enc->PutBit(tok.first, tok.second, "bit");
        break;
      }
      case Token::ABit: {
        enc->PutABit(tok.first, bin_symbol, "abit");
        break;
      }
      case Token::Bool: {
        enc->PutBool(tok.first, "bool");
        break;
      }
      case Token::ASymbol: {
        enc->PutASymbol(tok.second, asymbol, "asymbol");
        break;
      }
      case Token::Symbol: {
        enc->PutSymbol(tok.second, *(*dicts)[tok.first], "symbol");
        break;
      }
      case Token::Uniform: {
        enc->PutUValue(tok.first, tok.second, "U_value");
        break;
      }
      case Token::Signed: {
        enc->PutSUValue((int32_t)tok.first - ((1 << tok.second) >> 1),
                        tok.second, "S_value");
        break;
      }
      case Token::Range: {
        enc->PutRValue(tok.first, tok.second, "R_value");
        break;
      }
      case Token::MinMax: {
        enc->PutRange(tok.first, tok.second, tok.third, "minmax");
        break;
      }
      default:
        assert(false);
    }
  }
}

// Checks that an ANS buffer contains a list of tokens.
// 'codes' is the usual codes used for symbols in ANS: a vector of SymbolInfo
// per dictionary.
bool VerifyBuffer(const uint8_t* const buf, size_t size,
                  const std::vector<VectorNoCtor<ANSSymbolInfo>>* const codes,
                  ANSBinSymbol* const bin_symbol,
                  ANSAdaptiveSymbol* const asymbol,
                  const std::vector<Token>& tokens) {
  ExternalDataSource data_source(buf, size);
  ANSDec dec(&data_source);
  for (const Token& tok : tokens) {
    switch (tok.type) {
      case Token::Bit:
        if (dec.ReadBit(tok.second, "bit") != (uint32_t)tok.first) {
          printf("Bit value mismatch.");
          return false;
        }
        break;
      case Token::ABit:
        assert(bin_symbol != nullptr);
        if (dec.ReadABit(bin_symbol, "abit") != (uint32_t)tok.first) {
          printf("ABit value mismatch.");
          return false;
        }
        break;
      case Token::Bool:
        if (dec.ReadBool("bool") != (uint32_t)tok.first) {
          printf("Bool value mismatch.");
          return false;
        }
        break;
      case Token::ASymbol:
        assert(asymbol != nullptr);
        if (dec.ReadASymbol(asymbol, "asymbol") != (uint32_t)tok.second) {
          printf("ASymbol value mismatch.");
          return false;
        }
        break;
      case Token::Symbol:
        assert(codes != nullptr);
        if (dec.ReadSymbol((*codes)[tok.first].data(), ANS_LOG_TAB_SIZE,
                           "symbol") != (uint32_t)tok.second) {
          printf("Symbol value mismatch.");
          return false;
        }
        break;
      case Token::Uniform:
        if (dec.ReadUValue(tok.second, "U_value") != (uint32_t)tok.first) {
          printf("U value mismatch.");
          return false;
        }
        break;
      case Token::Signed:
        if (dec.ReadSUValue(tok.second, "S_value") + ((1 << tok.second) >> 1) !=
            tok.first) {
          printf("S value mismatch.");
          return false;
        }
        break;
      case Token::Range:
        if (dec.ReadRValue(tok.second, "R_value") != (uint32_t)tok.first) {
          printf("R value mismatch.");
          return false;
        }
        break;
      case Token::MinMax:
        if (dec.ReadRange(tok.second, tok.third, "minmax") != tok.first) {
          printf("MinMax value mismatch.");
          return false;
        }
        break;
      default:
        assert(0);
        break;
    }
  }
  return ((dec.GetStatus() == WP2_STATUS_OK) &&
          (data_source.GetNumNextBytes() == 0));
}

//------------------------------------------------------------------------------

class TestANS : public ::testing::Test {
  void SetUp() override {
    WP2DspReset();
    WP2MathInit();
    ANSInit();
  }
 public:
  static constexpr uint32_t kMaxRange = 1 << kANSMaxRangeBits;
};

// Test ANS for bits, U-value and R-values.
TEST_F(TestANS, Test0) {
  const uint32_t kNumTokens = 700 * kTestFactor;
  ANSEnc enc;
  UniformIntDistribution random(kBaseSeed);
  std::vector<Token> tokens;

  GenerateRandomEnc(
      kNumTokens, {{Token::Bit, Token::Uniform, Token::Signed, Token::Range}},
      /*max_asymbol=*/0, /*n_dict=*/0, /*max_symbol=*/0, /*max_range=*/0,
      /*bin_symbol=*/nullptr, /*asymbol=*/nullptr, /*dicts=*/nullptr, &random,
      &tokens, &enc);

  EXPECT_WP2_OK(enc.Assemble());
  const size_t size = enc.BufferSize();
  const uint8_t* buf = enc.Buffer();
  const int size_expected = (int)std::lround(enc.GetCostFull());
  if (kVerbose) {
    printf("%zu symbols -> actual bits: %d / expected bits:%d\n", tokens.size(),
           (int)size * 8, size_expected);
  }
  EXPECT_TRUE(VerifyBuffer(buf, size, /*codes=*/nullptr, /*bin_symbol=*/nullptr,
                           /*asymbol=*/nullptr, tokens));
}

//------------------------------------------------------------------------------

// Test HeaderEnc/HeaderDec.
TEST_F(TestANS, Test1) {
  const uint32_t kSize = 50 * kTestFactor;
  const uint32_t kMaxBits = 24;
  uint8_t buf[kSize];
  HeaderEnc enc(buf, kSize);
  HeaderDec dec(buf, 0);

  for (bool encoding : {true, false}) {
    UniformIntDistribution random(kBaseSeed);
    int N = 0;
    while (encoding ? enc.Ok() : dec.Ok()) {
      const uint32_t b = random.Get(1u, kMaxBits - 1u), sb = b - 1u;
      if (N + b > kSize * 8) break;

      const bool is_signed = random.FlipACoin();
      if (is_signed) {
        const int32_t v = random.Get(-(1 << sb), (1 << sb) - 1);
        if (encoding) {
          enc.PutSBits(v, b, "sbit");
        } else {
          EXPECT_EQ(dec.ReadSBits(b, "sbit"), v) << "Error at position N=" << N;
        }
      } else {
        const uint32_t v = random.Get(0u, (1u << b) - 1u);
        if (encoding) {
          enc.PutBits(v, b, "bit");
        } else {
          EXPECT_EQ(dec.ReadBits(b, "bit"), v) << "Error at position N=" << N;
        }
      }

      N += b;
    }

    if (encoding) {
      ASSERT_TRUE(enc.Ok()) << "Error encoding header.";
      dec = HeaderDec(buf, enc.Used());
    } else {
      ASSERT_TRUE(dec.Ok()) << "Error decoding header.";
    }

    if (kVerbose) {
      if (encoding) {
        printf("header encoding: ok: %d N=%d actual bytes: %d\n", (int)enc.Ok(),
               N, (int)enc.Used());
      } else {
        printf("header decoding: ok: %d N=%d\n", (int)dec.Ok(), N);
      }
    }
  }
}

//------------------------------------------------------------------------------

// Test different symbol inputs.
TEST_F(TestANS, Test2) {
  const uint32_t kNumSymbols = 700 * kTestFactor;
  UniformIntDistribution random(kBaseSeed);
  const uint32_t kMaxSymbol = 256;
  VectorNoCtor<ANSSymbolInfo> codes;
  EXPECT_TRUE(codes.resize(31));   // will be resized later

  for (size_t type : {0, 1, 2}) {
    ANSEnc enc;
    ANSDictionaries dicts;
    EXPECT_WP2_OK(dicts.Add(kMaxSymbol));
    ANSDictionary* const dict = dicts.back();
    // Generate the symbols and fill the encoder.
    std::vector<uint32_t> symbols(kNumSymbols);
    Vector_u32 counts;
    ASSERT_TRUE(counts.resize(kMaxSymbol));
    for (auto& c : counts) c = 0;
    for (uint32_t& s : symbols) {
      s = random.Get(0u, kMaxSymbol - 1u);
      ++counts[s];
      dict->RecordSymbol(s);
    }
    // Depending on the case, change counts.
    switch (type) {
      case 0: {
        // Nominal case, use the normal probabilities.
        const Vector_u32& counts0 = dict->Counts();
        EXPECT_EQ(dict->MaxSymbol(), kMaxSymbol);
        for (size_t i = 0; i < kMaxSymbol; ++i) {
          EXPECT_EQ(counts[i], counts0[i]) << "Count error at #" << i;
        }
        break;
      }
      case 1: {
        // Quantize using MaxFreq.
        const uint32_t max_freq =
            *std::max_element(std::begin(counts), std::end(counts)) / 2;
        EXPECT_GE(max_freq, 1u);
        EXPECT_TRUE(ANSCountsQuantize(false, max_freq, kMaxSymbol, &counts[0],
                                      nullptr));
        EXPECT_WP2_OK(dict->SetQuantizedCounts(counts));
        break;
      }
      case 2: {
        // Impose a pre-defined random set of probabilities.
        for (uint32_t& c : counts) {
          c = (c > 0u) ? random.Get(1u, 100u) : 0u;
        }
        EXPECT_WP2_OK(dict->SetQuantizedCounts(counts));
        break;
      }
      default:
        EXPECT_FALSE(true);
    }
    EXPECT_WP2_OK(dict->ToCodingTable());
    for (uint32_t s : symbols) {
      enc.PutSymbol(s, *dict, "id");
    }
    const float cost1 = enc.GetCost();
    const float cost2 = enc.GetCost(dicts);
    // Error is usually less than 0.1%.
    EXPECT_NEAR(cost1, cost2, 0.001 * cost1) << "type: " << type;

    EXPECT_WP2_OK(enc.Assemble());

    // Decode the buffer.
    const size_t size = enc.BufferSize();
    const uint8_t* buf = enc.Buffer();
    ExternalDataSource data_source(buf, size);
    ANSDec dec(&data_source);
    ASSERT_WP2_OK(ANSCountsToSpreadTable(&counts[0], kMaxSymbol,
                                         ANS_LOG_TAB_SIZE, codes));
    for (size_t i = 0; i < kNumSymbols; ++i) {
      const uint32_t s = dec.ReadSymbol(codes.data(), ANS_LOG_TAB_SIZE, "id");
      EXPECT_EQ(s, symbols[i]) << "Error case " << type << " at position " << i;
    }
  }
}

//------------------------------------------------------------------------------

// Test some convenience functions: Freeze, GetCounts, StoreVector, ReadVector.
TEST_F(TestANS, Test3) {
  const uint32_t kNumSymbols = 700 * kTestFactor;
  UniformIntDistribution random(kBaseSeed);
  const uint32_t kMaxSymbol = 256;

  // Generate the symbols and fill the dictionary.
  ANSEnc enc;
  ANSDictionaries dicts;
  EXPECT_WP2_OK(dicts.Add(kMaxSymbol));
  ANSDictionary* const dict = dicts.back();
  std::vector<uint32_t> symbols(kNumSymbols);
  for (uint32_t& i : symbols) {
    i = random.Get(0u, kMaxSymbol - 1u);
    dict->RecordSymbol(i);
  }

  // Write the stream.
  constexpr uint32_t vector_size = 10;
  constexpr uint32_t num_vectors = 10;
  uint32_t vector[num_vectors * vector_size];
  uint32_t nnzs[num_vectors];
  {
    const Vector_u32& counts = dict->Counts();
    const uint32_t size = dict->MaxSymbol();
    const uint32_t max_freq = *std::max_element(counts.begin(), counts.end());

    enc.PutRange(max_freq, 1, kANSMaxRange, "max_freq");
    enc.AddDebugPrefix("vector");
    std::vector<OptimizeArrayStorageStat> stats(size);
    const float cost =
        StoreVector(counts.data(), size, max_freq, stats.data(), &enc);
    if (kVerbose) printf("Storage cost: %f\n", cost);
    // Store random vectors.
    for (uint32_t max_nnz = 1; max_nnz <= num_vectors; ++max_nnz) {
      uint32_t* const sub_vector = &vector[vector_size * (max_nnz - 1)];
      std::fill(sub_vector, sub_vector + vector_size, 0u);
      for (uint32_t j = 0; j < max_nnz; ++j) {
        const uint32_t freq = random.Get(0, 100);
        sub_vector[j == max_nnz - 1 ? vector_size - 1
                                    : random.Get(0u, vector_size - 1)] = freq;
      }
      nnzs[max_nnz - 1] =
          vector_size - std::count(sub_vector, sub_vector + vector_size, 0);
      const uint32_t val_upper =
          *std::max_element(sub_vector, sub_vector + vector_size);
      StoreVectorNnz(sub_vector, vector_size, nnzs[max_nnz - 1], val_upper,
                     stats.data(), &enc);
    }
    EXPECT_WP2_OK(dict->ToCodingTable());
    for (const uint32_t& i : symbols) enc.PutSymbol(i, *dict, "id");
  }
  EXPECT_WP2_OK(enc.Assemble());

  // Decode the meta-info.
  VectorNoCtor<ANSSymbolInfo> codes;
  EXPECT_TRUE(codes.resize(535));
  const size_t size = enc.BufferSize();
  const uint8_t* buf = enc.Buffer();
  ExternalDataSource data_source(buf, size);
  ANSDec dec(&data_source);
  {
    const uint32_t max_freq = dec.ReadRange(1, kANSMaxRange, "max_freq");
    Vector_u32 counts;
    EXPECT_TRUE(counts.resize(kMaxSymbol));
    dec.AddDebugPrefix("vector");
    ReadVector(&dec, max_freq, counts);
    // Read random vectors.
    for (uint32_t max_nnz = 1; max_nnz <= num_vectors; ++max_nnz) {
      std::array<uint32_t, vector_size> vector_read;
      const uint32_t* const sub_vector = &vector[vector_size * (max_nnz - 1)];
      const uint32_t val_upper =
          *std::max_element(sub_vector, sub_vector + vector_size);
      ReadVectorNnz(&dec, nnzs[max_nnz - 1], val_upper, vector_read);
      EXPECT_TRUE(
          std::equal(sub_vector, sub_vector + vector_size, &vector_read[0]));
    }
    ASSERT_WP2_OK(ANSCountsToSpreadTable(&counts[0], kMaxSymbol,
                                         ANS_LOG_TAB_SIZE, codes));
  }

  // Decode the buffer and make sure it matches the input
  for (size_t i = 0; i < kNumSymbols; ++i) {
    const uint32_t s = dec.ReadSymbol(codes.data(), ANS_LOG_TAB_SIZE, "id");
    EXPECT_EQ(s, symbols[i]) << "Error at position #" << i;
  }
}

//------------------------------------------------------------------------------

struct TokenInfo {
  std::string name;
  Token::Type type;
  float err_threshold;  // Threshold for the estimation error, in %.
};

class TestANSCostEstimation
    : public ::testing::TestWithParam<std::tuple<uint32_t, TokenInfo>> {};

// Test the cost estimation method
TEST_P(TestANSCostEstimation, Simple) {
// Disable this test with WP2_ENC_DEC_MATCH as it considerably increases the
// size of the stream (by adding hashes) and therefore makes the thresholds
// useless.
#if !defined(WP2_ENC_DEC_MATCH)
  const uint32_t kMaxSymbol = 256;
  const uint32_t kMaxASymbol = 10;
  const uint32_t kNDict = 3;
  const uint32_t kNValues = 40000 * kTestFactor;
  const TokenInfo& info = std::get<1>(GetParam());

  WP2MathInit();
  ANSInit();

  int size = 0;
  ANSEnc enc0;
  ANSEncCounter enc1;
  for (ANSEncBase* const enc : {(ANSEncBase*)&enc0, (ANSEncBase*)&enc1}) {
    UniformIntDistribution random(/*seed=*/std::get<0>(GetParam()));
    ANSDictionaries sym_dicts;
    ANSBinSymbol bin_symbol;
    ANSAdaptiveSymbol asymbol;
    asymbol.SetAdaptationSpeed(ANSAdaptiveSymbol::Method::kAOM);
    asymbol.InitFromUniform(APROBA_MAX_SYMBOL);
    std::vector<Token> tokens;

    const uint32_t n_dict = (info.type == Token::Symbol) ? kNDict : 0;
    const uint32_t max_symbol = (info.type == Token::Symbol) ? kMaxSymbol : 0;
    const uint32_t max_asymbol =
        (info.type == Token::ASymbol) ? kMaxASymbol : 0;
    const uint32_t max_range = 0;
    ANSBinSymbol* const bin_symbol_stats =
        (info.type == Token::ABit) ? &bin_symbol : nullptr;
    ANSAdaptiveSymbol* const asymbol_stats =
        (info.type == Token::ASymbol) ? &asymbol : nullptr;
    ANSDictionaries* const dicts =
        (info.type == Token::Symbol) ? &sym_dicts : nullptr;

    GenerateRandomEnc(kNValues, {info.type}, max_asymbol, n_dict, max_symbol,
                      max_range, bin_symbol_stats, asymbol_stats, dicts,
                      &random, &tokens, enc);

    if (enc == &enc0) {
      EXPECT_WP2_OK(enc0.Assemble());
      size = (int)(8 * enc0.BufferSize());  // store for later comparison
      if (info.type != Token::Symbol) {     // omit the dictionary case for now
        ANSBinSymbol new_bin_symbol;
        ANSAdaptiveSymbol new_asymbol;
        new_asymbol.SetAdaptationSpeed(ANSAdaptiveSymbol::Method::kAOM);
        new_asymbol.InitFromUniform(APROBA_MAX_SYMBOL);
        EXPECT_TRUE(VerifyBuffer(enc0.Buffer(), enc0.BufferSize(),
                                 /*codes=*/nullptr, &new_bin_symbol,
                                 &new_asymbol, tokens));
      }
    }

    const int size_estimated = (int)std::lround(enc->GetCostFull(sym_dicts));
    const float error = 100.f * size / size_estimated - 100.f;
    EXPECT_LE(fabs((double)error), info.err_threshold);
    if (kVerbose) {
      printf("%s: actual bits: %d / estimated bits:%d. %1.3f%% bigger.\n",
             info.name.c_str(), size, size_estimated, error);
    }
  }
#endif  // WP2_ENC_DEC_MATCH
}

INSTANTIATE_TEST_SUITE_P(
    TestANSCostEstimationInstantiation, TestANSCostEstimation,
    ::testing::Combine(::testing::Values(kBaseSeed),
                       ::testing::ValuesIn(std::vector<TokenInfo>{
                           {"Bit    ", Token::Bit, 0.15f},
                           {"Bool   ", Token::Bool, 0.01f},
                           {"ASymbol", Token::ASymbol, 2.0f},
                           {"Symbol ", Token::Symbol, 0.02f},
                           {"Range  ", Token::Range, 0.05f},
                           {"MinMax ", Token::MinMax, 0.06f},
                           {"Uniform", Token::Uniform, 0.01f},
                           {"Signed ", Token::Signed, 0.01f},
                           {"ABit   ", Token::ABit, 0.5f}})));

//------------------------------------------------------------------------------

// Test the Quantizer.
TEST_F(TestANS, Test5) {
  UniformIntDistribution random(/*seed=*/kBaseSeed);
  const uint32_t kMaxSymbol = 256;
  const uint32_t kNumSymbols = 7000 * kTestFactor;

  // Create some stream with symbols.
  ANSEnc enc;
  ANSDictionaries dicts;
  EXPECT_WP2_OK(dicts.Add(kMaxSymbol));
  ANSDictionary* const dict = dicts.back();
  for (uint32_t i = 0; i < kNumSymbols; ++i) {
    dict->RecordSymbol(random.Get(0u, kMaxSymbol - 1u));
  }

  // Compute the sparse histogram.
  const Vector_u32& counts = dict->Counts();
  std::vector<uint32_t> histogram(kMaxSymbol);
  std::vector<uint16_t> mapping(kMaxSymbol);
  uint32_t size_sparse = 0;
  for (size_t i = 0; i < dict->MaxSymbol(); ++i) {
    if (counts[i] == 0) continue;
    mapping[size_sparse] = (uint16_t)i;
    histogram[size_sparse++] = counts[i];
  }

  // Optimize the storage of proba + symbols.
  Quantizer quantizer;
  EXPECT_WP2_OK(quantizer.Allocate(kMaxSymbol));

  Quantizer::Config* config = nullptr;
  quantizer.Quantize(histogram.data(), mapping.data(), size_sparse, kMaxSymbol,
                     /*max_count=*/16, /*speed=*/5, &config);

  // Display results.
  if (kVerbose) {
    printf("Cost: %f\n", config->cost_);
    if (config->next_ != nullptr) {
      printf("Histogram is recursive:\n");
    }
    while (config != nullptr) {
      if (config->param_.is_sparse_) {
        printf("Histogram is sparse:\n(index, count - 1):\n");
        for (size_t i = 0; i < config->size_to_write_; ++i) {
          printf("(%d,%d) ", mapping[i], config->histogram_to_write_[i]);
        }
        printf("\n");
      } else {
        printf("Histogram is not sparse:\n");
        for (size_t i = 0; i < config->size_to_write_; ++i) {
          printf("%d ", config->histogram_to_write_[i]);
        }
        printf("\n");
      }
      if (config->param_.type_ == Quantizer::ConfigType::Huffman) {
        printf("The counts above should be interpreted as exponents of 2.\n");
      }
      config = config->next_;
    }
  }
}

//------------------------------------------------------------------------------

/**
 * Test dedicated to encoding of integer in a range.
 */
TEST_F(TestANS, Test6a) {
  constexpr uint32_t kNumTokens = 700 * kTestFactor;
  // Verify behavior with randomized inputs.
  for (Token::Type type : {Token::Range, Token::MinMax}) {
    UniformIntDistribution random(kBaseSeed);
    ANSEnc encoder;
    std::vector<Token> tokens;

    GenerateRandomEnc(kNumTokens, {type}, /*max_asymbol=*/0, /*n_dict=*/0,
                      /*max_symbol=*/0, kMaxRange, /*bin_symbol=*/nullptr,
                      /*asymbol=*/nullptr, /*dicts=*/nullptr, &random, &tokens,
                      &encoder);
    EXPECT_WP2_OK(encoder.Assemble());
    EXPECT_TRUE(VerifyBuffer(encoder.Buffer(), encoder.BufferSize(),
                             /*codes=*/nullptr, /*bin_symbol=*/nullptr,
                             /*asymbol=*/nullptr, tokens));
  }
}

TEST_F(TestANS, Test6b) {
  // Exhaustive verification for small integer range going from 1 to
  // max_tested_range.
  const uint32_t max_tested_range = 256;
  ANSEnc encoder;
  std::vector<Token> tokens;
  tokens.reserve(max_tested_range * (max_tested_range - 1) / 2);
  for (uint16_t range = 1; range < max_tested_range; ++range) {
    for (uint16_t value = 0; value < range; ++value) {
      encoder.PutRValue(value, range, "R_value");
      tokens.emplace_back(Token::Range, value, range);
    }
  }
  EXPECT_WP2_OK(encoder.Assemble());
  EXPECT_TRUE(VerifyBuffer(encoder.Buffer(), encoder.BufferSize(),
                           /*codes=*/nullptr, /*bin_symbol=*/nullptr,
                           /*asymbol=*/nullptr, tokens));
}

TEST_F(TestANS, Test6c) {
  // Verification for large integer ranges going from min_tested_range to
  // kMaxRange.
  const uint32_t min_tested_range = kMaxRange - 255;
  ANSEnc encoder;
  std::vector<Token> tokens;
  tokens.reserve(kMaxRange * (kMaxRange - 1) / 2 -
                  (min_tested_range - 1) * (min_tested_range - 2) / 2);
  for (uint16_t range = min_tested_range; range < kMaxRange; ++range) {
    for (uint16_t value = range / 2; value < range; ++value) {
      encoder.PutRValue(value, range, "R_value");
      tokens.emplace_back(Token::Range, value, range);
    }
  }
  EXPECT_WP2_OK(encoder.Assemble());
  EXPECT_TRUE(VerifyBuffer(encoder.Buffer(), encoder.BufferSize(),
                           /*codes=*/nullptr, /*bin_symbol=*/nullptr,
                           /*asymbol=*/nullptr, tokens));
}

TEST_F(TestANS, Test6d) {
  // Sparse verification for min/max ranges.
  ANSEnc encoder;
  std::vector<Token> tokens;
  for (uint16_t hi = 0; hi < kMaxRange; hi += hi / 2 + 1) {
    for (uint16_t lo = 0; lo <= hi; lo += lo / 2 + 1) {
      for (uint16_t value = lo; value <= hi; ++value) {
        encoder.PutRange(value, lo, hi, "minmax");
        tokens.emplace_back(Token::MinMax, value, lo, hi);
      }
    }
  }
  EXPECT_WP2_OK(encoder.Assemble());
  EXPECT_TRUE(VerifyBuffer(encoder.Buffer(), encoder.BufferSize(),
                           /*codes=*/nullptr, /*bin_symbol=*/nullptr,
                           /*asymbol=*/nullptr, tokens));
}

//------------------------------------------------------------------------------

// Generate the codes for each dictionary, given a list of tokens.
void GetCodesFromTokens(const std::vector<uint32_t>& max_symbols,
                        const std::vector<Token>& tokens,
                        std::vector<VectorNoCtor<ANSSymbolInfo>>& codes) {
  // Compute the counts of symbols.
  const size_t n_dicts = max_symbols.size();
  std::vector<std::vector<uint32_t>> counts(n_dicts);
  for (size_t i = 0; i < n_dicts; ++i) counts[i].resize(max_symbols[i], 0);
  for (const Token& tok : tokens) {
    if (tok.type == Token::Symbol) ++counts[tok.first][tok.second];
  }

  codes.resize(n_dicts);
  for (size_t i = 0; i < n_dicts; ++i) {
    EXPECT_WP2_OK(ANSCountsToSpreadTable(
        &counts[i][0], max_symbols[i], ANS_LOG_TAB_SIZE, codes[i]));
  }
}

// Test for miscellaneous ANS functions
TEST_F(TestANS, Test7) {
  const uint32_t kNumElements = 700 * kTestFactor;
  WP2MathInit();
  ANSInit();
  // Test for adding new dictionaries after some symbols are inserted.
  {
    const std::vector<uint32_t> max_symbols = {17, 17, 13, 29, 101};
    std::vector<Token> tokens;
    ANSEnc enc;
    ANSDictionaries dicts;
    UniformIntDistribution random(kBaseSeed);
    for (uint32_t max_symbol : max_symbols) {
      GenerateRandomEnc(kNumElements, {Token::Symbol}, /*max_asymbol=*/0,
                        /*n_dict=*/1, max_symbol, /*max_range=*/0,
                        /*bin_symbol=*/nullptr, /*asymbol=*/nullptr, &dicts,
                        &random, &tokens, &enc);
    }
    std::vector<VectorNoCtor<ANSSymbolInfo>> codes;
    GetCodesFromTokens(max_symbols, tokens, codes);
    EXPECT_WP2_OK(enc.Assemble());
    EXPECT_TRUE(VerifyBuffer(enc.Buffer(), enc.BufferSize(), &codes,
                             /*bin_symbol=*/nullptr, /*asymbol=*/nullptr,
                             tokens));
  }
  // Test appending.
  {
    const std::vector<uint32_t> max_symbols[2] = {{19, 21, 3, 512},
                                                  {7, 5, 405}};
    UniformIntDistribution random(kBaseSeed);
    ANSEnc enc[2];
    ANSDictionaries dicts[2];
    std::vector<Token> tokens;
    for (size_t i = 0; i < 2; ++i) {
      for (uint32_t max_symbol : max_symbols[i]) {
        GenerateRandomEnc(
            kNumElements, {Token::Symbol}, /*max_asymbol=*/0, /*n_dict=*/1,
            max_symbol, /*max_range=*/0, /*bin_symbol=*/nullptr,
            /*asymbol=*/nullptr, &dicts[i], &random, &tokens, &enc[i]);
      }
    }
    EXPECT_WP2_OK(enc[0].Append(enc[1]));
    EXPECT_WP2_OK(dicts[0].AppendAndClear(&dicts[1]));
    EXPECT_WP2_OK(dicts[0].ToCodingTable());
    EXPECT_WP2_OK(enc[0].Assemble());
    // Update the ids in the tokens.
    size_t n_dict_0 = max_symbols[0].size();
    for (uint32_t i = kNumElements * n_dict_0; i < tokens.size(); ++i) {
      tokens[i].first += n_dict_0;
    }
    // Update the max_symbols.
    std::vector<uint32_t> max_symbols_new = max_symbols[0];
    max_symbols_new.insert(max_symbols_new.end(), max_symbols[1].begin(),
                           max_symbols[1].end());

    // Decode the buffer.
    std::vector<VectorNoCtor<ANSSymbolInfo>> codes;
    GetCodesFromTokens(max_symbols_new, tokens, codes);
    EXPECT_TRUE(VerifyBuffer(enc[0].Buffer(), enc[0].BufferSize(), &codes,
                             /*bin_symbol=*/nullptr, /*asymbol=*/nullptr,
                             tokens));
  }
}

//------------------------------------------------------------------------------

struct TokenSetInfo {
  std::string name;
  std::vector<Token::Type> types;
};

class TestANSBitCounts
    : public ::testing::TestWithParam<std::tuple<uint32_t, TokenSetInfo>> {};

// This test only works in WP2_BITTRACE mode!
TEST_P(TestANSBitCounts, Simple) {
#if defined(WP2_BITTRACE)
  const int kNumDicts = 3;
  const std::vector<uint32_t> max_symbols(kNumDicts, 256);

  WP2MathInit();
  ANSInit();

  const TokenSetInfo& info = std::get<1>(GetParam());
  ANSEnc enc;
  ANSDictionaries dicts;
  UniformIntDistribution random(/*seed=*/std::get<0>(GetParam()));
  std::vector<Token> tokens;
  GenerateRandomEnc(/*n_elements=*/1000 * kTestFactor, info.types,
                    /*max_asymbol=*/0, kNumDicts, /*max_symbol=*/256,
                    /*max_range=*/3263, /*bin_symbol=*/nullptr,
                    /*asymbol=*/nullptr, &dicts, &random, &tokens, &enc);
  std::vector<VectorNoCtor<ANSSymbolInfo>> codes(kNumDicts);
  if (std::find(info.types.begin(), info.types.end(), Token::Symbol) !=
      info.types.end()) {
    GetCodesFromTokens(max_symbols, tokens, codes);
  }
  EXPECT_WP2_OK(enc.Assemble());
  const size_t size = enc.BufferSize();
  const uint8_t* buf = enc.Buffer();
  const int size_expected = (int)std::lround(enc.GetCostFull(dicts));
  ExternalDataSource data_source(buf, size);
  ANSDec dec(&data_source);
  for (const Token& tok : tokens) {
    if (tok.type == Token::Bit) {
      (void)dec.ReadBit(tok.second, "bit");
    } else if (tok.type == Token::Bool) {
      (void)dec.ReadBool("bool");
    } else if (tok.type == Token::Symbol) {
      (void)dec.ReadSymbol(codes[tok.first].data(), ANS_LOG_TAB_SIZE, "symbol");
    } else if (tok.type == Token::Uniform) {
      (void)dec.ReadUValue(tok.second, "U_value");
    } else if (tok.type == Token::Range) {
      (void)dec.ReadRValue(tok.second, "R_value");
    } else if (tok.type == Token::MinMax) {
      (void)dec.ReadRange(tok.second, tok.third, "minmax");
    } else if (tok.type == Token::Signed) {
      (void)dec.ReadSUValue(tok.second, "S_value");
    } else {
      EXPECT_FALSE(true);
    }
  }
  EXPECT_WP2_OK(dec.GetStatus());
  const double bit_count = dec.GetBitCount();

  // difference between final size and bit-trace:
  const double error = (size - bit_count / 8.) / size;
  if (kVerbose) {
    printf("%s: %d elemt. Bits: actual: %d, reported: %.2lf (%.2lf%% off)\n",
           info.name.c_str(), (int)tokens.size(), (int)size * 8, bit_count,
           error * 100);
  }
  EXPECT_LE(std::abs(error), ANSDec::kBitCountAccuracy);

  // difference between GetCost() and bit-trace:
  const double internal_error = 1. - size_expected / bit_count;
  const double kInternalErrorTolerance = 0.001;
  if (kVerbose) {
    const auto report = dec.GetBitTraces();
    for (auto it = report.begin(); it != report.end(); ++it) {
      printf("%s : %.lf bits\n", it->first.c_str(), it->second.bits);
    }
    printf("-------------\n");
  }
  EXPECT_LE(std::abs(internal_error), kInternalErrorTolerance);

#endif
}

INSTANTIATE_TEST_SUITE_P(
    TestANSBitCountsInstantiation, TestANSBitCounts,
    ::testing::Combine(::testing::Values(kBaseSeed),
                       ::testing::ValuesIn(std::vector<TokenSetInfo>{
                           {"Bit     ", {Token::Bit}},
                           {"Bool    ", {Token::Bool}},
                           {"Symbol  ", {Token::Symbol}},
                           {"Range   ", {Token::Range}},
                           {"MinMax  ", {Token::MinMax}},
                           {"Uniform ", {Token::Uniform}},
                           {"Signed  ", {Token::Signed}},
                           {"All     ",
                            {Token::Bit, Token::Symbol, Token::Uniform,
                             Token::Range}}})));

// Make sure ANS does not spin forever when stats have a tail.
TEST_F(TestANS, TestSpin) {
  ANSDictionary dict;
  EXPECT_WP2_OK(dict.Init(20));
  dict.RecordSymbol(0, 2 * ANS_TAB_SIZE);
  for (uint32_t i = 1; i < 20; ++i) dict.RecordSymbol(i, 1);
  // This should not spin.
  EXPECT_WP2_OK(dict.ToCodingTable());
}

//------------------------------------------------------------------------------

class TestANSAdaptiveSymbol : public ::testing::TestWithParam<uint32_t> {};

TEST_P(TestANSAdaptiveSymbol, Simple) {
  const uint32_t kNumSymbols = 1000 * kTestFactor;
  UniformIntDistribution random(/*seed=*/GetParam());

  WP2MathInit();
  ANSInit();

  uint32_t rnd_range = (random.Get(0, 35) == 0)
                           ? (2u << APROBA_BITS) / APROBA_MAX_SYMBOL
                           : random.Get(1u, 30u);
  uint32_t cumul[APROBA_MAX_SYMBOL + 1];
  uint32_t real_cumul[APROBA_MAX_SYMBOL] = {0};
  uint32_t sum = 0;
  for (uint32_t& c : cumul) {
    sum = c = random.Get(sum, sum + rnd_range + 1u);
    rnd_range = rnd_range * 3 / 4;
  }

  std::vector<uint8_t> syms(kNumSymbols);
  for (uint8_t& sym : syms) {
    const uint32_t p = random.Get(0u, sum);
    uint32_t s = 0;
    for (s = 0; s < APROBA_MAX_SYMBOL - 1; ++s) {
      if (p < cumul[s + 1]) break;
    }
    sym = s;
    ++real_cumul[s];
  }

  if (random.Get(0, 19) == 0) {
    std::sort(syms.begin(), syms.end());  // we'll have to be very adaptive!
  }

  // Encode
  const uint32_t adapt_speed = random.Get(0u, APROBA_BITS - 1u);
  ANSAdaptiveSymbol w_dict;
  w_dict.SetAdaptationSpeed(ANSAdaptiveSymbol::Method::kConstant, adapt_speed);
  w_dict.InitFromUniform(APROBA_MAX_SYMBOL);
  ANSEnc enc;
  enc.Clear();
  for (const uint8_t& s : syms) enc.PutASymbol(s, &w_dict, "asym");
  EXPECT_WP2_OK(enc.Assemble());
  const size_t size = enc.BufferSize();
  const uint8_t* buf = enc.Buffer();

  if (kVerbose) {
    for (uint32_t s = 0; s < APROBA_MAX_SYMBOL; ++s) {
      const float real_p = real_cumul[s] * 100.f / syms.size();
      const float dict_p = w_dict.GetProba(s) * 100.f;
      printf(
          "#%d theory:%.2f%%  real:%.2f%%  adapt-dict:%.2f%% -> error=%.1f%%\n",
          s, (cumul[s + 1] - cumul[s]) * 100.f / sum, real_p, dict_p,
          fabs(1. - dict_p / real_p) * 100.);
    }
    printf("%d symbols, size=%d\n", (int)syms.size(), (int)size);
  }

  // Decode and check
  ExternalDataSource data_source(buf, size);
  ANSDec dec(&data_source);
  ANSAdaptiveSymbol r_dict;
  r_dict.SetAdaptationSpeed(ANSAdaptiveSymbol::Method::kConstant, adapt_speed);
  r_dict.InitFromUniform(APROBA_MAX_SYMBOL);
  for (const uint8_t& sym : syms) {
    const uint32_t s = dec.ReadASymbol(&r_dict, "asym");
    EXPECT_EQ(s, sym) << ": expected symbol " << sym << ", got " << s;
  }

  // Match ending probabilities
  for (uint32_t s = 0; s < APROBA_MAX_SYMBOL; ++s) {
    EXPECT_EQ(w_dict.GetProba(s), r_dict.GetProba(s))
        << "Symbol #" << s << ": Final proba error: expected "
        << w_dict.GetProba(s) << ", got " << r_dict.GetProba(s);
  }
}

TEST_P(TestANSAdaptiveSymbol, VsDictionary) {
// Disable this test with WP2_ENC_DEC_MATCH as it considerably increases the
// size of the stream (by adding hashes) and therefore makes the thresholds
// useless.
#if !defined(WP2_ENC_DEC_MATCH)
  UniformIntDistribution random(/*seed=*/GetParam());

  // generate a pdf[]
  uint32_t len = random.Get(1000u, 1999u) * kTestFactor;
  uint32_t pdf[APROBA_MAX_SYMBOL];
  for (uint32_t& i : pdf) i = random.Get(0u, len / 16u - 1u);
  const uint32_t speed = random.Get(0u, 0x3fffu);  // *must* be < 32768 in SSE2

  // remove some symbols
  uint32_t nsym =
      (random.Get(0, 19) == 0) ? APROBA_MAX_SYMBOL : random.Get(6u, 15u);
  while (nsym++ < APROBA_MAX_SYMBOL) {
    pdf[random.Get(0u, APROBA_MAX_SYMBOL - 1u)] = 0;
  }

  len = 0;
  for (uint32_t c : pdf) len += c;

  // generate a message
  std::vector<uint8_t> msg(len);
  for (uint32_t s = 0, k = 0; s < APROBA_MAX_SYMBOL; ++s) {
    for (uint32_t n = 0; n < pdf[s]; ++n) msg[k++] = (uint8_t)s;
  }
  // shuffle the message
  Shuffle(msg.begin(), msg.end(), /*seed=*/GetParam());

  // compute the theoretical bit length
  float total = 32.;  // ~32 bits for ANS padding
  for (uint32_t v : pdf) {
    if (v) total -= log2(1. * v / len) * v;
  }

  WP2MathInit();
  ANSInit();

  // populate dictionary stats
  ANSDictionaries dicts;
  EXPECT_WP2_OK(dicts.Add(APROBA_MAX_SYMBOL));
  ANSDictionary* const dict = dicts.back();
  for (uint8_t s : msg) dict->RecordSymbol(s);

  ANSEnc enc;
  enc.Clear();
  ANSAdaptiveSymbol asym;
  asym.SetAdaptationSpeed(ANSAdaptiveSymbol::Method::kConstant, speed);
  EXPECT_WP2_OK(asym.InitFromCounts(pdf, APROBA_MAX_SYMBOL));

  for (uint8_t s : msg) enc.PutASymbol(s, &asym, "test1");
  EXPECT_WP2_OK(enc.Assemble());
  const float size1 = 8.f * enc.BufferSize();

  enc.Clear();
  EXPECT_WP2_OK(dict->ToCodingTable());
  for (uint8_t s : msg) enc.PutSymbol(s, *dict, "test2");
  EXPECT_WP2_OK(enc.Assemble());
  const float size2 = 8.f * enc.BufferSize();

  // tolerance of ASym vs Dictionary
  const float kTolerance[16] = {
      350.f, 150.f, 75.f, 50.f, 15.f, 7.2f, 3.5f, 2.8f, 1.9f, 1.9f,
      // For some reasons, starting at speed=10, it becomes more unstable.
      6.0f, 20.0f, 50.0f, 70.f, 50.f,
      // Speed=15 is really non-adaptive, though, with low error.
      0.3f};
  const float kTolerance2 = 0.2f;  // dictionary tolerance (quite low!)
  const float err1 = fabs(100. * (size1 - size2) / total);
  const float err2 = fabs(100. * (size2 - total) / total);
  const int lspeed = 15 - (int)log2(1 + speed);
  EXPECT_LT(err1, kTolerance[lspeed]);
  EXPECT_LT(err2, kTolerance2);

  if (kVerbose || err1 > kTolerance[lspeed] || err2 > kTolerance2) {
    printf("#%d symbols, expected bits=%.1f ", len, total);
    printf("final bits1 = %.1f bits2=%.1f [speed=%d lspeed=%d %.3f]\n",
           size1, size2, speed, lspeed, err1);
  }

  // Write one symbol with sym and the next with asym.
  dict->ResetCounts();
  for (uint32_t i = 0; i < msg.size(); i += 2) dict->RecordSymbol(msg[i]);
  EXPECT_WP2_OK(dict->ToCodingTable());

  enc.Clear();
  EXPECT_WP2_OK(asym.InitFromCounts(pdf, APROBA_MAX_SYMBOL));

  for (uint32_t i = 0; i < msg.size(); ++i) {
    if (i % 2 == 0) {
      enc.PutSymbol(msg[i], *dict, "sym");
    } else {
      enc.PutASymbol(msg[i], &asym, "asym");
    }
  }
  EXPECT_WP2_OK(enc.Assemble());
  const float size = 8.f * enc.BufferSize();
  const float cost = enc.GetCost(dicts);

  // tolerance of ASym vs Dictionary
  const float err = (size == 0.f) ? 0.f : fabs(100. * (size - cost) / size);
  constexpr float kTolerance3 = 1.05f;
  EXPECT_LT(err, kTolerance3);

  if (kVerbose || err > kTolerance3) {
    printf("#%d symbols, expected cost=%.1f ", len, cost);
    printf("final bits = %.1f [speed=%d lspeed=%d %.3f]\n",
           size, speed, lspeed, err);
  }
#endif  // WP2_ENC_DEC_MATCH
}

INSTANTIATE_TEST_SUITE_P(TestANSAdaptiveSymbolInstantiation,
                         TestANSAdaptiveSymbol,
                         ::testing::Range(kBaseSeed,
                                          kBaseSeed + 30 * kTestFactor));

//------------------------------------------------------------------------------

void TestUniformProbability(ANSAdaptiveSymbol* const asym,
                            uint16_t max_symbol) {
  uint32_t symbol_to_num_occurrences[APROBA_MAX_SYMBOL]{0u};
  for (uint32_t proba = 0; proba < APROBA_MAX; ++proba) {
    const uint16_t symbol = asym->GetSymbol(proba).symbol;
    EXPECT_LT(symbol, max_symbol);
    ++symbol_to_num_occurrences[symbol];
  }
  // Make sure they have almost equal probabilities.
  uint32_t num_occurrences1 = symbol_to_num_occurrences[0];
  uint32_t num_occurrences2 = symbol_to_num_occurrences[0];
  for (uint16_t symbol = 1; symbol < max_symbol; ++symbol) {
    if (symbol_to_num_occurrences[symbol] != num_occurrences1) {
      if (num_occurrences1 == num_occurrences2) {
        num_occurrences2 = symbol_to_num_occurrences[symbol];
        EXPECT_EQ(num_occurrences1 + 1, num_occurrences2);
      } else {
        EXPECT_EQ(symbol_to_num_occurrences[symbol], num_occurrences2);
      }
    }
  }
}

TEST_F(TestANS, TestANSAdaptiveSymbol) {
  ANSAdaptiveSymbol asym;
  for (uint16_t max_symbol = 1; max_symbol < APROBA_MAX_SYMBOL; ++max_symbol) {
    asym.SetAdaptationSpeed(ANSAdaptiveSymbol::Method::kAOM);
    asym.InitFromUniform(max_symbol);
    TestUniformProbability(&asym, max_symbol);
  }
}

// Tests that symbols can be initialized from CDF's.
TEST_F(TestANS, TestANSAdaptiveSymbolFromCDF) {
  UniformIntDistribution random;
  for (uint16_t max_symbol = 1; max_symbol < APROBA_MAX_SYMBOL; ++max_symbol) {
    // Get random probabilities.
    uint16_t counts[APROBA_MAX_SYMBOL] = {0u};
    uint32_t sum = 0u;
    for (uint32_t s = 0; s < max_symbol; ++s) {
      counts[s] = random.Get((s == 0) ? 1u : 0u, 100u);
      sum += counts[s];
    }
    // Deduce a cdf.
    uint16_t cdf[APROBA_MAX_SYMBOL] = {0};
    for (uint32_t s = 1; s < max_symbol; ++s) {
      cdf[s] = cdf[s - 1] + counts[s - 1];
    }
    ANSAdaptiveSymbol asym;
    EXPECT_WP2_OK(asym.InitFromCDF(cdf, max_symbol, sum));
    // Verify probabilities are the same.
    for (uint32_t s = 0; s < max_symbol; ++s) {
      EXPECT_NEAR((double)asym.GetProba(s), (double)counts[s] / sum, 0.001);
    }
  }
}

//------------------------------------------------------------------------------

class TestANSMapping : public ::testing::TestWithParam<uint32_t> {};

TEST_P(TestANSMapping, Simple) {
  ANSInit();
  UniformIntDistribution random(kBaseSeed);
  const uint16_t range = random.Get<uint16_t>(1, 100);
  std::vector<OptimizeArrayStorageStat> stats(range);
  std::vector<uint16_t> mapping_ini(range);
  Vector_u16 mapping_res;
  const uint32_t pattern_type = GetParam();
  for (uint32_t n = 1; n < range; ++n) {
    if (pattern_type == 0) {
      mapping_ini[0] = random.Get(0u, range - n);
      for (uint32_t j = 1; j < n; ++j) {
        mapping_ini[j] =
            random.Get<uint16_t>(mapping_ini[j - 1] + 1, range - (n - j));
      }
    } else {
      for (uint32_t j = 0; j < n / 2; ++j) mapping_ini[j] = j;
      for (uint32_t j = n / 2; j < n; ++j) mapping_ini[j] = j + range - n;
    }
    for (uint32_t j = 1; j < n; ++j) {
      ASSERT_LT(mapping_ini[j - 1], mapping_ini[j])
          << "pos: " << j << " size:" << n << " range:" << range;
      ASSERT_LT(mapping_ini[j], range) << "pos: " << j;
    }
    // Store the mapping to an encoder.
    ANSEnc enc;
    StoreMapping(mapping_ini.data(), n, range, stats.data(), &enc);
    enc.PutRange(31, 0, 263, "verif-bit");
    EXPECT_WP2_OK(enc.Assemble());
    // Read back the mapping.
    const size_t size = enc.BufferSize();
    const uint8_t* buf = enc.Buffer();
    ExternalDataSource data_source(buf, size);
    ANSDec dec(&data_source);
    EXPECT_WP2_OK(LoadMapping(&dec, n, range, mapping_res));
    for (uint32_t i = 0; i < n; ++i) {
      EXPECT_EQ(mapping_ini[i], mapping_res[i]);
    }
    EXPECT_EQ(dec.ReadRange(0, 263, "verif-bit"), 31);
  }
}

INSTANTIATE_TEST_SUITE_P(TestANSMappingInstantiation, TestANSMapping,
                         ::testing::Range(0u, 2u));

//------------------------------------------------------------------------------

TEST_F(TestANS, DictionarySymbolCost) {
  ANSDictionary dict;
  EXPECT_NE(dict.Init(0), WP2_STATUS_OK);
  EXPECT_WP2_OK(dict.Init(/*max_symbol=*/4));

  EXPECT_EQ(dict.SymbolCost(0), 0);
  EXPECT_EQ(dict.SymbolCost(1), 0);
  EXPECT_EQ(dict.SymbolCost(2), 0);
  EXPECT_EQ(dict.SymbolCost(3), 0);

  dict.RecordSymbol(1, 2);

  EXPECT_EQ(dict.SymbolCost(0), 0);
  EXPECT_EQ(dict.SymbolCost(1), 0);
  EXPECT_EQ(dict.SymbolCost(2), 0);
  EXPECT_EQ(dict.SymbolCost(3), 0);

  dict.RecordSymbol(0, 2);

  EXPECT_EQ(dict.SymbolCost(0), 1);
  EXPECT_EQ(dict.SymbolCost(1), 1);
  EXPECT_EQ(dict.SymbolCost(2), 0);
  EXPECT_EQ(dict.SymbolCost(3), 0);

  dict.RecordSymbol(3);

  EXPECT_NEAR(dict.SymbolCost(0), 1.32, 0.01);
  EXPECT_NEAR(dict.SymbolCost(1), 1.32, 0.01);
  EXPECT_EQ(dict.SymbolCost(2), 0);
  EXPECT_NEAR(dict.SymbolCost(3), 2.32, 0.01);

  // Test max_symbol.
  EXPECT_EQ(dict.SymbolCost(0, /*max_symbol=*/2), 1);
  EXPECT_EQ(dict.SymbolCost(1, /*max_symbol=*/2), 1);
  EXPECT_EQ(dict.SymbolCost(2, /*max_symbol=*/2), 0);

  EXPECT_EQ(dict.SymbolCost(0, /*max_symbol=*/1), 1);
  EXPECT_EQ(dict.SymbolCost(1, /*max_symbol=*/1), 1);

  EXPECT_EQ(dict.SymbolCost(0, /*max_symbol=*/0), 0);

  // Test quantized dictionary.
  dict.ResetCounts();
  dict.RecordSymbol(0, 1);
  dict.RecordSymbol(2, 20);
  for (uint32_t i = 0; i < 2; ++i) {
    Vector_u32 quantized_counts;
    ASSERT_TRUE(quantized_counts.resize(4));
    quantized_counts[0] = 2;
    quantized_counts[1] = 0;
    quantized_counts[2] = 2;
    quantized_counts[3] = 0;
    EXPECT_WP2_OK(dict.SetQuantizedCounts(quantized_counts));
    EXPECT_EQ(dict.SymbolCost(0), 1);
    EXPECT_EQ(dict.SymbolCost(1), 0);
    EXPECT_EQ(dict.SymbolCost(2), 1);
    EXPECT_EQ(dict.SymbolCost(3), 0);
    // Even if no symbol has been recorded, the cost needs to be the one from
    // the quantized counts.
    dict.ResetCounts();
  }
}

//------------------------------------------------------------------------------

TEST_F(TestANS, LargeRange) {
  // List of <value, min, max>
  // Some edge cases.
  std::vector<std::tuple<uint32_t, uint32_t, uint32_t>> test_cases = {
      {100, 0, kANSMaxRange - 2},
      {100, 0, kANSMaxRange - 1},
      {100, 0, kANSMaxRange},
      {kANSMaxRange - 2, 0, kANSMaxRange - 2},
      {kANSMaxRange - 1, 0, kANSMaxRange - 1},
      {kANSMaxRange, 0, kANSMaxRange},
  };
  // More random test cases.
  UniformIntDistribution random(kBaseSeed);
  for (uint32_t i = 0; i < 1; ++i) {
    const uint32_t max = random.Get(1u, kANSMaxRange * 10);
    const uint32_t min = random.Get(0u, max - 1);
    const uint32_t value = random.Get(min, max);
    test_cases.push_back(std::make_tuple(value, min, max));
  }

  // Write.
  ANSEnc enc;
  for (const auto& t : test_cases) {
    const uint32_t value = std::get<0>(t);
    const uint32_t min = std::get<1>(t);
    const uint32_t max = std::get<2>(t);
    SCOPED_TRACE(SPrintf("value: %d min: %d max: %d", value, min, max));
    EXPECT_EQ(value, PutLargeRange(value, min, max, &enc, "label"));

    // Check cost.
    ANSEncCounter counter;
    PutLargeRange(value, min, max, &counter, "label");
    const uint32_t range = max - min + 1;
    EXPECT_NEAR(counter.GetCost(), std::log2(range), 0.01);
  }

  // Read.
  ASSERT_WP2_OK(enc.Assemble());
  const size_t size = enc.BufferSize();
  const uint8_t* buf = enc.Buffer();
  ExternalDataSource data_source(buf, size);
  ANSDec dec(&data_source);
  for (const auto& t : test_cases) {
    const uint32_t value = std::get<0>(t);
    const uint32_t min = std::get<1>(t);
    const uint32_t max = std::get<2>(t);
    SCOPED_TRACE(SPrintf("value: %d min: %d max: %d", value, min, max));
    EXPECT_EQ(value, ReadLargeRange(min, max, &dec, "label"));
  }
}

}  // namespace
}  // namespace WP2
