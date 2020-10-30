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

// Test SymbolWriter/SymbolReader.

#include <algorithm>
#include <limits>
#include <memory>

#include "examples/example_utils.h"
#include "include/helpers.h"
#include "src/common/symbols.h"
#include "src/dec/symbols_dec.h"
#include "src/enc/symbols_enc.h"
#include "src/utils/random.h"

namespace WP2 {
namespace {

//------------------------------------------------------------------------------

using StorageMethod = SymbolsInfo::StorageMethod;

// Expose SymbolsInfo::Init() for test purposes.
class SymbolsInfoTest : public SymbolsInfo {
 public:
  using SymbolsInfo::SetInfo;
};

//------------------------------------------------------------------------------

void SetupInfo(uint32_t num_symbol_types, uint32_t max_range,
               bool can_be_negative, uint32_t num_clusters,
               UniformIntDistribution* random, SymbolsInfoTest* info) {
  for (uint32_t symbol = 0; symbol < num_symbol_types; ++symbol) {
    const int32_t max = random->Get(0u, max_range - 1);
    info->SetInfo(symbol, can_be_negative ? -max : 0, max, num_clusters,
                  StorageMethod::kAuto);
  }
}

class DataGenerator {
 public:
  struct Symbol {
    uint32_t s;
    int32_t value;
    uint32_t cluster;
    bool can_be_negative;
  };
  // Enum to hint at the distribution of a symbol.
  enum class kHint {
    kUniform,  // uniform distribution
    kDict,     // random distribution
    kGolomb,   // exponential distribution
    kTrivial   // unique value
  };
  DataGenerator(const SymbolsInfoTest& symbols_info,
                UniformIntDistribution* random)
      : symbols_info_(&symbols_info), random_(random) {}

  // Adds data consisting of 'num_symbols', for the given symbol 'sym',
  // in a random range chosen in [0, 'max_range'].
  void AddData(uint32_t sym, uint32_t num_symbols, kHint distribution,
               SymbolRecorder* const recorder, std::vector<Symbol>* symbols) {
    ANSEncNoop enc;
    const uint32_t range = symbols_info_->GetMaxRange(sym);
    const uint32_t num_symbol_types = symbols_info_->Size();
    counts_.resize(num_symbol_types);
    const uint32_t num_clusters = symbols_info_->NumClusters(sym);
    counts_[sym].resize(num_clusters);
    for (uint32_t c = 0; c < num_clusters; ++c) {
      counts_[sym][c].resize(range, 0);
    }

    // Create some biased statistics by repeating values a certain number of
    // times and then uniformly picking from the vector.
    std::vector<uint32_t> probas;
    int32_t trivial_value = std::numeric_limits<int32_t>::max();
    switch (distribution) {
      case kHint::kDict:
        probas.reserve(num_symbol_types);
        for (uint32_t r = 0; r < range; ++r) {
          const uint32_t tmp = random_->Get(0u, 1000u);
          for (uint32_t j = 0; j <= tmp; ++j) probas.push_back(r);
        }
        break;
      case kHint::kGolomb:
        probas.reserve(num_symbol_types);
        for (uint32_t r = 0; r < range; ++r) {
          const double rand_max = 1000. * std::exp(-(double)r);
          const uint32_t tmp =
              random_->Get(0u, std::max((uint32_t)rand_max, 1u));
          for (uint32_t j = 0; j <= tmp; ++j) probas.push_back(r);
        }
        break;
      case kHint::kTrivial:
        trivial_value = random_->Get(0u, range);
        break;
      default:
        break;
    }

    // Create some stream using the biased statistics.
    const bool can_be_negative = (symbols_info_->Min(/*cluster=*/0, sym) < 0);
    bool has_at_least_one_negative = false;
    for (uint32_t i = 0; i < num_symbols; ++i) {
      const uint32_t cluster = random_->Get(0u, num_clusters - 1);
      int32_t value = std::numeric_limits<int32_t>::max();;
      switch (distribution) {
        case kHint::kDict:
        case kHint::kGolomb:
          // Other symbols have biased data for which a dictionary will be more
          // efficient.
          value = probas[random_->Get(0, (int)probas.size() - 1)];
          break;
        case kHint::kUniform:
          // Getting uniform probabilities for the symbol will force it to be
          // stored as a range by the SymbolWriter.
          value = random_->Get(0u, range - 1);
          break;
        case kHint::kTrivial:
          value = trivial_value;
          break;
      }
      Symbol s;
      s.s = sym;
      s.value = value + symbols_info_->Min(/*cluster=*/0, sym);
      has_at_least_one_negative |= (s.value < 0);
      s.can_be_negative = can_be_negative;
      s.cluster = cluster;
      symbols->push_back(s);
      ++(counts_)[sym][cluster][value];
      recorder->Process(cluster, sym, s.value, "label", &enc);
    }
    // Make sure there is at least one negative value if asked.
    EXPECT_TRUE(!can_be_negative || has_at_least_one_negative);
  }

  void WriteHeader(const SymbolsInfoTest& info, const SymbolRecorder& recorder,
                   uint32_t* max_nnz, ANSDictionaries* const dicts,
                   ANSEnc* const enc, SymbolWriter* const sw) {
    ASSERT_WP2_OK(sw->Init(info));
    ASSERT_WP2_OK(sw->Allocate());
    // Get the maximum number of non-zero values by aggregating over all
    // clusters.
    *max_nnz = 0;
    for (uint32_t s = 0; s < symbols_info_->Size(); ++s) {
      if (counts_[s].empty() > 0) continue;
      for (uint32_t c = 0; c < symbols_info_->NumClusters(s); ++c) {
        *max_nnz = std::max(*max_nnz, std::accumulate(counts_[s][c].begin(),
                                                      counts_[s][c].end(), 0u));
      }
    }
    // Write the headers.
    for (uint32_t s = 0; s < symbols_info_->Size(); ++s) {
      for (uint32_t c = 0; c < symbols_info_->NumClusters(s); ++c) {
        ASSERT_WP2_OK(
            sw->WriteHeader(c, *max_nnz, s, recorder, "counts", enc, dicts));
      }
    }
  }

  void ReadHeader(const SymbolsInfoTest& info, uint32_t max_nnz,
                  ANSDec* const dec, SymbolReader* const sr) {
    ASSERT_WP2_OK(sr->Init(info, dec));
    ASSERT_WP2_OK(sr->Allocate());
    for (uint32_t s = 0; s < symbols_info_->Size(); ++s) {
      for (uint32_t c = 0; c < symbols_info_->NumClusters(s); ++c) {
        ASSERT_WP2_OK(sr->ReadHeader(c, max_nnz, s, "counts"));
      }
    }
  }

 private:
  const SymbolsInfoTest* const symbols_info_;
  UniformIntDistribution* random_;
  // Per symbol, per cluster, per value.
  std::vector<std::vector<std::vector<uint32_t>>> counts_;
};

constexpr std::array<DataGenerator::kHint, 4> kDistributions = {
    DataGenerator::kHint::kUniform, DataGenerator::kHint::kDict,
    DataGenerator::kHint::kGolomb, DataGenerator::kHint::kTrivial};

//------------------------------------------------------------------------------

TEST(SymbolsInfo, BasicTest) {
  SymbolsInfoTest info;
  info.SetInfo(/*sym=*/0, /*min=*/-9, /*max=*/9, /*num_clusters=*/2,
               StorageMethod::kAuto);
  EXPECT_WP2_OK(
      info.SetMinMax(/*cluster=*/1, /*sym=*/0, /*min=*/-4, /*max=*/4));
  info.SetInfo(/*sym=*/1, /*min=*/0, /*max=*/1, /*num_clusters=*/3,
               StorageMethod::kAdaptiveBit);
  EXPECT_WP2_OK(info.SetStartingProba(/*cluster=*/0, /*sym=*/1, 3, 12));
  EXPECT_WP2_OK(info.SetStartingProba(/*cluster=*/2, /*sym=*/1, 8, 2));
  info.SetInfo(/*sym=*/3, /*min=*/0, /*max=*/3,
               /*num_clusters=*/3, StorageMethod::kAuto);
  EXPECT_WP2_OK(
      info.SetMinMax(/*cluster=*/2, /*sym=*/3, /*min=*/0, /*max=*/19));

  // Number of symbol. Even though we didn't set any info for symbol '2', it is
  // assumed to exist.
  EXPECT_EQ(info.Size(), 4u);
  EXPECT_EQ(info.GetMaxRange(), 19u);
  EXPECT_EQ(info.MaxRangeSum(), 41u);
  EXPECT_EQ(info.RangeSum(), 62u);

  EXPECT_EQ(info.Range(/*cluster=*/0, 0), 19u);
  EXPECT_EQ(info.Range(/*cluster=*/1, 0), 9u);
  EXPECT_EQ(info.GetMaxRange(0), 19u);
  EXPECT_EQ(info.Range(/*cluster=*/0, 1), 2u);
  EXPECT_EQ(info.Range(/*cluster=*/1, 1), 2u);
  EXPECT_EQ(info.Range(/*cluster=*/2, 1), 2u);
  EXPECT_EQ(info.GetMaxRange(1), 2u);
  EXPECT_EQ(info.GetMaxRange(2), 0u);
  EXPECT_EQ(info.Range(/*cluster=*/0, 3), 4u);
  EXPECT_EQ(info.Range(/*cluster=*/1, 3), 4u);
  EXPECT_EQ(info.Range(/*cluster=*/2, 3), 20u);
  EXPECT_EQ(info.GetMaxRange(3), 20u);

  EXPECT_EQ(info.NumClusters(0), 2u);
  EXPECT_EQ(info.NumClusters(1), 3u);
  EXPECT_EQ(info.NumClusters(2), 0u);
  EXPECT_EQ(info.NumClusters(3), 3u);

  EXPECT_EQ(info.Method(0), StorageMethod::kAuto);
  EXPECT_EQ(info.Method(1), StorageMethod::kAdaptiveBit);
  EXPECT_EQ(info.Method(2), StorageMethod::kAuto);
  EXPECT_EQ(info.Method(3), StorageMethod::kAuto);

  SymbolsInfoTest copy;
  ASSERT_WP2_OK(copy.CopyFrom(info));
  EXPECT_TRUE(copy == info);
  copy.SetInfo(/*sym=*/2, /*min=*/0, /*max=*/4, /*num_clusters=*/1,
               StorageMethod::kAuto);
  EXPECT_FALSE(copy == info);
  EXPECT_EQ(copy.Size(), 4u);
  EXPECT_EQ(copy.GetMaxRange(), 19u);
  EXPECT_EQ(copy.MaxRangeSum(), 46u);
  EXPECT_EQ(copy.RangeSum(), 67u);

  EXPECT_EQ(info.StartingProbaP0(/*cluster=*/0, /*sym=*/1), 3u);
  EXPECT_EQ(info.StartingProbaP1(/*cluster=*/0, /*sym=*/1), 12u);
  EXPECT_EQ(info.StartingProbaP0(/*cluster=*/1, /*sym=*/1), 1u);
  EXPECT_EQ(info.StartingProbaP1(/*cluster=*/1, /*sym=*/1), 1u);
  EXPECT_EQ(info.StartingProbaP0(/*cluster=*/2, /*sym=*/1), 8u);
  EXPECT_EQ(info.StartingProbaP1(/*cluster=*/2, /*sym=*/1), 2u);
}

TEST(SymbolsInfo, LossessSymbolsInfo) {
  WP2L::LosslessSymbolsInfo info;
  info.Init(/*has_alpha=*/true, WP2_Argb_32, /*num_clusters=*/3);
  info.SetCacheRange(1 << 2);
  EXPECT_EQ(info.Range(/*cluster=*/0, WP2L::kSymbolA), (uint32_t)(1 << 8));
  EXPECT_EQ(info.NumClusters(WP2L::kSymbolA), 3u);

  EXPECT_EQ(info.Range(/*cluster=*/0, WP2L::kSymbolCache), (uint32_t)(1 << 2));
  EXPECT_EQ(info.NumClusters(WP2L::kSymbolCache), 3u);

  info.SetCacheRange(1 << 4);
  EXPECT_EQ(info.Range(/*cluster=*/0, WP2L::kSymbolCache), (uint32_t)(1 << 4));
  EXPECT_EQ(info.NumClusters(WP2L::kSymbolCache), 3u);

  info.SetNumClusters(2);
  EXPECT_EQ(info.NumClusters(WP2L::kSymbolA), 2u);
  EXPECT_EQ(info.NumClusters(WP2L::kSymbolCache), 2u);
}

// Helper function for the following two tests.
void VerifyProcessRead(const SymbolsInfoTest& info, uint32_t num_symbols,
                       uint32_t num_symbol_types,
                       DataGenerator* const generator) {
  std::vector<DataGenerator::Symbol> symbols;
  SymbolRecorder recorder;
  ASSERT_WP2_OK(recorder.Allocate(info, num_symbols));

  for (uint32_t s = 0; s < num_symbol_types; ++s) {
    generator->AddData(
        s, num_symbols,
        kDistributions[s * kDistributions.size() / num_symbol_types], &recorder,
        &symbols);
  }
  Shuffle(symbols.begin(), symbols.end(), 0);

  // Create the symbol writer.
  std::unique_ptr<WP2::SymbolWriter>
      sw(new (WP2Allocable::nothrow) WP2::SymbolWriter);
  ASSERT_TRUE(sw != nullptr);
  ANSEnc enc;
  ANSDictionaries dicts;
  uint32_t max_nnz = 0;
  generator->WriteHeader(info, recorder, &max_nnz, &dicts, &enc, sw.get());
  const double header_cost = enc.GetCost();

  // Write the symbols.
  // We do not send max_nnz, the ranges, the clusters and the symbol types as
  // they do not require the SymbolWriter, they usually are just sent as raw
  // signal in the ANS encoder.
  double enc_cost = 0;
  for (const DataGenerator::Symbol& sym : symbols) {
    sw->ProcessWithCost(sym.cluster, sym.s, sym.value, "sym", &enc, &enc_cost);
  }
  EXPECT_WP2_OK(enc.Assemble());
  const size_t size = enc.BufferSize();
  const uint8_t* buf = enc.Buffer();

  // Decode the symbol reader statistics.
  ExternalDataSource data_source(buf, size);
  ANSDec dec(&data_source);
  SymbolReader sr;
  generator->ReadHeader(info, max_nnz, &dec, &sr);

  // Verify the symbols do match the original data.
  double dec_cost = 0;
  for (const DataGenerator::Symbol& sym : symbols) {
    ASSERT_EQ(sr.Read(sym.cluster, sym.s, "sym", &dec_cost), sym.value);
  }
  const double symbol_cost = enc.GetCost() - header_cost;
  EXPECT_NEAR(symbol_cost, enc_cost, 0.2f);
  EXPECT_NEAR(symbol_cost, dec_cost, 0.2f);
  ASSERT_WP2_OK(dec.GetStatus());
}

// Test basic functionality of the SymbolWriter/SymbolReader.
TEST(SymbolsTest, Basic) {
  static constexpr uint32_t kNumSymbol = 1000u;
  static constexpr uint32_t kNumSymbolTypes = 5u;
  static constexpr uint32_t kMaxRange = 100u;
  static constexpr uint32_t kNumClusters = 10u;
  UniformIntDistribution random(/*seed=*/0);

  static_assert(kNumSymbolTypes <= kSymbolNumMax, "Need less symbols.");

  // Generate the data.
  SymbolsInfoTest info;
  SetupInfo(kNumSymbolTypes, kMaxRange, /*can_be_negative=*/true, kNumClusters,
            &random, &info);
  DataGenerator generator(info, &random);

  VerifyProcessRead(info, kNumSymbol, kNumSymbolTypes, &generator);
}

// Test the SymbolWriter/SymbolReader when CDF is given.
TEST(SymbolsTest, CDF) {
  ANSInit();
  static constexpr uint32_t kNumSymbol = 10000u;
  static constexpr uint32_t kNumSymbolTypes = 2u;
  static constexpr uint32_t kMaxRange = APROBA_MAX_SYMBOL;
  static constexpr uint32_t kNumClusters = 10u;
  static constexpr uint32_t kMaxCount = 1000u;
  UniformIntDistribution random(/*seed=*/0);

  static_assert(kNumSymbolTypes <= kSymbolNumMax, "Need less symbols.");

  // Generate the data.
  SymbolsInfoTest info;
  SetupInfo(kNumSymbolTypes, kMaxRange, /*can_be_negative=*/false, kNumClusters,
            &random, &info);
  uint16_t cdfs[kNumSymbolTypes][APROBA_MAX_SYMBOL];
  for (uint32_t i = 0; i < kNumSymbolTypes; ++i) {
    cdfs[i][0] = 0;
    for (uint32_t j = 1; j < APROBA_MAX_SYMBOL; ++j) {
      const uint32_t c = random.Get(1u, kMaxCount);
      cdfs[i][j] = cdfs[i][j - 1] + c;
    }
    // Normalize.
    const uint32_t sum =
        cdfs[i][APROBA_MAX_SYMBOL - 1] + random.Get(1u, kMaxCount);
    for (uint32_t j = 1; j < APROBA_MAX_SYMBOL; ++j) {
      cdfs[i][j] = cdfs[i][j] * ANS_MAX_SYMBOLS / sum;
    }
  }
  for (uint32_t i = 0; i < kNumSymbolTypes; ++i) {
    info.SetInfo(i, /*min=*/0, /*max=*/info.Range(/*cluster=*/0, i) - 1,
                 info.NumClusters(i), StorageMethod::kAdaptiveSym);
    EXPECT_WP2_OK(info.SetInitialCDF(cdfs[i], ANS_MAX_SYMBOLS,
                                     /*cluster=*/0, i));
  }
  DataGenerator generator(info, &random);

  VerifyProcessRead(info, kNumSymbol, kNumSymbolTypes, &generator);
}

//------------------------------------------------------------------------------

// Test the SymbolWriter/SymbolReader for uniform probability.
TEST(SymbolsTest, Uniform) {
  static constexpr uint32_t kNumSymbol = 10000u;
  static constexpr uint32_t kMaxRange = 100u;
  UniformIntDistribution random(/*seed=*/0);

  // Generate the data.
  SymbolsInfoTest info;
  SetupInfo(/*num_symbol_types=*/1, kMaxRange, /*can_be_negative=*/true,
            /*num_clusters=*/1, &random, &info);
  DataGenerator generator(info, &random);

  SymbolRecorder recorder;
  ASSERT_WP2_OK(recorder.Allocate(info, kNumSymbol));

  std::vector<DataGenerator::Symbol> symbols;
  generator.AddData(0, kNumSymbol, DataGenerator::kHint::kUniform, &recorder,
                    &symbols);

  // Create the symbol writer.
  std::unique_ptr<WP2::SymbolWriter>
      sw(new (WP2Allocable::nothrow) WP2::SymbolWriter);
  ASSERT_TRUE(sw != nullptr);
  ANSEnc enc;
  ANSDictionaries dicts;
  uint32_t max_nnz = 0;
  generator.WriteHeader(info, recorder, &max_nnz, &dicts, &enc, sw.get());

  // Write the symbols.
  // We do not send max_nnz, the ranges, the clusters and the symbol types as
  // they do not require the SymbolWriter, they usually are just sent as raw
  // signal in the ANS encoder.
  int32_t min = 0, max = 0;
  for (const DataGenerator::Symbol& sym : symbols) {
    sw->Process(sym.cluster, sym.s, sym.value, "sym", &enc);
    min = std::min(min, sym.value);
    max = std::max(max, sym.value);
  }
  const uint32_t range = max - min + 1;
  EXPECT_GT(range, 0u);
  EXPECT_WP2_OK(enc.Assemble());
  const size_t size = enc.BufferSize();

  // With uniform data, the symbol writer should choose the range storage
  // method (there is a little overhead for storing the range).
  ASSERT_NEAR(symbols.size() * std::log2(range) / 8, size, size * 0.0015);
}

//------------------------------------------------------------------------------

// Test the fact that the SymbolWriter/SymbolReader can work when restricting
// the range of values.
TEST(SymbolsTest, MaxValue) {
  const uint32_t kNumSymbol = 1000u;
  const uint32_t kNumSymbolTypes = 5u;
  const uint32_t kMaxRange = 100u;
  const uint32_t kNumClusters = 1u;
  ASSERT_EQ(kNumClusters, 1u);
  UniformIntDistribution random(/*seed=*/0);

  static_assert(kNumSymbolTypes < kSymbolNumMax, "Need less symbols.");

  // Generate the data.
  SymbolsInfoTest info;
  SetupInfo(kNumSymbolTypes, kMaxRange, /*can_be_negative=*/true, kNumClusters,
            &random, &info);
  DataGenerator generator(info, &random);
  std::vector<DataGenerator::Symbol> symbols;

  SymbolRecorder recorder;
  ASSERT_WP2_OK(recorder.Allocate(info, kNumSymbol));

  for (uint32_t s = 0; s < kNumSymbolTypes; ++s) {
    generator.AddData(
        s, kNumSymbol,
        kDistributions[s * kDistributions.size() / kNumSymbolTypes], &recorder,
        &symbols);
  }
  Shuffle(symbols.begin(), symbols.end(), 0);

  // Create the symbol writer.
  std::unique_ptr<WP2::SymbolWriter>
      sw(new (WP2Allocable::nothrow) WP2::SymbolWriter);
  ASSERT_TRUE(sw != nullptr);
  ANSEnc enc;
  ANSDictionaries dicts;
  uint32_t max_nnz = 0;
  generator.WriteHeader(info, recorder, &max_nnz, &dicts, &enc, sw.get());

  // Write the data. We artificially cap to multiples of the value.
  for (float mul : {1.f, 1.2f, 2.f, 100.f}) {
    for (const DataGenerator::Symbol& sym : symbols) {
      sw->Process(sym.cluster, sym.s, sym.value, /*max_value=*/
                  std::abs(sym.value) * mul, "sym", &enc);
    }
  }
  EXPECT_WP2_OK(enc.Assemble());
  const size_t size = enc.BufferSize();
  const uint8_t* buf = enc.Buffer();

  // Decode the symbol reader statistics.
  ExternalDataSource data_source(buf, size);
  ANSDec dec(&data_source);
  SymbolReader sr;
  generator.ReadHeader(info, max_nnz, &dec, &sr);

  // Verify the symbols do match the original data.
  for (float mul : {1.f, 1.2f, 2.f, 100.f}) {
    for (const DataGenerator::Symbol& sym : symbols) {
      int32_t value;
      ASSERT_WP2_OK(
          sr.TryRead(sym.cluster, sym.s,
                     /*max_value=*/std::abs(sym.value) * mul, "sym", &value));
      ASSERT_EQ(std::abs(value), std::abs(sym.value)) << "mul " << mul;
    }
  }
  ASSERT_WP2_OK(dec.GetStatus());
}

// Test the SymbolWriter/SymbolReader for StorageMethod::kAdaptive.
TEST(SymbolsTest, StorageMethodAdaptive) {
  ANSInit();
  static constexpr uint32_t kNumSymbol = 1000u;
  static constexpr uint32_t kMaxRange = 10u;
  UniformIntDistribution random(/*seed=*/0);

  // Generate the data.
  SymbolsInfoTest info;
  info.SetInfo(/*sym=*/0, /*min=*/0, /*max=*/kMaxRange - 1, /*num_clusters=*/1,
               StorageMethod::kAdaptiveSym);
  DataGenerator generator(info, &random);

  SymbolRecorder recorder;
  ASSERT_WP2_OK(recorder.Allocate(info, kNumSymbol));

  std::vector<DataGenerator::Symbol> symbols;
  generator.AddData(0, kNumSymbol, DataGenerator::kHint::kUniform, &recorder,
                    &symbols);

  // Create the symbol writer.
  std::unique_ptr<WP2::SymbolWriter>
      sw(new (WP2Allocable::nothrow) WP2::SymbolWriter);
  ASSERT_TRUE(sw != nullptr);
  ANSEnc enc;
  ANSDictionaries dicts;
  uint32_t max_nnz = 0;
  generator.WriteHeader(info, recorder, &max_nnz, &dicts, &enc, sw.get());

  EXPECT_EQ(enc.GetCost(), 0);  // Header should be free.

  // Write the symbols.
  for (const DataGenerator::Symbol& sym : symbols) {
    sw->Process(sym.cluster, sym.s, sym.value, "sym", &enc);
  }
  EXPECT_WP2_OK(enc.Assemble());
  size_t size = enc.BufferSize();

  // Since the data is uniform, the adaptive bit shouldn't win much.
  EXPECT_NEAR(symbols.size() * std::log2(kMaxRange) / 8, size, size * 0.1);

  // Decode the symbol reader statistics.
  const uint8_t* buf = enc.Buffer();
  ExternalDataSource data_source(buf, size);
  ANSDec dec(&data_source);
  SymbolReader sr;
  generator.ReadHeader(info, max_nnz, &dec, &sr);

  // Verify the symbols do match the original data.
  for (const DataGenerator::Symbol& sym : symbols) {
    ASSERT_EQ(sr.Read(sym.cluster, sym.s, "sym"), sym.value);
  }
  ASSERT_WP2_OK(dec.GetStatus());

  // Try again with very predictable data.
  enc.Clear();
  for (uint32_t i = 0; i < kNumSymbol; ++i) {
    // Long string of zeros followed by a long string of ones.
    const uint32_t value = (i < kNumSymbol / 2) ? 0 : 1;
    sw->Process(/*cluster=*/0, /*sym=*/0, value, "sym", &enc);
  }
  EXPECT_WP2_OK(enc.Assemble());
  size = enc.BufferSize();
  // Since the data is very correlated it should be pretty cheap.
  EXPECT_LT(size, kNumSymbol * std::log2(kMaxRange) * 0.015);
}

//------------------------------------------------------------------------------

// Test the SymbolWriter/SymbolReader for StorageMethod::kAdaptiveWithAutoSpeed.
TEST(SymbolsTest, StorageMethodAdaptiveWithAutoSpeed) {
  ANSInit();
  static constexpr uint32_t kNumSymbol = 100u;
  static constexpr uint32_t kMaxRange = 10u;

  // Generate the data.
  SymbolsInfoTest info;
  info.SetInfo(/*sym=*/0, /*min=*/0, /*max=*/kMaxRange - 1, /*num_clusters=*/1,
               StorageMethod::kAdaptiveWithAutoSpeed);

  SymbolRecorder recorder;
  ASSERT_WP2_OK(recorder.Allocate(info, kNumSymbol));

  // Record symbols.
  ANSEncNoop noop;
  std::vector<DataGenerator::Symbol> symbols;
  for (uint32_t i = 0; i < kNumSymbol; ++i) {
    // Long string of zeros followed by a long string of ones.
    const int32_t value = (i < kNumSymbol / 2) ? 0 : 1;
    symbols.push_back(
        {/*s=*/0, value, /*cluster=*/0, /*can_be_negative=*/true});
    recorder.Process(/*cluster=*/0, /*sym=*/0, value, "sym", &noop);
  }

  // Create the symbol writer.
  std::unique_ptr<WP2::SymbolWriter>
      sw(new (WP2Allocable::nothrow) WP2::SymbolWriter);
  ASSERT_TRUE(sw != nullptr);
  ASSERT_WP2_OK(sw->Init(info));
  ASSERT_WP2_OK(sw->Allocate());
  ANSEnc enc;
  ANSDictionaries dicts;
  uint32_t max_nnz = kNumSymbol;
  ASSERT_WP2_OK(
      sw->WriteHeader(max_nnz, /*sym=*/0, recorder, "header", &enc, &dicts));

  for (const DataGenerator::Symbol& sym : symbols) {
    sw->Process(sym.cluster, sym.s, sym.value, "sym", &enc);
  }
  EXPECT_WP2_OK(enc.Assemble());
  size_t size = enc.BufferSize();

  // Since the data is very correlated it should be pretty cheap.
  EXPECT_LT(size, kNumSymbol * std::log2(kMaxRange) * 0.02);

  // Decode the symbol reader statistics.
  const uint8_t* buf = enc.Buffer();
  ExternalDataSource data_source(buf, size);
  ANSDec dec(&data_source);
  SymbolReader sr;
  ASSERT_WP2_OK(sr.Init(info, &dec));
  ASSERT_WP2_OK(sr.Allocate());
  ASSERT_WP2_OK(sr.ReadHeader(max_nnz, /*sym=*/0, "header"));

  // Verify the symbols do match the original data.
  for (const DataGenerator::Symbol& sym : symbols) {
    ASSERT_EQ(sr.Read(sym.cluster, sym.s, "sym"), sym.value);
  }
  ASSERT_WP2_OK(dec.GetStatus());
}

//------------------------------------------------------------------------------

// Test that SymbolCounter gives accurate measurements.
TEST(SymbolsTest, SymbolCounter) {
  ANSInit();
  static constexpr uint32_t kNumSymbol = 10000u;
  static constexpr uint32_t kMaxRange = 100u;
  static constexpr uint32_t kNumClusters = 10u;
  UniformIntDistribution random(/*seed=*/0);

  // Generate the data.
  SymbolsInfoTest info;
  SetupInfo(/*num_symbol_types=*/3, kMaxRange, /*can_be_negative=*/false,
            kNumClusters, &random, &info);
  // Adaptive bit.
  info.SetInfo(3, /*min=*/0, /*max=*/1, /*num_clusters=*/4,
               StorageMethod::kAdaptiveBit);
  EXPECT_WP2_OK(
      info.SetStartingProba(/*cluster=*/2, /*sym=*/3, /*p0=*/3, /*p1=*/7));
  EXPECT_WP2_OK(
      info.SetStartingProba(/*cluster=*/3, /*sym=*/3, /*p0=*/100, /*p1=*/1));
  // Adaptive symbol.
  info.SetInfo(info.Size(), /*min=*/0, /*max=*/9, /*num_clusters=*/3,
               StorageMethod::kAdaptiveSym);
  // Test that binary symbols are not only kAdaptiveBit but can also be
  // kAdaptiveSym.
  info.SetInfo(info.Size(), /*min=*/0, /*max=*/1, /*num_clusters=*/2,
               StorageMethod::kAdaptiveSym);
  const uint32_t num_symbol_types = info.Size();
  DataGenerator generator(info, &random);

  std::vector<DataGenerator::Symbol> symbols;

  SymbolRecorder recorder;
  ASSERT_WP2_OK(recorder.Allocate(info, kNumSymbol));

  for (uint32_t s = 0; s < num_symbol_types; ++s) {
    const DataGenerator::kHint distribution =
        (info.Method(s) == StorageMethod::kAdaptiveSym)
            ? DataGenerator::kHint::kDict
            : DataGenerator::kHint::kUniform;
    generator.AddData(s, kNumSymbol, distribution, &recorder, &symbols);
  }
  Shuffle(symbols.begin(), symbols.end(), 0);

  // Create the symbol writer and write the header.
  std::unique_ptr<WP2::SymbolWriter>
      sw(new (WP2Allocable::nothrow) WP2::SymbolWriter);
  ASSERT_TRUE(sw != nullptr);
  ANSEnc enc;
  ANSDictionaries dicts;
  uint32_t max_nnz = 0;
  generator.WriteHeader(info, recorder, &max_nnz, &dicts, &enc, sw.get());
  const float header_cost = enc.GetCost();

  // Write the symbols.
  // We do not send max_nnz, the ranges, the clusters and the symbol types as
  // they do not require the SymbolWriter, they usually are just sent as raw
  // signal in the ANS encoder.
  for (const DataGenerator::Symbol& sym : symbols) {
    sw->Process(sym.cluster, sym.s, sym.value, "sym", &enc);
  }
  const float real_cost = enc.GetCost() - header_cost;

  recorder.ResetRecord(/*reset_backup=*/true);
  ANSEncCounter enc_counter;
  SymbolCounter symbol_counter(&recorder);
  ASSERT_WP2_OK(symbol_counter.Allocate(/*syms=*/{3, 4, 5}));
  for (const DataGenerator::Symbol& sym : symbols) {
    symbol_counter.Process(sym.cluster, sym.s, sym.value, "sym", &enc_counter);
  }

  EXPECT_NEAR(enc_counter.GetCost(), real_cost, 9.1);
}

TEST(SymbolsTest, SymbolCounter_Dictionary) {
  WP2MathInit();
  // Generate the data.
  SymbolsInfoTest info;
  info.SetInfo(0, /*min=*/0, /*max=*/7, /*num_clusters=*/1,
               StorageMethod::kAuto);

  SymbolRecorder recorder;
  ASSERT_WP2_OK(recorder.Allocate(info, /*num_records=*/0));

  ANSEncCounter counter;
  SymbolCounter symbol_counter(&recorder);
  ASSERT_WP2_OK(symbol_counter.Allocate(/*syms=*/{}));
  symbol_counter.Process(0, 1, "test", &counter);
  // Before anything was recorded, all values are assumed to have the same cost
  // which is log2(range) = log2(8) = 3
  EXPECT_EQ(counter.GetCost(), 3);
  counter.Reset();
  symbol_counter.Process(0, 3, "test", &counter);
  EXPECT_EQ(counter.GetCost(), 3);

  ANSEncNoop noop;
  std::vector<uint8_t> data(400);
  for (uint32_t i = 0; i < 400; ++i) data[i] = (i < 100) ? 1 : 3;
  std::mt19937 gen(0);
  std::shuffle(data.begin(), data.end(), gen);
  for (uint8_t i : data) recorder.Process(/*sym=*/0, i, "test", &noop);

  counter.Reset();
  symbol_counter.Process(0, 1, "test", &counter);
  // After recording, costs are based on the recorded stats. -log2(1/4) = 2
  EXPECT_EQ(counter.GetCost(), 2);
  counter.Reset();
  symbol_counter.Process(0, 3, "test", &counter);
  // -log2(3/4)
  EXPECT_NEAR(counter.GetCost(), 0.415, 0.001);
  counter.Reset();
  symbol_counter.Process(0, 2, "test", &counter);
  // Value "2" was never recorded so we assume a cost above -1*log2(1/400).
  EXPECT_NEAR(counter.GetCost(), 8.643, 0.001);
  ASSERT_WP2_OK(recorder.MakeBackup());
  float costs[2];
  for (uint32_t reset_backup : {false, true}) {
    recorder.ResetRecord(reset_backup);
    counter.Reset();
    for (uint8_t i : data) {
      symbol_counter.Process(/*sym=*/0, i, "test", &counter);
    }
    costs[reset_backup ? 1 : 0] = counter.GetCost();
  }
  // If we use the recorded probabilities, we are much better at predicting the
  // cost.
  EXPECT_LT(costs[0], 0.5 * costs[1]);
}

// Test that SymbolCounter gives accurate measurements when given a
// SymbolRecorder that has already recorded a bunch of data.
TEST(SymbolsTest, SymbolCounter_Midway) {
  WP2MathInit();
  ANSInit();
  static constexpr uint32_t kNumSymbol = 10000u;
  UniformIntDistribution random(/*seed=*/0);

  // Here we only test adaptive bits and symbols, for which the counter should
  // know exactly how much space they'll take.
  SymbolsInfoTest info;
  // Adaptive bit.
  info.SetInfo(0, /*min=*/0, /*max=*/1, /*num_clusters=*/4,
               StorageMethod::kAdaptiveBit);
  EXPECT_WP2_OK(
      info.SetStartingProba(/*cluster=*/2, /*sym=*/0, /*p0=*/3, /*p1=*/7));
  EXPECT_WP2_OK(
      info.SetStartingProba(/*cluster=*/3, /*sym=*/0, /*p0=*/100, /*p1=*/1));
  // Adaptive symbol.
  info.SetInfo(1, /*min=*/0, /*max=*/9, /*num_clusters=*/3,
               StorageMethod::kAdaptiveSym);
  const uint32_t num_symbol_types = info.Size();
  DataGenerator generator(info, &random);

  std::vector<DataGenerator::Symbol> symbols;

  SymbolRecorder recorder;
  ASSERT_WP2_OK(recorder.Allocate(info, kNumSymbol));

  // Generate data.
  for (uint32_t s = 0; s < num_symbol_types; ++s) {
    generator.AddData(
        s, kNumSymbol,
        kDistributions[s * kDistributions.size() / num_symbol_types], &recorder,
        &symbols);
  }
  Shuffle(symbols.begin(), symbols.end(), 0);

  // Create the symbol writer and write the header.
  std::unique_ptr<WP2::SymbolWriter>
      sw(new (WP2Allocable::nothrow) WP2::SymbolWriter);
  ASSERT_TRUE(sw != nullptr);
  ANSEnc enc;
  ANSDictionaries dicts;
  uint32_t max_nnz = 0;
  generator.WriteHeader(info, recorder, &max_nnz, &dicts, &enc, sw.get());

  recorder.ResetRecord(/*reset_backup=*/true);

  // Write and record half of the symbols.
  // We do not send max_nnz, the ranges, the clusters and the symbol types as
  // they do not require the SymbolWriter, they usually are just sent as raw
  // signal in the ANS encoder.
  const uint32_t half = symbols.size() / 2;
  for (uint32_t i = 0; i < half; ++i) {
    const DataGenerator::Symbol& sym = symbols[i];
    sw->Process(sym.cluster, sym.s, sym.value, "sym", &enc);
    recorder.Process(sym.cluster, sym.s, sym.value, "sym", &enc);
  }

  const float half_cost = enc.GetCost();

  // Write the second half of symbols.
  for (uint32_t i = half; i < symbols.size(); ++i) {
    const DataGenerator::Symbol& sym = symbols[i];
    sw->Process(sym.cluster, sym.s, sym.value, "sym", &enc);
  }
  const float real_cost = enc.GetCost() - half_cost;

  ANSEncCounter enc_counter;
  SymbolCounter symbol_counter(&recorder);
  ASSERT_WP2_OK(symbol_counter.Allocate(/*syms=*/{0, 1}));
  for (uint32_t i = half; i < symbols.size(); ++i) {
    const DataGenerator::Symbol& sym = symbols[i];
    symbol_counter.Process(sym.cluster, sym.s, sym.value, "sym", &enc_counter);
  }

  EXPECT_NEAR(enc_counter.GetCost(), real_cost, 1);
}

struct GolombSymbol {
  uint32_t value;
  uint32_t range;
  uint32_t prefix_size;  // 1 or 2
};

TEST(SymbolsTest, SimpleGolomb) {
  WP2MathInit();
  static constexpr uint32_t kNumSymbol = 10000u;

  UniformIntDistribution random(/*seed=*/0);
  std::vector<GolombSymbol> symbols;
  symbols.reserve(kNumSymbol);
  for (uint32_t i = 0; i < kNumSymbol; ++i) {
    const uint32_t range = random.Get(1u, kANSMaxRange);
    const uint32_t value = random.Get(0u, range - 1);
    const uint32_t prefix_size = random.Get(1u, 2u);
    symbols.push_back({value, range, prefix_size});
  }

  ANSEnc enc;
  for (GolombSymbol& s : symbols) {
    WriteGolomb(s.value, /*min=*/0, /*max=*/s.range - 1, s.prefix_size, &enc,
                "test");
  }

  ASSERT_WP2_OK(enc.Assemble());
  ExternalDataSource data_source(enc.Buffer(), enc.BufferSize());
  ANSDec dec(&data_source);

  for (GolombSymbol& s : symbols) {
    ASSERT_EQ(s.value, ReadGolomb(/*min=*/0, /*max=*/s.range - 1, s.prefix_size,
                                  &dec, "test"));
  }
}

class SymbolWriterForTest : public SymbolWriter {
 public:
  using SymbolWriter::AddAdaptiveBit;
  using SymbolWriter::AddAdaptiveSymbol;
  using SymbolWriter::AddDict;
  using SymbolWriter::AddGolomb;
  using SymbolWriter::AddRange;
};

class SymbolReaderForTest : public SymbolReader {
 public:
  using SymbolReader::AddAdaptiveBit;
  using SymbolReader::AddAdaptiveSymbol;
  using SymbolReader::AddDict;
  using SymbolReader::AddGolomb;
  using SymbolReader::AddRange;
};

// bool: use_max
class SymbolCostTest : public ::testing::TestWithParam<bool> {
 protected:
  double Decode(SymbolReader* const sr, const std::vector<int32_t>& values,
                const std::vector<uint32_t>& max_values) {
    double dec_cost = 0;
    for (uint32_t i = 0; i < kNumValues; ++i) {
      SCOPED_TRACE(SPrintf("value  %d", i));
      if (!max_values.empty()) {
        int32_t v = 0;
        EXPECT_WP2_OK(sr->TryRead(kCluster, kSymbol, max_values[i], "label", &v,
                              &dec_cost));
        EXPECT_EQ(values[i], v);
      } else {
        EXPECT_EQ(values[i], sr->Read(kSymbol, "label", &dec_cost));
      }
    }
    return dec_cost;
  }

  UniformIntDistribution random_;

  static constexpr uint32_t kCluster = 0;
  static constexpr uint32_t kSymbol = 0;
  static constexpr uint32_t kNumValues = 1000;
  static constexpr uint32_t kRange = 10;
  // Allow 0.05% error per value.
  static constexpr double kCostError =
      (kNumValues * 0.05 / 100) < 0.01 ? 0.01 : (kNumValues * 0.05 / 100);
};

TEST_P(SymbolCostTest, TrivialCost) {
  const bool use_max = GetParam();

  const uint32_t kValue = 8;

  SymbolsInfo info;
  info.SetInfo(kSymbol, /*min=*/0, /*max=*/kRange - 1, /*num_clusters=*/1,
               StorageMethod::kAuto);

  SymbolWriterForTest sw;
  ASSERT_WP2_OK(sw.Init(info));
  ASSERT_WP2_OK(sw.Allocate());
  sw.AddTrivial(kCluster, kSymbol, kValue);

  std::vector<int32_t> values;
  std::vector<uint32_t> max_values;
  double enc_cost = 0;
  ANSEnc enc;
  for (uint32_t i = 0; i < kNumValues; ++i) {
    values.push_back(kValue);
    if (use_max) {
      const uint32_t max = random_.Get(kValue, kRange - 1);
      sw.ProcessWithCost(kCluster, kSymbol, kValue, max, "label", &enc,
                         &enc_cost);
      max_values.push_back(max);
    } else {
      sw.ProcessWithCost(kCluster, kSymbol, kValue, "label", &enc, &enc_cost);
    }
  }

  ASSERT_WP2_OK(enc.Assemble());
  ExternalDataSource data_source(enc.Buffer(), enc.BufferSize());
  ANSDec dec(&data_source);

  SymbolReaderForTest sr;
  ASSERT_WP2_OK(sr.Init(info, &dec));
  ASSERT_WP2_OK(sr.Allocate());

  sr.AddTrivial(kCluster, kSymbol, kValue);

  const double dec_cost = Decode(&sr, values, max_values);

  EXPECT_EQ(0., enc.GetCost());
  EXPECT_EQ(0., enc_cost);
  EXPECT_EQ(0., dec_cost);
}

TEST_P(SymbolCostTest, RangeCost) {
  const bool use_max = GetParam();

  SymbolsInfo info;
  info.SetInfo(kSymbol, /*min=*/0, /*max=*/kRange - 1, /*num_clusters=*/1,
               StorageMethod::kAuto);

  SymbolWriterForTest sw;
  ASSERT_WP2_OK(sw.Init(info));
  ASSERT_WP2_OK(sw.Allocate());
  sw.AddRange(kCluster, kSymbol, /*mapping=*/nullptr, /*size=*/0, kRange);

  std::vector<int32_t> values;
  std::vector<uint32_t> max_values;
  double enc_cost = 0;
  ANSEnc enc;
  for (uint32_t i = 0; i < kNumValues; ++i) {
    const uint32_t max =
        use_max ? random_.Get<uint32_t>(0, kRange - 1) : kRange - 1;
    const uint32_t v = random_.Get(0u, max);
    values.push_back(v);
    if (use_max) {
      max_values.push_back(max);
      sw.ProcessWithCost(kCluster, kSymbol, v, max, "label", &enc, &enc_cost);
    } else {
      sw.ProcessWithCost(kCluster, kSymbol, v, "label", &enc, &enc_cost);
    }
  }

  ASSERT_WP2_OK(enc.Assemble());
  ExternalDataSource data_source(enc.Buffer(), enc.BufferSize());
  ANSDec dec(&data_source);

  SymbolReaderForTest sr;
  ASSERT_WP2_OK(sr.Init(info, &dec));
  ASSERT_WP2_OK(sr.Allocate());

  sr.AddRange(kCluster, kSymbol, /*infos=*/nullptr, kRange);

  const double dec_cost = Decode(&sr, values, max_values);

  EXPECT_NEAR(enc.GetCost(), enc_cost, kCostError);
  EXPECT_NEAR(enc.GetCost(), dec_cost, kCostError);
}

TEST_P(SymbolCostTest, DictCost) {
  const bool use_max = GetParam();

  SymbolsInfo info;
  info.SetInfo(kSymbol, /*min=*/0, /*max=*/kRange - 1, /*num_clusters=*/1,
               StorageMethod::kAuto);

  std::vector<uint32_t> histogram = {10, 2, 50, 20, 5};
  std::vector<uint16_t> mapping = {0, 2, 3, 6, 7};
  const uint32_t symbol_size = histogram.size();
  VectorNoCtor<ANSSymbolInfo> infos;
  ASSERT_TRUE(infos.resize(symbol_size));
  for (uint32_t i = 0; i < symbol_size; ++i) {
    infos[i].freq = histogram[i];
    infos[i].symbol = mapping[i];
  }

  SymbolWriterForTest sw;
  ASSERT_WP2_OK(sw.Init(info));
  ASSERT_WP2_OK(sw.Allocate());

  ANSDictionaries dicts;
  ASSERT_WP2_OK(sw.AddDict(kCluster, kSymbol, histogram.data(),
                           histogram.data(), mapping.data(), symbol_size,
                           &dicts));

  std::vector<int32_t> values;
  std::vector<uint32_t> max_values;
  double enc_cost = 0;
  ANSEnc enc;
  for (uint32_t i = 0; i < kNumValues; ++i) {
    const uint32_t max_i = use_max
                               ? random_.Get<uint32_t>(0, mapping.size() - 1)
                               : mapping.size() - 1;
    const uint32_t max = mapping[max_i];
    const uint32_t v = mapping[random_.Get<uint32_t>(0, max_i)];
    values.push_back(v);
    if (use_max) {
      max_values.push_back(max);
      sw.ProcessWithCost(kCluster, kSymbol, v, max, "label", &enc, &enc_cost);
    } else {
      sw.ProcessWithCost(kCluster, kSymbol, v, "label", &enc, &enc_cost);
    }
  }

  ASSERT_WP2_OK(enc.Assemble());
  ExternalDataSource data_source(enc.Buffer(), enc.BufferSize());
  ANSDec dec(&data_source);

  SymbolReaderForTest sr;
  ASSERT_WP2_OK(sr.Init(info, &dec));
  ASSERT_WP2_OK(sr.Allocate());

  ASSERT_WP2_OK(sr.AddDict(kCluster, kSymbol, &infos));

  const double dec_cost = Decode(&sr, values, max_values);

  EXPECT_NEAR(enc.GetCost(), enc_cost, kCostError);
  EXPECT_NEAR(enc.GetCost(), dec_cost, kCostError);
}

TEST_P(SymbolCostTest, GolombCost) {
  const bool use_max = GetParam();

  SymbolsInfo info;
  info.SetInfo(kSymbol, /*min=*/0, /*max=*/kRange - 1, /*num_clusters=*/1,
               StorageMethod::kAuto);

  std::vector<uint16_t> mapping = {1, 2, 4, 6, 9};

  const uint32_t kGolombPrefixSize = 1;
  std::vector<uint32_t> golomb_histogram = {3, 7, 1, 2, 1};
  std::vector<uint16_t> golomb_mapping = {0, 1, 2, 3, 4};
  const uint32_t golomb_size = golomb_histogram.size();
  VectorNoCtor<ANSSymbolInfo> golomb_infos;
  ASSERT_TRUE(golomb_infos.resize(golomb_size));
  for (uint32_t i = 0; i < golomb_size; ++i) {
    golomb_infos[i].freq = golomb_histogram[i];
    golomb_infos[i].symbol = golomb_mapping[i];
  }

  SymbolWriterForTest sw;
  ASSERT_WP2_OK(sw.Init(info));
  ASSERT_WP2_OK(sw.Allocate());
  ANSDictionaries dicts;
  ASSERT_WP2_OK(sw.AddGolomb(kCluster, kSymbol, golomb_histogram.data(),
                             golomb_histogram.data(), golomb_mapping.data(),
                             golomb_size, kGolombPrefixSize, &dicts));

  std::vector<int32_t> values;
  std::vector<uint32_t> max_values;
  double enc_cost = 0;
  ANSEnc enc;
  for (uint32_t i = 0; i < kNumValues; ++i) {
    const uint32_t max_i = use_max
                               ? random_.Get<uint32_t>(0, mapping.size() - 1)
                               : mapping.size() - 1;
    const uint32_t max = mapping[max_i];
    const uint32_t v = mapping[random_.Get<uint32_t>(0, max_i)];
    values.push_back(v);
    if (use_max) {
      max_values.push_back(max);
      sw.ProcessWithCost(kCluster, kSymbol, v, max, "label", &enc, &enc_cost);
    } else {
      sw.ProcessWithCost(kCluster, kSymbol, v, "label", &enc, &enc_cost);
    }
  }

  ASSERT_WP2_OK(enc.Assemble());
  ExternalDataSource data_source(enc.Buffer(), enc.BufferSize());
  ANSDec dec(&data_source);

  SymbolReaderForTest sr;
  ASSERT_WP2_OK(sr.Init(info, &dec));
  ASSERT_WP2_OK(sr.Allocate());

  ASSERT_WP2_OK(
      sr.AddGolomb(kCluster, kSymbol, &golomb_infos, kGolombPrefixSize));

  const double dec_cost = Decode(&sr, values, max_values);

  EXPECT_NEAR(enc.GetCost(), enc_cost, kCostError);
  EXPECT_NEAR(enc.GetCost(), dec_cost, kCostError);
}

TEST_P(SymbolCostTest, AdaptiveBitCost) {
  WP2MathInit();
  const bool use_max = GetParam();

  SymbolsInfo info;
  info.SetInfo(kSymbol, /*min=*/0, /*max=*/1, /*num_clusters=*/1,
               StorageMethod::kAdaptiveBit);
  ANSDictionaries dicts;

  const uint32_t p0 = 3, p1 = 10;

  SymbolWriterForTest sw;
  ASSERT_WP2_OK(sw.Init(info));
  ASSERT_WP2_OK(sw.Allocate());
  ASSERT_WP2_OK(sw.AddAdaptiveBit(kCluster, kSymbol, p0, p1));

  std::vector<int32_t> values;
  std::vector<uint32_t> max_values;
  double enc_cost = 0;
  ANSEnc enc;
  for (uint32_t i = 0; i < kNumValues; ++i) {
    const uint32_t max = use_max ? random_.FlipACoin() : 1;
    const uint32_t v = random_.Get<uint32_t>(0, max);
    values.push_back(v);
    if (use_max) {
      max_values.push_back(max);
      sw.ProcessWithCost(kCluster, kSymbol, v, max, "label", &enc, &enc_cost);
    } else {
      sw.ProcessWithCost(kCluster, kSymbol, v, "label", &enc, &enc_cost);
    }
  }

  ASSERT_WP2_OK(enc.Assemble());
  ExternalDataSource data_source(enc.Buffer(), enc.BufferSize());
  ANSDec dec(&data_source);

  SymbolReaderForTest sr;
  ASSERT_WP2_OK(sr.Init(info, &dec));
  ASSERT_WP2_OK(sr.Allocate());

  ASSERT_WP2_OK(sr.AddAdaptiveBit(kCluster, kSymbol, p0, p1));

  const double dec_cost = Decode(&sr, values, max_values);

  EXPECT_NEAR(enc.GetCost(), enc_cost, kCostError);
  EXPECT_NEAR(enc.GetCost(), dec_cost, kCostError);
}

TEST_P(SymbolCostTest, AdaptiveSymCost) {
  ANSInit();
  const bool use_max = GetParam();

  for (int method = 0; method < (int)ANSAdaptiveSymbol::Method::kNum;
       ++method) {
    SCOPED_TRACE(SPrintf("method %d", method));
    const uint32_t speed = (method == (int)ANSAdaptiveSymbol::Method::kAOM)
                               ? kANSAProbaInvalidSpeed
                               : 42;

    SymbolsInfo info;
    info.SetInfo(kSymbol, /*min=*/0, /*max=*/kRange - 1, /*num_clusters=*/1,
                 StorageMethod::kAdaptiveSym);

    SymbolWriterForTest sw;
    ASSERT_WP2_OK(sw.Init(info));
    ASSERT_WP2_OK(sw.Allocate());
    ASSERT_WP2_OK(sw.AddAdaptiveSymbol(
        kCluster, kSymbol, (ANSAdaptiveSymbol::Method)method, speed));

    std::vector<int32_t> values;
    std::vector<uint32_t> max_values;
    double enc_cost = 0;
    ANSEnc enc;
    for (uint32_t i = 0; i < kNumValues; ++i) {
      const uint32_t max =
          use_max ? random_.Get<uint32_t>(0, kRange - 1) : kRange - 1;
      const uint32_t v = random_.Get<uint32_t>(0, max);
      values.push_back(v);
      if (use_max) {
        max_values.push_back(max);
        sw.ProcessWithCost(kCluster, kSymbol, v, max, "label", &enc, &enc_cost);
      } else {
        sw.ProcessWithCost(kCluster, kSymbol, v, "label", &enc, &enc_cost);
      }
    }

    ASSERT_WP2_OK(enc.Assemble());
    ExternalDataSource data_source(enc.Buffer(), enc.BufferSize());
    ANSDec dec(&data_source);

    SymbolReaderForTest sr;
    ASSERT_WP2_OK(sr.Init(info, &dec));
    ASSERT_WP2_OK(sr.Allocate());

    ASSERT_WP2_OK(sr.AddAdaptiveSymbol(
        kCluster, kSymbol, (ANSAdaptiveSymbol::Method)method, speed));

    const double dec_cost = Decode(&sr, values, max_values);

    EXPECT_NEAR(enc.GetCost(), enc_cost, kCostError);
    EXPECT_NEAR(enc.GetCost(), dec_cost, kCostError);
  }
}

INSTANTIATE_TEST_SUITE_P(SymbolReaderTestInstantiation, SymbolCostTest,
                         ::testing::Values(false, true));

//------------------------------------------------------------------------------

// Test an edge case with the capped SymbolWriter/SymbolReader.
#if !(defined(WP2_BITTRACE) || defined(WP2_TRACE) || defined(WP2_ENC_DEC_MATCH))
// Only test when there is no tracing as we will artificially set the state of
// the ANS decoder with a U value and read it as a symbol.
TEST(SymbolsTest, MaxValuePrecision) {
  // Define the faulty statistics.
  constexpr uint32_t kRange = 6;
  SymbolsInfoTest info;
  info.SetInfo(0, /*min=*/0, /*max=*/kRange - 1, /*num_clusters=*/1,
               StorageMethod::kAuto);
  constexpr uint32_t counts[kRange] = {4096, 2048, 8192, 512, 1024, 512};

  // Create the symbol writer.
  std::unique_ptr<WP2::SymbolWriter>
      sw(new (WP2Allocable::nothrow) WP2::SymbolWriter);
  ASSERT_TRUE(sw != nullptr);
  ANSEnc enc;
  ANSDictionaries dicts;
  ASSERT_WP2_OK(sw->Init(info));
  ASSERT_WP2_OK(sw->Allocate());
  constexpr uint32_t max_nnz = 10000;
  ASSERT_WP2_OK(sw->WriteHeader(/*cluster=*/0, max_nnz, /*sym=*/0,
                                counts, "counts", &enc, &dicts));

  // This value is actually the value of the ANS state that can crash the capped
  // decoding. It is right at the end of the last interval, which is not
  // properly defined: it goes up to 16382 if we use the usual roundings, but
  // we want it to go to the end of the tab size.
  enc.PutUValue(ANS_TAB_SIZE - 1, ANS_LOG_TAB_SIZE, "hack");
  // No need to write data as the stats are stored anyway. Just wrap up.
  EXPECT_WP2_OK(enc.Assemble());
  const size_t size = enc.BufferSize();
  const uint8_t* buf = enc.Buffer();

  // Decode the symbol reader statistics.
  ExternalDataSource data_source(buf, size);
  ANSDec dec(&data_source);
  SymbolReader sr;
  ASSERT_WP2_OK(sr.Init(info, &dec));
  ASSERT_WP2_OK(sr.Allocate());
  ASSERT_WP2_OK(sr.ReadHeader(/*cluster=*/0, max_nnz, /*sym=*/0, "counts"));

  // Check the decoder does not crash.
  int32_t value;
  ASSERT_WP2_OK(
      sr.TryRead(/*cluster=*/0, /*sym=*/0, /*max_value=*/1, "sym", &value));
}
#endif

}  // namespace
}  // namespace WP2
