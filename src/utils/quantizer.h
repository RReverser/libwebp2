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
// Defines utilities to quantize a distribution to better store it and the
// population it represents.
//
// Author: Vincent Rabaud (vrabaud@google.com)

#ifndef WP2_ENC_LOSSLESS_QUANTIZER_H_
#define WP2_ENC_LOSSLESS_QUANTIZER_H_

#include "src/utils/ans_utils.h"
#include "src/utils/vector.h"
#include "src/wp2/format_constants.h"

namespace WP2 {

class Quantizer {
 public:
  enum ConfigType { Raw, Huffman, HuffmanANS };
  // Config defining how a histogram is stored.
  // E.g.  [0,1,8,0,4]
  // If is_sparse_==true, it is stored as sparse, histogram_ contains [1,8,4]
  // and 'runs' the mappings differences minus 1: [2,1,2] (the first element is
  // a difference to -1) minus 1 hence [1,0,1]
  // histogram_to_write_ contains how histogram_ will be finally written, which
  // in Raw type is the same, but in Huffman is [0,3,2] (the powers of 2 used to
  // represent the histogram.
  // histogram_ represents a distribution of symbol and if we were to write
  // data with it, it would have an overal cost of cost_. Overall means the cost
  // of writing the symbols with probabilities defined by histogram_ but also
  // of storing histogram_.
  struct ConfigParam {
    // The following members differentiate the different configurations.
    ConfigType type_;
    bool is_sparse_;
    uint8_t max_freq_bits_;
  };
  struct HistogramSparse {
    uint32_t* counts;
    uint16_t* mapping;
    uint32_t nnz;  // number of non-zero counts
  };
  struct Config {
    ConfigParam param_;
    // The following members are the results of the above parameters.
    uint32_t* histogram_to_write_;
    size_t size_to_write_;
    HistogramSparse h_;
    float cost_;
    Config* next_;
  };

  Quantizer() {
    for (auto& c : configs_used_) c = false;
  }

  // range_max is the maximum range of the used symbols (maximum value + 1).
  WP2Status Allocate(uint32_t range_max);

  // Given a histogram of symbols, their range, as well as the bit size
  // of the number of pixels n_pixels_bits, find the best way to store the
  // symbols represented by 'histogram' with the data in 'histogram' by
  // minimizing their overall size. 'size_sparse' is an upper bound of the
  // number of non-zero counts. 'speed' is between 0 (fastest) and 9 (slowest).
  // 'max_count' is the maximum value a count can take.
  void Quantize(const uint32_t* const histogram, const uint16_t* const mapping,
                size_t size_sparse, uint32_t symbol_range, uint32_t max_count,
                int speed, Config** const config_best);

 private:
  // Free the config from the buffer to make it re-usable.
  void ResetConfig(const Config* c) {
    while (c != nullptr) {
      const size_t c_ind = c - &configs_[0];
      configs_used_[c_ind] = false;
      c = c->next_;
    }
  }

  // Implementation of the quantization that can call itself to compress its own
  // coefficients.
  // 'cost_max' is the value above which we can early exit.
  void QuantizeImpl(const uint32_t* const histogram,
                    const uint16_t* const mapping, size_t size_sparse,
                    uint32_t symbol_range, uint32_t n_pixels_bits, int speed,
                    size_t config_index,
                    float cost_max = std::numeric_limits<float>::max());

  // TODO(vrabaud) Investigate to see if it is worth having more levels.
  static constexpr size_t kConfigNbr = 2;
  Config configs_[kConfigNbr];
  bool configs_used_[kConfigNbr];
  // For each recursion level, each set of param can be Raw/Huffman, sparse or
  // not and have several levels of recursion.
  uint8_t bits_[kConfigNbr][kMaxFreqBits];
  Vector_u32 histogram_sub_[kConfigNbr];
  Vector_u16 mapping_sub_[kConfigNbr];

  uint32_t histogram_size_max_;
  uint32_t* histogram_quantized_[kConfigNbr];
  uint32_t* histogram_to_write_[kConfigNbr];
  Data buffer_;  // Big memory chunk storing the buffers above.

  VectorNoCtor<OptimizeArrayStorageStat> stats_buffer_;
};
}  // namespace WP2

#endif  // WP2_ENC_LOSSLESS_QUANTIZER_H_
