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

#include "src/utils/quantizer.h"

#include <algorithm>

#include "src/common/constants.h"
#include "src/dsp/math.h"
#include "src/utils/ans_utils.h"
#include "src/utils/utils.h"

namespace WP2 {

WP2Status Quantizer::Allocate(uint32_t range_max) {
  // We need at least kMaxFreqBits elements to store Huffman probabilities.
  histogram_size_max_ = std::max(kMaxFreqBits + 1u, range_max);
  const size_t histogram_size = 4 * sizeof(uint32_t) + sizeof(uint16_t);
  WP2_CHECK_STATUS(buffer_.Resize(
      kConfigNbr * histogram_size_max_ * histogram_size, /*keep_bytes=*/false));

  for (size_t i = 0; i < kConfigNbr; ++i) {
    // Make sure the recursion buffers are properly allocated.
    WP2_CHECK_ALLOC_OK(histogram_sub_[i].resize(kMaxFreqBits + 1));
    WP2_CHECK_ALLOC_OK(mapping_sub_[i].resize(kMaxFreqBits + 1));
  }

  WP2_CHECK_ALLOC_OK(stats_buffer_.resize(histogram_size_max_));
  return WP2_STATUS_OK;
}

void Quantizer::Quantize(const uint32_t* const histogram,
                         const uint16_t* const mapping, size_t size_sparse,
                         uint32_t symbol_range, uint32_t max_count, int speed,
                         Config** const config_best) {
  assert(speed >= 0 && speed <= 9);
  *config_best = &configs_[0];
  configs_used_[0] = true;
  for (size_t i = 1; i < kConfigNbr; ++i) configs_used_[i] = false;

  uint32_t* buffer = (uint32_t*)buffer_.bytes;
  for (size_t i = 0; i < kConfigNbr; ++i) {
    // Make sure the recursion buffers are properly allocated.
    histogram_quantized_[i] = buffer;
    buffer += histogram_size_max_;
    histogram_to_write_[i] = buffer;
    buffer += histogram_size_max_;
    configs_[i].histogram_to_write_ = buffer;
    buffer += histogram_size_max_;
    configs_[i].h_.counts = buffer;
    buffer += histogram_size_max_;
    configs_[i].h_.mapping = (uint16_t*)buffer;
    buffer += histogram_size_max_ * sizeof(uint16_t) / sizeof(uint32_t);
  }
  QuantizeImpl(histogram, mapping, size_sparse, symbol_range,
               1 + WP2Log2Floor(max_count), speed, 0);
}

// The quantization implementation uses the Config cache of the class during
// recursion.
void Quantizer::QuantizeImpl(const uint32_t* const histogram,
                             const uint16_t* const mapping, size_t size_sparse,
                             uint32_t symbol_range, uint32_t n_pixels_bits,
                             int speed, size_t config_index, float cost_max) {
  assert(config_index < kConfigNbr);
  Config* const config_best = &configs_[config_index];
  if (size_sparse == 0) {
    config_best->param_.type_ = Raw;
    config_best->param_.is_sparse_ = false;
    config_best->param_.max_freq_bits_ = 0;
    config_best->histogram_to_write_ = nullptr;
    config_best->size_to_write_ = 0;
    config_best->h_.counts = nullptr;
    config_best->h_.mapping = nullptr;
    config_best->h_.nnz = 0;
    config_best->cost_ = 0;
    config_best->next_ = nullptr;
    return;
  }

  // We simplify the mapping to be stored as the difference between
  // consecutive terms.
  config_best->h_.nnz = size_sparse;
  config_best->next_ = nullptr;
  std::copy(mapping, mapping + size_sparse, config_best->h_.mapping);

  // Pre-compute some quantization configurations for clarity.
  uint8_t* const bits = bits_[config_index];
  size_t n_max_freq_bits = 0;
  {
    uint32_t max_freq_bits_max =
        1 + WP2Log2Floor(*std::max_element(histogram, histogram + size_sparse));
    max_freq_bits_max = std::min(max_freq_bits_max, kMaxFreqBits);
    const uint32_t max_freq_bits_min =
        (1 * speed + max_freq_bits_max * (9 - speed)) / 9;
    for (uint32_t max_freq_bits = max_freq_bits_max;
         max_freq_bits >= max_freq_bits_min; --max_freq_bits) {
      bits[n_max_freq_bits++] = max_freq_bits;
    }
    // Start by doing some strong quantization to prune out the ones that would
    // barely change the bits cost.
    std::swap(bits[0], bits[(n_max_freq_bits - 1) / 2]);
    if (n_max_freq_bits > 1) {
      std::swap(bits[1], bits[(n_max_freq_bits - 1) / 4]);
    }
    if (n_max_freq_bits > 2) {
      std::swap(bits[2], bits[(n_max_freq_bits - 1) * 3 / 4]);
    }
  }

  uint32_t* histogram_quantized = histogram_quantized_[config_index];
  uint32_t* histogram_to_write = histogram_to_write_[config_index];

  // Compute some constant costs.
  float cost_probability_size, cost_mapping;
  if (symbol_range == size_sparse) {
    // If we use the whole range, actually force sparse usage. This seems
    // counter-intuitive but it's actually due to the optimization done for
    // sparse data: 1 is subtracted to all values.
    cost_probability_size = std::numeric_limits<float>::max();
    cost_mapping = 0.f;
  } else {
    cost_probability_size = WP2Log2(symbol_range + 1 - size_sparse);
    // Small speed improvements: if we use the whole range, no need to even
    // consider a potential mapping.
    cost_mapping = StoreMapping(config_best->h_.mapping, size_sparse,
                                symbol_range, stats_buffer_.data(), nullptr);
  }
  // Cost for the sparse bit and the histogram type.
  const float cost_sparse_bit_and_type = 1 + std::log2(3);
  n_pixels_bits = std::min(n_pixels_bits, kMaxFreqBits);

  // Go over the quantification configurations and pick the best one.
  config_best->cost_ = cost_max;
  for (size_t i = 0; i < n_max_freq_bits; ++i) {
    // No need to continue if any pre-computed cost we have is already bigger.
    if (cost_sparse_bit_and_type + cost_mapping > config_best->cost_ &&
        cost_sparse_bit_and_type + cost_probability_size > config_best->cost_) {
      break;
    }
    for (const ConfigType type : {Raw, Huffman}) {
      uint8_t max_freq_bits = bits[i];
      // Compute the cost of storing the data with the quantized probabilities.
      float cost_ini;

      if (type == Huffman) {
        ANSCountsQuantizeHuffman(size_sparse,
                                 histogram,
                                 histogram_quantized,
                                 max_freq_bits,
                                 &cost_ini);
      } else {
        const uint32_t freq_max = (1u << max_freq_bits) - 1;
        ANSCountsQuantize(false, freq_max, size_sparse, histogram,
                          histogram_quantized, &cost_ini);
      }

      // Add the cost for the sparse bit and the histogram type.
      cost_ini += cost_sparse_bit_and_type;

      for (const bool is_sparse : {false, true}) {
        // Reset variables.
        max_freq_bits = bits[i];
        float cost = cost_ini;
        // Add the cost for the mapping or the probability size.
        cost += is_sparse ? cost_mapping : cost_probability_size;

        // No need to continue if the bit cost without the histogram cost is
        // already bigger.
        if (cost >= config_best->cost_) continue;

        // Prepare the vector in which we will write the histogram.
        size_t histogram_to_write_size;
        if (is_sparse) {
          // Store the mapping and the probabilities: it's basically a sparse
          // vector.
          histogram_to_write_size = size_sparse;
        } else {
          // Store the full probabilities: the vector is not sparse and can
          // contain 0 values.
          histogram_to_write_size = mapping[size_sparse - 1] + 1;
          std::fill(histogram_to_write,
                    histogram_to_write + histogram_to_write_size, 0);
        }

        // Fill the histogram.
        Config* config_next = nullptr;
        float cost_histo;
        if (type == Huffman) {
          uint32_t max_n_bits = 0;
          for (uint32_t k = 0; k < size_sparse; ++k) {
            // +1 as we encode 0 as a frequency of 0.
            uint32_t n_bits =
                WP2Log2Floor(histogram_quantized[k]) + (is_sparse ? 0 : 1);
            if (is_sparse) {
              histogram_to_write[k] = n_bits;
            } else {
              histogram_to_write[mapping[k]] = n_bits;
            }
            max_n_bits = std::max(max_n_bits, n_bits);
          }
          max_freq_bits = 1 + WP2Log2Floor(max_n_bits);
          if (is_sparse) {
            cost_histo = StoreVector(
                histogram_to_write, histogram_to_write_size,
                (1 << max_freq_bits) - 1, stats_buffer_.data(), nullptr);
          } else {
            cost_histo = StoreVectorNnz(
                histogram_to_write, histogram_to_write_size, size_sparse,
                (1 << max_freq_bits) - 1, stats_buffer_.data(), nullptr);
          }
          cost_histo += WP2Log2Fast(1 + WP2Log2Floor(n_pixels_bits));
          // For Huffman, and huffman only, we try to compress the coefficients
          // again. We only do it for that mode as it has more compact data:
          // numbers in [0, kMaxFreqBits].
          for (size_t k = 0; k < kConfigNbr; ++k) {
            if (!configs_used_[k]) {
              config_next = &configs_[k];
              configs_used_[k] = true;
              break;
            }
          }
          if (config_next != nullptr) {
            // Define the histogram to use for recursion.
            uint32_t* const histogram_sub = histogram_sub_[config_index].data();
            uint16_t* const mapping_sub = mapping_sub_[config_index].data();
            const size_t histogram_sub_size_max = kMaxFreqBits + 1;

            std::fill(histogram_sub, histogram_sub + histogram_sub_size_max, 0);
            uint32_t nnz_left = size_sparse;
            for (uint32_t k = 0; k < histogram_to_write_size; ++k) {
              assert(histogram_to_write[k] < histogram_sub_size_max);
              if (!is_sparse && k + nnz_left == histogram_to_write_size) {
                // If we only have non-zeros left, we can store counts - 1.
                assert(histogram_to_write[k] > 0u);
                ++histogram_sub[histogram_to_write[k] - 1];
                --nnz_left;
              } else {
                ++histogram_sub[histogram_to_write[k]];
                if (histogram_to_write[k] > 0) --nnz_left;
              }
            }
            // Define the sparse version of the histogram.
            uint32_t ind = 0;
            for (size_t k = 0; k < histogram_sub_size_max; ++k) {
              if (histogram_sub[k] != 0) {
                assert(ind < histogram_sub_size_max);
                histogram_sub[ind] = histogram_sub[k];
                mapping_sub[ind] = k;
                ++ind;
              }
            }
            const float cost0 = WP2Log2Fast(1 + kMaxFreqBits);
            const float cost_max_tmp = cost_histo - cost0;
            const uint32_t histogram_to_write_size_bits =
                1 + WP2Log2Floor(histogram_to_write_size);
            QuantizeImpl(histogram_sub, mapping_sub, ind, 1 + kMaxFreqBits,
                         histogram_to_write_size_bits, speed,
                         config_next - configs_, cost_max_tmp);
            if (config_next->cost_ < cost_max_tmp) {
              // Also add to the cost the size of the array.
              config_next->cost_ += cost0;
              cost_histo = config_next->cost_;
            } else {
              // If we cannot do better, reset the used configurations.
              ResetConfig(config_next);
              config_next = nullptr;
            }
          }
        } else {
          assert(stats_buffer_.size() >= histogram_to_write_size);
          if (is_sparse) {
            for (uint32_t k = 0; k < size_sparse; ++k) {
              histogram_to_write[k] = histogram_quantized[k] - 1;
            }
            cost_histo = StoreVector(
                histogram_to_write, histogram_to_write_size,
                (1 << max_freq_bits) - 1, stats_buffer_.data(), nullptr);
          } else {
            for (uint32_t k = 0; k < size_sparse; ++k) {
              histogram_to_write[mapping[k]] = histogram_quantized[k];
            }
            cost_histo = StoreVectorNnz(
                histogram_to_write, histogram_to_write_size, size_sparse,
                (1 << max_freq_bits) - 1, stats_buffer_.data(), nullptr);
          }
          cost_histo += WP2Log2Fast(n_pixels_bits);
        }

        // Add the cost of the probabilities.
        cost += cost_histo;

        // Keep the best configuration.
        if (cost >= config_best->cost_) {
          ResetConfig(config_next);
          continue;
        }
        if (config_next == nullptr) {
          config_best->param_.type_ = type;
        } else {
          config_best->param_.type_ = HuffmanANS;
        }
        ResetConfig(config_best->next_);
        config_best->next_ = config_next;
        config_best->param_.is_sparse_ = is_sparse;
        config_best->param_.max_freq_bits_ = max_freq_bits;
        config_best->size_to_write_ = histogram_to_write_size;
        std::swap(config_best->histogram_to_write_, histogram_to_write);
        if (is_sparse) {
          std::swap(config_best->h_.counts, histogram_quantized);
        } else {
          std::copy(histogram_quantized,
                    histogram_quantized + histogram_to_write_size,
                    config_best->h_.counts);
        }
        config_best->cost_ = cost;
      }
    }
  }

  // Reset the common temporary variables.
  histogram_quantized_[config_index] = histogram_quantized;
  histogram_to_write_[config_index] = histogram_to_write;
}

}  // namespace WP2
