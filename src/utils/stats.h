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
//   Simple program-wide histogram collection (not MT-safe!)
//
// Usage:
//  static Stats my_stats;   // <- relies on dtor for printing stats
//  ...
//  my_stats.Add(key, value);
//
// Author: Skal (pascal.massimino@gmail.com)

#ifndef WP2_UTILS_STATS_H_
#define WP2_UTILS_STATS_H_

#include <algorithm>
#include <cstdio>
#include <functional>   // for std::less<>
#include <unordered_map>
#include <vector>

#include "src/utils/utils.h"  // for WP2Print()

namespace WP2 {

template<typename K = uint32_t> struct Stats {
  explicit Stats(const char* const name = nullptr,
                 const char* const format = "%lu",   // how to print type "K"
                 bool sort_by_count = false)   // just sort and print the counts
      : name_(name), format_(format), sort_by_count_(sort_by_count) {}
  void Add(K key, uint64_t value = 1) {  // TODO(skal): mutex?
    const auto i = stats_.find(key);
    if (i == stats_.end()) {
      stats_[key].total_ = value;
      stats_[key].counts_ = 1;
    } else {
      i->second.total_ += value;
      i->second.counts_ += 1;
    }
    ++total_;
  }

  ~Stats() {
    if (total_ > 0) {
      printf(" === STATISTICS (%s) count: %lu ===\n",
             (name_ == nullptr) ? "undef" : name_,
             (unsigned long)total_);  // NOLINT
      auto stats = GetStats();
      const double norm = 100. / total_;
      for (const auto& i : stats) {
        printf("# ");
        WP2Print(format_, i.first);
        const uint64_t count = i.second.counts_;
        printf(" : average: %.1f", count ? 1.f * i.second.total_ / count : 0);
        printf(" \tcount: %lu \t(%.2f%%)\n",
               (unsigned long)count, count * norm);  // NOLINT
      }
    }
  }

 protected:
  const char* const name_;
  const char* const format_;
  bool sort_by_count_;
  uint64_t total_ = 0;

  struct Record { uint64_t total_ = 0, counts_ = 0; };
  typedef std::pair<K, Record> Entry;
  std::unordered_map<K, Record> stats_;

  std::vector<Entry> GetStats() const {
    std::vector<Entry> stats(stats_.begin(), stats_.end());
    std::sort(stats.begin(), stats.end(), sort_by_count_ ? cmp : cmp_first);
    return stats;
  }

  static bool cmp(const Entry&a, const Entry&b) {
    return (a.second.counts_ > b.second.counts_);
  }
  static bool cmp_first(const Entry&a, const Entry&b) {
    return std::less<K>()(a.first, b.first);
  }
};

}    // namespace WP2

#endif  /* WP2_UTILS_STATS_H_ */
