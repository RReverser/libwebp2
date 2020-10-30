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
//  Partitioning of a source frame into blocks.
//
// Author: Skal (pascal.massimino@gmail.com)
//

#include <algorithm>
#include <cmath>

#include "src/common/integral.h"
#include "src/common/lossy/block_size_io.h"
#include "src/enc/analysis.h"
#include "src/enc/partition_score_func.h"
#include "src/enc/partitioner.h"
#include "src/enc/wp2_enc_i.h"

// Uncomment the following to collect and print block-type statistics
// #define BLOCK_STATS

#if defined(BLOCK_STATS)
#include "src/utils/stats.h"
#endif

namespace WP2 {

constexpr uint32_t kRDOptAreaSize = kMaxBlockSizePix;

//------------------------------------------------------------------------------

WP2Status ExtractBlockPartition(const EncoderConfig& config,
                                const GlobalParams& gparams,
                                const YUVPlane& yuv, const Rectangle& tile_rect,
                                VectorNoCtor<Block>* const blocks) {
  const PartitionMethod final_partition_method =
      FinalPartitionMethod(config, tile_rect.width, tile_rect.height);
  std::unique_ptr<PartitionScoreFunc> score_func;
  TileScoreFunc* tile_score_func = nullptr;
  MultiScoreFunc* multi_score_func = nullptr;
  switch (final_partition_method) {
    case MULTIPASS_PARTITIONING:
      multi_score_func = new (WP2Allocable::nothrow) MultiScoreFunc();
      score_func.reset(multi_score_func);
      break;
    case BLOCK_ENCODE_PARTITIONING:
    case SPLIT_RECURSE_PARTITIONING:
      score_func.reset(new (WP2Allocable::nothrow) BlockScoreFunc());
      break;
    case AREA_ENCODE_PARTITIONING:
      score_func.reset(new (WP2Allocable::nothrow)
                           AreaScoreFunc(kRDOptAreaSize, kRDOptAreaSize));
      break;
    case SUB_AREA_ENCODE_PARTITIONING:
      score_func.reset(new (WP2Allocable::nothrow) SubAreaScoreFunc(
          kRDOptAreaSize, kRDOptAreaSize));
      break;
    case TILE_ENCODE_PARTITIONING:
    case EXHAUSTIVE_PARTITIONING:
      tile_score_func = new (WP2Allocable::nothrow) TileScoreFunc();
      score_func.reset(tile_score_func);
      break;
    case ALL_4X4_PARTITIONING:
      score_func.reset(new (WP2Allocable::nothrow) FixedSizeScoreFunc(BLK_4x4));
      break;
    case ALL_8X8_PARTITIONING:
      score_func.reset(new (WP2Allocable::nothrow) FixedSizeScoreFunc(BLK_8x8));
      break;
    case ALL_16X16_PARTITIONING:
      score_func.reset(new (WP2Allocable::nothrow)
                           FixedSizeScoreFunc(BLK_16x16));
      break;
    case ALL_32X32_PARTITIONING:
      score_func.reset(new (WP2Allocable::nothrow)
                           FixedSizeScoreFunc(BLK_32x32));
      break;
    default:
      assert(false);
  }
  WP2_CHECK_ALLOC_OK(score_func != nullptr);

  std::unique_ptr<Partitioner> partitioner;
  switch (final_partition_method) {
    case MULTIPASS_PARTITIONING:
      partitioner.reset(new (WP2Allocable::nothrow)
                            MultiPassPartitioner(multi_score_func));
      break;
    case BLOCK_ENCODE_PARTITIONING:
    case TILE_ENCODE_PARTITIONING:
    case ALL_4X4_PARTITIONING:
    case ALL_8X8_PARTITIONING:
    case ALL_16X16_PARTITIONING:
    case ALL_32X32_PARTITIONING:
      partitioner.reset(new (WP2Allocable::nothrow) TopLeftBestPartitioner());
      break;
    case SPLIT_RECURSE_PARTITIONING:
      partitioner.reset(new (WP2Allocable::nothrow) SplitRecursePartitioner());
      break;
    case AREA_ENCODE_PARTITIONING:
      partitioner.reset(new (WP2Allocable::nothrow) AreaRDOptPartitioner(
          kRDOptAreaSize, kRDOptAreaSize));
      break;
    case SUB_AREA_ENCODE_PARTITIONING:
      partitioner.reset(new (WP2Allocable::nothrow) SubAreaRDOptPartitioner(
          kRDOptAreaSize, kRDOptAreaSize));
      break;
    case EXHAUSTIVE_PARTITIONING:
      assert(tile_score_func != nullptr);
      partitioner.reset(new (WP2Allocable::nothrow)
                            ExhaustivePartitioner(tile_score_func));
      break;
    default:
      assert(false);
  }
  WP2_CHECK_ALLOC_OK(partitioner != nullptr);

  WP2_CHECK_STATUS(score_func->Init(config, tile_rect, yuv, gparams));

  WP2_CHECK_STATUS(partitioner->Init(config, yuv, tile_rect, score_func.get()));
  WP2_CHECK_STATUS(partitioner->GetBestPartition(blocks));

#if defined(BLOCK_STATS)
  static WP2::Stats<const char*> blk_stats("block-type", "%s",
                                           /*sort_by_count=*/true);
  for (const Block& b : *out) blk_stats.Add(kDimNames[b.dim]);
#endif

  return WP2_STATUS_OK;
}

//------------------------------------------------------------------------------

}  // namespace WP2
