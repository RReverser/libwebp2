// Copyright 2020 Google LLC
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
// WebP encode.
//
// Author: Yannis Guyon (yguyon@google.com)

#include <cstdio>
#include <cmath>

#ifdef HAVE_CONFIG_H
#include "wp2/config.h"
#endif

#include "./image_enc.h"
#include "./imageio_util.h"
#include "src/dsp/math.h"
#include "src/utils/utils.h"
#include "src/wp2/format_constants.h"
#ifdef WP2_HAVE_WEBP
#include "webp/encode.h"
#endif

namespace WP2 {

// -----------------------------------------------------------------------------

#if defined(WP2_HAVE_WEBP)

static int WriterFunc(const uint8_t* data, size_t data_size,
                      const WebPPicture* picture) {
  Writer* const writer = reinterpret_cast<Writer*>(picture->custom_ptr);
  return writer->Append(data, data_size) ? 1 : 0;
}

WP2Status CompressWebP(const ArgbBuffer& buffer, const WebPConfig& config,
                       Writer* const output) {
  WP2_CHECK_OK(!buffer.IsEmpty(), WP2_STATUS_INVALID_PARAMETER);
  WP2_CHECK_OK(output != nullptr, WP2_STATUS_NULL_PARAMETER);
  WP2_CHECK_OK(WebPValidateConfig(&config), WP2_STATUS_INVALID_CONFIGURATION);
  // Convert WP2::ArgbBuffer to WebPPicture
  WebPPicture picture;
  WP2_CHECK_STATUS(ConvertToWebPPicture(buffer, &picture));
  picture.writer = WriterFunc;
  picture.custom_ptr = reinterpret_cast<void*>(output);
  WP2_CHECK_OK(WebPEncode(&config, &picture), WP2_STATUS_INVALID_PARAMETER);
  WebPPictureFree(&picture);
  return WP2_STATUS_OK;
}

WP2Status CompressWebP(const ArgbBuffer& buffer, float quality, int speed,
                       Writer* const output) {
  WebPConfig config;
  WP2_CHECK_OK(WebPConfigInit(&config), WP2_STATUS_VERSION_MISMATCH);
  // Map wp2 'quality' in [lossy:kMaxLossyQuality:kMaxQuality]
  // to WebP in [lossy:near-lossless:lossless].
  if (quality <= (float)kMaxLossyQuality) {
    config.quality = Clamp((quality / kMaxLossyQuality) * 100.f, 0.f, 100.f);
  } else {
    WP2_CHECK_OK(WebPConfigPreset(&config, WEBP_PRESET_DEFAULT, quality),
                 WP2_STATUS_VERSION_MISMATCH);
    config.lossless = 1;
    config.near_lossless =
        Clamp<int>(std::lround(((quality - kMaxLossyQuality) /
                                (kMaxQuality - kMaxLossyQuality)) *
                               100.f),
                   0, 100);
    config.exact = (config.near_lossless == 100);
    config.quality = Clamp(speed / 9.f * 100.f, 0.f, 100.f);
  }
  config.method = Clamp(DivRound(speed * 6, 9), 0, 6);
  WP2_CHECK_STATUS(CompressWebP(buffer, config, output));
  return WP2_STATUS_OK;
}

#else

WP2Status CompressWebP(const ArgbBuffer&, float, int, Writer*) {
  return WP2_STATUS_UNSUPPORTED_FEATURE;
}

#endif  // defined(WP2_HAVE_WEBP)

// -----------------------------------------------------------------------------

}  // namespace WP2
