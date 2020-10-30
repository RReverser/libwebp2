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
//  Simple WP2-to-SDL wrapper. Useful for emscripten.
//
// Author: Skal (pascal.massimino@gmail.com)

#ifdef HAVE_CONFIG_H
#include "src/wp2/config.h"
#endif

#if defined(WP2_HAVE_SDL)

#include "./wp2_to_sdl.h"

#include <cstdio>

#include "wp2/decode.h"

#if defined(WP2_HAVE_JUST_SDL_H)
#include <SDL.h>
#else
#include <SDL/SDL.h>
#endif

static bool init_ok = false;

int WP2ToSDL(const char* const data, unsigned int data_size) {
  int ok = 0;

  if (!WP2CheckVersion()) return 0;

  if (!init_ok) {
    SDL_Init(SDL_INIT_VIDEO);
    init_ok = true;
  }

  WP2::BitstreamFeatures bitstream;
  if (bitstream.Read((const uint8_t*)data, data_size) != WP2_STATUS_OK) {
    return 0;
  }

  SDL_Surface* const screen =
      SDL_SetVideoMode(bitstream.width, bitstream.height, 32, SDL_SWSURFACE);
  if (screen == NULL) {
    fprintf(stderr, "Unable to set video mode (32bpp %dx%d)!\n",
            bitstream.width, bitstream.height);
    return 0;
  }

  SDL_Surface* const surface =
      SDL_CreateRGBSurface(SDL_SWSURFACE,
                           bitstream.width, bitstream.height, 32,
                           0x000000ffu,   // R mask
                           0x0000ff00u,   // G mask
                           0x00ff0000u,   // B mask
                           0xff000000u);  // A mask
  if (surface == NULL) {
    fprintf(stderr, "Unable to create %dx%d RGBA surface!\n",
            bitstream.width, bitstream.height);
    return 0;
  }

  if (SDL_MUSTLOCK(surface)) SDL_LockSurface(surface);

#if SDL_BYTEORDER == SDL_BIG_ENDIAN
  const WP2SampleFormat format = WP2_BGRA_32;
#else
  const WP2SampleFormat format = WP2_RGBA_32;
#endif
  WP2Status status;
  WP2::ArgbBuffer output(format);
  status = output.SetExternal((uint32_t)surface->w, (uint32_t)surface->h,
                              (uint8_t*)surface->pixels,
                              (uint32_t)surface->pitch);
  if (status == WP2_STATUS_OK) {
    status = WP2::Decode((const uint8_t*)data, data_size, &output);
  }
  if (SDL_MUSTLOCK(surface)) SDL_UnlockSurface(surface);
  if (status != WP2_STATUS_OK) goto Error;
  if (SDL_BlitSurface(surface, NULL, screen, NULL) || SDL_Flip(screen)) {
    goto Error;
  }

  ok = 1;

 Error:
  SDL_FreeSurface(surface);
  SDL_FreeSurface(screen);
  return ok;
}

//------------------------------------------------------------------------------

#else

int WP2ToSDL(const char* const data, unsigned int data_size) {
  (void)data;
  (void)data_size;
  return 0;
}

#endif  // WP2_HAVE_SDL
