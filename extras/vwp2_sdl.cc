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
// Simple SDL-based WP2 file viewer.
// Does not support animation, just static images.
//
// Press 'q' to exit.
//
// Author: Skal (pascal.massimino@gmail.com)

#include <cstdio>

#ifdef HAVE_CONFIG_H
#include "wp2/config.h"
#endif

#if defined(WP2_HAVE_SDL)

#include "./wp2_to_sdl.h"
#include "imageio/imageio_util.h"
#include "wp2/decode.h"

#if defined(WP2_HAVE_JUST_SDL_H)
#include <SDL.h>
#else
#include <SDL/SDL.h>
#endif

static void ProcessEvents() {
  bool done = false;
  SDL_Event event;
  while (!done && SDL_WaitEvent(&event)) {
    if (event.type == SDL_KEYUP && event.key.keysym.sym == SDLK_q) {
      done = true;
      break;
    }
  }
}

int main(int argc, char* argv[]) {
  bool ok = false;
  for (int c = 1; c < argc; ++c) {
    if (!strcmp(argv[c], "-h")) {
      printf("Usage: %s [-h] image.wp2 [more_files.wp2...]\n", argv[0]);
      return 0;
    }
    const char* const file = argv[c];
    WP2::Data data;
    if (WP2::IoUtilReadFile(file, &data) != WP2_STATUS_OK) {
      fprintf(stderr, "Error opening file: %s\n", file);
      goto Error;
    }
    if (data.size != (size_t)(int)data.size) {
      fprintf(stderr, "File too large.\n");
      goto Error;
    }
    ok = !!WP2ToSDL((const char*)data.bytes, (int)data.size);
    if (!ok) {
      fprintf(stderr, "Error displaying file %s\n", file);
      goto Error;
    }
    data.Clear();   // free some memory
    ProcessEvents();
  }
  ok = 1;

 Error:
  SDL_Quit();
  return ok ? 0 : 1;
}

#else  // !WP2_HAVE_SDL

int main(int argc, const char *argv[]) {
  fprintf(stderr, "SDL support not enabled in %s.\n", argv[0]);
  (void)argc;
  return 0;
}

#endif
