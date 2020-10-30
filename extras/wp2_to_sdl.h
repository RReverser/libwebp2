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

#ifndef WP2_EXTRAS_WP2_TO_SDL_H_
#define WP2_EXTRAS_WP2_TO_SDL_H_

// Exports the method WP2ToSDL(const char* data, int data_size) which decodes
// a WP2 bitstream into an RGBA SDL surface.
// Return false on failure.
extern "C" int WP2ToSDL(const char* const data, unsigned int data_size);

#endif  // WP2_EXTRAS_WP2_TO_SDL_H_
