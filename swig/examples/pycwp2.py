#!/usr/bin/env python3
# Copyright 2020 Google LLC
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#      http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

# python3 (for gpylint)
"""Simple python3 wp2 encoder.

Encode images into wp2 format.

  Typical usage example:

  export PYTHONPATH=$PYTHONPATH:/path/to/libwebp2/module
  python3 examples/cwp2.py /path/to/input.png --output /path/to/output.wp2
                           --quality 80 --speed 7
"""

import argparse
import libwebp2 as wp2


def main():
  parser = argparse.ArgumentParser(description='Simple wp2 encoder.')

  parser.add_argument(
      '-v',
      '--version',
      action='store_true',
      help='Print libwebp2 version.',
      dest='version')

  parser.add_argument('input', nargs='?', help='Path to the image to encode.')
  parser.add_argument(
      '-f',
      '--frames',
      nargs='+',
      help='Paths to the images to encode as an animation, interleaved with '
      'durations in milliseconds. Usage: frame1.png 12 frame2.png 16 ...',
      dest='frames')
  parser.add_argument(
      '-o', '--output', help='Path to the encoded image.', dest='output')

  parser.add_argument(
      '-q',
      '--quality',
      default='75.0',
      type=float,
      help='Quality: from 0 (lossy) to 100 (lossless).',
      dest='quality')
  parser.add_argument(
      '-s',
      '--speed',
      default='6',
      type=int,
      help='Speed: from 0 (faster) to 9 (slower, better).',
      dest='speed')

  args = parser.parse_args()

  # ----------------------------------------------------------------------------

  if args.version:

    def format_version(v, n):  # Returns a string with n numbers: 'X.Y.Z'
      return '.'.join([str((v >> (8 * (n - i - 1))) & 0xff) for i in range(n)])

    print('libwebp2 version:     ', format_version(wp2.WP2GetVersion(), 3))
    print('libwebp2 ABI version: ', format_version(wp2.WP2GetABIVersion(), 2))
    exit(0)

  # ----------------------------------------------------------------------------

  if not wp2.WP2CheckVersion():
    print('error: version mismatch')
    exit(1)

  # ----------------------------------------------------------------------------

  config = wp2.EncoderConfig()
  config.quality = args.quality
  config.speed = args.speed
  if not config.IsValid():
    print('error: invalid configuration')
    exit(1)
  memory_writer = wp2.MemoryWriter()

  # ----------------------------------------------------------------------------

  if args.frames:
    if args.input:
      print('error: \'input\' argument must not be specified with --frames')
      exit(1)
    if (len(args.frames) % 2) != 0:
      print('error: --frames requires an even number of values')
      exit(1)
    num_frames = len(args.frames) // 2
    animation_encoder = wp2.AnimationEncoder()
    for f in range(num_frames):
      file_name = args.frames[f * 2]
      duration_ms = int(args.frames[f * 2 + 1])
      buffer_argb = wp2.ArgbBuffer()
      wp2.ReadImage(file_name, buffer_argb)
      animation_encoder.AddFrame(buffer_argb, duration_ms)
    animation_encoder.Encode(memory_writer, config)
  else:
    if not args.input:
      print('error: the following arguments are required: input')
      exit(1)
    original_image = wp2.ArgbBuffer()
    image_reader = wp2.ImageReader(args.input, original_image)
    [_, is_last, duration_ms] = image_reader.ReadFrame()
    if duration_ms == wp2.ImageReader.kInfiniteDuration:
      wp2.Encode(original_image, memory_writer, config)
    else:
      animation_encoder = wp2.AnimationEncoder()
      animation_encoder.AddFrame(original_image, duration_ms)
      while not is_last:
        [_, is_last, duration_ms] = image_reader.ReadFrame()
        animation_encoder.AddFrame(original_image, duration_ms)
      animation_encoder.Encode(memory_writer, config)

  # ----------------------------------------------------------------------------

  if args.output:
    wp2.IoUtilWriteFile(memory_writer.mem_, memory_writer.size_,
                        args.output, True)

  # ----------------------------------------------------------------------------


if __name__ == '__main__':
  main()
