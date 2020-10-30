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
"""Simple python3 wp2 decoder.

Decode wp2 format images.

  Typical usage example:

  export PYTHONPATH=$PYTHONPATH:/path/to/libwebp2/module
  python3 examples/dwp2.py /path/to/input.wp2 --output /path/to/output.png
"""

import argparse
import os
import libwebp2 as wp2


def main():
  parser = argparse.ArgumentParser(description='Simple wp2 decoder.')

  parser.add_argument(
      '-v',
      '--version',
      action='store_true',
      help='Print libwebp2 version.',
      dest='version')

  parser.add_argument('input', nargs='?', help='Path to the image to decode.')
  parser.add_argument(
      '-o', '--output', help='Path to the decoded image.', dest='output')
  parser.add_argument(
      '-f',
      '--frames_folder',
      type=str,
      help='Path to the folder where frames will be saved. '
      'Outputs: frame0_[duration]ms.png frame1_[duration]ms.png ...',
      dest='frames_folder')

  args = parser.parse_args()

  # ----------------------------------------------------------------------------

  if args.version:

    def format_version(v, n):  # Returns a string with n numbers: 'X.Y.Z'
      return '.'.join([str((v >> (8 * (n - i - 1))) & 0xff) for i in range(n)])

    print('libwebp2 version:     ' + format_version(wp2.WP2GetVersion(), 3))
    print('libwebp2 ABI version: ' + format_version(wp2.WP2GetABIVersion(), 2))
    exit(0)

  # ----------------------------------------------------------------------------

  if not args.input:
    print('error: the following arguments are required: input')
    exit(1)

  buffer_argb = wp2.ArgbBuffer()
  wp2.ReadImage(args.input, buffer_argb)

  # ----------------------------------------------------------------------------

  if args.output:
    wp2.SaveImage(buffer_argb, args.output, True)
  elif not args.frames_folder:
    print('Decoding success but no output file specified.')

  # ----------------------------------------------------------------------------

  if args.frames_folder:
    with open(args.input, 'rb') as input_file:
      encoded_bytes = input_file.read()
      decoder = wp2.ArrayDecoder()
      decoder.SetInput(encoded_bytes)
      frame_index = 0
      while decoder.ReadFrame():
        frame_file_name = 'frame{}_{}ms.png'.format(
            frame_index, decoder.GetFrameDurationMs())
        frame_file_path = os.path.join(args.frames_folder, frame_file_name)
        wp2.SaveImage(decoder.GetPixels(), frame_file_path, True)
        frame_index += 1
      decoder.GetStatus()  # Exception-triggering check

  # ----------------------------------------------------------------------------


if __name__ == '__main__':
  main()
