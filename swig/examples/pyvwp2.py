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
"""Simple python3 wp2 viewer.

Display a wp2 image from a given path or URL.

  Required dependencies:

    sudo apt install python3-tk python3-pil.imagetk

  Typical usage example:

  export PYTHONPATH=$PYTHONPATH:/path/to/libwebp2/module
  python3 examples/pyvwp2.py /path/to/input.wp2
"""

import argparse
import datetime
from PIL import Image
from PIL import ImageTk
import requests
import tkinter
import libwebp2 as wp2

# ------------------------------------------------------------------------------


class Viewer(object):
  """Main class.

  Call loop() after creation.
  """

  def __init__(self, file_name):
    self._file_name = file_name

    self._window = tkinter.Tk()
    self._window.title('Wp2 viewer')
    self._window.configure(background='grey')
    self._window.bind('<Escape>', lambda e: self._window.destroy())
    self._window.bind('<q>', lambda e: self._window.destroy())

    self._file_label = tkinter.Label(self._window, text=self._file_name)
    self._file_label.pack(side='top', fill='both', expand='no')

    self._info_label = tkinter.Label(self._window, text='Initialization')
    self._info_label.pack(side='top', fill='both', expand='no')

    self._img_label = tkinter.Label(self._window, text='loading')
    self._img_label.pack(side='bottom', fill='both', expand='yes')
    self._incr_dec_step_ms = 10  # 0 skips UI refresh, use >= 1 ms.

    # These might get garbage collected if not referenced.
    self._img = None
    self._img_tk = None
    self._checkerboard = None

    # Incremental decoding.
    self._input_bytes = None
    self._buffer_argb = None
    self._decoder = None
    self._last_status = wp2.WP2_STATUS_NOT_ENOUGH_DATA

    # Animations.
    self._expected_elapsed_ms = 0
    self._first_frame_timestamp = datetime.datetime.now()
    self._loop_index = 0

    # URL.
    self._request = None
    self._chunk_generator = None

  def _init_decoder(self):
    self._decoder = wp2.ArrayDecoder()
    self._decoder.SetInput(self._input_bytes)
    self._buffer_argb = self._decoder.GetPixels()
    self._expected_elapsed_ms = 0
    self._first_frame_timestamp = datetime.datetime.now()

  def _create_checkerboard(self, img_size):
    """Draws into the 'self._checkerboard' buffer."""
    bytearray_rgba = bytearray(img_size[0] * img_size[1] * 4)
    for y in range(img_size[1]):
      for x in range(img_size[0]):
        color = 255 if ((int(x) // 8) % 2 == (int(y) // 8) % 2) else 200
        bytearray_rgba[(y * img_size[0] + x) * 4 + 0] = color
        bytearray_rgba[(y * img_size[0] + x) * 4 + 1] = color
        bytearray_rgba[(y * img_size[0] + x) * 4 + 2] = color
        bytearray_rgba[(y * img_size[0] + x) * 4 + 3] = 255
    bytes_rgba = bytes(bytearray_rgba)
    self._checkerboard = Image.frombytes('RGBA', img_size, bytes_rgba)

  def _display_bytes(self, has_alpha):
    """Updates the window with the image stored in 'self._buffer_argb'."""
    if not self._buffer_argb.IsEmpty():
      # Note: converting the whole buffer to WP2_RGBA_32 at every chunk is not
      # efficient.
      buffer_rgba = wp2.ArgbBuffer(wp2.WP2_RGBA_32)
      buffer_rgba.ConvertFrom(self._buffer_argb)
      bytes_rgba = buffer_rgba.GetBytes()
      img_size = (buffer_rgba.width, buffer_rgba.height)

      self._img = Image.frombytes('RGBA', img_size, bytes_rgba)
      if has_alpha:
        if self._checkerboard is None:
          self._create_checkerboard(img_size)
        # Tkinter does not support transparency so the checkerboard is pasted
        # under the image.
        self._img = Image.alpha_composite(self._checkerboard, self._img)

      self._img_tk = ImageTk.PhotoImage(self._img)
      self._img_label.configure(image=self._img_tk)
      # Note: it would be possible to get rid of PIL dependency thanks to
      # tkinter.PhotoImage if it handled alpha.

  def _get_elapsed_ms_since_first_frame(self):
    timestamp_now = datetime.datetime.now()
    duration = (timestamp_now - self._first_frame_timestamp)
    return int(duration.total_seconds() * 1000)

  def _update(self):
    """Continues the incremental decoding, frame by frame if needed."""
    if self._get_elapsed_ms_since_first_frame() >= self._expected_elapsed_ms:
      read_frame = self._decoder.ReadFrame()
      if read_frame:
        self._expected_elapsed_ms += self._decoder.GetFrameDurationMs()
      features = self._decoder.TryGetDecodedFeatures()

      if read_frame and self._decoder.TryGetFrameDecodedFeatures(
          self._decoder.GetCurrentFrameIndex()).is_last:
        if not features.is_animation:
          return
        if features.loop_count > 0:
          self._loop_index += 1
          if self._loop_index >= features.loop_count:
            return
        self._init_decoder()  # Rewind animation to first frame.

      info = '{} bytes, {}x{}'.format(
          len(self._input_bytes),
          self._decoder.GetPixels().width, self._decoder.GetPixels().height)
      self._info_label.configure(text=info)

      has_alpha = features.has_alpha if features is not None else False
      self._display_bytes(has_alpha)

    wait_for = (
        self._expected_elapsed_ms - self._get_elapsed_ms_since_first_frame())
    self._window.after(max(wait_for, self._incr_dec_step_ms), self._update)

  def _receive_chunk(self):
    """Receives the next chunk of data requested by 'self._request'."""
    chunk = next(self._chunk_generator, None)
    if chunk:
      # It could be more efficient to decode chunk by chunk with StreamDecoder
      # but the decoded frames or the chunks would be stored for the next loop.
      self._input_bytes += chunk
      self._window.after(self._incr_dec_step_ms, self._receive_chunk)

  def _run(self):
    """Reads or fetches image, depending on '_file_name' being a file or URL."""
    if self._file_name.lower().startswith(('http://', 'https://')):
      self._input_bytes = bytes()
      self._request = requests.get(self._file_name, stream=True)
      if self._request.status_code != 200:
        print('error: status=', self._request.status_code)
        exit(1)
      elif self._request.headers['Content-type'] != 'image/wp2':
        print('error: Content-type=', self._request.headers['Content-type'])
        exit(1)
      self._chunk_generator = self._request.iter_content(
          chunk_size=None, decode_unicode=False)
      self._window.after(self._incr_dec_step_ms, self._receive_chunk)
    else:
      with open(self._file_name, 'rb') as input_file:
        self._input_bytes = input_file.read()
    self._init_decoder()
    self._update()

  def loop(self):
    self._window.after(self._incr_dec_step_ms, self._run)
    self._window.mainloop()


# ------------------------------------------------------------------------------


def main():
  parser = argparse.ArgumentParser(description='Simple wp2 viewer.')

  parser.add_argument(
      '-v',
      '--version',
      action='store_true',
      help='Print libwebp2 version.',
      dest='version')

  parser.add_argument('input', help='Path or URL of the image to display.')

  args = parser.parse_args()

  if args.version:

    def format_version(v, n):  # Returns a string with n numbers: 'X.Y.Z'
      return '.'.join([str((v >> (8 * (n - i - 1))) & 0xff) for i in range(n)])

    print('libwebp2 version:     ' + format_version(wp2.WP2GetVersion(), 3))
    print('libwebp2 ABI version: ' + format_version(wp2.WP2GetABIVersion(), 2))
    exit(0)

  viewer = Viewer(args.input)
  viewer.loop()


if __name__ == '__main__':
  main()
