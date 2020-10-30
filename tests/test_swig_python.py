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
"""Unit test for SWIG python interface.

Depends on libwebp2 module.
"""

import os
import unittest
import libwebp2 as wp2


def GetFilePath(file_name):
  return os.path.join('./testdata', file_name)


class Libwp2Test(unittest.TestCase):

  def Encode(self, quality, speed):
    file_path = GetFilePath('source1_64x48.png')

    original_image = wp2.ArgbBuffer()
    wp2.ReadImage(file_path, original_image)

    config = wp2.EncoderConfig()
    config.quality = quality
    config.speed = speed
    self.assertTrue(config.IsValid())

    memory_writer = wp2.MemoryWriter()
    wp2.Encode(original_image, memory_writer, config)

    return memory_writer.GetBytes()

  def Decode(self, encoded_bytes):
    decoded_argb = wp2.ArgbBuffer()
    wp2.Decode(encoded_bytes, decoded_argb)

  def testVersion(self):
    self.assertGreaterEqual(wp2.WP2GetVersion(), 1)

  def testEncodeDecode(self):
    # Test different combinations of quality, speed.
    self.Decode(self.Encode(0, 5))
    self.Decode(self.Encode(50, 0))
    self.Decode(self.Encode(75, 9))
    self.Decode(self.Encode(99, 7))
    self.Decode(self.Encode(100, 3))

  def testIncrementalDecoder(self):
    encoded_bytes = self.Encode(80, 4)
    self.assertGreaterEqual(len(encoded_bytes), 12)
    half_size = len(encoded_bytes) // 2  # Integer division.
    first_half = encoded_bytes[:half_size]
    second_half = encoded_bytes[half_size:]

    stream_decoder = wp2.StreamDecoder()
    stream_decoder.AppendInput(first_half)
    self.assertFalse(stream_decoder.ReadFrame())
    stream_decoder.AppendInput(second_half)
    self.assertTrue(stream_decoder.ReadFrame())
    stream_decoder.GetStatus()  # Exception-raising check

    array_decoder = wp2.ArrayDecoder()
    array_decoder.SetInput(first_half)
    self.assertFalse(array_decoder.ReadFrame())
    array_decoder.SetInput(encoded_bytes)
    self.assertTrue(array_decoder.ReadFrame())
    array_decoder.GetStatus()  # Exception-raising check

    decoded_area = array_decoder.GetDecodedArea()
    self.assertGreater(decoded_area.width, 0)
    self.assertGreater(decoded_area.height, 0)

  def testException(self):
    image = wp2.ArgbBuffer()
    self.assertRaises(RuntimeError, wp2.ReadImage, 'missing/file', image)
    self.assertRaises(RuntimeError, wp2.Decode, None, None)


if __name__ == '__main__':
  unittest.main()
