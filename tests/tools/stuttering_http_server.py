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

"""Simple HTTP server with limited bandwidth.

Sends the file requested with GET at a constant bitrate.

Typical usage example:
  python3 stuttering_http_server.py -f /path/to/the/folder/to/serve
Now resources can be accessed through:
  http://127.0.0.1:8080/test.png
Specify a bandwidth of 256 bytes sent every 1000 ms with:
  http://127.0.0.1:8080/test.png?bytes=256&ms=1000
Showcase (type in a browser):
  http://127.0.0.1:8080/test.png?display=true
"""

import argparse
import http.server
import os
import socket
import time
import urllib.parse

parser = argparse.ArgumentParser(description='Stuttering HTTP server.')

parser.add_argument(
    '-f',
    '--folder',
    type=str,
    default='.',
    help='Folder to serve.',
    dest='folder')

parser.add_argument(
    '-b', '--bytes', type=int, default=256, help='b bytes sent..', dest='bytes')

parser.add_argument(
    '-m', '--ms', type=float, default=1000, help='..every m ms.', dest='ms')

parser.add_argument(
    '-i', '--ip', type=str, default='localhost', help='IP address.', dest='ip')

parser.add_argument(
    '-p', '--port', type=int, default=8080, help='Port.', dest='port')

parser.add_argument(
    '--transfer',
    type=str,
    default='chunked',
    choices={'chunked', 'normal'},
    help='Transfer mode.',
    dest='transfer')

parser.add_argument(
    '--test', action='store_true', help='Quit after transfer.', dest='test')

args = parser.parse_args()

# ------------------------------------------------------------------------------


class HTTPRequestHandler(http.server.BaseHTTPRequestHandler):
  """Main class.

  Run it through http.server.HTTPServer() to answer requests.
  """

  # Use HTTP 1.1 as 1.0 doesn't support chunked encoding.
  if args.transfer == 'chunked':
    protocol_version = 'HTTP/1.1'

  def send_text(self, string):
    self.send_response(200)
    self.send_header('Content-type', 'text/html')
    self.send_header('Content-length', str(len(string)))
    self.end_headers()
    self.wfile.write(bytearray(string, 'utf8'))

  def do_GET(self):  # pylint: disable=invalid-name
    """Called when a GET request is received."""
    if self.path == '/favicon.ico':
      self.send_response(404)
      return

    # URL parsing.
    parsed_url = urllib.parse.urlparse(self.path)
    url_args = dict(urllib.parse.parse_qsl(parsed_url.query))
    if parsed_url.path == '/':  # Nothing requested, warn the user and return.
      url = 'http://' + args.ip + ':' + str(args.port) + '/test.png'
      self.send_text('Please pass a relative file path.<br>'
                     'Example:<br>'
                     '<a href=' + url + '>' + url + '</a>')
      return

    # File opening.
    file_path = os.path.normpath(args.folder + parsed_url.path)

    try:
      file = open(file_path, 'rb')
      data = file.read()
      file.close()

    except IOError:  # File not found, warn the user and return.
      self.send_text('File ' + file_path + ' does not exist.')
      return

    bandwidth_bytes = args.bytes
    bandwidth_ms = args.ms
    # Override bandwidth parameters if supplied in URL.
    url_args_bandwidth_bytes = url_args.get('bytes')
    if url_args_bandwidth_bytes:
      bandwidth_bytes = int(url_args_bandwidth_bytes)
    url_args_bandwidth_ms = url_args.get('ms')
    if url_args_bandwidth_ms:
      bandwidth_ms = int(url_args_bandwidth_ms)

    # Image display.
    url_args_display_image = url_args.get('display')
    if url_args_display_image == 'true':
      url = ('http://' + args.ip + ':' + str(args.port) + '' + parsed_url.path +
             '?bytes=' + str(bandwidth_bytes) + '&ms=' + str(bandwidth_ms))
      self.send_text('<img src="' + url + '" />')
      return

    self.send_response(200)

    # Sending the bytes.
    file_extension = os.path.splitext(file_path)[1][1:]
    self.send_header('Content-type', 'image/' + file_extension.lower())
    if args.transfer == 'chunked':
      self.send_header('Transfer-encoding', 'chunked')  # Needed for streaming.
    else:
      self.send_header('Content-length', str(len(data)))
    self.end_headers()

    start_index = 0
    while start_index < len(data):
      end_index = start_index + bandwidth_bytes
      if end_index > len(data):
        end_index = len(data)

      chunk = data[start_index:end_index]
      if args.transfer == 'chunked':
        # See 'Transfer-Encoding: chunked' documentation for syntax.
        chunk = b'%X\r\n%s\r\n' % (len(chunk), chunk)
      self.wfile.write(chunk)

      time.sleep(bandwidth_ms / 1000)
      start_index += bandwidth_bytes
    if args.transfer == 'chunked':
      # Send an empty chunk to notify it's the end.
      self.wfile.write(b'0\r\n\r\n')
    self.wfile.flush()

    if args.test:
      self.server.running = False


class HTTPServerV6(http.server.HTTPServer):
  address_family = socket.AF_INET6
  running = True

  def serve_until_stopped(self):
    while self.running:
      self.handle_request()


httpd = HTTPServerV6((args.ip, args.port), HTTPRequestHandler)
print('Serving at http://' + args.ip + ':' + str(args.port))
httpd.serve_until_stopped()
