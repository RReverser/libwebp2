# SWIG interface for webp2.

Currently supported languages:

*   python

## Build

To build the python libwebp2 SWIG interface, enable the CMake flag
`WP2_BUILD_SWIG_PY`. The libwebp2 module will be in the specified build folder, in
subdirectory `swig`.

## Python examples

### From Python3 command line

```python
import sys
sys.path.append('/path/to/libwebp2/module')
import libwebp2
print('libwebp2 version: {}'.format(libwebp2.WP2GetVersion()))
```

### From bash comand line

    export PYTHONPATH=$PYTHONPATH:/path/to/libwebp2/module

To compress:

    python3 examples/pycwp2.py /path/to/input.png -o /path/to/output.wp2 -q 80

To decompress:

    python3 examples/pydwp2.py /path/to/input.wp2 -o /path/to/output.png

To view a wp2 from a file:

    python3 examples/pyvwp2.py /path/to/input.wp2

To simulate a slow bandwidth and view incremental wp2 decoding, run in parallel:

    python3 ../tests/tools/stuttering_http_server.py
            --folder /path/to/folder

then

    python3 examples/pyvwp2.py
            http://localhost:8080/test.wp2?bytes=1024&ms=1000
