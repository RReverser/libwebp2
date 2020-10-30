# This makefile is a simpler alternative to the autoconf-based build
# system, for simple local building of the libraries and tools.
# It will not install the libraries system-wide, but just create the 'cwp2'
# and 'dwp2' tools in the examples/ directory, along with the static
# libraries 'src/libwebp2.a', 'src/libwebp2decoder.a'
#
# To build the library and examples, use 'make'.

#### Customizable part ####

VERSION=0.1
ARCHIVE_FILE=wp2-$(VERSION).tar.gz

# These flags assume you have libpng, libjpeg, libtiff, libgif and libwebp
# installed. If not, either follow the install instructions below or just
# comment out the next five lines.
EXTRA_FLAGS= -DWP2_HAVE_PNG -DWP2_HAVE_JPEG -DWP2_HAVE_TIFF -DWP2_HAVE_GIF \
             -DWP2_HAVE_WEBP
IMG_LIBS= -lpng -lz -ljpeg -ltiff -lgif \
          -lwebp -lwebpdemux -Wno-unused-command-line-argument

# point to Googletest include dir and objects.
GTEST_DIR=../googletest
GTEST_FLAGS  = -I$(GTEST_DIR)/googletest/include
GTEST_FLAGS += -I$(GTEST_DIR)/googlemock/include
GTEST_LIBS  = $(GTEST_DIR)/googletest/make/gtest-all.o
GTEST_LIBS += $(GTEST_DIR)/googlemock/make/gtest_main.o

ifeq ($(strip $(shell uname)), Darwin)
  # Work around a problem linking tables marked as common symbols,
  # cf., src/enc/yuv.[hc]
  # Failure observed with: gcc 4.2.1 and 4.0.1.
  EXTRA_FLAGS += -fno-common
  EXTRA_FLAGS += -DHAVE_GLUT_GLUT_H -DWP2_HAVE_OPENGL
  EXTRA_FLAGS += -Wno-deprecated-declarations
  EXTRA_FLAGS += -I/usr/local/include
  EXTRA_LIBS  += -L/usr/local/lib
  GL_LIBS = -framework GLUT -framework OpenGL
  # for now, we force linking again MacOS 10.14 sdk, since
  # OpenGL is broken with 10.15.
  SYSROOT=-isysroot /Applications/Xcode.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX10.14.sdk
  CFLAGS += $(SYSROOT)
  CXXFLAGS += $(SYSROOT)
  LDFLAGS += $(SYSROOT)
else
  EXTRA_FLAGS += -I/usr/local/include -DWP2_HAVE_OPENGL
  EXTRA_LIBS  += -L/usr/local/lib
  GL_LIBS = -lglut -lGL
endif

# SDL flags: use sdl-config if it exists
SDL_CONFIG = $(shell sdl-config --version 2> /dev/null)
ifneq ($(SDL_CONFIG),)
  SDL_LIBS = $(shell sdl-config --libs)
  SDL_FLAGS = $(shell sdl-config --cflags)
else
  # use best-guess
  SDL_LIBS = -lSDL
  SDL_FLAGS =
endif

# To install libraries on Mac OS X:
# 1. Install MacPorts (http://www.macports.org/install.php)
# 2. Run "sudo port install jpeg"
# 3. Run "sudo port install libpng"
# 4. Run "sudo port install tiff"
# 5. Run "sudo port install giflib"

# To install libraries on Linux:
# 1. Run "sudo apt-get install libjpeg62-dev"
# 2. Run "sudo apt-get install libpng12-dev"
# 3. Run "sudo apt-get install libtiff4-dev"
# 4. Run "sudo apt-get install libgif-dev"

# Uncomment for build for 32bit platform
# Alternatively, you can just use the command
# 'make -f makefile.unix EXTRA_FLAGS=-m32' to that effect.
# EXTRA_FLAGS += -m32

# Extra flags to enable byte swap for 16 bit colorspaces.
# EXTRA_FLAGS += -DWP2_SWAP_16BIT_CSP

# Extra flags to enable multi-threading
EXTRA_FLAGS += -DWP2_USE_THREAD
EXTRA_LIBS += -lpthread

# Control symbol visibility. Comment out if your compiler doesn't support it.
EXTRA_FLAGS += -fvisibility=hidden

# Extra flags to emulate C89 strictness with the full ANSI
EXTRA_FLAGS += -Wmissing-declarations
# EXTRA_FLAGS += -Wdeclaration-after-statement
EXTRA_FLAGS += -Wshadow
EXTRA_FLAGS += -Wformat-security
# EXTRA_FLAGS += -Wvla

# SSE4.2-specific flags. If unsure about your compiler's defines, try:
#    $(CC) -dM -E - < /dev/null | grep SSE
ifeq ($(HAVE_SSE), 1)
EXTRA_FLAGS += -DWP2_HAVE_SSE
endif
EXTRA_FLAGS += -msse4

# AVX2-specific flags:
ifeq ($(HAVE_AVX2), 1)
EXTRA_FLAGS += -DWP2_HAVE_AVX2
src/dsp/%_avx2.o: EXTRA_FLAGS += -mavx2
endif


# NEON-specific flags:
# EXTRA_FLAGS += -march=armv7-a -mfloat-abi=hard -mfpu=neon -mtune=cortex-a8
# -> seems to make the overall lib slower: -fno-split-wide-types

# MIPS (MSA) 32-bit build specific flags for mips32r5 (p5600):
# EXTRA_FLAGS += -mips32r5 -mabi=32 -mtune=p5600 -mmsa -mfp64
# EXTRA_FLAGS += -msched-weight -mload-store-pairs

# MIPS (MSA) 64-bit build specific flags for mips64r6 (i6400):
# EXTRA_FLAGS += -mips64r6 -mabi=64 -mtune=i6400 -mmsa -mfp64
# EXTRA_FLAGS += -msched-weight -mload-store-pairs

# specific flags for C/C++
EXTRA_C_FLAGS = $(EXTRA_FLAGS)
EXTRA_C_FLAGS += -Wextra -Wold-style-definition
EXTRA_C_FLAGS += -Wmissing-prototypes

EXTRA_CXX_FLAGS = $(EXTRA_FLAGS)
EXTRA_CXX_FLAGS += -std=c++11

#### Nothing should normally be changed below this line ####

AR = ar
ARFLAGS = r
CPPFLAGS = -Isrc/ -I. -Wall
ifeq ($(DEBUG), 1)
  CFLAGS += -g $(EXTRA_C_FLAGS)
  CXXFLAGS += -g $(EXTRA_CXX_FLAGS)
else
  CFLAGS += -O3 -DNDEBUG $(EXTRA_C_FLAGS)
  CXXFLAGS += -O3 -DNDEBUG $(EXTRA_CXX_FLAGS)
endif

ifeq ($(TRACE), 1)
  CFLAGS += -DWP2_TRACE
  CXXFLAGS += -DWP2_TRACE
endif

ifeq ($(ERROR_TRACE), 1)
  CFLAGS += -DWP2_ERROR_TRACE
  CXXFLAGS += -DWP2_ERROR_TRACE
endif

ifeq ($(BITTRACE), 1)
  CFLAGS += -DWP2_BITTRACE
  CXXFLAGS += -DWP2_BITTRACE
endif

# for AOM_ACCOUNTING/INSPECTION, we need a specially built libaom, with the
# flags CONFIG_ACCOUNTING/CONFIG_INSPECTION defined and access to inner
# headers av1/decoder/...
ifneq ($(AOM_PATH),)
  CFLAGS += -I${AOM_PATH} -I${AOM_PATH}/build -DWP2_HAVE_AOM_DBG
  CXXFLAGS += -I${AOM_PATH} -I${AOM_PATH}/build -DWP2_HAVE_AOM_DBG
  EXTRA_LIBS += -L${AOM_PATH}/build
  # we'll need AOM support compiled too, now:
  AOM=1
endif

ifeq ($(AOM), 1)
  CFLAGS += -DWP2_HAVE_AOM
  CXXFLAGS += -DWP2_HAVE_AOM
  EXTRA_LIBS += -laom
endif

# AVIF support through libavif
ifeq ($(AVIF), 1)
  CFLAGS += -DWP2_HAVE_AVIF
  CXXFLAGS += -DWP2_HAVE_AVIF
  EXTRA_LIBS += -lavif
endif

ifeq ($(SJPEG), 1)
  CFLAGS += -DWP2_HAVE_SJPEG
  CXXFLAGS += -DWP2_HAVE_SJPEG
  EXTRA_LIBS += -lsjpeg
endif

ifeq ($(ENC_DEC_MATCH), 1)
  CFLAGS += -DWP2_ENC_DEC_MATCH
  CXXFLAGS += -DWP2_ENC_DEC_MATCH
endif

CC = gcc -std=gnu99
CXX = g++
INSTALL = install
GROFF = /usr/bin/groff
COL = /usr/bin/col
LDFLAGS += $(EXTRA_LIBS) $(EXTRA_FLAGS) -lm

ANIM_UTIL_OBJS =

DEC_OBJS = \
    src/dec/alpha_dec.o \
    src/dec/anim_dec.o \
    src/dec/bypass_dec.o \
    src/dec/features_dec.o \
    src/dec/filters/alpha_filter.o \
    src/dec/filters/block_map_filter.o \
    src/dec/filters/deblocking_filter.o \
    src/dec/filters/directional_filter.o \
    src/dec/filters/grain_filter.o \
    src/dec/filters/intertile_filter.o \
    src/dec/filters/intratile_filter.o \
    src/dec/filters/restoration_filter.o \
    src/dec/incr/decoder_api.o \
    src/dec/incr/decoder_context.o \
    src/dec/incr/decoder_info.o \
    src/dec/incr/decoder_skip.o \
    src/dec/incr/decoder_stages.o \
    src/dec/incr/decoder_state.o \
    src/dec/incr/decoder_tiles.o \
    src/dec/lossless/group4_dec.o \
    src/dec/lossless/losslessi_dec.o \
    src/dec/lossless_dec.o \
    src/dec/lossy_dec.o \
    src/dec/main_dec.o \
    src/dec/residuals_dec.o \
    src/dec/residuals_dec_aom.o \
    src/dec/symbols_dec.o \
    src/dec/syntax_dec.o \
    src/dec/tile_dec.o \

DSP_DEC_OBJS = \
    src/dsp/alpha.o \
    src/dsp/ans_dsp.o \
    src/dsp/argb_converter.o \
    src/dsp/cpu.o \
    src/dsp/cdef_dsp.o \
    src/dsp/csp_dsp.o \
    src/dsp/deblocking_filter_dsp.o \
    src/dsp/dec_dsp.o \
    src/dsp/intra_c.o \
    src/dsp/quant_c.o \
    src/dsp/quant_sse.o \
    src/dsp/transf_c.o \
    src/dsp/transf_sse.o \
    src/dsp/lossless/lossless.o \
    src/dsp/grain_c.o \
    src/dsp/math.o \
    src/dsp/ssim_dsp.o \
    src/dsp/wiener_dsp.o \

ENC_OBJS = \
    src/enc/alpha_enc.o \
    src/enc/analysis_enc.o \
    src/enc/anim/anim_enc.o \
    src/enc/anim/anim_frame.o \
    src/enc/anim/anim_rectangle.o \
    src/enc/bypass_enc.o \
    src/enc/config_enc.o \
    src/enc/distortion_enc.o \
    src/enc/lossless/analysisl_enc.o \
    src/enc/lossless/backward_references_cost_enc.o \
    src/enc/lossless/backward_references_enc.o \
    src/enc/lossless/group4_enc.o \
    src/enc/lossless/histogram_enc.o \
    src/enc/lossless/palette.o \
    src/enc/lossless/predictorl_enc.o \
    src/enc/lossless/losslessi_enc.o \
    src/enc/symbols_enc.o \
    src/enc/lossless_enc.o \
    src/enc/lossy_enc.o \
    src/enc/main_enc.o \
    src/enc/nearlossless_enc.o \
    src/enc/partition_enc.o \
    src/enc/partition_score_func.o \
    src/enc/partitioner.o \
    src/enc/predictor_enc.o \
    src/enc/residuals_enc.o \
    src/enc/residuals_enc_aom.o \
    src/enc/screen_content/screen_enc.o \
    src/enc/segment_enc.o \
    src/enc/syntax_enc.o \
    src/enc/writer_enc.o \

DSP_ENC_OBJS = \
    src/dsp/enc_dsp.o \
    src/dsp/lossless/lossless_dsp_enc.o \
    src/dsp/score_dsp.o \

EX_FORMAT_DEC_OBJS = \
    imageio/anim_image_dec.o \
    imageio/bmpdec.o \
    imageio/gifdec.o \
    imageio/gifdec_metadata.o \
    imageio/image_dec.o \
    imageio/jpegdec.o \
    imageio/pngdec.o \
    imageio/pnmdec.o \
    imageio/tiffdec.o \
    imageio/webpdec.o \
    imageio/webpenc.o \
    imageio/wp2dec.o \

EX_FORMAT_ENC_OBJS = \
    imageio/bmpenc.o \
    imageio/image_enc.o \
    imageio/file_format.o \
    imageio/pngenc.o \
    imageio/pnmenc.o \
    imageio/tiffenc.o \

EX_UTIL_OBJS = \

IMAGE_UTIL_OBJS = \
    imageio/imageio_util.o \

ANS_OBJS = \
    src/utils/ans.o \
    src/utils/ans_cost_table.o \
    src/utils/ans_utils.o \

PREVIEW_OBJS = \
    src/common/preview/preview.o \
    src/common/preview/preview_rasterizer.o \
    src/dec/preview/preview_dec.o \
    src/enc/preview/preview_analysis.o \
    src/enc/preview/preview_color.o \
    src/enc/preview/preview_config.o \
    src/enc/preview/preview_enc.o \
    src/enc/preview/preview_opt.o \

UTILS_DEC_OBJS = \
    src/utils/argb_buffer.o \
    src/utils/csp.o \
    src/utils/data.o \
    src/utils/context_switch.o \
    src/utils/data_source.o \
    src/utils/data_source_context.o \
    src/utils/data_source_stream.o \
    src/utils/front_mgr.o \
    src/utils/memory.o \
    src/utils/orientation.o \
    src/utils/plane.o \
    src/utils/quantizer.o \
    src/utils/status.o \
    src/utils/thread_utils.o \
    src/utils/utils.o \
    src/utils/wiener.o \

UTILS_ENC_OBJS = \
    src/common/integral.o \
    src/common/global_params.o \

COMMON_DEC_OBJS = \
    src/common/filters/rstr_flt_params.o \
    src/common/lossless/color_cache.o \
    src/common/lossy/block.o \
    src/common/lossy/block_size.o \
    src/common/lossy/block_size_io.o \
    src/common/lossy/context.o \
    src/common/lossy/predictor.o \
    src/common/lossy/quant_mtx.o \
    src/common/lossy/residuals.o \
    src/common/lossy/residuals_aom.o \
    src/common/lossy/rnd_mtx.o \
    src/common/lossy/segment.o \
    src/common/lossy/transforms.o \
    src/common/progress_watcher.o \
    src/common/symbols.o \
    src/common/vdebug.o \

EXAMPLE_UTILS_OBJS = \
    examples/example_utils.o \

EXTRAS_OBJS = \
    extras/aom_utils.o \
    extras/ccsp_imageio.o \
    extras/wp2_dec_12b.o \
    extras/extras.o \
    extras/y4menc.o \
    extras/y4mdec.o \

LIBWP2ANS_OBJS = $(ANS_OBJS)
LIBWP2PREVIEW_OBJS = $(PREVIEW_OBJS)
LIBWP2DECODER_OBJS = $(DEC_OBJS) $(DSP_DEC_OBJS) $(UTILS_DEC_OBJS) \
                     $(LIBWP2ANS_OBJS) $(LIBWP2PREVIEW_OBJS) $(COMMON_DEC_OBJS)
LIBWP2_OBJS = $(LIBWP2DECODER_OBJS) $(ENC_OBJS) $(DSP_ENC_OBJS) \
               $(UTILS_ENC_OBJS)

HDRS_INSTALLED = \
    src/wp2/base.h \
    src/wp2/debug.h \
    src/wp2/decode.h \
    src/wp2/encode.h \
    src/wp2/format_constants.h \

HDRS = \
    src/common/color_precision.h \
    src/common/filters/rstr_flt_params.h \
    src/common/lossless/color_cache.h \
    src/common/lossy/block.h \
    src/common/lossy/block_size.h \
    src/common/lossy/block_size_io.h \
    src/common/lossy/context.h \
    src/common/lossy/predictor.h \
    src/common/lossy/quant_mtx.h \
    src/common/lossy/segment.h \
    src/common/lossy/transforms.h \
    src/common/lossy/residuals.h \
    src/common/lossy/residuals_aom.h \
    src/common/lossy/rnd_mtx.h \
    src/common/lossy/aom/array_2d.h \
    src/common/lossy/aom/cdfs.inc \
    src/common/integral.h \
    src/common/global_params.h \
    src/common/preview/preview.h \
    src/common/progress_watcher.h \
    src/common/symbols.h \
    src/dsp/dsp.h \
    src/dsp/lossless/lossless.h \
    src/dsp/lossless/lossless_common.h \
    src/dsp/math.h \
    src/dec/filters/alpha_filter.h \
    src/dec/filters/block_map_filter.h \
    src/dec/filters/deblocking_filter.h \
    src/dec/filters/directional_filter.h \
    src/dec/filters/grain_filter.h \
    src/dec/filters/intertile_filter.h \
    src/dec/filters/intratile_filter.h \
    src/dec/filters/restoration_filter.h \
    src/dec/incr/decoder_context.h \
    src/dec/incr/decoder_info.h \
    src/dec/incr/decoder_skip.h \
    src/dec/incr/decoder_state.h \
    src/dec/lossless/losslessi_dec.h \
    src/dec/preview/preview_dec.h \
    src/dec/residuals_dec_aom.h \
    src/dec/symbols_dec.h \
    src/dec/tile_dec.h \
    src/dec/wp2_dec_i.h \
    src/enc/analysis.h \
    src/enc/anim/anim_enc.h \
    src/enc/lossless/losslessi_enc.h \
    src/enc/partition_score_func.h \
    src/enc/partitioner.h \
    src/enc/preview/preview_enc.h \
    src/enc/residuals_enc_aom.h \
    src/enc/screen_content/screen_enc.h \
    src/enc/symbols_enc.h \
    src/enc/wp2_enc_i.h \
    src/utils/ans.h \
    src/utils/ans_utils.h \
    src/utils/context_switch.h \
    src/utils/csp.h \
    src/utils/data_source.h \
    src/utils/data_source_context.h \
    src/utils/data_source_stream.h \
    src/utils/front_mgr.h \
    src/utils/orientation.h \
    src/utils/plane.h \
    src/utils/random.h \
    src/utils/stats.h \
    src/utils/thread_utils.h \
    src/utils/utils.h \
    src/utils/vector.h \
    src/utils/wiener.h \
    src/wp2/format_constants.h \
    $(HDRS_INSTALLED) \

OUT_LIBS =
OUT_LIBS += imageio/libimageio_util.a
OUT_LIBS += imageio/libimagedec.a
OUT_LIBS += imageio/libimageenc.a
OUT_LIBS += extras/libextras.a
OUT_LIBS += examples/example_utils.a
OUT_LIBS += src/libwebp2decoder.a
OUT_LIBS += src/libwebp2ans.a
OUT_LIBS += src/libwebp2preview.a
OUT_LIBS += src/libwebp2.a

OUT_EXAMPLES =
OUT_EXAMPLES += examples/cwp2
OUT_EXAMPLES += examples/dwp2
OUT_EXAMPLES += examples/vwp2

EXTRAS_LIB = extras/libextras.a

EXTRA_EXAMPLES =
EXTRA_EXAMPLES += examples/get_disto
EXTRA_EXAMPLES += examples/rd_curve
EXTRA_EXAMPLES += extras/vwp2_sdl
EXTRA_EXAMPLES += extras/av1enc
EXTRA_EXAMPLES += extras/mk_preview
EXTRA_EXAMPLES += extras/ewp2

# libs for testing. Let's just include everything we could need
BASE_LIBS =
BASE_LIBS += extras/libextras.a
BASE_LIBS += imageio/libimageenc.a
BASE_LIBS += imageio/libimagedec.a
BASE_LIBS += imageio/libimageio_util.a
BASE_LIBS += examples/example_utils.a
BASE_LIBS += src/libwebp2.a

TEST_LIBS =
TEST_LIBS += tests/include/helpers_incr.o
TEST_LIBS += tests/include/helpers_filter.o
TEST_LIBS += tests/include/helpers.o
TEST_LIBS += $(BASE_LIBS)
TEST_LIBS += $(GTEST_LIBS)
TEST_CXX_FLAGS =
TEST_CXX_FLAGS += $(GTEST_FLAGS) -lpthread -Wno-unused-command-line-argument

OTHER_EXAMPLES =

OUTPUT = $(OUT_LIBS) $(OUT_EXAMPLES)
ifeq ($(MAKECMDGOALS),clean)
  OUTPUT += $(EXTRA_EXAMPLES) $(OTHER_EXAMPLES)
  OUTPUT += $(EXTRAS_LIB)
endif

TEST_BINARIES =
TEST_BINARIES += tests/test_alloc
TEST_BINARIES += tests/test_analysis
TEST_BINARIES += tests/test_anim
TEST_BINARIES += tests/test_anim_incr
TEST_BINARIES += tests/test_anim_rectangle
TEST_BINARIES += tests/test_alpha_filter
TEST_BINARIES += tests/test_imageio
TEST_BINARIES += tests/test_imageio_conversion
TEST_BINARIES += tests/test_imageio_util
TEST_BINARIES += tests/test_imageio_anim
TEST_BINARIES += tests/test_example_utils
TEST_BINARIES += tests/test_ans
TEST_BINARIES += tests/test_block_size
TEST_BINARIES += tests/test_buffer
TEST_BINARIES += tests/test_dec
TEST_BINARIES += tests/test_diffusion
TEST_BINARIES += tests/test_intertile_filter
TEST_BINARIES += tests/test_deblocking_filter
TEST_BINARIES += tests/test_directional_filter
TEST_BINARIES += tests/test_restoration_filter
TEST_BINARIES += tests/test_extras
TEST_BINARIES += tests/test_header
TEST_BINARIES += tests/test_global_params
TEST_BINARIES += tests/test_math
TEST_BINARIES += tests/test_malloc_fail_at
TEST_BINARIES += tests/test_plane
TEST_BINARIES += tests/test_pessimization
TEST_BINARIES += tests/test_custom_csp
TEST_BINARIES += tests/test_custom_csp_dec
TEST_BINARIES += tests/test_custom_csp_enc
TEST_BINARIES += tests/test_custom_csp_enc_dec
TEST_BINARIES += tests/test_predictor
TEST_BINARIES += tests/test_partitions
TEST_BINARIES += tests/test_partition_score_func
TEST_BINARIES += tests/test_preview_api
TEST_BINARIES += tests/test_api
TEST_BINARIES += tests/test_optim
TEST_BINARIES += tests/test_score
TEST_BINARIES += tests/test_simple_enc_dec
TEST_BINARIES += tests/test_simple_encode
TEST_BINARIES += tests/test_slow_enc
TEST_BINARIES += tests/test_segment
TEST_BINARIES += tests/test_transf
TEST_BINARIES += tests/test_utils
TEST_BINARIES += tests/test_uv_mode
TEST_BINARIES += tests/test_aom_residuals
TEST_BINARIES += tests/test_residuals
TEST_BINARIES += tests/test_anim_frame
TEST_BINARIES += tests/test_incr_dec
TEST_BINARIES += tests/test_quant_mtx
TEST_BINARIES += tests/test_histogram
TEST_BINARIES += tests/test_group4
TEST_BINARIES += tests/test_palette
TEST_BINARIES += tests/test_color_cache
TEST_BINARIES += tests/test_color_precision
TEST_BINARIES += tests/test_coded_block
TEST_BINARIES += tests/test_context
TEST_BINARIES += tests/test_csp
TEST_BINARIES += tests/test_dsp
TEST_BINARIES += tests/test_dsp_speed
TEST_BINARIES += tests/test_data_source
TEST_BINARIES += tests/test_data_source_stream
TEST_BINARIES += tests/test_distortion
TEST_BINARIES += tests/test_orientation
TEST_BINARIES += tests/test_front_mgr
TEST_BINARIES += tests/test_formats
TEST_BINARIES += tests/test_speeds
TEST_BINARIES += tests/test_thread
TEST_BINARIES += tests/test_tile_shapes
TEST_BINARIES += tests/test_syntax_writer
TEST_BINARIES += tests/test_version
TEST_BINARIES += tests/test_progress
TEST_BINARIES += tests/ansz
TEST_BINARIES += tests/test_assert_bin
TEST_BINARIES += tests/test_symbols
TEST_BINARIES += tests/test_vector

ex: $(OUT_EXAMPLES)
all: ex $(EXTRA_EXAMPLES) $(OTHER_EXAMPLES) $(TEST_BINARIES)
extras: $(EXTRAS_LIB)

%.o: %.c $(HDRS)
	$(CC) $(CFLAGS) $(CPPFLAGS) -c $< -o $@

%.o: %.cc $(HDRS)
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -c $< -o $@

imageio/libimagedec.a: $(EX_FORMAT_DEC_OBJS)
imageio/libimageenc.a: $(EX_FORMAT_ENC_OBJS)
imageio/libimageio_util.a: $(IMAGE_UTIL_OBJS)
examples/example_utils.a: $(EXAMPLE_UTILS_OBJS)
extras/libextras.a: $(EXTRAS_OBJS)
src/libwebp2decoder.a: $(LIBWP2DECODER_OBJS)
src/libwebp2ans.a: $(LIBWP2ANS_OBJS)
src/libwebp2preview.a: $(LIBWP2PREVIEW_OBJS)
src/libwebp2.a: $(LIBWP2_OBJS)

%.a:
	$(AR) $(ARFLAGS) $@ $^

examples/cwp2: examples/cwp2.o
examples/cwp2: $(BASE_LIBS)
examples/cwp2: override EXTRA_LIBS += $(IMG_LIBS)

examples/dwp2: examples/dwp2.o
examples/dwp2: $(BASE_LIBS)
examples/dwp2: override EXTRA_LIBS += $(IMG_LIBS)

examples/vwp2: examples/vwp2.o
examples/vwp2: $(BASE_LIBS)
examples/vwp2: extras/libextras.a
examples/vwp2: override EXTRA_LIBS += $(IMG_LIBS) $(GL_LIBS)

examples/get_disto: examples/get_disto.o
examples/get_disto: $(BASE_LIBS)
examples/get_disto: override EXTRA_LIBS += $(IMG_LIBS)

examples/rd_curve: examples/rd_curve.o
examples/rd_curve: $(BASE_LIBS)
examples/rd_curve: extras/libextras.a
examples/rd_curve: override EXTRA_LIBS += $(IMG_LIBS)

extras/vwp2_sdl: extras/vwp2_sdl.o
extras/vwp2_sdl: extras/wp2_to_sdl.o
extras/vwp2_sdl: $(BASE_LIBS)
extras/vwp2_sdl: override EXTRA_FLAGS += -DWP2_HAVE_SDL $(SDL_FLAGS)
extras/vwp2_sdl: override EXTRA_LIBS += $(SDL_LIBS) $(IMG_LIBS)

$(OUT_EXAMPLES) $(EXTRA_EXAMPLES) $(OTHER_EXAMPLES):
	$(CXX) -o $@ $^ $(LDFLAGS)

tests/test_alloc: tests/test_alloc.o
tests/test_alloc: $(TEST_LIBS)
tests/test_alloc: override EXTRA_LIBS += $(IMG_LIBS)
tests/test_alloc: override EXTRA_CXX_FLAGS += $(TEST_CXX_FLAGS)

tests/test_analysis: tests/test_analysis.o
tests/test_analysis: $(TEST_LIBS)
tests/test_analysis: override EXTRA_LIBS += $(IMG_LIBS)
tests/test_analysis: override EXTRA_CXX_FLAGS += $(TEST_CXX_FLAGS)

tests/test_progress: tests/test_progress.o
tests/test_progress: $(TEST_LIBS)
tests/test_progress: override EXTRA_LIBS += $(IMG_LIBS)
tests/test_progress: override EXTRA_CXX_FLAGS += $(TEST_CXX_FLAGS)

tests/test_speeds: tests/test_speeds.o
tests/test_speeds: $(TEST_LIBS)
tests/test_speeds: override EXTRA_LIBS += $(IMG_LIBS)
tests/test_speeds: override EXTRA_CXX_FLAGS += $(TEST_CXX_FLAGS)

tests/test_thread: tests/test_thread.o
tests/test_thread: $(TEST_LIBS)
tests/test_thread: override EXTRA_LIBS += $(IMG_LIBS)
tests/test_thread: override EXTRA_CXX_FLAGS += $(TEST_CXX_FLAGS)

tests/test_tile_shapes: tests/test_tile_shapes.o
tests/test_tile_shapes: $(TEST_LIBS)
tests/test_tile_shapes: override EXTRA_LIBS += $(IMG_LIBS)
tests/test_tile_shapes: override EXTRA_CXX_FLAGS += $(TEST_CXX_FLAGS)

tests/test_version: tests/test_version.o
tests/test_version: $(TEST_LIBS)
tests/test_version: override EXTRA_LIBS += $(IMG_LIBS)
tests/test_version: override EXTRA_CXX_FLAGS += $(TEST_CXX_FLAGS)

tests/test_syntax_writer: tests/test_syntax_writer.o
tests/test_syntax_writer: $(TEST_LIBS)
tests/test_syntax_writer: override EXTRA_LIBS += $(IMG_LIBS)
tests/test_syntax_writer: override EXTRA_CXX_FLAGS += $(TEST_CXX_FLAGS)

tests/test_alpha_filter: tests/test_alpha_filter.o
tests/test_alpha_filter: $(TEST_LIBS)
tests/test_alpha_filter: override EXTRA_LIBS += $(IMG_LIBS)
tests/test_alpha_filter: override EXTRA_CXX_FLAGS += $(TEST_CXX_FLAGS)

tests/test_anim: tests/test_anim.o
tests/test_anim: $(TEST_LIBS)
tests/test_anim: override EXTRA_LIBS += $(IMG_LIBS)
tests/test_anim: override EXTRA_CXX_FLAGS += $(TEST_CXX_FLAGS)

tests/test_anim_incr: tests/test_anim_incr.o
tests/test_anim_incr: $(TEST_LIBS)
tests/test_anim_incr: override EXTRA_LIBS += $(IMG_LIBS)
tests/test_anim_incr: override EXTRA_CXX_FLAGS += $(TEST_CXX_FLAGS)

tests/test_anim_rectangle: tests/test_anim_rectangle.o
tests/test_anim_rectangle: $(TEST_LIBS)
tests/test_anim_rectangle: override EXTRA_LIBS += $(IMG_LIBS)
tests/test_anim_rectangle: override EXTRA_CXX_FLAGS += $(TEST_CXX_FLAGS)

tests/test_imageio: tests/test_imageio.o
tests/test_imageio: $(TEST_LIBS)
tests/test_imageio: override EXTRA_LIBS += $(IMG_LIBS)
tests/test_imageio: override EXTRA_CXX_FLAGS += $(TEST_CXX_FLAGS)

tests/test_imageio_conversion: tests/test_imageio_conversion.o
tests/test_imageio_conversion: $(TEST_LIBS)
tests/test_imageio_conversion: override EXTRA_LIBS += $(IMG_LIBS)
tests/test_imageio_conversion: override EXTRA_CXX_FLAGS += $(TEST_CXX_FLAGS)

tests/test_imageio_util: tests/test_imageio_util.o
tests/test_imageio_util: $(TEST_LIBS)
tests/test_imageio_util: override EXTRA_LIBS += $(IMG_LIBS)
tests/test_imageio_util: override EXTRA_CXX_FLAGS += $(TEST_CXX_FLAGS)

tests/test_imageio_anim: tests/test_imageio_anim.o
tests/test_imageio_anim: $(TEST_LIBS)
tests/test_imageio_anim: override EXTRA_LIBS += $(IMG_LIBS)
tests/test_imageio_anim: override EXTRA_CXX_FLAGS += $(TEST_CXX_FLAGS)

tests/test_example_utils: tests/test_example_utils.o
tests/test_example_utils: examples/example_utils.a
tests/test_example_utils: $(TEST_LIBS)
tests/test_example_utils: override EXTRA_LIBS += $(IMG_LIBS)
tests/test_example_utils: override EXTRA_CXX_FLAGS += $(TEST_CXX_FLAGS)

tests/test_header: tests/test_header.o
tests/test_header: $(TEST_LIBS)
tests/test_header: override EXTRA_LIBS += $(IMG_LIBS)
tests/test_header: override EXTRA_CXX_FLAGS += $(TEST_CXX_FLAGS)

tests/test_global_params: tests/test_global_params.o
tests/test_global_params: $(TEST_LIBS)
tests/test_global_params: override EXTRA_LIBS += $(IMG_LIBS)
tests/test_global_params: override EXTRA_CXX_FLAGS += $(TEST_CXX_FLAGS)

tests/test_transf: tests/test_transf.o
tests/test_transf: $(TEST_LIBS)
tests/test_transf: override EXTRA_LIBS += $(IMG_LIBS)
tests/test_transf: override EXTRA_CXX_FLAGS += $(TEST_CXX_FLAGS)

tests/test_utils: tests/test_utils.o
tests/test_utils: $(TEST_LIBS)
tests/test_utils: override EXTRA_LIBS += $(IMG_LIBS)
tests/test_utils: override EXTRA_CXX_FLAGS += $(TEST_CXX_FLAGS)

tests/test_uv_mode: tests/test_uv_mode.o
tests/test_uv_mode: $(TEST_LIBS)
tests/test_uv_mode: override EXTRA_LIBS += $(IMG_LIBS)
tests/test_uv_mode: override EXTRA_CXX_FLAGS += $(TEST_CXX_FLAGS)

tests/test_aom_residuals: tests/test_aom_residuals.o
tests/test_aom_residuals: $(TEST_LIBS)
tests/test_aom_residuals: override EXTRA_LIBS += $(IMG_LIBS)
tests/test_aom_residuals: override EXTRA_CXX_FLAGS += $(TEST_CXX_FLAGS)

tests/test_residuals: tests/test_residuals.o
tests/test_residuals: $(TEST_LIBS)
tests/test_residuals: override EXTRA_LIBS += $(IMG_LIBS)
tests/test_residuals: override EXTRA_CXX_FLAGS += $(TEST_CXX_FLAGS)

tests/test_diffusion: tests/test_diffusion.o
tests/test_diffusion: $(TEST_LIBS)
tests/test_diffusion: override EXTRA_LIBS += $(IMG_LIBS)
tests/test_diffusion: override EXTRA_CXX_FLAGS += $(TEST_CXX_FLAGS)

tests/test_intertile_filter: tests/test_intertile_filter.o
tests/test_intertile_filter: $(TEST_LIBS)
tests/test_intertile_filter: override EXTRA_LIBS += $(IMG_LIBS)
tests/test_intertile_filter: override EXTRA_CXX_FLAGS += $(TEST_CXX_FLAGS)

tests/test_deblocking_filter: tests/test_deblocking_filter.o
tests/test_deblocking_filter: $(TEST_LIBS)
tests/test_deblocking_filter: override EXTRA_LIBS += $(IMG_LIBS)
tests/test_deblocking_filter: override EXTRA_CXX_FLAGS += $(TEST_CXX_FLAGS)

tests/test_directional_filter: tests/test_directional_filter.o
tests/test_directional_filter: $(TEST_LIBS)
tests/test_directional_filter: override EXTRA_LIBS += $(IMG_LIBS)
tests/test_directional_filter: override EXTRA_CXX_FLAGS += $(TEST_CXX_FLAGS)

tests/test_restoration_filter: tests/test_restoration_filter.o
tests/test_restoration_filter: $(TEST_LIBS)
tests/test_restoration_filter: override EXTRA_LIBS += $(IMG_LIBS)
tests/test_restoration_filter: override EXTRA_CXX_FLAGS += $(TEST_CXX_FLAGS)

tests/test_quant_mtx: tests/test_quant_mtx.o
tests/test_quant_mtx: $(TEST_LIBS)
tests/test_quant_mtx: override EXTRA_LIBS += $(IMG_LIBS)
tests/test_quant_mtx: override EXTRA_CXX_FLAGS += $(TEST_CXX_FLAGS)

tests/test_color_cache: tests/lossless/test_color_cache.o
tests/test_color_cache: $(TEST_LIBS)
tests/test_color_cache: override EXTRA_LIBS += $(IMG_LIBS)
tests/test_color_cache: override EXTRA_CXX_FLAGS += $(TEST_CXX_FLAGS)

tests/test_group4: tests/lossless/test_group4.o
tests/test_group4: $(TEST_LIBS)
tests/test_group4: override EXTRA_LIBS += $(IMG_LIBS)
tests/test_group4: override EXTRA_CXX_FLAGS += $(TEST_CXX_FLAGS)

tests/test_palette: tests/lossless/test_palette.o
tests/test_palette: $(TEST_LIBS)
tests/test_palette: override EXTRA_LIBS += $(IMG_LIBS)
tests/test_palette: override EXTRA_CXX_FLAGS += $(TEST_CXX_FLAGS)

tests/test_histogram: tests/lossless/test_histogram.o
tests/test_histogram: $(TEST_LIBS)
tests/test_histogram: override EXTRA_LIBS += $(IMG_LIBS)
tests/test_histogram: override EXTRA_CXX_FLAGS += $(TEST_CXX_FLAGS)

tests/test_color_precision: tests/test_color_precision.o
tests/test_color_precision: $(TEST_LIBS)
tests/test_color_precision: override EXTRA_LIBS += $(IMG_LIBS)
tests/test_color_precision: override EXTRA_CXX_FLAGS += $(TEST_CXX_FLAGS)

tests/test_coded_block: tests/test_coded_block.o
tests/test_coded_block: $(TEST_LIBS)
tests/test_coded_block: override EXTRA_LIBS += $(IMG_LIBS)
tests/test_coded_block: override EXTRA_CXX_FLAGS += $(TEST_CXX_FLAGS)

tests/test_segment: tests/test_segment.o
tests/test_segment: $(TEST_LIBS)
tests/test_segment: override EXTRA_LIBS += $(IMG_LIBS)
tests/test_segment: override EXTRA_CXX_FLAGS += $(TEST_CXX_FLAGS)

tests/test_simple_enc_dec: tests/test_simple_enc_dec.o
tests/test_simple_enc_dec: $(BASE_LIBS)  # doesn't use gtest
tests/test_simple_enc_dec: override EXTRA_LIBS += $(IMG_LIBS)
tests/test_simple_enc_dec: override EXTRA_CXX_FLAGS += $(TEST_CXX_FLAGS)

tests/test_simple_encode: tests/test_simple_encode.o
tests/test_simple_encode: $(TEST_LIBS)
tests/test_simple_encode: override EXTRA_LIBS += $(IMG_LIBS)
tests/test_simple_encode: override EXTRA_CXX_FLAGS += $(TEST_CXX_FLAGS)

tests/test_slow_enc: tests/test_slow_enc.o
tests/test_slow_enc: $(TEST_LIBS)
tests/test_slow_enc: override EXTRA_LIBS += $(IMG_LIBS)
tests/test_slow_enc: override EXTRA_CXX_FLAGS += $(TEST_CXX_FLAGS)

tests/test_context: tests/test_context.o
tests/test_context: $(TEST_LIBS)
tests/test_context: override EXTRA_LIBS += $(IMG_LIBS)
tests/test_context: override EXTRA_CXX_FLAGS += $(TEST_CXX_FLAGS)

tests/test_csp: tests/test_csp.o
tests/test_csp: $(TEST_LIBS)
tests/test_csp: override EXTRA_LIBS += $(IMG_LIBS)
tests/test_csp: override EXTRA_CXX_FLAGS += $(TEST_CXX_FLAGS)

tests/test_dsp: tests/test_dsp.o
tests/test_dsp: $(TEST_LIBS)
tests/test_dsp: override EXTRA_LIBS += $(IMG_LIBS)
tests/test_dsp: override EXTRA_CXX_FLAGS += $(TEST_CXX_FLAGS)

tests/test_dsp_speed: tests/test_dsp_speed.o
tests/test_dsp_speed: $(TEST_LIBS)
tests/test_dsp_speed: override EXTRA_LIBS += $(IMG_LIBS)
tests/test_dsp_speed: override EXTRA_CXX_FLAGS += $(TEST_CXX_FLAGS)

tests/test_dec: tests/test_dec.o
tests/test_dec: $(TEST_LIBS)
tests/test_dec: override EXTRA_LIBS += $(IMG_LIBS)
tests/test_dec: override EXTRA_CXX_FLAGS += $(TEST_CXX_FLAGS)

tests/test_data_source: tests/test_data_source.o
tests/test_data_source: $(TEST_LIBS)
tests/test_data_source: override EXTRA_LIBS += $(IMG_LIBS)
tests/test_data_source: override EXTRA_CXX_FLAGS += $(TEST_CXX_FLAGS)

tests/test_data_source_stream: tests/test_data_source_stream.o
tests/test_data_source_stream: $(TEST_LIBS)
tests/test_data_source_stream: override EXTRA_LIBS += $(IMG_LIBS)
tests/test_data_source_stream: override EXTRA_CXX_FLAGS += $(TEST_CXX_FLAGS)

tests/test_distortion: tests/test_distortion.o
tests/test_distortion: $(TEST_LIBS)
tests/test_distortion: override EXTRA_LIBS += $(IMG_LIBS)
tests/test_distortion: override EXTRA_CXX_FLAGS += $(TEST_CXX_FLAGS)

tests/test_orientation: tests/test_orientation.o
tests/test_orientation: $(TEST_LIBS)
tests/test_orientation: override EXTRA_LIBS += $(IMG_LIBS)
tests/test_orientation: override EXTRA_CXX_FLAGS += $(TEST_CXX_FLAGS)

tests/test_anim_frame: tests/test_anim_frame.o
tests/test_anim_frame: $(TEST_LIBS)
tests/test_anim_frame: override EXTRA_LIBS += $(IMG_LIBS)
tests/test_anim_frame: override EXTRA_CXX_FLAGS += $(TEST_CXX_FLAGS)

tests/test_incr_dec: tests/test_incr_dec.o
tests/test_incr_dec: $(TEST_LIBS)
tests/test_incr_dec: override EXTRA_LIBS += $(IMG_LIBS)
tests/test_incr_dec: override EXTRA_CXX_FLAGS += $(TEST_CXX_FLAGS)

tests/test_malloc_fail_at: tests/test_malloc_fail_at.o
tests/test_malloc_fail_at: $(TEST_LIBS)
tests/test_malloc_fail_at: override EXTRA_LIBS += $(IMG_LIBS)
tests/test_malloc_fail_at: override EXTRA_CXX_FLAGS += $(TEST_CXX_FLAGS)

tests/test_math: tests/test_math.o
tests/test_math: $(TEST_LIBS)
tests/test_math: override EXTRA_LIBS += $(IMG_LIBS)
tests/test_math: override EXTRA_CXX_FLAGS += $(TEST_CXX_FLAGS)

tests/test_ans: tests/test_ans.o
tests/test_ans: $(TEST_LIBS)
tests/test_ans: override EXTRA_LIBS += $(IMG_LIBS)
tests/test_ans: override EXTRA_CXX_FLAGS += $(TEST_CXX_FLAGS)

tests/test_buffer: tests/test_buffer.o
tests/test_buffer: $(TEST_LIBS)
tests/test_buffer: override EXTRA_LIBS += $(IMG_LIBS)
tests/test_buffer: override EXTRA_CXX_FLAGS += $(TEST_CXX_FLAGS)

tests/test_block_size: tests/test_block_size.o
tests/test_block_size: $(TEST_LIBS)
tests/test_block_size: override EXTRA_LIBS += $(IMG_LIBS)
tests/test_block_size: override EXTRA_CXX_FLAGS += $(TEST_CXX_FLAGS)

tests/test_vector: tests/test_vector.o
tests/test_vector: $(TEST_LIBS)
tests/test_vector: override EXTRA_LIBS += $(IMG_LIBS)
tests/test_vector: override EXTRA_CXX_FLAGS += $(TEST_CXX_FLAGS)

tests/test_plane: tests/test_plane.o
tests/test_plane: $(TEST_LIBS)
tests/test_plane: override EXTRA_LIBS += $(IMG_LIBS)
tests/test_plane: override EXTRA_CXX_FLAGS += $(TEST_CXX_FLAGS)

tests/test_pessimization: tests/test_pessimization.o
tests/test_pessimization: $(TEST_LIBS)
tests/test_pessimization: override EXTRA_LIBS += $(IMG_LIBS)
tests/test_pessimization: override EXTRA_CXX_FLAGS += $(TEST_CXX_FLAGS)

tests/test_custom_csp: tests/test_custom_csp.o
tests/test_custom_csp: $(TEST_LIBS)
tests/test_custom_csp: override EXTRA_LIBS += $(IMG_LIBS)
tests/test_custom_csp: override EXTRA_CXX_FLAGS += $(TEST_CXX_FLAGS)

tests/test_custom_csp_dec: tests/test_custom_csp_dec.o
tests/test_custom_csp_dec: $(TEST_LIBS)
tests/test_custom_csp_dec: override EXTRA_LIBS += $(IMG_LIBS)
tests/test_custom_csp_dec: override EXTRA_CXX_FLAGS += $(TEST_CXX_FLAGS)

tests/test_custom_csp_enc: tests/test_custom_csp_enc.o
tests/test_custom_csp_enc: $(TEST_LIBS)
tests/test_custom_csp_enc: override EXTRA_LIBS += $(IMG_LIBS)
tests/test_custom_csp_enc: override EXTRA_CXX_FLAGS += $(TEST_CXX_FLAGS)

tests/test_custom_csp_enc_dec: tests/test_custom_csp_enc_dec.o
tests/test_custom_csp_enc_dec: $(TEST_LIBS)
tests/test_custom_csp_enc_dec: override EXTRA_LIBS += $(IMG_LIBS)
tests/test_custom_csp_enc_dec: override EXTRA_CXX_FLAGS += $(TEST_CXX_FLAGS)

tests/test_predictor: tests/test_predictor.o
tests/test_predictor: $(TEST_LIBS)
tests/test_predictor: override EXTRA_LIBS += $(IMG_LIBS)
tests/test_predictor: override EXTRA_CXX_FLAGS += $(TEST_CXX_FLAGS)

tests/test_partitions: tests/test_partitions.o
tests/test_partitions: $(TEST_LIBS)
tests/test_partitions: override EXTRA_LIBS += $(IMG_LIBS)
tests/test_partitions: override EXTRA_CXX_FLAGS += $(TEST_CXX_FLAGS)

tests/test_partition_score_func: tests/test_partition_score_func.o
tests/test_partition_score_func: $(TEST_LIBS)
tests/test_partition_score_func: override EXTRA_LIBS += $(IMG_LIBS)
tests/test_partition_score_func: override EXTRA_CXX_FLAGS += $(TEST_CXX_FLAGS)

tests/test_api: tests/test_api.o
tests/test_api: $(TEST_LIBS)
tests/test_api: override EXTRA_LIBS += $(IMG_LIBS)
tests/test_api: override EXTRA_CXX_FLAGS += $(TEST_CXX_FLAGS)

tests/test_preview_api: tests/test_preview_api.o
tests/test_preview_api: $(TEST_LIBS)
tests/test_preview_api: override EXTRA_LIBS += $(IMG_LIBS)
tests/test_preview_api: override EXTRA_CXX_FLAGS += $(TEST_CXX_FLAGS)

tests/test_optim: tests/test_optim.o
tests/test_optim: $(TEST_LIBS)
tests/test_optim: override EXTRA_LIBS += $(IMG_LIBS)
tests/test_optim: override EXTRA_CXX_FLAGS += $(TEST_CXX_FLAGS)

tests/test_front_mgr: tests/test_front_mgr.o
tests/test_front_mgr: $(TEST_LIBS)
tests/test_front_mgr: override EXTRA_LIBS += $(IMG_LIBS)
tests/test_front_mgr: override EXTRA_CXX_FLAGS += $(TEST_CXX_FLAGS)

tests/test_formats: tests/test_formats.o
tests/test_formats: $(TEST_LIBS)
tests/test_formats: override EXTRA_LIBS += $(IMG_LIBS)
tests/test_formats: override EXTRA_CXX_FLAGS += $(TEST_CXX_FLAGS)

tests/test_extras: tests/test_extras.o
tests/test_extras: $(TEST_LIBS)
tests/test_extras: override EXTRA_LIBS += $(IMG_LIBS)
tests/test_extras: override EXTRA_CXX_FLAGS += $(TEST_CXX_FLAGS)

tests/test_symbols: tests/test_symbols.o
tests/test_symbols: $(TEST_LIBS)
tests/test_symbols: override EXTRA_LIBS += $(IMG_LIBS)
tests/test_symbols: override EXTRA_CXX_FLAGS += $(TEST_CXX_FLAGS)

tests/test_score: tests/test_score.o
tests/test_score: $(TEST_LIBS)
tests/test_score: override EXTRA_LIBS += $(IMG_LIBS)
tests/test_score: override EXTRA_CXX_FLAGS += $(TEST_CXX_FLAGS)

tests/test_assert_bin: tests/test_assert.o
tests/test_assert_bin: src/libwebp2.a

tests/ansz: tests/ansz.o
tests/ansz: $(BASE_LIBS)
tests/ansz: override EXTRA_LIBS += $(IMG_LIBS)

extras/av1enc: extras/av1enc.o
extras/av1enc: $(BASE_LIBS)
extras/av1enc: override EXTRA_LIBS += $(IMG_LIBS) -lpthread

extras/mk_preview: extras/mk_preview.o
extras/mk_preview: $(BASE_LIBS)
extras/mk_preview: override EXTRA_LIBS += $(IMG_LIBS) -lpthread

extras/ewp2: extras/ewp2.o
extras/ewp2: $(BASE_LIBS)
extras/ewp2: override EXTRA_LIBS += $(IMG_LIBS) -lpthread


$(TEST_BINARIES):
	$(CXX) -o $@ $^ $(LDFLAGS)

dist: DESTDIR := dist
dist: OUT_EXAMPLES += $(EXTRA_EXAMPLES)
dist: all
	$(INSTALL) -m755 -d $(DESTDIR)/include/wp2 \
	           $(DESTDIR)/bin $(DESTDIR)/doc $(DESTDIR)/lib
	$(INSTALL) -m755 -s $(OUT_EXAMPLES) $(DESTDIR)/bin
	$(INSTALL) -m644 $(HDRS_INSTALLED) $(DESTDIR)/include/wp2
	$(INSTALL) -m644 src/libwebp2.a $(DESTDIR)/lib
	umask 022; \
	for m in man/[cd]wp2.1; do \
	  basenam=$$(basename $$m .1); \
	  $(GROFF) -t -e -man -T utf8 $$m \
	    | $(COL) -bx >$(DESTDIR)/doc/$${basenam}.txt; \
	  $(GROFF) -t -e -man -T html $$m \
	    | $(COL) -bx >$(DESTDIR)/doc/$${basenam}.html; \
	done

test: test_io test_cmd test_enc test_alloc test_ans test_extras test_ansz \
	test_assert test_math test_transf test_utils test_front_mgr test_plane \
	test_preview_api test_optim test_csp test_dsp test_vector test_buffer \
	test_symbols test_segment test_quant_mtx test_incr_dec test_api \
	test_coded_block test_distortion test_predictor test_score \
	test_anim test_anim_incr test_imageio_anim test_example_utils \
	test_orientation test_partitions test_uv_mode test_aom_residuals \
	test_residuals test_speeds test_progress test_intertile_filter \
	test_directional_filter test_deblocking_filter test_restoration_filter \
	test_header test_custom_csp test_custom_csp_dec test_malloc_fail_at \
	test_header test_custom_csp_enc test_custom_csp_enc_dec test_dec \
	test_thread test_pessimization test_formats test_global_params \
	test_data_source test_data_source_stream test_alpha_filter \
	test_anim_rectangle test_block_size test_context test_color_precision \
	test_imageio test_imageio_conversion test_imageio_util test_dsp_speed \
	test_simple_enc_dec test_simple_encode test_slow_enc test_io_pnm \
	test_color_cache test_group4 test_histogram test_analysis \
	test_partition_score_func test_tile_shapes test_version \
	test_syntax_writer test_anim_frame test_diffusion test_palette

test_io: $(OUT_EXAMPLES) $(EXTRA_EXAMPLES) $(OTHER_EXAMPLES)
	cd tests && ./test_io.sh ../examples
test_io_pnm: $(OUT_EXAMPLES) $(EXTRA_EXAMPLES) $(OTHER_EXAMPLES)
	cd tests && ./test_io_pnm.sh ../examples
test_cmd: $(OUT_EXAMPLES) $(EXTRA_EXAMPLES) $(OTHER_EXAMPLES)
	cd tests && ./test_cmd.sh ../examples > /dev/null
test_enc: $(OUT_EXAMPLES) $(EXTRA_EXAMPLES) $(OTHER_EXAMPLES)
	cd tests && ./test_enc.sh ../examples
test_enc_dir: $(OUT_EXAMPLES) $(EXTRA_EXAMPLES) $(OTHER_EXAMPLES)
	cd tests && ./test_enc_dir.sh ../examples
test_enc_anim: $(OUT_EXAMPLES) $(EXTRA_EXAMPLES) $(OTHER_EXAMPLES)
	cd tests && ./test_enc_anim.sh ../examples

test_alloc: tests/test_alloc
	cd tests && ./test_alloc
test_analysis: tests/test_analysis
	cd tests && ./test_analysis
test_speeds: tests/test_speeds
	cd tests && ./test_speeds
test_thread: tests/test_thread
	cd tests && ./test_thread
test_tile_shapes: tests/test_tile_shapes
	cd tests && ./test_tile_shapes
test_version: tests/test_version
	cd tests && ./test_version
test_syntax_writer: tests/test_syntax_writer
	cd tests && ./test_syntax_writer
test_progress: tests/test_progress
	cd tests && ./test_progress
test_alpha_filter: tests/test_alpha_filter
	cd tests && ./test_alpha_filter
test_anim: tests/test_anim
	cd tests && ./test_anim
test_anim_rectangle: tests/test_anim_rectangle
	cd tests && ./test_anim_rectangle
test_anim_incr: tests/test_anim_incr
	cd tests && ./test_anim_incr
test_imageio: tests/test_imageio
	cd tests && ./test_imageio
test_imageio_conversion: tests/test_imageio_conversion
	cd tests && ./test_imageio_conversion
test_imageio_util: tests/test_imageio_util
	cd tests && ./test_imageio_util
test_imageio_anim: tests/test_imageio_anim
	cd tests && ./test_imageio_anim
test_example_utils: tests/test_example_utils
	cd tests && ./test_example_utils
test_header: tests/test_header
	cd tests && ./test_header
test_global_params: tests/test_global_params
	cd tests && ./test_global_params
test_slow_enc: tests/test_slow_enc
	cd tests && ./test_slow_enc
test_simple_encode: tests/test_simple_encode
	cd tests && ./test_simple_encode
test_simple_enc_dec: tests/test_simple_enc_dec
	cd tests && ./test_simple_enc_dec
test_segment: tests/test_segment
	cd tests && ./test_segment
test_utils: tests/test_utils
	cd tests && ./test_utils
test_transf: tests/test_transf
	cd tests && ./test_transf
test_anim_frame: tests/test_anim_frame
	cd tests && ./test_anim_frame
test_incr_dec: tests/test_incr_dec
	cd tests && ./test_incr_dec
test_quant_mtx: tests/test_quant_mtx
	cd tests && ./test_quant_mtx
test_diffusion: tests/test_diffusion
	cd tests && ./test_diffusion
test_intertile_filter: tests/test_intertile_filter
	cd tests && ./test_intertile_filter
test_deblocking_filter: tests/test_deblocking_filter
	cd tests && ./test_deblocking_filter
test_directional_filter: tests/test_directional_filter
	cd tests && ./test_directional_filter
test_restoration_filter: tests/test_restoration_filter
	cd tests && ./test_restoration_filter
test_color_cache: tests/test_color_cache
	cd tests && ./test_color_cache
test_group4: tests/test_group4
	cd tests && ./test_group4
test_palette: tests/test_palette
	cd tests && ./test_palette
test_histogram: tests/test_histogram
	cd tests && ./test_histogram
test_color_precision: tests/test_color_precision
	cd tests && ./test_color_precision
test_coded_block: tests/test_coded_block
	cd tests && ./test_coded_block
test_context: tests/test_context
	cd tests && ./test_context
test_csp: tests/test_csp
	cd tests && ./test_csp
test_dsp: tests/test_dsp
	cd tests && ./test_dsp
test_dsp_speed: tests/test_dsp_speed
	cd tests && ./test_dsp_speed
test_dec: tests/test_dec
	cd tests && ./test_dec
test_data_source: tests/test_data_source
	cd tests && ./test_data_source
test_data_source_stream: tests/test_data_source_stream
	cd tests && ./test_data_source_stream
test_distortion: tests/test_distortion
	cd tests && ./test_distortion
test_orientation: tests/test_orientation
	cd tests && ./test_orientation
test_math: tests/test_math
	cd tests && ./test_math
test_malloc_fail_at: tests/test_malloc_fail_at
	cd tests && ./test_malloc_fail_at
test_buffer: tests/test_buffer
	cd tests && ./test_buffer
test_block_size: tests/test_block_size
	cd tests && ./test_block_size

test_ans: tests/test_ans
	cd tests && ./test_ans
test_ansz: tests/ansz
	cd tests && ./test_ansz.sh
test_extras: tests/test_extras
	cd tests && ./test_extras
test_front_mgr: tests/test_front_mgr
	cd tests && ./test_front_mgr
test_formats: tests/test_formats
	cd tests && ./test_formats
test_plane: tests/test_plane
	cd tests && ./test_plane
test_pessimization: tests/test_pessimization
	cd tests && ./test_pessimization
test_custom_csp: tests/test_custom_csp
	cd tests && ./test_custom_csp
test_custom_csp_dec: tests/test_custom_csp_dec
	cd tests && ./test_custom_csp_dec
test_custom_csp_enc: tests/test_custom_csp_enc
	cd tests && ./test_custom_csp_enc
test_custom_csp_enc_dec: tests/test_custom_csp_enc_dec
	cd tests && ./test_custom_csp_enc_dec
test_predictor: tests/test_predictor
	cd tests && ./test_predictor
test_partitions: tests/test_partitions
	cd tests && ./test_partitions
test_partition_score_func: tests/test_partition_score_func
	cd tests && ./test_partition_score_func
test_api: tests/test_api
	cd tests && ./test_api
test_preview_api: tests/test_preview_api
	cd tests && ./test_preview_api
test_optim: tests/test_optim
	cd tests && ./test_optim source0.ppm
test_score: tests/test_score
	cd tests && ./test_score
test_assert: tests/test_assert_bin
	cd tests && ./test_assert.sh
test_vector: tests/test_vector
	cd tests && ./test_vector
test_symbols: tests/test_symbols
	cd tests && ./test_symbols
test_uv_mode: tests/test_uv_mode
	cd tests && ./test_uv_mode
test_aom_residuals: tests/test_aom_residuals
	cd tests && ./test_aom_residuals
test_residuals: tests/test_residuals
	cd tests && ./test_residuals

clean:
	$(RM) $(OUTPUT) $(TEST_BINARIES) *~
	$(RM) -rf Testing
	for d in examples extras imageio src src/dec src/dsp src/enc src/wp2 \
	         src/utils tests man doc swig cmake src/enc/lossless \
	         src/dec/filters src/dec/incr src/dec/lossless src/dsp/lossless \
	         src/utils/lossless src/common/lossless src/common/lossy \
	         src/enc/preview src/dec/preview src/common/preview \
	         src/common/filters/ src/common src/enc/anim src/enc/screen_content \
	         tests/include; \
	do \
           $(RM) -f $$d/*.o $$d/*~ $$d/*.a; \
        done
	$(RM) -f wp2_js/*.bc wp2_js/*~

DIST_FILES= \
            API.txt \
            Makefile \
            README.md \
            README.wp2_js \
            COPYING \
            CMakeLists.txt \
            cmake/compiler.cmake \
            cmake/cpu.cmake \
            cmake/thread.cmake \
            cmake/wp2Config.cmake.in \
            doc/REMOVE.ME \
            examples/cwp2.cc \
            examples/dwp2.cc \
            examples/vwp2.cc \
            examples/get_disto.cc \
            examples/rd_curve.cc \
            examples/example_utils.cc \
            examples/example_utils.h \
            examples/stopwatch.h \
            extras/aom_utils.h \
            extras/aom_utils.cc \
            extras/av1enc.cc \
            extras/mk_preview.cc \
            extras/ccsp_imageio.h \
            extras/ccsp_imageio.cc \
            extras/extras.h \
            extras/ewp2.cc \
            extras/extras.cc \
            extras/vwp2_sdl.cc \
            extras/wp2_dec_12b.cc \
            extras/wp2_to_sdl.h \
            extras/wp2_to_sdl.cc \
            extras/y4mdec.cc \
            extras/y4menc.cc \
            imageio/anim_image_dec.cc \
            imageio/anim_image_dec.h \
            imageio/bmpdec.cc \
            imageio/bmpenc.cc \
            imageio/file_format.cc \
            imageio/file_format.h \
            imageio/gifdec.cc \
            imageio/gifdec_metadata.cc \
            imageio/image_dec.cc \
            imageio/image_dec.h \
            imageio/image_enc.cc \
            imageio/image_enc.h \
            imageio/imageio_util.cc \
            imageio/imageio_util.h \
            imageio/jpegdec.cc \
            imageio/pngdec.cc \
            imageio/pngenc.cc \
            imageio/pnmdec.cc \
            imageio/pnmenc.cc \
            imageio/tiffdec.cc \
            imageio/tiffenc.cc \
            imageio/webpdec.cc \
            imageio/webpenc.cc \
            imageio/wp2dec.cc \
            man/dwp2.1 \
            man/cwp2.1 \
            man/vwp2.1 \
            src/common/color_precision.h \
            src/common/constants.h \
            src/common/filters/rstr_flt_params.h \
            src/common/filters/rstr_flt_params.cc \
            src/common/integral.h \
            src/common/integral.cc \
            src/common/global_params.h \
            src/common/global_params.cc \
            src/common/lossless/color_cache.h \
            src/common/lossless/color_cache.cc \
            src/common/lossy/aom/array_2d.h \
            src/common/lossy/aom/cdfs.inc \
            src/common/lossy/block.h \
            src/common/lossy/block.cc \
            src/common/lossy/block_size.h \
            src/common/lossy/block_size.cc \
            src/common/lossy/block_size_io.h \
            src/common/lossy/block_size_io.cc \
            src/common/lossy/context.h \
            src/common/lossy/context.cc \
            src/common/lossy/predictor.h \
            src/common/lossy/predictor.cc \
            src/common/lossy/quant_mtx.h \
            src/common/lossy/quant_mtx.cc \
            src/common/lossy/residuals.h \
            src/common/lossy/residuals.cc \
            src/common/lossy/residuals_aom.h \
            src/common/lossy/residuals_aom.cc \
            src/common/lossy/rnd_mtx.h \
            src/common/lossy/rnd_mtx.cc \
            src/common/lossy/segment.h \
            src/common/lossy/segment.cc \
            src/common/lossy/transforms.h \
            src/common/lossy/transforms.cc \
            src/common/preview/preview.h \
            src/common/preview/preview.cc \
            src/common/preview/preview_rasterizer.cc \
            src/common/progress_watcher.h \
            src/common/progress_watcher.cc \
            src/common/symbols.h \
            src/common/symbols.cc \
            src/common/vdebug.cc \
            src/dec/alpha_dec.cc \
            src/dec/anim_dec.cc \
            src/dec/bypass_dec.cc \
            src/dec/features_dec.cc \
            src/dec/filters/alpha_filter.h \
            src/dec/filters/alpha_filter.cc \
            src/dec/filters/block_map_filter.h \
            src/dec/filters/block_map_filter.cc \
            src/dec/filters/deblocking_filter.h \
            src/dec/filters/deblocking_filter.cc \
            src/dec/filters/directional_filter.h \
            src/dec/filters/directional_filter.cc \
            src/dec/filters/grain_filter.h \
            src/dec/filters/grain_filter.cc \
            src/dec/filters/intertile_filter.h \
            src/dec/filters/intertile_filter.cc \
            src/dec/filters/intratile_filter.h \
            src/dec/filters/intratile_filter.cc \
            src/dec/filters/restoration_filter.h \
            src/dec/filters/restoration_filter.cc \
            src/dec/incr/decoder_api.cc \
            src/dec/incr/decoder_context.h \
            src/dec/incr/decoder_context.cc \
            src/dec/incr/decoder_info.h \
            src/dec/incr/decoder_info.cc \
            src/dec/incr/decoder_skip.h \
            src/dec/incr/decoder_skip.cc \
            src/dec/incr/decoder_stages.cc \
            src/dec/incr/decoder_state.h \
            src/dec/incr/decoder_state.cc \
            src/dec/incr/decoder_tiles.cc \
            src/dec/lossless_dec.cc \
            src/dec/lossy_dec.cc \
            src/dec/main_dec.cc \
            src/dec/preview/preview_dec.h \
            src/dec/preview/preview_dec.cc \
            src/dec/residuals_dec.cc \
            src/dec/residuals_dec_aom.h \
            src/dec/residuals_dec_aom.cc \
            src/dec/symbols_dec.h \
            src/dec/symbols_dec.cc \
            src/dec/syntax_dec.cc \
            src/dec/tile_dec.h \
            src/dec/tile_dec.cc \
            src/dec/wp2_dec_i.h \
            src/dsp/alpha.cc \
            src/dsp/ans_dsp.cc \
            src/dsp/argb_converter.cc \
            src/dsp/cpu.cc \
            src/dsp/cdef_dsp.cc \
            src/dsp/csp_dsp.cc \
            src/dsp/deblocking_filter_dsp.cc \
            src/dsp/dec_dsp.cc \
            src/dsp/dsp.h \
            src/dsp/enc_dsp.cc \
            src/dsp/grain_c.cc \
            src/dsp/math.cc \
            src/dsp/math.h \
            src/dsp/ssim_dsp.cc \
            src/dsp/intra_c.cc  \
            src/dsp/quant_c.cc  \
            src/dsp/quant_sse.cc  \
            src/dsp/transf_c.cc  \
            src/dsp/transf_sse.cc  \
            src/dsp/wiener_dsp.cc \
            src/enc/alpha_enc.cc \
            src/enc/analysis.h \
            src/enc/analysis_enc.cc \
            src/enc/anim/anim_enc.h \
            src/enc/anim/anim_enc.cc \
            src/enc/anim/anim_frame.cc \
            src/enc/anim/anim_rectangle.cc \
            src/enc/bypass_enc.cc \
            src/enc/config_enc.cc \
            src/enc/distortion_enc.cc \
            src/enc/lossless_enc.cc \
            src/enc/lossy_enc.cc \
            src/enc/main_enc.cc \
            src/enc/nearlossless_enc.cc \
            src/enc/partition_enc.cc \
            src/enc/partition_score_func.h \
            src/enc/partition_score_func.cc \
            src/enc/partitioner.h \
            src/enc/partitioner.cc \
            src/enc/predictor_enc.cc \
            src/enc/syntax_enc.cc \
            src/enc/preview/preview_analysis.cc \
            src/enc/preview/preview_color.cc \
            src/enc/preview/preview_config.cc \
            src/enc/preview/preview_enc.h \
            src/enc/preview/preview_enc.cc \
            src/enc/preview/preview_opt.cc \
            src/enc/residuals_enc.cc \
            src/enc/residuals_enc_aom.h \
            src/enc/residuals_enc_aom.cc \
            src/enc/screen_content/screen_enc.h \
            src/enc/screen_content/screen_enc.cc \
            src/enc/segment_enc.cc \
            src/enc/symbols_enc.h \
            src/enc/symbols_enc.cc \
            src/enc/wp2_enc_i.h \
            src/enc/writer_enc.cc \
            src/utils/ans.h \
            src/utils/ans.cc \
            src/utils/ans_cost_table.cc \
            src/utils/ans_utils.h \
            src/utils/ans_utils.cc \
            src/utils/context_switch.h \
            src/utils/context_switch.cc \
            src/utils/csp.h \
            src/utils/csp.cc \
            src/utils/data.cc \
            src/utils/data_source.h \
            src/utils/data_source.cc \
            src/utils/data_source_context.h \
            src/utils/data_source_context.cc \
            src/utils/data_source_stream.h \
            src/utils/data_source_stream.cc \
            src/utils/argb_buffer.cc \
            src/utils/front_mgr.h \
            src/utils/front_mgr.cc \
            src/utils/memory.cc \
            src/utils/orientation.cc \
            src/utils/orientation.h \
            src/utils/plane.cc \
            src/utils/plane.h \
            src/utils/quantizer.cc \
            src/utils/quantizer.h \
            src/utils/random.h \
            src/utils/stats.h \
            src/utils/status.cc \
            src/utils/thread_utils.cc \
            src/utils/thread_utils.h \
            src/utils/utils.cc \
            src/utils/utils.h \
            src/utils/vector.h \
            src/utils/wiener.h \
            src/utils/wiener.cc \
            src/wp2/base.h \
            src/wp2/debug.h \
            src/wp2/decode.h \
            src/wp2/encode.h \
            src/wp2/format_constants.h \
            src/dec/lossless/group4_dec.cc \
            src/dec/lossless/losslessi_dec.cc \
            src/dec/lossless/losslessi_dec.h \
            src/dsp/lossless/lossless.cc \
            src/dsp/lossless/lossless.h \
            src/dsp/lossless/lossless_common.h \
            src/dsp/lossless/lossless_dsp_enc.cc \
            src/enc/lossless/analysisl_enc.cc \
            src/enc/lossless/histogram_enc.h \
            src/enc/lossless/predictorl_enc.cc \
            src/enc/lossless/losslessi_enc.cc \
            src/enc/symbols_enc.h \
            src/enc/symbols_enc.cc \
            src/enc/lossless/losslessi_enc.h \
            src/enc/lossless/backward_references_cost_enc.cc \
            src/enc/lossless/backward_references_enc.cc \
            src/enc/lossless/backward_references_enc.h \
            src/enc/lossless/group4_enc.cc \
            src/enc/lossless/histogram_enc.cc \
            src/enc/lossless/palette.h \
            src/enc/lossless/palette.cc \
            swig/examples/pycwp2.py \
            swig/examples/pydwp2.py \
            swig/examples/pyvwp2.py \
            swig/CMakeLists.txt \
            swig/libwebp2.i \
            swig/README.md \
            tests/fuzz/common_tags.dict \
            tests/fuzz/common_wp2_tags.dict \
            tests/fuzz/fuzz_anim.cc \
            tests/fuzz/fuzz_anim_and_dec_config.cc \
            tests/fuzz/fuzz_anim_and_enc_config.cc \
            tests/fuzz/fuzz_anim_enc_config.cc \
            tests/fuzz/fuzz_image.cc \
            tests/fuzz/fuzz_image_and_dec_config.cc \
            tests/fuzz/fuzz_image_and_enc_config.cc \
            tests/fuzz/fuzz_image_enc_config.cc \
            tests/fuzz/fuzz_imageio_anim_dec.cc \
            tests/fuzz/fuzz_imageio_image_dec_enc.cc \
            tests/fuzz/fuzz_imageio_image_enc.cc \
            tests/fuzz/fuzz_utils.cc \
            tests/fuzz/fuzz_utils.h \
            tests/fuzz/repro_fuzz.cc \
            tests/include/helpers.h \
            tests/include/helpers.cc \
            tests/include/helpers_filter.h \
            tests/include/helpers_filter.cc \
            tests/include/helpers_incr.h \
            tests/include/helpers_incr.cc \
            tests/CMakeLists.txt \
            tests/ansz.cc \
            tests/test_alloc.cc \
            tests/test_alpha_filter.cc \
            tests/test_analysis.cc \
            tests/test_anim.cc \
            tests/test_anim_frame.cc \
            tests/test_anim_incr.cc \
            tests/test_anim_rectangle.cc \
            tests/test_ans.cc \
            tests/test_ansz.sh \
            tests/test_aom_residuals.cc \
            tests/test_api.cc \
            tests/test_assert.cc \
            tests/test_assert.sh \
            tests/test_buffer.cc \
            tests/test_block_size.cc \
            tests/test_cmd.sh \
            tests/lossless/test_color_cache.cc \
            tests/lossless/test_histogram.cc \
            tests/lossless/test_group4.cc \
            tests/lossless/test_palette.cc \
            tests/test_coded_block.cc \
            tests/test_color_precision.cc \
            tests/test_context.cc \
            tests/test_csp.cc \
            tests/test_custom_csp.cc \
            tests/test_custom_csp_dec.cc \
            tests/test_custom_csp_enc.cc \
            tests/test_custom_csp_enc_dec.cc \
            tests/test_dsp.cc \
            tests/test_dsp_speed.cc \
            tests/test_data_source.cc \
            tests/test_data_source_stream.cc \
            tests/test_deblocking_filter.cc \
            tests/test_dec.cc \
            tests/test_diffusion.cc \
            tests/test_directional_filter.cc \
            tests/test_distortion.cc \
            tests/test_enc.sh \
            tests/test_enc_anim.sh \
            tests/test_enc_dir.sh \
            tests/test_enc_ewp2.sh \
            tests/test_example_utils.cc \
            tests/test_extras.cc \
            tests/test_formats.cc \
            tests/test_front_mgr.cc \
            tests/test_header.cc \
            tests/test_global_params.cc \
            tests/test_imageio.cc \
            tests/test_imageio_anim.cc \
            tests/test_imageio_conversion.cc \
            tests/test_imageio_util.cc \
            tests/test_incr_dec.cc \
            tests/test_intertile_filter.cc \
            tests/test_io.sh \
            tests/test_io_pnm.sh \
            tests/test_io_ewp2.sh \
            tests/test_math.cc \
            tests/test_malloc_fail_at.cc \
            tests/test_optim.cc \
            tests/test_orientation.cc \
            tests/test_partitions.cc \
            tests/test_partition_score_func.cc \
            tests/test_pessimization.cc \
            tests/test_plane.cc \
            tests/test_pessimization.cc \
            tests/test_predictor.cc \
            tests/test_preview_api.cc \
            tests/test_progress.cc \
            tests/test_score.cc \
            tests/test_quant_mtx.cc \
            tests/test_restoration_filter.cc \
            tests/test_residuals.cc \
            tests/test_segment.cc \
            tests/test_simple_enc_dec.cc \
            tests/test_simple_encode.cc \
            tests/test_slow_enc.cc \
            tests/test_stuttering_server.sh \
            tests/test_speeds.cc \
            tests/test_symbols.cc \
            tests/test_swig_examples.sh \
            tests/test_swig_python.py \
            tests/test_syntax_writer.cc \
            tests/test_thread.cc \
            tests/test_tile_shapes.cc \
            tests/test_transf.cc \
            tests/test_utils.cc \
            tests/test_uv_mode.cc \
            tests/test_vector.cc \
            tests/test_version.cc \
            tests/testdata \
            tests/testdata/alpha_ramp.lossy.webp \
            tests/testdata/alpha_ramp.pam \
            tests/testdata/alpha_ramp.png \
            tests/testdata/alpha_ramp.tiff \
            tests/testdata/alpha_ramp.webp \
            tests/testdata/animation0.gif \
            tests/testdata/animation0.webp \
            tests/testdata/animation1.webp \
            tests/testdata/source0.pgm \
            tests/testdata/source0.ppm \
            tests/testdata/source1.itl.png \
            tests/testdata/source1.png \
            tests/testdata/source1_1x1.png \
            tests/testdata/source1_1x48.png \
            tests/testdata/source1_32x32.png \
            tests/testdata/source1_4x4.png \
            tests/testdata/source1_64x1.png \
            tests/testdata/source1_64x48.png \
            tests/testdata/source2.tiff \
            tests/testdata/source3.jpg \
            tests/testdata/source4.ll.webp \
            tests/testdata/source4.webp \
            tests/testdata/specific \
            tests/testdata/specific/bug449_1001.pam \
            tests/testdata/specific/bug449_1001.pgm \
            tests/testdata/specific/bug449_1001.ppm \
            tests/testdata/specific/bug449_1001_g.pam \
            tests/testdata/specific/bug449_255.pam \
            tests/testdata/specific/bug449_255.pgm \
            tests/testdata/specific/bug449_255.ppm \
            tests/testdata/specific/bug449_255_g.pam \
            tests/testdata/specific/bug449_37.pam \
            tests/testdata/specific/bug449_37.pgm \
            tests/testdata/specific/bug449_37.ppm \
            tests/testdata/specific/bug449_37_g.pam \
            tests/testdata/specific/bug449_65535.pam \
            tests/testdata/specific/bug449_65535.pgm \
            tests/testdata/specific/bug449_65535.ppm \
            tests/testdata/specific/bug449_65535_g.pam \
            tests/testdata/specific/bug449_ref.png \
            tests/testdata/specific/bug449_ref_gray.png \
            tests/testdata/specific/bug449_ref_gray_alpha.png \
            tests/testdata/specific/bug449_ref_noalpha.png \
            tests/tools/stuttering_http_server.py \
            wp2_js/index.html \
            wp2_js/index_wasm.html \
            wp2_js/test_wp2_js.wp2 \
            wp2_js/test_wp2_wasm.wp2 \

pak: clean
ifeq ($(strip $(shell uname)), Darwin)
	COPYFILE_DISABLE=1 tar -s ,.*,wp2-$(VERSION)\/~,p -czf $(ARCHIVE_FILE) \
  $(DIST_FILES)
else
	tar --transform="s#^#wp2-$(VERSION)/#" -czf $(ARCHIVE_FILE) $(DIST_FILES)
endif
	@echo "GENERATED ARCHIVE: $(ARCHIVE_FILE)"

.PHONY: all clean dist ex pak
.SUFFIXES:
