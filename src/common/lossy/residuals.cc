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
//   Common code for residuals
//
// Author: Vincent Rabaud (vrabaud@google.com)

#include "src/common/lossy/residuals.h"

#include <cstdarg>

#include "src/common/lossy/residuals_aom.h"
#include "src/common/lossy/aom/cdfs.inc"
#include "src/common/symbols.h"
#include "src/utils/utils.h"

namespace WP2 {

//------------------------------------------------------------------------------
// Zigzag orders (similar to AV1)

static const uint16_t kZigzag_2x2[] = { 0, 1, 2, 3 };
static const uint16_t kZigzag_2x4[] = { 0, 1, 2, 3, 4, 5, 6, 7 };
static const uint16_t kZigzag_4x2[] = { 0, 4, 1, 5, 2, 6, 3, 7 };
static const uint16_t kZigzag_2x8[] = {
  0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15
};
static const uint16_t kZigzag_8x2[] = {
  0, 8, 1, 9, 2, 10, 3, 11, 4, 12, 5, 13, 6, 14, 7, 15
};
static const uint16_t kZigzag_4x4[] = {
  0, 1, 4, 8, 5, 2, 3, 6, 9, 12, 13, 10, 7, 11, 14, 15
};
static const uint16_t kZigzag_4x8[] = {
  0,  1,  4,  2,  5,  8,  3,  6,  9,  12, 7,  10, 13, 16, 11, 14,
  17, 20, 15, 18, 21, 24, 19, 22, 25, 28, 23, 26, 29, 27, 30, 31,
};

static const uint16_t kZigzag_8x4[] = {
  0,  8, 1,  16, 9,  2, 24, 17, 10, 3, 25, 18, 11, 4,  26, 19,
  12, 5, 27, 20, 13, 6, 28, 21, 14, 7, 29, 22, 15, 30, 23, 31,
};

static const uint16_t kZigzag_4x16[] = {
  0,  1,  4,  2,  5,  8,  3,  6,  9,  12, 7,  10, 13, 16, 11, 14,
  17, 20, 15, 18, 21, 24, 19, 22, 25, 28, 23, 26, 29, 32, 27, 30,
  33, 36, 31, 34, 37, 40, 35, 38, 41, 44, 39, 42, 45, 48, 43, 46,
  49, 52, 47, 50, 53, 56, 51, 54, 57, 60, 55, 58, 61, 59, 62, 63,
};

static const uint16_t kZigzag_16x4[] = {
  0,  16, 1,  32, 17, 2,  48, 33, 18, 3,  49, 34, 19, 4,  50, 35,
  20, 5,  51, 36, 21, 6,  52, 37, 22, 7,  53, 38, 23, 8,  54, 39,
  24, 9,  55, 40, 25, 10, 56, 41, 26, 11, 57, 42, 27, 12, 58, 43,
  28, 13, 59, 44, 29, 14, 60, 45, 30, 15, 61, 46, 31, 62, 47, 63,
};

static const uint16_t kZigzag_8x32[] = {
  0,   1,   8,   2,   9,   16,  3,   10,  17,  24,  4,   11,  18,  25,  32,
  5,   12,  19,  26,  33,  40,  6,   13,  20,  27,  34,  41,  48,  7,   14,
  21,  28,  35,  42,  49,  56,  15,  22,  29,  36,  43,  50,  57,  64,  23,
  30,  37,  44,  51,  58,  65,  72,  31,  38,  45,  52,  59,  66,  73,  80,
  39,  46,  53,  60,  67,  74,  81,  88,  47,  54,  61,  68,  75,  82,  89,
  96,  55,  62,  69,  76,  83,  90,  97,  104, 63,  70,  77,  84,  91,  98,
  105, 112, 71,  78,  85,  92,  99,  106, 113, 120, 79,  86,  93,  100, 107,
  114, 121, 128, 87,  94,  101, 108, 115, 122, 129, 136, 95,  102, 109, 116,
  123, 130, 137, 144, 103, 110, 117, 124, 131, 138, 145, 152, 111, 118, 125,
  132, 139, 146, 153, 160, 119, 126, 133, 140, 147, 154, 161, 168, 127, 134,
  141, 148, 155, 162, 169, 176, 135, 142, 149, 156, 163, 170, 177, 184, 143,
  150, 157, 164, 171, 178, 185, 192, 151, 158, 165, 172, 179, 186, 193, 200,
  159, 166, 173, 180, 187, 194, 201, 208, 167, 174, 181, 188, 195, 202, 209,
  216, 175, 182, 189, 196, 203, 210, 217, 224, 183, 190, 197, 204, 211, 218,
  225, 232, 191, 198, 205, 212, 219, 226, 233, 240, 199, 206, 213, 220, 227,
  234, 241, 248, 207, 214, 221, 228, 235, 242, 249, 215, 222, 229, 236, 243,
  250, 223, 230, 237, 244, 251, 231, 238, 245, 252, 239, 246, 253, 247, 254,
  255,
};

static const uint16_t kZigzag_32x8[] = {
  0,   32,  1,   64,  33,  2,   96,  65,  34,  3,   128, 97,  66,  35,  4,
  160, 129, 98,  67,  36,  5,   192, 161, 130, 99,  68,  37,  6,   224, 193,
  162, 131, 100, 69,  38,  7,   225, 194, 163, 132, 101, 70,  39,  8,   226,
  195, 164, 133, 102, 71,  40,  9,   227, 196, 165, 134, 103, 72,  41,  10,
  228, 197, 166, 135, 104, 73,  42,  11,  229, 198, 167, 136, 105, 74,  43,
  12,  230, 199, 168, 137, 106, 75,  44,  13,  231, 200, 169, 138, 107, 76,
  45,  14,  232, 201, 170, 139, 108, 77,  46,  15,  233, 202, 171, 140, 109,
  78,  47,  16,  234, 203, 172, 141, 110, 79,  48,  17,  235, 204, 173, 142,
  111, 80,  49,  18,  236, 205, 174, 143, 112, 81,  50,  19,  237, 206, 175,
  144, 113, 82,  51,  20,  238, 207, 176, 145, 114, 83,  52,  21,  239, 208,
  177, 146, 115, 84,  53,  22,  240, 209, 178, 147, 116, 85,  54,  23,  241,
  210, 179, 148, 117, 86,  55,  24,  242, 211, 180, 149, 118, 87,  56,  25,
  243, 212, 181, 150, 119, 88,  57,  26,  244, 213, 182, 151, 120, 89,  58,
  27,  245, 214, 183, 152, 121, 90,  59,  28,  246, 215, 184, 153, 122, 91,
  60,  29,  247, 216, 185, 154, 123, 92,  61,  30,  248, 217, 186, 155, 124,
  93,  62,  31,  249, 218, 187, 156, 125, 94,  63,  250, 219, 188, 157, 126,
  95,  251, 220, 189, 158, 127, 252, 221, 190, 159, 253, 222, 191, 254, 223,
  255,
};

static const uint16_t kZigzag_8x8[] = {
  0,  1,  8,  16, 9,  2,  3,  10, 17, 24, 32, 25, 18, 11, 4,  5,
  12, 19, 26, 33, 40, 48, 41, 34, 27, 20, 13, 6,  7,  14, 21, 28,
  35, 42, 49, 56, 57, 50, 43, 36, 29, 22, 15, 23, 30, 37, 44, 51,
  58, 59, 52, 45, 38, 31, 39, 46, 53, 60, 61, 54, 47, 55, 62, 63
};

static const uint16_t kZigzag_8x16[] = {
  0,   1,   8,   2,   9,   16,  3,   10,  17,  24,  4,   11,  18,  25,  32,
  5,   12,  19,  26,  33,  40,  6,   13,  20,  27,  34,  41,  48,  7,   14,
  21,  28,  35,  42,  49,  56,  15,  22,  29,  36,  43,  50,  57,  64,  23,
  30,  37,  44,  51,  58,  65,  72,  31,  38,  45,  52,  59,  66,  73,  80,
  39,  46,  53,  60,  67,  74,  81,  88,  47,  54,  61,  68,  75,  82,  89,
  96,  55,  62,  69,  76,  83,  90,  97,  104, 63,  70,  77,  84,  91,  98,
  105, 112, 71,  78,  85,  92,  99,  106, 113, 120, 79,  86,  93,  100, 107,
  114, 121, 87,  94,  101, 108, 115, 122, 95,  102, 109, 116, 123, 103, 110,
  117, 124, 111, 118, 125, 119, 126, 127,
};

static const uint16_t kZigzag_16x8[] = {
  0,  16,  1,   32, 17,  2,   48,  33,  18, 3,  64,  49,  34,  19,  4,   80,
  65, 50,  35,  20, 5,   96,  81,  66,  51, 36, 21,  6,   112, 97,  82,  67,
  52, 37,  22,  7,  113, 98,  83,  68,  53, 38, 23,  8,   114, 99,  84,  69,
  54, 39,  24,  9,  115, 100, 85,  70,  55, 40, 25,  10,  116, 101, 86,  71,
  56, 41,  26,  11, 117, 102, 87,  72,  57, 42, 27,  12,  118, 103, 88,  73,
  58, 43,  28,  13, 119, 104, 89,  74,  59, 44, 29,  14,  120, 105, 90,  75,
  60, 45,  30,  15, 121, 106, 91,  76,  61, 46, 31,  122, 107, 92,  77,  62,
  47, 123, 108, 93, 78,  63,  124, 109, 94, 79, 125, 110, 95,  126, 111, 127,
};


static const uint16_t kZigzag_16x32[] = {
  0,   1,   16,  2,   17,  32,  3,   18,  33,  48,  4,   19,  34,  49,  64,
  5,   20,  35,  50,  65,  80,  6,   21,  36,  51,  66,  81,  96,  7,   22,
  37,  52,  67,  82,  97,  112, 8,   23,  38,  53,  68,  83,  98,  113, 128,
  9,   24,  39,  54,  69,  84,  99,  114, 129, 144, 10,  25,  40,  55,  70,
  85,  100, 115, 130, 145, 160, 11,  26,  41,  56,  71,  86,  101, 116, 131,
  146, 161, 176, 12,  27,  42,  57,  72,  87,  102, 117, 132, 147, 162, 177,
  192, 13,  28,  43,  58,  73,  88,  103, 118, 133, 148, 163, 178, 193, 208,
  14,  29,  44,  59,  74,  89,  104, 119, 134, 149, 164, 179, 194, 209, 224,
  15,  30,  45,  60,  75,  90,  105, 120, 135, 150, 165, 180, 195, 210, 225,
  240, 31,  46,  61,  76,  91,  106, 121, 136, 151, 166, 181, 196, 211, 226,
  241, 256, 47,  62,  77,  92,  107, 122, 137, 152, 167, 182, 197, 212, 227,
  242, 257, 272, 63,  78,  93,  108, 123, 138, 153, 168, 183, 198, 213, 228,
  243, 258, 273, 288, 79,  94,  109, 124, 139, 154, 169, 184, 199, 214, 229,
  244, 259, 274, 289, 304, 95,  110, 125, 140, 155, 170, 185, 200, 215, 230,
  245, 260, 275, 290, 305, 320, 111, 126, 141, 156, 171, 186, 201, 216, 231,
  246, 261, 276, 291, 306, 321, 336, 127, 142, 157, 172, 187, 202, 217, 232,
  247, 262, 277, 292, 307, 322, 337, 352, 143, 158, 173, 188, 203, 218, 233,
  248, 263, 278, 293, 308, 323, 338, 353, 368, 159, 174, 189, 204, 219, 234,
  249, 264, 279, 294, 309, 324, 339, 354, 369, 384, 175, 190, 205, 220, 235,
  250, 265, 280, 295, 310, 325, 340, 355, 370, 385, 400, 191, 206, 221, 236,
  251, 266, 281, 296, 311, 326, 341, 356, 371, 386, 401, 416, 207, 222, 237,
  252, 267, 282, 297, 312, 327, 342, 357, 372, 387, 402, 417, 432, 223, 238,
  253, 268, 283, 298, 313, 328, 343, 358, 373, 388, 403, 418, 433, 448, 239,
  254, 269, 284, 299, 314, 329, 344, 359, 374, 389, 404, 419, 434, 449, 464,
  255, 270, 285, 300, 315, 330, 345, 360, 375, 390, 405, 420, 435, 450, 465,
  480, 271, 286, 301, 316, 331, 346, 361, 376, 391, 406, 421, 436, 451, 466,
  481, 496, 287, 302, 317, 332, 347, 362, 377, 392, 407, 422, 437, 452, 467,
  482, 497, 303, 318, 333, 348, 363, 378, 393, 408, 423, 438, 453, 468, 483,
  498, 319, 334, 349, 364, 379, 394, 409, 424, 439, 454, 469, 484, 499, 335,
  350, 365, 380, 395, 410, 425, 440, 455, 470, 485, 500, 351, 366, 381, 396,
  411, 426, 441, 456, 471, 486, 501, 367, 382, 397, 412, 427, 442, 457, 472,
  487, 502, 383, 398, 413, 428, 443, 458, 473, 488, 503, 399, 414, 429, 444,
  459, 474, 489, 504, 415, 430, 445, 460, 475, 490, 505, 431, 446, 461, 476,
  491, 506, 447, 462, 477, 492, 507, 463, 478, 493, 508, 479, 494, 509, 495,
  510, 511,
};

static const uint16_t kZigzag_32x16[] = {
  0,   32,  1,   64,  33,  2,   96,  65,  34,  3,   128, 97,  66,  35,  4,
  160, 129, 98,  67,  36,  5,   192, 161, 130, 99,  68,  37,  6,   224, 193,
  162, 131, 100, 69,  38,  7,   256, 225, 194, 163, 132, 101, 70,  39,  8,
  288, 257, 226, 195, 164, 133, 102, 71,  40,  9,   320, 289, 258, 227, 196,
  165, 134, 103, 72,  41,  10,  352, 321, 290, 259, 228, 197, 166, 135, 104,
  73,  42,  11,  384, 353, 322, 291, 260, 229, 198, 167, 136, 105, 74,  43,
  12,  416, 385, 354, 323, 292, 261, 230, 199, 168, 137, 106, 75,  44,  13,
  448, 417, 386, 355, 324, 293, 262, 231, 200, 169, 138, 107, 76,  45,  14,
  480, 449, 418, 387, 356, 325, 294, 263, 232, 201, 170, 139, 108, 77,  46,
  15,  481, 450, 419, 388, 357, 326, 295, 264, 233, 202, 171, 140, 109, 78,
  47,  16,  482, 451, 420, 389, 358, 327, 296, 265, 234, 203, 172, 141, 110,
  79,  48,  17,  483, 452, 421, 390, 359, 328, 297, 266, 235, 204, 173, 142,
  111, 80,  49,  18,  484, 453, 422, 391, 360, 329, 298, 267, 236, 205, 174,
  143, 112, 81,  50,  19,  485, 454, 423, 392, 361, 330, 299, 268, 237, 206,
  175, 144, 113, 82,  51,  20,  486, 455, 424, 393, 362, 331, 300, 269, 238,
  207, 176, 145, 114, 83,  52,  21,  487, 456, 425, 394, 363, 332, 301, 270,
  239, 208, 177, 146, 115, 84,  53,  22,  488, 457, 426, 395, 364, 333, 302,
  271, 240, 209, 178, 147, 116, 85,  54,  23,  489, 458, 427, 396, 365, 334,
  303, 272, 241, 210, 179, 148, 117, 86,  55,  24,  490, 459, 428, 397, 366,
  335, 304, 273, 242, 211, 180, 149, 118, 87,  56,  25,  491, 460, 429, 398,
  367, 336, 305, 274, 243, 212, 181, 150, 119, 88,  57,  26,  492, 461, 430,
  399, 368, 337, 306, 275, 244, 213, 182, 151, 120, 89,  58,  27,  493, 462,
  431, 400, 369, 338, 307, 276, 245, 214, 183, 152, 121, 90,  59,  28,  494,
  463, 432, 401, 370, 339, 308, 277, 246, 215, 184, 153, 122, 91,  60,  29,
  495, 464, 433, 402, 371, 340, 309, 278, 247, 216, 185, 154, 123, 92,  61,
  30,  496, 465, 434, 403, 372, 341, 310, 279, 248, 217, 186, 155, 124, 93,
  62,  31,  497, 466, 435, 404, 373, 342, 311, 280, 249, 218, 187, 156, 125,
  94,  63,  498, 467, 436, 405, 374, 343, 312, 281, 250, 219, 188, 157, 126,
  95,  499, 468, 437, 406, 375, 344, 313, 282, 251, 220, 189, 158, 127, 500,
  469, 438, 407, 376, 345, 314, 283, 252, 221, 190, 159, 501, 470, 439, 408,
  377, 346, 315, 284, 253, 222, 191, 502, 471, 440, 409, 378, 347, 316, 285,
  254, 223, 503, 472, 441, 410, 379, 348, 317, 286, 255, 504, 473, 442, 411,
  380, 349, 318, 287, 505, 474, 443, 412, 381, 350, 319, 506, 475, 444, 413,
  382, 351, 507, 476, 445, 414, 383, 508, 477, 446, 415, 509, 478, 447, 510,
  479, 511,
};

static const uint16_t kZigzag_16x16[] = {
  0,   1,   16,  32,  17,  2,   3,   18,  33,  48,  64,  49,  34,  19,  4,
  5,   20,  35,  50,  65,  80,  96,  81,  66,  51,  36,  21,  6,   7,   22,
  37,  52,  67,  82,  97,  112, 128, 113, 98,  83,  68,  53,  38,  23,  8,
  9,   24,  39,  54,  69,  84,  99,  114, 129, 144, 160, 145, 130, 115, 100,
  85,  70,  55,  40,  25,  10,  11,  26,  41,  56,  71,  86,  101, 116, 131,
  146, 161, 176, 192, 177, 162, 147, 132, 117, 102, 87,  72,  57,  42,  27,
  12,  13,  28,  43,  58,  73,  88,  103, 118, 133, 148, 163, 178, 193, 208,
  224, 209, 194, 179, 164, 149, 134, 119, 104, 89,  74,  59,  44,  29,  14,
  15,  30,  45,  60,  75,  90,  105, 120, 135, 150, 165, 180, 195, 210, 225,
  240, 241, 226, 211, 196, 181, 166, 151, 136, 121, 106, 91,  76,  61,  46,
  31,  47,  62,  77,  92,  107, 122, 137, 152, 167, 182, 197, 212, 227, 242,
  243, 228, 213, 198, 183, 168, 153, 138, 123, 108, 93,  78,  63,  79,  94,
  109, 124, 139, 154, 169, 184, 199, 214, 229, 244, 245, 230, 215, 200, 185,
  170, 155, 140, 125, 110, 95,  111, 126, 141, 156, 171, 186, 201, 216, 231,
  246, 247, 232, 217, 202, 187, 172, 157, 142, 127, 143, 158, 173, 188, 203,
  218, 233, 248, 249, 234, 219, 204, 189, 174, 159, 175, 190, 205, 220, 235,
  250, 251, 236, 221, 206, 191, 207, 222, 237, 252, 253, 238, 223, 239, 254,
  255
};

static const uint16_t kZigzag_32x32[] = {
  0,    1,    32,   64,   33,   2,   3,    34,   65,   96,   128,  97,  66,
  35,   4,    5,    36,   67,   98,  129,  160,  192,  161,  130,  99,  68,
  37,   6,    7,    38,   69,   100, 131,  162,  193,  224,  256,  225, 194,
  163,  132,  101,  70,   39,   8,   9,    40,   71,   102,  133,  164, 195,
  226,  257,  288,  320,  289,  258, 227,  196,  165,  134,  103,  72,  41,
  10,   11,   42,   73,   104,  135, 166,  197,  228,  259,  290,  321, 352,
  384,  353,  322,  291,  260,  229, 198,  167,  136,  105,  74,   43,  12,
  13,   44,   75,   106,  137,  168, 199,  230,  261,  292,  323,  354, 385,
  416,  448,  417,  386,  355,  324, 293,  262,  231,  200,  169,  138, 107,
  76,   45,   14,   15,   46,   77,  108,  139,  170,  201,  232,  263, 294,
  325,  356,  387,  418,  449,  480, 512,  481,  450,  419,  388,  357, 326,
  295,  264,  233,  202,  171,  140, 109,  78,   47,   16,   17,   48,  79,
  110,  141,  172,  203,  234,  265, 296,  327,  358,  389,  420,  451, 482,
  513,  544,  576,  545,  514,  483, 452,  421,  390,  359,  328,  297, 266,
  235,  204,  173,  142,  111,  80,  49,   18,   19,   50,   81,   112, 143,
  174,  205,  236,  267,  298,  329, 360,  391,  422,  453,  484,  515, 546,
  577,  608,  640,  609,  578,  547, 516,  485,  454,  423,  392,  361, 330,
  299,  268,  237,  206,  175,  144, 113,  82,   51,   20,   21,   52,  83,
  114,  145,  176,  207,  238,  269, 300,  331,  362,  393,  424,  455, 486,
  517,  548,  579,  610,  641,  672, 704,  673,  642,  611,  580,  549, 518,
  487,  456,  425,  394,  363,  332, 301,  270,  239,  208,  177,  146, 115,
  84,   53,   22,   23,   54,   85,  116,  147,  178,  209,  240,  271, 302,
  333,  364,  395,  426,  457,  488, 519,  550,  581,  612,  643,  674, 705,
  736,  768,  737,  706,  675,  644, 613,  582,  551,  520,  489,  458, 427,
  396,  365,  334,  303,  272,  241, 210,  179,  148,  117,  86,   55,  24,
  25,   56,   87,   118,  149,  180, 211,  242,  273,  304,  335,  366, 397,
  428,  459,  490,  521,  552,  583, 614,  645,  676,  707,  738,  769, 800,
  832,  801,  770,  739,  708,  677, 646,  615,  584,  553,  522,  491, 460,
  429,  398,  367,  336,  305,  274, 243,  212,  181,  150,  119,  88,  57,
  26,   27,   58,   89,   120,  151, 182,  213,  244,  275,  306,  337, 368,
  399,  430,  461,  492,  523,  554, 585,  616,  647,  678,  709,  740, 771,
  802,  833,  864,  896,  865,  834, 803,  772,  741,  710,  679,  648, 617,
  586,  555,  524,  493,  462,  431, 400,  369,  338,  307,  276,  245, 214,
  183,  152,  121,  90,   59,   28,  29,   60,   91,   122,  153,  184, 215,
  246,  277,  308,  339,  370,  401, 432,  463,  494,  525,  556,  587, 618,
  649,  680,  711,  742,  773,  804, 835,  866,  897,  928,  960,  929, 898,
  867,  836,  805,  774,  743,  712, 681,  650,  619,  588,  557,  526, 495,
  464,  433,  402,  371,  340,  309, 278,  247,  216,  185,  154,  123, 92,
  61,   30,   31,   62,   93,   124, 155,  186,  217,  248,  279,  310, 341,
  372,  403,  434,  465,  496,  527, 558,  589,  620,  651,  682,  713, 744,
  775,  806,  837,  868,  899,  930, 961,  992,  993,  962,  931,  900, 869,
  838,  807,  776,  745,  714,  683, 652,  621,  590,  559,  528,  497, 466,
  435,  404,  373,  342,  311,  280, 249,  218,  187,  156,  125,  94,  63,
  95,   126,  157,  188,  219,  250, 281,  312,  343,  374,  405,  436, 467,
  498,  529,  560,  591,  622,  653, 684,  715,  746,  777,  808,  839, 870,
  901,  932,  963,  994,  995,  964, 933,  902,  871,  840,  809,  778, 747,
  716,  685,  654,  623,  592,  561, 530,  499,  468,  437,  406,  375, 344,
  313,  282,  251,  220,  189,  158, 127,  159,  190,  221,  252,  283, 314,
  345,  376,  407,  438,  469,  500, 531,  562,  593,  624,  655,  686, 717,
  748,  779,  810,  841,  872,  903, 934,  965,  996,  997,  966,  935, 904,
  873,  842,  811,  780,  749,  718, 687,  656,  625,  594,  563,  532, 501,
  470,  439,  408,  377,  346,  315, 284,  253,  222,  191,  223,  254, 285,
  316,  347,  378,  409,  440,  471, 502,  533,  564,  595,  626,  657, 688,
  719,  750,  781,  812,  843,  874, 905,  936,  967,  998,  999,  968, 937,
  906,  875,  844,  813,  782,  751, 720,  689,  658,  627,  596,  565, 534,
  503,  472,  441,  410,  379,  348, 317,  286,  255,  287,  318,  349, 380,
  411,  442,  473,  504,  535,  566, 597,  628,  659,  690,  721,  752, 783,
  814,  845,  876,  907,  938,  969, 1000, 1001, 970,  939,  908,  877, 846,
  815,  784,  753,  722,  691,  660, 629,  598,  567,  536,  505,  474, 443,
  412,  381,  350,  319,  351,  382, 413,  444,  475,  506,  537,  568, 599,
  630,  661,  692,  723,  754,  785, 816,  847,  878,  909,  940,  971, 1002,
  1003, 972,  941,  910,  879,  848, 817,  786,  755,  724,  693,  662, 631,
  600,  569,  538,  507,  476,  445, 414,  383,  415,  446,  477,  508, 539,
  570,  601,  632,  663,  694,  725, 756,  787,  818,  849,  880,  911, 942,
  973,  1004, 1005, 974,  943,  912, 881,  850,  819,  788,  757,  726, 695,
  664,  633,  602,  571,  540,  509, 478,  447,  479,  510,  541,  572, 603,
  634,  665,  696,  727,  758,  789, 820,  851,  882,  913,  944,  975, 1006,
  1007, 976,  945,  914,  883,  852, 821,  790,  759,  728,  697,  666, 635,
  604,  573,  542,  511,  543,  574, 605,  636,  667,  698,  729,  760, 791,
  822,  853,  884,  915,  946,  977, 1008, 1009, 978,  947,  916,  885, 854,
  823,  792,  761,  730,  699,  668, 637,  606,  575,  607,  638,  669, 700,
  731,  762,  793,  824,  855,  886, 917,  948,  979,  1010, 1011, 980, 949,
  918,  887,  856,  825,  794,  763, 732,  701,  670,  639,  671,  702, 733,
  764,  795,  826,  857,  888,  919, 950,  981,  1012, 1013, 982,  951, 920,
  889,  858,  827,  796,  765,  734, 703,  735,  766,  797,  828,  859, 890,
  921,  952,  983,  1014, 1015, 984, 953,  922,  891,  860,  829,  798, 767,
  799,  830,  861,  892,  923,  954, 985,  1016, 1017, 986,  955,  924, 893,
  862,  831,  863,  894,  925,  956, 987,  1018, 1019, 988,  957,  926, 895,
  927,  958,  989,  1020, 1021, 990, 959,  991,  1022, 1023
};

static const uint16_t* const kZigzags[] = {
  // same ordering as TrfSize
  kZigzag_2x2, kZigzag_4x2,  kZigzag_8x2,
  kZigzag_2x4, kZigzag_4x4,  kZigzag_8x4,  kZigzag_16x4,
  kZigzag_2x8, kZigzag_4x8,  kZigzag_8x8,  kZigzag_16x8,  kZigzag_32x8,
               kZigzag_4x16, kZigzag_8x16, kZigzag_16x16, kZigzag_32x16,
                             kZigzag_8x32, kZigzag_16x32, kZigzag_32x32,
  nullptr   // sentinel
};
STATIC_ASSERT_ARRAY_SIZE(kZigzags, TRF_LAST + 1);

//------------------------------------------------------------------------------

// clang-format off
// Array containing the symmetric dimension for a given dimension.
constexpr TrfSize kTDimSym[TRF_LAST] = {
    TRF_2x2, TRF_2x4,  TRF_2x8,                          // *_2
    TRF_4x2, TRF_4x4,  TRF_4x8,  TRF_4x16,               // *_4
    TRF_8x2, TRF_8x4,  TRF_8x8,  TRF_8x16,  TRF_8x32,    // *_8
             TRF_16x4, TRF_16x8, TRF_16x16, TRF_16x32,   // *_16
                       TRF_32x8, TRF_32x16, TRF_32x32};  // *_32

// Define a no_range for ease of read.
#define NO_R {255, 255}
uint8_t ResidualBoxAnalyzer::bound_range_x[TRF_LAST][2] = {
      NO_R,   NO_R,  NO_R,                      // *_2
    {0, 0}, {0, 2}, {0, 5},    NO_R,            // *_4
    {0, 0}, {0, 2}, {0, 6}, {0, 14}, {0, 30},   // *_8
            {0, 2}, {0, 6}, {0, 14}, {0, 30},   // *_16
                    {0, 6}, {0, 14}, {0, 30}};  // *_32
// clang-format on
uint8_t ResidualBoxAnalyzer::bound_range_per_x
    [TRF_LAST - ResidualBoxAnalyzer::kIgnoreNum][kMaxBlockSizePix][2] = {
        {NO_R, NO_R, {4, 4}, {2, 5}, {3, 6}, {2, 6}, {3, 6}, NO_R,
         NO_R, NO_R, NO_R,   NO_R,   NO_R,   NO_R,   NO_R,   NO_R,
         NO_R, NO_R, NO_R,   NO_R,   NO_R,   NO_R,   NO_R,   NO_R,
         NO_R, NO_R, NO_R,   NO_R,   NO_R,   NO_R,   NO_R,   NO_R},
        {NO_R,   NO_R,   {4, 5}, {4, 6}, {4, 6}, {4, 6}, {4, 6}, {4, 6},
         {4, 6}, {4, 6}, {4, 6}, {4, 6}, {4, 6}, {4, 6}, {5, 6}, NO_R,
         NO_R,   NO_R,   NO_R,   NO_R,   NO_R,   NO_R,   NO_R,   NO_R,
         NO_R,   NO_R,   NO_R,   NO_R,   NO_R,   NO_R,   NO_R,   NO_R},
        {NO_R,   NO_R,   {4, 5}, {4, 6}, {4, 6}, {4, 6}, {4, 6}, {4, 6},
         {4, 6}, {4, 6}, {4, 6}, {4, 6}, {4, 6}, {4, 6}, {4, 6}, {4, 6},
         {4, 6}, {4, 6}, {4, 6}, {4, 6}, {4, 6}, {4, 6}, {4, 6}, {4, 6},
         {4, 6}, {4, 6}, {4, 6}, {4, 6}, {4, 6}, {4, 6}, {6, 6}, NO_R},
        {NO_R, NO_R, NO_R, NO_R, NO_R, NO_R, NO_R, NO_R, NO_R, NO_R, NO_R,
         NO_R, NO_R, NO_R, NO_R, NO_R, NO_R, NO_R, NO_R, NO_R, NO_R, NO_R,
         NO_R, NO_R, NO_R, NO_R, NO_R, NO_R, NO_R, NO_R, NO_R, NO_R},
        {NO_R, NO_R, NO_R, NO_R, {2, 13}, {2, 14}, {3, 14}, NO_R,
         NO_R, NO_R, NO_R, NO_R, NO_R,    NO_R,    NO_R,    NO_R,
         NO_R, NO_R, NO_R, NO_R, NO_R,    NO_R,    NO_R,    NO_R,
         NO_R, NO_R, NO_R, NO_R, NO_R,    NO_R,    NO_R,    NO_R},
        {NO_R,    NO_R,    NO_R,    {5, 13}, {3, 14}, {4, 14}, {3, 14}, {4, 14},
         {3, 14}, {4, 14}, {3, 14}, {4, 14}, {3, 14}, {4, 14}, {5, 14}, NO_R,
         NO_R,    NO_R,    NO_R,    NO_R,    NO_R,    NO_R,    NO_R,    NO_R,
         NO_R,    NO_R,    NO_R,    NO_R,    NO_R,    NO_R,    NO_R,    NO_R},
        {NO_R,    NO_R,    NO_R,    {4, 13}, {4, 14}, {4, 14}, {4, 14}, {4, 14},
         {4, 14}, {4, 14}, {4, 14}, {4, 14}, {4, 14}, {4, 14}, {4, 14}, {4, 14},
         {4, 14}, {4, 14}, {4, 14}, {4, 14}, {4, 14}, {4, 14}, {4, 14}, {4, 14},
         {4, 14}, {4, 14}, {4, 14}, {4, 14}, {4, 14}, {4, 14}, {6, 14}, NO_R},
        {NO_R, NO_R, NO_R, NO_R, {2, 29}, {2, 29}, {3, 30}, NO_R,
         NO_R, NO_R, NO_R, NO_R, NO_R,    NO_R,    NO_R,    NO_R,
         NO_R, NO_R, NO_R, NO_R, NO_R,    NO_R,    NO_R,    NO_R,
         NO_R, NO_R, NO_R, NO_R, NO_R,    NO_R,    NO_R,    NO_R},
        {NO_R,    NO_R,    NO_R,    NO_R,    {3, 29}, {3, 29}, {3, 30}, {3, 30},
         {3, 30}, {3, 30}, {3, 30}, {3, 30}, {3, 30}, {3, 30}, {4, 30}, NO_R,
         NO_R,    NO_R,    NO_R,    NO_R,    NO_R,    NO_R,    NO_R,    NO_R,
         NO_R,    NO_R,    NO_R,    NO_R,    NO_R,    NO_R,    NO_R,    NO_R},
        {NO_R,    NO_R,    NO_R,    {5, 29}, {3, 29}, {4, 29}, {3, 30}, {4, 30},
         {3, 30}, {4, 30}, {3, 30}, {4, 30}, {3, 30}, {4, 30}, {3, 30}, {4, 30},
         {3, 30}, {4, 30}, {3, 30}, {4, 30}, {3, 30}, {4, 30}, {3, 30}, {4, 30},
         {3, 30}, {4, 30}, {3, 30}, {4, 30}, {3, 30}, {4, 30}, {5, 30}, NO_R}};
// TODO(vrabaud) the values in bound_range_per_y are very close (by one unit or
// less) to the ones in bound_range_per_x. A merge should be tried. Plus, we
// could try to have the same range for each fixed value to reduce the size.
uint8_t ResidualBoxAnalyzer::bound_range_per_y
    [TRF_LAST - ResidualBoxAnalyzer::kIgnoreNum][kMaxBlockSizePix][2] = {
        {NO_R, NO_R, {3, 5}, {3, 6}, {2, 6}, {3, 6}, {4, 6}, NO_R,
         NO_R, NO_R, NO_R,   NO_R,   NO_R,   NO_R,   NO_R,   NO_R,
         NO_R, NO_R, NO_R,   NO_R,   NO_R,   NO_R,   NO_R,   NO_R,
         NO_R, NO_R, NO_R,   NO_R,   NO_R,   NO_R,   NO_R,   NO_R},
        {NO_R, NO_R, NO_R, NO_R, {2, 13}, {2, 14}, {3, 14}, NO_R,
         NO_R, NO_R, NO_R, NO_R, NO_R,    NO_R,    NO_R,    NO_R,
         NO_R, NO_R, NO_R, NO_R, NO_R,    NO_R,    NO_R,    NO_R,
         NO_R, NO_R, NO_R, NO_R, NO_R,    NO_R,    NO_R,    NO_R},
        {NO_R, NO_R, NO_R, NO_R, {2, 29}, {2, 29}, {3, 30}, NO_R,
         NO_R, NO_R, NO_R, NO_R, NO_R,    NO_R,    NO_R,    NO_R,
         NO_R, NO_R, NO_R, NO_R, NO_R,    NO_R,    NO_R,    NO_R,
         NO_R, NO_R, NO_R, NO_R, NO_R,    NO_R,    NO_R,    NO_R},
        {NO_R, NO_R, NO_R, NO_R, NO_R, NO_R, NO_R, NO_R, NO_R, NO_R, NO_R,
         NO_R, NO_R, NO_R, NO_R, NO_R, NO_R, NO_R, NO_R, NO_R, NO_R, NO_R,
         NO_R, NO_R, NO_R, NO_R, NO_R, NO_R, NO_R, NO_R, NO_R, NO_R},
        {NO_R,   NO_R,   {4, 5}, {4, 6}, {4, 6}, {4, 6}, {4, 6}, {4, 6},
         {4, 6}, {4, 6}, {4, 6}, {4, 6}, {4, 6}, {4, 6}, {5, 6}, NO_R,
         NO_R,   NO_R,   NO_R,   NO_R,   NO_R,   NO_R,   NO_R,   NO_R,
         NO_R,   NO_R,   NO_R,   NO_R,   NO_R,   NO_R,   NO_R,   NO_R},
        {NO_R,    NO_R,    NO_R,    {4, 12}, {4, 13}, {3, 14}, {4, 14}, {3, 14},
         {4, 14}, {3, 14}, {4, 14}, {3, 14}, {4, 14}, {3, 14}, {4, 14}, NO_R,
         NO_R,    NO_R,    NO_R,    NO_R,    NO_R,    NO_R,    NO_R,    NO_R,
         NO_R,    NO_R,    NO_R,    NO_R,    NO_R,    NO_R,    NO_R,    NO_R},
        {NO_R,    NO_R,    NO_R,    NO_R,    {3, 29}, {3, 29}, {3, 30}, {3, 30},
         {3, 30}, {3, 30}, {3, 30}, {3, 30}, {3, 30}, {3, 30}, {4, 30}, NO_R,
         NO_R,    NO_R,    NO_R,    NO_R,    NO_R,    NO_R,    NO_R,    NO_R,
         NO_R,    NO_R,    NO_R,    NO_R,    NO_R,    NO_R,    NO_R,    NO_R},
        {NO_R,   NO_R,   {4, 5}, {4, 6}, {4, 6}, {4, 6}, {4, 6}, {4, 6},
         {4, 6}, {4, 6}, {4, 6}, {4, 6}, {4, 6}, {4, 6}, {4, 6}, {4, 6},
         {4, 6}, {4, 6}, {4, 6}, {4, 6}, {4, 6}, {4, 6}, {4, 6}, {4, 6},
         {4, 6}, {4, 6}, {4, 6}, {4, 6}, {4, 6}, {4, 6}, {6, 6}, NO_R},
        {NO_R,    NO_R,    NO_R,    {4, 13}, {4, 14}, {4, 14}, {4, 14}, {4, 14},
         {4, 14}, {4, 14}, {4, 14}, {4, 14}, {4, 14}, {4, 14}, {4, 14}, {4, 14},
         {4, 14}, {4, 14}, {4, 14}, {4, 14}, {4, 14}, {4, 14}, {4, 14}, {4, 14},
         {4, 14}, {4, 14}, {4, 14}, {4, 14}, {4, 14}, {4, 14}, {6, 14}, NO_R},
        {NO_R,    NO_R,    NO_R,    {4, 28}, {4, 29}, {3, 30}, {4, 30}, {3, 30},
         {4, 30}, {3, 30}, {4, 30}, {3, 30}, {4, 30}, {3, 30}, {4, 30}, {3, 30},
         {4, 30}, {3, 30}, {4, 30}, {3, 30}, {4, 30}, {3, 30}, {4, 30}, {3, 30},
         {4, 30}, {3, 30}, {4, 30}, {3, 30}, {4, 30}, {3, 30}, {6, 30}, NO_R}};
#undef NO_R

void ResidualBoxAnalyzer::FindBoundingBox(const int16_t* const coeffs,
                                          uint32_t size, TrfSize dim,
                                          uint32_t* const max_x,
                                          uint32_t* const max_y) {
  const uint32_t bw = TrfWidth[dim];
  const uint32_t shift = TrfLog2[bw];
  *max_x = 0, *max_y = 0;
  for (uint32_t i = 0; i < size; ++i) {
    if (coeffs[i] != 0) {
      // TODO(vrabaud) pre-compute or cache those.
      const uint32_t x = i & (bw - 1);
      const uint32_t y = i >> shift;
      if (x > *max_x) *max_x = x;
      if (y > *max_y) *max_y = y;
    }
  }
  assert(*max_x < bw && *max_y < TrfHeight[dim]);
}

void ResidualBoxAnalyzer::ShouldUseBounds(TrfSize tdim,
                                          uint32_t last_zigzag_ind,
                                          uint32_t max_x, uint32_t max_y,
                                          bool* const use_bounds_x,
                                          bool* const use_bounds_y) {
  uint32_t min_zig_zag_ind_x, min_zig_zag_ind_y;
  FindBounds(tdim, max_x, max_y, &min_zig_zag_ind_x, &min_zig_zag_ind_y);
  *use_bounds_x = (last_zigzag_ind >= min_zig_zag_ind_x);
  *use_bounds_y = (last_zigzag_ind >= min_zig_zag_ind_y);
}

void ResidualBoxAnalyzer::CanUseBounds(TrfSize tdim,
                                       bool* const can_use_bounds_x,
                                       bool* const can_use_bounds_y) {
  *can_use_bounds_x = (bound_range_x[tdim][0] != 255);
  *can_use_bounds_y = (bound_range_x[kTDimSym[tdim]][0] != 255);
}

void ResidualBoxAnalyzer::FindBounds(TrfSize tdim, uint32_t max_x,
                                     uint32_t max_y,
                                     uint32_t* const min_zig_zag_ind_x,
                                     uint32_t* const min_zig_zag_ind_y
                                     ) {
  const uint8_t bw = TrfWidth[tdim];
  const uint8_t bh = TrfHeight[tdim];
  ResidualIterator iter(tdim);
  uint32_t remaining_indices = (max_x + 1) * (max_y + 1);
  uint32_t num_zeros_x = 0u, num_zeros_y = 0u;

  const uint32_t bw_log = TrfLog2[bw], bh_log = TrfLog2[bh];
  bool can_use_bounds_x = false;
  bool can_use_bounds_y = false;
  *min_zig_zag_ind_x = *min_zig_zag_ind_y = bw * bh;
  for (uint16_t zigzag_ind = 0; remaining_indices > 0; ++iter, ++zigzag_ind) {
    const uint32_t x = iter.x(), y = iter.y();
    if (x > max_x) {
      ++num_zeros_x;
      if (!can_use_bounds_x && num_zeros_x >= bw_log) {
        can_use_bounds_x = true;
        *min_zig_zag_ind_x = zigzag_ind;
      }
    }
    if (y > max_y) {
      ++num_zeros_y;
      if (!can_use_bounds_y && num_zeros_y >= bh_log) {
        can_use_bounds_y = true;
        *min_zig_zag_ind_y = zigzag_ind;
      }
    }
    if (x <= max_x && y <= max_y) --remaining_indices;
  }
}

void ResidualBoxAnalyzer::GetRangeX(TrfSize tdim, uint8_t* const min,
                                    uint8_t* const max) {
  *min = bound_range_x[tdim][0];
  *max = bound_range_x[tdim][1];
}

void ResidualBoxAnalyzer::GetRangeY(TrfSize tdim, uint8_t* const min,
                                    uint8_t* const max) {
  const TrfSize tdim_sym = kTDimSym[tdim];
  *min = bound_range_x[tdim_sym][0];
  *max = bound_range_x[tdim_sym][1];
}

void ResidualBoxAnalyzer::GetRangePerX(TrfSize tdim, uint8_t max_x,
                                       uint8_t* const min, uint8_t* const max) {
  if ((uint32_t)tdim < kIgnoreNum) {
    // Done to reduce memory usage on bound_range_x_.
    *min = *max = 255;
  } else {
    *min = bound_range_per_x[tdim - kIgnoreNum][max_x][0];
    *max = bound_range_per_x[tdim - kIgnoreNum][max_x][1];
  }
}

void ResidualBoxAnalyzer::GetRangePerY(TrfSize tdim, uint8_t max_y,
                                       uint8_t* const min, uint8_t* const max) {
  if ((uint32_t)tdim < kIgnoreNum) {
    // Done to reduce memory usage on bound_range_y_.
    *min = *max = 255;
  } else {
    *min = bound_range_per_y[tdim - kIgnoreNum][max_y][0];
    *max = bound_range_per_y[tdim - kIgnoreNum][max_y][1];
  }
}

// Uncomment the code below if you want to generate the ResidualBoxAnalyzer
// tables.
#if 0
template <class T>
inline void UpdateRange(const T v, T* const range) {
  range[0] = std::min(v, range[0]);
  range[1] = std::max(v, range[1]);
}

void ResidualBoxAnalyzer::FillData() {
  // The code below fills the static arrays. Careful, it is not thread-safe.
  static bool is_data_filled = false;
  if (is_data_filled) return;
  // per transform, per x/y, min/max.
  static uint8_t bound_range_per_x_t[TRF_LAST][kMaxBlockSizePix][2];
  static uint8_t bound_range_per_y_t[TRF_LAST][kMaxBlockSizePix][2];

  for (uint32_t i = 0; i < (uint32_t)TRF_LAST; ++i) {
    const auto tdim = (TrfSize)i;
    // Initialize the x stats.
    const uint8_t bw = TrfWidth[tdim];
    const uint8_t bh = TrfHeight[tdim];
    for (uint8_t y = 0; y < bh; ++y) {
      bound_range_per_y_t[i][y][0] = bw - 1;
      bound_range_per_y_t[i][y][1] = 0;
    }
    for (uint8_t y = bh; y < kMaxBlockSizePix; ++y) {
      bound_range_per_y_t[i][y][0] = bound_range_per_y_t[i][y][1] = 255;
    }
    bound_range_x[i][0] = bw - 1;
    bound_range_x[i][1] = 0;
    // Initialize the y stats.
    for (uint8_t x = 0; x < bw; ++x) {
      bound_range_per_x_t[i][x][0] = bh - 1;
      bound_range_per_x_t[i][x][1] = 0;
    }
    for (uint8_t x = bw; x < kMaxBlockSizePix; ++x) {
      bound_range_per_x_t[i][x][0] = bound_range_per_x_t[i][x][1] = 255;
    }
    bound_range_y[i][0] = bh - 1;
    bound_range_y[i][1] = 0;

    // Iterate over the block.
    for (uint8_t max_y = 0; max_y < bh; ++max_y) {
      for (uint8_t max_x = 0; max_x < bw; ++max_x) {
        uint32_t min_zig_zag_ind_x, min_zig_zag_ind_y;
        bool can_use_bounds_x, can_use_bounds_y;
        FindBounds(tdim, max_x, max_y, &min_zig_zag_ind_x, &min_zig_zag_ind_y,
                   &can_use_bounds_x, &can_use_bounds_y);

        // Update the range of the bounds.
        if (can_use_bounds_x) UpdateRange(max_x, bound_range_x[i]);
        if (can_use_bounds_y) UpdateRange(max_y, bound_range_y[i]);
        // Update the range of the bounds per x/y.
        if (can_use_bounds_x && can_use_bounds_y) {
          UpdateRange(max_y, bound_range_per_x_t[i][max_x]);
          UpdateRange(max_x, bound_range_per_y_t[i][max_y]);
        }
      }
    }
    for (uint8_t x = 0; x < bw; ++x) {
      if (bound_range_per_x_t[i][x][0] > bound_range_per_x_t[i][x][1]) {
        bound_range_per_x_t[i][x][0] = bound_range_per_x_t[i][x][1] = 255;
      }
    }
    for (uint8_t y = 0; y < bh; ++y) {
      if (bound_range_per_y_t[i][y][0] > bound_range_per_y_t[i][y][1]) {
        bound_range_per_y_t[i][y][0] = bound_range_per_y_t[i][y][1] = 255;
      }
    }
    if (bound_range_x[i][0] > bound_range_x[i][1]) {
      bound_range_x[i][0] = bound_range_x[i][1] = 255;
    }
    if (bound_range_y[i][0] > bound_range_y[i][1]) {
      bound_range_y[i][0] = bound_range_y[i][1] = 255;
    }
  }
  for (uint32_t i = 0; i < TRF_LAST; ++i) {
    for (uint32_t j = 0; j < kMaxBlockSizePix; ++j) {
      for (uint32_t k = 0; k < 2; ++k) {
        if (i < kIgnoreNum) {
          assert(bound_range_per_x_t[i][j][k] == 255);
          assert(bound_range_per_y_t[i][j][k] == 255);
        } else {
          bound_range_per_x[i - kIgnoreNum][j][k] =
              bound_range_per_x_t[i][j][k];
          bound_range_per_y[i - kIgnoreNum][j][k] =
              bound_range_per_y_t[i][j][k];
        }
      }
    }
  }
  is_data_filled = true;
}
#endif

//------------------------------------------------------------------------------

ResidualIterator::ResidualIterator(TrfSize dim)
    : max_i_(kNumCoeffs[dim]),
      zigzag_(kZigzags[dim]),
      dim_(dim),
      bw_(TrfWidth[dim_]),
      bw_shift_(TrfLog2[bw_]) {
  assert(((bw_ - 1) & bw_) == 0);  // is a power of 2
  Reset();
}

void ResidualIterator::Reset() {
  ind_ = 0u;
  x_ = y_ = 0u;
}

void ResidualIterator::operator++() {
  assert(ind_ < max_i_);
  ++ind_;
  if (ind_ >= max_i_) return;
  x_ = zigzag_[ind_] & (bw_ - 1);
  y_ = zigzag_[ind_] >> bw_shift_;
}

const uint16_t* ResidualIterator::GetZigzag(TrfSize dim) {
  return kZigzags[dim];
}

BoundedResidualIterator::BoundedResidualIterator(TrfSize dim, bool use_bounds_x,
                                                 bool use_bounds_y,
                                                 uint32_t max_x, uint32_t max_y)
    : ResidualIterator(dim),
      max_x_(max_x),
      max_y_(max_y),
      use_bounds_x_(use_bounds_x),
      use_bounds_y_(use_bounds_y) {
  Reset();
}

void BoundedResidualIterator::operator++() {
  assert(max_num_coeffs_left_ > 0);
  --max_num_coeffs_left_;
  do {
    ++ind_;
    if (ind_ >= max_i_) return;
    x_ = zigzag_[ind_] & (bw_ - 1);
    y_ = zigzag_[ind_] >> bw_shift_;
    // Get to the next coefficients if we are out of bounds.
  } while (x_ > max_x_ || y_ > max_y_);
}

void BoundedResidualIterator::Reset() {
  ResidualIterator::Reset();
  can_eob_max_x_ = !use_bounds_x_;
  can_eob_max_y_ = !use_bounds_y_;
  max_num_coeffs_left_ = (max_x_ + 1) * (max_y_ + 1);
}

void BoundedResidualIterator::SetAsNonZero() {
  // We must have reached both box boundaries before being allowed to write an
  // End Of Block. So remember if we have done so.
  if (x_ == max_x_) can_eob_max_x_ = true;
  if (y_ == max_y_) can_eob_max_y_ = true;
}

bool BoundedResidualIterator::CanEOB() const {
  return (can_eob_max_x_ && can_eob_max_y_);
}
//------------------------------------------------------------------------------

void ResidualIO::Init(bool use_aom_coeffs) {
  use_aom_coeffs_ = use_aom_coeffs;
}

uint32_t ResidualIO::GetSector(uint32_t x, uint32_t y, TrfSize dim) {
  uint32_t quadrant = 0;
  if (y >= TrfHeight[dim] / 2) quadrant += 2;
  if (x >= TrfWidth[dim] / 2) quadrant += 1;
  return quadrant;
}

static uint32_t DoGetCluster(uint32_t channel, uint32_t num_channels,
                             EncodingMethod method = EncodingMethod::kMethod0,
                             TrfSize dim = TRF_2x2, uint32_t sector = 0) {
  assert((uint32_t)method < kNumResidualStats);
  assert(channel < num_channels);

  uint32_t cluster = sector;
  cluster = cluster * TRF_LAST + dim;
  cluster = cluster * kNumResidualStats + ResidualIO::GetMethodIndex(method);
  cluster = cluster * num_channels + channel;
  return cluster;
}

uint32_t ResidualIO::GetCluster(Channel channel, uint32_t num_channels,
                                EncodingMethod method, TrfSize dim,
                                uint32_t sector) {
  return DoGetCluster(channel, num_channels, method, dim, sector);
}

uint32_t ResidualIO::GetClusterMergedUV(Channel channel, uint32_t num_channels,
                                        EncodingMethod method, TrfSize dim,
                                        uint32_t sector) {
  const uint32_t c = (channel == kYChannel)                           ? 0
                     : (channel == kUChannel || channel == kVChannel) ? 1
                                                                      : 2;
  if (num_channels >= 3) num_channels -= 1;
  return DoGetCluster(c, num_channels, method, dim, sector);
}

WP2Status InitializeAOMSymbol(Symbol sym, uint32_t quantizer_context,
                              SymbolsInfo* const info);

WP2Status ResidualIO::InitializeInfo(bool has_alpha, uint32_t quality_hint,
                                     bool use_aom_coeffs,
                                     SymbolsInfo* const info) {
  const uint32_t num_channels = (has_alpha ? 4 : 3);
  for (Channel channel : {kYChannel, kUChannel, kAChannel}) {
    if (channel == kAChannel && !has_alpha) continue;
    // For method 0, favor sparse data with low probabilities for 0.
    // For method 1, favor dense data.
    for (const TrfSize dim : kAllTrfSizes) {
      for (uint32_t sector = 0; sector < ResidualIO::kNumSectors; ++sector) {
        WP2_CHECK_STATUS(info->SetStartingProba(
            ResidualIO::GetClusterMergedUV(
                channel, num_channels, EncodingMethod::kMethod0, dim, sector),
            kSymbolResidualIsZero, 2, 8));
        WP2_CHECK_STATUS(info->SetStartingProba(
            ResidualIO::GetClusterMergedUV(
                channel, num_channels, EncodingMethod::kMethod1, dim, sector),
            kSymbolResidualIsZero, 8, 2));
      }
    }
  }
  assert(quality_hint <= kMaxLossyQualityHint);
  const int quantizer_context =
      3 - DivRound(quality_hint * 3u, kMaxLossyQualityHint);

  for (Symbol sym : kSymbolsForAOMCoeffs) {
    if (use_aom_coeffs) {
      WP2_CHECK_STATUS(InitializeAOMSymbol(sym, quantizer_context, info));
    } else {
      info->SetInfo(sym, /*min=*/SymbolsInfo::kInvalidBound,
                    /*max=*/SymbolsInfo::kInvalidBound, /*num_clusters=*/0,
                    SymbolsInfo::StorageMethod::kUnused);
    }
  }
  return WP2_STATUS_OK;
}

void ResidualIO::SetGeometry(uint32_t num_coeffs,
                             EncodingMethod* const encoding_method) {
  if (num_coeffs == 0) {
    *encoding_method = EncodingMethod::kAllZero;
  } else if (num_coeffs == 1) {
    *encoding_method = EncodingMethod::kDCOnly;
  } else {
    *encoding_method = EncodingMethod::kMethod0;
  }
}

//------------------------------------------------------------------------------
// AOM specific code.

// Fills the extents for the CDF of a symbol 'sym', as well as its number.
static void FillExtents(Symbol sym, uint32_t* const extents,
                        uint32_t* const num_extents) {
  uint32_t c[kCoordMaxSize + 1];
  *num_extents = 0;
  switch (sym) {
    case kAOMAllZero: {
      CoordinateFiller<kCoordMaxSize>::Do(libgav1::kDefaultAllZeroCdf, c,
                                          num_extents);
      break;
    }
    case kAOMEOBPT16: {
      CoordinateFiller<kCoordMaxSize>::Do(libgav1::kDefaultEobPt16Cdf, c,
                                          num_extents);
      break;
    }
    case kAOMEOBPT32: {
      CoordinateFiller<kCoordMaxSize>::Do(libgav1::kDefaultEobPt32Cdf, c,
                                          num_extents);
      break;
    }
    case kAOMEOBPT64: {
      CoordinateFiller<kCoordMaxSize>::Do(libgav1::kDefaultEobPt64Cdf, c,
                                          num_extents);
      break;
    }
    case kAOMEOBPT128: {
      CoordinateFiller<kCoordMaxSize>::Do(libgav1::kDefaultEobPt128Cdf, c,
                                          num_extents);
      break;
    }
    case kAOMEOBPT256: {
      CoordinateFiller<kCoordMaxSize>::Do(libgav1::kDefaultEobPt256Cdf, c,
                                          num_extents);
      break;
    }
    case kAOMEOBPT512: {
      CoordinateFiller<kCoordMaxSize>::Do(libgav1::kDefaultEobPt512Cdf, c,
                                          num_extents);
      break;
    }
    case kAOMEOBPT1024: {
      CoordinateFiller<kCoordMaxSize>::Do(libgav1::kDefaultEobPt1024Cdf, c,
                                          num_extents);
      break;
    }
    case kAOMEOBExtra: {
      CoordinateFiller<kCoordMaxSize>::Do(libgav1::kDefaultEobExtraCdf, c,
                                          num_extents);
      break;
    }
    case kAOMCoeffBaseEOB: {
      CoordinateFiller<kCoordMaxSize>::Do(libgav1::kDefaultCoeffBaseEobCdf, c,
                                          num_extents);
      break;
    }
    case kAOMCoeffBase: {
      CoordinateFiller<kCoordMaxSize>::Do(libgav1::kDefaultCoeffBaseCdf, c,
                                          num_extents);
      break;
    }
    case kAOMCoeffBaseRange: {
      CoordinateFiller<kCoordMaxSize>::Do(libgav1::kDefaultCoeffBaseRangeCdf, c,
                                          num_extents);
      break;
    }
    case kAOMDCSign: {
      CoordinateFiller<kCoordMaxSize>::Do(libgav1::kDefaultDcSignCdf, c,
                                          num_extents);
      break;
    }
    default:
      assert(false);
  }
  // Remove the first element that corresponds to
  // kCoefficientQuantizerContexts. The last one is the CDF itself.
  --*num_extents;
  std::copy(c + 1, c + 1 + *num_extents, extents);
}

// Returns the number of clusters for a given symbol as defined in AOM.
static uint32_t GetNumClusters(Symbol sym) {
  uint32_t extents[kCoordMaxSize] = {};
  uint32_t num_extents;
  FillExtents(sym, extents, &num_extents);
  uint32_t prod = extents[0];
  // The last element is the CDF itself.
  for (uint32_t i = 1; i < num_extents - 1; ++i) prod *= extents[i];

  return prod;
}

// Returns the CDF for a given symbol 'sym' in a given DecoderContext.
static const uint16_t* GetCDF(Symbol sym, uint32_t quantizer_context,
                              uint32_t cluster) {
  const uint32_t num_clusters = GetNumClusters(sym);
  cluster += quantizer_context * num_clusters;
  switch (sym) {
    case kAOMAllZero:
      return GetCDF(cluster, libgav1::kDefaultAllZeroCdf);
    case kAOMEOBPT16:
      return GetCDF(cluster, libgav1::kDefaultEobPt16Cdf);
    case kAOMEOBPT32:
      return GetCDF(cluster, libgav1::kDefaultEobPt32Cdf);
    case kAOMEOBPT64:
      return GetCDF(cluster, libgav1::kDefaultEobPt64Cdf);
    case kAOMEOBPT128:
      return GetCDF(cluster, libgav1::kDefaultEobPt128Cdf);
    case kAOMEOBPT256:
      return GetCDF(cluster, libgav1::kDefaultEobPt256Cdf);
    case kAOMEOBPT512:
      return GetCDF(cluster, libgav1::kDefaultEobPt512Cdf);
    case kAOMEOBPT1024:
      return GetCDF(cluster, libgav1::kDefaultEobPt1024Cdf);
    case kAOMEOBExtra:
      return GetCDF(cluster, libgav1::kDefaultEobExtraCdf);
    case kAOMCoeffBaseEOB:
      return GetCDF(cluster, libgav1::kDefaultCoeffBaseEobCdf);
    case kAOMCoeffBase:
      return GetCDF(cluster, libgav1::kDefaultCoeffBaseCdf);
    case kAOMCoeffBaseRange:
      return GetCDF(cluster, libgav1::kDefaultCoeffBaseRangeCdf);
    case kAOMDCSign:
      return GetCDF(cluster, libgav1::kDefaultDcSignCdf);
    default:
      assert(false);
  }
  return nullptr;
}

uint32_t ResidualIO::GetMethodIndex(EncodingMethod method) {
  assert(method == EncodingMethod::kMethod0 ||
         method == EncodingMethod::kMethod1);
  const uint32_t res = (method == EncodingMethod::kMethod0) ? 0 : 1;
  assert(res < kNumResidualStats);
  return res;
}

// Initializes the CDFs for symbol 'sym' in 'info', based on the default for the
// given quantizer context.
WP2Status InitializeAOMSymbol(Symbol sym, uint32_t quantizer_context,
                              SymbolsInfo* const info) {
  const uint32_t num_clusters = GetNumClusters(sym);

  uint32_t extents[kCoordMaxSize];
  uint32_t num_extents;
  FillExtents(sym, extents, &num_extents);
  const uint32_t range = extents[num_extents - 1];
  info->SetInfo(sym, /*min=*/0, /*max=*/range - 1, num_clusters,
                SymbolsInfo::StorageMethod::kAdaptiveSym);
  for (uint32_t i = 0; i < num_clusters; ++i) {
    const uint16_t* cdf = GetCDF(sym, quantizer_context, i);
    assert(cdf[0] == 0);
    // TODO(vrabaud) all CDFs in the code should have the same max_proba.
    WP2_CHECK_STATUS(info->SetInitialCDF(cdf, /*max_proba=*/1u << 14, i, sym));
  }
  return WP2_STATUS_OK;
}

}  // namespace WP2
