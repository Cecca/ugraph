#pragma once

#include "prelude.hpp"

#define COLOR_STREAM_FN(name, code)\
    inline std::ostream& name(std::ostream& stream)\
    { stream << code; return stream; }

namespace termcolor {

    COLOR_STREAM_FN(reset,   "\033[00m")
    COLOR_STREAM_FN(grey,    "\033[30m")
    COLOR_STREAM_FN(red,     "\033[31m")
    COLOR_STREAM_FN(green,   "\033[32m")
    COLOR_STREAM_FN(yellow,  "\033[33m")
    COLOR_STREAM_FN(blue,    "\033[34m")
    COLOR_STREAM_FN(magenta, "\033[35m")
    COLOR_STREAM_FN(cyan,    "\033[36m")
    COLOR_STREAM_FN(white,   "\033[37m")

}
