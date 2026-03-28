cmake_minimum_required(VERSION 3.16)
project(ff_poly VERSION 2.0.0 LANGUAGES C)

# --- ff_poly core sources ---
set(FF_POLY_SOURCES
    ff.c
    ff2k.c
    ffext.c
    ffpoly.c
    ffpolyalloc.c
    ffpolybig.c
    ffpolyfromroots.c
    ffpolymultipoint.c
    ffpolysmall.c
    polyparse.c
)

# --- zn_poly sources ---
set(ZN_POLY_SOURCES
    zn_poly/array.c
    zn_poly/ks_support.c
    zn_poly/misc.c
    zn_poly/mpn_mulmid.c
    zn_poly/mul.c
    zn_poly/mul_fft.c
    zn_poly/mul_ks.c
    zn_poly/mulmid.c
    zn_poly/mulmid_ks.c
    zn_poly/nuss.c
    zn_poly/pack.c
    zn_poly/pmf.c
    zn_poly/pmfvec_fft.c
    zn_poly/tuning.c
    zn_poly/zn_mod.c
)

add_library(ff_poly STATIC ${FF_POLY_SOURCES} ${ZN_POLY_SOURCES})

target_include_directories(ff_poly
    # Private: ff_poly's own sources use bare includes like "ff.h"
    PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}
    PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/zn_poly
    # Public: consumers get the symlink-based include layout
    PUBLIC ${CMAKE_CURRENT_BINARY_DIR}/include
)

# Create an include layout so that #include "ff_poly.h" and #include "ff_poly/ff.h" both work
# ff_poly.h -> include/ff_poly.h
# ff_poly/ff.h -> include/ff_poly/ff.h  (symlink to source dir)
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/include)
file(CREATE_LINK
    ${CMAKE_CURRENT_SOURCE_DIR}/ff_poly.h
    ${CMAKE_CURRENT_BINARY_DIR}/include/ff_poly.h
    SYMBOLIC
)
file(CREATE_LINK
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_CURRENT_BINARY_DIR}/include/ff_poly
    SYMBOLIC
)

include(CheckCCompilerFlag)
check_c_compiler_flag("-march=native" FF_POLY_SUPPORTS_MARCH_NATIVE)
if(FF_POLY_SUPPORTS_MARCH_NATIVE)
    set(FF_POLY_ARCH_FLAG "-march=native")
else()
    set(FF_POLY_ARCH_FLAG "")
endif()

target_compile_options(ff_poly PRIVATE
    -O3
    -fomit-frame-pointer
    -funroll-loops
    ${FF_POLY_ARCH_FLAG}
    -pedantic
    -std=gnu99
)

# On aarch64, unsigned long and uint64_t are distinct types (unsigned long long),
# producing harmless format and pointer warnings throughout the codebase
if(CMAKE_SYSTEM_PROCESSOR MATCHES "aarch64|arm64|ARM64")
    target_compile_options(ff_poly PRIVATE
        -Wno-format
        -Wno-incompatible-pointer-types
        -Wno-strict-prototypes
        -Wno-empty-body
        -Wno-return-type
    )
endif()

# Find GMP (static)
find_library(GMP_STATIC_LIBRARY libgmp.a HINTS
    /opt/homebrew/lib
    /opt/homebrew/opt/gmp/lib
    /usr/local/lib
    /usr/lib
    /usr/lib/x86_64-linux-gnu
    /usr/lib/aarch64-linux-gnu
)
find_path(GMP_INCLUDE_DIR gmp.h HINTS
    /opt/homebrew/include
    /opt/homebrew/opt/gmp/include
)

if(NOT GMP_STATIC_LIBRARY)
    message(FATAL_ERROR "GMP static library (libgmp.a) not found. Install libgmp-dev.")
endif()
if(NOT GMP_INCLUDE_DIR)
    message(FATAL_ERROR "GMP headers not found. Install libgmp-dev.")
endif()

target_include_directories(ff_poly PUBLIC ${GMP_INCLUDE_DIR})
target_link_libraries(ff_poly PUBLIC ${GMP_STATIC_LIBRARY} m)
