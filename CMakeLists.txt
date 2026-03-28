cmake_minimum_required(VERSION 3.16)
project(classpoly VERSION 1.0.3 LANGUAGES C)

# Require 64-bit
if(NOT CMAKE_SIZEOF_VOID_P EQUAL 8)
    message(FATAL_ERROR "classpoly requires a 64-bit platform")
endif()

# Detect native architecture for optimized builds
include(CheckCCompilerFlag)
check_c_compiler_flag("-march=native" COMPILER_SUPPORTS_MARCH_NATIVE)
if(COMPILER_SUPPORTS_MARCH_NATIVE)
    set(ARCH_FLAG "-march=native")
else()
    set(ARCH_FLAG "")
endif()

# Build ff_poly as a static library
add_subdirectory(ff_poly_v2.0.0)

# --- classpoly sources ---
set(CLASSPOLY_SOURCES
    src/classpoly_program.c
    src/bipoly.c
    src/class_inv.c
    src/classpoly.c
    src/classpoly_crt.c
    src/classpoly_inv.c
    src/compute_classpoly.c
    src/compute_classpoly_parallel.c
    src/phi_eval.c
    src/phi_gcd.c
    src/phi_poly.c
    src/polycosts.c
    src/velu.c
    src/crt.c
    src/cstd.c
    src/ecurve.c
    src/ecurve_ladic.c
    src/evec.c
    src/findcurve.c
    src/mpzutil.c
    src/ntutil.c
    src/prime.c
    src/pickprimes.c
    src/iqclass.c
    src/qform.c
    src/table.c
    src/tecurve.c
)

add_executable(classpoly ${CLASSPOLY_SOURCES})

target_include_directories(classpoly PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src)

target_compile_options(classpoly PRIVATE
    -O3
    -fomit-frame-pointer
    -funroll-loops
    ${ARCH_FLAG}
    -pedantic
    -std=gnu99
)

# Link ff_poly statically (brings in GMP and math library transitively)
target_link_libraries(classpoly PRIVATE ff_poly)

# Place the binary in the build directory root
set_target_properties(classpoly PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}
)
