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

# On aarch64, unsigned long and uint64_t are distinct types (unsigned long long),
# producing harmless format and pointer warnings throughout the codebase
if(CMAKE_SYSTEM_PROCESSOR MATCHES "aarch64|arm64|ARM64")
    target_compile_options(classpoly PRIVATE
        -Wno-format
        -Wno-incompatible-pointer-types
        -Wno-strict-prototypes
        -Wno-empty-body
        -Wno-return-type
    )
endif()

# Compile-time default for modular polynomial directory
set(CLASSPOLY_PHI_DIR "${CMAKE_INSTALL_PREFIX}/share/classpoly/phi_files"
    CACHE STRING "Default directory for modular polynomial files")
target_compile_definitions(classpoly PRIVATE
    CLASSPOLY_DEFAULT_PHI_DIR="${CLASSPOLY_PHI_DIR}"
)

# Link ff_poly statically (brings in GMP and math library transitively)
target_link_libraries(classpoly PRIVATE ff_poly)

# Place the binary in the build directory root
set_target_properties(classpoly PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}
)

# --- test target ---
# Same sources as classpoly but replace main with test harness
set(CLASSPOLY_LIB_SOURCES ${CLASSPOLY_SOURCES})
list(REMOVE_ITEM CLASSPOLY_LIB_SOURCES src/classpoly_program.c)

enable_testing()

add_executable(test_classpoly tests/test_classpoly.c ${CLASSPOLY_LIB_SOURCES})
target_include_directories(test_classpoly PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src)
target_compile_options(test_classpoly PRIVATE
    -O3
    -fomit-frame-pointer
    -funroll-loops
    ${ARCH_FLAG}
    -pedantic
    -std=gnu99
)
if(CMAKE_SYSTEM_PROCESSOR MATCHES "aarch64|arm64|ARM64")
    target_compile_options(test_classpoly PRIVATE
        -Wno-format
        -Wno-incompatible-pointer-types
        -Wno-strict-prototypes
        -Wno-empty-body
        -Wno-return-type
    )
endif()
target_link_libraries(test_classpoly PRIVATE ff_poly)
set_target_properties(test_classpoly PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}
)

add_test(NAME classpoly_correctness
    COMMAND test_classpoly ${CMAKE_SOURCE_DIR}/tests/data
)

# --- install targets ---
install(TARGETS classpoly RUNTIME DESTINATION bin)

# Install modular polynomial files by extracting tarballs at install time
install(CODE "
    file(MAKE_DIRECTORY \"${CLASSPOLY_PHI_DIR}\")
    execute_process(
        COMMAND sh -c \"ls '${CMAKE_SOURCE_DIR}/phi_polys'/*.tar.gz | xargs -P 4 -I {} tar xzf {} -C '${CLASSPOLY_PHI_DIR}'\"
        RESULT_VARIABLE _phi_result)
    if(NOT _phi_result EQUAL 0)
        message(FATAL_ERROR \"Failed to extract modular polynomials\")
    endif()
    message(STATUS \"Extracted modular polynomials to ${CLASSPOLY_PHI_DIR}\")
")
