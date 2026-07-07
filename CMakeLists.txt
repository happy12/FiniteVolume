cmake_minimum_required(VERSION 3.20)
project(FiniteVolume LANGUAGES CXX RC)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

set(CMAKE_CXX_FLAGS_DEBUG "-g -O0")
set(CMAKE_CXX_FLAGS_RELEASE "-O2 -ffast-math")

# Single output location regardless of config, matching the Makefile. Ninja
# Multi-Config appends a per-config subdirectory unless the per-config
# variables are set explicitly too.
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/build/bin)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_DEBUG ${CMAKE_SOURCE_DIR}/build/bin)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_RELEASE ${CMAKE_SOURCE_DIR}/build/bin)

set(GENERATED_DIR ${CMAKE_SOURCE_DIR}/build/generated)
set(VERSION_H ${GENERATED_DIR}/Version.h)

file(GLOB FV_SOURCES CONFIGURE_DEPENDS ${CMAKE_SOURCE_DIR}/src/*.cpp)

# Regenerates build/generated/Version.h from version/Version.h.in on every
# build, bumping version/build_number.txt each time -- same behavior as the
# Makefile's version_header target.
add_custom_target(version_header ALL
    COMMAND ${CMAKE_COMMAND} -E make_directory ${GENERATED_DIR}
    COMMAND ${CMAKE_COMMAND}
        -DBUILD_NUMBER_FILE=${CMAKE_SOURCE_DIR}/version/build_number.txt
        -DVERSION_H_IN=${CMAKE_SOURCE_DIR}/version/Version.h.in
        -DVERSION_H_OUT=${VERSION_H}
        -P ${CMAKE_SOURCE_DIR}/version/bump_build_number.cmake
    BYPRODUCTS ${VERSION_H}
    COMMENT "Bumping FiniteVolume build number"
)

add_executable(FiniteVolume ${FV_SOURCES} src/FiniteVolume.rc)
add_dependencies(FiniteVolume version_header)

target_include_directories(FiniteVolume PRIVATE ${CMAKE_SOURCE_DIR}/include ${GENERATED_DIR})
target_compile_options(FiniteVolume PRIVATE -Wall -Wextra -fopenmp)

# Release-only cross-TU inlining -- e.g. GradientReconstruction.cpp's
# face_gradient() called from AdvectionDiffusionFVMSolver.cpp/
# NavierStokesFVMSolver.cpp's hot per-face loops, which per-TU compilation
# can't inline across. INTERPROCEDURAL_OPTIMIZATION_RELEASE (rather than a
# plain -flto compile/link flag) lets CMake pair the compile+link flags
# itself and only if the toolchain actually supports it.
include(CheckIPOSupported)
check_ipo_supported(RESULT ipo_supported)
if(ipo_supported)
    set_target_properties(FiniteVolume PROPERTIES INTERPROCEDURAL_OPTIMIZATION_RELEASE ON)
endif()

# -I lets windres find the generated Version.h; OBJECT_DEPENDS makes the
# resource recompile whenever the icon or the version header changes.
set_source_files_properties(src/FiniteVolume.rc PROPERTIES
    COMPILE_FLAGS "-I${GENERATED_DIR}"
    OBJECT_DEPENDS "${CMAKE_SOURCE_DIR}/icons/favicon256.ico;${VERSION_H}"
)

target_link_options(FiniteVolume PRIVATE -static -fopenmp)
