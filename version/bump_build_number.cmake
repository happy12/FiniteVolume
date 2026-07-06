# Runs at build time via `cmake -P`, from the version_header custom target
# in CMakeLists.txt. Increments version/build_number.txt and regenerates
# Version.h from Version.h.in, mirroring the Makefile's version_header
# target.
#
# Required -D variables (set by the parent CMakeLists.txt):
#   BUILD_NUMBER_FILE   version/build_number.txt
#   VERSION_H_IN        version/Version.h.in
#   VERSION_H_OUT       build/generated/Version.h

file(READ "${BUILD_NUMBER_FILE}" _current_build)
string(STRIP "${_current_build}" _current_build)
if(NOT _current_build MATCHES "^[0-9]+$")
    set(_current_build 0)
endif()
math(EXPR BUILD_NUMBER "${_current_build} + 1")
file(WRITE "${BUILD_NUMBER_FILE}" "${BUILD_NUMBER}\n")

configure_file("${VERSION_H_IN}" "${VERSION_H_OUT}" @ONLY)
