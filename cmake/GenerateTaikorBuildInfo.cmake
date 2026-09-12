if(NOT DEFINED TAIKOR_BUILD_INFO_DIR OR NOT DEFINED TAIKOR_BUILD_VERSION)
    message(FATAL_ERROR "Taikor build identity requires an output directory and version")
endif()

# A UTC timestamp distinguishes builds across machines and build directories,
# unlike a local counter that could label unrelated binaries with the same ID.
string(TIMESTAMP TAIKOR_BUILD_NUMBER "%Y%m%d.%H%M%S" UTC)
file(MAKE_DIRECTORY "${TAIKOR_BUILD_INFO_DIR}")
file(WRITE "${TAIKOR_BUILD_INFO_DIR}/TaikorBuildInfo.h"
    "// Generated at build time; do not edit.\n"
    "#pragma once\n"
    "#define TAIKOR_BUILD_VERSION \"${TAIKOR_BUILD_VERSION}\"\n"
    "#define TAIKOR_BUILD_NUMBER \"${TAIKOR_BUILD_NUMBER}\"\n")
message(STATUS "Taikor ${TAIKOR_BUILD_VERSION} / build ${TAIKOR_BUILD_NUMBER} (UTC)")
