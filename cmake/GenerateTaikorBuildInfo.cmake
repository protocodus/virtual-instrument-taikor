if(NOT DEFINED TAIKOR_BUILD_INFO_DIR OR NOT DEFINED TAIKOR_BUILD_VERSION)
    message(FATAL_ERROR "Taikor build identity requires an output directory and version")
endif()

# The visible identity is X.Y.Z-build.<run>, matching the package names; a
# local build has no CI run and says so. The UTC timestamp still tells local
# builds apart from each other.
if(DEFINED ENV{GITHUB_RUN_NUMBER} AND "$ENV{GITHUB_RUN_NUMBER}" MATCHES "^[1-9][0-9]*$")
    set(TAIKOR_BUILD_NUMBER "$ENV{GITHUB_RUN_NUMBER}")
else()
    set(TAIKOR_BUILD_NUMBER "local")
endif()
string(TIMESTAMP TAIKOR_BUILD_TIME "%Y-%m-%d %H:%M:%S" UTC)
file(MAKE_DIRECTORY "${TAIKOR_BUILD_INFO_DIR}")
file(WRITE "${TAIKOR_BUILD_INFO_DIR}/TaikorBuildInfo.h"
    "// Generated at build time; do not edit.\n"
    "#pragma once\n"
    "#define TAIKOR_BUILD_VERSION \"${TAIKOR_BUILD_VERSION}\"\n"
    "#define TAIKOR_BUILD_NUMBER \"${TAIKOR_BUILD_NUMBER}\"\n"
    "#define TAIKOR_BUILD_TIME \"${TAIKOR_BUILD_TIME}\"\n")
message(STATUS "Taikor ${TAIKOR_BUILD_VERSION}-build.${TAIKOR_BUILD_NUMBER} (${TAIKOR_BUILD_TIME} UTC)")
