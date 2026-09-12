#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
BUILD_DIR="${BUILD_DIR:-${PROJECT_DIR}/build-linux}"
CONFIG="${CONFIG:-Release}"
artefacts="${BUILD_DIR}/Taikor_artefacts/${CONFIG}"
dist="${BUILD_DIR}/dist"
version="$(sed -n 's/^CMAKE_PROJECT_VERSION:STATIC=//p' "${BUILD_DIR}/CMakeCache.txt")"
prefix="$(python3 "${SCRIPT_DIR}/release_metadata.py" --version "${version}")"

test -d "${artefacts}/VST3/Taikor.vst3"
test -s "${artefacts}/CLAP/Taikor.clap"
test -x "${artefacts}/Standalone/Taikor"
mkdir -p "${dist}"
archive="${dist}/${prefix}-Linux-x64.tar.gz"
contents=("VST3/Taikor.vst3" "CLAP/Taikor.clap" "Standalone/Taikor")
for notice in "${PROJECT_DIR}/LICENSE" "${PROJECT_DIR}/THIRD_PARTY_NOTICES.md" \
              "${PROJECT_DIR}"/ThirdParty/*-LICENSE.md \
              "${PROJECT_DIR}/ThirdParty/ECHOTHIEF-LICENSE.pdf" \
              "${PROJECT_DIR}/ThirdParty/VOXENGO-IMPULSES-LICENSE.txt"; do
    test -s "${notice}"
    cp "${notice}" "${artefacts}/$(basename "${notice}")"
    contents+=("$(basename "${notice}")")
done
# Only remove this platform's earlier release archives from the staging tree.
rm -f "${dist}/Taikor-"*-Linux-x64.tar.gz
tar -czf "${archive}" -C "${artefacts}" "${contents[@]}"
test -s "${archive}"
tar -tzf "${archive}" >/dev/null
printf 'Packaged %s\n' "${archive}"
