#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
BUILD_DIR="${BUILD_DIR:-${PROJECT_DIR}/build-linux}"
CONFIG="${CONFIG:-Release}"
artefacts="${BUILD_DIR}/Taikor_artefacts/${CONFIG}"
dist="${BUILD_DIR}/dist"
version_lines="$(sed -n 's/^CMAKE_PROJECT_VERSION:STATIC=//p' "${BUILD_DIR}/CMakeCache.txt")"
version_count="$(printf '%s\n' "${version_lines}" | awk 'NF { count += 1 } END { print count + 0 }')"
if (( version_count != 1 )); then
    echo "error: expected exactly one CMake project version" >&2
    exit 1
fi
version="${version_lines}"
prefix="$(python3 "${SCRIPT_DIR}/release_metadata.py" --version "${version}")"

required=(
    "VST3/Taikor.vst3/Contents/x86_64-linux/Taikor.so"
    "VST3/Taikor.vst3/Contents/Resources/moduleinfo.json"
    "CLAP/Taikor.clap"
    "Standalone/Taikor"
)
for required_file in "${required[@]}"; do
    if [[ ! -s "${artefacts}/${required_file}" ]]; then
        echo "error: missing or empty distribution file: ${artefacts}/${required_file}" >&2
        exit 1
    fi
done
if [[ ! -x "${artefacts}/Standalone/Taikor" ]]; then
    echo "error: standalone artifact is not executable" >&2
    exit 1
fi
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
archive_listing="$(tar -tzf "${archive}")"
for required_file in "${required[@]}"; do
    if ! printf '%s\n' "${archive_listing}" | grep -Fxq "${required_file}"; then
        echo "error: archive is missing required entry: ${required_file}" >&2
        exit 1
    fi
done
printf 'Packaged %s\n' "${archive}"
