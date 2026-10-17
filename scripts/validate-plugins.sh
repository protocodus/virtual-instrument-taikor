#!/usr/bin/env bash
# Host validation of a Release build: pluginval at strictness 10 on the VST3
# on every platform (pluginval 1.0.4 has no CLAP host) and auval on macOS.
# Shared by main-build.yml and release.yml so both validate identically.
set -euo pipefail

BUILD_DIR="${1:?usage: validate-plugins.sh build-dir [config]}"
CONFIG="${2:-Release}"
PLUGINVAL_VERSION=v1.0.4
artefacts="${BUILD_DIR}/Taikor_artefacts/${CONFIG}"
tools="${BUILD_DIR}/pluginval"
mkdir -p "${tools}"

case "$(uname -s)" in
    Darwin)
        archive=pluginval_macOS.zip
        sha256=3c4c533bda0c5059eea3ddaea752d757ee2025041f0f47e6bcb0e87f6082b29f
        pluginval="${tools}/pluginval.app/Contents/MacOS/pluginval" ;;
    Linux)
        archive=pluginval_Linux.zip
        sha256=c01c49d8063965c4c2dea8324468336768f5c9139e0b1caebde14c2400b55352
        pluginval="${tools}/pluginval" ;;
    MINGW* | MSYS* | CYGWIN*)
        archive=pluginval_Windows.zip
        sha256=c08e61ce3b96db41636f8ec7e76f4c7e2c13ebdac7fa1b5a1f52b4f32ec715ab
        pluginval="${tools}/pluginval.exe" ;;
    *)
        echo "error: unsupported platform $(uname -s)" >&2
        exit 1 ;;
esac

curl -sSfL --retry 3 -o "${tools}/${archive}" \
    "https://github.com/Tracktion/pluginval/releases/download/${PLUGINVAL_VERSION}/${archive}"
if command -v sha256sum >/dev/null 2>&1; then
    echo "${sha256}  ${tools}/${archive}" | sha256sum -c -
else
    echo "${sha256}  ${tools}/${archive}" | shasum -a 256 -c -
fi
case "${archive}" in
    pluginval_macOS.zip) ditto -x -k "${tools}/${archive}" "${tools}" ;;
    pluginval_Linux.zip) unzip -q -o "${tools}/${archive}" -d "${tools}" ;;
    *) 7z x -y -o"${tools}" "${tools}/${archive}" >/dev/null ;;
esac

validate=("${pluginval}" --strictness-level 10 --verbose
          --output-dir "${BUILD_DIR}/pluginval-logs"
          --validate "${artefacts}/VST3/Taikor.vst3")
if [[ "$(uname -s)" == "Linux" && -z "${DISPLAY:-}" ]]; then
    xvfb-run -a "${validate[@]}"
else
    "${validate[@]}"
fi

if [[ "$(uname -s)" == "Darwin" ]]; then
    components="${HOME}/Library/Audio/Plug-Ins/Components"
    mkdir -p "${components}"
    rm -rf "${components}/Taikor.component"
    ditto "${artefacts}/AU/Taikor.component" "${components}/Taikor.component"
    # The registrar caches the component list; restarting it rescans.
    killall -9 AudioComponentRegistrar 2>/dev/null || true
    auval -strict -v aumu Tko1 Tkor
fi
