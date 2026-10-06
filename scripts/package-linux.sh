#!/usr/bin/env bash
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
python3 "${SCRIPT_DIR}/release_metadata.py" --package --platform Linux --arch "${ARTIFACT_ARCH:-x64}" \
    --build-dir "${BUILD_DIR:-${PROJECT_DIR}/build-linux}" --config "${CONFIG:-Release}"
