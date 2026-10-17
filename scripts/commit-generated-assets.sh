#!/usr/bin/env bash
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

if [[ $# -lt 2 || $# -gt 4 ]]; then
    echo "usage: $0 source_sha generated_asset_dir [remote] [distribution_url]" >&2
    exit 1
fi
SOURCE_SHA="$1"
ASSET_DIR="$(cd "$2" && pwd)"
REMOTE="${3:-origin}"
DISTRIBUTION_URL="${4:-}"
BUILD_NUMBER="$(python3 "${SCRIPT_DIR}/release_metadata.py")"
REPO_ROOT="$(git rev-parse --show-toplevel)"
cd "${REPO_ROOT}"

if [[ "$(git rev-parse HEAD)" != "${SOURCE_SHA}" ]]; then
    echo "error: checkout HEAD does not match source SHA ${SOURCE_SHA}" >&2
    exit 1
fi
if [[ -n "$(git status --porcelain --untracked-files=no)" ]]; then
    echo "error: tracked checkout must be clean before publishing generated assets" >&2
    exit 1
fi
if [[ "${ASSET_DIR}" == "${REPO_ROOT}" ]]; then
    echo "error: generated assets must be staged separately from the checkout" >&2
    exit 1
fi

README_COPY="$(mktemp "${TMPDIR:-/tmp}/taikor-readme.XXXXXX")"
trap 'rm -f "${README_COPY}"' EXIT
python3 - "${ASSET_DIR}" "${README_COPY}" <<'PY'
import pathlib
import re
import sys
import wave

root = pathlib.Path(sys.argv[1])
readme_bytes = (root / "README.md").read_bytes()
readme = readme_bytes.decode("utf-8")
begin = "<!-- peaks-table-begin:"
end = "<!-- peaks-table-end -->"
if readme.count(begin) != 1 or readme.count(end) != 1 or readme.index(begin) >= readme.index(end):
    sys.exit("error: generated README has missing or invalid peaks table markers")
wavs = list((root / "Docs/audio").glob("*.wav"))
if len(wavs) != 27 or any(not re.fullmatch(r"[0-9][0-9]-.*\.wav", p.name) for p in wavs):
    sys.exit("error: expected exactly 27 renderer-owned top-level WAV files")
for wav in wavs:
    with wav.open("rb") as audio:
        header = audio.read(12)
    if header[:4] != b"RIFF" or header[8:12] != b"WAVE" or int.from_bytes(header[4:8], "little") + 8 != wav.stat().st_size:
        sys.exit(f"error: invalid or empty rendered WAV: {wav}")
    with wave.open(str(wav), "rb") as audio:
        if audio.getnchannels() != 2 or audio.getsampwidth() != 2 or audio.getnframes() == 0:
            sys.exit(f"error: expected nonempty stereo 16-bit rendered WAV: {wav}")
        audio.setpos(audio.getnframes() - 1)
        if len(audio.readframes(1)) != 4:
            sys.exit(f"error: incomplete rendered WAV: {wav}")
# The preset demos: demos.json is the renderer's own list, so every MP3 it
# names must be beside it, and nothing else numbered may be.
import json
manifest = json.loads((root / "Docs/audio/demos.json").read_text(encoding="utf-8"))
named = sorted(demo["file"] for demo in manifest["demos"])
mp3s = sorted(p.name for p in (root / "Docs/audio").glob("[0-9][0-9]-*.mp3"))
if not named or named != mp3s:
    sys.exit(f"error: preset MP3s {mp3s} do not match the renderer's manifest {named}")
for mp3 in mp3s:
    data = (root / "Docs/audio" / mp3).read_bytes()
    if len(data) < 4096 or not (data[:3] == b"ID3" or data[0] == 0xFF):
        sys.exit(f"error: {mp3} is not a complete MP3")
for marker in ("<!-- audio-demos-begin", "<!-- audio-demos-end -->"):
    if readme.count(marker) != 1:
        sys.exit("error: generated README has missing or repeated audio demo markers")

png = root / "Docs/screenshots/taikor.png"
with png.open("rb") as screenshot:
    header = screenshot.read(8)
    screenshot.seek(-12, 2)
    trailer = screenshot.read()
if header != b"\x89PNG\r\n\x1a\n" or trailer != b"\0\0\0\0IEND\xaeB`\x82":
    sys.exit("error: screenshot is not a complete PNG")

# This temporary copy is only installed after the source/main guard below.
pathlib.Path(sys.argv[2]).write_bytes(readme_bytes)
PY
if [[ -n "${DISTRIBUTION_URL}" ]]; then
    python3 "${SCRIPT_DIR}/../Tools/update_readme_blocks.py" download --readme "${README_COPY}" \
        --url "${DISTRIBUTION_URL}" --build "${BUILD_NUMBER}" --commit "${SOURCE_SHA}" \
        --repository "${GITHUB_REPOSITORY:-protocodus/virtual-instrument-taikor}"
fi

# Fetch immediately before staging so an older run cannot overwrite newer work.
git fetch --no-tags "${REMOTE}" refs/heads/main
if [[ "$(git rev-parse FETCH_HEAD)" != "${SOURCE_SHA}" ]]; then
    echo "::notice::main advanced past ${SOURCE_SHA}; skipping generated assets."
    exit 0
fi

# Only numbered renderer outputs are owned here. Historical preview directories
# and manually named top-level WAVs must remain untouched.
git rm -q --ignore-unmatch -- ':(top,glob)Docs/audio/[0-9][0-9]-*.wav' ':(top,glob)Docs/audio/[0-9][0-9]-*.mp3'
mkdir -p Docs/audio Docs/screenshots
cp "${README_COPY}" README.md
cp "${ASSET_DIR}/Docs/screenshots/taikor.png" Docs/screenshots/taikor.png
generated_paths=()
for wav in "${ASSET_DIR}"/Docs/audio/[0-9][0-9]-*.wav; do
    destination="Docs/audio/${wav##*/}"
    cp "${wav}" "${destination}"
    generated_paths+=("${destination}")
done
for mp3 in "${ASSET_DIR}"/Docs/audio/[0-9][0-9]-*.mp3; do
    destination="Docs/audio/${mp3##*/}"
    cp "${mp3}" "${destination}"
    generated_paths+=("${destination}")
done
git add -- README.md Docs/screenshots/taikor.png "${generated_paths[@]}"
if git diff --cached --quiet; then
    echo "::notice::Generated assets are unchanged; no commit needed."
    exit 0
fi

git -c user.name='github-actions[bot]' \
    -c user.email='41898282+github-actions[bot]@users.noreply.github.com' \
    commit -m "chore: refresh generated media for ${SOURCE_SHA:0:12}"
if git push "${REMOTE}" HEAD:refs/heads/main; then
    echo "::notice::Published generated assets for ${SOURCE_SHA}."
    exit 0
fi

# A concurrent merge may win between the fetch and push. Never force an update
# or apply generated media to a different source revision.
git fetch --no-tags "${REMOTE}" refs/heads/main
if [[ "$(git rev-parse FETCH_HEAD)" != "${SOURCE_SHA}" ]]; then
    echo "::notice::main advanced during publication; skipping generated assets."
    exit 0
fi
echo "error: generated asset push failed while main still matches ${SOURCE_SHA}" >&2
exit 1
