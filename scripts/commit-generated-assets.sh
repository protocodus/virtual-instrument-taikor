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
python3 - "${ASSET_DIR}" "${README_COPY}" "${SOURCE_SHA}" "${DISTRIBUTION_URL}" "$#" "${BUILD_NUMBER}" <<'PY'
import os
import hashlib
import json
import math
import pathlib
import re
import struct
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

# The suite is independently owned; never glob-copy its large working masters
# or trust a passing QA report from a different source revision.
suite = root / "Docs/audio/suite"
suite_names = {"taikor-suite-96k24.wav", "taikor-suite.mp3", "taikor-suite.mid",
               "taikor-suite-score.csv", "taikor-suite-qa.json"}
if not suite.is_dir() or {p.name for p in suite.iterdir()} != suite_names:
    sys.exit("error: expected exactly the five suite delivery files")
if any(not (suite / name).is_file() or (suite / name).is_symlink()
       or (suite / name).stat().st_size == 0 for name in suite_names):
    sys.exit("error: suite deliveries must be nonempty regular files")

def digest(path):
    result = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            result.update(block)
    return result.hexdigest()

def reject(message):
    sys.exit(f"error: suite {message}")

try:
    qa = json.loads((suite / "taikor-suite-qa.json").read_text())
    if qa["schema_version"] != 1 or qa["status"] != "passed" or qa["failures"] != []:
        reject("QA did not pass")
    checkout = pathlib.Path.cwd()
    sources = {str(path.relative_to(checkout)): digest(path)
               for path in sorted((checkout / "Source/DSP").glob("*")) if path.is_file()}
    sources["Tools/RenderSuite.cpp"] = digest(checkout / "Tools/RenderSuite.cpp")
    # The producer records this legacy alternative when it exists, too.
    alternative = checkout / "Tools/RenderTaikoSuite.cpp"
    if alternative.is_file():
        sources["Tools/RenderTaikoSuite.cpp"] = digest(alternative)
    if qa["source_files"] != sources:
        reject("source manifest does not match the publishing checkout")
    if qa["production_script"]["sha256"] != digest(checkout / "Tools/ProduceSuite.py"):
        reject("production script does not match the publishing checkout")
    processing = qa["processing"]
    if processing["room_ir_sha256"] != digest(checkout / "Assets/ImpulseResponses/theater.wav"):
        reject("room impulse does not match the publishing checkout")
    if any(processing[field] != 384000 for field in
           ("engine_sample_rate_hz", "master_sample_rate_hz")):
        reject("was not produced at the required 384 kHz quality")
    if qa["source"]["sample_rate_hz"] != 384000 or qa["source"]["frames"] != 23040000:
        reject("source is not a one-minute 384 kHz render")
    for name, rate in (("taikor-suite-96k24.wav", 96000), ("taikor-suite.mp3", 48000)):
        artifact = qa["artifacts"][name]
        path = suite / name
        if artifact["sha256"] != digest(path) or artifact["bytes"] != path.stat().st_size:
            reject(f"delivery does not match its QA report: {name}")
        if (artifact["sample_rate_hz"] != rate or artifact["channels"] != 2
                or artifact["frames"] != rate * 60 or artifact["duration_seconds"] != 60.0
                or artifact["clipped_sample_count"] != 0 or artifact["nonfinite_sample_count"] != 0):
            reject(f"invalid delivery measurements: {name}")
        for field in ("rms_dbfs", "estimated_true_peak_dbtp"):
            value = artifact[field]
            if not isinstance(value, (float, int)) or not math.isfinite(value):
                reject(f"missing finite {field}: {name}")
        if artifact["rms_dbfs"] < -100 or artifact["estimated_true_peak_dbtp"] >= 0:
            reject(f"silent or clipping delivery: {name}")
    for name in ("taikor-suite.mid", "taikor-suite-score.csv"):
        artifact = qa["performance_files"][name]
        if artifact["sha256"] != digest(suite / name) or artifact["bytes"] != (suite / name).stat().st_size:
            reject(f"performance does not match its QA report: {name}")
except (KeyError, TypeError, ValueError, OSError) as error:
    reject(f"invalid QA or missing provenance: {error}")

# Read RIFF chunks directly so PCM24 with either ordinary PCM or the extensible
# PCM header emitted by FFmpeg works on older Python versions as well.
with (suite / "taikor-suite-96k24.wav").open("rb") as audio:
    header = audio.read(12)
    length = (suite / "taikor-suite-96k24.wav").stat().st_size
    if header[:4] != b"RIFF" or header[8:] != b"WAVE" or int.from_bytes(header[4:8], "little") + 8 != length:
        reject("listening WAV has an invalid RIFF header")
    fmt, data_bytes = None, None
    while audio.tell() < length:
        chunk = audio.read(8)
        if len(chunk) != 8:
            reject("listening WAV has a truncated chunk")
        size = int.from_bytes(chunk[4:], "little")
        if audio.tell() + size + (size & 1) > length:
            reject("listening WAV has a truncated payload")
        if chunk[:4] == b"fmt ":
            if fmt is not None or size > 1024:
                reject("listening WAV has an invalid format chunk")
            fmt = audio.read(size)
        else:
            if chunk[:4] == b"data":
                if data_bytes is not None:
                    reject("listening WAV has duplicate audio data")
                data_bytes = size
            audio.seek(size, 1)
        audio.seek(size & 1, 1)
    if fmt is None or len(fmt) < 16:
        reject("listening WAV has no PCM format")
    tag, channels, rate, byte_rate, align, bits = struct.unpack_from("<HHIIHH", fmt)
    pcm_guid = bytes.fromhex("0100000000001000800000aa00389b71")
    pcm = tag == 1 or (tag == 65534 and len(fmt) >= 40
                      and int.from_bytes(fmt[18:20], "little") == 24 and fmt[24:40] == pcm_guid)
    if not pcm or (channels, rate, byte_rate, align, bits, data_bytes) != (2, 96000, 576000, 6, 24, 34560000):
        reject("listening WAV must be exactly 60 seconds of stereo 96 kHz PCM24")

png = root / "Docs/screenshots/taikor-standalone.png"
with png.open("rb") as screenshot:
    header = screenshot.read(8)
    screenshot.seek(-12, 2)
    trailer = screenshot.read()
if header != b"\x89PNG\r\n\x1a\n" or trailer != b"\0\0\0\0IEND\xaeB`\x82":
    sys.exit("error: screenshot is not a complete PNG")

if sys.argv[5] == "4":
    source_sha, url = sys.argv[3:5]
    repository = os.environ.get("GITHUB_REPOSITORY", "protocodus/virtual-instrument-taikor")
    pattern = rf"https://github\.com/{re.escape(repository)}/actions/runs/[0-9]+/artifacts/[0-9]+"
    if not re.fullmatch(pattern, url):
        sys.exit("error: distribution URL must be a GitHub Actions artifact URL for this repository")
    link_begin = b"<!-- distribution-link-begin -->"
    link_end = b"<!-- distribution-link-end -->"
    if readme_bytes.count(link_begin) != 1 or readme_bytes.count(link_end) != 1:
        sys.exit("error: generated README must have exactly one pair of distribution link markers")
    start = readme_bytes.index(link_begin) + len(link_begin)
    stop = readme_bytes.index(link_end)
    peak_start = readme_bytes.index(begin.encode())
    peak_stop = readme_bytes.index(end.encode()) + len(end)
    if start > stop or (start < peak_stop and stop > peak_start):
        sys.exit("error: distribution link markers are reversed or overlap the peaks table")
    link = (f"\n**[Download latest distribution]({url})** — build {sys.argv[6]}, built from "
            f"[`{source_sha[:12]}`](https://github.com/{repository}/commit/{source_sha}).\n")
    readme_bytes = readme_bytes[:start] + link.encode("utf-8") + readme_bytes[stop:]

# Keep every byte outside the link block, including the rendered peaks table.
# This temporary copy is only installed after the source/main guard below.
pathlib.Path(sys.argv[2]).write_bytes(readme_bytes)
PY

# Fetch immediately before staging so an older run cannot overwrite newer work.
git fetch --no-tags "${REMOTE}" refs/heads/main
if [[ "$(git rev-parse FETCH_HEAD)" != "${SOURCE_SHA}" ]]; then
    echo "::notice::main advanced past ${SOURCE_SHA}; skipping generated assets."
    exit 0
fi

# Only numbered diagnostics and five named suite deliveries are owned here.
# Historical previews and manually named files must remain untouched.
git rm -q --ignore-unmatch -- ':(top,glob)Docs/audio/[0-9][0-9]-*.wav'
mkdir -p Docs/audio/suite Docs/screenshots
cp "${README_COPY}" README.md
cp "${ASSET_DIR}/Docs/screenshots/taikor-standalone.png" Docs/screenshots/taikor-standalone.png
generated_paths=()
for wav in "${ASSET_DIR}"/Docs/audio/[0-9][0-9]-*.wav; do
    destination="Docs/audio/${wav##*/}"
    cp "${wav}" "${destination}"
    generated_paths+=("${destination}")
done
for name in taikor-suite-96k24.wav taikor-suite.mp3 taikor-suite.mid \
            taikor-suite-score.csv taikor-suite-qa.json; do
    destination="Docs/audio/suite/${name}"
    cp "${ASSET_DIR}/Docs/audio/suite/${name}" "${destination}"
    generated_paths+=("${destination}")
done
git add -- README.md Docs/screenshots/taikor-standalone.png "${generated_paths[@]}"
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
