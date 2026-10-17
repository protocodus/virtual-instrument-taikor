#!/usr/bin/env python3
"""Render the README demo set twice and require the expected count and identical bytes.

usage: check_audio_demos.py <TaikorRenderDemos> <work-dir> <expected-count>

The renderer itself rejects non-finite, silent and limiter-clipped takes. This
check adds what only a second render can show: the shipping engine renders the
same bytes from the same score, so an unchanged engine commits nothing.
"""

import pathlib
import shutil
import subprocess
import sys


def render(renderer: str, directory: pathlib.Path) -> dict[str, bytes]:
    shutil.rmtree(directory, ignore_errors=True)
    directory.mkdir(parents=True)
    subprocess.run([renderer, str(directory)], check=True)
    return {path.name: path.read_bytes() for path in sorted(directory.glob("*.wav"))}


def main() -> int:
    renderer, work, expected = sys.argv[1], pathlib.Path(sys.argv[2]), int(sys.argv[3])
    first = render(renderer, work / "first")
    if len(first) != expected:
        print(f"error: rendered {len(first)} demo WAVs, expected {expected}", file=sys.stderr)
        return 1
    second = render(renderer, work / "second")
    differing = sorted(name for name in first if first[name] != second.get(name))
    if differing or first.keys() != second.keys():
        print("error: a second render differs: " + ", ".join(differing or ["file set"]),
              file=sys.stderr)
        return 1
    print(f"Taikor: {expected} demo WAVs rendered twice with identical bytes")
    return 0


if __name__ == "__main__":
    sys.exit(main())
