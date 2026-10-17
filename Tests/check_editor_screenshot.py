#!/usr/bin/env python3
"""Run the processor suite for the editor screenshot and check the image it wrote.

usage: check_editor_screenshot.py <TaikorPluginProcessorTests> <output.png> <committed.png>

The image must be a complete PNG, have the committed screenshot's pixel size,
and not be blank. Pixel equality is not required: font rasterisation differs
between machines even when the editor source is unchanged.
"""

import os
import pathlib
import struct
import subprocess
import sys
import zlib

SIGNATURE = b"\x89PNG\r\n\x1a\n"
TRAILER = b"\0\0\0\0IEND\xaeB`\x82"


def chunks(data: bytes):
    offset = len(SIGNATURE)
    while offset < len(data):
        length, kind = struct.unpack(">I4s", data[offset:offset + 8])
        yield kind, data[offset + 8:offset + 8 + length]
        offset += 12 + length


def size(data: bytes) -> tuple[int, int]:
    return struct.unpack(">II", data[16:24])


def main() -> int:
    tests, output, committed = sys.argv[1], pathlib.Path(sys.argv[2]), pathlib.Path(sys.argv[3])
    output.parent.mkdir(parents=True, exist_ok=True)
    output.unlink(missing_ok=True)
    subprocess.run([tests], check=True, env=dict(os.environ, TAIKOR_EDITOR_SNAPSHOT=str(output)))

    if not output.is_file():
        print(f"error: the suite did not render {output}", file=sys.stderr)
        return 1
    data = output.read_bytes()
    if not data.startswith(SIGNATURE) or not data.endswith(TRAILER):
        print(f"error: {output} is not a complete PNG", file=sys.stderr)
        return 1
    if size(data) != size(committed.read_bytes()):
        print(f"error: {output} is {size(data)}, the committed screenshot is "
              f"{size(committed.read_bytes())}", file=sys.stderr)
        return 1
    # A blank editor decompresses to a handful of distinct byte values; a real
    # one, with artwork, text and controls, uses nearly all of them.
    pixels = zlib.decompress(b"".join(body for kind, body in chunks(data) if kind == b"IDAT"))
    if len(set(pixels)) < 64:
        print(f"error: {output} looks blank", file=sys.stderr)
        return 1
    print(f"Rendered {output} ({size(data)[0]}x{size(data)[1]}, {len(data)} bytes)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
