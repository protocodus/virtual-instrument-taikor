#!/usr/bin/env python3
"""Package a Release Windows build, failing if a plug-in or notice is missing."""

import argparse
from pathlib import Path
import re
from zipfile import ZIP_DEFLATED, ZipFile

from release_metadata import package_prefix


def package(build_dir: Path, config: str) -> Path:
    project_dir = Path(__file__).resolve().parent.parent
    artefacts = build_dir / "Taikor_artefacts" / config
    versions = re.findall(
        r"^CMAKE_PROJECT_VERSION:STATIC=(.+)$",
        (build_dir / "CMakeCache.txt").read_text(encoding="utf-8"),
        re.MULTILINE,
    )
    if len(versions) != 1 or not re.fullmatch(r"\d+(?:\.\d+){1,3}", versions[0]):
        raise ValueError("Expected one valid CMake project version")
    prefix = package_prefix(versions[0])

    required = (
        "VST3/Taikor.vst3/Contents/x86_64-win/Taikor.vst3",
        "CLAP/Taikor.clap",
        "Standalone/Taikor.exe",
    )
    notices = (
        "LICENSE",
        "THIRD_PARTY_NOTICES.md",
        "ThirdParty/JUCE-LICENSE.md",
        "ThirdParty/CLAP-JUCE-EXTENSIONS-LICENSE.md",
        "ThirdParty/CLAP-LICENSE.md",
        "ThirdParty/CLAP-HELPERS-LICENSE.md",
        "ThirdParty/ECHOTHIEF-LICENSE.pdf",
        "ThirdParty/VOXENGO-IMPULSES-LICENSE.txt",
    )
    for source in [*(artefacts / name for name in required),
                   *(project_dir / name for name in notices)]:
        if not source.is_file() or source.stat().st_size == 0:
            raise FileNotFoundError(f"Missing or empty distribution file: {source}")

    dist = build_dir / "dist"
    dist.mkdir(parents=True, exist_ok=True)
    archive_path = dist / f"{prefix}-Windows-x64.zip"
    # Wildcard-driven publication must not upload a previous build as well.
    for stale in dist.glob("Taikor-*-Windows-x64.zip"):
        stale.unlink()
    with ZipFile(archive_path, "w", compression=ZIP_DEFLATED) as archive:
        # Include the complete VST3 bundle, including its generated module info.
        for directory in ("VST3/Taikor.vst3", "CLAP", "Standalone"):
            for source in sorted((artefacts / directory).rglob("*")):
                if source.is_file():
                    archive.write(source, source.relative_to(artefacts).as_posix())
        for name in notices:
            source = project_dir / name
            archive.write(source, source.name)

    with ZipFile(archive_path) as archive:
        if archive.testzip() is not None:
            raise ValueError(f"Corrupt archive: {archive_path}")
        for name in (*required, *(Path(name).name for name in notices)):
            if archive.getinfo(name).file_size == 0:
                raise ValueError(f"Empty archive entry: {name}")

    return archive_path


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--build-dir", type=Path, default=Path("build-win"))
    parser.add_argument("--config", default="Release")
    args = parser.parse_args()
    print(f"Packaged {package(args.build_dir, args.config)}")
