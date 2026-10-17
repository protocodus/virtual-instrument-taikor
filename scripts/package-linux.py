#!/usr/bin/env python3
"""Package a Release Linux build, failing if a plug-in or notice is missing."""

import argparse
from pathlib import Path
import re
import tarfile

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

    if not (artefacts / "VST3/Taikor.vst3").is_dir():
        raise FileNotFoundError(f"Missing VST3 bundle in {artefacts}")
    required = ("CLAP/Taikor.clap", "Standalone/Taikor")
    notices = (
        project_dir / "LICENSE",
        project_dir / "THIRD_PARTY_NOTICES.md",
        *sorted((project_dir / "ThirdParty").glob("*-LICENSE.md")),
        project_dir / "ThirdParty/ECHOTHIEF-LICENSE.pdf",
        project_dir / "ThirdParty/VOXENGO-IMPULSES-LICENSE.txt",
    )
    for source in [*(artefacts / name for name in required), *notices]:
        if not source.is_file() or source.stat().st_size == 0:
            raise FileNotFoundError(f"Missing or empty distribution file: {source}")

    dist = build_dir / "dist"
    dist.mkdir(parents=True, exist_ok=True)
    archive_path = dist / f"{prefix}-Linux-x64.tar.gz"
    # Wildcard-driven publication must not upload a previous build as well.
    for stale in dist.glob("Taikor-*-Linux-x64.tar.gz"):
        stale.unlink()
    with tarfile.open(archive_path, "w:gz") as archive:
        for name in ("VST3/Taikor.vst3", *required):
            archive.add(artefacts / name, name)
        for source in notices:
            archive.add(source, source.name)

    with tarfile.open(archive_path) as archive:
        for name in (*required, *(source.name for source in notices)):
            if archive.getmember(name).size == 0:
                raise ValueError(f"Empty archive entry: {name}")

    return archive_path


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--build-dir", type=Path, default=Path("build-linux"))
    parser.add_argument("--config", default="Release")
    args = parser.parse_args()
    print(f"Packaged {package(args.build_dir, args.config)}")
