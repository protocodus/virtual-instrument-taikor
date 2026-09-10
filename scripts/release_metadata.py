#!/usr/bin/env python3
"""Shared release filenames: CI run number, or an explicit local-build label."""

import argparse
import os
import re


def build_number() -> str:
    value = os.environ.get("GITHUB_RUN_NUMBER")
    if value is None:
        return "local"
    if not re.fullmatch(r"[1-9][0-9]*", value):
        raise ValueError("GITHUB_RUN_NUMBER must be a positive integer when set")
    return value


def package_prefix(version: str) -> str:
    if not re.fullmatch(r"[0-9]+(?:\.[0-9]+){1,3}", version):
        raise ValueError("Expected a valid project version")
    return f"Taikor-{version}-build-{build_number()}"


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--version", help="Print the package prefix for this version")
    args = parser.parse_args()
    try:
        print(package_prefix(args.version) if args.version else build_number())
    except ValueError as error:
        parser.exit(1, f"error: {error}\n")
