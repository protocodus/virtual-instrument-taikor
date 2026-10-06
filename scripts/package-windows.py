#!/usr/bin/env python3
"""Package the configured Windows binaries without changing their identity."""
import argparse
from pathlib import Path
from release_metadata import package as package_release

def package(build_dir, config):
    return package_release(build_dir, config, 'Windows', 'x64')

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--build-dir', type=Path, default=Path('build-win'))
    parser.add_argument('--config', default='Release')
    args = parser.parse_args()
    print(package(args.build_dir, args.config))
