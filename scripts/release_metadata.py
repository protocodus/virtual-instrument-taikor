#!/usr/bin/env python3
"""Package compiled identities into immutable version/build/platform directories."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import plistlib
import re
import shutil
import subprocess
import tarfile
import tempfile
import zipfile

ROOT = Path(__file__).resolve().parents[1]
PRODUCT = 'Taikor'
BUILD_PATTERN = r'[1-9][0-9]*(?:\.[1-9][0-9]*)?'


def build_number():
    value = os.environ.get('BUILD_NUMBER')
    if value is None:
        run = os.environ.get('GITHUB_RUN_ID')
        value = f"{run}.{os.environ.get('GITHUB_RUN_ATTEMPT', '1')}" if run else '1'
    if not re.fullmatch(BUILD_PATTERN, value):
        raise ValueError('Build number must be a positive integer or positive run.attempt pair')
    return value


def package_prefix(version, number=None):
    if not re.fullmatch(r'[0-9]+\.[0-9]+\.[0-9]+', version):
        raise ValueError('Expected a three-part semantic version')
    return f'{PRODUCT}-{version}-build.{number or build_number()}'


def load_identity(build_dir):
    identity = json.loads((Path(build_dir) / 'build-release.json').read_text())
    if identity['product'] != PRODUCT or not re.fullmatch(BUILD_PATTERN, identity['build_number']):
        raise ValueError('Invalid compiled product/build identity')
    package_prefix(identity['version'], identity['build_number'])
    if 'BUILD_NUMBER' in os.environ and build_number() != identity['build_number']:
        raise ValueError('BUILD_NUMBER disagrees with the configured build; rebuild before packaging')
    if 'GITHUB_RUN_ID' in os.environ and 'BUILD_NUMBER' not in os.environ and build_number() != identity['build_number']:
        raise ValueError('CI run identity disagrees with the configured build; rebuild before packaging')
    if os.environ.get('VERSION') and os.environ['VERSION'] != identity['version']:
        raise ValueError('VERSION disagrees with the configured build')
    return identity


def stable_dist_root(project_dir=ROOT):
    if os.environ.get('DIST_ROOT'):
        return Path(os.environ['DIST_ROOT']).expanduser().resolve()
    try:
        top = Path(subprocess.check_output(['git', '-C', str(project_dir), 'rev-parse', '--show-toplevel'], text=True, stderr=subprocess.DEVNULL).strip()).resolve()
        common = Path(subprocess.check_output(['git', '-C', str(project_dir), 'rev-parse', '--path-format=absolute', '--git-common-dir'], text=True, stderr=subprocess.DEVNULL).strip()).resolve()
        return common.parent / project_dir.resolve().relative_to(top) / 'dist'
    except (subprocess.CalledProcessError, FileNotFoundError, ValueError):
        return project_dir.resolve() / 'dist'


def distribution(identity, platform, arch):
    if platform not in ('macOS', 'Windows', 'Linux') or arch not in ('universal', 'arm64', 'x64'):
        raise ValueError('Unsupported platform or architecture')
    return stable_dist_root() / identity['version'] / f"build.{identity['build_number']}" / f'{platform}-{arch}'


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def binary_architectures(data):
    magic = data[:4]
    if magic in (b'\xca\xfe\xba\xbe', b'\xca\xfe\xba\xbf'):
        count = int.from_bytes(data[4:8], 'big'); stride = 32 if magic[-1] == 191 else 20
        cpus = [int.from_bytes(data[8+i*stride:12+i*stride], 'big') for i in range(count)]
    elif magic in (b'\xcf\xfa\xed\xfe', b'\xce\xfa\xed\xfe'):
        cpus = [int.from_bytes(data[4:8], 'little')]
    elif magic == b'\x7fELF':
        return {62: ['x64'], 183: ['arm64']}.get(int.from_bytes(data[18:20], 'little' if data[5] == 1 else 'big'), [])
    elif magic[:2] == b'MZ':
        pe = int.from_bytes(data[60:64], 'little')
        return {0x8664: ['x64'], 0xaa64: ['arm64']}.get(int.from_bytes(data[pe+4:pe+6], 'little'), [])
    else:
        return []
    return sorted({{0x1000007: 'x64', 0x100000c: 'arm64'}.get(c, 'unknown') for c in cpus})


def binary_slices(data):
    if data[:4] not in (b'\xca\xfe\xba\xbe', b'\xca\xfe\xba\xbf'):
        return [data]
    wide = data[3] == 191
    stride, width = (32, 8) if wide else (20, 4)
    count = int.from_bytes(data[4:8], 'big')
    if count < 1 or 8 + count * stride > len(data):
        raise ValueError('Invalid universal binary header')
    slices = []
    for index in range(count):
        start = 8 + index * stride + 8
        offset = int.from_bytes(data[start:start+width], 'big')
        size = int.from_bytes(data[start+width:start+2*width], 'big')
        if size == 0 or offset < 8 + count * stride or offset + size > len(data):
            raise ValueError('Invalid universal binary slice')
        slices.append(data[offset:offset+size])
    return slices


def require_empty_destination(dest):
    if dest.exists() and (not dest.is_dir() or any(dest.iterdir())):
        raise ValueError(f'Distribution already exists; choose a new build number: {dest}')


def verify_binaries(build_dir, config, platform, arch, identity=None):
    identity = identity or load_identity(build_dir)
    artefacts = Path(build_dir) / f'{PRODUCT}_artefacts' / config
    formats = ['VST3', 'Standalone'] + (['AU'] if platform == 'macOS' else [])
    if PRODUCT != 'Septum': formats.append('CLAP')
    expected_archs = {'arm64', 'x64'} if arch == 'universal' else {arch}
    records = []
    for fmt in formats:
        suffix = {'VST3': '.vst3', 'AU': '.component', 'CLAP': '.clap', 'Standalone': '.app' if platform == 'macOS' else '.exe' if platform == 'Windows' else ''}[fmt]
        bundle = artefacts / fmt / (PRODUCT + suffix)
        if platform == 'macOS':
            binary = bundle / 'Contents/MacOS' / PRODUCT
            plist = plistlib.loads((bundle / 'Contents/Info.plist').read_bytes())
            if (plist.get('CFBundleShortVersionString'), plist.get('CFBundleVersion')) != (identity['version'], identity['build_number']):
                raise ValueError(f'Bundle identity disagrees with configured identity: {bundle}')
        elif fmt == 'VST3':
            binary = bundle / 'Contents' / ('x86_64-win' if platform == 'Windows' else 'x86_64-linux') / (PRODUCT + ('.vst3' if platform == 'Windows' else '.so'))
        else:
            binary = bundle
        data = binary.read_bytes()
        if any(identity['binary_identity'].encode('ascii') not in part for part in binary_slices(data)):
            raise ValueError(f'Binary has a stale or missing compiled build identity: {binary}')
        actual_archs = binary_architectures(data)
        if set(actual_archs) != expected_archs:
            raise ValueError(f'Unexpected binary architectures {actual_archs}: {binary}')
        records.append({'format': fmt, 'bundle': str(bundle), 'binary': str(binary), 'binary_sha256': hashlib.sha256(data).hexdigest(), 'architectures': actual_archs})
    return records


def finalize(build_dir, config, platform, arch, files):
    identity = load_identity(build_dir)
    binaries = verify_binaries(build_dir, config, platform, arch, identity)
    dest = distribution(identity, platform, arch)
    stem = package_prefix(identity['version'], identity['build_number']) + f'-{platform}-{arch}'
    packages = []
    for item in files:
        path = Path(item).resolve()
        if path.parent != dest.resolve() or not path.name.startswith(stem + '.') or not path.is_file() or path.stat().st_size == 0:
            raise ValueError(f'Unexpected or missing distribution: {path}')
        if path.suffix == '.zip':
            with zipfile.ZipFile(path) as archive:
                if archive.testzip(): raise ValueError(f'Corrupt ZIP: {path}')
        elif path.name.endswith('.tar.gz'):
            with tarfile.open(path, 'r:gz') as archive:
                if not archive.getmembers(): raise ValueError(f'Empty tarball: {path}')
        packages.append({'path': path.name, 'sha256': digest(path), 'bytes': path.stat().st_size})
    manifest = dest / (stem + '.manifest.json')
    manifest.write_text(json.dumps({**identity, 'platform': platform, 'architecture': arch, 'configuration': config, 'binaries': binaries, 'packages': packages}, indent=2) + '\n')
    sums = [f"{p['sha256']}  {p['path']}" for p in packages]
    sums.append(f'{digest(manifest)}  {manifest.name}')
    (dest / (stem + '.sha256')).write_text('\n'.join(sums) + '\n')
    return dest


def package(build_dir, config, platform, arch='x64'):
    identity = load_identity(build_dir)
    binaries = verify_binaries(build_dir, config, platform, arch, identity)
    dest = distribution(identity, platform, arch)
    require_empty_destination(dest)
    dest.mkdir(parents=True, exist_ok=True)
    stem = package_prefix(identity['version'], identity['build_number']) + f'-{platform}-{arch}'
    extension = '.zip' if platform == 'Windows' else '.tar.gz'
    output = dest / (stem + extension)
    for required in (ROOT/'LICENSE', ROOT/'THIRD_PARTY_NOTICES.md', ROOT/'ThirdParty/JUCE-LICENSE.md'):
        if not required.is_file() or required.stat().st_size == 0:
            raise FileNotFoundError(f'Missing or empty distribution notice: {required}')
    with tempfile.TemporaryDirectory(prefix=f'{PRODUCT}-package-') as folder:
        staging = Path(folder)
        for record in binaries:
            source = Path(record['bundle']); target = staging / record['format'] / source.name
            target.parent.mkdir(parents=True, exist_ok=True)
            if source.is_dir(): shutil.copytree(source, target)
            else: shutil.copy2(source, target)
        for notice in [ROOT/'LICENSE', ROOT/'THIRD_PARTY_NOTICES.md', *sorted((ROOT/'ThirdParty').glob('*'))]:
            if notice.is_file() and (notice.name == 'LICENSE' or notice.suffix.lower() in ('.md', '.txt', '.pdf')):
                target = staging / notice.relative_to(ROOT); target.parent.mkdir(parents=True, exist_ok=True); shutil.copy2(notice, target)
        (staging/'build-release.json').write_text(json.dumps(identity, indent=2)+'\n')
        temporary = dest / (stem + '.pending')
        try:
            if platform == 'Windows':
                with zipfile.ZipFile(temporary, 'w', compression=zipfile.ZIP_DEFLATED) as archive:
                    for path in sorted(staging.rglob('*')):
                        if path.is_file(): archive.write(path, path.relative_to(staging))
            else:
                with tarfile.open(temporary, 'w:gz') as archive:
                    for path in sorted(staging.iterdir()): archive.add(path, arcname=path.name)
            temporary.replace(output)
        finally:
            temporary.unlink(missing_ok=True)
    finalize(build_dir, config, platform, arch, [output])
    return output


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--version')
    parser.add_argument('--build-dir', type=Path)
    parser.add_argument('--config', default='Release')
    parser.add_argument('--field')
    parser.add_argument('--platform', choices=['macOS', 'Windows', 'Linux'])
    parser.add_argument('--arch', default='x64')
    parser.add_argument('--dist', action='store_true')
    parser.add_argument('--verify', action='store_true')
    parser.add_argument('--require-empty', action='store_true')
    parser.add_argument('--package', action='store_true')
    parser.add_argument('--finalize', nargs='+')
    args = parser.parse_args()
    try:
        if args.build_dir is None:
            print(package_prefix(args.version) if args.version else build_number()); return
        identity = load_identity(args.build_dir)
        if args.require_empty: require_empty_destination(distribution(identity, args.platform, args.arch))
        if args.field: print(identity[args.field])
        elif args.dist: print(distribution(identity, args.platform, args.arch))
        elif args.verify: print(json.dumps(verify_binaries(args.build_dir, args.config, args.platform, args.arch)))
        elif args.package: print(package(args.build_dir, args.config, args.platform, args.arch))
        elif args.finalize: print(finalize(args.build_dir, args.config, args.platform, args.arch, args.finalize))
        else: print(json.dumps(identity, indent=2))
    except (ValueError, KeyError, OSError, json.JSONDecodeError) as error:
        parser.exit(1, f'error: {error}\n')


if __name__ == '__main__': main()
