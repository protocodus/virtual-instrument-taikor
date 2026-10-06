#!/usr/bin/env python3
"""Packaging regressions: compiled identity, provenance, worktrees and history."""
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import plistlib
import subprocess
import tempfile
import unittest
from unittest.mock import patch
import zipfile

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location('release_metadata_under_test', ROOT/'scripts/release_metadata.py')
metadata = importlib.util.module_from_spec(spec)
spec.loader.exec_module(metadata)


def write_fixture(build, platform='Windows', version='0.3.2', number='42'):
    build = Path(build); build.mkdir(parents=True, exist_ok=True)
    product = metadata.PRODUCT
    identity = {'product':product, 'version':version, 'build_number':number,
                'source_revision':'a'*40, 'source_dirty':False,
                'binary_identity':f'PROTOCODUS_BUILD_IDENTITY={product}|{version}|{number}|'+('a'*40)}
    (build/'build-release.json').write_text(json.dumps(identity))
    if platform == 'Windows':
        header=bytearray(128);header[:2]=b'MZ';header[60:64]=(64).to_bytes(4,'little');header[64:68]=b'PE\0\0';header[68:70]=(0x8664).to_bytes(2,'little')
    elif platform == 'Linux':
        header=bytearray(64);header[:6]=b'\x7fELF\x02\x01';header[18:20]=(62).to_bytes(2,'little')
    else:
        header=bytearray(48);header[:4]=b'\xca\xfe\xba\xbe';header[4:8]=(2).to_bytes(4,'big')
        slices=[]
        for cpu in [0x1000007,0x100000c]:
            thin=bytearray(32);thin[:4]=b'\xcf\xfa\xed\xfe';thin[4:8]=cpu.to_bytes(4,'little')
            slices.append(bytes(thin)+identity['binary_identity'].encode('ascii'))
        offset=48
        for index,(cpu,thin) in enumerate(zip([0x1000007,0x100000c],slices)):
            start=8+index*20;header[start:start+4]=cpu.to_bytes(4,'big');header[start+8:start+12]=offset.to_bytes(4,'big');header[start+12:start+16]=len(thin).to_bytes(4,'big');offset+=len(thin)
        header.extend(b''.join(slices))
    content=bytes(header)+identity['binary_identity'].encode('ascii')
    formats=['VST3','Standalone']+(['AU'] if platform=='macOS' else [])+(['CLAP'] if product!='Septum' else [])
    for fmt in formats:
        suffix={'VST3':'.vst3','AU':'.component','CLAP':'.clap','Standalone':'.app' if platform=='macOS' else '.exe' if platform=='Windows' else ''}[fmt]
        bundle=build/f'{product}_artefacts/Release'/fmt/(product+suffix)
        if platform=='macOS':
            binary=bundle/'Contents/MacOS'/product
            binary.parent.mkdir(parents=True,exist_ok=True)
            (bundle/'Contents/Info.plist').write_bytes(plistlib.dumps({'CFBundleShortVersionString':version,'CFBundleVersion':number,'CFBundlePackageType':'APPL' if fmt=='Standalone' else 'BNDL'}))
        elif fmt=='VST3':
            binary=bundle/'Contents'/('x86_64-win' if platform=='Windows' else 'x86_64-linux')/(product+('.vst3' if platform=='Windows' else '.so'))
        else: binary=bundle
        binary.parent.mkdir(parents=True,exist_ok=True);binary.write_bytes(content);binary.chmod(0o755)
    return identity


class ReleaseMetadataTests(unittest.TestCase):
    def setUp(self):
        self.folder=tempfile.TemporaryDirectory();self.addCleanup(self.folder.cleanup)
        self.temp=Path(self.folder.name).resolve()
        env=patch.dict(os.environ,{},clear=True);env.start();self.addCleanup(env.stop)

    def test_build_identity_defaults_and_ci_attempts(self):
        self.assertEqual(metadata.build_number(),'1')
        for value in ['1','100','100.2']:
            with patch.dict(os.environ,{'BUILD_NUMBER':value}): self.assertEqual(metadata.build_number(),value)
        for value in ['0','1.0','01','1.2.3','../2','','2\n']:
            with patch.dict(os.environ,{'BUILD_NUMBER':value}):
                with self.assertRaises(ValueError):metadata.build_number()
        with patch.dict(os.environ,{'GITHUB_RUN_ID':'987654','GITHUB_RUN_ATTEMPT':'3'}): self.assertEqual(metadata.build_number(),'987654.3')

    def test_packaging_cannot_relabel_configured_or_compiled_binaries(self):
        build=self.temp/'build';identity=write_fixture(build)
        with patch.dict(os.environ,{'BUILD_NUMBER':'43'}):
            with self.assertRaisesRegex(ValueError,'disagrees'):metadata.load_identity(build)
        identity['build_number']='43';identity['binary_identity']=identity['binary_identity'].replace('|42|','|43|')
        (build/'build-release.json').write_text(json.dumps(identity))
        with self.assertRaisesRegex(ValueError,'stale or missing'):metadata.verify_binaries(build,'Release','Windows','x64')

    def test_architecture_is_checked_from_binary(self):
        build=self.temp/'build';write_fixture(build)
        with self.assertRaisesRegex(ValueError,'architectures'):metadata.verify_binaries(build,'Release','Windows','arm64')

    def test_versioned_packages_preserve_older_builds_and_include_receipts(self):
        build=self.temp/'build';write_fixture(build)
        with patch.dict(os.environ,{'DIST_ROOT':str(self.temp/'dist')}):
            old=self.temp/'dist/0.3.2/build.41/Windows-x64/previous.zip';old.parent.mkdir(parents=True);old.write_bytes(b'previous')
            archive=metadata.package(build,'Release','Windows','x64')
            stem=f'{metadata.PRODUCT}-0.3.2-build.42-Windows-x64'
            self.assertEqual(archive.name,stem+'.zip');self.assertEqual(archive.parent,self.temp/'dist/0.3.2/build.42/Windows-x64')
            self.assertEqual(old.read_bytes(),b'previous')
            with zipfile.ZipFile(archive) as contents:
                self.assertEqual(contents.read('build-release.json').decode().strip(),json.dumps(metadata.load_identity(build),indent=2))
                self.assertGreater(contents.getinfo('LICENSE').file_size,0)
            manifest=json.loads((archive.parent/(stem+'.manifest.json')).read_text())
            self.assertEqual(manifest['packages'][0]['sha256'],hashlib.sha256(archive.read_bytes()).hexdigest())
            for line in (archive.parent/(stem+'.sha256')).read_text().splitlines():
                digest,path=line.split('  ',1);self.assertEqual(digest,hashlib.sha256((archive.parent/path).read_bytes()).hexdigest())

    def test_same_build_cannot_overwrite_a_previous_distribution(self):
        build=self.temp/'build';write_fixture(build)
        with patch.dict(os.environ,{'DIST_ROOT':str(self.temp/'dist')}):
            archive=metadata.package(build,'Release','Windows','x64');before=archive.read_bytes()
            with self.assertRaisesRegex(ValueError,'already exists'):metadata.package(build,'Release','Windows','x64')
            self.assertEqual(before,archive.read_bytes())

    def test_every_universal_slice_must_have_the_frozen_identity(self):
        build=self.temp/'build';write_fixture(build,'macOS')
        metadata.verify_binaries(build,'Release','macOS','universal')
        binary=build/f'{metadata.PRODUCT}_artefacts/Release/VST3/{metadata.PRODUCT}.vst3/Contents/MacOS/{metadata.PRODUCT}'
        data=binary.read_bytes();marker=metadata.load_identity(build)['binary_identity'].encode()
        binary.write_bytes(data.replace(marker,marker.replace(b'|42|',b'|41|'),1))
        with self.assertRaisesRegex(ValueError,'stale or missing'):metadata.verify_binaries(build,'Release','macOS','universal')

    def test_primary_checkout_dist_is_shared_by_nested_project_worktree(self):
        repo=self.temp/'repo';worktree=self.temp/'worktree';repo.mkdir()
        def git(*args):return subprocess.run(['git',*map(str,args)],check=True,capture_output=True,text=True)
        git('init','-b','main',repo);(repo/'fixture').write_text('fixture')
        git('-C',repo,'add','fixture');git('-C',repo,'-c','user.name=Test','-c','user.email=test@example.invalid','commit','-m','fixture')
        git('-C',repo,'worktree','add','--detach',worktree)
        project=worktree/'instrument';project.mkdir()
        self.assertEqual(metadata.stable_dist_root(project),repo/'instrument/dist')
        with patch.dict(os.environ,{'DIST_ROOT':str(self.temp/'override')}):self.assertEqual(metadata.stable_dist_root(project),self.temp/'override')

    def test_no_git_falls_back_to_project_dist(self):
        self.assertEqual(metadata.stable_dist_root(self.temp),self.temp/'dist')

if __name__=='__main__':unittest.main()
