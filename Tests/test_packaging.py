"""Exercise release naming, archive contents and final-artifact README publishing."""

import importlib.util
import os
from pathlib import Path
import subprocess
import sys
import tarfile
import tempfile
import unittest
from unittest.mock import patch
import wave
from zipfile import ZipFile


ROOT = Path(__file__).resolve().parents[1]
SCRIPTS = ROOT / "scripts"
IR_LICENSES = ("ECHOTHIEF-LICENSE.pdf", "VOXENGO-IMPULSES-LICENSE.txt")
sys.path.insert(0, str(SCRIPTS))
from release_metadata import build_number, package_prefix

spec = importlib.util.spec_from_file_location("windows_package", SCRIPTS / "package-windows.py")
windows = importlib.util.module_from_spec(spec)
spec.loader.exec_module(windows)


def run(*args, cwd=None, env=None, check=True):
    return subprocess.run(args, cwd=cwd, env=env, check=check,
                          text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)


class PackagingTests(unittest.TestCase):
    def test_build_labels_and_invalid_input(self):
        with patch.dict(os.environ, {}, clear=True):
            self.assertEqual(build_number(), "local")
            self.assertEqual(package_prefix("0.3.2"), "Taikor-0.3.2-build-local")
        for value in ("1", "42", "100001"):
            with patch.dict(os.environ, {"GITHUB_RUN_NUMBER": value}):
                self.assertEqual(package_prefix("0.3.2"), f"Taikor-0.3.2-build-{value}")
        for value in ("", "0", "../123", "42/filename", " 42", "42\n"):
            with patch.dict(os.environ, {"GITHUB_RUN_NUMBER": value}):
                with self.assertRaises(ValueError):
                    package_prefix("0.3.2")

    def test_windows_archive_names_contents_and_stale_cleanup(self):
        with tempfile.TemporaryDirectory() as directory:
            build = Path(directory)
            (build / "CMakeCache.txt").write_text("CMAKE_PROJECT_VERSION:STATIC=0.3.2\n")
            artefacts = build / "Taikor_artefacts/Release"
            names = ("VST3/Taikor.vst3/Contents/x86_64-win/Taikor.vst3",
                     "VST3/Taikor.vst3/Contents/Resources/moduleinfo.json",
                     "CLAP/Taikor.clap", "Standalone/Taikor.exe")
            for name in names:
                path = artefacts / name
                path.parent.mkdir(parents=True, exist_ok=True)
                path.write_bytes(b"fixture binary")
            (build / "dist").mkdir()
            (build / "dist/Taikor-0.3.2-build-41-Windows-x64.zip").write_bytes(b"stale")
            with patch.dict(os.environ, {"GITHUB_RUN_NUMBER": "42"}):
                archive = windows.package(build, "Release")
            self.assertEqual(archive.name, "Taikor-0.3.2-build-42-Windows-x64.zip")
            self.assertEqual(list((build / "dist").iterdir()), [archive])
            with ZipFile(archive) as contents:
                for name in (*names, "LICENSE", "THIRD_PARTY_NOTICES.md"):
                    self.assertGreater(contents.getinfo(name).file_size, 0)
                for name in IR_LICENSES:
                    self.assertEqual(contents.read(name), (ROOT / "ThirdParty" / name).read_bytes())

    def test_linux_archive_uses_same_naming_and_contains_binaries(self):
        with tempfile.TemporaryDirectory() as directory:
            build = Path(directory)
            (build / "CMakeCache.txt").write_text("CMAKE_PROJECT_VERSION:STATIC=0.3.2\n")
            artefacts = build / "Taikor_artefacts/Release"
            names = ("VST3/Taikor.vst3/Contents/x86_64-linux/Taikor.so",
                     "CLAP/Taikor.clap", "Standalone/Taikor")
            for name in names:
                path = artefacts / name
                path.parent.mkdir(parents=True, exist_ok=True)
                path.write_bytes(b"fixture binary")
                path.chmod(0o755)
            environment = dict(os.environ, GITHUB_RUN_NUMBER="42", BUILD_DIR=str(build))
            run("bash", str(SCRIPTS / "package-linux.sh"), env=environment)
            archive = build / "dist/Taikor-0.3.2-build-42-Linux-x64.tar.gz"
            with tarfile.open(archive) as contents:
                for name in (*names, "LICENSE", "THIRD_PARTY_NOTICES.md"):
                    self.assertGreater(contents.getmember(name).size, 0)
                for name in IR_LICENSES:
                    self.assertEqual(contents.extractfile(name).read(),
                                     (ROOT / "ThirdParty" / name).read_bytes())

    def test_macos_license_list_includes_unaltered_ir_licenses(self):
        # Exercise the actual declaration without signing or packaging binaries.
        script = (SCRIPTS / "sign-and-package-macos.sh").read_text()
        declaration = script.split("THIRD_PARTY_LICENSES=(", 1)[1].split("\n)", 1)[0]
        command = 'THIRD_PARTY_LICENSES=(' + declaration + '\n)\nprintf "%s\\n" "${THIRD_PARTY_LICENSES[@]}"'
        result = run("bash", "-c", command,
                     env=dict(os.environ, PROJECT_DIR=str(ROOT)))
        paths = {Path(line) for line in result.stdout.splitlines()}
        for name in IR_LICENSES:
            source = ROOT / "ThirdParty" / name
            self.assertIn(source, paths)
            self.assertGreater(source.stat().st_size, 0)

    def test_cross_platform_manifest_rejects_different_build(self):
        workflow = (ROOT / ".github/workflows/nightly.yml").read_text()
        section = workflow.split("      - name: Check the four-file manifest\n", 1)[1]
        body = section.split("        run: |\n", 1)[1].split("\n      - name:", 1)[0]
        script = "\n".join(line[10:] for line in body.splitlines())
        with tempfile.TemporaryDirectory() as directory:
            dist = Path(directory) / "dist"
            dist.mkdir()
            for suffix in ("macOS-universal.zip", "macOS-universal.pkg",
                           "Linux-x64.tar.gz", "Windows-x64.zip"):
                (dist / f"Taikor-0.3.2-build-42-{suffix}").write_bytes(b"package")
            environment = dict(os.environ, GITHUB_RUN_NUMBER="42")
            run("bash", "-c", script, cwd=directory, env=environment)
            correct = dist / "Taikor-0.3.2-build-42-Windows-x64.zip"
            correct.rename(dist / "Taikor-0.3.2-build-41-Windows-x64.zip")
            self.assertNotEqual(run("bash", "-c", script, cwd=directory,
                                   env=environment, check=False).returncode, 0)

    def test_generated_readme_links_final_uploaded_artifact_with_build_provenance(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            repo, remote, media = root / "repo", root / "remote.git", root / "media"
            repo.mkdir()
            run("git", "init", "-b", "main", str(repo))
            run("git", "init", "--bare", str(remote))
            run("git", "config", "user.name", "Packaging test", cwd=repo)
            run("git", "config", "user.email", "packaging@example.invalid", cwd=repo)
            original = ("# Fixture\n<!-- distribution-link-begin -->\nOld link\n"
                        "<!-- distribution-link-end -->\n"
                        "<!-- peaks-table-begin: fixture -->\nunchanged peaks\n"
                        "<!-- peaks-table-end -->\nFooter\n")
            (repo / "README.md").write_text(original)
            run("git", "add", "README.md", cwd=repo)
            run("git", "commit", "-m", "fixture", cwd=repo)
            source = run("git", "rev-parse", "HEAD", cwd=repo).stdout.strip()
            run("git", "remote", "add", "origin", str(remote), cwd=repo)
            run("git", "push", "origin", "main", cwd=repo)
            (media / "Docs/audio").mkdir(parents=True)
            (media / "Docs/screenshots").mkdir()
            (media / "README.md").write_text(original)
            (media / "Docs/screenshots/taikor-standalone.png").write_bytes(
                (ROOT / "Docs/screenshots/taikor-standalone.png").read_bytes())
            for number in range(1, 28):
                with wave.open(str(media / f"Docs/audio/{number:02}-fixture.wav"), "wb") as audio:
                    audio.setparams((2, 2, 48000, 0, "NONE", "not compressed"))
                    audio.writeframes(b"\0\0\0\0")
            url = "https://github.com/example/taikor/actions/runs/123456/artifacts/789012"
            environment = dict(os.environ, GITHUB_RUN_NUMBER="42", GITHUB_REPOSITORY="example/taikor")
            # A failed metadata validation must not publish the generated assets.
            rejected = run("bash", str(SCRIPTS / "commit-generated-assets.sh"), source,
                           str(media), "origin", url, cwd=repo,
                           env=dict(environment, GITHUB_RUN_NUMBER="../42"), check=False)
            self.assertNotEqual(rejected.returncode, 0)
            self.assertEqual(run("git", "rev-parse", "HEAD", cwd=repo).stdout.strip(), source)
            run("bash", str(SCRIPTS / "commit-generated-assets.sh"), source,
                str(media), "origin", url, cwd=repo, env=environment)
            published = run("git", "--git-dir", str(remote), "show", "main:README.md").stdout
            self.assertIn(f"[Download latest distribution]({url})", published)
            self.assertIn("build 42, built from", published)
            self.assertIn(f"/commit/{source}", published)
            self.assertTrue(published.endswith(original.split("<!-- distribution-link-end -->", 1)[1]))


if __name__ == "__main__":
    unittest.main()
