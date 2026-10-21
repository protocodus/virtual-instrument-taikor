"""Exercise release naming, archive contents and final-artifact README publishing."""

import importlib.util
import hashlib
import json
import os
from pathlib import Path
import subprocess
import struct
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


def suite_fixture(repo, media):
    """Small source fixtures plus a sparse, correctly sized PCM24 delivery."""
    for name in ("Source/DSP/TaikoEngine.cpp", "Source/DSP/TaikoEngine.h",
                 "Tools/RenderSuite.cpp", "Tools/ProduceSuite.py",
                 "Assets/ImpulseResponses/theater.wav"):
        path = repo / name
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(f"fixture source: {name}\n".encode())
    suite = media / "Docs/audio/suite"
    suite.mkdir(parents=True)
    path = suite / "taikor-suite-96k24.wav"
    # FFmpeg uses WAVE_FORMAT_EXTENSIBLE for this delivery. Include an unknown
    # odd-sized chunk to exercise RIFF chunk padding as well as the PCM subtype.
    fmt = struct.pack("<HHIIHHHHI", 65534, 2, 96000, 576000, 6, 24, 22, 24, 3)
    fmt += bytes.fromhex("0100000000001000800000aa00389b71")
    data_bytes = 96000 * 60 * 2 * 3
    chunks = b"JUNK" + struct.pack("<I", 1) + b"x\0"
    chunks += b"fmt " + struct.pack("<I", len(fmt)) + fmt
    chunks += b"data" + struct.pack("<I", data_bytes)
    with path.open("wb") as audio:
        audio.write(b"RIFF" + struct.pack("<I", 4 + len(chunks) + data_bytes) + b"WAVE" + chunks)
        audio.seek(data_bytes - 1, 1)
        audio.write(b"\0")
    (suite / "taikor-suite.mp3").write_bytes(b"ID3 fixture MP3 validated by production QA")
    (suite / "taikor-suite.mid").write_bytes(b"MThd\0\0\0\6\0\1\0\5\7\x80")
    (suite / "taikor-suite-score.csv").write_text("seconds,instance,event\n0,1,note\n")

    def digest(path):
        return hashlib.sha256(path.read_bytes()).hexdigest()

    qa = {
        "schema_version": 1, "status": "passed", "failures": [],
        "source_files": {str(path.relative_to(repo)): digest(path)
                         for path in sorted((repo / "Source/DSP").glob("*"))},
        "production_script": {"sha256": digest(repo / "Tools/ProduceSuite.py")},
        "processing": {"engine_sample_rate_hz": 384000, "master_sample_rate_hz": 384000,
                       "room_ir_sha256": digest(repo / "Assets/ImpulseResponses/theater.wav")},
        "source": {"sample_rate_hz": 384000, "frames": 23040000},
        "artifacts": {}, "performance_files": {},
    }
    qa["source_files"]["Tools/RenderSuite.cpp"] = digest(repo / "Tools/RenderSuite.cpp")
    for name, rate in (("taikor-suite-96k24.wav", 96000), ("taikor-suite.mp3", 48000)):
        path = suite / name
        qa["artifacts"][name] = {
            "sha256": digest(path), "bytes": path.stat().st_size,
            "sample_rate_hz": rate, "channels": 2, "frames": 60 * rate,
            "duration_seconds": 60.0, "clipped_sample_count": 0, "nonfinite_sample_count": 0,
            "rms_dbfs": -20.0, "estimated_true_peak_dbtp": -1.2,
        }
    for name in ("taikor-suite.mid", "taikor-suite-score.csv"):
        path = suite / name
        qa["performance_files"][name] = {"sha256": digest(path), "bytes": path.stat().st_size}
    (suite / "taikor-suite-qa.json").write_text(json.dumps(qa))
    return suite, qa


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
                     "VST3/Taikor.vst3/Contents/Resources/moduleinfo.json",
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

    def test_linux_package_rejects_ambiguous_project_version(self):
        with tempfile.TemporaryDirectory() as directory:
            build = Path(directory)
            (build / "CMakeCache.txt").write_text(
                "CMAKE_PROJECT_VERSION:STATIC=0.3.2\n"
                "CMAKE_PROJECT_VERSION:STATIC=0.3.3\n")
            result = run("bash", str(SCRIPTS / "package-linux.sh"),
                         env=dict(os.environ, BUILD_DIR=str(build)), check=False)
            self.assertNotEqual(result.returncode, 0)
            self.assertIn("exactly one CMake project version", result.stderr)

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

    def test_nightly_stages_only_the_five_suite_deliveries(self):
        workflow = (ROOT / ".github/workflows/nightly.yml").read_text()
        section = workflow.split("      - name: Produce and stage the one-minute suite\n", 1)[1]
        body = section.split("        run: |\n", 1)[1].split("\n      - name:", 1)[0]
        script = "\n".join(line[10:] for line in body.splitlines())
        names = {"taikor-suite-96k24.wav", "taikor-suite.mp3", "taikor-suite.mid",
                 "taikor-suite-score.csv", "taikor-suite-qa.json"}
        masters = {"taikor-suite-dry.wav", "taikor-suite-master-384k-float.wav"}
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            (root / "Tools").mkdir()
            # Exercise the workflow's actual shell staging with a lightweight
            # producer standing in for the expensive physical-model render.
            producer = root / "Tools/ProduceSuite.py"
            producer.write_text(
                "from pathlib import Path\nimport sys\n"
                "output = Path(sys.argv[sys.argv.index('--output') + 1])\n"
                "output.mkdir(parents=True, exist_ok=True)\n"
                f"for name in {sorted(names | masters)!r}:\n"
                "    (output / name).write_bytes(b'fresh output')\n")
            run("bash", "-c", script, cwd=root)
            suite = root / "build-dsp/generated-media/Docs/audio/suite"
            self.assertEqual({path.name for path in suite.iterdir()}, names)
            for name in masters:
                self.assertTrue((root / "build-dsp/suite-master" / name).is_file())
            # A failed renderer must not allow any stale result to be staged.
            producer.write_text("raise SystemExit(1)\n")
            for path in suite.iterdir():
                path.unlink()
            self.assertNotEqual(run("bash", "-c", script, cwd=root, check=False).returncode, 0)
            self.assertEqual(list(suite.iterdir()), [])

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
            suite, original_qa = suite_fixture(repo, media)
            preserved = ("Docs/audio/manual.wav", "Docs/audio/historical/preview.wav",
                         "Docs/audio/suite/listening-notes.txt")
            for name in preserved:
                path = repo / name
                path.parent.mkdir(parents=True, exist_ok=True)
                path.write_bytes(b"keep this manually authored file")
            obsolete = repo / "Docs/audio/00-obsolete.wav"
            obsolete.write_bytes(b"old numbered diagnostic")
            run("git", "add", ".", cwd=repo)
            run("git", "commit", "-m", "fixture", cwd=repo)
            source = run("git", "rev-parse", "HEAD", cwd=repo).stdout.strip()
            run("git", "remote", "add", "origin", str(remote), cwd=repo)
            run("git", "push", "origin", "main", cwd=repo)
            (media / "Docs/audio").mkdir(parents=True, exist_ok=True)
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

            def publish(check=True):
                return run("bash", str(SCRIPTS / "commit-generated-assets.sh"), source,
                           str(media), "origin", url, cwd=repo, env=environment, check=check)

            # Provenance failures must be caught before the checkout or remote
            # changes, even when the source-SHA publication guard would pass.
            mutations = {
                "failed QA": lambda qa: qa.update(status="failed"),
                "reported failure": lambda qa: qa.update(failures=["clipping"]),
                "different DSP": lambda qa: qa["source_files"].update({"Source/DSP/TaikoEngine.cpp": "0" * 64}),
                "missing DSP": lambda qa: qa["source_files"].pop("Source/DSP/TaikoEngine.h"),
                "different producer": lambda qa: qa["production_script"].update(sha256="0" * 64),
                "different room": lambda qa: qa["processing"].update(room_ir_sha256="0" * 64),
                "draft rate": lambda qa: qa["processing"].update(engine_sample_rate_hz=48000),
                "corrupt WAV": lambda qa: qa["artifacts"]["taikor-suite-96k24.wav"].update(sha256="0" * 64),
                "short playback": lambda qa: qa["artifacts"]["taikor-suite.mp3"].update(frames=48),
                "nonfinite peak": lambda qa: qa["artifacts"]["taikor-suite.mp3"].update(estimated_true_peak_dbtp=float("nan")),
                "clipping": lambda qa: qa["artifacts"]["taikor-suite.mp3"].update(estimated_true_peak_dbtp=0.1),
                "different score": lambda qa: qa["performance_files"]["taikor-suite-score.csv"].update(sha256="0" * 64),
            }
            qa_path = suite / "taikor-suite-qa.json"
            for label, mutate in mutations.items():
                with self.subTest(suite_validation=label):
                    qa = json.loads(json.dumps(original_qa))
                    mutate(qa)
                    qa_path.write_text(json.dumps(qa))
                    result = publish(check=False)
                    self.assertNotEqual(result.returncode, 0, result.stdout)
                    self.assertIn("error: suite", result.stderr)
                    self.assertEqual(run("git", "rev-parse", "HEAD", cwd=repo).stdout.strip(), source)
                    self.assertEqual(run("git", "status", "--porcelain", cwd=repo).stdout, "")
            qa_path.write_text(json.dumps(original_qa))
            master = suite / "taikor-suite-dry.wav"
            master.write_bytes(b"large working masters do not belong in publication")
            result = publish(check=False)
            self.assertNotEqual(result.returncode, 0)
            self.assertIn("exactly the five suite delivery files", result.stderr)
            master.unlink()
            score = suite / "taikor-suite-score.csv"
            score_bytes = score.read_bytes()
            score.unlink()
            self.assertNotEqual(publish(check=False).returncode, 0)
            score.write_bytes(score_bytes)

            # A forged format claim in QA cannot make a malformed WAV publish.
            wav = suite / "taikor-suite-96k24.wav"
            with wav.open("r+b") as audio:
                original_header = audio.read(4)
                audio.seek(0)
                audio.write(b"NOPE")
            qa = json.loads(json.dumps(original_qa))
            qa["artifacts"][wav.name]["sha256"] = hashlib.sha256(wav.read_bytes()).hexdigest()
            qa_path.write_text(json.dumps(qa))
            result = publish(check=False)
            self.assertNotEqual(result.returncode, 0)
            self.assertIn("invalid RIFF header", result.stderr)
            with wav.open("r+b") as audio:
                audio.write(original_header)
            qa_path.write_text(json.dumps(original_qa))

            publish()
            published = run("git", "--git-dir", str(remote), "show", "main:README.md").stdout
            self.assertIn(f"[Download latest distribution]({url})", published)
            self.assertIn("build 42, built from", published)
            self.assertIn(f"/commit/{source}", published)
            self.assertTrue(published.endswith(original.split("<!-- distribution-link-end -->", 1)[1]))
            for name in preserved:
                self.assertEqual((repo / name).read_bytes(), b"keep this manually authored file")
            self.assertFalse(obsolete.exists())
            self.assertEqual(len(list((repo / "Docs/audio").glob("[0-9][0-9]-*.wav"))), 27)
            for name in original_qa["artifacts"] | original_qa["performance_files"] | {"taikor-suite-qa.json": {}}:
                self.assertEqual((repo / "Docs/audio/suite" / name).read_bytes(), (suite / name).read_bytes())
                run("git", "--git-dir", str(remote), "cat-file", "-e", f"main:Docs/audio/suite/{name}")
            self.assertFalse((repo / "Docs/audio/suite/taikor-suite-dry.wav").exists())


if __name__ == "__main__":
    unittest.main()
