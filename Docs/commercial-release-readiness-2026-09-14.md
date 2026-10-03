# Commercial release readiness

Status: engineering validation is in progress. This document is a release
checklist, not a legal opinion or a claim that the product is commercially
cleared.

## Current build and package facts

- Project version is `1.0.0` in `CMakeLists.txt`; package names add the CI run
  number or `build-local` for local packaging.
- Supported formats are macOS VST3, Audio Unit, CLAP and standalone; Windows
  x64 VST3, CLAP and standalone; and Linux x64 VST3, CLAP and standalone.
- The normal macOS script builds arm64 + x86_64 when `BUILD_UNIVERSAL=ON`.
  Current CI packaging uses ad-hoc bundle signing, creates an unsigned PKG
  unless an installer identity is supplied, and notarizes only when
  `NOTARY_PROFILE` is supplied with both Developer ID identities.
- The offline ensemble pool is intentionally offline-only. Live Apple audio
  helper-thread integration is a separate host concern: any future live helper
  thread must join the host-provided `AudioWorkgroup`. The generic pool must be
  exercised through `prepareOfflineRendering()` and
  `getOfflineWorkerCount()` rather than used by live audio callbacks.
- No CPU improvement, benchmark result or release-performance target is
  asserted here. Existing comparison/reference locations are
  `build-cpu-sep12/native` and `build-release-sep14/reference`; record measured
  results only after the central build run supplies them.

## Blocking gates before a commercial release

- **Third-party asset rights:** obtain and archive written permission or an
  approved distribution interpretation for the two EchoThief impulse responses.
  Review and satisfy every Voxengo redistribution condition for the opera
  response. Confirm that the exact embedded WAV bytes and required notices are
  shipped together. `THIRD_PARTY_NOTICES.md` and `ThirdParty/` document
  provenance and supplied terms; they do not certify legal clearance.
- **JUCE licensing:** choose and document either a commercial JUCE licence or
  complete AGPLv3 compliance for the combined distributed work. The vendored
  notice is not proof that either path has been completed.
- **Apple signing and notarization:** use the intended Developer ID Application
  and Installer identities, enable hardened runtime where required, submit the
  final PKG to Apple notarization, staple it, and verify the final ZIP/PKG and
  every independently installable bundle after signing. Replace the current
  ad-hoc/unsigned CI path for the release build.
- **Host matrix:** test installation, scanning, MIDI note/CC behavior, state
  save/restore, audio layout and automation in representative current hosts on
  macOS 11+ arm64 and x86_64, Windows x64 and Linux x64. The build matrix proves
  artifacts can compile and package; it does not prove host compatibility.
  Run the official `pluginval` 1.0.4 release against the rebuilt VST3 on each
  applicable platform and run macOS `auval` against the rebuilt Audio Unit.
  Keep those tools and test outputs isolated from installed plug-in locations;
  both results remain pending until the central validation run supplies them.
- **Release artifact review:** verify that all four final assets share one
  version and CI build number, contain the required plug-in binaries and notice
  files, have matching SHA-256 records, and come from the intended source
  commit. Do not publish a partial platform set as the commercial release.

## Reproducible validation commands

JUCE-free DSP build, tests and offline benchmark smoke test:

```bash
cmake -S . -B build-dsp -G Ninja -DCMAKE_BUILD_TYPE=Release \
  -DTAIKOR_BUILD_PLUGIN=OFF -DBUILD_TESTING=ON
cmake --build build-dsp --parallel
ctest --test-dir build-dsp --output-on-failure
./build-dsp/TaikorBenchmarkEnsembleParallel
```

ASan + UBSan validation:

```bash
cmake -S . -B build-sanitized -G Ninja -DCMAKE_BUILD_TYPE=RelWithDebInfo \
  -DTAIKOR_BUILD_PLUGIN=OFF -DBUILD_TESTING=ON \
  -DTAIKOR_ENABLE_ASAN_UBSAN=ON
cmake --build build-sanitized --parallel
ASAN_OPTIONS=halt_on_error=1:detect_leaks=1 \
UBSAN_OPTIONS=halt_on_error=1:print_stacktrace=1 \
ctest --test-dir build-sanitized --output-on-failure
```

CTest keeps the automatic EngineFuzz corpus bounded: one default run, a
`0x12345678` seed at 200 iterations, and a `0x87654321` seed at 100 iterations.
An optional longer campaign can be run manually when time permits:

```bash
./build-dsp/TaikorEngineFuzzTests --seed 0x12345678 --iterations 10000
```

That long campaign is not an automatic 120-second CTest gate and must not be
reported as a CI result unless it is explicitly run and recorded.

The mutually exclusive optional ThreadSanitizer configuration is
`-DTAIKOR_ENABLE_TSAN=ON`; use it only on a supported Clang/GNU host and do not
combine it with `TAIKOR_ENABLE_ASAN_UBSAN=ON`.

Packaging regression checks:

```bash
python3 -m unittest discover -s Tests -p 'test_packaging.py' -v
```

Platform package commands are documented in the README and the scripts under
`scripts/`. The macOS signing/notarization command must be run with the actual
release identities and keychain profile before treating its output as a
commercial candidate.
