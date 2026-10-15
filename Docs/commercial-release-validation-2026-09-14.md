# Commercial release candidate: validation record

Date: 2026-09-14. Platform tested: local macOS arm64.

Engineering hardening is implemented and the native candidate is rebuilt.
**This is not clearance to sell or publish the instrument.** Asset rights,
JUCE licensing, signing/notarization, and the remaining host/platform checks
in `commercial-release-readiness-2026-09-14.md` are still release gates.
Installed plug-ins were not replaced. Nothing was committed, pushed, merged,
or published during this work.

Nine specialist agents contributed, with at most six running concurrently.
The coordinator implemented the offline ensemble pool and ran integration
validation and performance comparisons.

## Implemented hardening

- MIDI preserves the host's insertion order at equal sample timestamps. UI
  auditions precede host events at sample zero. Negative and out-of-block
  timestamps are discarded rather than moved to an unintended sample.
- Supported MIDI packets require exactly three bytes and seven-bit data.
  Note-offs and velocity-zero note-ons preserve this percussion instrument's
  natural tails. The bounded UI queue drops overflow without waiting; the
  4,096-event callback budget panics and discards the remaining burst.
- UI audition generations prevent old events crossing panic and lifecycle
  boundaries. Invalid channel/block layouts are handled without dereferencing
  missing channels. Every callback exit consumes input MIDI.
- State loading bounds input to 1 MiB, XML nesting to eight levels and element
  count to 512. Unsafe XML/encoding is rejected before parsing. Known scalar
  parameters are rebuilt canonically; malformed values default and finite
  values clamp/snap. Decimal exponent parsing no longer uses the vulnerable
  signed-integer accumulation path. Existing IDs and parameter order remain.
- Sample rates, controllers, parameter conversions, articulation mapping,
  filter input and output boundaries are hardened against invalid/nonfinite
  values. Reverb history cannot retain injected NaN/infinity. Final host
  output has a finite-value guard.
- Measurement scratch-allocation failures return finite unavailable/default
  readouts and remain retryable, instead of terminating through `noexcept`.
- Editor callbacks use lifetime-safe references; missing parameter controls,
  display values, build metadata and prepared-state readouts are guarded.
- CMake includes the new tests, portable thread linkage, optional mutually
  exclusive ASan/UBSan and TSan builds, bounded fuzz cases and CI sanitizer
  coverage. Packaging tests cover version ambiguity and manifest consistency;
  the Linux packaging script remains compatible with macOS Bash 3.2 tests.

## Multicore scope

The pool is **offline-only**. Up to three helpers are created during preparation
when the host reports offline rendering. The caller also renders jobs. Workers
have separate per-player output storage, inherit the caller's floating-point
environment, and finish before the original ascending-member mix runs.
Small segments and insufficient active companions use the serial path.
Worker creation failure also falls back to serial operation.

There is no worker creation in the audio callback. Real-time callbacks never
use this generic pool. A mode change without offline preparation safely stays
serial. Reset, preparation and release retain the normal host requirement that
they do not overlap rendering of the same instance.

Live helper threads require suitable real-time scheduling and host audio
workgroup integration. Offline timings are not a live deadline guarantee:
[Apple Audio Workgroups](https://developer.apple.com/documentation/audiotoolbox/understanding-audio-workgroups).

## Functional and sanitizer results

| Check | Result |
| --- | --- |
| Native build | Standalone, VST3, AU, CLAP, tests and tools built successfully |
| Native regression suite | All 34 tests pass across the full run and targeted corrected-test rerun |
| Packaging | All seven Python regression tests pass |
| ASan/UBSan DSP tests | CPU meter, modal advance, ensemble robustness, offline stress, high-pass and limiter pass |
| ASan/UBSan engine fuzz | Seed `0xC0FFEE`, 200 iterations replayed twice; 400 stream operations, 512 enum probes, deterministic digest and partition checks pass |
| Allocation-failure injection | Three measurement allocation failures, recovery and allocation-free reuse pass |
| Extended ASan/UBSan offline stress | Two independent host threads, eight iterations each; 412,224 stereo sample pairs, zero failures |
| Extended TSan offline stress | Two independent host threads, 12 iterations each; 618,336 stereo sample pairs, zero reported races/failures |
| ASan/UBSan plug-in tests | Processor/editor regression, malformed-state fuzz and MIDI-ordering tests all pass |
| Utility seeded sweeps | 54,272 limiter/high-pass/CPU-meter cases, including invalid-value injections, pass |
| VST3 pluginval | Version 1.0.4, strictness 10, seed 1413564747; 44.1/48/96 kHz and 64/256/1024-sample buffers: SUCCESS |
| VST3 signature integrity | `codesign --verify --deep --strict` succeeds; this does not establish Developer ID trust or notarization |

The first full native run passed 33/34 tests. Two assertions in the newly added
MIDI test incorrectly required different sound when a smoothed pitch/damping
target and a note swapped order without any intervening sample. The test was
corrected, not the production smoothing. It now checks equality for that case,
effects of a 256-sample delay, unchanged preceding audio, and block partition
equivalence. The corrected test passes natively and under ASan/UBSan. Tests
for note-captured CC16/17/18 retain their intentional order-difference checks.

The separate full-plugin sanitizer build required a retry after overlapping
build jobs interrupted object generation. The successful retry and final test
results are retained separately; the failed build log is not hidden.

Local ASan runs used `detect_leaks=0`. These are bounded sanitizer/fuzz results,
not proof of crash freedom, complete leak coverage, or whole-plug-in race
freedom. TSan coverage here is the independent-instance/offline-worker path.
The optional 10,000-iteration campaign is not claimed as completed.

## Sound preservation

The saved pre-hardening paired-renderer executables are in
`build-release-sep14/reference`. The full DSP and plug-in-path matrix compared
14 DSP and ten plug-in/effect workloads at 48 kHz/256 and 96 kHz/64, with two
alternating baseline/candidate rounds. **All 96 audio captures are bit-exact.**
The independent offline benchmark and concurrency tests also compare serial
and parallel output bit-for-bit.

No mode-count reduction, float-state conversion, lower sample rate, shortened
tail, replacement IR, or reduced-quality sound mode was introduced. No new
human listening study was performed; equality establishes identical audio for
the compared corpus, not an exhaustive statement about every possible input.

## Offline performance

Three alternating rounds after warmup; 0.4 seconds of four-drum audio per run.
Numbers below are medians. CPU time is process-wide, including helpers.

| Players | Rate/block | Serial wall ms | Parallel wall ms | Less bounce time | Serial CPU ms | Parallel CPU ms |
| --- | --- | ---: | ---: | ---: | ---: | ---: |
| 4 | 48 kHz / 256 | 481.7204 | 257.7760 | 46.5% | 475.7030 | 492.5600 |
| 4 | 96 kHz / 64 | 923.3042 | 518.7230 | 43.8% | 914.2580 | 957.2450 |
| 8 | 48 kHz / 256 | 1003.5994 | 405.6994 | 59.6% | 968.1370 | 994.7820 |
| 8 | 96 kHz / 64 | 2375.1698 | 840.8312 | 64.6% | 1973.2640 | 1962.0360 |

This is about 2.5-2.8 times the measured eight-player bounce throughput, not a
60-65% reduction in total CPU work or heat. Scheduling affects wall time; the
CSV also records p95/worst block times. Those outliers reinforce why this pool
is not enabled for real-time use.

## Real-time-path CPU impact

The paired modal option remains enabled in the tested native build. Relative
to the saved pre-hardening paired build, results are mixed and generally small.
CPU below means elapsed rendering time divided by the audio deadline, not
whole-machine utilization. Values over 100% exceed that deadline.

| Workload | Rate/block | Before CPU | After CPU | Relative change |
| --- | --- | ---: | ---: | --- |
| Eight players, four drums | 48 kHz / 256 | 252.4763% | 248.3225% | 1.6% less |
| Eight players, four drums | 96 kHz / 64 | 485.0208% | 485.2262% | Essentially unchanged |
| Dense plug-in path | 48 kHz / 256 | 31.0802% | 31.6558% | 1.9% more |
| Automated plug-in path | 48 kHz / 256 | 33.9443% | 34.1822% | 0.7% more |
| Dense plug-in path | 96 kHz / 64 | 59.1615% | 59.2334% | 0.1% more |

Idle safety-check overhead is small in absolute terms: the DSP idle case at
48 kHz rises from 0.0273% to 0.0328% of the audio deadline. Do not advertise a
general live CPU reduction from this hardening pass. The heavy eight-player,
four-drum workload still needs substantial further work for live operation.

## Remaining release gates and exact artifact identity

- The AU built, but direct-path pluginval scanning found no registered AU.
  No successful `auval` result is claimed. The installed AU/VST was not replaced
  merely to make the validator find this candidate.
- Windows/Linux builds and representative DAW installation/session tests were
  not run locally. Their CI and host matrix remain required.
- The published Voxengo terms distinguish commercial use from redistribution
  and prohibit distribution profit. The repository also flags missing separate
  EchoThief redistribution permission. Preserve written clearance or resolve
  the assets before a paid release; do not silently replace sound assets.
  [Voxengo terms](https://www.voxengo.com/impulses/).
- JUCE licensing, Developer ID signing/notarization and final installer checks
  remain unconfirmed. No release was published.

Validated native VST3 bundle:
`build-cpu-sep12/native/Taikor_artefacts/Release/VST3/Taikor.vst3`.
Bundle version and short version: `1.0.0`.
SHA-256 of `Contents/MacOS/Taikor`:
`d683cdb99f18fcb6200f900a1d4627aec2fcc44454a1c5e281f84298f8816def`.

## Evidence files

All paths below are relative to the repository root.

- `build-release-sep14/native-tests.log`
- `build-release-sep14/native-rerun.log`
- `build-release-sep14/packaging-tests-final.log`
- `build-release-sep14/asan-tests.log`
- `build-release-sep14/asan-engine-fuzz.log`
- `build-release-sep14/asan-offline-stress.log`
- `build-release-sep14/tsan-offline-stress.log`
- `build-release-sep14/asan-plugin-tests-final.log`
- `build-release-sep14/pluginval-vst3-final.log`
- `build-release-sep14/pluginval-au.log`
- `build-release-sep14/offline-benchmark.csv`
- `build-release-sep14/audio-cpu-matrix/summary.json`
- `build-release-sep14/audio-cpu-matrix.log`
