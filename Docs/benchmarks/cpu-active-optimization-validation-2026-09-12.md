# Active-playback CPU optimization validation

Date: 2026-09-12. Host: Apple M1 Max, macOS, native arm64 Release build.

These are incremental gains against the already optimized implementation from
`cpu-optimization-validation-2026-09-12.md`, not against the original engine.
The new reference executables and source snapshots are preserved in
`build-cpu-active/reference/`.

## Changes

1. Cache the strain-contributing mode indices instead of scanning the complete
   bank every sample. Eligibility is exactly the former membrane/entry predicate.
2. Reuse that sparse traversal for the control-rate batter and rear strain
   calculations, retaining their original double-precision accumulation.
3. Cache the nonzero shell-boundary ports for both velocity gathering and impulse
   application, avoiding two complete bank scans per sample.

The fixed-capacity index arrays preserve ascending bank order and stop at the
active-mode boundary. They are rebuilt after construction, lifetime sorting and
boundary reconfiguration; silence clears their counts. No allocation occurs in
these render traversals. Every participating mode and every arithmetic operation
within its contribution is retained. No sample-rate, precision, resonator-count,
tail, reverb, or quality-setting reduction was introduced.

## CPU results

The main matrix used all 16 DSP/plugin/reverb scenarios at 48 and 96 kHz, with
64- and 256-sample blocks. Each render lasted 0.6 seconds including 0.2 seconds
of tail, with 0.1 seconds of warmup, three repetitions, and two alternating
reference/candidate rounds. Builds and tests did not overlap timing runs.

Across all four configurations, active engine reductions were:

| Engine workload | Relative CPU reduction |
| --- | ---: |
| Solo | 15.3-25.1% |
| Dense | 14.7-24.7% |
| Four-player ensemble | 12.1-21.9% |
| Eight-player ensemble | 15.5-18.9% |
| Control automation | 12.7-19.1% |

Representative 48 kHz / 256-sample results:

| Workload | Reference budget used | Candidate budget used | Reduction |
| --- | ---: | ---: | ---: |
| Solo engine | 11.592% | 9.521% | 17.9% |
| Dense engine | 67.959% | 56.972% | 16.2% |
| Four-player ensemble | 280.325% | 218.802% | 21.9% |
| Eight-player ensemble | 532.055% | 439.714% | 17.4% |
| Dense full plugin | 72.573% | 64.295% | 11.4% |

Budget used means measured render time divided by audio duration, not whole-machine
CPU utilization. Values above 100% still exceed the real-time budget; these changes
do not make every ensemble configuration real-time safe. The engine's measured
p95 block times also improved in all active matrix configurations.

### Longer plugin confirmation

The short matrix contained a dense-plugin timing regression at 96 kHz / 64 samples
(-8.5% reduction), and only a 5.9% gain for plugin automation at 48 kHz / 256 samples.
Both cases were repeated after the regression suite, with 1.2-second renders,
0.4-second tails, 0.25-second warmup, seven repetitions and three alternating
rounds. These are separate fixtures/durations, not replacements for the raw matrix:

| Workload | Reference budget used | Candidate budget used | Reduction |
| --- | ---: | ---: | ---: |
| Dense plugin, 96 kHz / 64 | 103.872% | 89.610% | 13.7% |
| Plugin automation, 48 kHz / 256 | 45.429% | 40.375% | 11.1% |

The longer measurements did not reproduce the short-run regression. This desktop
was not an isolated benchmark host; timing variability remains visible in the raw
results. No new idle or reverb-only improvement is claimed.

## Audio and regression checks

- All 136 reference/candidate float32 stereo capture comparisons were bit-identical:
  128 in the matrix, two in the initial dense check, and six in the longer checks.
- Every captured sample was finite, with zero sample error against the reference.
  The tested output is identical, rather than merely below an audibility threshold.
- The extended strain test passed 576 bank/state combinations, including shared
  and legacy cavity paths, all four drum families, 8/48/192 kHz, reversed bank
  ordering and partial/empty active banks. It also compares sparse shell updates
  against the original full-bank calculation, including untouched inactive modes.
- The complete CTest suite passed: 23/23 tests, including physical-model,
  shell-boundary passivity/reciprocity, automation, reverb and plugin regressions.
- The complete native build passed, including AU, VST3, CLAP and Standalone.

Bit identity is established for the tested renders, not an exhaustive proof for
every possible host/input sequence. No separate subjective listening claim is made.

## Reproduction and artifacts

Preserve the reference before rebuilding; do not replace it with candidate binaries.
Each comparator output directory must be new.

```sh
cmake --build build-cpu-sep12/native --parallel 6
python3 Tools/BenchmarkCompare.py \
  --baseline build-cpu-active/reference --candidate build-cpu-sep12/native \
  --output build-cpu-active/matrix-repeat --program all \
  --seconds .6 --tail .2 --warmup .1 --repeats 3 --rounds 2
ctest --test-dir build-cpu-sep12/native --output-on-failure --parallel 4
python3 Tools/BenchmarkCompare.py \
  --baseline build-cpu-active/reference --candidate build-cpu-sep12/native \
  --output build-cpu-active/plugin-dense-repeat --program path \
  --case plugin-dense --config 96000:64 \
  --seconds 1.2 --tail .4 --warmup .25 --repeats 7 --rounds 3
python3 Tools/BenchmarkCompare.py \
  --baseline build-cpu-active/reference --candidate build-cpu-sep12/native \
  --output build-cpu-active/plugin-controls-repeat --program path \
  --case plugin-controls --config 48000:256 \
  --seconds 1.2 --tail .4 --warmup .25 --repeats 7 --rounds 3
```

Machine-readable results and captures are in `build-cpu-active/matrix/`,
`build-cpu-active/sparse-dense/`, `build-cpu-active/plugin-dense-confirm/` and
`build-cpu-active/plugin-controls-confirm/`. Each contains `summary.json`, CSV
timings, stderr logs and the uncompressed reference/candidate audio captures.
The build and full regression logs are `build-cpu-active/full-build.log` and
`build-cpu-active/tests.log`.
