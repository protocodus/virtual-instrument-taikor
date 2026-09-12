# CPU optimization and audio validation, 2026-09-12

The largest measured savings are in idle processing. Solo engine processing
improves by 2.6-7.6%; dense ensembles remain approximately unchanged. The final
implementation produces bit-identical audio in all 160 comparisons against the
original reference, covering 80 scenario/configuration combinations. All 23
CTest tests pass.

## Measurement conditions

- Apple M1 Max, native arm64, Release builds with matching compiler options.
- CPU denotes processing wall time divided by rendered audio duration. It is
  a percentage of one real-time audio thread's budget, not total application
  CPU including the editor. Values above 100% exceed that budget in these
  offline stress cases, in both reference and candidate.
- Main matrix: 48/96 kHz, 64/256-sample blocks, one-second scores including
  0.3 seconds of tail, 0.1 seconds of warmup, two repetitions per invocation.
- Two rounds reverse executable order: reference/candidate, then
  candidate/reference. Builds and tests do not run concurrently with benchmarks.
- The reference is the working tree at the start of this task, including its
  existing uncommitted changes. It is not the repository's HEAD revision.
- The starting processor contained a meter shortcut calling an unavailable
  ensemble method. The reference omits that shortcut so it can compile; its
  audio sample generation is unchanged. The final shortcut checks silence
  before rendering and preserves the original meter recurrence.

## Representative results

CPU columns are percentages of the real-time thread budget. Negative reduction
means the candidate measured slower. Changes around 1% and isolated outliers
are not evidence of a reliable speedup or slowdown.

| Workload | Rate / buffer | Reference CPU | Final CPU | Reduction |
| --- | --- | ---: | ---: | ---: |
| Idle plugin callback | 48 kHz / 64 | 0.1550% | 0.0510% | 67.1% |
| Idle plugin callback | 48 kHz / 256 | 0.1510% | 0.0353% | 76.6% |
| Idle plugin callback | 96 kHz / 64 | 0.3121% | 0.1053% | 66.3% |
| Idle plugin callback | 96 kHz / 256 | 0.3043% | 0.0706% | 76.8% |
| Solo engine | 48 kHz / 64 | 10.2523% | 9.4697% | 7.6% |
| Solo engine | 48 kHz / 256 | 9.5163% | 9.2520% | 2.8% |
| Solo engine | 96 kHz / 64 | 19.0842% | 18.2786% | 4.2% |
| Solo engine | 96 kHz / 256 | 18.7509% | 18.2557% | 2.6% |
| Dense plugin callback | 48 kHz / 64 | 53.5788% | 52.8655% | 1.3% |
| Dense plugin callback | 48 kHz / 256 | 53.6348% | 54.0560% | -0.8% |
| Dense plugin callback | 96 kHz / 64 | 106.7481% | 103.8432% | 2.7% |
| Dense plugin callback | 96 kHz / 256 | 104.9393% | 105.7329% | -0.8% |
| Eight-player engine | 48 kHz / 256 | 447.5191% | 447.7338% | -0.05% |
| Eight-player engine | 96 kHz / 256 | 868.3859% | 875.1851% | -0.8% |

Active reverb has no substantial, consistent improvement. Short measurements
contained large outliers, so Opera and room switching were repeated with
two-second scores, seven repetitions and three alternating rounds. Opera's
reduction ranged from -1.0% to +1.0%; room switching ranged from -0.3% to +1.5%.
The larger short-run excursions were not reproduced. All original measurements
are retained rather than replacing the unfavorable observations.

A second implementation pass removed unnecessary activity refreshes in settled
reverb. Compared directly with the first pass, silent-reverb CPU improved
2.3-4.6%, with identical output. Its absolute CPU cost was already below 0.01%.

## Disposition of the ten ideas

1. **Active voice lists:** `TaikoEngine` keeps compact, slot-ordered contact and
   physical-bank views. Rendering and the contact solver iterate these views;
   block entry and retirement refresh them without allocation or reordered sums.
2. **Settled controls:** test whether the next floating-point smoothing step
   actually changes the state. Moving controls retain their original sample
   recurrence and pitch-cache invalidation timing; no epsilon snapping is added.
3. **Frozen output:** settled silent engines fill zeros without a per-sample
   housekeeping loop. Moving controls and nonzero display envelopes still advance.
4. **Reverb activity cache:** retain the active-room indices and smoothing state;
   refresh them when parameters or fades change, rather than every settled chunk.
5. **Unchanged reverb parameters:** the processor skips repeated room/mix updates,
   and the reverb setter also exits for unchanged sanitized values. Prepare and
   state restoration retain their parameter initialization behavior.
6. **Meter optimization:** preserve the exact sample peak/release recurrence.
   A plain block maximum is not equivalent to that recurrence. Fully silent
   segments skip sample reads, and zero meters need no decay loop. Silence is
   established before rendering, including pending ensemble hits and wet state.
7. **Settled dry reverb:** a cached bypass returns before chunk processing;
   history is cleared on the transition instead of repeatedly resetting it.
8. **Wet copies:** convolution uses `ProcessContextNonReplacing`, reading the
   original dry block and writing each wet buffer without the external copy pair.
9. **Block invariants:** smoothing/activity decisions and the convolution input
   block are established outside sample loops. Public pointer/preparation guards
   were already at block entry and remain there.
10. **Benchmarking:** `TaikorBenchmarkAudioPath` supplements `TaikorBenchmarkCPU`
    with full-callback and isolated-room cases, process/event block counters,
    timing distributions and deterministic captures. `Tools/BenchmarkCompare.py`
    alternates builds, records CPU metrics, checks every sample for finiteness,
    compares capture bytes, and exits unsuccessfully if audio differs.

The physical model, full stereo IRs, sample rate, latency and arithmetic
precision are preserved.

## Audio and regression coverage

- 128 comparisons cover the final DSP/callback/reverb matrix: idle, solo, dense,
  four/eight performers, live controls, all rooms, dry bypass, fades and silence.
- Two comparisons render a complete 14-second solo take at 48 and 96 kHz,
  including physical retirement and the transition to frozen output.
- Four comparisons render three-second room impulse tails across both rates
  and buffer sizes, extending past the declared IR tail bound.
- Two comparisons exercise plugin controller/room/panic automation with seed 42.
- 24 comparisons repeat Opera and room switching with longer scores and more
  repetitions. These bring the final-versus-original total to 160.
- Every capture has zero sample error, including signed-zero bits. Each timed
  invocation also verifies deterministic audio between its own repetitions.
- Eight additional comparisons validate the reverb refinement against pass one.
- `Taikor.CPUFastPath` tests exact meter behavior through MIDI slicing, panic,
  reverb and a 13-second block that ends frozen. It also tests repeated reverb
  setters, reprepare at different rates, reset and cached bypass state.
- All 23 CTest tests pass, including engine physics, ensemble, IR reverb,
  processor and the two benchmark smoke tests. AU, VST3, CLAP and Standalone
  targets build successfully for native arm64.

The existing processor round-trip test expected nondefault reverb values
without setting them before saving. Its setup now sets those values; its
restoration assertions remain intact.

## Required procedure for future CPU changes

Preserve a same-options Release reference before editing audio code. Benchmark
the change, compare deterministic full-level audio captures, and run the
relevant tests. Report absolute CPU cost as well as relative reduction. Use
longer repeated runs to investigate apparent regressions; do not accept timing
outliers or lower-quality rendering as an optimization.

Example using fresh native build/reference/output directories:

```sh
cmake -S . -B build-cpu -DCMAKE_BUILD_TYPE=Release -DTAIKOR_BUILD_UNIVERSAL=OFF
cmake --build build-cpu --target TaikorBenchmarkCPU TaikorBenchmarkAudioPath --parallel 6
mkdir build-cpu-reference
cp build-cpu/TaikorBenchmarkCPU build-cpu/TaikorBenchmarkAudioPath build-cpu-reference/

# After the implementation change:
cmake --build build-cpu --parallel 6
python3 Tools/BenchmarkCompare.py \
  --baseline build-cpu-reference --candidate build-cpu \
  --output build-cpu-comparison --repeats 3 --rounds 2
ctest --test-dir build-cpu --output-on-failure --parallel 4
```

Use `--program dsp` for a plugin-free build. Use a new output directory per run;
the comparator refuses to mix new measurements with existing captures.

Local evidence for this run is under `build-cpu-sep12/`: `baseline/` contains
the reference executables and `reference/` contains the initial source snapshot.
`matrix/summary.json` supplies the final engine results; `path-final/summary.json`
supplies the final processor/reverb results. Other reports are in `drum-tails/`,
`room-tails/`, `seed42/`, `refinement/`, `opera-confirm/` and `switch-confirm/`.
Each report retains raw timing CSVs and interleaved stereo float32 captures.
`final-tests.log` records the passing suite.
