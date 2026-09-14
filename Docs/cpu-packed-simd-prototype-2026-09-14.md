# Packed SIMD Prototype Validation

## Decision

**Not accepted for installation.** The prototype preserves the tested audio, but the broader paired benchmarks show roughly 10-12% more CPU in the demanding DSP fixtures. An initial isolated comparison suggested a 6.8% improvement; the larger matrix and long-tail comparison did not reproduce that benefit.

The existing installed VST was not replaced. The experimental implementation remains in the local source and build. `TAIKOR_ENABLE_PACKED_MODAL_BANK` currently defaults to `ON`; disabling that option and rebuilding selects the established scalar storage/rendering path. Disabling the experiment and continuing with a lower-overhead layout has been proposed, but not performed pending the follow-up decision.

## Implementation

The prototype adds a single-source-of-state, structure-of-arrays modal bank for Apple silicon. Five contiguous double-precision arrays hold the recurrence coefficients and current/previous states. Existing contact, damping, boundary, and observation code accesses those same arrays through bound resonator views rather than a periodically synchronized copy.

The recurrence processes two independent modes per Neon vector. Its fused arithmetic follows the original Apple Clang release assembly. The nonlinear correction, rim feedback reduction order, and stereo summation retain their previous logical ordering. Other platforms retain the scalar engine path, and the CMake option provides an explicit scalar fallback on Apple silicon.

Copy construction creates independent state; assignment copies values without replacing the destination slot's bindings. This is necessary for the existing insertion sort and geometry rebuilds. The implementation uses no per-render allocation or reduced mode count, lower precision, tail truncation, or lower simulation rate.

The adapter introduces reference indirection and changes metadata layout. These are plausible costs that can outweigh the recurrence improvement, but their individual contributions have not been isolated by a new hardware-counter profile. The results do not establish that SIMD itself is unsuitable, only that this particular integration is not a useful optimization.

## Correctness and sound

- All 26 CTest tests passed; total suite time was 158.01 seconds.
- The expanded modal-advance test passed 960 bit-exact cases across all four drum families and 8, 48, 96, 192, and 384 kHz, including odd/partial/empty active prefixes and rim/nonlinear combinations.
- The dedicated packed-bank test passed all 18 active-prefix lengths, long recurrence runs, sorting, copying, move construction, and direct coefficient/state changes.
- AddressSanitizer and UndefinedBehaviorSanitizer passed the packed-bank ownership/recurrence test.
- The initial ensemble comparison contributed 2 bit-identical audio comparisons.
- The full DSP/plugin benchmark matrix contributed 96 bit-identical audio comparisons.
- A control-automation capture with 14 seconds of stimulus and a 12-second tail contributed 1 bit-identical comparison.

**Total: 99 bit-identical before/after audio comparisons.** This establishes equality for the tested corpus, not a universal proof over every possible input. No subjective listening claim is needed to explain equality of those captures, and no new human listening study was performed.

## Performance

The baseline executables were rebuilt from the current source before the experiment and preserved separately. Reference and candidate runs were alternated in A/B and B/A order. The main matrix used three repetitions per run, two rounds, 0.6 seconds of stimulus, 0.2 seconds of tail, and 0.1 seconds of warmup. Builds and CTest completed before the matrix began.

Percentages below mean processing time divided by audio duration, not total-machine CPU percentage. Negative savings indicate a regression.

| DSP fixture | Configuration | Baseline | Prototype | CPU reduction |
|---|---|---:|---:|---:|
| Solo | 48 kHz / 256 | 5.7354% | 5.9441% | -3.6% |
| Dense | 48 kHz / 256 | 32.5265% | 35.9846% | -10.6% |
| Four-player ensemble | 48 kHz / 256 | 131.2211% | 142.7103% | -8.8% |
| Eight-player ensemble | 48 kHz / 256 | 258.5696% | 284.9580% | -10.2% |
| Controls | 48 kHz / 256 | 82.9881% | 92.5832% | -11.6% |
| Dense | 96 kHz / 64 | 63.5143% | 70.0018% | -10.2% |
| Four-player ensemble | 96 kHz / 64 | 253.6992% | 278.2028% | -9.7% |
| Eight-player ensemble | 96 kHz / 64 | 504.8135% | 557.7849% | -10.5% |
| Controls | 96 kHz / 64 | 165.4483% | 181.7200% | -9.8% |

Eight players on a single drum family also regressed: approximately 0.5-2.2% at 48 kHz / 256, and 2.5-6.2% at 96 kHz / 64 in this matrix. The long control/tail comparison measured 97.9130% to 109.2504%, an 11.6% regression; it used one repetition and is primarily a correctness/stability check.

The full-plugin and unchanged reverb-only fixtures showed substantial variability, including large differences in code not modified by the prototype. Those results should not be attributed to the modal change, and they reinforce the need to avoid selecting a candidate from one favorable measurement. No general speedup claim is supported.

The first isolated eight-player comparison measured 371.8807% to 346.7257%, or 6.8% lower CPU. Its much higher absolute baseline than the subsequent matrix and the contradictory broader results make it unsuitable as the acceptance result.

## Next implementation direction

Retain the established renderer as the default before proceeding. A subsequent prototype should avoid paying a bound-reference adapter cost throughout every hot physical-state reader. Candidates include a more direct compact state-access design or a smaller kernel-only layout change whose full-engine benefit can be isolated before broad integration.

Keep precision and physics unchanged for that comparison. Reuse the new ownership/odd-prefix coverage, compare full-engine workloads, and do not combine layout changes with mode retirement or lower-rate tails in the same experiment. This result supplies no progress toward the 90% reduction target.

## Local evidence

- [Main benchmark summary](/Users/vojta/Dev/virtual-instrument-taikor/build-cpu-sep14-simd/matrix/summary.json)
- [Main benchmark log](/Users/vojta/Dev/virtual-instrument-taikor/build-cpu-sep14-simd/matrix.log)
- [Initial ensemble summary](/Users/vojta/Dev/virtual-instrument-taikor/build-cpu-sep14-simd/ensemble-first/summary.json)
- [Long-tail summary](/Users/vojta/Dev/virtual-instrument-taikor/build-cpu-sep14-simd/long-tail/summary.json)
- [CTest log](/Users/vojta/Dev/virtual-instrument-taikor/build-cpu-sep14-simd/tests.log)
- [Build log](/Users/vojta/Dev/virtual-instrument-taikor/build-cpu-sep14-simd/build-all.log)
- [Packed bank implementation](/Users/vojta/Dev/virtual-instrument-taikor/Source/DSP/PackedModalBank.h)
- [Packed bank tests](/Users/vojta/Dev/virtual-instrument-taikor/Tests/PackedModalBankTests.cpp)

The pre-experiment DSP sources and benchmark executables are preserved in `build-cpu-sep14-simd/reference/`. No commit, push, or merge was performed.
