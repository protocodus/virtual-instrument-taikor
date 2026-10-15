# Taikor CPU Optimization Research

## Executive recommendation

The strongest next investment is a compact, vectorized implementation of the existing physical model, followed by an adaptive tail renderer only where its assumptions can be demonstrated. This sequence first removes implementation overhead without removing acoustic detail, then investigates reducing the amount of simulation that is genuinely necessary. Multicore rendering should be treated as a separate playback-reliability feature, not as evidence of lower total CPU consumption.

Three priorities follow:

1. **Retain every mode and improve the execution layout.** Replace the large, strided mode representation in the recurrence hot path with compact structure-of-arrays or small array-of-structure-of-arrays tiles. Prototype SIMD across modes and, separately, across players. Preserve precision, event timing, feedback ordering, and ordered output accumulation initially.
2. **Establish which state is actually eligible for cheaper evolution.** Distinguish contact, nonlinear ringing, coupled ringing, genuinely free decay, and silence. A ringing drum is not necessarily a linear oscillator bank merely because its stick has detached.
3. **Investigate perceptually transparent tail simplification.** Candidates include selective multirate evolution, bounded mode retirement, and reduced coupling models. These require state-aware transitions, restriction of approximation error, and listening tests, not just a quieter difference signal.

These are recommendations, not measured new speedups. The last available Taikor measurements already establish a 12-24% incremental DSP improvement with 99 bit-identical audio comparisons. They do not establish the performance or sound quality of any proposal in this report.[^1]

**A further 90% reduction in total CPU is an unproven research target.** No source reviewed demonstrates a drop-in, sound-equivalent 10x acceleration of Taikor's particular nonlinear, coupled model on an M1 Max. Large accelerations reported for other modal synthesizers are useful architectural evidence, but their workloads, latency budgets, reference implementations, and approximations differ.

The proposed delivery is a sequence of independently benchmarked changes, not a single rewrite that changes layout, precision, physics, and scheduling simultaneously. The scalar reference must remain available throughout development.

## Existing performance and the actual target

### Last measured reference

The following values are from the September 13 local validation on an Apple M1 Max, after the preceding optimizations. They are historical benchmark measurements, not new measurements taken for this report. The eight-player, all-drum fixture is a deliberate stress case; it should not be described as the cost of every eight-player musical passage.[^1]

| Workload | 48 kHz / 256 samples | 96 kHz / 64 samples |
|---|---:|---:|
| Eight players, O-daiko | 60.7% | 110.3% |
| Eight players, Nagado-daiko | 62.6% | 116.1% |
| Eight players, Okedo-daiko | 53.1% | 102.6% |
| Eight players, Shime-daiko | 51.5% | 99.4% |
| Eight players, all four drum families | 258.8% | 482.5% |
| Full-plugin dense fixture, a different workload | 31.7% | 60.5% |

The percentage is processing time divided by the corresponding audio duration. It is not a percentage of all CPU cores. The DSP stress fixture and the full-plugin dense fixture are not interchangeable measurements.

At 48 kHz with a 256-sample buffer, the nominal deadline is 5.333 ms. At 96 kHz with a 64-sample buffer, it is only 0.667 ms. The previous all-drum stress measurements also had p95 processing times of approximately 16.918 ms and 3.861 ms respectively, so reducing the average alone will not resolve playback reliability.[^1]

### Throughput needed for useful headroom

An 80% deadline-use target is an engineering planning margin, not a guarantee against dropouts. Host overhead and other instruments consume time too. Applying that target to the prior average DSP measurements gives:

| All-drum eight-player configuration | Required overall speedup to reach 80% | Required reduction in elapsed rendering time |
|---|---:|---:|
| 48 kHz / 256 | 3.23x | 69.1% |
| 96 kHz / 64 | 6.03x | 83.4% |

These are calculations from the benchmark, not forecasts. A 2x total speedup would still leave the 48 kHz all-drum stress fixture at about 129% of its deadline. A 90% reduction relative to the September 13 implementation means a 10x speedup, reducing those averages to approximately 25.9% and 48.2%.

There must be two independently reported success measures:

- **Deadline performance:** callback elapsed time, high percentiles, observed maximum, and missed deadlines.
- **Total work:** CPU time summed across audio and worker threads per rendered audio second; energy measurements where practical.

Multithreading can improve the first while leaving the second unchanged or making it worse. GPU offload can similarly lower CPU accounting while moving the work elsewhere.

### Why one fast kernel cannot guarantee 90%

For a fraction `p` of execution accelerated by a factor `s`, the ideal overall speedup is `1 / ((1 - p) + p / s)`. For example, making a stage that occupies 55% of runtime four times faster reduces overall time by only 41.25%. Eliminating that stage entirely would still save no more than 55%.

The earlier sampling profile suggested substantial modal-render work, with additional strain, boundary, and control-update costs. Those sampled shares came from the older implementation under profiling and must not be used as precise shares of the current build. The calculation explains the architectural constraint; it is not a prediction based on a freshly measured hotspot breakdown.

## Why drums have different costs

The prior default-parameter diagnostics found the same counts for all four families: 296 modes, 74 strain modes, and 54 active shell-boundary ports. At those settings, the difference is therefore not simply that one drum has a larger allocated mode bank.[^1]

| Drum | Model diameter | Nominal Don contact at velocity 0.85 | Model tail diagnostic |
|---|---:|---:|---:|
| O-daiko | 150 cm | 0.451 ms | 5.856 s |
| Nagado-daiko | 78 cm | 0.264 ms | 2.877 s |
| Okedo-daiko | 40 cm | 0.140 ms | 1.730 s |
| Shime-daiko | 30 cm | 0.145 ms | 2.280 s |

These are internal model diagnostics, not acoustic measurements of real instruments. Contact duration, damping, and state-dependent processing provide plausible explanations for differing costs, but the earlier work did not isolate their individual causal contributions. The slightly higher Nagado result in the eight-player matrix also shows why diameter alone is not a reliable CPU predictor.

The scaling problem is more direct. Each player owns independent physical drum state, and the current rendering is serial. Eight players with four active families imply `8 * 4 * 296 = 9,472` resonators. Updating them once per output sample produces roughly 455 million modal recurrence updates per second at 48 kHz and 909 million at 96 kHz, before the other processing stages.

A useful future diagnostic should report cost by both drum and state: contact-active samples, nonlinear-active samples, coupled-tail samples, free-tail samples, and mode-samples processed. That can distinguish a more expensive update from a drum simply remaining active longer. Do not infer that a 0.451 ms contact explains a multi-second rendering burden by itself.

## Optimization opportunities

The priorities below combine expected applicability, potential scope, and sound risk. No row has a validated Taikor speedup yet. Opportunities overlap and their hypothetical savings must not be added together.

| Priority | Candidate | Primary benefit | Sound classification | Recommendation |
|---|---|---|---|---|
| 1 | Compact recurrence storage and SIMD across modes | Less data movement and more work per instruction | Exact-result candidate | First prototype |
| 2 | SIMD across players | Vectorized processing without reordering each player's mode sum | Exact-result candidate | Compare against priority 1 |
| 3 | Dependency-aware loop specialization and fusion | Fewer branches and repeated loads | Exact-result candidate | Follow measured bottlenecks |
| 4 | Immutable coefficient sharing and dependency-based updates | Less repeated setup and coefficient work | Exact-result candidate | Audit incremental opportunities |
| 5 | Profile-guided optimization and link-time optimization | Better generated code | Exact-result candidate | Low-intrusion experiment |
| 6 | Exact zero/disconnected-state skipping | Avoid mathematically unnecessary work | Exact-result candidate | Narrow eligibility only |
| 7 | Selective precision and arithmetic changes | Denser vectors and cheaper arithmetic | Numerically changed | After exact SIMD |
| 8 | State-aware adaptive multirate tails | Fewer mode updates | Perceptual approximation | Highest-priority architecture experiment |
| 9 | Error-budgeted mode retirement | Fewer active modes | Perceptual approximation | Start with late, uncoupled tails |
| 10 | Reduced or factored coupling models | Less boundary/strain work | Exact algebra or approximation, depending on method | Profile and derive first |
| 11 | Shared linear tail banks | Fewer duplicated bank evolutions | Conditional mathematical equivalence | Restricted research case |
| 12 | Frequency-domain late-tail synthesis | Amortized synthesis of many modes | Perceptual approximation | Secondary tail experiment |
| 13 | Real-time multicore player rendering | Shorter critical path | Sound-preserving candidate | Separate reliability track |
| 14 | GPU or learned surrogate | Different compute architecture/model | High integration or sound risk | Defer |

### 1. Compact recurrence storage and SIMD across modes

Arm documents contiguous structure-of-arrays layouts as favorable for SIMD. Neon supports four 32-bit or two 64-bit floating-point elements per 128-bit vector; lane width is not an application speedup guarantee.[^2][^3] LLVM provides vectorization diagnostics that explain why particular loops do or do not vectorize.[^4]

**Taikor application.** The September 13 implementation places frequently used fields in a 128-byte prefix of a 384-byte aligned `Mode`. That was an improvement over the previous arrangement, but it still leaves the recurrence fields separated by a large stride. For the 9,472-mode stress case, complete mode records occupy about 3.47 MiB, while the hot prefixes alone account for about 1.16 MiB. These are storage calculations, not measurements of bytes fetched on every sample.

Use separate contiguous arrays for recurrence state and coefficients, with other arrays for observation, coupling, and infrequently changed metadata. As an illustration, four double-precision scalars per mode occupy about 0.29 MiB for the same population; this is not a claim that Taikor's entire hot model can fit in four scalars. A small tiled layout may offer better locality across consecutive stages than either the current large records or one enormous global collection of arrays.

Prototype on one drum first. Keep its scalar implementation, preserve its existing arithmetic precision, and compare recurrence state as well as audio. Vectorize independent modes at the same sample instant. Do not accidentally advance each mode through a whole block while other modes remain behind: per-sample coupling can make that transformation incorrect.

Initially preserve the original ordered stereo sum and the order of boundary reductions. A changed horizontal reduction can alter rounding even when every individual oscillator update is identical. If ordered mixing becomes the bottleneck, evaluate it as a separate change rather than weakening the entire build's numerical rules.

**Stop condition:** the compact/vectorized path must improve full-engine benchmarks, not merely an isolated recurrence loop. Conversion, gathering, sparse corrections, and output mixing may absorb the kernel gain.

### 2. SIMD across players

**Taikor-specific proposal.** Assign each SIMD lane to one player's copy of a drum, rather than to a different mode. Iterate through mode indices in their original order. Each lane can then accumulate its own output in the same sequence as the scalar implementation, before the final player outputs are combined in the existing player order.

This potentially exposes parallelism in observation and coupling calculations that would otherwise involve cross-lane reductions. It does not require the players to have identical coefficients, provided each lane loads the correct values. It does require compatible traversal shapes or explicit handling of divergence.

Preserve player-specific dynamics, pan, humanization, and event order. Inactive lanes, different articulations, and different sparse-mode sets can waste vector capacity. Repacking eight players every sample would likely undermine the design, so establish stable groups during setup or at appropriate existing event boundaries.

The useful comparison is mode-lane SIMD versus player-lane SIMD versus a small tiled hybrid, measured across solo, sparse ensembles, and all-eight-player workloads. Do not commit the entire engine to a layout before that comparison.

### 3. Dependency-aware specialization and loop fusion

**Taikor-specific proposal.** Generate a small number of kernels for meaningfully distinct states, such as active contact, nonlinear ringing, and proven free evolution. Keep the sparse nonlinear correction already introduced; do not count its existing benefit again. Separate coefficient and observation data by update frequency so the free recurrence does not repeatedly visit unrelated fields.

Before fusing loops, write down the dependency sequence. In the known implementation, the free modal step, nonlinear correction, boundary update, and observation have meaningful ordering. Fusion is safe only where the same values are consumed at the same logical sample. A smaller number of loops is not necessarily faster if it introduces dependencies, increases register pressure, or prevents vectorization.

Preserve the order of each floating-point reduction. Where two passes genuinely consume the same unchanged state, retaining an intermediate result can be preferable to recalculating it. Where one pass changes state used by another, retaining an old value silently changes the model.

**Stop condition:** reject a specialization whose branch-selection or maintenance overhead outweighs its gain in realistic mixed-state workloads.

### 4. Immutable coefficients and dependency-based updates

**Taikor-specific proposal.** Share immutable data for players with genuinely identical physical parameter sets, without sharing their dynamic resonator state. Organize derived coefficients by dependency: geometry, sample rate, tension, damping, strike position, microphone placement, and state-dependent nonlinear terms must not be treated as one interchangeable cache key.

The previous work already optimized settled smoothers and cached UI measurement results. The remaining question is whether repeated coefficient construction or repeated identical immutable data still contributes significantly. Do not advertise already-cached work as a new optimization.

Where expensive preparation can happen outside the callback, publish completed immutable objects through a bounded ownership protocol. Avoid allocations, object destruction, and reference-count cleanup on the audio thread. Preserve the defined sample timing of audible parameter changes; delayed publication is not automatically sound-equivalent.

Distinguish an exact cache from a lookup approximation. Reusing a coefficient for unchanged inputs is one category. Quantizing a moving control, interpolating a table, or updating it less often is a different category requiring sound tests.

### 5. Compiler-guided improvements

LLVM's loop-vectorization reports identify blocked vectorization, and ThinLTO enables cross-module optimization.[^4][^5] These are tools for improving the implementation, not evidence that another compiler flag will provide a particular percentage gain.

Experiment with optimization remarks, appropriate release settings, link-time optimization, and representative profile-guided optimization. Training workloads should include all drum families, dense ensembles, long tails, and parameter changes. Benchmark non-instrumented release binaries after training, and retain small-workload cases to detect specialization regressions.

Do not enable global `-ffast-math` as the first step. Clang documents that it permits reassociation, reciprocal transformations, and assumptions about exceptional values; contraction into fused operations also affects rounding.[^6] Preserve the reference's contraction behavior for exact-result work. Neither globally enabling FMA nor globally disabling an FMA already used by the reference guarantees bit identity.

Available Apple Clang behavior must be established for the installed toolchain. Current upstream LLVM documentation is not proof that every documented option behaves identically in that toolchain.

### 6. Exact zero and disconnected-state skipping

**Taikor-specific proposal.** A mode can be skipped exactly only when its state and all relevant future inputs satisfy conditions that make the skipped evolution redundant. Zero microphone gain is insufficient: the mode may still influence a stick contact, strain estimate, shell port, rear head, or another observable state.

A safe initial case is a zero-state subsystem with no excitation, no coupling input, and no hidden update obligation until a known event. More general elimination requires proving that the removed state cannot affect any retained state or future output under the supported controls. A current zero crossing is not zero energy.

If eligible modes are sparse, maintain bounded index sets and rebuild them at existing state-change boundaries. The bookkeeping itself needs benchmarking. The already implemented idle fast paths mean this opportunity may be modest; it should not lead the 90% strategy.

### 7. Selective precision and arithmetic changes

**Taikor-specific proposal.** After a same-precision SIMD reference exists, consider float32 for well-conditioned subsets while retaining higher precision for sensitive state and reductions. Do not begin by converting the entire engine to float32. Larger lane counts provide a hardware opportunity, not a quality guarantee.[^3]

The sensitivity concern follows directly from a conventional resonator coefficient `a = 2*r*cos(theta)`: a small coefficient perturbation produces approximately `delta_theta = -delta_a / (2*r*sin(theta))`, with `r` held fixed. Near zero frequency, the denominator is small. Long decays and coupled feedback can also accumulate errors over many updates. The exact consequences depend on Taikor's realization and parameter range.

Test explicit FMA, reciprocal approximations, or a different resonator realization as separate candidates. Preserve original denominators and safety constraints in the initial SIMD path. An alternative can be more numerically accurate yet still differ from the accepted instrument sound.

**Acceptance:** long-duration state stability, pitch and decay checks, dense retriggering, parameter extremes, and the perceptual protocol below. A good average signal-to-error ratio alone is insufficient.

### 8. State-aware adaptive multirate tails

Zambon's thesis presents full-rate contact processing followed by multirate free evolution, including state transformations. It also describes interpolation delay, transition artifacts, and limitations of its single-object MATLAB demonstration. Its larger tabulated speedups assume a particular mode distribution and exclude interpolation cost; they are not end-to-end Taikor measurements.[^7] A separate Microsoft Research implementation also uses lower-rate processing for lower-frequency modal components.[^8]

**Taikor application.** Introduce an eligibility state machine, not a fixed timer after note-on. Contact completion is necessary but may not be sufficient: tension modulation, shell coupling, rear-head interactions, or external control changes can keep a bank unsuitable for independent lower-rate evolution.

Start with a deliberately restricted experimental case: one drum, fixed controls, no active contact, and a subsystem proven independent of full-rate feedback. Retain full-rate evolution for every ineligible component. Preserve pole behavior when changing rate, map state consistently, and account for the interpolation filters' history.

Require an event to promote the relevant state back to full rate before its first affected sample. A new strike on a partly decayed drum is a critical case. Resetting the bank, replaying a sampled tail, or starting the attack from zero is not equivalent to continuing the physical state.

Do not hide conversion delay. Align the attack and tail paths, report any added plugin latency, and reject transitions that smear the transient or comb-filter the overlap. At 96 kHz / 64 samples, even a small additional delay is material. The initial experiment should not quietly turn every 96 kHz project into a 48 kHz internal simulation.

**Stop condition:** if realistic playing leaves little eligible time, or if preserving coupling and conversion costs removes the benefit, abandon this branch before a full-engine integration. Measure eligibility and net cost before building a complicated architecture around an assumed long free tail.

### 9. Error-budgeted mode retirement

Human-subject work by van den Doel and colleagues supports perceptually informed selection of modal components for tested contact sounds. It does not establish a universal number of modes for drums or prove that feedback states may be removed safely.[^9]

**Taikor application.** Begin with absolute, conservative error budgeting in genuinely uncoupled late tails, not cross-instrument masking. Derive a bound on each candidate mode's future output contribution using both state components and both microphone projections. Require that the combined omitted contribution stays below the selected budget; evaluating modes independently can underestimate the sum.

For deterministic output bounds, the triangle inequality gives `abs(error[n]) <= sum(bound_k[n])`. An energy-summing rule can be less conservative, but needs assumptions about correlation. Eight equal, coherent errors can add about 18 dB in amplitude relative to one, and 32 can add about 30 dB. The correct population for a mode-level bound may be much larger than the player count.

Use hysteresis and smooth output transitions where necessary, without confusing output fading with correct internal-state handling. A mode removed from evolution cannot necessarily be restored faithfully after a later control change. Protect components involved in feedback, beating pairs, stereo differentiation, or position-dependent excitation unless their approximation has its own demonstrated bound.

Cross-player masking is a later option because a drum may be soloed, panned, or exposed when another player stops. Do not depend on unrelated DAW tracks masking an error; the plugin does not control that context.

### 10. Coupling structure and reduced models

Ducceschi, Bilbao, and Webb describe nonlinear modal networks with an energy-based formulation and a rank-one update that can be evaluated in linear complexity using the Sherman-Morrison identity. This is evidence that solver structure can matter more than superficial code tuning, not evidence that Taikor currently uses a dense iterative solve.[^10]

**Taikor application.** First identify whether any remaining boundary or strain operation contains redundant projections or reusable low-rank structure. Exact algebraic factorization and a reduced physical model are different proposals. The existing sparse strain and boundary traversals already exploit some structure and must remain the comparison point.

If an approximation is necessary, preserve the effective behavior at the physical input and output ports, not just the average spectrum of one recorded hit. A mode that radiates weakly may still be dynamically important. Fit and evaluate across articulations, velocities, tensions, and control trajectories, keeping some combinations held out.

Preserve passivity or an appropriate energy bound for unforced, fixed-parameter operation. An audible model that behaves well for one second but gains energy during long tails is not an optimization. Deliberate parameter changes can exchange energy, so energy tests must specify their boundary conditions.

### 11. Shared linear tail banks

**Conditional mathematical opportunity.** For identical linear dynamics,

```text
x_i[n+1] = A*x_i[n] + B*u_i[n]
y_i[n]   = C*x_i[n]

X[n]     = sum_i x_i[n]
X[n+1]   = A*X[n] + B*sum_i u_i[n]
sum_i y_i[n] = C*X[n]
```

Thus, if only the summed output is needed, several identical linear systems may be represented by an aggregated state. Fixed output weights can sometimes be handled by separate weighted aggregate states, for example for left and right channels. This derivation is mathematical; changing the order of floating-point operations need not be bit-identical.

**Taikor restriction.** Individual nonlinear feedback, varying player tuning, independent damping, changing panning, humanization, or the need to retrieve one player's physical state can invalidate aggregation. Once individual states are discarded, a future per-player edit cannot generally recover them from the sum. Maintaining full-rate shadow copies would consume the work intended to be saved.

Investigate only a strictly defined linear tail cohort with a correct plan for future controls and retriggers. Eight players do not automatically imply an 8x reducible workload. This is not recommended as the first architecture change.

### 12. Frequency-domain late-tail synthesis

Bonneel and colleagues report 5-8x acceleration of modal summation in their test scenes using sparse frequency-domain synthesis, with quality limitations particularly for fast decays and high frequencies. Their implementation used 1024-point transforms and 512-sample reconstruction steps at 44.1 kHz. Those conditions differ materially from a 64-sample, 96 kHz instrument callback.[^11]

**Taikor application.** Consider this only for late components whose behavior is sufficiently predictable. Keep attacks and interactive nonlinear state in the time domain. A shared stereo spectral tail bus might amortize transforms, but changes to individual players, their spatial weights, or coupling can undermine that advantage.

Do not adopt audiovisual impact-delay tricks from game rendering for live musical notes. Preserve MIDI timing. Benchmark transform, accumulation, transition, and conversion costs at the required latency, rather than quoting a historical synthesis-kernel speedup.

The alternative is a separate offline-rendering optimization, but that would not solve live eight-player playback.

### 13. Real-time multicore player rendering

Apple's Audio Workgroups coordinate real-time threads that share an audio deadline. JUCE exposes workgroup integration and a recommended parallel-thread count, while noting that this recommendation does not account for current system load.[^12][^13] Host and wrapper support must be established for the actual VST3 deployment; documentation for Audio Units is not proof that every VST3 host provides the same context.[^14]

**Taikor application.** Use a persistent, preallocated worker arrangement over independent player state. Build per-player event streams deterministically, render into preallocated buffers, and mix in the original player order. Avoid a generic task pool, per-callback thread creation, logging, and allocation. Preserve shared random-number/event sequencing when moving work out of the original serial loop.

Dispatch only when workload size justifies it. Small blocks, one active player, multiple plugin instances, and host oversubscription can reverse a parallel speedup. Include the calling audio thread in the resource budget rather than allocating one worker per logical CPU.

Correctness needs an explicit ownership/completion protocol. If a worker is late, the audio thread cannot safely render the same mutable player state concurrently as a supposed fallback. Change scheduling between completed jobs; do not introduce a race to conceal a missed deadline.

Measure worker wake-up time, longest player job, join time, callback percentiles, and total worker CPU. Fewer dropouts with the same total work is still valuable, but it is a different result from a 90% CPU reduction.

### 14. GPU and learned alternatives

GPU modal synthesis has been demonstrated for cymbals, but the cited implementation used an NVIDIA 1080 Ti and primarily 256-sample buffers at 44.1 kHz. It also approximated nonlinear behavior rather than reproducing Taikor's solver.[^15] Apple's Metal guidance warns against synchronous completion waits because they introduce pipeline latency; its illustrative timing should not be treated as a measured M1 Max audio bound.[^16]

**Assessment.** GPU throughput is not a sufficient reason to place a GPU completion dependency inside a 0.667 ms callback. Device contention, latency, precision support, state synchronization, and total energy would all need evaluation. Defer this until CPU vectorization and tail experiments have been quantified.

A 2025 DAFx paper combines analytic linear modal evolution with learned nonlinear dynamics for a proof-of-concept string model.[^17] It is an interesting longer-term research direction, not evidence of transparent substitution for all Taikor articulations and control ranges. Training coverage, stability, retrigger behavior, and worst-case inference time would be substantial additional obligations.

## Changes that should not be presented as transparent optimizations

- Cutting the default mode count to an arbitrary value such as 32 or 64.
- Truncating every tail after a fixed time regardless of its state.
- Replacing each player's independent physical response with one cloned player output.
- Disabling shell, rear-head, strain, or contact behavior without perceptual validation.
- Enabling global fast-math or converting every state to float32 without separate evaluation.
- Delaying or staggering MIDI events to flatten CPU peaks.
- Increasing the audio buffer and calling the result lower total CPU without measuring it.
- Reporting a lower callback meter after adding workers as a reduction in whole-process CPU.
- Replacing live interactive physics with static samples and testing only isolated note-ons.

A modal piano implementation provides useful precedent for making deliberate, acoustically informed simplifications while preserving musical interactions.[^18] It does not justify assuming the same simplifications are valid for a two-headed drum. The distinction is between an explicitly tested model change and an accidental loss of behavior hidden inside a performance patch.

## Sound-quality acceptance

### Exact-result track

For layout, scheduling, and unchanged-arithmetic candidates, aim for identical output and state on each supported target relative to that target's scalar reference. Cross-platform bit identity is a different requirement and is not implied. Fix random seeds, input events, initial state, parameter trajectories, compiler settings, and reference build identity.

Retain the 99 prior audio comparisons and 25 tests as historical regression coverage, not proof for new code. Expand the corpus before introducing approximation or new scheduling. Include long renders, because small state divergence can be hidden by a short attack comparison.

Required cases include all four families, every supported articulation, quiet and loud hits, center and rim positions, damping/chokes, controls moving during ringing, repeated hits on an already vibrating head, simultaneous and staggered ensembles, stereo and mono output, reverb off and on, and variable block sizes. Compare dry output as well as the complete plugin path.

An exact mismatch is a failure of the exact-result claim, not automatically evidence of audible degradation. It should trigger diagnosis and an explicit decision about whether the candidate belongs in the approximation track.

### Approximation track

Use objective checks to identify problems and choose difficult listening excerpts. Useful measurements include residual peak and RMS, time-local error around attacks and transitions, multiresolution spectra, decay slopes, modal frequency/beat behavior, and stereo correlation. Examine quiet tails after a dense passage as well as the loudest moments.

No universal residual threshold establishes inaudibility for every sound, gain setting, or listener. A long-term error average can hide a short click. Independent loudness normalization can conceal an incorrect gain change, and automatic time alignment can conceal added latency. Retain unmodified-output comparisons alongside any diagnostic alignment or level matching.

ITU-R BS.1116-3 is designed for small audio impairments and is the closer methodological fit to a transparency target. BS.1534-3/MUSHRA addresses intermediate-quality evaluation and is better suited to comparing visibly different research variants than to proving an optimized build is transparent.[^19][^20] Both remained listed as in force when consulted.

**Proposed listening protocol:** randomized, blinded reference/candidate trials, with repeatable switching, hidden references, and appropriate reference-quality listening equipment. Include trained listeners and musicians familiar with percussion. Use BS.1116-style impairment judgments for near-transparent candidates; ABX identification can provide an additional discrimination check.

Set the smallest unacceptable degradation and the statistical analysis before collecting results. Failure to find a statistically significant difference is not proof of equivalence. Report listener count, trial count, uncertainty, difficult excerpts, and any repeatable preference for the reference. Do not label an informal listening session as formal ITU compliance.

PEAQ, standardized in ITU-R BS.1387-2, can provide an additional objective indicator where an appropriate implementation is available.[^21] It does not replace listening or stability tests for an interactive nonlinear synthesizer. A candidate must also preserve the instrument's response to controls, not just match a selection of finished recordings.

### Release rule

Ship an exact candidate only after its correctness and real-workload performance gates pass. Ship an approximate candidate only after its numerical, musical-interaction, and perceptual gates also pass. Keep the established renderer available during evaluation and define fallback eligibility before integrating the new path.

Do not silently vary the physics according to instantaneous CPU load. If adaptive rendering is introduced, it needs deterministic eligibility and an explicit error policy so the same performance remains musically consistent under different background load.

## Benchmark design

### Isolate diagnosis from acceptance timing

Apple's CPU profiling guidance distinguishes call-stack profiling from investigation of pipeline and memory bottlenecks.[^22] Use hardware counters and compiler reports to explain a candidate's performance, then use separate, uninstrumented runs for the reported speedup. A profile collected while a benchmark runs is not the clean performance comparison.

Record build identity, architecture, compiler and flags, sample rate, actual block sizes, workload seed, and active feature settings. Native arm64 execution should be distinguished from translated execution when comparing host configurations; Apple documents that x86_64 applications on Apple silicon use Rosetta.[^23] The previous local benchmark was already arm64, so native execution is not an unclaimed optimization waiting to be counted again.

### Workload matrix

| Dimension | Minimum coverage |
|---|---|
| Audio configuration | Existing 48 kHz / 256 and 96 kHz / 64 anchors; expand to smaller/larger and variable blocks |
| Population | 1, 2, 4, and 8 players; one family and all four |
| Timing | Aligned hits, flams, rolls, sparse passages, and attacks after long silence |
| Dynamics | Soft through maximum velocity; damping and parameter extremes |
| Interaction | Restrikes, controls during ringing, and transitions into/out of every fast path |
| Output | Dry, wet, mono compatibility, and representative spatial settings |
| Environment | Standalone, actual VST3 host, editor open/closed, several plugin instances |
| Duration | Short reproducible fixtures plus sustained runs long enough to expose thermal and scheduling behavior |

Pair reference and candidate runs in alternating or randomized order. Use enough repetitions to estimate variation, and report the distribution rather than the single fastest result. Keep profiling, builds, and other heavy background activity out of acceptance measurements. Use an independent held-out workload set for PGO and approximation tuning.

### Required output

- Mean rendering time and CPU-seconds per audio second.
- p95, p99, p99.9 where sample counts support them, observed maximum, and deadline misses.
- Per-drum and per-state work counts, with sampling overhead separately characterized.
- SIMD occupancy or effective lane utilization for sparse and divergent workloads.
- For tail approximations: eligible mode-samples, transition counts, conversion cost, and net savings.
- For parallel rendering: main-thread elapsed time, worker CPU sum, scheduling overhead, and host behavior.
- Audio comparison results and an explicit exact/approximate classification.

Do not choose a candidate on aggregate average alone if it regresses the worst drum or causes rare attack spikes. A practical milestone is reducing all-drum eight-player deadline use with margin on the actual machine and host, followed by the lower total-CPU target.

## Implementation handoff

This is a proposed future change sequence. No production DSP, build settings, installed plugin, or tests are changed by this report.

### Phase 1: Compact exact renderer

1. Preserve the September 13 scalar renderer and deterministic corpus as the comparison baseline.
2. Prototype a compact recurrence representation for one drum; keep the original coefficients, precision, event ordering, and mode population.
3. Implement mode-lane SIMD with ordered reductions and separately compare player-lane SIMD on the ensemble fixtures.
4. Integrate the winning layout through nonlinear correction, boundary, and observation stages without changing their logical order.
5. Benchmark the full engine and plugin, not only the recurrence microbenchmark; publish exact audio/state comparisons.

**Likely edit scope:** [TaikoEngine.h](/Users/vojta/Dev/virtual-instrument-taikor/Source/DSP/TaikoEngine.h), [TaikoEngine.cpp](/Users/vojta/Dev/virtual-instrument-taikor/Source/DSP/TaikoEngine.cpp), [ModalAdvanceTests.cpp](/Users/vojta/Dev/virtual-instrument-taikor/Tests/ModalAdvanceTests.cpp), and [BenchmarkCPU.cpp](/Users/vojta/Dev/virtual-instrument-taikor/Tools/BenchmarkCPU.cpp). Treat these as architectural entry points, not instructions to replace every file.

**Exit gate:** a reproducible improvement beyond run-to-run noise on the expensive fixtures, no unacceptable sparse-case regression, and passing exact-result checks. The desired gain is not assumed in advance.

### Phase 2: State cost and safe-tail eligibility

1. Instrument contact, nonlinear, boundary, and free-evolution eligibility in a diagnostic build.
2. Measure the fraction of total time and mode-samples that can actually use an independent tail renderer.
3. Derive and test eligibility invariants and event-driven promotion, including retriggers and control movement.
4. Decide whether the measured opportunity justifies multirate or pruning work before implementing either globally.

**Likely edit scope:** engine state/traversal code and benchmark reporting. Retain the existing production CPU meter's lightweight behavior; detailed diagnostic counters need not become permanent per-sample overhead.

**Exit gate:** quantified potential savings after transition/conversion overhead, plus a correct state-continuation design. If eligible tails are rare, redirect effort to coupling and vectorization rather than weakening the eligibility rules to force a result.

### Phase 3: One bounded approximation

1. Choose either selective multirate free tails or bounded late-mode retirement, not both initially.
2. Keep the full renderer as reference and as the path for ineligible states.
3. Evaluate state behavior, audio differences, latency, and the expanded held-out corpus.
4. Conduct blind listening tests before enabling the candidate by default.

**Likely edit scope:** a separately isolated tail-rendering implementation, explicit transitions in the engine, and new dedicated quality/transition tests. A new module is preferable to scattering approximation conditions throughout unrelated code.

**Exit gate:** demonstrated total-work savings and no unacceptable audible or interactive degradation under the predefined protocol. A strong speedup does not waive the quality gate.

### Separate track: Eight-player scheduling

Prototype persistent real-time workers only after the state ownership and deterministic mix contract are written down. The relevant integration points are [PluginProcessor.cpp](/Users/vojta/Dev/virtual-instrument-taikor/Source/PluginProcessor.cpp), [PluginProcessor.h](/Users/vojta/Dev/virtual-instrument-taikor/Source/PluginProcessor.h), and the engine's player rendering boundary. Changes to [AudioCpuMeter.h](/Users/vojta/Dev/virtual-instrument-taikor/Source/DSP/AudioCpuMeter.h) should distinguish callback load from any added aggregate worker measurement.

**Exit gate:** lower high-percentile callback time in the actual host, bounded ownership and shutdown behavior, no race-based fallback, and separately reported total CPU. This track can improve usability before a 90% total-work reduction is attainable.

## Evidence limits and decision

The local performance and structural facts are tied to the September 13 implementation and its preserved validation results. There is no new profile, prototype speedup, or listening-test result for the proposals above. The earlier benchmark matrix used relatively short repeatable fixtures; sustained host performance still requires separate evaluation.

Several important research papers are old because the core modal-synthesis techniques are established. Their age does not invalidate the mathematics, but their measured speedups on historical CPUs, consoles, or GPUs are not transferable to an M1 Max. Current Apple, Arm, LLVM, and JUCE documentation supports implementation choices, while the exact project dependency revision and host capabilities must be checked during implementation.

The most consequential unresolved question is how much of Taikor's ringing state can evolve cheaply without altering feedback. Until that is measured, it is not responsible to assign a predicted 70%, 80%, or 90% gain to a hybrid renderer.

**Decision: implement and benchmark the compact exact SIMD renderer first. Then let measured state costs determine whether the next major investment is adaptive tails, coupling reduction, or multicore scheduling.** This offers the strongest route to substantial savings while protecting the accepted sound and the instrument's response to playing.

## Sources

External documentation and source status consulted September 14, 2026. Conference dates and document title pages take precedence over search-engine crawl or upload dates. Local evidence is not a public URL.

[^1]: Taikor project. [CPU ensemble and meter validation report](/Users/vojta/Dev/virtual-instrument-taikor/Docs/cpu-ensemble-and-meter-2026-09-13.md), September 13, 2026. Local implementation and benchmark record. Associated [matrix summary](/Users/vojta/Dev/virtual-instrument-taikor/build-cpu-sep13/matrix/summary.json), [drum diagnostics](/Users/vojta/Dev/virtual-instrument-taikor/build-cpu-sep13/drum-physics.csv), and [test log](/Users/vojta/Dev/virtual-instrument-taikor/build-cpu-sep13/tests.log). Historical results, not rerun for this report.

[^2]: Arm. [Structure of arrays: Optimize SIMD code with vectorization-friendly data layout](https://learn.arm.com/learning-paths/cross-platform/vectorization-friendly-data-layout/a-more-complex-problem-revisited/). Living learning-path documentation; publication date not stated. Supports contiguous layout and vector occupancy considerations.

[^3]: Arm. [Neon](https://www.arm.com/technologies/neon). Living architecture overview; publication date not stated. Floating-point vector lane widths and available implementation approaches.

[^4]: LLVM Project. [Auto-Vectorization in LLVM](https://llvm.org/docs/Vectorizers.html). Living documentation; version/date not fixed by this URL. Loop vectorization, optimization remarks, and floating-point reduction constraints.

[^5]: LLVM Project / Clang. [ThinLTO](https://clang.llvm.org/docs/ThinLTO.html). Living documentation; version/date not fixed by this URL. Cross-module optimization and build architecture.

[^6]: LLVM Project / Clang. [Clang Compiler User's Manual](https://clang.llvm.org/docs/UsersManual.html), sections on floating-point behavior and profile-guided optimization. Living documentation; version/date not fixed by this URL. Fast-math, contraction, and PGO semantics.

[^7]: Stefano Zambon. [Accurate Sound Synthesis of 3D Object Collisions in Interactive Virtual Scenarios](https://www.di.univr.it/documenti/AllegatiOA/allegatooa_18655.pdf). PhD thesis, University of Verona, April 23, 2012. Chapter 5, especially sections 5.2-5.4 and Table 5.1. Multirate/state-switching methods and limitations.

[^8]: D. Brandon Lloyd, Nikunj Raghuvanshi, and Naga K. Govindaraju. [Sound Synthesis for Impact Sounds in Video Games](https://www.microsoft.com/en-us/research/wp-content/uploads/2016/10/5.pdf). ACM I3D, February 18-20, 2011, pp. 55-62. Sections 3.4.3-3.4.4. Multirate modal mixing and SIMD in a different, primarily feedforward synthesis architecture.

[^9]: K. van den Doel, D. K. Pai, T. Adam, L. Kortchmar, and K. Pichora-Fuller. [Measurements of Perceptual Quality of Contact Sound Models](https://www.icad.org/Proceedings/2002/vandenDoelPai2002.pdf). International Conference on Auditory Display, July 2-5, 2002. Human-subject evidence for modal selection on the tested objects.

[^10]: Michele Ducceschi, Stefan Bilbao, and Craig Webb. [Real-Time Modal Synthesis of Nonlinearly Interconnected Networks](https://www.pure.ed.ac.uk/ws/portalfiles/portal/377936620/Bilbao2023DAFxRealTime.pdf). DAFx, September 4-7, 2023. Energy quadratisation, structured updates, and real-time implementation. Also indexed in the [DAFx archive](https://www.dafx.de/paper-archive/details/QBDKQDp2t_iVFXWNnffSrg).

[^11]: Nicolas Bonneel, George Drettakis, Nicolas Tsingos, Isabelle Viaud-Delmon, and Doug James. [Fast Modal Sounds with Scalable Frequency-Domain Synthesis](https://www-sop.inria.fr/reves/Basilic/2008/BDTVJ08/FastModalSounds.pdf). ACM Transactions on Graphics / SIGGRAPH, 2008. Sections 4 and 7. Frequency-domain acceleration, quality limitations, and evaluation configuration.

[^12]: Apple. [Understanding Audio Workgroups](https://developer.apple.com/documentation/audiotoolbox/understanding-audio-workgroups?changes=_9) and [Meet Audio Workgroups](https://developer.apple.com/videos/play/wwdc2020/10224/), WWDC 2020. Common-deadline real-time scheduling and auxiliary rendering threads.

[^13]: JUCE. [juce::AudioWorkgroup Class Reference](https://docs.juce.com/master/classjuce_1_1AudioWorkgroup.html). Living `master` documentation. Thread membership and `getMaxParallelThreadCount`, including the system-load caveat.

[^14]: JUCE. [juce::AudioProcessor Class Reference](https://docs.juce.com/master/classjuce_1_1AudioProcessor.html), `audioWorkgroupContextChanged`, block-size, and latency contracts. Living `master` documentation. Apple. [Porting your audio code to Apple silicon](https://developer.apple.com/documentation/apple-silicon/porting-your-audio-code-to-apple-silicon?changes=l_1_5). Living documentation. Audio Unit workgroup integration; not a guarantee for a particular VST3 host.

[^15]: Travis Skare and Jonathan Abel. [Real-Time Modal Synthesis of Crash Cymbals with Nonlinear Approximations, using a GPU](https://www.dafx.de/paper-archive/2019/DAFx2019_paper_48.pdf). DAFx, September 2-6, 2019. GPU implementation, hardware/buffer configuration, and approximated nonlinear behavior.

[^16]: Apple. [Tuning Hints](https://developer.apple.com/documentation/metalperformanceshaders/tuning-hints?changes=_3_1&language=objc). Metal Performance Shaders living documentation; publication date not stated. CPU/GPU pipeline and synchronous-wait cautions, not a current-device real-time latency guarantee.

[^17]: Victor Zheleznov, Stefan Bilbao, Alec Wright, and Simon King. [Learning Nonlinear Dynamics in Physical Modelling Synthesis using Neural Ordinary Differential Equations](https://www.dafx.de/paper-archive/2025/DAFx25_paper_37.pdf). DAFx, 2025. [Author demonstration page](https://victorzheleznov.github.io/dafx25/). Proof of concept for a nonlinear string, not a validated Taikor replacement.

[^18]: Balazs Bank, Stefano Zambon, and Federico Fontana. [A Modal-Based Real-Time Piano Synthesizer](https://home.mit.bme.hu/~bank/publist/taslp10.pdf). IEEE Transactions on Audio, Speech, and Language Processing, 18(4), pp. 809-821, May 2010. Modal instrument implementation and model/interaction tradeoffs.

[^19]: International Telecommunication Union. [ITU-R BS.1116-3: Methods for the subjective assessment of small impairments in audio systems](https://www.itu.int/rec/R-REC-BS.1116-3-201502-I/en). February 2015; listed as in force. Near-transparency subjective assessment.

[^20]: International Telecommunication Union. [ITU-R BS.1534-3: Method for the subjective assessment of intermediate quality level of audio systems](https://www.itu.int/rec/R-REC-BS.1534-3-201510-I/en). October 2015; listed as in force. MUSHRA's intended quality range.

[^21]: International Telecommunication Union. [ITU-R BS.1387: Method for objective measurements of perceived audio quality](https://www.itu.int/rec/R-REC-BS.1387). Current listed revision BS.1387-2, May 2023, in force. Objective perceived-quality measurement.

[^22]: Apple. [Optimize CPU performance with Instruments](https://developer.apple.com/videos/play/wwdc2025/308/). WWDC 2025. Profiling, CPU counters, and bottleneck-directed optimization.

[^23]: Apple. [Porting your macOS apps to Apple silicon](https://developer.apple.com/documentation/Apple-Silicon/porting-your-macos-apps-to-apple-silicon?changes=lates_7_6). Living documentation; publication date not stated. Native/universal binaries and Rosetta execution.
