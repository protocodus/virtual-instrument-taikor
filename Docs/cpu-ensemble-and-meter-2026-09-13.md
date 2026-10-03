# Ensemble CPU, live meter, and the 90% target

## Result

This pass reduces active DSP render time by approximately 12-24% against the
September 12 optimized build. It does **not** achieve a 90% reduction, and the
worst eight-player workload still misses the real-time budget on this M1 Max.
All 99 before/after audio comparisons are bit-identical. No modes, players,
precision, sample rate, tails, or physical interactions were removed.

## Implemented

- Place each mode's audio-rate fields in an aligned 128-byte prefix. Cold
  construction data no longer shares all of those hot cache lines. `Mode` grows
  from 328 to 384 bytes; this trades some memory for better traversal locality.
- Advance every free recurrence in a simple first pass, then apply the original
  nonlinear division only to the cached parametric subset. The division consumes
  the exact stored free displacement, without advancing history twice. The shell
  coupling and stereo observation order remain unchanged.
- Cache complete per-drum panel measurements by parameters, live controllers,
  sample rate and realism-feature mask. The editor previously reconstructed a
  physical modal bank on every timer refresh, even with unchanged controls.
- Add an audio-callback CPU indicator, smoothed over approximately 250 ms, with a
  recent decaying peak. Peaks at or above 100% are red. The accessible description
  includes the number of callbacks exceeding their deadlines since preparation.
  Measurement uses two clock reads per callback and lock-free publication, not
  per-sample timers, logging, or locks.
- Add reproducible individual-drum and individual-drum/eight-player fixtures to
  `TaikorBenchmarkCPU`, plus modal-equivalence, CPU-meter and readout-cache tests.

The percentage displayed is elapsed callback time divided by its audio duration.
It is **not** total-computer CPU utilization. A 100% reading consumes the entire
audio deadline. Average load below 100% does not guarantee glitch-free peaks.

## Why some drums cost more

At the measured default parameters and 48 kHz, all four drums contain 296 modes,
74 strain-contributing modes, and 54 active boundary ports. The difference here
is not a different mode count. The instruments have different physical dynamics:

| Drum | Diameter | Nominal Don contact | Model readout tail |
| --- | ---: | ---: | ---: |
| O-daiko | 150 cm | 0.451 ms | 5.856 s |
| Nagado-daiko | 78 cm | 0.264 ms | 2.877 s |
| Okedo-daiko | 40 cm | 0.140 ms | 1.730 s |
| Shime-daiko | 30 cm | 0.145 ms | 2.280 s |

These are model diagnostics at default parameters and velocity 0.85, not
measurements of real acoustic drums. Contact duration, envelope decay, and
state-dependent tension/continuum updates change which expensive work stays
active. The profile identifies modal rendering as the dominant cost, with strain,
boundary coupling, and retuning also substantial; it does not isolate an exact
causal percentage for each difference between drum families.

Each ensemble player owns independent physical drums and currently renders
serially on the same audio thread. Eight players on one drum mean 2,368 resonators.
With all four drums ringing, there are 32 physical banks and 9,472 resonators:
about 455 million recurrence updates per second at 48 kHz, before the other DSP.

## Matched benchmark results

Native arm64 Release, Apple M1 Max. The matrix used 48 kHz / 256 samples and
96 kHz / 64 samples, 0.6-second renders including a 0.2-second tail, 0.1-second
warmup, three repetitions, and two alternating reference/candidate rounds.
Builds and tests did not overlap the clean timing runs.

48 kHz / 256 samples, percentage of one real-time audio-thread budget:

| Workload | Before | After | Relative reduction |
| --- | ---: | ---: | ---: |
| Dense full plugin | 41.14% | 31.65% | 23.1% |
| Four players, all drums | 165.60% | 129.70% | 21.7% |
| Eight players, all drums | 335.79% | 258.77% | 22.9% |
| Eight players, O-daiko | 69.11% | 60.71% | 12.2% |
| Eight players, Nagado-daiko | 71.22% | 62.56% | 12.2% |
| Eight players, Okedo-daiko | 61.61% | 53.06% | 13.9% |
| Eight players, Shime-daiko | 60.75% | 51.49% | 15.3% |

At 96 kHz / 64 samples, eight-player individual-drum averages still range from
99.4% to 116.1%; all four drums cost 482.5%, down from 628.7%. Their p95 block
times still exceed the 0.667 ms deadline. The instrument is not generally
eight-player real-time safe at that configuration. No idle or reverb-only CPU
improvement is claimed, and these desktop timings are not scheduling guarantees.

## Validation and artifacts

- 96 bit-identical capture comparisons across the expanded matrix.
- Two bit-identical initial eight-player stress comparisons.
- One bit-identical 14-second automation render with a 12-second tail.
- 288 exact modal-step comparisons, covering nonlinear/linear, historical rim,
  partial/empty banks and shared-cavity indexing at 8/48/192 kHz.
- 576 existing strain/boundary state comparisons passed.
- All 25 CTest tests passed. The processor tests were also run with editor
  snapshot generation, and the resulting CPU indicator was visually inspected.
- Complete native AU, VST3, CLAP and Standalone build passed.

The reference DSP source/library and binaries are in `build-cpu-sep13/reference/`.
The reference CPU harness was rebuilt with the new drum fixtures against the
preserved old DSP library and headers. Plugin comparisons use the preserved
pre-change executable. Results and captures are in `build-cpu-sep13/matrix/`,
`build-cpu-sep13/ensemble8-first/` and `build-cpu-sep13/long-tail/`.
The screenshot is `build-cpu-sep13/cpu-indicator.png`; logs and profile data are in
the same `build-cpu-sep13/` directory.

Example repeat, using a fresh output directory:

```sh
python3 Tools/BenchmarkCompare.py \
  --baseline build-cpu-sep13/reference --candidate build-cpu-sep12/native \
  --output build-cpu-sep13/repeat-ensemble8 --program dsp --case ensemble8 \
  --config 48000:256 --seconds 1 --tail .3 --warmup .1 --repeats 3 --rounds 2
```

Individual-drum cases are `drum0` through `drum3`, or `drum0-8` through `drum3-8`
for eight independent players. The ordinary `ensemble8` case plays all four drums.

## What the 90% target requires next

A 90% total-CPU reduction is a 10x throughput improvement, not another small
loop cleanup. The current profile does not justify promising that from caching.

For glitch-free eight-player playback, the next architecture to investigate is
preallocated real-time worker rendering of independent players, with deterministic
mix order and proper host/audio-workgroup scheduling. This can reduce callback
completion time by using more cores, but does not reduce total CPU work by 90%.
It is not implemented in this pass.

For substantially lower total CPU and heat, investigate a more deeply vectorized
modal representation or a separate reduced/hybrid ensemble model. The latter
would change the model and needs explicit quality acceptance, full dynamic
control/late-restrike comparisons and listening tests before becoming a default.
No such quality tradeoff has been silently enabled here.
