# Taikor — decision log

Directions chosen by ear, recorded per the A–Z listening-test convention in
the repository's `CLAUDE.md`. A choice made by ear is recorded as made by ear,
never written up as though a measurement had settled it.

## 2026-09-05 — remove the additional airborne click

**The question.** The attack still sounded artificial. Does it improve when
the additional differentiated-force click is removed?

| Letter | What it was | Outcome |
| --- | --- | --- |
| A | Current engine, including the previously preferred continuum weighting | baseline |
| C | Same engine with only the extra normal-force-derivative click removed | **chosen** |

**The comparison.** Identical 12.4-second, 32-hit scores across the four
drums, with four Don and four Edge strokes per drum at velocities 0.4, 0.6,
0.8 and 1.0. Humanise was zero, Stick Noise was 0.35, Drive was zero, Low
Cut was off, Performer was P1 and Ensemble Size was 1. Both used the same
reset state, 44.1 kHz sample rate and 256-sample blocks. Whole-file stereo
RMS was matched to −30 dBFS using constant gain only: +15.870748 dB for A
and +15.928435 dB for C. Sources, audio and key are retained in the working
directory `build/hit-realism-review/`.

**Verdict: C, by ear.** The user chose C when asked which attack sounded
more natural. This is a listening preference, not a measured taiko match.

**What it licensed.** Remove the separate normal-force-derivative pressure
layer and its unsubstantiated contact-patch low-pass. Keep the contact
mechanics, resonant head and shell, accepted continuum weighting, and tack
source with its existing propagation intact. The final second of the A/C
renders is sample-identical. This does not claim that real bachi produce no
acceleration sound; that source needs its own geometry and radiation model.
The force-history continuum candidate is a separate comparison, with C
as its accepted baseline.

## 2026-09-05 — statistical upper spectrum follows impulse displacement

**The question.** Does the derived modal impulse-displacement weighting
improve the upper spectrum compared with the previous fixed rising tilt?

| Letter | What it was | Outcome |
| --- | --- | --- |
| A | Shipping engine at `1e20b319`, unchanged | baseline |
| B | Band levels from Weyl modal count and squared impulse displacement under tension and bending | **tentatively preferred** |

**The comparison.** Identical 12.4-second, 32-hit scores across all four
drums: four Don and four Edge strokes per drum, at velocities 0.4, 0.6,
0.8 and 1.0. Humanise, Stick Noise and Drive were zero, Low Cut was off,
Performer was P1 and Ensemble Size was 1. Both used the same reset state,
44.1 kHz sample rate and 256-sample blocks, with natural tail overlap.
Whole-file stereo RMS was matched to −30 dBFS using constant gain only:
+14.536871 dB for A and +15.971877 dB for B. The WAVs, raw renders,
render settings and key remain in `build/realism-review/continuum/`.

**Verdict: B, tentatively, by ear.** The user said “B is better i guess.”
This records a tentative listening preference, not a measured match to real
taiko or a comparison with other products.

**What it licensed.** Retain B's derived relative band levels in the live
engine, with the existing first-band anchor and 2HP/7LP filters. No fitted
constant is introduced by this choice. The experimental shared-cavity
solver was inactive in both renders and is outside this verdict. This is
a separate comparison from the repeated-stroke variation choice below.

## 2026-09-05 — repeated strokes retain contact variation

**The question.** Should Humanise 0 produce identical strokes, or retain
subtle differences in the physical stick/head contact?

| Letter | What it was | Outcome |
| --- | --- | --- |
| A | Previous engine, with exact repetition at Humanise 0 | baseline |
| B | Always-on impact-speed and actual contact-stiffness variation, with fresh contact texture | **chosen** |

**The comparison.** Identical 12.4-second, 32-hit scores across all four
drums, with Don and Rimshot at velocity 0.80. Humanise and Stick Noise were
zero, Drive was zero, Low Cut was off, and Performer was P1. Both started
from the same reset state at 44.1 kHz with 256-sample render blocks. Natural
tails overlapped. Whole-file stereo RMS was matched to −26 dBFS using only
constant gain: +8.425475 dB for A and +8.376360 dB for B. Review files and
their key remain in the working directory `build/stroke-variation-review/`.

**Verdict: B, by ear.** The user chose “B is better.” This records a listening
preference, not a measurement of a player's natural variation limits.

**What it licensed.** Retain the current B implementation: at Humanise 0,
incoming speed varies by up to ±0.5%, and a ±0.4% Hertz contact-time factor
changes the stiffness used by both the duration estimate and the live
collision. Humanise adds wider variation and position scatter. The authored
position at zero, MIDI onset timing and resting drum tuning stay fixed.
Successive hits advance the variation through silence and Panic; reset
replays the same performance. Stick Noise was disabled in the comparison,
so this choice specifically supports the mechanical contact variation.

## 2026-08-19 — the body moves with the pair

**The question.** Mic Distance moved the head's near field and the head's
continuum, but left the wooden shell at one fixed level — so backing the pair
off thinned the drum around a body that never receded. What should replace the
fixed level?

**The set.** Six candidates rendered against the shipping engine: Don Rim at
three Mic Distances on the ō-daiko and the okedo, level-matched on the
factory-position sections, letters carrying no hint of the mechanism.

| Letter | What it was | Outcome |
| --- | --- | --- |
| A | shipping engine | baseline |
| B | perspective, normalised per drum | **chosen** |
| C | `radiationEfficiency` on the level, re-pinned 29.5× | not chosen |
| D | perspective, normalised against one fixed 15.95 cm | not chosen; peak 0.985, past the suite's 0.95 clause |
| E | `radiationEfficiency` on the level, un-re-pinned | not chosen; the body nearly disappears |
| F | perspective, near field only | not chosen |

**Verdict: B, by ear.** The measurements rule out C, D and E on their own
terms — C needs a 29.5× fitted re-pin, D breaks a suite clause, E is
uncalibrated — but nothing measurable separates B from F, which differ only in
whether the propagating share the head carries reaches the body as well. B
keeps the construction identical to the head's; that it also sounds right is
the listener's finding, not a derivation.

**What it licensed.** A ring mode is now read the way a membrane mode is: an
evanescent term at its own circumferential wavenumber `n/R` over the
wall-to-capsule path, plus the same proximity lift and propagating share, taken
as a ratio against each drum's own capsule distance at the factory Mic
Distance. Every factory preset renders exactly as it did. At 40 cm the
ō-daiko's lowest ring mode is 3.6 dB down and the okedo's 16.0 dB, while the
top of either bank moves under a decibel.

Guarded by `testTheBodyMovesWithThePair`, which holds four clauses recomputed
from the resolved drum: the factory Mic Distance moves nothing on any drum or
ring mode; every other distance moves it monotonically, closer being louder;
the highest ring mode falls off with distance at less than half the rate of the
lowest, which is the evanescent signature rather than a taper; and a rim shot's
wooden bank recedes as the pair backs off while a head-only stroke still has no
wooden bank at all.

Note that gating the shell's *level* on `radiationEfficiency` — once named as
the missing term — is the wrong term in the wrong place: it is a far-field
power law and the pair is in the shell's near field.

## 2026-09-08 — five realism mechanisms, one switch each

**The question.** Five mechanisms were added across all four drums, each
derivable and each audible, and the model cannot say which of them improve
the instrument rather than merely change it. Which ship on?

| Letter | What it was | Outcome |
| --- | --- | --- |
| A | Shipping engine, every review switch off — sample-identical to the previous release | baseline |
| B | Extended bank: 76 resolved entries as 152 resonators, the 56 new ones driven one-way under the unchanged continuum | verdict pending |
| C | Felt contact law: exponent 1.5 → 2.5 below Bachi Hardness 50 %, pinned at the neutral stroke | verdict pending |
| D | Rear-head path: the far head's share reaches the pair down the body and round the rim, later and more spread | verdict pending |
| E | Per-entry tension ripple: each resolved entry stiffened by its own squared slope | verdict pending |
| F | Distinct ensemble drums: companions get their own tension, hide, shell and damping | verdict pending |
| G | All five together | verdict pending |

**The comparison.** Identical 12.9-second scores rendered through the
shipping `EnsembleEngine` path at 48 kHz in 256-sample blocks from one reset
state: every drum struck Don at velocities 0.5 and 1.0 and Edge at 0.85 with
the factory bachi; then Don at 0.6 and 1.0 on the ō-daiko, nagadō and shime
with Bachi Hardness 0.3, which is the only passage that exercises C; then a
short four-player phrase at Variation 0.4, the only passage that exercises
F. Humanise was 0, Stick Noise 0.35, Drive 0, Low Cut off, Performer P1.
Whole-file stereo RMS was matched to −30 dBFS with constant gain only:
+13.1759 dB for A, +12.9714 for B, +13.1538 for C, +12.9440 for D,
+13.1725 for E, +13.1937 for F and +12.7587 for G. The set and its key
were rendered to the session's working directory and sent to the user; it
is not committed.

**Verdict: pending.** The five ship on behind `TaikoEngine::setRealismFeatures`
until the user has listened; whichever letters are not preferred are to be
switched off there, not removed, and this entry updated with the choice.

**What the measurements already settled.** A uniform Berger ripple across
the bank was tried before E and rejected by measurement, not by ear: it
pumped the nagadō-daiko's (0,1) and (2,1) modes parametrically from the
(1,1) mode and moved the pad's heard pitch from 119 to 171 Hz. Solving the
contact against the full 76-entry bank was likewise measured and rejected: it
takes the ō-daiko's momentum transfer from 13 % to 3 % of the stick's and
folds the velocity law, because the calibrated contact is a modal-truncation
artefact. Neither is a candidate here.
