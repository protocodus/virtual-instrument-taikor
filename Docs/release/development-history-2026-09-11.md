# Taikor development history to 2026-09-11

Moved verbatim from the README's version history on 2026-10-17, when the README
ledger was shortened to the portfolio format. None of these dated changes was
ever a tagged release; all of them are part of the unreleased 1.0.0.

## 2026-09-11

- Added build numbers to every platform package and the final distribution
  artifact; the download link records the matching build and source commit.
- Added independent family skin and bachi profiles, fixed uneven-head modal
  axes, the shared head/air cavity, and reciprocal shell boundary coupling.
- Added rear-head strokes through **CC18** (0–63 front, 64–127 rear), captured
  per hit including delayed ensemble companions. **CC121** restores the front.
- Updated physical readouts and live state transfer for the shared cavity;
  retained the vintage Japanese editor and existing host parameter IDs.
- Aligned the residual-noise filters with their physical excitation frequencies,
  preventing misplaced noise from overwhelming a drum's body partial.

## 2026-09-09

- **Redesigned the instrument in indigo, ivory and vermilion**, with original
  Japanese woodblock artwork, illustrated drum selectors, a high-contrast
  playing grid and separate body, stroke and output panels. Performance controls remain
  on one screen, with editable values, keyboard focus and host automation preserved.

## 2026-09-08

- **Resolved the head to seventy-six modes** — 152 resonators against the
  former 40 — with the bachi contact still solved against the calibrated
  twenty entries and the further fifty-six driven one-way beneath the
  unchanged statistical continuum. Their share sits 20–25 dB under a Don on
  the ō-daiko; the sample peak of a full Don moves by 0.7 dB and the 500 ms
  RMS by under 0.1 dB.
- **Gave soft bachi a felt contact law.** Below Bachi Hardness 50 % the
  contact exponent rises from the Hertz 1.5 to 2.5, with the stiffness pinned
  so a neutral stroke keeps its contact time; a felt beater's contact now
  shortens as *v^−0.43* rather than *v^−0.2*. The hard bachi the family is
  played with is untouched.
- **The far head reaches the pair by its own path.** On every cavity-coupled
  axisymmetric branch the resonant head's share arrives down the body and
  round the rim, later and more spread, carried as the mode's quadrature
  residue.
- **Each resolved entry carries its own tension ripple** beside the shared
  follower, so a mode is stiffened by its own motion. A uniform ripple was
  measured to pump the nagadō-daiko's (0,1) and (2,1) modes parametrically
  and was not kept.
- **Ensemble companions are distinct drums.** Above Variation 0 each player
  gets its own head tension, hide, shell and damping offsets, hashed from its
  seat.
- **Every one of the five sits behind a review switch**
  (`TaikoEngine::setRealismFeatures`); with all five off the engine renders
  sample-identically to the previous release. They ship on: the listening
  test recorded in `Docs/decisions.md` chose all five together.
- **The panel is readable and shorter.** The control decks sit on a
  translucent washi wash, every label and readout is larger, Velocity Curve
  is a three-position Curve switch, Ensemble Size an eight-position Players
  switch, and Resonant Head, Air Coupling and Shell Resonance stay host
  parameters without a knob. Each grid pad shows its strike map beside its
  note name, both sized to the pad.
- **Added two switched candidates for the wood and the hide**, both off: the
  head's boundary shear at the rim driving the shell on every stroke, and a
  per-drum hide thickness variation replacing the single 5.1-cent split. The
  first fails its own screen against real captures and the reason is recorded;
  neither has been to a listening test.
- **Two candidates for the head's noise character were rendered and rejected
  by ear**, and stay in `TaikoEngine::setRealismFeatures` switched off: the
  continuum handed off above the whole resolved bank, and the bachi contact
  solved against every resolved entry. See `Docs/decisions.md` for what the
  measurements said and why the listening test overruled them.

## 2026-09-06

- **Added CLAP for macOS, Windows and Linux**, alongside VST3, the macOS Audio
  Unit and the standalone instrument.

## 2026-09-05

- **The upper spectrum follows the solved strike force.** Causal modal-response
  integration replaces nominal Hertz spectral shading and positive-only noise
  injection. Sustained contact can now cancel high-mode excitation. The accepted
  tension/bending weighting remains, while the attack and upper tail can change
  in existing sessions. The reference impulse still anchors overall voicing.

- **Removed the extra airborne click.** The hit now reaches the output through
  the contact-driven drum, its statistical upper spectrum and the applicable
  tack source. The separate differentiated-force layer and its unsupported
  patch-size filter were removed after the user preferred the comparison
  without them. This changes the attack of existing sessions.

- **Corrected the high-mode statistical impulse response.** Relative band
  levels now integrate modal population and force-to-displacement response
  using each head's tension/bending dispersion. The original crossover band,
  filter topology and output protection remain in place; higher bands change
  level, so existing sessions can sound different above the resolved modes.
  This is a model correction, not a recording-derived calibration.

- **Added Ensemble Size and Ensemble Variation.** Each played row can trigger
  1–8 independent drum models with up to 30 ms of companion timing spread and
  bounded strike-placement scatter. The lead stays on the MIDI timestamp,
  all members share one protected output stage, and Size 1 retains the solo
  behavior. The two controls append host parameters without moving existing
  automation; older sessions load with Size 1.

- **Kept repeated hits alive at Humanise 0.** A small impact-speed and
  contact-duration floor now varies successive strokes even with Stick Noise
  at zero. Performer identities also vary this tightest gesture, contact noise
  is reseeded on every hit, and the sequence uses both halves of its 64-bit
  stroke counter. Reset still reproduces the same performance; Panic and
  natural silence preserve the advancing sequence. Existing sessions that
  relied on machine-identical Humanise-zero hits now receive this variation.
  Selected as **B by ear** in the level-matched
  [repeated-stroke comparison](../decisions.md#2026-09-05--repeated-strokes-retain-contact-variation).

## 2026-08-19

- **Gave the wooden shell a perspective.** Mic Distance moved the head's near
  field and the head's continuum but left the body at one fixed level, so
  backing the pair off thinned the drum around a shell that never receded. A
  ring mode is now read the way a membrane mode is — an evanescent term at its
  own circumferential wavenumber `n/R` over the wall-to-capsule path, plus the
  same proximity lift and propagating share — taken as a ratio against each
  drum's own capsule distance at the factory Mic Distance, so every factory
  preset renders exactly as it did. At 40 cm the ō-daiko's lowest ring mode is
  3.6 dB down and the okedo's 16.0 dB, while the top of either bank moves under
  a decibel. Chosen by ear from a six-way listening test; see
  [Docs/decisions.md](../decisions.md).
- **Adjudicated fifteen proposed acoustic and performance mechanisms** against
  the shipping engine; none was released. Two README claims were corrected in
  the process: Body Depth *does* move the decay of the axisymmetric pair, so
  what the model lacks is a loss belonging to the enclosed *air* rather than
  authority over the tail; and that missing loss factor is 2.8e-4 to 6.6e-4 by
  Kirchhoff's boundary-layer result, worth at worst 1.49 % of the decay over
  four octaves, eleven Body Depths and three Air Couplings.

## 2026-08-17

- **Removed the unsupported one-way shell copy from head-only Don, Ka and Tsu.**
  Their normal contact force already loses energy through the resolved
  shell-dependent boundary, but it was also being sent after the reciprocal
  solve into six shell oscillators that contributed no sensing or compliance.
  On the light okedo that fixed 191 Hz ladder overwhelmed the head by 26.8 dB
  on Ka. The wooden bank now belongs to Don Rim alone; okedo Ka falls by about
  20 dB over 5–30 ms and Tsu by about 10 dB. No replacement EQ, radiation
  scalar or fitted gain was added.
- **Preserved the unresolved head's calibrated RMS** when live Pitch, Tension
  Mod or a structural rebuild moves a continuum filter. Ordinary attack stretch
  moves only 0.009–0.016 dB, while two-octave automation no longer inherits up
  to several decibels of passband gain.
- **Made the zero-azimuth analytic observer rank the detuned cosine pole** the
  renderer actually builds rather than its unsplit parent. The four sounding
  pitches now close on exact heard octaves at 59.747 / 119.495 / 238.990 /
  477.979 Hz.

## 2026-08-16

- **Corrected the analytic `sin^1.5` Hertz reference impulse** from an
  unrelated `sqrt(sin)` integral to 1.7480383695280799, then bounded each
  contact's stochastic-continuum observation by the direct rigid-target
  Hunt–Crossley squared-force integral. A factory chū Don remains uncapped at
  about 0.953 of the limit; the hostile shime Don's roughly 446× request is
  bounded, while 80–500 ms shime tails and dense-gesture onsets stay intact.
- **Stopped treating the two heads as coincident radiators.** The axisymmetric
  loss now includes their finite body-depth separation through
  `b² + r² + 2br·sinc(ωL/c)`, preserving the former zero-depth limit and a
  non-negative passive power. The factory shime's opposing lower branch moves
  from 1.685 to 0.851 s T60 without changing its first-80-ms hit level.
- **Gave each axisymmetric membrane pair its own finite-column cavity factor**
  instead of reusing the (0,1) result. The factory higher upper branches move
  down by 12.2 / 4.8 / 1.3 cents while the tuned pair is unchanged, and the
  higher solves run once on the octave search's winning drum, adding no
  per-sample work.
- **Gave each unresolved-head octave a wavelength-dependent microphone-distance
  law.** All five bands keep their exact factory-position level, while a
  3–40 cm move now attenuates the factory ō-daiko's band targets by roughly
  9.1 / 14.2 / 18.6 / 21.1 / 22.1 dB from low to high instead of one flat gain.
- **Removed the scalar observer's false omnidirectional floor** from every
  non-axisymmetric head mode: the angular factor now multiplies the propagating
  term as well as the local field, so all 32 resolved multipole resonators and
  the pitch readout go to zero on their nodal azimuths.
- **Added explicit polar strike placement.** A ±180° Strike Azimuth parameter
  reaches modal drive, reciprocal contact and the two-microphone observer, with
  CC16/CC17 for sample-accurate absolute azimuth and radial overrides and CC121
  to clear both. At the 0° default the former angle-jitter expression is
  retained exactly.
- **Removed Strike Position's clamped dead travel:** each articulation now
  spends its full bipolar range moving continuously between its own written
  position and the same centre/rim endpoints it already reached.
- **Made structural automation continuous in physical two-head coordinates.**
  Axisymmetric rebuilds preserve batter/rear displacement, velocity and pending
  force while redistributing a live Tsu loss in the new eigenbasis; a bachi
  still touching the head receives matching sensing and force projections.
- **Removed two synthetic-noise shortcuts from the hit.** The unresolved head
  now uses the same Hertz contact spectrum as the resolved modes instead of a
  one-pole brightness law, and the direct airborne click differentiates only
  the solved normal force rather than differentiating roughness a second time
  into a flat Nyquist shelf.
- **Removed three sources of cheap hit noise** without lowering the drum's
  resonant object: rope-laced okedo and shime no longer emit byō tack chatter,
  the remaining tacked drums confine it to 2.6–9 kHz with exact host-rate
  normalisation, and Humanise scatters equal bounded radial and tangential
  distances instead of swinging edge strokes through metre-scale arcs.
- **Added four stable Performer identities** for honest multi-instance
  ensembles. P2–P4 salt only the variation Humanise already owns — contact
  point, speed and duration — so layered instances change the performed
  resonant hit rather than merely adding noise. P1 preserves the established
  sequence at introduction; no identity hides a
  pitch, drum-size or gain offset.
- **Appended the automatable Velocity Curve** as the 25th host slot, without
  moving any established parameter ID or automation index. `pow(v, exp2(c))`
  gives Soft 100 / Linear / Hard 100 controller response; Linear takes an exact
  identity branch and is the default. It calibrates the player and controller
  before Velocity Depth rather than inventing an acoustic effect.
- **Added a dependency-free controlled-capture inventory preflight.** Its
  strict TSV contract checks synchronized calibrated channels, one complete
  3×3×2 mobility core, the low-speed Ka/Tsu/Rim strike matrix, held-palm and
  simultaneous head/hoop metadata, repeated measured-speed bins, matching
  fixture state, raw traction maps, and bachi identity, mass and tip
  metrology. It validates acquisition coverage only — it neither reads samples
  nor blesses the current contact discontinuities.
