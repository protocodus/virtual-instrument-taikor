# Historical engine notes through 2026-09-08

This is the earlier engine description and its measurement record. The skin, bachi, cavity, rear-head and shell sections are superseded by [the September 9 implementation](physical-realism-2026-09-09.md). Numbers below refer to the earlier revision, not the current plugin.

#### The sound engine

##### The head

A circular membrane of radius *a* under tension *T* with areal density *σ* has
modes at *f(m,n) = c·λ(m,n) / 2πa*, where *c = √(T/σ)* and *λ(m,n)* is the *n*-th
zero of the Bessel function *J(m)*. Taikor resolves seventy-six such modes —
eight axisymmetric and sixty-eight with a circumferential order of up to
nineteen — as a hundred and fifty-two resonators, every degenerate pair and
every cavity-coupled axisymmetric pair taking two.

The first twenty entries — four axisymmetric, sixteen with an order — are the
bank the instrument was calibrated on, and they are the only ones the bachi
contact is solved against. The further fifty-six are driven one-way by the
solved contact force and never load the stick. That split is deliberate.
Solving the contact against the whole bank is the more complete physics, and
it was measured: it moves the ō-daiko from passing 13 % of the stick's
momentum into the head to passing 3 %, folds the velocity law and takes the
attack glide with it — because the contact as calibrated is a
modal-truncation artefact. The head under the stick is exactly as compliant
as twenty modes make it, and every constant that was pinned by ear was pinned
against that. So the extended entries add what a noise band cannot carry, the
individually audible partials, underneath a statistical continuum that stays
where it was voiced, and the contact keeps the load it was tuned with. Their
own contribution sits 20–25 dB under a Don on the ō-daiko. What closing the
split properly would need is under Known gaps.

Each of the five mechanisms added in the 2026-09-08 release — this extension,
the felt contact law, the far head's own path to the pair, the per-mode
tension ripple and the ensemble's distinct drums — sits behind a review
switch in `TaikoEngine::setRealismFeatures`, and with every switch off the
engine renders sample-identically to the release before. They ship on. Two
further switches were rendered as candidates and rejected by ear, and ship
off: the continuum handed off above the whole resolved bank rather than above
the calibrated twenty entries, and the bachi contact solved against every
resolved entry. Removing the noise bed measured closer to real captures on
every partial-balance descriptor and sounded like a struck string, because
the resolved bank underneath it is too sparse to stand alone. The listening
test and its numbers are recorded in `Docs/decisions.md`.

A treated cowhide head can also resist bending. A
[Japanese-drum diaphragm study](https://doi.org/10.1250/ast.30.348) measured a
Young's modulus around 3.5 GPa and investigated a *stretched plate* model.
Taikor includes that mechanism, with effective thickness and tension selected
for its family profiles. These are model parameters, not measurements of four
owned instruments.
A plate resists bending as well as stretching, and bending adds a term in the
fourth power of the wavenumber to *ω²*: *ω² = (T k² + D k⁴)/σ*, with the
flexural rigidity *D = E h³ / 12(1 − ν²)*. So the ratios between the modes are
not constants of the geometry. They open out with the mode's order, and they
open out further the smaller the drum and the thicker its hide.

Within this model, that term changes the spacing of the family's partials.
The following numbers are model predictions, not recorded-drum measurements:
the top of the resolved bank sits 30 cents above
where an ideal membrane would put it on the ō-daiko, 33 on the nagadō-daiko, 13 on
the okedo and 6 on the shime — because *B = D/(Ta²)* falls as the hide gets
thinner (as the cube of its thickness) and as it is pulled tighter, and both of
those go the same way up the family. Head Material moves it as hard again: a
thin synthetic film is an ideal membrane to within half a cent, and a thick hide
stretches its top mode by well over a semitone.

The stretch is taken relative to the *(0,1)* mode rather than applied
absolutely, because a drum is tuned by the pitch it sounds. A player brings the
fundamental back where it belongs with the ropes or the tacks, and what
stiffness leaves behind afterwards is the spread above it. That is also what
keeps an octave an octave: the stiffness parameter falls as the tension rises
and as the square of the radius, so the two Drum Layout constructions reach the
same tuning by different physical routes and an absolute shift would put the
keyboard out of tune with itself.

No stroke lands on the geometric centre, because every mode with a
circumferential order has *J(m)(0) = 0* and a strike at radius zero drives the
four axisymmetric modes and nothing else — a note with an attack and no body
behind it. A real taiko does the same thing if you manage to hit its exact
middle, which is why players do not: a full Don lands a hand's width in, close
enough to keep the fundamental and far enough out to wake the rest of the head.

Modes with a circumferential order come in degenerate pairs, the same shape
rotated by a quarter of its own period. A real head is never quite uniform, so
the pair sits a fraction of a percent apart and beats. That asymmetry belongs to
the hide rather than to the stroke, so it is seeded from a fixed constant: the
same drum splits the same way every time it is hit.

##### Published frequency comparison

The stiffness hypothesis has an external check in
[`Tools/BenchmarkPublishedModes.py`](../Tools/BenchmarkPublishedModes.py).
It compares eight measured frequencies from Suzuki et al.'s 2009 single-headed
nagadō-daiko experiment against the published plate theories and Taikor's
single-head frequency equation using the paper's material inputs. After
removing one common tuning offset, RMS spacing error is 44.27 cents for an
ideal membrane, 75.79–78.95 cents for the published plate theories, and
85.70–89.09 cents for Taikor's equation at the two inferred tensions.
The simpler membrane hypothesis fits this specimen better. The paper itself
reports that its measured low-mode spacing contradicted the predicted bending
effect. The open body, floor loading and actual rim condition are unresolved;
these figures assess equation hypotheses, not the full plugin or its presets.

Run `python3 Tools/BenchmarkPublishedModes.py` to reproduce the comparison.
The tool retains the paper's mode labels, setup and distinction between measured
frequencies and inferred inputs. Its tests check analytic limits and invariance
to a common tuning change. This is an external reference, not a fitted factory
profile or evidence of superiority to other instruments.

##### The air on the head

The air a mode has to move rides along with it as added mass, and lowers it. How
much depends on how much air the mode actually displaces, so the fundamental is
loaded far more than the high modes, and a light synthetic head is loaded far
more than a heavy hide. That is why a thin head sounds lower than its tension
alone predicts.

The released approximation currently uses that added mass in the frequency
solve. Propagating the same approximate coefficient through modal force,
contact mobility, cavity coupling and damping is a larger coordinate migration,
not a harmless cleanup: it materially changes the drum's attack and which
partial is heard. That complete migration is therefore kept behind the
capture/revoicing boundary described under "Known gaps" rather than
silently treating an unmeasured loading law as a calibrated mechanical mass.

##### The air inside the body

A taiko is a closed drum, and the enclosed air is a spring between its two heads.
The current pressure approximation couples only axisymmetric head modes,
weighted by their net volume displacement. Other head modes leave the total
volume unchanged, but a real cavity can still carry nonuniform pressure fields;
those transverse acoustic modes are outside this approximation.

The result is that each axisymmetric mode splits in two: a **breathing** mode
where both heads move outward together, lifted well above its uncoupled
frequency by the air spring, and a volume-preserving mode that is left roughly
where it was. On the default drum the pair lands at about 32.7 Hz and 61.4 Hz. The
breathing mode has a net volume source; other branches can radiate as spatial
multipoles. The microphone perspective and branch losses determine which one
is heard most strongly.

The spring is a column and not an infinite one. *ρc²/L* is what a cavity is
worth only while the wavelength runs away from the body, and this instrument
leaves that limit inside its own range: the body's first axial resonance is
134.5 Hz on the ō-daiko, 183.2 on the nagadō-daiko, 343.5 on the okedo and 817.7 on the
shime, every one of them inside the resolved bank of the drum it belongs to.
What each head actually drives is a rigidly terminated column of length *L/2* —
the volume-changing motion is symmetric about the midplane, so that plane
behaves like a wall — and its stiffness is *x cot x* times the lumped value,
with *x = ωL/2c*. For the lowest pair it is the same number at low frequency
and less than it as the body gets deep against the wavelength: 0.82 on the
ō-daiko, 0.77 on the
nagadō-daiko, 0.49 on the okedo — the longest body in the family relative to its
head — and 0.65 on the shime, falling to 0.05 there at full Body Depth. It has
to be solved for rather than computed, because the stiffness depends on the
frequency it sets. Each axisymmetric pair now converges on its own
factor; the three higher solves are deferred until the octave search has chosen
the final drum. In this independent-pair approximation their cavity factors
do not change the lowest pair used by that search; this does not guarantee
that higher-mode changes are perceptually pitch-neutral. The audio loop never
sees the iteration. The model reports the lowest pair's
factor for the same reason: it is an answer the drum has to converge on rather
than an expression anything can write down.

Musically this is what stops a long-bodied drum being an air spring with a hide
attached, and it is why the four drums split their lowest pair so differently:
1.88, 1.39, 1.07 and 1.08 times the fundamental going up the family. The
keyboard is an octave in the pitch each drum is heard at and it is not an octave
in the breathing branch, which steps 747 / 1730 / 1206 cents — because the branch
above the fundamental is lifted by a column whose length is a property of each
instrument rather than of a scaling.

Where the column passes its own quarter-wave the stiffness reaches zero and the
two heads stop being tied together at all. Above that the air is mass-like
rather than stiff, which is a real thing this model has nowhere to put, so the
answer there is the one an open body already gets: one axisymmetric mode,
reported twice. It takes a body longer than half its head's own wavelength to
reach, which is the same thing as the half column passing its quarter-wave. The
lowest pair of each taiko-proportioned drum stays below it, but higher radial
modes already cross it: the factory ō-daiko's four cached factors are about
0.823 / 0.610 / 0 / 0. The zeroes are a continuous truncation of this scalar
spring model, not a claim that the higher acoustic field disappears.

Radiation, mounting and material loss jointly determine the decay. On the
factory ō-daiko, the breathing branch has a T60 of about 1.02 seconds, the
lower branch about 0.57 seconds, and the slowest higher radial branch about
3.81 seconds. The editor therefore reports the longest-lived audible
axisymmetric branch rather than treating the lowest pair as the whole tail.

Turn Air Coupling all the way down and there is no split at all: the two heads
are independent, and a stroke on the batter head cannot reach the far one. The
editor then shows the same figure for both, because an open body has one
axisymmetric mode you can hear rather than two. It matters which one is named:
with the resonant head slack, the far head's mode is the *lower* of the pair,
so reading off the lower frequency reported a mode that nothing was driving.

##### The shell

The wooden body's ring modes come from the standard thin-cylinder result, so
the shell material moves their frequencies, their spacing and their Q together.
Rimshot drives those modes through a force-over-modal-mass path because the
bachi catches the hoop and body along with the head: a heavy carved log refuses
to move while a light laminated shell rings. Don, Edge and Muted hit only the
membrane. They retain the shell-dependent boundary loss, but do not feed an
unmodelled one-way copy of their normal contact force into the wooden bank.
That distinction matters most on the okedo: its stave shell takes two and a
half times as much out of the head at the rim as the ō-daiko's solid zelkova
does, which is why it is the driest of the four, without becoming a fixed tone
laid over every head stroke.

The body's Q is low, because a drum shell is a thick, short piece of wood
clamped at both ends by the hoops rather than a free bar. That matters more than
it sounds: when the shell rang longer than the head, it put a wooden pitch on
top of the drum where the body should only have been adding weight, and it left
a hand laid on the head unable to damp anything anyone could still hear — since
a hand on the head does not touch the body.

The body moves with the microphones, and is read the way the head is. A ring
mode's shape around the shell is *cos(nθ)*, so its spatial wavenumber is *n/R*,
and the same evanescent law that carries the head's shape to a close capsule
carries the wall's — over the path from the wall to the capsule rather than from
the head to it, because that is where the wood is. The same proximity lift and
propagating share sit on top of it. Mic Distance used to move the head's near
field and the head's continuum and leave the shell nailed at one level, so
backing the pair off thinned the drum around a body that never receded.

It is taken as a ratio against each drum's own capsule distance at the factory
Mic Distance, exactly as the continuum's perspective law is — per drum, because
the pair sits proportionally closer to the small heads. Every factory preset
therefore renders as it did. Moving the pair out to 40 cm now takes the
ō-daiko's lowest ring mode down 3.6 dB and the okedo's down 16.0 dB, while the
top of either bank barely moves: the evanescent rate is
*√((n/R)² − (ω/c)²)*, and a ring mode's frequency climbs as about *n²* where its
wavenumber climbs as *n*, so the high ring modes are already propagating and do
not care where the pair stands. The light stave okedo moves most, which is the
drum whose body is loudest against its own head.

What the shell still does not have is a radiation-aware *decay* — its Q is a
single lumped figure that already contains whatever the body radiates, so a
derived radiation loss cannot be added on top of it without counting it twice —
or an angular doublet. See the 2026-08-19 shell-perspective entry in
[Docs/decisions.md](decisions.md), including why gating the shell's *level*
on `radiationEfficiency` is the wrong term in the wrong place.

##### Above the modes: the head's continuum

A modal bank can only resolve so far. The mode table runs to the Bessel zeros
around *λ = 13*, which on a large drum puts the highest resolved mode a couple
of hundred hertz up — and a real head goes on having modes for another five
octaves above that. Taikor approximates this unresolved region with five
overlapping statistical bands, each carrying the head's loss law and driven
by the same contact as the resolved modes. The handoff is a selected model
boundary: some low bands on the large drums begin before genuine modal overlap.

Their relative levels now follow the statistical **impulse displacement** of a
head under both tension and bending. With dimensionless squared wavenumber
`q = (k a)^2`, Weyl's leading modal count is `N = q/4`; the head's dispersion
is `omega^2 = C q (1 + B q)`. Each mode's displacement after a force impulse
scales as `1/omega`. Integrating squared displacement over the modes in each
band therefore gives
`W = log1p((q_high - q_low) / (q_low * (1 + B*q_high)))`.
The band's relative RMS is `sqrt(W/W_first)`. A tension-dominated membrane
has flat octave displacement RMS; in the bending limit it falls as
`1/sqrt(f)`. The first-band level remains the existing reference.

This replaces a fixed rising tilt that conflated membrane and plate mode
density and omitted the frequency dependence of impulse displacement. The
[modal impulse-response formulation](https://www.dafx.de/paper-archive/2009/papers/paper_77.pdf)
is checked here against independently enumerated membrane/plate modes and
their summed squared displacement amplitudes. Statistical input/output coupling is still
approximated; this is not a measured microphone-pressure spectrum or a
reciprocal mechanical residual. The integral uses nominal band limits,
whereas the actual filters overlap.

Each octave is a serial two-pole high-pass followed by a seven-pole low-pass.
Its lower skirt therefore falls at twelve decibels per octave and its upper one
at forty-two, keeping the loud crossover region out of the four octaves above
it. The exact stationary variance of that nine-state filter is solved once when
the stroke is built, so its target RMS does not inherit the filter geometry or
where the band lies against Nyquist.

That calibration now survives a live pole move too. Pitch, Tension Mod and a
structural rebuild can widen or narrow a band's filter without creating or
destroying unresolved-head level: the engine carries the stored RMS envelope
and filter memory into the new variance coordinate, and converts later flam or
roll contacts before they join the same ringing field. A full ordinary stroke
changes by only 0.009--0.016 dB; stepped pitch automation had accumulated as
much as several decibels of filter-bandwidth gain. This is an automation
correctness fix, not a substitute for the measured high-mode mobility still
needed to finish the soft-hit attack.

The continuum supplies broadband attack energy above the resolved modes.
Earlier documentation claimed agreement within a few decibels with two real
recordings, but did not identify the recordings or retain a reproducible
comparison. That claim is unverified. Current regression tests establish
internal model behavior and stability, not a measured match to acoustic taiko.

Its level reference is the strongest resolved modal displacement response to a
unit impulse, observed at the factory microphone position. The first statistical
band retains that level anchor. Each octave also follows its own wavelength-
dependent distance gain. That is a property of the drum and its placement, and
of nothing else. Regression tests isolate the uppermost statistical band at
44.1, 48, 96 and 192 kHz and require it to remain within 2 dB; the complete
4–10 kHz response remains within 1.5 dB. A separate invariant checks all five
near-to-far gains against the coherent-patch equation and requires them to rise
strictly with frequency. Moving Mic Distance under an already-ringing head also
remaps the five stored envelopes by those gains without clearing their filter
history. Energy arriving later from a bachi that is still in contact is converted
from its trigger-time perspective to the live one as well. Every higher octave
is measured in isolation, so a lower band's skirt cannot masquerade as the whole
statistical tail.

It has to stay in its place, though, and its place is much smaller than it
looks. Left too loud it does not sit above the resolved bank, it buries it: the
sustain becomes a bed of noise with the drum's pitched ring somewhere
underneath, which measures as a body and does not sound like one. In the model
voicing, the sustained low-mid is primarily resolved modal energy; the
continuum supplements it. This balance still needs controlled acoustic captures.

The upper spectrum now follows the **actual contact-force history**. Each band
uses 64 quadrature frequencies distributed through its Weyl-weighted modal
population. A complex response at each frequency obeys
`dz/dt = (-sigma + i*omega)*z + F(t)`, integrated exactly for the force held
over one sample. The band's squared amplitude follows the weighted sum of
`|z|^2`. Later force can reinforce or cancel earlier excitation, just as it can
in a forced mode; the former positive-only `F^2` envelope could only add to it.
This is a causal approximation of the
[forced modal response](https://www.dafx.de/paper-archive/2009/papers/paper_77.pdf),
so the attack starts during contact rather than waiting for the stick to leave.

That distinction matters when a moving head keeps the bachi in contact longer
than the nominal Hertz estimate. The previous model timed its noise injection
from the solved collision but shaded its spectrum using the estimated pulse.
Those two descriptions could disagree by several milliseconds. The new model
uses one force history for both timing and spectral cancellation. Its absolute
level still uses the nominal reference pulse's impulse as a compatibility
anchor; it is not a measured pressure calibration. The finite quadrature
approximates a band rather than resolving all its physical modes. For 24
captured factory Don/Edge contacts, 64 nodes agreed with dense integration to
within 0.001 dB in bands within 30 dB of each attack's strongest continuum band.
That checks numerical convergence for those contacts, not acoustic realism or
every possible control setting.

Separate strikes still contribute independent statistical powers to one shared
field. Its existing position weighting and microphone correlation remain:
shorter wavelengths become less correlated between the two microphones.
Distance, palm damping and filter-variance changes scale both the shared field
and the active contacts' response states, keeping subsequent cancellation in
the same amplitude coordinates.

The attack comes from the contact-driven drum and its statistical upper modes.
The former additional differentiated-force click has been removed following a
[listening comparison](decisions.md#2026-09-05--remove-the-additional-airborne-click).
Its nominal 3.5 kHz low-pass was not derived from a millimetre-sized contact
patch: that patch is acoustically compact at these frequencies. A separate
bachi [acceleration-noise model](https://www.cs.cornell.edu/projects/Sound/impact/)
would need the stick's radiating geometry and motion. Stick/hide roughness still
excites the resonant drum, and the tack source retains its airborne path.

##### Where the body comes from

Four things take energy out of a struck head, and which of them dominates
decides whether the drum has a body at all.

**Radiation** is the largest, and it is the one that separates the modes. How
fast a mode loses energy to the air goes as the square of the air it actually
moves, and a mode with nodal circles moves very little: its annuli alternate in
sign and cancel before the sound has left the head. Integrating *J(0)(λr/a)*
over the disc gives a net volume of *2·J(1)(λ)/λ*, and the same *J(1)(λ)²*
appears in the modal mass, so the two cancel and what is left is a bare
*4/λ²* — the identical weighting the cavity coupling carries, for the identical
reason. Both are net-volume couplings.

The two heads are not placed at the same point. For an axisymmetric pair the
radiation loss therefore uses the coherent two-source power
*b² + r² + 2br·sinc(ωL/c)*, where *b* and *r* are the two head shares and *L*
is the body depth. It reduces exactly to the former *(b+r)²* at zero depth and
is non-negative at every separation. An opposing pair has no monopole moment,
but across a real body it still radiates as a dipole; treating the heads as
coincident had left the shime's lower branch ringing for 1.69 seconds instead
of 0.85 and made the body sound like a second synthetic resonator. This
finite-separation term changes damping only. The front microphone observer
remains the bounded scalar approximation described below until the owned
near/far pressure captures can identify rear-head diffraction and phase.

In the model, the fundamental's more piston-like motion couples strongly to
radiation, while higher modes can ring longer. Earlier numerical comparisons
to unnamed real recordings are not retained as validation: their acquisition
conditions and analysis cannot be checked from this repository.

**The hide's own loss** is viscoelastic, and that is two terms rather than one.
The hysteretic part has a frequency-independent loss angle and damps as *ω*; the
viscous part follows the rate of strain and damps as *ω²*. Only the pair works
here, because this model resolves the low modes individually and treats
the region above its selected crossover as a continuum, and no single power of *ω*
serves both: set for the body, the continuum rings for the best part of a second
as a bed of noise behind the drum; set for the continuum, the body is gone
before it is heard.

**The rim** takes the rest. It is the only term that does not scale with
frequency, which makes it the ceiling on how long anything can ring — a mode
cannot outlast *6.9/edgeLoss* however little else touches it — and measured
decays would be needed to identify it for a particular drum. The factory value
is a model choice, not a retained recording fit. Head Damping scales it from
almost nothing to a great deal, so
the long ō-daiko boom is still there at the bottom of the control. A mode with a
circumferential order pays more of it, because those shapes are pressed against
the boundary rather than spread across the head.

**The mounting** takes what is left, and only at the very bottom. The lowest
modes of a drum do not stay in the head: they move the shell, the hoops and
whatever the drum is stood on. The term is steeply low-pass, because a mode has
to be long enough to move the whole instrument before any of this applies — and
that steepness is measured, not chosen. A gentler skirt reaches far enough up to
cost the body most of what makes it a body.

Where that shelf begins is a property of the drum, not an absolute pitch: a mode
moves the shell when its wavelength is on the order of the instrument's own
size, which is a comparison and therefore scale-invariant. The corner tracks the
radius, as every other frequency in the model already does. Pinned at a fixed
55 Hz it did the opposite of what it describes — a bigger drum slid its whole
modal set down through a shelf that did not move, so the stand ate more of the
instrument the larger the instrument got, and the ō-daiko end of the keyboard
came out both the quietest and the shortest thing on it. Measured across the
four drums: 60 / 119 / 239 / 477 Hz of sounding pitch against 4.2 / 3.0 /
1.5 / 1.7 s of tail, and a full-velocity Don's 20–63 Hz band running 0 / −13.9 /
−28.0 / −36.9 dB against the ō-daiko's. The big drum is bigger in every way that
matters, rather than only in name — and the shime outlasts the okedo despite
being smaller, because its dense carved body absorbs a third of what the okedo's
staves do.

##### The stick

A bachi is a moving mass, not a prescribed force envelope. MIDI velocity gives
it an incoming speed; its compression against the moving head produces a
Hertzian *δ^1.5* force with constrained Hunt–Crossley loss, and contact ends
when stick and hide actually separate. Duration, peak force and rebound all
emerge from that coupled motion.

Below Bachi Hardness 50 % the contact stops being Hertzian. A felt or
wound beater compresses as a power law whose exponent sits well above the
1.5 of two elastic solids — [Chaigne and Doutaut's](https://doi.org/10.1121/1.418117)
xylophone mallets and [Stulov's](https://doi.org/10.1121/1.412299) piano
hammers both fall between 2 and 3 — so the exponent rises linearly from 1.5
at the half-way point to 2.5 at the softest setting. The stiffness for each
exponent is pinned so that a neutral stroke at 72 % velocity keeps the
contact time the Hertz law gave it, and the duration then goes as
*v^((1−α)/(1+α))*: −0.2 for Hertz, −0.43 for felt. A soft beater's contact
therefore shortens far more between a ghost note and a full blow, which is
the brightening a felt bachi actually has, while the hard bachi the family
is played with is untouched. Everything downstream — the pulse's impulse
integral and mean square, the momentum bookkeeping — is evaluated for the
exponent in force rather than for 1.5.

The analytic reference pulse is `sin(πt/τ)^1.5`. Its exact impulse integral
is `sqrt(π)·Γ(5/4)/Γ(7/4) = 1.7480383695280799`; the former 2.3963 is
the integral of `sqrt(sin(x))` and understated the reference peak force. This
reference anchors the statistical residual's level; the moving-head contact
force itself still comes from the coupled solve.

The contact is advanced with a discrete-gradient IMP-2 scheme whose free modal
poles exactly match the drum's existing resonators. The same spatial projection
senses head displacement and spreads force back into the modes, so the coupling
is reciprocal and cannot pull on the hide. Simultaneous bachi contacts are solved
together through one small symmetric compliance system, rather than in an
audibly order-dependent sequence.

The stick's mass scales with the drum being played, because nobody hits a
shime-daiko with an odaiko club. Leaving it fixed made the smallest drums about
twenty-five decibels louder than the largest — a property of the wrong stick
rather than of the instrument.

What it is not is a contact *patch*. The bachi meets the head at a point, and a
real one meets it over a small disc whose radius grows with the force, which
would low-pass the modal bank spatially by *2J₁(k·a_c)/(k·a_c)*. That is left
out for two reasons, both measured. The engine carries only the Hertz product
*K = (4/3)E\*√R*, never the tip radius on its own, so *a_c* cannot be formed
without drawing one — and halving the assumed tip halves the answer. And with
the 12 mm dowel radius of the retired stick model standing in for the tip's
curvature, the factor is worth **0.12 dB** at most on the four family
instruments at full velocity, at the top of the seventy-six-entry bank under
a Rimshot on the shime (0.009 dB for a Don on the ō-daiko, 0.045 dB for a Don
on the shime; it was 0.033 dB when the bank stopped at twenty entries). It
reaches 23 dB only on a 3 cm head struck with the softest beater, whose
highest entries have wavelengths of a few millimetres on that head — which the
controls reach and no taiko is.
Note also which way it runs: a *softer* stick makes a *larger* patch, because
contact stiffness falls far faster than peak force rises, so the mechanism would
dull soft strokes and leave hard ones alone rather than the other way about.
`testTheContactPatchWouldNotBeAudibleOnTheResolvedBank` keeps all of that
recomputed from the engine's own contact solve.

##### The attack pitch glide

A membrane clamped at its rim cannot move without getting longer, and a longer
head is a tighter one. The tension it gains goes as the square of its
displacement — the von Kármán / Berger term — and the pitch as the square root
of the tension, so a struck head starts sharp and settles. Taikor estimates that
strain from the area-mean squared slopes of every resolved batter-head mode,
holds its peaks with a 40 ms release follower, and retunes the shared head from
the result. The follower smooths the reduced-order estimate; it is not a
scripted note envelope, and the head's own modal decay still determines what
keeps it alive.

The resulting bend depends on the strike and the drum. A hard stroke bends
further than a light one because it pushes the head further. A slack head bends
more than a tight one, because the tension a given displacement adds is measured
against the tension already there. An Edge stroke can bend
further than a Don even with less average displacement, because its edge contact
excites shorter, high-gradient shapes and strain depends on slope rather than on
mean height. The depth is computed after the model's one output-level
calibration has been divided out, so that constant cannot reach the drum's
pitch.

It is a first-order expansion, so it is applied through a form that agrees with
it exactly while the displacement is small and saturates where the expansion
stops describing the head. That matters at the edge of the controls rather than
in the middle of them: the fractional tension rise goes as the fourth inverse
power of the radius, so the smallest head at no tension reached fifteen
semitones of bend before it was bounded, and the factory drum reaches a tenth of
a tension at full velocity.

The shared batter-head tension increment now moves each resolved mode according
to the fraction of its stiffness contributed by batter-head tension. Cavity,
resonant-head and bending stiffness are held fixed during this attack strain;
an air-dominated branch consequently glides less than a head-dominated one.
This uses a first-order projection onto the existing modes, checked against an
independent perturbed two-head eigensolve. It follows the distinction between
tensile and bending terms in [Avanzini, Bank and Borin's nonlinear percussion
model](https://home.mit.bme.hu/~bank/publist/jasa12.pdf).

This remains a quasi-static, frozen-mode approximation. It changes the poles
without transferring energy between modes or including the unresolved
continuum in the mechanical strain. Strong nonlinear intermodal coupling and
absolute acoustic calibration require measured head displacement and pitch
glide; improved equations alone do not establish perceptual equivalence to a
recorded taiko.

One shared tension rise is the Berger approximation: the strain is averaged
over the head and every mode feels the same stretch. The head's real strain
is local, and a mode is stiffened most by its own displacement. Each resolved
entry therefore also carries its own ripple — the squared slope of its own
motion, smoothed with the follower's constant, applied as an implicit
per-sample stiffness factor on that entry's resonators so it cannot run away.
A uniform ripple was tried first and rejected by measurement rather than by
ear: on the nagadō-daiko the (0,1) and (2,1) modes sit almost exactly at
twice the (1,1) mode, and a common tension ripple at the (1,1) frequency
pumped them parametrically until the heard pitch of the pad climbed from
119 to 171 Hz. A mode's own ripple has no such partner, and the pad still
measures where the engine says it sounds.

The glide carries the head's continuum with it as well as its resolved modes,
because the continuum is the same head. That is the whole of what Tension Mod
does above a kilohertz, and it is worth about a decibel across the control —
a bend, not a brightness control. Rewriting a running resonator's coefficients
under its own state is not exactly energy-conserving, and it is easy to assume
the difference is heard as spray; measured with the continuum silenced, over
30–80 ms with a settled high-pass and no analysis window involved, what the
rewrite leaves above 1.2 kHz sits 102 dB under the stroke that made it and does
not move with Tension Mod at all.

##### The tack line

A nagado-daiko is *byō-uchi*: the head is not roped on, it is nailed to the
shell with a ring of iron tacks. Each of them holds down the head's tension
over its share of the circumference — a few hundred newtons on the factory drum
— and a stroke that catches the hoop has to beat that before it lifts the head
at all. Past it the tacks chatter against the wood, which is the metallic edge a
firm rim shot has and a light one has no trace of at all.

It is a threshold rather than a level, so it does not fade in: below the preload
there is nothing. Raising Head Tension raises the preload because each tack
carries one fixed spacing of membrane tension; changing diameter at the same
tension does not. Stick Noise owns the level, since this is contact noise. And
it is the one part of the instrument that does not scale with the drum: a byō is
a nail, and the same nails go into a nagadō-daiko and an ō-daiko, so the rattle
keeps its own band across those tacked drums. The rope-laced Okedo and
cord-laced Shime have no tack source at all. Only Rimshot beats the preload at
ordinary velocities; an Edge stroke reaches the hoop with a third of the force and leaves
the tacks alone.

The chatter is confined to its declared 2.6–9 kHz metal-on-wood band by two
serial high-pass and seven low-pass stages. Its exact filter-state variance is
normalised once at the host rate, while only the released filter's energy inside
that physical band is retained. The former one-pole pair put about 55% of its
power above 9 kHz at 48 kHz and moved its centroid with the host clock; that
cheap Nyquist shelf is no longer part of the hit.

##### One shared ringing state per drum

A stroke lands on whatever the head is already doing. Each bachi senses the
surface velocity made by all resolved membrane modes at its own strike point,
and its force is returned through exactly the same shapes and modal masses. A
centre stroke therefore couples strongly to the boom; an edge stroke largely
leaves that boom alone and works on the circumferential modes instead. A second
stroke can shorten, reinforce or delay its contact according to the phase of the
already-moving hide—none of those interactions is scripted.

The nonlinear solve is passive in the resolved head-plus-stick coordinates: it
never creates an adhesive force, and without player input their discrete total
energy cannot rise. Modal displacement stays continuous through every contact.
The live poles used by the solve follow attack glide, automation and pitch-wheel
retuning, so the contact never exchanges energy with a different recurrence from
the one that is actually rendered.

Structural automation can rotate the two normal modes made by an axisymmetric
head pair. A rebuild therefore passes through physical batter- and rear-head
displacement and velocity instead of copying the old upper/lower labels. Pending
force increments make the same transform, and a bachi still in contact receives
fresh reciprocal projections in the new basis.

The Muted stroke goes further. Its free hand remains for 180 ms as a 55 mm-radius damping
patch on motion already ringing on that drum. The same five-point radial area
projection reaches every membrane mode; patch area and modal mass make the same
palm bite harder on a small head than on a five-shaku one, while the angular
orientation of each degenerate pair remains phase averaged.
Control updates set a continuous extra pole loss without stepping either modal
displacement or physical velocity. The continuum receives the corresponding
phase-averaged RMS loss, while every other drum remains bit-identical.

Each playable drum now has one canonical ringing state. MIDI notes allocate
contact transients, project their forces into stable physical-mode IDs, and are
summed before the bank advances once; they do not own another rendered head.
The residual continuum is likewise one persistent field per drum, with new
contacts adding energy in quadrature. Contact slots still reuse the old `Voice`
storage type and therefore carry some unused arrays; that is storage debt, not a
second physical drum.

##### Two microphones, and where the stereo comes from

Taikor's output is a **close stereo pair**, and its image is a consequence of the
model rather than a widener bolted onto it.

A mode whose pattern on the head is finer than the sound it makes cannot
radiate: its field is evanescent and dies as *e^(−√(k(s)² − k²)·d)* above the
surface. Right on the head the two microphones therefore read the *shape* of the
membrane under them, and because every mode with a circumferential order reaches
two different points with a different sign and amplitude, the pair genuinely
decorrelates. A hand's width back, only what the drum radiates survives, and the
image closes towards mono. That narrowing is a mechanism rather than a width
control, and it is the whole of what Mic Distance does to the image.

The two heads are two monopoles a body apart, and the pair stands in front
of one of them. On every cavity-coupled axisymmetric branch the batter head's
propagating share arrives across the capsule distance as it always did, while
the far head's has to come down the body and round the rim: its path is the
body depth plus the diagonal from the rim to the capsule, so it reaches the
pair later and more spread, through the same one-pole spreading law. The
delay puts the two shares out of phase by the mode's frequency times the
extra path, and the mode's quadrature residue carries that phase exactly, in
the convention the complex observer established. What that changes is not a
level but a comb across the breathing branches that moves with Body Depth and
Mic Distance — which is what a real far head does to a close pair. The
continuum's anchor is still measured with both heads in phase at the factory
position, so the statistical layer keeps its level.

It is also a physically bounded softening above the resolved bank. The dense
high modes no longer inherit one loud partial's flat distance gain. Each of the
continuum's five octaves is levelled at the 15.95 cm factory position, then
given a coherent patch radius of *c(head)/(2f)*. A long-wavelength patch keeps
some near-field reach; a short one tends towards the pressure law of a compact
source, so backing the pair from 3 to 40 cm attenuates progressively more of the
top instead of making the distant drum brighter. On the factory ō-daiko the
five target-level drops rise from about 9.1 to 22.1 dB. The factory voicing is
unchanged because that position is the exact unity point of the new law.

The tack source has an additional airborne path with separate propagation
distances and arrival times at the two microphones. The former standalone
stick click used that path too; removing it leaves the head's stereo image
to its modal observations and statistical correlation model.

Fully opened, the two sit about fifty degrees of arc apart, which is what a
close pair over one head actually is. It is worth being strict about that: at a
hundred and twenty-six degrees the capsules sit either side of the nodal
diameter of every mode of order one, and the edge strokes — the ones that drive
those modes hardest — come out of phase.

Up to and including the default 50 % width — everything the microphones actually
captured — no stroke ever inverts, anywhere in the microphone range: the worst
case over the four strokes on all four drums, both microphone controls swept in
tenths, is a correlation of +0.19 — an Edge stroke on the shime with the pair right down
on the head and fully opened — which is a decorrelated pair rather than an
out-of-phase one. The regression suite sweeps that space on the ō-daiko, at
every stroke and at both ends of both controls, and it reaches +0.31 there. Past 50 % the
width control exaggerates the side signal beyond the measurement, and with the
pair close in and fully opened that can push strokes out of phase — the same
thing that happens when a real wide spaced pair is pushed through a widener, and
worth a phase check if the mix has to fold down.

##### Drum Layout

**Drum Layout** is a two-position physical choice:

- **1 Drum** uses the ō-daiko design described by the controls on every keyboard
  row and reaches the four pitches by retuning its head.
- **4 Drums** — the default — resolves each row as its own ō-daiko, nagadō-daiko,
  okedo-daiko or shime-daiko, with an independent diameter, body, hide and shell.

Both layouts land on the same four sounding pitches — 59.75 / 119.49 / 238.99 /
477.98 Hz — but they do not sound alike. Air load, cavity stiffness, radiation
and modal density depend on physical size and cannot be reproduced by pitch
shifting one body. The model therefore solves the required tension/size rather
than transposing a recording.

“1 Drum” describes one physical *design*, not monophony. Each row still owns an
independent ringing state, so a chord can hold four differently tuned copies and
a new row never steals another row's tail. “4 Drums” instead gives those four
states the four family geometries in the table.

This used to be a continuous morph called Octave Body. Intermediate geometry
was not physically useful: one necessary modal-identity handover changed head
size and decay abruptly even though the target pitch stayed fixed. Old projects
keep the same parameter ID and restore safely, with values below the midpoint
mapping to 1 Drum and values at or above it mapping to 4 Drums.

The mode each row is tuned by is latched to the one that instrument is heard at,
rather than re-selected from whichever peak happens to win under automation.
At the factory voicing, rendered Don and Muted strokes across the rows are within
7 cents of exact octaves. Far from the calibrated family a different partial can
become dominant; the tuning itself remains continuous, but a listener may then
name another partial as the drum's pitch.

### Known gaps

**Shared-stage acoustics are not modeled.** Repeated hits share the ringing
state of their own drum, and that drum's heads couple through its enclosed
air. Drum rows and ensemble members do not transfer energy to one another. Each uses its own
local close-pair observation; their outputs are summed without a shared room,
stage coordinates, floor coupling, microphone spill or sympathetic excitation
of unstruck drums. Performer changes gestures, not acoustic connections between
instances.

A real stage has shared acoustic paths. Kodo's sound engineer describes room
microphones, reverberation and spill from a nearby ōdaiko into another
instrument's microphone in [Capturing Taiko](https://www.kodo.or.jp/archives/kodobeat/kb79.pdf).
Sympathetic vibration is a separate mechanical question from hearing spill or
room reflections; its strength needs measurement. A future stage model needs
source/microphone positions and measured room paths, then passive inter-drum
coupling validated with unstruck-head captures. A common reverb alone would not
establish that the drums mechanically excite one another. The output limiter
prevents sample overload; it does not remove acoustic spill.

The room, the player's body, the stand, and the far head's own radiation into
the space behind the drum. The enclosed air carries the stiffness of a finite
column but not its mass or its own resonances: above the body's first axial
resonance — 134.5 Hz on the ō-daiko, 817.7 on the shime, and between 87.9 and 285.8 Hz
across Body Depth on the ō-daiko — the column is treated as absent rather than
as the mass it becomes.
Nothing anywhere associates a loss with the enclosed air either. Body Depth still
has authority over the tail — it moves the breathing branch's T60 on the factory
ō-daiko from 0.893 s at the shallowest body to 1.016 s in the middle and 0.941 s
at the deepest, and it is not a monotone control, because radiation falls as the
branch comes down in frequency while the mounting loss rises towards its corner
— but every bit of that comes from where the branch lands, not from the air.

The missing term is an omission with a number behind it rather than an
oversight. Wall thermal exchange gives an enclosed volume a loss factor
*(γ−1)·δ_t·S/2V*, with *δ_t = √(2α/ω)* the thermal boundary layer; on the
resolved drums that runs 2.8e-4 to 6.6e-4, and it only reaches the mode through
the share of that mode's stiffness the cavity actually holds. Taken over four
octaves, eleven Body Depths and three Air Couplings, the largest decay it would
add anywhere is **1.49 %** of the branch's own — the shallowest small drum at
full coupling, T60 0.3078 s against 0.3032 s. On the factory ō-daiko's breathing
branch it is 0.57 %, and on the branch below it, which barely touches the
cavity, 0.0001 %. The README used to compare that figure against radiation;
on the lower branch radiation is only 0.7 % of the loss and the mounting is
92.6 %, so the comparison was against the wrong term as well as being the wrong
size. `testTheEnclosedAirIsLosslessOnlyWhereThatIsInaudible` recomputes all of
this from the resolved drum rather than remembering it.

Each axisymmetric pair gets its own positive-stiffness column solve. On the
factory ō-daiko their factors are 0.823, 0.610, 0 and 0; relative to reusing the
lowest pair's factor, the upper branches of the second, third and fourth pairs
move down by 12.2, 4.8 and 1.3 cents while the pair the keyboard is tuned by is
unchanged. This is still a scalar one-dimensional correction. Higher radial
modes excite a spatial pressure field, and at the quarter-wave the current
model continuously floors the stiffness to zero rather than introducing the
cavity mass and poles required beyond it.

An experimental replacement is implemented in
[`Source/DSP/CoupledCavity.cpp`](../Source/DSP/CoupledCavity.cpp). It uses eight
physical head coordinates and two shared axial air coordinates by default,
with up to four available for convergence checks. Integrating the air's kinetic
and compression energies gives symmetric mass and stiffness matrices; the
generalized eigensolve retains both cavity inertia and acoustic resonances.
This follows the importance of internal-air mass and compression identified
in [Suzuki and Hwang's Japanese-drum study](https://doi.org/10.1250/ast.29.215),
using an explicit cylindrical, plane-wave reduction rather than claiming the
paper's barrel geometry or measured pressure transfer.

`build/TaikorCavityStudy` compares this solver with the current four factory
drums. Increasing from two to four axial coordinates moves the selected
fundamental-head pairs by less than 0.61 cents in that comparison. Tests also
check the energy integrals against spatial quadrature, positive energy,
mass-orthogonal eigenvectors, the zero-coupling limit, rigid-head air poles,
and convergence against the one-dimensional wave equation. The two-air-mode
solve takes about 10 microseconds on the development machine.

**The experimental cavity solver is not in the plugin's audio path.** It
changes each mode into a mixture of several head and air coordinates; live
contact, damping, microphone observation and retuning must carry that complete
state before activation. The study holds the existing exterior-air frequency
approximation fixed and does not establish calibrated taiko pressure or decay.
An independent-pair inertia shortcut was checked and rejected: it differs from
the shared model by up to about 86 cents on the factory pairs.

Neither the bank's extension nor the statistical continuum is a mechanical
load in the bachi solve. Contact sees the forty coordinates of the calibrated
twenty entries; afterwards its force history drives the further hundred and
twelve resonators and the higher stochastic bands, and none of them exerts a
reciprocal force on the stick. Letting the full bank load the contact is the
right physics, and it was measured: it takes the ō-daiko's momentum transfer
from 13 % to 3 % of the stick's, folds the velocity law and removes the attack
glide, because the shipping contact stiffness and every level pinned by ear
were pinned against a head only as compliant as twenty modes make it. The
ideal membrane's real driving-point mobility is *ω/4T*, and the truncated
solve behaves instead like the resistive *1/(8√(Tσ))*. Closing the split
needs the contact re-calibrated against the full bank with driving-point
mobility captures, and then every by-ear constant re-pinned — a listening
exercise rather than a derivation. Causal force-history integration fixes
temporal and spectral consistency in that observation; it does not complete
the missing mechanical
coupling. The former positive-only `∫F²/Z dt` exposure ceiling no longer gates
the upper spectrum, since truncating the later force would also remove its
physical phase cancellation. The retained exposure calculation is diagnostic
only, and its per-unit-length `Z` does not make it a mechanical energy measure.
Output protection remains downstream of the complete instrument.

**The body is silent on three of the four strokes.** Don, Edge and Muted all
carry no shell gain, so only a rim shot rings the wood. On the factory ō-daiko
that matters more than it sounds: all 152 membrane resonators sit between 33
and 447 Hz, so above 447 Hz an ordinary stroke is the statistical bed and
nothing else, while the six shell ring modes that could fill 0.5-1.7 kHz are
switched off. Resolving that region with modes instead is not an option — the
Weyl count puts about 3100 resonators between 447 Hz and 2 kHz on that drum —
so the wood is the only discrete source available there. A head-to-shell path
driven by the membrane's own boundary shear at the rim is implemented behind
a review switch and is off: it selects modes correctly and is clearly audible,
but being one-way it lets an undiminished head pump the shell, which stretches
the fundamental band's T60 from 2.24 s to 3.51 against a real 1.90. Closing it
needs the coupling loss that a one-way path omits. See `Docs/decisions.md`.

The hide is likewise more uniform than a real one. Every degenerate pair is
split by one fixed constant worth 5.1 cents, beating once every 3.4 seconds at
100 Hz, and the fundamental band's envelope ripple measures 1.7 dB where stage
captures hold 4.9. A per-drum thickness variation is implemented behind a
second switch, also off: at 2.5 % it moves the ripple to 2.1 dB, which is the
right direction and not far enough, so the depth still needs a screened range
and a listening test.

The continuum is also carrying more of the drum than it should. On the factory
ō-daiko its lowest bands are louder than the head's own modes in the pitch
region: the strongest partial below 500 Hz measures 119 Hz with the bands
present and 59.5 Hz without them, while the model's own readout puts that
drum's sounding pitch at 59.7 Hz and real ō-daiko-class captures sit at 61 Hz.
Handing the bands off above the whole resolved bank was rendered as a
candidate and measured: it moves spectral flatness, tail balance and partial
count onto the real distribution, and it was rejected by ear because what it
leaves behind is a sparse set of long-ringing partials over empty space, which
reads as a struck string rather than a drum. Both facts point the same way.
The bands are standing in for resolved modes that are too few and too quiet,
so the level of the extended bank and the losses that thin it are what need
the captures, not the handoff frequency. See `Docs/decisions.md`.

The continuum's absolute level and spatial weighting still need controlled
driving-point mobility and pressure captures. Its first bands on the large
drums can begin before statistical overlap, and its finite frequency nodes can
miss fine interference across a band. Closing these gaps needs a passive
dynamic residual with sufficient frequency resolution, rather than treating
the omitted modes as an uncalibrated resistance.

Resolved modes still reach the capsules through real scalar observations. They
carry mode-dependent near-field cancellation, and every non-axisymmetric mode
keeps its nodal `cos(m theta)` or `sin(m theta)` sign through both the local and
propagating terms rather than acquiring a false mono floor. They still lack a
complex propagation phase and a complete finite-disk radiation pattern. The
continuum now has a frequency-dependent distance law, but its coherent-patch
radius is a bounded physical approximation rather than a solved spatial field;
it does not yet include capsule directivity or the room. The far head now
reaches the pair by its own longer path on the axisymmetric branches, but only
as a delayed monopole: it has no directivity of its own and no path around
the outside of the shell. A complex Rayleigh observer exists as a tested
prototype, but stays out of the release signal path until controlled force,
head-motion and near/far
pressure captures can calibrate it without silently changing the keyboard's
tuning anchor.

Only the struck body rings. The contact force acts equally and oppositely on
both bodies, but on all four strokes the stick that struck the drum is a force
and not a sounding object: nothing here rings the bachi. The engine used to
carry a free-free bar model for one, reached by the stick-against-stick stroke,
and that stroke is no longer on the grid — so the model went with it rather than
sitting unreachable. Building it back for the striker is cheap in modes and
expensive in level, because what a bachi is worth against the drum it is hitting
was only ever pinned by how the stick-against-stick stroke sounded.

Several coefficients remain voicing choices rather than identified physical
parameters: radiation damping strength, shell and tack observation levels,
the statistical upper spectrum's level relative to the resolved head, and the
attack pitch glide's shape factor. Their values need controlled captures of
the mounted drum and microphone arrangement. The glide factor stands in for
the difference between the
modal states the engine has and the mean square slope the tension rise depends
on. Everything about how those terms *vary* with size, material, position and
stroke is computed.
