# Physical realism implementation — September 9, 2026

This revision implements the four mechanisms proposed in the
[construction and acoustics research](realism-research-2026-09-09.md).
It preserves the vintage Japanese editor, sixteen-note articulation map and
existing serialized host parameter IDs. It is a reduced physical model with
explicit construction priors; it is not a measured reproduction of a named
maker's drum.

## Bachi and skin

`BachiModel.h` separates stick mass, tip radius, elastic compliance and covering.
The four families have distinct wooden reference sticks. Diameter changes no
longer resize the striking mass. Hardness trims wood/contact compliance; a soft
wooden tip remains Hertzian with exponent 1.5. A felt covering is an explicit
internal profile type, not an automatic result of turning Hardness down.

Force and surface sensing use the same finite circular footprint. A circular
mode with wavenumber k receives the disk-average factor 2J1(kr)/(kr), including
its sign. At the rim the patch shrinks to stay inside the membrane. This is a
uniform circular contact approximation, not a measured pressure distribution
or a changing Hertz contact radius during impact.

`PhysicalDrumProfile` gives each head its own areal mass, thickness, effective
elastic modulus, loss angle and viscous loss. Rear tension has an independent
family ratio. All four factory profiles are natural hide, including a thick,
high-tension tsuke-shime design. HEAD trims the family profile around 75%; its
lowest fifth blends toward an explicit synthetic-film endpoint. The shime
head's retail folded hem is not used as membrane thickness.

Fixed modal imperfections split squared frequencies about the same mean and
rotate nodal axes. One basis is used for strike projection, reciprocal sensing,
left/right observation and live state transfer. The axes belong to the drum;
Humanise 0 still preserves the player's authored point. These perturbations
are engineering priors, not estimates of variation among all handmade drums.

## Shared heads and air

The first four radial shapes on each head share one axial air field, represented
by two sine coordinates. Integrating kinetic and compression energies gives
symmetric generalized mass and stiffness matrices. The eigensolve returns
unit-mass modes, each potentially containing several head and air coordinates.

The production engine now carries that full basis through contact force and
sensing, microphone projection, skin/radiation/thermal loss, finite palm
projection, strain and structural state remapping. Automation reconstructs
physical head and air displacement/velocity before projecting into the new
basis. In-band pure-air coordinates remain stored even at zero coupling so a
structural change does not silently discard their state. Their force and
observation projections naturally vanish when decoupled.

The shared model replaces the first four independent radial pairs. The four
higher radial pairs retain the previous bounded column approximation. This
is still a cylinder with plane-wave internal air: transverse pressure modes,
barrel flare, leakage and detailed thermoviscous boundary layers are omitted.
The thermal decay used by the engine is a small, explicit wall-loss prior.

Factory tuning uses explicit head-tension priors at the declared dimensions;
it no longer changes drum size to follow a different observed air resonance.
The Nagado profile has a 1.75 mounting-loss scale, shortening the competing
low breathing mode without moving its frequency. This is a voicing prior,
not an identified support impedance from a measured instrument. The single
drum layout's upper octaves use virtual tensions beyond a validated hide range.

Live Pitch and wheel changes project the change in head stiffness into the
current shared basis. The enclosed air does not transpose independently as
though it were a stretched skin. A subsequent strike at the new tuning solves
the full cavity again and remaps the ringing physical state into that basis.
During wheel movement alone the fixed-basis projection omits changing mode
shapes and off-diagonal stiffness; it avoids a full eigensolve for every bend.

The statistical batter spectrum follows the skin's own spatial dispersion.
Adding an arbitrarily weak air resonance therefore cannot jump its band
cutoff. Its digital filters now peak at the same frequencies used for force
integration; the previous asymmetric filter cascade shifted their peaks down
into the resolved spectrum. Variance normalization and live retuning use the
same corrected coefficients. The old cavity column factor remains available
as a legacy diagnostic;
the active readout separately identifies the shared solver and internal-air
restoring-energy participation.

## Play either head

MIDI CC18 values 0–63 select the front head; 64–127 select the rear. Selection is
captured by each note, including delayed ensemble companions. A later control
change cannot move an existing contact to the other skin. CC121 and reset
restore front strikes. No keyswitch steals a note from the playing grid.

Both skins have independent angular banks in addition to the common radial
cavity coordinates. A rear strike drives and senses the rear skin's physical
coordinates. Mean strain and each entry's own tension ripple use the selected
skin's displacement and elastic coefficient. A muted rear articulation lays
its patch on that skin; CC1 remains a palm on the front head.

The microphones stay in front of the batter. The rear angular observer uses
a reduced effective path around the rim with propagation phase and spreading;
it is not a full diffraction solution for a wooden barrel. The five unresolved
statistical bands currently belong only to the batter. Rear strokes excite the
resolved spectrum and do not feed an artificial copy of the front noise layer.

## Wooden shell and boundary

Six orthotropic ring modes represent carved, stave and tensioned-hoop
construction, with per-mode physical masses, frequencies and losses. A signed
near-rim transformer links the mechanical head modes to the shell through a
reciprocal viscous bearing-edge port. The same pre-impulse relative velocity
acts on all connected coordinates. Its exact exponential step cannot increase
the physical modal kinetic energy; displacement is unchanged by that step.

Only the reciprocally solved head coordinates join the port. Observation-only
upper modes must not form an indirect unaccounted contact load through the
shell. A rimshot also has a direct external shell force. Head-only strokes do
not duplicate that force into the shell: their shell motion comes from the
boundary transfer.

Radiation gain, family voicing and Shell Resonance belong to microphone
observation, separate from physical mass and momentum. The output was adjusted
for the newly exposed lightweight shell masses rather than imposing a fictitious
minimum moving mass. Equations, assumptions and a level audit are in
[the shell implementation note](shell-boundary-implementation-2026-09-09.md).

The boundary is an effective lossy projection, not an equality between axial
head velocity and radial shell velocity. It uses one averaged two-end port and
one orientation per ring order, with no axial shell modes or elastic bearing
spring. Measured directional mobility would be needed to identify these terms.

## Earlier review switches

The upstream `headToShellPath` and `hideInhomogeneity` candidates retain their
original switch bits and remain disabled. The first selects the earlier
one-way shear experiment in place of the reciprocal boundary; the second
selects the stronger geometry-seeded split while retaining consistent angular
projections. Their historical screening notes remain in `Docs/decisions.md`.
They are not the production realism mechanisms. The shared cavity uses its
own bit so a prior review mask cannot accidentally enable it.

## Verification and remaining calibration

Focused tests exercise independent skins, fixed angular bases, finite contact,
generalized mass normalization and reciprocal compliance, head/air state
transfer, selected-side mute persistence, rear-head contact and strain,
ensemble event capture, shell energy loss and bidirectional runtime transfer.
Block-partition replay is checked at 8, 48 and 192 kHz, alongside the existing
contact, output, host-state and automation regressions.

The original twenty spatial entries per struck head plus the shared cavity
remain the reciprocal contact load. Higher resolved entries are driven by the
solved force without feeding back; a complete broadband mechanical termination
still needs contact-mobility calibration. The self-strain approximation omits
cross-entry parametric pumping and updates cavity tension in a fixed basis.
These are deliberate limits, not claims of a full nonlinear finite-element drum.

A controlled dataset is still needed to calibrate acoustic accuracy: identified
drums and both skins, documented mounting and microphones, known bachi,
repeated velocities and positions, and dry recordings. Numerical consistency
and stable rendering do not by themselves establish perceptual realism.
