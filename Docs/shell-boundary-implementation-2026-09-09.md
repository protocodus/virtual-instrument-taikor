# Shell boundary prototype — updated 11 September 2026

This implements the fourth proposal in `realism-research-2026-09-09.md`. The source audit in that report describes the earlier shipped engine. The current candidate has reciprocal shell/head coupling; its construction and microphone gains remain engineering priors requiring measured and listening calibration.

## Mechanical reduction

`Source/DSP/ShellBoundary.cpp` supplies six circumferential ring coordinates, orders 2–7. Wood has distinct axial and transverse Young's moduli and reciprocal Poisson ratios. The transverse flexural rigidity sets ring stiffness; both radial and tangential displacement contribute modal mass. Effective carved-barrel, thin-stave and tensioned-ring properties differ by family. The host's 1 Drum / 4 Drums setting selects the reference or family construction; the internal profile representation supports interpolation. Shell Material is a trim around each family's nominal construction, so the existing family control offsets do not count the same material difference twice.

This is a ring reduction, not the seven-mode fitted barrel model in [Hwang and Suzuki (2016)](https://www.jstage.jst.go.jp/article/ast/37/3/37_E1526/_article). It retains one directional coordinate per order and omits its orthogonal partner, axial shell bending, resolved stave joints and individual hoops. Its family densities, wall thicknesses, elastic constants and damping are unmeasured priors. The shime's effective added ring mass does not resolve a ring's separate motion.

Each ring couples to the matching order of both heads. A narrow near-rim sampling band uses the same fixed rotated membrane basis as excitation and observation. The two heads contribute half-weighted projections; the wall contributes the opposite sign. This is an effective mechanical transformer between axial membrane motion and radial ring motion, not a measured rim mobility or an assertion that those orthogonal motions are equal. The coupling currently adds dissipation and reciprocal transfer; it does not add an elastic rim spring or solve a flexible boundary eigenproblem.

## Passive discrete step

For signed generalized port coefficients `b`, inverse modal masses `M^-1`, damping `c >= 0` and step `h`, define `u = b^T v` and `g = b^T M^-1 b`. After each free resonator step, apply

```
J = -(1 - exp(-c*g*h)) * u / g
v <- v + M^-1 b J
```

Current displacement stays fixed. The kinetic energy change is exactly `-(1-exp(-2*c*g*h))*u*u/(2*g)`, which is nonpositive for every positive step. Opposite ports exchange energy in both directions. A head strike has zero direct shell contact projection; its shell response arises through this boundary. Don Rim retains its direct hoop-force path.

Velocity is recovered from the cached damped-pole quadrature. A velocity impulse updates only the stored previous displacement. Impulse coefficients are cached when configuring modes; the audio step is two linear passes and six scalar accumulators. Missing or filtered-out rings disable their corresponding boundary instead of grounding the head against an absent wall.

This proof covers the boundary step and its composition with exact unforced damped poles in their physical modal energy. It does not by itself prove passivity of the complete nonlinear contact discretization, which uses a separate midpoint energy metric and has its own strict tests.

## Radiation and level calibration

Shell coordinates now use their actual per-mode mass. The previous minimum-mass substitute and acoustic multiplier inside displacement are removed. Family microphone voicing is applied only in observation, together with radius squared, the ring-index taper, Shell Resonance and the existing microphone-distance ratio. This keeps output gain separate from force, energy and head/shell transfer.

The observer calibration is 1100 / 1000 / 150 / 1500 for ō-daiko / nagadō / okedō / shime, interpolated by the internal family mix. These are explicit musical level choices, not measured radiation efficiencies. During the initial 9 September prototype, retaining the old common gain of 4200 produced a +7.6 dBFS pre-limiter peak on a hard okedō rim accent. That historical figure used the earlier geometry/tension mapping; it is not a measurement of the current candidate. Physical shell masses and observer gains were not changed for the audit below.

The 11 September audit uses the current fixed dimensions and explicit family tension profiles, 48 kHz, default 4 Drums parameters except Humanise and Tension Modulation off, one front-head stroke at velocity 1, and one second of raw stereo output. Raw samples are scaled by the default output gain of 0.075 and the engine's −4 dB reference calibration. Values precede DC filtering, Drive, output high-pass filtering and limiting. They are the larger channel's peak, not a loudness match or evidence of listening preference.

The shell-only pass mutes both heads' microphone coefficients, the statistical continuum's output and the tack observation while retaining the same mechanical contact and head–shell transfer. Its result is independent of residual-layer voicing. Full output includes those layers and was rerun after correcting the residual filter's actual peak to match its nominal frequency. The correction left every shell-only measurement unchanged.

| Family | Full Edge peak, dBFS | Shell-only Edge peak, dBFS | Full Don Rim peak, dBFS | Shell-only Don Rim peak, dBFS |
| --- | ---: | ---: | ---: | ---: |
| Ō-daiko | −19.27 | −84.95 | −13.96 | −16.38 |
| Nagadō | −25.52 | −98.75 | −19.40 | −25.86 |
| Okedō | −29.32 | −105.51 | −17.27 | −18.28 |
| Shime | −37.22 | −101.87 | −22.95 | −22.96 |

Velocity 0.45 and 0.9 were also checked; every audited full and shell-only peak remained below 0 dBFS before the limiter. At velocity 1, the normal edge stroke's shell-only peak lies 65–76 dB below the complete stroke with this voicing. The intended effect is boundary-dependent decay and a controlled rim accent, not a loud wooden overlay. [Ono et al. (2009)](https://www.jstage.jst.go.jp/article/ast/30/6/30_6_410/_article) supports prioritizing the membrane's boundary over the shell's own emitted sound, but does not establish these numeric gains or damping values.

## Verification

`ShellBoundaryTests` checks unequal-mass and distributed signed-port passivity, both transfer directions, exact dissipation, time-step composition, damped-pole conversion, construction differences and geometric scaling. `ShellBoundaryIntegrationTests` checks the real builder/render path, absent-ring bypass, shell motion after an ordinary head strike and head motion after an initially moving shell. It also tests unforced energy decay and the boundary step on states reached during active front/rear head and rim contacts at 8, 48 and 192 kHz.

The integration suite also changes Shell Resonance on a ringing drum and verifies that its sound changes while its pole coefficients and physical state remain identical. `AngularRebuildIntegrationTests` preserves both heads' physical displacement, velocity and pending input when a structural change rotates a split angular pair while bringing its missing partner below Nyquist.

Existing distance tests now reflect the thin okedō wall's subcoincidence modes: higher ring orders lose their evanescent near field more quickly with distance. Direct head-to-shell force remains forbidden, while Shell Resonance may now change the subtle boundary-transferred response. Global nonlinear-contact energy tolerances are unchanged.
