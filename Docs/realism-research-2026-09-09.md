**Taikor realism research — 9 September 2026**

Research against `d645d4a6d0942deab741fbf9b40ec6aba7064bb8`. The four current voices were studied in parallel: ō-daiko, nagadō-daiko, okedō-daiko and tsuke-shime-daiko. This is a proposal and source audit; it does not change the DSP. Expected audible benefits below are hypotheses to test, not improvements established by listening.

The resulting changes and their limits are recorded in the
[implementation note](physical-realism-2026-09-09.md). Descriptions of the
existing engine below refer to the audited revision above.

The strongest next steps are to give the sticks and skins more independent physical properties, complete the shared head/air model, and represent the mechanical boundary between head and shell. Taikor already has nonlinear contact, repeated strikes on a shared ringing drum, polar strike position, extended membrane modes, pitch modulation, and a rear-head radiation path. Those existing mechanisms should be retained and improved.

| Voice | Physical distinctions that matter | Implication for Taikor |
| --- | --- | --- |
| Ō-daiko | A size category; Taikor chooses a large carved, tacked-head barrel. Kodo describes large carved hardwood bodies and whole cowhide heads. | Identify a particular drum and bachi combination. A large head does not uniquely determine stick weight, shell dimensions or skin properties. [Kodo](https://www.kodo.or.jp/en/about_en/faq) |
| Nagadō-daiko | Carved body with heads nailed in place; tension is normally set during manufacture/reheading. | Prioritize actual modal spacing and rim support, including stable differences around the head. [TaikoMasa](https://www.taikomasa.co.jp/english/products/nagadoudaiko.html) |
| Okedō-daiko | Stave/rope construction varies. Miyamoto's carried model uses cedar and usually horsehide, with cowhide optional; seam placement changes the active striking area. | Give natural-hide material, active span, light bachi and mounting their own identities. Carried and standing okedō should be distinct reference configurations. [Miyamoto](https://www.miyamoto-unosuke.co.jp/products/13030402014004) |
| Tsuke-shime-daiko | Cowhide stitched over metal rings, tensioned by rope or bolts. Available grades vary skin thickness and body height; high pitch does not imply a thin synthetic-like head. | Separate hide thickness from material and tension. Choose a named grade before calibrating it. [Miyamoto](https://www.miyamoto-unosuke.co.jp/collections/shimedaiko) |

These names do not exhaust Japanese drum traditions. The current shime voice specifically represents tsuke-shime; a theatrical shime with a central skin patch would need its own mass distribution. [NSW Arts Unit's instrument guide](https://artsunit.nsw.edu.au/de_modules/aie/the-beauty-of-8/taiko.html) Likewise, a shallow hira/ōhira is a different geometry from the current large barrel. Catalog leather diameter is not necessarily the freely vibrating head diameter, and shime rim/hem thickness is not the vibrating skin thickness.

**1. Give each family a real bachi profile and model its finite contact area.**

The present `drumContactTerms()` scales one 0.19 kg stick directly with drum radius. Before collision effective-mass calculation, this assigns approximately 518 g to the nominal 1.50 m ō-daiko and 138 g to the nominal 40 cm okedō. It has no independent tip footprint or stick geometry. The lower half of Bachi Hardness becomes felt contact, which cannot separately represent a light, softer piece of bare wood.

Miyamoto distinguishes bachi by use, wood, length and shape. Performer Tsuyoshi Maeda lists a magnolia katsugi bachi at 365 × 19 mm and approximately 50 g, and describes adjusting the tip's chamfer for touch. These provide concrete reference configurations. Static stick mass alone does not determine effective impact mass: grip and motion also matter. [Miyamoto bachi guide](https://www.miyamoto-unosuke.co.jp/pages/wadaiko-bachi-erabikata), [Maeda's specification](https://www.tsuyoshimaeda.com/blank-1)

Start with separate effective mass, bare-wood compliance, tip size and contact-loss parameters for large hinoki, nagadō wood, light okedō and shime sticks. Keep felt as a separate option. Project force over a small contact patch and use the same projection when sensing head motion; this preserves reciprocal contact. Distributed mallet contact is established in general percussion models, though its parameters still need taiko measurements. [Laird's percussion-model thesis, §6.2](https://slab.org/software/laird/joel_laird_thesis.pdf)

Then revisit the higher modes' mechanical load. Today only the original 20 modal entries react back on the stick; the extended bank and statistical upper spectrum are driven afterward. The previous full-bank experiment reduced ō-daiko momentum transfer and damaged the velocity response, so simply enabling that flag repeats a rejected experiment. Finite-area contact and a passive high-frequency residual should be calibrated together before changing the spectral balance.

Expected benefit: more distinct thump/crack, clearer soft strokes and plausible rebound across the four drums. Profile separation is a relatively small implementation; completing the high-frequency contact mechanics is larger. This is my first prototype choice.

**2. Give each head its own measured skin, tension and spatial imperfections.**

`headMaterial` currently moves areal density, volumetric density, bending stiffness and loss together along a synthetic-film-to-cowhide interpolation. The shime profile selects 0.18 on that axis to obtain a thin, high voice. This makes several independent physical choices inseparable. A real high-tension, thick-hide shime is a useful alternative target: Miyamoto's four-chō model explicitly uses thicker cowhide than its three-chō version. Its listed roughly 26 mm hem must not be used as membrane thickness. [Four-chō shime](https://www.miyamoto-unosuke.co.jp/collections/shimedaiko/products/13020101040601)

Separate areal mass, effective bending stiffness, damping, active radius and tension in the internal family description. Allow a thin natural hide without blending its material toward polyester. Calibrate each head independently; the current rear density is always 92% of the batter's density, with rear tension limited to 85–115% of the batter's.

Use a smooth, fixed density/tension map to generate slightly altered frequencies and mode shapes. At present every drum uses the same small modal-pair split by entry index and ideal Bessel shapes. An offline spectral solve could give each instrument persistent asymmetric nodal patterns; force and microphone projections must both use those patterns. This technique is supported by studies of nonuniform tabla/tom heads, so applying it to taiko is an inference requiring validation. [Samejima and Fukuda, 2016](https://www.jstage.jst.go.jp/article/ast/37/6/37_E1556/_article/-char/en)

The published nagadō evidence argues against blindly adding stiffness. Suzuki et al. found that low-mode spacing on their measured single-headed specimen did not show the theoretically predicted bending effect. I reran `python3 Tools/BenchmarkPublishedModes.py`: after removing one common tuning offset, the ideal-membrane hypothesis gives 44.27 cents RMS spacing error, while Taikor's single-head equation at the paper's inferred tensions gives 85.70–89.09 cents. These compare equation hypotheses for one open, floor-influenced specimen, not the full plugin's four presets. [Suzuki et al., 2009](https://www.jstage.jst.go.jp/article/ast/30/5/30_5_348/_pdf)

Expected benefit: a more convincing pitched body, stable position-dependent differences, and natural beating without random retuning each hit. Fit modal decay and upper-bank strength alongside frequencies; the existing decision log records that reducing the continuum alone exposed an unconvincing sparse ring. Preserve the current MIDI layout while evaluating measured physical profiles; the present exact-octave family tuning is a musical design constraint, not evidence of real drum dimensions.

**3. Complete the shared two-head/cavity model, then support strikes on either head.**

Production still couples each axisymmetric radial pair independently and floors its air stiffness beyond the positive branch. `Source/DSP/CoupledCavity.cpp` is a tested offline prototype, not an active audio-path feature. It supplies shared air inertia and axial resonances that the released model omits. Suzuki–Hwang provides taiko-specific theoretical support for both air inertia and compression and for the effects of unequal heads, using an idealized 48 cm rigid-body drum. [Suzuki and Hwang, 2008](https://www.jstage.jst.go.jp/article/ast/29/3/29_3_215/_pdf/-char/en)

I reran `build-ci-macos/TaikorCavityStudy`. These are selected model eigenfrequency pairs, not measured sounding pitches:

| Factory geometry | Current independent pair, Hz | Shared cavity, two axial coordinates, Hz |
| --- | --- | --- |
| Ō-daiko | 32.650 / 61.372 | 26.961 / 58.835 |
| Nagadō | 68.048 / 94.506 | 56.578 / 93.625 |
| Okedō | 238.989 / 256.636 | 197.370 / 255.899 |
| Shime | 477.979 / 515.119 | 428.289 / 513.781 |

Increasing two axial coordinates to four changes these selected pairs by at most 0.61 cents. That demonstrates convergence for this reduced model, not acoustic accuracy. The prototype still omits transverse pressure modes, barrel flare, damping and complete radiation.

Integration must transform excitation, contact sensing, damping, microphone observation and automation into the full mixed head/air basis. Replacing frequencies while retaining old participation factors would be inconsistent. After that, a side selector can strike either physical head while preserving one shared ringing state. Asano's specialized Kanade provides a documented reference with different head thicknesses and independently adjustable tension, pitch and timbre; that independent tuning hardware is not universal to all okedō. Merely triggering two independent pitched voices would omit their interaction. [Asano Kanade](https://www.asano.jp/files/kanade/)

Expected benefit: more plausible development of the low/mid body resonance, head-to-head energy transfer and front/back contrast. Strong physical rationale, substantial implementation scope, and audible gain still to establish.

**4. Make shell construction affect the head through a physical boundary.**

The current wooden bank uses six isotropic cylinder-ring modes. Direct hoop strikes excite it; ordinary head strikes do not transfer force into audible shell modes. That omission is deliberate: an earlier duplicated force path made the light okedō shell dominate. Restoring the duplicated force would recreate the problem.

A Japanese barrel study matched seven measured low modes with a model that accounts for different wood stiffness along and across the grain; cross-grain stiffness was crucial. This supports a small offline-derived carved-shell model. Extending the method to stave bodies and tensioned hoops requires separate calibration. [Hwang and Suzuki, 2016](https://www.jstage.jst.go.jp/article/ast/37/3/37_E1526/_article)

Connect shell/hoop coordinates to the head through a passive mechanical boundary. Fit compliance and losses at the rim before increasing shell radiation. Ono et al.'s wadaiko prototypes showed shell-material effects through membrane tension and boundary conditions, while the shell's own sound did not explain the overall difference. That is evidence for prioritizing the boundary, not making every drum audibly ring like wood. [Ono et al., 2009](https://www.jstage.jst.go.jp/article/ast/30/6/30_6_410/_article)

Expected benefit: more characteristic edge/hoop response and decay differences between heavy carved barrels, light staves and tensioned hoops. A few reduced shell modes are feasible; the strength of head-to-shell coupling needs measurement.

**How to decide whether a candidate is better**

Use an identified mounted drum from each family, with measured active diameter, depth, both heads and bachi. Capture quiet, medium and hard strokes at several radii and angles, with repeated hits and rolls; record both heads for the okedō comparison. Log temperature/humidity, mounting, microphone geometry and processing. Maker descriptions establish construction, but do not provide the controlled modal/force dataset needed to tune these mechanisms.

For each isolated candidate, compare contact duration/rebound, modal frequencies and decay, attack-to-tail spectral balance, velocity response and near/far pressure. Keep calibration strokes separate from held-out validation strokes. Use consistent playback gain across a velocity series so level normalization does not conceal a broken dynamic law, followed by blind listening at matched reference loudness. Audit energy/passivity, live retuning and sample-rate stability before putting a candidate into a release.

The first practical experiment should separate bare-wood bachi properties and skin properties, with one reversible candidate per mechanism. The shared cavity follows as the larger structural change. Full-bank contact and continuum handoff stay experimental until their coupled response is convincing. Room acoustics, environmental drift and extra shell loudness have lower priority for this pass.
