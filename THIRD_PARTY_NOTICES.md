# Third-party notices

Taikor's original source code is covered by the project's `LICENSE` file.
The build downloads or consumes the separate dependencies below; those
dependencies are not relicensed by Taikor.

## JUCE 8.0.14

JUCE is Copyright (c) Raw Material Software Limited.

Repository: <https://github.com/juce-framework/JUCE/tree/8.0.14>

A copy of the JUCE licence text is vendored at
[`ThirdParty/JUCE-LICENSE.md`](ThirdParty/JUCE-LICENSE.md) so that it can be
included in every distributed bundle and installer package alongside Taikor's
own licence, as the packaging script does.

License: JUCE framework modules are dual-licensed under the GNU Affero General
Public License version 3 (AGPLv3) and the commercial JUCE licence. Building or
distributing Taikor with JUCE therefore requires either compliance with the
AGPLv3 for the complete combined work or an appropriate commercial JUCE
licence. Review the JUCE 8 licence terms before distribution:
<https://github.com/juce-framework/JUCE/blob/8.0.14/LICENSE.md>

JUCE includes or interfaces with additional third-party components. Their
copyright notices and licence terms are listed in JUCE's own `LICENSE.md` and
in the JUCE source tree. The VST3 SDK portions used through JUCE are identified
there as MIT-licensed; Apple's Audio Unit frameworks are supplied by the macOS
SDK and remain subject to Apple's terms.

## CLAP support

The CLAP plug-in uses the MIT-licensed `clap-juce-extensions` wrapper, pinned
to commit `c1a5ad025f95d01e03267857fa8276ebeed16500`:
<https://github.com/free-audio/clap-juce-extensions/tree/c1a5ad025f95d01e03267857fa8276ebeed16500>.
It is Copyright 2019-2020, Paul Walker. The complete upstream licence is
included in [`ThirdParty/CLAP-JUCE-EXTENSIONS-LICENSE.md`](ThirdParty/CLAP-JUCE-EXTENSIONS-LICENSE.md).

The wrapper brings in two dependencies at the revisions recorded by its
submodules:

- CLAP 1.2.7, commit `29ffcc273be7c7c651f6c9953b99e69700e2387a`:
  <https://github.com/free-audio/clap/tree/29ffcc273be7c7c651f6c9953b99e69700e2387a>.
  Its MIT licence is included in [`ThirdParty/CLAP-LICENSE.md`](ThirdParty/CLAP-LICENSE.md).
- CLAP helpers, commit `a61bcdf0ecc2c8db1e80bfe8bf9cb7e8d9fd2bbc`:
  <https://github.com/free-audio/clap-helpers/tree/a61bcdf0ecc2c8db1e80bfe8bf9cb7e8d9fd2bbc>.
  Its MIT licence is included in [`ThirdParty/CLAP-HELPERS-LICENSE.md`](ThirdParty/CLAP-HELPERS-LICENSE.md).

Both dependencies are Copyright (c) 2021 Alexandre BIQUE. These licence files
are included in the distribution archives and macOS installer packages.

## Embedded room impulse responses

The room reverb uses the following third-party impulse responses. The WAV
contents are unaltered copies of the supplied originals; only their filenames
have changed. These assets retain their authors' copyrights and are not
relicensed under Taikor's source-code licence.

| Embedded file | Original file | Source |
| --- | --- | --- |
| `Assets/ImpulseResponses/hall.wav` | `SteinmanHall.wav` | EchoThief Impulse Response Library, `Venues` |
| `Assets/ImpulseResponses/theater.wav` | `MillsGreekTheater.wav` | EchoThief Impulse Response Library, `Venues` |
| `Assets/ImpulseResponses/opera.wav` | `Scala Milan Opera Hall.wav` | Voxengo IM Reverbs Pack |

### EchoThief: Hall and Theater

EchoThief Impulse Response Library is copyright Dr. Chris Warren.
Source: <https://www.echothief.com/downloads/>.

The two embedded WAV files match the corresponding originals in the official
[EchoThief library archive](https://www.echothief.com/wp-content/uploads/2024/07/EchoThiefImpulseResponseLibrary.zip).
The archive's complete, unaltered licence is included at
[`ThirdParty/ECHOTHIEF-LICENSE.pdf`](ThirdParty/ECHOTHIEF-LICENSE.pdf).

That licence permits derivative artistic work, including using convolution to
create reverberation. Other uses are subject to contacting the author at
`cwarren@sdsu.edu`; the published licence does not grant permission to
redistribute the original impulse responses inside another product. Separate
permission is therefore needed before distributing these bundled originals.
The licence also excludes AI training.

### Voxengo: Opera

The Scala Milan Opera Hall response was created with Voxengo Impulse Modeler.
All copyrights and intellectual property rights in this impulse response
remain exclusively owned by Aleksey Vaneev.
Source and published terms: <https://www.voxengo.com/impulses/>.

The complete licence supplied alongside the original WAV is preserved at
[`ThirdParty/VOXENGO-IMPULSES-LICENSE.txt`](ThirdParty/VOXENGO-IMPULSES-LICENSE.txt).
It grants royalty-free use, including commercial use, but redistribution has
additional conditions: copies must remain complete and unaltered, include the
copyright notice and all conditions, carry no charge, and generate no direct
or indirect distribution profit. Both distributor and recipient must
acknowledge Aleksey Vaneev's continuing exclusive ownership. Distribution of
this asset is expressly subject to that acknowledgment and the complete
licence; include the licence and this notice with every distributed copy.
