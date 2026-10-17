#!/usr/bin/env python3
"""Write the committed MIDI phrases for the README demos.

usage: write_preset_midi.py [output-dir]   (default: Docs/audio/midi)

Writes one phrase per factory preset into the output directory and one per
More demos showcase into its more/ subdirectory.

Writes Docs/audio/midi/<NN>-<slug>.mid, one per entry of
Source/DSP/FactoryPresets.cpp in the same order; the slug is the preset name
lowercased with every other character run turned into one hyphen, which is
how Tools/PresetDemoRenderer.h finds the file. The files are committed; run
this only when a phrase is changed on purpose, and commit that change on its
own. Every phrase is written in kumi-daiko idiom for what its preset is for.

Taikor's keyboard: each drum is one row of four keys, C = Don (open centre),
C# = Edge, D = Muted (palm on the head), D# = Rimshot. C3 (48) is the o-daiko,
C4 (60) the nagado-daiko, C5 (72) the okedo and C6 (84) the shime. CC1 lays a
hand on the head.
"""

import pathlib
import struct
import sys

PPQ = 480
ODAIKO, NAGADO, OKEDO, SHIME = 48, 60, 72, 84
DON, EDGE, MUTED, RIM = 0, 1, 2, 3


class Phrase:
    def __init__(self, bpm: float):
        self.bpm = bpm
        self.events: list[tuple[int, int, bytes]] = []  # (tick, order, message)

    def hit(self, beat: float, drum: int, stroke: int, velocity: int) -> None:
        tick = round(beat * PPQ)
        note = drum + stroke
        self.events.append((tick, 1, bytes([0x90, note, velocity])))
        # Drums ignore note length; a short note keeps the file tidy in a DAW.
        self.events.append((tick + PPQ // 8, 0, bytes([0x80, note, 0])))

    def pattern(self, start: float, step: float, drum: int, strokes: str, velocities: list[int]) -> float:
        """One character per step: D Don, K Edge, T Muted, R Rimshot, . rest."""
        codes = {"D": DON, "K": EDGE, "T": MUTED, "R": RIM}
        for index, symbol in enumerate(strokes):
            if symbol != ".":
                self.hit(start + index * step, drum, codes[symbol], velocities[index % len(velocities)])
        return start + len(strokes) * step

    def hand(self, beat: float, amount: int) -> None:
        self.cc(beat, 1, amount)

    def cc(self, beat: float, controller: int, value: int) -> None:
        # Controllers sort before notes on the same tick, so a stroke uses them.
        self.events.append((round(beat * PPQ), 0, bytes([0xB0, controller, value])))

    def bend(self, beat: float, amount: float) -> None:
        """amount in -1..1; 0 is the centred wheel."""
        value = max(0, min(16383, 8192 + round(amount * 8191)))
        self.events.append((round(beat * PPQ), 0, bytes([0xE0, value & 0x7F, value >> 7])))

    def end(self, beat: float) -> None:
        self.events.append((round(beat * PPQ), 2, b""))

    def smf(self) -> bytes:
        tempo = round(60_000_000 / self.bpm)
        track = bytearray(b"\x00\xff\x51\x03" + tempo.to_bytes(3, "big"))
        now = 0
        for tick, _, message in sorted(self.events):
            track += variable_length(tick - now)
            now = tick
            track += message if message else b"\xff\x2f\x00"
        if not track.endswith(b"\xff\x2f\x00"):
            track += b"\x00\xff\x2f\x00"
        return (b"MThd" + struct.pack(">IHHH", 6, 0, 1, PPQ)
                + b"MTrk" + struct.pack(">I", len(track)) + bytes(track))


def variable_length(value: int) -> bytes:
    out = [value & 0x7F]
    value >>= 7
    while value:
        out.append(0x80 | (value & 0x7F))
        value >>= 7
    return bytes(reversed(out))


# Each phrase ends with rests so the last stroke's tail is part of the demo.

def kumi_daiko_stage() -> Phrase:
    # The whole family: shime ji-uchi underneath, nagado main pattern, o-daiko
    # on the downbeats, okedo answering, and a unison finish.
    p = Phrase(120)
    for bar in range(6):
        b = bar * 4
        p.pattern(b, 0.5, SHIME, "DKDKDKDK", [78, 52])
        p.pattern(b, 1.0, ODAIKO, "D..." if bar % 2 else "D.D.", [112, 90])
        if bar >= 1:
            p.pattern(b, 0.5, NAGADO, "D.DKD.RK", [105, 70, 95, 60])
        if bar >= 3:
            p.pattern(b + 0.25, 0.5, OKEDO, ".K.K.KDK", [64, 58, 70, 86])
    for drum in (ODAIKO, NAGADO, OKEDO, SHIME):
        p.hit(24, drum, DON, 120)
    p.hit(25, NAGADO, RIM, 110)
    p.hit(26, ODAIKO, DON, 127)
    p.end(32)
    return p


def festival_ensemble() -> Phrase:
    # A matsuri pattern played in unison by the group: don-don, kara-ka.
    p = Phrase(132)
    for bar in range(8):
        b = bar * 4
        p.pattern(b, 0.5, NAGADO, "DD.KKKD." if bar % 4 != 3 else "DDKKRKRD", [112, 100, 70, 66])
        p.pattern(b, 2.0, ODAIKO, "DD", [118, 96])
        p.pattern(b, 0.5, SHIME, "DKDKDKDK", [72, 48])
    p.hit(32, ODAIKO, DON, 127)
    p.hit(32, NAGADO, RIM, 124)
    p.end(40)
    return p


def close_studio() -> Phrase:
    # Nagado solo with the detail a close pair shows: ghost muted strokes,
    # edge-of-head accents and a crescendo into rimshots.
    p = Phrase(96)
    for bar in range(4):
        b = bar * 4
        p.pattern(b, 0.25, NAGADO, "D.T.K.TDD.TKR.T.", [100, 30, 62, 34, 92, 88, 32, 70, 110, 28])
    velocities = [40 + 6 * index for index in range(12)]
    p.pattern(16, 0.125, NAGADO, "KTKTKTKTKTKT", velocities)
    p.hit(17.5, NAGADO, RIM, 120)
    p.hit(18.5, NAGADO, DON, 112)
    p.end(24)
    return p


def tight_ji_uchi() -> Phrase:
    # The base rhythm that keeps a group together: don-ko don-ko on the shime,
    # okedo doubling the offbeats, accented every bar.
    p = Phrase(140)
    for bar in range(8):
        b = bar * 4
        p.pattern(b, 0.25, SHIME, "DKDKDKDKDKDKDKDK", [100, 54, 80, 54])
        p.pattern(b, 0.5, OKEDO, ".D.D.D.R" if bar % 2 else ".D.D.D.D", [84])
        p.hit(b, NAGADO, DON, 96)
    p.hit(32, SHIME, RIM, 118)
    p.hit(32, OKEDO, RIM, 110)
    p.end(36)
    return p


def soft_bachi() -> Phrase:
    # A quiet introduction: soft open strokes, a swelling roll and a single
    # full stroke, the way a piece often opens before the group comes in.
    p = Phrase(72)
    p.pattern(0, 1.0, ODAIKO, "D.D.", [36, 30])
    p.pattern(4, 1.0, NAGADO, "D.DD", [40, 32, 46])
    roll = [18 + 3 * index for index in range(24)]
    p.pattern(8, 0.125, NAGADO, "DT" * 12, roll)
    p.hit(11, ODAIKO, DON, 88)
    p.pattern(13, 0.5, NAGADO, "K.K.D...", [34, 30, 44])
    p.end(18)
    return p


def single_odaiko() -> Phrase:
    # One drum retuned four ways: a slow call walking down the rows, with a
    # hand laid on the head to stop a stroke short.
    p = Phrase(66)
    for index, drum in enumerate((SHIME, OKEDO, NAGADO, ODAIKO)):
        p.hit(index * 1.5, drum, DON, 104)
    p.hit(6.5, ODAIKO, DON, 116)
    p.hand(8.5, 110)
    p.hand(9.0, 0)
    p.hit(9.0, NAGADO, RIM, 108)
    p.hit(10.0, ODAIKO, DON, 124)
    p.end(16)
    return p


def synthetic_heads() -> Phrase:
    # A practice-room groove on okedo and nagado, edge strokes in between.
    p = Phrase(116)
    for bar in range(6):
        b = bar * 4
        p.pattern(b, 0.5, OKEDO, "DKKDKKDK", [96, 58, 52])
        p.pattern(b, 1.0, NAGADO, "D.TR" if bar % 2 else "D.D.", [104, 50, 112])
    p.hit(24, NAGADO, DON, 120)
    p.hit(24, OKEDO, RIM, 116)
    p.end(28)
    return p


def distant_pair() -> Phrase:
    # An o-daiko call and the whole group's answer, heard from the front of a
    # small stage with the tails left to ring.
    p = Phrase(84)
    p.pattern(0, 1.0, ODAIKO, "D.D.DD..", [120, 100, 110, 124])
    p.pattern(8, 0.5, NAGADO, "DDKKD.R.", [110, 104, 72, 70, 112, 118])
    p.pattern(8, 0.5, SHIME, "DKDKDKDK", [80, 54])
    for drum in (ODAIKO, NAGADO, OKEDO, SHIME):
        p.hit(12, drum, DON, 124)
    p.end(20)
    return p


# More demos: playing techniques the presets do not feature. Each renders
# through the factory preset named beside it in Tools/RenderPresetDemos.cpp.

def rolls_and_flams() -> Phrase:
    # Flams on the nagado, a press roll swelling into a rimshot, then the
    # shime's single-stroke roll accelerating under an o-daiko call.
    p = Phrase(90)
    for bar in range(2):
        b = bar * 4
        for beat in (0, 1, 2, 3):
            p.hit(b + beat, NAGADO, DON, 70)          # grace note
            p.hit(b + beat + 0.06, NAGADO, DON, 112)  # main stroke
    swell = [30 + 4 * index for index in range(24)]
    p.pattern(8, 1 / 6, NAGADO, "DT" * 12, swell)
    p.hit(12, NAGADO, RIM, 124)
    t, step, hand = 13.0, 0.25, 0
    while step > 0.07:
        p.hit(t, SHIME, (DON, EDGE)[hand], 60 + int((0.25 - step) * 300))
        t, step, hand = t + step, step * 0.93, 1 - hand
    p.hit(t, ODAIKO, DON, 124)
    p.hit(t, SHIME, RIM, 120)
    return p


def dynamics() -> Phrase:
    # One open stroke from a ghost note to a full-arm hit and back down, on
    # the nagado and then the o-daiko.
    p = Phrase(100)
    rise = [8, 16, 26, 38, 50, 62, 74, 86, 98, 110, 120, 127]
    p.pattern(0, 0.75, NAGADO, "D" * 12, rise)
    p.pattern(9, 0.75, ODAIKO, "D" * 12, list(reversed(rise)))
    return p


def hand_damping() -> Phrase:
    # The o-daiko left to ring, then a hand laid on it (CC1) a little, more,
    # and fully, so each stroke is cut shorter; then accents choked on the beat.
    p = Phrase(80)
    for index, amount in enumerate((0, 40, 80, 127)):
        b = index * 3
        p.hand(b, 0)
        p.hit(b, ODAIKO, DON, 116)
        p.hand(b + 1, amount)
    p.hand(12, 0)
    for beat in range(4):
        p.hit(12 + beat, ODAIKO, DON, 120)
        p.hand(12 + beat + 0.5, 127)
        p.hand(12 + beat + 0.99, 0)
    p.hand(16, 0)
    return p


def pitch_wheel() -> Phrase:
    # The wheel presses the head: a held o-daiko note swept sharp and back,
    # nagado strokes on a pressed head, then the wheel dipped flat.
    p = Phrase(90)
    p.hit(0, ODAIKO, DON, 120)
    for step in range(25):
        p.bend(0.5 + step * 0.12, (1 - abs(step - 12) / 12.0))
    p.bend(3.6, 0)
    for index, amount in enumerate((0.0, 0.33, 0.66, 1.0, 0.66, 0.33, 0.0)):
        p.bend(5 + index, amount)
        p.hit(5 + index, NAGADO, DON, 104)
    p.hit(13, ODAIKO, DON, 120)
    for step in range(13):
        p.bend(13.5 + step * 0.15, -step / 12.0)
    p.bend(16, 0)
    return p


def strike_placement() -> Phrase:
    # CC17 walks the stroke from the centre to the rim, CC16 then carries it
    # once around the head, and CC18 answers on the rear head.
    p = Phrase(96)
    for index in range(9):
        p.cc(index * 0.75, 17, min(127, index * 16))
        p.hit(index * 0.75, NAGADO, DON, 108)
    p.cc(7.5, 17, 96)
    for index in range(9):
        p.cc(7.5 + index * 0.5, 16, min(127, index * 16))
        p.hit(7.5 + index * 0.5, NAGADO, DON, 104)
    p.cc(13, 121, 0)
    for index in range(8):
        p.cc(13 + index * 0.75, 18, 127 if index % 2 else 0)
        p.hit(13 + index * 0.75, ODAIKO, DON, 112)
    p.cc(19, 18, 0)
    return p


MORE = [
    ("rolls-and-flams", rolls_and_flams),
    ("dynamics", dynamics),
    ("hand-damping", hand_damping),
    ("pitch-wheel", pitch_wheel),
    ("strike-placement", strike_placement),
]


# Same order as Source/DSP/FactoryPresets.cpp.
PHRASES = [
    ("kumi-daiko-stage", kumi_daiko_stage),
    ("festival-ensemble", festival_ensemble),
    ("close-studio", close_studio),
    ("tight-ji-uchi", tight_ji_uchi),
    ("soft-bachi", soft_bachi),
    ("single-o-daiko", single_odaiko),
    ("synthetic-heads", synthetic_heads),
    ("distant-pair", distant_pair),
]


def main() -> int:
    out = pathlib.Path(sys.argv[1] if len(sys.argv) > 1 else "Docs/audio/midi")
    out.mkdir(parents=True, exist_ok=True)
    for directory, phrases in ((out, PHRASES), (out / "more", MORE)):
        directory.mkdir(parents=True, exist_ok=True)
        for number, (slug, phrase) in enumerate(phrases, start=1):
            path = directory / f"{number:02}-{slug}.mid"
            path.write_bytes(phrase().smf())
            print(f"Wrote {path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
