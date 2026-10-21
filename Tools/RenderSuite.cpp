// HIBIKI / Resonance — a one-minute performance of the shipping Taikor model.
// Four independent instances let palms, head choice and pitch pressure belong
// to one player while the others keep ringing. No alternate realism switches,
// samples, equalisation or compression are used here. ProduceSuite.py adds the
// documented room and delivery mastering. The default is the host's maximum
// supported 384 kHz, with unquantised float32 output and conservative gain.
//
// Score: 120 BPM, 28 bars in 4/4, final downbeat at 55 s, tail through 60 s.
// MIDI reproduces notes/CCs; CSV also records host-only ensemble automation.

#include "DSP/EnsembleEngine.h"

#include <algorithm>
#include <array>
#include <bit>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <future>
#include <iomanip>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

namespace
{
using taikor::Articulation;
constexpr auto don = Articulation::Don;
constexpr auto edge = Articulation::Ka;
constexpr auto muted = Articulation::Tsu;
constexpr auto rim = Articulation::DonRim;
constexpr double seconds = 60.0;
constexpr int blockSize = 256;
constexpr int ppq = 1920;
constexpr double pi = 3.14159265358979323846;
constexpr std::array<const char*, 4> names { "O-daiko", "Nagado-daiko", "Okedo", "Shime" };
constexpr std::array<float, 4> stage { -0.08f, -0.32f, 0.32f, 0.12f };
constexpr std::array<float, 4> faders { 1.12f, 1.0f, 0.83f, 0.64f };

enum class Kind { hit, damping, bend, players };
struct Event
{
    int tick = 0;
    Kind kind = Kind::hit;
    Articulation articulation = don;
    int value = 0, position = 64, azimuth = 64;
    bool rear = false;
    std::string technique;
};
using Part = std::vector<Event>;
using Score = std::array<Part, 4>;

int encodeBipolar (double value)
{
    return std::clamp (static_cast<int> (std::lround (64 + value * (value < 0 ? 64 : 63))), 0, 127);
}
float decodeBipolar (int value)
{
    return (value - 64) / (value < 64 ? 64.0f : 63.0f);
}

// Quantise authored controllers and velocities exactly as the accompanying
// MIDI does. Deterministic microtiming never moves structural downbeats.
void hit (Score& score, int drum, double beat, Articulation articulation,
          float velocity, const char* technique = "groove", bool rear = false,
          float position = 0.0f, bool tight = false)
{
    const auto index = static_cast<int> (score[drum].size());
    const int hand = index % 2 == 0 ? 1 : -1;
    const double offset = tight ? 0.0 : 0.005 * std::sin (index * 2.399 + drum * 1.7);
    Event event;
    event.tick = std::max (0, static_cast<int> (std::lround ((beat + offset * 2.0) * ppq)));
    event.articulation = articulation;
    event.value = std::clamp (static_cast<int> (std::lround (velocity * 127)), 1, 127);
    event.position = encodeBipolar (position);
    event.azimuth = encodeBipolar (hand * 0.48 / pi);
    event.rear = rear;
    event.technique = technique;
    score[drum].push_back (event);
}

void control (Score& score, int drum, double beat, Kind kind, int value, const char* label)
{
    Event event;
    event.tick = static_cast<int> (std::lround (beat * ppq));
    event.kind = kind;
    event.value = value;
    event.technique = label;
    score[drum].push_back (event);
}

void flam (Score& score, int drum, double beat, float velocity, Articulation stroke = don)
{
    hit (score, drum, beat - 0.060, don, velocity * 0.50f, "flam grace", false, 0.04f, true);
    hit (score, drum, beat, stroke, velocity, "flam accent", false, 0.0f, true);
}

void press (Score& score, int drum, double beat, float velocity, bool rear = false)
{
    double gap = 0.043;
    for (int bounce = 0; bounce < 7; ++bounce)
    {
        hit (score, drum, beat, muted, velocity, "press contacts", rear, 0.08f, true);
        beat += gap;
        gap *= 0.85;
        velocity *= 0.87f;
    }
}

void theme (Score& score, int drum, double beat, float velocity, bool answer = false)
{
    // A recurring 3+3+2 eighth-note phrase. Its silhouette survives every
    // orchestration change; ghosts and the answering cadence are subordinate.
    hit (score, drum, beat, don, velocity, "theme", false, 0, true);
    hit (score, drum, beat + 0.75, muted, velocity * 0.47f, "ghost");
    hit (score, drum, beat + 1.5, don, velocity * 0.89f, "theme");
    hit (score, drum, beat + 2.5, edge, velocity * 0.66f, "edge answer");
    hit (score, drum, beat + 3, answer ? rim : don,
         velocity * (answer ? 0.70f : 0.95f), "theme cadence");
    if (answer)
        hit (score, drum, beat + 3.5, edge, velocity * 0.63f, "turnaround");
}

void pulse (Score& score, int drum, double beat, float velocity, int count = 8)
{
    for (int note = 0; note < count; ++note)
        hit (score, drum, beat + note * 0.5,
             note % 4 == 3 ? muted : edge,
             velocity * (note % 2 == 0 ? 1.0f : 0.76f), "alternating hands", false, 0, note == 0);
}

Score compose()
{
    Score score;
    // 0–8 s: a recognizable hook with air around the opening low-drum voice.
    hit (score, 0, 0.04, don, 0.97f, "opening", false, 0, true);
    hit (score, 1, 1.5, don, 0.81f, "answer");
    hit (score, 1, 3, edge, 0.68f, "answer");
    hit (score, 0, 4, don, 0.90f, "low anchor", false, 0, true);
    theme (score, 1, 4, 0.84f);
    hit (score, 2, 7.5, muted, 0.55f, "pickup");
    for (int bar = 2; bar < 4; ++bar)
    {
        const double beat = bar * 4.0;
        hit (score, 0, beat, don, 0.91f, "low anchor", false, 0, true);
        theme (score, 1, beat, 0.87f, bar == 3);
        pulse (score, 3, beat, 0.58f);
        hit (score, 2, beat + 2, don, 0.69f, "cross answer");
    }
    flam (score, 2, 15.5, 0.84f);

    // 8–16 s: interlocking groove, shaped accents and audible muted ghosts.
    for (int bar = 4; bar < 8; ++bar)
    {
        const double beat = bar * 4.0;
        hit (score, 0, beat, don, 0.92f, "low anchor", false, 0, true);
        hit (score, 0, beat + 3, bar % 2 == 0 ? muted : don, 0.69f, "low answer");
        theme (score, 1, beat, 0.88f, bar % 2 == 1);
        pulse (score, 3, beat, 0.64f);
        hit (score, 2, beat + 1, don, 0.76f, "interlock");
        hit (score, 2, beat + 2.75, muted, 0.44f, "ghost");
        hit (score, 2, beat + 3.5, edge, 0.66f, "interlock");
    }
    flam (score, 2, 31.75, 0.82f, rim);

    // 16–24 s: the high drums become soloists. Rear responses are played
    // against sparse low anchors so the changed resonating head is exposed.
    for (int bar = 8; bar < 12; ++bar)
    {
        const double beat = bar * 4.0;
        hit (score, 0, beat, don, 0.81f, "dialogue anchor", false, 0, true);
        hit (score, 1, beat + 2, muted, 0.58f, "muted support");
        for (int note = 0; note < 6; ++note)
        {
            const bool rear = note >= 3;
            hit (score, 2, beat + note * 0.5,
                 note % 3 == 2 ? edge : don, rear ? 0.88f : 0.82f,
                 rear ? "rear-head answer" : "front-head call", rear, 0.08f);
        }
        hit (score, 3, beat + 3, don, 0.76f, "shime answer");
        if (bar % 2 == 0)
            press (score, 3, beat + 3.5, 0.66f);
        else
            flam (score, 3, beat + 3.5, 0.68f, rim);
    }

    // 24–32 s: an intimate break. Damping belongs only to the O-daiko
    // instance; all the other players' decays remain physical and continuous.
    hit (score, 0, 48, don, 0.88f, "open before palm", false, 0, true);
    control (score, 0, 48.85, Kind::damping, 112, "palm choke");
    hit (score, 0, 49.5, don, 0.63f, "palm held");
    control (score, 0, 50.25, Kind::damping, 0, "palm release");
    hit (score, 0, 51, don, 0.75f, "open after release");
    for (int note = 0; note < 7; ++note)
        hit (score, 2, 52 + note * 0.5, note == 5 ? muted : don,
             0.32f + note * 0.055f, "position and velocity phrase", false,
             -0.30f + note * 0.105f);
    hit (score, 1, 56, don, 0.85f, "pitch pressure", false, 0, true);
    for (int step = 0; step <= 24; ++step)
    {
        const float bend = 0.24f * static_cast<float> (std::sin (pi * step / 24.0));
        control (score, 1, 56.18 + step * 0.06, Kind::bend,
                 static_cast<int> (std::lround (8192 + bend * 8191)), "pressure arch");
    }
    control (score, 1, 58, Kind::bend, 8192, "pressure released");
    hit (score, 1, 58.5, don, 0.58f, "soft answer");
    press (score, 2, 59.5, 0.68f);
    for (int note = 0; note < 8; ++note)
        hit (score, 3, 60 + note * 0.5, note % 2 ? muted : edge,
             0.40f + note * 0.037f, "pickup crescendo");

    // 32–40 s: the theme returns. A measured accelerando grows from eighths
    // through triplets to sixteenths; strike position travels outwards with it.
    for (int drum = 0; drum < 4; ++drum)
        control (score, drum, 64, Kind::players, 2, "two players per part");
    for (int bar = 16; bar < 20; ++bar)
    {
        const double beat = bar * 4.0;
        hit (score, 0, beat, don, 0.87f, "rebuild anchor", false, 0, true);
        theme (score, 1, beat, 0.82f + (bar - 16) * 0.025f, bar % 2 == 1);
        pulse (score, 3, beat, 0.59f + (bar - 16) * 0.025f);
        const int notes = bar == 16 ? 8 : bar < 19 ? 12 : 16;
        for (int note = 0; note < notes; ++note)
        {
            const float progress = static_cast<float> (bar - 16 + note / static_cast<double> (notes)) / 4;
            hit (score, 2, beat + note * 4.0 / notes, note % 4 == 3 ? edge : don,
                 0.48f + progress * 0.34f, "accelerating roll", false,
                 -0.08f + progress * 0.38f);
        }
    }

    // 40–48 s: the four-part ensemble opens out. Each of the 16 grid voices
    // now has a musical role; rimshots are phrase punctuation, not every beat.
    for (int drum = 0; drum < 4; ++drum)
        control (score, drum, 80, Kind::players, drum == 3 ? 2 : 4, "tutti players");
    for (int bar = 20; bar < 24; ++bar)
    {
        const double beat = bar * 4.0;
        hit (score, 0, beat, don, 0.97f, "tutti anchor", false, 0, true);
        hit (score, 0, beat + 1.5, edge, 0.80f, "low edge answer");
        hit (score, 0, beat + 2.5, muted, 0.70f, "low muted answer");
        hit (score, 0, beat + 3, bar % 2 ? rim : don,
             bar % 2 ? 0.61f : 0.90f, "low cadence");
        theme (score, 1, beat, 0.94f, bar % 2 == 1);
        theme (score, 2, beat + 0.5, 0.82f, bar % 2 == 0);
        pulse (score, 3, beat, 0.70f, 6);
        hit (score, 3, beat + 3, don, 0.87f, "high accent");
        if (bar % 2)
            hit (score, 3, beat + 3.5, rim, 0.64f, "high rim cadence");
    }

    // 48–56 s: compressed reprise, answered rolls, two stop-time gestures,
    // a breath, then one final flam/unison. Last stroke at exactly 55 seconds.
    theme (score, 1, 96, 0.95f, true);
    theme (score, 2, 96.5, 0.87f, true);
    hit (score, 0, 96, don, 0.98f, "coda anchor", false, 0, true);
    pulse (score, 3, 96, 0.71f);
    hit (score, 0, 99, rim, 0.66f, "coda accent");
    press (score, 2, 99.5, 0.74f);
    hit (score, 0, 100, don, 0.93f, "roll anchor", false, 0, true);
    for (int note = 0; note < 24; ++note)
    {
        const double beat = 100 + note * (note < 16 ? 0.25 : 0.0);
        const double time = note < 16 ? beat : 104 + (note - 16) * 0.125;
        const int drum = note % 4 < 2 ? 1 : 2;
        hit (score, drum, time, note % 4 == 3 ? edge : don,
             0.56f + note * 0.015f, "final rising roll", false, 0.05f, true);
    }
    for (const double beat : { 105.0, 106.5 })
        for (int drum = 0; drum < 4; ++drum)
            hit (score, drum, beat, drum == 1 ? rim : don,
                 drum == 1 ? 0.67f : 0.93f, "stop time", false, 0, true);
    for (int drum = 0; drum < 4; ++drum)
    {
        flam (score, drum, 110, drum == 3 ? 0.86f : 0.99f);
        control (score, drum, 112, Kind::bend, 8192, "neutral tail");
    }
    for (auto& part : score)
    {
        std::stable_sort (part.begin(), part.end(), [] (const Event& a, const Event& b)
        { return a.tick < b.tick; });
        int hand = 1;
        for (auto& event : part)
            if (event.kind == Kind::hit)
            {
                event.azimuth = encodeBipolar (hand * 0.48 / pi);
                hand = -hand;
            }
    }
    return score;
}

void validate (const Score& score)
{
    std::array<int, 4> controllers {};
    int rearHits = 0;
    for (int drum = 0; drum < 4; ++drum)
    {
        std::array<int, 4> strokes {};
        for (const auto& event : score[drum])
        {
            if (event.tick < 0 || event.tick >= 120 * ppq)
                throw std::runtime_error ("score event outside one-minute duration");
            ++controllers[static_cast<int> (event.kind)];
            if (event.kind == Kind::hit)
            {
                ++strokes[static_cast<int> (event.articulation)];
                rearHits += event.rear ? 1 : 0;
                if (event.tick > 110 * ppq)
                    throw std::runtime_error ("last hit leaves insufficient natural tail");
            }
        }
        if (std::any_of (strokes.begin(), strokes.end(), [] (int count) { return count == 0; }))
            throw std::runtime_error ("score does not expose every grid voice");
        std::cout << names[drum] << ": Don " << strokes[0] << ", Edge " << strokes[1]
                  << ", Muted " << strokes[2] << ", Rimshot " << strokes[3] << '\n';
    }
    if (rearHits == 0 || controllers[1] == 0 || controllers[2] == 0 || controllers[3] == 0)
        throw std::runtime_error ("missing expressive performance controls");
}

struct Audio { std::vector<float> interleaved; };

Audio renderPart (const Part& score, int drum, int rate)
{
    auto engine = std::make_unique<taikor::EnsembleEngine>();
    taikor::EngineParameters parameters;
    parameters.outputGain = 0.025f;
    parameters.performer = drum;
    parameters.humanise = 0.22f;
    parameters.ensembleVariation = 0.32f;
    // The original microphone image remains intact before modest stage panning.
    engine->prepare (rate, blockSize);
    engine->setParameters (parameters);
    engine->reset();
    const auto pan = taikor::StereoPan::atPosition (stage[drum]);
    Audio audio;
    audio.interleaved.resize (static_cast<std::size_t> (seconds * rate) * 2);
    std::array<float, blockSize> left {}, right {};
    std::size_t cursor = 0, eventIndex = 0;
    const auto length = audio.interleaved.size() / 2;
    const auto sampleAt = [rate] (const Event& e)
    { return static_cast<std::size_t> (std::llround (e.tick * 0.5 / ppq * rate)); };
    while (cursor < length)
    {
        while (eventIndex < score.size() && sampleAt (score[eventIndex]) <= cursor)
        {
            const auto& e = score[eventIndex++];
            switch (e.kind)
            {
                case Kind::hit:
                    engine->setStrikePositionOverride (decodeBipolar (e.position));
                    engine->setStrikeAzimuthOverride (decodeBipolar (e.azimuth) * static_cast<float> (pi));
                    engine->setRearHeadStrike (e.rear);
                    engine->trigger (e.articulation, drum, e.value / 127.0f);
                    break;
                case Kind::damping: engine->setHandDamping (e.value / 127.0f); break;
                case Kind::bend:
                    engine->setPitchBend ((e.value - 8192) / 8192.0f);
                    break;
                case Kind::players:
                    parameters.ensembleSize = e.value;
                    engine->setParameters (parameters);
                    break;
            }
        }
        auto end = std::min (length, cursor + blockSize);
        if (eventIndex < score.size())
            end = std::min (end, sampleAt (score[eventIndex]));
        const auto count = static_cast<int> (end - cursor);
        engine->process (left.data(), right.data(), count);
        for (int frame = 0; frame < count; ++frame)
        {
            if (! std::isfinite (left[frame]) || ! std::isfinite (right[frame])
                || std::max (std::abs (left[frame]), std::abs (right[frame])) >= 0.88f)
                throw std::runtime_error ("invalid sample or insufficient engine limiter headroom");
            pan.apply (left[frame], right[frame]);
            audio.interleaved[2 * (cursor + frame)] = left[frame] * faders[drum];
            audio.interleaved[2 * (cursor + frame) + 1] = right[frame] * faders[drum];
        }
        cursor = end;
    }
    return audio;
}

void little (std::ostream& output, std::uint32_t value, int bytes)
{
    for (int i = 0; i < bytes; ++i)
        output.put (static_cast<char> ((value >> (8 * i)) & 255));
}

void writeAudio (const std::filesystem::path& path, const Audio& audio, int rate)
{
    std::ofstream file (path, std::ios::binary);
    const auto bytes = static_cast<std::uint32_t> (audio.interleaved.size() * sizeof (float));
    file.write ("RIFF", 4); little (file, bytes + 48, 4); file.write ("WAVEfmt ", 8);
    little (file, 16, 4); little (file, 3, 2); little (file, 2, 2);
    little (file, rate, 4); little (file, rate * 8, 4);
    little (file, 8, 2); little (file, 32, 2);
    file.write ("fact", 4); little (file, 4, 4);
    little (file, static_cast<std::uint32_t> (audio.interleaved.size() / 2), 4);
    file.write ("data", 4); little (file, bytes, 4);
    if constexpr (std::endian::native == std::endian::little)
        file.write (reinterpret_cast<const char*> (audio.interleaved.data()), bytes);
    else
        for (const float value : audio.interleaved)
            little (file, std::bit_cast<std::uint32_t> (value), 4);
    file.close();
    if (! file)
        throw std::runtime_error ("failed to write float WAV");
}

void writeScore (const std::filesystem::path& directory, const Score& score)
{
    std::ofstream csv (directory / "taikor-suite-score.csv");
    csv << "seconds,instance,drum,event,stroke,value,position_cc17,azimuth_cc16,rear,technique\n";
    constexpr std::array<const char*, 4> kinds { "note", "CC1", "pitch_wheel", "ensemble_size" };
    for (int drum = 0; drum < 4; ++drum)
        for (const auto& e : score[drum])
            csv << std::fixed << std::setprecision (6) << e.tick * 0.5 / ppq << ',' << drum + 1
                << ',' << names[drum] << ',' << kinds[static_cast<int> (e.kind)] << ','
                << (e.kind == Kind::hit ? taikor::getArticulationDisplayName (e.articulation) : "")
                << ',' << e.value << ',' << e.position << ',' << e.azimuth << ','
                << (e.rear ? 1 : 0) << ',' << e.technique << '\n';
    csv.close();
    if (! csv)
        throw std::runtime_error ("failed to write score CSV");
}

using Bytes = std::vector<std::uint8_t>;
void variable (Bytes& out, std::uint32_t value)
{
    std::array<std::uint8_t, 4> buffer {};
    int count = 0;
    buffer[count++] = static_cast<std::uint8_t> (value & 127);
    while ((value >>= 7) != 0)
        buffer[count++] = static_cast<std::uint8_t> ((value & 127) | 128);
    while (count > 0)
        out.push_back (buffer[--count]);
}
void big (std::ostream& out, std::uint32_t value, int count)
{
    for (int i = count - 1; i >= 0; --i)
        out.put (static_cast<char> ((value >> (i * 8)) & 255));
}
void midiChunk (std::ostream& file, const Bytes& bytes)
{
    file.write ("MTrk", 4); big (file, static_cast<std::uint32_t> (bytes.size()), 4);
    file.write (reinterpret_cast<const char*> (bytes.data()), static_cast<std::streamsize> (bytes.size()));
}
Bytes meta (int type, const std::string& message)
{
    Bytes bytes { 0xff, static_cast<std::uint8_t> (type) };
    variable (bytes, static_cast<std::uint32_t> (message.size()));
    bytes.insert (bytes.end(), message.begin(), message.end());
    return bytes;
}
void writeMidi (const std::filesystem::path& directory, const Score& score)
{
    std::ofstream file (directory / "taikor-suite.mid", std::ios::binary);
    file.write ("MThd", 4); big (file, 6, 4); big (file, 1, 2); big (file, 5, 2); big (file, ppq, 2);
    Bytes conductor { 0, 0xff, 0x51, 3, 7, 0xa1, 0x20, 0, 0xff, 0x58, 4, 4, 2, 24, 8 };
    const std::array<const char*, 8> sections {
        "HIBIKI - Opening", "Interlocking groove", "Front and rear dialogue", "Intimate break",
        "Rebuild", "Tutti", "Coda", "Resonance" };
    for (int i = 0; i < 8; ++i)
    {
        variable (conductor, i == 0 ? 0 : 16 * ppq);
        const auto bytes = meta (6, sections[i]); conductor.insert (conductor.end(), bytes.begin(), bytes.end());
    }
    variable (conductor, 8 * ppq); conductor.insert (conductor.end(), { 0xff, 0x2f, 0 });
    midiChunk (file, conductor);
    for (int drum = 0; drum < 4; ++drum)
    {
        struct Message { int tick; Bytes bytes; };
        std::vector<Message> messages;
        messages.push_back ({ 0, meta (3, names[drum]) });
        messages.push_back ({ 0, meta (1, "Separate Taikor instance; P" + std::to_string (drum + 1)
            + "; Humanise 22%; Variation 32%; outputGain .025; other defaults. Ensemble automation in CSV.") });
        const auto cc = [&messages, drum] (int tick, int number, int value)
        { messages.push_back ({ tick, { static_cast<std::uint8_t> (0xb0 + drum),
                       static_cast<std::uint8_t> (number), static_cast<std::uint8_t> (value) } }); };
        cc (0, 121, 0);
        for (const auto& e : score[drum])
        {
            if (e.kind == Kind::hit)
            {
                cc (e.tick, 17, e.position); cc (e.tick, 16, e.azimuth); cc (e.tick, 18, e.rear ? 127 : 0);
                const auto note = static_cast<std::uint8_t> (taikor::midiNoteFor (e.articulation, drum));
                messages.push_back ({ e.tick, { static_cast<std::uint8_t> (0x90 + drum), note,
                                               static_cast<std::uint8_t> (e.value) } });
                messages.push_back ({ e.tick + 6, { static_cast<std::uint8_t> (0x80 + drum), note, 0 } });
            }
            else if (e.kind == Kind::damping)
                cc (e.tick, 1, e.value);
            else if (e.kind == Kind::bend)
                messages.push_back ({ e.tick, { static_cast<std::uint8_t> (0xe0 + drum),
                    static_cast<std::uint8_t> (e.value & 127), static_cast<std::uint8_t> (e.value >> 7) } });
            else
                messages.push_back ({ e.tick, meta (1, "Set Ensemble Size = " + std::to_string (e.value)) });
        }
        std::stable_sort (messages.begin(), messages.end(), [] (const Message& a, const Message& b)
        { return a.tick < b.tick; });
        Bytes track;
        int last = 0;
        for (const auto& message : messages)
        {
            variable (track, message.tick - last);
            track.insert (track.end(), message.bytes.begin(), message.bytes.end());
            last = message.tick;
        }
        variable (track, 120 * ppq - last); track.insert (track.end(), { 0xff, 0x2f, 0 });
        midiChunk (file, track);
    }
    file.close();
    if (! file) throw std::runtime_error ("failed to write MIDI");
}
} // namespace

int main (int argc, char** argv)
{
    try
    {
        int rate = 384000;
        bool scoreOnly = false;
        std::filesystem::path directory = "build-demo-suite";
        for (int arg = 1; arg < argc; ++arg)
        {
            const std::string value = argv[arg];
            if (value == "--rate" && arg + 1 < argc) rate = std::stoi (argv[++arg]);
            else if (value == "--score-only") scoreOnly = true;
            else if (value == "--help")
            {
                std::cout << "TaikorRenderSuite [--rate 8000..384000] [--score-only] [output-directory]\n";
                return 0;
            }
            else if (value.starts_with ("--")) throw std::runtime_error ("unknown or incomplete option: " + value);
            else directory = value;
        }
        if (rate < 8000 || rate > 384000) throw std::runtime_error ("rate outside engine's supported range");
        const auto score = compose();
        validate (score);
        std::filesystem::create_directories (directory);
        writeScore (directory, score);
        writeMidi (directory, score);
        if (scoreOnly) return 0;
        std::cout << "Rendering four independent instances at " << rate << " Hz..." << std::endl;
        std::array<std::future<Audio>, 4> tasks;
        for (int drum = 0; drum < 4; ++drum)
            tasks[drum] = std::async (std::launch::async, [&, drum] { return renderPart (score[drum], drum, rate); });
        auto mix = tasks[0].get();
        std::cout << names[0] << " rendered." << std::endl;
        for (int drum = 1; drum < 4; ++drum)
        {
            const auto part = tasks[drum].get();
            for (std::size_t i = 0; i < mix.interleaved.size(); ++i)
                mix.interleaved[i] += part.interleaved[i];
            std::cout << names[drum] << " rendered." << std::endl;
        }
        float peak = 0;
        for (const auto sample : mix.interleaved)
        {
            if (! std::isfinite (sample)) throw std::runtime_error ("nonfinite mix");
            peak = std::max (peak, std::abs (sample));
        }
        if (peak < 0.0001f || peak >= 0.88f) throw std::runtime_error ("mix silent or insufficient headroom");
        writeAudio (directory / "taikor-suite-dry.wav", mix, rate);
        std::cout << "60.000 s stereo float32 dry master; peak " << 20 * std::log10 (peak) << " dBFS.\n";
        return 0;
    }
    catch (const std::exception& error)
    {
        std::cerr << "Taikor suite: " << error.what() << '\n';
        return 1;
    }
}
