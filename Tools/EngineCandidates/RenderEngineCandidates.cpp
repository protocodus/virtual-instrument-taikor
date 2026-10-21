// Renders one stroke of one drum through the shipping engine (letter A) and
// through alternative drum engines (B, C, D, ...) for a blind A-Z listening
// test, following the repository's CLAUDE.md protocol:
//
// - A is the shipping engine, unchanged, rendered through EnsembleEngine at
//   Ensemble Size 1 exactly as the plug-in renders it.
// - Every letter plays the same MIDI note, velocity, reset state, sample rate
//   (48 kHz), block size (256), length and controls. The candidates are given
//   the drum, stroke and capsule positions the shipping engine resolves for
//   those controls (ShippingSpec.cpp), including the stroke's seeded speed and
//   contact-stiffness variation.
// - Each candidate's microphone pressure is trimmed by one constant so its
//   whole-take stereo RMS equals the shipping engine's raw pair, then passed
//   through the shipping output stage (TaikoEngine::processWithEnsemble on an
//   idle engine, verified here to reproduce A bit for bit from A's own raw
//   pair). The limiter must not engage for any letter; the run fails if it
//   does.
// - The finished files are level-matched on whole-file stereo RMS (all
//   samples of both channels pooled) with one constant gain each, to a common
//   target low enough that every peak stays under -1 dBFS, and written as
//   24-bit PCM named by letter only. key.md records which engine each letter
//   is, every trim, and each engine's diagnostics; it is written for the
//   listener to read only after deciding.
//
// Usage:
//   TaikorRenderEngineCandidates --out DIR --letters B=modal,C=berger,D=vonkarman
//       [--note 60] [--velocity 127] [--seconds 4.5] [--date YYYY-MM-DD]
//       [--reference R=take.f32 [--reference-note provenance.md]]
//   TaikorRenderEngineCandidates --out DIR --dev ENGINE [--seconds S]
//       renders one engine alone as dev-ENGINE.wav (development only).
// Engine ids: berger, vonkarman, modal.

#include "CandidateEngine.h"
#include "FdtdContracts.h"
#include "ShippingSpec.h"

#include "DSP/OutputLimiter.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

namespace
{
using namespace taikor::candidates;
namespace fs = std::filesystem;

struct Letter
{
    char letter = 'A';
    std::string engineId;          // "shipping" for A
    std::string description;
    std::string diagnostics;
    Take raw;                      // microphone pair before the output stage
    Take out;                      // after the shipping output stage
    double rawTrimDb = 0.0;
    double prepareSeconds = 0.0, renderSeconds = 0.0;
    double finalGainDb = 0.0;
};

double rms (const Take& take)
{
    double sum = 0.0;
    for (std::size_t i = 0; i < take.left.size(); ++i)
        sum += double (take.left[i]) * take.left[i] + double (take.right[i]) * take.right[i];
    return std::sqrt (sum / (2.0 * static_cast<double> (take.left.size())));
}

double peak (const Take& take)
{
    double result = 0.0;
    for (std::size_t i = 0; i < take.left.size(); ++i)
        result = std::max ({ result, std::abs (double (take.left[i])), std::abs (double (take.right[i])) });
    return result;
}

bool finite (const Take& take)
{
    for (std::size_t i = 0; i < take.left.size(); ++i)
        if (! std::isfinite (take.left[i]) || ! std::isfinite (take.right[i]))
            return false;
    return true;
}

double db (double gain) { return 20.0 * std::log10 (gain); }

void scale (Take& take, double gain)
{
    for (auto& v : take.left) v = static_cast<float> (v * gain);
    for (auto& v : take.right) v = static_cast<float> (v * gain);
}

void appendLittleEndian (std::vector<std::uint8_t>& bytes, std::uint32_t value, int count)
{
    for (int i = 0; i < count; ++i)
        bytes.push_back (static_cast<std::uint8_t> ((value >> (8 * i)) & 0xffu));
}

// 24-bit PCM stereo. Samples are rounded, never clamped silently: the caller
// guarantees |x| < 1.
bool writeWav24 (const fs::path& path, const Take& take, double sampleRate, double gain)
{
    const auto frames = static_cast<std::uint32_t> (take.left.size());
    const std::uint32_t dataBytes = frames * 2u * 3u;
    std::vector<std::uint8_t> bytes;
    bytes.reserve (44u + dataBytes);
    const auto tag = [&bytes] (const char* text) { bytes.insert (bytes.end(), text, text + 4); };
    tag ("RIFF");
    appendLittleEndian (bytes, 36u + dataBytes, 4);
    tag ("WAVE");
    tag ("fmt ");
    appendLittleEndian (bytes, 16u, 4);
    appendLittleEndian (bytes, 1u, 2);
    appendLittleEndian (bytes, 2u, 2);
    appendLittleEndian (bytes, static_cast<std::uint32_t> (sampleRate), 4);
    appendLittleEndian (bytes, static_cast<std::uint32_t> (sampleRate) * 6u, 4);
    appendLittleEndian (bytes, 6u, 2);
    appendLittleEndian (bytes, 24u, 2);
    tag ("data");
    appendLittleEndian (bytes, dataBytes, 4);
    const auto encode = [gain] (float value)
    {
        const double x = double (value) * gain;
        const auto q = static_cast<std::int32_t> (std::lround (x * 8388607.0));
        return static_cast<std::uint32_t> (q) & 0xffffffu;
    };
    for (std::size_t i = 0; i < take.left.size(); ++i)
    {
        appendLittleEndian (bytes, encode (take.left[i]), 3);
        appendLittleEndian (bytes, encode (take.right[i]), 3);
    }
    std::ofstream file (path, std::ios::binary);
    file.write (reinterpret_cast<const char*> (bytes.data()), static_cast<std::streamsize> (bytes.size()));
    return static_cast<bool> (file);
}

bool writeFloat32 (const fs::path& path, const Take& take)
{
    std::ofstream file (path, std::ios::binary);
    for (std::size_t i = 0; i < take.left.size(); ++i)
    {
        file.write (reinterpret_cast<const char*> (&take.left[i]), 4);
        file.write (reinterpret_cast<const char*> (&take.right[i]), 4);
    }
    return static_cast<bool> (file);
}

// Each engine is compiled in only when its library is linked, so one engine's
// development build never depends on another's (see CMakeLists.txt).
std::unique_ptr<Engine> makeEngine ([[maybe_unused]] const std::string& id)
{
#if TAIKOR_WITH_FDTD
    if (id == "berger")
        return makeFdtdEngine (makeBergerPotential());
#if TAIKOR_WITH_VONKARMAN
    if (id == "vonkarman")
        if (auto potential = makeVonKarmanPotential())
            return makeFdtdEngine (std::move (potential));
#endif
#endif
#if TAIKOR_WITH_MODAL
    if (id == "modal")
        return makeModalEngine();
#endif
    return nullptr;
}

// Renders a candidate through the same schedule as A: pre-roll, strike on
// strikeSample, blockSize blocks.
bool renderCandidate (Letter& letter, const ShippingStroke& stroke,
                      const ResolvedShipping& shipping, std::string& error)
{
    auto engine = makeEngine (letter.engineId);
    if (! engine)
    {
        error = "engine '" + letter.engineId + "' is not built";
        return false;
    }
    letter.description = engine->description();
    const auto t0 = std::chrono::steady_clock::now();
    if (! engine->prepare (shipping.drum, shipping.mics, stroke.sampleRate, error))
        return false;
    const auto t1 = std::chrono::steady_clock::now();
    letter.raw.left.assign (static_cast<std::size_t> (stroke.totalSamples), 0.0f);
    letter.raw.right.assign (static_cast<std::size_t> (stroke.totalSamples), 0.0f);
    int position = 0;
    const auto run = [&] (int until)
    {
        while (position < until)
        {
            const int count = std::min (stroke.blockSize, until - position);
            engine->process (letter.raw.left.data() + position, letter.raw.right.data() + position, count);
            position += count;
        }
    };
    run (stroke.strikeSample);
    engine->strike (shipping.stroke);
    run (stroke.totalSamples);
    const auto t2 = std::chrono::steady_clock::now();
    letter.prepareSeconds = std::chrono::duration<double> (t1 - t0).count();
    letter.renderSeconds = std::chrono::duration<double> (t2 - t1).count();
    letter.diagnostics = engine->diagnostics();
    if (! finite (letter.raw))
    {
        error = "non-finite microphone output";
        return false;
    }
    if (rms (letter.raw) <= 0.0)
    {
        error = "silent microphone output";
        return false;
    }
    return true;
}

std::map<std::string, std::string> parseArguments (int argc, char** argv)
{
    std::map<std::string, std::string> result;
    for (int i = 1; i + 1 < argc; i += 2)
        result[argv[i]] = argv[i + 1];
    return result;
}

} // namespace

int main (int argc, char** argv)
{
    auto args = parseArguments (argc, argv);
    if (! args.count ("--out") || (! args.count ("--letters") && ! args.count ("--dev")))
    {
        std::fprintf (stderr, "usage: %s --out DIR (--letters B=id,C=id,... | --dev id) "
                              "[--note 60] [--velocity 127] [--seconds 4.5] [--date D]\n", argv[0]);
        return 2;
    }
    const fs::path outDir = args["--out"];
    fs::create_directories (outDir);

    ShippingStroke stroke;
    stroke.parameters = listeningTestParameters();
    stroke.midiNote = args.count ("--note") ? std::stoi (args["--note"]) : 60;
    stroke.velocity = (args.count ("--velocity") ? std::stof (args["--velocity"]) : 127.0f) / 127.0f;
    const double seconds = args.count ("--seconds") ? std::stod (args["--seconds"]) : 4.5;
    stroke.totalSamples = static_cast<int> (std::lround (seconds * stroke.sampleRate));
    stroke.strikeSample = static_cast<int> (std::lround (0.1 * stroke.sampleRate));

    ResolvedShipping shipping;
    try { shipping = resolveShipping (stroke); }
    catch (const std::exception& e)
    {
        std::fprintf (stderr, "resolve failed: %s\n", e.what());
        return 1;
    }
    std::printf ("%s", shipping.report.c_str());

    // A, and the proof that the output stage the candidates go through is the
    // shipping one: A's own raw pair must come out of it bit for bit.
    Letter a;
    a.letter = 'A';
    a.engineId = "shipping";
    a.description = "Shipping engine (TaikoEngine via EnsembleEngine), unchanged.";
    {
        const auto t0 = std::chrono::steady_clock::now();
        a.out = renderShipping (stroke);
        a.renderSeconds = std::chrono::duration<double> (std::chrono::steady_clock::now() - t0).count();
        a.raw = renderShippingRaw (stroke);
        const Take check = applyShippingOutputStage (stroke, a.raw);
        std::size_t differing = 0;
        for (std::size_t i = 0; i < check.left.size(); ++i)
            differing += (check.left[i] != a.out.left[i]) + (check.right[i] != a.out.right[i]);
        std::printf ("output-stage check: %zu of %zu samples differ from the shipping render\n",
                     differing, 2 * check.left.size());
        if (differing != 0)
        {
            std::fprintf (stderr, "the injected output stage does not reproduce the shipping engine\n");
            return 1;
        }
    }
    const double rawTarget = rms (a.raw);

    std::vector<Letter> letters;
    if (args.count ("--dev"))
    {
        Letter dev;
        dev.letter = '?';
        dev.engineId = args["--dev"];
        letters.push_back (std::move (dev));
    }
    else
    {
        letters.push_back (std::move (a));
        std::stringstream list (args["--letters"]);
        std::string item;
        while (std::getline (list, item, ','))
        {
            if (item.size() < 3 || item[1] != '=')
            {
                std::fprintf (stderr, "bad letter assignment '%s'\n", item.c_str());
                return 2;
            }
            Letter letter;
            letter.letter = item[0];
            letter.engineId = item.substr (2);
            letters.push_back (std::move (letter));
        }
    }

    // An optional real recording, labelled rather than blind (--reference
    // R=file.f32: interleaved stereo float32 at 48 kHz, totalSamples frames,
    // onset on strikeSample). It gets the same raw trim, output stage and
    // final RMS match as every engine, and no other processing.
    std::string referenceNote;
    if (args.count ("--reference") && ! args.count ("--dev"))
    {
        const std::string item = args["--reference"];
        if (item.size() < 3 || item[1] != '=')
        {
            std::fprintf (stderr, "bad reference '%s'\n", item.c_str());
            return 2;
        }
        std::ifstream file (item.substr (2), std::ios::binary);
        std::vector<float> interleaved (2 * static_cast<std::size_t> (stroke.totalSamples));
        file.read (reinterpret_cast<char*> (interleaved.data()),
                   static_cast<std::streamsize> (interleaved.size() * sizeof (float)));
        if (file.gcount() != static_cast<std::streamsize> (interleaved.size() * sizeof (float)))
        {
            std::fprintf (stderr, "reference must hold exactly %d stereo float32 frames\n", stroke.totalSamples);
            return 2;
        }
        Letter reference;
        reference.letter = item[0];
        reference.engineId = "reference";
        reference.description = "Real recording (labelled reference; see below).";
        for (std::size_t i = 0; i < interleaved.size(); i += 2)
        {
            reference.raw.left.push_back (interleaved[i]);
            reference.raw.right.push_back (interleaved[i + 1]);
        }
        if (args.count ("--reference-note"))
        {
            std::ifstream note (args["--reference-note"]);
            referenceNote.assign (std::istreambuf_iterator<char> (note), {});
        }
        letters.push_back (std::move (reference));
    }

    for (auto& letter : letters)
    {
        if (letter.engineId == "shipping")
            continue;
        if (letter.engineId == "reference")
        {
            const double trim = rawTarget / rms (letter.raw);
            letter.rawTrimDb = db (trim);
            Take trimmed = letter.raw;
            scale (trimmed, trim);
            letter.out = applyShippingOutputStage (stroke, trimmed);
            continue;
        }
        std::string error;
        std::printf ("rendering %c (%s)...\n", letter.letter, letter.engineId.c_str());
        std::fflush (stdout);
        if (! renderCandidate (letter, stroke, shipping, error))
        {
            std::fprintf (stderr, "%c (%s) failed: %s\n", letter.letter, letter.engineId.c_str(), error.c_str());
            return 1;
        }
        const double trim = rawTarget / rms (letter.raw);
        letter.rawTrimDb = db (trim);
        Take trimmed = letter.raw;
        scale (trimmed, trim);
        letter.out = applyShippingOutputStage (stroke, trimmed);
        std::printf ("  prepared in %.1f s, rendered in %.1f s; raw trim %+.4f dB\n%s\n",
                     letter.prepareSeconds, letter.renderSeconds, letter.rawTrimDb,
                     letter.diagnostics.c_str());
    }

    // The limiter passes samples under its ceiling unchanged, so a take whose
    // output never reaches it has been through a linear stage.
    for (const auto& letter : letters)
    {
        const double p = peak (letter.out);
        if (p >= 0.99 * taikor::OutputLimiter::ceiling)
        {
            std::fprintf (stderr, "%c reached the output limiter (peak %.4f)\n", letter.letter, p);
            return 1;
        }
    }

    // Common whole-file RMS target: -30 dBFS unless a take's crest factor
    // would then put its peak above -1 dBFS.
    double target = -30.0;
    for (const auto& letter : letters)
        target = std::min (target, -1.0 - db (peak (letter.out) / rms (letter.out)));
    for (auto& letter : letters)
        letter.finalGainDb = target - db (rms (letter.out));

    if (args.count ("--dev"))
    {
        const auto& dev = letters.front();
        writeWav24 (outDir / ("dev-" + dev.engineId + ".wav"), dev.out, stroke.sampleRate,
                    std::pow (10.0, dev.finalGainDb / 20.0));
        writeFloat32 (outDir / ("dev-" + dev.engineId + "-raw.f32"), dev.raw);
        writeFloat32 (outDir / "dev-shipping-raw.f32", a.raw);
        std::printf ("wrote %s\n", (outDir / ("dev-" + dev.engineId + ".wav")).string().c_str());
        return 0;
    }

    const fs::path analysis = outDir / "analysis";
    fs::create_directories (analysis);
    for (const auto& letter : letters)
    {
        const std::string name (1, letter.letter);
        if (! writeWav24 (outDir / (name + ".wav"), letter.out, stroke.sampleRate,
                          std::pow (10.0, letter.finalGainDb / 20.0)))
        {
            std::fprintf (stderr, "could not write %s.wav\n", name.c_str());
            return 1;
        }
        writeFloat32 (analysis / (name + "-raw.f32"), letter.raw);
        writeFloat32 (analysis / (name + "-out.f32"), letter.out);
    }

    std::ofstream key (outDir / "key.md");
    key << "# Key: engine-candidate listening test"
        << (args.count ("--date") ? " (" + args["--date"] + ")" : std::string()) << "\n\n"
        << "Read this only after choosing a letter by ear.\n\n"
        << "## Material\n\n"
        << "One stroke: MIDI note " << stroke.midiNote << " (nagado-daiko Don) at velocity "
        << std::lround (stroke.velocity * 127.0f) << "/127, struck "
        << stroke.strikeSample / stroke.sampleRate << " s into a "
        << stroke.totalSamples / stroke.sampleRate << " s take. 48 kHz, 256-sample blocks, "
        << "one reset state. Factory controls (Mic Distance 16 cm, Output -22.5 dB, Stick Noise "
        << stroke.parameters.strikeNoise << ", Drive 0, Low Cut off, Width 50 %, Performer P1, "
        << "Ensemble Size 1) with Humanise 0. The plug-in's reverb (default Off) is not part of "
        << "the engine and is not applied.\n\n"
        << "## What each letter is\n\n| Letter | Engine |\n| --- | --- |\n";
    for (const auto& letter : letters)
        key << "| " << letter.letter << " | " << letter.description << " |\n";
    key << "\n## Level matching\n\n"
        << "Each candidate's raw microphone pair was first trimmed by one constant to the shipping "
        << "engine's raw whole-take stereo RMS, then passed through the shipping output stage "
        << "(verified to reproduce A bit for bit from A's own raw pair); no take reached the "
        << "limiter. The finished takes were then matched on whole-file stereo RMS (both channels "
        << "pooled) to " << target << " dBFS with one constant gain each.\n\n"
        << "| Letter | Raw trim (dB) | Final gain (dB) | Final peak (dBFS) | Prepare (s) | Render (s) |\n"
        << "| --- | ---: | ---: | ---: | ---: | ---: |\n";
    key.precision (4);
    key << std::fixed;
    for (const auto& letter : letters)
        key << "| " << letter.letter << " | " << letter.rawTrimDb << " | " << letter.finalGainDb
            << " | " << db (peak (letter.out)) + letter.finalGainDb << " | " << letter.prepareSeconds
            << " | " << letter.renderSeconds << " |\n";
    if (! referenceNote.empty())
        key << "\n## The reference recording\n\n" << referenceNote << "\n";
    key << "\n## The drum and stroke every letter was given\n\n" << shipping.report
        << "\n## Engine diagnostics\n\n";
    for (const auto& letter : letters)
        if (! letter.diagnostics.empty())
            key << "### " << letter.letter << "\n\n" << letter.diagnostics << "\n\n";
    std::printf ("wrote %zu takes and key.md to %s (target %.2f dBFS RMS)\n", letters.size(),
                 outDir.string().c_str(), target);
    return 0;
}
