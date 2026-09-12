#!/usr/bin/env python3
"""Alternate reference/candidate runs and require identical float32 captures.

Build both executables in Release with the same compiler/options, and preserve
reference binaries before editing DSP. Timing runs are deliberately sequential.
"""

import argparse
import array
import csv
import json
import math
from pathlib import Path
import statistics
import subprocess
import sys


def compare_audio(reference, candidate):
    reference_files = {p.name: p for p in reference.glob("*.f32le")}
    candidate_files = {p.name: p for p in candidate.glob("*.f32le")}
    if not reference_files or reference_files.keys() != candidate_files.keys():
        raise ValueError("Capture inventory differs or is empty")
    results = {}
    for name, reference_path in sorted(reference_files.items()):
        before = reference_path.read_bytes()
        after = candidate_files[name].read_bytes()
        if len(before) != len(after) or len(before) % 8:
            raise ValueError(f"Invalid stereo capture length: {name}")
        samples_a, samples_b = array.array("f"), array.array("f")
        samples_a.frombytes(before)
        samples_b.frombytes(after)
        if sys.byteorder != "little":
            samples_a.byteswap()
            samples_b.byteswap()
        if not all(math.isfinite(v) for v in samples_a) or not all(
            math.isfinite(v) for v in samples_b
        ):
            raise ValueError(f"Non-finite audio: {name}")
        exact = before == after
        peak_error = 0.0
        snr = "infinite"
        if not exact:
            errors = [float(b) - float(a) for a, b in zip(samples_a, samples_b)]
            peak_error = max(map(abs, errors), default=0.0)
            error_energy = math.fsum(e * e for e in errors)
            reference_energy = math.fsum(float(a) * a for a in samples_a)
            if error_energy:
                snr = 10.0 * math.log10(reference_energy / error_energy) if reference_energy else "-infinite"
        results[name] = {
            "bit_exact": exact,
            "frames": len(before) // 8,
            "peak_error": peak_error,
            "peak_error_dbfs": 20.0 * math.log10(peak_error) if peak_error else "-infinite",
            "snr_db": snr,
        }
    return results


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--baseline", type=Path, required=True, help="Directory of preserved executables")
    parser.add_argument("--candidate", type=Path, required=True, help="Directory of candidate executables")
    parser.add_argument("--output", type=Path, required=True, help="New results directory")
    parser.add_argument("--config", action="append", choices=["48000:64", "48000:256", "96000:64", "96000:256"])
    parser.add_argument("--program", choices=["all", "dsp", "path"], default="all")
    parser.add_argument("--case", default="all")
    parser.add_argument("--seconds", type=float, default=3.0)
    parser.add_argument("--tail", type=float, default=1.0)
    parser.add_argument("--warmup", type=float, default=0.25)
    parser.add_argument("--repeats", type=int, default=3)
    parser.add_argument("--rounds", type=int, default=2)
    parser.add_argument("--seed", type=int, default=1)
    args = parser.parse_args()
    if args.rounds < 1 or args.rounds > 10:
        parser.error("--rounds must be 1..10")
    args.output.mkdir(parents=True, exist_ok=False)
    configs = args.config or ["48000:64", "48000:256", "96000:64", "96000:256"]
    programs = {"dsp": "TaikorBenchmarkCPU", "path": "TaikorBenchmarkAudioPath"}
    if args.program != "all":
        programs = {args.program: programs[args.program]}
    measurements = {}
    audio_results = {}
    for config in configs:
        rate, block = config.split(":")
        for program in programs.values():
            for round_index in range(args.rounds):
                # A/B, then B/A: reduce bias from warmup, clock and temperature.
                order = ["baseline", "candidate"] if round_index % 2 == 0 else ["candidate", "baseline"]
                captures = {}
                for version in order:
                    label = f"{program}-{rate}-b{block}-r{round_index}-{version}"
                    capture = args.output / label
                    executable = getattr(args, version).resolve() / program
                    if sys.platform == "win32":
                        executable = executable.with_suffix(".exe")
                    command = [str(executable), "--rate", rate, "--block", block,
                               "--case", args.case, "--seconds", str(args.seconds),
                               "--tail", str(args.tail), "--repeats", str(args.repeats),
                               "--warmup", str(args.warmup), "--seed", str(args.seed),
                               "--capture", str(capture.resolve())]
                    print(f"Running {label}", flush=True)
                    result = subprocess.run(command, capture_output=True, text=True, check=True)
                    (args.output / f"{label}.csv").write_text(result.stdout)
                    (args.output / f"{label}.stderr").write_text(result.stderr)
                    rows = list(csv.DictReader(result.stdout.splitlines()))
                    if not rows:
                        raise ValueError(f"No benchmark rows: {label}")
                    for row in rows:
                        key = f"{row['case']}:{rate}:{block}"
                        entry = measurements.setdefault(key, {"baseline": [], "candidate": []})
                        entry[version].append(row)
                    captures[version] = capture
                comparisons = compare_audio(captures["baseline"], captures["candidate"])
                audio_results[f"{program}:{config}:round{round_index}"] = comparisons
                print(f"  Audio: {sum(r['bit_exact'] for r in comparisons.values())}/{len(comparisons)} bit-exact", flush=True)
    summaries = {}
    for key, versions in measurements.items():
        cpu = {v: statistics.median(float(r["median_cpu_percent"]) for r in rows)
               for v, rows in versions.items()}
        reduction = (1.0 - cpu["candidate"] / cpu["baseline"]) * 100.0 if cpu["baseline"] else 0.0
        summaries[key] = {
            "baseline_cpu_percent": cpu["baseline"],
            "candidate_cpu_percent": cpu["candidate"],
            "reduction_percent": reduction,
            "baseline_p95_block_us": statistics.median(float(r["p95_block_us"]) for r in versions["baseline"]),
            "candidate_p95_block_us": statistics.median(float(r["p95_block_us"]) for r in versions["candidate"]),
        }
        print(f"{key:32} {cpu['baseline']:9.4f}% -> {cpu['candidate']:9.4f}% CPU ({reduction:+.1f}% reduction)")
    exact = all(r["bit_exact"] for group in audio_results.values() for r in group.values())
    report = {"all_audio_bit_exact": exact, "options": vars(args),
              "cpu": summaries, "audio": audio_results}
    (args.output / "summary.json").write_text(json.dumps(report, indent=2, default=str) + "\n")
    if not exact:
        print("FAIL: audio changed; inspect summary.json before accepting the optimization", file=sys.stderr)
        return 1
    print("PASS: every baseline/candidate audio capture is bit-exact")
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except (ValueError, OSError, subprocess.CalledProcessError) as error:
        print(f"BenchmarkCompare: {error}", file=sys.stderr)
        if isinstance(error, subprocess.CalledProcessError):
            print(error.stderr, file=sys.stderr)
        sys.exit(1)
