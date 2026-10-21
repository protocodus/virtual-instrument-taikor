#!/usr/bin/env python3
"""Render the one-minute Taikor suite with dynamic and listening masters.

Requires NumPy, SciPy and FFmpeg with libmp3lame. Example:
  python3 Tools/ProduceSuite.py --renderer build/TaikorRenderSuite

The physical-model performance and floating-point master run at 384 kHz. The
bundled theater convolution is calculated at 96 kHz and interpolated back into
the master; the report records this separately from the engine sample rate.
The 384 kHz master keeps the full performance dynamics. The listening editions
use a documented stereo-linked lookahead limiter for a louder presentation.
"""

import argparse
import csv
from datetime import datetime, timezone
import hashlib
import json
import math
from pathlib import Path
import re
import shutil
import struct
import subprocess
import sys
import tempfile

import numpy as np
import scipy
from scipy import ndimage, signal
from scipy.io import wavfile


ROOT = Path(__file__).resolve().parents[1]
MASTER_RATE = 384000
LISTENING_RATE = 96000
MP3_RATE = 48000
DURATION = 60.0
TARGET_TRUE_PEAK_DB = -1.2
CHUNK_FRAMES = 262144
RESAMPLE_WINDOW = ("kaiser", 12.0)
SUITE_SECTIONS = (("Opening", 0, 8), ("Groove", 8, 16),
                  ("Dialogue", 16, 24), ("Intimate break", 24, 32),
                  ("Rebuild", 32, 40), ("Tutti", 40, 48),
                  ("Coda", 48, 56), ("Resonance", 56, 60))


def db(value):
    return 20.0 * math.log10(value) if value > 0.0 else None


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def run(command):
    result = subprocess.run(command, capture_output=True, text=True)
    if result.returncode:
        raise RuntimeError(
            f"Command failed ({result.returncode}): {' '.join(map(str, command))}\n"
            f"{result.stdout[-4000:]}\n{result.stderr[-8000:]}"
        )
    return result


def read_audio(path, mmap=False):
    try:
        rate, samples = wavfile.read(path, mmap=mmap)
    except ValueError:
        # Three-byte PCM cannot be mapped by scipy.io.wavfile.
        rate, samples = wavfile.read(path)
    if samples.ndim == 1:
        samples = samples[:, None]
    if np.issubdtype(samples.dtype, np.integer):
        if samples.dtype == np.uint8:
            samples = (samples.astype(np.float32) - 128.0) / 128.0
        else:
            samples = samples.astype(np.float32) / float(
                2 ** (np.iinfo(samples.dtype).bits - 1)
            )
    elif samples.dtype != np.float32:
        samples = samples.astype(np.float32)
    return rate, samples


def wav_metadata(path):
    with Path(path).open("rb") as stream:
        header = stream.read(12)
        if header[:4] != b"RIFF" or header[8:] != b"WAVE":
            raise ValueError(f"Expected a RIFF WAVE file: {path}")
        while chunk := stream.read(8):
            if len(chunk) != 8:
                break
            tag, length = struct.unpack("<4sI", chunk)
            if tag == b"fmt ":
                payload = stream.read(length)
                code, channels, rate, byte_rate, alignment, bits = struct.unpack(
                    "<HHIIHH", payload[:16])
                if code == 0xFFFE and len(payload) >= 40:
                    code = struct.unpack("<H", payload[24:26])[0]
                return {"encoding": {1: "PCM", 3: "IEEE float"}.get(code, str(code)),
                        "bits_per_sample": bits, "channels": channels,
                        "sample_rate_hz": rate, "byte_rate": byte_rate,
                        "block_alignment_bytes": alignment}
            stream.seek(length + (length & 1), 1)
    raise ValueError(f"Missing WAV format chunk: {path}")


def resample(samples, source_rate, target_rate):
    if source_rate == target_rate:
        return np.array(samples, dtype=np.float32, copy=True)
    common = math.gcd(source_rate, target_rate)
    up, down = target_rate // common, source_rate // common
    factor = max(up, down)
    # 64 taps per polyphase branch, with a high-rejection Kaiser window.
    kernel = signal.firwin(64 * factor + 1, 1.0 / factor,
                           window=RESAMPLE_WINDOW).astype(np.float32)
    return signal.resample_poly(samples, up, down, axis=0, window=kernel,
                                padtype="constant").astype(np.float32, copy=False)


def audio_metrics(samples, rate):
    count, channels = samples.shape
    totals = np.zeros(channels, dtype=np.float64)
    squares = np.zeros(channels, dtype=np.float64)
    cross = 0.0
    peaks = np.zeros(channels, dtype=np.float64)
    clipped = nonfinite = 0
    for start in range(0, count, CHUNK_FRAMES):
        block = np.asarray(samples[start:start + CHUNK_FRAMES], dtype=np.float64)
        nonfinite += int(np.count_nonzero(~np.isfinite(block)))
        if not np.isfinite(block).all():
            block = np.nan_to_num(block, nan=0.0, posinf=0.0, neginf=0.0)
        totals += np.sum(block, axis=0)
        squares += np.einsum("ij,ij->j", block, block)
        peaks = np.maximum(peaks, np.max(np.abs(block), axis=0))
        clipped += int(np.count_nonzero(np.abs(block) >= 1.0))
        if channels == 2:
            cross += float(np.dot(block[:, 0], block[:, 1]))
    means = totals / max(count, 1)
    rms = np.sqrt(squares / max(count, 1))
    overall_rms = float(np.sqrt(np.sum(squares) / max(count * channels, 1)))
    correlation = None
    if channels == 2 and count:
        variance = np.maximum(squares - totals * means, 0.0)
        denominator = math.sqrt(float(variance[0] * variance[1]))
        if denominator > 0:
            correlation = float(np.clip((cross - totals[0] * means[1])
                                        / denominator, -1.0, 1.0))
    peak = float(np.max(peaks))
    return {
        "sample_rate_hz": rate,
        "channels": channels,
        "frames": count,
        "duration_seconds": count / rate,
        "sample_peak": peak,
        "sample_peak_dbfs": db(peak),
        "channel_peak_dbfs": [db(float(value)) for value in peaks],
        "rms_dbfs": db(overall_rms),
        "channel_rms_dbfs": [db(float(value)) for value in rms],
        "crest_factor_db": db(peak / overall_rms) if overall_rms else None,
        "channel_dc_offset": means.tolist(),
        "stereo_correlation": correlation,
        "clipped_sample_count": clipped,
        "nonfinite_sample_count": nonfinite,
        "first_frame": [float(value) if np.isfinite(value) else None
                        for value in samples[0]] if count else [],
        "last_frame": [float(value) if np.isfinite(value) else None
                       for value in samples[-1]] if count else [],
    }


def estimated_true_peak(samples, factor=8):
    """Oversample overlapping blocks; ignore artificial internal boundaries."""
    kernel = signal.firwin(64 * factor + 1, 1.0 / factor,
                           window=RESAMPLE_WINDOW).astype(np.float32)
    overlap = 128  # Exceeds the filter's 32 input-frame half width.
    maximum = 0.0
    for start in range(0, len(samples), 65536):
        stop = min(start + 65536, len(samples))
        lo, hi = max(0, start - overlap), min(len(samples), stop + overlap)
        for channel in range(samples.shape[1]):
            interpolated = signal.resample_poly(
                samples[lo:hi, channel], factor, 1, window=kernel,
                padtype="constant")
            core = interpolated[(start - lo) * factor:(stop - lo) * factor]
            maximum = max(maximum, float(np.max(np.abs(core))))
    return maximum


def limit_listening(samples, rate):
    """Apply one linked, anticipatory gain envelope without clipping samples.

    The quiet passages receive the same fixed drive as everything else. Only
    peaks requiring headroom change the envelope; neither loudness nor score
    sections control the gain. The final interpolated-peak check catches the
    small residual overshoot possible with an exponentially smoothed attack.
    """
    drive_db, ceiling_db = 8.0, TARGET_TRUE_PEAK_DB
    lookahead_seconds, attack_seconds, release_seconds = 0.003, 0.00025, 0.070
    peaks = np.max(np.abs(samples), axis=1)
    required = np.maximum(0.0, 20.0 * np.log10(np.maximum(peaks, 1.0e-12))
                          + drive_db - ceiling_db).astype(np.float32)
    del peaks
    lookahead_frames = 2 * round(lookahead_seconds * rate / 2)
    future = ndimage.maximum_filter1d(
        required, size=lookahead_frames + 1, origin=-(lookahead_frames // 2),
        mode="constant", cval=0.0)
    del required
    envelope = np.empty(len(samples), dtype=np.float32)
    attack = math.exp(-1.0 / (attack_seconds * rate))
    release = math.exp(-1.0 / (release_seconds * rate))
    previous = float(future[0])  # An opening transient gets a primed detector.
    for index, value in enumerate(future):
        target = float(value)
        coefficient = attack if target > previous else release
        previous = coefficient * previous + (1.0 - coefficient) * target
        envelope[index] = previous
    del future
    for start in range(0, len(samples), CHUNK_FRAMES):
        attenuation = envelope[start:start + CHUNK_FRAMES]
        gain = np.power(10.0, (drive_db - attenuation) / 20.0).astype(np.float32)
        samples[start:start + CHUNK_FRAMES] *= gain[:, None]
    peak = estimated_true_peak(samples)
    correction = min(1.0, 10.0 ** (ceiling_db / 20.0) / max(peak, 1.0e-30))
    samples *= np.float32(correction)
    maximum_reduction = float(np.max(envelope))
    stats = {
        "method": "Stereo-linked lookahead peak limiter; identical gain for both channels; no hard clipping",
        "presentation_gain_db": drive_db,
        "lookahead_seconds": lookahead_frames / rate,
        "attack_time_constant_seconds": attack_seconds,
        "release_time_constant_seconds": release_seconds,
        "detector": "Maximum absolute L/R sample, future-window maximum, asymmetric one-pole attenuation in dB",
        "nominal_maximum_reduction_db": 8.0,
        "maximum_limiter_reduction_db": maximum_reduction,
        "mean_limiter_reduction_db": float(np.mean(envelope, dtype=np.float64)),
        "active_fraction_above_0_1_db_reduction": float(np.mean(envelope > 0.1)),
        "final_true_peak_trim_db": db(correction),
        "maximum_total_reduction_including_peak_trim_db": maximum_reduction - db(correction),
        "target_true_peak_dbtp": ceiling_db,
        "section_gain_reduction": [],
    }
    for label, start, stop in SUITE_SECTIONS:
        section = envelope[round(start * rate):round(stop * rate)]
        if len(section):
            stats["section_gain_reduction"].append({
                "name": label, "maximum_db": float(np.max(section)),
                "mean_db": float(np.mean(section, dtype=np.float64)),
                "active_fraction_above_0_1_db": float(np.mean(section > 0.1))})
    if stats["maximum_total_reduction_including_peak_trim_db"] > 8.1:
        raise ValueError("Listening limiter exceeded its 8 dB reduction budget plus 0.1 dB peak-trim allowance")
    return stats


def ebur128(ffmpeg, path):
    result = run([ffmpeg, "-hide_banner", "-nostdin", "-i", str(path),
                  "-af", "ebur128=peak=true", "-f", "null", "-"])
    summary = result.stderr.rsplit("Summary:", 1)[-1]
    patterns = {
        "integrated_loudness_lufs": r"Integrated loudness:\s+I:\s+(-?[\d.]+)",
        "loudness_range_lu": r"Loudness range:\s+LRA:\s+(-?[\d.]+)",
        "true_peak_dbtp": r"True peak:\s+Peak:\s+(-?[\d.]+)",
    }
    values = {}
    for name, pattern in patterns.items():
        match = re.search(pattern, summary)
        if not match:
            raise RuntimeError(f"FFmpeg did not report {name} for {path}: {summary}")
        values[name] = float(match.group(1))
    values["measurement"] = "FFmpeg ebur128=peak=true; summary precision 0.1 dB"
    return values


def score_sections(path):
    if path.exists():
        with path.open(newline="") as stream:
            rows = list(csv.DictReader(stream))
        boundaries = []
        for row in rows:
            label = row.get("section", row.get("section_name", ""))
            time_text = next((row[key] for key in
                              ("seconds", "time_s", "start_seconds", "time_seconds", "time")
                              if row.get(key)), None)
            if label and time_text is not None:
                seconds = float(time_text)
                if not boundaries or boundaries[-1][1] != label:
                    boundaries.append((max(0.0, seconds), label))
        if boundaries:
            boundaries[0] = (0.0, boundaries[0][1])
            return [(label, start, boundaries[index + 1][0]
                     if index + 1 < len(boundaries) else DURATION)
                    for index, (start, label) in enumerate(boundaries)
                    if start < DURATION]
    return [(label, float(start), min(float(stop), DURATION))
            for label, start, stop in SUITE_SECTIONS if start < DURATION]


def failures_for(name, metrics):
    failures = []
    if metrics["nonfinite_sample_count"]:
        failures.append(f"{name}: nonfinite samples")
    if metrics["clipped_sample_count"]:
        failures.append(f"{name}: clipped samples")
    if metrics["channels"] != 2:
        failures.append(f"{name}: expected two channels")
    if abs(metrics["duration_seconds"] - DURATION) > 1.0 / metrics["sample_rate_hz"]:
        failures.append(f"{name}: expected exactly 60 seconds")
    if metrics["rms_dbfs"] is None or metrics["rms_dbfs"] < -100.0:
        failures.append(f"{name}: silent or effectively silent audio")
    return failures


def produce(args, report):
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=True)
    ffmpeg = shutil.which(args.ffmpeg)
    if not ffmpeg:
        raise RuntimeError(f"FFmpeg executable not found: {args.ffmpeg}")
    dry_path = output / "taikor-suite-dry.wav"
    if not args.skip_render:
        if args.renderer is None:
            raise ValueError("Provide --renderer, or use --skip-render with an existing dry render")
        renderer = args.renderer.resolve()
        command = [str(renderer), "--rate", str(MASTER_RATE), str(output)]
        print("Rendering physical-model performance at 384 kHz…", flush=True)
        rendered = run(command)
        report["renderer"] = {"path": str(renderer), "sha256": sha256(renderer),
                              "command": command, "stdout": rendered.stdout,
                              "stderr": rendered.stderr}
    else:
        report["renderer"] = {"skipped": True, "path": str(args.renderer.resolve())
                              if args.renderer else None}
        if args.renderer and args.renderer.is_file():
            report["renderer"]["sha256"] = sha256(args.renderer)
    dry_format = wav_metadata(dry_path)
    if dry_format["encoding"] != "IEEE float" or dry_format["bits_per_sample"] != 32:
        raise ValueError("Dry render must be unquantised IEEE float32 WAV")
    rate, dry = read_audio(dry_path, mmap=True)
    if rate != MASTER_RATE:
        raise ValueError(f"Dry render must be native {MASTER_RATE} Hz, got {rate}")
    source_metrics = audio_metrics(dry, rate)
    report["source"] = {"file": str(dry_path), "sha256": sha256(dry_path),
                        "format": dry_format,
                        **source_metrics}
    failures = failures_for("dry render", source_metrics)
    if failures:
        raise ValueError("; ".join(failures))

    ir_path = args.impulse_response.resolve()
    ir_rate, ir = read_audio(ir_path)
    if ir.shape[1] == 1:
        ir = np.repeat(ir, 2, axis=1)
    if ir.shape[1] != 2 or not len(ir) or not np.isfinite(ir).all():
        raise ValueError("The room impulse must contain finite mono or stereo audio")
    original_onset = int(np.argmax(np.max(np.abs(ir), axis=1)
                                   > 0.01 * np.max(np.abs(ir)))) / ir_rate
    ir = resample(ir, ir_rate, LISTENING_RATE)
    ir -= np.mean(ir, axis=0, dtype=np.float64).astype(np.float32)
    high_pass = signal.butter(2, 80.0, "highpass", fs=LISTENING_RATE,
                              output="sos").astype(np.float32)
    ir = signal.sosfilt(high_pass, ir, axis=0).astype(np.float32, copy=False)
    ir_energy = math.sqrt(float(np.einsum("ij,ij->", ir, ir, dtype=np.float64)) / 2.0)
    if ir_energy < 1.0e-10:
        raise ValueError("The room impulse has no usable energy")
    ir /= ir_energy  # One common factor preserves the IR's channel balance.
    predelay_seconds = max(0.0, 0.020 - original_onset)
    delay_frames = round(predelay_seconds * LISTENING_RATE)
    report["processing"] = {
        "engine_sample_rate_hz": MASTER_RATE,
        "master_sample_rate_hz": MASTER_RATE,
        "room_processing_sample_rate_hz": LISTENING_RATE,
        "room_ir": str(ir_path), "room_ir_sha256": sha256(ir_path),
        "room_ir_original_sample_rate_hz": ir_rate,
        "room_ir_original_onset_seconds": original_onset,
        "room_added_predelay_seconds": delay_frames / LISTENING_RATE,
        "room_highpass_hz": 80.0, "room_highpass_order": 2,
        "room_ir_energy_normalization_factor": ir_energy,
        "wet_to_dry_rms_ratio_before_master_fades": args.room_ratio,
        "dry_gain_before_mastering": 1.0,
        "resampling": "SciPy polyphase FIR; 64 taps per branch; Kaiser beta 12; float32",
        "entrance_fade_seconds": 0.005, "tail_fade_seconds": 0.300,
        "fade_shape": "half cosine; first and last master frames are zero",
        "mastering": "384 kHz float master: one constant linear gain; no compression, limiting, or loudness normalization",
        "target_true_peak_dbtp": TARGET_TRUE_PEAK_DB,
        "true_peak_estimation": "8x FIR interpolation of 96 kHz delivery with overlapping blocks",
        "wav_delivery": "96 kHz stereo PCM 24-bit; triangular dither at 24-bit precision",
        "mp3_delivery": "48 kHz stereo; libmp3lame 320 kbps CBR; compression_level 0",
    }

    with tempfile.TemporaryDirectory(prefix=".taikor-master-", dir=output) as temporary:
        temp = Path(temporary)
        frames96 = round(DURATION * LISTENING_RATE)
        wet = np.lib.format.open_memmap(temp / "wet.npy", mode="w+",
                                       dtype=np.float32, shape=(frames96, 2))
        dry_energy = wet_energy = 0.0
        print("Adding the theater room and preserving the performance dynamics…", flush=True)
        for channel in range(2):
            dry96 = resample(dry[:, channel], MASTER_RATE, LISTENING_RATE)
            dry_energy += float(np.einsum("i,i->", dry96, dry96, dtype=np.float64))
            convolution = signal.fftconvolve(dry96, ir[:, channel], mode="full")
            wet[:, channel] = 0.0
            usable = min(len(convolution), frames96 - delay_frames)
            wet[delay_frames:delay_frames + usable, channel] = convolution[:usable]
            wet_energy += float(np.einsum("i,i->", wet[:, channel], wet[:, channel],
                                         dtype=np.float64))
            del dry96, convolution
        wet_gain = args.room_ratio * math.sqrt(dry_energy / max(wet_energy, 1.0e-30))
        report["processing"]["wet_gain_after_ir_normalization"] = wet_gain
        master = np.lib.format.open_memmap(temp / "master.npy", mode="w+",
                                          dtype=np.float32, shape=dry.shape)
        for channel in range(2):
            interpolated = resample(wet[:, channel], LISTENING_RATE, MASTER_RATE)
            interpolated *= np.float32(wet_gain)
            np.add(dry[:, channel], interpolated[:len(dry)], out=master[:, channel])
            del interpolated
        del wet
        entrance_frames, tail_frames = round(0.005 * MASTER_RATE), round(0.300 * MASTER_RATE)
        entrance = (0.5 - 0.5 * np.cos(np.linspace(0.0, math.pi, entrance_frames))).astype(np.float32)
        tail = (0.5 + 0.5 * np.cos(np.linspace(0.0, math.pi, tail_frames))).astype(np.float32)
        master[:entrance_frames] *= entrance[:, None]
        master[-tail_frames:] *= tail[:, None]
        delivery = np.empty((frames96, 2), dtype=np.float32)
        for channel in range(2):
            delivery[:, channel] = resample(master[:, channel], MASTER_RATE, LISTENING_RATE)
        before_peak = max(audio_metrics(master, MASTER_RATE)["sample_peak"],
                          estimated_true_peak(delivery))
        if before_peak <= 1.0e-10:
            raise ValueError("The mixed master is silent")
        gain = 10.0 ** (TARGET_TRUE_PEAK_DB / 20.0) / before_peak
        report["processing"]["master_linear_gain"] = gain
        report["processing"]["master_gain_db"] = db(gain)
        report["processing"]["pre_master_peak_dbfs_or_dbtp"] = db(before_peak)
        for start in range(0, len(master), CHUNK_FRAMES):
            master[start:start + CHUNK_FRAMES] *= np.float32(gain)
        delivery *= np.float32(gain)
        master_path = output / "taikor-suite-master-384k-float.wav"
        wavfile.write(master_path, MASTER_RATE, master)
        del master
        print("Preparing the louder listening editions with stereo-linked lookahead limiting…", flush=True)
        report["processing"]["listening_mastering"] = limit_listening(delivery, LISTENING_RATE)
        listening_float = temp / "listening-96k-float.wav"
        wavfile.write(listening_float, LISTENING_RATE, delivery)
        listening_path = output / "taikor-suite-96k24.wav"
        run([ffmpeg, "-hide_banner", "-loglevel", "error", "-nostdin", "-y",
             "-i", str(listening_float), "-af",
             "aresample=96000:osf=s32:dither_method=triangular:output_sample_bits=24",
             "-c:a", "pcm_s24le", "-ar", str(LISTENING_RATE),
             "-metadata", "title=Taikor — HIBIKI / Resonance",
             "-metadata", "comment=Physical-model performance; theater room; stereo-linked listening limiter",
             str(listening_path)])
        mp3_float = temp / "listening-48k-float.wav"
        mp3_samples = resample(delivery, LISTENING_RATE, MP3_RATE)
        del delivery
        wavfile.write(mp3_float, MP3_RATE, mp3_samples)
        del mp3_samples
        mp3_path = output / "taikor-suite.mp3"
        decoded_path = temp / "decoded-mp3.wav"
        mp3_gain = 1.0
        # A codec can introduce overshoots. Keep any codec-only attenuation
        # explicit; the lossless master always retains its original target.
        for attempt in range(3):
            command = [ffmpeg, "-hide_banner", "-loglevel", "error", "-nostdin", "-y",
                       "-i", str(mp3_float), "-af", f"volume={mp3_gain:.12g}",
                       "-c:a", "libmp3lame", "-b:a", "320k", "-compression_level", "0",
                       "-ar", str(MP3_RATE), "-metadata", "title=Taikor — HIBIKI / Resonance",
                       "-metadata", "album=Taikor Demonstrations", str(mp3_path)]
            run(command)
            run([ffmpeg, "-hide_banner", "-loglevel", "error", "-nostdin", "-y",
                 "-i", str(mp3_path), "-c:a", "pcm_f32le", str(decoded_path)])
            decoded_rate, decoded = read_audio(decoded_path, mmap=True)
            mp3_peak = estimated_true_peak(decoded)
            if mp3_peak <= 10.0 ** (-0.3 / 20.0):
                break
            if attempt < 2:
                del decoded
                mp3_gain *= 10.0 ** (-0.5 / 20.0) / mp3_peak
        report["processing"]["mp3_only_linear_gain"] = mp3_gain
        report["processing"]["mp3_only_gain_db"] = db(mp3_gain)
        print("Measuring the masters, encoded playback, and section dynamics…", flush=True)
        artifacts = {}
        failures = []
        for path in (master_path, listening_path, mp3_path):
            if path == mp3_path:
                artifact_rate, samples = read_audio(decoded_path, mmap=True)
            else:
                artifact_rate, samples = read_audio(path, mmap=True)
            metrics = audio_metrics(samples, artifact_rate)
            metrics.update({"file": str(path), "bytes": path.stat().st_size,
                            "sha256": sha256(path), "ebur128": ebur128(ffmpeg, path)})
            if path.suffix == ".wav":
                metrics["format"] = wav_metadata(path)
            if metrics["ebur128"]["true_peak_dbtp"] >= 0.0:
                failures.append(f"{path.name}: FFmpeg measured intersample clipping")
            if path != master_path:
                metrics["estimated_true_peak_dbtp"] = db(estimated_true_peak(samples))
                metrics["true_peak_oversampling_factor"] = 8
                if metrics["estimated_true_peak_dbtp"] is not None and metrics["estimated_true_peak_dbtp"] >= 0.0:
                    failures.append(f"{path.name}: estimated intersample clipping")
            metrics["sections"] = []
            for label, start, stop in score_sections(output / "taikor-suite-score.csv"):
                section = samples[round(start * artifact_rate):round(stop * artifact_rate)]
                if len(section):
                    section_metrics = audio_metrics(section, artifact_rate)
                    metrics["sections"].append({"name": label, "start_seconds": start,
                                                "end_seconds": stop,
                                                "rms_dbfs": section_metrics["rms_dbfs"],
                                                "sample_peak_dbfs": section_metrics["sample_peak_dbfs"],
                                                "crest_factor_db": section_metrics["crest_factor_db"]})
            failures.extend(failures_for(path.name, metrics))
            artifacts[path.name] = metrics
            del samples
        report["artifacts"] = artifacts
        report["master_sha256"] = artifacts[master_path.name]["sha256"]
        listening_stats = report["processing"]["listening_mastering"]
        before_loudness = artifacts[master_path.name]["ebur128"]
        after_loudness = artifacts[listening_path.name]["ebur128"]
        listening_stats.update({
            "input_integrated_loudness_lufs": before_loudness["integrated_loudness_lufs"],
            "output_integrated_loudness_lufs": after_loudness["integrated_loudness_lufs"],
            "measured_loudness_increase_lu": after_loudness["integrated_loudness_lufs"]
                - before_loudness["integrated_loudness_lufs"],
            "input_loudness_range_lu": before_loudness["loudness_range_lu"],
            "output_loudness_range_lu": after_loudness["loudness_range_lu"],
        })
        report["failures"] = failures
        report["status"] = "passed" if not failures else "failed"
        report["qa_scope"] = (
            "Automated format, duration, finite-sample, silence, clipping, loudness and stereo measurements. "
            "These measurements do not certify musical or subjective listening quality. "
            "MP3 duration and sample metrics are measured after gapless FFmpeg decoding."
        )
        del decoded
    for name in ("taikor-suite.mid", "taikor-suite-score.csv"):
        path = output / name
        if path.exists():
            report.setdefault("performance_files", {})[name] = {
                "file": str(path), "sha256": sha256(path), "bytes": path.stat().st_size}
    report["software"]["ffmpeg"] = run([ffmpeg, "-version"]).stdout.splitlines()[0]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--renderer", type=Path, help="Path to the TaikorRenderSuite executable")
    parser.add_argument("--output", type=Path, default=Path("build-demo-suite"))
    parser.add_argument("--skip-render", action="store_true", help="Master an existing 384 kHz dry render")
    parser.add_argument("--room-ratio", type=float, default=0.15,
                        help="Wet/dry RMS ratio before fades and mastering (default: 0.15)")
    parser.add_argument("--impulse-response", type=Path,
                        default=ROOT / "Assets/ImpulseResponses/theater.wav")
    parser.add_argument("--ffmpeg", default="ffmpeg")
    args = parser.parse_args()
    if not 0.0 <= args.room_ratio <= 0.5:
        parser.error("--room-ratio must be between 0 and 0.5")
    if not args.skip_render and args.renderer is None:
        parser.error("--renderer is required unless --skip-render is used")
    args.output.mkdir(parents=True, exist_ok=True)
    report = {
        "schema_version": 1, "created_utc": datetime.now(timezone.utc).isoformat(),
        "status": "failed", "software": {"python": sys.version.split()[0],
                                           "numpy": np.__version__, "scipy": scipy.__version__},
        "production_script": {"path": str(Path(__file__).resolve()), "sha256": sha256(__file__)},
        "source_files": {str(path.relative_to(ROOT)): sha256(path)
                         for path in sorted((ROOT / "Source/DSP").glob("*")) if path.is_file()},
    }
    for path in (ROOT / "Tools/RenderSuite.cpp", ROOT / "Tools/RenderTaikoSuite.cpp"):
        if path.exists():
            report["source_files"][str(path.relative_to(ROOT))] = sha256(path)
    try:
        produce(args, report)
    except (OSError, ValueError, RuntimeError) as error:
        report["status"] = "failed"
        report.setdefault("failures", []).append(str(error))
        print(f"Production failed: {error}", file=sys.stderr)
    report_path = args.output.resolve() / "taikor-suite-qa.json"
    report_path.write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")
    print(f"QA {report['status']}: {report_path}", flush=True)
    return 0 if report["status"] == "passed" else 1


if __name__ == "__main__":
    raise SystemExit(main())
