#!/usr/bin/env python3
"""Reproduce expression profiles, simulations, and temporal-to-spatial landmarks."""
from pathlib import Path
import argparse
import shutil
import subprocess
import sys

ROOT = Path(__file__).resolve().parent
MODEL = ROOT / "analysis/intron_exon"
MODEL_STAGES = ("verify", "empirical", "sweep", "fits", "numerics", "assessment",
                "figures", "manuscript-figures", "movie", "all")


def model_run(stage, output, fresh=False):
    """Run the original model in a disposable copy, preserving recorded results."""
    if stage == "verify":
        working = MODEL
    else:
        from seqimaging.data import output_directory
        working = output_directory(output) / "model"
        if not working.exists():
            shutil.copytree(MODEL, working,
                            ignore=shutil.ignore_patterns("__pycache__", "derived", "figures"))
    command = [sys.executable, str(working / "reproduce_analysis.py"), "--stage", stage]
    if fresh:
        command.append("--fresh")
    subprocess.run(command, cwd=working, check=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    sub.add_parser("verify", help="Check source data, figure assignments, and saved model results")
    profiles = sub.add_parser("profiles", help="Recreate the spatial and temporal numerical panels")
    profiles.add_argument("--output", type=Path, default=ROOT / "outputs")
    model = sub.add_parser("model", help="Run an intron/exon model stage in the output directory")
    model.add_argument("--stage", choices=MODEL_STAGES, default="verify")
    model.add_argument("--output", type=Path, default=ROOT / "outputs")
    model.add_argument("--fresh", action="store_true", help="Fit without archived warm starts")
    peaks = sub.add_parser("peaks", help="Rebuild the second-revision peak-order analysis and Fig. 7")
    peaks.add_argument("--output", type=Path, default=ROOT / "outputs/peak_order")
    peaks.add_argument("--bootstraps", type=int, default=20000)
    peaks.add_argument("--skip-figures", action="store_true")
    peaks.add_argument("--apply-confirmed-corrections", action="store_true",
                       help="Separate sensitivity run with the recovered late KR-cohort eve background")
    landmarks = sub.add_parser("landmarks", help="Rebuild the current landmark analysis and optionally Fig. 7")
    landmarks.add_argument("--output", type=Path, default=ROOT / "outputs/landmarks")
    landmarks.add_argument("--figure", action="store_true", help="Also rebuild the six-panel Fig. 7; requires Inkscape")
    args = parser.parse_args()
    if args.command == "verify":
        from seqimaging.validation import verify
        verify()
        sys.stdout.flush()
        model_run("verify", ROOT / "outputs")
    elif args.command == "profiles":
        from seqimaging.figures import reproduce_profiles
        reproduce_profiles(args.output)
    elif args.command == "landmarks":
        working = args.output.resolve()
        source = ROOT / "analysis/landmarks"
        if working == source or source in working.parents:
            parser.error("Choose an output directory outside analysis/landmarks")
        if working.exists():
            parser.error("Landmark output directory already exists; choose a fresh --output path")
        shutil.copytree(source, working,
                        ignore=shutil.ignore_patterns("__pycache__", "outputs", "*.tmp", "*.writing"))
        command = [sys.executable, str(working / "reproduce.py")]
        if args.figure:
            command.append("--figure")
        subprocess.run(command, cwd=working, check=True)
    elif args.command == "peaks":
        command = [sys.executable, str(ROOT / "analysis/peak_order/run_analysis.py"),
                   "--output-dir", str(args.output), "--bootstraps", str(args.bootstraps)]
        if args.skip_figures:
            command.append("--skip-figures")
        if args.apply_confirmed_corrections:
            command.append("--apply-confirmed-corrections")
        subprocess.run(command, cwd=ROOT, check=True)
    else:
        model_run(args.stage, args.output, args.fresh)


if __name__ == "__main__":
    main()
