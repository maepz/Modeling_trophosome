#!/usr/bin/env python3
"""Run the historical slide-13 and matched-neutral scalability regimes."""

from __future__ import annotations

import argparse
import os
import shlex
import subprocess
import sys
from pathlib import Path

PROFILE_NAMES = (
    "historical-slide13",
    "matched-neutral-capacity",
    "matched-neutral-mutation",
)


def _profile_arguments(profile: str) -> list[str]:
    common = [
        "--baseline-hosts",
        "100",
        "--baseline-workers",
        "10",
        "--bottleneck",
        "10",
        "--growth-factor",
        "1.2",
        "--escape-fraction",
        "0.01",
        "--initial-strains",
        "100",
    ]
    if profile == "historical-slide13":
        return common + [
            "--models",
            "v1_3,v3_1",
            "--dimensions",
            "within_host_population_size",
            "--carrying-capacities",
            "10000,100000,1000000,10000000,100000000",
            "--baseline-mutation-rate",
            "1e-12",
            "--steady-generations",
            "0",
            "--historical-v1-3-steady-generations",
            "50",
            "--v3-steady-generations",
            "0",
            "--host-generations",
            "20",
            "--regime-label",
            "slide 13 historical replication",
        ]
    if profile == "matched-neutral-capacity":
        return common + [
            "--models",
            "legacy,v3_1,exact",
            "--dimensions",
            "within_host_population_size",
            "--carrying-capacities",
            "10000,100000,1000000,10000000,100000000",
            "--baseline-mutation-rate",
            "1e-12",
            "--steady-generations",
            "0",
            "--host-generations",
            "1",
            "--regime-label",
            "matched neutral low mutation",
        ]
    if profile == "matched-neutral-mutation":
        return common + [
            "--models",
            "legacy,v3_1,exact",
            "--dimensions",
            "mutation_rate",
            "--mutation-rates",
            "1e-12,1e-9,1e-6,1e-4,1e-3",
            "--baseline-carrying-capacity",
            "100000",
            "--steady-generations",
            "0",
            "--host-generations",
            "1",
            "--regime-label",
            "matched neutral mutation supply",
        ]
    raise ValueError(f"unknown profile: {profile}")


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Run a faithful slide-13 benchmark and two matched-neutral controls. "
            "The historical profile is intentionally expensive."
        )
    )
    parser.add_argument(
        "--profiles",
        default=",".join(PROFILE_NAMES),
        help=f"comma-separated subset of: {', '.join(PROFILE_NAMES)}",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("results/benchmarks/regime-comparison"),
    )
    parser.add_argument("--repeats", type=int, default=3)
    parser.add_argument("--historical-repeats", type=int, default=1)
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args()

    profiles = tuple(item.strip() for item in args.profiles.split(",") if item)
    unknown = set(profiles) - set(PROFILE_NAMES)
    if unknown:
        parser.error(f"unknown profiles: {', '.join(sorted(unknown))}")
    if args.repeats < 1 or args.historical_repeats < 1:
        parser.error("repeat counts must be positive")

    repository = Path(__file__).resolve().parents[1]
    benchmark = repository / "scripts" / "benchmark_v1_3_vs_exact.py"
    output_dir = args.output_dir
    if not output_dir.is_absolute():
        output_dir = repository / output_dir

    environment = os.environ.copy()
    source_path = str(repository / "src")
    existing_pythonpath = environment.get("PYTHONPATH")
    environment["PYTHONPATH"] = (
        f"{source_path}{os.pathsep}{existing_pythonpath}"
        if existing_pythonpath
        else source_path
    )

    for profile in profiles:
        repeats = (
            args.historical_repeats if profile == "historical-slide13" else args.repeats
        )
        output = output_dir / f"{profile}.csv"
        command = [
            sys.executable,
            str(benchmark),
            "--output",
            str(output),
            "--repeats",
            str(repeats),
            *_profile_arguments(profile),
        ]
        print(shlex.join(command), flush=True)
        if args.dry_run:
            continue
        output_dir.mkdir(parents=True, exist_ok=True)
        subprocess.run(command, cwd=repository, env=environment, check=True)


if __name__ == "__main__":
    main()
