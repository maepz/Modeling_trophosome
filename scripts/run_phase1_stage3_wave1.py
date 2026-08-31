#!/usr/bin/env python3
"""Run the frozen Phase 1 Stage 3 first wave, starting with three safety runs."""

from __future__ import annotations

import argparse
import csv
import json
import os
import shutil
import signal
import subprocess
import sys
import threading
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import UTC, datetime
from pathlib import Path

from prepare_phase1_stage3_wave1 import (
    CELLS,
    EXPECTED_RUNS,
    EXPERIMENT_ID,
    HOST_GENERATIONS,
    MODEL_SPEC_VERSION,
    OUTPUT_SCHEMA_VERSION,
    SEED_BLOCKS,
    SMOKE_CELLS,
    SOFTWARE_VERSION,
    STAGE_DIRECTORY,
    verify_files,
)
from run_phase1_first_pilot import _atomic_json, _sha256
from run_phase1_first_pilot_v2_1 import (
    _prepare_scratch,
    _require_frozen_source,
    _run_one,
    _verify_runtime,
)
from run_phase1_second_pilot import _resolved_config_sha256

from trophosome.config import load_config
from trophosome.simulation import _output_fields


def load_rows(repository: Path) -> tuple[Path, Path, list[dict[str, str]]]:
    work = repository / "experiments/work/trophosome"
    layout = work / "layout.local.json"
    if not layout.is_file():
        raise RuntimeError(
            "Create layout.local.json using scripts/hpc/README.md first."
        )
    scratch = Path(json.loads(layout.read_text())["scratch"])
    if not scratch.is_absolute():
        raise RuntimeError("the machine-local scratch path must be absolute")
    with (
        work / "p01-neutral-feedback/manifests" / f"{EXPERIMENT_ID}-runs.tsv"
    ).open() as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    return work, scratch, rows


def completion_issues(
    rows: list[dict[str, str]], *, work: Path, scratch: Path
) -> list[str]:
    """Read-only audit, also used before a completed run is skipped."""
    issues = []
    for row in rows:
        label = row["run_id"]
        output = scratch / row["scratch_relative_path"]
        config_path = work / row["config_path"]
        try:
            if _sha256(config_path) != row["config_sha256"]:
                raise ValueError("raw configuration hash differs from manifest")
            config = load_config(config_path)
            expected = json.loads(json.dumps(config.to_dict()))
            resolved = json.loads((output / "resolved_config.json").read_text())
            completed = json.loads((output / "completion.json").read_text())
            if resolved != expected:
                raise ValueError("resolved configuration differs from frozen TOML")
            if completed.get("complete") is not True:
                raise ValueError("completion marker is not committed")
            for field, value in (
                ("config_sha256", _resolved_config_sha256(expected)),
                ("model_spec_version", MODEL_SPEC_VERSION),
                ("output_schema_version", OUTPUT_SCHEMA_VERSION),
                ("software_version", SOFTWARE_VERSION),
            ):
                if completed.get(field) != value:
                    raise ValueError(f"completed {field} differs")
            sizes = completed.get("output_sizes", {})
            if set(sizes) != set(_output_fields(config)):
                raise ValueError("completion table inventory differs")
            for name, size in sizes.items():
                if (output / name).stat().st_size != size:
                    raise ValueError(f"committed size differs for {name}")
            final_name = "final_environment_rep000.npz"
            if completed.get("final_environment_sha256", {}).get(final_name) != _sha256(
                output / final_name
            ):
                raise ValueError("final environmental state checksum differs")
            execution = json.loads((output / "execution-summary.json").read_text())
            if (
                execution.get("run_id") != label
                or execution.get("status") != "complete"
            ):
                raise ValueError("execution record is not complete for this run")
            provenance = json.loads((output / "provenance.json").read_text())
            if not completed.get("source_sha256") or completed[
                "source_sha256"
            ] != provenance.get("source_sha256"):
                raise ValueError("completion and provenance source hashes differ")
        except (OSError, ValueError, KeyError, TypeError) as exc:
            issues.append(f"{label}: {exc}")
    return issues


def smoke_assessment(rows: list[dict[str, str]], *, work: Path, scratch: Path) -> dict:
    selected = [
        r for r in rows if r["cell"] in SMOKE_CELLS and r["seed_block_id"] == "sb0001"
    ]
    issues = completion_issues(selected, work=work, scratch=scratch)
    if len(selected) != 3:
        issues.append("the smoke subset must contain exactly three populations")
    if issues:
        return {"passed": False, "issues": issues}
    measurements = {}
    for row in selected:
        output = scratch / row["scratch_relative_path"]
        record = json.loads((output / "execution-summary.json").read_text())
        # A resumed attempt's elapsed_seconds covers only that attempt and cannot
        # safely estimate total runtime. Review it rather than underestimating.
        if record.get("resumed"):
            issues.append(
                f"{row['cell']}: resumed smoke run needs manual resource review"
            )
        hours = float(record["elapsed_seconds"]) / 3600
        if hours > 48 or hours <= 0:
            issues.append(
                f"{row['cell']}: runtime is outside the 0-48 hour safety range"
            )
        measurements[row["cell"]] = {
            "hours": hours,
            "gib": sum(p.stat().st_size for p in output.rglob("*") if p.is_file())
            / 1024**3,
        }
    projected_gib = projected_hours = maximum_hours = 0.0
    for cell in CELLS:
        if cell.mutation_probability == "0":
            prediction = {
                key: value * cell.hosts / 10000
                for key, value in measurements["c0049"].items()
            }
        else:
            # Both ends of the H range are measured at maximal feedback. Include
            # small-H overhead by interpolating, rather than scaling it to zero.
            weight = (cell.hosts - 100) / (10000 - 100)
            prediction = {
                key: (1 - weight) * measurements["c0034"][key]
                + weight * measurements["c0050"][key]
                for key in ("gib", "hours")
            }
        # Conservative 2x margin; feedback and mutation may scale nonlinearly.
        projected_gib += 2 * prediction["gib"] * len(SEED_BLOCKS)
        hours = 2 * prediction["hours"]
        projected_hours += hours * len(SEED_BLOCKS)
        maximum_hours = max(maximum_hours, hours)
    if projected_gib > 350:
        issues.append("projected full-wave storage exceeds 350 GiB (70% of 500 GiB)")
    if maximum_hours > 48:
        issues.append("projected runtime exceeds 48 hours for at least one population")
    available = shutil.disk_usage(scratch).free / 1024**3
    if available < projected_gib:
        issues.append(
            "scratch has less free space than the conservative wave projection"
        )
    return {
        "passed": not issues,
        "issues": issues,
        "measurements": measurements,
        "projected_output_gib": projected_gib,
        "projected_population_hours": projected_hours,
        "maximum_projected_population_hours": maximum_hours,
        "free_scratch_gib": available,
        "method": (
            "maximal-feedback smoke runs: scale mutation-free H, interpolate "
            "mutation-enabled H between 100 and 10000; 2x safety factor; "
            "not a guarantee"
        ),
    }


def build_report(repository: Path, python: str) -> int:
    return subprocess.run(
        [
            python,
            str(repository / "scripts/build_phase1_stage3_wave1_report.py"),
            "--repository",
            str(repository),
        ],
        cwd=repository,
        check=False,
    ).returncode


def verify_source_location(python: str, repository: Path) -> None:
    result = subprocess.run(
        [python, "-c", "import trophosome; print(trophosome.__file__)"],
        capture_output=True,
        text=True,
        check=True,
    )
    imported = Path(result.stdout.strip()).resolve().parent
    if imported != (repository / "src/trophosome").resolve():
        raise RuntimeError(
            "The selected Python is not importing this frozen repository. "
            "Activate trophosome and run: python -m pip install -e '.[report]'"
        )


def main() -> int:
    parser = argparse.ArgumentParser(
        description=__doc__,
        epilog=(
            "First use --prepare-only, then --smoke-only (3 included populations), "
            "then --check-smoke. Run without a mode to finish all 288 new runs. "
            "--report-only never launches simulations."
        ),
    )
    parser.add_argument(
        "--repository", type=Path, default=Path(__file__).resolve().parents[1]
    )
    parser.add_argument("--python", default=sys.executable)
    parser.add_argument(
        "--jobs",
        type=int,
        default=8,
        help="simultaneous populations; two host workers each",
    )
    parser.add_argument("--monitor-interval", type=float, default=2)
    parser.add_argument(
        "--smoke-only",
        action="store_true",
        help="run c0034, c0049, c0050 with sb0001 only",
    )
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument(
        "--prepare-only",
        action="store_true",
        help="check frozen inputs and create scratch directories only",
    )
    mode.add_argument(
        "--dry-run",
        action="store_true",
        help="print the selected population IDs without writing or simulating",
    )
    mode.add_argument(
        "--report-only",
        action="store_true",
        help="audit all 288 new results and reused references, then rebuild the report",
    )
    mode.add_argument(
        "--check-smoke",
        action="store_true",
        help="assess the three safety runs without launching anything",
    )
    parser.add_argument(
        "--no-report",
        action="store_true",
        help="defer reporting; preserve all simulation outputs",
    )
    args = parser.parse_args()
    if args.jobs < 1 or args.monitor_interval <= 0:
        parser.error("jobs and monitor interval must be positive")
    if args.report_only and (args.smoke_only or args.no_report):
        parser.error(
            "--report-only cannot be combined with --smoke-only or --no-report"
        )
    repository = args.repository.resolve()
    differences = verify_files(repository)
    if differences:
        raise SystemExit("Frozen Stage 3 inputs differ:\n" + "\n".join(differences))
    work, scratch, rows = load_rows(repository)
    _verify_runtime(Path(os.path.abspath(args.python)))
    selected = [
        r
        for r in rows
        if not args.smoke_only
        or (r["cell"] in SMOKE_CELLS and r["seed_block_id"] == "sb0001")
    ]
    print(
        f"Preflight passed: {len(selected)} new populations; "
        f"{HOST_GENERATIONS} passages; m=0.1; "
        f"H <= 10000; {args.jobs} concurrent populations, "
        f"{2 * args.jobs} host workers.",
        flush=True,
    )
    if len(rows) != EXPECTED_RUNS:
        raise SystemExit(
            "The frozen run manifest must contain exactly 288 new populations"
        )
    if args.dry_run:
        for row in selected:
            print(f"{row['run_id']} -> {scratch / row['scratch_relative_path']}")
        return 0
    if args.report_only:
        return build_report(repository, args.python)
    if args.check_smoke or (not args.smoke_only and not args.prepare_only):
        assessment = smoke_assessment(rows, work=work, scratch=scratch)
        print(json.dumps(assessment, indent=2), flush=True)
        if not assessment["passed"]:
            print("Full wave not launched. Complete/review the three smoke runs first.")
            return 1
        if args.check_smoke:
            return 0
    existing = [
        r
        for r in selected
        if (scratch / r["scratch_relative_path"] / "completion.json").exists()
    ]
    issues = completion_issues(existing, work=work, scratch=scratch)
    if issues:
        raise SystemExit("Existing completions failed the audit:\n" + "\n".join(issues))
    if not args.prepare_only:
        verify_source_location(args.python, repository)
        _require_frozen_source(repository, False)
    _prepare_scratch(selected, work=work, scratch=scratch)
    if args.prepare_only:
        print("Prepared isolated Stage 3 directories; no simulations launched.")
        return 0
    stop = threading.Event()
    previous = {sig: signal.getsignal(sig) for sig in (signal.SIGINT, signal.SIGTERM)}

    def request_stop(signum, _frame):
        print(
            f"{signal.Signals(signum).name}: stopping; "
            "latest checkpoints are retained.",
            flush=True,
        )
        stop.set()

    results = []
    try:
        for sig in previous:
            signal.signal(sig, request_stop)
        with ThreadPoolExecutor(max_workers=args.jobs) as pool:
            futures = {
                pool.submit(
                    _run_one,
                    row,
                    repository=repository,
                    work=work,
                    scratch=scratch,
                    python=Path(args.python),
                    monitor_interval=args.monitor_interval,
                    stop_event=stop,
                ): row
                for row in selected
            }
            for future in as_completed(futures):
                try:
                    result = future.result()
                except Exception as exc:
                    result = {
                        "run_id": futures[future]["run_id"],
                        "status": "launcher-error",
                        "error": str(exc),
                    }
                results.append(result)
                if result["status"] not in {"complete", "already-complete"}:
                    stop.set()
    finally:
        for sig, handler in previous.items():
            signal.signal(sig, handler)
    stamp = datetime.now(UTC).strftime("%Y%m%dT%H%M%S%fZ")
    launcher_directory = (
        scratch / "p01-neutral-feedback" / STAGE_DIRECTORY / "_launcher"
    )
    launcher_directory.mkdir(parents=True, exist_ok=True)
    _atomic_json(
        launcher_directory / f"summary-{stamp}.json",
        {"experiment_id": EXPERIMENT_ID, "runs": results},
    )
    if len(results) != len(selected) or any(
        r["status"] not in {"complete", "already-complete"} for r in results
    ):
        print(
            "Some runs are incomplete. Review run.err, "
            "then repeat this command to resume."
        )
        return 1
    if args.smoke_only:
        print(
            "Three smoke runs complete. Use --check-smoke "
            "before launching the full wave."
        )
        return 0
    return 0 if args.no_report else build_report(repository, args.python)


if __name__ == "__main__":
    raise SystemExit(main())
