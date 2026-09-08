#!/usr/bin/env python3
"""Run the frozen 16-cell Phase 1 Stage 3 bridge experiment."""

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
from typing import Any

from prepare_phase1_stage3_wave3 import (
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
    VARIANT_TAG,
    verify_files,
)
from run_phase1_first_pilot import _atomic_json, _directory_size, _sha256
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
            "Create layout.local.json using scripts/hpc/README.md before preflight."
        )
    scratch = Path(json.loads(layout.read_text(encoding="utf-8"))["scratch"])
    if not scratch.is_absolute():
        raise RuntimeError("the machine-local scratch path must be absolute")
    manifest = work / "p01-neutral-feedback/manifests" / f"{EXPERIMENT_ID}-runs.tsv"
    with manifest.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    return work, scratch, rows


def _verify_manifest(rows: list[dict[str, str]], *, work: Path) -> None:
    if len(rows) != EXPECTED_RUNS:
        raise RuntimeError(f"expected {EXPECTED_RUNS} Wave 3 runs, found {len(rows)}")
    if len({row["run_id"] for row in rows}) != EXPECTED_RUNS:
        raise RuntimeError("Wave 3 run IDs are not unique")
    if len({row["scratch_relative_path"] for row in rows}) != EXPECTED_RUNS:
        raise RuntimeError("Wave 3 scratch paths are not unique")
    expected_seeds = dict(SEED_BLOCKS)
    for row in rows:
        if row["experiment_id"] != EXPERIMENT_ID:
            raise RuntimeError(f"unexpected experiment ID in {row['run_id']}")
        if row["variant_id"] != VARIANT_TAG:
            raise RuntimeError(f"unexpected variant ID in {row['run_id']}")
        path = work / row["config_path"]
        if not path.is_file() or _sha256(path) != row["config_sha256"]:
            raise RuntimeError(f"configuration checksum differs: {path}")
        if int(row["master_seed"]) != expected_seeds[row["seed_block_id"]]:
            raise RuntimeError(f"seed differs in {row['run_id']}")
        config = load_config(path)
        if config.host.host_generations != HOST_GENERATIONS:
            raise RuntimeError(f"passage endpoint differs in {path}")
        if config.migration.mode != "fixed_regional_pool":
            raise RuntimeError(f"regional migration mode differs in {path}")
        if config.migration.regional_counts != config.environment.initial_counts:
            raise RuntimeError(f"regional composition differs from focal start: {path}")
        if config.evolution.mutation_probability != 0:
            raise RuntimeError(f"mutation is active in {path}")
        if (
            config.evolution.within_host_selection
            or config.evolution.free_living_selection
        ):
            raise RuntimeError(f"selection is active in {path}")
        if config.output.environment_counts_mode != "all":
            raise RuntimeError(f"environmental trajectories are incomplete in {path}")


def completion_issues(
    rows: list[dict[str, str]], *, work: Path, scratch: Path
) -> list[str]:
    """Audit committed Wave 3 outputs without modifying them."""
    issues: list[str] = []
    for row in rows:
        label = row["run_id"]
        output = scratch / row["scratch_relative_path"]
        config_path = work / row["config_path"]
        try:
            if _sha256(config_path) != row["config_sha256"]:
                raise ValueError("raw configuration hash differs from manifest")
            config = load_config(config_path)
            expected = json.loads(json.dumps(config.to_dict()))
            resolved = json.loads(
                (output / "resolved_config.json").read_text(encoding="utf-8")
            )
            completed = json.loads(
                (output / "completion.json").read_text(encoding="utf-8")
            )
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
                if (output / name).stat().st_size != int(size):
                    raise ValueError(f"committed size differs for {name}")
            final_name = "final_environment_rep000.npz"
            if completed.get("final_environment_sha256", {}).get(final_name) != _sha256(
                output / final_name
            ):
                raise ValueError("final environmental state checksum differs")
            execution = json.loads(
                (output / "execution-summary.json").read_text(encoding="utf-8")
            )
            if execution.get("run_id") != label or execution.get("status") not in {
                "complete",
                "already-complete",
            }:
                raise ValueError("execution record is not complete for this run")
            provenance = json.loads(
                (output / "provenance.json").read_text(encoding="utf-8")
            )
            if not completed.get("source_sha256") or completed[
                "source_sha256"
            ] != provenance.get("source_sha256"):
                raise ValueError("completion and provenance source hashes differ")
        except (OSError, ValueError, KeyError, TypeError) as exc:
            issues.append(f"{label}: {exc}")
    return issues


def smoke_assessment(
    rows: list[dict[str, str]], *, work: Path, scratch: Path
) -> dict[str, Any]:
    selected = [
        row
        for row in rows
        if row["cell"] in SMOKE_CELLS and row["seed_block_id"] == "sb0001"
    ]
    issues = completion_issues(selected, work=work, scratch=scratch)
    if len(selected) != 3:
        issues.append("the smoke subset must contain exactly three populations")
    if issues:
        return {"passed": False, "issues": issues}

    measurements: dict[int, dict[str, float]] = {}
    cells_by_id = {cell.short_id: cell for cell in CELLS}
    for row in selected:
        output = scratch / row["scratch_relative_path"]
        record = json.loads(
            (output / "execution-summary.json").read_text(encoding="utf-8")
        )
        if record.get("resumed"):
            issues.append(
                f"{row['cell']}: resumed smoke run needs manual resource review"
            )
        hours = float(record["elapsed_seconds"]) / 3600
        if hours <= 0 or hours > 48:
            issues.append(
                f"{row['cell']}: runtime is outside the 0-48 hour safety range"
            )
        measurements[cells_by_id[row["cell"]].hosts] = {
            "hours": hours,
            "gib": _directory_size(output) / 1024**3,
        }

    projected_gib = 0.0
    projected_population_hours = 0.0
    maximum_hours = 0.0
    for cell in CELLS:
        measurement = measurements[cell.hosts]
        # These bridge cells span extreme feedback and migration values. A 2x
        # margin is deliberately retained because resource use may be nonlinear.
        predicted_hours = 2 * measurement["hours"]
        predicted_gib = 2 * measurement["gib"]
        projected_population_hours += predicted_hours * len(SEED_BLOCKS)
        projected_gib += predicted_gib * len(SEED_BLOCKS)
        maximum_hours = max(maximum_hours, predicted_hours)
    if projected_gib > 350:
        issues.append("projected Wave 3 storage exceeds 350 GiB (70% of 500 GiB)")
    if maximum_hours > 48:
        issues.append("projected runtime exceeds 48 hours for at least one population")
    available = shutil.disk_usage(scratch).free / 1024**3
    if available < projected_gib:
        issues.append("scratch has less free space than the conservative projection")
    return {
        "passed": not issues,
        "issues": issues,
        "measurements_by_H": measurements,
        "projected_output_gib": projected_gib,
        "projected_population_hours": projected_population_hours,
        "maximum_projected_population_hours": maximum_hours,
        "free_scratch_gib": available,
        "method": "one smoke population per H stratum, multiplied by a 2x margin",
    }


def _compile_community(repository: Path, python: str, *, mode: str) -> int:
    command = [
        os.path.abspath(python),
        str(repository / "scripts/compile_phase1_stage3_dbrda_inputs.py"),
        "--repository",
        str(repository),
        "--through-wave",
        "3",
    ]
    if mode == "endpoint":
        command.append("--endpoint-only")
    elif mode == "prc":
        command.append("--prc-only")
    return subprocess.run(command, check=False).returncode


def main() -> int:
    parser = argparse.ArgumentParser(
        description=__doc__,
        epilog=(
            "Use --prepare-only, --smoke-only and --check-smoke before the full "
            "launch. After completion, --community-only builds the updated "
            "Wave 1+2+3 X, Y, Y-prime and PRC tables."
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
        help="simultaneous populations; each uses two host workers",
    )
    parser.add_argument("--monitor-interval", type=float, default=2)
    parser.add_argument("--cell", action="append", help="restrict to c0091-style cell")
    parser.add_argument(
        "--seed-block", action="append", help="restrict to sb0001-style seed"
    )
    parser.add_argument(
        "--smoke-only",
        action="store_true",
        help="run c0096, c0098 and c0100 with sb0001 only",
    )
    parser.add_argument(
        "--allow-dirty-source",
        action="store_true",
        help="allow an explicitly reviewed uncommitted model source",
    )
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument(
        "--prepare-only",
        action="store_true",
        help="verify inputs and create isolated scratch folders; do not simulate",
    )
    mode.add_argument(
        "--dry-run",
        action="store_true",
        help="list selected populations and destinations without writing",
    )
    mode.add_argument(
        "--check-smoke",
        action="store_true",
        help="audit the three safety runs and estimate Wave 3 resources",
    )
    mode.add_argument(
        "--community-only",
        action="store_true",
        help="compile updated Wave 1+2+3 endpoint and PRC community tables",
    )
    mode.add_argument(
        "--endpoint-only",
        action="store_true",
        help="compile updated Wave 1+2+3 passage-100 X/Y/TV tables only",
    )
    mode.add_argument(
        "--prc-only",
        action="store_true",
        help="compile the updated Wave 1+2+3 passage 0-100 PRC table only",
    )
    args = parser.parse_args()
    if args.jobs < 1 or args.monitor_interval <= 0:
        parser.error("jobs and monitor interval must be positive")
    if args.smoke_only and (args.community_only or args.endpoint_only or args.prc_only):
        parser.error("community-table modes cannot be combined with --smoke-only")

    repository = args.repository.resolve()
    if args.community_only:
        return _compile_community(repository, args.python, mode="community")
    if args.endpoint_only:
        return _compile_community(repository, args.python, mode="endpoint")
    if args.prc_only:
        return _compile_community(repository, args.python, mode="prc")

    differences = verify_files(repository)
    if differences:
        raise SystemExit("Frozen Wave 3 inputs differ:\n" + "\n".join(differences))
    work, scratch, rows = load_rows(repository)
    _verify_runtime(Path(os.path.abspath(args.python)))
    _verify_manifest(rows, work=work)
    selected = [
        row
        for row in rows
        if (not args.cell or row["cell"] in set(args.cell))
        and (not args.seed_block or row["seed_block_id"] in set(args.seed_block))
        and (
            not args.smoke_only
            or (row["cell"] in SMOKE_CELLS and row["seed_block_id"] == "sb0001")
        )
    ]
    if not selected:
        raise SystemExit("The requested Wave 3 selection contains no populations")
    print(
        f"Preflight passed: {len(selected)} populations; {HOST_GENERATIONS} passages; "
        f"16 bridge cells; up to {args.jobs} populations and {2 * args.jobs} "
        "host workers concurrently.",
        flush=True,
    )
    if args.dry_run:
        for row in selected:
            print(f"{row['run_id']} -> {scratch / row['scratch_relative_path']}")
        return 0
    if args.check_smoke or (
        not args.smoke_only
        and not args.prepare_only
        and not args.cell
        and not args.seed_block
    ):
        assessment = smoke_assessment(rows, work=work, scratch=scratch)
        print(json.dumps(assessment, indent=2), flush=True)
        if not assessment["passed"]:
            print("Full Wave 3 not launched. Complete/review the smoke runs first.")
            return 1
        if args.check_smoke:
            return 0

    existing = [
        row
        for row in selected
        if (scratch / row["scratch_relative_path"] / "completion.json").is_file()
    ]
    issues = completion_issues(existing, work=work, scratch=scratch)
    if issues:
        raise SystemExit("Existing completions failed the audit:\n" + "\n".join(issues))
    _prepare_scratch(selected, work=work, scratch=scratch)
    if args.prepare_only:
        print("Prepared isolated Wave 3 directories; no simulations launched.")
        return 0
    if len(existing) != len(selected):
        _require_frozen_source(repository, args.allow_dirty_source)

    stop = threading.Event()
    previous_handlers = {
        sig: signal.getsignal(sig) for sig in (signal.SIGINT, signal.SIGTERM)
    }

    def request_stop(signum: int, _frame: Any) -> None:
        print(
            f"{signal.Signals(signum).name}: stopping; checkpoints are retained.",
            flush=True,
        )
        stop.set()

    results: list[dict[str, Any]] = []
    try:
        for sig in previous_handlers:
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
        for sig, handler in previous_handlers.items():
            signal.signal(sig, handler)

    launcher = scratch / "p01-neutral-feedback" / STAGE_DIRECTORY / "_launcher"
    launcher.mkdir(parents=True, exist_ok=True)
    stamp = datetime.now(UTC).strftime("%Y%m%dT%H%M%S%fZ")
    _atomic_json(
        launcher / f"summary-{stamp}.json",
        {"experiment_id": EXPERIMENT_ID, "runs": results},
    )
    if len(results) != len(selected) or any(
        result["status"] not in {"complete", "already-complete"} for result in results
    ):
        print("Some runs are incomplete. Review run.err, then repeat to resume.")
        return 1
    issues = completion_issues(selected, work=work, scratch=scratch)
    if issues:
        raise SystemExit("Wave 3 completion audit failed:\n" + "\n".join(issues))
    if args.smoke_only:
        print("Three smoke populations complete; run --check-smoke next.")
    else:
        print(
            "Wave 3 is complete. Run --community-only to compile the updated "
            "db-RDA, Hellinger-RDA and PRC input tables."
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
