#!/usr/bin/env python3
"""Run the Phase 1 equilibrium-and-precision second pilot on an HPC node."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import signal
import subprocess
import sys
import threading
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import UTC, datetime
from pathlib import Path
from typing import Any

from prepare_phase1_second_pilot import (
    CELLS,
    DESIGN_STEM,
    EXPERIMENT_ID,
    HOST_GENERATIONS,
    MODEL_SPEC_VERSION,
    OUTPUT_SCHEMA_VERSION,
    SEED_BLOCKS,
    SOFTWARE_VERSION,
    STAGE_DIRECTORY,
    VARIANT_TAG,
)
from run_phase1_first_pilot import _atomic_json, _sha256
from run_phase1_first_pilot_v2_1 import (
    _prepare_scratch,
    _require_frozen_source,
    _run_one,
    _verify_runtime,
)

from trophosome.config import load_config

EXPECTED_CELLS = {cell.cell_id for cell in CELLS}
EXPECTED_SEED_BLOCKS = {seed_block_id for seed_block_id, _ in SEED_BLOCKS}
EXPECTED_RUNS = len(EXPECTED_CELLS) * len(EXPECTED_SEED_BLOCKS)
MANIFEST_NAME = f"{DESIGN_STEM}-runs.tsv"


def _resolved_config_sha256(payload: dict[str, Any]) -> str:
    serialized = json.dumps(payload, sort_keys=True, separators=(",", ":")).encode(
        "utf-8"
    )
    return hashlib.sha256(serialized).hexdigest()


def _normalise_cell_id(value: str) -> str:
    if value.startswith("p01-s02-c"):
        return value
    if value.startswith("c") and value[1:].isdigit():
        return f"p01-s02-c{int(value[1:]):04d}"
    if value.isdigit():
        return f"p01-s02-c{int(value):04d}"
    raise argparse.ArgumentTypeError("cell must look like c0021, 21, or p01-s02-c0021")


def _load_rows(repository: Path) -> tuple[Path, Path, list[dict[str, str]]]:
    work = repository / "experiments/work/trophosome"
    layout_path = work / "layout.local.json"
    if not layout_path.is_file():
        raise RuntimeError(
            f"machine-local layout is missing: {layout_path}. Create it using "
            "the instructions under 'Create the machine-local storage layout' "
            "in scripts/hpc/README.md, then rerun this command."
        )
    layout = json.loads(layout_path.read_text(encoding="utf-8"))
    scratch = Path(layout["scratch"])
    manifest = work / "p01-neutral-feedback/manifests" / MANIFEST_NAME
    with manifest.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    return work, scratch, rows


def _select_rows(
    rows: list[dict[str, str]],
    *,
    cells: set[str] | None,
    seed_blocks: set[str] | None,
) -> list[dict[str, str]]:
    selected_cells = cells or EXPECTED_CELLS
    selected_seeds = seed_blocks or EXPECTED_SEED_BLOCKS
    order = {cell.cell_id: index for index, cell in enumerate(CELLS)}
    selected = [
        row
        for row in rows
        if row["cell_id"] in selected_cells and row["seed_block_id"] in selected_seeds
    ]
    selected.sort(key=lambda row: (order[row["cell_id"]], row["seed_block_id"]))
    return selected


def _verify_manifest(
    rows: list[dict[str, str]],
    selected: list[dict[str, str]],
    *,
    work: Path,
) -> None:
    if len(rows) != EXPECTED_RUNS:
        raise RuntimeError(f"expected {EXPECTED_RUNS} manifest rows, found {len(rows)}")
    if len({row["run_id"] for row in rows}) != EXPECTED_RUNS:
        raise RuntimeError("run IDs are not unique")
    if len({row["scratch_relative_path"] for row in rows}) != EXPECTED_RUNS:
        raise RuntimeError("scratch output paths are not unique")
    if {row["cell_id"] for row in rows} != EXPECTED_CELLS:
        raise RuntimeError("manifest does not contain the six frozen sentinel cells")
    if {row["seed_block_id"] for row in rows} != EXPECTED_SEED_BLOCKS:
        raise RuntimeError(
            f"manifest does not contain the {len(EXPECTED_SEED_BLOCKS)} frozen "
            "seed blocks"
        )
    if not selected:
        raise RuntimeError("the requested cell/seed selection contains no runs")

    for row in selected:
        if row["experiment_id"] != EXPERIMENT_ID:
            raise RuntimeError(f"unexpected experiment ID in {row['run_id']}")
        if row["variant_id"] != VARIANT_TAG:
            raise RuntimeError(f"unexpected variant ID in {row['run_id']}")
        config_path = work / row["config_path"]
        if not config_path.is_file() or _sha256(config_path) != row["config_sha256"]:
            raise RuntimeError(f"configuration checksum differs: {config_path}")
        config = load_config(config_path)
        if config.host.host_generations != HOST_GENERATIONS:
            raise RuntimeError(f"host horizon is not 250 in {config_path}")
        if config.output.environment_counts_mode != "all":
            raise RuntimeError(
                f"complete environmental trajectories are disabled in {config_path}"
            )
        if config.migration.mode != "fixed_regional_pool":
            raise RuntimeError(f"migration is disabled in {config_path}")
        if config.migration.fraction != 0.1:
            raise RuntimeError(f"migration fraction is not 0.1 in {config_path}")
        if config.migration.regional_counts != config.environment.initial_counts:
            raise RuntimeError(
                f"regional and focal starting compositions differ in {config_path}"
            )
        if (
            config.evolution.within_host_selection
            or config.evolution.free_living_selection
        ):
            raise RuntimeError(f"selection is active in Phase 1 config {config_path}")


def _write_launcher_summary(
    results: list[dict[str, Any]], *, scratch: Path, started_at: str
) -> Path:
    directory = scratch / "p01-neutral-feedback" / STAGE_DIRECTORY / "_launcher"
    directory.mkdir(parents=True, exist_ok=True)
    stamp = datetime.now(UTC).strftime("%Y%m%dT%H%M%SZ")
    path = directory / f"launcher-summary-{stamp}.json"
    counts: dict[str, int] = {}
    for result in results:
        status = str(result["status"])
        counts[status] = counts.get(status, 0) + 1
    _atomic_json(
        path,
        {
            "launcher_summary_schema_version": "1.0.0",
            "experiment_id": EXPERIMENT_ID,
            "variant_id": VARIANT_TAG,
            "started_at": started_at,
            "completed_at": datetime.now(UTC).isoformat(),
            "status_counts": counts,
            "runs": results,
        },
    )
    return path


def _report_readiness_issues(
    rows: list[dict[str, str]], *, work: Path, scratch: Path
) -> list[str]:
    issues: list[str] = []
    if len(rows) != EXPECTED_RUNS:
        issues.append(f"manifest contains {len(rows)} runs; expected {EXPECTED_RUNS}")
    for row in rows:
        run_directory = scratch / row["scratch_relative_path"]
        completion_path = run_directory / "completion.json"
        if not completion_path.is_file():
            issues.append(f"{row['run_id']}: completion.json is missing")
            continue
        try:
            completion = json.loads(completion_path.read_text(encoding="utf-8"))
        except (OSError, ValueError) as exc:
            issues.append(f"{row['run_id']}: completion.json is unreadable ({exc})")
            continue
        if completion.get("complete") is not True:
            issues.append(f"{row['run_id']}: completion marker is not committed")
        config = load_config(work / row["config_path"])
        expected_resolved = json.loads(json.dumps(config.to_dict()))
        resolved_path = run_directory / "resolved_config.json"
        try:
            observed_resolved = json.loads(resolved_path.read_text(encoding="utf-8"))
        except (OSError, ValueError) as exc:
            issues.append(
                f"{row['run_id']}: resolved_config.json is unreadable ({exc})"
            )
            continue
        if observed_resolved != expected_resolved:
            issues.append(
                f"{row['run_id']}: resolved configuration differs from frozen TOML"
            )
        expected_resolved_hash = _resolved_config_sha256(expected_resolved)
        expected_versions = {
            "model_spec_version": MODEL_SPEC_VERSION,
            "output_schema_version": OUTPUT_SCHEMA_VERSION,
            "software_version": SOFTWARE_VERSION,
        }
        for field, expected in expected_versions.items():
            if completion.get(field) != expected:
                issues.append(
                    f"{row['run_id']}: {field}={completion.get(field)!r}; "
                    f"expected {expected!r}"
                )
        if completion.get("config_sha256") != expected_resolved_hash:
            issues.append(
                f"{row['run_id']}: completed resolved-configuration hash differs"
            )
        output_sizes = completion.get("output_sizes")
        if not isinstance(output_sizes, dict):
            issues.append(f"{row['run_id']}: output-size audit is missing")
            continue
        for name in ("environment_counts.csv", "host_generation_summary.csv"):
            path = run_directory / name
            if name not in output_sizes or not path.is_file():
                issues.append(
                    f"{row['run_id']}: required trajectory table {name} missing"
                )
            elif path.stat().st_size != output_sizes[name]:
                issues.append(f"{row['run_id']}: committed size differs for {name}")
    return issues


def _build_report(repository: Path, python: Path) -> int:
    command = [
        str(python),
        str(repository / "scripts/build_phase1_second_pilot_report.py"),
        "--repository",
        str(repository),
    ]
    print(
        f"All {EXPECTED_RUNS} populations passed the completion gate; "
        "analysing equilibrium and knitting the report.",
        flush=True,
    )
    result = subprocess.run(command, cwd=repository, check=False)
    if result.returncode:
        print(
            "Report generation failed, but completed simulation outputs were not "
            "modified. Retry with --report-only after resolving the reported issue.",
            flush=True,
        )
    return result.returncode


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--repository", type=Path, default=Path(__file__).resolve().parents[1]
    )
    parser.add_argument(
        "--python",
        type=Path,
        default=Path(sys.executable),
        help="Python interpreter in the trophosome environment",
    )
    parser.add_argument(
        "--jobs",
        type=int,
        default=8,
        help="independent populations to run concurrently (default: 8)",
    )
    parser.add_argument(
        "--cell",
        action="append",
        type=_normalise_cell_id,
        help="restrict to a sentinel cell; repeat for several",
    )
    parser.add_argument(
        "--seed-block",
        action="append",
        choices=sorted(EXPECTED_SEED_BLOCKS),
        help="restrict to a seed block; repeat for several",
    )
    parser.add_argument("--monitor-interval", type=float, default=2.0)
    parser.add_argument("--continue-on-error", action="store_true")
    parser.add_argument(
        "--allow-dirty-source",
        action="store_true",
        help="permit launching from an uncommitted source tree",
    )
    parser.add_argument(
        "--prepare-only",
        action="store_true",
        help="verify inputs and create run directories without simulating",
    )
    parser.add_argument(
        "--report-only",
        action="store_true",
        help=(
            f"audit all {EXPECTED_RUNS} completed populations, rerun the analysis, "
            "and rebuild the PDF and Markdown without simulating"
        ),
    )
    parser.add_argument(
        "--no-report",
        action="store_true",
        help=(
            "do not automatically build the report after all "
            f"{EXPECTED_RUNS} runs complete"
        ),
    )
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args()
    if args.jobs < 1:
        parser.error("--jobs must be at least 1")
    if args.monitor_interval <= 0:
        parser.error("--monitor-interval must be positive")
    if sum((args.prepare_only, args.report_only, args.dry_run)) > 1:
        parser.error("choose only one of --prepare-only, --report-only, or --dry-run")
    if args.report_only and args.no_report:
        parser.error("--report-only and --no-report cannot be combined")

    repository = args.repository.resolve()
    python = Path(os.path.abspath(args.python))
    if not python.is_file():
        parser.error(f"Python interpreter does not exist: {python}")

    work, scratch, rows = _load_rows(repository)
    selected = _select_rows(
        rows,
        cells=set(args.cell) if args.cell else None,
        seed_blocks=set(args.seed_block) if args.seed_block else None,
    )
    _verify_runtime(python)
    _verify_manifest(rows, selected, work=work)

    populations = len(selected)
    cell_count = len({row["cell_id"] for row in selected})
    selected_seed_count = len({row["seed_block_id"] for row in selected})
    print(
        f"Preflight passed: {populations} populations from {cell_count} sentinel "
        f"cells; 250 passages; {selected_seed_count} selected of "
        f"{len(EXPECTED_SEED_BLOCKS)} frozen seed blocks; m=0.1; up to "
        f"{args.jobs} populations and {args.jobs * 2} host workers concurrently.",
        flush=True,
    )
    if args.dry_run:
        for row in selected:
            print(
                f"{row['run_id']} -> {scratch / row['scratch_relative_path']}",
                flush=True,
            )
        return 0
    if args.report_only:
        issues = _report_readiness_issues(rows, work=work, scratch=scratch)
        if issues:
            print("The report completion gate failed:\n" + "\n".join(issues))
            return 1
        return _build_report(repository, python)

    _prepare_scratch(selected, work=work, scratch=scratch)
    if args.prepare_only:
        print(f"Prepared run directories below {scratch}", flush=True)
        return 0

    _require_frozen_source(repository, args.allow_dirty_source)
    started_at = datetime.now(UTC).isoformat()
    stop_event = threading.Event()
    results: list[dict[str, Any]] = []
    previous_sigint = signal.getsignal(signal.SIGINT)
    previous_sigterm = signal.getsignal(signal.SIGTERM)

    def request_stop(signum: int, _frame: object) -> None:
        if stop_event.is_set():
            raise KeyboardInterrupt
        print(
            f"{signal.Signals(signum).name} received; stopping active simulations "
            "and retaining their latest valid checkpoints...",
            flush=True,
        )
        stop_event.set()

    signal.signal(signal.SIGINT, request_stop)
    signal.signal(signal.SIGTERM, request_stop)
    try:
        with ThreadPoolExecutor(max_workers=args.jobs) as executor:
            futures = {
                executor.submit(
                    _run_one,
                    row,
                    repository=repository,
                    work=work,
                    scratch=scratch,
                    python=python,
                    monitor_interval=args.monitor_interval,
                    stop_event=stop_event,
                ): row
                for row in selected
            }
            for future in as_completed(futures):
                row = futures[future]
                try:
                    result = future.result()
                except Exception as exc:
                    result = {
                        "run_id": row["run_id"],
                        "status": "launcher-error",
                        "error": f"{type(exc).__name__}: {exc}",
                    }
                    print(f"[{row['run_id']}] launcher error: {exc}", flush=True)
                results.append(result)
                if result["status"] in {"failed", "launcher-error"}:
                    if not args.continue_on_error:
                        stop_event.set()
    except KeyboardInterrupt:
        stop_event.set()
    finally:
        signal.signal(signal.SIGINT, previous_sigint)
        signal.signal(signal.SIGTERM, previous_sigterm)

    summary_path = _write_launcher_summary(
        results, scratch=scratch, started_at=started_at
    )
    successful = sum(
        result["status"] in {"complete", "already-complete"} for result in results
    )
    print(
        f"Launcher finished: {successful}/{populations} selected populations are "
        f"complete. Summary: {summary_path}",
        flush=True,
    )
    if successful != populations:
        return 1
    if args.no_report:
        return 0
    readiness_issues = _report_readiness_issues(rows, work=work, scratch=scratch)
    if readiness_issues:
        print(
            f"Automatic report deferred because the complete {EXPECTED_RUNS}-run "
            "set is not yet available:\n" + "\n".join(readiness_issues[:20]),
            flush=True,
        )
        if len(readiness_issues) > 20:
            print(f"... and {len(readiness_issues) - 20} more issues", flush=True)
        return 0
    return _build_report(repository, python)


if __name__ == "__main__":
    raise SystemExit(main())
