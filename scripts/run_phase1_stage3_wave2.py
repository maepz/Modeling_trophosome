#!/usr/bin/env python3
"""Run Wave 2 to passage 100, then only frozen adaptive continuations."""

from __future__ import annotations

import argparse
import csv
import fcntl
import json
import os
import resource
import shutil
import signal
import subprocess
import sys
import threading
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import UTC, datetime
from pathlib import Path
from typing import Any

from assess_phase1_stage3_wave2_horizon import decision_path, write_or_verify
from prepare_phase1_stage3_wave2 import (
    EXPECTED_RUNS,
    EXPERIMENT_ID,
    INITIAL_HORIZON,
    INTERMEDIATE_HORIZON,
    MAXIMUM_HORIZON,
    MODEL_SPEC_VERSION,
    OUTPUT_SCHEMA_VERSION,
    SMOKE_CELLS,
    SOFTWARE_VERSION,
    STAGE_DIRECTORY,
    VARIANT_TAG,
    verify_files,
)
from run_phase1_first_pilot import (
    _archive_unrestartable_attempt,
    _atomic_json,
    _directory_size,
    _process_tree_rss_kib,
)
from run_phase1_first_pilot import (
    _sha256 as _sha256_file,
)
from run_phase1_first_pilot_v2_1 import (
    _prepare_scratch,
    _require_frozen_source,
    _verify_runtime,
)
from run_phase1_second_pilot import _resolved_config_sha256

from trophosome.checkpointing import load_recovery_checkpoint
from trophosome.config import load_config
from trophosome.simulation import _output_fields

HORIZONS = (INITIAL_HORIZON, INTERMEDIATE_HORIZON, MAXIMUM_HORIZON)


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
    path = work / "p01-neutral-feedback/manifests" / f"{EXPERIMENT_ID}-runs.tsv"
    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    return work, scratch, rows


def selected_cells_for_horizon(work: Path, horizon: int) -> set[str] | None:
    if horizon == INITIAL_HORIZON:
        return None
    previous_horizon = (
        INITIAL_HORIZON if horizon == INTERMEDIATE_HORIZON else INTERMEDIATE_HORIZON
    )
    path = decision_path(work, previous_horizon)
    if not path.is_file():
        raise RuntimeError(
            f"the frozen passage-{previous_horizon} adaptive decision is missing; "
            f"run --assess-only --horizon {previous_horizon} first"
        )
    payload = json.loads(path.read_text(encoding="utf-8"))
    if payload.get("continuation_horizon") != horizon:
        raise RuntimeError(f"adaptive decision does not authorize passage {horizon}")
    return set(payload["selected_cell_ids"])


def _verify_manifest(rows: list[dict[str, str]], *, work: Path) -> None:
    if len(rows) != EXPECTED_RUNS:
        raise RuntimeError(f"expected {EXPECTED_RUNS} new runs, found {len(rows)}")
    if len({row["run_id"] for row in rows}) != EXPECTED_RUNS:
        raise RuntimeError("Wave 2 run IDs are not unique")
    if len({row["scratch_relative_path"] for row in rows}) != EXPECTED_RUNS:
        raise RuntimeError("Wave 2 scratch paths are not unique")
    for row in rows:
        path = work / row["config_path"]
        if _sha256_file(path) != row["config_sha256"]:
            raise RuntimeError(f"configuration checksum differs: {path}")
        config = load_config(path)
        if config.host.host_generations != MAXIMUM_HORIZON:
            raise RuntimeError(f"maximum horizon differs in {path}")
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


def _state_generation(output: Path) -> int | None:
    if (output / "completion.json").is_file():
        return MAXIMUM_HORIZON
    pause = output / "pause.json"
    if pause.is_file():
        return int(
            json.loads(pause.read_text(encoding="utf-8"))["last_completed_generation"]
        )
    return None


def state_issues(
    rows: list[dict[str, str]], *, work: Path, scratch: Path, horizon: int
) -> list[str]:
    """Audit a paused or completed state at or beyond the requested horizon."""

    issues: list[str] = []
    for row in rows:
        output = scratch / row["scratch_relative_path"]
        label = row["run_id"]
        try:
            config_path = work / row["config_path"]
            if _sha256_file(config_path) != row["config_sha256"]:
                raise ValueError("raw configuration hash differs")
            config = load_config(config_path)
            expected = json.loads(json.dumps(config.to_dict()))
            if json.loads((output / "resolved_config.json").read_text()) != expected:
                raise ValueError("resolved configuration differs")
            state_generation = _state_generation(output)
            if state_generation is None or state_generation < horizon:
                raise ValueError(
                    f"trajectory reached {state_generation}, not passage {horizon}"
                )
            completion = output / "completion.json"
            if completion.is_file():
                record = json.loads(completion.read_text())
                if record.get("complete") is not True:
                    raise ValueError("completion marker is not committed")
                size_field = "output_sizes"
            else:
                pause_path = output / "pause.json"
                record = json.loads(pause_path.read_text())
                if record.get("status") != "paused":
                    raise ValueError("pause marker is invalid")
                checkpoint_path = output / "checkpoints" / record["checkpoint"]
                if _sha256_file(checkpoint_path) != record["checkpoint_sha256"]:
                    raise ValueError("paused checkpoint checksum differs")
                checkpoint = load_recovery_checkpoint(checkpoint_path)
                if checkpoint.metadata["last_completed_generation"] != state_generation:
                    raise ValueError("pause marker and checkpoint generation differ")
                size_field = "output_sizes"
            for field, value in (
                ("config_sha256", _resolved_config_sha256(expected)),
                ("model_spec_version", MODEL_SPEC_VERSION),
                ("output_schema_version", OUTPUT_SCHEMA_VERSION),
                ("software_version", SOFTWARE_VERSION),
            ):
                if record.get(field) != value:
                    raise ValueError(f"{field} differs")
            sizes = record.get(size_field, {})
            if set(sizes) != set(_output_fields(config)):
                raise ValueError("output table inventory differs")
            for name, size in sizes.items():
                if (output / name).stat().st_size != int(size):
                    raise ValueError(f"committed size differs for {name}")
            with (output / "host_generation_summary.csv").open(
                newline="", encoding="utf-8"
            ) as handle:
                summaries = list(csv.DictReader(handle))
            if len(summaries) != state_generation:
                raise ValueError("generation summary length differs")
            execution = json.loads((output / "execution-summary.json").read_text())
            if execution.get("run_id") != label or execution.get("status") not in {
                "paused",
                "complete",
                "already-paused",
                "already-complete",
                "already-beyond-target",
            }:
                raise ValueError("execution summary does not record a valid state")
        except (OSError, ValueError, KeyError, TypeError) as exc:
            issues.append(f"{label}: {exc}")
    return issues


def _run_one(
    row: dict[str, str],
    *,
    horizon: int,
    repository: Path,
    work: Path,
    scratch: Path,
    python: Path,
    monitor_interval: float,
    stop_event: threading.Event,
) -> dict[str, Any]:
    if stop_event.is_set():
        return {"run_id": row["run_id"], "status": "cancelled-before-start"}
    config = work / row["config_path"]
    output = scratch / row["scratch_relative_path"]
    output.mkdir(parents=True, exist_ok=True)
    lock_directory = scratch / "_launch-locks" / VARIANT_TAG
    lock_directory.mkdir(parents=True, exist_ok=True)
    lock_path = lock_directory / f"{row['run_id']}.lock"

    with lock_path.open("a+", encoding="utf-8") as lock:
        try:
            fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError:
            return {"run_id": row["run_id"], "status": "locked"}

        reached = _state_generation(output)
        if reached is not None and reached >= horizon:
            status = (
                "already-complete"
                if reached == MAXIMUM_HORIZON
                else "already-paused"
                if reached == horizon
                else "already-beyond-target"
            )
            summary = {
                "execution_summary_schema_version": "1.2.0",
                "experiment_id": EXPERIMENT_ID,
                "variant_id": VARIANT_TAG,
                "run_id": row["run_id"],
                "cell_id": row["cell_id"],
                "seed_block_id": row["seed_block_id"],
                "status": status,
                "target_horizon": horizon,
                "reached_horizon": reached,
                "output_bytes": _directory_size(output),
            }
            _atomic_json(output / "execution-summary.json", summary)
            return summary

        run_artifacts = any(
            (output / name).exists()
            for name in ("resolved_config.json", "provenance.json", "checkpoints")
        )
        checkpoint_exists = any(
            (output / "checkpoints").glob("checkpoint-rep*-gen*.npz")
        )
        if run_artifacts and not checkpoint_exists:
            archive = _archive_unrestartable_attempt(output, scratch, row["run_id"])
            print(f"[{row['run_id']}] archived incomplete attempt at {archive}")
            run_artifacts = False
        resume = run_artifacts
        command = [
            str(python),
            "-m",
            "trophosome.cli",
            "run",
            str(config),
            "--output",
            str(output),
            "--repository",
            str(repository),
            "--pause-after-generation",
            str(horizon),
        ]
        if resume:
            command.append("--resume")
        summary: dict[str, Any] = {
            "execution_summary_schema_version": "1.2.0",
            "experiment_id": EXPERIMENT_ID,
            "variant_id": VARIANT_TAG,
            "run_id": row["run_id"],
            "cell_id": row["cell_id"],
            "seed_block_id": row["seed_block_id"],
            "started_at": datetime.now(UTC).isoformat(),
            "status": "running",
            "target_horizon": horizon,
            "resumed": resume,
            "command": command,
        }
        _atomic_json(output / "execution-summary.json", summary)
        print(f"[{row['run_id']}] {'resuming' if resume else 'starting'} to {horizon}")
        started = time.monotonic()
        peak_rss_kib = 0
        measurements = 0
        interrupt_started: float | None = None
        with (
            (output / "run.out").open("a", encoding="utf-8") as stdout,
            (output / "run.err").open("a", encoding="utf-8") as stderr,
        ):
            process = subprocess.Popen(
                command,
                cwd=repository,
                stdout=stdout,
                stderr=stderr,
                text=True,
                start_new_session=True,
                env={**os.environ, "PYTHONUNBUFFERED": "1"},
            )
            while process.poll() is None:
                rss = _process_tree_rss_kib(process.pid)
                if rss is not None:
                    peak_rss_kib = max(peak_rss_kib, rss)
                    measurements += 1
                if stop_event.is_set():
                    if interrupt_started is None:
                        os.killpg(process.pid, signal.SIGINT)
                        interrupt_started = time.monotonic()
                    elif time.monotonic() - interrupt_started > 30:
                        os.killpg(process.pid, signal.SIGTERM)
                time.sleep(monitor_interval)
        elapsed = time.monotonic() - started
        if not measurements:
            child_peak = int(resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss)
            peak_rss_kib = (
                child_peak // 1024 if sys.platform == "darwin" else child_peak
            )
        reached = _state_generation(output)
        successful = (
            process.returncode == 0 and reached is not None and reached >= horizon
        )
        if successful:
            status = "complete" if reached == MAXIMUM_HORIZON else "paused"
        elif stop_event.is_set():
            status = "interrupted"
        else:
            status = "failed"
        summary.update(
            {
                "completed_at": datetime.now(UTC).isoformat(),
                "elapsed_seconds": elapsed,
                "peak_process_tree_rss_kib": peak_rss_kib or None,
                "memory_measurements": measurements,
                "return_code": process.returncode,
                "reached_horizon": reached,
                "output_bytes": _directory_size(output),
                "status": status,
            }
        )
        _atomic_json(output / "execution-summary.json", summary)
        print(
            f"[{row['run_id']}] {status} at passage {reached}; {elapsed / 60:.1f} min"
        )
        return summary


def smoke_assessment(rows: list[dict[str, str]], *, work: Path, scratch: Path) -> dict:
    selected = [
        row
        for row in rows
        if row["cell"] in SMOKE_CELLS and row["seed_block_id"] == "sb0001"
    ]
    issues = state_issues(selected, work=work, scratch=scratch, horizon=INITIAL_HORIZON)
    measurements: dict[str, dict[str, float]] = {}
    for row in selected:
        try:
            output = scratch / row["scratch_relative_path"]
            execution = json.loads((output / "execution-summary.json").read_text())
            if execution.get("resumed"):
                issues.append(
                    f"{row['cell']}: resumed smoke timing needs manual review"
                )
            hours = float(execution["elapsed_seconds"]) / 3600
            if not 0 < hours <= 48:
                issues.append(f"{row['cell']}: runtime is outside 0-48 hours")
            measurements[row["cell"]] = {
                "hours": hours,
                "gib": _directory_size(output) / 1024**3,
            }
        except (OSError, KeyError, TypeError, ValueError) as exc:
            issues.append(f"{row['run_id']}: cannot assess resources ({exc})")
    projected_gib = (
        sum(value["gib"] for value in measurements.values())
        * (EXPECTED_RUNS / max(1, len(measurements)))
        * 2
    )
    if projected_gib > 350:
        issues.append("conservative initial-wave projection exceeds 350 GiB")
    free_gib = shutil.disk_usage(scratch).free / 1024**3
    if free_gib < projected_gib:
        issues.append("scratch has less free space than the conservative projection")
    return {
        "passed": not issues,
        "issues": issues,
        "measurements": measurements,
        "projected_initial_output_gib": projected_gib,
        "free_scratch_gib": free_gib,
        "method": (
            "three included passage-100 runs; mean extrapolation with 2x allowance"
        ),
    }


def main() -> int:
    parser = argparse.ArgumentParser(
        description=__doc__,
        epilog=(
            "Run --prepare-only, --smoke-only, and --check-smoke before the full "
            "passage-100 batch. Later horizons require the frozen decision from "
            "the preceding assessment. No extension changes the passage-100 endpoint."
        ),
    )
    parser.add_argument(
        "--repository",
        type=Path,
        default=Path(__file__).resolve().parents[1],
        help="repository containing the frozen Wave 2 files",
    )
    parser.add_argument(
        "--python",
        default=sys.executable,
        help="Python executable from the trophosome environment",
    )
    parser.add_argument(
        "--jobs",
        type=int,
        default=8,
        help="number of independent populations to run at the same time",
    )
    parser.add_argument(
        "--monitor-interval",
        type=float,
        default=2,
        help="seconds between memory and interruption checks",
    )
    parser.add_argument(
        "--horizon",
        type=int,
        choices=HORIZONS,
        default=INITIAL_HORIZON,
        help=(
            "passage boundary to reach; 500 and 1000 use only conditions "
            "authorized by the preceding frozen decision"
        ),
    )
    parser.add_argument("--cell", action="append", help="restrict to c0051-style cell")
    parser.add_argument(
        "--seed-block", action="append", help="restrict to sb0001-style seed"
    )
    parser.add_argument(
        "--smoke-only",
        action="store_true",
        help="run the three included passage-100 resource checks only",
    )
    parser.add_argument(
        "--no-assessment",
        action="store_true",
        help="finish the selected boundary without creating its adaptive decision",
    )
    parser.add_argument(
        "--allow-dirty-source",
        action="store_true",
        help="allow an explicitly reviewed uncommitted source for the initial run",
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
        help="list the authorized populations and destinations; do not write",
    )
    mode.add_argument(
        "--check-smoke",
        action="store_true",
        help="audit the three safety runs and estimate initial-batch resources",
    )
    mode.add_argument(
        "--assess-only",
        action="store_true",
        help="freeze the adaptive decision from existing results; do not simulate",
    )
    args = parser.parse_args()
    if args.jobs < 1 or args.monitor_interval <= 0:
        parser.error("jobs and monitor interval must be positive")
    if args.smoke_only and args.horizon != INITIAL_HORIZON:
        parser.error("smoke-only applies only to the passage-100 batch")
    repository = args.repository.resolve()
    differences = verify_files(repository)
    if differences:
        raise SystemExit("Frozen Wave 2 inputs differ:\n" + "\n".join(differences))
    work, scratch, rows = load_rows(repository)
    _verify_runtime(Path(os.path.abspath(args.python)))
    _verify_manifest(rows, work=work)
    authorized_cells = selected_cells_for_horizon(work, args.horizon)
    selected = [
        row
        for row in rows
        if (authorized_cells is None or row["cell_id"] in authorized_cells)
        and (not args.cell or row["cell"] in set(args.cell))
        and (not args.seed_block or row["seed_block_id"] in set(args.seed_block))
        and (
            not args.smoke_only
            or (row["cell"] in SMOKE_CELLS and row["seed_block_id"] == "sb0001")
        )
    ]
    if not selected and not args.assess_only:
        raise SystemExit("The requested Wave 2 selection contains no new runs")
    print(
        f"Preflight passed: {len(selected)} populations toward passage "
        f"{args.horizon}; up to {args.jobs} populations and {2 * args.jobs} "
        "host workers concurrently.",
        flush=True,
    )
    if args.dry_run:
        for row in selected:
            print(f"{row['run_id']} -> {scratch / row['scratch_relative_path']}")
        return 0
    if args.assess_only:
        if args.horizon == MAXIMUM_HORIZON:
            parser.error("there is no extension decision after passage 1000")
        write_or_verify(repository, args.horizon, verify=False)
        return 0
    if args.check_smoke or (
        args.horizon == INITIAL_HORIZON
        and not args.smoke_only
        and not args.prepare_only
        and not args.cell
        and not args.seed_block
    ):
        assessment = smoke_assessment(rows, work=work, scratch=scratch)
        print(json.dumps(assessment, indent=2), flush=True)
        if not assessment["passed"]:
            print(
                "Full Wave 2 batch not launched; complete/review the smoke runs first."
            )
            return 1
        if args.check_smoke:
            return 0
    _prepare_scratch(selected, work=work, scratch=scratch)
    if args.prepare_only:
        print("Prepared isolated Wave 2 directories; no simulations launched.")
        return 0
    needs_computation = any(
        (_state_generation(scratch / row["scratch_relative_path"]) or 0) < args.horizon
        for row in selected
    )
    if needs_computation and args.horizon == INITIAL_HORIZON:
        _require_frozen_source(repository, args.allow_dirty_source)
    elif needs_computation:
        # Adaptive decisions are derived, intentionally uncommitted files.  A
        # full Git-clean check would therefore block every continuation.  The
        # model itself validates the exact src/trophosome + pyproject checksum
        # against each recovery checkpoint before appending any row.
        print(
            "Continuation will verify the exact model-source checksum stored "
            "in every checkpoint.",
            flush=True,
        )
    stop = threading.Event()
    previous_handlers = {
        sig: signal.getsignal(sig) for sig in (signal.SIGINT, signal.SIGTERM)
    }

    def request_stop(signum, _frame):
        print(f"{signal.Signals(signum).name}: stopping; checkpoints are retained.")
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
                    horizon=args.horizon,
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
                if result["status"] not in {
                    "paused",
                    "complete",
                    "already-paused",
                    "already-complete",
                    "already-beyond-target",
                }:
                    stop.set()
    finally:
        for sig, handler in previous_handlers.items():
            signal.signal(sig, handler)
    launcher = scratch / "p01-neutral-feedback" / STAGE_DIRECTORY / "_launcher"
    launcher.mkdir(parents=True, exist_ok=True)
    stamp = datetime.now(UTC).strftime("%Y%m%dT%H%M%S%fZ")
    _atomic_json(
        launcher / f"summary-g{args.horizon}-{stamp}.json",
        {
            "experiment_id": EXPERIMENT_ID,
            "target_horizon": args.horizon,
            "runs": results,
        },
    )
    successful = {
        "paused",
        "complete",
        "already-paused",
        "already-complete",
        "already-beyond-target",
    }
    if len(results) != len(selected) or any(
        r["status"] not in successful for r in results
    ):
        print("Some runs are incomplete. Review run.err, then repeat to resume.")
        return 1
    issues = state_issues(selected, work=work, scratch=scratch, horizon=args.horizon)
    if issues:
        raise SystemExit("Wave 2 state audit failed:\n" + "\n".join(issues))
    if args.smoke_only:
        print(
            "Three included smoke populations reached passage 100; run --check-smoke."
        )
        return 0
    complete_batch = not args.cell and not args.seed_block
    if complete_batch and not args.no_assessment and args.horizon < MAXIMUM_HORIZON:
        write_or_verify(repository, args.horizon, verify=False)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
