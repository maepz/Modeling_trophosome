#!/usr/bin/env python3
"""Run the model-2.1 fixed-pool first pilot safely on an HPC node.

The launcher reads the frozen 60-run manifest, prepares isolated scratch
directories, runs several independent populations concurrently, records time,
memory and logs per population, and resumes from valid model checkpoints after
an interruption. Once all 60 populations pass the completion gate, it invokes
the dedicated migration-aware analysis and report-knitting workflow.
"""

from __future__ import annotations

import argparse
import csv
import fcntl
import json
import os
import resource
import signal
import subprocess
import sys
import threading
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import UTC, datetime
from pathlib import Path
from typing import Any

from prepare_phase1_first_pilot_v2_1 import (
    EXPERIMENT_ID,
    MODEL_SPEC_VERSION,
    OUTPUT_SCHEMA_VERSION,
    SOFTWARE_VERSION,
    VARIANT_TAG,
)
from run_phase1_first_pilot import (
    CORE_ORDER,
    EXTENSION_ORDER,
    PILOT_ORDER,
    _archive_unrestartable_attempt,
    _atomic_json,
    _directory_size,
    _normalise_cell_id,
    _process_tree_rss_kib,
    _sha256,
)

from trophosome import (
    MODEL_SPEC_VERSION as INSTALLED_MODEL_SPEC_VERSION,
)
from trophosome import (
    OUTPUT_SCHEMA_VERSION as INSTALLED_OUTPUT_SCHEMA_VERSION,
)
from trophosome import __version__ as INSTALLED_SOFTWARE_VERSION
from trophosome.config import load_config

MANIFEST_NAME = f"phase1-first-pilot-{VARIANT_TAG}-runs.tsv"
EXPECTED_SEED_BLOCKS = {"sb0001", "sb0002", "sb0003"}


def _load_rows(repository: Path) -> tuple[Path, Path, list[dict[str, str]]]:
    work = repository / "experiments" / "work" / "trophosome"
    layout_path = work / "layout.local.json"
    if not layout_path.is_file():
        raise RuntimeError(
            f"machine-local layout is missing: {layout_path}. "
            "Create it using the instructions under 'Create the machine-local "
            "storage layout' in scripts/hpc/README.md, then rerun this command."
        )
    layout = json.loads(layout_path.read_text(encoding="utf-8"))
    scratch = Path(layout["scratch"])
    manifest = work / "p01-neutral-feedback" / "manifests" / MANIFEST_NAME
    with manifest.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    return work, scratch, rows


def _select_rows(
    rows: list[dict[str, str]],
    *,
    tier: str,
    cells: set[str] | None,
    seed_blocks: set[str] | None,
) -> list[dict[str, str]]:
    if tier == "core":
        tier_cells = set(CORE_ORDER)
    elif tier == "extension":
        tier_cells = set(EXTENSION_ORDER)
    else:
        tier_cells = set(PILOT_ORDER)
    if cells is not None:
        tier_cells &= cells
    selected_seeds = seed_blocks or EXPECTED_SEED_BLOCKS
    order = {cell_id: index for index, cell_id in enumerate(PILOT_ORDER)}
    selected = [
        row
        for row in rows
        if row["cell_id"] in tier_cells and row["seed_block_id"] in selected_seeds
    ]
    selected.sort(
        key=lambda row: (
            order[row["cell_id"]],
            row["seed_block_id"],
        )
    )
    return selected


def _verify_runtime(python: Path) -> None:
    current = {
        "model_spec_version": INSTALLED_MODEL_SPEC_VERSION,
        "output_schema_version": INSTALLED_OUTPUT_SCHEMA_VERSION,
        "software_version": INSTALLED_SOFTWARE_VERSION,
    }
    expected = {
        "model_spec_version": MODEL_SPEC_VERSION,
        "output_schema_version": OUTPUT_SCHEMA_VERSION,
        "software_version": SOFTWARE_VERSION,
    }
    if current != expected:
        raise RuntimeError(
            f"launcher imported the wrong trophosome version: {current}; "
            f"expected {expected}"
        )

    command = [
        str(python),
        "-c",
        (
            "import json, trophosome; "
            "print(json.dumps({'model_spec_version': "
            "trophosome.MODEL_SPEC_VERSION, 'output_schema_version': "
            "trophosome.OUTPUT_SCHEMA_VERSION, 'software_version': "
            "trophosome.__version__}))"
        ),
    ]
    result = subprocess.run(command, check=True, capture_output=True, text=True)
    selected = json.loads(result.stdout.strip())
    if selected != expected:
        raise RuntimeError(
            f"selected Python contains the wrong trophosome version: {selected}; "
            f"expected {expected}"
        )


def _verify_manifest(
    rows: list[dict[str, str]],
    selected: list[dict[str, str]],
    *,
    work: Path,
) -> None:
    if len(rows) != 60:
        raise RuntimeError(f"expected 60 manifest rows, found {len(rows)}")
    if len({row["run_id"] for row in rows}) != 60:
        raise RuntimeError("run IDs are not unique")
    if len({row["scratch_relative_path"] for row in rows}) != 60:
        raise RuntimeError("scratch output paths are not unique")
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
        if config.migration.mode != "fixed_regional_pool":
            raise RuntimeError(f"migration is disabled in {config_path}")
        if config.migration.fraction != 0.1:
            raise RuntimeError(f"migration fraction is not 0.1 in {config_path}")
        if config.migration.regional_counts != config.environment.initial_counts:
            raise RuntimeError(
                f"regional and initial focal compositions differ in {config_path}"
            )
        if (
            config.evolution.within_host_selection
            or config.evolution.free_living_selection
        ):
            raise RuntimeError(f"selection is active in Phase 1 config {config_path}")


def _scratch_manifest(
    row: dict[str, str], work: Path, run_directory: Path
) -> dict[str, object]:
    return {
        "scratch_manifest_schema_version": "1.1.0",
        "experiment_id": row["experiment_id"],
        "variant_id": row["variant_id"],
        "run_id": row["run_id"],
        "cell_id": row["cell_id"],
        "pilot_tier": row.get("pilot_tier", "second-pilot"),
        "seed_block_id": row["seed_block_id"],
        "master_seed": int(row["master_seed"]),
        "within_run_replicate_index": int(row["within_run_replicate_index"]),
        "config_path": str(work / row["config_path"]),
        "config_sha256": row["config_sha256"],
        "output_path": str(run_directory),
        "status": "prepared-not-launched",
    }


def _prepare_scratch(
    selected: list[dict[str, str]], *, work: Path, scratch: Path
) -> None:
    for row in selected:
        run_directory = scratch / row["scratch_relative_path"]
        run_directory.mkdir(parents=True, exist_ok=True)
        path = run_directory / "run.json"
        content = (
            json.dumps(
                _scratch_manifest(row, work, run_directory),
                indent=2,
                sort_keys=True,
            )
            + "\n"
        )
        if path.exists() and path.read_text(encoding="utf-8") != content:
            other_artifacts = [
                item for item in run_directory.iterdir() if item.name != "run.json"
            ]
            if other_artifacts:
                raise RuntimeError(
                    f"refusing to replace a different run manifest after launch: {path}"
                )
        path.write_text(content, encoding="utf-8")


def _require_frozen_source(repository: Path, allow_dirty_source: bool) -> None:
    result = subprocess.run(
        ["git", "status", "--porcelain"],
        cwd=repository,
        capture_output=True,
        text=True,
        check=False,
    )
    if result.returncode != 0:
        print(
            "Warning: no Git worktree was detected; the model will still record "
            "its source checksum.",
            flush=True,
        )
        return
    changed = [line for line in result.stdout.splitlines() if line.strip()]
    if changed and not allow_dirty_source:
        preview = "\n".join(changed[:12])
        suffix = "\n..." if len(changed) > 12 else ""
        raise RuntimeError(
            "the source tree is not frozen. Commit the model and generated pilot "
            "files before launching, or deliberately pass --allow-dirty-source.\n"
            f"{preview}{suffix}"
        )


def _run_one(
    row: dict[str, str],
    *,
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
    lock_directory = scratch / "_launch-locks" / row["variant_id"]
    lock_directory.mkdir(parents=True, exist_ok=True)
    lock_path = lock_directory / f"{row['run_id']}.lock"

    with lock_path.open("a+", encoding="utf-8") as lock:
        try:
            fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError:
            print(f"[{row['run_id']}] already controlled by another launcher")
            return {"run_id": row["run_id"], "status": "locked"}

        completion = output / "completion.json"
        if completion.is_file():
            print(f"[{row['run_id']}] already complete", flush=True)
            return {
                "run_id": row["run_id"],
                "status": "already-complete",
                "output_bytes": _directory_size(output),
            }

        run_artifacts_exist = any(
            (output / name).exists()
            for name in ("resolved_config.json", "provenance.json", "checkpoints")
        )
        checkpoint_exists = any(
            (output / "checkpoints").glob("checkpoint-rep*-gen*.npz")
        )
        if run_artifacts_exist and not checkpoint_exists:
            archive = _archive_unrestartable_attempt(output, scratch, row["run_id"])
            print(
                f"[{row['run_id']}] preserved pre-checkpoint partial attempt at "
                f"{archive}; restarting cleanly",
                flush=True,
            )
            run_artifacts_exist = False
        resume = run_artifacts_exist

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
        ]
        if resume:
            command.append("--resume")

        summary_path = output / "execution-summary.json"
        summary: dict[str, Any] = {
            "execution_summary_schema_version": "1.1.0",
            "experiment_id": row["experiment_id"],
            "variant_id": row["variant_id"],
            "run_id": row["run_id"],
            "cell_id": row["cell_id"],
            "seed_block_id": row["seed_block_id"],
            "started_at": datetime.now(UTC).isoformat(),
            "status": "running",
            "resumed": resume,
            "command": command,
        }
        _atomic_json(summary_path, summary)
        print(
            f"[{row['run_id']}] {'resuming' if resume else 'starting'}",
            flush=True,
        )

        started = time.monotonic()
        peak_rss_kib = 0
        memory_measurements = 0
        interrupt_started: float | None = None
        return_code = 1
        stdout_path = output / "run.out"
        stderr_path = output / "run.err"
        with (
            stdout_path.open("a", encoding="utf-8") as stdout,
            stderr_path.open("a", encoding="utf-8") as stderr,
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
                rss_kib = _process_tree_rss_kib(process.pid)
                if rss_kib is not None:
                    peak_rss_kib = max(peak_rss_kib, rss_kib)
                    memory_measurements += 1
                if stop_event.is_set():
                    if interrupt_started is None:
                        os.killpg(process.pid, signal.SIGINT)
                        interrupt_started = time.monotonic()
                    elif time.monotonic() - interrupt_started > 30:
                        os.killpg(process.pid, signal.SIGTERM)
                time.sleep(monitor_interval)
            return_code = process.returncode

        elapsed = time.monotonic() - started
        if not memory_measurements:
            child_peak = int(resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss)
            peak_rss_kib = (
                child_peak // 1024 if sys.platform == "darwin" else child_peak
            )
            memory_mode = "cumulative-child-peak"
        else:
            memory_mode = "process-tree-polling"

        successful = return_code == 0 and completion.is_file()
        if successful:
            status = "complete"
        elif stop_event.is_set():
            status = "interrupted"
        else:
            status = "failed"
        summary.update(
            {
                "completed_at": datetime.now(UTC).isoformat(),
                "elapsed_seconds": elapsed,
                "peak_process_tree_rss_kib": peak_rss_kib or None,
                "memory_measurement_mode": memory_mode,
                "memory_measurements": memory_measurements,
                "return_code": return_code,
                "output_bytes": _directory_size(output),
                "status": status,
            }
        )
        _atomic_json(summary_path, summary)
        print(
            f"[{row['run_id']}] {status} in {elapsed / 60:.1f} min; "
            f"peak RSS {peak_rss_kib / 1024:.0f} MiB; "
            f"output {summary['output_bytes'] / 1024**2:.1f} MiB",
            flush=True,
        )
        return summary


def _write_launcher_summary(
    results: list[dict[str, Any]], *, scratch: Path, started_at: str
) -> Path:
    directory = (
        scratch / "p01-neutral-feedback" / f"s01-pilot-{VARIANT_TAG}" / "_launcher"
    )
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


def _report_readiness_issues(rows: list[dict[str, str]], *, scratch: Path) -> list[str]:
    issues: list[str] = []
    if len(rows) != 60:
        issues.append(f"manifest contains {len(rows)} runs; expected 60")
    for row in rows:
        completion_path = scratch / row["scratch_relative_path"] / "completion.json"
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
    return issues


def _build_report(repository: Path, python: Path) -> int:
    command = [
        str(python),
        str(repository / "scripts" / "build_phase1_first_pilot_v2_1_report.py"),
        "--repository",
        str(repository),
    ]
    print("All 60 populations passed the completion gate; knitting report.", flush=True)
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
        "--repository",
        type=Path,
        default=Path(__file__).resolve().parents[1],
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
        "--tier",
        choices=("core", "extension", "all"),
        default="all",
        help="run the 12 core cells, 8 extension cells, or all 20",
    )
    parser.add_argument(
        "--cell",
        action="append",
        type=_normalise_cell_id,
        help="restrict to a cell; repeat for several",
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
        help="audit all completed populations and knit the report without simulating",
    )
    parser.add_argument(
        "--no-report",
        action="store_true",
        help="do not automatically knit a report when all 60 populations are complete",
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
        tier=args.tier,
        cells=set(args.cell) if args.cell else None,
        seed_blocks=set(args.seed_block) if args.seed_block else None,
    )
    _verify_runtime(python)
    _verify_manifest(rows, selected, work=work)

    populations = len(selected)
    cell_count = len({row["cell_id"] for row in selected})
    print(
        f"Preflight passed: {populations} populations from {cell_count} cells; "
        f"m=0.1; up to {args.jobs} populations and {args.jobs * 2} host workers "
        "concurrently.",
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
        issues = _report_readiness_issues(rows, scratch=scratch)
        if issues:
            print(
                "The report completion gate failed:\n" + "\n".join(issues),
                flush=True,
            )
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
        signal_name = signal.Signals(signum).name
        print(
            f"{signal_name} received; stopping active simulations and retaining "
            "their latest valid checkpoints...",
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
                except Exception as exc:  # preserve all other independent outputs
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
    readiness_issues = _report_readiness_issues(rows, scratch=scratch)
    if readiness_issues:
        print(
            "Automatic report deferred because the complete 60-run set is not yet "
            "available:\n" + "\n".join(readiness_issues[:20]),
            flush=True,
        )
        if len(readiness_issues) > 20:
            print(f"... and {len(readiness_issues) - 20} more issues", flush=True)
        return 0
    return _build_report(repository, python)


if __name__ == "__main__":
    raise SystemExit(main())
