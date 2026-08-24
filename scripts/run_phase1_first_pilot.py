#!/usr/bin/env python3
"""Run the staged Phase 1 first pilot and record resource use.

With no cell or seed selectors, the launcher performs the complete workflow:
run or verify the core 12 cells, analyse and report them, apply the expansion
safety gates, run the eight extension cells, then overwrite the report with the
20-cell result. Supplying ``--cell`` or ``--seed-block`` retains the low-level,
scheduler-independent mode for selected populations.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import resource
import shutil
import signal
import subprocess
import sys
import time
from datetime import UTC, datetime
from pathlib import Path
from typing import Any

EXPERIMENT_ID = "phase1-first-pilot-20cell"
CORE_ORDER = (
    "p01-s01-c0001",
    "p01-s01-c0002",
    "p01-s01-c0003",
    "p01-s01-c0007",
    "p01-s01-c0008",
    "p01-s01-c0009",
    "p01-s01-c0010",
    "p01-s01-c0011",
    "p01-s01-c0004",
    "p01-s01-c0012",
    "p01-s01-c0005",
    "p01-s01-c0006",
)
EXTENSION_ORDER = (
    "p01-s01-c0018",
    "p01-s01-c0019",
    "p01-s01-c0020",
    "p01-s01-c0016",
    "p01-s01-c0013",
    "p01-s01-c0017",
    "p01-s01-c0014",
    "p01-s01-c0015",
)
PILOT_ORDER = (*CORE_ORDER, *EXTENSION_ORDER)
SEED_BLOCKS = ("sb0001", "sb0002", "sb0003")


def _experiment_id(cell_id: str) -> str:
    number = int(cell_id.rsplit("c", 1)[1])
    if number <= 12:
        return "phase1-first-pilot-core12"
    if number <= 17:
        return "phase1-first-pilot-extension5"
    return "phase1-first-pilot-alpha-extension3"


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _atomic_json(path: Path, payload: dict[str, Any]) -> None:
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    os.replace(temporary, path)


def _directory_size(path: Path) -> int:
    return sum(item.stat().st_size for item in path.rglob("*") if item.is_file())


def _archive_unrestartable_attempt(
    output: Path, scratch: Path, run_id: str
) -> Path:
    """Preserve a partial pre-checkpoint attempt and recreate a clean run folder."""
    timestamp = datetime.now(UTC).strftime("%Y%m%dT%H%M%SZ")
    archive_parent = scratch / "_failed-attempts" / run_id
    archive_parent.mkdir(parents=True, exist_ok=True)
    archive = archive_parent / timestamp
    suffix = 1
    while archive.exists():
        archive = archive_parent / f"{timestamp}-{suffix:02d}"
        suffix += 1
    shutil.move(str(output), str(archive))
    output.mkdir(parents=True, exist_ok=True)
    archived_manifest = archive / "run.json"
    if archived_manifest.is_file():
        shutil.copy2(archived_manifest, output / "run.json")
    return archive


def _process_tree_rss_kib(root_pid: int) -> int | None:
    """Return summed resident memory for a process and its descendants."""

    try:
        result = subprocess.run(
            ["ps", "-axo", "pid=,ppid=,rss="],
            check=True,
            capture_output=True,
            text=True,
        )
    except (OSError, subprocess.SubprocessError):
        return None
    rows: list[tuple[int, int, int]] = []
    for line in result.stdout.splitlines():
        fields = line.split()
        if len(fields) != 3:
            continue
        try:
            rows.append(tuple(int(field) for field in fields))
        except ValueError:
            continue
    descendants = {root_pid}
    changed = True
    while changed:
        changed = False
        for pid, parent, _rss in rows:
            if parent in descendants and pid not in descendants:
                descendants.add(pid)
                changed = True
    return sum(rss for pid, _parent, rss in rows if pid in descendants)


def _normalise_cell_id(value: str) -> str:
    if value.startswith("p01-s01-c"):
        return value
    if value.startswith("c") and value[1:].isdigit():
        return f"p01-s01-c{int(value[1:]):04d}"
    if value.isdigit():
        return f"p01-s01-c{int(value):04d}"
    raise argparse.ArgumentTypeError("cell must look like c0001, 1, or p01-s01-c0001")


def _load_rows(repository: Path) -> tuple[Path, list[dict[str, str]]]:
    work = repository / "experiments" / "work" / "trophosome"
    layout = json.loads((work / "layout.local.json").read_text(encoding="utf-8"))
    scratch = Path(layout["scratch"])
    manifest = (
        work / "p01-neutral-feedback" / "manifests" / "phase1-first-pilot-runs.tsv"
    )
    with manifest.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    return scratch, rows


def _select_rows(
    rows: list[dict[str, str]],
    seed_blocks: set[str],
    cells: set[str] | None,
) -> list[dict[str, str]]:
    order = {cell: index for index, cell in enumerate(PILOT_ORDER)}
    selected = [
        row
        for row in rows
        if row["seed_block_id"] in seed_blocks
        and (cells is None or row["cell_id"] in cells)
    ]
    selected.sort(
        key=lambda row: (
            row["seed_block_id"],
            order.get(row["cell_id"], len(order)),
        )
    )
    return selected


def _run_rows(
    selected: list[dict[str, str]],
    *,
    repository: Path,
    work: Path,
    scratch: Path,
    python: Path,
    monitor_interval: float,
    continue_on_error: bool,
) -> tuple[list[dict[str, Any]], int]:
    results = []
    exit_code = 0
    for row in selected:
        result = _run_one(
            row,
            repository,
            work,
            scratch,
            python,
            monitor_interval,
        )
        results.append(result)
        if result["status"] == "failed":
            exit_code = 1
            if not continue_on_error:
                break
    return results, exit_code


def _checked(command: list[str], repository: Path) -> None:
    print("Workflow step: " + " ".join(command), flush=True)
    subprocess.run(command, cwd=repository, check=True)


def _run_full_workflow(
    *,
    repository: Path,
    work: Path,
    scratch: Path,
    rows: list[dict[str, str]],
    python: Path,
    monitor_interval: float,
    continue_on_error: bool,
    core_only: bool,
) -> int:
    phase = work / "p01-neutral-feedback"
    analysis_script = phase / "analysis" / "analyse_first_pilot.py"
    core_analysis = phase / "analysis" / "s01-pilot-core-derived"
    full_analysis = phase / "analysis" / "s01-pilot-derived"
    design = phase / "design" / "phase1-first-pilot-cells.tsv"
    core_design = phase / "design" / "phase1-first-pilot-core-cells.tsv"
    report = repository / "output" / "pdf" / "phase1-first-pilot-report.pdf"
    markdown = repository / "docs" / "phase1-first-pilot-report.md"
    assets = repository / "docs" / "figures" / "phase1-pilot-report"
    safety = phase / "analysis" / "phase1-first-pilot-expansion-safety.json"

    _checked(
        [str(python), "scripts/prepare_phase1_first_pilot.py", "--verify"],
        repository,
    )
    _checked(
        [str(python), "scripts/prepare_phase1_pilot_scratch.py", "--write"],
        repository,
    )

    core_rows = _select_rows(rows, set(SEED_BLOCKS), set(CORE_ORDER))
    _, exit_code = _run_rows(
        core_rows,
        repository=repository,
        work=work,
        scratch=scratch,
        python=python,
        monitor_interval=monitor_interval,
        continue_on_error=continue_on_error,
    )
    if exit_code:
        return exit_code
    _checked([str(python), str(analysis_script), "--tier", "core"], repository)
    _checked(
        [
            str(python),
            "-m",
            "trophosome.cli",
            "report",
            "--analysis",
            str(core_analysis),
            "--design",
            str(core_design),
            "--output",
            str(report),
            "--markdown",
            str(markdown),
            "--assets",
            str(assets),
            "--title",
            "Phase 1 first-pilot report",
        ],
        repository,
    )
    _checked(
        [
            str(python),
            "scripts/assess_phase1_pilot_expansion.py",
            "--analysis",
            str(core_analysis),
            "--design",
            str(design),
            "--core-report",
            str(report),
            "--output",
            str(safety),
        ],
        repository,
    )
    if core_only:
        print("Core-only workflow complete; extension was not launched.", flush=True)
        return 0

    extension_rows = _select_rows(rows, set(SEED_BLOCKS), set(EXTENSION_ORDER))
    _, exit_code = _run_rows(
        extension_rows,
        repository=repository,
        work=work,
        scratch=scratch,
        python=python,
        monitor_interval=monitor_interval,
        continue_on_error=continue_on_error,
    )
    if exit_code:
        return exit_code
    _checked([str(python), str(analysis_script), "--tier", "all"], repository)
    _checked(
        [
            str(python),
            "-m",
            "trophosome.cli",
            "report",
            "--analysis",
            str(full_analysis),
            "--design",
            str(design),
            "--output",
            str(report),
            "--markdown",
            str(markdown),
            "--assets",
            str(assets),
            "--title",
            "Phase 1 first-pilot report",
        ],
        repository,
    )
    print(f"Full 20-cell workflow complete. Updated report: {report}", flush=True)
    return 0


def _run_one(
    row: dict[str, str],
    repository: Path,
    work: Path,
    scratch: Path,
    python: Path,
    monitor_interval: float,
) -> dict[str, Any]:
    config = work / row["config_path"]
    output = scratch / row["scratch_relative_path"]
    output.mkdir(parents=True, exist_ok=True)
    if _sha256(config) != row["config_sha256"]:
        raise RuntimeError(f"configuration checksum differs: {config}")

    completion = output / "completion.json"
    run_artifacts_exist = any(
        (output / name).exists()
        for name in ("resolved_config.json", "provenance.json", "checkpoints")
    )
    checkpoint_exists = any(
        (output / "checkpoints").glob("checkpoint-rep*-gen*.npz")
    )
    if not completion.exists() and run_artifacts_exist and not checkpoint_exists:
        archive = _archive_unrestartable_attempt(output, scratch, row["run_id"])
        print(
            f"[{row['run_id']}] partial attempt preceded its first checkpoint; "
            f"preserved at {archive} and restarting cleanly",
            flush=True,
        )
        run_artifacts_exist = False
    resume = not completion.exists() and run_artifacts_exist
    if completion.exists():
        status = "already-complete"
        print(f"[{row['run_id']}] {status}", flush=True)
        return {
            "run_id": row["run_id"],
            "status": status,
            "output_bytes": _directory_size(output),
        }

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

    started_at = datetime.now(UTC).isoformat()
    summary_path = output / "execution-summary.json"
    summary: dict[str, Any] = {
        "execution_summary_schema_version": "1.0.0",
        "experiment_id": _experiment_id(row["cell_id"]),
        "run_id": row["run_id"],
        "cell_id": row["cell_id"],
        "seed_block_id": row["seed_block_id"],
        "started_at": started_at,
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
    stdout_path = output / "run.out"
    stderr_path = output / "run.err"
    process: subprocess.Popen[str] | None = None
    try:
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
                time.sleep(monitor_interval)
            return_code = process.returncode
    except KeyboardInterrupt:
        if process is not None and process.poll() is None:
            os.killpg(process.pid, signal.SIGINT)
            try:
                return_code = process.wait(timeout=30)
            except subprocess.TimeoutExpired:
                os.killpg(process.pid, signal.SIGTERM)
                return_code = process.wait(timeout=30)
        else:
            return_code = 130
        summary["status"] = "interrupted"
        raise
    finally:
        elapsed = time.monotonic() - started
        if not memory_measurements:
            child_peak = int(resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss)
            # ru_maxrss is bytes on macOS and KiB on Linux.
            peak_rss_kib = (
                child_peak // 1024 if sys.platform == "darwin" else child_peak
            )
            memory_mode = "cumulative-child-peak"
        else:
            memory_mode = "process-tree-polling"
        summary.update(
            {
                "completed_at": datetime.now(UTC).isoformat(),
                "elapsed_seconds": elapsed,
                "peak_process_tree_rss_kib": (peak_rss_kib if peak_rss_kib else None),
                "memory_measurement_mode": memory_mode,
                "memory_measurements": memory_measurements,
                "output_bytes": _directory_size(output),
            }
        )
        if "return_code" not in summary and "return_code" in locals():
            summary["return_code"] = return_code
        _atomic_json(summary_path, summary)

    successful = return_code == 0 and completion.is_file()
    summary["status"] = "complete" if successful else "failed"
    summary["return_code"] = return_code
    summary["output_bytes"] = _directory_size(output)
    _atomic_json(summary_path, summary)
    print(
        f"[{row['run_id']}] {summary['status']} in {elapsed / 60:.1f} min; "
        f"peak RSS {peak_rss_kib / 1024:.0f} MiB; "
        f"output {summary['output_bytes'] / 1024**2:.1f} MiB",
        flush=True,
    )
    return summary


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--repository",
        type=Path,
        default=Path(__file__).resolve().parents[1],
    )
    parser.add_argument(
        "--seed-block",
        action="append",
        default=None,
        help="selected-run mode: seed block to run; repeat for several",
    )
    parser.add_argument(
        "--cell",
        action="append",
        type=_normalise_cell_id,
        help="selected-run mode: cell to run; repeat for several",
    )
    parser.add_argument(
        "--python",
        type=Path,
        help="Python interpreter containing trophosome (default: repository .venv)",
    )
    parser.add_argument("--monitor-interval", type=float, default=2.0)
    parser.add_argument("--continue-on-error", action="store_true")
    parser.add_argument(
        "--core-only",
        action="store_true",
        help="run, report and assess the core but do not launch the extension",
    )
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args()

    repository = args.repository.resolve()
    # Keep a virtual-environment interpreter path as written.  Resolving its
    # symlink to the base interpreter would bypass the environment's packages.
    python = Path(
        os.path.abspath(args.python or repository / ".venv" / "bin" / "python")
    )
    if not python.is_file():
        parser.error(f"Python interpreter does not exist: {python}")
    scratch, rows = _load_rows(repository)
    selected_mode = bool(args.seed_block or args.cell)
    if not selected_mode and args.dry_run:
        print(
            "Staged workflow: core 12 x 3 seed blocks -> core analysis/report -> "
            "safety gate -> extension 8 x 3 seed blocks -> 20-cell analysis/report"
        )
        return 0
    work = repository / "experiments" / "work" / "trophosome"
    if not selected_mode:
        return _run_full_workflow(
            repository=repository,
            work=work,
            scratch=scratch,
            rows=rows,
            python=python,
            monitor_interval=args.monitor_interval,
            continue_on_error=args.continue_on_error,
            core_only=args.core_only,
        )

    seed_blocks = set(args.seed_block or ["sb0001"])
    cells = set(args.cell) if args.cell else set(CORE_ORDER)
    selected = _select_rows(rows, seed_blocks, cells)
    if not selected:
        parser.error("no runs match the requested cells and seed blocks")
    known_seed_blocks = {row["seed_block_id"] for row in rows}
    unknown = seed_blocks - known_seed_blocks
    if unknown:
        parser.error(f"unknown seed block(s): {', '.join(sorted(unknown))}")
    print(
        f"Selected {len(selected)} population(s): "
        + ", ".join(row["run_id"] for row in selected),
        flush=True,
    )
    if args.dry_run:
        return 0

    results, exit_code = _run_rows(
        selected,
        repository=repository,
        work=work,
        scratch=scratch,
        python=python,
        monitor_interval=args.monitor_interval,
        continue_on_error=args.continue_on_error,
    )

    aggregate_path = (
        scratch
        / "p01-neutral-feedback"
        / "s01-pilot"
        / ("execution-summary-" + "-".join(sorted(seed_blocks)) + ".json")
    )
    aggregate_path.parent.mkdir(parents=True, exist_ok=True)
    _atomic_json(
        aggregate_path,
        {
            "execution_summary_schema_version": "1.0.0",
            "experiment_id": EXPERIMENT_ID,
            "created_at": datetime.now(UTC).isoformat(),
            "seed_blocks": sorted(seed_blocks),
            "results": results,
        },
    )
    print(f"Wrote aggregate summary: {aggregate_path}", flush=True)
    return exit_code


if __name__ == "__main__":
    raise SystemExit(main())
