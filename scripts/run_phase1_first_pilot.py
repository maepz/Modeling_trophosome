#!/usr/bin/env python3
"""Run selected Phase 1 first-pilot populations and record resource use.

The launcher is deliberately scheduler-independent.  It runs populations
sequentially by default, which is appropriate for the first Mac pilot.  On an
HPC, separate invocations can safely run different cell/seed-block pairs.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import signal
import subprocess
import time
from datetime import UTC, datetime
from pathlib import Path
from typing import Any

EXPERIMENT_ID = "phase1-first-pilot-core12"
PILOT_ORDER = (
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
    raise argparse.ArgumentTypeError(
        "cell must look like c0001, 1, or p01-s01-c0001"
    )


def _load_rows(repository: Path) -> tuple[Path, list[dict[str, str]]]:
    work = repository / "experiments" / "work" / "trophosome"
    layout = json.loads((work / "layout.local.json").read_text(encoding="utf-8"))
    scratch = Path(layout["scratch"])
    manifest = (
        work
        / "p01-neutral-feedback"
        / "manifests"
        / "phase1-first-pilot-runs.tsv"
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
        "experiment_id": EXPERIMENT_ID,
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
        with stdout_path.open("a", encoding="utf-8") as stdout, stderr_path.open(
            "a", encoding="utf-8"
        ) as stderr:
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
        summary.update(
            {
                "completed_at": datetime.now(UTC).isoformat(),
                "elapsed_seconds": elapsed,
                "peak_process_tree_rss_kib": (
                    peak_rss_kib if memory_measurements else None
                ),
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
        help="seed block to run; repeat for several (default: sb0001)",
    )
    parser.add_argument(
        "--cell",
        action="append",
        type=_normalise_cell_id,
        help="cell to run; repeat for several (default: all 12)",
    )
    parser.add_argument(
        "--python",
        type=Path,
        help="Python interpreter containing trophosome (default: repository .venv)",
    )
    parser.add_argument("--monitor-interval", type=float, default=2.0)
    parser.add_argument("--continue-on-error", action="store_true")
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
    seed_blocks = set(args.seed_block or ["sb0001"])
    cells = set(args.cell) if args.cell else None
    scratch, rows = _load_rows(repository)
    selected = _select_rows(rows, seed_blocks, cells)
    if not selected:
        parser.error("no runs match the requested cells and seed blocks")
    known_seed_blocks = {row["seed_block_id"] for row in rows}
    unknown = seed_blocks - known_seed_blocks
    if unknown:
        parser.error(f"unknown seed block(s): {', '.join(sorted(unknown))}")
    work = repository / "experiments" / "work" / "trophosome"

    print(
        f"Selected {len(selected)} population(s): "
        + ", ".join(row["run_id"] for row in selected),
        flush=True,
    )
    if args.dry_run:
        return 0

    results = []
    exit_code = 0
    for row in selected:
        result = _run_one(
            row,
            repository,
            work,
            scratch,
            python,
            args.monitor_interval,
        )
        results.append(result)
        if result["status"] == "failed":
            exit_code = 1
            if not args.continue_on_error:
                break

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
