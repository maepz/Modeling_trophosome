#!/usr/bin/env python3
"""Create or verify machine-local scratch folders for Phase 1 pilot runs."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path


def expected_runs(
    repository: Path,
) -> tuple[Path, list[tuple[Path, dict[str, object]]]]:
    work = repository / "experiments" / "work" / "trophosome"
    layout = json.loads((work / "layout.local.json").read_text(encoding="utf-8"))
    scratch = Path(layout["scratch"])
    runs_path = (
        work
        / "p01-neutral-feedback"
        / "manifests"
        / "phase1-first-pilot-runs.tsv"
    )
    with runs_path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    expected = []
    for row in rows:
        run_directory = scratch / row["scratch_relative_path"]
        metadata: dict[str, object] = {
            "scratch_manifest_schema_version": "1.0.0",
            "experiment_id": "phase1-first-pilot-core12",
            "run_id": row["run_id"],
            "cell_id": row["cell_id"],
            "seed_block_id": row["seed_block_id"],
            "master_seed": int(row["master_seed"]),
            "within_run_replicate_index": int(
                row["within_run_replicate_index"]
            ),
            "config_path": str(work / row["config_path"]),
            "config_sha256": row["config_sha256"],
            "output_path": str(run_directory),
            "status": "prepared-not-launched",
        }
        expected.append((run_directory, metadata))
    return scratch, expected


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--repository",
        type=Path,
        default=Path(__file__).resolve().parents[1],
    )
    parser.add_argument("--write", action="store_true")
    parser.add_argument("--verify", action="store_true")
    args = parser.parse_args()
    if args.write == args.verify:
        parser.error("choose exactly one of --write or --verify")
    scratch, runs = expected_runs(args.repository.resolve())
    if args.write:
        for run_directory, metadata in runs:
            run_directory.mkdir(parents=True, exist_ok=True)
            path = run_directory / "run.json"
            content = json.dumps(metadata, indent=2, sort_keys=True) + "\n"
            if path.exists() and path.read_text(encoding="utf-8") != content:
                raise SystemExit(f"refusing to replace different run manifest: {path}")
            path.write_text(content, encoding="utf-8")
        print(f"Prepared {len(runs)} run directories below {scratch}")
        return 0
    mismatches = []
    for run_directory, metadata in runs:
        path = run_directory / "run.json"
        content = json.dumps(metadata, indent=2, sort_keys=True) + "\n"
        if not path.is_file() or path.read_text(encoding="utf-8") != content:
            mismatches.append(str(path))
    if mismatches:
        raise SystemExit("scratch manifests differ:\n" + "\n".join(mismatches))
    print(f"Verified {len(runs)} run directories below {scratch}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
