#!/usr/bin/env python3
"""Audit 288 new populations plus reused references and report; never simulate."""

from __future__ import annotations

import argparse
import hashlib
from datetime import UTC, datetime
from pathlib import Path

from analyse_phase1_stage3_wave1 import analyse
from prepare_phase1_stage3_wave1 import EXPERIMENT_ID, STAGE_DIRECTORY
from run_phase1_first_pilot import _atomic_json, _sha256

from trophosome.stage3_report import generate_report


def build(repository: Path) -> list[Path]:
    repository = repository.resolve()
    phase = repository / "experiments/work/trophosome/p01-neutral-feedback"
    derived = phase / "analysis" / f"{STAGE_DIRECTORY}-derived"
    derived.mkdir(parents=True, exist_ok=True)
    completion = derived / "report-completion.json"
    # Never treat a stale report as current after a failed rebuild.
    _atomic_json(
        completion,
        {"complete": False, "status": "auditing", "experiment_id": EXPERIMENT_ID},
    )
    analyse(repository)
    artifacts = generate_report(
        analysis=derived,
        design=phase / "design" / f"{EXPERIMENT_ID}-cells.tsv",
        pdf=repository / "output/pdf" / f"{EXPERIMENT_ID}-report.pdf",
        markdown=repository / "docs" / f"{EXPERIMENT_ID}-report.md",
        assets=repository / "docs/figures" / f"{EXPERIMENT_ID}-report",
    )
    inputs = [
        phase / "manifests" / f"{EXPERIMENT_ID}-manifest.json",
        phase / "design" / f"{EXPERIMENT_ID}-reference-endpoints.tsv",
        repository / "scripts/analyse_phase1_stage3_wave1.py",
        repository / "scripts/prepare_phase1_stage3_wave1.py",
        repository / "scripts/run_phase1_stage3_wave1.py",
        phase / "analysis/analyse_first_pilot.py",
        phase / "analysis/analyse_second_pilot.py",
        repository / "src/trophosome/stage3_report.py",
        repository / "src/trophosome/second_pilot_report.py",
        Path(__file__).resolve(),
    ]
    artifacts.extend(sorted(derived.glob("*.tsv")))
    artifacts.extend(
        [derived / "analysis-summary.json", derived / "analysis-audit.json"]
    )
    digest = hashlib.sha256()
    for path in inputs + artifacts:
        digest.update(str(path.relative_to(repository)).encode())
        digest.update(_sha256(path).encode())
    _atomic_json(
        completion,
        {
            "complete": True,
            "experiment_id": EXPERIMENT_ID,
            "completed_at": datetime.now(UTC).isoformat(),
            "analysis_audit": "PASS",
            "input_and_output_sha256": digest.hexdigest(),
            "inputs": [
                {"path": str(path.relative_to(repository)), "sha256": _sha256(path)}
                for path in inputs
            ],
            "outputs": [
                {"path": str(path.relative_to(repository)), "sha256": _sha256(path)}
                for path in artifacts
            ],
        },
    )
    print(f"Report complete: {artifacts[0]}")
    print(f"Editable report: {artifacts[1]}")
    return artifacts


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--repository", type=Path, default=Path(__file__).resolve().parents[1]
    )
    args = parser.parse_args()
    build(args.repository.resolve())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
