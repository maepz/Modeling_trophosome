#!/usr/bin/env python3
"""Audit, analyse, and build the completed Phase 1 second-pilot report."""

from __future__ import annotations

import argparse
import csv
import hashlib
import importlib.util
import json
import subprocess
import sys
from datetime import UTC, datetime
from pathlib import Path
from typing import Any

VARIANT_TAG = "v210-m010-g250"
EXPERIMENT_ID = f"phase1-second-pilot-{VARIANT_TAG}"
STAGE_DIRECTORY = f"s02-equilibrium-precision-{VARIANT_TAG}"


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _atomic_json(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    temporary.replace(path)


def _paths(repository: Path) -> dict[str, Path]:
    work = repository / "experiments/work/trophosome"
    phase = work / "p01-neutral-feedback"
    derived = phase / "analysis" / f"{STAGE_DIRECTORY}-derived"
    return {
        "work": work,
        "manifest": phase / "manifests" / f"{EXPERIMENT_ID}-runs.tsv",
        "design": phase / "design" / f"{EXPERIMENT_ID}-cells.tsv",
        "analysis_script": phase / "analysis" / "analyse_second_pilot.py",
        "derived": derived,
        "pdf": repository / "output/pdf" / f"{EXPERIMENT_ID}-report.pdf",
        "markdown": repository / "docs" / f"{EXPERIMENT_ID}-report.md",
        "assets": repository / "docs/figures" / f"{EXPERIMENT_ID}-report",
        "completion": derived / "report-completion.json",
    }


def _input_fingerprint(repository: Path, paths: dict[str, Path]) -> str:
    layout_path = paths["work"] / "layout.local.json"
    layout = json.loads(layout_path.read_text(encoding="utf-8"))
    scratch = Path(layout["scratch"])
    with paths["manifest"].open(newline="", encoding="utf-8") as handle:
        runs = list(csv.DictReader(handle, delimiter="\t"))
    digest = hashlib.sha256()
    source_paths = (
        paths["manifest"],
        paths["design"],
        paths["analysis_script"],
        repository / "src/trophosome/second_pilot_report.py",
        Path(__file__).resolve(),
    )
    for path in source_paths:
        digest.update(str(path.relative_to(repository)).encode())
        digest.update(b"\0")
        digest.update(path.read_bytes())
        digest.update(b"\0")
    for run in sorted(runs, key=lambda row: row["run_id"]):
        completion = scratch / run["scratch_relative_path"] / "completion.json"
        digest.update(run["run_id"].encode())
        digest.update(b"\0")
        digest.update(completion.read_bytes() if completion.is_file() else b"MISSING")
        digest.update(b"\0")
    return digest.hexdigest()


def _completion_is_current(
    completion_path: Path, fingerprint: str
) -> tuple[bool, dict[str, Any] | None]:
    if not completion_path.is_file():
        return False, None
    try:
        payload = json.loads(completion_path.read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return False, None
    if payload.get("complete") is not True:
        return False, payload
    if payload.get("input_fingerprint_sha256") != fingerprint:
        return False, payload
    for item in payload.get("outputs", []):
        path = Path(item["path"])
        if not path.is_file() or _sha256(path) != item["sha256"]:
            return False, payload
    return True, payload


def _check_dependencies() -> None:
    missing = [
        name
        for name in ("matplotlib", "numpy", "reportlab")
        if importlib.util.find_spec(name) is None
    ]
    if missing:
        raise RuntimeError(
            "reporting dependencies are missing: "
            + ", ".join(missing)
            + ". Install the trophosome report extra in the active environment."
        )


def _run(command: list[str], repository: Path) -> None:
    print("Report step: " + " ".join(command), flush=True)
    subprocess.run(command, cwd=repository, check=True)


def main() -> int:
    parser = argparse.ArgumentParser(
        description=__doc__,
        epilog=(
            "This is the stand-alone report path: it never launches simulations. "
            "It stops unless all 72 committed populations pass the audit."
        ),
    )
    parser.add_argument(
        "--repository", type=Path, default=Path(__file__).resolve().parents[1]
    )
    parser.add_argument("--report-date", metavar="YYYY-MM-DD")
    parser.add_argument(
        "--force", action="store_true", help="rebuild even when inputs are unchanged"
    )
    args = parser.parse_args()
    repository = args.repository.resolve()
    paths = _paths(repository)
    fingerprint = _input_fingerprint(repository, paths)
    current, payload = _completion_is_current(paths["completion"], fingerprint)
    if current and not args.force:
        print(f"Report is already current: {payload['pdf_path']}", flush=True)
        return 0

    _check_dependencies()
    _run(
        [
            sys.executable,
            str(paths["analysis_script"]),
            "--repository",
            str(repository),
        ],
        repository,
    )
    summary = json.loads(
        (paths["derived"] / "analysis-summary.json").read_text(encoding="utf-8")
    )
    if summary.get("audit_status") != "PASS":
        raise RuntimeError("analysis audit did not pass; refusing to build the report")

    from trophosome.second_pilot_report import generate_second_pilot_report

    generate_second_pilot_report(
        analysis=paths["derived"],
        design=paths["design"],
        output=paths["pdf"],
        markdown=paths["markdown"],
        assets=paths["assets"],
        report_date=args.report_date,
    )
    required = (
        paths["pdf"],
        paths["markdown"],
        paths["derived"] / "analysis-summary.json",
        paths["derived"] / "analysis-audit.json",
        paths["derived"] / "environment-trajectories.tsv",
        paths["derived"] / "run-window-summaries.tsv",
        paths["derived"] / "stationarity-screen.tsv",
        paths["derived"] / "stationarity-history.tsv",
        paths["derived"] / "precision-recommendations.tsv",
        paths["assets"] / "environmental-trajectories.png",
        paths["assets"] / "stationarity-screen.png",
        paths["assets"] / "continuing-fluctuations.png",
        paths["assets"] / "precision-recommendations.png",
    )
    missing = [str(path) for path in required if not path.is_file()]
    if missing:
        raise RuntimeError("report output is missing: " + ", ".join(missing))
    if not paths["pdf"].read_bytes().startswith(b"%PDF-"):
        raise RuntimeError("the report does not have a valid PDF header")
    markdown = paths["markdown"].read_text(encoding="utf-8").lower()
    for required_text in (
        "stationarity",
        "not definitive equilibrium",
        "audit: `pass`",
    ):
        if required_text not in markdown:
            raise RuntimeError(f"report is missing required text: {required_text}")

    completion = {
        "report_completion_schema_version": "1.0.0",
        "complete": True,
        "completed_at": datetime.now(UTC).isoformat(),
        "experiment_id": EXPERIMENT_ID,
        "variant_id": VARIANT_TAG,
        "input_fingerprint_sha256": fingerprint,
        "analysis_audit": summary["audit_status"],
        "pdf_path": str(paths["pdf"]),
        "markdown_path": str(paths["markdown"]),
        "assets_path": str(paths["assets"]),
        "outputs": [
            {
                "path": str(path),
                "sha256": _sha256(path),
                "bytes": path.stat().st_size,
            }
            for path in required
        ],
    }
    _atomic_json(paths["completion"], completion)
    print(f"Built self-contained report: {paths['pdf']}", flush=True)
    print(f"Editable companion: {paths['markdown']}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
