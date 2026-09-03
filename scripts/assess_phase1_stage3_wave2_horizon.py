#!/usr/bin/env python3
"""Freeze a Wave 2 decision to continue selected trajectories to 500 or 1000."""

from __future__ import annotations

import argparse
import csv
import json
import math
import statistics
from collections import defaultdict
from pathlib import Path

from prepare_phase1_first_pilot import _sha256, _tsv_text
from prepare_phase1_stage3_wave2 import (
    ADAPTIVE_MIGRATION_LEVELS,
    CELLS,
    EXPERIMENT_ID,
    INITIAL_HORIZON,
    INTERMEDIATE_HORIZON,
    MAXIMUM_HORIZON,
    SEED_BLOCKS,
    STAGE_DIRECTORY,
)

T_CRITICAL_90_DF11 = 1.795884819
ABSOLUTE_TV_FLOOR = 0.002
RELATIVE_TV_MARGIN = 0.25
WINDOWS = {
    INITIAL_HORIZON: ((51, 75), (76, 100)),
    INTERMEDIATE_HORIZON: ((401, 450), (451, 500)),
}


def _work_and_scratch(repository: Path) -> tuple[Path, Path]:
    work = repository / "experiments/work/trophosome"
    layout = work / "layout.local.json"
    if not layout.is_file():
        raise RuntimeError("machine-local layout.local.json is missing")
    scratch = Path(json.loads(layout.read_text(encoding="utf-8"))["scratch"])
    return work, scratch


def _new_tv_trajectory(path: Path, initial_counts: tuple[int, ...]) -> dict[int, float]:
    initial_total = sum(initial_counts)
    initial = {
        index: count / initial_total for index, count in enumerate(initial_counts)
    }
    differences: dict[int, float] = defaultdict(lambda: 1.0)
    with path.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            if int(row["replicate"]) != 0:
                raise ValueError(f"unexpected replicate in {path}")
            generation = int(row["generation"])
            strain = int(row["strain_id"])
            frequency = int(row["count"]) / initial_total
            baseline = initial.get(strain, 0.0)
            if strain in initial:
                differences[generation] += abs(frequency - baseline) - baseline
            else:
                differences[generation] += frequency
    return {generation: 0.5 * value for generation, value in differences.items()}


def load_trajectories(
    repository: Path, horizon: int
) -> tuple[dict[tuple[str, str, int], float], list[dict[str, str]]]:
    if horizon not in WINDOWS:
        raise ValueError(f"horizon must be one of {sorted(WINDOWS)}")
    work, scratch = _work_and_scratch(repository)
    phase = work / "p01-neutral-feedback"
    manifest_path = phase / "manifests" / f"{EXPERIMENT_ID}-runs.tsv"
    with manifest_path.open(newline="", encoding="utf-8") as handle:
        runs = list(csv.DictReader(handle, delimiter="\t"))
    initial_payload = json.loads(
        (work / "common/initial-populations/ip001-fisher100.json").read_text(
            encoding="utf-8"
        )
    )
    initial_counts = tuple(int(value) for value in initial_payload["scaled_counts"])
    values: dict[tuple[str, str, int], float] = {}
    fingerprints: list[dict[str, str]] = []

    required_cells = {cell.cell_id for cell in CELLS}
    if horizon == INTERMEDIATE_HORIZON:
        previous = decision_path(work, INITIAL_HORIZON)
        decision = json.loads(previous.read_text(encoding="utf-8"))
        required_cells = set(decision["selected_cell_ids"])
        fingerprints.append(
            {
                "kind": "previous-decision",
                "path": str(previous.relative_to(work)),
                "sha256": _sha256(previous.read_text(encoding="utf-8")),
            }
        )

    reference_path = phase / "design" / f"{EXPERIMENT_ID}-reused-trajectories.tsv"
    if horizon == INITIAL_HORIZON:
        with reference_path.open(newline="", encoding="utf-8") as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                generation = int(row["generation"])
                values[(row["cell_id"], row["seed_block_id"], generation)] = float(
                    row["TV"]
                )
        fingerprints.append(
            {
                "kind": "frozen-reuse",
                "path": str(reference_path.relative_to(work)),
                "sha256": _sha256(reference_path.read_text(encoding="utf-8")),
            }
        )

    for row in runs:
        if row["cell_id"] not in required_cells:
            continue
        output = scratch / row["scratch_relative_path"]
        state_path = output / (
            "completion.json"
            if (output / "completion.json").is_file()
            else "pause.json"
        )
        if not state_path.is_file():
            raise RuntimeError(f"missing paused/completed state for {row['run_id']}")
        state = json.loads(state_path.read_text(encoding="utf-8"))
        observed_horizon = (
            MAXIMUM_HORIZON
            if state_path.name == "completion.json"
            else int(state["last_completed_generation"])
        )
        if observed_horizon < horizon:
            raise RuntimeError(
                f"{row['run_id']} has reached {observed_horizon}, not {horizon}"
            )
        environment_path = output / "environment_counts.csv"
        trajectory = _new_tv_trajectory(environment_path, initial_counts)
        expected_generations = set(range(horizon + 1))
        if not expected_generations <= set(trajectory):
            raise RuntimeError(f"incomplete environmental trajectory: {row['run_id']}")
        for generation in expected_generations:
            values[(row["cell_id"], row["seed_block_id"], generation)] = trajectory[
                generation
            ]
        # Fingerprint only the assessed prefix.  A later authorized resume
        # appends rows and replaces pause.json, but it cannot change passages
        # already committed in the checkpointed trajectory.
        prefix = json.dumps(
            [
                format(trajectory[generation], ".17g")
                for generation in range(horizon + 1)
            ],
            separators=(",", ":"),
        )
        fingerprints.append(
            {
                "kind": "TV-trajectory-prefix",
                "run_id": row["run_id"],
                "through_generation": str(horizon),
                "sha256": _sha256(prefix),
            }
        )

    expected_keys = {
        (cell_id, seed, generation)
        for cell_id in required_cells
        for seed, _master in SEED_BLOCKS
        for generation in range(horizon + 1)
    }
    missing = expected_keys - set(values)
    if missing:
        preview = ", ".join(map(str, sorted(missing)[:3]))
        raise RuntimeError(f"missing {len(missing)} trajectory values ({preview})")
    extra_cells = {key[0] for key in values} - required_cells
    if extra_cells:
        for key in [key for key in values if key[0] in extra_cells]:
            del values[key]
    return values, fingerprints


def _window_change(
    values: dict[tuple[str, str, int], float],
    cell_id: str,
    horizon: int,
    *,
    control_id: str | None = None,
) -> dict[str, object]:
    first, second = WINDOWS[horizon]
    changes: list[float] = []
    late_means: list[float] = []
    for seed, _master in SEED_BLOCKS:

        def observed(generation: int, selected_seed: str = seed) -> float:
            value = values[(cell_id, selected_seed, generation)]
            if control_id is not None:
                value -= values[(control_id, selected_seed, generation)]
            return value

        early = statistics.fmean(observed(g) for g in range(first[0], first[1] + 1))
        late = statistics.fmean(observed(g) for g in range(second[0], second[1] + 1))
        changes.append(late - early)
        late_means.append(late)
    mean_change = statistics.fmean(changes)
    se = statistics.stdev(changes) / math.sqrt(len(changes))
    lower = mean_change - T_CRITICAL_90_DF11 * se
    upper = mean_change + T_CRITICAL_90_DF11 * se
    late_mean = statistics.fmean(late_means)
    margin = max(ABSOLUTE_TV_FLOOR, RELATIVE_TV_MARGIN * abs(late_mean))
    return {
        "mean_change": mean_change,
        "ci90_lower": lower,
        "ci90_upper": upper,
        "late_mean": late_mean,
        "margin": margin,
        "stable": lower >= -margin and upper <= margin,
    }


def build_decision(
    repository: Path, horizon: int
) -> tuple[dict, list[dict[str, object]]]:
    values, fingerprints = load_trajectories(repository, horizon)
    by_alpha_m = {
        (cell.alpha_target, cell.migration_fraction): cell
        for cell in CELLS
        if cell.panel == "alpha-by-m"
    }
    eligible = [
        cell
        for cell in CELLS
        if cell.panel == "alpha-by-m"
        and cell.migration_fraction in ADAPTIVE_MIGRATION_LEVELS
    ]
    if horizon == INTERMEDIATE_HORIZON:
        work, _scratch = _work_and_scratch(repository)
        previous = json.loads(decision_path(work, INITIAL_HORIZON).read_text())
        eligible = [
            cell for cell in eligible if cell.cell_id in previous["selected_cell_ids"]
        ]

    diagnostics: list[dict[str, object]] = []
    diagnostics_by_cell: dict[str, dict[str, dict[str, object]]] = defaultdict(dict)
    for cell in eligible:
        raw = _window_change(values, cell.cell_id, horizon)
        diagnostics_by_cell[cell.cell_id]["raw-TV"] = raw
        diagnostics.append(
            {
                "cell_id": cell.cell_id,
                "cell": cell.short_id,
                "alpha": cell.alpha_target,
                "m": cell.migration_fraction,
                "diagnostic": "raw-TV",
                **raw,
            }
        )
        if cell.alpha_target != "0":
            control = by_alpha_m[("0", cell.migration_fraction)]
            effect = _window_change(
                values, cell.cell_id, horizon, control_id=control.cell_id
            )
            diagnostics_by_cell[cell.cell_id]["host-induced-TV"] = effect
            diagnostics.append(
                {
                    "cell_id": cell.cell_id,
                    "cell": cell.short_id,
                    "alpha": cell.alpha_target,
                    "m": cell.migration_fraction,
                    "diagnostic": "host-induced-TV",
                    **effect,
                }
            )

    if horizon == INITIAL_HORIZON:
        target = INTERMEDIATE_HORIZON
        mandatory = {
            by_alpha_m[(alpha, migration)].cell_id
            for alpha in ("0", "0.1")
            for migration in ADAPTIVE_MIGRATION_LEVELS
        }
    else:
        target = MAXIMUM_HORIZON
        mandatory = {
            by_alpha_m[(alpha, migration)].cell_id
            for alpha in ("0", "0.1")
            for migration in ("0", "0.001")
        }

    selected = set(mandatory)
    reasons: dict[str, list[str]] = defaultdict(list)
    for cell_id in mandatory:
        reasons[cell_id].append("pre-specified-anchor")

    migrations_with_unstable_controls = {
        cell.migration_fraction
        for cell in eligible
        if cell.alpha_target == "0"
        and not diagnostics_by_cell[cell.cell_id]["raw-TV"]["stable"]
    }
    for cell in eligible:
        if cell.alpha_target == "0":
            continue
        diagnostics_for_cell = diagnostics_by_cell[cell.cell_id]
        failures = [
            name
            for name, record in diagnostics_for_cell.items()
            if not record["stable"]
        ]
        if failures or cell.migration_fraction in migrations_with_unstable_controls:
            selected.add(cell.cell_id)
            control = by_alpha_m[("0", cell.migration_fraction)]
            selected.add(control.cell_id)
            for name in failures:
                reasons[cell.cell_id].append(f"unresolved-{name}")
            if cell.migration_fraction in migrations_with_unstable_controls:
                reasons[cell.cell_id].append("unresolved-same-m-control")
                reasons[control.cell_id].append("unresolved-raw-TV")

    selected_cells = sorted(
        (cell for cell in CELLS if cell.cell_id in selected),
        key=lambda cell: cell.number,
    )
    decision = {
        "adaptive_horizon_decision_schema_version": "1.0.0",
        "experiment_id": EXPERIMENT_ID,
        "status": "frozen-outcome-dependent-decision",
        "assessed_horizon": horizon,
        "continuation_horizon": target,
        "selected_cell_ids": [cell.cell_id for cell in selected_cells],
        "selected_cells": [cell.short_id for cell in selected_cells],
        "selected_populations": len(selected_cells) * len(SEED_BLOCKS),
        "reasons": {cell.cell_id: reasons[cell.cell_id] for cell in selected_cells},
        "windows": [list(window) for window in WINDOWS[horizon]],
        "stability_margin": (
            f"max({ABSOLUTE_TV_FLOOR} TV, "
            f"{RELATIVE_TV_MARGIN} * absolute late-window mean)"
        ),
        "uncertainty": "paired 90% Student-t interval across 12 seed blocks",
        "input_fingerprints": sorted(
            fingerprints, key=lambda row: (row.get("kind", ""), row.get("run_id", ""))
        ),
        "interpretation": (
            "This decision controls exploratory time-horizon diagnostics only; "
            "passage 100 remains the primary Wave 2 endpoint."
        ),
    }
    return decision, diagnostics


def analysis_directory(work: Path) -> Path:
    return work / "p01-neutral-feedback" / "analysis" / f"{STAGE_DIRECTORY}-derived"


def decision_path(work: Path, horizon: int) -> Path:
    return analysis_directory(work) / f"adaptive-horizon-decision-g{horizon}.json"


def diagnostic_path(work: Path, horizon: int) -> Path:
    return analysis_directory(work) / f"adaptive-horizon-diagnostics-g{horizon}.tsv"


def write_or_verify(repository: Path, horizon: int, *, verify: bool) -> None:
    work, _scratch = _work_and_scratch(repository)
    decision, diagnostics = build_decision(repository, horizon)
    decision_text = json.dumps(decision, indent=2) + "\n"
    diagnostic_text = _tsv_text(diagnostics, list(diagnostics[0]))
    outputs = {
        decision_path(work, horizon): decision_text,
        diagnostic_path(work, horizon): diagnostic_text,
    }
    mismatches = [
        str(path)
        for path, content in outputs.items()
        if not path.is_file() or path.read_text(encoding="utf-8") != content
    ]
    if verify:
        if mismatches:
            raise RuntimeError(
                "adaptive decision differs or is missing: " + ", ".join(mismatches)
            )
        return
    existing = [path for path in outputs if path.is_file()]
    if existing and mismatches:
        raise RuntimeError(
            "an existing adaptive decision is immutable but no longer matches "
            "its inputs"
        )
    for path, content in outputs.items():
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(content, encoding="utf-8")
    print(
        f"Frozen passage-{horizon} decision: continue "
        f"{len(decision['selected_cell_ids'])} cells "
        f"({decision['selected_populations']} populations) to passage "
        f"{decision['continuation_horizon']}."
    )


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--repository", type=Path, default=Path(__file__).resolve().parents[1]
    )
    parser.add_argument("--horizon", type=int, choices=sorted(WINDOWS), required=True)
    parser.add_argument("--verify", action="store_true")
    args = parser.parse_args()
    write_or_verify(args.repository.resolve(), args.horizon, verify=args.verify)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
