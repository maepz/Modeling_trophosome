#!/usr/bin/env python3
"""Audit Wave 2 passage-100 outputs and compile portable analysis tables."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import statistics
import sys
from collections import defaultdict
from collections.abc import Iterable
from pathlib import Path
from typing import Any

import numpy as np
from prepare_phase1_stage3_wave2 import (
    ALPHA_LEVELS,
    CELLS,
    EXPECTED_RUNS,
    EXPERIMENT_ID,
    INITIAL_HORIZON,
    MIGRATION_LEVELS,
    REUSED_CELLS,
    SEED_BLOCKS,
    STAGE_DIRECTORY,
    K,
    verify_files,
)
from run_phase1_first_pilot import _atomic_json, _sha256
from run_phase1_stage3_wave2 import (
    _state_generation,
    _verify_manifest,
    load_rows,
    state_issues,
)

ANALYSIS = (
    Path(__file__).resolve().parents[1]
    / "experiments/work/trophosome/p01-neutral-feedback/analysis"
)
sys.path.insert(0, str(ANALYSIS))
from analyse_first_pilot import (  # noqa: E402
    composition_metrics,
    diversity_metrics,
)

T_CRITICAL_90_DF11 = 1.795884819
TAIL_START = 51
EARLY_WINDOW = (51, 75)
LATE_WINDOW = (76, 100)
ENDPOINT_RESPONSES = ("D0", "D1", "D2", "evenness", "TV")
ANALYSIS_RESPONSES = (
    "D0_g100",
    "D1_g100",
    "D2_g100",
    "evenness_g100",
    "TV_g100",
    "TV_tail_mean",
    "TV_tail_sd",
)
RELATIVE_RESPONSES = {"D0_g100", "D1_g100", "D2_g100"}


def _read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def _atomic_tsv(
    path: Path, rows: list[dict[str, Any]], fields: Iterable[str]
) -> None:
    if not rows:
        raise ValueError(f"cannot write an empty analysis table: {path.name}")
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    with temporary.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=list(fields),
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(rows)
    temporary.replace(path)


def _text_sha256(text: str) -> str:
    return hashlib.sha256(text.encode("utf-8")).hexdigest()


def _mean_ci90(values: list[float]) -> tuple[float, float, float, float]:
    if len(values) != len(SEED_BLOCKS) or not all(
        math.isfinite(value) for value in values
    ):
        raise ValueError("summary requires 12 finite matched-seed values")
    mean = statistics.fmean(values)
    sd = statistics.stdev(values)
    half_width = T_CRITICAL_90_DF11 * sd / math.sqrt(len(values))
    return mean, sd, mean - half_width, mean + half_width


def _read_environment_prefix(
    path: Path, initial: dict[int, int]
) -> tuple[dict[int, dict[int, int]], str]:
    states: dict[int, dict[int, int]] = defaultdict(dict)
    with path.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            if int(row["replicate"]) != 0:
                raise ValueError("environment table contains an unexpected replicate")
            generation = int(row["generation"])
            if generation > INITIAL_HORIZON:
                continue
            strain = int(row["strain_id"])
            count = int(row["count"])
            if strain < 0 or count <= 0:
                raise ValueError("environment table contains an invalid strain count")
            if strain in states[generation]:
                raise ValueError("duplicate environment generation/strain row")
            states[generation][strain] = count
    expected = set(range(INITIAL_HORIZON + 1))
    if set(states) != expected:
        raise ValueError("environment table lacks part of the passage-100 prefix")
    if states[0] != initial:
        raise ValueError("initial environment differs from the frozen population")
    if any(sum(counts.values()) != K for counts in states.values()):
        raise ValueError("environmental capacity differs from one billion cells")
    if any(set(counts) - set(initial) for counts in states.values()):
        raise ValueError("mutation labels occur in a mutation-free Wave 2 run")
    canonical = json.dumps(
        [
            [generation, strain, count]
            for generation in range(INITIAL_HORIZON + 1)
            for strain, count in sorted(states[generation].items())
        ],
        separators=(",", ":"),
    )
    return dict(states), _text_sha256(canonical)


def _read_host_summary_prefix(
    path: Path,
) -> tuple[dict[int, dict[str, str]], str]:
    rows: dict[int, dict[str, str]] = {}
    with path.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            if int(row["replicate"]) != 0:
                raise ValueError("host summary contains an unexpected replicate")
            generation = int(row["host_generation"])
            if generation > INITIAL_HORIZON:
                continue
            if generation in rows:
                raise ValueError("duplicate host-generation summary row")
            rows[generation] = row
    if set(rows) != set(range(1, INITIAL_HORIZON + 1)):
        raise ValueError("host summary lacks part of the passage-100 prefix")
    canonical = json.dumps(
        [rows[generation] for generation in range(1, INITIAL_HORIZON + 1)],
        sort_keys=True,
        separators=(",", ":"),
    )
    return rows, _text_sha256(canonical)


def _base_row(
    *,
    run_id: str,
    cell: dict[str, str],
    seed: str,
    generation: int,
    source_role: str,
    source_run_id: str,
) -> dict[str, Any]:
    return {
        "run_id": run_id,
        "cell_id": cell["cell_id"],
        "cell": cell["cell"],
        "seed_block_id": seed,
        "panel": cell["panel"],
        "H": int(cell["H"]),
        "B": int(cell["B"]),
        "HB": int(cell["H"]) * int(cell["B"]),
        "alpha_target": float(cell["alpha_target"]),
        "alpha": float(cell["alpha"]),
        "m": float(cell["m"]),
        "generation": generation,
        "source_role": source_role,
        "source_run_id": source_run_id,
    }


def _new_run_rows(
    run: dict[str, str],
    cell: dict[str, str],
    output: Path,
    initial: dict[int, int],
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    states, environment_sha = _read_environment_prefix(
        output / "environment_counts.csv", initial
    )
    host_summaries, host_sha = _read_host_summary_prefix(
        output / "host_generation_summary.csv"
    )
    rows: list[dict[str, Any]] = []
    previous = states[0]
    for generation in range(INITIAL_HORIZON + 1):
        counts = states[generation]
        diversity = diversity_metrics(counts)
        row = _base_row(
            run_id=run["run_id"],
            cell=cell,
            seed=run["seed_block_id"],
            generation=generation,
            source_role="new-simulation",
            source_run_id=run["run_id"],
        )
        row.update(
            {
                "D0": diversity["richness"],
                "shannon": diversity["shannon"],
                "simpson": diversity["simpson"],
                "D1": diversity["hill_q1"],
                "D2": diversity["hill_q2"],
                "evenness": diversity["evenness"],
                "TV": composition_metrics(counts, initial)["total_variation"],
                "turnover": 0.0
                if generation == 0
                else composition_metrics(counts, previous)["total_variation"],
                "realized_host_feedback": "",
                "mean_adult_richness": "",
                "mean_adult_gene_diversity": "",
            }
        )
        if generation:
            summary = host_summaries[generation]
            row.update(
                {
                    "realized_host_feedback": float(
                        summary["realized_host_feedback"]
                    ),
                    "mean_adult_richness": float(summary["mean_adult_richness"]),
                    "mean_adult_gene_diversity": float(
                        summary["mean_adult_gene_diversity"]
                    ),
                }
            )
        rows.append(row)
        previous = counts
    return rows, {
        "run_id": run["run_id"],
        "cell_id": run["cell_id"],
        "seed_block_id": run["seed_block_id"],
        "source_role": "new-simulation",
        "source_run_id": run["run_id"],
        "reached_generation": _state_generation(output),
        "environment_prefix_sha256": environment_sha,
        "host_summary_prefix_sha256": host_sha,
    }


def _reused_rows(
    phase: Path, cells: dict[str, dict[str, str]]
) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    path = phase / "design" / f"{EXPERIMENT_ID}-reused-trajectories.tsv"
    frozen_sha = _sha256(path)
    source = _read_tsv(path)
    rows: list[dict[str, Any]] = []
    inputs: dict[tuple[str, str], dict[str, Any]] = {}
    for record in source:
        cell = cells[record["cell_id"]]
        seed = record["seed_block_id"]
        generation = int(record["generation"])
        d1 = float(record["D1"])
        d2 = float(record["D2"])
        run_id = f"{record['cell_id']}-v210-adaptive-g1000-{seed}-reused"
        row = _base_row(
            run_id=run_id,
            cell=cell,
            seed=seed,
            generation=generation,
            source_role="reused-passage-100",
            source_run_id=record["source_run_id"],
        )
        key = (record["cell_id"], seed)
        tv = float(record["TV"])
        row.update(
            {
                "D0": int(float(record["D0"])),
                "shannon": math.log(d1),
                "simpson": 1.0 - 1.0 / d2,
                "D1": d1,
                "D2": d2,
                "evenness": float(record["evenness"]),
                "TV": tv,
                "turnover": "",
                "realized_host_feedback": ""
                if generation == 0
                else float(cell["alpha"]),
                "mean_adult_richness": "",
                "mean_adult_gene_diversity": "",
            }
        )
        rows.append(row)
        inputs.setdefault(
            key,
            {
                "run_id": run_id,
                "cell_id": record["cell_id"],
                "seed_block_id": seed,
                "source_role": "reused-passage-100",
                "source_run_id": record["source_run_id"],
                "reached_generation": INITIAL_HORIZON,
                "environment_prefix_sha256": frozen_sha,
                "host_summary_prefix_sha256": "not-frozen-for-reuse",
            },
        )
    expected = len(REUSED_CELLS) * len(SEED_BLOCKS) * (INITIAL_HORIZON + 1)
    keys = {
        (row["cell_id"], row["seed_block_id"], row["generation"]) for row in rows
    }
    if len(rows) != expected or len(keys) != expected:
        raise ValueError(
            f"expected {expected} unique reused trajectory rows, found "
            f"{len(rows)} ({len(keys)} unique)"
        )
    return rows, list(inputs.values())


def _run_summaries(
    trajectories: list[dict[str, Any]],
) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    grouped: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for row in trajectories:
        grouped[row["run_id"]].append(row)
    endpoints: list[dict[str, Any]] = []
    tails: list[dict[str, Any]] = []
    for run_id, values in grouped.items():
        values.sort(key=lambda row: row["generation"])
        if [row["generation"] for row in values] != list(
            range(INITIAL_HORIZON + 1)
        ):
            raise ValueError(f"{run_id}: incomplete or duplicated trajectory")
        endpoint = dict(values[-1])
        endpoints.append(endpoint)
        tail = {
            key: endpoint[key]
            for key in (
                "run_id",
                "cell_id",
                "cell",
                "seed_block_id",
                "panel",
                "H",
                "B",
                "HB",
                "alpha_target",
                "alpha",
                "m",
                "source_role",
                "source_run_id",
            )
        }
        selected = [
            row for row in values if TAIL_START <= row["generation"] <= 100
        ]
        for response in ENDPOINT_RESPONSES:
            observed = [float(row[response]) for row in selected]
            tail[f"{response}_tail_mean"] = statistics.fmean(observed)
            tail[f"{response}_tail_sd"] = statistics.stdev(observed)
        early = [
            float(row["TV"])
            for row in values
            if EARLY_WINDOW[0] <= row["generation"] <= EARLY_WINDOW[1]
        ]
        late = [
            float(row["TV"])
            for row in values
            if LATE_WINDOW[0] <= row["generation"] <= LATE_WINDOW[1]
        ]
        tail["TV_51_75_mean"] = statistics.fmean(early)
        tail["TV_76_100_mean"] = statistics.fmean(late)
        tail["TV_late_change"] = tail["TV_76_100_mean"] - tail["TV_51_75_mean"]
        tails.append(tail)
    endpoints.sort(key=lambda row: (row["cell_id"], row["seed_block_id"]))
    tails.sort(key=lambda row: (row["cell_id"], row["seed_block_id"]))
    return endpoints, tails


def _analysis_records(
    endpoints: list[dict[str, Any]], tails: list[dict[str, Any]]
) -> dict[tuple[str, str], dict[str, float]]:
    tail_index = {(row["cell_id"], row["seed_block_id"]): row for row in tails}
    records: dict[tuple[str, str], dict[str, float]] = {}
    for endpoint in endpoints:
        key = (endpoint["cell_id"], endpoint["seed_block_id"])
        tail = tail_index[key]
        records[key] = {
            "D0_g100": float(endpoint["D0"]),
            "D1_g100": float(endpoint["D1"]),
            "D2_g100": float(endpoint["D2"]),
            "evenness_g100": float(endpoint["evenness"]),
            "TV_g100": float(endpoint["TV"]),
            "TV_tail_mean": float(tail["TV_tail_mean"]),
            "TV_tail_sd": float(tail["TV_tail_sd"]),
        }
    return records


def _cell_summaries(
    records: dict[tuple[str, str], dict[str, float]],
    matrix: list[dict[str, str]],
) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for cell in matrix:
        for response in ANALYSIS_RESPONSES:
            values = [
                records[cell["cell_id"], seed][response]
                for seed, _master in SEED_BLOCKS
            ]
            mean, sd, low, high = _mean_ci90(values)
            rows.append(
                {
                    "cell_id": cell["cell_id"],
                    "cell": cell["cell"],
                    "panel": cell["panel"],
                    "H": int(cell["H"]),
                    "B": int(cell["B"]),
                    "HB": int(cell["H"]) * int(cell["B"]),
                    "alpha_target": float(cell["alpha_target"]),
                    "alpha": float(cell["alpha"]),
                    "m": float(cell["m"]),
                    "response": response,
                    "n": len(values),
                    "mean": mean,
                    "sd": sd,
                    "ci90_low": low,
                    "ci90_high": high,
                }
            )
    return rows


def _paired_values(
    records: dict[tuple[str, str], dict[str, float]],
    treatment: str,
    reference: str,
    response: str,
) -> tuple[list[float], str]:
    relative = response in RELATIVE_RESPONSES
    values = []
    for seed, _master in SEED_BLOCKS:
        observed = records[treatment, seed][response]
        baseline = records[reference, seed][response]
        if relative:
            if baseline == 0:
                raise ValueError(f"zero reference for relative contrast: {reference}")
            values.append(observed / baseline - 1.0)
        else:
            values.append(observed - baseline)
    return values, "relative change" if relative else "absolute difference"


def _contrast_row(
    *,
    label: str,
    treatment: str,
    reference: str,
    response: str,
    values: list[float],
    scale: str,
) -> dict[str, Any]:
    mean, sd, low, high = _mean_ci90(values)
    margin = 0.02 if response == "evenness_g100" else 0.05
    status = (
        "increase"
        if low > margin
        else "decrease"
        if high < -margin
        else "equivalent"
        if low >= -margin and high <= margin
        else "uncertain"
    )
    return {
        "contrast": label,
        "treatment": treatment,
        "reference": reference,
        "response": response,
        "n_pairs": len(values),
        "scale": scale,
        "mean": mean,
        "sd": sd,
        "ci90_low": low,
        "ci90_high": high,
        "equivalence_margin": margin,
        "status": status,
    }


def _h_by_b_contrasts(
    records: dict[tuple[str, str], dict[str, float]],
    matrix: list[dict[str, str]],
) -> list[dict[str, Any]]:
    grouped: dict[int, list[dict[str, str]]] = defaultdict(list)
    for cell in matrix:
        if cell["panel"] == "H-by-B":
            grouped[int(cell["H"]) * int(cell["B"])].append(cell)
    rows: list[dict[str, Any]] = []
    for total_founders, cells in sorted(grouped.items()):
        if len(cells) < 2:
            continue
        cells.sort(key=lambda cell: int(cell["H"]))
        for reference, treatment in zip(cells[:-1], cells[1:], strict=True):
            for response in ANALYSIS_RESPONSES:
                values, scale = _paired_values(
                    records,
                    treatment["cell_id"],
                    reference["cell_id"],
                    response,
                )
                row = _contrast_row(
                    label=(
                        f"fixed-HB-{total_founders}-H{treatment['H']}-vs-"
                        f"H{reference['H']}"
                    ),
                    treatment=treatment["cell_id"],
                    reference=reference["cell_id"],
                    response=response,
                    values=values,
                    scale=scale,
                )
                row.update(
                    {
                        "HB": total_founders,
                        "H_treatment": int(treatment["H"]),
                        "B_treatment": int(treatment["B"]),
                        "H_reference": int(reference["H"]),
                        "B_reference": int(reference["B"]),
                    }
                )
                rows.append(row)
    return rows


def _alpha_by_m_contrasts(
    records: dict[tuple[str, str], dict[str, float]],
    matrix: list[dict[str, str]],
) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    cells = {
        (cell["alpha_target"], cell["m"]): cell
        for cell in matrix
        if cell["panel"] == "alpha-by-m"
    }
    effects: list[dict[str, Any]] = []
    effect_values: dict[tuple[str, str, str], list[float]] = {}
    for alpha in ALPHA_LEVELS[1:]:
        for migration in MIGRATION_LEVELS:
            treatment = cells[alpha, migration]
            control = cells["0", migration]
            for response in ANALYSIS_RESPONSES:
                values, scale = _paired_values(
                    records,
                    treatment["cell_id"],
                    control["cell_id"],
                    response,
                )
                effect_values[alpha, migration, response] = values
                row = _contrast_row(
                    label=f"host-effect-alpha-{alpha}-m-{migration}",
                    treatment=treatment["cell_id"],
                    reference=control["cell_id"],
                    response=response,
                    values=values,
                    scale=scale,
                )
                row.update({"alpha_target": float(alpha), "m": float(migration)})
                effects.append(row)

    interactions: list[dict[str, Any]] = []
    for alpha in ALPHA_LEVELS[1:]:
        for migration in MIGRATION_LEVELS[1:]:
            for response in ANALYSIS_RESPONSES:
                at_m = effect_values[alpha, migration, response]
                at_zero = effect_values[alpha, "0", response]
                values = [
                    observed - baseline
                    for observed, baseline in zip(at_m, at_zero, strict=True)
                ]
                row = _contrast_row(
                    label=f"migration-modification-alpha-{alpha}-m-{migration}-vs-0",
                    treatment=cells[alpha, migration]["cell_id"],
                    reference=cells[alpha, "0"]["cell_id"],
                    response=response,
                    values=values,
                    scale="difference in paired host effects",
                )
                row.update(
                    {
                        "alpha_target": float(alpha),
                        "m_treatment": float(migration),
                        "m_reference": 0.0,
                        "interpretation": (
                            "Negative values mean migration reduces the host effect."
                        ),
                    }
                )
                interactions.append(row)
    return effects, interactions


def _fit_ols(x: np.ndarray, y: np.ndarray) -> tuple[np.ndarray, float, float, float]:
    coefficients, _residuals, rank, _singular = np.linalg.lstsq(x, y, rcond=None)
    if rank != x.shape[1]:
        raise ValueError("H-by-B model matrix is rank deficient")
    residual = y - x @ coefficients
    rss = float(residual @ residual)
    tss = float((y - y.mean()) @ (y - y.mean()))
    r_squared = 1.0 - rss / tss if tss else math.nan
    aic = len(y) * math.log(max(rss / len(y), 1e-300)) + 2 * x.shape[1]
    return coefficients, rss, r_squared, aic


def _h_by_b_models(
    records: dict[tuple[str, str], dict[str, float]],
    matrix: list[dict[str, str]],
) -> list[dict[str, Any]]:
    cells = [cell for cell in matrix if cell["panel"] == "H-by-B"]
    log_h = np.asarray([math.log(float(cell["H"])) for cell in cells])
    log_b = np.asarray([math.log(float(cell["B"])) for cell in cells])
    centred_h = log_h - log_h.mean()
    centred_b = log_b - log_b.mean()
    centred_hb = log_h + log_b - (log_h + log_b).mean()
    specifications = {
        "HB-only": (
            np.column_stack([np.ones(len(cells)), centred_hb]),
            ("intercept", "log_HB_centered"),
        ),
        "separate-H-B-interaction": (
            np.column_stack(
                [np.ones(len(cells)), centred_h, centred_b, centred_h * centred_b]
            ),
            ("intercept", "log_H_centered", "log_B_centered", "log_H_x_log_B"),
        ),
    }
    by_model: dict[str, dict[str, list[float]]] = {
        model: defaultdict(list) for model in specifications
    }
    for seed, _master in SEED_BLOCKS:
        observed = np.asarray(
            [records[cell["cell_id"], seed]["TV_g100"] for cell in cells]
        )
        if np.any(observed <= 0):
            raise ValueError("H-by-B log(TV) model requires positive passage-100 TV")
        response = np.log(observed)
        for model, (x, terms) in specifications.items():
            coefficients, rss, r_squared, aic = _fit_ols(x, response)
            for term, value in zip(terms, coefficients, strict=True):
                by_model[model][f"coefficient:{term}"].append(float(value))
            by_model[model]["fit:AIC"].append(aic)
            by_model[model]["fit:R_squared"].append(r_squared)
            by_model[model]["fit:RSS"].append(rss)

    rows: list[dict[str, Any]] = []
    for model, quantities in by_model.items():
        for quantity, values in quantities.items():
            mean, sd, low, high = _mean_ci90(values)
            rows.append(
                {
                    "model": model,
                    "quantity": quantity,
                    "n_seed_blocks": len(values),
                    "mean": mean,
                    "sd": sd,
                    "ci90_low": low,
                    "ci90_high": high,
                    "interpretation": (
                        "Per-seed descriptive fit to log passage-100 TV; "
                        "interval summarizes variation across matched seed blocks."
                    ),
                }
            )
    delta = [
        separate - hb
        for separate, hb in zip(
            by_model["separate-H-B-interaction"]["fit:AIC"],
            by_model["HB-only"]["fit:AIC"],
            strict=True,
        )
    ]
    mean, sd, low, high = _mean_ci90(delta)
    rows.append(
        {
            "model": "separate-minus-HB-only",
            "quantity": "comparison:delta_AIC",
            "n_seed_blocks": len(delta),
            "mean": mean,
            "sd": sd,
            "ci90_low": low,
            "ci90_high": high,
            "interpretation": (
                "Negative values favor separate H, B and interaction terms; "
                "positive values favor the simpler HB-only model."
            ),
        }
    )
    return rows


def _derive_tables(
    trajectories: list[dict[str, Any]], matrix: list[dict[str, str]]
) -> dict[str, list[dict[str, Any]]]:
    endpoints, tails = _run_summaries(trajectories)
    records = _analysis_records(endpoints, tails)
    effects, interactions = _alpha_by_m_contrasts(records, matrix)
    return {
        "environment-trajectories-g100": trajectories,
        "run-endpoints-g100": endpoints,
        "run-tail-summaries-g100": tails,
        "cell-summaries-g100": _cell_summaries(records, matrix),
        "h-by-b-paired-contrasts": _h_by_b_contrasts(records, matrix),
        "h-by-b-model-comparison": _h_by_b_models(records, matrix),
        "alpha-by-m-contrasts": effects,
        "alpha-by-m-interactions": interactions,
    }


def analyse(repository: Path) -> Path:
    repository = repository.resolve()
    work, scratch, runs = load_rows(repository)
    phase = work / "p01-neutral-feedback"
    derived = phase / "analysis" / f"{STAGE_DIRECTORY}-derived"
    derived.mkdir(parents=True, exist_ok=True)
    audit_path = derived / "analysis-audit-g100.json"

    issues = verify_files(repository)
    try:
        _verify_manifest(runs, work=work)
    except (OSError, RuntimeError, ValueError) as exc:
        issues.append(str(exc))
    if len(runs) != EXPECTED_RUNS:
        issues.append(f"manifest has {len(runs)} rather than {EXPECTED_RUNS} new runs")
    issues.extend(
        state_issues(runs, work=work, scratch=scratch, horizon=INITIAL_HORIZON)
    )
    if issues:
        _atomic_json(
            audit_path,
            {
                "status": "FAIL",
                "scope": "Wave 2 passage-100 primary analysis",
                "issues": issues,
            },
        )
        raise RuntimeError("Wave 2 summary gate failed:\n" + "\n".join(issues))

    matrix_path = phase / "design" / f"{EXPERIMENT_ID}-cells.tsv"
    matrix = _read_tsv(matrix_path)
    cells = {row["cell_id"]: row for row in matrix}
    if len(matrix) != len(CELLS) or set(cells) != {cell.cell_id for cell in CELLS}:
        raise RuntimeError("Wave 2 design table does not contain the frozen 40 cells")
    initial_payload = json.loads(
        (work / "common/initial-populations/ip001-fisher100.json").read_text(
            encoding="utf-8"
        )
    )
    initial = {
        index: int(count)
        for index, count in enumerate(initial_payload["scaled_counts"])
    }

    trajectories: list[dict[str, Any]] = []
    inputs: list[dict[str, Any]] = []
    processing_issues: list[str] = []
    for number, run in enumerate(runs, start=1):
        if number == 1 or number % 25 == 0 or number == len(runs):
            print(f"Summarising population {number}/{len(runs)}", flush=True)
        output = scratch / run["scratch_relative_path"]
        try:
            rows, input_record = _new_run_rows(
                run, cells[run["cell_id"]], output, initial
            )
            trajectories.extend(rows)
            inputs.append(input_record)
        except (KeyError, OSError, TypeError, ValueError) as exc:
            processing_issues.append(f"{run['run_id']}: {exc}")
    try:
        reused, reused_inputs = _reused_rows(phase, cells)
        trajectories.extend(reused)
        inputs.extend(reused_inputs)
    except (KeyError, OSError, TypeError, ValueError) as exc:
        processing_issues.append(f"reused passage-100 populations: {exc}")
    if processing_issues:
        _atomic_json(
            audit_path,
            {
                "status": "FAIL",
                "scope": "Wave 2 passage-100 primary analysis",
                "issues": processing_issues,
            },
        )
        raise RuntimeError(
            "Wave 2 biological summary failed:\n" + "\n".join(processing_issues)
        )

    trajectories.sort(
        key=lambda row: (row["cell_id"], row["seed_block_id"], row["generation"])
    )
    expected_populations = len(CELLS) * len(SEED_BLOCKS)
    expected_trajectories = expected_populations * (INITIAL_HORIZON + 1)
    if (
        len(trajectories) != expected_trajectories
        or len(inputs) != expected_populations
    ):
        raise RuntimeError(
            "compiled population count differs from the frozen Wave 2 design"
        )
    tables = _derive_tables(trajectories, matrix)
    tables["analysis-inputs-g100"] = sorted(
        inputs, key=lambda row: (row["cell_id"], row["seed_block_id"])
    )
    for name, rows in tables.items():
        _atomic_tsv(derived / f"{name}.tsv", rows, list(rows[0]))

    table_records = [
        {
            "path": str((derived / f"{name}.tsv").relative_to(work)),
            "rows": len(rows),
            "sha256": _sha256(derived / f"{name}.tsv"),
        }
        for name, rows in tables.items()
    ]
    summary = {
        "experiment_id": EXPERIMENT_ID,
        "status": "PASS",
        "scope": "complete passage-100 primary Wave 2 analysis",
        "new_populations": len(runs),
        "reused_populations": len(REUSED_CELLS) * len(SEED_BLOCKS),
        "total_populations": expected_populations,
        "cells": len(CELLS),
        "seed_blocks_per_cell": len(SEED_BLOCKS),
        "environmental_states": expected_trajectories,
        "primary_endpoint": INITIAL_HORIZON,
        "tail_window": [TAIL_START, INITIAL_HORIZON],
        "uncertainty": "paired 90% Student-t intervals across 12 seed blocks",
        "interpretation": (
            "Passage 100 is the complete exploratory primary endpoint. Adaptive "
            "continuations do not alter these tables."
        ),
        "tables": table_records,
    }
    _atomic_json(derived / "analysis-summary-g100.json", summary)
    _atomic_json(
        audit_path,
        {
            "status": "PASS",
            "scope": summary["scope"],
            "issues": [],
            "new_populations": len(runs),
            "reused_populations": len(REUSED_CELLS) * len(SEED_BLOCKS),
            "total_populations": expected_populations,
            "input_records": len(inputs),
            "environmental_states": expected_trajectories,
            "design_sha256": _sha256(matrix_path),
            "tables": table_records,
        },
    )
    print(
        f"Compiled {expected_populations} populations and "
        f"{expected_trajectories} environmental states in {derived}."
    )
    return derived


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--repository", type=Path, default=Path(__file__).resolve().parents[1]
    )
    args = parser.parse_args()
    analyse(args.repository)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
