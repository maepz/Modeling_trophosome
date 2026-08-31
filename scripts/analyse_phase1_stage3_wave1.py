#!/usr/bin/env python3
"""Audit and summarise the Stage 3 first wave as finite-time exploratory results."""

from __future__ import annotations

import argparse
import csv
import json
import math
import statistics
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np
from prepare_phase1_stage3_wave1 import (
    CELLS,
    CONTROL_ID,
    EXPECTED_RUNS,
    EXPERIMENT_ID,
    FEEDBACK_LEVELS,
    HOST_GENERATIONS,
    HOST_LEVELS,
    MUTATION_LEVELS,
    SEED_BLOCKS,
    STAGE_DIRECTORY,
    TAIL_START,
    K,
    verify_files,
)
from run_phase1_first_pilot import _atomic_json, _sha256
from run_phase1_stage3_wave1 import completion_issues, load_rows

from trophosome.second_pilot_report import _configure_matplotlib

_configure_matplotlib()
# Reuse the audited biological metrics and memory-bounded lineage tracing.
ANALYSIS = (
    Path(__file__).resolve().parents[1]
    / "experiments/work/trophosome/p01-neutral-feedback/analysis"
)
sys.path.insert(0, str(ANALYSIS))
from analyse_first_pilot import (  # noqa: E402
    composition_metrics,
    diversity_metrics,
    read_lineages,
    read_tsv,
    trace_roots,
    write_tsv,
)
from analyse_second_pilot import _root_collapsed_counts  # noqa: E402

RESPONSES = ("D0", "D1", "D2", "evenness", "TV")
DISPLAY_RESPONSES = (*RESPONSES, "shannon", "simpson")
CELL_LOOKUP = {
    (c.hosts, c.alpha_target, c.mutation_probability): c.cell_id for c in CELLS
}


def contrast_definitions() -> list[tuple[str, str, str]]:
    result = []
    for alpha in FEEDBACK_LEVELS:
        for u in MUTATION_LEVELS:
            for small, large in zip(HOST_LEVELS[:-1], HOST_LEVELS[1:], strict=True):
                result.append(
                    (
                        f"pooling-H{large}-vs-{small}-a{alpha}-u{u}",
                        CELL_LOOKUP[large, alpha, u],
                        CELL_LOOKUP[small, alpha, u],
                    )
                )
    for h in HOST_LEVELS:
        for u in MUTATION_LEVELS:
            for low, high in zip(
                FEEDBACK_LEVELS[:-1], FEEDBACK_LEVELS[1:], strict=True
            ):
                result.append(
                    (
                        f"feedback-H{h}-a{high}-vs-{low}-u{u}",
                        CELL_LOOKUP[h, high, u],
                        CELL_LOOKUP[h, low, u],
                    )
                )
        for alpha in FEEDBACK_LEVELS:
            result.append(
                (
                    f"mutation-H{h}-a{alpha}",
                    CELL_LOOKUP[h, alpha, "1e-10"],
                    CELL_LOOKUP[h, alpha, "0"],
                )
            )
    return result


def mean_ci90(values: list[float]) -> tuple[float, float, float]:
    """Two-sided Student t interval; twelve independent seed blocks, df=11."""
    if len(values) != len(SEED_BLOCKS) or not all(math.isfinite(x) for x in values):
        raise ValueError("interval requires twelve finite population values")
    mean = statistics.fmean(values)
    half_width = 1.7958848187036691 * statistics.stdev(values) / math.sqrt(len(values))
    return mean, mean - half_width, mean + half_width


def read_trajectory(path: Path) -> dict[int, dict[int, int]]:
    result: dict[int, dict[int, int]] = defaultdict(dict)
    with path.open(newline="") as handle:
        for row in csv.DictReader(handle):
            generation, strain, count = (
                int(row["generation"]),
                int(row["strain_id"]),
                int(row["count"]),
            )
            if int(row["replicate"]) != 0 or strain < 0 or count <= 0:
                raise ValueError("invalid environmental row")
            if strain in result[generation]:
                raise ValueError("duplicate environmental generation/strain row")
            result[generation][strain] = count
    if set(result) != set(range(HOST_GENERATIONS + 1)):
        raise ValueError("environment trajectory must contain generations 0-100")
    if any(sum(counts.values()) != K for counts in result.values()):
        raise ValueError("environment capacity differs from one billion cells")
    return dict(result)


def summarise_cells(endpoints: list[dict], matrix: list[dict]) -> list[dict]:
    rows = []
    for cell in matrix:
        selected = [r for r in endpoints if r["cell_id"] == cell["cell_id"]]
        for response in DISPLAY_RESPONSES:
            mean, low, high = mean_ci90([float(r[response]) for r in selected])
            margin = (
                0.05 * abs(mean)
                if response in ("D0", "D1", "D2")
                else (
                    0.02
                    if response == "evenness"
                    else 0.05
                    if response == "TV"
                    else None
                )
            )
            rows.append(
                {
                    "cell_id": cell["cell_id"],
                    "cell": cell["cell"],
                    "response": response,
                    "n": len(selected),
                    "mean": mean,
                    "ci90_low": low,
                    "ci90_high": high,
                    "sd": statistics.stdev(float(r[response]) for r in selected),
                    "precision_margin": margin if margin is not None else "",
                    "more_seeds_flag": margin is not None and (high - low) / 2 > margin,
                }
            )
    return rows


def paired_contrasts(
    endpoints: list[dict], references: list[dict], matrix: list[dict]
) -> list[dict]:
    indexed = {(r["cell_id"], r["seed_block_id"]): r for r in endpoints + references}
    if len(indexed) != len(endpoints) + len(references):
        raise ValueError("duplicate cell/seed endpoint")
    contrasts = contrast_definitions() + [
        (f"host-return-{c['cell']}", c["cell_id"], CONTROL_ID)
        for c in matrix
        if c["cell_id"] != CONTROL_ID
    ]
    rows = []
    for label, treatment, reference in contrasts:
        for response in RESPONSES:
            relative = response in ("D0", "D1", "D2")
            values = []
            for seed, _ in SEED_BLOCKS:
                a, b = (
                    float(indexed[treatment, seed][response]),
                    float(indexed[reference, seed][response]),
                )
                values.append(a / b - 1 if relative else a - b)
            mean, low, high = mean_ci90(values)
            margin = 0.05 if relative or response == "TV" else 0.02
            status = (
                "increase"
                if low > margin
                else "decrease"
                if high < -margin
                else (
                    "equivalent" if low >= -margin and high <= margin else "uncertain"
                )
            )
            rows.append(
                {
                    "contrast": label,
                    "treatment": treatment,
                    "reference": reference,
                    "response": response,
                    "n_pairs": len(SEED_BLOCKS),
                    "scale": "relative change" if relative else "absolute difference",
                    "mean": mean,
                    "ci90_low": low,
                    "ci90_high": high,
                    "margin": margin,
                    "status": status,
                }
            )
    return rows


def tv_interactions(endpoints: list[dict]) -> list[dict]:
    """Difference of feedback effects across H, calculated within each seed."""
    indexed = {(r["cell_id"], r["seed_block_id"]): float(r["TV"]) for r in endpoints}
    result = []
    for small, large in zip(HOST_LEVELS[:-1], HOST_LEVELS[1:], strict=True):
        for low, high in zip(FEEDBACK_LEVELS[:-1], FEEDBACK_LEVELS[1:], strict=True):
            by_u = {}
            for u in MUTATION_LEVELS:

                def tv(h, a, seed, mutation=u):
                    return indexed[CELL_LOOKUP[h, a, mutation], seed]

                by_u[u] = [
                    (tv(large, high, seed) - tv(large, low, seed))
                    - (tv(small, high, seed) - tv(small, low, seed))
                    for seed, _ in SEED_BLOCKS
                ]
            for u, values in (
                *by_u.items(),
                (
                    "mutation-modification",
                    [a - b for a, b in zip(by_u["1e-10"], by_u["0"], strict=True)],
                ),
            ):
                mean, low_ci, high_ci = mean_ci90(values)
                result.append(
                    {
                        "H_low": small,
                        "H_high": large,
                        "alpha_low_target": low,
                        "alpha_high_target": high,
                        "mutation": u,
                        "response": "TV",
                        "n_pairs": len(SEED_BLOCKS),
                        "mean": mean,
                        "ci90_low": low_ci,
                        "ci90_high": high_ci,
                        "interpretation": "Negative: smaller feedback effect at high H"
                        if u != "mutation-modification"
                        else "Change in H x alpha interaction caused by mutation",
                    }
                )
    return result


def audit_migration(path: Path) -> None:
    totals = defaultdict(lambda: [0, 0])
    with path.open(newline="") as handle:
        for row in csv.DictReader(handle):
            if int(row["replicate"]) != 0:
                raise ValueError("unexpected migration replicate")
            e, i = int(row["emigrant_count"]), int(row["immigrant_count"])
            if min(e, i) < 0:
                raise ValueError("negative migration count")
            totals[int(row["generation"])][0] += e
            totals[int(row["generation"])][1] += i
    if set(totals) != set(range(1, HOST_GENERATIONS + 1)) or any(
        pair != [K // 10, K // 10] for pair in totals.values()
    ):
        raise ValueError("migration generations or replacement totals differ")


def classify(contrasts: list[dict], matrix: list[dict]) -> list[dict]:
    result = []
    for cell in matrix:
        if cell["cell_id"] == CONTROL_ID:
            continue
        responses = {
            r["response"]: r["status"]
            for r in contrasts
            if r["contrast"] == f"host-return-{cell['cell']}"
        }
        if responses["D0"] == "increase" and responses["evenness"] == "decrease":
            category = "richness_increase_evenness_decrease"
        elif responses["D1"] == responses["D2"] == "decrease":
            category = "effective_diversity_decrease"
        elif responses["D1"] == responses["D2"] == "increase":
            category = "effective_diversity_increase"
        elif all(s == "equivalent" for s in responses.values()):
            category = "negligible_for_five_measured_statistics"
        else:
            category = "mixed_or_unresolved"
        result.append(
            {
                "cell_id": cell["cell_id"],
                "cell": cell["cell"],
                "classification": category,
                **responses,
            }
        )
    return result


def analyse(repository: Path) -> Path:
    work, scratch, runs = load_rows(repository)
    phase = work / "p01-neutral-feedback"
    derived = phase / "analysis" / f"{STAGE_DIRECTORY}-derived"
    derived.mkdir(parents=True, exist_ok=True)
    audit_path = derived / "analysis-audit.json"
    issues = verify_files(repository)
    if len(runs) != EXPECTED_RUNS:
        issues.append("manifest does not contain all 288 new populations")
    issues.extend(completion_issues(runs, work=work, scratch=scratch))
    if issues:
        _atomic_json(audit_path, {"status": "FAIL", "issues": issues})
        raise RuntimeError("Report completion gate failed:\n" + "\n".join(issues))
    matrix = read_tsv(phase / "design" / f"{EXPERIMENT_ID}-cells.tsv")
    by_cell = {r["cell_id"]: r for r in matrix}
    references = read_tsv(phase / "design" / f"{EXPERIMENT_ID}-reference-endpoints.tsv")
    initial = json.loads(
        (work / "common/initial-populations/ip001-fisher100.json").read_text()
    )["scaled_counts"]
    initial_counts = dict(enumerate(initial))
    endpoints, trajectories, tails, resource_rows = [], [], [], []
    source_hashes = set()
    for run in runs:
        print(f"Analysing {run['run_id']}", flush=True)
        output = scratch / run["scratch_relative_path"]
        try:
            states = read_trajectory(output / "environment_counts.csv")
            if states[0] != initial_counts:
                raise ValueError("initial environment differs from frozen population")
            with np.load(
                output / "final_environment_rep000.npz", allow_pickle=False
            ) as final:
                observed = {
                    int(k): int(v)
                    for k, v in zip(final["genotype_ids"], final["counts"], strict=True)
                }
            if states[HOST_GENERATIONS] != observed:
                raise ValueError("final archive differs from passage-100 counts")
            with (output / "host_generation_summary.csv").open() as handle:
                summaries = list(csv.DictReader(handle))
            if len(summaries) != HOST_GENERATIONS or {
                int(r["host_generation"]) for r in summaries
            } != set(range(1, HOST_GENERATIONS + 1)):
                raise ValueError("host-generation summary is incomplete or duplicated")
            audit_migration(output / "migration_counts.csv")
            ids = set().union(*(counts.keys() for counts in states.values()))
            parents, lineages, max_mutants = read_lineages(
                output / "strain_lineage_events.csv", ids
            )
            if any(child <= parent for child, parent in parents.items()):
                raise ValueError(
                    "lineage identifiers do not follow parent-before-child order"
                )
            if max_mutants > 100000:
                raise ValueError(
                    "a transition exceeds the frozen mutant materialization cap"
                )
            roots = trace_roots(parents, ids)
            if by_cell[run["cell_id"]]["u"] == "0" and (
                lineages or any(i >= 100 for i in ids)
            ):
                raise ValueError("new mutation labels exist in a mutation-free cell")
            local = []
            previous = states[0]
            for generation, counts in sorted(states.items()):
                metrics = diversity_metrics(counts)
                collapsed = diversity_metrics(_root_collapsed_counts(counts, roots))
                row = {
                    "run_id": run["run_id"],
                    "cell_id": run["cell_id"],
                    "cell": run["cell"],
                    "seed_block_id": run["seed_block_id"],
                    "generation": generation,
                    "D0": metrics["richness"],
                    "shannon": metrics["shannon"],
                    "simpson": metrics["simpson"],
                    "D1": metrics["hill_q1"],
                    "D2": metrics["hill_q2"],
                    "evenness": metrics["evenness"],
                    "TV": composition_metrics(counts, initial_counts)[
                        "total_variation"
                    ],
                    "turnover": composition_metrics(counts, previous)[
                        "total_variation"
                    ],
                    "root_D0": collapsed["richness"],
                    "root_D1": collapsed["hill_q1"],
                    "root_D2": collapsed["hill_q2"],
                    "mutant_richness": sum(i >= 100 for i in counts),
                    "mutant_abundance_fraction": sum(
                        n for i, n in counts.items() if i >= 100
                    )
                    / K,
                }
                local.append(row)
                previous = counts
            trajectories.extend(local)
            endpoints.append(local[-1])
            tail = {
                k: local[-1][k] for k in ("run_id", "cell_id", "cell", "seed_block_id")
            }
            for response in RESPONSES:
                values = [
                    float(r[response])
                    for r in local
                    if TAIL_START <= r["generation"] <= HOST_GENERATIONS
                ]
                tail[f"{response}_mean"] = statistics.fmean(values)
                tail[f"{response}_sd"] = statistics.stdev(values)
            tail["TV_late_change"] = statistics.fmean(
                r["TV"] for r in local if 76 <= r["generation"] <= 100
            ) - statistics.fmean(r["TV"] for r in local if 51 <= r["generation"] <= 75)
            tails.append(tail)
            completed = json.loads((output / "completion.json").read_text())
            source_hashes.add(completed["source_sha256"])
            execution = json.loads((output / "execution-summary.json").read_text())
            resource_rows.append(
                {
                    "run_id": run["run_id"],
                    "generated_lineages": lineages,
                    "maximum_transition_mutants": max_mutants,
                    "last_attempt_hours": float(execution["elapsed_seconds"]) / 3600,
                    "resumed": execution.get("resumed", False),
                    "output_gib": sum(
                        p.stat().st_size for p in output.rglob("*") if p.is_file()
                    )
                    / 1024**3,
                    "source_sha256": completed["source_sha256"],
                    "completion_sha256": _sha256(output / "completion.json"),
                }
            )
        except (ValueError, OSError, KeyError, RuntimeError) as exc:
            issues.append(f"{run['run_id']}: {exc}")
    if len(source_hashes) != 1:
        issues.append(
            "Stage 3 populations were not all produced by one frozen model source"
        )
    if issues:
        _atomic_json(audit_path, {"status": "FAIL", "issues": issues})
        raise RuntimeError("Biological/output audit failed:\n" + "\n".join(issues))
    controls = [r for r in references if r["cell_id"] == CONTROL_ID]
    primary = endpoints + controls
    summaries = summarise_cells(primary, matrix)
    supplementary = [r for r in references if r["cell_id"] != CONTROL_ID]
    supplementary_cells = list({r["cell_id"]: r for r in supplementary}.values())
    reference_summaries = summarise_cells(supplementary, supplementary_cells)
    contrasts = paired_contrasts(endpoints, references, matrix)
    interactions = tv_interactions(endpoints)
    categories = classify(contrasts, matrix)
    drift = []
    for cell in CELLS:
        values = [r["TV_late_change"] for r in tails if r["cell_id"] == cell.cell_id]
        mean, low, high = mean_ci90(values)
        drift.append(
            {
                "cell_id": cell.cell_id,
                "cell": cell.short_id,
                "n": len(values),
                "mean": mean,
                "ci90_low": low,
                "ci90_high": high,
                "margin": 0.05,
                "review_longer_runs": not (-0.05 <= low and high <= 0.05),
                "scope": "Exploratory mean-drift check; not proof of stationarity",
            }
        )
    for name, rows in (
        ("environment-trajectories", trajectories),
        ("run-endpoints", endpoints),
        ("reused-reference-endpoints", references),
        ("supplementary-reference-summaries", reference_summaries),
        ("run-tail-summaries", tails),
        ("cell-summaries", summaries),
        ("paired-contrasts", contrasts),
        ("tv-interactions", interactions),
        ("late-tv-drift", drift),
        ("classifications", categories),
        ("resource-summary", resource_rows),
    ):
        write_tsv(derived / f"{name}.tsv", rows, list(rows[0]))
    _atomic_json(
        audit_path,
        {
            "status": "PASS",
            "issues": [],
            "populations": len(runs),
            "primary_populations_including_reused_control": len(primary),
            "source_sha256": next(iter(source_hashes)),
        },
    )
    _atomic_json(
        derived / "analysis-summary.json",
        {
            "experiment_id": EXPERIMENT_ID,
            "audit_status": "PASS",
            "data_origin": "synthetic-test-fixture"
            if source_hashes == {"synthetic-test-source"}
            else "simulation",
            "populations": len(runs),
            "cells": len(matrix),
            "seeds_per_cell": len(SEED_BLOCKS),
            "primary_endpoint": HOST_GENERATIONS,
            "primary_response": "TV",
            "primary_populations_including_reused_control": len(primary),
            "new_cells": len(CELLS),
            "reused_control_populations": len(controls),
            "supplementary_reference_populations": len(supplementary),
            "interpretation": (
                "Exploratory finite-time outcomes, not equilibrium claims."
            ),
            "reference_populations": len(references),
            "total_output_gib": sum(r["output_gib"] for r in resource_rows),
            "more_seeds_flags": sum(bool(r["more_seeds_flag"]) for r in summaries),
            "uncertainty": (
                "90% Student-t intervals across twelve independent seed blocks; "
                "not multiplicity-adjusted"
            ),
            "TV_contrast": (
                "Difference in distance from the initial community, not distance "
                "between mutant IDs from different runs."
            ),
            "next_step": (
                "Review paired H x feedback TV interactions, precision, and late drift "
                "before choosing additional seeds, intermediate cells, "
                "or longer validation runs. "
                "No automatic launches."
            ),
        },
    )
    return derived


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--repository", type=Path, default=Path(__file__).resolve().parents[1]
    )
    args = parser.parse_args()
    print(analyse(args.repository.resolve()))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
