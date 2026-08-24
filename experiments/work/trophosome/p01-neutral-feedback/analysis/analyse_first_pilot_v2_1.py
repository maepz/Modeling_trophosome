#!/usr/bin/env python3
"""Audit and summarize the completed model-2.1 fixed-pool first pilot.

The script is read-only with respect to raw scratch outputs. It requires all 60
populations to be complete and internally consistent before publishing the
standardized analysis tables consumed by ``trophosome report``.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
from collections import defaultdict
from datetime import UTC, datetime
from pathlib import Path
from typing import Any

from analyse_first_pilot import (
    INITIAL_IDS,
    INITIAL_RICHNESS,
    SEED_ORDER,
    composition_metrics,
    describe,
    diversity_metrics,
    pooling_metrics,
    read_environment,
    read_lineages,
    read_tsv,
    sha256,
    trace_roots,
    unique_mutant_ids,
    write_tsv,
)

VARIANT_TAG = "v210-m010"
EXPERIMENT_ID = f"phase1-first-pilot-{VARIANT_TAG}"
MODEL_SPEC_VERSION = "2.1.0"
OUTPUT_SCHEMA_VERSION = "2.3.0"
SOFTWARE_VERSION = "0.7.0"
MIGRATION_FRACTION = 0.1
MIGRATION_REPLACEMENT_COUNT = 100_000_000
ENVIRONMENT_CAPACITY = 1_000_000_000
HOST_GENERATIONS = 5
CELL_ORDER = tuple(f"p01-s01-c{index:04d}" for index in range(1, 21))
REQUIRED_RAW_FILES = (
    "completion.json",
    "environment_counts.csv",
    "execution-summary.json",
    "final_environment_rep000.npz",
    "host_adult_counts.csv",
    "host_generation_summary.csv",
    "infection_counts.csv",
    "migration_counts.csv",
    "pooled_host_counts_and_occupancy.csv",
    "provenance.json",
    "release_counts.csv",
    "strain_lineage_events.csv",
    "strain_origins.csv",
)


def _atomic_json(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    temporary.replace(path)


def _audit_payload(status: str, issues: list[str]) -> dict[str, Any]:
    return {
        "analysis_audit_schema_version": "1.0.0",
        "experiment_id": EXPERIMENT_ID,
        "variant_id": VARIANT_TAG,
        "checked_at": datetime.now(UTC).isoformat(),
        "status": status,
        "expected_populations": 60,
        "issues": issues,
    }


def completion_gate_issues(
    run_rows: list[dict[str, str]], *, work: Path, scratch: Path
) -> list[str]:
    """Return every reason the 60-run set is not ready for analysis."""

    issues: list[str] = []
    if len(run_rows) != 60:
        issues.append(f"manifest contains {len(run_rows)} runs; expected 60")
    if len({row.get("run_id") for row in run_rows}) != len(run_rows):
        issues.append("manifest run IDs are not unique")
    if len({row.get("scratch_relative_path") for row in run_rows}) != len(run_rows):
        issues.append("manifest scratch paths are not unique")

    for run in run_rows:
        run_id = run.get("run_id", "unknown-run")
        config = work / run.get("config_path", "")
        if not config.is_file():
            issues.append(f"{run_id}: configuration is missing")
        elif sha256(config) != run.get("config_sha256"):
            issues.append(f"{run_id}: configuration checksum differs from manifest")

        output = scratch / run.get("scratch_relative_path", "")
        missing = [name for name in REQUIRED_RAW_FILES if not (output / name).is_file()]
        if missing:
            issues.append(f"{run_id}: missing {', '.join(missing)}")
            continue

        try:
            completion = json.loads(
                (output / "completion.json").read_text(encoding="utf-8")
            )
        except (OSError, ValueError) as exc:
            issues.append(f"{run_id}: unreadable completion.json ({exc})")
            continue
        if completion.get("complete") is not True:
            issues.append(f"{run_id}: completion marker is not committed")
        expected_versions = {
            "model_spec_version": MODEL_SPEC_VERSION,
            "output_schema_version": OUTPUT_SCHEMA_VERSION,
            "software_version": SOFTWARE_VERSION,
        }
        for field, expected in expected_versions.items():
            if completion.get(field) != expected:
                issues.append(
                    f"{run_id}: {field}={completion.get(field)!r}; "
                    f"expected {expected!r}"
                )
        final_path = output / "final_environment_rep000.npz"
        recorded_final = completion.get("final_environment_sha256", {}).get(
            final_path.name
        )
        if not recorded_final or sha256(final_path) != recorded_final:
            issues.append(f"{run_id}: final-environment checksum differs")
        for name, expected_size in completion.get("output_sizes", {}).items():
            path = output / name
            if not path.is_file() or path.stat().st_size != int(expected_size):
                issues.append(f"{run_id}: committed size differs for {name}")
    return issues


def _migration_metrics(
    path: Path, run_id: str, audit_issues: list[str]
) -> dict[str, float | int]:
    by_generation: dict[int, dict[str, dict[int, int]]] = defaultdict(
        lambda: {"emigrant": {}, "immigrant": {}}
    )
    with path.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            generation = int(row["generation"])
            strain_id = int(row["strain_id"])
            by_generation[generation]["emigrant"][strain_id] = int(
                row["emigrant_count"]
            )
            by_generation[generation]["immigrant"][strain_id] = int(
                row["immigrant_count"]
            )

    if set(by_generation) != set(range(1, HOST_GENERATIONS + 1)):
        audit_issues.append(
            f"{run_id}: migration table does not contain generations "
            f"1-{HOST_GENERATIONS}"
        )
    for generation, values in sorted(by_generation.items()):
        emigrants = sum(values["emigrant"].values())
        immigrants = sum(values["immigrant"].values())
        if emigrants != MIGRATION_REPLACEMENT_COUNT:
            audit_issues.append(
                f"{run_id}: generation {generation} has {emigrants} emigrants"
            )
        if immigrants != MIGRATION_REPLACEMENT_COUNT:
            audit_issues.append(
                f"{run_id}: generation {generation} has {immigrants} immigrants"
            )

    final = by_generation.get(HOST_GENERATIONS, {"emigrant": {}, "immigrant": {}})
    keys = set(final["emigrant"]) | set(final["immigrant"])
    total = sum(final["emigrant"].values())
    sample_tv = (
        0.5
        * sum(
            abs(final["emigrant"].get(key, 0) - final["immigrant"].get(key, 0))
            for key in keys
        )
        / total
        if total
        else math.nan
    )
    return {
        "migration_total_emigrants": sum(
            sum(values["emigrant"].values()) for values in by_generation.values()
        ),
        "migration_total_immigrants": sum(
            sum(values["immigrant"].values()) for values in by_generation.values()
        ),
        "migration_emigrant_richness_g5": sum(
            count > 0 for count in final["emigrant"].values()
        ),
        "migration_immigrant_richness_g5": sum(
            count > 0 for count in final["immigrant"].values()
        ),
        "migration_sample_tv_g5": sample_tv,
    }


def _load_cells(path: Path) -> dict[str, dict[str, Any]]:
    cells: dict[str, dict[str, Any]] = {}
    for row in read_tsv(path):
        cells[row["cell_id"]] = {
            **row,
            "H": int(row["H"]),
            "f": float(row["f"]),
            "e": int(row["e"]),
            "R": int(row["R"]),
            "alpha": float(row["alpha"]),
            "u": float(row["u"]),
            "m": float(row["m"]),
        }
    if set(cells) != set(CELL_ORDER):
        raise RuntimeError("v2.1 design does not contain the expected 20 cells")
    return cells


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--repository",
        type=Path,
        default=Path(__file__).resolve().parents[5],
    )
    args = parser.parse_args()
    repository = args.repository.resolve()
    work = repository / "experiments" / "work" / "trophosome"
    phase = work / "p01-neutral-feedback"
    layout = json.loads((work / "layout.local.json").read_text(encoding="utf-8"))
    scratch = Path(layout["scratch"])
    design_path = phase / "design" / f"phase1-first-pilot-{VARIANT_TAG}-cells.tsv"
    manifest_path = phase / "manifests" / f"phase1-first-pilot-{VARIANT_TAG}-runs.tsv"
    derived_dir = phase / "analysis" / f"s01-pilot-{VARIANT_TAG}-derived"
    audit_path = derived_dir / "analysis-audit.json"
    cells = _load_cells(design_path)
    run_rows = read_tsv(manifest_path)

    preflight_issues = completion_gate_issues(run_rows, work=work, scratch=scratch)
    if preflight_issues:
        _atomic_json(audit_path, _audit_payload("FAIL", preflight_issues))
        raise SystemExit(
            "v2.1 pilot is not ready for reporting:\n" + "\n".join(preflight_issues)
        )

    initial_payload = json.loads(
        (work / "common" / "initial-populations" / "ip001-fisher100.json").read_text(
            encoding="utf-8"
        )
    )
    initial_counts = {
        strain_id: int(count)
        for strain_id, count in enumerate(initial_payload["scaled_counts"])
    }
    initial_metrics = diversity_metrics(initial_counts)
    endpoints: list[dict[str, Any]] = []
    trajectories: list[dict[str, Any]] = []
    audit_issues: list[str] = []
    commits: set[str] = set()
    source_hashes: set[str] = set()
    platforms: set[str] = set()

    for run in run_rows:
        run_id = run["run_id"]
        cell_id = run["cell_id"]
        cell = cells[cell_id]
        output = scratch / run["scratch_relative_path"]
        final_counts = read_environment(output / "environment_counts.csv")
        if sum(final_counts.values()) != ENVIRONMENT_CAPACITY:
            audit_issues.append(f"{run_id}: final environmental capacity differs")
        final_metrics = diversity_metrics(final_counts)
        composition = composition_metrics(final_counts, initial_counts)
        parents, generated_lineages, maximum_transition_mutants = read_lineages(
            output / "strain_lineage_events.csv", final_counts
        )
        roots = trace_roots(parents, final_counts)
        root_counts: dict[int, int] = defaultdict(int)
        for strain_id, count in final_counts.items():
            root_counts[roots[strain_id]] += count
        root_metrics = diversity_metrics(dict(root_counts))
        root_composition = composition_metrics(dict(root_counts), initial_counts)
        adult_mutants = unique_mutant_ids(output / "host_adult_counts.csv")
        released_mutants = unique_mutant_ids(output / "release_counts.csv")
        final_mutants = {
            strain_id for strain_id in final_counts if strain_id >= INITIAL_RICHNESS
        }
        mutant_abundance = sum(final_counts[strain_id] for strain_id in final_mutants)
        pooling = pooling_metrics(output, cell["H"])
        migration = _migration_metrics(
            output / "migration_counts.csv", run_id, audit_issues
        )
        execution = json.loads(
            (output / "execution-summary.json").read_text(encoding="utf-8")
        )
        provenance = json.loads(
            (output / "provenance.json").read_text(encoding="utf-8")
        )
        commits.add(str(provenance.get("git_commit")))
        source_hashes.add(str(provenance["source_sha256"]))
        platforms.add(str(provenance["platform"]))
        for field, expected in (
            ("model_spec_version", MODEL_SPEC_VERSION),
            ("output_schema_version", OUTPUT_SCHEMA_VERSION),
            ("software_version", SOFTWARE_VERSION),
        ):
            if provenance.get(field) != expected:
                audit_issues.append(f"{run_id}: provenance {field} differs")

        with (output / "host_generation_summary.csv").open(
            newline="", encoding="utf-8"
        ) as handle:
            summary_rows = list(csv.DictReader(handle))
        if len(summary_rows) != HOST_GENERATIONS:
            audit_issues.append(
                f"{run_id}: expected {HOST_GENERATIONS} host-generation summaries"
            )
            continue
        for summary in summary_rows:
            generation = int(summary["host_generation"])
            replacement_count = int(summary["migration_replacement_count"])
            realized_migration = float(summary["realized_migration_fraction"])
            expected_feedback = float(summary["expected_host_feedback_after_migration"])
            realized_feedback = float(summary["realized_host_feedback"])
            if replacement_count != MIGRATION_REPLACEMENT_COUNT:
                audit_issues.append(
                    f"{run_id}: generation {generation} replacement count differs"
                )
            if not math.isclose(realized_migration, MIGRATION_FRACTION):
                audit_issues.append(
                    f"{run_id}: generation {generation} realized m differs"
                )
            if not math.isclose(
                expected_feedback,
                realized_feedback * (1.0 - MIGRATION_FRACTION),
            ):
                audit_issues.append(
                    f"{run_id}: generation {generation} post-migration feedback differs"
                )
            if int(summary["post_migration_richness"]) != int(
                summary["environment_richness"]
            ):
                audit_issues.append(
                    f"{run_id}: generation {generation} selection-neutral "
                    "richness differs"
                )

        final_summary = summary_rows[-1]
        peak_rss_kib = execution.get("peak_process_tree_rss_kib")
        if peak_rss_kib is None:
            audit_issues.append(f"{run_id}: peak process-tree memory is missing")
            peak_rss_kib = 0
        endpoint: dict[str, Any] = {
            "run_id": run_id,
            "cell_id": cell_id,
            "cell": cell_id.rsplit("-", 1)[-1],
            "seed_block_id": run["seed_block_id"],
            "master_seed": int(run["master_seed"]),
            "H": cell["H"],
            "f": cell["f"],
            "e": cell["e"],
            "R": cell["R"],
            "alpha": cell["alpha"],
            "u": cell["u"],
            "m": cell["m"],
            **final_metrics,
            "richness_change_pct": 100.0
            * (final_metrics["richness"] / initial_metrics["richness"] - 1.0),
            "hill_q1_change_pct": 100.0
            * (final_metrics["hill_q1"] / initial_metrics["hill_q1"] - 1.0),
            "hill_q2_change_pct": 100.0
            * (final_metrics["hill_q2"] / initial_metrics["hill_q2"] - 1.0),
            "evenness_change": final_metrics["evenness"] - initial_metrics["evenness"],
            **composition,
            "original_strains_lost": sum(
                strain_id not in final_counts for strain_id in INITIAL_IDS
            ),
            "original_strains_retained": sum(
                strain_id in final_counts for strain_id in INITIAL_IDS
            ),
            "final_mutant_richness": len(final_mutants),
            "mutant_abundance_fraction": mutant_abundance / sum(final_counts.values()),
            "generated_mutant_lineages": generated_lineages,
            "maximum_mutants_in_one_transition": maximum_transition_mutants,
            "adult_mutant_lineages": len(adult_mutants),
            "released_mutant_lineages": len(released_mutants),
            "root_richness": root_metrics["richness"],
            "root_hill_q1": root_metrics["hill_q1"],
            "root_hill_q2": root_metrics["hill_q2"],
            "root_evenness": root_metrics["evenness"],
            "root_total_variation": root_composition["total_variation"],
            **pooling,
            **migration,
            "post_return_richness_g5": int(final_summary["post_return_richness"]),
            "post_migration_richness_g5": int(final_summary["post_migration_richness"]),
            "migration_richness_change_g5": int(
                final_summary["post_migration_richness"]
            )
            - int(final_summary["post_return_richness"]),
            "realized_host_feedback_g5": float(final_summary["realized_host_feedback"]),
            "expected_host_feedback_after_migration_g5": float(
                final_summary["expected_host_feedback_after_migration"]
            ),
            "realized_migration_fraction": float(
                final_summary["realized_migration_fraction"]
            ),
            "elapsed_minutes": float(execution["elapsed_seconds"]) / 60.0,
            "output_mib": float(execution["output_bytes"]) / 1024**2,
            "peak_rss_mib": float(peak_rss_kib) / 1024.0,
        }
        endpoints.append(endpoint)

        trajectories.append(
            {
                "cell_id": cell_id,
                "cell": cell_id.rsplit("-", 1)[-1],
                "seed_block_id": run["seed_block_id"],
                "generation": 0,
                "environment_richness": int(initial_metrics["richness"]),
                "environment_gene_diversity": initial_metrics["gene_diversity"],
                "post_return_richness": "",
                "post_migration_richness": "",
                "realized_host_feedback": "",
                "expected_host_feedback_after_migration": "",
                "migration_replacement_count": "",
                "realized_migration_fraction": "",
                "mean_adult_richness": "",
                "mean_adult_gene_diversity": "",
            }
        )
        for summary in summary_rows:
            trajectories.append(
                {
                    "cell_id": cell_id,
                    "cell": cell_id.rsplit("-", 1)[-1],
                    "seed_block_id": run["seed_block_id"],
                    "generation": int(summary["host_generation"]),
                    "environment_richness": int(summary["environment_richness"]),
                    "environment_gene_diversity": float(
                        summary["environment_gene_diversity"]
                    ),
                    "post_return_richness": int(summary["post_return_richness"]),
                    "post_migration_richness": int(summary["post_migration_richness"]),
                    "realized_host_feedback": float(summary["realized_host_feedback"]),
                    "expected_host_feedback_after_migration": float(
                        summary["expected_host_feedback_after_migration"]
                    ),
                    "migration_replacement_count": int(
                        summary["migration_replacement_count"]
                    ),
                    "realized_migration_fraction": float(
                        summary["realized_migration_fraction"]
                    ),
                    "mean_adult_richness": float(summary["mean_adult_richness"]),
                    "mean_adult_gene_diversity": float(
                        summary["mean_adult_gene_diversity"]
                    ),
                }
            )

    if audit_issues:
        _atomic_json(audit_path, _audit_payload("FAIL", audit_issues))
        raise SystemExit("v2.1 analysis audit failed:\n" + "\n".join(audit_issues))

    endpoints.sort(
        key=lambda row: (
            CELL_ORDER.index(row["cell_id"]),
            SEED_ORDER.index(row["seed_block_id"]),
        )
    )
    trajectories.sort(
        key=lambda row: (
            CELL_ORDER.index(row["cell_id"]),
            SEED_ORDER.index(row["seed_block_id"]),
            row["generation"],
        )
    )
    write_tsv(derived_dir / "run-endpoints.tsv", endpoints, list(endpoints[0]))
    write_tsv(
        derived_dir / "environment-trajectories.tsv",
        trajectories,
        list(trajectories[0]),
    )

    summary_measures = (
        "richness",
        "shannon",
        "simpson",
        "hill_q1",
        "hill_q2",
        "evenness",
        "total_variation",
        "original_strains_lost",
        "final_mutant_richness",
        "mutant_abundance_fraction",
        "root_richness",
        "root_total_variation",
        "infection_richness_g5",
        "pooled_adult_richness_g5",
        "release_richness_g5",
        "post_return_richness_g5",
        "post_migration_richness_g5",
        "migration_richness_change_g5",
        "migration_sample_tv_g5",
        "realized_host_feedback_g5",
        "expected_host_feedback_after_migration_g5",
        "realized_migration_fraction",
        "elapsed_minutes",
        "output_mib",
        "peak_rss_mib",
    )
    cell_summaries: list[dict[str, Any]] = []
    for cell_id in CELL_ORDER:
        rows = [row for row in endpoints if row["cell_id"] == cell_id]
        cell = cells[cell_id]
        cell_summary: dict[str, Any] = {
            "cell_id": cell_id,
            "cell": cell_id.rsplit("-", 1)[-1],
            "label": cell["label"],
            "H": cell["H"],
            "f": cell["f"],
            "R": cell["R"],
            "alpha": cell["alpha"],
            "u": cell["u"],
            "m": cell["m"],
        }
        for measure in summary_measures:
            for statistic, value in describe(
                [float(row[measure]) for row in rows]
            ).items():
                cell_summary[f"{measure}_{statistic}"] = value
        cell_summaries.append(cell_summary)
    write_tsv(
        derived_dir / "cell-summaries.tsv",
        cell_summaries,
        list(cell_summaries[0]),
    )

    resources_by_seed: dict[str, dict[str, float | int]] = {}
    for seed_id in SEED_ORDER:
        rows = [row for row in endpoints if row["seed_block_id"] == seed_id]
        resources_by_seed[seed_id] = {
            "populations": len(rows),
            "elapsed_minutes": sum(row["elapsed_minutes"] for row in rows),
            "output_mib": sum(row["output_mib"] for row in rows),
            "peak_rss_mib": max(row["peak_rss_mib"] for row in rows),
        }
    analysis_summary = {
        "analysis_schema_version": "1.1.0",
        "experiment_id": EXPERIMENT_ID,
        "variant_id": VARIANT_TAG,
        "pilot_tier": "all",
        "analysis_scope": "exploratory feasibility and calibration pilot",
        "raw_scratch": str(scratch),
        "populations": len(endpoints),
        "cells": len(cells),
        "seed_blocks": list(SEED_ORDER),
        "model_spec_version": MODEL_SPEC_VERSION,
        "output_schema_version": OUTPUT_SCHEMA_VERSION,
        "software_version": SOFTWARE_VERSION,
        "software_git_commits": sorted(commits),
        "source_sha256": sorted(source_hashes),
        "benchmark_platforms": sorted(platforms),
        "benchmark_hardware_note": "HPC run: " + "; ".join(sorted(platforms)),
        "initial_metrics": initial_metrics,
        "migration": {
            "mode": "fixed_regional_pool",
            "fraction_m": MIGRATION_FRACTION,
            "replacement_count_per_passage": MIGRATION_REPLACEMENT_COUNT,
            "regional_pool_id": "rp001-fisher100-fixed",
            "regional_composition": "identical to initial ip001-fisher100",
            "regional_pool_depleted": False,
            "stage": "after host return and environmental capacity regulation",
        },
        "resources_by_seed_block": resources_by_seed,
        "overall_resources": {
            "elapsed_hours": sum(row["elapsed_minutes"] for row in endpoints) / 60.0,
            "output_mib": sum(row["output_mib"] for row in endpoints),
            "peak_rss_mib": max(row["peak_rss_mib"] for row in endpoints),
        },
        "mutation_materialization": {
            "configured_limit_per_transition": 100000,
            "largest_realized_count": max(
                row["maximum_mutants_in_one_transition"] for row in endpoints
            ),
            "largest_fraction_of_limit": max(
                row["maximum_mutants_in_one_transition"] for row in endpoints
            )
            / 100000,
        },
        "audit_status": "PASS",
        "audit_issues": [],
        "limitations": [
            "Three independent populations per cell support descriptive pilot "
            "inference only.",
            "Five host passages do not test equilibrium.",
            "Migration is fixed at m=0.1, so this pilot does not estimate the "
            "effect of m.",
            "The regional pool begins with the focal composition; alternative "
            "regional compositions are not tested.",
            "Complete labelled environmental composition is retained only at "
            "generation 5.",
        ],
    }
    _atomic_json(derived_dir / "analysis-summary.json", analysis_summary)
    _atomic_json(audit_path, _audit_payload("PASS", []))
    print(f"Analysed {len(endpoints)} populations across {len(cells)} cells")
    print(f"Derived tables: {derived_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
