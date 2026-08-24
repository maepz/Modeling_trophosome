#!/usr/bin/env python3
"""Generate or verify the model-2.1 fixed-pool rerun of the first pilot.

This is an isolated variant of the completed 20-cell Phase 1 pilot.  It keeps
the original biological cell matrix and matched seeds, enables fixed-regional-
pool migration at m=0.1, and writes to version-specific config, manifest and
scratch paths so that historical results cannot be overwritten.
"""

from __future__ import annotations

import argparse
import json
from dataclasses import asdict
from pathlib import Path

from prepare_phase1_first_pilot import (
    CAPACITY_RATIO,
    CELLS,
    ENVIRONMENT_CAPACITY,
    GROWTH_FACTOR,
    HOST_GENERATIONS,
    SEED_BLOCKS,
    STEADY_GENERATIONS,
    B,
    K,
    PilotCell,
    _array,
    _sha256,
    _tsv_text,
)

MODEL_SPEC_VERSION = "2.1.0"
OUTPUT_SCHEMA_VERSION = "2.3.0"
SOFTWARE_VERSION = "0.7.0"
MIGRATION_FRACTION = "0.1"
MIGRATION_REPLACEMENT_COUNT = 100_000_000
VARIANT_TAG = "v210-m010"
EXPERIMENT_ID = f"phase1-first-pilot-{VARIANT_TAG}"
REGIONAL_POOL_ID = "rp001-fisher100-fixed"


def _config_text(
    cell: PilotCell,
    counts: list[int],
    checksum: str,
    *,
    seed: int,
    run_id: str | None = None,
) -> str:
    fitness = [1.0] * len(counts)
    host_mode = "full" if cell.hosts <= 100 else "panel"
    run_comment = "" if run_id is None else f"# Run: {run_id}\n"
    return f'''# Phase 1 first-pilot rerun: {cell.cell_id} ({VARIANT_TAG})
# {cell.label}
{run_comment}# Scientific model {MODEL_SPEC_VERSION}; fixed regional-pool migration.
# Initial focal and regional populations use ip001-fisher100 ({checksum}).
# The regional source is fixed and is not depleted by immigration.
model = "wright_fisher_counts"
seed = {seed}
replicates = 1

[environment]
initial_counts = {_array(counts)}
initial_within_host_fitness = {_array(fitness)}
initial_free_living_fitness = {_array(fitness)}
sampling_mode = "reservoir"
capacity_ratio = {CAPACITY_RATIO}

[migration]
mode = "fixed_regional_pool"
fraction = {MIGRATION_FRACTION}
regional_counts = {_array(counts)}

[host]
population_size = {cell.hosts}
infection_bottleneck = {B}
carrying_capacity = {K}
growth_factor = {GROWTH_FACTOR}
steady_generations = {STEADY_GENERATIONS}
host_generations = {HOST_GENERATIONS}
escape_fraction = {cell.escape_fraction}

[evolution]
mutation_probability = {cell.mutation_probability}
mutation_effect_mean = 0.0
mutation_effect_sd = 0.0
within_host_selection = false
free_living_selection = false
fitness_floor = 0.0
max_materialized_mutants = 100000

[output]
snapshot_interval = 100
checkpoint_interval = "1h"
checkpoint_keep = 2
retain_host_histories = false
environment_counts_mode = "final"
host_counts_mode = "{host_mode}"
host_panel_size = 100

[execution]
workers = 2
host_batch_size = 8
in_flight_batches_per_worker = 1
'''


def _run_id(cell: PilotCell, seed_block_id: str) -> str:
    return f"{cell.cell_id}-{VARIANT_TAG}-{seed_block_id}"


def _matrix_row(cell: PilotCell) -> dict[str, object]:
    expected_post_migration_feedback = cell.feedback_alpha * (
        1.0 - float(MIGRATION_FRACTION)
    )
    return {
        "cell_id": cell.cell_id,
        "variant_id": VARIANT_TAG,
        "label": cell.label,
        "experimental_group": cell.group,
        "comparison_sets": "|".join(cell.comparison_sets),
        "pilot_tier": cell.pilot_tier,
        "H": cell.hosts,
        "f": cell.escape_fraction,
        "e": cell.escape_cells_per_host,
        "R": cell.total_return,
        "alpha": format(cell.feedback_alpha, ".12g"),
        "u": cell.mutation_probability,
        "migration_mode": "fixed_regional_pool",
        "m": MIGRATION_FRACTION,
        "migration_replacement_count": MIGRATION_REPLACEMENT_COUNT,
        "expected_host_feedback_after_migration": format(
            expected_post_migration_feedback, ".12g"
        ),
        "regional_pool_id": REGIONAL_POOL_ID,
        "host_counts_mode": "full" if cell.hosts <= 100 else "panel",
        "model_spec_version": MODEL_SPEC_VERSION,
        "status": "planned",
    }


def build_files(repository: Path) -> dict[Path, str]:
    work = repository / "experiments" / "work" / "trophosome"
    population_path = work / "common" / "initial-populations" / "ip001-fisher100.json"
    population = json.loads(population_path.read_text(encoding="utf-8"))
    counts = [int(value) for value in population["scaled_counts"]]
    count_checksum = str(population["scaled_counts_sha256"])
    if len(counts) != 100 or sum(counts) != ENVIRONMENT_CAPACITY:
        raise ValueError("ip001-fisher100 does not contain 100 counts summing to N_E")

    phase = work / "p01-neutral-feedback"
    config_directory = phase / "configs" / f"s01-pilot-{VARIANT_TAG}"
    design_directory = phase / "design"
    manifest_directory = phase / "manifests"
    files: dict[Path, str] = {}
    config_records: list[dict[str, object]] = []
    run_records: list[dict[str, object]] = []
    matrix_rows = [_matrix_row(cell) for cell in CELLS]
    seed_rows = [
        {
            "seed_block_id": seed_block_id,
            "master_seed": master_seed,
            "within_run_replicate_index": 0,
            "use": "matched exploratory replicate; reused from the original pilot",
        }
        for seed_block_id, master_seed in SEED_BLOCKS
    ]

    for cell in CELLS:
        config_path = config_directory / f"{cell.cell_id}-{VARIANT_TAG}.toml"
        config_text = _config_text(
            cell,
            counts,
            count_checksum,
            seed=SEED_BLOCKS[0][1],
        )
        files[config_path] = config_text
        config_records.append(
            {
                "cell_id": cell.cell_id,
                "variant_id": VARIANT_TAG,
                "config_path": str(config_path.relative_to(work)),
                "config_sha256": _sha256(config_text),
            }
        )
        for seed_block_id, master_seed in SEED_BLOCKS:
            run_id = _run_id(cell, seed_block_id)
            run_config_path = config_directory / "runs" / f"{run_id}.toml"
            run_text = _config_text(
                cell,
                counts,
                count_checksum,
                seed=master_seed,
                run_id=run_id,
            )
            files[run_config_path] = run_text
            run_records.append(
                {
                    "experiment_id": EXPERIMENT_ID,
                    "run_id": run_id,
                    "cell_id": cell.cell_id,
                    "variant_id": VARIANT_TAG,
                    "pilot_tier": cell.pilot_tier,
                    "seed_block_id": seed_block_id,
                    "master_seed": master_seed,
                    "within_run_replicate_index": 0,
                    "config_path": str(run_config_path.relative_to(work)),
                    "config_sha256": _sha256(run_text),
                    "scratch_relative_path": str(
                        Path("p01-neutral-feedback")
                        / f"s01-pilot-{VARIANT_TAG}"
                        / cell.cell_id
                        / seed_block_id
                    ),
                    "status": "prepared",
                }
            )

    matrix_fields = list(matrix_rows[0])
    design_stem = f"phase1-first-pilot-{VARIANT_TAG}"
    files[design_directory / f"{design_stem}-cells.tsv"] = _tsv_text(
        matrix_rows, matrix_fields
    )
    files[design_directory / f"{design_stem}-core-cells.tsv"] = _tsv_text(
        [row for row in matrix_rows if row["pilot_tier"] == "core"], matrix_fields
    )
    files[design_directory / f"{design_stem}-extension-cells.tsv"] = _tsv_text(
        [row for row in matrix_rows if row["pilot_tier"] == "extension"],
        matrix_fields,
    )

    files[manifest_directory / f"{design_stem}-seed-blocks.tsv"] = _tsv_text(
        seed_rows,
        ["seed_block_id", "master_seed", "within_run_replicate_index", "use"],
    )
    run_fields = list(run_records[0])
    files[manifest_directory / f"{design_stem}-runs.tsv"] = _tsv_text(
        run_records, run_fields
    )

    manifest = {
        "experiment_manifest_schema_version": "1.1.0",
        "experiment_id": EXPERIMENT_ID,
        "design_id": "phase1-first-pilot-20cell",
        "variant_id": VARIANT_TAG,
        "execution_strategy": (
            "run all 20 cells with three matched seed blocks; the core and "
            "extension labels are retained for monitoring; after all 60 runs pass "
            "the completion audit, automatically knit the migration-aware report"
        ),
        "status": "prepared-not-launched",
        "confirmatory": False,
        "model_family": "wright_fisher_counts",
        "model_spec_version": MODEL_SPEC_VERSION,
        "output_schema_version": OUTPUT_SCHEMA_VERSION,
        "software_version": SOFTWARE_VERSION,
        "initial_population_id": population["initial_population_id"],
        "initial_counts_sha256": count_checksum,
        "regional_pool": {
            "regional_pool_id": REGIONAL_POOL_ID,
            "mode": "fixed_regional_pool",
            "composition": "identical to ip001-fisher100",
            "counts_sha256": count_checksum,
            "depleted_by_immigration": False,
            "fraction_m": float(MIGRATION_FRACTION),
            "replacement_count_per_passage": MIGRATION_REPLACEMENT_COUNT,
        },
        "neutral_profiles": {
            "architecture_profile_id": "arch-fixed-regional-pool-v1",
            "selection_profile_id": "sel-neutral-v1",
            "fitness_profile_id": "fit-neutral-v1",
        },
        "fixed_parameters": {
            "B": B,
            "K": K,
            "growth_factor": float(GROWTH_FACTOR),
            "steady_generations": STEADY_GENERATIONS,
            "host_generations": HOST_GENERATIONS,
            "c": float(CAPACITY_RATIO),
            "N_E": ENVIRONMENT_CAPACITY,
            "sampling_mode": "reservoir",
            "migration_mode": "fixed_regional_pool",
            "m": float(MIGRATION_FRACTION),
            "regional_pool_id": REGIONAL_POOL_ID,
            "within_host_selection": False,
            "free_living_selection": False,
            "mutation_effect_mean": 0.0,
            "mutation_effect_sd": 0.0,
            "environment_counts_mode": "final",
            "replicates_per_run": 1,
            "workers_per_population": 2,
        },
        "seed_blocks": seed_rows,
        "runs": run_records,
        "cells": [
            {
                **asdict(cell),
                "cell_id": cell.cell_id,
                "pilot_tier": cell.pilot_tier,
                "escape_cells_per_host": cell.escape_cells_per_host,
                "total_return": cell.total_return,
                "feedback_alpha": cell.feedback_alpha,
                "config_path": record["config_path"],
                "config_sha256": record["config_sha256"],
            }
            for cell, record in zip(CELLS, config_records, strict=True)
        ],
        "source_freeze": {
            "required_before_launch": True,
            "note": (
                "Transfer one frozen model-2.1.0 source revision to the HPC "
                "before launching the first population."
            ),
        },
        "automatic_reporting": {
            "completion_gate": "all 60 populations complete and internally valid",
            "analysis_directory": (
                f"p01-neutral-feedback/analysis/s01-pilot-{VARIANT_TAG}-derived"
            ),
            "pdf": f"output/pdf/phase1-first-pilot-{VARIANT_TAG}-report.pdf",
            "markdown": f"docs/phase1-first-pilot-{VARIANT_TAG}-report.md",
            "migration_aware": True,
        },
    }
    files[manifest_directory / f"{design_stem}-manifest.json"] = (
        json.dumps(manifest, indent=2) + "\n"
    )
    return files


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

    repository = args.repository.resolve()
    files = build_files(repository)
    if args.verify:
        mismatches = [
            str(path.relative_to(repository))
            for path, expected in files.items()
            if not path.is_file() or path.read_text(encoding="utf-8") != expected
        ]
        if mismatches:
            raise SystemExit("model-2.1 pilot files differ:\n" + "\n".join(mismatches))
        print(f"Verified {len(files)} model-2.1 first-pilot files")
        return 0

    for path, content in files.items():
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(content, encoding="utf-8")
    print(f"Wrote {len(files)} model-2.1 first-pilot files")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
