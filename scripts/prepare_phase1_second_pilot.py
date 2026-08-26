#!/usr/bin/env python3
"""Generate or verify the Phase 1 equilibrium-and-precision second pilot."""

from __future__ import annotations

import argparse
import csv
import json
import math
import statistics
from dataclasses import asdict, dataclass
from decimal import Decimal
from pathlib import Path

from prepare_phase1_first_pilot import (
    CAPACITY_RATIO,
    ENVIRONMENT_CAPACITY,
    GROWTH_FACTOR,
    STEADY_GENERATIONS,
    B,
    K,
    _array,
    _sha256,
    _tsv_text,
)

MODEL_SPEC_VERSION = "2.1.0"
OUTPUT_SCHEMA_VERSION = "2.3.0"
SOFTWARE_VERSION = "0.7.0"
MIGRATION_FRACTION = "0.1"
MIGRATION_REPLACEMENT_COUNT = 100_000_000
HOST_GENERATIONS = 250
VARIANT_TAG = "v210-m010-g250"
EXPERIMENT_ID = f"phase1-second-pilot-{VARIANT_TAG}"
DESIGN_STEM = f"phase1-second-pilot-{VARIANT_TAG}"
STAGE_DIRECTORY = f"s02-equilibrium-precision-{VARIANT_TAG}"
REGIONAL_POOL_ID = "rp001-fisher100-fixed"
SEED_BLOCKS = tuple((f"sb{number:04d}", 665 + number) for number in range(1, 13))


@dataclass(frozen=True)
class SentinelCell:
    number: int
    source_first_pilot_cell: str
    label: str
    sentinel_role: str
    group: str
    comparison_sets: tuple[str, ...]
    hosts: int
    escape_fraction: str
    mutation_probability: str

    @property
    def cell_id(self) -> str:
        return f"p01-s02-c{self.number:04d}"

    @property
    def short_id(self) -> str:
        return f"c{self.number:04d}"

    @property
    def escape_cells_per_host(self) -> int:
        return int(Decimal(self.escape_fraction) * K)

    @property
    def total_return(self) -> int:
        return self.hosts * self.escape_cells_per_host

    @property
    def feedback_alpha(self) -> float:
        return self.total_return / (ENVIRONMENT_CAPACITY + self.total_return)

    @property
    def feedback_exposure_rate(self) -> float:
        if self.feedback_alpha <= 0:
            return 0.0
        return -math.log1p(-self.feedback_alpha)

    @property
    def diagnostic_window_generations(self) -> int:
        if self.feedback_exposure_rate == 0:
            return 50
        return max(20, math.ceil(5 / self.feedback_exposure_rate))

    @property
    def host_counts_mode(self) -> str:
        return "full" if self.hosts <= 100 else "panel"


CELLS = (
    SentinelCell(
        21,
        "p01-s01-c0001",
        "Migration-only no-return control",
        "no return",
        "EQ-NR",
        ("host-passage", "migration-baseline"),
        100,
        "0",
        "0",
    ),
    SentinelCell(
        22,
        "p01-s01-c0003",
        "Mutation-free baseline feedback, H=100",
        "baseline host passage",
        "EQ-B",
        ("host-passage", "mutation", "pooling", "feedback"),
        100,
        "1e-2",
        "0",
    ),
    SentinelCell(
        23,
        "p01-s01-c0009",
        "Mutation-enabled baseline feedback, H=100",
        "within-host mutation",
        "EQ-M",
        ("mutation",),
        100,
        "1e-2",
        "1e-10",
    ),
    SentinelCell(
        24,
        "p01-s01-c0005",
        "Many-host fixed-return extreme, H=10000",
        "host pooling",
        "EQ-H",
        ("pooling",),
        10_000,
        "1e-4",
        "0",
    ),
    SentinelCell(
        25,
        "p01-s01-c0016",
        "Weak feedback, alpha about 0.091",
        "weak feedback",
        "EQ-F",
        ("feedback",),
        100,
        "1e-3",
        "0",
    ),
    SentinelCell(
        26,
        "p01-s01-c0018",
        "Strong feedback, alpha about 0.909",
        "strong feedback",
        "EQ-F",
        ("feedback", "escape-range-sensitivity"),
        100,
        "1e-1",
        "0",
    ),
)


def _run_id(cell: SentinelCell, seed_block_id: str) -> str:
    return f"{cell.cell_id}-{VARIANT_TAG}-{seed_block_id}"


def _config_text(
    cell: SentinelCell,
    counts: list[int],
    checksum: str,
    *,
    seed: int,
    run_id: str | None = None,
) -> str:
    fitness = [1.0] * len(counts)
    run_comment = "" if run_id is None else f"# Run: {run_id}\n"
    return f'''# Phase 1 second pilot: {cell.cell_id} ({VARIANT_TAG})
# {cell.label}; derived from {cell.source_first_pilot_cell}.
{run_comment}# Equilibrium-and-precision screen; 250 fixed host passages.
# Scientific model {MODEL_SPEC_VERSION}; fixed regional-pool migration.
# Initial focal and regional populations use ip001-fisher100 ({checksum}).
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
environment_counts_mode = "all"
host_counts_mode = "{cell.host_counts_mode}"
host_panel_size = 100

[execution]
workers = 2
host_batch_size = 8
in_flight_batches_per_worker = 1
'''


def _first_pilot_costs(repository: Path) -> dict[str, dict[str, float]]:
    path = (
        repository
        / "experiments/work/trophosome/p01-neutral-feedback/analysis"
        / "s01-pilot-v210-m010-derived/run-endpoints.tsv"
    )
    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    costs: dict[str, dict[str, float]] = {}
    for source_cell in {cell.source_first_pilot_cell for cell in CELLS}:
        short_id = source_cell.rsplit("-", 1)[-1]
        selected = [row for row in rows if row["cell"] == short_id]
        if len(selected) != 3:
            raise ValueError(f"expected three first-pilot results for {source_cell}")
        costs[source_cell] = {
            "runtime_minutes_5g": statistics.median(
                float(row["elapsed_minutes"]) for row in selected
            ),
            "output_mib_5g": statistics.median(
                float(row["output_mib"]) for row in selected
            ),
            "peak_rss_mib": statistics.median(
                float(row["peak_rss_mib"]) for row in selected
            ),
        }
    return costs


def _matrix_row(
    cell: SentinelCell, costs: dict[str, dict[str, float]]
) -> dict[str, object]:
    source_cost = costs[cell.source_first_pilot_cell]
    scale = HOST_GENERATIONS / 5
    return {
        "cell_id": cell.cell_id,
        "cell": cell.short_id,
        "source_first_pilot_cell": cell.source_first_pilot_cell,
        "label": cell.label,
        "sentinel_role": cell.sentinel_role,
        "experimental_group": cell.group,
        "comparison_sets": "|".join(cell.comparison_sets),
        "H": cell.hosts,
        "f": cell.escape_fraction,
        "e": cell.escape_cells_per_host,
        "R": cell.total_return,
        "alpha": format(cell.feedback_alpha, ".12g"),
        "u": cell.mutation_probability,
        "m": MIGRATION_FRACTION,
        "host_generations": HOST_GENERATIONS,
        "feedback_exposure_at_horizon": format(
            HOST_GENERATIONS * cell.feedback_exposure_rate, ".6g"
        ),
        "diagnostic_window_generations": cell.diagnostic_window_generations,
        "diagnostic_windows": 4,
        "host_counts_mode": cell.host_counts_mode,
        "environment_counts_mode": "all",
        "first_pilot_runtime_minutes_5g": format(
            source_cost["runtime_minutes_5g"], ".8g"
        ),
        "projected_runtime_hours_per_run": format(
            source_cost["runtime_minutes_5g"] * scale / 60, ".6g"
        ),
        "projected_output_gib_per_run": format(
            source_cost["output_mib_5g"] * scale / 1024, ".6g"
        ),
        "first_pilot_peak_rss_mib": format(source_cost["peak_rss_mib"], ".6g"),
        "model_spec_version": MODEL_SPEC_VERSION,
        "status": "planned",
    }


def _registry_cell_row(cell: SentinelCell) -> dict[str, str]:
    mnemonic = (
        f"h{cell.hosts}-f{cell.escape_fraction.replace('-', 'm')}-"
        f"u{cell.mutation_probability.replace('-', 'm')}-c1-b10-k{K}-"
        f"ts{STEADY_GENERATIONS}-hg{HOST_GENERATIONS}-archfixedpool-"
        "selneutral-fitneutral"
    )
    return {
        "cell_id": cell.cell_id,
        "phase_id": "p01",
        "stage_id": "s02",
        "label": cell.label,
        "mnemonic": mnemonic,
        "cell_dirname": cell.cell_id,
        "experimental_group": cell.group,
        "comparison_set": "|".join(cell.comparison_sets),
        "confirmatory": "false",
        "architecture_profile_id": "arch-fixed-regional-pool-v1",
        "selection_profile_id": "sel-neutral-v1",
        "fitness_profile_id": "fit-neutral-v1",
        "initial_population_id": "ip001-fisher100",
        "config_path": (
            f"p01-neutral-feedback/configs/{STAGE_DIRECTORY}/"
            f"{cell.cell_id}-{VARIANT_TAG}.toml"
        ),
        "status": "prepared",
        "notes": (
            "Exploratory equilibrium-and-precision sentinel derived from "
            f"{cell.source_first_pilot_cell}; 12 seed blocks and 250 passages."
        ),
    }


def _parameter_rows(cell: SentinelCell) -> list[dict[str, str]]:
    values = (
        ("H", str(cell.hosts), "integer", "hosts", "input"),
        ("f", cell.escape_fraction, "float", "fraction", "input"),
        ("e", str(cell.escape_cells_per_host), "integer", "cells_per_host", "derived"),
        ("R", str(cell.total_return), "integer", "cells", "derived"),
        ("alpha", format(cell.feedback_alpha, ".12g"), "float", "fraction", "derived"),
        (
            "u",
            cell.mutation_probability,
            "float",
            "probability_per_genome_per_bacterial_generation",
            "input",
        ),
        ("m", MIGRATION_FRACTION, "float", "fraction", "input"),
        ("regional_pool_id", REGIONAL_POOL_ID, "string", "", "nuisance"),
        ("host_counts_mode", cell.host_counts_mode, "string", "", "technical"),
        ("host_panel_size", "100", "integer", "hosts", "technical"),
        ("B", str(B), "integer", "cells_per_infection", "nuisance"),
        ("K", str(K), "integer", "cells_per_host", "nuisance"),
        ("growth_factor", GROWTH_FACTOR, "float", "ratio", "nuisance"),
        (
            "steady_generations",
            str(STEADY_GENERATIONS),
            "integer",
            "bacterial_generations",
            "nuisance",
        ),
        (
            "host_generations",
            str(HOST_GENERATIONS),
            "integer",
            "host_passages",
            "technical",
        ),
        ("c", CAPACITY_RATIO, "float", "ratio", "nuisance"),
        ("N_E", str(ENVIRONMENT_CAPACITY), "integer", "cells", "derived"),
        ("sampling_mode", "reservoir", "string", "", "nuisance"),
        ("within_host_selection", "false", "boolean", "", "nuisance"),
        ("free_living_selection", "false", "boolean", "", "nuisance"),
        ("mutation_effect_mean", "0.0", "float", "fitness_units", "nuisance"),
        ("mutation_effect_sd", "0.0", "float", "fitness_units", "nuisance"),
        ("environment_counts_mode", "all", "string", "", "technical"),
        ("planned_seed_blocks", "12", "integer", "seed_blocks", "technical"),
        (
            "seed_block_set",
            "phase1-second-pilot-sb0001-sb0012",
            "string",
            "",
            "technical",
        ),
        (
            "source_first_pilot_cell",
            cell.source_first_pilot_cell,
            "string",
            "",
            "design",
        ),
        (
            "diagnostic_window_generations",
            str(cell.diagnostic_window_generations),
            "integer",
            "host_passages",
            "analysis",
        ),
    )
    return [
        {
            "cell_id": cell.cell_id,
            "parameter_name": name,
            "value": value,
            "value_type": value_type,
            "unit": unit,
            "role": role,
        }
        for name, value, value_type, unit, role in values
    ]


def build_files(repository: Path) -> dict[Path, str]:
    work = repository / "experiments" / "work" / "trophosome"
    population_path = work / "common/initial-populations/ip001-fisher100.json"
    population = json.loads(population_path.read_text(encoding="utf-8"))
    counts = [int(value) for value in population["scaled_counts"]]
    count_checksum = str(population["scaled_counts_sha256"])
    if len(counts) != 100 or sum(counts) != ENVIRONMENT_CAPACITY:
        raise ValueError("ip001-fisher100 does not contain 100 counts summing to N_E")

    phase = work / "p01-neutral-feedback"
    config_directory = phase / "configs" / STAGE_DIRECTORY
    design_directory = phase / "design"
    manifest_directory = phase / "manifests"
    costs = _first_pilot_costs(repository)
    matrix_rows = [_matrix_row(cell, costs) for cell in CELLS]
    files: dict[Path, str] = {}
    run_records: list[dict[str, object]] = []
    config_records: list[dict[str, object]] = []

    for cell in CELLS:
        base_path = config_directory / f"{cell.cell_id}-{VARIANT_TAG}.toml"
        base_text = _config_text(cell, counts, count_checksum, seed=SEED_BLOCKS[0][1])
        files[base_path] = base_text
        config_records.append(
            {
                "cell_id": cell.cell_id,
                "config_path": str(base_path.relative_to(work)),
                "config_sha256": _sha256(base_text),
            }
        )
        for seed_block_id, master_seed in SEED_BLOCKS:
            run_id = _run_id(cell, seed_block_id)
            run_path = config_directory / "runs" / f"{run_id}.toml"
            run_text = _config_text(
                cell,
                counts,
                count_checksum,
                seed=master_seed,
                run_id=run_id,
            )
            files[run_path] = run_text
            run_records.append(
                {
                    "experiment_id": EXPERIMENT_ID,
                    "run_id": run_id,
                    "cell_id": cell.cell_id,
                    "cell": cell.short_id,
                    "variant_id": VARIANT_TAG,
                    "seed_block_id": seed_block_id,
                    "master_seed": master_seed,
                    "within_run_replicate_index": 0,
                    "config_path": str(run_path.relative_to(work)),
                    "config_sha256": _sha256(run_text),
                    "scratch_relative_path": str(
                        Path("p01-neutral-feedback")
                        / STAGE_DIRECTORY
                        / cell.cell_id
                        / seed_block_id
                    ),
                    "status": "prepared",
                }
            )

    files[design_directory / f"{DESIGN_STEM}-cells.tsv"] = _tsv_text(
        matrix_rows, list(matrix_rows[0])
    )
    seed_rows = [
        {
            "seed_block_id": seed_block_id,
            "master_seed": master_seed,
            "within_run_replicate_index": 0,
            "use": (
                "matched exploratory precision replicate; first three continue "
                "the seed series and nine are new"
            ),
        }
        for seed_block_id, master_seed in SEED_BLOCKS
    ]
    files[manifest_directory / f"{DESIGN_STEM}-seed-blocks.tsv"] = _tsv_text(
        seed_rows, list(seed_rows[0])
    )
    files[manifest_directory / f"{DESIGN_STEM}-runs.tsv"] = _tsv_text(
        run_records, list(run_records[0])
    )

    projected_runtime_hours = sum(
        float(row["projected_runtime_hours_per_run"]) * len(SEED_BLOCKS)
        for row in matrix_rows
    )
    projected_output_gib = sum(
        float(row["projected_output_gib_per_run"]) * len(SEED_BLOCKS)
        for row in matrix_rows
    )
    manifest = {
        "experiment_manifest_schema_version": "1.2.0",
        "experiment_id": EXPERIMENT_ID,
        "design_id": "phase1-equilibrium-precision-six-sentinel",
        "variant_id": VARIANT_TAG,
        "stage": "s02-equilibrium-precision",
        "status": "prepared-not-launched",
        "confirmatory": False,
        "purpose": (
            "estimate equilibrium-screen behavior, continuing fluctuations, "
            "required host-passage horizon and confirmatory replicate precision"
        ),
        "model_family": "wright_fisher_counts",
        "model_spec_version": MODEL_SPEC_VERSION,
        "output_schema_version": OUTPUT_SCHEMA_VERSION,
        "software_version": SOFTWARE_VERSION,
        "initial_population_id": population["initial_population_id"],
        "initial_counts_sha256": count_checksum,
        "sentinel_selection_source": (
            "phase1-first-pilot-v210-m010 analysis committed in 0998f87"
        ),
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
            "environment_counts_mode": "all",
            "replicates_per_run": 1,
            "workers_per_population": 2,
        },
        "generation_rule": {
            "formula": "min(5000, max(250, ceil(20 / -log(1-alpha))))",
            "selected_horizon": HOST_GENERATIONS,
            "note": (
                "All selected positive-feedback sentinels reach the 250-passage "
                "minimum; the no-return control uses the same fixed horizon."
            ),
        },
        "equilibrium_screen": {
            "windows": 4,
            "target_feedback_exposure_per_window": 5,
            "minimum_window_generations": 20,
            "no_return_window_generations": 50,
            "responses": ["D0", "D1", "D2", "evenness", "TV"],
            "equivalence_margins": {
                "D0_relative": 0.05,
                "D1_relative": 0.05,
                "D2_relative": 0.05,
                "evenness_absolute": 0.02,
                "TV_absolute": 0.05,
            },
            "full_equilibrium_limit": (
                "Alternative-initial-condition convergence is not tested in this "
                "stage, so passing is a stationarity screen, not definitive "
                "equilibrium."
            ),
        },
        "resource_projection_from_first_pilot": {
            "method": "linear scaling from 5 to 250 host passages",
            "aggregate_runtime_hours": projected_runtime_hours,
            "aggregate_output_gib": projected_output_gib,
            "runtime_budget_hours_per_run": 48,
            "storage_budget_gib": 100,
            "expansion_fraction_limit": 0.70,
        },
        "seed_blocks": seed_rows,
        "cells": [
            {
                **asdict(cell),
                **row,
                "config_path": config_record["config_path"],
                "config_sha256": config_record["config_sha256"],
            }
            for cell, row, config_record in zip(
                CELLS, matrix_rows, config_records, strict=True
            )
        ],
        "runs": run_records,
        "automatic_reporting": {
            "completion_gate": "all 72 populations complete and internally valid",
            "analysis_directory": (
                f"p01-neutral-feedback/analysis/{STAGE_DIRECTORY}-derived"
            ),
            "pdf": f"output/pdf/{EXPERIMENT_ID}-report.pdf",
            "markdown": f"docs/{EXPERIMENT_ID}-report.md",
            "report_only_option": "--report-only",
            "standalone_builder": "scripts/build_phase1_second_pilot_report.py",
        },
        "source_freeze": {"required_before_launch": True},
    }
    files[manifest_directory / f"{DESIGN_STEM}-manifest.json"] = (
        json.dumps(manifest, indent=2) + "\n"
    )
    return files


def _registry_expectations() -> tuple[list[dict[str, str]], list[dict[str, str]]]:
    return (
        [_registry_cell_row(cell) for cell in CELLS],
        [row for cell in CELLS for row in _parameter_rows(cell)],
    )


def _sync_registry(
    path: Path, expected: list[dict[str, str]], *, write: bool
) -> list[str]:
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        fields = list(reader.fieldnames or ())
        rows = list(reader)
    by_key = {(row["cell_id"], row.get("parameter_name", "")): row for row in rows}
    mismatches: list[str] = []
    additions: list[dict[str, str]] = []
    for row in expected:
        key = (row["cell_id"], row.get("parameter_name", ""))
        current = by_key.get(key)
        if current is None:
            additions.append(row)
        elif current != row:
            mismatches.append(f"registry row differs: {path.name} {key}")
    if write:
        managed_cell_ids = {row["cell_id"] for row in expected}
        preserved = [row for row in rows if row["cell_id"] not in managed_cell_ids]
        with path.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=fields, lineterminator="\n")
            writer.writeheader()
            writer.writerows(preserved + expected)
        return []
    if additions:
        mismatches.extend(
            f"registry row is missing: {path.name} "
            f"{(row['cell_id'], row.get('parameter_name', ''))}"
            for row in additions
        )
    return mismatches


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--repository", type=Path, default=Path(__file__).resolve().parents[1]
    )
    parser.add_argument("--write", action="store_true")
    parser.add_argument("--verify", action="store_true")
    args = parser.parse_args()
    if args.write == args.verify:
        parser.error("choose exactly one of --write or --verify")

    repository = args.repository.resolve()
    files = build_files(repository)
    work = repository / "experiments/work/trophosome"
    expected_cells, expected_parameters = _registry_expectations()
    registry_mismatches = _sync_registry(
        work / "registry/cells.csv", expected_cells, write=args.write
    )
    registry_mismatches.extend(
        _sync_registry(
            work / "registry/cell_parameters.csv",
            expected_parameters,
            write=args.write,
        )
    )
    if registry_mismatches:
        raise SystemExit(
            "second-pilot registry differs:\n" + "\n".join(registry_mismatches)
        )

    if args.verify:
        mismatches = [
            str(path.relative_to(repository))
            for path, expected in files.items()
            if not path.is_file() or path.read_text(encoding="utf-8") != expected
        ]
        if mismatches:
            raise SystemExit("second-pilot files differ:\n" + "\n".join(mismatches))
        print(f"Verified {len(files)} second-pilot files and both registries")
        return 0

    for path, content in files.items():
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(content, encoding="utf-8")
    print(f"Wrote {len(files)} second-pilot files and updated both registries")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
