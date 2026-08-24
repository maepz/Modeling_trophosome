#!/usr/bin/env python3
"""Generate or verify the staged 12-core plus 8-extension pilot files."""

from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
from dataclasses import asdict, dataclass
from decimal import Decimal
from pathlib import Path

SEED_BLOCKS = (("sb0001", 666), ("sb0002", 667), ("sb0003", 668))
K = 1_000_000_000
B = 10
GROWTH_FACTOR = "1.2"
STEADY_GENERATIONS = 500
HOST_GENERATIONS = 5
CAPACITY_RATIO = "1.0"
ENVIRONMENT_CAPACITY = 1_000_000_000


@dataclass(frozen=True)
class PilotCell:
    number: int
    label: str
    group: str
    comparison_sets: tuple[str, ...]
    hosts: int
    escape_fraction: str
    mutation_probability: str
    extension: bool = False

    @property
    def cell_id(self) -> str:
        return f"p01-s01-c{self.number:04d}"

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
    def pilot_tier(self) -> str:
        return "extension" if self.extension else "core"


CELLS = (
    PilotCell(1, "No return, mutation off", "NR", ("no-return",), 100, "0", "0"),
    PilotCell(
        2,
        "No return, within-host mutation on",
        "NRM",
        ("no-return", "mutation-gating"),
        100,
        "0",
        "1e-10",
    ),
    PilotCell(
        3,
        "Matched return R=1e9, H=100",
        "MR0",
        ("matched-r1e9-u0", "mutation-bracket", "feedback-h100"),
        100,
        "1e-2",
        "0",
    ),
    PilotCell(
        4,
        "Matched return R=1e9, H=1000",
        "MR0",
        ("matched-r1e9-u0", "feedback-h1000"),
        1_000,
        "1e-3",
        "0",
    ),
    PilotCell(
        5,
        "Matched return R=1e9, H=10000",
        "MR0",
        ("matched-r1e9-u0",),
        10_000,
        "1e-4",
        "0",
    ),
    PilotCell(
        6,
        "Matched return R=1e9, H=100000",
        "MR0",
        ("matched-r1e9-u0",),
        100_000,
        "1e-5",
        "0",
    ),
    PilotCell(
        7,
        "Mutation bracket u=1e-12",
        "MUT",
        ("mutation-bracket",),
        100,
        "1e-2",
        "1e-12",
    ),
    PilotCell(
        8,
        "Mutation bracket u=1e-11",
        "MUT",
        ("mutation-bracket",),
        100,
        "1e-2",
        "1e-11",
    ),
    PilotCell(
        9,
        "Mutation bracket u=1e-10",
        "MUT",
        ("mutation-bracket", "matched-r1e9-u1e-10"),
        100,
        "1e-2",
        "1e-10",
    ),
    PilotCell(
        10,
        "Mutation bracket u=1e-9",
        "MUT",
        ("mutation-bracket",),
        100,
        "1e-2",
        "1e-9",
    ),
    PilotCell(
        11,
        "Very weak feedback boundary",
        "FB0",
        ("feedback-h100",),
        100,
        "1e-5",
        "0",
    ),
    PilotCell(
        12,
        "Strong feedback boundary",
        "FB0",
        ("feedback-h1000",),
        1_000,
        "1e-2",
        "0",
    ),
    PilotCell(
        13,
        "Matched return R=1e9, H=1000, u=1e-10",
        "MRM",
        ("matched-r1e9-u1e-10",),
        1_000,
        "1e-3",
        "1e-10",
        True,
    ),
    PilotCell(
        14,
        "Matched return R=1e9, H=10000, u=1e-10",
        "MRM",
        ("matched-r1e9-u1e-10",),
        10_000,
        "1e-4",
        "1e-10",
        True,
    ),
    PilotCell(
        15,
        "Matched return R=1e9, H=100000, u=1e-10",
        "MRM",
        ("matched-r1e9-u1e-10",),
        100_000,
        "1e-5",
        "1e-10",
        True,
    ),
    PilotCell(
        16,
        "Matched return R=1e8, H=100",
        "MRLOW",
        ("matched-r1e8-u0", "feedback-h100"),
        100,
        "1e-3",
        "0",
        True,
    ),
    PilotCell(
        17,
        "Matched return R=1e8, H=10000",
        "MRLOW",
        ("matched-r1e8-u0",),
        10_000,
        "1e-5",
        "0",
        True,
    ),
    PilotCell(
        18,
        "Alpha series H=100, alpha=0.909; f=0.1 sensitivity",
        "FBA",
        ("feedback-h100", "escape-range-sensitivity"),
        100,
        "1e-1",
        "0",
        True,
    ),
    PilotCell(
        19,
        "Alpha series H=1000, alpha=0.001; f=1e-6 sensitivity",
        "FBA",
        ("feedback-h1000", "escape-range-sensitivity"),
        1_000,
        "1e-6",
        "0",
        True,
    ),
    PilotCell(
        20,
        "Alpha series H=1000, alpha=0.091",
        "FBA",
        ("feedback-h1000", "matched-r1e8-u0"),
        1_000,
        "1e-4",
        "0",
        True,
    ),
)


def _array(values: list[int] | list[float]) -> str:
    lines = []
    for start in range(0, len(values), 10):
        chunk = ", ".join(str(value) for value in values[start : start + 10])
        lines.append(f"  {chunk},")
    return "[\n" + "\n".join(lines) + "\n]"


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
    return f'''# Phase 1 first pilot: {cell.cell_id}
# {cell.label}
{run_comment}# One run contains one stochastic population. Child random streams
# come from this master seed for every host, passage, and model process.
# Initial population: ip001-fisher100 ({checksum})
model = "wright_fisher_counts"
seed = {seed}
replicates = 1

[environment]
initial_counts = {_array(counts)}
initial_within_host_fitness = {_array(fitness)}
initial_free_living_fitness = {_array(fitness)}
sampling_mode = "reservoir"
capacity_ratio = {CAPACITY_RATIO}

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


def _tsv_text(rows: list[dict[str, object]], fields: list[str]) -> str:
    stream = io.StringIO()
    writer = csv.DictWriter(
        stream,
        fieldnames=fields,
        delimiter="\t",
        lineterminator="\n",
    )
    writer.writeheader()
    writer.writerows(rows)
    return stream.getvalue()


def _csv_text(rows: list[dict[str, object]], fields: list[str]) -> str:
    stream = io.StringIO()
    writer = csv.DictWriter(stream, fieldnames=fields, lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows)
    return stream.getvalue()


def _sha256(text: str) -> str:
    return hashlib.sha256(text.encode()).hexdigest()


def _mnemonic_number(value: str) -> str:
    return value.lower().replace("e-", "em").replace("e+", "e").replace(".", "p")


def _cell_mnemonic(cell: PilotCell) -> str:
    return "-".join(
        (
            f"h{cell.hosts}",
            f"f{_mnemonic_number(cell.escape_fraction)}",
            f"u{_mnemonic_number(cell.mutation_probability)}",
            "c1",
            "b10",
            "k1000000000",
            "ts500",
            "archpanmictic",
            "selneutral",
            "fitneutral",
        )
    )


def _parameter_rows(cell: PilotCell) -> list[dict[str, object]]:
    host_mode = "full" if cell.hosts <= 100 else "panel"
    values = (
        ("H", cell.hosts, "integer", "hosts", "input"),
        (
            "f",
            cell.escape_fraction,
            "float" if "e" in cell.escape_fraction else "integer",
            "fraction",
            "input",
        ),
        ("e", cell.escape_cells_per_host, "integer", "cells_per_host", "derived"),
        ("R", cell.total_return, "integer", "cells", "derived"),
        (
            "alpha",
            format(cell.feedback_alpha, ".12g"),
            "float" if cell.feedback_alpha else "integer",
            "fraction",
            "derived",
        ),
        (
            "u",
            cell.mutation_probability,
            "float" if "e" in cell.mutation_probability else "integer",
            "probability_per_genome_per_bacterial_generation",
            "input",
        ),
        ("host_counts_mode", host_mode, "string", "", "technical"),
        ("B", B, "integer", "cells_per_infection", "nuisance"),
        ("K", K, "integer", "cells_per_host", "nuisance"),
        ("growth_factor", GROWTH_FACTOR, "float", "ratio", "nuisance"),
        (
            "steady_generations",
            STEADY_GENERATIONS,
            "integer",
            "bacterial_generations",
            "nuisance",
        ),
        ("host_generations", HOST_GENERATIONS, "integer", "host_passages", "technical"),
        ("c", CAPACITY_RATIO, "float", "ratio", "nuisance"),
        ("N_E", ENVIRONMENT_CAPACITY, "integer", "cells", "derived"),
        ("sampling_mode", "reservoir", "string", "", "nuisance"),
        ("within_host_selection", "false", "boolean", "", "nuisance"),
        ("free_living_selection", "false", "boolean", "", "nuisance"),
        ("mutation_effect_mean", "0.0", "float", "fitness_units", "nuisance"),
        ("mutation_effect_sd", "0.0", "float", "fitness_units", "nuisance"),
        ("environment_counts_mode", "final", "string", "", "technical"),
        (
            "planned_seed_blocks",
            len(SEED_BLOCKS),
            "integer",
            "seed_blocks",
            "technical",
        ),
        (
            "seed_block_set",
            "phase1-first-pilot-sb0001-sb0003",
            "string",
            "",
            "technical",
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
    population_path = work / "common" / "initial-populations" / "ip001-fisher100.json"
    population = json.loads(population_path.read_text(encoding="utf-8"))
    counts = [int(value) for value in population["scaled_counts"]]
    count_checksum = str(population["scaled_counts_sha256"])
    if len(counts) != 100 or sum(counts) != ENVIRONMENT_CAPACITY:
        raise ValueError("ip001-fisher100 does not contain 100 counts summing to N_E")

    config_directory = work / "p01-neutral-feedback" / "configs" / "s01-pilot"
    design_directory = work / "p01-neutral-feedback" / "design"
    manifest_directory = work / "p01-neutral-feedback" / "manifests"
    files: dict[Path, str] = {}
    config_records = []
    run_records = []
    matrix_rows = []
    seed_rows = [
        {
            "seed_block_id": seed_block_id,
            "master_seed": master_seed,
            "within_run_replicate_index": 0,
            "use": "exploratory first-pilot replicate",
        }
        for seed_block_id, master_seed in SEED_BLOCKS
    ]
    registry_cell_rows: list[dict[str, object]] = []
    registry_parameter_rows: list[dict[str, object]] = []
    for cell in CELLS:
        config_path = config_directory / f"{cell.cell_id}.toml"
        text = _config_text(
            cell,
            counts,
            count_checksum,
            seed=SEED_BLOCKS[0][1],
        )
        files[config_path] = text
        config_records.append(
            {
                "cell_id": cell.cell_id,
                "config_path": str(config_path.relative_to(work)),
                "config_sha256": _sha256(text),
            }
        )
        for seed_block_id, master_seed in SEED_BLOCKS:
            run_id = f"{cell.cell_id}-{seed_block_id}"
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
                    "run_id": run_id,
                    "cell_id": cell.cell_id,
                    "seed_block_id": seed_block_id,
                    "master_seed": master_seed,
                    "within_run_replicate_index": 0,
                    "config_path": str(run_config_path.relative_to(work)),
                    "config_sha256": _sha256(run_text),
                    "scratch_relative_path": str(
                        Path("p01-neutral-feedback")
                        / "s01-pilot"
                        / cell.cell_id
                        / seed_block_id
                    ),
                    "status": "prepared",
                }
            )
        matrix_rows.append(
            {
                "cell_id": cell.cell_id,
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
                "host_counts_mode": "full" if cell.hosts <= 100 else "panel",
                "status": "planned",
            }
        )
        registry_cell_rows.append(
            {
                "cell_id": cell.cell_id,
                "phase_id": "p01",
                "stage_id": "s01",
                "label": cell.label,
                "mnemonic": _cell_mnemonic(cell),
                "cell_dirname": cell.cell_id,
                "experimental_group": cell.group,
                "comparison_set": "|".join(cell.comparison_sets),
                "confirmatory": "false",
                "architecture_profile_id": "arch-panmictic-v1",
                "selection_profile_id": "sel-neutral-v1",
                "fitness_profile_id": "fit-neutral-v1",
                "initial_population_id": "ip001-fisher100",
                "config_path": str(config_path.relative_to(work)),
                "status": "prepared",
                "notes": (
                    "Exploratory first-pilot "
                    f"{cell.pilot_tier} cell; not confirmatory; "
                    "not yet launched."
                ),
            }
        )
        registry_parameter_rows.extend(_parameter_rows(cell))

    matrix_fields = [
        "cell_id",
        "label",
        "experimental_group",
        "comparison_sets",
        "pilot_tier",
        "H",
        "f",
        "e",
        "R",
        "alpha",
        "u",
        "host_counts_mode",
        "status",
    ]
    files[design_directory / "phase1-first-pilot-cells.tsv"] = _tsv_text(
        matrix_rows, matrix_fields
    )
    files[design_directory / "phase1-first-pilot-core-cells.tsv"] = _tsv_text(
        [row for row in matrix_rows if row["pilot_tier"] == "core"], matrix_fields
    )
    files[design_directory / "phase1-first-pilot-extension-cells.tsv"] = _tsv_text(
        [row for row in matrix_rows if row["pilot_tier"] == "extension"],
        matrix_fields,
    )
    files[work / "registry" / "cells.csv"] = _csv_text(
        registry_cell_rows,
        [
            "cell_id",
            "phase_id",
            "stage_id",
            "label",
            "mnemonic",
            "cell_dirname",
            "experimental_group",
            "comparison_set",
            "confirmatory",
            "architecture_profile_id",
            "selection_profile_id",
            "fitness_profile_id",
            "initial_population_id",
            "config_path",
            "status",
            "notes",
        ],
    )
    files[work / "registry" / "cell_parameters.csv"] = _csv_text(
        registry_parameter_rows,
        ["cell_id", "parameter_name", "value", "value_type", "unit", "role"],
    )

    files[manifest_directory / "phase1-first-pilot-seed-blocks.tsv"] = _tsv_text(
        seed_rows,
        ["seed_block_id", "master_seed", "within_run_replicate_index", "use"],
    )
    files[manifest_directory / "phase1-first-pilot-runs.tsv"] = _tsv_text(
        run_records,
        [
            "run_id",
            "cell_id",
            "seed_block_id",
            "master_seed",
            "within_run_replicate_index",
            "config_path",
            "config_sha256",
            "scratch_relative_path",
            "status",
        ],
    )

    manifest = {
        "experiment_manifest_schema_version": "1.0.0",
        "experiment_id": "phase1-first-pilot-core12",
        "design_id": "phase1-first-pilot-20cell",
        "execution_strategy": "run core 12, assess safety, then run extension 8",
        "status": "prepared-not-launched",
        "confirmatory": False,
        "model_family": "wright_fisher_counts",
        "model_spec_version": "2.0.0",
        "output_schema_version": "2.2.0",
        "software_version": "0.6.0",
        "initial_population_id": population["initial_population_id"],
        "initial_counts_sha256": count_checksum,
        "neutral_profiles": {
            "architecture_profile_id": "arch-panmictic-v1",
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
            "within_host_selection": False,
            "free_living_selection": False,
            "mutation_effect_mean": 0.0,
            "mutation_effect_sd": 0.0,
            "environment_counts_mode": "final",
            "replicates_per_run": 1,
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
            "note": "Commit the prepared code and files before starting sb0001.",
        },
    }
    files[manifest_directory / "phase1-first-pilot-manifest.json"] = (
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
            raise SystemExit("pilot files differ:\n" + "\n".join(mismatches))
        print(f"Verified {len(files)} Phase 1 first-pilot files")
        return 0
    for path, content in files.items():
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(content, encoding="utf-8")
    print(f"Wrote {len(files)} Phase 1 first-pilot files")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
