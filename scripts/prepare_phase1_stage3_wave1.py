#!/usr/bin/env python3
"""Freeze Stage 3 part one: 24 new cells plus one reused no-return control."""

from __future__ import annotations

import argparse
import csv
import json
import math
import statistics
from dataclasses import dataclass
from decimal import ROUND_HALF_EVEN, Decimal
from itertools import product
from pathlib import Path

from prepare_phase1_first_pilot import _array, _sha256, _tsv_text
from prepare_phase1_second_pilot import _sync_registry

from trophosome.config import HostConfig
from trophosome.count_model import population_size_schedule

MODEL_SPEC_VERSION = "2.1.0"
OUTPUT_SCHEMA_VERSION = "2.3.0"
SOFTWARE_VERSION = "0.7.0"
VARIANT_TAG = "v210-m010-g100"
EXPERIMENT_ID = f"phase1-stage3-wave1-{VARIANT_TAG}"
STAGE_DIRECTORY = f"s03-parameter-map-wave1-{VARIANT_TAG}"
HOST_GENERATIONS = 100
TAIL_START = 51
K = 1_000_000_000
B = 10
SEED_BLOCKS = tuple((f"sb{n:04d}", 665 + n) for n in range(1, 13))
HOST_LEVELS = (100, 1000, 10000)
FEEDBACK_LEVELS = ("0.001", "0.01", "0.1", "0.99")
MUTATION_LEVELS = ("0", "1e-10")
RETURN_QUANTUM = math.lcm(*HOST_LEVELS)
CONTROL_ID = "p01-s02-c0021"
SMOKE_CELLS = ("c0034", "c0049", "c0050")
EXPECTED_RUNS = 288
TOTAL_OFFSPRING = sum(
    population_size_schedule(
        HostConfig(
            population_size=100,
            infection_bottleneck=B,
            carrying_capacity=K,
            growth_factor=1.2,
            steady_generations=500,
            host_generations=HOST_GENERATIONS,
            escape_fraction=0.01,
        )
    )
)


@dataclass(frozen=True)
class MappingCell:
    number: int
    hosts: int
    alpha_target: str
    mutation_probability: str

    @property
    def purpose(self) -> str:
        return "Host abundance x feedback; " + (
            "mutation-free" if self.mutation_probability == "0" else "mutation-enabled"
        )

    @property
    def short_id(self) -> str:
        return f"c{self.number:04d}"

    @property
    def cell_id(self) -> str:
        return f"p01-s03-{self.short_id}"

    @property
    def escape_cells(self) -> int:
        return self.total_return // self.hosts

    @property
    def escape_fraction(self) -> str:
        return format(Decimal(self.escape_cells) / K, "g")

    @property
    def total_return(self) -> int:
        target = Decimal(self.alpha_target)
        quanta = (K * target / (1 - target) / RETURN_QUANTUM).to_integral_value(
            rounding=ROUND_HALF_EVEN
        )
        return int(quanta) * RETURN_QUANTUM

    @property
    def alpha(self) -> float:
        return self.total_return / (K + self.total_return)

    @property
    def host_counts_mode(self) -> str:
        return "full" if self.hosts <= 100 else "panel"

    @property
    def escape_range(self) -> str:
        fraction = Decimal(self.escape_fraction)
        return (
            "below-primary-range"
            if fraction < Decimal("1e-5")
            else "above-primary-range"
            if fraction > Decimal("1e-2")
            else "primary-range"
        )


CELLS = tuple(
    MappingCell(number, hosts, alpha, mutation)
    for number, (hosts, alpha, mutation) in enumerate(
        product(HOST_LEVELS, FEEDBACK_LEVELS, MUTATION_LEVELS), start=27
    )
)


def _config_text(cell: MappingCell, counts: list[int], seed: int) -> str:
    fitness = _array([1.0] * len(counts))
    return f'''# {EXPERIMENT_ID}: {cell.cell_id}; {cell.purpose}
# Exploratory finite-time outcome at 100 passages, not assumed equilibrium.
# Requested alpha={cell.alpha_target}; realized alpha={cell.alpha:.15g}.
# R={cell.total_return} exactly matched across H; {cell.escape_range}.
model = "wright_fisher_counts"
seed = {seed}
replicates = 1

[environment]
initial_counts = {_array(counts)}
initial_within_host_fitness = {fitness}
initial_free_living_fitness = {fitness}
sampling_mode = "reservoir"
capacity_ratio = 1.0

[migration]
mode = "fixed_regional_pool"
fraction = 0.1
regional_counts = {_array(counts)}

[host]
population_size = {cell.hosts}
infection_bottleneck = 10
carrying_capacity = {K}
growth_factor = 1.2
steady_generations = 500
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


def matrix_rows() -> list[dict[str, object]]:
    rows = [
        {
            "cell_id": cell.cell_id,
            "cell": cell.short_id,
            "stratum": "mutation-free"
            if cell.mutation_probability == "0"
            else "mutation-enabled",
            "H": cell.hosts,
            "f": cell.escape_fraction,
            "e": cell.escape_cells,
            "R": cell.total_return,
            "alpha": format(cell.alpha, ".12g"),
            "alpha_target": cell.alpha_target,
            "escape_range": cell.escape_range,
            "source_role": "new-grid-cell",
            "u": cell.mutation_probability,
            "U": float(cell.mutation_probability) * TOTAL_OFFSPRING,
            "m": 0.1,
            "host_generations": HOST_GENERATIONS,
            "host_counts_mode": cell.host_counts_mode,
            "purpose": cell.purpose,
        }
        for cell in CELLS
    ]
    rows.insert(
        0,
        {
            "cell_id": CONTROL_ID,
            "cell": "c0021",
            "stratum": "no-return-control",
            "H": 100,
            "f": "0",
            "e": 0,
            "R": 0,
            "alpha": "0",
            "alpha_target": "0",
            "escape_range": "no-return-control",
            "source_role": "reused-control",
            "u": "0",
            "U": 0.0,
            "m": 0.1,
            "host_generations": HOST_GENERATIONS,
            "host_counts_mode": "full",
            "purpose": "Shared no-return control; extract Stage 2 passage 100",
        },
    )
    return rows


def _reference_files(work: Path) -> tuple[list[dict[str, str]], dict[str, object]]:
    """Freeze 12 matched passage-100 endpoints and tail summaries per Stage 2 cell."""
    phase = work / "p01-neutral-feedback"
    frozen = phase / "design" / f"{EXPERIMENT_ID}-reference-endpoints.tsv"
    manifest = phase / "manifests" / f"{EXPERIMENT_ID}-manifest.json"
    if frozen.is_file() and manifest.is_file():
        provenance = json.loads(manifest.read_text())["references"].copy()
        checksum = provenance.pop("frozen_endpoints_sha256")
        if _sha256(frozen.read_text()) != checksum:
            raise ValueError("frozen Stage 2 reference endpoint checksum differs")
        with frozen.open(newline="") as handle:
            return list(csv.DictReader(handle, delimiter="\t")), provenance
    derived = phase / "analysis/s02-equilibrium-precision-v210-m010-g250-derived"
    audit = derived / "analysis-audit.json"
    if json.loads(audit.read_text())["status"] != "PASS":
        raise ValueError("Stage 2 reference analysis has not passed its audit")
    source = derived / "environment-trajectories.tsv"
    completion_path = derived / "report-completion.json"
    completion = json.loads(completion_path.read_text())
    source_record = next(
        r
        for r in completion["outputs"]
        if r["path"].endswith("/environment-trajectories.tsv")
    )
    if not completion.get("complete") or source_record["sha256"] != _sha256(
        source.read_text()
    ):
        raise ValueError("Stage 2 trajectory checksum differs from completed report")
    with source.open(newline="") as handle:
        trajectory = [
            r
            for r in csv.DictReader(handle, delimiter="\t")
            if TAIL_START <= int(r["generation"]) <= HOST_GENERATIONS
            and r["seed_block_id"] in dict(SEED_BLOCKS)
        ]
    rows = [r for r in trajectory if int(r["generation"]) == HOST_GENERATIONS]
    expected = {
        (f"p01-s02-c{n:04d}", seed) for n in range(21, 27) for seed, _ in SEED_BLOCKS
    }
    if (
        len(rows) != 72
        or {(r["cell_id"], r["seed_block_id"]) for r in rows} != expected
    ):
        raise ValueError(
            "Stage 2 references do not contain six cells x twelve endpoints"
        )
    fields = (
        "run_id",
        "cell_id",
        "cell",
        "seed_block_id",
        "generation",
        "D0",
        "shannon",
        "simpson",
        "D1",
        "D2",
        "evenness",
        "TV",
        "root_D0",
        "root_D1",
        "root_D2",
        "mutant_richness",
        "mutant_abundance_fraction",
    )
    with (
        phase / "design/phase1-second-pilot-v210-m010-g250-cells.tsv"
    ).open() as handle:
        reference_cells = {
            r["cell_id"]: r for r in csv.DictReader(handle, delimiter="\t")
        }
    # Simpson (1 - sum p_i^2) is exactly 1 - 1/D2 for the stored Hill definition.
    for row in rows:
        row["simpson"] = str(1 - 1 / float(row["D2"]))
    rows = [{field: row[field] for field in fields} for row in rows]
    for row in rows:
        cell = reference_cells[row["cell_id"]]
        row.update({key: cell[key] for key in ("H", "f", "alpha", "u", "m")})
        row["source_role"] = (
            "reused-control"
            if row["cell_id"] == CONTROL_ID
            else "supplementary-reference"
        )
        tail = [
            r
            for r in trajectory
            if (r["cell_id"], r["seed_block_id"])
            == (row["cell_id"], row["seed_block_id"])
        ]
        if len(tail) != 50 or {int(r["generation"]) for r in tail} != set(
            range(51, 101)
        ):
            raise ValueError("Stage 2 reference tail is incomplete or duplicated")
        for response in ("D0", "D1", "D2", "evenness", "TV"):
            row[f"{response}_tail_mean"] = str(
                statistics.fmean(float(r[response]) for r in tail)
            )
            row[f"{response}_tail_sd"] = str(
                statistics.stdev(float(r[response]) for r in tail)
            )
    rows.sort(key=lambda r: (r["cell_id"], r["seed_block_id"]))
    return rows, {
        "source": str(source.relative_to(work)),
        "sha256": _sha256(source.read_text()),
        "audit": str(audit.relative_to(work)),
        "audit_sha256": _sha256(audit.read_text()),
        "completion_sha256": _sha256(completion_path.read_text()),
        "note": (
            "Frozen Stage 2 passage-100 endpoints and passages 51-100 summaries; "
            "c0021 is the shared control, c0022-c0026 are off-grid "
            "supplementary references. These are derived summaries, "
            "not copied passage-100 labelled count archives."
        ),
    }


def build_files(repository: Path) -> dict[Path, str]:
    work = repository / "experiments/work/trophosome"
    phase = work / "p01-neutral-feedback"
    initial = json.loads(
        (work / "common/initial-populations/ip001-fisher100.json").read_text()
    )
    counts = initial["scaled_counts"]
    if len(counts) != 100 or sum(counts) != K:
        raise ValueError("invalid frozen initial population")
    files: dict[Path, str] = {}
    runs: list[dict[str, object]] = []
    for cell in CELLS:
        directory = phase / "configs" / STAGE_DIRECTORY
        files[directory / f"{cell.cell_id}-{VARIANT_TAG}.toml"] = _config_text(
            cell, counts, SEED_BLOCKS[0][1]
        )
        for seed_block, seed in SEED_BLOCKS:
            run_id = f"{cell.cell_id}-{VARIANT_TAG}-{seed_block}"
            path = directory / "runs" / f"{run_id}.toml"
            content = _config_text(cell, counts, seed)
            files[path] = content
            runs.append(
                {
                    "experiment_id": EXPERIMENT_ID,
                    "run_id": run_id,
                    "cell_id": cell.cell_id,
                    "cell": cell.short_id,
                    "variant_id": VARIANT_TAG,
                    "pilot_tier": "stage3-wave1",
                    "seed_block_id": seed_block,
                    "master_seed": seed,
                    "within_run_replicate_index": 0,
                    "config_path": str(path.relative_to(work)),
                    "config_sha256": _sha256(content),
                    "scratch_relative_path": str(
                        Path("p01-neutral-feedback")
                        / STAGE_DIRECTORY
                        / cell.cell_id
                        / seed_block
                    ),
                    "status": "prepared",
                }
            )
    matrix = matrix_rows()
    references, provenance = _reference_files(work)
    reference_text = _tsv_text(references, list(references[0]))
    files[phase / "design" / f"{EXPERIMENT_ID}-reference-endpoints.tsv"] = (
        reference_text
    )
    files[phase / "design" / f"{EXPERIMENT_ID}-cells.tsv"] = _tsv_text(
        matrix, list(matrix[0])
    )
    seeds = [
        {
            "seed_block_id": block,
            "master_seed": seed,
            "use": "exploratory; matched to Stage 2; not held-out confirmation",
        }
        for block, seed in SEED_BLOCKS
    ]
    for suffix, rows in (("runs", runs), ("seed-blocks", seeds)):
        files[phase / "manifests" / f"{EXPERIMENT_ID}-{suffix}.tsv"] = _tsv_text(
            rows, list(rows[0])
        )
    manifest = {
        "experiment_manifest_schema_version": "1.0.0",
        "experiment_id": EXPERIMENT_ID,
        "model_spec_version": MODEL_SPEC_VERSION,
        "software_version": SOFTWARE_VERSION,
        "output_schema_version": OUTPUT_SCHEMA_VERSION,
        "status": "prepared-not-launched",
        "confirmatory": False,
        "cells": matrix,
        "seed_blocks": seeds,
        "runs": runs,
        "initial_population_id": "ip001-fisher100",
        "initial_counts_sha256": initial["scaled_counts_sha256"],
        "total_offspring_per_host_lifetime": TOTAL_OFFSPRING,
        "estimand": "TV at passage 100; not assumed equilibrium",
        "new_cells": len(CELLS),
        "new_populations": EXPECTED_RUNS,
        "primary_cells_including_reused_control": len(matrix),
        "primary_populations_including_reused_control": 300,
        "shared_control_cell_id": CONTROL_ID,
        "return_rounding": (
            "Round target N_E*alpha/(1-alpha) to nearest multiple of 10000, "
            "ties to even; then divide R exactly by H. All H have identical R."
        ),
        "references": {
            **provenance,
            "frozen_endpoints_sha256": _sha256(reference_text),
        },
        "smoke_test": {
            "cells": list(SMOKE_CELLS),
            "seed_block": "sb0001",
            "included_in_new_populations": True,
            "runtime_limit_hours": 48,
            "projected_storage_limit_gib": 350,
            "storage_budget_gib": 500,
        },
        "analysis": {
            "primary_response": "TV",
            "primary_endpoint": HOST_GENERATIONS,
            "descriptive_tail": [TAIL_START, HOST_GENERATIONS],
            "late_drift_check": "Paired TV mean: passages 76-100 minus 51-75",
            "uncertainty": "90% Student-t intervals across twelve populations",
            "relative_diversity_margin": 0.05,
            "evenness_margin": 0.02,
            "TV_margin": 0.05,
            "interactions": "Paired TV H x alpha differences-in-differences, by u",
            "note": (
                "Exploratory contrasts; no multiplicity adjustment "
                "or automatic adaptive launches."
            ),
        },
        "report_only_option": "--report-only",
        "retention": (
            "No automatic deletion; retain raw outputs until audited and archived."
        ),
    }
    files[phase / "manifests" / f"{EXPERIMENT_ID}-manifest.json"] = (
        json.dumps(manifest, indent=2) + "\n"
    )
    return files


def registry_rows() -> tuple[list[dict[str, str]], list[dict[str, str]]]:
    cells, parameters = [], []
    for cell in CELLS:
        path = (
            f"p01-neutral-feedback/configs/{STAGE_DIRECTORY}/"
            f"{cell.cell_id}-{VARIANT_TAG}.toml"
        )
        cells.append(
            {
                "cell_id": cell.cell_id,
                "phase_id": "p01",
                "stage_id": "s03",
                "label": cell.purpose,
                "mnemonic": (
                    f"h{cell.hosts}-f{cell.escape_fraction}"
                    f"-u{cell.mutation_probability}-g100-m010"
                ),
                "cell_dirname": cell.cell_id,
                "experimental_group": "MAP-U0"
                if cell.mutation_probability == "0"
                else "MAP-U+",
                "comparison_set": "stage3-wave1",
                "confirmatory": "false",
                "architecture_profile_id": "arch-fixed-regional-pool-v1",
                "selection_profile_id": "sel-neutral-v1",
                "fitness_profile_id": "fit-neutral-v1",
                "initial_population_id": "ip001-fisher100",
                "config_path": path,
                "status": "prepared",
                "notes": (
                    "Exploratory finite-time map; twelve matched seeds; not launched. "
                    + cell.escape_range
                ),
            }
        )
        for name, value, kind, unit, role in (
            ("H", cell.hosts, "integer", "hosts", "input"),
            ("f", cell.escape_fraction, "float", "fraction", "input"),
            (
                "u",
                cell.mutation_probability,
                "float",
                "per_genome_per_bacterial_generation",
                "input",
            ),
            ("R", cell.total_return, "integer", "cells", "derived"),
            ("alpha", format(cell.alpha, ".12g"), "float", "fraction", "derived"),
            ("alpha_target", cell.alpha_target, "float", "fraction", "input"),
            ("escape_range", cell.escape_range, "string", "", "derived"),
            (
                "U",
                float(cell.mutation_probability) * TOTAL_OFFSPRING,
                "float",
                "events_per_host_lifetime",
                "derived",
            ),
            ("B", B, "integer", "cells", "nuisance"),
            ("K", K, "integer", "cells", "nuisance"),
            ("N_E", K, "integer", "cells", "nuisance"),
            ("m", 0.1, "float", "fraction", "nuisance"),
            (
                "host_generations",
                HOST_GENERATIONS,
                "integer",
                "host_passages",
                "technical",
            ),
            (
                "planned_seed_blocks",
                len(SEED_BLOCKS),
                "integer",
                "seed_blocks",
                "technical",
            ),
            ("host_counts_mode", cell.host_counts_mode, "string", "", "technical"),
            ("within_host_selection", "false", "boolean", "", "nuisance"),
            ("free_living_selection", "false", "boolean", "", "nuisance"),
        ):
            parameters.append(
                {
                    "cell_id": cell.cell_id,
                    "parameter_name": name,
                    "value": str(value),
                    "value_type": kind,
                    "unit": unit,
                    "role": role,
                }
            )
    return cells, parameters


def verify_files(repository: Path) -> list[str]:
    return [
        str(path.relative_to(repository))
        for path, content in build_files(repository).items()
        if not path.is_file() or path.read_text(encoding="utf-8") != content
    ]


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--repository", type=Path, default=Path(__file__).resolve().parents[1]
    )
    action = parser.add_mutually_exclusive_group(required=True)
    action.add_argument("--write", action="store_true")
    action.add_argument("--verify", action="store_true")
    args = parser.parse_args()
    repository = args.repository.resolve()
    if args.write:
        layout = repository / "experiments/work/trophosome/layout.local.json"
        if layout.is_file():
            stage = (
                Path(json.loads(layout.read_text())["scratch"])
                / "p01-neutral-feedback"
                / STAGE_DIRECTORY
            )
            if stage.is_dir() and any(stage.rglob("provenance.json")):
                raise SystemExit(
                    "Stage 3 execution has begun; frozen inputs cannot be rewritten"
                )
    files = build_files(repository)
    issues = []
    for name, rows in zip(
        ("cells.csv", "cell_parameters.csv"), registry_rows(), strict=True
    ):
        issues.extend(
            _sync_registry(
                repository / "experiments/work/trophosome/registry" / name,
                rows,
                write=args.write,
            )
        )
    if args.verify:
        issues.extend(verify_files(repository))
        if issues:
            raise SystemExit("Stage 3 frozen inputs differ:\n" + "\n".join(issues))
        print(
            f"Verified {len(files)} Stage 3 files and both registries; "
            f"{EXPECTED_RUNS} new populations."
        )
    else:
        for path, content in files.items():
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_text(content, encoding="utf-8")
        print(f"Wrote {len(files)} Stage 3 files; no simulations launched.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
