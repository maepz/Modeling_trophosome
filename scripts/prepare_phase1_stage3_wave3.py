#!/usr/bin/env python3
"""Freeze the 16-cell Phase 1 Stage 3 bridge experiment."""

from __future__ import annotations

import argparse
import json
from dataclasses import dataclass
from decimal import Decimal
from pathlib import Path

from prepare_phase1_first_pilot import _array, _sha256, _tsv_text
from prepare_phase1_second_pilot import _sync_registry

MODEL_SPEC_VERSION = "2.1.0"
OUTPUT_SCHEMA_VERSION = "2.3.0"
SOFTWARE_VERSION = "0.7.0"
VARIANT_TAG = "v210-bridge-g100"
EXPERIMENT_ID = f"phase1-stage3-wave3-{VARIANT_TAG}"
STAGE_DIRECTORY = f"s03-parameter-map-wave3-{VARIANT_TAG}"
HOST_GENERATIONS = 100
K = 1_000_000_000
U = "0"
SEED_BLOCKS = tuple((f"sb{number:04d}", 665 + number) for number in range(1, 7))
RETURN_BY_ALPHA_TARGET = {
    "0.01": 10_100_000,
    "0.1": 111_110_000,
    "0.99": 99_000_000_000,
}
EXPECTED_CELLS = 16
EXPECTED_RUNS = EXPECTED_CELLS * len(SEED_BLOCKS)
SMOKE_CELLS = ("c0096", "c0098", "c0100")


@dataclass(frozen=True)
class BridgeCell:
    number: int
    bridge_group: str
    bridge_code: str
    hosts: int
    alpha_target: str
    bottleneck: int
    migration_fraction: str
    purpose: str
    extension: bool = False

    @property
    def short_id(self) -> str:
        return f"c{self.number:04d}"

    @property
    def cell_id(self) -> str:
        return f"p01-s03-{self.short_id}"

    @property
    def total_return(self) -> int:
        total = RETURN_BY_ALPHA_TARGET[self.alpha_target]
        if total % self.hosts:
            raise ValueError("frozen return is not divisible by host abundance")
        return total

    @property
    def escape_cells(self) -> int:
        return self.total_return // self.hosts

    @property
    def escape_fraction(self) -> str:
        return format(Decimal(self.escape_cells) / K, "g")

    @property
    def alpha(self) -> float:
        return self.total_return / (K + self.total_return)

    @property
    def host_counts_mode(self) -> str:
        return "full" if self.hosts <= 100 else "panel"

    @property
    def panel(self) -> str:
        names = {
            "A": "bridge-A-alpha-by-B",
            "B": "bridge-B-H-by-m",
            "C": "bridge-C-B-by-m",
            "D": "bridge-D-HB-by-alpha",
        }
        return names[self.bridge_group]

    @property
    def design_tier(self) -> str:
        return "two-cell-extension" if self.extension else "bridge-core"


CELLS = (
    BridgeCell(
        91, "A", "A1", 1000, "0.01", 1, "0.1", "Weak feedback, severe bottleneck"
    ),
    BridgeCell(
        92, "A", "A2", 1000, "0.01", 50, "0.1", "Weak feedback, wide bottleneck"
    ),
    BridgeCell(
        93, "A", "A3", 1000, "0.1", 1, "0.1", "Intermediate feedback, severe bottleneck"
    ),
    BridgeCell(
        94, "A", "A4", 1000, "0.1", 50, "0.1", "Intermediate feedback, wide bottleneck"
    ),
    BridgeCell(
        95, "A", "A5", 1000, "0.99", 1, "0.1", "Strong feedback, severe bottleneck"
    ),
    BridgeCell(
        96, "A", "A6", 1000, "0.99", 50, "0.1", "Strong feedback, wide bottleneck"
    ),
    BridgeCell(97, "B", "B1", 100, "0.1", 10, "0.01", "Few hosts, weak immigration"),
    BridgeCell(98, "B", "B2", 100, "0.1", 10, "0.9", "Few hosts, strong immigration"),
    BridgeCell(99, "B", "B3", 10000, "0.1", 10, "0.01", "Many hosts, weak immigration"),
    BridgeCell(
        100, "B", "B4", 10000, "0.1", 10, "0.9", "Many hosts, strong immigration"
    ),
    BridgeCell(
        101, "C", "C1", 1000, "0.1", 1, "0.01", "Severe bottleneck, weak immigration"
    ),
    BridgeCell(
        102, "C", "C2", 1000, "0.1", 1, "0.9", "Severe bottleneck, strong immigration"
    ),
    BridgeCell(
        103, "C", "C3", 1000, "0.1", 50, "0.01", "Wide bottleneck, weak immigration"
    ),
    BridgeCell(
        104, "C", "C4", 1000, "0.1", 50, "0.9", "Wide bottleneck, strong immigration"
    ),
    BridgeCell(
        105,
        "D",
        "D1",
        10000,
        "0.01",
        1,
        "0.1",
        "Matched-HB extension at weak feedback",
        True,
    ),
    BridgeCell(
        106,
        "D",
        "D2",
        10000,
        "0.99",
        1,
        "0.1",
        "Matched-HB extension at strong feedback",
        True,
    ),
)


def _config_text(cell: BridgeCell, counts: list[int], seed: int) -> str:
    fitness = _array([1.0] * len(counts))
    return f'''# {EXPERIMENT_ID}: {cell.cell_id}; bridge {cell.bridge_code}
# {cell.purpose}. Exploratory passage-100 bridge design.
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
fraction = {cell.migration_fraction}
regional_counts = {_array(counts)}

[host]
population_size = {cell.hosts}
infection_bottleneck = {cell.bottleneck}
carrying_capacity = {K}
growth_factor = 1.2
steady_generations = 500
host_generations = {HOST_GENERATIONS}
escape_fraction = {cell.escape_fraction}

[evolution]
mutation_probability = {U}
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
    return [
        {
            "cell_id": cell.cell_id,
            "cell": cell.short_id,
            "bridge_code": cell.bridge_code,
            "bridge_group": cell.bridge_group,
            "panel": cell.panel,
            "design_tier": cell.design_tier,
            "H": cell.hosts,
            "B": cell.bottleneck,
            "HB": cell.hosts * cell.bottleneck,
            "f": cell.escape_fraction,
            "e": cell.escape_cells,
            "R": cell.total_return,
            "alpha": format(cell.alpha, ".12g"),
            "alpha_target": cell.alpha_target,
            "m": cell.migration_fraction,
            "u": U,
            "host_generations": HOST_GENERATIONS,
            "host_counts_mode": cell.host_counts_mode,
            "purpose": cell.purpose,
        }
        for cell in CELLS
    ]


def build_files(repository: Path) -> dict[Path, str]:
    work = repository / "experiments/work/trophosome"
    phase = work / "p01-neutral-feedback"
    initial = json.loads(
        (work / "common/initial-populations/ip001-fisher100.json").read_text(
            encoding="utf-8"
        )
    )
    counts = [int(value) for value in initial["scaled_counts"]]
    if len(counts) != 100 or sum(counts) != K:
        raise ValueError("invalid frozen initial population")

    files: dict[Path, str] = {}
    runs: list[dict[str, object]] = []
    configs = phase / "configs" / STAGE_DIRECTORY
    for cell in CELLS:
        base = configs / f"{cell.cell_id}-{VARIANT_TAG}.toml"
        files[base] = _config_text(cell, counts, SEED_BLOCKS[0][1])
        for seed_block, seed in SEED_BLOCKS:
            run_id = f"{cell.cell_id}-{VARIANT_TAG}-{seed_block}"
            path = configs / "runs" / f"{run_id}.toml"
            content = _config_text(cell, counts, seed)
            files[path] = content
            runs.append(
                {
                    "experiment_id": EXPERIMENT_ID,
                    "run_id": run_id,
                    "cell_id": cell.cell_id,
                    "cell": cell.short_id,
                    "bridge_code": cell.bridge_code,
                    "bridge_group": cell.bridge_group,
                    "variant_id": VARIANT_TAG,
                    "pilot_tier": "stage3-wave3",
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
    seed_rows = [
        {
            "seed_block_id": block,
            "master_seed": seed,
            "use": "exploratory matched population; shared with Waves 1 and 2",
        }
        for block, seed in SEED_BLOCKS
    ]
    files[phase / "design" / f"{EXPERIMENT_ID}-cells.tsv"] = _tsv_text(
        matrix, list(matrix[0])
    )
    files[phase / "manifests" / f"{EXPERIMENT_ID}-seed-blocks.tsv"] = _tsv_text(
        seed_rows, list(seed_rows[0])
    )
    files[phase / "manifests" / f"{EXPERIMENT_ID}-runs.tsv"] = _tsv_text(
        runs, list(runs[0])
    )
    manifest = {
        "experiment_manifest_schema_version": "1.0.0",
        "experiment_id": EXPERIMENT_ID,
        "variant_id": VARIANT_TAG,
        "status": "prepared-not-launched",
        "confirmatory": False,
        "model_family": "wright_fisher_counts",
        "model_spec_version": MODEL_SPEC_VERSION,
        "software_version": SOFTWARE_VERSION,
        "output_schema_version": OUTPUT_SCHEMA_VERSION,
        "initial_population_id": "ip001-fisher100",
        "initial_counts_sha256": initial["scaled_counts_sha256"],
        "new_cells": len(CELLS),
        "core_bridge_cells": sum(not cell.extension for cell in CELLS),
        "two_cell_extension": sum(cell.extension for cell in CELLS),
        "seed_blocks": seed_rows,
        "new_populations": len(runs),
        "cells": matrix,
        "runs": runs,
        "primary_endpoint": HOST_GENERATIONS,
        "fixed_parameters": {
            "u": 0.0,
            "K": K,
            "N_E": K,
            "capacity_ratio": 1.0,
            "growth_factor": 1.2,
            "steady_within_host_generations": 500,
            "sampling_mode": "reservoir",
            "migration_mode": "fixed_regional_pool",
            "regional_pool": "same frozen composition as the focal start",
            "within_host_selection": False,
            "free_living_selection": False,
        },
        "target_interactions": ["alpha:B", "H:m", "B:m"],
        "target_model": (
            "H + alpha + B + m + H:alpha + H:B + alpha:m + alpha:B + H:m + B:m"
        ),
        "interpretation": (
            "Biologically constrained bridge augmentation; not a claim that an "
            "unconstrained algorithm found the globally D-optimal 16-cell subset."
        ),
        "smoke_test": {
            "cells": list(SMOKE_CELLS),
            "seed_block": "sb0001",
            "included_in_new_populations": True,
            "runtime_limit_hours": 48,
            "projected_storage_limit_gib": 350,
            "storage_budget_gib": 500,
        },
        "community_compilation": {
            "launcher_option": "--community-only",
            "through_wave": 3,
            "includes": "Waves 1, 2 and 3 at passage 100, plus PRC trajectories",
        },
        "retention": "No automatic deletion; retain raw outputs until audited.",
    }
    files[phase / "manifests" / f"{EXPERIMENT_ID}-manifest.json"] = (
        json.dumps(manifest, indent=2) + "\n"
    )
    return files


def registry_rows() -> tuple[list[dict[str, str]], list[dict[str, str]]]:
    cells: list[dict[str, str]] = []
    parameters: list[dict[str, str]] = []
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
                "label": f"Bridge {cell.bridge_code}: {cell.purpose}",
                "mnemonic": (
                    f"h{cell.hosts}-b{cell.bottleneck}-a{cell.alpha_target}"
                    f"-m{cell.migration_fraction}-u0-g100-bridge"
                ),
                "cell_dirname": cell.cell_id,
                "experimental_group": cell.panel,
                "comparison_set": "stage3-wave3-bridge",
                "confirmatory": "false",
                "architecture_profile_id": "arch-fixed-regional-pool-v1",
                "selection_profile_id": "sel-neutral-v1",
                "fitness_profile_id": "fit-neutral-v1",
                "initial_population_id": "ip001-fisher100",
                "config_path": path,
                "status": "prepared",
                "notes": (
                    f"Exploratory Wave 3 {cell.design_tier}; six matched seeds; "
                    "passage-100 endpoint."
                ),
            }
        )
        values = (
            ("bridge_code", cell.bridge_code, "string", "", "design"),
            ("bridge_group", cell.bridge_group, "string", "", "design"),
            ("design_tier", cell.design_tier, "string", "", "design"),
            ("H", cell.hosts, "integer", "hosts", "input"),
            ("B", cell.bottleneck, "integer", "cells", "input"),
            ("HB", cell.hosts * cell.bottleneck, "integer", "cells", "derived"),
            ("f", cell.escape_fraction, "float", "fraction", "input"),
            ("e", cell.escape_cells, "integer", "cells_per_host", "derived"),
            ("R", cell.total_return, "integer", "cells", "derived"),
            ("alpha", format(cell.alpha, ".12g"), "float", "fraction", "derived"),
            ("alpha_target", cell.alpha_target, "float", "fraction", "input"),
            ("m", cell.migration_fraction, "float", "fraction", "input"),
            ("u", U, "float", "per_genome_per_bacterial_generation", "nuisance"),
            ("K", K, "integer", "cells", "nuisance"),
            ("N_E", K, "integer", "cells", "nuisance"),
            (
                "host_generations",
                HOST_GENERATIONS,
                "integer",
                "host_passages",
                "analysis",
            ),
            ("host_counts_mode", cell.host_counts_mode, "string", "", "technical"),
            ("within_host_selection", "false", "boolean", "", "nuisance"),
            ("free_living_selection", "false", "boolean", "", "nuisance"),
            (
                "planned_seed_blocks",
                len(SEED_BLOCKS),
                "integer",
                "seed_blocks",
                "technical",
            ),
        )
        parameters.extend(
            {
                "cell_id": cell.cell_id,
                "parameter_name": name,
                "value": str(value),
                "value_type": value_type,
                "unit": unit,
                "role": role,
            }
            for name, value, value_type, unit, role in values
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
            scratch = Path(json.loads(layout.read_text(encoding="utf-8"))["scratch"])
            stage = scratch / "p01-neutral-feedback" / STAGE_DIRECTORY
            if stage.is_dir() and any(stage.rglob("provenance.json")):
                raise SystemExit(
                    "Wave 3 execution has begun; frozen inputs cannot be rewritten"
                )

    files = build_files(repository)
    issues: list[str] = []
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
            raise SystemExit("Wave 3 frozen inputs differ:\n" + "\n".join(issues))
        print(
            f"Verified {len(files)} Wave 3 files and both registries; "
            f"{len(CELLS)} cells and {EXPECTED_RUNS} populations."
        )
    else:
        for path, content in files.items():
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_text(content, encoding="utf-8")
        print(
            f"Wrote {len(files)} Wave 3 files; {len(CELLS)} cells and "
            f"{EXPECTED_RUNS} populations; no simulations launched."
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
