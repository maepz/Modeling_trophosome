#!/usr/bin/env python3
"""Freeze the Phase 1 Stage 3 Wave 2 equivalence experiments."""

from __future__ import annotations

import argparse
import csv
import json
from dataclasses import dataclass
from decimal import Decimal
from itertools import product
from pathlib import Path

from prepare_phase1_first_pilot import _array, _sha256, _tsv_text
from prepare_phase1_second_pilot import _sync_registry

MODEL_SPEC_VERSION = "2.1.0"
OUTPUT_SCHEMA_VERSION = "2.3.0"
SOFTWARE_VERSION = "0.7.0"
VARIANT_TAG = "v210-adaptive-g1000"
EXPERIMENT_ID = f"phase1-stage3-wave2-{VARIANT_TAG}"
STAGE_DIRECTORY = f"s03-parameter-map-wave2-{VARIANT_TAG}"
INITIAL_HORIZON = 100
INTERMEDIATE_HORIZON = 500
MAXIMUM_HORIZON = 1000
K = 1_000_000_000
U = "0"
SEED_BLOCKS = tuple((f"sb{number:04d}", 665 + number) for number in range(1, 13))
H_LEVELS = (100, 1000, 10000)
B_LEVELS = (1, 5, 10, 50)
ALPHA_LEVELS = ("0", "0.01", "0.1", "0.99")
MIGRATION_LEVELS = ("0", "0.001", "0.01", "0.1", "0.5", "0.9", "0.99")
ADAPTIVE_MIGRATION_LEVELS = ("0", "0.001", "0.01")
RETURN_BY_ALPHA_TARGET = {
    "0": 0,
    "0.01": 10_100_000,
    "0.1": 111_110_000,
    "0.5": 1_000_000_000,
    "0.99": 99_000_000_000,
}
EXPECTED_NEW_CELLS = 34
EXPECTED_REUSED_CELLS = 6
EXPECTED_RUNS = EXPECTED_NEW_CELLS * len(SEED_BLOCKS)
SMOKE_CELLS = ("c0051", "c0062", "c0084")


@dataclass(frozen=True)
class Wave2Cell:
    number: int
    panel: str
    hosts: int
    bottleneck: int
    alpha_target: str
    migration_fraction: str
    reused_source: str = ""

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
    def label(self) -> str:
        if self.panel == "H-by-B":
            return f"Host pooling: H={self.hosts}, B={self.bottleneck}"
        return (
            f"Feedback by regional exchange: alpha={self.alpha_target}, "
            f"m={self.migration_fraction}"
        )

    @property
    def initial_source_role(self) -> str:
        return "reused-passage-100" if self.reused_source else "new-simulation"


_reuse = {
    ("H-by-B", 100, 10, "0.5", "0.1"): "p01-s02-c0022",
    ("H-by-B", 10000, 10, "0.5", "0.1"): "p01-s02-c0024",
    ("alpha-by-m", 1000, 10, "0", "0.1"): "p01-s02-c0021",
    ("alpha-by-m", 1000, 10, "0.01", "0.1"): "p01-s03-c0037",
    ("alpha-by-m", 1000, 10, "0.1", "0.1"): "p01-s03-c0039",
    ("alpha-by-m", 1000, 10, "0.99", "0.1"): "p01-s03-c0041",
}

_cells: list[Wave2Cell] = []
number = 51
for hosts, bottleneck in product(H_LEVELS, B_LEVELS):
    key = ("H-by-B", hosts, bottleneck, "0.5", "0.1")
    _cells.append(Wave2Cell(number, *key[:-2], key[-2], key[-1], _reuse.get(key, "")))
    number += 1
for alpha, migration in product(ALPHA_LEVELS, MIGRATION_LEVELS):
    key = ("alpha-by-m", 1000, 10, alpha, migration)
    _cells.append(Wave2Cell(number, *key[:-2], key[-2], key[-1], _reuse.get(key, "")))
    number += 1
CELLS = tuple(_cells)
NEW_CELLS = tuple(cell for cell in CELLS if not cell.reused_source)
REUSED_CELLS = tuple(cell for cell in CELLS if cell.reused_source)


def _config_text(cell: Wave2Cell, counts: list[int], seed: int) -> str:
    fitness = _array([1.0] * len(counts))
    return f'''# {EXPERIMENT_ID}: {cell.cell_id}; {cell.label}
# Primary endpoint: passage 100. The trajectory may continue, without rerunning,
# only under the frozen adaptive time-horizon rule (maximum passage 1000).
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
host_generations = {MAXIMUM_HORIZON}
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
            "panel": cell.panel,
            "H": cell.hosts,
            "B": cell.bottleneck,
            "f": cell.escape_fraction,
            "e": cell.escape_cells,
            "R": cell.total_return,
            "alpha": format(cell.alpha, ".12g"),
            "alpha_target": cell.alpha_target,
            "m": cell.migration_fraction,
            "u": U,
            "initial_horizon": INITIAL_HORIZON,
            "maximum_horizon": MAXIMUM_HORIZON,
            "host_counts_mode": cell.host_counts_mode,
            "initial_source_role": cell.initial_source_role,
            "reused_source_cell_id": cell.reused_source,
            "adaptive_eligible": str(
                cell.panel == "alpha-by-m"
                and cell.migration_fraction in ADAPTIVE_MIGRATION_LEVELS
            ).lower(),
            "purpose": cell.label,
        }
        for cell in CELLS
    ]


def _reference_rows(
    work: Path,
) -> tuple[list[dict[str, str]], dict[str, object]]:
    phase = work / "p01-neutral-feedback"
    frozen = phase / "design" / f"{EXPERIMENT_ID}-reused-trajectories.tsv"
    manifest = phase / "manifests" / f"{EXPERIMENT_ID}-manifest.json"
    if frozen.is_file() and manifest.is_file():
        payload = json.loads(manifest.read_text(encoding="utf-8"))
        provenance = dict(payload["reused_trajectory_provenance"])
        expected = provenance.pop("frozen_sha256")
        provenance.pop("note", None)
        if _sha256(frozen.read_text(encoding="utf-8")) != expected:
            raise ValueError("frozen Wave 2 reused trajectory checksum differs")
        with frozen.open(newline="", encoding="utf-8") as handle:
            return list(csv.DictReader(handle, delimiter="\t")), provenance

    sources = {
        "stage2": phase
        / (
            "analysis/s02-equilibrium-precision-v210-m010-g250-derived/"
            "environment-trajectories.tsv"
        ),
        "wave1": phase
        / (
            "analysis/s03-parameter-map-wave1-v210-m010-g100-derived/"
            "environment-trajectories.tsv"
        ),
    }
    audits = {
        "stage2": phase
        / (
            "analysis/s02-equilibrium-precision-v210-m010-g250-derived/"
            "analysis-audit.json"
        ),
        "wave1": phase
        / "analysis/s03-parameter-map-wave1-v210-m010-g100-derived/analysis-audit.json",
    }
    for audit in audits.values():
        if json.loads(audit.read_text(encoding="utf-8"))["status"] != "PASS":
            raise ValueError(f"reference analysis did not pass: {audit}")

    by_source = {cell.reused_source: cell for cell in REUSED_CELLS}
    collected: list[dict[str, str]] = []
    for source_name, path in sources.items():
        with path.open(newline="", encoding="utf-8") as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                cell = by_source.get(row["cell_id"])
                if cell is None or row["seed_block_id"] not in dict(SEED_BLOCKS):
                    continue
                generation = int(row["generation"])
                if generation > INITIAL_HORIZON:
                    continue
                collected.append(
                    {
                        "cell_id": cell.cell_id,
                        "cell": cell.short_id,
                        "source_cell_id": row["cell_id"],
                        "source_run_id": row["run_id"],
                        "source_analysis": source_name,
                        "seed_block_id": row["seed_block_id"],
                        "generation": str(generation),
                        "D0": row["D0"],
                        "D1": row["D1"],
                        "D2": row["D2"],
                        "evenness": row["evenness"],
                        "TV": row["TV"],
                    }
                )
    expected = len(REUSED_CELLS) * len(SEED_BLOCKS) * (INITIAL_HORIZON + 1)
    keys = {
        (row["cell_id"], row["seed_block_id"], row["generation"]) for row in collected
    }
    if len(collected) != expected or len(keys) != expected:
        raise ValueError(
            f"expected {expected} unique reused trajectory rows, found "
            f"{len(collected)} ({len(keys)} unique)"
        )
    collected.sort(
        key=lambda row: (
            row["cell_id"],
            row["seed_block_id"],
            int(row["generation"]),
        )
    )
    provenance = {
        name: {
            "path": str(path.relative_to(work)),
            "sha256": _sha256(path.read_text(encoding="utf-8")),
            "audit_path": str(audits[name].relative_to(work)),
            "audit_sha256": _sha256(audits[name].read_text(encoding="utf-8")),
        }
        for name, path in sources.items()
    }
    return collected, provenance


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
        if cell.reused_source:
            continue
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
                    "panel": cell.panel,
                    "variant_id": VARIANT_TAG,
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
    references, source_provenance = _reference_rows(work)
    reference_text = _tsv_text(references, list(references[0]))
    files[phase / "design" / f"{EXPERIMENT_ID}-cells.tsv"] = _tsv_text(
        matrix, list(matrix[0])
    )
    files[phase / "design" / f"{EXPERIMENT_ID}-reused-trajectories.tsv"] = (
        reference_text
    )
    seed_rows = [
        {
            "seed_block_id": block,
            "master_seed": seed,
            "use": "exploratory matched population; shared with Stage 2 and Wave 1",
        }
        for block, seed in SEED_BLOCKS
    ]
    files[phase / "manifests" / f"{EXPERIMENT_ID}-seed-blocks.tsv"] = _tsv_text(
        seed_rows, list(seed_rows[0])
    )
    files[phase / "manifests" / f"{EXPERIMENT_ID}-runs.tsv"] = _tsv_text(
        runs, list(runs[0])
    )

    manifest = {
        "experiment_manifest_schema_version": "1.3.0",
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
        "primary_endpoint": INITIAL_HORIZON,
        "maximum_configured_horizon": MAXIMUM_HORIZON,
        "new_cells": len(NEW_CELLS),
        "reused_cells": len(REUSED_CELLS),
        "total_analysis_cells": len(CELLS),
        "new_populations": len(runs),
        "reused_populations": len(REUSED_CELLS) * len(SEED_BLOCKS),
        "panels": {
            "H_by_B": {
                "question": (
                    "Can many severe bottlenecks collectively provide a "
                    "representative return?"
                ),
                "H": list(H_LEVELS),
                "B": list(B_LEVELS),
                "alpha": 0.5,
                "m": 0.1,
                "matched_HB": [1000, 5000, 10000, 50000],
            },
            "alpha_by_m": {
                "question": (
                    "Does regional immigration erase or stabilize host-induced change?"
                ),
                "H": 1000,
                "B": 10,
                "alpha": [float(value) for value in ALPHA_LEVELS],
                "m": [float(value) for value in MIGRATION_LEVELS],
                "immediate_host_signal": "alpha * (1 - m)",
            },
        },
        "fixed_parameters": {
            "u": 0.0,
            "K": K,
            "N_E": K,
            "sampling_mode": "reservoir",
            "migration_mode": "fixed_regional_pool",
            "regional_pool": "same frozen composition as the focal start",
            "within_host_selection": False,
            "free_living_selection": False,
            "environment_counts_mode": "all",
            "replicates_per_run": 1,
            "workers_per_population": 2,
        },
        "adaptive_time_horizon": {
            "primary_inference_remains_at": INITIAL_HORIZON,
            "eligible_panel": "alpha-by-m only",
            "eligible_m": [float(value) for value in ADAPTIVE_MIGRATION_LEVELS],
            "assessment_horizons": [INITIAL_HORIZON, INTERMEDIATE_HORIZON],
            "possible_continuations": [INTERMEDIATE_HORIZON, MAXIMUM_HORIZON],
            "required_to_500": {
                "alpha": [0.0, 0.1],
                "m": [0.0, 0.001, 0.01],
                "reason": "pre-specified control and central-feedback anchors",
            },
            "required_to_1000": {
                "alpha": [0.0, 0.1],
                "m": [0.0, 0.001],
                "reason": (
                    "no-immigration diagnostic and m=0.001 memory half-life "
                    "of approximately 693 passages"
                ),
            },
            "additional_extension_rule": (
                "extend a positive-feedback eligible cell, together with its "
                "same-m alpha=0 control, when either its raw TV or its paired "
                "host-induced TV change fails the frozen stability rule"
            ),
            "windows_at_100": [[51, 75], [76, 100]],
            "windows_at_500": [[401, 450], [451, 500]],
            "stability_rule": {
                "uncertainty": "paired 90% Student-t interval across 12 seed blocks",
                "margin": "max(0.002 TV, 0.25 * absolute late-window mean)",
                "pass": "the complete interval lies inside plus/minus the margin",
            },
            "trajectory_integrity": (
                "all new TOMLs have host_generations=1000; execution pauses at "
                "100 or 500 and resumes the same verified checkpoint"
            ),
            "scope_warning": (
                "adaptive extensions diagnose time dependence; they do not "
                "change the pre-specified passage-100 primary comparisons"
            ),
        },
        "equivalence": {
            "biological_TV_margin": 0.05,
            "contour_ratio_interval": [0.8, 1.25],
            "uncertainty": "paired seed-block log-ratio intervals",
        },
        "analysis": {
            "exploratory": True,
            "seed_block_pairing": (
                "all contrasts and uncertainty preserve the 12 matched seeds"
            ),
            "primary_endpoint": {
                "generation": INITIAL_HORIZON,
                "response": "TV distance from the shared initial composition",
            },
            "H_by_B": {
                "models": [
                    "log(TV) ~ log(H*B) + seed-block intercept",
                    "log(TV) ~ log(H) * log(B) + seed-block intercept",
                ],
                "comparison": (
                    "compare HB-only fit with separate H, B and interaction; "
                    "report beta_H and beta_B"
                ),
                "exact_HB_contrasts": [1000, 5000, 10000, 50000],
            },
            "alpha_by_m": {
                "paired_host_effect": ("TV(alpha,m) - TV(alpha=0,m) within seed block"),
                "candidate_models": [
                    "immediate signal alpha*(1-m) only",
                    "separate alpha, m and alpha-by-m interaction",
                    "finite-time environmental-memory model",
                ],
                "erasure": (
                    "lower mean host-induced TV relative to the same-m control"
                ),
                "stabilization": (
                    "lower passages-51-100 temporal SD and/or lower among-seed "
                    "variation, reported separately from erasure"
                ),
            },
            "registered_responses": [
                "TV at passage 100",
                "mean TV during passages 51-100",
                "within-population temporal SD of TV during passages 51-100",
                "among-seed TV variation",
                "D1",
                "D2",
                "evenness",
            ],
            "future_grid_additions": (
                "B=2 or B=20 require a separate reviewed design amendment if "
                "the frozen H-by-B surface shows important curvature"
            ),
        },
        "cells": matrix,
        "seed_blocks": seed_rows,
        "runs": runs,
        "reused_trajectory_provenance": {
            **source_provenance,
            "frozen_sha256": _sha256(reference_text),
            "note": (
                "Complete environmental generations 0-100 for five exact prior "
                "conditions plus the environmentally equivalent no-return c0021 "
                "control; 12 matched seed blocks per condition"
            ),
        },
        "smoke_test": {
            "cells": list(SMOKE_CELLS),
            "seed_block": "sb0001",
            "included_in_new_populations": True,
        },
        "launch_policy": (
            "freeze and commit before launch; no simulation is run by this generator"
        ),
    }
    files[phase / "manifests" / f"{EXPERIMENT_ID}-manifest.json"] = (
        json.dumps(manifest, indent=2) + "\n"
    )
    return files


def registry_rows() -> tuple[list[dict[str, str]], list[dict[str, str]]]:
    cells: list[dict[str, str]] = []
    parameters: list[dict[str, str]] = []
    for cell in CELLS:
        config_path = (
            f"p01-neutral-feedback/configs/{STAGE_DIRECTORY}/"
            f"{cell.cell_id}-{VARIANT_TAG}.toml"
        )
        cells.append(
            {
                "cell_id": cell.cell_id,
                "phase_id": "p01",
                "stage_id": "s03",
                "label": cell.label,
                "mnemonic": (
                    f"h{cell.hosts}-b{cell.bottleneck}-a{cell.alpha_target}"
                    f"-m{cell.migration_fraction}-u0-g100-adaptive1000"
                ),
                "cell_dirname": cell.cell_id,
                "experimental_group": cell.panel,
                "comparison_set": "stage3-wave2",
                "confirmatory": "false",
                "architecture_profile_id": "arch-fixed-regional-pool-v1",
                "selection_profile_id": "sel-neutral-v1",
                "fitness_profile_id": "fit-neutral-v1",
                "initial_population_id": "ip001-fisher100",
                "config_path": config_path,
                "status": (
                    "prepared-reused-at-p100" if cell.reused_source else "prepared"
                ),
                "notes": (
                    "Exploratory Wave 2; primary endpoint passage 100; "
                    "frozen checkpoint-based adaptive horizon to passage 1000."
                ),
            }
        )
        values = (
            ("panel", cell.panel, "string", "", "design"),
            ("H", cell.hosts, "integer", "hosts", "input"),
            ("B", cell.bottleneck, "integer", "cells", "input"),
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
                "primary_horizon",
                INITIAL_HORIZON,
                "integer",
                "host_passages",
                "analysis",
            ),
            (
                "maximum_horizon",
                MAXIMUM_HORIZON,
                "integer",
                "host_passages",
                "technical",
            ),
            (
                "adaptive_eligible",
                str(
                    cell.panel == "alpha-by-m"
                    and cell.migration_fraction in ADAPTIVE_MIGRATION_LEVELS
                ).lower(),
                "boolean",
                "",
                "analysis",
            ),
            ("reused_source_cell_id", cell.reused_source, "string", "", "design"),
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
                    "Wave 2 execution has begun; frozen inputs cannot be rewritten"
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
            raise SystemExit("Wave 2 frozen inputs differ:\n" + "\n".join(issues))
        print(
            f"Verified {len(files)} Wave 2 files and both registries; "
            f"{len(CELLS)} cells ({len(NEW_CELLS)} new, "
            f"{len(REUSED_CELLS)} reused) and "
            f"{len(runs_from_files(files))} new populations."
        )
    else:
        for path, content in files.items():
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_text(content, encoding="utf-8")
        print(
            f"Wrote {len(files)} Wave 2 files; {len(CELLS)} cells, "
            f"{EXPECTED_RUNS} new populations; no simulations launched."
        )
    return 0


def runs_from_files(files: dict[Path, str]) -> list[str]:
    """Return generated per-population TOMLs (used only for concise reporting)."""

    return [
        str(path)
        for path in files
        if path.parent.name == "runs" and path.suffix == ".toml"
    ]


if __name__ == "__main__":
    raise SystemExit(main())
