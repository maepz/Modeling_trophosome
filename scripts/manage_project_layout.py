#!/usr/bin/env python3
"""Create and query the manifest-driven trophosome project layout.

The script deliberately uses only the Python standard library so that it can be
run on a laptop or an HPC login node before the model environment is installed.

Examples
--------
Initialize the default layout below a home/project root::

    python scripts/manage_project_layout.py init --root "$HOME"

Register a parameter cell and create its stable-ID scratch directory::

    python scripts/manage_project_layout.py register-cell \
        --root "$HOME" \
        --cell-id p01-s03-c0042 \
        --label "Baseline mutation mapping" \
        --group M \
        --architecture arch-panmictic-v1 \
        --param H=10000 \
        --param B=10 \
        --param K=1e9 \
        --param steady_generations=500 \
        --param f=1e-4 \
        --param u=1e-10 \
        --param c=1

Review it later using only the stable cell ID::

    python scripts/manage_project_layout.py show-cell \
        --root "$HOME" p01-s03-c0042
"""

from __future__ import annotations

import argparse
import csv
import json
import re
import sys
from collections.abc import Iterable
from dataclasses import dataclass
from decimal import Decimal, InvalidOperation
from pathlib import Path
from typing import Any

LAYOUT_SCHEMA_VERSION = "1.2.0"
DEFAULT_PROJECT = "trophosome"

DEFAULT_PHASES: tuple[tuple[str, str], ...] = (
    ("p01", "neutral-feedback"),
    ("p02", "selection"),
    ("p03", "architecture"),
)

DEFAULT_INITIAL_PHASE_IDS: tuple[str, ...] = ("p01",)

DEFAULT_STAGES: tuple[tuple[str, str], ...] = (
    ("s01", "pilot"),
    ("s02", "equilibrium-precision"),
    ("s03", "parameter-map"),
    ("s04", "confirmatory"),
    ("s05", "sensitivities"),
)

CELL_FIELDS: tuple[str, ...] = (
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
)

PARAMETER_FIELDS: tuple[str, ...] = (
    "cell_id",
    "parameter_name",
    "value",
    "value_type",
    "unit",
    "role",
)

PREFERRED_PARAMETER_ORDER: tuple[str, ...] = (
    "H",
    "host_abundance_H",
    "population_size",
    "B",
    "infection_bottleneck_B",
    "infection_bottleneck",
    "K",
    "carrying_capacity_K",
    "carrying_capacity",
    "growth_factor",
    "steady_generations",
    "host_generations",
    "f",
    "escape_fraction",
    "e",
    "escape_cells_per_host_e",
    "R",
    "total_return_R",
    "u",
    "mutation_probability_u",
    "mutation_probability",
    "c",
    "capacity_ratio_c",
    "capacity_ratio",
    "N_E",
    "environmental_capacity_NE",
    "alpha",
    "feedback_alpha",
    "architecture_mode",
    "lobe_count",
)

MNEMONIC_PARAMETERS: tuple[tuple[tuple[str, ...], str], ...] = (
    (("H", "host_abundance_H", "population_size"), "h"),
    (("f", "escape_fraction"), "f"),
    (("u", "mutation_probability_u", "mutation_probability"), "u"),
    (("c", "capacity_ratio_c", "capacity_ratio"), "c"),
    (("B", "infection_bottleneck_B", "infection_bottleneck"), "b"),
    (("K", "carrying_capacity_K", "carrying_capacity"), "k"),
    (("steady_generations",), "ts"),
)

CELL_ID_RE = re.compile(r"(?P<phase>p\d{2})-(?P<stage>s\d{2})-c(?P<number>\d{4,})")
SHORT_CELL_ID_RE = re.compile(r"c\d{4,}")


class LayoutError(ValueError):
    """Raised when a layout or registry operation is invalid."""


@dataclass(frozen=True)
class LayoutPaths:
    base_root: Path
    project: str
    work: Path
    scratch: Path
    data: Path

    @property
    def registry(self) -> Path:
        return self.work / "registry"

    @property
    def cells_registry(self) -> Path:
        return self.registry / "cells.csv"

    @property
    def parameters_registry(self) -> Path:
        return self.registry / "cell_parameters.csv"

    @property
    def layout_file(self) -> Path:
        return self.work / "layout.local.json"


def _expanded(path: str | Path) -> Path:
    return Path(path).expanduser().resolve()


def resolve_layout(args: argparse.Namespace) -> LayoutPaths:
    base_root = _expanded(args.root)
    project = args.project
    work_parent = _expanded(args.work_root) if args.work_root else base_root / "work"
    scratch_parent = (
        _expanded(args.scratch_root) if args.scratch_root else base_root / "scratch"
    )
    data_parent = _expanded(args.data_root) if args.data_root else base_root / "data"
    return LayoutPaths(
        base_root=base_root,
        project=project,
        work=work_parent / project,
        scratch=scratch_parent / project,
        data=data_parent / project,
    )


def phase_dirname(phase_id: str) -> str:
    for candidate, label in DEFAULT_PHASES:
        if candidate == phase_id:
            return f"{candidate}-{label}"
    raise LayoutError(
        f"unknown phase {phase_id!r}; add it to DEFAULT_PHASES before registering cells"
    )


def stage_dirname(stage_id: str) -> str:
    for candidate, label in DEFAULT_STAGES:
        if candidate == stage_id:
            return f"{candidate}-{label}"
    raise LayoutError(
        f"unknown stage {stage_id!r}; add it to DEFAULT_STAGES before registering cells"
    )


def parse_cell_id(value: str) -> tuple[str, str, str]:
    """Return canonical cell ID, phase ID, and stage ID.

    The input may be a bare cell ID, a cell directory path, or a complete run
    ID. This makes copy/paste from scheduler output convenient.
    """

    match = CELL_ID_RE.search(value)
    if match is None:
        raise LayoutError(
            f"could not find a cell ID in {value!r}; expected pNN-sNN-cNNNN"
        )
    canonical = match.group(0)
    return canonical, match.group("phase"), match.group("stage")


def ensure_csv(path: Path, fields: Iterable[str]) -> None:
    if path.exists():
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("x", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=list(fields),
            lineterminator="\n",
        )
        writer.writeheader()


def write_text_if_missing(path: Path, content: str) -> None:
    if path.exists():
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("x", encoding="utf-8") as handle:
        handle.write(content)


def initialize_layout(
    layout: LayoutPaths,
    phase_ids: Iterable[str] = DEFAULT_INITIAL_PHASE_IDS,
) -> None:
    requested_phase_ids = tuple(dict.fromkeys(phase_ids))
    for phase_id in requested_phase_ids:
        phase_dirname(phase_id)

    common_directories = (
        layout.work / "common" / "architecture-profiles",
        layout.work / "common" / "fitness-profiles",
        layout.work / "common" / "initial-populations",
        layout.work / "common" / "environments",
        layout.work / "common" / "shared-analysis",
        layout.work / "registry",
        layout.work / "environments",
        layout.work / "inventories",
        layout.scratch / "staging",
        layout.data / "releases",
    )
    for directory in common_directories:
        directory.mkdir(parents=True, exist_ok=True)

    for phase_id in requested_phase_ids:
        phase_name = phase_dirname(phase_id)
        work_phase = layout.work / phase_name
        scratch_phase = layout.scratch / phase_name
        for directory in (
            work_phase / "design",
            work_phase / "manifests" / "cells",
            work_phase / "manifests" / "runs",
            work_phase / "wrappers" / "local",
            work_phase / "wrappers" / "slurm",
            work_phase / "analysis",
            work_phase / "inventories",
        ):
            directory.mkdir(parents=True, exist_ok=True)
        for stage_id, stage_label in DEFAULT_STAGES:
            stage_name = f"{stage_id}-{stage_label}"
            (work_phase / "configs" / stage_name).mkdir(parents=True, exist_ok=True)
            (work_phase / "logs" / stage_name).mkdir(parents=True, exist_ok=True)
            (scratch_phase / stage_name).mkdir(parents=True, exist_ok=True)

    ensure_csv(layout.cells_registry, CELL_FIELDS)
    ensure_csv(layout.parameters_registry, PARAMETER_FIELDS)
    write_text_if_missing(
        layout.registry / "parameter-dictionary.md",
        """# Cell parameter dictionary

`cells.csv` stores stable cell identity and experimental metadata.

`cell_parameters.csv` stores an open-ended, long-form parameter registry. New
parameters may be added without changing its columns. The authoritative model
input remains the versioned TOML configuration and each run's
`resolved_config.json`.

One cell is one scientific parameter combination. One `(cell_id, replicate_id)`
pair is one stochastic run. Reusing replicate IDs across cells creates matched
seed blocks. Never reuse a cell ID after changing a scientific parameter.
""",
    )
    initialized_phases = [
        (phase_id, label)
        for phase_id, label in DEFAULT_PHASES
        if (layout.work / f"{phase_id}-{label}").is_dir()
        or (layout.scratch / f"{phase_id}-{label}").is_dir()
    ]
    payload = {
        "layout_schema_version": LAYOUT_SCHEMA_VERSION,
        "project": layout.project,
        "base_root": str(layout.base_root),
        "work": str(layout.work),
        "scratch": str(layout.scratch),
        "data": str(layout.data),
        "initialized_phases": [
            {"phase_id": phase_id, "directory": f"{phase_id}-{label}"}
            for phase_id, label in initialized_phases
        ],
        "known_phases": [
            {"phase_id": phase_id, "directory": f"{phase_id}-{label}"}
            for phase_id, label in DEFAULT_PHASES
        ],
        "stages": [
            {"stage_id": stage_id, "directory": f"{stage_id}-{label}"}
            for stage_id, label in DEFAULT_STAGES
        ],
    }
    with layout.layout_file.open("w", encoding="utf-8") as handle:
        json.dump(payload, handle, indent=2)
        handle.write("\n")


def read_csv(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        raise LayoutError(f"registry not found: {path}; run the init command first")
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def append_csv(path: Path, fields: Iterable[str], row: dict[str, Any]) -> None:
    with path.open("a", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=list(fields),
            extrasaction="ignore",
            lineterminator="\n",
        )
        writer.writerow({field: row.get(field, "") for field in fields})


def parse_assignments(values: list[str], option: str) -> dict[str, str]:
    parsed: dict[str, str] = {}
    for value in values:
        if "=" not in value:
            raise LayoutError(f"{option} expects NAME=VALUE, received {value!r}")
        name, assigned = value.split("=", 1)
        name = name.strip()
        assigned = assigned.strip()
        if not name or not assigned:
            raise LayoutError(f"{option} expects non-empty NAME=VALUE")
        if name in parsed:
            raise LayoutError(f"duplicate {option} name: {name}")
        parsed[name] = assigned
    return parsed


def infer_value_type(value: str) -> str:
    lowered = value.lower()
    if lowered in {"true", "false"}:
        return "boolean"
    try:
        integer = int(value)
    except ValueError:
        integer = None
    if integer is not None and str(integer) == value:
        return "integer"
    try:
        Decimal(value)
    except InvalidOperation:
        return "string"
    return "float"


def slug(value: str) -> str:
    value = value.strip().lower()
    value = re.sub(r"[^a-z0-9]+", "-", value)
    return value.strip("-") or "na"


def compact_number(value: str) -> str:
    """Return a filename-safe compact representation of a numeric value."""

    try:
        number = Decimal(value)
    except InvalidOperation:
        return slug(value)
    if not number.is_finite():
        return slug(value)
    if number == number.to_integral_value():
        return str(number.quantize(Decimal(1)))
    normalized = format(number.normalize(), "E").lower()
    coefficient, exponent = normalized.split("e", 1)
    coefficient = coefficient.rstrip("0").rstrip(".")
    exponent_number = int(exponent)
    exponent_text = (
        f"m{abs(exponent_number)}"
        if exponent_number < 0
        else f"p{exponent_number}"
    )
    coefficient = coefficient.replace("-", "m").replace(".", "p")
    return f"{coefficient}e{exponent_text}"


def _first_parameter(parameters: dict[str, str], names: Iterable[str]) -> str | None:
    for name in names:
        if name in parameters:
            return parameters[name]
    return None


def profile_mnemonic(value: str, prefix: str) -> str:
    cleaned = value
    for removable in ("-v1", "-v2", "profile", "arch-", "sel-", "fit-"):
        cleaned = cleaned.replace(removable, "")
    return prefix + slug(cleaned).replace("-", "")[:14]


def make_mnemonic(
    parameters: dict[str, str],
    architecture_profile: str,
    selection_profile: str,
    fitness_profile: str,
) -> str:
    parts: list[str] = []
    for names, prefix in MNEMONIC_PARAMETERS:
        value = _first_parameter(parameters, names)
        if value is not None:
            parts.append(prefix + compact_number(value))
    if architecture_profile:
        parts.append(profile_mnemonic(architecture_profile, "arch"))
    if selection_profile:
        parts.append(profile_mnemonic(selection_profile, "sel"))
    if fitness_profile:
        parts.append(profile_mnemonic(fitness_profile, "fit"))
    return "-".join(parts) if parts else "cell"


def config_relative_path(phase_id: str, stage_id: str, cell_id: str) -> str:
    return str(
        Path(phase_dirname(phase_id))
        / "configs"
        / stage_dirname(stage_id)
        / f"{cell_id}.toml"
    )


def register_cell(args: argparse.Namespace, layout: LayoutPaths) -> dict[str, str]:
    cell_id, phase_id, stage_id = parse_cell_id(args.cell_id)
    if cell_id != args.cell_id:
        raise LayoutError("register-cell requires a bare canonical cell ID")
    initialize_layout(layout, phase_ids=(phase_id,))
    existing = read_csv(layout.cells_registry)
    if any(row["cell_id"] == cell_id for row in existing):
        raise LayoutError(f"cell ID is already registered: {cell_id}")

    parameters = parse_assignments(args.param, "--param")
    units = parse_assignments(args.unit, "--unit")
    roles = parse_assignments(args.role, "--role")
    types = parse_assignments(args.value_type, "--value-type")
    unknown_metadata = (set(units) | set(roles) | set(types)) - set(parameters)
    if unknown_metadata:
        raise LayoutError(
            "parameter metadata supplied for missing parameters: "
            + ", ".join(sorted(unknown_metadata))
        )

    mnemonic = args.mnemonic or make_mnemonic(
        parameters,
        args.architecture,
        args.selection,
        args.fitness,
    )
    mnemonic = slug(mnemonic)
    cell_dirname = cell_id
    relative_config = config_relative_path(phase_id, stage_id, cell_id)
    scratch_cell = (
        layout.scratch
        / phase_dirname(phase_id)
        / stage_dirname(stage_id)
        / cell_dirname
    )
    if scratch_cell.exists():
        raise LayoutError(
            "scratch directory already exists for the proposed cell: "
            f"{scratch_cell}"
        )
    row = {
        "cell_id": cell_id,
        "phase_id": phase_id,
        "stage_id": stage_id,
        "label": args.label,
        "mnemonic": mnemonic,
        "cell_dirname": cell_dirname,
        "experimental_group": args.group,
        "comparison_set": args.comparison_set,
        "confirmatory": str(bool(args.confirmatory)).lower(),
        "architecture_profile_id": args.architecture,
        "selection_profile_id": args.selection,
        "fitness_profile_id": args.fitness,
        "initial_population_id": args.initial_population,
        "config_path": relative_config,
        "status": args.status,
        "notes": args.notes,
    }
    append_csv(layout.cells_registry, CELL_FIELDS, row)
    for name, value in parameters.items():
        append_csv(
            layout.parameters_registry,
            PARAMETER_FIELDS,
            {
                "cell_id": cell_id,
                "parameter_name": name,
                "value": value,
                "value_type": types.get(name, infer_value_type(value)),
                "unit": units.get(name, ""),
                "role": roles.get(name, "input"),
            },
        )

    scratch_cell.mkdir(parents=True, exist_ok=False)
    cell_metadata = {
        **row,
        "parameters": parameters,
        "scratch_path": str(scratch_cell),
    }
    with (scratch_cell / "cell.json").open("x", encoding="utf-8") as handle:
        json.dump(cell_metadata, handle, indent=2, sort_keys=True)
        handle.write("\n")
    return row


def parameter_sort_key(row: dict[str, str]) -> tuple[int, str]:
    try:
        priority = PREFERRED_PARAMETER_ORDER.index(row["parameter_name"])
    except ValueError:
        priority = len(PREFERRED_PARAMETER_ORDER)
    return priority, row["parameter_name"].lower()


def get_cell(layout: LayoutPaths, requested: str) -> dict[str, Any]:
    cell_rows = read_csv(layout.cells_registry)
    short_reference = requested.strip().lower()
    if SHORT_CELL_ID_RE.fullmatch(short_reference):
        matches = [
            row
            for row in cell_rows
            if row["cell_id"].endswith(f"-{short_reference}")
        ]
        if not matches:
            raise LayoutError(f"cell is not registered: {short_reference}")
        if len(matches) > 1:
            matching_ids = ", ".join(row["cell_id"] for row in matches)
            raise LayoutError(
                f"cell number {short_reference} is ambiguous; use one of: "
                f"{matching_ids}"
            )
        canonical = matches[0]["cell_id"]
    else:
        canonical, _, _ = parse_cell_id(requested)
        matches = [row for row in cell_rows if row["cell_id"] == canonical]
        if not matches:
            raise LayoutError(f"cell is not registered: {canonical}")
    parameter_rows = [
        row
        for row in read_csv(layout.parameters_registry)
        if row["cell_id"] == canonical
    ]
    parameter_rows.sort(key=parameter_sort_key)
    cell = matches[0]
    scratch_path = (
        layout.scratch
        / phase_dirname(cell["phase_id"])
        / stage_dirname(cell["stage_id"])
        / cell["cell_dirname"]
    )
    return {
        "cell": cell,
        "parameters": parameter_rows,
        "paths": {
            "scratch": str(scratch_path),
            "config": str(layout.work / cell["config_path"]),
        },
    }


def get_mnemonic(layout: LayoutPaths, requested: str) -> str:
    """Return the registered mnemonic for a cell reference."""

    return str(get_cell(layout, requested)["cell"]["mnemonic"])


def get_cell_ids_from_mnemonic(
    layout: LayoutPaths,
    requested: str,
    phase_id: str | None = None,
    stage_id: str | None = None,
) -> list[str]:
    """Return all cell IDs matching an exact mnemonic and optional filters."""

    normalized = slug(requested)
    matches = [
        row
        for row in read_csv(layout.cells_registry)
        if row["mnemonic"] == normalized
        and (phase_id is None or row["phase_id"] == phase_id)
        and (stage_id is None or row["stage_id"] == stage_id)
    ]
    if not matches:
        filters = ""
        if phase_id or stage_id:
            filters = (
                f" with phase={phase_id or '*'} and stage={stage_id or '*'}"
            )
        raise LayoutError(f"mnemonic is not registered{filters}: {normalized}")
    return [row["cell_id"] for row in matches]


def print_cell_human(payload: dict[str, Any]) -> None:
    cell = payload["cell"]
    print(f"Cell:       {cell['cell_id']}")
    print(f"Mnemonic:   {cell['mnemonic']}")
    print(f"Directory:  {cell['cell_dirname']}")
    print(f"Label:      {cell['label'] or '-'}")
    print(f"Phase:      {cell['phase_id']}")
    print(f"Stage:      {cell['stage_id']}")
    print(f"Group:      {cell['experimental_group'] or '-'}")
    print(f"Status:     {cell['status'] or '-'}")
    print(f"Confirmatory: {cell['confirmatory']}")
    print(f"Architecture: {cell['architecture_profile_id'] or '-'}")
    print(f"Selection:    {cell['selection_profile_id'] or '-'}")
    print(f"Fitness:      {cell['fitness_profile_id'] or '-'}")
    print(f"Initial pop.: {cell['initial_population_id'] or '-'}")
    print(f"Config:      {payload['paths']['config']}")
    print(f"Scratch:     {payload['paths']['scratch']}")
    print("Parameters:")
    if not payload["parameters"]:
        print("  (none registered)")
    for parameter in payload["parameters"]:
        suffix = f" {parameter['unit']}" if parameter["unit"] else ""
        role = f" [{parameter['role']}]" if parameter["role"] else ""
        print(
            f"  {parameter['parameter_name']:<30} "
            f"{parameter['value']}{suffix}{role}"
        )
    if cell["notes"]:
        print(f"Notes: {cell['notes']}")


def list_cells(args: argparse.Namespace, layout: LayoutPaths) -> None:
    rows = read_csv(layout.cells_registry)
    for row in rows:
        if args.phase and row["phase_id"] != args.phase:
            continue
        if args.stage and row["stage_id"] != args.stage:
            continue
        if args.status and row["status"] != args.status:
            continue
        haystack = " ".join(row.values()).lower()
        if args.contains and args.contains.lower() not in haystack:
            continue
        print(
            f"{row['cell_id']:<16} {row['mnemonic']:<55} "
            f"{row['status']:<12} {row['label']}"
        )


def add_layout_arguments(parser: argparse.ArgumentParser) -> None:
    parser.add_argument(
        "--root",
        required=True,
        help="Base path containing work/, scratch/, and data/ (for example $HOME)",
    )
    parser.add_argument("--project", default=DEFAULT_PROJECT)
    parser.add_argument(
        "--work-root",
        help="Optional parent for the project work directory; defaults to ROOT/work",
    )
    parser.add_argument(
        "--scratch-root",
        help="Optional parent for scratch results; defaults to ROOT/scratch",
    )
    parser.add_argument(
        "--data-root",
        help="Optional parent for retained data; defaults to ROOT/data",
    )


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Initialize and query a trophosome experiment project layout"
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    init_parser = subparsers.add_parser(
        "init", help="Create the work/scratch/data hierarchy and empty registries"
    )
    add_layout_arguments(init_parser)
    init_parser.add_argument(
        "--phase",
        action="append",
        choices=[phase_id for phase_id, _ in DEFAULT_PHASES],
        help="Phase to initialize; may be repeated and defaults to p01 only",
    )

    register_parser = subparsers.add_parser(
        "register-cell", help="Register a scientific parameter cell"
    )
    add_layout_arguments(register_parser)
    register_parser.add_argument("--cell-id", required=True)
    register_parser.add_argument("--label", default="")
    register_parser.add_argument("--mnemonic", default="")
    register_parser.add_argument("--group", default="")
    register_parser.add_argument("--comparison-set", default="")
    register_parser.add_argument("--confirmatory", action="store_true")
    register_parser.add_argument("--architecture", default="")
    register_parser.add_argument("--selection", default="")
    register_parser.add_argument("--fitness", default="")
    register_parser.add_argument("--initial-population", default="")
    register_parser.add_argument("--status", default="planned")
    register_parser.add_argument("--notes", default="")
    register_parser.add_argument(
        "--param",
        action="append",
        default=[],
        metavar="NAME=VALUE",
        help="Cell parameter; may be repeated and may use future parameter names",
    )
    register_parser.add_argument(
        "--unit",
        action="append",
        default=[],
        metavar="NAME=UNIT",
        help="Optional unit for a registered parameter",
    )
    register_parser.add_argument(
        "--role",
        action="append",
        default=[],
        metavar="NAME=ROLE",
        help="Optional role, such as input, derived, nuisance, or technical",
    )
    register_parser.add_argument(
        "--value-type",
        action="append",
        default=[],
        metavar="NAME=TYPE",
        help="Optional explicit value type; otherwise inferred",
    )

    show_parser = subparsers.add_parser(
        "show-cell", help="Print all registered metadata and parameters for a cell"
    )
    add_layout_arguments(show_parser)
    show_parser.add_argument(
        "cell_id",
        help="Cell number, full cell ID, cell directory path, or run ID",
    )
    show_parser.add_argument("--json", action="store_true", dest="as_json")

    mnemonic_parser = subparsers.add_parser(
        "mnemonic", help="Return the mnemonic registered for a cell ID"
    )
    add_layout_arguments(mnemonic_parser)
    mnemonic_parser.add_argument(
        "cell_id",
        help="Cell number, full cell ID, cell directory path, or run ID",
    )

    cell_id_parser = subparsers.add_parser(
        "cell-id", help="Return cell IDs registered for an exact mnemonic"
    )
    add_layout_arguments(cell_id_parser)
    cell_id_parser.add_argument("mnemonic")
    cell_id_parser.add_argument(
        "--phase", help="Optional phase filter, for example p01"
    )
    cell_id_parser.add_argument(
        "--stage", help="Optional stage filter, for example s03"
    )

    list_parser = subparsers.add_parser(
        "list-cells", help="List and filter registered cells"
    )
    add_layout_arguments(list_parser)
    list_parser.add_argument("--phase")
    list_parser.add_argument("--stage")
    list_parser.add_argument("--status")
    list_parser.add_argument("--contains")
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    layout = resolve_layout(args)
    try:
        if args.command == "init":
            initialize_layout(
                layout,
                phase_ids=(
                    tuple(args.phase) if args.phase else DEFAULT_INITIAL_PHASE_IDS
                ),
            )
            print(f"Work:    {layout.work}")
            print(f"Scratch: {layout.scratch}")
            print(f"Data:    {layout.data}")
        elif args.command == "register-cell":
            row = register_cell(args, layout)
            print(f"Registered {row['cell_id']}")
            print(f"Mnemonic:  {row['mnemonic']}")
            print(f"Directory: {row['cell_dirname']}")
        elif args.command == "show-cell":
            payload = get_cell(layout, args.cell_id)
            if args.as_json:
                json.dump(payload, sys.stdout, indent=2)
                sys.stdout.write("\n")
            else:
                print_cell_human(payload)
        elif args.command == "mnemonic":
            print(get_mnemonic(layout, args.cell_id))
        elif args.command == "cell-id":
            cell_ids = get_cell_ids_from_mnemonic(
                layout,
                args.mnemonic,
                phase_id=args.phase,
                stage_id=args.stage,
            )
            for cell_id in cell_ids:
                print(cell_id)
        elif args.command == "list-cells":
            list_cells(args, layout)
        else:  # pragma: no cover - argparse prevents this branch
            parser.error(f"unknown command: {args.command}")
    except (LayoutError, OSError, csv.Error) as exc:
        parser.exit(2, f"error: {exc}\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
