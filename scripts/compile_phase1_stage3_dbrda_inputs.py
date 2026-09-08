#!/usr/bin/env python3
"""Compile Stage 3 community-analysis inputs for db-RDA, RDA and PRC."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import io
import json
import math
import os
import time
from collections import defaultdict
from collections.abc import Iterable, Iterator, Sequence
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Literal

import numpy as np

PASSAGE = 100
TRAJECTORY_PASSAGES = tuple(range(PASSAGE + 1))
ANCESTRAL_LINEAGES = tuple(range(100))
ANALYSIS_DIRECTORY = "s03-parameter-map-community-g100-derived"
WAVE3_ANALYSIS_DIRECTORY = "s03-parameter-map-community-wave3-g100-derived"
WAVE1_EXPERIMENT = "phase1-stage3-wave1-v210-m010-g100"
WAVE2_EXPERIMENT = "phase1-stage3-wave2-v210-adaptive-g1000"
WAVE3_EXPERIMENT = "phase1-stage3-wave3-v210-bridge-g100"
STAGE2_EXPERIMENT = "phase1-second-pilot-v210-m010-g250"
SEED_BLOCKS = tuple(f"sb{number:04d}" for number in range(1, 13))
WAVE3_SEED_BLOCKS = SEED_BLOCKS[:6]
CACHE_SCHEMA_VERSION = "1.0.0"
DEFAULT_WORKERS = 4
DEFAULT_PROGRESS_EVERY = 10
DEFAULT_GZIP_LEVEL = 1
CompileMode = Literal["community", "endpoint", "prc"]

X_FIELDS = (
    "sample_id",
    "analysis_set",
    "panel",
    "analysis_role",
    "include_primary_dbrda",
    "cell_id",
    "cell",
    "seed_block_id",
    "master_seed",
    "source_role",
    "source_cell_id",
    "source_run_id",
    "source_alias_count",
    "passage",
    "H",
    "B",
    "HB",
    "log10_H",
    "log10_B",
    "log10_HB",
    "f",
    "e",
    "R",
    "alpha",
    "alpha_target",
    "m",
    "u",
    "mutation_enabled",
    "host_return_enabled",
    "immediate_host_signal",
)

PRC_FIELDS = (
    "trajectory_sample_id",
    "population_sample_id",
    *X_FIELDS[1:],
)


@dataclass(frozen=True)
class Sample:
    analysis_set: str
    panel: str
    analysis_role: str
    include_primary_dbrda: bool
    cell: dict[str, str]
    seed_block_id: str
    source_role: str
    source_cell_id: str
    source_run: dict[str, str]

    @property
    def sample_id(self) -> str:
        return (
            f"{self.analysis_set}__{self.cell['cell']}__"
            f"{self.seed_block_id}__g{PASSAGE}"
        )


@dataclass(frozen=True)
class SourceTask:
    sample: Sample
    work: Path
    scratch: Path
    capacity: int
    kind: Literal["endpoint", "trajectory"]
    cache_directory: Path | None
    rebuild_cache: bool


def _read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _format_value(value: Any) -> Any:
    if isinstance(value, bool):
        return str(value).lower()
    if isinstance(value, float):
        if not math.isfinite(value):
            raise ValueError("analysis tables cannot contain non-finite values")
        return format(value, ".15g")
    return value


def _atomic_tsv(
    path: Path, rows: Iterable[Sequence[Any]], fields: Sequence[str]
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    with temporary.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(fields)
        for row in rows:
            writer.writerow([_format_value(value) for value in row])
    temporary.replace(path)


def _atomic_json(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    temporary.replace(path)


def _manifest_index(path: Path) -> dict[tuple[str, str], dict[str, str]]:
    rows = _read_tsv(path)
    index = {(row["cell_id"], row["seed_block_id"]): row for row in rows}
    if len(index) != len(rows):
        raise ValueError(f"duplicate cell/seed rows in {path.name}")
    return index


def _resolve_samples(
    phase: Path, *, through_wave: Literal[2, 3] = 2
) -> dict[str, list[Sample]]:
    if through_wave not in {2, 3}:
        raise ValueError("through_wave must be 2 or 3")
    design = phase / "design"
    manifests = phase / "manifests"
    wave1_cells = _read_tsv(design / f"{WAVE1_EXPERIMENT}-cells.tsv")
    wave2_cells = _read_tsv(design / f"{WAVE2_EXPERIMENT}-cells.tsv")
    wave1_runs = _manifest_index(manifests / f"{WAVE1_EXPERIMENT}-runs.tsv")
    wave2_runs = _manifest_index(manifests / f"{WAVE2_EXPERIMENT}-runs.tsv")
    stage2_runs = _manifest_index(manifests / f"{STAGE2_EXPERIMENT}-runs.tsv")

    result: dict[str, list[Sample]] = defaultdict(list)
    for cell in wave1_cells:
        for seed in SEED_BLOCKS:
            if cell["source_role"] == "new-grid-cell":
                source = wave1_runs[cell["cell_id"], seed]
                source_role = "wave1-simulation"
                analysis_role = "factorial-treatment"
                include_primary = True
            else:
                source = stage2_runs[cell["cell_id"], seed]
                source_role = "reused-stage2-no-return-control"
                analysis_role = "paired-no-return-control"
                include_primary = False
            result["wave1_h_alpha_u"].append(
                Sample(
                    analysis_set="wave1_h_alpha_u",
                    panel="H-by-alpha-by-u",
                    analysis_role=analysis_role,
                    include_primary_dbrda=include_primary,
                    cell=cell,
                    seed_block_id=seed,
                    source_role=source_role,
                    source_cell_id=cell["cell_id"],
                    source_run=source,
                )
            )

    for cell in wave2_cells:
        analysis_set = (
            "wave2a_h_by_b" if cell["panel"] == "H-by-B" else "wave2b_alpha_by_m"
        )
        source_cell_id = cell["reused_source_cell_id"] or cell["cell_id"]
        for seed in SEED_BLOCKS:
            if not cell["reused_source_cell_id"]:
                source = wave2_runs[source_cell_id, seed]
                role = "wave2-simulation"
            elif source_cell_id.startswith("p01-s02-"):
                source = stage2_runs[source_cell_id, seed]
                role = "reused-stage2-simulation"
            elif source_cell_id.startswith("p01-s03-"):
                source = wave1_runs[source_cell_id, seed]
                role = "reused-wave1-simulation"
            else:
                raise ValueError(f"unknown reused source cell: {source_cell_id}")
            result[analysis_set].append(
                Sample(
                    analysis_set=analysis_set,
                    panel=cell["panel"],
                    analysis_role="factorial-treatment",
                    include_primary_dbrda=True,
                    cell=cell,
                    seed_block_id=seed,
                    source_role=role,
                    source_cell_id=source_cell_id,
                    source_run=source,
                )
            )

    if through_wave == 3:
        wave3_cells = _read_tsv(design / f"{WAVE3_EXPERIMENT}-cells.tsv")
        wave3_runs = _manifest_index(manifests / f"{WAVE3_EXPERIMENT}-runs.tsv")
        for cell in wave3_cells:
            for seed in WAVE3_SEED_BLOCKS:
                source = wave3_runs[cell["cell_id"], seed]
                result["wave3_bridge"].append(
                    Sample(
                        analysis_set="wave3_bridge",
                        panel=cell["panel"],
                        analysis_role=(
                            "bridge-extension-treatment"
                            if cell["design_tier"] == "two-cell-extension"
                            else "bridge-core-treatment"
                        ),
                        include_primary_dbrda=True,
                        cell=cell,
                        seed_block_id=seed,
                        source_role="wave3-simulation",
                        source_cell_id=cell["cell_id"],
                        source_run=source,
                    )
                )

    expected = {
        "wave1_h_alpha_u": 25 * len(SEED_BLOCKS),
        "wave2a_h_by_b": 12 * len(SEED_BLOCKS),
        "wave2b_alpha_by_m": 28 * len(SEED_BLOCKS),
    }
    if through_wave == 3:
        expected["wave3_bridge"] = 16 * len(WAVE3_SEED_BLOCKS)
    observed = {name: len(samples) for name, samples in result.items()}
    if observed != expected:
        raise ValueError(
            f"unexpected analysis-set sizes: {observed}, expected {expected}"
        )
    for name, samples in result.items():
        identifiers = [sample.sample_id for sample in samples]
        if len(set(identifiers)) != len(identifiers):
            raise ValueError(f"duplicate sample IDs in {name}")
    return dict(result)


def _required_csv_indices(
    header: Sequence[str], required: Sequence[str], *, description: str
) -> tuple[int, ...]:
    try:
        return tuple(header.index(field) for field in required)
    except ValueError as exc:
        raise ValueError(f"{description} has an unexpected header") from exc


def _environment_at_passage(path: Path, passage: int) -> dict[int, int]:
    result: dict[int, int] = {}
    previous_generation = -1
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.reader(handle)
        try:
            header = next(reader)
        except StopIteration as exc:
            raise ValueError(f"environmental table is empty: {path}") from exc
        replicate_index, generation_index, strain_index, count_index = (
            _required_csv_indices(
                header,
                ("replicate", "generation", "strain_id", "count"),
                description="environmental table",
            )
        )
        maximum_index = max(
            replicate_index, generation_index, strain_index, count_index
        )
        for row in reader:
            if not row:
                continue
            if len(row) <= maximum_index:
                raise ValueError("environmental table contains a truncated row")
            if int(row[replicate_index]) != 0:
                raise ValueError("environmental table contains an unexpected replicate")
            generation = int(row[generation_index])
            if generation < previous_generation:
                raise ValueError("environmental rows are not ordered by generation")
            previous_generation = generation
            if generation < passage:
                continue
            if generation > passage:
                break
            strain = int(row[strain_index])
            count = int(row[count_index])
            if strain < 0 or count <= 0:
                raise ValueError("environmental table contains an invalid strain count")
            if strain in result:
                raise ValueError("duplicate strain at the requested passage")
            result[strain] = count
    if not result:
        raise ValueError(f"environmental table does not contain passage {passage}")
    return result


def _environment_trajectory(
    path: Path, *, first_passage: int = 0, last_passage: int = PASSAGE
) -> dict[int, dict[int, int]]:
    """Read one complete, ordered environmental trajectory from a source run."""
    expected = tuple(range(first_passage, last_passage + 1))
    result: dict[int, dict[int, int]] = {passage: {} for passage in expected}
    previous_generation = -1
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.reader(handle)
        try:
            header = next(reader)
        except StopIteration as exc:
            raise ValueError(f"environmental table is empty: {path}") from exc
        replicate_index, generation_index, strain_index, count_index = (
            _required_csv_indices(
                header,
                ("replicate", "generation", "strain_id", "count"),
                description="environmental table",
            )
        )
        maximum_index = max(
            replicate_index, generation_index, strain_index, count_index
        )
        for row in reader:
            if not row:
                continue
            if len(row) <= maximum_index:
                raise ValueError("environmental table contains a truncated row")
            if int(row[replicate_index]) != 0:
                raise ValueError("environmental table contains an unexpected replicate")
            generation = int(row[generation_index])
            if generation < previous_generation:
                raise ValueError("environmental rows are not ordered by generation")
            previous_generation = generation
            if generation < first_passage:
                continue
            if generation > last_passage:
                break
            strain = int(row[strain_index])
            count = int(row[count_index])
            if strain < 0 or count <= 0:
                raise ValueError("environmental table contains an invalid strain count")
            counts = result[generation]
            if strain in counts:
                raise ValueError(
                    f"duplicate strain at environmental passage {generation}"
                )
            counts[strain] = count
    missing = [passage for passage in expected if not result[passage]]
    if missing:
        preview = ", ".join(str(value) for value in missing[:10])
        suffix = " ..." if len(missing) > 10 else ""
        raise ValueError(
            "environmental trajectory is incomplete; missing passage(s) "
            f"{preview}{suffix}"
        )
    return result


def _reversed_binary_lines(
    path: Path, block_size: int = 1024 * 1024
) -> Iterator[bytes]:
    """Yield non-empty data lines from a numeric CSV in reverse file order."""
    with path.open("rb") as handle:
        handle.seek(0, 2)
        position = handle.tell()
        remainder = b""
        while position > 0:
            amount = min(block_size, position)
            position -= amount
            handle.seek(position)
            chunk = handle.read(amount) + remainder
            lines = chunk.split(b"\n")
            remainder = lines[0]
            for line in reversed(lines[1:]):
                if line:
                    yield line.rstrip(b"\r")
        if remainder:
            yield remainder.rstrip(b"\r")


def _retained_parents(path: Path, final_ids: Iterable[int]) -> dict[int, int]:
    needed = {strain for strain in final_ids if strain not in ANCESTRAL_LINEAGES}
    if not needed:
        return {}
    with path.open("rb") as handle:
        header = handle.readline().decode("utf-8").rstrip("\r\n").split(",")
    try:
        child_index = header.index("strain_id")
        parent_index = header.index("parent_strain_id")
    except ValueError as exc:
        raise ValueError("lineage-event table has an unexpected header") from exc
    parents: dict[int, int] = {}
    for raw_line in _reversed_binary_lines(path):
        fields = raw_line.split(b",")
        try:
            child = int(fields[child_index])
        except (IndexError, ValueError):
            continue
        if child not in needed:
            continue
        try:
            parent = int(fields[parent_index])
        except (IndexError, ValueError) as exc:
            raise ValueError(f"invalid parent for retained strain {child}") from exc
        parents[child] = parent
        if parent not in ANCESTRAL_LINEAGES:
            needed.add(parent)
    missing = needed - set(parents) - set(ANCESTRAL_LINEAGES)
    if missing:
        example = min(missing)
        raise ValueError(f"lineage parentage is missing for retained strain {example}")
    return parents


def _root(strain: int, parents: dict[int, int], cache: dict[int, int]) -> int:
    if strain in cache:
        return cache[strain]
    path: list[int] = []
    current = strain
    while current not in cache:
        path.append(current)
        if current not in parents:
            raise ValueError(
                f"lineage parentage is missing for retained strain {current}"
            )
        current = parents[current]
    root = cache[current]
    for descendant in path:
        cache[descendant] = root
    return root


def _committed_output_check(output: Path, filenames: Iterable[str]) -> str:
    marker = output / "completion.json"
    if not marker.is_file():
        marker = output / "pause.json"
    if not marker.is_file():
        raise ValueError("completion/pause marker is missing")
    payload = json.loads(marker.read_text(encoding="utf-8"))
    sizes = payload.get("output_sizes", {})
    for filename in filenames:
        path = output / filename
        if not path.is_file():
            raise ValueError(f"required source table is missing: {filename}")
        if filename in sizes and path.stat().st_size != int(sizes[filename]):
            raise ValueError(f"committed size differs for {filename}")
    return marker.name


def _ancestral_trajectory(
    sample: Sample, *, work: Path, scratch: Path, capacity: int
) -> tuple[np.ndarray, dict[str, Any]]:
    """Return passages 0--100 collapsed to the common ancestral lineages."""
    source = sample.source_run
    config = work / source["config_path"]
    if _sha256(config) != source["config_sha256"]:
        raise ValueError("source configuration checksum differs from its manifest")
    output = scratch / source["scratch_relative_path"]
    marker = _committed_output_check(
        output, ("environment_counts.csv", "strain_lineage_events.csv")
    )
    counts_by_passage = _environment_trajectory(
        output / "environment_counts.csv",
        first_passage=TRAJECTORY_PASSAGES[0],
        last_passage=TRAJECTORY_PASSAGES[-1],
    )
    retained_strains = {
        strain for counts in counts_by_passage.values() for strain in counts
    }
    parents = _retained_parents(output / "strain_lineage_events.csv", retained_strains)
    root_cache = {lineage: lineage for lineage in ANCESTRAL_LINEAGES}
    collapsed = np.zeros(
        (len(TRAJECTORY_PASSAGES), len(ANCESTRAL_LINEAGES)), dtype=np.int64
    )
    for row_index, passage in enumerate(TRAJECTORY_PASSAGES):
        counts = counts_by_passage[passage]
        observed_capacity = sum(counts.values())
        if observed_capacity != capacity:
            raise ValueError(
                f"passage-{passage} capacity is {observed_capacity}, "
                f"expected {capacity}"
            )
        for strain, count in counts.items():
            root = _root(strain, parents, root_cache)
            if root not in ANCESTRAL_LINEAGES:
                raise ValueError(
                    f"strain {strain} traces to invalid ancestral root {root}"
                )
            collapsed[row_index, root] += count
        if int(collapsed[row_index].sum()) != capacity:
            raise ValueError(
                f"ancestral collapse does not preserve capacity at passage {passage}"
            )
    frequencies = collapsed.astype(np.float64) / capacity
    if not np.allclose(frequencies.sum(axis=1), 1.0, rtol=0, atol=1e-12):
        raise ValueError("ancestral trajectory frequencies do not sum to one")
    endpoint_counts = counts_by_passage[PASSAGE]
    endpoint_collapsed = collapsed[PASSAGE]
    return frequencies, {
        "source_run_id": source["run_id"],
        "source_cell_id": source["cell_id"],
        "seed_block_id": source["seed_block_id"],
        "scratch_relative_path": source["scratch_relative_path"],
        "state_marker": marker,
        "first_passage": TRAJECTORY_PASSAGES[0],
        "last_passage": TRAJECTORY_PASSAGES[-1],
        "passages": len(TRAJECTORY_PASSAGES),
        "raw_strain_richness_g100": len(endpoint_counts),
        "ancestral_lineage_richness_g100": int(np.count_nonzero(endpoint_collapsed)),
        "retained_mutant_strains_g100": sum(
            strain not in ANCESTRAL_LINEAGES for strain in endpoint_counts
        ),
        "configuration_sha256": source["config_sha256"],
    }


def _ancestral_frequencies(
    sample: Sample, *, work: Path, scratch: Path, capacity: int
) -> tuple[np.ndarray, dict[str, Any]]:
    source = sample.source_run
    config = work / source["config_path"]
    if _sha256(config) != source["config_sha256"]:
        raise ValueError("source configuration checksum differs from its manifest")
    output = scratch / source["scratch_relative_path"]
    marker = _committed_output_check(
        output, ("environment_counts.csv", "strain_lineage_events.csv")
    )
    counts = _environment_at_passage(output / "environment_counts.csv", PASSAGE)
    observed_capacity = sum(counts.values())
    if observed_capacity != capacity:
        raise ValueError(
            f"passage-{PASSAGE} capacity is {observed_capacity}, expected {capacity}"
        )
    parents = _retained_parents(output / "strain_lineage_events.csv", counts)
    root_cache = {lineage: lineage for lineage in ANCESTRAL_LINEAGES}
    collapsed = np.zeros(len(ANCESTRAL_LINEAGES), dtype=np.int64)
    for strain, count in counts.items():
        root = _root(strain, parents, root_cache)
        if root not in ANCESTRAL_LINEAGES:
            raise ValueError(f"strain {strain} traces to invalid ancestral root {root}")
        collapsed[root] += count
    if int(collapsed.sum()) != capacity:
        raise ValueError("ancestral collapse does not preserve environmental capacity")
    frequencies = collapsed.astype(np.float64) / capacity
    if not math.isclose(float(frequencies.sum()), 1.0, abs_tol=1e-12):
        raise ValueError("ancestral frequencies do not sum to one")
    return frequencies, {
        "source_run_id": source["run_id"],
        "source_cell_id": source["cell_id"],
        "seed_block_id": source["seed_block_id"],
        "scratch_relative_path": source["scratch_relative_path"],
        "state_marker": marker,
        "first_passage": PASSAGE,
        "last_passage": PASSAGE,
        "passages": 1,
        "raw_strain_richness_g100": len(counts),
        "ancestral_lineage_richness_g100": int(np.count_nonzero(collapsed)),
        "retained_mutant_strains_g100": sum(
            strain not in ANCESTRAL_LINEAGES for strain in counts
        ),
        "configuration_sha256": source["config_sha256"],
    }


def _source_signature(
    sample: Sample, *, work: Path, scratch: Path, capacity: int
) -> dict[str, Any]:
    """Build a cheap signature that invalidates stale per-source caches."""
    source = sample.source_run
    config = work / source["config_path"]
    if _sha256(config) != source["config_sha256"]:
        raise ValueError("source configuration checksum differs from its manifest")
    output = scratch / source["scratch_relative_path"]
    filenames = ("environment_counts.csv", "strain_lineage_events.csv")
    marker_name = _committed_output_check(output, filenames)
    environment_stat = (output / filenames[0]).stat()
    lineage_stat = (output / filenames[1]).stat()
    return {
        "cache_schema_version": CACHE_SCHEMA_VERSION,
        "source_run_id": source["run_id"],
        "configuration_sha256": source["config_sha256"],
        "state_marker": marker_name,
        "state_marker_sha256": _sha256(output / marker_name),
        "environment_counts_size": environment_stat.st_size,
        "environment_counts_mtime_ns": environment_stat.st_mtime_ns,
        "strain_lineage_events_size": lineage_stat.st_size,
        "strain_lineage_events_mtime_ns": lineage_stat.st_mtime_ns,
        "capacity": capacity,
        "endpoint_passage": PASSAGE,
        "trajectory_first_passage": TRAJECTORY_PASSAGES[0],
        "trajectory_last_passage": TRAJECTORY_PASSAGES[-1],
        "ancestral_lineages": len(ANCESTRAL_LINEAGES),
    }


def _source_cache_path(
    cache_directory: Path, run_id: str, kind: Literal["endpoint", "trajectory"]
) -> Path:
    digest = hashlib.sha256(run_id.encode("utf-8")).hexdigest()
    return cache_directory / kind / f"{digest}.npz"


def _load_source_cache(
    path: Path,
    *,
    signature: dict[str, Any],
    kind: Literal["endpoint", "trajectory"],
) -> tuple[np.ndarray, dict[str, Any]] | None:
    if not path.is_file():
        return None
    try:
        with np.load(path, allow_pickle=False) as archive:
            frequencies = np.asarray(archive["frequencies"], dtype=np.float64)
            metadata = json.loads(str(archive["metadata_json"].item()))
        if metadata.get("kind") != kind:
            return None
        if metadata.get("signature") != signature:
            return None
        provenance = metadata["provenance"]
        expected_shape = (
            (len(ANCESTRAL_LINEAGES),)
            if kind == "endpoint"
            else (len(TRAJECTORY_PASSAGES), len(ANCESTRAL_LINEAGES))
        )
        if frequencies.shape != expected_shape:
            return None
        rows = frequencies.reshape(1, -1) if kind == "endpoint" else frequencies
        if np.any(rows < 0) or not np.allclose(
            rows.sum(axis=1), 1.0, rtol=0, atol=1e-12
        ):
            return None
        return frequencies, provenance
    except (OSError, ValueError, KeyError, json.JSONDecodeError):
        return None


def _write_source_cache(
    path: Path,
    frequencies: np.ndarray,
    provenance: dict[str, Any],
    *,
    signature: dict[str, Any],
    kind: Literal["endpoint", "trajectory"],
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(f".{os.getpid()}.tmp")
    metadata = json.dumps(
        {"kind": kind, "signature": signature, "provenance": provenance},
        sort_keys=True,
        separators=(",", ":"),
    )
    try:
        with temporary.open("wb") as handle:
            # The cache is intentionally uncompressed: it is small compared with
            # the raw CSVs and much faster to load and save after interruption.
            np.savez(handle, frequencies=frequencies, metadata_json=metadata)
        temporary.replace(path)
    finally:
        if temporary.exists():
            temporary.unlink()


def _process_source_task(
    task: SourceTask,
) -> tuple[str, np.ndarray, dict[str, Any], str]:
    sample = task.sample
    run_id = sample.source_run["run_id"]
    signature = _source_signature(
        sample, work=task.work, scratch=task.scratch, capacity=task.capacity
    )
    if task.cache_directory is not None and not task.rebuild_cache:
        requested_path = _source_cache_path(task.cache_directory, run_id, task.kind)
        cached = _load_source_cache(requested_path, signature=signature, kind=task.kind)
        if cached is not None:
            frequencies, provenance = cached
            return run_id, frequencies, provenance, f"{task.kind}-hit"
        if task.kind == "endpoint":
            trajectory_path = _source_cache_path(
                task.cache_directory, run_id, "trajectory"
            )
            trajectory_cache = _load_source_cache(
                trajectory_path, signature=signature, kind="trajectory"
            )
            if trajectory_cache is not None:
                trajectory, provenance = trajectory_cache
                return run_id, trajectory[PASSAGE], provenance, "trajectory-hit"

    if task.kind == "endpoint":
        frequencies, provenance = _ancestral_frequencies(
            sample,
            work=task.work,
            scratch=task.scratch,
            capacity=task.capacity,
        )
    else:
        frequencies, provenance = _ancestral_trajectory(
            sample,
            work=task.work,
            scratch=task.scratch,
            capacity=task.capacity,
        )
    if task.cache_directory is not None:
        cache_path = _source_cache_path(task.cache_directory, run_id, task.kind)
        _write_source_cache(
            cache_path,
            frequencies,
            provenance,
            signature=signature,
            kind=task.kind,
        )
        cache_status = f"{task.kind}-written"
    else:
        cache_status = "disabled"
    return run_id, frequencies, provenance, cache_status


def _x_record(
    sample: Sample, source_alias_counts: dict[str, int], *, passage: int = PASSAGE
) -> dict[str, Any]:
    cell = sample.cell
    h = int(cell["H"])
    b = int(cell.get("B") or 10)
    alpha = float(cell["alpha"])
    migration = float(cell["m"])
    mutation = float(cell["u"])
    return {
        "sample_id": sample.sample_id,
        "analysis_set": sample.analysis_set,
        "panel": sample.panel,
        "analysis_role": sample.analysis_role,
        "include_primary_dbrda": sample.include_primary_dbrda,
        "cell_id": cell["cell_id"],
        "cell": cell["cell"],
        "seed_block_id": sample.seed_block_id,
        "master_seed": int(sample.source_run["master_seed"]),
        "source_role": sample.source_role,
        "source_cell_id": sample.source_cell_id,
        "source_run_id": sample.source_run["run_id"],
        "source_alias_count": source_alias_counts[sample.source_run["run_id"]],
        "passage": passage,
        "H": h,
        "B": b,
        "HB": h * b,
        "log10_H": math.log10(h),
        "log10_B": math.log10(b),
        "log10_HB": math.log10(h * b),
        "f": float(cell["f"]),
        "e": int(cell["e"]),
        "R": int(cell["R"]),
        "alpha": alpha,
        "alpha_target": float(cell["alpha_target"]),
        "m": migration,
        "u": mutation,
        "mutation_enabled": mutation > 0,
        "host_return_enabled": alpha > 0,
        "immediate_host_signal": alpha * (1 - migration),
    }


def _x_row(sample: Sample, source_alias_counts: dict[str, int]) -> list[Any]:
    record = _x_record(sample, source_alias_counts)
    return [record[field] for field in X_FIELDS]


def total_variation_matrix(frequencies: np.ndarray) -> np.ndarray:
    if frequencies.ndim != 2 or frequencies.shape[1] != len(ANCESTRAL_LINEAGES):
        raise ValueError("Y must have one column for every ancestral lineage")
    if np.any(frequencies < 0) or not np.allclose(
        frequencies.sum(axis=1), 1.0, rtol=0, atol=1e-12
    ):
        raise ValueError("every row of Y must be a non-negative frequency vector")
    size = frequencies.shape[0]
    distances = np.zeros((size, size), dtype=np.float64)
    for index in range(size - 1):
        values = 0.5 * np.abs(frequencies[index + 1 :] - frequencies[index]).sum(axis=1)
        distances[index, index + 1 :] = values
        distances[index + 1 :, index] = values
    if np.any(distances < 0) or np.any(distances > 1 + 1e-12):
        raise ValueError("TV distances fall outside zero to one")
    return distances


def _write_master_triplet(
    output: Path,
    samples: list[Sample],
    frequencies: np.ndarray,
) -> dict[str, Any]:
    sample_ids = [sample.sample_id for sample in samples]
    if frequencies.shape != (len(samples), len(ANCESTRAL_LINEAGES)):
        raise ValueError("Y dimensions differ from X")
    distances = total_variation_matrix(frequencies)
    source_alias_counts: dict[str, int] = defaultdict(int)
    for sample in samples:
        source_alias_counts[sample.source_run["run_id"]] += 1
    for analysis_set in {sample.analysis_set for sample in samples}:
        selected = [sample for sample in samples if sample.analysis_set == analysis_set]
        source_ids = [sample.source_run["run_id"] for sample in selected]
        if len(source_ids) != len(set(source_ids)):
            raise ValueError(f"a source population is duplicated within {analysis_set}")
    x_path = output / "x-explanatory-g100.tsv"
    y_path = output / "y-ancestral-frequencies-g100.tsv"
    d_path = output / "yprime-tv-g100.tsv"
    _atomic_tsv(
        x_path,
        (_x_row(sample, source_alias_counts) for sample in samples),
        X_FIELDS,
    )
    _atomic_tsv(
        y_path,
        (
            [sample_id, *frequencies[index].tolist()]
            for index, sample_id in enumerate(sample_ids)
        ),
        ("sample_id", *(f"ancestral_{lineage:03d}" for lineage in ANCESTRAL_LINEAGES)),
    )
    _atomic_tsv(
        d_path,
        (
            [sample_id, *distances[index].tolist()]
            for index, sample_id in enumerate(sample_ids)
        ),
        ("sample_id", *sample_ids),
    )
    return {
        "samples": len(samples),
        "ancestral_lineages": len(ANCESTRAL_LINEAGES),
        "x": {"path": x_path.name, "sha256": _sha256(x_path)},
        "y": {"path": y_path.name, "sha256": _sha256(y_path)},
        "yprime": {"path": d_path.name, "sha256": _sha256(d_path)},
        "maximum_tv": float(distances.max()),
        "reused_analysis_aliases": sum(
            count - 1 for count in source_alias_counts.values()
        ),
    }


def _write_prc_trajectory(
    output: Path,
    samples: list[Sample],
    trajectories: dict[str, np.ndarray],
    *,
    gzip_level: int = DEFAULT_GZIP_LEVEL,
    progress_every: int = 50,
) -> dict[str, Any]:
    """Write a gzip-compressed long table without materializing it in memory."""
    source_alias_counts: dict[str, int] = defaultdict(int)
    for sample in samples:
        source_alias_counts[sample.source_run["run_id"]] += 1
    path = output / "prc-ancestral-trajectories-g0-g100.tsv.gz"
    temporary = path.with_suffix(path.suffix + ".tmp")
    lineage_fields = tuple(f"ancestral_{lineage:03d}" for lineage in ANCESTRAL_LINEAGES)
    row_count = 0
    started = time.monotonic()
    with temporary.open("wb") as raw_handle:
        with gzip.GzipFile(
            filename="",
            mode="wb",
            fileobj=raw_handle,
            compresslevel=gzip_level,
            mtime=0,
        ) as gzip_handle:
            with io.TextIOWrapper(
                gzip_handle, encoding="utf-8", newline=""
            ) as text_handle:
                writer = csv.writer(text_handle, delimiter="\t", lineterminator="\n")
                writer.writerow((*PRC_FIELDS, *lineage_fields))
                for sample_index, sample in enumerate(samples, start=1):
                    run_id = sample.source_run["run_id"]
                    if run_id not in trajectories:
                        raise ValueError(f"PRC trajectory is missing for {run_id}")
                    trajectory = trajectories[run_id]
                    expected_shape = (
                        len(TRAJECTORY_PASSAGES),
                        len(ANCESTRAL_LINEAGES),
                    )
                    if trajectory.shape != expected_shape:
                        raise ValueError(
                            f"PRC trajectory for {run_id} has shape "
                            f"{trajectory.shape}, expected {expected_shape}"
                        )
                    if np.any(trajectory < 0) or not np.allclose(
                        trajectory.sum(axis=1), 1.0, rtol=0, atol=1e-12
                    ):
                        raise ValueError(
                            f"PRC trajectory for {run_id} is not compositional"
                        )
                    for passage in TRAJECTORY_PASSAGES:
                        x_record = _x_record(
                            sample, source_alias_counts, passage=passage
                        )
                        metadata = {
                            "trajectory_sample_id": (
                                f"{sample.analysis_set}__{sample.cell['cell']}__"
                                f"{sample.seed_block_id}__g{passage}"
                            ),
                            "population_sample_id": sample.sample_id,
                            **{field: x_record[field] for field in X_FIELDS[1:]},
                        }
                        writer.writerow(
                            [_format_value(metadata[field]) for field in PRC_FIELDS]
                            + [
                                _format_value(value)
                                for value in trajectory[passage].tolist()
                            ]
                        )
                        row_count += 1
                    if (
                        sample_index == len(samples)
                        or sample_index % progress_every == 0
                    ):
                        elapsed = time.monotonic() - started
                        rate = sample_index / elapsed if elapsed > 0 else 0.0
                        eta = (len(samples) - sample_index) / rate if rate > 0 else 0.0
                        print(
                            f"Wrote {sample_index}/{len(samples)} PRC population "
                            f"trajectories (ETA {_format_duration(eta)}).",
                            flush=True,
                        )
    temporary.replace(path)
    expected_rows = len(samples) * len(TRAJECTORY_PASSAGES)
    if row_count != expected_rows:
        raise ValueError(f"PRC table has {row_count} rows, expected {expected_rows}")
    return {
        "path": path.name,
        "sha256": _sha256(path),
        "rows": row_count,
        "population_trajectories": len(samples),
        "unique_source_trajectories": len(trajectories),
        "passages_per_trajectory": len(TRAJECTORY_PASSAGES),
        "first_passage": TRAJECTORY_PASSAGES[0],
        "last_passage": TRAJECTORY_PASSAGES[-1],
        "ancestral_lineages": len(ANCESTRAL_LINEAGES),
        "gzip_compression_level": gzip_level,
    }


def _format_duration(seconds: float) -> str:
    seconds = max(0, round(seconds))
    hours, remainder = divmod(seconds, 3600)
    minutes, seconds = divmod(remainder, 60)
    if hours:
        return f"{hours:d}h {minutes:02d}m"
    if minutes:
        return f"{minutes:d}m {seconds:02d}s"
    return f"{seconds:d}s"


def _unique_source_samples(samples: Iterable[Sample]) -> list[Sample]:
    unique: dict[str, Sample] = {}
    for sample in samples:
        unique.setdefault(sample.source_run["run_id"], sample)
    return list(unique.values())


def _compile_source_populations(
    samples: list[Sample],
    *,
    work: Path,
    scratch: Path,
    capacity: int,
    kind: Literal["endpoint", "trajectory"],
    workers: int,
    cache_directory: Path | None,
    rebuild_cache: bool,
    progress_every: int,
) -> tuple[
    dict[str, np.ndarray],
    dict[str, dict[str, Any]],
    dict[str, int],
    list[str],
]:
    source_samples = _unique_source_samples(samples)
    tasks = [
        SourceTask(
            sample=sample,
            work=work,
            scratch=scratch,
            capacity=capacity,
            kind=kind,
            cache_directory=cache_directory,
            rebuild_cache=rebuild_cache,
        )
        for sample in source_samples
    ]
    total = len(tasks)
    frequencies: dict[str, np.ndarray] = {}
    provenance: dict[str, dict[str, Any]] = {}
    cache_counts: dict[str, int] = defaultdict(int)
    issues: list[str] = []
    started = time.monotonic()

    def record_result(
        expected_run_id: str,
        result: tuple[str, np.ndarray, dict[str, Any], str] | None,
        error: BaseException | None,
        completed: int,
    ) -> None:
        if error is not None:
            issues.append(f"{expected_run_id}: {type(error).__name__}: {error}")
        elif result is not None:
            run_id, values, source, cache_status = result
            if run_id != expected_run_id:
                issues.append(
                    f"{expected_run_id}: worker returned unexpected source {run_id}"
                )
            else:
                frequencies[run_id] = values
                provenance[run_id] = source
                cache_counts[cache_status] += 1
        if completed == 1 or completed == total or completed % progress_every == 0:
            elapsed = time.monotonic() - started
            rate = completed / elapsed if elapsed > 0 else 0.0
            eta = (total - completed) / rate if rate > 0 else 0.0
            hits = sum(
                count
                for status, count in cache_counts.items()
                if status.endswith("-hit")
            )
            print(
                f"Processed {completed}/{total} source populations "
                f"({rate:.2f}/s; ETA {_format_duration(eta)}; "
                f"cache hits {hits})",
                flush=True,
            )

    if workers == 1:
        for completed, task in enumerate(tasks, start=1):
            run_id = task.sample.source_run["run_id"]
            try:
                result = _process_source_task(task)
            except Exception as exc:  # worker boundary: report every source failure
                record_result(run_id, None, exc, completed)
            else:
                record_result(run_id, result, None, completed)
    else:
        with ProcessPoolExecutor(max_workers=workers) as executor:
            future_to_run = {
                executor.submit(_process_source_task, task): task.sample.source_run[
                    "run_id"
                ]
                for task in tasks
            }
            for completed, future in enumerate(as_completed(future_to_run), start=1):
                run_id = future_to_run[future]
                try:
                    result = future.result()
                except Exception as exc:  # worker boundary: report every source failure
                    record_result(run_id, None, exc, completed)
                else:
                    record_result(run_id, result, None, completed)
    return frequencies, provenance, dict(cache_counts), issues


def compile_inputs(
    repository: Path,
    output: Path | None = None,
    *,
    through_wave: Literal[2, 3] = 2,
    mode: CompileMode = "community",
    workers: int = DEFAULT_WORKERS,
    cache_directory: Path | None = None,
    use_cache: bool = True,
    rebuild_cache: bool = False,
    progress_every: int = DEFAULT_PROGRESS_EVERY,
    gzip_level: int = DEFAULT_GZIP_LEVEL,
) -> Path:
    if through_wave not in {2, 3}:
        raise ValueError("through_wave must be 2 or 3")
    if mode not in {"community", "endpoint", "prc"}:
        raise ValueError(f"unknown compilation mode: {mode}")
    if workers < 1:
        raise ValueError("workers must be at least one")
    if progress_every < 1:
        raise ValueError("progress_every must be at least one")
    if not 0 <= gzip_level <= 9:
        raise ValueError("gzip_level must be between zero and nine")
    if rebuild_cache and not use_cache:
        raise ValueError("rebuild_cache cannot be combined with use_cache=False")
    compilation_started = time.monotonic()
    work = repository / "experiments/work/trophosome"
    phase = work / "p01-neutral-feedback"
    layout = work / "layout.local.json"
    if not layout.is_file():
        raise RuntimeError(
            "machine-local layout is missing; create layout.local.json using "
            "scripts/hpc/README.md"
        )
    scratch = Path(json.loads(layout.read_text(encoding="utf-8"))["scratch"])
    if not scratch.is_absolute():
        raise RuntimeError("the machine-local scratch path must be absolute")
    initial_payload = json.loads(
        (work / "common/initial-populations/ip001-fisher100.json").read_text(
            encoding="utf-8"
        )
    )
    initial_counts = [int(value) for value in initial_payload["scaled_counts"]]
    if len(initial_counts) != len(ANCESTRAL_LINEAGES):
        raise RuntimeError("the frozen starting population does not have 100 lineages")
    capacity = sum(initial_counts)
    default_directory = (
        WAVE3_ANALYSIS_DIRECTORY if through_wave == 3 else ANALYSIS_DIRECTORY
    )
    output = output or phase / "analysis" / default_directory
    output.mkdir(parents=True, exist_ok=True)
    if use_cache:
        cache_directory = cache_directory or (
            scratch
            / "p01-neutral-feedback"
            / "analysis-cache"
            / "s03-community-g100-v1"
        )
        cache_directory.mkdir(parents=True, exist_ok=True)
    else:
        cache_directory = None

    analysis_sets = _resolve_samples(phase, through_wave=through_wave)
    all_samples = [sample for samples in analysis_sets.values() for sample in samples]
    source_kind: Literal["endpoint", "trajectory"] = (
        "endpoint" if mode == "endpoint" else "trajectory"
    )
    unique_sources = _unique_source_samples(all_samples)
    cache_description = str(cache_directory) if cache_directory else "disabled"
    print(
        f"Compiling {mode} inputs from {len(unique_sources)} unique source "
        f"populations with {workers} worker(s).",
        flush=True,
    )
    print(f"Resumable cache: {cache_description}", flush=True)
    source_values, provenance, cache_counts, issues = _compile_source_populations(
        all_samples,
        work=work,
        scratch=scratch,
        capacity=capacity,
        kind=source_kind,
        workers=workers,
        cache_directory=cache_directory,
        rebuild_cache=rebuild_cache,
        progress_every=progress_every,
    )

    audit_path = output / "community-input-audit-g100.json"
    if issues:
        _atomic_json(
            audit_path,
            {
                "status": "FAIL",
                "mode": mode,
                "through_wave": through_wave,
                "passage": PASSAGE,
                "issues": issues,
                "cache_directory": cache_description,
                "note": "No complete requested community-table release was produced.",
            },
        )
        raise RuntimeError(
            "community-input audit failed:\n"
            + "\n".join(issues[:30])
            + (f"\n... and {len(issues) - 30} more" if len(issues) > 30 else "")
        )

    if len(source_values) != len(unique_sources):
        raise RuntimeError("internal source/result count differs after audit")
    triplet: dict[str, Any] | None = None
    if mode in {"community", "endpoint"}:
        print("Writing endpoint X, Y and Y-prime tables.", flush=True)
        frequency_rows = []
        for sample in all_samples:
            values = source_values[sample.source_run["run_id"]]
            frequency_rows.append(
                values if source_kind == "endpoint" else values[PASSAGE]
            )
        triplet = _write_master_triplet(output, all_samples, np.vstack(frequency_rows))
    prc_trajectory: dict[str, Any] | None = None
    if mode in {"community", "prc"}:
        print(
            f"Writing compressed PRC trajectory table at gzip level {gzip_level}.",
            flush=True,
        )
        prc_trajectory = _write_prc_trajectory(
            output,
            all_samples,
            source_values,
            gzip_level=gzip_level,
            progress_every=max(progress_every, 1),
        )

    source_path = output / "community-source-provenance-g100.tsv"
    source_fields = (
        "source_run_id",
        "source_cell_id",
        "seed_block_id",
        "scratch_relative_path",
        "state_marker",
        "first_passage",
        "last_passage",
        "passages",
        "raw_strain_richness_g100",
        "ancestral_lineage_richness_g100",
        "retained_mutant_strains_g100",
        "configuration_sha256",
    )
    _atomic_tsv(
        source_path,
        (
            [provenance[sample.source_run["run_id"]][field] for field in source_fields]
            for sample in unique_sources
        ),
        source_fields,
    )
    cache_hits = sum(
        count for status, count in cache_counts.items() if status.endswith("-hit")
    )
    cache_writes = sum(
        count for status, count in cache_counts.items() if status.endswith("-written")
    )
    payload = {
        "status": "PASS",
        "schema_version": "3.0.0",
        "mode": mode,
        "through_wave": through_wave,
        "passage": PASSAGE,
        "environmental_capacity": capacity,
        "response_definition": (
            "Y contains relative abundances after mutant descendants are collapsed "
            "to the 100 frozen ancestral lineages."
        ),
        "distance_definition": "D_TV(i,j) = 0.5 * sum_k(abs(Y[i,k] - Y[j,k])).",
        "bray_curtis_equivalence": (
            "Because every row of Y sums to one, D_TV equals Bray-Curtis dissimilarity."
        ),
        "master_triplet": triplet,
        "prc_trajectory": prc_trajectory,
        "performance": {
            "workers": workers,
            "elapsed_seconds": time.monotonic() - compilation_started,
            "source_kind": source_kind,
            "source_populations": len(unique_sources),
            "cache_directory": cache_description,
            "cache_rebuild_requested": rebuild_cache,
            "cache_hits": cache_hits,
            "cache_writes": cache_writes,
            "cache_status_counts": cache_counts,
        },
        "analysis_sets": [
            {
                "analysis_set": name,
                "samples_in_master": len(samples),
                "primary_dbrda_samples": sum(
                    sample.include_primary_dbrda for sample in samples
                ),
            }
            for name, samples in analysis_sets.items()
        ],
        "unique_source_populations": len(provenance),
        "mandatory_analysis_rule": (
            "Subset X by analysis_set and include_primary_dbrda, then use the "
            "same sample IDs to subset Y and both axes of Yprime. Do not fit a "
            "single db-RDA to the unfiltered master matrix."
        ),
        "wave3_bridge_rule": (
            "For the predeclared bridge-augmented model, retain sb0001-sb0006, "
            "u=0 and primary rows from Waves 1-3, deduplicate source_run_id, "
            "then subset Y and both axes of Yprime by the retained sample_id."
            if through_wave == 3
            else None
        ),
        "mandatory_prc_rule": (
            "Subset the trajectory table to one biologically coherent contrast "
            "with one declared reference treatment. Retain every passage and "
            "permute whole population trajectories within matched seed blocks."
        ),
        "source_provenance": {
            "path": source_path.name,
            "sha256": _sha256(source_path),
        },
    }
    _atomic_json(audit_path, payload)
    elapsed = time.monotonic() - compilation_started
    print(
        f"Compiled {mode} inputs for {len(all_samples)} analysis populations "
        f"from {len(provenance)} unique sources in {_format_duration(elapsed)}.\n"
        f"Portable outputs: {output}",
        flush=True,
    )
    return output


def main() -> int:
    parser = argparse.ArgumentParser(
        description=__doc__,
        epilog=(
            "This command reads existing HPC scratch outputs and never launches "
            "or changes a simulation."
        ),
    )
    parser.add_argument(
        "--repository", type=Path, default=Path(__file__).resolve().parents[1]
    )
    parser.add_argument(
        "--output",
        type=Path,
        help=(
            "optional output directory; defaults to the portable Phase 1 analysis tree"
        ),
    )
    parser.add_argument(
        "--through-wave",
        type=int,
        choices=(2, 3),
        default=2,
        help=(
            "latest completed Stage 3 wave included in the master tables; "
            "default: %(default)s"
        ),
    )
    requested = parser.add_mutually_exclusive_group()
    requested.add_argument(
        "--endpoint-only",
        action="store_true",
        help="compile passage-100 X, Y and Y-prime without the PRC trajectory",
    )
    requested.add_argument(
        "--prc-only",
        action="store_true",
        help="compile the passage 0-100 PRC trajectory without endpoint matrices",
    )
    parser.add_argument(
        "--workers",
        type=int,
        default=os.environ.get("TROPHOSOME_COMMUNITY_WORKERS", str(DEFAULT_WORKERS)),
        help=(
            "source populations processed concurrently; default: %(default)s "
            "(or TROPHOSOME_COMMUNITY_WORKERS)"
        ),
    )
    parser.add_argument(
        "--cache-directory",
        type=Path,
        help="optional machine-local per-source cache location",
    )
    cache = parser.add_mutually_exclusive_group()
    cache.add_argument(
        "--no-cache",
        action="store_true",
        help="disable resumable per-source caching",
    )
    cache.add_argument(
        "--rebuild-cache",
        action="store_true",
        help="ignore valid cached sources and replace them after recomputation",
    )
    parser.add_argument(
        "--progress-every",
        type=int,
        default=DEFAULT_PROGRESS_EVERY,
        metavar="N",
        help="print progress after every N source populations; default: %(default)s",
    )
    parser.add_argument(
        "--gzip-level",
        type=int,
        default=os.environ.get(
            "TROPHOSOME_COMMUNITY_GZIP_LEVEL", str(DEFAULT_GZIP_LEVEL)
        ),
        help=(
            "PRC gzip level from 0 (none) to 9 (smallest); default: %(default)s "
            "(or TROPHOSOME_COMMUNITY_GZIP_LEVEL)"
        ),
    )
    args = parser.parse_args()
    if args.workers < 1:
        parser.error("workers must be at least one")
    if args.progress_every < 1:
        parser.error("progress-every must be at least one")
    if not 0 <= args.gzip_level <= 9:
        parser.error("gzip-level must be between zero and nine")
    mode: CompileMode = (
        "endpoint" if args.endpoint_only else "prc" if args.prc_only else "community"
    )
    compile_inputs(
        args.repository.resolve(),
        args.output,
        through_wave=args.through_wave,
        mode=mode,
        workers=args.workers,
        cache_directory=(
            args.cache_directory.resolve() if args.cache_directory else None
        ),
        use_cache=not args.no_cache,
        rebuild_cache=args.rebuild_cache,
        progress_every=args.progress_every,
        gzip_level=args.gzip_level,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
