"""Streaming runner for the main exact-count Wright--Fisher prototype."""

from __future__ import annotations

import csv
import hashlib
import json
import os
import platform
import subprocess
import sys
import tempfile
import time
from collections import deque
from collections.abc import Iterable, Iterator
from concurrent.futures import Executor, ProcessPoolExecutor, ThreadPoolExecutor
from dataclasses import asdict, dataclass
from datetime import UTC, datetime
from pathlib import Path
from typing import TextIO

import numpy as np

from trophosome import (
    MODEL_FAMILY,
    MODEL_SPEC_VERSION,
    OUTPUT_SCHEMA_VERSION,
    __version__,
)
from trophosome.checkpointing import (
    CheckpointError,
    CheckpointIntegrityError,
    CheckpointMismatchError,
    RecoveryCheckpoint,
    checkpoint_candidates,
    load_recovery_checkpoint,
    remove_recovery_checkpoints,
    write_recovery_checkpoint,
)
from trophosome.config import (
    EvolutionConfig,
    HostConfig,
    ModelConfig,
    parse_duration_seconds,
)
from trophosome.count_model import (
    IdAllocator,
    LineageEvent,
    LineageRecorder,
    PopulationState,
    free_living_selection_step,
    merge_populations,
    population_size_schedule,
    proportional_rescale,
    sample_population,
    simulate_within_host,
    subtract_population,
)
from trophosome.metrics import population_metrics

SEED_COORDINATE_SCHEME_VERSION = "1.0.0"


@dataclass(frozen=True)
class HostGenerationSummary:
    replicate: int
    host_generation: int
    environment_capacity: int
    environment_size: int
    environment_richness: int
    release_pool_size: int
    post_mix_size: int
    realized_host_feedback: float
    mean_adult_richness: float
    mean_adult_gene_diversity: float
    environment_gene_diversity: float
    environment_mean_within_host_fitness: float
    environment_mean_free_living_fitness: float


@dataclass(frozen=True)
class _HostTask:
    replicate: int
    host_generation: int
    host_id: int
    founders: PopulationState
    founder_depths: np.ndarray
    host: HostConfig
    evolution: EvolutionConfig
    size_schedule: tuple[int, ...]
    escape_count: int
    id_start: int
    id_stop: int
    master_seed: int


@dataclass(frozen=True)
class _HostResult:
    host_id: int
    adult: PopulationState
    escapees: PopulationState | None
    lineage_events: tuple[LineageEvent, ...]
    adult_richness: int
    adult_gene_diversity: float
    adult_mean_within_host_fitness: float
    adult_mean_free_living_fitness: float


class _PopulationAccumulator:
    """Accumulate many sparse populations once, without repeated full merges."""

    def __init__(self) -> None:
        self.counts: dict[int, int] = {}
        self.within_host_fitness: dict[int, float] = {}
        self.free_living_fitness: dict[int, float] = {}

    def add(self, population: PopulationState) -> None:
        for genotype_id, count, host_value, environment_value in zip(
            population.genotype_ids.tolist(),
            population.counts.tolist(),
            population.within_host_fitness.tolist(),
            population.free_living_fitness.tolist(),
            strict=True,
        ):
            genotype_id = int(genotype_id)
            host_value = float(host_value)
            environment_value = float(environment_value)
            previous_host = self.within_host_fitness.get(genotype_id)
            previous_environment = self.free_living_fitness.get(genotype_id)
            if previous_host is not None and previous_host != host_value:
                raise ValueError(
                    f"conflicting within-host fitness for strain {genotype_id}"
                )
            if (
                previous_environment is not None
                and previous_environment != environment_value
            ):
                raise ValueError(
                    f"conflicting free-living fitness for strain {genotype_id}"
                )
            self.within_host_fitness[genotype_id] = host_value
            self.free_living_fitness[genotype_id] = environment_value
            self.counts[genotype_id] = self.counts.get(genotype_id, 0) + int(count)

    @property
    def size(self) -> int:
        return sum(self.counts.values())

    def to_population(self) -> PopulationState | None:
        if not self.counts:
            return None
        ids = np.asarray(sorted(self.counts), dtype=np.int64)
        return PopulationState._trusted(
            ids,
            np.asarray([self.counts[int(i)] for i in ids], dtype=np.int64),
            np.asarray(
                [self.within_host_fitness[int(i)] for i in ids], dtype=np.float64
            ),
            np.asarray(
                [self.free_living_fitness[int(i)] for i in ids], dtype=np.float64
            ),
        )


def _git_commit(repository: Path) -> str | None:
    result = subprocess.run(
        ["git", "rev-parse", "HEAD"],
        cwd=repository,
        text=True,
        capture_output=True,
        check=False,
    )
    return result.stdout.strip() if result.returncode == 0 else None


def _config_hash(config: ModelConfig) -> str:
    serialized = json.dumps(
        config.to_dict(), sort_keys=True, separators=(",", ":")
    ).encode("utf-8")
    return hashlib.sha256(serialized).hexdigest()


def _restart_config_hash(config: ModelConfig) -> str:
    """Hash state-affecting inputs while allowing execution-only resume changes."""

    values = config.to_dict()
    values.pop("execution", None)
    output = values.get("output")
    if isinstance(output, dict):
        output.pop("checkpoint_interval", None)
        output.pop("checkpoint_keep", None)
    serialized = json.dumps(values, sort_keys=True, separators=(",", ":")).encode(
        "utf-8"
    )
    return hashlib.sha256(serialized).hexdigest()


def _source_hash(repository: Path) -> str:
    """Hash maintained source files so dirty or untracked code is represented."""

    root = repository
    if not (root / "src" / "trophosome").is_dir():
        root = Path(__file__).resolve().parents[2]
    paths = sorted((root / "src" / "trophosome").rglob("*.py"))
    if (root / "pyproject.toml").is_file():
        paths.append(root / "pyproject.toml")
    if not paths:
        raise RuntimeError("cannot locate maintained trophosome source files")
    digest = hashlib.sha256()
    for path in sorted(paths):
        relative = path.relative_to(root).as_posix()
        digest.update(relative.encode("utf-8"))
        digest.update(b"\0")
        digest.update(path.read_bytes())
        digest.update(b"\0")
    return digest.hexdigest()


def _sync_directory(path: Path) -> None:
    flags = os.O_RDONLY
    if hasattr(os, "O_DIRECTORY"):
        flags |= os.O_DIRECTORY
    try:
        descriptor = os.open(path, flags)
    except OSError:
        return
    try:
        os.fsync(descriptor)
    finally:
        os.close(descriptor)


def _atomic_json(path: Path, value: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(
        "w", dir=path.parent, delete=False, encoding="utf-8"
    ) as handle:
        json.dump(value, handle, indent=2, sort_keys=True)
        handle.write("\n")
        handle.flush()
        os.fsync(handle.fileno())
        temporary = Path(handle.name)
    os.replace(temporary, path)
    _sync_directory(path.parent)


def _save_population(path: Path, population: PopulationState) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(
        "wb", suffix=".npz", dir=path.parent, delete=False
    ) as handle:
        np.savez_compressed(
            handle,
            genotype_ids=population.genotype_ids,
            counts=population.counts,
            within_host_fitness=population.within_host_fitness,
            free_living_fitness=population.free_living_fitness,
        )
        handle.flush()
        os.fsync(handle.fileno())
        temporary = Path(handle.name)
    os.replace(temporary, path)
    _sync_directory(path.parent)


def _rng(
    master_seed: int,
    replicate: int,
    host_generation: int,
    host_id: int,
    stage: int,
) -> np.random.Generator:
    """Return an order-independent random stream for one stochastic stage."""

    seed = np.random.SeedSequence(
        [master_seed, replicate, host_generation, host_id, stage]
    )
    return np.random.default_rng(seed)


def _reservoir_founders(
    environment: PopulationState,
    cumulative_probabilities: np.ndarray,
    sample_size: int,
    rng: np.random.Generator,
) -> PopulationState:
    """Draw a small bottleneck without scanning every strain in a multinomial."""

    indices = np.searchsorted(
        cumulative_probabilities, rng.random(sample_size), side="right"
    )
    indices = np.minimum(indices, environment.richness - 1)
    unique, counts = np.unique(indices, return_counts=True)
    return PopulationState._trusted(
        environment.genotype_ids[unique].copy(),
        counts.astype(np.int64),
        environment.within_host_fitness[unique].copy(),
        environment.free_living_fitness[unique].copy(),
    )


def _run_host(task: _HostTask) -> _HostResult:
    recorder = LineageRecorder.from_founders(
        task.founders.genotype_ids, task.founder_depths
    )
    adult, _ = simulate_within_host(
        task.founders,
        task.host,
        task.evolution,
        _rng(
            task.master_seed,
            task.replicate,
            task.host_generation,
            task.host_id,
            1,
        ),
        IdAllocator(task.id_start, task.id_stop),
        record_summaries=False,
        size_schedule=task.size_schedule,
        lineage_recorder=recorder,
        free_living_fitness_rng=_rng(
            task.master_seed,
            task.replicate,
            task.host_generation,
            task.host_id,
            4,
        ),
    )
    escapees = None
    if task.escape_count:
        escapees = sample_population(
            adult,
            task.escape_count,
            _rng(
                task.master_seed,
                task.replicate,
                task.host_generation,
                task.host_id,
                2,
            ),
            replace=False,
        )
    metrics = population_metrics(adult)
    return _HostResult(
        host_id=task.host_id,
        adult=adult,
        escapees=escapees,
        lineage_events=tuple(recorder.events),
        adult_richness=adult.richness,
        adult_gene_diversity=metrics.gene_diversity,
        adult_mean_within_host_fitness=metrics.mean_within_host_fitness,
        adult_mean_free_living_fitness=metrics.mean_free_living_fitness,
    )


def _run_host_batch(tasks: tuple[_HostTask, ...]) -> tuple[_HostResult, ...]:
    return tuple(_run_host(task) for task in tasks)


def _batches(values: Iterable[_HostTask], size: int) -> Iterator[tuple[_HostTask, ...]]:
    batch: list[_HostTask] = []
    for value in values:
        batch.append(value)
        if len(batch) == size:
            yield tuple(batch)
            batch.clear()
    if batch:
        yield tuple(batch)


def _bounded_host_map(
    executor: Executor,
    task_batches: Iterable[tuple[_HostTask, ...]],
    max_in_flight: int,
) -> Iterator[tuple[_HostResult, ...]]:
    """Run batches with bounded submission while retaining deterministic order.

    ``Executor.map`` eagerly consumes its input on supported Python versions.
    For very large host populations that quietly retains one future and its
    serialized arguments per batch. This small FIFO window keeps memory tied to
    worker count while producing results in host order.
    """

    if max_in_flight < 1:
        raise ValueError("max_in_flight must be positive")
    pending = deque()
    for batch in task_batches:
        pending.append(executor.submit(_run_host_batch, batch))
        if len(pending) >= max_in_flight:
            yield pending.popleft().result()
    while pending:
        yield pending.popleft().result()


@dataclass
class _CsvOutput:
    path: Path
    handle: TextIO
    writer: csv.DictWriter


def _output_fields(config: ModelConfig) -> dict[str, list[str]]:
    fields = {
        "environment_counts.csv": [
            "replicate",
            "generation",
            "strain_id",
            "count",
            "frequency",
            "within_host_fitness",
            "free_living_fitness",
        ],
        "infection_counts.csv": [
            "replicate",
            "generation",
            "host_id",
            "strain_id",
            "count",
        ],
        "host_adult_summaries.csv": [
            "replicate",
            "generation",
            "host_id",
            "richness",
            "gene_diversity",
            "mean_within_host_fitness",
            "mean_free_living_fitness",
        ],
        "pooled_host_counts_and_occupancy.csv": [
            "replicate",
            "generation",
            "strain_id",
            "pooled_count",
            "occupied_hosts",
        ],
        "release_counts.csv": [
            "replicate",
            "generation",
            "host_id",
            "strain_id",
            "count",
        ],
        "strain_lineage_events.csv": [
            "replicate",
            "generation",
            "host_id",
            "within_host_generation",
            "strain_id",
            "parent_strain_id",
            "mutational_depth",
            "within_host_fitness",
            "free_living_fitness",
        ],
        "host_generation_summary.csv": list(
            HostGenerationSummary.__dataclass_fields__
        ),
    }
    if config.output.host_counts_mode != "summary":
        fields["host_adult_counts.csv"] = [
            "replicate",
            "generation",
            "host_id",
            "strain_id",
            "count",
        ]
    return fields


def _open_csv_outputs(
    output: Path, fields_by_name: dict[str, list[str]], append: bool
) -> dict[str, _CsvOutput]:
    opened: dict[str, _CsvOutput] = {}
    try:
        for name, fields in fields_by_name.items():
            path = output / name
            if append:
                with path.open(newline="", encoding="utf-8") as existing:
                    header = next(csv.reader(existing), None)
                if header != fields:
                    raise CheckpointMismatchError(
                        f"output header does not match the checkpoint: {name}"
                    )
            handle = path.open("a" if append else "w", newline="", encoding="utf-8")
            writer = csv.DictWriter(handle, fieldnames=fields)
            if not append:
                writer.writeheader()
            opened[name] = _CsvOutput(path=path, handle=handle, writer=writer)
    except Exception:
        for item in opened.values():
            item.handle.close()
        raise
    return opened


def _sync_csv_outputs(outputs: dict[str, _CsvOutput]) -> dict[str, int]:
    offsets: dict[str, int] = {}
    for name, item in outputs.items():
        item.handle.flush()
        os.fsync(item.handle.fileno())
        offsets[name] = os.fstat(item.handle.fileno()).st_size
    return offsets


def _truncate_outputs(
    output: Path,
    fields_by_name: dict[str, list[str]],
    offsets: dict[str, int],
) -> None:
    if set(offsets) != set(fields_by_name):
        raise CheckpointMismatchError(
            "checkpoint output table set does not match this configuration"
        )
    for name, offset in offsets.items():
        if not isinstance(offset, int) or offset < 0:
            raise CheckpointIntegrityError("checkpoint contains an invalid file offset")
        path = output / name
        if not path.is_file() or path.stat().st_size < offset:
            raise CheckpointIntegrityError(
                f"output file is missing or shorter than its safe offset: {name}"
            )
        with path.open("r+b") as handle:
            handle.truncate(offset)
            handle.flush()
            os.fsync(handle.fileno())
    _sync_directory(output)


def _summary_from_dict(values: dict[str, object]) -> HostGenerationSummary:
    integer_fields = {
        "replicate",
        "host_generation",
        "environment_capacity",
        "environment_size",
        "environment_richness",
        "release_pool_size",
        "post_mix_size",
    }
    converted = {
        name: int(value) if name in integer_fields else float(value)
        for name, value in values.items()
    }
    return HostGenerationSummary(**converted)


def _load_summaries(path: Path) -> list[HostGenerationSummary]:
    with path.open(newline="", encoding="utf-8") as handle:
        return [_summary_from_dict(row) for row in csv.DictReader(handle)]


def _validate_checkpoint_for_run(
    checkpoint: RecoveryCheckpoint,
    config: ModelConfig,
    restart_hash: str,
    source_hash: str,
    environment_capacity: int,
    id_block_size: int,
    root_id_ceiling: int,
    fields_by_name: dict[str, list[str]],
    output: Path,
) -> None:
    metadata = checkpoint.metadata
    expected = {
        "restart_config_sha256": restart_hash,
        "source_sha256": source_hash,
        "model_family": MODEL_FAMILY,
        "model_spec_version": MODEL_SPEC_VERSION,
        "output_schema_version": OUTPUT_SCHEMA_VERSION,
        "software_version": __version__,
        "seed_coordinate_scheme_version": SEED_COORDINATE_SCHEME_VERSION,
        "master_seed": config.seed,
        "environment_capacity": environment_capacity,
        "id_block_size": id_block_size,
        "root_id_ceiling": root_id_ceiling,
    }
    for name, value in expected.items():
        if metadata.get(name) != value:
            raise CheckpointMismatchError(
                f"checkpoint {name} does not match the requested run"
            )
    replicate = metadata.get("replicate")
    generation = metadata.get("last_completed_generation")
    if not isinstance(replicate, int) or not 0 <= replicate < config.replicates:
        raise CheckpointIntegrityError("checkpoint replicate index is invalid")
    if not isinstance(generation, int) or not 1 <= generation <= (
        config.host.host_generations
    ):
        raise CheckpointIntegrityError("checkpoint host generation is invalid")
    expected_summary_count = replicate * config.host.host_generations + generation
    summaries = metadata.get("summaries")
    if not isinstance(summaries, list) or len(summaries) != expected_summary_count:
        raise CheckpointIntegrityError("checkpoint summary history is incomplete")
    if checkpoint.environment.size != environment_capacity:
        raise CheckpointIntegrityError("checkpoint environment has the wrong capacity")
    offsets = metadata.get("output_offsets")
    if not isinstance(offsets, dict):
        raise CheckpointIntegrityError("checkpoint has no committed output offsets")
    _truncate_outputs(output, fields_by_name, offsets)


def _find_resume_checkpoint(
    output: Path,
    config: ModelConfig,
    restart_hash: str,
    source_hash: str,
    environment_capacity: int,
    id_block_size: int,
    root_id_ceiling: int,
    fields_by_name: dict[str, list[str]],
) -> RecoveryCheckpoint:
    candidates = checkpoint_candidates(output / "checkpoints")
    if not candidates:
        raise CheckpointError("no recovery checkpoint exists in the output directory")
    integrity_failures: list[str] = []
    for path in candidates:
        try:
            checkpoint = load_recovery_checkpoint(path)
            _validate_checkpoint_for_run(
                checkpoint,
                config,
                restart_hash,
                source_hash,
                environment_capacity,
                id_block_size,
                root_id_ceiling,
                fields_by_name,
                output,
            )
            return checkpoint
        except CheckpointMismatchError:
            raise
        except CheckpointIntegrityError as exc:
            integrity_failures.append(f"{path.name}: {exc}")
    details = "; ".join(integrity_failures)
    raise CheckpointError(f"no valid recovery checkpoint was found ({details})")


def _checkpoint_metadata(
    config: ModelConfig,
    repository: Path,
    replicate: int,
    generation: int,
    environment_capacity: int,
    id_block_size: int,
    root_id_ceiling: int,
    all_summaries: list[HostGenerationSummary],
    output_offsets: dict[str, int],
    config_hash: str,
    restart_hash: str,
    source_hash: str,
) -> dict[str, object]:
    at_replicate_end = generation == config.host.host_generations
    return {
        "replicate": replicate,
        "last_completed_generation": generation,
        "next_replicate": replicate + 1 if at_replicate_end else replicate,
        "next_generation": 1 if at_replicate_end else generation + 1,
        "environment_capacity": environment_capacity,
        "master_seed": config.seed,
        "seed_coordinate_scheme_version": SEED_COORDINATE_SCHEME_VERSION,
        "config_sha256": config_hash,
        "restart_config_sha256": restart_hash,
        "resolved_config": config.to_dict(),
        "git_commit": _git_commit(repository),
        "source_sha256": source_hash,
        "software_version": __version__,
        "model_family": MODEL_FAMILY,
        "model_spec_version": MODEL_SPEC_VERSION,
        "output_schema_version": OUTPUT_SCHEMA_VERSION,
        "id_block_size": id_block_size,
        "root_id_ceiling": root_id_ceiling,
        "summaries": [asdict(summary) for summary in all_summaries],
        "output_offsets": output_offsets,
        "created_at": datetime.now(UTC).isoformat(),
    }


def _sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _verify_final_environment(path: Path, expected_size: int) -> str:
    try:
        with np.load(path, allow_pickle=False) as payload:
            state = PopulationState(
                np.asarray(payload["genotype_ids"], dtype=np.int64),
                np.asarray(payload["counts"], dtype=np.int64),
                np.asarray(payload["within_host_fitness"], dtype=np.float64),
                np.asarray(payload["free_living_fitness"], dtype=np.float64),
            )
    except (OSError, ValueError, KeyError) as exc:
        raise RuntimeError(f"invalid final environmental state {path.name}") from exc
    if state.size != expected_size:
        raise RuntimeError(f"final environmental state has wrong size: {path.name}")
    return _sha256_file(path)


def _verify_completed_run(
    output: Path,
    config: ModelConfig,
    fields_by_name: dict[str, list[str]],
    restart_hash: str,
    source_hash: str,
) -> tuple[HostGenerationSummary, ...]:
    completion_path = output / "completion.json"
    try:
        completion = json.loads(completion_path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise CheckpointError("completion record is unreadable") from exc
    if completion.get("complete") is not True:
        raise CheckpointError("completion record is not committed")
    expected = {
        "restart_config_sha256": restart_hash,
        "source_sha256": source_hash,
        "model_spec_version": MODEL_SPEC_VERSION,
        "output_schema_version": OUTPUT_SCHEMA_VERSION,
        "software_version": __version__,
    }
    for name, value in expected.items():
        if completion.get(name) != value:
            raise CheckpointMismatchError(
                f"completed run {name} does not match the requested run"
            )
    recorded_sizes = completion.get("output_sizes")
    if not isinstance(recorded_sizes, dict) or set(recorded_sizes) != set(
        fields_by_name
    ):
        raise CheckpointError("completion record has an invalid output table set")
    for name, size in recorded_sizes.items():
        path = output / name
        if not path.is_file() or path.stat().st_size != size:
            raise CheckpointError(f"completed output size does not match: {name}")
    summaries = _load_summaries(output / "host_generation_summary.csv")
    expected_count = config.replicates * config.host.host_generations
    if len(summaries) != expected_count:
        raise CheckpointError("completed generation summary is incomplete")
    environment_capacity = round(
        config.environment.capacity_ratio * config.host.carrying_capacity
    )
    final_hashes = completion.get("final_environment_sha256")
    if not isinstance(final_hashes, dict):
        raise CheckpointError("completion record has no final environment checksums")
    for replicate in range(config.replicates):
        name = f"final_environment_rep{replicate:03d}.npz"
        observed = _verify_final_environment(output / name, environment_capacity)
        if final_hashes.get(name) != observed:
            raise CheckpointError(f"final environment checksum does not match: {name}")
    remove_recovery_checkpoints(output / "checkpoints")
    return tuple(summaries)


def _record_resume(
    output: Path,
    checkpoint: RecoveryCheckpoint,
    config: ModelConfig,
) -> None:
    path = output / "provenance.json"
    try:
        provenance = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise CheckpointError("cannot read run provenance during resume") from exc
    resumes = provenance.setdefault("resumes", [])
    if not isinstance(resumes, list):
        raise CheckpointError("run provenance has an invalid resume history")
    resumes.append(
        {
            "checkpoint": checkpoint.path.name,
            "resumed_at": datetime.now(UTC).isoformat(),
            "execution": asdict(config.execution),
            "config_sha256": _config_hash(config),
        }
    )
    _atomic_json(path, provenance)


def run_simulation(
    config: ModelConfig,
    output_directory: str | Path,
    repository: str | Path,
    *,
    resume: bool = False,
) -> tuple[HostGenerationSummary, ...]:
    """Run or safely resume exact-count replicates at host-generation boundaries."""

    config.validate()
    output = Path(output_directory)
    output.mkdir(parents=True, exist_ok=True)
    repository_path = Path(repository)
    fields_by_name = _output_fields(config)
    config_hash = _config_hash(config)
    restart_hash = _restart_config_hash(config)
    source_hash = _source_hash(repository_path)
    environment_capacity = round(
        config.environment.capacity_ratio * config.host.carrying_capacity
    )
    schedule = population_size_schedule(config.host)
    transitions = max(1, len(schedule))
    id_block_size = transitions * config.evolution.max_materialized_mutants
    root_id_ceiling = len(config.environment.initial_counts)
    total_host_lifetimes = (
        config.replicates * config.host.host_generations * config.host.population_size
    )
    if root_id_ceiling + total_host_lifetimes * id_block_size > np.iinfo(np.int64).max:
        raise ValueError("deterministic strain-ID reservations exceed int64 capacity")

    completion_path = output / "completion.json"
    if completion_path.exists():
        if not resume:
            raise FileExistsError(
                "the output directory already contains a completed run; choose a "
                "new directory"
            )
        return _verify_completed_run(
            output, config, fields_by_name, restart_hash, source_hash
        )

    checkpoint: RecoveryCheckpoint | None = None
    if resume:
        checkpoint = _find_resume_checkpoint(
            output,
            config,
            restart_hash,
            source_hash,
            environment_capacity,
            id_block_size,
            root_id_ceiling,
            fields_by_name,
        )
        all_summaries = [
            _summary_from_dict(item) for item in checkpoint.metadata["summaries"]
        ]
        _record_resume(output, checkpoint, config)
        csv_outputs = _open_csv_outputs(output, fields_by_name, append=True)
    else:
        run_artifacts = [
            output / "resolved_config.json",
            output / "provenance.json",
            output / "checkpoints",
            *[output / name for name in fields_by_name],
            *output.glob("final_environment_rep*.npz"),
        ]
        if any(path.exists() for path in run_artifacts):
            raise FileExistsError(
                "the output directory already contains run files; use --resume "
                "or choose a new directory"
            )
        _atomic_json(output / "resolved_config.json", config.to_dict())
        _atomic_json(
            output / "provenance.json",
            {
                "created_at": datetime.now(UTC).isoformat(),
                "config_sha256": config_hash,
                "restart_config_sha256": restart_hash,
                "git_commit": _git_commit(repository_path),
                "source_sha256": source_hash,
                "model": config.model,
                "model_family": MODEL_FAMILY,
                "model_spec_version": MODEL_SPEC_VERSION,
                "numpy_version": np.__version__,
                "output_schema_version": OUTPUT_SCHEMA_VERSION,
                "platform": platform.platform(),
                "python": sys.version,
                "seed": config.seed,
                "seed_coordinate_scheme_version": (
                    SEED_COORDINATE_SCHEME_VERSION
                ),
                "software_version": __version__,
                "resumes": [],
            },
        )
        all_summaries: list[HostGenerationSummary] = []
        csv_outputs = _open_csv_outputs(output, fields_by_name, append=False)

    environment_writer = csv_outputs["environment_counts.csv"].writer
    infection_writer = csv_outputs["infection_counts.csv"].writer
    adult_summary_writer = csv_outputs["host_adult_summaries.csv"].writer
    pooled_writer = csv_outputs["pooled_host_counts_and_occupancy.csv"].writer
    release_writer = csv_outputs["release_counts.csv"].writer
    lineage_writer = csv_outputs["strain_lineage_events.csv"].writer
    summary_writer = csv_outputs["host_generation_summary.csv"].writer
    adult_counts_output = csv_outputs.get("host_adult_counts.csv")
    adult_counts_writer = (
        None if adult_counts_output is None else adult_counts_output.writer
    )

    checkpoint_interval = parse_duration_seconds(config.output.checkpoint_interval)
    last_checkpoint_time = time.monotonic()
    executor: Executor | None = None
    if config.execution.workers > 1:
        try:
            executor = ProcessPoolExecutor(max_workers=config.execution.workers)
        except (OSError, PermissionError):
            executor = ThreadPoolExecutor(max_workers=config.execution.workers)

    if checkpoint is None:
        start_replicate = 0
        resume_generation = 0
    else:
        checkpoint_replicate = int(checkpoint.metadata["replicate"])
        checkpoint_generation = int(
            checkpoint.metadata["last_completed_generation"]
        )
        if checkpoint_generation == config.host.host_generations:
            start_replicate = checkpoint_replicate + 1
            resume_generation = 0
        else:
            start_replicate = checkpoint_replicate
            resume_generation = checkpoint_generation

    try:
        for replicate in range(start_replicate, config.replicates):
            if (
                checkpoint is not None
                and replicate == start_replicate
                and resume_generation > 0
            ):
                environment = checkpoint.environment
                active_depths = {
                    int(genotype_id): int(depth)
                    for genotype_id, depth in zip(
                        environment.genotype_ids.tolist(),
                        checkpoint.active_depths.tolist(),
                        strict=True,
                    )
                }
                first_generation = resume_generation + 1
            else:
                environment = proportional_rescale(
                    PopulationState.from_counts(
                        config.environment.initial_counts,
                        config.environment.initial_within_host_fitness,
                        config.environment.initial_free_living_fitness,
                    ),
                    environment_capacity,
                    _rng(config.seed, replicate, 0, 0, 3),
                )
                active_depths = {
                    int(genotype_id): 0
                    for genotype_id in environment.genotype_ids.tolist()
                }
                first_generation = 1

            def write_environment(generation: int) -> None:
                for genotype_id, count, host_value, environment_value in zip(
                    environment.genotype_ids.tolist(),
                    environment.counts.tolist(),
                    environment.within_host_fitness.tolist(),
                    environment.free_living_fitness.tolist(),
                    strict=True,
                ):
                    environment_writer.writerow(
                        {
                            "replicate": replicate,
                            "generation": generation,
                            "strain_id": genotype_id,
                            "count": count,
                            "frequency": count / environment.size,
                            "within_host_fitness": host_value,
                            "free_living_fitness": environment_value,
                        }
                    )

            if (
                first_generation == 1
                and config.output.environment_counts_mode == "all"
            ):
                write_environment(0)
            for host_generation in range(
                first_generation, config.host.host_generations + 1
            ):
                source = environment
                remaining: PopulationState | None = environment
                cumulative = np.cumsum(source.counts, dtype=np.float64)
                cumulative /= cumulative[-1]
                escape_count = int(
                    round(config.host.carrying_capacity * config.host.escape_fraction)
                )
                if config.output.host_counts_mode == "full":
                    retained_hosts = set(range(config.host.population_size))
                elif config.output.host_counts_mode == "panel":
                    panel_size = min(
                        config.output.host_panel_size, config.host.population_size
                    )
                    panel_rng = _rng(config.seed, replicate, host_generation, 0, 9)
                    retained_hosts = set(
                        panel_rng.choice(
                            config.host.population_size,
                            size=panel_size,
                            replace=False,
                        ).tolist()
                    )
                else:
                    retained_hosts = set()

                def tasks() -> Iterator[_HostTask]:
                    nonlocal remaining
                    for host_id in range(config.host.population_size):
                        if config.environment.sampling_mode == "finite":
                            if remaining is None:
                                raise RuntimeError(
                                    "finite environmental population was depleted"
                                )
                            founders = sample_population(
                                remaining,
                                config.host.infection_bottleneck,
                                _rng(
                                    config.seed,
                                    replicate,
                                    host_generation,
                                    host_id,
                                    0,
                                ),
                                replace=False,
                            )
                            remaining = subtract_population(remaining, founders)
                        else:
                            founders = _reservoir_founders(
                                source,
                                cumulative,
                                config.host.infection_bottleneck,
                                _rng(
                                    config.seed,
                                    replicate,
                                    host_generation,
                                    host_id,
                                    0,
                                ),
                            )
                        for genotype_id, count in zip(
                            founders.genotype_ids.tolist(),
                            founders.counts.tolist(),
                            strict=True,
                        ):
                            infection_writer.writerow(
                                {
                                    "replicate": replicate,
                                    "generation": host_generation,
                                    "host_id": host_id,
                                    "strain_id": genotype_id,
                                    "count": count,
                                }
                            )
                        ordinal = (
                            replicate
                            * config.host.host_generations
                            * config.host.population_size
                            + (host_generation - 1) * config.host.population_size
                            + host_id
                        )
                        id_start = root_id_ceiling + ordinal * id_block_size
                        yield _HostTask(
                            replicate=replicate,
                            host_generation=host_generation,
                            host_id=host_id,
                            founders=founders,
                            founder_depths=np.asarray(
                                [
                                    active_depths[int(genotype_id)]
                                    for genotype_id in founders.genotype_ids
                                ],
                                dtype=np.int64,
                            ),
                            host=config.host,
                            evolution=config.evolution,
                            size_schedule=schedule,
                            escape_count=escape_count,
                            id_start=id_start,
                            id_stop=id_start + id_block_size,
                            master_seed=config.seed,
                        )

                task_batches = _batches(tasks(), config.execution.host_batch_size)
                result_batches: Iterable[tuple[_HostResult, ...]]
                if executor is None:
                    result_batches = map(_run_host_batch, task_batches)
                else:
                    result_batches = _bounded_host_map(
                        executor,
                        task_batches,
                        max_in_flight=(
                            config.execution.workers
                            * config.execution.in_flight_batches_per_worker
                        ),
                    )

                release_accumulator = _PopulationAccumulator()
                adult_accumulator = _PopulationAccumulator()
                occupancy: dict[int, int] = {}
                adult_richness = 0.0
                adult_gene_diversity = 0.0
                for result_batch in result_batches:
                    for result in result_batch:
                        adult_richness += result.adult_richness
                        adult_gene_diversity += result.adult_gene_diversity
                        adult_accumulator.add(result.adult)
                        for genotype_id in result.adult.genotype_ids.tolist():
                            genotype_id = int(genotype_id)
                            occupancy[genotype_id] = occupancy.get(genotype_id, 0) + 1
                        adult_summary_writer.writerow(
                            {
                                "replicate": replicate,
                                "generation": host_generation,
                                "host_id": result.host_id,
                                "richness": result.adult_richness,
                                "gene_diversity": result.adult_gene_diversity,
                                "mean_within_host_fitness": (
                                    result.adult_mean_within_host_fitness
                                ),
                                "mean_free_living_fitness": (
                                    result.adult_mean_free_living_fitness
                                ),
                            }
                        )
                        if (
                            adult_counts_writer is not None
                            and result.host_id in retained_hosts
                        ):
                            for genotype_id, count in zip(
                                result.adult.genotype_ids.tolist(),
                                result.adult.counts.tolist(),
                                strict=True,
                            ):
                                adult_counts_writer.writerow(
                                    {
                                        "replicate": replicate,
                                        "generation": host_generation,
                                        "host_id": result.host_id,
                                        "strain_id": genotype_id,
                                        "count": count,
                                    }
                                )
                        for event in result.lineage_events:
                            active_depths[event.strain_id] = event.mutational_depth
                            lineage_writer.writerow(
                                {
                                    "replicate": replicate,
                                    "generation": host_generation,
                                    "host_id": result.host_id,
                                    "within_host_generation": (
                                        event.within_host_generation
                                    ),
                                    "strain_id": event.strain_id,
                                    "parent_strain_id": event.parent_strain_id,
                                    "mutational_depth": event.mutational_depth,
                                    "within_host_fitness": event.within_host_fitness,
                                    "free_living_fitness": event.free_living_fitness,
                                }
                            )
                        if result.escapees is not None:
                            release_accumulator.add(result.escapees)
                            for genotype_id, count in zip(
                                result.escapees.genotype_ids.tolist(),
                                result.escapees.counts.tolist(),
                                strict=True,
                            ):
                                release_writer.writerow(
                                    {
                                        "replicate": replicate,
                                        "generation": host_generation,
                                        "host_id": result.host_id,
                                        "strain_id": genotype_id,
                                        "count": count,
                                    }
                                )

                pooled_adults = adult_accumulator.to_population()
                assert pooled_adults is not None
                for genotype_id, count in zip(
                    pooled_adults.genotype_ids.tolist(),
                    pooled_adults.counts.tolist(),
                    strict=True,
                ):
                    pooled_writer.writerow(
                        {
                            "replicate": replicate,
                            "generation": host_generation,
                            "strain_id": genotype_id,
                            "pooled_count": count,
                            "occupied_hosts": occupancy[int(genotype_id)],
                        }
                    )

                released = release_accumulator.to_population()
                base = (
                    remaining
                    if config.environment.sampling_mode == "finite"
                    else source
                )
                if base is None and released is None:
                    raise RuntimeError(
                        "environment went extinct: all cells infected hosts and "
                        "none returned"
                    )
                if base is None:
                    assert released is not None
                    post_mix = released
                elif released is None:
                    post_mix = base
                else:
                    post_mix = merge_populations(base, released)
                release_pool_size = 0 if released is None else released.size
                realized_feedback = release_pool_size / post_mix.size
                if config.evolution.free_living_selection:
                    environment = free_living_selection_step(
                        post_mix,
                        environment_capacity,
                        _rng(config.seed, replicate, host_generation, 0, 5),
                    )
                else:
                    environment = proportional_rescale(
                        post_mix,
                        environment_capacity,
                        _rng(config.seed, replicate, host_generation, 0, 3),
                    )
                active_depths = {
                    int(genotype_id): active_depths[int(genotype_id)]
                    for genotype_id in environment.genotype_ids.tolist()
                }
                environment_metrics = population_metrics(environment)
                summary = HostGenerationSummary(
                    replicate=replicate,
                    host_generation=host_generation,
                    environment_capacity=environment_capacity,
                    environment_size=environment.size,
                    environment_richness=environment.richness,
                    release_pool_size=release_pool_size,
                    post_mix_size=post_mix.size,
                    realized_host_feedback=realized_feedback,
                    mean_adult_richness=(adult_richness / config.host.population_size),
                    mean_adult_gene_diversity=(
                        adult_gene_diversity / config.host.population_size
                    ),
                    environment_gene_diversity=environment_metrics.gene_diversity,
                    environment_mean_within_host_fitness=(
                        environment_metrics.mean_within_host_fitness
                    ),
                    environment_mean_free_living_fitness=(
                        environment_metrics.mean_free_living_fitness
                    ),
                )
                all_summaries.append(summary)
                summary_writer.writerow(asdict(summary))
                if (
                    config.output.environment_counts_mode == "all"
                    or host_generation == config.host.host_generations
                ):
                    write_environment(host_generation)

                at_replicate_end = (
                    host_generation == config.host.host_generations
                )
                if at_replicate_end:
                    _save_population(
                        output / f"final_environment_rep{replicate:03d}.npz",
                        environment,
                    )
                checkpoint_due = (
                    time.monotonic() - last_checkpoint_time >= checkpoint_interval
                )
                if checkpoint_due or at_replicate_end:
                    offsets = _sync_csv_outputs(csv_outputs)
                    depths = np.asarray(
                        [
                            active_depths[int(genotype_id)]
                            for genotype_id in environment.genotype_ids
                        ],
                        dtype=np.int64,
                    )
                    metadata = _checkpoint_metadata(
                        config,
                        repository_path,
                        replicate,
                        host_generation,
                        environment_capacity,
                        id_block_size,
                        root_id_ceiling,
                        all_summaries,
                        offsets,
                        config_hash,
                        restart_hash,
                        source_hash,
                    )
                    write_recovery_checkpoint(
                        output / "checkpoints",
                        metadata,
                        environment,
                        depths,
                        config.output.checkpoint_keep,
                    )
                    last_checkpoint_time = time.monotonic()

        output_sizes = _sync_csv_outputs(csv_outputs)
        expected_summary_count = config.replicates * config.host.host_generations
        if len(all_summaries) != expected_summary_count:
            raise RuntimeError("generation summary is incomplete at finalization")
        final_hashes = {}
        for replicate in range(config.replicates):
            name = f"final_environment_rep{replicate:03d}.npz"
            final_hashes[name] = _verify_final_environment(
                output / name, environment_capacity
            )
        _atomic_json(
            completion_path,
            {
                "complete": True,
                "completed_at": datetime.now(UTC).isoformat(),
                "config_sha256": config_hash,
                "restart_config_sha256": restart_hash,
                "source_sha256": source_hash,
                "software_version": __version__,
                "model_spec_version": MODEL_SPEC_VERSION,
                "output_schema_version": OUTPUT_SCHEMA_VERSION,
                "output_sizes": output_sizes,
                "final_environment_sha256": final_hashes,
            },
        )
        remove_recovery_checkpoints(output / "checkpoints")
        return tuple(all_summaries)
    finally:
        if executor is not None:
            executor.shutdown(wait=True, cancel_futures=True)
        for item in csv_outputs.values():
            item.handle.close()
