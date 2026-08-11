"""Validated atomic recovery checkpoints for generation-boundary restart."""

from __future__ import annotations

import hashlib
import json
import os
import re
import uuid
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np

from trophosome.count_model import PopulationState

CHECKPOINT_SCHEMA_VERSION = "1.0.0"
_CHECKPOINT_NAME = re.compile(
    r"^checkpoint-rep(?P<replicate>\d+)-gen(?P<generation>\d+)\.npz$"
)


class CheckpointError(RuntimeError):
    """Base class for checkpoint or restart failures."""


class CheckpointIntegrityError(CheckpointError):
    """Raised when a checkpoint is incomplete, unreadable, or corrupted."""


class CheckpointMismatchError(CheckpointError):
    """Raised when a valid checkpoint belongs to a different run."""


@dataclass(frozen=True)
class RecoveryCheckpoint:
    """One validated recovery point at a completed host-generation boundary."""

    path: Path
    metadata: dict[str, Any]
    environment: PopulationState
    active_depths: np.ndarray


def _canonical_json(value: object) -> bytes:
    return json.dumps(
        value,
        allow_nan=False,
        separators=(",", ":"),
        sort_keys=True,
    ).encode("utf-8")


def _payload_checksum(metadata: dict[str, Any], arrays: dict[str, np.ndarray]) -> str:
    digest = hashlib.sha256()
    digest.update(_canonical_json(metadata))
    for name in sorted(arrays):
        array = np.ascontiguousarray(arrays[name])
        digest.update(name.encode("utf-8"))
        digest.update(array.dtype.str.encode("ascii"))
        digest.update(_canonical_json(array.shape))
        digest.update(array.tobytes())
    return digest.hexdigest()


def _sync_directory(path: Path) -> None:
    """Best-effort directory synchronization after an atomic rename."""

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


def checkpoint_candidates(directory: Path) -> tuple[Path, ...]:
    """Return recovery checkpoint files from newest logical boundary to oldest."""

    if not directory.is_dir():
        return ()
    candidates = []
    for path in directory.glob("checkpoint-rep*-gen*.npz"):
        match = _CHECKPOINT_NAME.fullmatch(path.name)
        if match is not None:
            candidates.append(
                (
                    int(match.group("replicate")),
                    int(match.group("generation")),
                    path,
                )
            )
    return tuple(path for _, _, path in sorted(candidates, reverse=True))


def load_recovery_checkpoint(path: str | Path) -> RecoveryCheckpoint:
    """Open and fully validate one recovery checkpoint."""

    checkpoint_path = Path(path)
    try:
        with np.load(checkpoint_path, allow_pickle=False) as payload:
            required = {
                "metadata",
                "genotype_ids",
                "counts",
                "within_host_fitness",
                "free_living_fitness",
                "active_depths",
            }
            missing = required.difference(payload.files)
            if missing:
                raise CheckpointIntegrityError(
                    f"checkpoint is missing arrays: {', '.join(sorted(missing))}"
                )
            metadata_text = str(payload["metadata"].item())
            metadata = json.loads(metadata_text)
            arrays = {
                "genotype_ids": np.asarray(payload["genotype_ids"], dtype=np.int64),
                "counts": np.asarray(payload["counts"], dtype=np.int64),
                "within_host_fitness": np.asarray(
                    payload["within_host_fitness"], dtype=np.float64
                ),
                "free_living_fitness": np.asarray(
                    payload["free_living_fitness"], dtype=np.float64
                ),
                "active_depths": np.asarray(payload["active_depths"], dtype=np.int64),
            }
    except CheckpointIntegrityError:
        raise
    except (OSError, ValueError, KeyError, json.JSONDecodeError) as exc:
        raise CheckpointIntegrityError(
            f"cannot read recovery checkpoint {checkpoint_path.name}: {exc}"
        ) from exc

    if not isinstance(metadata, dict):
        raise CheckpointIntegrityError("checkpoint metadata must be a JSON object")
    if metadata.get("checkpoint_schema_version") != CHECKPOINT_SCHEMA_VERSION:
        raise CheckpointIntegrityError("unsupported checkpoint schema version")
    if metadata.get("complete") is not True:
        raise CheckpointIntegrityError("checkpoint was not completely committed")
    recorded_checksum = metadata.get("payload_sha256")
    if not isinstance(recorded_checksum, str):
        raise CheckpointIntegrityError("checkpoint has no payload checksum")
    checksum_metadata = dict(metadata)
    del checksum_metadata["payload_sha256"]
    observed_checksum = _payload_checksum(checksum_metadata, arrays)
    if observed_checksum != recorded_checksum:
        raise CheckpointIntegrityError("checkpoint payload checksum does not match")

    try:
        environment = PopulationState(
            arrays["genotype_ids"],
            arrays["counts"],
            arrays["within_host_fitness"],
            arrays["free_living_fitness"],
        )
    except ValueError as exc:
        raise CheckpointIntegrityError(f"invalid environmental state: {exc}") from exc
    active_depths = arrays["active_depths"]
    if len(active_depths) != environment.richness or np.any(active_depths < 0):
        raise CheckpointIntegrityError(
            "active mutational depths must be non-negative and align with strains"
        )
    return RecoveryCheckpoint(
        path=checkpoint_path,
        metadata=metadata,
        environment=environment,
        active_depths=active_depths,
    )


def write_recovery_checkpoint(
    directory: str | Path,
    metadata: dict[str, Any],
    environment: PopulationState,
    active_depths: np.ndarray,
    keep: int,
) -> Path:
    """Atomically write, reopen, validate, and retain the newest checkpoints."""

    if keep < 2:
        raise ValueError("at least two recovery checkpoints must be retained")
    depth_array = np.asarray(active_depths, dtype=np.int64)
    if len(depth_array) != environment.richness or np.any(depth_array < 0):
        raise ValueError("active depths must align with the environmental strains")

    checkpoint_directory = Path(directory)
    checkpoint_directory.mkdir(parents=True, exist_ok=True)
    committed_metadata = dict(metadata)
    committed_metadata.update(
        {
            "checkpoint_schema_version": CHECKPOINT_SCHEMA_VERSION,
            "complete": True,
        }
    )
    arrays = {
        "genotype_ids": environment.genotype_ids,
        "counts": environment.counts,
        "within_host_fitness": environment.within_host_fitness,
        "free_living_fitness": environment.free_living_fitness,
        "active_depths": depth_array,
    }
    committed_metadata["payload_sha256"] = _payload_checksum(
        committed_metadata, arrays
    )

    replicate = int(committed_metadata["replicate"])
    generation = int(committed_metadata["last_completed_generation"])
    destination = checkpoint_directory / (
        f"checkpoint-rep{replicate:03d}-gen{generation:06d}.npz"
    )
    temporary = checkpoint_directory / (
        f".{destination.name}.tmp-{uuid.uuid4().hex}"
    )
    try:
        with temporary.open("wb") as handle:
            np.savez_compressed(
                handle,
                metadata=np.asarray(_canonical_json(committed_metadata).decode()),
                **arrays,
            )
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary, destination)
        _sync_directory(checkpoint_directory)
        load_recovery_checkpoint(destination)
    finally:
        temporary.unlink(missing_ok=True)

    candidates = checkpoint_candidates(checkpoint_directory)
    for obsolete in candidates[keep:]:
        obsolete.unlink(missing_ok=True)
    _sync_directory(checkpoint_directory)
    return destination


def remove_recovery_checkpoints(directory: str | Path) -> None:
    """Remove recovery-only files after a completed run has been verified."""

    checkpoint_directory = Path(directory)
    if not checkpoint_directory.is_dir():
        return
    patterns = ("checkpoint-rep*-gen*.npz", ".checkpoint-*.tmp-*", ".*.npz.tmp-*")
    for pattern in patterns:
        for path in checkpoint_directory.glob(pattern):
            if path.is_file():
                path.unlink()
    try:
        checkpoint_directory.rmdir()
    except OSError:
        _sync_directory(checkpoint_directory)
