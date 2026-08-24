"""Typed, validated configuration for reproducible simulations."""

from __future__ import annotations

import math
import re
import tomllib
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any, Literal


class ConfigurationError(ValueError):
    """Raised when a configuration is internally inconsistent."""


_DURATION_PATTERN = re.compile(
    r"^(?P<value>(?:\d+(?:\.\d*)?|\.\d+))(?P<unit>ms|s|m|h)$"
)
_DURATION_MULTIPLIERS = {"ms": 0.001, "s": 1.0, "m": 60.0, "h": 3600.0}


def parse_duration_seconds(value: str) -> float:
    """Convert a compact wall-clock duration such as ``30m`` or ``1h``."""

    if not isinstance(value, str):
        raise ConfigurationError(
            "output.checkpoint_interval must be a duration string such as "
            "'30m', '1h', or '2h'"
        )
    match = _DURATION_PATTERN.fullmatch(value.strip().lower())
    if match is None:
        raise ConfigurationError(
            "output.checkpoint_interval must be a duration string such as "
            "'30m', '1h', or '2h'"
        )
    seconds = float(match.group("value")) * _DURATION_MULTIPLIERS[match.group("unit")]
    if not math.isfinite(seconds) or seconds <= 0:
        raise ConfigurationError("output.checkpoint_interval must be positive")
    return seconds


@dataclass(frozen=True)
class EnvironmentConfig:
    initial_counts: tuple[int, ...]
    initial_within_host_fitness: tuple[float, ...]
    initial_free_living_fitness: tuple[float, ...]
    sampling_mode: Literal["reservoir", "finite"] = "reservoir"
    capacity_ratio: float = 1.0


@dataclass(frozen=True)
class MigrationConfig:
    """Exchange between the focal environment and a fixed regional source."""

    mode: Literal["none", "fixed_regional_pool"] = "none"
    fraction: float = 0.0
    regional_counts: tuple[int, ...] = ()


@dataclass(frozen=True)
class HostConfig:
    population_size: int
    infection_bottleneck: int
    carrying_capacity: int
    growth_factor: float
    steady_generations: int
    host_generations: int
    escape_fraction: float


@dataclass(frozen=True)
class EvolutionConfig:
    mutation_probability: float
    mutation_effect_mean: float = -0.01
    mutation_effect_sd: float = 0.01
    within_host_selection: bool = True
    free_living_selection: bool = False
    fitness_floor: float = 0.0
    max_materialized_mutants: int = 100_000


@dataclass(frozen=True)
class OutputConfig:
    snapshot_interval: int = 1
    checkpoint_interval: str = "1h"
    checkpoint_keep: int = 2
    retain_host_histories: bool = False
    environment_counts_mode: Literal["all", "final"] = "all"
    host_counts_mode: Literal["summary", "panel", "full"] = "summary"
    host_panel_size: int = 100


@dataclass(frozen=True)
class ExecutionConfig:
    workers: int = 1
    host_batch_size: int = 32
    in_flight_batches_per_worker: int = 2


@dataclass(frozen=True)
class ModelConfig:
    seed: int
    replicates: int
    environment: EnvironmentConfig
    host: HostConfig
    evolution: EvolutionConfig
    migration: MigrationConfig = MigrationConfig()
    output: OutputConfig = OutputConfig()
    execution: ExecutionConfig = ExecutionConfig()
    model: Literal["wright_fisher_counts"] = "wright_fisher_counts"

    def validate(self) -> None:
        errors: list[str] = []
        if self.seed < 0:
            errors.append("seed must be non-negative")
        if self.replicates < 1:
            errors.append("replicates must be at least 1")
        if not self.environment.initial_counts:
            errors.append("environment.initial_counts must not be empty")
        if len(self.environment.initial_counts) != len(
            self.environment.initial_within_host_fitness
        ):
            errors.append(
                "initial_counts and initial_within_host_fitness must have equal length"
            )
        if len(self.environment.initial_counts) != len(
            self.environment.initial_free_living_fitness
        ):
            errors.append(
                "initial_counts and initial_free_living_fitness must have equal length"
            )
        if any(value < 0 for value in self.environment.initial_counts):
            errors.append("all initial counts must be non-negative")
        if not any(value > 0 for value in self.environment.initial_counts):
            errors.append("at least one initial focal count must be positive")
        if any(
            not math.isfinite(value) or value < 0
            for value in self.environment.initial_within_host_fitness
        ):
            errors.append(
                "all initial within-host fitness values must be finite and non-negative"
            )
        if any(
            not math.isfinite(value) or value < 0
            for value in self.environment.initial_free_living_fitness
        ):
            errors.append(
                "all initial free-living fitness values must be finite and non-negative"
            )
        if self.environment.initial_within_host_fitness and not any(
            count > 0 and value > 0
            for count, value in zip(
                self.environment.initial_counts,
                self.environment.initial_within_host_fitness,
                strict=True,
            )
        ):
            errors.append(
                "at least one initially focal genotype must have positive "
                "within-host fitness"
            )
        if self.environment.initial_free_living_fitness and not any(
            count > 0 and value > 0
            for count, value in zip(
                self.environment.initial_counts,
                self.environment.initial_free_living_fitness,
                strict=True,
            )
        ):
            errors.append(
                "at least one initially focal genotype must have positive "
                "free-living fitness"
            )
        if (
            not math.isfinite(self.environment.capacity_ratio)
            or self.environment.capacity_ratio <= 0
        ):
            errors.append("environment.capacity_ratio must be finite and positive")

        migration = self.migration
        if migration.mode not in {"none", "fixed_regional_pool"}:
            errors.append("migration.mode must be 'none' or 'fixed_regional_pool'")
        if not math.isfinite(migration.fraction) or not 0 <= migration.fraction <= 1:
            errors.append("migration.fraction must be finite and between 0 and 1")
        if migration.mode == "none":
            if migration.fraction != 0:
                errors.append("migration.fraction must be 0 when migration is disabled")
            if migration.regional_counts:
                errors.append(
                    "migration.regional_counts must be empty when migration is disabled"
                )
            if any(value == 0 for value in self.environment.initial_counts):
                errors.append(
                    "all initial focal counts must be positive when migration "
                    "is disabled"
                )
        else:
            if len(migration.regional_counts) != len(self.environment.initial_counts):
                errors.append(
                    "migration.regional_counts and environment.initial_counts "
                    "must have equal length"
                )
            if any(value < 0 for value in migration.regional_counts):
                errors.append("all regional counts must be non-negative")
            if not any(value > 0 for value in migration.regional_counts):
                errors.append("at least one regional count must be positive")
            if len(migration.regional_counts) == len(
                self.environment.initial_counts
            ) and any(
                focal == 0 and regional == 0
                for focal, regional in zip(
                    self.environment.initial_counts,
                    migration.regional_counts,
                    strict=True,
                )
            ):
                errors.append(
                    "every declared strain must occur in the focal or regional pool"
                )

        host = self.host
        if host.population_size < 1:
            errors.append("host.population_size must be at least 1")
        if host.infection_bottleneck < 1:
            errors.append("host.infection_bottleneck must be at least 1")
        if host.carrying_capacity < host.infection_bottleneck:
            errors.append("carrying_capacity must be at least the infection bottleneck")
        if not math.isfinite(host.growth_factor) or host.growth_factor < 1:
            errors.append("growth_factor must be finite and at least 1")
        if (
            host.growth_factor == 1
            and host.carrying_capacity > host.infection_bottleneck
        ):
            errors.append(
                "growth_factor=1 cannot grow a bottleneck to carrying capacity"
            )
        if host.steady_generations < 0 or host.host_generations < 1:
            errors.append("steady_generations must be >=0 and host_generations >=1")
        if not 0 <= host.escape_fraction <= 1:
            errors.append("escape_fraction must be between 0 and 1")

        evolution = self.evolution
        if not 0 <= evolution.mutation_probability <= 1:
            errors.append("mutation_probability must be between 0 and 1")
        if evolution.mutation_effect_sd < 0:
            errors.append("mutation_effect_sd must be non-negative")
        if evolution.fitness_floor < 0:
            errors.append("fitness_floor must be non-negative")
        if evolution.max_materialized_mutants < 1:
            errors.append("max_materialized_mutants must be at least 1")
        if self.output.snapshot_interval < 1:
            errors.append("snapshot_interval must be at least 1")
        try:
            parse_duration_seconds(self.output.checkpoint_interval)
        except ConfigurationError as exc:
            errors.append(str(exc))
        if self.output.checkpoint_keep < 2:
            errors.append("checkpoint_keep must be at least 2")
        if self.output.retain_host_histories:
            errors.append(
                "retain_host_histories is intentionally unsupported in the "
                "scalable path; "
                "write summaries or sampled lineages instead"
            )
        if self.output.environment_counts_mode not in {"all", "final"}:
            errors.append("environment_counts_mode must be 'all' or 'final'")
        if self.output.host_counts_mode not in {"summary", "panel", "full"}:
            errors.append("host_counts_mode must be 'summary', 'panel', or 'full'")
        if self.output.host_panel_size < 1:
            errors.append("host_panel_size must be at least 1")
        if self.execution.workers < 1:
            errors.append("execution.workers must be at least 1")
        if self.execution.host_batch_size < 1:
            errors.append("execution.host_batch_size must be at least 1")
        if self.execution.in_flight_batches_per_worker < 1:
            errors.append("execution.in_flight_batches_per_worker must be at least 1")

        capacity = round(self.environment.capacity_ratio * host.carrying_capacity)
        if self.environment.sampling_mode == "finite":
            demand = host.population_size * host.infection_bottleneck
            if capacity < demand:
                errors.append(
                    "finite environmental sampling needs effective capacity >= "
                    "population_size * infection_bottleneck"
                )
        if capacity < 1:
            errors.append("effective environmental capacity must be at least 1")
        if errors:
            raise ConfigurationError("; ".join(errors))

    def feasibility_warnings(self) -> tuple[str, ...]:
        expected = self.host.carrying_capacity * self.evolution.mutation_probability
        warnings: list[str] = []
        if expected > self.evolution.max_materialized_mutants:
            warnings.append(
                "Expected new mutants per adult bacterial generation "
                f"({expected:,.0f}) exceed max_materialized_mutants "
                f"({self.evolution.max_materialized_mutants:,}); use an aggregate, "
                "rare-variant, or sample-first model."
            )
        if self.host.population_size >= 1_000 and self.output.snapshot_interval == 1:
            warnings.append(
                "Per-generation snapshots across thousands of hosts are expensive; "
                "prefer endpoint summaries and sampled-host trajectories."
            )
        growth_steps = 0
        size = self.host.infection_bottleneck
        while size < self.host.carrying_capacity:
            size = min(
                self.host.carrying_capacity,
                max(size + 1, math.ceil(size * self.host.growth_factor)),
            )
            growth_steps += 1
        transitions = (
            self.host.population_size
            * self.host.host_generations
            * (growth_steps + self.host.steady_generations)
            * self.replicates
        )
        if transitions > 10_000_000:
            warnings.append(
                f"The configuration requires about {transitions:,.0f} within-host "
                "population transitions; stage a smaller timing pilot before the "
                "full run or use an endpoint transition approximation."
            )
        return tuple(warnings)

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


def _section(data: dict[str, Any], name: str) -> dict[str, Any]:
    value = data.get(name)
    if not isinstance(value, dict):
        raise ConfigurationError(f"missing [{name}] section")
    return value


def load_config(path: str | Path) -> ModelConfig:
    """Load and validate a TOML configuration."""

    config_path = Path(path)
    with config_path.open("rb") as handle:
        data = tomllib.load(handle)

    environment_data = _section(data, "environment")
    host_data = _section(data, "host")
    evolution_data = _section(data, "evolution")
    migration_data = data.get("migration", {})
    if not isinstance(migration_data, dict):
        raise ConfigurationError("[migration] must be a TOML table")
    output_data = data.get("output", {})
    if not isinstance(output_data, dict):
        raise ConfigurationError("[output] must be a TOML table")
    execution_data = data.get("execution", {})
    if not isinstance(execution_data, dict):
        raise ConfigurationError("[execution] must be a TOML table")

    try:
        model_name = data.get("model", "wright_fisher_counts")
        if model_name == "exact_counts":
            model_name = "wright_fisher_counts"
        initial_counts = tuple(int(x) for x in environment_data["initial_counts"])
        within_host_values = environment_data.get(
            "initial_within_host_fitness", environment_data.get("initial_fitness")
        )
        if within_host_values is None:
            raise KeyError("initial_within_host_fitness")
        free_living_values = environment_data.get("initial_free_living_fitness")
        if free_living_values is None:
            # Backward-compatible Phase 1 configs remain environmentally neutral.
            free_living_values = [1.0] * len(initial_counts)
        config = ModelConfig(
            seed=int(data.get("seed", 666)),
            replicates=int(data.get("replicates", 1)),
            model=model_name,
            environment=EnvironmentConfig(
                initial_counts=initial_counts,
                initial_within_host_fitness=tuple(float(x) for x in within_host_values),
                initial_free_living_fitness=tuple(float(x) for x in free_living_values),
                sampling_mode=environment_data.get("sampling_mode", "reservoir"),
                capacity_ratio=float(environment_data.get("capacity_ratio", 1.0)),
            ),
            host=HostConfig(**host_data),
            evolution=EvolutionConfig(**evolution_data),
            migration=MigrationConfig(
                mode=migration_data.get("mode", "none"),
                fraction=float(migration_data.get("fraction", 0.0)),
                regional_counts=tuple(
                    int(x) for x in migration_data.get("regional_counts", ())
                ),
            ),
            output=OutputConfig(**output_data),
            execution=ExecutionConfig(**execution_data),
        )
    except (KeyError, TypeError, ValueError) as exc:
        raise ConfigurationError(f"invalid configuration: {exc}") from exc
    if config.model != "wright_fisher_counts":
        raise ConfigurationError(f"unsupported model: {config.model!r}")
    if config.environment.sampling_mode not in {"reservoir", "finite"}:
        raise ConfigurationError(
            "environment.sampling_mode must be 'reservoir' or 'finite'"
        )
    config.validate()
    return config
