"""Scalable exact-count Wright--Fisher dynamics.

The kernel is exact for the declared haploid, discrete-generation,
infinite-alleles model. Census size is represented by integer counts and work
scales with extant strain richness and materialised mutation events, not with
the number of bacterial cells.

Mutable buffers are used only inside a host simulation. Public population
snapshots remain immutable arrays, and mutation ancestry is written to a
separate append-only event stream rather than copied with every generation.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field

import numpy as np
from numpy.typing import NDArray

from trophosome.config import EvolutionConfig, HostConfig

IntArray = NDArray[np.int64]
FloatArray = NDArray[np.float64]


class MutationMaterializationError(RuntimeError):
    """Raised when an exact infinite-alleles step would create too many strains."""


@dataclass(frozen=True)
class PopulationState:
    """An immutable snapshot of extant strain counts.

    Only fields used by population transitions are retained here. Parentage and
    mutational depth live in :class:`LineageRecorder` and the lineage-event
    output table.
    """

    genotype_ids: IntArray
    counts: IntArray
    within_host_fitness: FloatArray
    free_living_fitness: FloatArray

    def __post_init__(self) -> None:
        arrays = (
            self.genotype_ids,
            self.counts,
            self.within_host_fitness,
            self.free_living_fitness,
        )
        lengths = {len(array) for array in arrays}
        if len(lengths) != 1:
            raise ValueError("population arrays must have equal length")
        if len(self.counts) == 0:
            raise ValueError("population must contain at least one extant genotype")
        if np.any(self.counts <= 0):
            raise ValueError("extant genotype counts must be positive")
        if len(np.unique(self.genotype_ids)) != len(self.genotype_ids):
            raise ValueError("genotype IDs must be unique")
        for name, values in (
            ("within-host fitness", self.within_host_fitness),
            ("free-living fitness", self.free_living_fitness),
        ):
            if np.any(~np.isfinite(values)) or np.any(values < 0):
                raise ValueError(f"{name} must be finite and non-negative")

    @classmethod
    def from_counts(
        cls,
        counts: list[int] | tuple[int, ...],
        within_host_fitness: list[float] | tuple[float, ...],
        free_living_fitness: list[float] | tuple[float, ...] | None = None,
    ) -> PopulationState:
        count_array = np.asarray(counts, dtype=np.int64)
        within_host_array = np.asarray(within_host_fitness, dtype=np.float64)
        if free_living_fitness is None:
            free_living_array = np.ones(len(count_array), dtype=np.float64)
        else:
            free_living_array = np.asarray(free_living_fitness, dtype=np.float64)
        ids = np.arange(len(count_array), dtype=np.int64)
        if not (len(count_array) == len(within_host_array) == len(free_living_array)):
            raise ValueError("population arrays must have equal length")
        if np.any(count_array < 0):
            raise ValueError("declared genotype counts must be non-negative")
        keep = count_array > 0
        if not np.any(keep):
            raise ValueError("population must contain at least one extant genotype")
        return cls(
            ids[keep],
            count_array[keep],
            within_host_array[keep],
            free_living_array[keep],
        )

    @classmethod
    def _trusted(
        cls,
        genotype_ids: IntArray,
        counts: IntArray,
        within_host_fitness: FloatArray,
        free_living_fitness: FloatArray,
    ) -> PopulationState:
        """Construct an internally validated snapshot without repeated scans."""

        instance = object.__new__(cls)
        object.__setattr__(instance, "genotype_ids", genotype_ids)
        object.__setattr__(instance, "counts", counts)
        object.__setattr__(instance, "within_host_fitness", within_host_fitness)
        object.__setattr__(instance, "free_living_fitness", free_living_fitness)
        return instance

    @property
    def fitness(self) -> FloatArray:
        """Backward-compatible alias for within-host fitness."""

        return self.within_host_fitness

    @property
    def size(self) -> int:
        return int(self.counts.sum())

    @property
    def richness(self) -> int:
        return len(self.counts)


@dataclass
class IdAllocator:
    """Allocate deterministic strain IDs, optionally within a reserved block."""

    next_id: int
    stop_id: int | None = None

    def take(self, count: int) -> IntArray:
        if count < 0:
            raise ValueError("ID count must be non-negative")
        stop = self.next_id + count
        if self.stop_id is not None and stop > self.stop_id:
            raise MutationMaterializationError(
                "host exceeded its deterministic strain-ID reservation"
            )
        values = np.arange(self.next_id, stop, dtype=np.int64)
        self.next_id = stop
        return values


@dataclass(frozen=True)
class LineageEvent:
    """One strain-changing mutation event within a host."""

    strain_id: int
    parent_strain_id: int
    mutational_depth: int
    within_host_generation: int
    within_host_fitness: float
    free_living_fitness: float


@dataclass
class LineageRecorder:
    """Append mutation events without burdening the active population state."""

    depth_by_id: dict[int, int]
    events: list[LineageEvent] = field(default_factory=list)

    @classmethod
    def from_founders(
        cls, founder_ids: IntArray, founder_depths: IntArray | None = None
    ) -> LineageRecorder:
        if founder_depths is None:
            founder_depths = np.zeros(len(founder_ids), dtype=np.int64)
        if len(founder_ids) != len(founder_depths):
            raise ValueError("founder IDs and depths must have equal length")
        return cls(
            {
                int(genotype_id): int(depth)
                for genotype_id, depth in zip(
                    founder_ids.tolist(), founder_depths.tolist(), strict=True
                )
            }
        )

    def record(
        self,
        strain_ids: IntArray,
        parent_ids: IntArray,
        within_host_fitness: FloatArray,
        free_living_fitness: FloatArray,
        generation: int,
    ) -> None:
        for strain_id, parent_id, host_value, environment_value in zip(
            strain_ids.tolist(),
            parent_ids.tolist(),
            within_host_fitness.tolist(),
            free_living_fitness.tolist(),
            strict=True,
        ):
            parent_depth = self.depth_by_id[int(parent_id)]
            depth = parent_depth + 1
            self.depth_by_id[int(strain_id)] = depth
            self.events.append(
                LineageEvent(
                    strain_id=int(strain_id),
                    parent_strain_id=int(parent_id),
                    mutational_depth=depth,
                    within_host_generation=generation,
                    within_host_fitness=float(host_value),
                    free_living_fitness=float(environment_value),
                )
            )


@dataclass(frozen=True)
class GenerationSummary:
    generation: int
    population_size: int
    richness: int
    mean_within_host_fitness: float
    mean_free_living_fitness: float
    new_mutants: int


def summarize(
    state: PopulationState, generation: int, new_mutants: int = 0
) -> GenerationSummary:
    return GenerationSummary(
        generation=generation,
        population_size=state.size,
        richness=state.richness,
        mean_within_host_fitness=float(
            np.average(state.within_host_fitness, weights=state.counts)
        ),
        mean_free_living_fitness=float(
            np.average(state.free_living_fitness, weights=state.counts)
        ),
        new_mutants=new_mutants,
    )


class _PopulationBuffer:
    """Capacity-managed mutable storage used by one host simulation."""

    def __init__(self, state: PopulationState) -> None:
        self.length = state.richness
        capacity = max(8, self.length)
        self.genotype_ids = np.empty(capacity, dtype=np.int64)
        self.counts = np.empty(capacity, dtype=np.int64)
        self.within_host_fitness = np.empty(capacity, dtype=np.float64)
        self.free_living_fitness = np.empty(capacity, dtype=np.float64)
        self.genotype_ids[: self.length] = state.genotype_ids
        self.counts[: self.length] = state.counts
        self.within_host_fitness[: self.length] = state.within_host_fitness
        self.free_living_fitness[: self.length] = state.free_living_fitness

    def _ensure_capacity(self, required: int) -> None:
        if required <= len(self.counts):
            return
        capacity = max(required, len(self.counts) * 2)
        ids = np.empty(capacity, dtype=np.int64)
        counts = np.empty(capacity, dtype=np.int64)
        within_host_fitness = np.empty(capacity, dtype=np.float64)
        free_living_fitness = np.empty(capacity, dtype=np.float64)
        ids[: self.length] = self.genotype_ids[: self.length]
        counts[: self.length] = self.counts[: self.length]
        within_host_fitness[: self.length] = self.within_host_fitness[: self.length]
        free_living_fitness[: self.length] = self.free_living_fitness[: self.length]
        self.genotype_ids = ids
        self.counts = counts
        self.within_host_fitness = within_host_fitness
        self.free_living_fitness = free_living_fitness

    @property
    def size(self) -> int:
        return int(self.counts[: self.length].sum())

    def snapshot(self) -> PopulationState:
        return PopulationState._trusted(
            self.genotype_ids[: self.length].copy(),
            self.counts[: self.length].copy(),
            self.within_host_fitness[: self.length].copy(),
            self.free_living_fitness[: self.length].copy(),
        )

    def step(
        self,
        target_size: int,
        evolution: EvolutionConfig,
        rng: np.random.Generator,
        id_allocator: IdAllocator,
        *,
        generation: int,
        lineage_recorder: LineageRecorder | None,
        free_living_fitness_rng: np.random.Generator | None,
    ) -> int:
        if target_size < 1:
            raise ValueError("target_size must be positive")

        ids = self.genotype_ids[: self.length]
        counts = self.counts[: self.length]
        within_host_fitness = self.within_host_fitness[: self.length]
        free_living_fitness = self.free_living_fitness[: self.length]
        if evolution.within_host_selection:
            weights = counts.astype(np.float64) * within_host_fitness
        else:
            weights = counts.astype(np.float64)
        total_weight = float(weights.sum())
        if not math.isfinite(total_weight) or total_weight <= 0:
            raise ValueError("population has no positive reproductive weight")

        offspring = rng.multinomial(target_size, weights / total_weight).astype(
            np.int64
        )
        if evolution.mutation_probability == 0:
            mutant_counts = np.zeros(self.length, dtype=np.int64)
            total_mutants = 0
        else:
            mutant_counts = rng.binomial(
                offspring, evolution.mutation_probability
            ).astype(np.int64)
            total_mutants = int(mutant_counts.sum())
        if total_mutants > evolution.max_materialized_mutants:
            raise MutationMaterializationError(
                f"step produced {total_mutants:,} mutants; the exact infinite-alleles "
                f"limit is {evolution.max_materialized_mutants:,}"
            )

        retained = offspring - mutant_counts
        keep = retained > 0
        survivor_ids = ids[keep].copy()
        survivor_counts = retained[keep]
        survivor_within_host_fitness = within_host_fitness[keep].copy()
        survivor_free_living_fitness = free_living_fitness[keep].copy()

        if total_mutants:
            parent_index = np.repeat(
                np.arange(self.length, dtype=np.int64), mutant_counts
            )
            parent_ids = ids[parent_index].copy()
            parent_within_host_fitness = within_host_fitness[parent_index].copy()
            parent_free_living_fitness = free_living_fitness[parent_index].copy()
            mutant_within_host_fitness = parent_within_host_fitness + rng.normal(
                evolution.mutation_effect_mean,
                evolution.mutation_effect_sd,
                size=total_mutants,
            )
            mutant_within_host_fitness = np.maximum(
                mutant_within_host_fitness, evolution.fitness_floor
            ).astype(np.float64)
            if free_living_fitness_rng is None:
                free_living_fitness_rng = rng
            mutant_free_living_fitness = parent_free_living_fitness + (
                free_living_fitness_rng.normal(
                    evolution.mutation_effect_mean,
                    evolution.mutation_effect_sd,
                    size=total_mutants,
                )
            )
            mutant_free_living_fitness = np.maximum(
                mutant_free_living_fitness, evolution.fitness_floor
            ).astype(np.float64)
            mutant_ids = id_allocator.take(total_mutants)
        else:
            parent_ids = np.empty(0, dtype=np.int64)
            mutant_within_host_fitness = np.empty(0, dtype=np.float64)
            mutant_free_living_fitness = np.empty(0, dtype=np.float64)
            mutant_ids = np.empty(0, dtype=np.int64)

        survivors = len(survivor_counts)
        new_length = survivors + total_mutants
        self._ensure_capacity(new_length)
        self.genotype_ids[:survivors] = survivor_ids
        self.counts[:survivors] = survivor_counts
        self.within_host_fitness[:survivors] = survivor_within_host_fitness
        self.free_living_fitness[:survivors] = survivor_free_living_fitness
        if total_mutants:
            mutant_slice = slice(survivors, new_length)
            self.genotype_ids[mutant_slice] = mutant_ids
            self.counts[mutant_slice] = 1
            self.within_host_fitness[mutant_slice] = mutant_within_host_fitness
            self.free_living_fitness[mutant_slice] = mutant_free_living_fitness
            if lineage_recorder is not None:
                lineage_recorder.record(
                    mutant_ids,
                    parent_ids,
                    mutant_within_host_fitness,
                    mutant_free_living_fitness,
                    generation,
                )
        self.length = new_length
        return total_mutants


def wright_fisher_step(
    state: PopulationState,
    target_size: int,
    evolution: EvolutionConfig,
    rng: np.random.Generator,
    id_allocator: IdAllocator,
    *,
    generation: int = 1,
    lineage_recorder: LineageRecorder | None = None,
    free_living_fitness_rng: np.random.Generator | None = None,
) -> tuple[PopulationState, int]:
    """Draw one exact haploid Wright--Fisher generation at count resolution."""

    buffer = _PopulationBuffer(state)
    mutants = buffer.step(
        target_size,
        evolution,
        rng,
        id_allocator,
        generation=generation,
        lineage_recorder=lineage_recorder,
        free_living_fitness_rng=free_living_fitness_rng,
    )
    return buffer.snapshot(), mutants


def sample_population(
    state: PopulationState,
    sample_size: int,
    rng: np.random.Generator,
    *,
    replace: bool,
) -> PopulationState:
    """Sample cells while retaining strain identity and both fitness traits."""

    if sample_size < 1:
        raise ValueError("sample_size must be positive")
    if not replace and sample_size > state.size:
        raise ValueError("sample size exceeds a finite population")
    if replace:
        sampled = rng.multinomial(sample_size, state.counts / state.size)
    elif state.size >= 1_000_000_000:
        # NumPy's fast multivariate implementation rejects source totals at
        # or above 1e9.  Sequential conditional hypergeometric draws are an
        # exact construction of the same multivariate distribution and scale
        # with strain richness rather than the number of bacterial cells.
        sampled = np.zeros(len(state.counts), dtype=np.int64)
        remaining_population = state.size
        remaining_sample = sample_size
        for index, count_value in enumerate(state.counts[:-1]):
            count = int(count_value)
            if remaining_sample == 0:
                break
            if remaining_sample == remaining_population:
                sampled[index:] = state.counts[index:]
                remaining_sample = 0
                break
            drawn = int(
                rng.hypergeometric(
                    count,
                    remaining_population - count,
                    remaining_sample,
                )
            )
            sampled[index] = drawn
            remaining_sample -= drawn
            remaining_population -= count
        if remaining_sample:
            sampled[-1] = remaining_sample
    else:
        sampled = rng.multivariate_hypergeometric(state.counts, sample_size)
    keep = sampled > 0
    return PopulationState._trusted(
        state.genotype_ids[keep].copy(),
        sampled[keep].astype(np.int64),
        state.within_host_fitness[keep].copy(),
        state.free_living_fitness[keep].copy(),
    )


def subtract_population(
    population: PopulationState, sample: PopulationState
) -> PopulationState | None:
    """Remove a without-replacement sample from a finite population."""

    removed = dict(
        zip(
            sample.genotype_ids.tolist(),
            sample.counts.tolist(),
            strict=True,
        )
    )
    remaining = population.counts.copy()
    for index, genotype_id in enumerate(population.genotype_ids.tolist()):
        remaining[index] -= removed.get(genotype_id, 0)
    if np.any(remaining < 0):
        raise ValueError("sample contains more cells than the source population")
    keep = remaining > 0
    if not np.any(keep):
        return None
    return PopulationState._trusted(
        population.genotype_ids[keep].copy(),
        remaining[keep],
        population.within_host_fitness[keep].copy(),
        population.free_living_fitness[keep].copy(),
    )


def merge_populations(*populations: PopulationState) -> PopulationState:
    """Merge populations by a single sort-and-reduce over strain IDs."""

    if not populations:
        raise ValueError("at least one population is required")
    if len(populations) == 1:
        population = populations[0]
        return PopulationState._trusted(
            population.genotype_ids.copy(),
            population.counts.copy(),
            population.within_host_fitness.copy(),
            population.free_living_fitness.copy(),
        )

    ids = np.concatenate([population.genotype_ids for population in populations])
    counts = np.concatenate([population.counts for population in populations])
    within_host_fitness = np.concatenate(
        [population.within_host_fitness for population in populations]
    )
    free_living_fitness = np.concatenate(
        [population.free_living_fitness for population in populations]
    )
    order = np.argsort(ids, kind="stable")
    ids = ids[order]
    counts = counts[order]
    within_host_fitness = within_host_fitness[order]
    free_living_fitness = free_living_fitness[order]
    duplicate = ids[1:] == ids[:-1]
    conflicting_within_host = within_host_fitness[1:] != within_host_fitness[:-1]
    conflicting_free_living = free_living_fitness[1:] != free_living_fitness[:-1]
    if np.any(duplicate & (conflicting_within_host | conflicting_free_living)):
        raise ValueError("conflicting fitness metadata for a shared genotype")
    starts = np.r_[0, np.flatnonzero(~duplicate) + 1]
    return PopulationState._trusted(
        ids[starts].copy(),
        np.add.reduceat(counts, starts).astype(np.int64),
        within_host_fitness[starts].copy(),
        free_living_fitness[starts].copy(),
    )


def fixed_pool_migration_step(
    focal: PopulationState,
    regional_pool: PopulationState,
    replacement_count: int,
    emigration_rng: np.random.Generator,
    immigration_rng: np.random.Generator,
) -> tuple[PopulationState, PopulationState | None, PopulationState | None]:
    """Exchange focal cells with a fixed, non-depleting regional pool.

    ``replacement_count`` focal cells emigrate as an exact without-replacement
    sample. The same number of immigrants is independently sampled with
    replacement from the fixed regional source. The returned population has the
    same census size as ``focal``. The two optional population results record
    the realized emigrant and immigrant strain counts for output and validation.
    """

    if replacement_count < 0 or replacement_count > focal.size:
        raise ValueError("migration replacement count must be between 0 and N_E")
    if replacement_count == 0:
        return merge_populations(focal), None, None

    emigrants = sample_population(
        focal,
        replacement_count,
        emigration_rng,
        replace=False,
    )
    residents = subtract_population(focal, emigrants)
    immigrants = sample_population(
        regional_pool,
        replacement_count,
        immigration_rng,
        replace=True,
    )
    migrated = (
        immigrants if residents is None else merge_populations(residents, immigrants)
    )
    if migrated.size != focal.size:
        raise RuntimeError("migration did not preserve focal population capacity")
    return migrated, emigrants, immigrants


def proportional_rescale(
    population: PopulationState,
    target_size: int,
    rng: np.random.Generator | None = None,
) -> PopulationState:
    """Apply strain-neutral Hamilton capacity regulation.

    Hamilton apportionment preserves the mixed frequencies as closely as integer
    counts permit. When multiple strains have the same fractional remainder at
    the allocation boundary, ``rng`` selects among them without reference to
    strain ID. Supplying ``rng`` is therefore required for label-exchangeable
    stochastic simulations; omitting it retains deterministic behavior for
    backwards-compatible utility calls with no tied boundary.
    """

    if target_size < 1:
        raise ValueError("target size must be positive")
    if population.size == target_size:
        return merge_populations(population)
    scaled = population.counts.astype(np.float64) * (target_size / population.size)
    counts = np.floor(scaled).astype(np.int64)
    remainder = target_size - int(counts.sum())
    if remainder:
        fractions = scaled - counts
        order = np.argsort(-fractions, kind="stable")
        cutoff = fractions[order[remainder - 1]]
        certain = np.flatnonzero(fractions > cutoff)
        tied = np.flatnonzero(fractions == cutoff)
        needed = remainder - len(certain)
        counts[certain] += 1
        if needed == len(tied) or rng is None:
            counts[tied[:needed]] += 1
        else:
            selected = rng.choice(tied, size=needed, replace=False)
            counts[selected] += 1
    keep = counts > 0
    return PopulationState._trusted(
        population.genotype_ids[keep].copy(),
        counts[keep],
        population.within_host_fitness[keep].copy(),
        population.free_living_fitness[keep].copy(),
    )


def free_living_selection_step(
    population: PopulationState,
    target_size: int,
    rng: np.random.Generator,
) -> PopulationState:
    """Apply one exact free-living Wright--Fisher selection transition.

    Mutations remain confined to the within-host phase. The free-living step
    samples ``target_size`` descendants with probabilities proportional to
    strain abundance multiplied by free-living fitness.
    """

    if target_size < 1:
        raise ValueError("target size must be positive")
    weights = population.counts.astype(np.float64) * population.free_living_fitness
    total_weight = float(weights.sum())
    if not math.isfinite(total_weight) or total_weight <= 0:
        raise ValueError("population has no positive free-living reproductive weight")
    offspring = rng.multinomial(target_size, weights / total_weight).astype(np.int64)
    keep = offspring > 0
    return PopulationState._trusted(
        population.genotype_ids[keep].copy(),
        offspring[keep],
        population.within_host_fitness[keep].copy(),
        population.free_living_fitness[keep].copy(),
    )


def population_size_schedule(
    host: HostConfig, *, initial_size: int | None = None
) -> tuple[int, ...]:
    """Precompute the target size of every within-host transition."""

    size = host.infection_bottleneck if initial_size is None else initial_size
    if size > host.carrying_capacity:
        raise ValueError("founder population exceeds carrying capacity")
    targets: list[int] = []
    while size < host.carrying_capacity:
        size = min(
            host.carrying_capacity,
            max(size + 1, math.ceil(size * host.growth_factor)),
        )
        targets.append(size)
    targets.extend([host.carrying_capacity] * host.steady_generations)
    return tuple(targets)


def simulate_within_host(
    founders: PopulationState,
    host: HostConfig,
    evolution: EvolutionConfig,
    rng: np.random.Generator,
    id_allocator: IdAllocator | None = None,
    *,
    snapshot_interval: int = 1,
    record_summaries: bool = True,
    size_schedule: tuple[int, ...] | None = None,
    lineage_recorder: LineageRecorder | None = None,
    free_living_fitness_rng: np.random.Generator | None = None,
) -> tuple[PopulationState, tuple[GenerationSummary, ...]]:
    """Grow founders to capacity and simulate the exact adult steady phase."""

    if id_allocator is None:
        id_allocator = IdAllocator(int(founders.genotype_ids.max()) + 1)
    if snapshot_interval < 1:
        raise ValueError("snapshot interval must be positive")
    if size_schedule is None:
        size_schedule = population_size_schedule(host, initial_size=founders.size)

    buffer = _PopulationBuffer(founders)
    summaries: list[GenerationSummary] = []
    if record_summaries:
        summaries.append(summarize(founders, 0))
    for generation, target_size in enumerate(size_schedule, start=1):
        mutants = buffer.step(
            target_size,
            evolution,
            rng,
            id_allocator,
            generation=generation,
            lineage_recorder=lineage_recorder,
            free_living_fitness_rng=free_living_fitness_rng,
        )
        if record_summaries and (
            generation % snapshot_interval == 0 or generation == len(size_schedule)
        ):
            summaries.append(summarize(buffer.snapshot(), generation, mutants))
    return buffer.snapshot(), tuple(summaries)
