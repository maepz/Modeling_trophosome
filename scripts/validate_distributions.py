#!/usr/bin/env python3
"""Validate exact-count transitions against known probability distributions.

This is a scientific validation harness, not a benchmark.  It repeatedly calls
the maintained count kernel and compares empirical distributions with analytic
Binomial/Hypergeometric laws and with a deliberately simple cell-level
Wright--Fisher reference implementation.
"""

from __future__ import annotations

import argparse
import json
import math
import platform
from dataclasses import asdict, dataclass
from datetime import UTC, datetime
from pathlib import Path

import numpy as np

from trophosome import MODEL_SPEC_VERSION, OUTPUT_SCHEMA_VERSION, __version__
from trophosome.config import EvolutionConfig, HostConfig
from trophosome.count_model import (
    IdAllocator,
    LineageRecorder,
    PopulationState,
    fixed_pool_migration_step,
    free_living_selection_step,
    merge_populations,
    population_size_schedule,
    proportional_rescale,
    sample_population,
    simulate_within_host,
    wright_fisher_step,
)
from trophosome.simulation import _reservoir_founders


@dataclass(frozen=True)
class ValidationResult:
    """One declared distributional comparison and its acceptance decision."""

    name: str
    passed: bool
    observations: int
    metrics: dict[str, float]
    criteria: dict[str, float]
    interpretation: str


def _count(state: PopulationState, genotype_id: int) -> int:
    matches = np.flatnonzero(state.genotype_ids == genotype_id)
    return 0 if len(matches) == 0 else int(state.counts[int(matches[0])])


def _mean_variance_result(
    name: str,
    values: np.ndarray,
    expected_mean: float,
    expected_variance: float,
    interpretation: str,
    *,
    z_limit: float = 6.0,
) -> ValidationResult:
    observations = len(values)
    observed_mean = float(np.mean(values))
    observed_variance = float(np.var(values, ddof=1))
    mean_se = math.sqrt(expected_variance / observations)
    mean_z = abs(observed_mean - expected_mean) / mean_se
    # This normal approximation is intentionally conservative.  The six-SE
    # threshold avoids making the release check sensitive to harmless RNG
    # variation while still detecting material implementation errors.
    variance_se = expected_variance * math.sqrt(2.0 / (observations - 1))
    variance_z = abs(observed_variance - expected_variance) / variance_se
    return ValidationResult(
        name=name,
        passed=mean_z <= z_limit and variance_z <= z_limit,
        observations=observations,
        metrics={
            "observed_mean": observed_mean,
            "expected_mean": expected_mean,
            "mean_z": mean_z,
            "observed_variance": observed_variance,
            "expected_variance": expected_variance,
            "variance_z": variance_z,
        },
        criteria={"maximum_absolute_z": z_limit},
        interpretation=interpretation,
    )


def _neutral_step(repetitions: int, seed: int) -> ValidationResult:
    state = PopulationState.from_counts([30, 70], [1.0, 1.0])
    evolution = EvolutionConfig(
        mutation_probability=0.0,
        within_host_selection=False,
    )
    rng = np.random.default_rng(seed)
    values = np.empty(repetitions, dtype=np.float64)
    for index in range(repetitions):
        observed, _ = wright_fisher_step(state, 80, evolution, rng, IdAllocator(2))
        values[index] = _count(observed, 0)
    probability = 0.3
    return _mean_variance_result(
        "neutral Wright-Fisher drift",
        values,
        80 * probability,
        80 * probability * (1.0 - probability),
        "Offspring count follows Binomial(N, parental frequency).",
    )


def _selected_step(repetitions: int, seed: int) -> ValidationResult:
    state = PopulationState.from_counts([30, 70], [2.0, 1.0])
    evolution = EvolutionConfig(
        mutation_probability=0.0,
        within_host_selection=True,
    )
    rng = np.random.default_rng(seed)
    values = np.empty(repetitions, dtype=np.float64)
    for index in range(repetitions):
        observed, _ = wright_fisher_step(state, 60, evolution, rng, IdAllocator(2))
        values[index] = _count(observed, 0)
    probability = (30 * 2.0) / (30 * 2.0 + 70)
    return _mean_variance_result(
        "fitness-weighted Wright-Fisher reproduction",
        values,
        60 * probability,
        60 * probability * (1.0 - probability),
        "Selection uses abundance multiplied by relative fitness.",
    )


def _free_living_selected_step(repetitions: int, seed: int) -> ValidationResult:
    state = PopulationState.from_counts(
        [30, 70],
        [1.0, 2.0],
        [2.0, 1.0],
    )
    rng = np.random.default_rng(seed)
    values = np.empty(repetitions, dtype=np.float64)
    for index in range(repetitions):
        observed = free_living_selection_step(state, 60, rng)
        values[index] = _count(observed, 0)
    probability = (30 * 2.0) / (30 * 2.0 + 70)
    return _mean_variance_result(
        "free-living fitness-weighted reproduction",
        values,
        60 * probability,
        60 * probability * (1.0 - probability),
        (
            "Environmental selection uses abundance multiplied by free-living "
            "fitness, independently of within-host fitness."
        ),
    )


def _dual_fitness_effects(repetitions: int, seed: int) -> ValidationResult:
    effect_mean = -0.01
    effect_sd = 0.02
    state = PopulationState.from_counts([1], [1.0], [1.0])
    evolution = EvolutionConfig(
        mutation_probability=1.0,
        mutation_effect_mean=effect_mean,
        mutation_effect_sd=effect_sd,
        within_host_selection=False,
        max_materialized_mutants=repetitions,
    )
    observed, mutants = wright_fisher_step(
        state,
        repetitions,
        evolution,
        np.random.default_rng(seed),
        IdAllocator(1),
        free_living_fitness_rng=np.random.default_rng(seed + 1),
    )
    if mutants != repetitions:
        raise AssertionError("all offspring should mutate in this validation")
    host_effect = observed.within_host_fitness - 1.0
    environment_effect = observed.free_living_fitness - 1.0
    expected_variance = effect_sd**2
    mean_se = effect_sd / math.sqrt(repetitions)
    variance_se = expected_variance * math.sqrt(2.0 / (repetitions - 1))
    host_mean_z = abs(float(np.mean(host_effect)) - effect_mean) / mean_se
    environment_mean_z = abs(float(np.mean(environment_effect)) - effect_mean) / mean_se
    host_variance = float(np.var(host_effect, ddof=1))
    environment_variance = float(np.var(environment_effect, ddof=1))
    host_variance_z = abs(host_variance - expected_variance) / variance_se
    environment_variance_z = abs(environment_variance - expected_variance) / variance_se
    correlation = float(np.corrcoef(host_effect, environment_effect)[0, 1])
    correlation_z = abs(np.arctanh(correlation)) * math.sqrt(repetitions - 3)
    maximum_z = max(
        host_mean_z,
        environment_mean_z,
        host_variance_z,
        environment_variance_z,
        correlation_z,
    )
    z_limit = 6.0
    return ValidationResult(
        name="independent dual-habitat mutation fitness effects",
        passed=bool(maximum_z <= z_limit),
        observations=repetitions,
        metrics={
            "within_host_effect_mean": float(np.mean(host_effect)),
            "free_living_effect_mean": float(np.mean(environment_effect)),
            "within_host_effect_variance": host_variance,
            "free_living_effect_variance": environment_variance,
            "effect_correlation": correlation,
            "maximum_z": maximum_z,
        },
        criteria={"maximum_absolute_z": z_limit},
        interpretation=(
            "New strains inherit two independently sampled fitness effects from "
            "the same declared normal distribution."
        ),
    )


def _mutation_count(repetitions: int, seed: int) -> tuple[ValidationResult, ...]:
    state = PopulationState.from_counts([40, 60], [1.0, 1.0])
    mutation_probability = 0.1
    evolution = EvolutionConfig(
        mutation_probability=mutation_probability,
        mutation_effect_mean=0.0,
        mutation_effect_sd=0.0,
        within_host_selection=False,
        max_materialized_mutants=100,
    )
    rng = np.random.default_rng(seed)
    total = np.empty(repetitions, dtype=np.float64)
    parent_zero = np.empty(repetitions, dtype=np.float64)
    for index in range(repetitions):
        recorder = LineageRecorder.from_founders(state.genotype_ids)
        _, total_mutants = wright_fisher_step(
            state,
            50,
            evolution,
            rng,
            IdAllocator(2),
            lineage_recorder=recorder,
        )
        total[index] = total_mutants
        parent_zero[index] = sum(
            event.parent_strain_id == 0 for event in recorder.events
        )
    total_probability = mutation_probability
    parent_zero_probability = 0.4 * mutation_probability
    return (
        _mean_variance_result(
            "mutation-event count",
            total,
            50 * total_probability,
            50 * total_probability * (1.0 - total_probability),
            "Every offspring mutates independently with the declared probability.",
        ),
        _mean_variance_result(
            "mutation-parent assignment",
            parent_zero,
            50 * parent_zero_probability,
            50 * parent_zero_probability * (1.0 - parent_zero_probability),
            "Mutant parent identities follow reproductive parental weights.",
        ),
    )


def _sampling_with_replacement(repetitions: int, seed: int) -> ValidationResult:
    state = PopulationState.from_counts([20, 80], [1.0, 1.0])
    rng = np.random.default_rng(seed)
    values = np.empty(repetitions, dtype=np.float64)
    for index in range(repetitions):
        values[index] = _count(sample_population(state, 25, rng, replace=True), 0)
    probability = 0.2
    return _mean_variance_result(
        "with-replacement population sampling",
        values,
        25 * probability,
        25 * probability * (1.0 - probability),
        "Reservoir-style sampling follows the Multinomial/Binomial law.",
    )


def _sampling_without_replacement(repetitions: int, seed: int) -> ValidationResult:
    state = PopulationState.from_counts([20, 80], [1.0, 1.0])
    rng = np.random.default_rng(seed)
    sample_size = 25
    values = np.empty(repetitions, dtype=np.float64)
    for index in range(repetitions):
        values[index] = _count(
            sample_population(state, sample_size, rng, replace=False), 0
        )
    probability = 0.2
    expected_variance = (
        sample_size
        * probability
        * (1.0 - probability)
        * (100 - sample_size)
        / (100 - 1)
    )
    return _mean_variance_result(
        "without-replacement population sampling",
        values,
        sample_size * probability,
        expected_variance,
        "Finite-population sampling follows the Hypergeometric law.",
    )


def _optimized_reservoir_founders(repetitions: int, seed: int) -> ValidationResult:
    state = PopulationState.from_counts([20, 30, 50], [1.0, 1.0, 1.0])
    cumulative = np.cumsum(state.counts, dtype=np.float64)
    cumulative /= cumulative[-1]
    rng = np.random.default_rng(seed)
    sample_size = 8
    values = np.empty(repetitions, dtype=np.float64)
    for index in range(repetitions):
        founders = _reservoir_founders(state, cumulative, sample_size, rng)
        values[index] = _count(founders, 0)
    probability = 0.2
    return _mean_variance_result(
        "optimized reservoir-founder sampling",
        values,
        sample_size * probability,
        sample_size * probability * (1.0 - probability),
        "The inverse-CDF optimization preserves reservoir founder sampling.",
    )


def _environmental_apportionment(repetitions: int, seed: int) -> ValidationResult:
    """Check that exact fractional ties do not privilege a strain label."""

    environment = PopulationState.from_counts([30, 70], [1.0, 1.0])
    rng = np.random.default_rng(seed)
    values = np.empty(repetitions, dtype=np.float64)
    for index in range(repetitions):
        release_counts = rng.multinomial(100, [0.3, 0.7])
        release = PopulationState.from_counts(release_counts.tolist(), [1.0, 1.0])
        updated = proportional_rescale(
            merge_populations(environment, release), 100, rng
        )
        values[index] = _count(updated, 0)
    # Let X be the returned strain-0 count. Before rounding, the updated count
    # is (30 + X) / 2. Odd X produces an exact half-cell tie and the seeded
    # tie-break contributes conditional variance 1/4 without changing the mean.
    probability_odd = (1.0 - (1.0 - 2.0 * 0.3) ** 100) / 2.0
    expected_variance = 100 * 0.3 * 0.7 / 4.0 + probability_odd / 4.0
    return _mean_variance_result(
        "environmental apportionment label neutrality",
        values,
        30.0,
        expected_variance,
        "Capacity-regulation ties are exchangeable with respect to strain ID.",
    )


def _fixed_pool_migration(
    repetitions: int, seed: int
) -> tuple[ValidationResult, ValidationResult]:
    focal = PopulationState.from_counts([30, 70], [1.0, 1.0])
    regional = PopulationState.from_counts([20, 80], [1.0, 1.0])
    sample_size = 20
    emigrant_values = np.empty(repetitions, dtype=np.float64)
    immigrant_values = np.empty(repetitions, dtype=np.float64)
    emigration_rng = np.random.default_rng(seed)
    immigration_rng = np.random.default_rng(seed + 1)
    for index in range(repetitions):
        migrated, emigrants, immigrants = fixed_pool_migration_step(
            focal,
            regional,
            sample_size,
            emigration_rng,
            immigration_rng,
        )
        if migrated.size != focal.size or emigrants is None or immigrants is None:
            raise AssertionError("migration must exchange equal non-zero samples")
        emigrant_values[index] = _count(emigrants, 0)
        immigrant_values[index] = _count(immigrants, 0)
    focal_probability = 0.3
    regional_probability = 0.2
    emigrant_variance = (
        sample_size
        * focal_probability
        * (1.0 - focal_probability)
        * (focal.size - sample_size)
        / (focal.size - 1)
    )
    return (
        _mean_variance_result(
            "fixed-pool focal emigration",
            emigrant_values,
            sample_size * focal_probability,
            emigrant_variance,
            (
                "Focal emigrants are sampled without replacement and follow the "
                "Hypergeometric law."
            ),
        ),
        _mean_variance_result(
            "fixed-pool regional immigration",
            immigrant_values,
            sample_size * regional_probability,
            sample_size * regional_probability * (1.0 - regional_probability),
            (
                "Immigrants are sampled with replacement from the unchanged "
                "regional composition and follow the Binomial law."
            ),
        ),
    )


def _binomial_probability(size: int, count: int, probability: float) -> float:
    if probability == 0.0:
        return float(count == 0)
    if probability == 1.0:
        return float(count == size)
    return (
        math.comb(size, count)
        * probability**count
        * (1.0 - probability) ** (size - count)
    )


def _neutral_endpoint_distribution(
    initial_count: int, initial_size: int, targets: tuple[int, ...]
) -> np.ndarray:
    distribution = np.zeros(initial_size + 1, dtype=np.float64)
    distribution[initial_count] = 1.0
    current_size = initial_size
    for target_size in targets:
        updated = np.zeros(target_size + 1, dtype=np.float64)
        for parent_count, parent_probability in enumerate(distribution):
            if parent_probability == 0.0:
                continue
            probability = parent_count / current_size
            for child_count in range(target_size + 1):
                updated[child_count] += parent_probability * _binomial_probability(
                    target_size, child_count, probability
                )
        distribution = updated
        current_size = target_size
    return distribution


def _multigeneration_drift(repetitions: int, seed: int) -> ValidationResult:
    state = PopulationState.from_counts([3, 7], [1.0, 1.0])
    evolution = EvolutionConfig(
        mutation_probability=0.0,
        within_host_selection=False,
    )
    targets = (10, 10, 10)
    expected = _neutral_endpoint_distribution(3, 10, targets)
    rng = np.random.default_rng(seed)
    observed = np.zeros(11, dtype=np.float64)
    for _ in range(repetitions):
        current = state
        for target in targets:
            current, _ = wright_fisher_step(
                current, target, evolution, rng, IdAllocator(2)
            )
        observed[_count(current, 0)] += 1
    observed /= repetitions
    total_variation = 0.5 * float(np.abs(observed - expected).sum())
    eligible = expected * repetitions >= 20
    standardized = np.abs(observed[eligible] - expected[eligible]) / np.sqrt(
        expected[eligible] * (1.0 - expected[eligible]) / repetitions
    )
    maximum_z = float(standardized.max(initial=0.0))
    tv_limit = max(0.015, 1.5 / math.sqrt(repetitions))
    z_limit = 6.0
    return ValidationResult(
        name="multi-generation neutral drift distribution",
        passed=total_variation <= tv_limit and maximum_z <= z_limit,
        observations=repetitions,
        metrics={
            "total_variation_distance": total_variation,
            "maximum_bin_z": maximum_z,
        },
        criteria={
            "maximum_total_variation_distance": tv_limit,
            "maximum_absolute_z": z_limit,
        },
        interpretation=(
            "Three exact count transitions match the recursively calculated "
            "Wright-Fisher endpoint distribution."
        ),
    )


def _individual_reference_features(
    host: HostConfig,
    mutation_probability: float,
    repetitions: int,
    seed: int,
) -> np.ndarray:
    rng = np.random.default_rng(seed)
    schedule = population_size_schedule(host)
    features = np.empty((repetitions, 3), dtype=np.float64)
    for replicate in range(repetitions):
        cells = np.zeros(host.infection_bottleneck, dtype=np.int64)
        next_id = 1
        for target_size in schedule:
            cells = cells[rng.integers(0, len(cells), size=target_size)].copy()
            mutated = rng.random(target_size) < mutation_probability
            mutation_count = int(mutated.sum())
            cells[mutated] = np.arange(next_id, next_id + mutation_count)
            next_id += mutation_count
        root_count = int(np.count_nonzero(cells == 0))
        non_root = cells[cells != 0]
        if len(non_root):
            _, non_root_counts = np.unique(non_root, return_counts=True)
            largest_mutant_clone = int(non_root_counts.max())
        else:
            largest_mutant_clone = 0
        features[replicate] = (
            len(np.unique(cells)),
            root_count,
            largest_mutant_clone,
        )
    return features


def _count_kernel_features(
    host: HostConfig,
    mutation_probability: float,
    repetitions: int,
    seed: int,
) -> np.ndarray:
    rng = np.random.default_rng(seed)
    schedule = population_size_schedule(host)
    evolution = EvolutionConfig(
        mutation_probability=mutation_probability,
        mutation_effect_mean=0.0,
        mutation_effect_sd=0.0,
        within_host_selection=False,
        max_materialized_mutants=host.carrying_capacity,
    )
    founders = PopulationState.from_counts([host.infection_bottleneck], [1.0])
    features = np.empty((repetitions, 3), dtype=np.float64)
    for replicate in range(repetitions):
        adult, _ = simulate_within_host(
            founders,
            host,
            evolution,
            rng,
            IdAllocator(1),
            size_schedule=schedule,
            record_summaries=False,
        )
        root_count = _count(adult, 0)
        mutant_counts = adult.counts[adult.genotype_ids != 0]
        largest_mutant_clone = int(mutant_counts.max(initial=0))
        features[replicate] = (
            adult.richness,
            root_count,
            largest_mutant_clone,
        )
    return features


def _cell_reference_comparison(repetitions: int, seed: int) -> ValidationResult:
    reference_repetitions = max(2_000, repetitions // 2)
    host = HostConfig(
        population_size=1,
        infection_bottleneck=2,
        carrying_capacity=16,
        growth_factor=2.0,
        steady_generations=2,
        host_generations=1,
        escape_fraction=0.0,
    )
    mutation_probability = 0.04
    count_features = _count_kernel_features(
        host, mutation_probability, reference_repetitions, seed
    )
    cell_features = _individual_reference_features(
        host, mutation_probability, reference_repetitions, seed + 1
    )
    labels = ("richness", "founder_cells", "largest_mutant_clone")
    metrics: dict[str, float] = {}
    maximum_z = 0.0
    for column, label in enumerate(labels):
        count_values = count_features[:, column]
        cell_values = cell_features[:, column]
        difference = float(np.mean(count_values) - np.mean(cell_values))
        standard_error = math.sqrt(
            np.var(count_values, ddof=1) / reference_repetitions
            + np.var(cell_values, ddof=1) / reference_repetitions
        )
        z_score = abs(difference) / standard_error
        metrics[f"{label}_mean_difference"] = difference
        metrics[f"{label}_z"] = z_score
        maximum_z = max(maximum_z, z_score)
    richness_bins = np.arange(1, host.carrying_capacity + 2)
    count_histogram, _ = np.histogram(count_features[:, 0], bins=richness_bins)
    cell_histogram, _ = np.histogram(cell_features[:, 0], bins=richness_bins)
    total_variation = 0.5 * float(
        np.abs(
            count_histogram / reference_repetitions
            - cell_histogram / reference_repetitions
        ).sum()
    )
    tv_limit = max(0.03, 2.0 / math.sqrt(reference_repetitions))
    z_limit = 6.0
    metrics["richness_total_variation_distance"] = total_variation
    metrics["maximum_mean_z"] = maximum_z
    return ValidationResult(
        name="cell-level reference trajectory",
        passed=maximum_z <= z_limit and total_variation <= tv_limit,
        observations=2 * reference_repetitions,
        metrics=metrics,
        criteria={
            "maximum_absolute_z": z_limit,
            "maximum_richness_total_variation_distance": tv_limit,
        },
        interpretation=(
            "Growth, mutation timing, drift, and mutant jackpot clones agree "
            "with an explicit individual-cell implementation."
        ),
    )


def run_validation(
    repetitions: int = 20_000, seed: int = 20260810
) -> list[ValidationResult]:
    """Run every release-gating distributional validation."""

    if repetitions < 2_000:
        raise ValueError("at least 2,000 repetitions are required")
    results = [
        _neutral_step(repetitions, seed),
        _selected_step(repetitions, seed + 1),
        _free_living_selected_step(repetitions, seed + 2),
        _dual_fitness_effects(repetitions, seed + 3),
        *_mutation_count(repetitions, seed + 4),
        _sampling_with_replacement(repetitions, seed + 5),
        _sampling_without_replacement(repetitions, seed + 6),
        _optimized_reservoir_founders(repetitions, seed + 7),
        _environmental_apportionment(repetitions, seed + 8),
        *_fixed_pool_migration(repetitions, seed + 9),
        _multigeneration_drift(repetitions, seed + 10),
        _cell_reference_comparison(repetitions, seed + 11),
    ]
    return results


def _payload(results: list[ValidationResult], seed: int) -> dict[str, object]:
    return {
        "created_at": datetime.now(UTC).isoformat(),
        "model_family": "wright_fisher_counts",
        "model_spec_version": MODEL_SPEC_VERSION,
        "output_schema_version": OUTPUT_SCHEMA_VERSION,
        "software_version": __version__,
        "seed": seed,
        "python": platform.python_version(),
        "numpy": np.__version__,
        "passed": all(result.passed for result in results),
        "checks": [asdict(result) for result in results],
    }


def _markdown(payload: dict[str, object]) -> str:
    checks = payload["checks"]
    assert isinstance(checks, list)
    lines = [
        "# Distributional validation report",
        "",
        f"- Overall result: **{'PASS' if payload['passed'] else 'FAIL'}**",
        f"- Model family: `{payload['model_family']}`",
        f"- Model specification: `{payload['model_spec_version']}`",
        f"- Software: `{payload['software_version']}`",
        f"- Output schema: `{payload['output_schema_version']}`",
        f"- Seed: `{payload['seed']}`",
        f"- Runtime: Python {payload['python']}, NumPy {payload['numpy']}",
        f"- Generated: {payload['created_at']}",
        "",
        "| Check | Result | Draws | Main diagnostic |",
        "|---|---:|---:|---|",
    ]
    for check in checks:
        assert isinstance(check, dict)
        metrics = check["metrics"]
        assert isinstance(metrics, dict)
        diagnostic = ", ".join(
            f"{key}={float(value):.4g}" for key, value in metrics.items()
        )
        lines.append(
            f"| {check['name']} | {'PASS' if check['passed'] else 'FAIL'} | "
            f"{check['observations']:,} | {diagnostic} |"
        )
    lines.extend(
        [
            "",
            "## Scope",
            "",
            "These checks validate the declared stochastic law: neutral and selected "
            "within-host Wright–Fisher reproduction, independent dual-habitat fitness "
            "effects, free-living selection, Bernoulli infinite-alleles mutation, "
            "mutation parentage, reservoir founder sampling, finite escape sampling, "
            "fixed-pool emigration and immigration, multi-generation drift, and "
            "mutation-timing jackpot effects.",
            "",
            "Hamilton environmental capacity regulation is checked for exact capacity, "
            "seeded label-neutral tie resolution, and no-return invariance when "
            "free-living selection is disabled. It is not "
            "yet compared with a biological equilibrium distribution; that belongs to "
            "the subsequent Phase 1 experiment.",
            "",
            "Acceptance limits were declared in the validation program and are stored "
            "with every result. They are deliberately conservative six-standard-error "
            "release gates, supplemented by total-variation limits for complete "
            "distributions.",
            "",
        ]
    )
    return "\n".join(lines)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repetitions", type=int, default=20_000)
    parser.add_argument("--seed", type=int, default=20260810)
    parser.add_argument("--json", type=Path)
    parser.add_argument("--markdown", type=Path)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    results = run_validation(args.repetitions, args.seed)
    payload = _payload(results, args.seed)
    if args.json is not None:
        args.json.parent.mkdir(parents=True, exist_ok=True)
        args.json.write_text(json.dumps(payload, indent=2) + "\n")
    if args.markdown is not None:
        args.markdown.parent.mkdir(parents=True, exist_ok=True)
        args.markdown.write_text(_markdown(payload))
    print(_markdown(payload))
    return 0 if payload["passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
