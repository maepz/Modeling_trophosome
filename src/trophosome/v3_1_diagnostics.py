"""Pre-run scaling diagnostics for the legacy V3.1 endpoint prototype.

This module does not run or alter V3.1.  It makes the quantities that control
the current implementation's cost explicit before a cluster job is submitted.
"""

from __future__ import annotations

import math
from dataclasses import asdict, dataclass


@dataclass(frozen=True)
class V31LoadEstimate:
    """Expected mutation workload under the assumptions currently used by V3.1."""

    infection_size: int
    carrying_capacity: int
    growth_factor: float
    growth_generations: int
    steady_generations: int
    lineage_generations: int
    mutation_probability: float
    mutation_carrier_probability: float
    expected_unmutated_cells: float
    expected_mutant_cells: float
    mutant_cells_for_screening: int
    legacy_ewens_theta: float
    expected_strain_richness_lower: float
    expected_strain_richness_upper: float
    warnings: tuple[str, ...]

    def to_dict(self) -> dict[str, object]:
        return asdict(self)


def growth_generations(
    infection_size: int, carrying_capacity: int, growth_factor: float
) -> int:
    """Return the integer growth duration used by the legacy V3.1 formula."""

    if infection_size < 1:
        raise ValueError("infection_size must be at least 1")
    if carrying_capacity < infection_size:
        raise ValueError("carrying_capacity must be at least infection_size")
    if not math.isfinite(growth_factor) or growth_factor <= 0:
        raise ValueError("growth_factor must be finite and positive")
    if carrying_capacity == infection_size:
        return 0
    if growth_factor <= 1:
        raise ValueError("growth_factor must be greater than 1 when growth is required")
    raw = math.log(carrying_capacity / infection_size) / math.log(growth_factor)
    return math.ceil(raw)


def mutation_carrier_probability(
    mutation_probability: float, lineage_generations: int
) -> float:
    """Probability of at least one mutation on a lineage of fixed length."""

    if not 0 <= mutation_probability <= 1:
        raise ValueError("mutation_probability must be between 0 and 1")
    if lineage_generations < 0:
        raise ValueError("lineage_generations must be non-negative")
    if mutation_probability == 0 or lineage_generations == 0:
        return 0.0
    if mutation_probability == 1:
        return 1.0
    return -math.expm1(lineage_generations * math.log1p(-mutation_probability))


def ewens_expected_richness_bounds(
    sample_size: int, theta: float
) -> tuple[float, float]:
    """Bounds for E[K_n] under ESF(theta), without looping over ``sample_size``.

    The exact expectation is ``sum(theta / (theta + i), i=0..n-1)``.  The
    integral bounds returned here remain inexpensive when ``n`` is enormous.
    """

    if sample_size < 0:
        raise ValueError("sample_size must be non-negative")
    if not math.isfinite(theta) or theta < 0:
        raise ValueError("theta must be finite and non-negative")
    if sample_size == 0:
        return 0.0, 0.0
    if theta == 0:
        return 1.0, 1.0
    lower = theta * math.log1p(sample_size / theta)
    upper = 1.0 + theta * math.log1p((sample_size - 1) / theta)
    return min(lower, float(sample_size)), min(upper, float(sample_size))


def estimate_v3_1_load(
    *,
    infection_size: int,
    carrying_capacity: int,
    growth_factor: float,
    steady_generations: int,
    mutation_probability: float,
    crp_cell_warning: int = 100_000,
    strain_warning: int = 100_000,
) -> V31LoadEstimate:
    """Estimate the two output-size limits in the current V3.1 implementation.

    V3.1 uses ``theta = carrying_capacity * mutation_probability``.  This
    function reports that legacy value without endorsing it as a calibrated
    population-scaled mutation parameter.
    """

    if steady_generations < 0:
        raise ValueError("steady_generations must be non-negative")
    growth = growth_generations(infection_size, carrying_capacity, growth_factor)
    lineage_total = growth + steady_generations
    mutated_probability = mutation_carrier_probability(
        mutation_probability, lineage_total
    )
    expected_mutants = carrying_capacity * mutated_probability
    screening_mutants = math.ceil(expected_mutants)
    theta = carrying_capacity * mutation_probability
    richness_lower, richness_upper = ewens_expected_richness_bounds(
        screening_mutants, theta
    )

    warnings: list[str] = []
    if screening_mutants == 0:
        warnings.append(
            "No mutant cells are expected. The legacy V3.1 implementation does "
            "not currently handle the zero-mutant branch."
        )
    if screening_mutants >= crp_cell_warning:
        warnings.append(
            "The current CRP visits approximately one inferred mutant cell per "
            f"iteration ({screening_mutants:,}) and rebuilds probabilities over "
            "all current strains. Use an efficient Ewens sampler; for extremely "
            "large values, use a bounded output sample or aggregate model."
        )
    if richness_upper >= strain_warning:
        warnings.append(
            "The expected Ewens richness may require materializing at least "
            f"about {richness_lower:,.0f}-{richness_upper:,.0f} strains. No data "
            "structure can store a complete graph in sublinear space relative "
            "to the number of strains requested."
        )

    return V31LoadEstimate(
        infection_size=infection_size,
        carrying_capacity=carrying_capacity,
        growth_factor=growth_factor,
        growth_generations=growth,
        steady_generations=steady_generations,
        lineage_generations=lineage_total,
        mutation_probability=mutation_probability,
        mutation_carrier_probability=mutated_probability,
        expected_unmutated_cells=carrying_capacity - expected_mutants,
        expected_mutant_cells=expected_mutants,
        mutant_cells_for_screening=screening_mutants,
        legacy_ewens_theta=theta,
        expected_strain_richness_lower=richness_lower,
        expected_strain_richness_upper=richness_upper,
        warnings=tuple(warnings),
    )
