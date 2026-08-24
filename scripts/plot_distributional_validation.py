#!/usr/bin/env python3
"""Create publication-ready figures for the distributional validation report."""

from __future__ import annotations

import argparse
import math
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.axes import Axes
from matplotlib.figure import Figure
from matplotlib.lines import Line2D
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch
from validate_distributions import (
    _count,
    _count_kernel_features,
    _individual_reference_features,
    _neutral_endpoint_distribution,
)

from trophosome.config import EvolutionConfig, HostConfig
from trophosome.count_model import (
    IdAllocator,
    LineageRecorder,
    PopulationState,
    fixed_pool_migration_step,
    free_living_selection_step,
    merge_populations,
    proportional_rescale,
    sample_population,
    wright_fisher_step,
)
from trophosome.simulation import _reservoir_founders

OBSERVED = "#2C7FB8"
EXPECTED = "#D95F0E"
REFERENCE = "#31A354"
NEUTRAL = "#5F6B76"


def _style() -> None:
    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 10,
            "axes.titlesize": 11,
            "axes.labelsize": 10,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "legend.frameon": False,
            "figure.dpi": 140,
            "savefig.dpi": 300,
        }
    )


def _save(figure: Figure, output: Path, stem: str) -> None:
    output.mkdir(parents=True, exist_ok=True)
    figure.savefig(output / f"{stem}.svg", bbox_inches="tight")
    figure.savefig(output / f"{stem}.png", bbox_inches="tight")
    plt.close(figure)


def _cycle_figure(output: Path) -> None:
    figure, axis = plt.subplots(figsize=(12, 4.8))
    axis.set_xlim(0, 12)
    axis.set_ylim(0, 5)
    axis.axis("off")

    nodes = [
        (0.4, 2.0, 1.7, 1.15, "Effective\nreservoir", "Environment"),
        (2.7, 2.0, 1.7, 1.15, "Founder\nsampling", "Infection"),
        (5.0, 2.0, 2.0, 1.15, "Growth, drift\nand mutation", "Within host"),
        (7.6, 2.0, 1.7, 1.15, "Escape\nsampling", "Host death"),
        (
            9.7,
            2.0,
            2.1,
            1.15,
            "Mixing, capacity\nand optional\nselection",
            "Environment",
        ),
    ]
    for index, (x, y, width, height, label, stage) in enumerate(nodes):
        color = OBSERVED if index in {1, 2, 3} else REFERENCE
        box = FancyBboxPatch(
            (x, y),
            width,
            height,
            boxstyle="round,pad=0.06,rounding_size=0.08",
            facecolor="white",
            edgecolor=color,
            linewidth=2,
        )
        axis.add_patch(box)
        axis.text(
            x + width / 2,
            y + height * 0.57,
            label,
            ha="center",
            va="center",
            weight="bold",
        )
        axis.text(
            x + width / 2,
            y - 0.24,
            stage,
            ha="center",
            va="top",
            color=NEUTRAL,
        )
        if index < len(nodes) - 1:
            next_x = nodes[index + 1][0]
            arrow = FancyArrowPatch(
                (x + width + 0.08, y + height / 2),
                (next_x - 0.08, y + height / 2),
                arrowstyle="-|>",
                mutation_scale=13,
                linewidth=1.5,
                color=NEUTRAL,
            )
            axis.add_patch(arrow)
    return_arrow = FancyArrowPatch(
        (10.75, 3.25),
        (1.25, 3.25),
        connectionstyle="arc3,rad=0.18",
        arrowstyle="-|>",
        mutation_scale=13,
        linewidth=1.5,
        color=NEUTRAL,
    )
    axis.add_patch(return_arrow)

    annotations = [
        (1.25, 0.7, "Reservoir founder\nsampling"),
        (6.0, 0.7, "Within-host reproduction,\nmutation and jackpots"),
        (8.45, 0.08, "Finite escape\nsampling"),
        (10.75, 0.7, "Neutral regulation or\nfree-living selection"),
    ]
    for x, y, text in annotations:
        axis.text(
            x,
            y,
            text,
            ha="center",
            va="center",
            color=NEUTRAL,
            fontsize=9,
        )
    axis.set_title(
        "Where each distributional validation acts in the biological cycle",
        pad=14,
        weight="bold",
    )
    _save(figure, output, "validation-cycle")


def _binomial_pmf(size: int, probability: float) -> tuple[np.ndarray, np.ndarray]:
    x = np.arange(size + 1)
    probability_values = np.asarray(
        [
            math.comb(size, int(value))
            * probability ** int(value)
            * (1.0 - probability) ** (size - int(value))
            for value in x
        ]
    )
    return x, probability_values


def _hypergeometric_pmf(
    population_size: int, focal_count: int, sample_size: int
) -> tuple[np.ndarray, np.ndarray]:
    lower = max(0, sample_size - (population_size - focal_count))
    upper = min(sample_size, focal_count)
    x = np.arange(lower, upper + 1)
    denominator = math.comb(population_size, sample_size)
    probabilities = np.asarray(
        [
            math.comb(focal_count, int(value))
            * math.comb(population_size - focal_count, sample_size - int(value))
            / denominator
            for value in x
        ]
    )
    return x, probabilities


def _histogram(values: np.ndarray, x: np.ndarray) -> np.ndarray:
    return np.asarray([(values == value).mean() for value in x])


def _distribution_panel(
    axis: Axes,
    x: np.ndarray,
    observed: np.ndarray,
    expected: np.ndarray,
    title: str,
    x_label: str,
) -> None:
    axis.bar(
        x,
        observed,
        width=0.82,
        color=OBSERVED,
        alpha=0.68,
        label="Exact-count observations",
    )
    axis.plot(
        x,
        expected,
        color=EXPECTED,
        marker="o",
        markersize=3.2,
        linewidth=1.6,
        label="Expected probability",
    )
    axis.set_title(title, loc="left", weight="bold")
    axis.set_xlabel(x_label)
    axis.set_ylabel("Probability")
    axis.grid(axis="y", color="#D9D9D9", linewidth=0.6, alpha=0.6)


def _analytical_samples(
    repetitions: int, seed: int
) -> list[tuple[str, str, np.ndarray, np.ndarray, np.ndarray]]:
    rng = np.random.default_rng(seed)
    panels: list[tuple[str, str, np.ndarray, np.ndarray, np.ndarray]] = []

    neutral_state = PopulationState.from_counts([30, 70], [1.0, 1.0])
    neutral_evolution = EvolutionConfig(
        mutation_probability=0.0, within_host_selection=False
    )
    neutral = np.empty(repetitions, dtype=np.int64)
    for index in range(repetitions):
        state, _ = wright_fisher_step(
            neutral_state, 80, neutral_evolution, rng, IdAllocator(2)
        )
        neutral[index] = _count(state, 0)
    x, expected = _binomial_pmf(80, 0.3)
    keep = expected > 0.00015
    panels.append(
        (
            "A  Neutral reproduction",
            "Offspring cells of strain A",
            x[keep],
            _histogram(neutral, x[keep]),
            expected[keep],
        )
    )

    selected_state = PopulationState.from_counts([30, 70], [2.0, 1.0])
    selected_evolution = EvolutionConfig(
        mutation_probability=0.0, within_host_selection=True
    )
    selected = np.empty(repetitions, dtype=np.int64)
    for index in range(repetitions):
        state, _ = wright_fisher_step(
            selected_state, 60, selected_evolution, rng, IdAllocator(2)
        )
        selected[index] = _count(state, 0)
    selected_probability = 60 / 130
    x, expected = _binomial_pmf(60, selected_probability)
    keep = expected > 0.00015
    panels.append(
        (
            "B  Fitness-weighted reproduction",
            "Offspring cells of fitter strain A",
            x[keep],
            _histogram(selected, x[keep]),
            expected[keep],
        )
    )

    free_living_state = PopulationState.from_counts([30, 70], [1.0, 2.0], [2.0, 1.0])
    free_living_selected = np.empty(repetitions, dtype=np.int64)
    for index in range(repetitions):
        state = free_living_selection_step(free_living_state, 60, rng)
        free_living_selected[index] = _count(state, 0)
    x, expected = _binomial_pmf(60, selected_probability)
    keep = expected > 0.00015
    panels.append(
        (
            "C  Free-living selection",
            "Environmental descendants of fitter strain A",
            x[keep],
            _histogram(free_living_selected, x[keep]),
            expected[keep],
        )
    )

    mutation_state = PopulationState.from_counts([40, 60], [1.0, 1.0])
    mutation_evolution = EvolutionConfig(
        mutation_probability=0.1,
        mutation_effect_mean=0.0,
        mutation_effect_sd=0.0,
        within_host_selection=False,
        max_materialized_mutants=100,
    )
    mutations = np.empty(repetitions, dtype=np.int64)
    for index in range(repetitions):
        recorder = LineageRecorder.from_founders(mutation_state.genotype_ids)
        _, mutations[index] = wright_fisher_step(
            mutation_state,
            50,
            mutation_evolution,
            rng,
            IdAllocator(2),
            lineage_recorder=recorder,
        )
    x, expected = _binomial_pmf(50, 0.1)
    keep = expected > 0.00015
    panels.append(
        (
            "D  Mutation events",
            "New mutations in one generation",
            x[keep],
            _histogram(mutations, x[keep]),
            expected[keep],
        )
    )

    reservoir = PopulationState.from_counts([20, 30, 50], [1.0, 1.0, 1.0])
    cumulative = np.cumsum(reservoir.counts, dtype=np.float64)
    cumulative /= cumulative[-1]
    founders = np.empty(repetitions, dtype=np.int64)
    for index in range(repetitions):
        sampled = _reservoir_founders(reservoir, cumulative, 8, rng)
        founders[index] = _count(sampled, 0)
    x, expected = _binomial_pmf(8, 0.2)
    panels.append(
        (
            "E  Infection founder sample",
            "Strain A cells among 8 founders",
            x,
            _histogram(founders, x),
            expected,
        )
    )

    finite_state = PopulationState.from_counts([20, 80], [1.0, 1.0])
    escape = np.empty(repetitions, dtype=np.int64)
    for index in range(repetitions):
        sampled = sample_population(finite_state, 25, rng, replace=False)
        escape[index] = _count(sampled, 0)
    x, expected = _hypergeometric_pmf(100, 20, 25)
    panels.append(
        (
            "F  Escape without replacement",
            "Strain A cells among 25 escapees",
            x,
            _histogram(escape, x),
            expected,
        )
    )

    environment = PopulationState.from_counts([30, 70], [1.0, 1.0])
    updated_counts = np.empty(repetitions, dtype=np.int64)
    for index in range(repetitions):
        release = PopulationState.from_counts(
            rng.multinomial(100, [0.3, 0.7]).tolist(), [1.0, 1.0]
        )
        updated = proportional_rescale(
            merge_populations(environment, release), 100, rng
        )
        updated_counts[index] = _count(updated, 0)
    release_x, release_probability = _binomial_pmf(100, 0.3)
    expected_environment = np.zeros(101, dtype=np.float64)
    for returned, probability in zip(
        release_x.tolist(), release_probability.tolist(), strict=True
    ):
        unrounded = (30 + returned) / 2
        lower = math.floor(unrounded)
        upper = math.ceil(unrounded)
        if lower == upper:
            expected_environment[lower] += probability
        else:
            expected_environment[lower] += probability / 2
            expected_environment[upper] += probability / 2
    x = np.flatnonzero(expected_environment > 0.00015)
    panels.append(
        (
            "G  Neutral environmental return",
            "Strain A cells after mixing and capacity",
            x,
            _histogram(updated_counts, x),
            expected_environment[x],
        )
    )
    return panels


def _analytical_figure(output: Path, repetitions: int, seed: int) -> None:
    panels = _analytical_samples(repetitions, seed)
    figure, axes = plt.subplots(2, 4, figsize=(16, 8.2), constrained_layout=True)
    for axis, (title, label, x, observed, expected) in zip(
        axes.flat, panels, strict=False
    ):
        _distribution_panel(axis, x, observed, expected, title, label)
    for axis in axes.flat[len(panels) :]:
        axis.axis("off")
    handles, labels = axes.flat[0].get_legend_handles_labels()
    figure.legend(
        handles,
        labels,
        loc="upper center",
        ncol=2,
        bbox_to_anchor=(0.5, 1.035),
    )
    figure.suptitle(
        "Exact-count observations follow the expected sampling laws "
        f"(n={repetitions:,})",
        y=1.075,
        fontsize=14,
        weight="bold",
    )
    _save(figure, output, "validation-analytical-distributions")


def _migration_figure(output: Path, repetitions: int, seed: int) -> None:
    focal = PopulationState.from_counts([30, 70], [1.0, 1.0])
    regional = PopulationState.from_counts([20, 80], [1.0, 1.0])
    emigrant_values = np.empty(repetitions, dtype=np.int64)
    immigrant_values = np.empty(repetitions, dtype=np.int64)
    emigration_rng = np.random.default_rng(seed)
    immigration_rng = np.random.default_rng(seed + 1)
    for index in range(repetitions):
        _, emigrants, immigrants = fixed_pool_migration_step(
            focal,
            regional,
            20,
            emigration_rng,
            immigration_rng,
        )
        assert emigrants is not None
        assert immigrants is not None
        emigrant_values[index] = _count(emigrants, 0)
        immigrant_values[index] = _count(immigrants, 0)

    emigrant_x, emigrant_expected = _hypergeometric_pmf(100, 30, 20)
    immigrant_x, immigrant_expected = _binomial_pmf(20, 0.2)
    figure, axes = plt.subplots(1, 2, figsize=(11.5, 4.5), constrained_layout=True)
    _distribution_panel(
        axes[0],
        emigrant_x,
        _histogram(emigrant_values, emigrant_x),
        emigrant_expected,
        "A  Focal emigration",
        "Strain A cells among 20 emigrants",
    )
    _distribution_panel(
        axes[1],
        immigrant_x,
        _histogram(immigrant_values, immigrant_x),
        immigrant_expected,
        "B  Immigration from fixed source",
        "Strain A cells among 20 immigrants",
    )
    for axis in axes:
        axis.legend(loc="upper right", fontsize=9)
    figure.suptitle(
        "Fixed-pool exchange follows the declared sampling laws",
        y=1.03,
        fontsize=14,
        weight="bold",
    )
    _save(figure, output, "validation-fixed-pool-migration")


def _dual_fitness_figure(output: Path, repetitions: int, seed: int) -> None:
    state = PopulationState.from_counts([1], [1.0], [1.0])
    evolution = EvolutionConfig(
        mutation_probability=1.0,
        mutation_effect_mean=-0.01,
        mutation_effect_sd=0.02,
        within_host_selection=False,
        max_materialized_mutants=repetitions,
    )
    mutants, _ = wright_fisher_step(
        state,
        repetitions,
        evolution,
        np.random.default_rng(seed),
        IdAllocator(1),
        free_living_fitness_rng=np.random.default_rng(seed + 1),
    )
    host_effects = mutants.within_host_fitness - 1.0
    environment_effects = mutants.free_living_fitness - 1.0
    correlation = float(np.corrcoef(host_effects, environment_effects)[0, 1])

    figure, axes = plt.subplots(1, 2, figsize=(11.5, 4.5), constrained_layout=True)
    bins = np.linspace(-0.09, 0.07, 45)
    axes[0].hist(
        host_effects,
        bins=bins,
        density=True,
        alpha=0.55,
        color=OBSERVED,
        label="Within-host effect",
    )
    axes[0].hist(
        environment_effects,
        bins=bins,
        density=True,
        alpha=0.45,
        color=REFERENCE,
        label="Free-living effect",
    )
    x = np.linspace(bins[0], bins[-1], 300)
    normal = np.exp(-0.5 * ((x + 0.01) / 0.02) ** 2) / (0.02 * math.sqrt(2 * math.pi))
    axes[0].plot(x, normal, color=EXPECTED, linewidth=2, label="Configured normal")
    axes[0].set_title("A  Both habitats use the same effect distribution", loc="left")
    axes[0].set_xlabel("Mutant fitness minus parental fitness")
    axes[0].set_ylabel("Density")
    axes[0].legend()

    plot_count = min(3_000, repetitions)
    axes[1].scatter(
        host_effects[:plot_count],
        environment_effects[:plot_count],
        s=9,
        alpha=0.28,
        color=OBSERVED,
        edgecolors="none",
    )
    axes[1].axhline(-0.01, color=NEUTRAL, linewidth=1, linestyle="--")
    axes[1].axvline(-0.01, color=NEUTRAL, linewidth=1, linestyle="--")
    axes[1].text(
        0.03,
        0.95,
        f"correlation = {correlation:.3f}",
        transform=axes[1].transAxes,
        va="top",
    )
    axes[1].set_title("B  Effects are sampled independently", loc="left")
    axes[1].set_xlabel("Within-host fitness effect")
    axes[1].set_ylabel("Free-living fitness effect")
    figure.suptitle(
        f"Independent dual-habitat fitness effects (n={repetitions:,} mutants)",
        fontsize=14,
        weight="bold",
    )
    _save(figure, output, "validation-dual-fitness-effects")


def _trajectory_figure(output: Path, repetitions: int, seed: int) -> None:
    figure, axes = plt.subplots(2, 2, figsize=(11.5, 8.5), constrained_layout=True)

    state = PopulationState.from_counts([3, 7], [1.0, 1.0])
    evolution = EvolutionConfig(mutation_probability=0.0, within_host_selection=False)
    expected = _neutral_endpoint_distribution(3, 10, (10, 10, 10))
    rng = np.random.default_rng(seed)
    observed_values = np.empty(repetitions, dtype=np.int64)
    for replicate in range(repetitions):
        current = state
        for target in (10, 10, 10):
            current, _ = wright_fisher_step(
                current, target, evolution, rng, IdAllocator(2)
            )
        observed_values[replicate] = _count(current, 0)
    x = np.arange(11)
    _distribution_panel(
        axes[0, 0],
        x,
        _histogram(observed_values, x),
        expected,
        "A  Drift accumulated across three generations",
        "Strain A cells at the endpoint",
    )

    host = HostConfig(
        population_size=1,
        infection_bottleneck=2,
        carrying_capacity=16,
        growth_factor=2.0,
        steady_generations=2,
        host_generations=1,
        escape_fraction=0.0,
    )
    reference_repetitions = max(2_000, repetitions // 2)
    count_features = _count_kernel_features(host, 0.04, reference_repetitions, seed + 1)
    cell_features = _individual_reference_features(
        host, 0.04, reference_repetitions, seed + 2
    )
    feature_details = [
        (0, "B  Adult strain richness", "Number of strains"),
        (1, "C  Founder-lineage abundance", "Cells retaining founder strain"),
        (2, "D  Largest mutant jackpot clone", "Cells in largest mutant clone"),
    ]
    for axis, (column, title, label) in zip(
        axes.flat[1:], feature_details, strict=True
    ):
        lower = int(
            min(count_features[:, column].min(), cell_features[:, column].min())
        )
        upper = int(
            max(count_features[:, column].max(), cell_features[:, column].max())
        )
        x = np.arange(lower, upper + 1)
        axis.plot(
            x,
            _histogram(count_features[:, column], x),
            color=OBSERVED,
            marker="o",
            markersize=3.5,
            linewidth=1.6,
            label="Exact-count model",
        )
        axis.plot(
            x,
            _histogram(cell_features[:, column], x),
            color=REFERENCE,
            marker="s",
            markersize=3.2,
            linewidth=1.5,
            label="Cell-by-cell reference",
        )
        axis.set_title(title, loc="left", weight="bold")
        axis.set_xlabel(label)
        axis.set_ylabel("Probability")
        axis.grid(axis="y", color="#D9D9D9", linewidth=0.6, alpha=0.6)
    handles = [
        Line2D(
            [0],
            [0],
            color=EXPECTED,
            marker="o",
            markersize=4,
            linewidth=1.6,
            label="Analytical expectation",
        ),
        Line2D(
            [0],
            [0],
            color=OBSERVED,
            marker="o",
            markersize=4,
            linewidth=1.6,
            label="Exact-count model",
        ),
        Line2D(
            [0],
            [0],
            color=REFERENCE,
            marker="s",
            markersize=4,
            linewidth=1.5,
            label="Cell-by-cell reference",
        ),
    ]
    figure.legend(
        handles=handles,
        loc="upper center",
        ncol=3,
        bbox_to_anchor=(0.5, 1.02),
    )
    figure.suptitle(
        "Multi-generation dynamics and early-mutation jackpots match "
        "independent references",
        y=1.06,
        fontsize=14,
        weight="bold",
    )
    _save(figure, output, "validation-trajectories-and-jackpots")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repetitions", type=int, default=30_000)
    parser.add_argument("--seed", type=int, default=20260810)
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("docs/figures"),
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    if args.repetitions < 2_000:
        raise ValueError("at least 2,000 repetitions are required")
    _style()
    _cycle_figure(args.output)
    _analytical_figure(args.output, args.repetitions, args.seed)
    _migration_figure(args.output, args.repetitions, args.seed + 9)
    _dual_fitness_figure(args.output, args.repetitions, args.seed + 3)
    _trajectory_figure(args.output, args.repetitions, args.seed + 20)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
