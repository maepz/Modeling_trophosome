#!/usr/bin/env python3
"""Summarise and plot the completed Phase 1 first pilot.

The script reads the machine-local scratch location from ``layout.local.json``.
It never modifies raw simulation results. Derived TSV/JSON files are written
beside this script and report figures are written below ``docs/figures``.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from collections import defaultdict
from collections.abc import Iterable
from pathlib import Path
from typing import Any

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch  # noqa: E402

INITIAL_RICHNESS = 100
INITIAL_IDS = tuple(range(INITIAL_RICHNESS))
CORE_CELL_ORDER = tuple(f"p01-s01-c{index:04d}" for index in range(1, 13))
FULL_CELL_ORDER = tuple(f"p01-s01-c{index:04d}" for index in range(1, 21))
CELL_ORDER = CORE_CELL_ORDER
SEED_ORDER = ("sb0001", "sb0002", "sb0003")
GROUP_COLOURS = {
    "NR": "#5F6B76",
    "NRM": "#009E73",
    "MR0": "#0072B2",
    "MUT": "#D55E00",
    "FB0": "#8C6BB1",
    "MRM": "#CC79A7",
    "MRLOW": "#009E73",
    "FBA": "#6A51A3",
}
CELL_MARKERS = (
    "o",
    "s",
    "^",
    "D",
    "P",
    "v",
    "<",
    ">",
    "X",
    "h",
    "*",
    "d",
    "p",
    "8",
    "H",
    "+",
    "x",
)
SEED_OFFSETS = (-0.18, 0.0, 0.18)


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, Any]], fields: Iterable[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=list(fields),
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(rows)


def cell_short(cell_id: str) -> str:
    return cell_id.rsplit("-", 1)[-1]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def diversity_metrics(counts: dict[int, int]) -> dict[str, float | int]:
    positive = np.asarray(
        [value for value in counts.values() if value > 0], dtype=float
    )
    total = float(positive.sum())
    frequencies = positive / total
    shannon = float(-np.sum(frequencies * np.log(frequencies)))
    richness = int(len(positive))
    return {
        "richness": richness,
        "shannon": shannon,
        "simpson": float(1.0 - np.sum(frequencies**2)),
        "hill_q1": float(math.exp(shannon)),
        "hill_q2": float(1.0 / np.sum(frequencies**2)),
        "evenness": float(shannon / math.log(richness)) if richness > 1 else math.nan,
        "gene_diversity": float(1.0 - np.sum(frequencies**2)),
        "top_frequency": float(frequencies.max()),
    }


def aligned_frequencies(
    first: dict[int, int], second: dict[int, int]
) -> tuple[np.ndarray, np.ndarray]:
    keys = sorted(set(first) | set(second))
    first_total = sum(first.values())
    second_total = sum(second.values())
    p = np.asarray([first.get(key, 0) / first_total for key in keys], dtype=float)
    q = np.asarray([second.get(key, 0) / second_total for key in keys], dtype=float)
    return p, q


def composition_metrics(
    first: dict[int, int], second: dict[int, int]
) -> dict[str, float]:
    p, q = aligned_frequencies(first, second)
    midpoint = 0.5 * (p + q)

    def kl(values: np.ndarray, reference: np.ndarray) -> float:
        keep = values > 0
        return float(np.sum(values[keep] * np.log(values[keep] / reference[keep])))

    return {
        "total_variation": float(0.5 * np.abs(p - q).sum()),
        "jensen_shannon": float(0.5 * kl(p, midpoint) + 0.5 * kl(q, midpoint)),
    }


def read_environment(path: Path) -> dict[int, int]:
    result: dict[int, int] = {}
    with path.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            result[int(row["strain_id"])] = int(row["count"])
    return result


def reversed_binary_lines(path: Path, block_size: int = 1024 * 1024) -> Iterable[bytes]:
    """Yield non-empty lines from a file in reverse without loading it into RAM."""
    with path.open("rb") as handle:
        handle.seek(0, 2)
        position = handle.tell()
        remainder = b""
        while position > 0:
            amount = min(block_size, position)
            position -= amount
            handle.seek(position)
            remainder = handle.read(amount) + remainder
            lines = remainder.split(b"\n")
            remainder = lines[0]
            for line in reversed(lines[1:]):
                if line:
                    yield line
        if remainder:
            yield remainder


def read_lineages(
    path: Path, final_ids: Iterable[int]
) -> tuple[dict[int, int], int, int]:
    """Summarize lineages and retain only ancestry needed by the final population.

    Mutation-enabled high-host cells can produce tens of millions of lineage rows.
    A forward pass computes event totals and the largest one-transition mutation
    count. A reverse pass follows only ancestors of lineages that remain in the
    final environment, keeping memory use proportional to retained ancestry rather
    than to every transient mutation.
    """
    generated_lineages = 0
    maximum = 0
    current_transition: tuple[int, int, int] | None = None
    current_transition_count = 0
    with path.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            generated_lineages += 1
            transition = (
                int(row["generation"]),
                int(row["host_id"]),
                int(row["within_host_generation"]),
            )
            if transition == current_transition:
                current_transition_count += 1
            else:
                maximum = max(maximum, current_transition_count)
                current_transition = transition
                current_transition_count = 1
    maximum = max(maximum, current_transition_count)

    needed = {
        int(strain_id) for strain_id in final_ids if int(strain_id) >= INITIAL_RICHNESS
    }
    parents: dict[int, int] = {}
    if needed:
        with path.open("rb") as handle:
            header = handle.readline().decode("utf-8").rstrip("\r\n").split(",")
        child_index = header.index("strain_id")
        parent_index = header.index("parent_strain_id")
        for raw_line in reversed_binary_lines(path):
            fields = raw_line.split(b",")
            try:
                child = int(fields[child_index])
            except (IndexError, ValueError):
                continue
            if child not in needed:
                continue
            parent = int(fields[parent_index])
            parents[child] = parent
            if parent >= INITIAL_RICHNESS:
                needed.add(parent)
    return parents, generated_lineages, maximum


def trace_roots(parents: dict[int, int], final_ids: Iterable[int]) -> dict[int, int]:
    cache = {strain_id: strain_id for strain_id in INITIAL_IDS}

    def root(strain_id: int) -> int:
        if strain_id in cache:
            return cache[strain_id]
        path: list[int] = []
        current = strain_id
        while current not in cache:
            path.append(current)
            if current not in parents:
                raise RuntimeError(f"lineage parent is missing for strain {current}")
            current = parents[current]
        result = cache[current]
        for descendant in path:
            cache[descendant] = result
        return result

    return {strain_id: root(strain_id) for strain_id in final_ids}


def unique_mutant_ids(path: Path) -> set[int]:
    result: set[int] = set()
    with path.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            strain_id = int(row["strain_id"])
            if strain_id >= INITIAL_RICHNESS:
                result.add(strain_id)
    return result


def pooling_metrics(output: Path, host_number: int) -> dict[str, float | int]:
    pooled_by_generation: dict[int, list[tuple[int, int]]] = defaultdict(list)
    with (output / "pooled_host_counts_and_occupancy.csv").open(
        newline="", encoding="utf-8"
    ) as handle:
        for row in csv.DictReader(handle):
            pooled_by_generation[int(row["generation"])].append(
                (int(row["strain_id"]), int(row["occupied_hosts"]))
            )

    infection_ids: dict[int, set[int]] = defaultdict(set)
    with (output / "infection_counts.csv").open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            infection_ids[int(row["generation"])].add(int(row["strain_id"]))

    release_ids: dict[int, set[int]] = defaultdict(set)
    with (output / "release_counts.csv").open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            release_ids[int(row["generation"])].add(int(row["strain_id"]))

    final_generation = 5
    pooled = pooled_by_generation.get(final_generation, [])
    original_occupancy_hosts = [
        occupied for strain_id, occupied in pooled if strain_id < INITIAL_RICHNESS
    ]
    return {
        "infection_richness_g5": len(infection_ids.get(final_generation, set())),
        "pooled_adult_richness_g5": len(pooled),
        "release_richness_g5": len(release_ids.get(final_generation, set())),
        "median_original_occupancy_hosts_g5": (
            float(np.median(original_occupancy_hosts))
            if original_occupancy_hosts
            else 0.0
        ),
        "median_original_occupancy_fraction_g5": (
            float(np.median(original_occupancy_hosts)) / host_number
            if original_occupancy_hosts
            else 0.0
        ),
        "original_strains_in_pooled_adults_g5": sum(
            strain_id < INITIAL_RICHNESS for strain_id, _occupied in pooled
        ),
    }


def describe(values: list[float]) -> dict[str, float]:
    array = np.asarray(values, dtype=float)
    return {
        "mean": float(np.mean(array)),
        "median": float(np.median(array)),
        "minimum": float(np.min(array)),
        "maximum": float(np.max(array)),
        "sd": float(np.std(array, ddof=1)) if len(array) > 1 else 0.0,
    }


def configure_plot_style() -> None:
    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 9.5,
            "axes.titlesize": 11,
            "axes.labelsize": 10,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "axes.grid": True,
            "axes.grid.axis": "y",
            "grid.alpha": 0.2,
            "legend.frameon": False,
            "figure.facecolor": "white",
            "savefig.facecolor": "white",
        }
    )


def save_figure(figure: plt.Figure, figure_dir: Path, stem: str) -> None:
    figure_dir.mkdir(parents=True, exist_ok=True)
    figure.savefig(figure_dir / f"{stem}.png", dpi=220, bbox_inches="tight")
    figure.savefig(figure_dir / f"{stem}.svg", bbox_inches="tight")
    plt.close(figure)


def plot_design(cells: dict[str, dict[str, Any]], figure_dir: Path) -> None:
    figure = plt.figure(figsize=(13, 8.6))
    grid = figure.add_gridspec(2, 1, height_ratios=(1.0, 1.8), hspace=0.30)
    life = figure.add_subplot(grid[0])
    life.set_axis_off()
    labels = (
        "Environmental\nreservoir",
        "Infection\n(B = 10)",
        "Within-host growth\n(K = 10⁹; 500 generations)",
        "Escape and\npooling",
        "Hamilton regulation\n(Nₑ = 10⁹)",
    )
    positions = np.linspace(0.08, 0.92, len(labels))
    for index, (x_value, label) in enumerate(zip(positions, labels, strict=True)):
        box = FancyBboxPatch(
            (x_value - 0.082, 0.38),
            0.164,
            0.28,
            boxstyle="round,pad=0.02,rounding_size=0.02",
            facecolor="#EAF3F8",
            edgecolor="#2C7FB8",
            linewidth=1.5,
            transform=life.transAxes,
        )
        life.add_patch(box)
        life.text(
            x_value, 0.52, label, ha="center", va="center", transform=life.transAxes
        )
        if index < len(labels) - 1:
            life.add_patch(
                FancyArrowPatch(
                    (x_value + 0.085, 0.52),
                    (positions[index + 1] - 0.085, 0.52),
                    arrowstyle="-|>",
                    mutation_scale=12,
                    color="#5F6B76",
                    transform=life.transAxes,
                )
            )
    life.add_patch(
        FancyArrowPatch(
            (0.92, 0.31),
            (0.08, 0.31),
            connectionstyle="arc3,rad=-0.10",
            arrowstyle="-|>",
            mutation_scale=12,
            color="#5F6B76",
            transform=life.transAxes,
        )
    )
    life.text(
        0.50,
        0.16,
        "Pooled escapees return to the environmental reservoir",
        ha="center",
        color="#5F6B76",
        transform=life.transAxes,
    )
    life.text(
        0.01,
        0.92,
        "A  Modelled host–environment cycle",
        weight="bold",
        transform=life.transAxes,
    )

    design = figure.add_subplot(grid[1])
    design.set_axis_off()
    rows = [
        (
            "No return",
            "Tests whether within-host events can affect a reservoir without return",
            ["p01-s01-c0001", "p01-s01-c0002"],
        ),
        (
            "Matched return",
            "R = 10⁹ and α = 0.5; H varies from 10² to 10⁵",
            ["p01-s01-c0003", "p01-s01-c0004", "p01-s01-c0005", "p01-s01-c0006"],
        ),
        (
            "Mutation bracket",
            "H = 100 and α = 0.5; u varies from 0 to 10⁻⁹",
            [
                "p01-s01-c0003",
                "p01-s01-c0007",
                "p01-s01-c0008",
                "p01-s01-c0009",
                "p01-s01-c0010",
            ],
        ),
        (
            "Feedback bracket",
            "Mutation off; α spans 0.001, 0.5 and 0.909",
            ["p01-s01-c0011", "p01-s01-c0003", "p01-s01-c0012"],
        ),
    ]
    y_values = (0.79, 0.56, 0.33, 0.10)
    for (title, subtitle, identifiers), y_value in zip(rows, y_values, strict=True):
        design.text(
            0.01,
            y_value + 0.115,
            title,
            weight="bold",
            va="center",
            transform=design.transAxes,
        )
        design.text(
            0.16,
            y_value + 0.115,
            subtitle,
            color="#5F6B76",
            va="center",
            transform=design.transAxes,
        )
        start = 0.17
        step = 0.15
        for index, cell_id in enumerate(identifiers):
            cell = cells[cell_id]
            x_value = start + index * step
            colour = GROUP_COLOURS[cell["experimental_group"]]
            design.add_patch(
                FancyBboxPatch(
                    (x_value, y_value - 0.055),
                    0.12,
                    0.11,
                    boxstyle="round,pad=0.012,rounding_size=0.018",
                    facecolor=colour,
                    edgecolor="white",
                    transform=design.transAxes,
                )
            )
            detail = cell_short(cell_id)
            if title == "Matched return":
                detail += f"\nH={cell['H']:.0e}"
            elif title == "Mutation bracket":
                detail += f"\nu={cell['u']:.0e}" if cell["u"] else "\nu=0"
            elif title == "Feedback bracket":
                detail += f"\nα={cell['alpha']:.3g}"
            elif cell_id.endswith("c0002"):
                detail += "\nu=10⁻¹⁰"
            else:
                detail += "\nu=0"
            design.text(
                x_value + 0.06,
                y_value,
                detail,
                ha="center",
                va="center",
                color="white",
                weight="bold",
                fontsize=8.5,
                transform=design.transAxes,
            )
    design.text(
        0.01,
        0.98,
        "B  Twelve-cell first-pilot design",
        weight="bold",
        transform=design.transAxes,
    )
    figure.suptitle(
        "Phase 1 first pilot: biological cycle and planned comparisons",
        fontsize=15,
        weight="bold",
    )
    save_figure(figure, figure_dir, "fig01-pilot-design")


def plot_feasibility(
    endpoints: list[dict[str, Any]], cells: dict[str, dict[str, Any]], figure_dir: Path
) -> None:
    figure, axes = plt.subplots(1, 3, figsize=(14, 4.6))
    measures = (
        ("elapsed_minutes", "Runtime per population (min)", True),
        ("output_mib", "Output per population (MiB)", True),
        ("peak_rss_mib", "Peak process-tree memory (MiB)", False),
    )
    for axis, (measure, label, log_y) in zip(axes, measures, strict=True):
        for index, cell_id in enumerate(CELL_ORDER):
            rows = [row for row in endpoints if row["cell_id"] == cell_id]
            colour = GROUP_COLOURS[cells[cell_id]["experimental_group"]]
            for seed_index, row in enumerate(rows):
                axis.scatter(
                    cells[cell_id]["H"],
                    row[measure],
                    s=42,
                    marker=CELL_MARKERS[index % len(CELL_MARKERS)],
                    color=colour,
                    alpha=0.82,
                    edgecolor="white",
                    linewidth=0.5,
                    label=cell_short(cell_id) if seed_index == 0 else None,
                )
        axis.set_xscale("log")
        if log_y:
            axis.set_yscale("log")
        axis.set_xlabel("Host number, H")
        axis.set_ylabel(label)
    axes[0].text(
        -0.2, 1.08, "A", transform=axes[0].transAxes, weight="bold", fontsize=13
    )
    axes[1].text(
        -0.2, 1.08, "B", transform=axes[1].transAxes, weight="bold", fontsize=13
    )
    axes[2].text(
        -0.2, 1.08, "C", transform=axes[2].transAxes, weight="bold", fontsize=13
    )
    handles, labels = axes[2].get_legend_handles_labels()
    figure.legend(
        handles,
        labels,
        ncol=6,
        loc="lower center",
        bbox_to_anchor=(0.5, -0.04),
        fontsize=8,
    )
    figure.suptitle(
        f"Computational feasibility across all {len(endpoints)} simulated populations",
        fontsize=14,
        weight="bold",
    )
    figure.tight_layout(rect=(0, 0.08, 1, 0.94))
    save_figure(figure, figure_dir, "fig02-computational-feasibility")


def plot_trajectories(
    trajectories: list[dict[str, Any]],
    cells: dict[str, dict[str, Any]],
    figure_dir: Path,
) -> None:
    themes = (
        (
            "Controls and feedback",
            (
                "p01-s01-c0001",
                "p01-s01-c0002",
                "p01-s01-c0011",
                "p01-s01-c0003",
                "p01-s01-c0012",
            ),
        ),
        (
            "Matched total return",
            ("p01-s01-c0003", "p01-s01-c0004", "p01-s01-c0005", "p01-s01-c0006"),
        ),
        (
            "Mutation bracket",
            (
                "p01-s01-c0003",
                "p01-s01-c0007",
                "p01-s01-c0008",
                "p01-s01-c0009",
                "p01-s01-c0010",
            ),
        ),
    )
    palettes = (
        ("#5F6B76", "#009E73", "#9467BD", "#0072B2", "#CC79A7"),
        ("#08306B", "#2171B5", "#4292C6", "#9ECAE1"),
        ("#0072B2", "#F6C85F", "#F28E2B", "#E66101", "#A63603"),
    )
    figure, axes = plt.subplots(2, 3, figsize=(15, 8), sharex=True)
    for column, (title, identifiers) in enumerate(themes):
        for cell_id, colour in zip(identifiers, palettes[column], strict=True):
            cell_rows = [row for row in trajectories if row["cell_id"] == cell_id]
            for seed_id in SEED_ORDER:
                seed_rows = sorted(
                    [row for row in cell_rows if row["seed_block_id"] == seed_id],
                    key=lambda item: item["generation"],
                )
                generations = [row["generation"] for row in seed_rows]
                axes[0, column].plot(
                    generations,
                    [row["environment_richness"] for row in seed_rows],
                    color=colour,
                    alpha=0.22,
                    linewidth=1,
                )
                axes[1, column].plot(
                    generations,
                    [row["environment_gene_diversity"] for row in seed_rows],
                    color=colour,
                    alpha=0.22,
                    linewidth=1,
                )
            generations = sorted({row["generation"] for row in cell_rows})
            med_richness = [
                float(
                    np.median(
                        [
                            row["environment_richness"]
                            for row in cell_rows
                            if row["generation"] == generation
                        ]
                    )
                )
                for generation in generations
            ]
            med_gene = [
                float(
                    np.median(
                        [
                            row["environment_gene_diversity"]
                            for row in cell_rows
                            if row["generation"] == generation
                        ]
                    )
                )
                for generation in generations
            ]
            axes[0, column].plot(
                generations,
                med_richness,
                color=colour,
                linewidth=2.2,
                label=cell_short(cell_id),
            )
            axes[1, column].plot(
                generations,
                med_gene,
                color=colour,
                linewidth=2.2,
                label=cell_short(cell_id),
            )
        axes[0, column].set_title(title, weight="bold")
        axes[0, column].legend(fontsize=8, ncol=2)
        axes[1, column].set_xlabel("Host passage")
        axes[1, column].set_xticks(range(0, 6))
    axes[0, 0].set_ylabel("Environmental richness")
    axes[1, 0].set_ylabel("Environmental gene diversity")
    axes[0, 0].text(
        -0.22, 1.08, "A", transform=axes[0, 0].transAxes, weight="bold", fontsize=13
    )
    axes[1, 0].text(
        -0.22, 1.08, "B", transform=axes[1, 0].transAxes, weight="bold", fontsize=13
    )
    figure.suptitle(
        "Short-term environmental trajectories over five host passages",
        fontsize=14,
        weight="bold",
    )
    figure.tight_layout(rect=(0, 0, 1, 0.95))
    save_figure(figure, figure_dir, "fig03-environmental-trajectories")


def plot_endpoints(
    endpoints: list[dict[str, Any]], cells: dict[str, dict[str, Any]], figure_dir: Path
) -> None:
    figure, axes = plt.subplots(2, 3, figsize=(15, 8.5))
    measures = (
        ("richness_change_pct", "Richness change (%)", (-5, 5)),
        ("hill_q1_change_pct", "D₁ change (%)", (-5, 5)),
        ("hill_q2_change_pct", "D₂ change (%)", (-5, 5)),
        ("evenness_change", "Evenness change", (-0.02, 0.02)),
        ("total_variation", "Total-variation distance", (0, 0.05)),
    )
    labels = [cell_short(cell_id) for cell_id in CELL_ORDER]
    for axis, (measure, ylabel, margin) in zip(axes.flat, measures, strict=False):
        for index, cell_id in enumerate(CELL_ORDER):
            rows = [row for row in endpoints if row["cell_id"] == cell_id]
            colour = GROUP_COLOURS[cells[cell_id]["experimental_group"]]
            values = [row[measure] for row in rows]
            for seed_index, value in enumerate(values):
                axis.scatter(
                    index + SEED_OFFSETS[seed_index],
                    value,
                    color=colour,
                    alpha=0.72,
                    s=30,
                    edgecolor="white",
                    linewidth=0.4,
                )
            axis.scatter(
                index,
                np.median(values),
                marker="_",
                s=190,
                linewidth=2.2,
                color="#172B3A",
                zorder=5,
            )
        axis.axhline(margin[0], color="#7F8C8D", linestyle="--", linewidth=1)
        if margin[1] != margin[0]:
            axis.axhline(margin[1], color="#7F8C8D", linestyle="--", linewidth=1)
        axis.set_ylabel(ylabel)
        axis.set_xticks(range(len(labels)))
        axis.set_xticklabels(labels, rotation=55, ha="right")
    axes[1, 2].set_axis_off()
    axes[1, 2].text(
        0.02,
        0.92,
        "Points: individual seed blocks\n"
        "Black bars: cell medians\n"
        "Dashed lines: agreed biological\n"
        "relevance margins",
        va="top",
        fontsize=11,
        transform=axes[1, 2].transAxes,
    )
    for label, axis in zip("ABCDE", list(axes.flat)[:5], strict=True):
        axis.text(
            -0.17, 1.07, label, transform=axis.transAxes, weight="bold", fontsize=13
        )
    figure.suptitle(
        "Endpoint environmental effects relative to the frozen starting reservoir",
        fontsize=14,
        weight="bold",
    )
    figure.tight_layout(rect=(0, 0, 1, 0.95))
    save_figure(figure, figure_dir, "fig04-endpoint-effects")


def plot_matched_return(
    endpoints: list[dict[str, Any]], cells: dict[str, dict[str, Any]], figure_dir: Path
) -> None:
    identifiers = ("p01-s01-c0003", "p01-s01-c0004", "p01-s01-c0005", "p01-s01-c0006")
    measures = (
        ("original_strains_lost", "Original strains lost"),
        ("root_total_variation", "Root-collapsed composition distance"),
        ("release_richness_g5", "Release-pool richness at passage 5"),
        (
            "median_original_occupancy_hosts_g5",
            "Median hosts occupied per original strain",
        ),
    )
    figure, axes = plt.subplots(2, 2, figsize=(11, 8))
    for axis, (measure, ylabel) in zip(axes.flat, measures, strict=True):
        for seed_id in SEED_ORDER:
            values = []
            hosts = []
            for cell_id in identifiers:
                row = next(
                    item
                    for item in endpoints
                    if item["cell_id"] == cell_id and item["seed_block_id"] == seed_id
                )
                hosts.append(cells[cell_id]["H"])
                values.append(row[measure])
            axis.plot(hosts, values, color="#0072B2", alpha=0.28, linewidth=1.2)
            axis.scatter(hosts, values, color="#0072B2", alpha=0.65, s=34)
        medians = [
            np.median([row[measure] for row in endpoints if row["cell_id"] == cell_id])
            for cell_id in identifiers
        ]
        axis.plot(
            [cells[cell_id]["H"] for cell_id in identifiers],
            medians,
            color="#172B3A",
            linewidth=2.3,
            marker="o",
        )
        axis.set_xscale("log")
        axis.set_xlabel("Host number, H (R = 10⁹ and α = 0.5)")
        axis.set_ylabel(ylabel)
    for label, axis in zip("ABCD", axes.flat, strict=True):
        axis.text(
            -0.18, 1.07, label, transform=axis.transAxes, weight="bold", fontsize=13
        )
    figure.suptitle(
        "Does pooling across more hosts matter when total return is fixed?",
        fontsize=14,
        weight="bold",
    )
    figure.tight_layout(rect=(0, 0, 1, 0.95))
    save_figure(figure, figure_dir, "fig05-matched-return-pooling")


def plot_mutation_bracket(
    endpoints: list[dict[str, Any]], cells: dict[str, dict[str, Any]], figure_dir: Path
) -> None:
    identifiers = (
        "p01-s01-c0003",
        "p01-s01-c0007",
        "p01-s01-c0008",
        "p01-s01-c0009",
        "p01-s01-c0010",
    )
    x_values = np.arange(len(identifiers))
    x_labels = ("0", "10⁻¹²", "10⁻¹¹", "10⁻¹⁰", "10⁻⁹")
    measures = (
        (
            "final_mutant_richness",
            "Mutant strains in final environment (+1; log scale)",
        ),
        ("mutant_abundance_fraction", "Final abundance derived from mutants"),
        ("richness", "Final labelled richness"),
        ("evenness", "Final Pielou evenness"),
    )
    figure, axes = plt.subplots(2, 2, figsize=(11, 8))
    for axis, (measure, ylabel) in zip(axes.flat, measures, strict=True):
        for seed_id in SEED_ORDER:
            values = [
                next(
                    row
                    for row in endpoints
                    if row["cell_id"] == cell_id and row["seed_block_id"] == seed_id
                )[measure]
                for cell_id in identifiers
            ]
            plotted_values = (
                np.asarray(values, dtype=float) + 1
                if measure == "final_mutant_richness"
                else values
            )
            axis.plot(
                x_values,
                plotted_values,
                color="#D55E00",
                alpha=0.25,
                linewidth=1.2,
            )
            axis.scatter(x_values, plotted_values, color="#D55E00", alpha=0.65, s=34)
        medians = [
            np.median([row[measure] for row in endpoints if row["cell_id"] == cell_id])
            for cell_id in identifiers
        ]
        plotted_medians = (
            np.asarray(medians, dtype=float) + 1
            if measure == "final_mutant_richness"
            else medians
        )
        axis.plot(
            x_values,
            plotted_medians,
            color="#172B3A",
            linewidth=2.3,
            marker="o",
        )
        axis.set_xticks(x_values)
        axis.set_xticklabels(x_labels)
        axis.set_xlabel("Whole-genome mutation probability, u")
        axis.set_ylabel(ylabel)
    axes[0, 0].set_yscale("log")
    for label, axis in zip("ABCD", axes.flat, strict=True):
        axis.text(
            -0.18, 1.07, label, transform=axis.transAxes, weight="bold", fontsize=13
        )
    figure.suptitle(
        "Mutation supply spans undetectable to strong environmental novelty",
        fontsize=14,
        weight="bold",
    )
    figure.tight_layout(rect=(0, 0, 1, 0.95))
    save_figure(figure, figure_dir, "fig06-mutation-bracket")


def plot_mutant_fate(endpoints: list[dict[str, Any]], figure_dir: Path) -> None:
    identifiers = (
        "p01-s01-c0002",
        "p01-s01-c0007",
        "p01-s01-c0008",
        "p01-s01-c0009",
        "p01-s01-c0010",
    )
    labels = {
        "p01-s01-c0002": "no return; u=10⁻¹⁰",
        "p01-s01-c0007": "u=10⁻¹²",
        "p01-s01-c0008": "u=10⁻¹¹",
        "p01-s01-c0009": "u=10⁻¹⁰",
        "p01-s01-c0010": "u=10⁻⁹",
    }
    stages = (
        "generated_mutant_lineages",
        "adult_mutant_lineages",
        "released_mutant_lineages",
        "final_mutant_richness",
    )
    stage_labels = ("Generated", "Reached adult", "Released", "Final environment")
    colours = ("#009E73", "#E69F00", "#56B4E9", "#D55E00", "#CC79A7")
    figure, axes = plt.subplots(1, 2, figsize=(13, 5.2))
    for cell_id, colour in zip(identifiers, colours, strict=True):
        medians = []
        for stage in stages:
            medians.append(
                np.median(
                    [row[stage] for row in endpoints if row["cell_id"] == cell_id]
                )
            )
        axes[0].plot(
            range(len(stages)),
            np.asarray(medians) + 1,
            marker="o",
            linewidth=2,
            color=colour,
            label=labels[cell_id],
        )
        for seed_id in SEED_ORDER:
            row = next(
                item
                for item in endpoints
                if item["cell_id"] == cell_id and item["seed_block_id"] == seed_id
            )
            axes[0].plot(
                range(len(stages)),
                np.asarray([row[stage] for stage in stages]) + 1,
                color=colour,
                alpha=0.13,
                linewidth=0.8,
            )
    axes[0].set_yscale("log")
    axes[0].set_xticks(range(len(stages)))
    axes[0].set_xticklabels(stage_labels, rotation=20, ha="right")
    axes[0].set_ylabel("Lineages + 1 (log scale)")
    axes[0].legend(fontsize=8)

    mutation_identifiers = identifiers[1:]
    x_values = np.arange(len(mutation_identifiers))
    x_labels = ("10⁻¹²", "10⁻¹¹", "10⁻¹⁰", "10⁻⁹")
    for seed_id in SEED_ORDER:
        labelled = []
        root = []
        for cell_id in mutation_identifiers:
            row = next(
                item
                for item in endpoints
                if item["cell_id"] == cell_id and item["seed_block_id"] == seed_id
            )
            labelled.append(row["richness"])
            root.append(row["root_richness"])
        axes[1].plot(x_values, labelled, color="#D55E00", alpha=0.2)
        axes[1].plot(x_values, root, color="#0072B2", alpha=0.2)
    labelled_median = [
        np.median([row["richness"] for row in endpoints if row["cell_id"] == cell_id])
        for cell_id in mutation_identifiers
    ]
    root_median = [
        np.median(
            [row["root_richness"] for row in endpoints if row["cell_id"] == cell_id]
        )
        for cell_id in mutation_identifiers
    ]
    axes[1].plot(
        x_values,
        labelled_median,
        color="#D55E00",
        marker="o",
        linewidth=2.3,
        label="Labelled strains",
    )
    axes[1].plot(
        x_values,
        root_median,
        color="#0072B2",
        marker="o",
        linewidth=2.3,
        label="Initial roots",
    )
    axes[1].set_xticks(x_values)
    axes[1].set_xticklabels(x_labels)
    axes[1].set_xlabel("Whole-genome mutation probability, u")
    axes[1].set_ylabel("Final richness")
    axes[1].legend()
    axes[0].text(
        -0.14, 1.07, "A", transform=axes[0].transAxes, weight="bold", fontsize=13
    )
    axes[1].text(
        -0.14, 1.07, "B", transform=axes[1].transAxes, weight="bold", fontsize=13
    )
    figure.suptitle(
        "Mechanistic fate of within-host mutation", fontsize=14, weight="bold"
    )
    figure.tight_layout(rect=(0, 0, 1, 0.95))
    save_figure(figure, figure_dir, "fig07-mutant-fate")


def main() -> int:
    global CELL_ORDER

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--repository",
        type=Path,
        default=Path(__file__).resolve().parents[5],
    )
    parser.add_argument(
        "--tier",
        choices=("core", "all"),
        default="all",
        help="analyse the completed 12-cell core or require all 20 cells",
    )
    args = parser.parse_args()
    repository = args.repository.resolve()
    work = repository / "experiments" / "work" / "trophosome"
    phase = work / "p01-neutral-feedback"
    layout = json.loads((work / "layout.local.json").read_text(encoding="utf-8"))
    scratch = Path(layout["scratch"])
    CELL_ORDER = CORE_CELL_ORDER if args.tier == "core" else FULL_CELL_ORDER
    derived_dir = (
        phase
        / "analysis"
        / ("s01-pilot-core-derived" if args.tier == "core" else "s01-pilot-derived")
    )
    figure_dir = repository / "docs" / "figures" / "phase1-pilot"

    cell_rows = read_tsv(phase / "design" / "phase1-first-pilot-cells.tsv")
    cells: dict[str, dict[str, Any]] = {}
    for row in cell_rows:
        if row["cell_id"] not in CELL_ORDER:
            continue
        cells[row["cell_id"]] = {
            **row,
            "H": int(row["H"]),
            "f": float(row["f"]),
            "e": int(row["e"]),
            "R": int(row["R"]),
            "alpha": float(row["alpha"]),
            "u": float(row["u"]),
        }
    run_rows = [
        row
        for row in read_tsv(phase / "manifests" / "phase1-first-pilot-runs.tsv")
        if row["cell_id"] in CELL_ORDER
    ]
    initial_payload = json.loads(
        (work / "common" / "initial-populations" / "ip001-fisher100.json").read_text(
            encoding="utf-8"
        )
    )
    initial_counts = {
        strain_id: int(count)
        for strain_id, count in enumerate(initial_payload["scaled_counts"])
    }
    initial_metrics = diversity_metrics(initial_counts)

    endpoints: list[dict[str, Any]] = []
    trajectories: list[dict[str, Any]] = []
    audit_issues: list[str] = []
    commits: set[str] = set()
    source_hashes: set[str] = set()
    platforms: set[str] = set()

    for run in run_rows:
        cell_id = run["cell_id"]
        cell = cells[cell_id]
        output = scratch / run["scratch_relative_path"]
        completion = json.loads(
            (output / "completion.json").read_text(encoding="utf-8")
        )
        if completion.get("complete") is not True:
            audit_issues.append(f"{run['run_id']}: completion is not committed")
        final_path = output / "final_environment_rep000.npz"
        if (
            sha256(final_path)
            != completion["final_environment_sha256"][final_path.name]
        ):
            audit_issues.append(f"{run['run_id']}: final checksum differs")
        final_counts = read_environment(output / "environment_counts.csv")
        if sum(final_counts.values()) != 1_000_000_000:
            audit_issues.append(f"{run['run_id']}: final capacity differs")
        final_metrics = diversity_metrics(final_counts)
        composition = composition_metrics(final_counts, initial_counts)
        parents, generated_lineages, maximum_transition_mutants = read_lineages(
            output / "strain_lineage_events.csv", final_counts
        )
        roots = trace_roots(parents, final_counts)
        root_counts: dict[int, int] = defaultdict(int)
        for strain_id, count in final_counts.items():
            root_counts[roots[strain_id]] += count
        root_metrics = diversity_metrics(dict(root_counts))
        root_composition = composition_metrics(dict(root_counts), initial_counts)
        adult_mutants = unique_mutant_ids(output / "host_adult_counts.csv")
        released_mutants = unique_mutant_ids(output / "release_counts.csv")
        final_mutants = {
            strain_id for strain_id in final_counts if strain_id >= INITIAL_RICHNESS
        }
        mutant_abundance = sum(final_counts[strain_id] for strain_id in final_mutants)
        pooling = pooling_metrics(output, cell["H"])
        execution = json.loads(
            (output / "execution-summary.json").read_text(encoding="utf-8")
        )
        provenance = json.loads(
            (output / "provenance.json").read_text(encoding="utf-8")
        )
        commits.add(provenance["git_commit"])
        source_hashes.add(provenance["source_sha256"])
        platforms.add(provenance["platform"])

        endpoint: dict[str, Any] = {
            "run_id": run["run_id"],
            "cell_id": cell_id,
            "cell": cell_short(cell_id),
            "seed_block_id": run["seed_block_id"],
            "master_seed": int(run["master_seed"]),
            "H": cell["H"],
            "f": cell["f"],
            "e": cell["e"],
            "R": cell["R"],
            "alpha": cell["alpha"],
            "u": cell["u"],
            **final_metrics,
            "richness_change_pct": 100.0
            * (final_metrics["richness"] / initial_metrics["richness"] - 1.0),
            "hill_q1_change_pct": 100.0
            * (final_metrics["hill_q1"] / initial_metrics["hill_q1"] - 1.0),
            "hill_q2_change_pct": 100.0
            * (final_metrics["hill_q2"] / initial_metrics["hill_q2"] - 1.0),
            "evenness_change": final_metrics["evenness"] - initial_metrics["evenness"],
            **composition,
            "original_strains_lost": sum(
                strain_id not in final_counts for strain_id in INITIAL_IDS
            ),
            "original_strains_retained": sum(
                strain_id in final_counts for strain_id in INITIAL_IDS
            ),
            "final_mutant_richness": len(final_mutants),
            "mutant_abundance_fraction": mutant_abundance / sum(final_counts.values()),
            "generated_mutant_lineages": generated_lineages,
            "maximum_mutants_in_one_transition": maximum_transition_mutants,
            "adult_mutant_lineages": len(adult_mutants),
            "released_mutant_lineages": len(released_mutants),
            "root_richness": root_metrics["richness"],
            "root_hill_q1": root_metrics["hill_q1"],
            "root_hill_q2": root_metrics["hill_q2"],
            "root_evenness": root_metrics["evenness"],
            "root_total_variation": root_composition["total_variation"],
            **pooling,
            "elapsed_minutes": execution["elapsed_seconds"] / 60.0,
            "output_mib": execution["output_bytes"] / 1024**2,
            "peak_rss_mib": execution["peak_process_tree_rss_kib"] / 1024.0,
        }
        endpoints.append(endpoint)

        trajectories.append(
            {
                "cell_id": cell_id,
                "cell": cell_short(cell_id),
                "seed_block_id": run["seed_block_id"],
                "generation": 0,
                "environment_richness": int(initial_metrics["richness"]),
                "environment_gene_diversity": initial_metrics["gene_diversity"],
                "mean_adult_richness": "",
                "mean_adult_gene_diversity": "",
            }
        )
        with (output / "host_generation_summary.csv").open(
            newline="", encoding="utf-8"
        ) as handle:
            summary_rows = list(csv.DictReader(handle))
        if len(summary_rows) != 5:
            audit_issues.append(f"{run['run_id']}: expected five generation summaries")
        for summary in summary_rows:
            trajectories.append(
                {
                    "cell_id": cell_id,
                    "cell": cell_short(cell_id),
                    "seed_block_id": run["seed_block_id"],
                    "generation": int(summary["host_generation"]),
                    "environment_richness": int(summary["environment_richness"]),
                    "environment_gene_diversity": float(
                        summary["environment_gene_diversity"]
                    ),
                    "mean_adult_richness": float(summary["mean_adult_richness"]),
                    "mean_adult_gene_diversity": float(
                        summary["mean_adult_gene_diversity"]
                    ),
                }
            )

    endpoints.sort(
        key=lambda row: (
            CELL_ORDER.index(row["cell_id"]),
            SEED_ORDER.index(row["seed_block_id"]),
        )
    )
    trajectories.sort(
        key=lambda row: (
            CELL_ORDER.index(row["cell_id"]),
            SEED_ORDER.index(row["seed_block_id"]),
            row["generation"],
        )
    )

    endpoint_fields = list(endpoints[0])
    write_tsv(derived_dir / "run-endpoints.tsv", endpoints, endpoint_fields)
    trajectory_fields = list(trajectories[0])
    write_tsv(
        derived_dir / "environment-trajectories.tsv", trajectories, trajectory_fields
    )

    summary_measures = (
        "richness",
        "shannon",
        "simpson",
        "hill_q1",
        "hill_q2",
        "evenness",
        "total_variation",
        "original_strains_lost",
        "final_mutant_richness",
        "mutant_abundance_fraction",
        "generated_mutant_lineages",
        "maximum_mutants_in_one_transition",
        "adult_mutant_lineages",
        "released_mutant_lineages",
        "root_richness",
        "root_total_variation",
        "infection_richness_g5",
        "pooled_adult_richness_g5",
        "release_richness_g5",
        "median_original_occupancy_hosts_g5",
        "median_original_occupancy_fraction_g5",
        "elapsed_minutes",
        "output_mib",
        "peak_rss_mib",
    )
    cell_summaries: list[dict[str, Any]] = []
    for cell_id in CELL_ORDER:
        rows = [row for row in endpoints if row["cell_id"] == cell_id]
        summary: dict[str, Any] = {
            "cell_id": cell_id,
            "cell": cell_short(cell_id),
            "label": cells[cell_id]["label"],
            "H": cells[cell_id]["H"],
            "f": cells[cell_id]["f"],
            "R": cells[cell_id]["R"],
            "alpha": cells[cell_id]["alpha"],
            "u": cells[cell_id]["u"],
        }
        for measure in summary_measures:
            statistics = describe([float(row[measure]) for row in rows])
            for statistic, value in statistics.items():
                summary[f"{measure}_{statistic}"] = value
        cell_summaries.append(summary)
    write_tsv(
        derived_dir / "cell-summaries.tsv", cell_summaries, list(cell_summaries[0])
    )

    resources_by_seed = {}
    for seed_id in SEED_ORDER:
        rows = [row for row in endpoints if row["seed_block_id"] == seed_id]
        resources_by_seed[seed_id] = {
            "populations": len(rows),
            "elapsed_minutes": sum(row["elapsed_minutes"] for row in rows),
            "output_mib": sum(row["output_mib"] for row in rows),
            "peak_rss_mib": max(row["peak_rss_mib"] for row in rows),
        }
    analysis_summary = {
        "analysis_schema_version": "1.0.0",
        "experiment_id": (
            "phase1-first-pilot-core12"
            if args.tier == "core"
            else "phase1-first-pilot-20cell"
        ),
        "pilot_tier": args.tier,
        "analysis_scope": "exploratory feasibility and calibration pilot",
        "raw_scratch": str(scratch),
        "populations": len(endpoints),
        "cells": len(cells),
        "seed_blocks": list(SEED_ORDER),
        "software_git_commits": sorted(commits),
        "source_sha256": sorted(source_hashes),
        "benchmark_platforms": sorted(platforms),
        "benchmark_hardware_note": (
            "Local Mac: 1.2 GHz quad-core Intel Core i7 with 8 GB RAM"
        ),
        "initial_metrics": initial_metrics,
        "resources_by_seed_block": resources_by_seed,
        "overall_resources": {
            "elapsed_hours": sum(row["elapsed_minutes"] for row in endpoints) / 60.0,
            "output_mib": sum(row["output_mib"] for row in endpoints),
            "peak_rss_mib": max(row["peak_rss_mib"] for row in endpoints),
        },
        "mutation_materialization": {
            "configured_limit_per_transition": 100000,
            "largest_realized_count": max(
                row["maximum_mutants_in_one_transition"] for row in endpoints
            ),
            "largest_fraction_of_limit": max(
                row["maximum_mutants_in_one_transition"] for row in endpoints
            )
            / 100000,
        },
        "audit_status": "PASS" if not audit_issues else "FAIL",
        "audit_issues": audit_issues,
        "limitations": [
            "Three independent populations per cell support descriptive "
            "pilot inference only.",
            "Five host passages do not test equilibrium.",
            "Complete labelled environmental composition is retained only "
            "at generation 5.",
        ],
    }
    derived_dir.mkdir(parents=True, exist_ok=True)
    (derived_dir / "analysis-summary.json").write_text(
        json.dumps(analysis_summary, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    if audit_issues:
        raise SystemExit("analysis audit failed:\n" + "\n".join(audit_issues))

    configure_plot_style()
    plot_design(cells, figure_dir)
    plot_feasibility(endpoints, cells, figure_dir)
    plot_trajectories(trajectories, cells, figure_dir)
    plot_endpoints(endpoints, cells, figure_dir)
    plot_matched_return(endpoints, cells, figure_dir)
    plot_mutation_bracket(endpoints, cells, figure_dir)
    plot_mutant_fate(endpoints, figure_dir)
    print(f"Analysed {len(endpoints)} populations across {len(cells)} cells")
    print(f"Derived tables: {derived_dir}")
    print(f"Figures: {figure_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
