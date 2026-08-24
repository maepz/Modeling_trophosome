"""Build a concise, biologist-readable PDF report from standardized pilot results."""

from __future__ import annotations

import csv
import html
import json
import math
from collections import defaultdict
from dataclasses import dataclass
from datetime import date
from pathlib import Path
from typing import Any


@dataclass(frozen=True)
class ReportArtifacts:
    """Paths created by :func:`generate_pilot_report`."""

    pdf: Path
    markdown: Path
    assets: Path


SLIDE28_MAIN_PURPOSES = {
    "c0001": "No-return baseline",
    "c0002": "Mutation occurs in hosts but nothing returns",
    "c0003": "Fixed-return, feedback comparison",
    "c0004": "Fixed-return, pooling comparison",
    "c0005": "Fixed-return, pooling comparison",
    "c0006": "Fixed-return, pooling comparison",
    "c0007": "Lowest positive mutation level",
    "c0008": "Low mutation",
    "c0009": "Intermediate mutation, pooling comparison",
    "c0010": "High mutation",
    "c0011": "feedback comparison",
    "c0012": "Strong feedback boundary",
    "c0013": "Intermediate mutation, pooling comparison",
    "c0014": "Intermediate mutation, pooling comparison",
    "c0015": "Intermediate mutation, pooling comparison",
    "c0016": "Weaker fixed-return, pooling comparison; feedback comparison",
    "c0017": "Weaker fixed-return, pooling comparison",
    "c0018": "feedback comparison",
    "c0019": "More hosts feedback comparison",
    "c0020": (
        "Weaker fixed-return, pooling comparison; More hosts, feedback comparison"
    ),
}

# Resolved from the Office theme and tint/shade transforms used by slide 28.
SLIDE28_PALETTE = {
    "black": "#000000",
    "white": "#FFFFFF",
    "grey": "#BFBFBF",
    "purple_light": "#ECD5E9",
    "purple_medium": "#C06AB9",
    "purple_dark": "#78206E",
    "blue_grey": "#D0DFE6",
    "blue_light": "#CFECF7",
    "blue_medium": "#6FC5E6",
    "green_light": "#DCEED5",
}


def _read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def _number(value: Any, default: float = 0.0) -> float:
    if value in (None, ""):
        return default
    return float(value)


def _median(values: list[float]) -> float:
    ordered = sorted(values)
    middle = len(ordered) // 2
    if len(ordered) % 2:
        return ordered[middle]
    return 0.5 * (ordered[middle - 1] + ordered[middle])


def _short_cell(identifier: str) -> str:
    return identifier.rsplit("-", 1)[-1]


def _format_number(value: float, digits: int = 3) -> str:
    if value == 0:
        return "0"
    absolute = abs(value)
    if absolute >= 10000 or absolute < 0.001:
        return f"{value:.{digits - 1}e}"
    return f"{value:.{digits}f}".rstrip("0").rstrip(".")


def _format_integer(value: float) -> str:
    return f"{int(round(value)):,}"


def _normalise_endpoints(
    endpoints: list[dict[str, str]], initial: dict[str, Any]
) -> None:
    initial_shannon = _number(initial.get("shannon"))
    initial_simpson = _number(
        initial.get("simpson"), 1.0 - 1.0 / _number(initial.get("hill_q2"), 1.0)
    )
    initial["simpson"] = initial_simpson
    for row in endpoints:
        if row.get("shannon", "") == "":
            row["shannon"] = str(math.log(_number(row["hill_q1"])))
        if row.get("simpson", "") == "":
            row["simpson"] = str(1.0 - 1.0 / _number(row["hill_q2"], 1.0))
        if row.get("shannon_change_pct", "") == "":
            row["shannon_change_pct"] = str(
                100.0 * (_number(row["shannon"]) / initial_shannon - 1.0)
            )
        if row.get("simpson_change", "") == "":
            row["simpson_change"] = str(_number(row["simpson"]) - initial_simpson)


def _cell_rows(endpoints: list[dict[str, str]], cell_id: str) -> list[dict[str, str]]:
    return [row for row in endpoints if row["cell_id"] == cell_id]


def _cell_median(endpoints: list[dict[str, str]], cell_id: str, measure: str) -> float:
    return _median([_number(row[measure]) for row in _cell_rows(endpoints, cell_id)])


def _load_report_data(
    analysis_dir: Path, design_path: Path
) -> tuple[
    dict[str, Any],
    list[dict[str, str]],
    list[dict[str, str]],
    list[dict[str, str]],
]:
    required = (
        analysis_dir / "analysis-summary.json",
        analysis_dir / "run-endpoints.tsv",
        analysis_dir / "environment-trajectories.tsv",
        design_path,
    )
    missing = [str(path) for path in required if not path.is_file()]
    if missing:
        raise FileNotFoundError("Missing report input(s): " + ", ".join(missing))
    summary = json.loads(required[0].read_text(encoding="utf-8"))
    endpoints = _read_tsv(required[1])
    trajectories = _read_tsv(required[2])
    design = _read_tsv(design_path)
    if not endpoints:
        raise ValueError("run-endpoints.tsv contains no populations")
    if not design:
        raise ValueError("design table contains no experimental cells")
    _normalise_endpoints(endpoints, summary["initial_metrics"])
    return summary, endpoints, trajectories, design


def _design_cell_id(row: dict[str, str]) -> str:
    return row["cell_id"]


def _comparison_sets(row: dict[str, str]) -> set[str]:
    return {value for value in row.get("comparison_sets", "").split("|") if value}


def _migration_fraction(design: list[dict[str, str]]) -> float:
    values = {_number(row.get("m")) for row in design if row.get("m", "") != ""}
    return next(iter(values)) if len(values) == 1 else 0.0


def _matched_return_group(
    design: list[dict[str, str]],
) -> list[dict[str, str]]:
    groups: dict[tuple[float, float, float], list[dict[str, str]]] = defaultdict(list)
    for row in design:
        return_size = _number(row.get("R"))
        mutation = _number(row.get("u"))
        if return_size > 0 and mutation == 0:
            key = (return_size, _number(row.get("alpha")), mutation)
            groups[key].append(row)
    candidates = [
        rows
        for rows in groups.values()
        if len({_number(row.get("H")) for row in rows}) >= 2
    ]
    if not candidates:
        return []
    return sorted(max(candidates, key=len), key=lambda row: _number(row.get("H")))


def _all_return_pooling_groups(
    design: list[dict[str, str]],
) -> list[list[dict[str, str]]]:
    groups: dict[tuple[float, float, float], list[dict[str, str]]] = defaultdict(list)
    for row in design:
        return_size = _number(row.get("R"))
        if return_size <= 0:
            continue
        key = (
            return_size,
            _number(row.get("alpha")),
            _number(row.get("u")),
        )
        groups[key].append(row)
    candidates = [
        sorted(rows, key=lambda row: _number(row.get("H")))
        for rows in groups.values()
        if len({_number(row.get("H")) for row in rows}) >= 2
    ]
    return sorted(
        candidates,
        key=lambda rows: (
            _number(rows[0].get("u")),
            -_number(rows[0].get("R")),
            -len(rows),
        ),
    )


def _extension_return_groups(
    design: list[dict[str, str]],
) -> list[list[dict[str, str]]]:
    primary = {_design_cell_id(row) for row in _matched_return_group(design)}
    return [
        rows
        for rows in _all_return_pooling_groups(design)
        if {_design_cell_id(row) for row in rows} != primary
        and any(
            value.startswith("matched-")
            for value in set.intersection(*(_comparison_sets(row) for row in rows))
        )
    ]


def _feedback_groups_by_host(
    design: list[dict[str, str]],
) -> list[list[dict[str, str]]]:
    """Return H-fixed mutation-free alpha series with at least three levels."""

    by_host: dict[int, list[dict[str, str]]] = defaultdict(list)
    for row in design:
        if _number(row.get("R")) > 0 and _number(row.get("u")) == 0:
            by_host[int(_number(row.get("H")))].append(row)
    groups = [
        sorted(rows, key=lambda row: _number(row.get("alpha")))
        for rows in by_host.values()
        if len({_number(row.get("alpha")) for row in rows}) >= 3
    ]
    return sorted(groups, key=lambda rows: _number(rows[0].get("H")))


def _figure2_pooling_groups(
    design: list[dict[str, str]],
) -> list[list[dict[str, str]]]:
    """Return the alpha=0.5 and alpha~=0.09 pooling series for Figure 2."""

    requested_series = {
        (1e9, 0.5, 0.0),
        (1e9, 0.5, 1e-10),
        (1e8, 1 / 11, 0.0),
    }
    groups = [
        group
        for group in _all_return_pooling_groups(design)
        if any(
            math.isclose(_number(group[0].get("R")), return_size)
            and math.isclose(_number(group[0].get("alpha")), alpha)
            and math.isclose(_number(group[0].get("u")), mutation)
            for return_size, alpha, mutation in requested_series
        )
    ]
    return sorted(
        groups,
        key=lambda group: (
            -_number(group[0].get("alpha")),
            _number(group[0].get("u")),
        ),
    )


def _mutation_group(design: list[dict[str, str]]) -> list[dict[str, str]]:
    groups: dict[tuple[float, float, float], list[dict[str, str]]] = defaultdict(list)
    for row in design:
        if _number(row.get("R")) <= 0:
            continue
        key = (
            _number(row.get("H")),
            _number(row.get("R")),
            _number(row.get("alpha")),
        )
        groups[key].append(row)
    candidates = [
        rows
        for rows in groups.values()
        if len({_number(row.get("u")) for row in rows}) >= 2
        and any(_number(row.get("u")) > 0 for row in rows)
    ]
    if not candidates:
        return []
    return sorted(max(candidates, key=len), key=lambda row: _number(row.get("u")))


def _feedback_group(design: list[dict[str, str]]) -> list[dict[str, str]]:
    """Return the first fixed-H alpha series, or the legacy broad bracket.

    The 20-cell design contains matched alpha series at H=100 and H=1,000.
    The 12- and 17-cell designs retain the earlier broad bracket, in which H
    also changes, so old pilot reports can still be regenerated honestly.
    """

    fixed_host_groups = _feedback_groups_by_host(design)
    if fixed_host_groups:
        return fixed_host_groups[0]

    by_alpha: dict[float, list[dict[str, str]]] = defaultdict(list)
    for row in design:
        if _number(row.get("R")) > 0 and _number(row.get("u")) == 0:
            by_alpha[_number(row.get("alpha"))].append(row)
    if len(by_alpha) < 3:
        return []
    representatives = [
        min(rows, key=lambda row: _number(row.get("H")))
        for _, rows in sorted(by_alpha.items())
    ]
    return sorted(representatives, key=lambda row: _number(row.get("alpha")))


def _seed_block(row: dict[str, str]) -> str:
    return row.get("seed_block_id", row.get("seed_block", ""))


def _paired_group_values(
    endpoints: list[dict[str, str]],
    group: list[dict[str, str]],
    measure: str,
) -> dict[str, list[float]]:
    """Collect one endpoint per cell in each complete seed block."""

    values_by_cell_seed: dict[tuple[str, str], float] = {}
    seed_sets: list[set[str]] = []
    for design_row in group:
        cell_id = _design_cell_id(design_row)
        cell_seeds: set[str] = set()
        for endpoint in _cell_rows(endpoints, cell_id):
            seed = _seed_block(endpoint)
            if seed:
                values_by_cell_seed[(cell_id, seed)] = _number(endpoint[measure])
                cell_seeds.add(seed)
        seed_sets.append(cell_seeds)
    complete_seeds = sorted(set.intersection(*seed_sets)) if seed_sets else []
    return {
        seed: [
            values_by_cell_seed[(_design_cell_id(design_row), seed)]
            for design_row in group
        ]
        for seed in complete_seeds
    }


def _log_log_slope(x_values: list[float], y_values: list[float]) -> float:
    x_logs = [math.log10(value) for value in x_values]
    y_logs = [math.log10(value) for value in y_values]
    x_mean = sum(x_logs) / len(x_logs)
    y_mean = sum(y_logs) / len(y_logs)
    denominator = sum((value - x_mean) ** 2 for value in x_logs)
    if denominator == 0:
        return 0.0
    return (
        sum(
            (x_value - x_mean) * (y_value - y_mean)
            for x_value, y_value in zip(x_logs, y_logs, strict=True)
        )
        / denominator
    )


def _cell_order(design: list[dict[str, str]]) -> list[str]:
    return [_design_cell_id(row) for row in design]


def _configure_plots() -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 8.5,
            "axes.titlesize": 10,
            "axes.labelsize": 9,
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


def _save_figure(figure: Any, path: Path) -> None:
    import matplotlib.pyplot as plt

    figure.savefig(path, dpi=220, bbox_inches="tight")
    plt.close(figure)


def _plot_endpoint_summary(
    endpoints: list[dict[str, str]], design: list[dict[str, str]], path: Path
) -> None:
    import matplotlib.pyplot as plt

    measures = (
        ("richness_change_pct", "Richness change (%)", (-5.0, 5.0)),
        ("shannon_change_pct", "Shannon entropy change (%)", None),
        ("simpson_change", "Simpson change", None),
        ("hill_q1_change_pct", "D1 change (%)", (-5.0, 5.0)),
        ("hill_q2_change_pct", "D2 change (%)", (-5.0, 5.0)),
        ("evenness_change", "Pielou evenness change", (-0.02, 0.02)),
        ("total_variation", "Composition distance", (0.05,)),
    )
    order = _cell_order(design)
    labels = [_short_cell(cell_id) for cell_id in order]
    figure, axes = plt.subplots(2, 4, figsize=(14, 6.8))
    for axis, (measure, ylabel, margins) in zip(
        list(axes.flat)[:7], measures, strict=True
    ):
        for index, cell_id in enumerate(order):
            values = [_number(row[measure]) for row in _cell_rows(endpoints, cell_id)]
            offsets = (-0.16, 0.0, 0.16)
            for replicate, value in enumerate(values):
                offset = offsets[replicate] if replicate < len(offsets) else 0.0
                axis.scatter(index + offset, value, color="#2C7FB8", alpha=0.65, s=24)
            axis.scatter(
                index,
                _median(values),
                marker="_",
                s=160,
                linewidth=2.0,
                color="#172B3A",
                zorder=5,
            )
        if margins:
            for margin in margins:
                axis.axhline(margin, color="#7F8C8D", linestyle="--", linewidth=0.9)
        axis.set_ylabel(ylabel)
        axis.set_xticks(range(len(labels)))
        axis.set_xticklabels(labels, rotation=55, ha="right")
    axes[1, 3].set_axis_off()
    axes[1, 3].text(
        0.04,
        0.92,
        "Dots: independent populations\n"
        "Black bars: cell medians\n"
        "Dashed lines: relevance margins",
        va="top",
        transform=axes[1, 3].transAxes,
        fontsize=10,
    )
    figure.suptitle("Environmental diversity after the pilot passages", weight="bold")
    figure.tight_layout(rect=(0, 0, 1, 0.95))
    _save_figure(figure, path)


def _plot_feasibility(
    endpoints: list[dict[str, str]], design: list[dict[str, str]], path: Path
) -> None:
    import matplotlib.pyplot as plt

    matched = _matched_return_group(design)
    mutation = _mutation_group(design)
    figure, axes = plt.subplots(2, 3, figsize=(12.5, 6.5))
    measures = (
        ("elapsed_minutes", "Runtime per population (min)", True),
        ("output_mib", "Diagnostic output (MiB)", True),
        ("peak_rss_mib", "Peak memory (MiB)", False),
    )
    seed_colors = ("#0072B2", "#E69F00", "#009E73", "#CC79A7")
    for column, (measure, ylabel, log_y) in enumerate(measures):
        axis = axes[0, column]
        hosts = [_number(row.get("H")) for row in matched]
        paired = _paired_group_values(endpoints, matched, measure)
        for index, (seed, values) in enumerate(paired.items()):
            axis.plot(
                hosts,
                values,
                color=seed_colors[index % len(seed_colors)],
                marker="o",
                linewidth=1.2,
                alpha=0.75,
                label=seed,
            )
        medians = [
            _cell_median(endpoints, _design_cell_id(row), measure) for row in matched
        ]
        axis.plot(hosts, medians, color="#172B3A", marker="o", linewidth=2.2)
        axis.set_xscale("log")
        if log_y:
            axis.set_yscale("log")
        axis.set_xlabel("# hosts H (u=0; total return fixed)")
        axis.set_ylabel(ylabel)
        if column == 0 and paired:
            axis.legend(title="Seed block", fontsize=7, title_fontsize=7)

        axis = axes[1, column]
        labels = [
            "0" if _number(row.get("u")) == 0 else f"{_number(row.get('u')):.0e}"
            for row in mutation
        ]
        x_values = list(range(len(mutation)))
        paired = _paired_group_values(endpoints, mutation, measure)
        for index, values in enumerate(paired.values()):
            axis.plot(
                x_values,
                values,
                color=seed_colors[index % len(seed_colors)],
                marker="o",
                linewidth=1.2,
                alpha=0.75,
            )
        medians = [
            _cell_median(endpoints, _design_cell_id(row), measure) for row in mutation
        ]
        axis.plot(x_values, medians, color="#172B3A", marker="o", linewidth=2.2)
        if log_y:
            axis.set_yscale("log")
        axis.set_xticks(x_values)
        axis.set_xticklabels(labels)
        axis.set_xlabel("Mutation rate u (H=100)")
        axis.set_ylabel(ylabel)
    figure.suptitle("Observed computational scaling in the pilot", weight="bold")
    figure.tight_layout(rect=(0, 0, 1, 0.93))
    _save_figure(figure, path)


def _plot_matched_return(
    endpoints: list[dict[str, str]],
    groups: list[list[dict[str, str]]],
    path: Path,
) -> None:
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D

    figure, axes = plt.subplots(1, 2, figsize=(10.5, 3.8))
    measures = (
        ("total_variation", "Composition distance"),
        ("release_richness_g5", "Release-pool richness at final passage"),
    )
    seed_colors = ("#0072B2", "#E69F00", "#009E73", "#CC79A7")
    mutation_styles = {
        0.0: {"marker": "o", "label": "u=0"},
        1e-10: {"marker": "s", "label": "u=1e-10"},
    }

    def alpha_style(alpha: float) -> dict[str, Any]:
        if math.isclose(alpha, 0.5):
            return {"linestyle": "-", "filled": True, "label": "alpha=0.5"}
        return {
            "linestyle": "--",
            "filled": False,
            "label": "alpha=0.09",
        }

    seed_names = sorted(
        {
            seed
            for group in groups
            for seed in _paired_group_values(endpoints, group, measures[0][0])
        }
    )
    colours_by_seed = {
        seed: seed_colors[index % len(seed_colors)]
        for index, seed in enumerate(seed_names)
    }
    for axis, (measure, ylabel) in zip(axes, measures, strict=True):
        for group in groups:
            hosts = [_number(row.get("H")) for row in group]
            mutation = _number(group[0].get("u"))
            alpha = _number(group[0].get("alpha"))
            mutation_style = mutation_styles[mutation]
            feedback_style = alpha_style(alpha)
            medians = [
                _cell_median(endpoints, _design_cell_id(row), measure) for row in group
            ]
            axis.plot(
                hosts,
                medians,
                color="#172B3A",
                marker=mutation_style["marker"],
                linestyle=feedback_style["linestyle"],
                markerfacecolor=("#172B3A" if feedback_style["filled"] else "none"),
                markeredgecolor="#172B3A",
                linewidth=2.8,
                alpha=0.45,
                zorder=1,
            )
            paired = _paired_group_values(endpoints, group, measure)
            for seed, values in paired.items():
                axis.plot(
                    hosts,
                    values,
                    color=colours_by_seed[seed],
                    marker=mutation_style["marker"],
                    linestyle=feedback_style["linestyle"],
                    markerfacecolor=(
                        colours_by_seed[seed] if feedback_style["filled"] else "none"
                    ),
                    markeredgecolor=colours_by_seed[seed],
                    markeredgewidth=1.15,
                    linewidth=1.35,
                    markersize=5,
                    alpha=0.9,
                    zorder=2,
                )
        axis.set_xscale("log")
        axis.set_xlabel("# hosts H")
        axis.set_ylabel(ylabel)
    seed_handles = [
        Line2D([0], [0], color=colours_by_seed[seed], linewidth=2, label=seed)
        for seed in seed_names
    ]
    if seed_handles:
        axes[0].legend(
            handles=seed_handles,
            title="Seed block (colour)",
            fontsize=7,
            title_fontsize=7,
        )
    mutation_rates = sorted({_number(group[0].get("u")) for group in groups})
    mutation_handles = [
        Line2D(
            [0],
            [0],
            color="#172B3A",
            marker=mutation_styles[mutation]["marker"],
            linestyle="none",
            markerfacecolor="#172B3A",
            label=mutation_styles[mutation]["label"],
        )
        for mutation in mutation_rates
    ]
    feedback_levels = sorted(
        {_number(group[0].get("alpha")) for group in groups}, reverse=True
    )
    feedback_handles = [
        Line2D(
            [0],
            [0],
            color="#172B3A",
            marker="o",
            linestyle=alpha_style(alpha)["linestyle"],
            markerfacecolor=("#172B3A" if alpha_style(alpha)["filled"] else "none"),
            markeredgecolor="#172B3A",
            linewidth=1.7,
            label=alpha_style(alpha)["label"],
        )
        for alpha in feedback_levels
    ]
    encoding_handles = mutation_handles + feedback_handles
    encoding_handles.append(
        Line2D(
            [0],
            [0],
            color="#172B3A",
            linewidth=3,
            alpha=0.45,
            label="Cell median",
        )
    )
    axes[1].legend(
        handles=encoding_handles,
        title="Encoding",
        fontsize=7,
        title_fontsize=7,
    )
    figure.suptitle(
        "Pooling across hosts within fixed-return feedback levels", weight="bold"
    )
    figure.tight_layout(rect=(0, 0, 1, 0.92))
    _save_figure(figure, path)


def _plot_mutation(
    endpoints: list[dict[str, str]], group: list[dict[str, str]], path: Path
) -> None:
    import matplotlib.pyplot as plt

    figure, axes = plt.subplots(1, 3, figsize=(12.5, 3.8))
    measures = (
        ("richness", "Labelled richness"),
        ("evenness", "Pielou evenness"),
        ("mutant_abundance_fraction", "Mutant abundance fraction"),
    )
    labels = [
        "0" if _number(row.get("u")) == 0 else f"{_number(row.get('u')):.0e}"
        for row in group
    ]
    x_values = list(range(len(group)))
    seed_colors = ("#0072B2", "#E69F00", "#009E73", "#CC79A7")
    for axis, (measure, ylabel) in zip(axes, measures, strict=True):
        paired = _paired_group_values(endpoints, group, measure)
        for index, (seed, values) in enumerate(paired.items()):
            axis.plot(
                x_values,
                values,
                color=seed_colors[index % len(seed_colors)],
                marker="o",
                linewidth=1.3,
                alpha=0.78,
                label=seed,
            )
        medians = [
            _cell_median(endpoints, _design_cell_id(row), measure) for row in group
        ]
        axis.plot(x_values, medians, color="#172B3A", marker="o", linewidth=2)
        axis.set_xticks(x_values)
        axis.set_xticklabels(labels)
        axis.set_xlabel("Mutation probability, u")
        axis.set_ylabel(ylabel)
    if _paired_group_values(endpoints, group, measures[0][0]):
        axes[0].legend(title="Seed block", fontsize=7, title_fontsize=7)
    figure.suptitle("Within-host mutation and environmental novelty", weight="bold")
    figure.tight_layout(rect=(0, 0, 1, 0.92))
    _save_figure(figure, path)


def _plot_feedback(
    endpoints: list[dict[str, str]],
    groups: list[list[dict[str, str]]],
    path: Path,
) -> None:
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D

    figure, axes = plt.subplots(1, 3, figsize=(12.5, 3.8))
    measures = (
        ("total_variation", "Composition distance"),
        ("hill_q1", "Hill D1"),
        ("evenness", "Pielou evenness"),
    )
    seed_colors = ("#0072B2", "#E69F00", "#009E73", "#CC79A7")
    fixed_host_comparison = all(
        len({int(_number(row.get("H"))) for row in group}) == 1 for group in groups
    )
    if not fixed_host_comparison:
        group = groups[0]
        labels = [_format_number(_number(row.get("alpha")), 3) for row in group]
        x_values = list(range(len(group)))
        for axis, (measure, ylabel) in zip(axes, measures, strict=True):
            paired = _paired_group_values(endpoints, group, measure)
            for index, (seed, values) in enumerate(paired.items()):
                axis.plot(
                    x_values,
                    values,
                    color=seed_colors[index % len(seed_colors)],
                    marker="o",
                    linewidth=1.3,
                    alpha=0.78,
                    label=seed,
                )
            medians = [
                _cell_median(endpoints, _design_cell_id(row), measure) for row in group
            ]
            axis.plot(x_values, medians, color="#172B3A", marker="o", linewidth=2)
            axis.set_xticks(x_values)
            axis.set_xticklabels(labels)
            axis.set_xlabel("Host-derived fraction alpha")
            axis.set_ylabel(ylabel)
        if _paired_group_values(endpoints, group, measures[0][0]):
            axes[0].legend(title="Seed block", fontsize=7, title_fontsize=7)
        figure.suptitle(
            "Broad feedback bracket (H and total return also change)", weight="bold"
        )
        figure.tight_layout(rect=(0, 0, 1, 0.92))
        _save_figure(figure, path)
        return

    host_styles = (
        {"marker": "o", "linestyle": "-"},
        {"marker": "s", "linestyle": "--"},
        {"marker": "^", "linestyle": ":"},
    )
    seed_names = sorted(
        {
            seed
            for group in groups
            for seed in _paired_group_values(endpoints, group, measures[0][0])
        }
    )
    colours_by_seed = {
        seed: seed_colors[index % len(seed_colors)]
        for index, seed in enumerate(seed_names)
    }
    for axis, (measure, ylabel) in zip(axes, measures, strict=True):
        for group_index, group in enumerate(groups):
            style = host_styles[group_index % len(host_styles)]
            alpha_values = [_number(row.get("alpha")) for row in group]
            medians = [
                _cell_median(endpoints, _design_cell_id(row), measure) for row in group
            ]
            axis.plot(
                alpha_values,
                medians,
                color="#172B3A",
                marker=style["marker"],
                linestyle=style["linestyle"],
                linewidth=2.8,
                alpha=0.45,
                zorder=1,
            )
            for seed, values in _paired_group_values(endpoints, group, measure).items():
                axis.plot(
                    alpha_values,
                    values,
                    color=colours_by_seed[seed],
                    marker=style["marker"],
                    linestyle=style["linestyle"],
                    linewidth=1.35,
                    markersize=5,
                    alpha=0.9,
                    zorder=2,
                )
        axis.set_xscale("log")
        axis.set_xlabel("Host-derived fraction alpha")
        axis.set_ylabel(ylabel)
    seed_handles = [
        Line2D([0], [0], color=colours_by_seed[seed], linewidth=2, label=seed)
        for seed in seed_names
    ]
    axes[0].legend(
        handles=seed_handles,
        title="Seed block (colour)",
        fontsize=7,
        title_fontsize=7,
    )
    host_handles = [
        Line2D(
            [0],
            [0],
            color="#172B3A",
            marker=host_styles[index % len(host_styles)]["marker"],
            linestyle=host_styles[index % len(host_styles)]["linestyle"],
            linewidth=1.7,
            label=f"H={int(_number(group[0].get('H'))):,}",
        )
        for index, group in enumerate(groups)
    ]
    host_handles.append(
        Line2D(
            [0],
            [0],
            color="#172B3A",
            linewidth=3,
            alpha=0.45,
            label="Cell median",
        )
    )
    axes[2].legend(
        handles=host_handles,
        title="Host abundance (shape / line)",
        fontsize=7,
        title_fontsize=7,
    )
    figure.suptitle(
        "Effect of feedback strength at fixed host abundance", weight="bold"
    )
    figure.tight_layout(rect=(0, 0, 1, 0.92))
    _save_figure(figure, path)


def _plot_extension_return_groups(
    endpoints: list[dict[str, str]],
    groups: list[list[dict[str, str]]],
    path: Path,
) -> None:
    import matplotlib.pyplot as plt

    measures = (
        ("total_variation", "Composition distance"),
        ("richness", "Labelled richness D0"),
        ("evenness", "Pielou evenness"),
    )
    figure, axes = plt.subplots(len(groups), 3, figsize=(12.5, 3.3 * len(groups)))
    if len(groups) == 1:
        axes = [axes]
    seed_colors = ("#0072B2", "#E69F00", "#009E73", "#CC79A7")
    for row_axes, group in zip(axes, groups, strict=True):
        hosts = [_number(row.get("H")) for row in group]
        return_size = _number(group[0].get("R"))
        mutation = _number(group[0].get("u"))
        for column, (measure, ylabel) in enumerate(measures):
            axis = row_axes[column]
            paired = _paired_group_values(endpoints, group, measure)
            for index, (seed, values) in enumerate(paired.items()):
                axis.plot(
                    hosts,
                    values,
                    color=seed_colors[index % len(seed_colors)],
                    marker="o",
                    linewidth=1.3,
                    alpha=0.78,
                    label=seed,
                )
            medians = [
                _cell_median(endpoints, _design_cell_id(row), measure) for row in group
            ]
            axis.plot(hosts, medians, color="#172B3A", marker="o", linewidth=2)
            axis.set_xscale("log")
            axis.set_xlabel("# hosts H")
            axis.set_ylabel(ylabel)
            axis.set_title(f"R={return_size:.0e}; u={mutation:.0e}", fontsize=9)
        if _paired_group_values(endpoints, group, measures[0][0]):
            row_axes[0].legend(title="Seed block", fontsize=7, title_fontsize=7)
    figure.suptitle("Extension: paired pooling comparisons", weight="bold")
    figure.tight_layout(rect=(0, 0, 1, 0.96))
    _save_figure(figure, path)


def _plot_migration(
    trajectories: list[dict[str, str]],
    design: list[dict[str, str]],
    path: Path,
) -> None:
    import matplotlib.pyplot as plt

    order = _cell_order(design)
    labels = [_short_cell(cell_id) for cell_id in order]
    final_rows = [
        row for row in trajectories if int(_number(row.get("generation"))) == 5
    ]
    figure, axes = plt.subplots(1, 2, figsize=(13, 4.3))

    before = []
    after = []
    richness_changes: list[list[float]] = []
    for cell_id in order:
        rows = [row for row in final_rows if row["cell_id"] == cell_id]
        before.append(_median([_number(row["realized_host_feedback"]) for row in rows]))
        after.append(
            _median(
                [_number(row["expected_host_feedback_after_migration"]) for row in rows]
            )
        )
        richness_changes.append(
            [
                _number(row["post_migration_richness"])
                - _number(row["post_return_richness"])
                for row in rows
            ]
        )

    x_values = list(range(len(order)))
    axes[0].plot(
        x_values,
        before,
        color="#D55E00",
        marker="o",
        linewidth=1.8,
        label="After host return and regulation",
    )
    axes[0].plot(
        x_values,
        after,
        color="#0072B2",
        marker="s",
        linewidth=1.8,
        label="After regional exchange",
    )
    axes[0].set_ylabel("Expected host-derived fraction")
    axes[0].set_xticks(x_values)
    axes[0].set_xticklabels(labels, rotation=55, ha="right")
    axes[0].legend(fontsize=7)

    offsets = (-0.16, 0.0, 0.16)
    for index, values in enumerate(richness_changes):
        for replicate, value in enumerate(values):
            axes[1].scatter(
                index + offsets[replicate],
                value,
                color="#2C7FB8",
                alpha=0.7,
                s=28,
            )
        axes[1].scatter(
            index,
            _median(values),
            marker="_",
            s=170,
            linewidth=2,
            color="#172B3A",
            zorder=5,
        )
    axes[1].axhline(0, color="#7F8C8D", linewidth=0.9, linestyle="--")
    axes[1].set_ylabel("Richness change caused by migration")
    axes[1].set_xticks(x_values)
    axes[1].set_xticklabels(labels, rotation=55, ha="right")
    axes[0].text(-0.14, 1.06, "A", transform=axes[0].transAxes, weight="bold")
    axes[1].text(-0.14, 1.06, "B", transform=axes[1].transAxes, weight="bold")
    figure.suptitle(
        "Effect of fixed regional-pool exchange at passage 5", weight="bold"
    )
    figure.tight_layout(rect=(0, 0, 1, 0.93))
    _save_figure(figure, path)


def _make_figures(
    endpoints: list[dict[str, str]],
    trajectories: list[dict[str, str]],
    design: list[dict[str, str]],
    assets: Path,
) -> dict[str, Path]:
    assets.mkdir(parents=True, exist_ok=True)
    _configure_plots()
    figures = {
        "endpoint": assets / "endpoint-diversity.png",
        "feasibility": assets / "computational-feasibility.png",
    }
    _plot_endpoint_summary(endpoints, design, figures["endpoint"])
    _plot_feasibility(endpoints, design, figures["feasibility"])
    if (
        _migration_fraction(design) > 0
        and trajectories
        and "post_migration_richness" in trajectories[0]
    ):
        figures["migration"] = assets / "fixed-pool-migration.png"
        _plot_migration(trajectories, design, figures["migration"])
    matched = _matched_return_group(design)
    if matched:
        figures["matched"] = assets / "matched-return.png"
        _plot_matched_return(
            endpoints, _figure2_pooling_groups(design), figures["matched"]
        )
    mutation = _mutation_group(design)
    if mutation:
        figures["mutation"] = assets / "mutation-bracket.png"
        _plot_mutation(endpoints, mutation, figures["mutation"])
    feedback = _feedback_group(design)
    if feedback:
        figures["feedback"] = assets / "feedback-bracket.png"
        feedback_groups = _feedback_groups_by_host(design)
        _plot_feedback(
            endpoints,
            feedback_groups if feedback_groups else [feedback],
            figures["feedback"],
        )
    extension_groups = _extension_return_groups(design)
    if extension_groups:
        figures["extension"] = assets / "extension-return-comparisons.png"
        _plot_extension_return_groups(endpoints, extension_groups, figures["extension"])
    return figures


def _design_table_rows(design: list[dict[str, str]]) -> list[list[str]]:
    header = [
        "Cell",
        "Main purpose",
        "# hosts H",
        "Escape fraction f",
        "Escapees/host e",
        "Total return R",
        "Host fraction alpha",
        "Mutation rate u",
    ]
    has_migration = any(row.get("m", "") != "" for row in design)
    if has_migration:
        header.append("Migration m")
    rows = [header]
    for row in design:
        cell = _short_cell(_design_cell_id(row)).lower()
        values = [
            cell,
            SLIDE28_MAIN_PURPOSES.get(
                cell, row.get("label", row.get("experimental_group", ""))
            ),
            _format_integer(_number(row.get("H"))),
            _format_number(_number(row.get("f"))),
            _format_integer(_number(row.get("e"))),
            _format_number(_number(row.get("R"))),
            _format_number(_number(row.get("alpha"))),
            _format_number(_number(row.get("u"))),
        ]
        if has_migration:
            values.append(_format_number(_number(row.get("m"))))
        rows.append(values)
    return rows


def _metric_table_rows(
    endpoints: list[dict[str, str]], design: list[dict[str, str]]
) -> list[list[str]]:
    rows = [["Cell", "D0", "Shannon H'", "Simpson", "D1", "D2", "Evenness", "TV"]]
    for cell_id in _cell_order(design):
        rows.append(
            [
                _short_cell(cell_id),
                _format_integer(_cell_median(endpoints, cell_id, "richness")),
                _format_number(_cell_median(endpoints, cell_id, "shannon"), 4),
                _format_number(_cell_median(endpoints, cell_id, "simpson"), 5),
                _format_number(_cell_median(endpoints, cell_id, "hill_q1"), 3),
                _format_number(_cell_median(endpoints, cell_id, "hill_q2"), 3),
                _format_number(_cell_median(endpoints, cell_id, "evenness"), 4),
                _format_number(_cell_median(endpoints, cell_id, "total_variation"), 4),
            ]
        )
    return rows


def _cell_label_list(rows: list[dict[str, str]]) -> str:
    return ", ".join(_short_cell(_design_cell_id(row)) for row in rows)


def _effect_summary_rows(
    endpoints: list[dict[str, str]], design: list[dict[str, str]]
) -> list[list[str]]:
    """Build cautious biological summaries from the detected comparison series."""

    rows = [
        [
            "Comparison",
            "Diversity / composition",
            "Richness D0",
            "Evenness",
            "Main interpretation",
        ]
    ]
    migration_fraction = _migration_fraction(design)
    no_return = [row for row in design if _number(row.get("R")) == 0]
    if no_return:
        if migration_fraction > 0:
            no_return_ids = [_design_cell_id(row) for row in no_return]
            median_tv = _median(
                [
                    _number(endpoint["total_variation"])
                    for cell_id in no_return_ids
                    for endpoint in _cell_rows(endpoints, cell_id)
                ]
            )
            richness_values = [
                _cell_median(endpoints, cell_id, "richness")
                for cell_id in no_return_ids
            ]
            evenness_values = [
                _cell_median(endpoints, cell_id, "evenness")
                for cell_id in no_return_ids
            ]
            rows.append(
                [
                    f"Migration-only baseline ({_cell_label_list(no_return)})",
                    f"Median composition distance {median_tv:.4f}",
                    "Median D0 "
                    + "-".join(
                        _format_integer(value) for value in sorted(set(richness_values))
                    ),
                    "Median "
                    + "-".join(
                        f"{value:.3f}" for value in sorted(set(evenness_values))
                    ),
                    "With f=0, changes arise from stochastic exchange with the fixed regional pool; within-host events cannot return.",
                ]
            )
        else:
            rows.append(
                [
                    f"No return ({_cell_label_list(no_return)})",
                    "No change",
                    "No change",
                    "No change",
                    "With f=0, within-host events cannot affect the environment.",
                ]
            )

    matched = _matched_return_group(design)
    if matched:
        smallest = matched[0]
        largest = matched[-1]
        tv_low = _cell_median(endpoints, _design_cell_id(smallest), "total_variation")
        tv_high = _cell_median(endpoints, _design_cell_id(largest), "total_variation")
        rows.append(
            [
                f"Host return, u=0 ({_short_cell(_design_cell_id(smallest))})",
                f"Redistribution (median TV {tv_low:.3f}); abundance-weighted change varies by seed block.",
                "No change",
                "Median decrease at H=100; seed-block dependent.",
                "A small number of independent host samples can amplify an unrepresentative infection pool.",
            ]
        )
        rows.append(
            [
                f"Increase H, fixed R, u=0 ({_cell_label_list(matched)})",
                f"Redistribution falls (median TV {tv_low:.3f} to {tv_high:.3f}).",
                "No change in this mutation-free series",
                "Returns toward the initial value in two of three seed blocks.",
                "At R=1e9, H controls averaging across independent infections; u was fixed at zero.",
            ]
        )

    mutation = _mutation_group(design)
    if mutation:
        baseline = mutation[0]
        highest = mutation[-1]
        richness_start = _cell_median(endpoints, _design_cell_id(baseline), "richness")
        richness_end = _cell_median(endpoints, _design_cell_id(highest), "richness")
        rows.append(
            [
                f"Increase u ({_cell_label_list(mutation)})",
                "Composition changes; median D1 and D2 increase, but direction is seed-block dependent.",
                f"Strong increase ({_format_integer(richness_start)} to {_format_integer(richness_end)})",
                "Decreases as rare mutant labels accumulate",
                "The u dose-response was measured at H=100; the extension adds an exploratory H by u comparison.",
            ]
        )

    feedback_groups = _feedback_groups_by_host(design)
    if feedback_groups:
        for group in feedback_groups:
            host_number = int(_number(group[0].get("H")))
            first_id = _design_cell_id(group[0])
            last_id = _design_cell_id(group[-1])
            first_tv = _cell_median(endpoints, first_id, "total_variation")
            last_tv = _cell_median(endpoints, last_id, "total_variation")
            first_d1 = _cell_median(endpoints, first_id, "hill_q1")
            last_d1 = _cell_median(endpoints, last_id, "hill_q1")
            richness_values = [
                _cell_median(endpoints, _design_cell_id(row), "richness")
                for row in group
            ]
            first_evenness = _cell_median(endpoints, first_id, "evenness")
            last_evenness = _cell_median(endpoints, last_id, "evenness")
            sensitivity = (
                "The alpha=0.909 endpoint uses f=0.1, outside the primary range."
                if host_number == 100
                else "The alpha=0.001 endpoint uses f=1e-6, outside the primary range."
            )
            rows.append(
                [
                    f"Increase alpha, H={host_number:,} ({_cell_label_list(group)})",
                    "Redistribution increases overall (median TV "
                    f"{_format_number(first_tv, 4)} to {_format_number(last_tv, 4)}); "
                    f"median D1 {first_d1:.1f} to {last_d1:.1f}.",
                    (
                        f"No change (median D0 {_format_integer(richness_values[0])} "
                        "throughout)."
                        if min(richness_values) == max(richness_values)
                        else f"Median D0 {_format_integer(min(richness_values))}-"
                        f"{_format_integer(max(richness_values))}."
                    ),
                    f"Median {first_evenness:.3f} to {last_evenness:.3f}.",
                    f"H and u are fixed; f, R and alpha increase together. {sensitivity}",
                ]
            )
    else:
        feedback = _feedback_group(design)
        if feedback:
            rows.append(
                [
                    f"Feedback bracket ({_cell_label_list(feedback)})",
                    "Weak feedback is negligible; stronger cells redistribute composition, but not monotonically.",
                    "No change in mutation-free cells",
                    "Small and seed-block dependent changes",
                    "Exploratory bracket only: alpha, H and R are confounded.",
                ]
            )
    for group in _extension_return_groups(design):
        first = group[0]
        last = group[-1]
        first_id = _design_cell_id(first)
        last_id = _design_cell_id(last)
        return_size = _number(first.get("R"))
        mutation = _number(first.get("u"))
        first_tv = _cell_median(endpoints, first_id, "total_variation")
        last_tv = _cell_median(endpoints, last_id, "total_variation")
        first_richness = _cell_median(endpoints, first_id, "richness")
        last_richness = _cell_median(endpoints, last_id, "richness")
        first_evenness = _cell_median(endpoints, first_id, "evenness")
        last_evenness = _cell_median(endpoints, last_id, "evenness")
        richness_values = [
            _cell_median(endpoints, _design_cell_id(row), "richness") for row in group
        ]
        evenness_values = [
            _cell_median(endpoints, _design_cell_id(row), "evenness") for row in group
        ]
        mutation_enabled = mutation > 0
        label = (
            "Mutation-enabled pooling" if mutation_enabled else "Lower-return pooling"
        )
        if mutation_enabled:
            richness_result = (
                "Moderate, non-monotonic increase "
                f"(median D0 {_format_integer(min(richness_values))}-"
                f"{_format_integer(max(richness_values))}; "
                f"{_format_integer(last_richness)} at highest H)."
            )
            evenness_result = (
                "Small, non-monotonic change "
                f"(median {min(evenness_values):.3f}-{max(evenness_values):.3f})."
            )
            interpretation = (
                "Alongside the mutation-free R=1e9 series, this explores H by u: "
                "redistribution falls in both, while rare mutant labels remain."
            )
        else:
            richness_result = (
                f"No change (median D0 {_format_integer(first_richness)} throughout)."
                if min(richness_values) == max(richness_values)
                else f"Median D0 {_format_integer(first_richness)} to "
                f"{_format_integer(last_richness)}."
            )
            evenness_result = (
                "Negligible net change "
                f"(median {first_evenness:.3f} to {last_evenness:.3f})."
            )
            interpretation = (
                f"At R={return_size:.0e} and u=0, pooling across more hosts still "
                f"reduces redistribution; {len(group)} H levels were tested."
            )
        rows.append(
            [
                f"{label} ({_cell_label_list(group)})",
                f"Redistribution falls (median TV {first_tv:.3f} to {last_tv:.3f}).",
                richness_result,
                evenness_result,
                interpretation,
            ]
        )
    return rows


def _resource_scaling_rows(
    endpoints: list[dict[str, str]], design: list[dict[str, str]]
) -> list[list[str]]:
    rows = [["# hosts H", "Median runtime", "Median output", "Median peak RAM"]]
    for design_row in _matched_return_group(design):
        cell_id = _design_cell_id(design_row)
        rows.append(
            [
                _format_integer(_number(design_row.get("H"))),
                f"{_cell_median(endpoints, cell_id, 'elapsed_minutes'):.2f} min",
                f"{_cell_median(endpoints, cell_id, 'output_mib'):.2f} MiB",
                f"{_cell_median(endpoints, cell_id, 'peak_rss_mib'):.1f} MiB",
            ]
        )
    return rows


def _resource_scaling_summary(
    endpoints: list[dict[str, str]], design: list[dict[str, str]]
) -> dict[str, float]:
    matched = _matched_return_group(design)
    hosts = [_number(row.get("H")) for row in matched]
    summary: dict[str, float] = {}
    for measure in ("elapsed_minutes", "output_mib", "peak_rss_mib"):
        medians = [
            _cell_median(endpoints, _design_cell_id(row), measure) for row in matched
        ]
        summary[f"{measure}_slope"] = _log_log_slope(hosts, medians)
        summary[f"{measure}_tenfold"] = 10 ** summary[f"{measure}_slope"]
    mutation = _mutation_group(design)
    if mutation:
        baseline = _design_cell_id(mutation[0])
        highest = _design_cell_id(mutation[-1])
        for measure in ("elapsed_minutes", "output_mib", "peak_rss_mib"):
            summary[f"{measure}_mutation_ratio"] = _cell_median(
                endpoints, highest, measure
            ) / _cell_median(endpoints, baseline, measure)
    return summary


def _key_findings(
    summary: dict[str, Any],
    endpoints: list[dict[str, str]],
    design: list[dict[str, str]],
) -> list[str]:
    findings: list[str] = []
    migration_fraction = _migration_fraction(design)
    if migration_fraction > 0:
        findings.append(
            f"Every cell used fixed regional-pool exchange at m={migration_fraction:.2f}. "
            f"This replaces {100 * migration_fraction:.0f}% of the focal environment "
            "after host return, so the expected host-derived fraction carried into the "
            f"next passage is multiplied by {1.0 - migration_fraction:.2f}."
        )
    no_return = [row for row in design if _number(row.get("R")) == 0]
    if no_return:
        maximum_tv = max(
            _number(endpoint["total_variation"])
            for row in no_return
            for endpoint in _cell_rows(endpoints, _design_cell_id(row))
        )
        if migration_fraction > 0:
            median_tv = _median(
                [
                    _number(endpoint["total_variation"])
                    for row in no_return
                    for endpoint in _cell_rows(endpoints, _design_cell_id(row))
                ]
            )
            findings.append(
                "No-return controls provide the migration-only baseline: their median "
                f"composition distance was {_format_number(median_tv, 4)} after five "
                "passages. Host-return effects should be interpreted relative to this "
                "background exchange."
            )
        elif maximum_tv == 0:
            findings.append(
                "No-return controls left the environmental composition exactly unchanged."
            )
    matched = _matched_return_group(design)
    if matched:
        smallest = matched[0]
        largest = matched[-1]
        findings.append(
            "At fixed total return, increasing host number from "
            f"{_format_integer(_number(smallest.get('H')))} to "
            f"{_format_integer(_number(largest.get('H')))} changed the median "
            "composition distance from "
            f"{_format_number(_cell_median(endpoints, _design_cell_id(smallest), 'total_variation'), 4)} "
            "to "
            f"{_format_number(_cell_median(endpoints, _design_cell_id(largest), 'total_variation'), 4)}. "
            "Total return alone therefore did not determine the short-term response."
        )
    mutation = _mutation_group(design)
    positive = [row for row in mutation if _number(row.get("u")) > 0]
    if positive:
        highest = positive[-1]
        cell_id = _design_cell_id(highest)
        findings.append(
            f"At the highest mutation level ({_number(highest.get('u')):.0e}), "
            f"median labelled richness was {_format_integer(_cell_median(endpoints, cell_id, 'richness'))}, "
            "but the median mutant abundance fraction was only "
            f"{_cell_median(endpoints, cell_id, 'mutant_abundance_fraction'):.2e}. "
            "Mutation created many rare labels rather than abundant new types."
        )
    feedback_groups = _feedback_groups_by_host(design)
    if feedback_groups:
        changes = []
        for group in feedback_groups:
            first = _design_cell_id(group[0])
            last = _design_cell_id(group[-1])
            changes.append(
                f"H={_format_integer(_number(group[0].get('H')))}: median TV "
                f"{_format_number(_cell_median(endpoints, first, 'total_variation'), 4)} to "
                f"{_format_number(_cell_median(endpoints, last, 'total_variation'), 4)}, median D1 "
                f"{_cell_median(endpoints, first, 'hill_q1'):.1f} to "
                f"{_cell_median(endpoints, last, 'hill_q1'):.1f}"
            )
        findings.append(
            "Across the fixed-H feedback series, increasing alpha from about 0.001 "
            "to 0.909 produced these endpoint changes (" + "; ".join(changes) + ")."
        )
    return findings


def _relative_asset(markdown: Path, asset: Path) -> str:
    import os

    return Path(os.path.relpath(asset, markdown.parent)).as_posix()


def _markdown_table(rows: list[list[str]]) -> str:
    widths = len(rows[0])
    lines = ["| " + " | ".join(rows[0]) + " |"]
    lines.append("|" + "|".join(["---"] * widths) + "|")
    lines.extend("| " + " | ".join(row) + " |" for row in rows[1:])
    return "\n".join(lines)


def _write_markdown(
    path: Path,
    title: str,
    report_date: str,
    summary: dict[str, Any],
    endpoints: list[dict[str, str]],
    design: list[dict[str, str]],
    figures: dict[str, Path],
) -> None:
    resources = summary["overall_resources"]
    findings = _key_findings(summary, endpoints, design)
    matched = _matched_return_group(design)
    mutation = _mutation_group(design)
    feedback = _feedback_group(design)
    feedback_groups = _feedback_groups_by_host(design)
    extension_groups = _extension_return_groups(design)
    scaling = _resource_scaling_summary(endpoints, design)
    migration_fraction = _migration_fraction(design)
    migration_scope = (
        "All cells include fixed regional-pool exchange after host return and "
        f"environmental regulation (m={migration_fraction:.2f}). The no-return cells "
        "therefore measure migration alone, rather than an unchanged environment."
        if migration_fraction > 0
        else "The environmental reservoir has no regional migration step in this pilot."
    )
    lines = [
        f"# {title}",
        "",
        "**Status:** exploratory pilot report; not a confirmatory analysis  ",
        f"**Report date:** {report_date}  ",
        f"**Experiment:** `{summary.get('experiment_id', 'pilot')}`",
        "",
        "## Purpose and scope",
        "",
        "This pilot tested whether the simulation is computationally feasible and whether "
        "the selected parameter levels produce informative biological responses. One complete "
        "simulated population is the independent replicate. The pilot does not establish "
        "equilibrium or provide confirmatory tests.",
        "",
        migration_scope,
        "",
        "## Pilot design",
        "",
        _markdown_table(_design_table_rows(design)),
        "",
        "**Colour key in the PDF:** grey marks baseline/reference purposes; purple "
        "shades distinguish the host-number H comparisons; blue shades mark the "
        "feedback-alpha comparisons; and green marks the mutation-rate u comparison. "
        "Coloured text identifies cells that also belong to another comparison series. "
        "The wording and colour assignments reproduce slide 28 of the project update.",
        "",
        "## Biological results",
        "",
    ]
    lines.extend(f"- {finding}" for finding in findings)
    lines.extend(
        [
            "",
            "### Summary of biological effects",
            "",
            _markdown_table(_effect_summary_rows(endpoints, design)),
            "",
            "The three seed blocks are matched across cells. Comparisons therefore use "
            "within-seed changes: each coloured line in the comparison figures joins the "
            "same seed block across parameter levels. Black lines show cell medians. With "
            "only three blocks, these are descriptive paired patterns, not hypothesis tests.",
            "",
            f"![Endpoint diversity]({_relative_asset(path, figures['endpoint'])})",
            "",
            "### Endpoint diversity indices",
            "",
            "Values are cell medians across independent population replicates.",
            "",
            _markdown_table(_metric_table_rows(endpoints, design)),
            "",
        ]
    )
    if "migration" in figures:
        lines.extend(
            [
                "### Fixed regional-pool migration",
                "",
                f"![Fixed-pool migration]({_relative_asset(path, figures['migration'])})",
                "",
                "Panel A shows how exchange reduces the expected host-derived fraction "
                "after return and regulation. Panel B isolates the richness change caused "
                "by replacing emigrants with immigrants during passage 5. All cells use "
                "the same m, so this pilot controls for migration but does not estimate its "
                "dose-response.",
                "",
            ]
        )
    if matched and "matched" in figures:
        figure2_groups = _figure2_pooling_groups(design)
        figure2_alphas = {_number(group[0].get("alpha")) for group in figure2_groups}
        figure2_mutations = {_number(group[0].get("u")) for group in figure2_groups}
        if len(figure2_alphas) > 1:
            matched_caption = (
                "Colours identify matched seed blocks. Circles show u=0 and squares "
                "show u=1e-10. Filled markers with solid lines show alpha=0.5; hollow "
                "markers with dashed lines show alpha=0.09. The thicker black "
                "underlays show cell medians."
            )
        elif len(figure2_mutations) > 1:
            matched_caption = (
                "Colours identify matched seed blocks. Circles show u=0 and squares "
                "show u=1e-10. The thicker black underlays show cell medians."
            )
        else:
            matched_caption = (
                "Coloured lines connect the same seed block across H; the black line "
                "is the cell median."
            )
        lines.extend(
            [
                "### Host pooling within fixed-return series",
                "",
                f"![Matched return]({_relative_asset(path, figures['matched'])})",
                "",
                matched_caption,
                "",
            ]
        )
    if mutation and "mutation" in figures:
        lines.extend(
            [
                "### Within-host mutation",
                "",
                f"![Mutation bracket]({_relative_asset(path, figures['mutation'])})",
                "",
                "Coloured lines connect the same seed block across u; the black line is "
                "the cell median.",
                "",
            ]
        )
    if feedback and "feedback" in figures:
        if feedback_groups:
            feedback_heading = "### Feedback strength at fixed host abundance"
            feedback_text = (
                "The same four alpha levels are evaluated separately at H=100 and "
                "H=1,000. Within each series H and u are fixed, while f, total return R "
                "and alpha increase together. Lines connect the same seed block. The "
                "extreme f values used for c0018 and c0019 are sensitivity levels outside "
                "the primary escape-rate range."
            )
        else:
            feedback_heading = "### Broad feedback bracket"
            feedback_text = (
                "This is not an isolated alpha comparison because H and total return also "
                "change. It is retained as an exploratory boundary check."
            )
        lines.extend(
            [
                feedback_heading,
                "",
                feedback_text,
                "",
                f"![Feedback bracket]({_relative_asset(path, figures['feedback'])})",
                "",
            ]
        )
    if extension_groups and "extension" in figures:
        lines.extend(
            [
                "### Extension pooling comparisons",
                "",
                "The extension completes the mutation-enabled R=10^9 series and the "
                "mutation-free R=10^8 series. Lines connect the same seed block.",
                "",
                f"![Extension comparisons]({_relative_asset(path, figures['extension'])})",
                "",
            ]
        )
    lines.extend(
        [
            "## Quality control and computational feasibility",
            "",
            f"The analysis included {summary.get('populations')} populations across "
            f"{summary.get('cells')} cells. The endpoint audit result was "
            f"**{summary.get('audit_status')}**. Summed runtime was "
            f"{resources['elapsed_hours']:.2f} hours, diagnostic output was "
            f"{resources['output_mib']:.1f} MiB, and maximum measured process-tree memory "
            f"was {resources['peak_rss_mib']:.1f} MiB. Every pilot run used two worker "
            "processes, five host passages and 500 steady within-host bacterial generations.",
            "",
            *(
                [
                    f"Each passage exchanged {100 * migration_fraction:.0f}% of the "
                    "one-billion-cell focal environment with the fixed regional source. "
                    "The completion audit verified equal emigrant and immigrant totals "
                    "for every population and passage.",
                    "",
                ]
                if migration_fraction > 0
                else []
            ),
            f"Benchmark machine: {summary.get('benchmark_hardware_note', 'not recorded')}. "
            "Timings should be recalibrated with a small batch on the target HPC before "
            "the larger experiment is launched.",
            "",
            f"![Computational feasibility]({_relative_asset(path, figures['feasibility'])})",
            "",
            _markdown_table(_resource_scaling_rows(endpoints, design)),
            "",
            f"Within the mutation-free matched-return series, the fitted host-number exponent "
            f"was {scaling['elapsed_minutes_slope']:.2f} for runtime and "
            f"{scaling['output_mib_slope']:.2f} for diagnostic output. Thus, a tenfold increase in H multiplied "
            f"median runtime by about {scaling['elapsed_minutes_tenfold']:.1f} and output by "
            f"about {scaling['output_mib_tenfold']:.1f}. Peak RAM was nearly independent of "
            "H because host tasks are streamed through a bounded queue.",
            "",
            f"At H=100, raising u from zero to the highest pilot level multiplied median "
            f"runtime by {scaling['elapsed_minutes_mutation_ratio']:.1f}, output by "
            f"{scaling['output_mib_mutation_ratio']:.1f}, and peak RAM by "
            f"{scaling['peak_rss_mib_mutation_ratio']:.2f}. This cost comes from materialized "
            "mutant lineages and their records.",
            "",
            "For this implementation, a useful work model is O(P H [G S + M]), where P is "
            "the number of host passages, G is the number of within-host transitions, S is "
            "extant strain richness and M is the number of materialized mutation events. "
            "The code works with strain counts, not individual bacterial cells. Therefore "
            "K, e and R do not create one operation per cell. The independent effect of e "
            "cannot be estimated from this pilot because e changes inversely with H in the "
            "matched-return series.",
            "",
            "**HPC starting allocation:** request two CPUs and 1 GiB RAM for each concurrent "
            "population. Begin with 8-16 concurrent populations, inspect CPU use and file-system "
            "throughput, then increase toward the CPU limit. With 100 available CPUs and two "
            "workers per population, 50 concurrent populations is the CPU ceiling. Rebenchmark "
            "before increasing workers, u, retained host detail, or the number of passages. "
            "At unchanged output settings, runtime and host-level output should grow roughly "
            "linearly with the number of host passages.",
            "",
            "## Interpretation and next step",
            "",
            "The pilot supports using a longer equilibrium and precision pilot before any "
            "confirmatory experiment. Three replicates per cell are sufficient for calibration "
            "and graphical ranges, but not for stable confidence intervals or equivalence tests.",
            "",
            "## Reproducibility",
            "",
            f"- Source commit: `{', '.join(summary.get('software_git_commits', []))}`",
            f"- Recorded source-tree SHA-256 values: "
            f"`{', '.join(summary.get('source_sha256', []))}`",
            f"- Seed blocks: `{', '.join(summary.get('seed_blocks', []))}`",
            f"- Analysis audit: `{summary.get('audit_status')}`",
            *(
                [
                    "- Multiple source-tree hashes reflect staged reporting, launcher and "
                    "pilot-design additions. The biological transition modules were unchanged; "
                    "each population records the exact hash and resolved configuration used."
                ]
                if len(summary.get("source_sha256", [])) > 1
                else []
            ),
            "- Complete numerical results remain in the machine-readable analysis tables.",
            "",
            "## Appendix: glossary and diversity measures",
            "",
            "- **Relative abundance, pᵢ:** fraction of environmental cells assigned to strain i.",
            "- **Labelled richness, D0:** number of distinct strain labels. Every label counts "
            "once, even if it contains only one cell.",
            "- **Shannon entropy, H':** H' = -Σᵢ pᵢ ln(pᵢ). Higher values mean that strain "
            "identity is less predictable. It is an index, not an effective number of strains.",
            "- **Simpson diversity:** 1 - Σᵢ pᵢ². This is the probability that two "
            "independently sampled cells have different strain labels.",
            "- **Hill D1:** D₁ = exp(H'). The number of equally abundant strains that would give "
            "the observed Shannon entropy.",
            "- **Hill D2:** D₂ = 1 / Σᵢ pᵢ². The number of equally abundant strains that would "
            "give the observed probability of sampling the same strain twice. It gives more "
            "weight to abundant strains than D1.",
            "- **Pielou evenness:** J = H' / ln(D₀). It compares the observed Shannon entropy "
            "with the maximum possible value for the observed richness.",
            "- **Composition distance (TV):** TV(A,B) = ½ Σᵢ |pᵢ(A) - pᵢ(B)|. A value of 0.20 means "
            "that 20% of abundance must be reassigned among labels to make two compositions "
            "identical.",
            "- **Cell:** one defined combination of experimental parameter values.",
            "- **Seed block:** one independently seeded population replicate reused as a "
            "reproducible block across cells.",
            "- **H:** number of hosts; **f:** escape fraction per host; **e:** escaping cells "
            "per host (fK); **R:** total returned cells (He); **alpha:** fraction of the "
            "pre-regulation mixture derived from hosts; **u:** whole-genome mutation "
            "probability per bacterial generation; **m:** fraction of the focal "
            "environment replaced from the fixed regional pool after each host return.",
            "",
        ]
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines), encoding="utf-8")


def _reportlab_imports() -> dict[str, Any]:
    try:
        from reportlab.lib import colors as report_colors
        from reportlab.lib.pagesizes import A4 as report_a4
        from reportlab.lib.styles import ParagraphStyle as ReportParagraphStyle
        from reportlab.lib.styles import getSampleStyleSheet as report_styles
        from reportlab.lib.units import mm as report_mm
        from reportlab.pdfbase import pdfmetrics as report_pdfmetrics
        from reportlab.pdfbase.ttfonts import TTFont as ReportTTFont
        from reportlab.platypus import (
            Image as ReportImage,
        )
        from reportlab.platypus import (
            KeepTogether as ReportKeepTogether,
        )
        from reportlab.platypus import (
            PageBreak as ReportPageBreak,
        )
        from reportlab.platypus import (
            Paragraph as ReportParagraph,
        )
        from reportlab.platypus import (
            SimpleDocTemplate as ReportDocument,
        )
        from reportlab.platypus import (
            Spacer as ReportSpacer,
        )
        from reportlab.platypus import (
            Table as ReportTable,
        )
        from reportlab.platypus import (
            TableStyle as ReportTableStyle,
        )
    except ImportError as error:
        raise RuntimeError(
            "PDF reporting requires the optional report dependencies. "
            "Install with: pip install 'trophosome-model[report]'"
        ) from error
    return {
        "colors": report_colors,
        "A4": report_a4,
        "ParagraphStyle": ReportParagraphStyle,
        "getSampleStyleSheet": report_styles,
        "mm": report_mm,
        "pdfmetrics": report_pdfmetrics,
        "TTFont": ReportTTFont,
        "Image": ReportImage,
        "KeepTogether": ReportKeepTogether,
        "PageBreak": ReportPageBreak,
        "Paragraph": ReportParagraph,
        "SimpleDocTemplate": ReportDocument,
        "Spacer": ReportSpacer,
        "Table": ReportTable,
        "TableStyle": ReportTableStyle,
    }


def _paragraph(text: str, style: Any, paragraph: Any) -> Any:
    return paragraph(html.escape(text).replace("\n", "<br/>"), style)


def _pdf_table(
    rows: list[list[str]],
    widths: list[float],
    imports: dict[str, Any],
    font_size: float = 7.5,
    *,
    extra_style: list[tuple[Any, ...]] | None = None,
    font_name: str = "Helvetica",
    paragraph_text_colours: dict[tuple[int, int], str] | None = None,
) -> Any:
    colors = imports["colors"]
    Paragraph = imports["Paragraph"]
    ParagraphStyle = imports["ParagraphStyle"]
    Table = imports["Table"]
    TableStyle = imports["TableStyle"]
    cell_style = ParagraphStyle(
        "table-cell",
        fontName=font_name,
        fontSize=font_size,
        leading=font_size + 2,
    )
    paragraph_text_colours = paragraph_text_colours or {}
    table_data = []
    for row_index, row in enumerate(rows):
        table_row = []
        for column_index, value in enumerate(row):
            text_colour = paragraph_text_colours.get(
                (column_index, row_index),
                "#16394F" if row_index == 0 else "#000000",
            )
            paragraph_style = ParagraphStyle(
                f"table-cell-{row_index}-{column_index}",
                parent=cell_style,
                fontName="Helvetica-Bold" if row_index == 0 else font_name,
                textColor=colors.HexColor(text_colour),
            )
            table_row.append(Paragraph(html.escape(value), paragraph_style))
        table_data.append(table_row)
    table = Table(table_data, colWidths=widths, repeatRows=1, hAlign="LEFT")
    commands: list[tuple[Any, ...]] = [
        ("BACKGROUND", (0, 0), (-1, 0), colors.HexColor("#DDEBF3")),
        ("TEXTCOLOR", (0, 0), (-1, 0), colors.HexColor("#16394F")),
        ("FONTNAME", (0, 0), (-1, 0), "Helvetica-Bold"),
        ("GRID", (0, 0), (-1, -1), 0.35, colors.HexColor("#BFCBD3")),
        ("VALIGN", (0, 0), (-1, -1), "TOP"),
        ("LEFTPADDING", (0, 0), (-1, -1), 4),
        ("RIGHTPADDING", (0, 0), (-1, -1), 4),
        ("TOPPADDING", (0, 0), (-1, -1), 4),
        ("BOTTOMPADDING", (0, 0), (-1, -1), 4),
    ]
    if extra_style:
        commands.extend(extra_style)
    table.setStyle(TableStyle(commands))
    return table


def _design_paragraph_text_colours(
    design: list[dict[str, str]],
) -> dict[tuple[int, int], str]:
    column_count = len(_design_table_rows(design)[0])
    colours = {(column, 0): SLIDE28_PALETTE["white"] for column in range(column_count)}
    colours[(2, 0)] = SLIDE28_PALETTE["purple_light"]
    colours[(6, 0)] = SLIDE28_PALETTE["blue_light"]
    colours[(7, 0)] = SLIDE28_PALETTE["green_light"]

    host_dark = {"c0009", "c0013", "c0014", "c0015"}
    host_blue_text = {"c0004", "c0012", "c0019", "c0020"}
    feedback_purple_text = {"c0016", "c0017", "c0020"}
    mutation_purple_text = {"c0009", "c0013", "c0014", "c0015"}
    for row_index, design_row in enumerate(design, start=1):
        cell_id = _short_cell(_design_cell_id(design_row)).lower()
        if cell_id in host_dark:
            colours[(2, row_index)] = SLIDE28_PALETTE["white"]
        if cell_id in host_blue_text:
            colours[(2, row_index)] = SLIDE28_PALETTE["blue_medium"]
        if cell_id in feedback_purple_text:
            colours[(6, row_index)] = SLIDE28_PALETTE["purple_medium"]
        if cell_id in mutation_purple_text:
            colours[(7, row_index)] = SLIDE28_PALETTE["purple_dark"]
    return colours


def _design_highlight_style(
    design: list[dict[str, str]], imports: dict[str, Any]
) -> list[tuple[Any, ...]]:
    colors = imports["colors"]
    palette = {name: colors.HexColor(value) for name, value in SLIDE28_PALETTE.items()}
    commands: list[tuple[Any, ...]] = [
        ("BACKGROUND", (0, 0), (-1, 0), palette["black"]),
        ("TEXTCOLOR", (0, 0), (-1, 0), palette["white"]),
        ("TEXTCOLOR", (2, 0), (2, 0), palette["purple_light"]),
        ("TEXTCOLOR", (6, 0), (6, 0), palette["blue_light"]),
        ("TEXTCOLOR", (7, 0), (7, 0), palette["green_light"]),
    ]
    grey_purpose = {"c0001", "c0002", "c0003"}
    host_light = {"c0003", "c0004", "c0005", "c0006"}
    host_dark = {"c0009", "c0013", "c0014", "c0015"}
    host_medium = {"c0016", "c0017", "c0020"}
    host_blue_text = {"c0004", "c0012", "c0019", "c0020"}
    feedback_blue_grey = {"c0003", "c0011"}
    feedback_light = {"c0016", "c0018"}
    feedback_medium = {"c0004", "c0012", "c0020"}
    feedback_purple_text = {"c0016", "c0017", "c0020"}
    mutation_green = {"c0003", "c0007", "c0008", "c0009", "c0010"}
    mutation_purple_text = {"c0009", "c0013", "c0014", "c0015"}

    for row_index, design_row in enumerate(design, start=1):
        cell_id = _short_cell(_design_cell_id(design_row)).lower()
        if cell_id in grey_purpose:
            commands.append(
                ("BACKGROUND", (1, row_index), (1, row_index), palette["grey"])
            )
        if cell_id in host_light:
            commands.append(
                (
                    "BACKGROUND",
                    (2, row_index),
                    (2, row_index),
                    palette["purple_light"],
                )
            )
        if cell_id in host_dark:
            commands.append(
                (
                    "BACKGROUND",
                    (2, row_index),
                    (2, row_index),
                    palette["purple_dark"],
                )
            )
            commands.append(
                ("TEXTCOLOR", (2, row_index), (2, row_index), palette["white"])
            )
        if cell_id in host_medium:
            commands.append(
                (
                    "BACKGROUND",
                    (2, row_index),
                    (2, row_index),
                    palette["purple_medium"],
                )
            )
        if cell_id in host_blue_text:
            commands.append(
                (
                    "TEXTCOLOR",
                    (2, row_index),
                    (2, row_index),
                    palette["blue_medium"],
                )
            )
        if cell_id in feedback_blue_grey:
            commands.append(
                (
                    "BACKGROUND",
                    (6, row_index),
                    (6, row_index),
                    palette["blue_grey"],
                )
            )
        if cell_id in feedback_light:
            commands.append(
                (
                    "BACKGROUND",
                    (6, row_index),
                    (6, row_index),
                    palette["blue_light"],
                )
            )
        if cell_id in feedback_medium:
            commands.append(
                (
                    "BACKGROUND",
                    (6, row_index),
                    (6, row_index),
                    palette["blue_medium"],
                )
            )
        if cell_id in feedback_purple_text:
            commands.append(
                (
                    "TEXTCOLOR",
                    (6, row_index),
                    (6, row_index),
                    palette["purple_medium"],
                )
            )
        if cell_id in mutation_green:
            commands.append(
                (
                    "BACKGROUND",
                    (7, row_index),
                    (7, row_index),
                    palette["green_light"],
                )
            )
        if cell_id in mutation_purple_text:
            commands.append(
                (
                    "TEXTCOLOR",
                    (7, row_index),
                    (7, row_index),
                    palette["purple_dark"],
                )
            )
    return commands


def _design_legend(
    _design: list[dict[str, str]], imports: dict[str, Any], font_size: float = 6.5
) -> Any:
    colors = imports["colors"]
    Paragraph = imports["Paragraph"]
    ParagraphStyle = imports["ParagraphStyle"]
    Table = imports["Table"]
    TableStyle = imports["TableStyle"]
    base_style = ParagraphStyle(
        "legend-cell",
        fontName="Helvetica",
        fontSize=font_size,
        leading=font_size + 2,
        alignment=1,
    )
    legend_items = [
        ("Baseline / reference", SLIDE28_PALETTE["grey"], "#000000"),
        ("H: R=1e9, u=0", SLIDE28_PALETTE["purple_light"], "#000000"),
        ("H: R=1e9, u=1e-10", SLIDE28_PALETTE["purple_dark"], "#FFFFFF"),
        ("H: R=1e8, u=0", SLIDE28_PALETTE["purple_medium"], "#000000"),
        ("Feedback alpha: blue shades", SLIDE28_PALETTE["blue_light"], "#000000"),
        ("Mutation u: green", SLIDE28_PALETTE["green_light"], "#000000"),
    ]
    labels = [label for label, _fill, _text in legend_items]
    fills = [fill for _label, fill, _text in legend_items]
    text_colours = [text for _label, _fill, text in legend_items]
    cell_width = 464 / len(labels)
    table = Table(
        [
            [
                Paragraph(
                    label,
                    ParagraphStyle(
                        f"legend-cell-{index}",
                        parent=base_style,
                        textColor=colors.HexColor(text_colours[index]),
                    ),
                )
                for index, label in enumerate(labels)
            ]
        ],
        colWidths=[cell_width] * len(labels),
        hAlign="LEFT",
    )
    commands: list[tuple[Any, ...]] = [
        ("GRID", (0, 0), (-1, -1), 0.35, colors.HexColor("#BFCBD3")),
        ("VALIGN", (0, 0), (-1, -1), "MIDDLE"),
        ("TOPPADDING", (0, 0), (-1, -1), 4),
        ("BOTTOMPADDING", (0, 0), (-1, -1), 4),
    ]
    for column, fill in enumerate(fills):
        commands.append(("BACKGROUND", (column, 0), (column, 0), colors.HexColor(fill)))
        commands.append(
            (
                "TEXTCOLOR",
                (column, 0),
                (column, 0),
                colors.HexColor(text_colours[column]),
            )
        )
    table.setStyle(TableStyle(commands))
    return table


def _scaled_image(path: Path, width: float, imports: dict[str, Any]) -> Any:
    Image = imports["Image"]
    image = Image(str(path))
    ratio = image.imageHeight / image.imageWidth
    image.drawWidth = width
    image.drawHeight = width * ratio
    return image


def _write_pdf(
    path: Path,
    title: str,
    report_date: str,
    summary: dict[str, Any],
    endpoints: list[dict[str, str]],
    design: list[dict[str, str]],
    figures: dict[str, Path],
) -> None:
    imports = _reportlab_imports()
    colors = imports["colors"]
    A4 = imports["A4"]
    getSampleStyleSheet = imports["getSampleStyleSheet"]
    ParagraphStyle = imports["ParagraphStyle"]
    Paragraph = imports["Paragraph"]
    SimpleDocTemplate = imports["SimpleDocTemplate"]
    Spacer = imports["Spacer"]
    PageBreak = imports["PageBreak"]
    KeepTogether = imports["KeepTogether"]
    mm = imports["mm"]
    pdfmetrics = imports["pdfmetrics"]
    TTFont = imports["TTFont"]

    from matplotlib import font_manager

    if "DejaVuSans" not in pdfmetrics.getRegisteredFontNames():
        pdfmetrics.registerFont(
            TTFont("DejaVuSans", font_manager.findfont("DejaVu Sans"))
        )

    path.parent.mkdir(parents=True, exist_ok=True)
    document = SimpleDocTemplate(
        str(path),
        pagesize=A4,
        leftMargin=18 * mm,
        rightMargin=18 * mm,
        topMargin=18 * mm,
        bottomMargin=18 * mm,
        title=title,
        author="Trophosome model reporting workflow",
    )
    base = getSampleStyleSheet()
    styles = {
        "title": ParagraphStyle(
            "report-title",
            parent=base["Title"],
            fontName="Helvetica-Bold",
            fontSize=20,
            leading=24,
            textColor=colors.HexColor("#16394F"),
            spaceAfter=12,
        ),
        "subtitle": ParagraphStyle(
            "subtitle",
            parent=base["Normal"],
            fontSize=9,
            leading=13,
            textColor=colors.HexColor("#5F6B76"),
            spaceAfter=6,
        ),
        "h1": ParagraphStyle(
            "h1",
            parent=base["Heading1"],
            fontName="Helvetica-Bold",
            fontSize=14,
            leading=17,
            textColor=colors.HexColor("#16394F"),
            spaceBefore=12,
            spaceAfter=7,
            keepWithNext=True,
        ),
        "h2": ParagraphStyle(
            "h2",
            parent=base["Heading2"],
            fontName="Helvetica-Bold",
            fontSize=11,
            leading=14,
            textColor=colors.HexColor("#1C5D86"),
            spaceBefore=9,
            spaceAfter=5,
            keepWithNext=True,
        ),
        "body": ParagraphStyle(
            "body",
            parent=base["BodyText"],
            fontName="Helvetica",
            fontSize=9.2,
            leading=13.2,
            textColor=colors.HexColor("#20303B"),
            spaceAfter=7,
        ),
        "caption": ParagraphStyle(
            "caption",
            parent=base["BodyText"],
            fontName="Helvetica-Oblique",
            fontSize=7.8,
            leading=10.2,
            textColor=colors.HexColor("#5F6B76"),
            spaceAfter=8,
        ),
        "bullet": ParagraphStyle(
            "bullet",
            parent=base["BodyText"],
            fontName="Helvetica",
            fontSize=9.2,
            leading=13.2,
            leftIndent=12,
            firstLineIndent=-7,
            bulletIndent=0,
            spaceAfter=4,
        ),
    }

    story: list[Any] = []
    migration_fraction = _migration_fraction(design)
    story.append(_paragraph(title, styles["title"], Paragraph))
    story.append(
        _paragraph(
            "Exploratory pilot report - not a confirmatory analysis",
            styles["subtitle"],
            Paragraph,
        )
    )
    story.append(
        _paragraph(
            f"Report date: {report_date} | Experiment: "
            f"{summary.get('experiment_id', 'pilot')} | "
            f"{summary.get('cells')} cells | {summary.get('populations')} populations",
            styles["subtitle"],
            Paragraph,
        )
    )
    story.append(Spacer(1, 6))
    story.append(_paragraph("Purpose and scope", styles["h1"], Paragraph))
    story.append(
        _paragraph(
            "This pilot tested computational feasibility, early biological effects and "
            "whether the selected parameter levels bracket informative responses. One "
            "complete simulated population is the independent replicate. The pilot does "
            "not establish equilibrium and does not provide confirmatory hypothesis tests.",
            styles["body"],
            Paragraph,
        )
    )
    if migration_fraction > 0:
        story.append(
            _paragraph(
                "All cells include fixed regional-pool exchange after host return and "
                f"environmental regulation (m={migration_fraction:.2f}). The no-return "
                "cells therefore provide a migration-only baseline rather than an "
                "unchanging environment.",
                styles["body"],
                Paragraph,
            )
        )
    story.append(_paragraph("Pilot design", styles["h1"], Paragraph))
    design_widths = (
        [28, 91, 38, 43, 48, 43, 49, 50, 40]
        if migration_fraction > 0
        else [32, 112, 44, 51, 55, 51, 57, 61]
    )
    story.append(
        _pdf_table(
            _design_table_rows(design),
            design_widths,
            imports,
            font_size=5.6 if migration_fraction > 0 else 5.9,
            extra_style=_design_highlight_style(design, imports),
            paragraph_text_colours=_design_paragraph_text_colours(design),
        )
    )
    story.append(Spacer(1, 5))
    story.append(_design_legend(design, imports))
    story.append(Spacer(1, 4))
    story.append(
        _paragraph(
            "Purpose wording and colour assignments follow slide 28. Grey marks baseline "
            "or reference purposes; purple shades distinguish the H series; blue shades "
            "mark feedback; green marks mutation; and coloured text shows overlaps between "
            "comparison series. Here e=fK and R=He.",
            styles["caption"],
            Paragraph,
        )
    )
    story.append(PageBreak())
    story.append(_paragraph("Biological results", styles["h1"], Paragraph))
    for finding in _key_findings(summary, endpoints, design):
        story.append(
            Paragraph(
                "- " + html.escape(finding),
                styles["bullet"],
            )
        )
    story.append(_paragraph("Summary of biological effects", styles["h2"], Paragraph))
    story.append(
        _pdf_table(
            _effect_summary_rows(endpoints, design),
            [75, 105, 80, 86, 128],
            imports,
            font_size=6.3,
        )
    )
    story.append(Spacer(1, 6))
    story.append(
        _paragraph(
            "Seed blocks are matched across cells. Each coloured line in the comparison "
            "figures therefore joins the same stochastic population across parameter levels; "
            "the black line is the cell median. With only three blocks, these are descriptive "
            "paired patterns, not hypothesis tests.",
            styles["body"],
            Paragraph,
        )
    )
    story.append(_scaled_image(figures["endpoint"], 174 * mm, imports))
    story.append(
        _paragraph(
            "Figure 1. Endpoint responses. Dots are independent simulated populations; "
            "black bars are cell medians. Shannon entropy and Simpson diversity are "
            "shown separately from their effective-number transformations D1 and D2.",
            styles["caption"],
            Paragraph,
        )
    )
    if "migration" in figures:
        story.append(
            _paragraph("Fixed regional-pool migration", styles["h2"], Paragraph)
        )
        story.append(
            _paragraph(
                "Regional exchange occurs after host return and environmental capacity "
                "regulation. It both attenuates the host-derived contribution and can "
                "change richness through the stochastic replacement of focal emigrants "
                "with regional immigrants.",
                styles["body"],
                Paragraph,
            )
        )
        story.append(_scaled_image(figures["migration"], 172 * mm, imports))
        story.append(
            _paragraph(
                "Migration diagnostic. Panel A compares the expected host-derived "
                "fraction immediately before and after exchange at passage 5. Panel B "
                "shows the richness change caused specifically by that exchange; dots "
                "are seed blocks and black bars are cell medians.",
                styles["caption"],
                Paragraph,
            )
        )

    matched = _matched_return_group(design)
    if matched and "matched" in figures:
        figure2_groups = _figure2_pooling_groups(design)
        feedback_overlay = (
            len({_number(group[0].get("alpha")) for group in figure2_groups}) > 1
        )
        mutation_overlay = (
            len({_number(group[0].get("u")) for group in figure2_groups}) > 1
        )
        story.append(
            _paragraph(
                "Host pooling within fixed-return series", styles["h2"], Paragraph
            )
        )
        story.append(
            _paragraph(
                (
                    "Within each feedback level, total return is held fixed while it is "
                    "distributed across different numbers of hosts: R=1e9 at alpha=0.5 "
                    "and R=1e8 at alpha=0.09. Marker shape additionally compares u=0 "
                    "with u=1e-10 where both mutation levels were simulated."
                    if feedback_overlay
                    else "Both mutation series returned R=1e9 bacteria but distributed "
                    "that return across different numbers of hosts. Comparing u=0 with "
                    "u=1e-10 shows whether the averaging effect of many independent "
                    "infections changes when within-host mutation is active."
                    if mutation_overlay
                    else "These cells returned the same total number of bacteria but "
                    "distributed that return across different numbers of hosts. This "
                    "separates the averaging effect of many independent infections from "
                    "total return size."
                ),
                styles["body"],
                Paragraph,
            )
        )
        story.append(_scaled_image(figures["matched"], 168 * mm, imports))
        story.append(
            _paragraph(
                (
                    "Figure 2. Host pooling within fixed-return series. Colours identify "
                    "seed blocks; circles show u=0 and squares show u=1e-10. Filled "
                    "markers with solid lines show alpha=0.5; hollow markers with dashed "
                    "lines show alpha=0.09. Thicker black underlays show cell medians."
                    if feedback_overlay
                    else "Figure 2. Matched-return comparison at u=0 and u=1e-10. "
                    "Colours identify seed blocks; circles show u=0 and squares show "
                    "u=1e-10. Thicker black underlays show cell medians."
                    if mutation_overlay
                    else "Figure 2. Matched-return comparison. Coloured lines join the "
                    "same seed block across H; the black line shows cell medians."
                ),
                styles["caption"],
                Paragraph,
            )
        )

    mutation = _mutation_group(design)
    if mutation and "mutation" in figures:
        story.append(_paragraph("Within-host mutation", styles["h2"], Paragraph))
        story.append(
            _paragraph(
                "Mutation adds new strain labels. Richness must therefore be interpreted "
                "with abundance-weighted indices and the total abundance carried by mutants.",
                styles["body"],
                Paragraph,
            )
        )
        story.append(_scaled_image(figures["mutation"], 172 * mm, imports))
        story.append(
            _paragraph(
                "Figure 3. Mutation comparison. Coloured lines join the same seed block "
                "across u; the black line shows cell medians.",
                styles["caption"],
                Paragraph,
            )
        )

    feedback = _feedback_group(design)
    feedback_groups = _feedback_groups_by_host(design)
    if feedback and "feedback" in figures:
        if feedback_groups:
            feedback_heading = "Feedback strength at fixed host abundance"
            feedback_text = (
                "The same four alpha levels are evaluated separately at H=100 and "
                "H=1,000. Within each series H and u are fixed; f, total return R and "
                "alpha increase together. Extreme escape fractions in c0018 and c0019 "
                "are sensitivity levels outside the primary range."
            )
            feedback_caption = (
                "Figure 4. Feedback-strength series paired by seed block. Solid circles "
                "show H=100; dashed squares show H=1,000."
            )
        else:
            feedback_heading = "Broad feedback bracket"
            feedback_text = (
                "This boundary check is not an isolated alpha comparison: H and total "
                "return also change. It should be read as an exploratory bracket, not as "
                "a dose-response curve for alpha."
            )
            feedback_caption = (
                "Figure 4. Feedback bracket paired by seed block. Alpha, H and R are "
                "confounded across these cells."
            )
        story.append(_paragraph(feedback_heading, styles["h2"], Paragraph))
        story.append(_paragraph(feedback_text, styles["body"], Paragraph))
        story.append(_scaled_image(figures["feedback"], 172 * mm, imports))
        story.append(
            _paragraph(
                feedback_caption,
                styles["caption"],
                Paragraph,
            )
        )

    extension_groups = _extension_return_groups(design)
    if extension_groups and "extension" in figures:
        story.append(PageBreak())
        story.append(
            _paragraph("Extension pooling comparisons", styles["h2"], Paragraph)
        )
        story.append(
            _paragraph(
                "The extension completes the mutation-enabled R=10^9 pooling series and "
                "the mutation-free R=10^8 series. Each coloured line joins the "
                "same seed block across H.",
                styles["body"],
                Paragraph,
            )
        )
        story.append(_scaled_image(figures["extension"], 172 * mm, imports))
        story.append(
            _paragraph(
                "Figure 5. Extension pooling comparisons paired by seed block. The "
                "black lines show cell medians.",
                styles["caption"],
                Paragraph,
            )
        )

    story.append(PageBreak())
    story.append(_paragraph("Endpoint diversity indices", styles["h1"], Paragraph))
    story.append(
        _paragraph(
            "Values are medians across independent populations. D0 is labelled richness; "
            "TV is total-variation composition distance from the frozen starting state.",
            styles["body"],
            Paragraph,
        )
    )
    if migration_fraction > 0:
        story.append(
            _paragraph(
                f"The migration audit confirmed that every passage exchanged exactly "
                f"{100 * migration_fraction:.0f}% of the one-billion-cell focal "
                "environment, with equal emigrant and immigrant totals.",
                styles["body"],
                Paragraph,
            )
        )
    story.append(
        _pdf_table(
            _metric_table_rows(endpoints, design),
            [31, 35, 61, 68, 45, 45, 55, 42],
            imports,
            font_size=6.7,
        )
    )

    story.append(
        _paragraph(
            "Quality control and computational feasibility", styles["h1"], Paragraph
        )
    )
    resources = summary["overall_resources"]
    story.append(
        _paragraph(
            f"All {summary.get('populations')} populations passed the endpoint checksum "
            f"and environmental-capacity audit. Summed runtime was "
            f"{resources['elapsed_hours']:.2f} hours, diagnostic output was "
            f"{resources['output_mib']:.1f} MiB, and the largest measured process-tree "
            f"memory was {resources['peak_rss_mib']:.1f} MiB. Every run used two worker "
            "processes, five host passages and 500 steady within-host bacterial generations.",
            styles["body"],
            Paragraph,
        )
    )
    benchmark_note = summary.get("benchmark_hardware_note")
    if benchmark_note:
        story.append(
            _paragraph(
                f"Benchmark machine: {benchmark_note}. Timings should be recalibrated with "
                "a small batch on the target HPC before the larger experiment is launched.",
                styles["body"],
                Paragraph,
            )
        )
    material = summary.get("mutation_materialization", {})
    if material:
        story.append(
            _paragraph(
                "Mutation materialization remained far below its safety limit: the largest "
                "within-host transition created "
                f"{material.get('largest_realized_count')} mutant lineages, compared with "
                f"the configured limit of {material.get('configured_limit_per_transition'):,}.",
                styles["body"],
                Paragraph,
            )
        )
    story.append(
        KeepTogether(
            [
                _scaled_image(figures["feasibility"], 145 * mm, imports),
                _paragraph(
                    "Figure 6. Computational cost paired by seed block. The upper row "
                    "varies H at fixed total return and u=0; the lower row varies u at "
                    "H=100. Black lines show cell medians.",
                    styles["caption"],
                    Paragraph,
                ),
            ]
        )
    )
    story.append(PageBreak())
    story.append(
        _paragraph("Observed scaling with host number", styles["h2"], Paragraph)
    )
    story.append(
        _pdf_table(
            _resource_scaling_rows(endpoints, design),
            [95, 115, 115, 115],
            imports,
            font_size=7.2,
        )
    )
    scaling = _resource_scaling_summary(endpoints, design)
    story.append(Spacer(1, 6))
    story.append(
        _paragraph(
            f"Across the mutation-free matched-return cells, the fitted host-number exponent "
            f"was {scaling['elapsed_minutes_slope']:.2f} for runtime and "
            f"{scaling['output_mib_slope']:.2f} for diagnostic output. A tenfold increase in H therefore "
            f"multiplied median runtime by about {scaling['elapsed_minutes_tenfold']:.1f} "
            f"and output by {scaling['output_mib_tenfold']:.1f}. Peak RAM was nearly flat "
            "with H because host tasks are streamed through a bounded queue.",
            styles["body"],
            Paragraph,
        )
    )
    story.append(
        _paragraph(
            f"At H=100, increasing u from zero to the highest pilot value multiplied median "
            f"runtime by {scaling['elapsed_minutes_mutation_ratio']:.1f}, output by "
            f"{scaling['output_mib_mutation_ratio']:.1f}, and peak RAM by "
            f"{scaling['peak_rss_mib_mutation_ratio']:.2f}. The extra cost comes from "
            "materialized mutant lineages and their records.",
            styles["body"],
            Paragraph,
        )
    )
    story.append(_paragraph("How computational work scales", styles["h2"], Paragraph))
    story.append(
        _paragraph(
            "In big-O notation, a useful approximation is O(P H [G S + M]): P is the number of host passages; "
            "H is host number; G is the number of within-host transitions; S is extant "
            "strain richness; and M is the number of materialized mutation events. The "
            "prototype operates on strain counts, not individual bacterial cells. Thus K, "
            "e and R do not create one operation per cell. Larger e can still alter realised "
            "richness and output. Its independent effect cannot be estimated here because e "
            "decreases as H increases in the matched-return series.",
            styles["body"],
            Paragraph,
        )
    )
    story.append(_paragraph("HPC starting allocation", styles["h2"], Paragraph))
    for recommendation in (
        "Request two CPUs and 1 GiB RAM per concurrently running population.",
        "Begin with 8-16 concurrent populations and inspect CPU utilisation and file-system throughput before increasing concurrency.",
        "With 100 CPUs and two workers per population, 50 populations is the CPU ceiling; RAM is unlikely to be limiting at this pilot scale.",
        "Rebenchmark if workers, u, retained host detail, host lifespan or passage number increase. At unchanged output settings, runtime and host-level output should grow roughly linearly with passage number.",
    ):
        story.append(Paragraph("- " + html.escape(recommendation), styles["bullet"]))
    story.append(_paragraph("Interpretation and next step", styles["h1"], Paragraph))
    story.append(
        _paragraph(
            "The pilot supports a longer equilibrium and precision pilot before any "
            "confirmatory experiment. Three populations per cell are suitable for "
            "calibration and graphical ranges, but not for stable confidence intervals "
            "or equivalence tests. Whole populations, not individual hosts or passages, "
            "must be used as independent replicates.",
            styles["body"],
            Paragraph,
        )
    )
    for limitation in summary.get("limitations", []):
        story.append(Paragraph("- " + html.escape(limitation), styles["bullet"]))

    story.append(_paragraph("Reproducibility", styles["h1"], Paragraph))
    source_hashes = summary.get("source_sha256", [])
    source_hash_display = ", ".join(f"{value[:12]}..." for value in source_hashes)
    if len(source_hashes) > 1:
        source_hash_display += (
            " (staged report/launcher/design additions; biological transition "
            "modules unchanged)"
        )
    provenance_rows = [
        ["Item", "Recorded value"],
        ["Experiment", str(summary.get("experiment_id", "pilot"))],
        ["Analysis audit", str(summary.get("audit_status"))],
        ["Source commit", ", ".join(summary.get("software_git_commits", []))],
        ["Source-tree SHA-256", source_hash_display],
        ["Seed blocks", ", ".join(summary.get("seed_blocks", []))],
        ["Analysis schema", str(summary.get("analysis_schema_version", ""))],
    ]
    story.append(_pdf_table(provenance_rows, [100, 340], imports, font_size=7.5))

    story.append(PageBreak())
    story.append(
        _paragraph("Appendix: glossary and diversity measures", styles["h1"], Paragraph)
    )
    glossary = (
        (
            "Relative abundance, pᵢ",
            "The fraction of environmental cells assigned to strain i.",
        ),
        (
            "Labelled richness, D0",
            "The number of distinct strain labels. Every label counts once, even if it "
            "contains only one cell.",
        ),
        (
            "Shannon entropy, H'",
            "H' = -Σᵢ pᵢ ln(pᵢ). Higher values mean that the identity of a sampled "
            "cell is less predictable. H' is an index, not a number of strains.",
        ),
        (
            "Simpson diversity",
            "1 - Σᵢ pᵢ². This is the probability that two independently sampled cells "
            "have different strain labels.",
        ),
        (
            "Hill D1",
            "D₁ = exp(H'). The number of equally abundant strains needed to reproduce the "
            "observed Shannon entropy. It discounts rare strains, but less strongly than D2.",
        ),
        (
            "Hill D2",
            "D₂ = 1 / Σᵢ pᵢ². The number of equally abundant strains needed to reproduce "
            "the observed probability of sampling the same strain twice. It emphasizes "
            "abundant strains.",
        ),
        (
            "Pielou evenness",
            "J = H' / ln(D₀). It compares observed Shannon entropy with the maximum possible "
            "entropy at the observed richness. It ranges from near zero to one.",
        ),
        (
            "Composition distance, TV",
            "TV(A,B) = ½ Σᵢ |pᵢ(A) - pᵢ(B)|. A distance of 0.20 means that 20% of abundance must "
            "be reassigned among strain labels to make the compositions identical.",
        ),
        (
            "Root-collapsed diversity",
            "Diversity after each mutant descendant is merged back into its original "
            "ancestral strain. This separates new labels from redistribution among roots.",
        ),
        (
            "Experimental cell",
            "One defined combination of biological parameter values.",
        ),
        (
            "Seed block",
            "One independently seeded population replicate reused as a reproducible block "
            "across experimental cells.",
        ),
        (
            "H, f, e, R, alpha, u and m",
            "H is host number; f is escape fraction per host; e is escaping cells per host "
            "(fK); R is total returned cells (He); alpha is the host-derived fraction before "
            "regulation; u is whole-genome mutation probability per bacterial generation; "
            "m is the fraction of the focal environment replaced from the fixed regional "
            "pool after each host return.",
        ),
    )
    story.append(
        _pdf_table(
            [["Term", "Meaning"], *[list(row) for row in glossary]],
            [125, 315],
            imports,
            font_size=7.8,
            font_name="DejaVuSans",
        )
    )

    def page_footer(canvas: Any, doc: Any) -> None:
        canvas.saveState()
        canvas.setFont("Helvetica", 7)
        canvas.setFillColor(colors.HexColor("#5F6B76"))
        canvas.drawString(18 * mm, 10 * mm, "Trophosome pilot report")
        canvas.drawRightString(A4[0] - 18 * mm, 10 * mm, f"Page {doc.page}")
        canvas.restoreState()

    document.build(story, onFirstPage=page_footer, onLaterPages=page_footer)


def generate_pilot_report(
    analysis_dir: Path,
    design_path: Path,
    output_pdf: Path,
    *,
    markdown_path: Path | None = None,
    assets_dir: Path | None = None,
    title: str = "Pilot report: host feedback and environmental symbiont diversity",
    report_date: str | None = None,
) -> ReportArtifacts:
    """Generate a PDF and transparent Markdown companion from pilot summary tables."""

    analysis_dir = analysis_dir.resolve()
    design_path = design_path.resolve()
    output_pdf = output_pdf.resolve()
    markdown = (
        markdown_path.resolve()
        if markdown_path is not None
        else output_pdf.with_suffix(".md")
    )
    assets = (
        assets_dir.resolve()
        if assets_dir is not None
        else markdown.parent / f"{markdown.stem}-assets"
    )
    report_date = report_date or date.today().isoformat()
    summary, endpoints, trajectories, design = _load_report_data(
        analysis_dir, design_path
    )
    figures = _make_figures(endpoints, trajectories, design, assets)
    _write_markdown(markdown, title, report_date, summary, endpoints, design, figures)
    _write_pdf(output_pdf, title, report_date, summary, endpoints, design, figures)
    return ReportArtifacts(pdf=output_pdf, markdown=markdown, assets=assets)
