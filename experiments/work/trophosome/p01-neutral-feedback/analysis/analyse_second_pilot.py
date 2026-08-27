#!/usr/bin/env python3
"""Audit and analyse the Phase 1 equilibrium-and-precision second pilot."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import statistics
from collections import defaultdict
from datetime import UTC, datetime
from pathlib import Path
from statistics import NormalDist
from typing import Any

import numpy as np
from analyse_first_pilot import (
    INITIAL_RICHNESS,
    composition_metrics,
    diversity_metrics,
    read_lineages,
    read_tsv,
    sha256,
    trace_roots,
    write_tsv,
)

from trophosome.config import load_config

VARIANT_TAG = "v210-m010-g250"
EXPERIMENT_ID = f"phase1-second-pilot-{VARIANT_TAG}"
DESIGN_STEM = f"phase1-second-pilot-{VARIANT_TAG}"
STAGE_DIRECTORY = f"s02-equilibrium-precision-{VARIANT_TAG}"
MODEL_SPEC_VERSION = "2.1.0"
OUTPUT_SCHEMA_VERSION = "2.3.0"
SOFTWARE_VERSION = "0.7.0"
HOST_GENERATIONS = 250
EXPECTED_SEED_BLOCKS = tuple(f"sb{number:04d}" for number in range(1, 21))
EXPECTED_RUNS = 6 * len(EXPECTED_SEED_BLOCKS)
MIGRATION_REPLACEMENT_COUNT = 100_000_000
RESPONSES = ("D0", "D1", "D2", "evenness", "TV")
SEPARATED_STABILITY_DIAGNOSTIC_VERSION = "1.0.0"
STABILITY_BOOTSTRAP_RESAMPLES = 5_000
STABILITY_BOOTSTRAP_SEED = 20_260_827
REQUIRED_RAW_FILES = (
    "completion.json",
    "environment_counts.csv",
    "execution-summary.json",
    "final_environment_rep000.npz",
    "host_adult_counts.csv",
    "host_adult_summaries.csv",
    "host_generation_summary.csv",
    "infection_counts.csv",
    "migration_counts.csv",
    "pooled_host_counts_and_occupancy.csv",
    "provenance.json",
    "release_counts.csv",
    "resolved_config.json",
    "strain_lineage_events.csv",
    "strain_origins.csv",
)
CONTRASTS = (
    ("host-passage", "p01-s02-c0022", "p01-s02-c0021"),
    ("mutation", "p01-s02-c0023", "p01-s02-c0022"),
    ("many-vs-few-hosts", "p01-s02-c0024", "p01-s02-c0022"),
    ("weak-vs-baseline-feedback", "p01-s02-c0025", "p01-s02-c0022"),
    ("strong-vs-baseline-feedback", "p01-s02-c0026", "p01-s02-c0022"),
)


def _atomic_json(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    temporary.replace(path)


def _resolved_config_sha256(payload: dict[str, Any]) -> str:
    serialized = json.dumps(payload, sort_keys=True, separators=(",", ":")).encode(
        "utf-8"
    )
    return hashlib.sha256(serialized).hexdigest()


def _audit_payload(status: str, issues: list[str]) -> dict[str, Any]:
    return {
        "analysis_audit_schema_version": "1.0.0",
        "experiment_id": EXPERIMENT_ID,
        "variant_id": VARIANT_TAG,
        "checked_at": datetime.now(UTC).isoformat(),
        "status": status,
        "expected_populations": EXPECTED_RUNS,
        "issues": issues,
    }


def completion_gate_issues(
    run_rows: list[dict[str, str]], *, work: Path, scratch: Path
) -> list[str]:
    issues: list[str] = []
    if len(run_rows) != EXPECTED_RUNS:
        issues.append(
            f"manifest contains {len(run_rows)} runs; expected {EXPECTED_RUNS}"
        )
    if len({row.get("run_id") for row in run_rows}) != len(run_rows):
        issues.append("manifest run IDs are not unique")
    if len({row.get("scratch_relative_path") for row in run_rows}) != len(run_rows):
        issues.append("manifest scratch paths are not unique")

    for run in run_rows:
        run_id = run.get("run_id", "unknown-run")
        config = work / run.get("config_path", "")
        frozen_resolved: dict[str, Any] | None = None
        if not config.is_file():
            issues.append(f"{run_id}: configuration is missing")
        elif sha256(config) != run.get("config_sha256"):
            issues.append(f"{run_id}: configuration checksum differs from manifest")
        else:
            frozen_resolved = json.loads(json.dumps(load_config(config).to_dict()))

        output = scratch / run.get("scratch_relative_path", "")
        missing = [name for name in REQUIRED_RAW_FILES if not (output / name).is_file()]
        if missing:
            issues.append(f"{run_id}: missing {', '.join(missing)}")
            continue
        try:
            completion = json.loads(
                (output / "completion.json").read_text(encoding="utf-8")
            )
        except (OSError, ValueError) as exc:
            issues.append(f"{run_id}: unreadable completion.json ({exc})")
            continue
        if completion.get("complete") is not True:
            issues.append(f"{run_id}: completion marker is not committed")
        for field, expected in (
            ("model_spec_version", MODEL_SPEC_VERSION),
            ("output_schema_version", OUTPUT_SCHEMA_VERSION),
            ("software_version", SOFTWARE_VERSION),
        ):
            if completion.get(field) != expected:
                issues.append(
                    f"{run_id}: {field}={completion.get(field)!r}; "
                    f"expected {expected!r}"
                )
        resolved_path = output / "resolved_config.json"
        try:
            observed_resolved = json.loads(resolved_path.read_text(encoding="utf-8"))
        except (OSError, ValueError) as exc:
            issues.append(f"{run_id}: unreadable resolved_config.json ({exc})")
            continue
        if frozen_resolved is None:
            continue
        if observed_resolved != frozen_resolved:
            issues.append(f"{run_id}: resolved configuration differs from frozen TOML")
        if completion.get("config_sha256") != _resolved_config_sha256(frozen_resolved):
            issues.append(f"{run_id}: completed resolved-configuration hash differs")
        final_path = output / "final_environment_rep000.npz"
        recorded_final = completion.get("final_environment_sha256", {}).get(
            final_path.name
        )
        if not recorded_final or sha256(final_path) != recorded_final:
            issues.append(f"{run_id}: final-environment checksum differs")
        for name, expected_size in completion.get("output_sizes", {}).items():
            path = output / name
            if not path.is_file() or path.stat().st_size != int(expected_size):
                issues.append(f"{run_id}: committed size differs for {name}")
    return issues


def read_environment_trajectory(path: Path) -> dict[int, dict[int, int]]:
    trajectory: dict[int, dict[int, int]] = defaultdict(dict)
    with path.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            trajectory[int(row["generation"])][int(row["strain_id"])] = int(
                row["count"]
            )
    return dict(trajectory)


def _root_collapsed_counts(
    counts: dict[int, int], roots: dict[int, int]
) -> dict[int, int]:
    collapsed: dict[int, int] = defaultdict(int)
    for strain_id, count in counts.items():
        collapsed[roots[strain_id]] += count
    return dict(collapsed)


def _lag_one(values: np.ndarray) -> float:
    if len(values) < 3 or np.all(values == values[0]):
        return 0.0
    return float(np.corrcoef(values[:-1], values[1:])[0, 1])


def integrated_autocorrelation_time(values: np.ndarray) -> float:
    values = np.asarray(values, dtype=float)
    if len(values) < 3 or np.all(values == values[0]):
        return 1.0
    centered = values - values.mean()
    denominator = float(np.dot(centered, centered))
    if denominator <= 0:
        return 1.0
    correlations: list[float] = []
    for lag in range(1, len(values)):
        correlation = float(np.dot(centered[:-lag], centered[lag:]) / denominator)
        correlations.append(correlation)
    positive_sum = 0.0
    for start in range(0, len(correlations), 2):
        pair = correlations[start : start + 2]
        if sum(pair) <= 0:
            break
        positive_sum += sum(pair)
    return max(1.0, 1.0 + 2.0 * positive_sum)


def _rank_normalize(values: np.ndarray) -> np.ndarray:
    flat = values.ravel()
    order = np.argsort(flat, kind="mergesort")
    ranks = np.empty(len(flat), dtype=float)
    start = 0
    while start < len(flat):
        stop = start + 1
        while stop < len(flat) and flat[order[stop]] == flat[order[start]]:
            stop += 1
        ranks[order[start:stop]] = 0.5 * (start + stop - 1) + 1
        start = stop
    probabilities = (ranks - 0.375) / (len(flat) + 0.25)
    normal = NormalDist()
    transformed = np.asarray([normal.inv_cdf(value) for value in probabilities])
    return transformed.reshape(values.shape)


def rank_normalized_split_rhat(chains: np.ndarray) -> float:
    chains = np.asarray(chains, dtype=float)
    if chains.ndim != 2 or chains.shape[0] < 2 or chains.shape[1] < 4:
        return math.nan
    half = chains.shape[1] // 2
    split = np.concatenate((chains[:, :half], chains[:, -half:]), axis=0)
    ranked = _rank_normalize(split)
    within = float(np.mean(np.var(ranked, axis=1, ddof=1)))
    between = float(half * np.var(np.mean(ranked, axis=1), ddof=1))
    if within == 0:
        return 1.0 if between == 0 else math.inf
    variance = ((half - 1) / half) * within + between / half
    return float(math.sqrt(variance / within))


def approximate_combined_ess(chains: np.ndarray) -> float:
    chains = np.asarray(chains, dtype=float)
    if chains.ndim != 2 or chains.size == 0:
        return math.nan
    taus = [integrated_autocorrelation_time(chain) for chain in chains]
    return float(sum(len(chain) / tau for chain, tau in zip(chains, taus, strict=True)))


def _mean_ci90(values: list[float]) -> tuple[float, float, float]:
    mean = statistics.fmean(values)
    if len(values) < 2:
        return mean, math.nan, math.nan
    standard_error = statistics.stdev(values) / math.sqrt(len(values))
    degrees_freedom = len(values) - 1
    z = NormalDist().inv_cdf(0.95)
    critical = (
        z
        + (z**3 + z) / (4 * degrees_freedom)
        + (5 * z**5 + 16 * z**3 + 3 * z) / (96 * degrees_freedom**2)
        + (3 * z**7 + 19 * z**5 + 17 * z**3 - 15 * z)
        / (384 * degrees_freedom**3)
    )
    return mean, mean - critical * standard_error, mean + critical * standard_error


def _margin(response: str, reference: float) -> float:
    if response in {"D0", "D1", "D2"}:
        return 0.05 * abs(reference)
    if response == "evenness":
        return 0.02
    if response == "TV":
        return 0.05
    raise ValueError(f"unknown response: {response}")


def _inside_margin(lower: float, upper: float, margin: float) -> bool:
    return math.isfinite(lower) and lower >= -margin and upper <= margin


def _assessment_pass(
    window_values: dict[str, list[np.ndarray]],
    response: str,
    window_ids: tuple[int, int, int],
) -> tuple[bool, float, float]:
    checks: list[bool] = []
    largest_ci_ratio = 0.0
    for window_id in window_ids:
        slopes: list[float] = []
        references: list[float] = []
        for values in window_values[response]:
            segment = values[window_id]
            x = np.arange(len(segment), dtype=float)
            slope = float(np.polyfit(x, segment, 1)[0] * len(segment))
            slopes.append(slope)
            references.append(float(np.mean(segment)))
        _, lower, upper = _mean_ci90(slopes)
        margin = _margin(response, statistics.fmean(references))
        checks.append(_inside_margin(lower, upper, margin))
        largest_ci_ratio = max(largest_ci_ratio, max(abs(lower), abs(upper)) / margin)

    for first, second in zip(window_ids[:-1], window_ids[1:], strict=True):
        differences = [
            float(np.mean(values[second]) - np.mean(values[first]))
            for values in window_values[response]
        ]
        references = [
            float(np.mean(values[first])) for values in window_values[response]
        ]
        _, lower, upper = _mean_ci90(differences)
        margin = _margin(response, statistics.fmean(references))
        checks.append(_inside_margin(lower, upper, margin))
        largest_ci_ratio = max(largest_ci_ratio, max(abs(lower), abs(upper)) / margin)
    return (
        all(checks),
        largest_ci_ratio,
        _margin(
            response,
            statistics.fmean(
                float(np.mean(values[window_ids[-1]]))
                for values in window_values[response]
            ),
        ),
    )


def _ci_state(lower: float, upper: float, margin: float) -> str:
    """Classify an interval relative to a symmetric biological margin."""

    if _inside_margin(lower, upper, margin):
        return "equivalent"
    if math.isfinite(lower) and lower > margin:
        return "increasing"
    if math.isfinite(upper) and upper < -margin:
        return "decreasing"
    return "inconclusive"


def _bootstrap_rng(*parts: object) -> np.random.Generator:
    token = ":".join(str(part) for part in (STABILITY_BOOTSTRAP_SEED, *parts))
    digest = hashlib.sha256(token.encode("utf-8")).digest()
    return np.random.default_rng(int.from_bytes(digest[:8], "big"))


def _bootstrap_sd_interval(
    values: np.ndarray,
    *,
    key: tuple[object, ...],
    resamples: int,
) -> tuple[float, float, float]:
    values = np.asarray(values, dtype=float)
    observed = float(np.std(values, ddof=1))
    indices = _bootstrap_rng(*key).integers(
        0, len(values), size=(resamples, len(values))
    )
    estimates = np.std(values[indices], axis=1, ddof=1)
    return (
        observed,
        float(np.quantile(estimates, 0.05)),
        float(np.quantile(estimates, 0.95)),
    )


def _bootstrap_sd_change_interval(
    first: np.ndarray,
    second: np.ndarray,
    *,
    key: tuple[object, ...],
    resamples: int,
) -> tuple[float, float, float]:
    first = np.asarray(first, dtype=float)
    second = np.asarray(second, dtype=float)
    if len(first) != len(second):
        raise ValueError("paired spread comparison has unequal seed counts")
    observed = float(np.std(second, ddof=1) - np.std(first, ddof=1))
    indices = _bootstrap_rng(*key).integers(
        0, len(first), size=(resamples, len(first))
    )
    estimates = np.std(second[indices], axis=1, ddof=1) - np.std(
        first[indices], axis=1, ddof=1
    )
    return (
        observed,
        float(np.quantile(estimates, 0.05)),
        float(np.quantile(estimates, 0.95)),
    )


def _stability_status(states: list[str], recent_states: list[str]) -> str:
    if states and all(state == "equivalent" for state in states):
        return "stable"
    if any(state in {"increasing", "decreasing"} for state in recent_states):
        return "continuing_change"
    return "inconclusive"


def _direction(states: list[str]) -> str:
    directions = {
        state for state in states if state in {"increasing", "decreasing"}
    }
    if not directions:
        return "none"
    if len(directions) == 1:
        return next(iter(directions))
    return "mixed"


def separated_stability_diagnostic(
    cell_id: str,
    cell_trajectories: list[list[dict[str, Any]]],
    *,
    window_length: int,
    end_generation: int = HOST_GENERATIONS,
    bootstrap_resamples: int = STABILITY_BOOTSTRAP_RESAMPLES,
) -> list[dict[str, Any]]:
    """Separate temporal settling, distribution stability, and seed heterogeneity.

    This is a post-hoc exploratory diagnostic. Slopes are first calculated inside
    each independently seeded population and then summarized across seed blocks.
    Changes in the location and spread of the seed-block distribution are assessed
    separately. R-hat and ESS are retained as secondary diagnostics and do not
    determine the biological classification.
    """

    if len(cell_trajectories) < 12:
        raise ValueError(f"{cell_id} needs at least 12 trajectories")
    if bootstrap_resamples < 100:
        raise ValueError("at least 100 bootstrap resamples are required")
    start = end_generation - 4 * window_length + 1
    if start < 1:
        raise ValueError(f"four diagnostic windows do not fit for {cell_id}")

    window_values: dict[str, np.ndarray] = {}
    for response in RESPONSES:
        seed_windows: list[list[np.ndarray]] = []
        for trajectory in cell_trajectories:
            by_generation = {int(row["generation"]): row for row in trajectory}
            windows: list[np.ndarray] = []
            for window_id in range(4):
                lower = start + window_id * window_length
                upper = lower + window_length
                windows.append(
                    np.asarray(
                        [
                            float(by_generation[generation][response])
                            for generation in range(lower, upper)
                        ],
                        dtype=float,
                    )
                )
            seed_windows.append(windows)
        window_values[response] = np.asarray(seed_windows)

    rows: list[dict[str, Any]] = []
    for response, values in window_values.items():
        window_means = np.mean(values, axis=2)
        trend_states: list[str] = []
        trend_ratios: list[float] = []
        for window_id in range(4):
            slopes = [
                float(
                    np.polyfit(np.arange(window_length), segment, 1)[0]
                    * window_length
                )
                for segment in values[:, window_id, :]
            ]
            references = [
                float(np.mean(segment)) for segment in values[:, window_id, :]
            ]
            _mean, lower, upper = _mean_ci90(slopes)
            margin = _margin(response, statistics.fmean(references))
            trend_states.append(_ci_state(lower, upper, margin))
            trend_ratios.append(max(abs(lower), abs(upper)) / margin)
        within_seed_trend_status = _stability_status(
            trend_states, trend_states[1:]
        )

        location_states: list[str] = []
        spread_states: list[str] = []
        spread_ratios: list[float] = []
        for first_window, second_window in ((0, 1), (1, 2), (2, 3)):
            first = window_means[:, first_window]
            second = window_means[:, second_window]
            differences = [float(b - a) for a, b in zip(first, second, strict=True)]
            _mean, lower, upper = _mean_ci90(differences)
            margin = _margin(response, float(np.mean(first)))
            location_states.append(_ci_state(lower, upper, margin))

            _spread_change, spread_lower, spread_upper = (
                _bootstrap_sd_change_interval(
                    first,
                    second,
                    key=(cell_id, response, first_window, second_window, "spread"),
                    resamples=bootstrap_resamples,
                )
            )
            spread_states.append(_ci_state(spread_lower, spread_upper, margin))
            spread_ratios.append(
                max(abs(spread_lower), abs(spread_upper)) / margin
            )

        location_status = _stability_status(
            location_states, location_states[1:]
        )
        spread_status = _stability_status(spread_states, spread_states[1:])
        distribution_states = location_states + spread_states
        recent_distribution_states = location_states[1:] + spread_states[1:]
        between_seed_distribution_status = _stability_status(
            distribution_states, recent_distribution_states
        )

        final_means = window_means[:, 3]
        final_reference = float(np.mean(final_means))
        heterogeneity_margin = _margin(response, final_reference)
        between_sd, between_sd_lower, between_sd_upper = _bootstrap_sd_interval(
            final_means,
            key=(cell_id, response, "final-between-seed-sd"),
            resamples=bootstrap_resamples,
        )
        if between_sd_lower > heterogeneity_margin:
            heterogeneity_status = "meaningful"
        elif between_sd_upper <= heterogeneity_margin:
            heterogeneity_status = "negligible"
        else:
            heterogeneity_status = "uncertain"

        chains = np.asarray(
            [np.concatenate((seed[1], seed[2], seed[3])) for seed in values]
        )
        rhat = rank_normalized_split_rhat(chains)
        ess = approximate_combined_ess(chains)

        if (
            within_seed_trend_status == "stable"
            and between_seed_distribution_status == "stable"
        ):
            if heterogeneity_status == "meaningful":
                classification = "stable_with_persistent_heterogeneity"
            elif heterogeneity_status == "negligible":
                classification = "stable_with_negligible_heterogeneity"
            else:
                classification = "stable_heterogeneity_uncertain"
        elif (
            within_seed_trend_status == "continuing_change"
            or between_seed_distribution_status == "continuing_change"
        ):
            classification = "continuing_change_detected"
        else:
            classification = "stability_unresolved"

        rows.append(
            {
                "cell_id": cell_id,
                "response": response,
                "assessment_end_generation": end_generation,
                "window_generations": window_length,
                "diagnostic_version": SEPARATED_STABILITY_DIAGNOSTIC_VERSION,
                "diagnostic_scope": "post_hoc_exploratory",
                "within_seed_trend_status": within_seed_trend_status,
                "within_seed_trend_direction": _direction(trend_states[1:]),
                "largest_trend_ci_to_margin_ratio": max(trend_ratios),
                "between_seed_location_status": location_status,
                "between_seed_location_direction": _direction(location_states[1:]),
                "between_seed_spread_status": spread_status,
                "between_seed_spread_direction": _direction(spread_states[1:]),
                "largest_spread_ci_to_margin_ratio": max(spread_ratios),
                "between_seed_distribution_status": (
                    between_seed_distribution_status
                ),
                "between_seed_sd_final_window": between_sd,
                "between_seed_sd_ci90_lower": between_sd_lower,
                "between_seed_sd_ci90_upper": between_sd_upper,
                "heterogeneity_margin": heterogeneity_margin,
                "between_seed_heterogeneity_status": heterogeneity_status,
                "rank_normalized_split_rhat_secondary": rhat,
                "approximate_combined_ess_secondary": ess,
                "final_classification": classification,
            }
        )
    return rows


def _matched_seed_blocks(
    trajectories_by_cell_seed: dict[tuple[str, str], list[dict[str, Any]]],
    cell_ids: tuple[str, ...],
) -> tuple[str, ...]:
    seed_sets = {
        cell_id: {
            seed_block
            for observed_cell, seed_block in trajectories_by_cell_seed
            if observed_cell == cell_id
        }
        for cell_id in cell_ids
    }
    if any(len(seed_blocks) < 12 for seed_blocks in seed_sets.values()):
        counts = ", ".join(
            f"{cell_id}={len(seed_blocks)}"
            for cell_id, seed_blocks in seed_sets.items()
        )
        raise ValueError(f"at least 12 seed blocks per cell are required ({counts})")
    unique_sets = {frozenset(seed_blocks) for seed_blocks in seed_sets.values()}
    if len(unique_sets) != 1:
        raise ValueError("cells do not contain the same matched seed blocks")
    return tuple(sorted(next(iter(unique_sets))))


def separated_stability_rows(
    trajectories_by_cell_seed: dict[tuple[str, str], list[dict[str, Any]]],
    design: dict[str, dict[str, str]],
) -> list[dict[str, Any]]:
    seed_blocks = _matched_seed_blocks(trajectories_by_cell_seed, tuple(design))
    rows: list[dict[str, Any]] = []
    for cell_id in design:
        rows.extend(
            separated_stability_diagnostic(
                cell_id,
                [
                    trajectories_by_cell_seed[(cell_id, seed_block)]
                    for seed_block in seed_blocks
                ],
                window_length=int(design[cell_id]["diagnostic_window_generations"]),
            )
        )
    return rows


def _separated_stability_summary(rows: list[dict[str, Any]]) -> dict[str, Any]:
    counts: dict[str, int] = defaultdict(int)
    for row in rows:
        counts[str(row["final_classification"])] += 1
    return {
        "diagnostic_version": SEPARATED_STABILITY_DIAGNOSTIC_VERSION,
        "scope": "post-hoc exploratory",
        "bootstrap_resamples": STABILITY_BOOTSTRAP_RESAMPLES,
        "classification_counts": dict(sorted(counts.items())),
        "rhat_and_ess_role": "secondary diagnostics; not classification criteria",
    }


def equilibrium_screen(
    cell_id: str,
    cell_trajectories: list[list[dict[str, Any]]],
    *,
    window_length: int,
    end_generation: int = HOST_GENERATIONS,
) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    if len(cell_trajectories) < 12:
        raise ValueError(f"{cell_id} needs at least 12 trajectories")
    start = end_generation - 4 * window_length + 1
    if start < 1:
        raise ValueError(f"four diagnostic windows do not fit for {cell_id}")

    window_values: dict[str, list[np.ndarray]] = {name: [] for name in RESPONSES}
    run_rows: list[dict[str, Any]] = []
    for trajectory in cell_trajectories:
        by_generation = {int(row["generation"]): row for row in trajectory}
        seed_block = str(trajectory[0]["seed_block_id"])
        for response in RESPONSES:
            segments: list[np.ndarray] = []
            for window_id in range(4):
                lower = start + window_id * window_length
                upper = lower + window_length
                segment = np.asarray(
                    [
                        float(by_generation[generation][response])
                        for generation in range(lower, upper)
                    ],
                    dtype=float,
                )
                segments.append(segment)
                mean = float(np.mean(segment))
                sd = float(np.std(segment, ddof=1)) if len(segment) > 1 else 0.0
                run_rows.append(
                    {
                        "cell_id": cell_id,
                        "seed_block_id": seed_block,
                        "response": response,
                        "window": window_id + 1,
                        "start_generation": lower,
                        "end_generation": upper - 1,
                        "mean": mean,
                        "median": float(np.median(segment)),
                        "sd": sd,
                        "cv": sd / abs(mean) if mean else math.nan,
                        "q05": float(np.quantile(segment, 0.05)),
                        "q95": float(np.quantile(segment, 0.95)),
                        "lag1_autocorrelation": _lag_one(segment),
                        "integrated_autocorrelation_time": (
                            integrated_autocorrelation_time(segment)
                        ),
                    }
                )
            window_values[response].append(np.asarray(segments))

    screen_rows: list[dict[str, Any]] = []
    for response in RESPONSES:
        previous_pass, previous_ratio, margin = _assessment_pass(
            window_values, response, (0, 1, 2)
        )
        current_pass, current_ratio, _ = _assessment_pass(
            window_values, response, (1, 2, 3)
        )
        chains = np.asarray(
            [
                np.concatenate((values[1], values[2], values[3]))
                for values in window_values[response]
            ]
        )
        rhat = rank_normalized_split_rhat(chains)
        ess = approximate_combined_ess(chains)
        diagnostic_pass = (
            previous_pass
            and current_pass
            and math.isfinite(rhat)
            and rhat < 1.05
            and ess >= 400
        )
        screen_rows.append(
            {
                "cell_id": cell_id,
                "response": response,
                "assessment_end_generation": end_generation,
                "window_generations": window_length,
                "previous_assessment_equivalent": previous_pass,
                "current_assessment_equivalent": current_pass,
                "largest_previous_ci_to_margin_ratio": previous_ratio,
                "largest_current_ci_to_margin_ratio": current_ratio,
                "equivalence_margin": margin,
                "rank_normalized_split_rhat": rhat,
                "approximate_combined_ess": ess,
                "stationarity_screen_pass": diagnostic_pass,
                "full_equilibrium_status": (
                    "not established: contrasting initial conditions not tested"
                ),
            }
        )
    return run_rows, screen_rows


def _assessment_endpoints(window_length: int) -> list[int]:
    endpoints = list(range(4 * window_length, HOST_GENERATIONS + 1, window_length))
    if not endpoints or endpoints[-1] != HOST_GENERATIONS:
        endpoints.append(HOST_GENERATIONS)
    return endpoints


def _persistent_passing_generation(rows: list[dict[str, Any]]) -> int | None:
    ordered = sorted(rows, key=lambda row: int(row["assessment_end_generation"]))
    for index, row in enumerate(ordered):
        if all(bool(later["stationarity_screen_pass"]) for later in ordered[index:]):
            return int(row["assessment_end_generation"])
    return None


def _recommended_replicates(raw: int) -> tuple[int, bool]:
    if raw <= 20:
        return 20, False
    batched = 20 + math.ceil((raw - 20) / 8) * 8
    return min(100, batched), batched > 100


def precision_rows(
    trajectories_by_cell_seed: dict[tuple[str, str], list[dict[str, Any]]],
    design: dict[str, dict[str, str]],
) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for contrast, treatment, reference in CONTRASTS:
        seed_blocks = _matched_seed_blocks(
            trajectories_by_cell_seed, (treatment, reference)
        )
        window_length = int(design[treatment]["diagnostic_window_generations"])
        reference_window = int(design[reference]["diagnostic_window_generations"])
        for response in ("D1", "D2", "TV"):
            differences: list[float] = []
            for seed in seed_blocks:
                treatment_rows = trajectories_by_cell_seed[(treatment, seed)]
                reference_rows = trajectories_by_cell_seed[(reference, seed)]
                treatment_value = statistics.fmean(
                    float(row[response])
                    for row in treatment_rows
                    if int(row["generation"]) > HOST_GENERATIONS - 3 * window_length
                )
                reference_value = statistics.fmean(
                    float(row[response])
                    for row in reference_rows
                    if int(row["generation"]) > HOST_GENERATIONS - 3 * reference_window
                )
                if response in {"D1", "D2"}:
                    differences.append(
                        (treatment_value - reference_value) / reference_value
                    )
                else:
                    differences.append(treatment_value - reference_value)
            standard_deviation = statistics.stdev(differences)
            desired_half_width = 0.05
            raw = max(
                1, math.ceil((1.96 * standard_deviation / desired_half_width) ** 2)
            )
            recommended, exceeds_cap = _recommended_replicates(raw)
            mean, lower, upper = _mean_ci90(differences)
            rows.append(
                {
                    "contrast": contrast,
                    "treatment_cell": treatment,
                    "reference_cell": reference,
                    "response": response,
                    "scale": "relative" if response in {"D1", "D2"} else "absolute",
                    "pilot_mean_difference": mean,
                    "pilot_sd_difference": standard_deviation,
                    "pilot_ci90_lower": lower,
                    "pilot_ci90_upper": upper,
                    "desired_ci_half_width": desired_half_width,
                    "formula_required_replicates": raw,
                    "recommended_replicates": recommended,
                    "exceeds_predeclared_cap": exceeds_cap,
                }
            )
    return rows


def _load_design(path: Path) -> dict[str, dict[str, str]]:
    rows = read_tsv(path)
    design = {row["cell_id"]: row for row in rows}
    if len(design) != 6:
        raise RuntimeError("second-pilot design does not contain six sentinel cells")
    return design


def _audit_migration(path: Path, run_id: str, issues: list[str]) -> None:
    totals: dict[int, list[int]] = defaultdict(lambda: [0, 0])
    with path.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            generation = int(row["generation"])
            totals[generation][0] += int(row["emigrant_count"])
            totals[generation][1] += int(row["immigrant_count"])
    if set(totals) != set(range(1, HOST_GENERATIONS + 1)):
        issues.append(f"{run_id}: migration generations are incomplete")
    for generation, (emigrants, immigrants) in totals.items():
        if emigrants != MIGRATION_REPLACEMENT_COUNT:
            issues.append(f"{run_id}: generation {generation} emigrants differ")
        if immigrants != MIGRATION_REPLACEMENT_COUNT:
            issues.append(f"{run_id}: generation {generation} immigrants differ")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--repository", type=Path, default=Path(__file__).resolve().parents[5]
    )
    parser.add_argument(
        "--diagnostics-from-derived",
        action="store_true",
        help=(
            "recompute only the post-hoc separated stability diagnostic from the "
            "audited environment-trajectories.tsv; do not inspect raw simulations"
        ),
    )
    args = parser.parse_args()
    repository = args.repository.resolve()
    work = repository / "experiments/work/trophosome"
    phase = work / "p01-neutral-feedback"
    layout = json.loads((work / "layout.local.json").read_text(encoding="utf-8"))
    scratch = Path(layout["scratch"])
    design_path = phase / "design" / f"{DESIGN_STEM}-cells.tsv"
    manifest_path = phase / "manifests" / f"{DESIGN_STEM}-runs.tsv"
    derived = phase / "analysis" / f"{STAGE_DIRECTORY}-derived"
    audit_path = derived / "analysis-audit.json"
    design = _load_design(design_path)

    if args.diagnostics_from_derived:
        trajectory_path = derived / "environment-trajectories.tsv"
        summary_path = derived / "analysis-summary.json"
        if not trajectory_path.is_file() or not summary_path.is_file():
            raise SystemExit(
                "derived-only diagnostics require environment-trajectories.tsv "
                "and analysis-summary.json"
            )
        trajectory_rows = read_tsv(trajectory_path)
        trajectories_by_cell_seed: dict[
            tuple[str, str], list[dict[str, Any]]
        ] = defaultdict(list)
        for row in trajectory_rows:
            trajectories_by_cell_seed[(row["cell_id"], row["seed_block_id"])].append(
                row
            )
        try:
            derived_seed_blocks = _matched_seed_blocks(
                dict(trajectories_by_cell_seed), tuple(design)
            )
        except ValueError as exc:
            raise SystemExit(f"derived trajectories are incomplete: {exc}") from exc
        expected_keys = {
            (cell_id, seed_block)
            for cell_id in design
            for seed_block in derived_seed_blocks
        }
        if set(trajectories_by_cell_seed) != expected_keys:
            raise SystemExit(
                "derived trajectories do not form a complete matched cell-by-seed "
                "matrix"
            )
        if any(
            len(rows) != HOST_GENERATIONS + 1
            for rows in trajectories_by_cell_seed.values()
        ):
            raise SystemExit("a derived trajectory does not contain 251 states")
        separated_rows = separated_stability_rows(
            dict(trajectories_by_cell_seed), design
        )
        write_tsv(
            derived / "separated-stability-diagnostic.tsv",
            separated_rows,
            list(separated_rows[0]),
        )
        summary = json.loads(summary_path.read_text(encoding="utf-8"))
        if summary.get("audit_status") != "PASS":
            raise SystemExit("the existing analysis audit is not PASS")
        summary["analysis_schema_version"] = "2.1.0"
        summary["separated_stability_diagnostic"] = (
            _separated_stability_summary(separated_rows)
        )
        _atomic_json(summary_path, summary)
        print(
            "Separated stability diagnostic written from audited derived trajectories: "
            f"{len(separated_rows)} cell-response classifications"
        )
        return 0

    run_rows = read_tsv(manifest_path)

    preflight_issues = completion_gate_issues(run_rows, work=work, scratch=scratch)
    if preflight_issues:
        _atomic_json(audit_path, _audit_payload("FAIL", preflight_issues))
        raise SystemExit(
            "second pilot is not ready for reporting:\n" + "\n".join(preflight_issues)
        )

    initial_payload = json.loads(
        (work / "common/initial-populations/ip001-fisher100.json").read_text(
            encoding="utf-8"
        )
    )
    initial_counts = {
        strain_id: int(count)
        for strain_id, count in enumerate(initial_payload["scaled_counts"])
    }
    trajectories: list[dict[str, Any]] = []
    trajectories_by_cell_seed: dict[tuple[str, str], list[dict[str, Any]]] = {}
    audit_issues: list[str] = []
    commits: set[str] = set()
    source_hashes: set[str] = set()
    platforms: set[str] = set()
    total_elapsed_seconds = 0.0
    total_output_bytes = 0
    maximum_rss_kib = 0
    total_generated_lineages = 0
    maximum_transition_mutants = 0

    for run in run_rows:
        run_id = run["run_id"]
        cell_id = run["cell_id"]
        seed_block = run["seed_block_id"]
        output = scratch / run["scratch_relative_path"]
        environment = read_environment_trajectory(output / "environment_counts.csv")
        expected_generations = set(range(HOST_GENERATIONS + 1))
        if set(environment) != expected_generations:
            audit_issues.append(f"{run_id}: environmental trajectory is incomplete")
            continue
        for generation, counts in environment.items():
            if sum(counts.values()) != 1_000_000_000:
                audit_issues.append(
                    f"{run_id}: environment capacity differs at generation {generation}"
                )
        all_environment_ids = set().union(
            *(counts.keys() for counts in environment.values())
        )
        parents, generated_lineages, maximum_mutants = read_lineages(
            output / "strain_lineage_events.csv", all_environment_ids
        )
        roots = trace_roots(parents, all_environment_ids)
        total_generated_lineages += generated_lineages
        maximum_transition_mutants = max(maximum_transition_mutants, maximum_mutants)

        with (output / "host_generation_summary.csv").open(
            newline="", encoding="utf-8"
        ) as handle:
            summaries = {
                int(row["host_generation"]): row for row in csv.DictReader(handle)
            }
        if set(summaries) != set(range(1, HOST_GENERATIONS + 1)):
            audit_issues.append(f"{run_id}: host-generation summary is incomplete")
            continue
        _audit_migration(output / "migration_counts.csv", run_id, audit_issues)

        run_trajectory: list[dict[str, Any]] = []
        previous_counts: dict[int, int] | None = None
        for generation in range(HOST_GENERATIONS + 1):
            counts = environment[generation]
            metrics = diversity_metrics(counts)
            root_counts = _root_collapsed_counts(counts, roots)
            root_metrics = diversity_metrics(root_counts)
            initial_composition = composition_metrics(counts, initial_counts)
            turnover = (
                0.0
                if previous_counts is None
                else composition_metrics(counts, previous_counts)["total_variation"]
            )
            mutant_abundance = sum(
                count
                for strain_id, count in counts.items()
                if strain_id >= INITIAL_RICHNESS
            )
            summary = summaries.get(generation)
            row: dict[str, Any] = {
                "run_id": run_id,
                "cell_id": cell_id,
                "cell": run["cell"],
                "seed_block_id": seed_block,
                "generation": generation,
                "feedback_exposure": (
                    generation * -math.log1p(-float(design[cell_id]["alpha"]))
                    if float(design[cell_id]["alpha"]) > 0
                    else 0.0
                ),
                "D0": metrics["richness"],
                "shannon": metrics["shannon"],
                "D1": metrics["hill_q1"],
                "D2": metrics["hill_q2"],
                "evenness": metrics["evenness"],
                "top_frequency": metrics["top_frequency"],
                "TV": initial_composition["total_variation"],
                "jensen_shannon": initial_composition["jensen_shannon"],
                "turnover": turnover,
                "mutant_richness": sum(
                    strain_id >= INITIAL_RICHNESS for strain_id in counts
                ),
                "mutant_abundance_fraction": mutant_abundance / sum(counts.values()),
                "root_D0": root_metrics["richness"],
                "root_D1": root_metrics["hill_q1"],
                "root_D2": root_metrics["hill_q2"],
                "root_evenness": root_metrics["evenness"],
                "mean_adult_richness": (
                    None if summary is None else float(summary["mean_adult_richness"])
                ),
                "mean_adult_gene_diversity": (
                    None
                    if summary is None
                    else float(summary["mean_adult_gene_diversity"])
                ),
                "realized_host_feedback": (
                    None
                    if summary is None
                    else float(summary["realized_host_feedback"])
                ),
                "expected_host_feedback_after_migration": (
                    None
                    if summary is None
                    else float(summary["expected_host_feedback_after_migration"])
                ),
            }
            trajectories.append(row)
            run_trajectory.append(row)
            previous_counts = counts
        trajectories_by_cell_seed[(cell_id, seed_block)] = run_trajectory

        execution = json.loads(
            (output / "execution-summary.json").read_text(encoding="utf-8")
        )
        provenance = json.loads(
            (output / "provenance.json").read_text(encoding="utf-8")
        )
        total_elapsed_seconds += float(execution.get("elapsed_seconds") or 0)
        total_output_bytes += int(execution.get("output_bytes") or 0)
        maximum_rss_kib = max(
            maximum_rss_kib, int(execution.get("peak_process_tree_rss_kib") or 0)
        )
        commits.add(str(provenance.get("git_commit")))
        source_hashes.add(str(provenance.get("source_sha256")))
        platforms.add(str(provenance.get("platform")))

    if audit_issues:
        _atomic_json(audit_path, _audit_payload("FAIL", audit_issues))
        raise SystemExit("second-pilot audit failed:\n" + "\n".join(audit_issues))

    trajectory_fields = list(trajectories[0])
    write_tsv(derived / "environment-trajectories.tsv", trajectories, trajectory_fields)

    run_window_rows: list[dict[str, Any]] = []
    screen_rows: list[dict[str, Any]] = []
    screen_history_rows: list[dict[str, Any]] = []
    seed_blocks = _matched_seed_blocks(trajectories_by_cell_seed, tuple(design))
    if seed_blocks != EXPECTED_SEED_BLOCKS:
        raise SystemExit(
            "completed trajectories do not contain the 20 frozen seed blocks"
        )
    for cell_id in design:
        cell_trajectories = [
            trajectories_by_cell_seed[(cell_id, seed_block)]
            for seed_block in seed_blocks
        ]
        run_windows, screens = equilibrium_screen(
            cell_id,
            cell_trajectories,
            window_length=int(design[cell_id]["diagnostic_window_generations"]),
        )
        window_length = int(design[cell_id]["diagnostic_window_generations"])
        cell_history: list[dict[str, Any]] = []
        for endpoint in _assessment_endpoints(window_length):
            _unused_windows, endpoint_screens = equilibrium_screen(
                cell_id,
                cell_trajectories,
                window_length=window_length,
                end_generation=endpoint,
            )
            cell_history.extend(endpoint_screens)
        for screen in screens:
            response_history = [
                row for row in cell_history if row["response"] == screen["response"]
            ]
            passing_endpoints = [
                int(row["assessment_end_generation"])
                for row in response_history
                if bool(row["stationarity_screen_pass"])
            ]
            screen["first_passing_assessment_generation"] = (
                min(passing_endpoints) if passing_endpoints else None
            )
            screen["persistent_stationarity_generation"] = (
                _persistent_passing_generation(response_history)
            )
        run_window_rows.extend(run_windows)
        screen_rows.extend(screens)
        screen_history_rows.extend(cell_history)
    write_tsv(
        derived / "run-window-summaries.tsv",
        run_window_rows,
        list(run_window_rows[0]),
    )
    write_tsv(derived / "stationarity-screen.tsv", screen_rows, list(screen_rows[0]))
    write_tsv(
        derived / "stationarity-history.tsv",
        screen_history_rows,
        list(screen_history_rows[0]),
    )

    separated_rows = separated_stability_rows(trajectories_by_cell_seed, design)
    write_tsv(
        derived / "separated-stability-diagnostic.tsv",
        separated_rows,
        list(separated_rows[0]),
    )

    precision = precision_rows(trajectories_by_cell_seed, design)
    write_tsv(derived / "precision-recommendations.tsv", precision, list(precision[0]))

    passed_responses = sum(bool(row["stationarity_screen_pass"]) for row in screen_rows)
    cells_passing_all = sum(
        all(
            bool(row["stationarity_screen_pass"])
            for row in screen_rows
            if row["cell_id"] == cell_id
        )
        for cell_id in design
    )
    cell_stationarity_generation: dict[str, int | None] = {}
    for cell_id in design:
        generations = [
            row["persistent_stationarity_generation"]
            for row in screen_rows
            if row["cell_id"] == cell_id
        ]
        cell_stationarity_generation[cell_id] = (
            max(int(value) for value in generations)
            if generations and all(value is not None for value in generations)
            else None
        )
    summary = {
        "analysis_schema_version": "2.1.0",
        "analysis_scope": "equilibrium-and-precision pilot; exploratory",
        "experiment_id": EXPERIMENT_ID,
        "variant_id": VARIANT_TAG,
        "audit_status": "PASS",
        "populations": EXPECTED_RUNS,
        "cells": len(design),
        "seed_blocks": len(seed_blocks),
        "host_generations": HOST_GENERATIONS,
        "stationarity_responses_passing": passed_responses,
        "stationarity_responses_total": len(screen_rows),
        "cells_passing_all_stationarity_responses": cells_passing_all,
        "estimated_persistent_stationarity_generation_by_cell": (
            cell_stationarity_generation
        ),
        "full_equilibrium_status": (
            "not established because contrasting initial conditions were not tested"
        ),
        "separated_stability_diagnostic": _separated_stability_summary(
            separated_rows
        ),
        "precision_recommendation_maximum": max(
            int(row["recommended_replicates"]) for row in precision
        ),
        "precision_contrasts_exceeding_cap": sum(
            bool(row["exceeds_predeclared_cap"]) for row in precision
        ),
        "mutation_materialization": {
            "generated_lineages": total_generated_lineages,
            "largest_realized_count_per_transition": maximum_transition_mutants,
            "configured_limit_per_transition": 100_000,
        },
        "resources": {
            "summed_elapsed_hours": total_elapsed_seconds / 3600,
            "summed_output_gib": total_output_bytes / 1024**3,
            "maximum_peak_rss_mib": maximum_rss_kib / 1024,
        },
        "software_git_commits": sorted(commits),
        "source_sha256": sorted(source_hashes),
        "benchmark_platforms": sorted(platforms),
        "limitations": [
            "Stationarity diagnostics are an initial screen, not proof of equilibrium.",
            "Convergence from contrasting initial strain compositions is not tested.",
            (
                "Linear first-pilot cost projections are planning estimates, "
                "not guarantees."
            ),
            "Confirmatory seed blocks must remain held out until final design freeze.",
        ],
    }
    _atomic_json(derived / "analysis-summary.json", summary)
    _atomic_json(audit_path, _audit_payload("PASS", []))
    print(
        f"Second-pilot audit PASS: {EXPECTED_RUNS} populations; "
        f"{passed_responses}/{len(screen_rows)} response screens passed; "
        f"derived tables written to {derived}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
