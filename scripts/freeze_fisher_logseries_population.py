#!/usr/bin/env python3
"""Generate or verify the frozen Phase 1 Fisher log-series population."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path

import numpy as np

DEFAULT_OUTPUT = Path(
    "experiments/work/trophosome/common/initial-populations/"
    "ip001-fisher100.json"
)


def _scaled_counts(raw_draws: np.ndarray, capacity: int) -> list[int]:
    """Scale positive weights exactly using deterministic Hamilton apportionment."""

    total = int(raw_draws.sum())
    floors = [(int(value) * capacity) // total for value in raw_draws]
    remainders = [(int(value) * capacity) % total for value in raw_draws]
    remaining = capacity - sum(floors)
    order = sorted(range(len(floors)), key=lambda index: (-remainders[index], index))
    for index in order[:remaining]:
        floors[index] += 1
    return floors


def _count_checksum(counts: list[int]) -> str:
    canonical = (",".join(str(value) for value in counts) + "\n").encode()
    return hashlib.sha256(canonical).hexdigest()


def build_payload() -> dict[str, object]:
    seed = 666
    fisher_a = 0.995
    richness = 100
    capacity = 1_000_000_000
    raw_draws = np.random.RandomState(seed).logseries(fisher_a, richness)
    counts = _scaled_counts(raw_draws, capacity)
    frequencies = np.asarray(counts, dtype=np.float64) / capacity
    shannon = float(-np.sum(frequencies * np.log(frequencies)))
    squared_sum = float(np.sum(frequencies * frequencies))
    return {
        "initial_population_schema_version": "1.0.0",
        "initial_population_id": "ip001-fisher100",
        "description": (
            "Frozen neutral 100-strain Fisher log-series realization for "
            "Phase 1 experiments"
        ),
        "generator": {
            "distribution": "numpy.random.RandomState.logseries",
            "random_number_generator": "MT19937 RandomState",
            "seed": seed,
            "fisher_logseries_a": fisher_a,
            "draw_count": richness,
        },
        "scaling": {
            "environmental_capacity": capacity,
            "method": "Hamilton largest remainder",
            "tie_break": "ascending zero-based strain_id",
            "canonical_checksum_encoding": "comma-separated counts plus LF",
        },
        "strain_id_order": "zero-based draw order, 0 through 99",
        "raw_logseries_draws": [int(value) for value in raw_draws],
        "scaled_counts": counts,
        "scaled_counts_sha256": _count_checksum(counts),
        "neutral_fitness": {
            "within_host": 1.0,
            "free_living": 1.0,
            "mutation_effect_mean": 0.0,
            "mutation_effect_sd": 0.0,
        },
        "summary": {
            "capacity": sum(counts),
            "richness": richness,
            "minimum_count": min(counts),
            "maximum_count": max(counts),
            "hill_q0": richness,
            "hill_q1": math.exp(shannon),
            "hill_q2": 1.0 / squared_sum,
            "shannon": shannon,
            "pielou_evenness": shannon / math.log(richness),
            "gene_diversity": 1.0 - squared_sum,
        },
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument(
        "--verify",
        action="store_true",
        help="Compare an existing artifact with a fresh deterministic realization",
    )
    args = parser.parse_args()
    expected = build_payload()
    if args.verify:
        observed = json.loads(args.output.read_text(encoding="utf-8"))
        if observed != expected:
            raise SystemExit(f"frozen population does not match: {args.output}")
        print(f"Verified {args.output}")
        print(f"Counts SHA-256: {expected['scaled_counts_sha256']}")
        return 0

    if args.output.exists():
        raise SystemExit(f"refusing to overwrite existing artifact: {args.output}")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(expected, indent=2) + "\n", encoding="utf-8")
    print(f"Wrote {args.output}")
    print(f"Counts SHA-256: {expected['scaled_counts_sha256']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
