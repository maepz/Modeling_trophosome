#!/usr/bin/env python3
"""Benchmark the exact count kernel without retaining trajectories."""

from __future__ import annotations

import argparse
import csv
import time
import tracemalloc
from pathlib import Path

import numpy as np

from trophosome.config import EvolutionConfig
from trophosome.count_model import IdAllocator, PopulationState, wright_fisher_step


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--richness", default="10,100,1000")
    parser.add_argument("--population-size", type=int, default=1_000_000_000)
    parser.add_argument("--generations", type=int, default=20)
    parser.add_argument("--mutation-probability", type=float, default=0.0)
    parser.add_argument("--max-materialized-mutants", type=int, default=100_000)
    parser.add_argument("--seed", type=int, default=666)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    rows = []
    for richness in (int(value) for value in args.richness.split(",")):
        base, remainder = divmod(args.population_size, richness)
        counts = np.full(richness, base, dtype=np.int64)
        counts[:remainder] += 1
        population = PopulationState.from_counts(
            counts.tolist(), np.ones(richness).tolist()
        )
        allocator = IdAllocator(richness)
        evolution = EvolutionConfig(
            mutation_probability=args.mutation_probability,
            max_materialized_mutants=args.max_materialized_mutants,
        )
        rng = np.random.default_rng(args.seed)
        tracemalloc.start()
        started = time.perf_counter()
        new_mutants = 0
        for _ in range(args.generations):
            population, mutants = wright_fisher_step(
                population, args.population_size, evolution, rng, allocator
            )
            new_mutants += mutants
        elapsed = time.perf_counter() - started
        _, peak = tracemalloc.get_traced_memory()
        tracemalloc.stop()
        rows.append(
            {
                "initial_richness": richness,
                "final_richness": population.richness,
                "population_size": population.size,
                "generations": args.generations,
                "mutation_probability": args.mutation_probability,
                "new_mutants": new_mutants,
                "elapsed_seconds": elapsed,
                "peak_tracemalloc_mb": peak / 1024**2,
            }
        )

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


if __name__ == "__main__":
    main()
