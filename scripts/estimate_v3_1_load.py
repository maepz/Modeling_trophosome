#!/usr/bin/env python3
"""Print a pre-run mutation-load estimate for the V3.1 prototype."""

from __future__ import annotations

import argparse
import json

from trophosome.v3_1_diagnostics import estimate_v3_1_load


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Estimate mutant-cell and strain-richness workloads in legacy V3.1 "
            "without running a simulation"
        )
    )
    parser.add_argument("--infection-size", type=int, required=True)
    parser.add_argument("--carrying-capacity", type=int, required=True)
    parser.add_argument("--growth-factor", type=float, required=True)
    parser.add_argument("--steady-generations", type=int, required=True)
    parser.add_argument("--mutation-probability", type=float, required=True)
    return parser


def main() -> int:
    args = build_parser().parse_args()
    estimate = estimate_v3_1_load(
        infection_size=args.infection_size,
        carrying_capacity=args.carrying_capacity,
        growth_factor=args.growth_factor,
        steady_generations=args.steady_generations,
        mutation_probability=args.mutation_probability,
    )
    print(json.dumps(estimate.to_dict(), indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
