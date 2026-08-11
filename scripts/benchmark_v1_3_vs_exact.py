#!/usr/bin/env python3
"""Compute benchmark for historical V1.3, neutral V1.4, V3.1, and exact counts.

The legacy V1.3 function always applies selection. Its V1.4 compatibility path
uses the same graph-based Wright--Fisher implementation with selection disabled,
so that path is labelled ``V1.3 neutral (legacy)`` in the benchmark.

All timed paths perform infection, within-host growth and steady generations,
escape sampling, host aggregation, and environmental return. File output is
excluded. Worker-pool startup is included because it is part of the legacy
host-generation implementation. The default remains the original one-generation
toy grid; command-line options expose the slide-13 and matched-neutral regimes.
"""

from __future__ import annotations

import argparse
import csv
import gc
import math
import time
from collections.abc import Callable, Iterable
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from trophosome.config import EvolutionConfig, HostConfig
from trophosome.count_model import (
    PopulationState,
    merge_populations,
    population_size_schedule,
    proportional_rescale,
)
from trophosome.simulation import (
    _batches,
    _bounded_host_map,
    _HostTask,
    _PopulationAccumulator,
    _reservoir_founders,
    _rng,
    _run_host_batch,
)

HISTORICAL_V1_3_LABEL = "V1.3 selected (historical)"
LEGACY_LABEL = "V1.3 neutral (legacy V1.4 switch)"
V3_1_LABEL = "V3.1 endpoint"
EXACT_LABEL = "Exact-count prototype"


@dataclass(frozen=True)
class Case:
    workers: int
    hosts: int
    mutation_probability: float
    carrying_capacity: int


def _legacy_environment(
    hosts: int,
    bottleneck: int,
    *,
    initial_strains: int = 1,
    seed: int = 666,
) -> object:
    import networkx as nx

    if initial_strains > 1:
        from project_package.generate_pop import generate_random_fisherlog_pop_unlinked

        return generate_random_fisherlog_pop_unlinked(i=initial_strains, seed=seed)

    population = nx.Graph()
    population.add_node(
        "0.0.0",
        abundance=max(1_000_000, hosts * bottleneck * 100),
        fitness=1.0,
    )
    return population


def run_legacy_case(
    case: Case,
    *,
    bottleneck: int,
    growth_factor: float,
    steady_generations: int,
    escape_fraction: float,
    host_generations: int = 1,
    initial_strains: int = 1,
    seed: int = 666,
) -> int:
    # Lazy import is important on spawn-based systems: exact-model workers must
    # not pay to import the legacy NetworkX, SciPy and plotting dependency tree.
    from project_package.run_model import run_generation_of_host_pop_v1_4

    environment = _legacy_environment(
        case.hosts,
        bottleneck,
        initial_strains=initial_strains,
        seed=seed,
    )
    touched = 0
    for host_generation in range(1, host_generations + 1):
        adults, environment, histories = run_generation_of_host_pop_v1_4(
            environment,
            n_worms=case.hosts,
            infection_sym_count=bottleneck,
            host_pop_gen=host_generation,
            escape_rate=escape_fraction,
            mutation_rate=case.mutation_probability,
            steady_state_runtime=steady_generations,
            max_runtime=10_000,
            growth_factor=growth_factor,
            pop_size_thr=case.carrying_capacity,
            stop_when_fixed=False,
            simplify=1,
            verbose=0,
            sampling_rate=1,
            nthreads=case.workers,
            intra_host_selection=False,
        )
        touched += len(adults) + len(environment)
        touched += sum(len(item) for item in histories)
    return touched


def run_historical_v1_3_case(
    case: Case,
    *,
    bottleneck: int,
    growth_factor: float,
    steady_generations: int,
    escape_fraction: float,
    host_generations: int = 20,
    initial_strains: int = 100,
    seed: int = 666,
) -> int:
    """Run the selected V1.3 path used by the slide-13 server benchmark."""

    from project_package.run_model import run_generation_of_host_pop_v1_3

    environment = _legacy_environment(
        case.hosts,
        bottleneck,
        initial_strains=initial_strains,
        seed=seed,
    )
    touched = 0
    for host_generation in range(1, host_generations + 1):
        adults, environment, histories = run_generation_of_host_pop_v1_3(
            environment,
            n_worms=case.hosts,
            infection_sym_count=bottleneck,
            host_pop_gen=host_generation,
            escape_rate=escape_fraction,
            mutation_rate=case.mutation_probability,
            steady_state_runtime=steady_generations,
            max_runtime=np.inf,
            growth_factor=growth_factor,
            pop_size_thr=case.carrying_capacity,
            stop_when_fixed=True,
            simplify=1,
            verbose=0,
            t=0,
            sampling_rate=1,
            nthreads=case.workers,
        )
        touched += len(adults) + len(environment)
        touched += sum(len(item) for item in histories)
    return touched


def run_v3_1_case(
    case: Case,
    *,
    bottleneck: int,
    growth_factor: float,
    steady_generations: int,
    escape_fraction: float,
    host_generations: int = 1,
    initial_strains: int = 1,
    seed: int = 666,
) -> int:
    from project_package.run_model import run_generation_of_host_pop_v3_1

    environment = _legacy_environment(
        case.hosts,
        bottleneck,
        initial_strains=initial_strains,
        seed=seed,
    )
    touched = 0
    for host_generation in range(1, host_generations + 1):
        adults, environment = run_generation_of_host_pop_v3_1(
            environment,
            n_worms=case.hosts,
            infection_sym_count=bottleneck,
            host_pop_gen=host_generation,
            escape_rate=escape_fraction,
            mutation_rate=case.mutation_probability,
            steady_state_runtime=steady_generations,
            growth_factor=growth_factor,
            pop_size_thr=case.carrying_capacity,
            nthreads=case.workers,
        )
        touched += len(adults) + len(environment)
    return touched


def run_exact_case(
    case: Case,
    *,
    bottleneck: int,
    growth_factor: float,
    steady_generations: int,
    escape_fraction: float,
    seed: int,
    host_generations: int = 1,
    initial_strains: int = 1,
) -> int:
    if host_generations != 1:
        raise ValueError("the exact comparison currently supports one host generation")
    host = HostConfig(
        population_size=case.hosts,
        infection_bottleneck=bottleneck,
        carrying_capacity=case.carrying_capacity,
        growth_factor=growth_factor,
        steady_generations=steady_generations,
        host_generations=1,
        escape_fraction=escape_fraction,
    )
    max_mutants = max(
        1_000,
        math.ceil(case.carrying_capacity * case.mutation_probability * 20),
    )
    evolution = EvolutionConfig(
        mutation_probability=case.mutation_probability,
        mutation_effect_mean=-0.01,
        mutation_effect_sd=0.01,
        within_host_selection=False,
        max_materialized_mutants=max_mutants,
    )
    legacy_environment = _legacy_environment(
        case.hosts,
        bottleneck,
        initial_strains=initial_strains,
        seed=seed,
    )
    environment_rows = [
        (attributes["abundance"], attributes["fitness"])
        for _, attributes in legacy_environment.nodes(data=True)
        if attributes["abundance"] > 0
    ]
    environment = PopulationState.from_counts(
        [row[0] for row in environment_rows],
        [row[1] for row in environment_rows],
    )
    cumulative = np.cumsum(environment.counts, dtype=np.float64)
    cumulative /= cumulative[-1]
    schedule = population_size_schedule(host)
    transitions = max(1, len(schedule))
    id_block_size = transitions * max_mutants
    first_new_strain_id = int(environment.genotype_ids.max()) + 1
    escape_count = round(case.carrying_capacity * escape_fraction)

    def tasks() -> Iterable[_HostTask]:
        for host_id in range(case.hosts):
            founders = _reservoir_founders(
                environment,
                cumulative,
                bottleneck,
                _rng(seed, 0, 1, host_id, 0),
            )
            id_start = first_new_strain_id + host_id * id_block_size
            yield _HostTask(
                replicate=0,
                host_generation=1,
                host_id=host_id,
                founders=founders,
                founder_depths=np.zeros(founders.richness, dtype=np.int64),
                host=host,
                evolution=evolution,
                size_schedule=schedule,
                escape_count=escape_count,
                id_start=id_start,
                id_stop=id_start + id_block_size,
                master_seed=seed,
            )

    task_batches = _batches(tasks(), max(1, min(8, case.hosts)))
    if case.workers == 1:
        result_batches = map(_run_host_batch, task_batches)
        executor = None
    else:
        try:
            executor = ProcessPoolExecutor(max_workers=case.workers)
        except (OSError, PermissionError):
            # Match the maintained runner's restricted-environment fallback.
            executor = ThreadPoolExecutor(max_workers=case.workers)
        result_batches = _bounded_host_map(
            executor,
            task_batches,
            max_in_flight=case.workers * 2,
        )

    adults = _PopulationAccumulator()
    releases = _PopulationAccumulator()
    event_count = 0
    try:
        for result_batch in result_batches:
            for result in result_batch:
                adults.add(result.adult)
                event_count += len(result.lineage_events)
                if result.escapees is not None:
                    releases.add(result.escapees)
    finally:
        if executor is not None:
            executor.shutdown(wait=True, cancel_futures=True)

    adult_population = adults.to_population()
    release_population = releases.to_population()
    assert adult_population is not None
    if release_population is None:
        updated_environment = environment
    else:
        updated_environment = proportional_rescale(
            merge_populations(environment, release_population),
            case.carrying_capacity,
        )
    return adult_population.richness + updated_environment.richness + event_count


def _comma_ints(value: str) -> tuple[int, ...]:
    return tuple(int(item) for item in value.split(","))


def _comma_floats(value: str) -> tuple[float, ...]:
    return tuple(float(item) for item in value.split(","))


def _measure(call: Callable[[], int]) -> float:
    gc.collect()
    started = time.perf_counter()
    touched = call()
    elapsed = time.perf_counter() - started
    if touched < 1:
        raise RuntimeError("benchmark produced an empty result")
    return elapsed


def _expected_mutation_workload(
    case: Case,
    *,
    model: str,
    bottleneck: int,
    growth_factor: float,
    steady_generations: int,
    escape_fraction: float,
    host_generations: int = 1,
) -> tuple[float, str]:
    host = HostConfig(
        population_size=case.hosts,
        infection_bottleneck=bottleneck,
        carrying_capacity=case.carrying_capacity,
        growth_factor=growth_factor,
        steady_generations=steady_generations,
        host_generations=1,
        escape_fraction=escape_fraction,
    )
    transitions = population_size_schedule(host)
    if model == V3_1_LABEL:
        probability = 1.0 - (1.0 - case.mutation_probability) ** len(transitions)
        return (
            host_generations * case.hosts * case.carrying_capacity * probability,
            "mutation-derived endpoint cells",
        )
    return (
        host_generations * case.hosts * sum(transitions) * case.mutation_probability,
        "forward mutation events",
    )


def _endpoint_mutation_load(
    case: Case,
    *,
    bottleneck: int,
    growth_factor: float,
    steady_generations: int,
    escape_fraction: float,
    host_generations: int = 1,
) -> float:
    """Return the shared endpoint mutation-supply axis used across models."""

    workload, _ = _expected_mutation_workload(
        case,
        model=V3_1_LABEL,
        bottleneck=bottleneck,
        growth_factor=growth_factor,
        steady_generations=steady_generations,
        escape_fraction=escape_fraction,
        host_generations=host_generations,
    )
    return workload


def _write_rows(path: Path, rows: list[dict[str, object]]) -> None:
    """Checkpoint the small raw benchmark table after every completed timing."""

    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--repeats", type=int, default=3)
    parser.add_argument("--models", default="legacy,v3_1,exact")
    parser.add_argument(
        "--dimensions",
        default=(
            "cpu_workers,individual_hosts,mutation_rate,within_host_population_size"
        ),
    )
    parser.add_argument("--workers", default="1,2,4")
    parser.add_argument("--hosts", default="1,2,4,8,16")
    parser.add_argument("--mutation-rates", default="0,0.0001,0.0005,0.001,0.005")
    parser.add_argument("--carrying-capacities", default="250,1000,4000,16000")
    parser.add_argument("--baseline-workers", type=int, default=1)
    parser.add_argument("--baseline-hosts", type=int, default=8)
    parser.add_argument("--baseline-mutation-rate", type=float, default=0.001)
    parser.add_argument("--baseline-carrying-capacity", type=int, default=4000)
    parser.add_argument("--bottleneck", type=int, default=10)
    parser.add_argument("--growth-factor", type=float, default=1.5)
    parser.add_argument("--steady-generations", type=int, default=5)
    parser.add_argument("--historical-v1-3-steady-generations", type=int)
    parser.add_argument("--v3-steady-generations", type=int)
    parser.add_argument("--escape-fraction", type=float, default=0.01)
    parser.add_argument("--host-generations", type=int, default=1)
    parser.add_argument("--initial-strains", type=int, default=1)
    parser.add_argument("--regime-label", default="one-generation toy")
    parser.add_argument("--seed", type=int, default=666)
    parser.add_argument("--max-expected-forward-mutations", type=float, default=1e6)
    parser.add_argument("--max-v3-mutation-derived-cells", type=float, default=1e5)
    args = parser.parse_args()
    if args.repeats < 1:
        parser.error("--repeats must be positive")
    if args.host_generations < 1:
        parser.error("--host-generations must be positive")
    if args.initial_strains < 1:
        parser.error("--initial-strains must be positive")

    all_dimensions = {
        "cpu_workers": _comma_ints(args.workers),
        "individual_hosts": _comma_ints(args.hosts),
        "mutation_rate": _comma_floats(args.mutation_rates),
        "within_host_population_size": _comma_ints(args.carrying_capacities),
    }
    selected_dimensions = tuple(args.dimensions.split(","))
    unknown = set(selected_dimensions) - set(all_dimensions)
    if unknown:
        parser.error(f"unknown dimensions: {', '.join(sorted(unknown))}")
    dimensions = {name: all_dimensions[name] for name in selected_dimensions}
    model_functions = {
        "v1_3": HISTORICAL_V1_3_LABEL,
        "legacy": LEGACY_LABEL,
        "v3_1": V3_1_LABEL,
        "exact": EXACT_LABEL,
    }
    selected_model_keys = tuple(args.models.split(","))
    unknown_models = set(selected_model_keys) - set(model_functions)
    if unknown_models:
        parser.error(f"unknown models: {', '.join(sorted(unknown_models))}")
    selected_models = tuple(model_functions[key] for key in selected_model_keys)
    if EXACT_LABEL in selected_models and args.host_generations != 1:
        parser.error("the exact comparison currently supports --host-generations=1")

    def steady_for(model: str) -> int:
        if (
            model == HISTORICAL_V1_3_LABEL
            and args.historical_v1_3_steady_generations is not None
        ):
            return args.historical_v1_3_steady_generations
        if model == V3_1_LABEL and args.v3_steady_generations is not None:
            return args.v3_steady_generations
        return args.steady_generations

    def case_for(dimension: str, value: int | float) -> Case:
        values: dict[str, int | float] = {
            "cpu_workers": args.baseline_workers,
            "individual_hosts": args.baseline_hosts,
            "mutation_rate": args.baseline_mutation_rate,
            "within_host_population_size": args.baseline_carrying_capacity,
        }
        values[dimension] = value
        return Case(
            workers=int(values["cpu_workers"]),
            hosts=int(values["individual_hosts"]),
            mutation_probability=float(values["mutation_rate"]),
            carrying_capacity=int(values["within_host_population_size"]),
        )

    # Warm imports and NumPy dispatch outside the recorded timings.
    warm = Case(workers=1, hosts=1, mutation_probability=0.001, carrying_capacity=50)
    if EXACT_LABEL in selected_models:
        run_exact_case(
            warm,
            bottleneck=min(args.bottleneck, 10),
            growth_factor=args.growth_factor,
            steady_generations=1,
            escape_fraction=args.escape_fraction,
            seed=args.seed,
        )
    if LEGACY_LABEL in selected_models:
        run_legacy_case(
            warm,
            bottleneck=min(args.bottleneck, 10),
            growth_factor=args.growth_factor,
            steady_generations=1,
            escape_fraction=args.escape_fraction,
        )
    if HISTORICAL_V1_3_LABEL in selected_models:
        run_historical_v1_3_case(
            warm,
            bottleneck=min(args.bottleneck, 10),
            growth_factor=args.growth_factor,
            steady_generations=1,
            escape_fraction=args.escape_fraction,
            host_generations=1,
            initial_strains=1,
            seed=args.seed,
        )
    if V3_1_LABEL in selected_models:
        run_v3_1_case(
            warm,
            bottleneck=min(args.bottleneck, 10),
            growth_factor=args.growth_factor,
            steady_generations=1,
            escape_fraction=args.escape_fraction,
        )

    rows: list[dict[str, object]] = []
    for dimension, values in dimensions.items():
        for value in values:
            case = case_for(dimension, value)
            for repeat in range(args.repeats):
                offset = repeat % len(selected_models)
                ordered_models = selected_models[offset:] + selected_models[:offset]
                for model in ordered_models:
                    workload, workload_kind = _expected_mutation_workload(
                        case,
                        model=model,
                        bottleneck=args.bottleneck,
                        growth_factor=args.growth_factor,
                        steady_generations=steady_for(model),
                        escape_fraction=args.escape_fraction,
                        host_generations=args.host_generations,
                    )
                    endpoint_mutation_load = _endpoint_mutation_load(
                        case,
                        bottleneck=args.bottleneck,
                        growth_factor=args.growth_factor,
                        steady_generations=steady_for(model),
                        escape_fraction=args.escape_fraction,
                        host_generations=args.host_generations,
                    )
                    if model == HISTORICAL_V1_3_LABEL:
                        call = lambda: run_historical_v1_3_case(
                            case,
                            bottleneck=args.bottleneck,
                            growth_factor=args.growth_factor,
                            steady_generations=steady_for(model),
                            escape_fraction=args.escape_fraction,
                            host_generations=args.host_generations,
                            initial_strains=args.initial_strains,
                            seed=args.seed + repeat,
                        )
                    elif model == LEGACY_LABEL:
                        call = lambda: run_legacy_case(
                            case,
                            bottleneck=args.bottleneck,
                            growth_factor=args.growth_factor,
                            steady_generations=steady_for(model),
                            escape_fraction=args.escape_fraction,
                            host_generations=args.host_generations,
                            initial_strains=args.initial_strains,
                            seed=args.seed + repeat,
                        )
                    elif model == V3_1_LABEL:
                        call = lambda: run_v3_1_case(
                            case,
                            bottleneck=args.bottleneck,
                            growth_factor=args.growth_factor,
                            steady_generations=steady_for(model),
                            escape_fraction=args.escape_fraction,
                            host_generations=args.host_generations,
                            initial_strains=args.initial_strains,
                            seed=args.seed + repeat,
                        )
                    else:
                        call = lambda: run_exact_case(
                            case,
                            bottleneck=args.bottleneck,
                            growth_factor=args.growth_factor,
                            steady_generations=steady_for(model),
                            escape_fraction=args.escape_fraction,
                            seed=args.seed + repeat,
                            host_generations=args.host_generations,
                            initial_strains=args.initial_strains,
                        )
                    status = "ok"
                    error = ""
                    workload_limit = (
                        args.max_v3_mutation_derived_cells
                        if model == V3_1_LABEL
                        else args.max_expected_forward_mutations
                    )
                    if workload > workload_limit:
                        elapsed = ""
                        status = "preflight_rejected"
                        error = (
                            f"estimated {workload_kind} {workload:,.0f} exceed "
                            f"configured limit {workload_limit:,.0f}"
                        )
                    else:
                        try:
                            elapsed = _measure(call)
                        except Exception as exc:
                            elapsed = ""
                            status = "failed"
                            error = f"{type(exc).__name__}: {str(exc).splitlines()[0]}"
                    rows.append(
                        {
                            "dimension": dimension,
                            "x_value": value,
                            "model": model,
                            "repeat": repeat,
                            "runtime_seconds": elapsed,
                            "status": status,
                            "error": error,
                            "regime": args.regime_label,
                            "estimated_mutation_workload": workload,
                            "mutation_workload_kind": workload_kind,
                            "endpoint_mutation_load": endpoint_mutation_load,
                            "workers": case.workers,
                            "individual_hosts": case.hosts,
                            "mutation_rate": case.mutation_probability,
                            "within_host_population_size": case.carrying_capacity,
                            "infection_bottleneck": args.bottleneck,
                            "growth_factor": args.growth_factor,
                            "steady_generations": steady_for(model),
                            "host_generations": args.host_generations,
                            "initial_strains": args.initial_strains,
                            "escape_fraction": args.escape_fraction,
                        }
                    )
                    _write_rows(args.output, rows)
                    if status == "ok":
                        print(
                            f"{dimension}={value} | {model} | repeat {repeat + 1}: "
                            f"{elapsed:.6f} s",
                            flush=True,
                        )
                    else:
                        print(
                            f"{dimension}={value} | {model} | repeat {repeat + 1}: "
                            f"{status.upper()} ({error})",
                            flush=True,
                        )

    _write_rows(args.output, rows)


if __name__ == "__main__":
    main()
