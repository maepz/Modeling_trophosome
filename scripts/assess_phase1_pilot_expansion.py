#!/usr/bin/env python3
"""Assess whether the eight-cell Phase 1 pilot extension can run safely."""

from __future__ import annotations

import argparse
import csv
import json
import statistics
from pathlib import Path
from typing import Any


def _read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def _median(rows: list[dict[str, str]], measure: str) -> float:
    return statistics.median(float(row[measure]) for row in rows)


def _cell_rows(endpoints: list[dict[str, str]], cell_id: str) -> list[dict[str, str]]:
    return [row for row in endpoints if row["cell_id"] == cell_id]


def _gate(name: str, passed: bool, evidence: str) -> dict[str, Any]:
    return {"name": name, "passed": passed, "evidence": evidence}


def assess(
    analysis: Path,
    full_design: Path,
    core_report: Path,
    *,
    runtime_budget_hours: float,
    storage_budget_gib: float,
    budget_fraction: float,
) -> dict[str, Any]:
    summary = json.loads((analysis / "analysis-summary.json").read_text())
    endpoints = _read_tsv(analysis / "run-endpoints.tsv")
    design = _read_tsv(full_design)
    design_by_cell = {row["cell_id"]: row for row in design}
    gates: list[dict[str, Any]] = []

    gates.append(
        _gate(
            "Core analysis and output audit",
            summary.get("audit_status") == "PASS"
            and int(summary.get("cells", 0)) == 12
            and int(summary.get("populations", 0)) == 36,
            f"audit={summary.get('audit_status')}; cells={summary.get('cells')}; "
            f"populations={summary.get('populations')}",
        )
    )
    gates.append(
        _gate(
            "Core PDF report generated",
            core_report.is_file() and core_report.stat().st_size > 0,
            str(core_report),
        )
    )

    no_return = [row for row in endpoints if float(row["R"]) == 0]
    maximum_no_return_tv = max(float(row["total_variation"]) for row in no_return)
    gates.append(
        _gate(
            "No-return environment unchanged",
            bool(no_return) and maximum_no_return_tv == 0,
            f"maximum total-variation distance={maximum_no_return_tv:g}",
        )
    )

    mutation_fraction = float(
        summary["mutation_materialization"]["largest_fraction_of_limit"]
    )
    gates.append(
        _gate(
            "Mutation materialization below 25% of limit",
            mutation_fraction < 0.25,
            f"largest fraction={mutation_fraction:.6%}",
        )
    )

    matched_ids = [f"p01-s01-c{number:04d}" for number in range(3, 7)]
    matched_values = {
        (
            int(design_by_cell[cell_id]["R"]),
            float(design_by_cell[cell_id]["alpha"]),
        )
        for cell_id in matched_ids
    }
    gates.append(
        _gate(
            "Matched-return core has equal R and alpha",
            matched_values == {(1_000_000_000, 0.5)},
            f"observed={sorted(matched_values)}",
        )
    )

    mutation_ids = [f"p01-s01-c{number:04d}" for number in range(7, 11)]
    mutant_richness = [
        _median(_cell_rows(endpoints, cell_id), "final_mutant_richness")
        for cell_id in mutation_ids
    ]
    maximum_mutant_fraction = max(
        _median(_cell_rows(endpoints, cell_id), "mutant_abundance_fraction")
        for cell_id in mutation_ids
    )
    informative_mutation = (
        min(mutant_richness) <= 5
        and max(mutant_richness) >= 10
        and maximum_mutant_fraction < 0.01
    )
    gates.append(
        _gate(
            "Mutation levels bracket informative novelty",
            informative_mutation,
            f"median final mutant richness={mutant_richness}; maximum median mutant "
            f"abundance fraction={maximum_mutant_fraction:.3g}",
        )
    )

    baseline_by_h: dict[int, dict[str, float]] = {}
    for cell_id in matched_ids:
        rows = _cell_rows(endpoints, cell_id)
        host_number = int(float(rows[0]["H"]))
        baseline_by_h[host_number] = {
            measure: _median(rows, measure)
            for measure in ("elapsed_minutes", "output_mib", "peak_rss_mib")
        }
    baseline_rows = _cell_rows(endpoints, "p01-s01-c0003")
    informative_rows = _cell_rows(endpoints, "p01-s01-c0009")
    mutation_multipliers = {
        measure: _median(informative_rows, measure) / _median(baseline_rows, measure)
        for measure in ("elapsed_minutes", "output_mib", "peak_rss_mib")
    }

    projections: list[dict[str, Any]] = []
    for number in range(13, 21):
        cell_id = f"p01-s01-c{number:04d}"
        row = design_by_cell[cell_id]
        host_number = int(row["H"])
        mutation_enabled = float(row["u"]) > 0
        multiplier = (
            mutation_multipliers
            if mutation_enabled
            else {
                measure: 1.0
                for measure in ("elapsed_minutes", "output_mib", "peak_rss_mib")
            }
        )
        baseline = baseline_by_h[host_number]
        projections.append(
            {
                "cell_id": cell_id,
                "H": host_number,
                "u": float(row["u"]),
                "populations": 3,
                "projected_elapsed_minutes": (
                    3 * baseline["elapsed_minutes"] * multiplier["elapsed_minutes"]
                ),
                "projected_output_mib": (
                    3 * baseline["output_mib"] * multiplier["output_mib"]
                ),
                "projected_peak_rss_mib_per_population": (
                    baseline["peak_rss_mib"] * multiplier["peak_rss_mib"]
                ),
            }
        )

    projected_hours = sum(row["projected_elapsed_minutes"] for row in projections) / 60
    projected_gib = sum(row["projected_output_mib"] for row in projections) / 1024
    projected_peak_rss_mib = max(
        row["projected_peak_rss_mib_per_population"] for row in projections
    )
    runtime_limit = runtime_budget_hours * budget_fraction
    storage_limit = storage_budget_gib * budget_fraction
    gates.append(
        _gate(
            "Projected extension runtime below budget gate",
            projected_hours <= runtime_limit,
            f"projected sequential runtime={projected_hours:.2f} h; "
            f"gate={runtime_limit:.2f} h",
        )
    )
    gates.append(
        _gate(
            "Projected extension storage below budget gate",
            projected_gib <= storage_limit,
            f"projected diagnostic output={projected_gib:.2f} GiB; "
            f"gate={storage_limit:.2f} GiB",
        )
    )

    return {
        "assessment_schema_version": "1.0.0",
        "core_experiment_id": summary.get("experiment_id"),
        "extension_cells": [f"p01-s01-c{number:04d}" for number in range(13, 21)],
        "budgets": {
            "runtime_hours": runtime_budget_hours,
            "storage_gib": storage_budget_gib,
            "maximum_fraction": budget_fraction,
        },
        "projection": {
            "sequential_runtime_hours": projected_hours,
            "diagnostic_output_gib": projected_gib,
            "peak_rss_mib_per_population": projected_peak_rss_mib,
            "mutation_cost_multipliers": mutation_multipliers,
            "cells": projections,
        },
        "gates": gates,
        "safe_to_expand": all(gate["passed"] for gate in gates),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--analysis", type=Path, required=True)
    parser.add_argument("--design", type=Path, required=True)
    parser.add_argument("--core-report", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--runtime-budget-hours", type=float, default=48.0)
    parser.add_argument("--storage-budget-gib", type=float, default=100.0)
    parser.add_argument("--budget-fraction", type=float, default=0.70)
    args = parser.parse_args()

    result = assess(
        args.analysis.resolve(),
        args.design.resolve(),
        args.core_report.resolve(),
        runtime_budget_hours=args.runtime_budget_hours,
        storage_budget_gib=args.storage_budget_gib,
        budget_fraction=args.budget_fraction,
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    for gate in result["gates"]:
        status = "PASS" if gate["passed"] else "FAIL"
        print(f"[{status}] {gate['name']}: {gate['evidence']}")
    print(
        "Extension decision: "
        + ("SAFE TO RUN" if result["safe_to_expand"] else "DO NOT RUN")
    )
    return 0 if result["safe_to_expand"] else 2


if __name__ == "__main__":
    raise SystemExit(main())
