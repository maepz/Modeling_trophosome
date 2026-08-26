from __future__ import annotations

import csv
import importlib.util
import json
import tempfile
import unittest
from pathlib import Path

from trophosome.second_pilot_report import generate_second_pilot_report

REPOSITORY = Path(__file__).resolve().parents[1]
DESIGN = (
    REPOSITORY
    / "experiments/work/trophosome/p01-neutral-feedback/design"
    / "phase1-second-pilot-v210-m010-g250-cells.tsv"
)
REPORT_DEPS = all(
    importlib.util.find_spec(module) is not None
    for module in ("matplotlib", "numpy", "reportlab")
)


def _write_tsv(path: Path, rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)


def _fixture(analysis: Path) -> None:
    cell_ids = [f"p01-s02-c{number:04d}" for number in range(21, 27)]
    trajectories: list[dict[str, object]] = []
    for cell_number, cell_id in enumerate(cell_ids):
        for seed in range(1, 13):
            for generation in range(251):
                wave = 0.15 * ((generation + seed) % 7 - 3)
                trajectories.append(
                    {
                        "run_id": f"{cell_id}-sb{seed:04d}",
                        "cell_id": cell_id,
                        "cell": cell_id.rsplit("-", 1)[-1],
                        "seed_block_id": f"sb{seed:04d}",
                        "generation": generation,
                        "feedback_exposure": generation * 0.69,
                        "D0": 100 + (1 if cell_number == 2 and generation > 100 else 0),
                        "shannon": 3.4,
                        "D1": 30 + cell_number + wave,
                        "D2": 20 + 0.5 * cell_number + wave,
                        "evenness": 0.74 + 0.005 * cell_number,
                        "top_frequency": 0.1,
                        "TV": 0.02 + 0.01 * cell_number + 0.0001 * generation,
                        "jensen_shannon": 0.001,
                        "turnover": 0.002,
                        "mutant_richness": 0,
                        "mutant_abundance_fraction": 0,
                        "root_D0": 100,
                        "root_D1": 30,
                        "root_D2": 20,
                        "root_evenness": 0.74,
                        "mean_adult_richness": 8,
                        "mean_adult_gene_diversity": 0.7,
                        "realized_host_feedback": 0.5,
                        "expected_host_feedback_after_migration": 0.45,
                    }
                )
    windows: list[dict[str, object]] = []
    for cell_id in cell_ids:
        for seed in range(1, 13):
            for response in ("D0", "D1", "D2", "evenness", "TV"):
                for window in range(1, 5):
                    windows.append(
                        {
                            "cell_id": cell_id,
                            "seed_block_id": f"sb{seed:04d}",
                            "response": response,
                            "window": window,
                            "start_generation": 171 + 20 * (window - 1),
                            "end_generation": 190 + 20 * (window - 1),
                            "mean": 1,
                            "median": 1,
                            "sd": 0.02,
                            "cv": 0.01 + 0.002 * (int(cell_id[-2:]) - 21),
                            "q05": 0.98,
                            "q95": 1.02,
                            "lag1_autocorrelation": 0.1,
                            "integrated_autorrelation_time": 1.2,
                        }
                    )
    stationarity: list[dict[str, object]] = []
    for cell_number, cell_id in enumerate(cell_ids):
        for response_number, response in enumerate(
            ("D0", "D1", "D2", "evenness", "TV")
        ):
            passed = not (cell_number == 5 and response_number == 4)
            stationarity.append(
                {
                    "cell_id": cell_id,
                    "response": response,
                    "assessment_end_generation": 250,
                    "window_generations": 20,
                    "previous_assessment_equivalent": passed,
                    "current_assessment_equivalent": passed,
                    "largest_previous_ci_to_margin_ratio": 0.5,
                    "largest_current_ci_to_margin_ratio": 0.6,
                    "equivalence_margin": 0.05,
                    "rank_normalized_split_rhat": 1.01 if passed else 1.08,
                    "approximate_combined_ess": 650 if passed else 320,
                    "stationarity_screen_pass": passed,
                    "first_passing_assessment_generation": 200 if passed else "",
                    "persistent_stationarity_generation": 200 if passed else "",
                    "full_equilibrium_status": "not established",
                }
            )
    contrasts = (
        "host-passage",
        "mutation",
        "many-vs-few-hosts",
        "weak-vs-baseline-feedback",
        "strong-vs-baseline-feedback",
    )
    precision: list[dict[str, object]] = []
    for contrast_number, contrast in enumerate(contrasts):
        for response in ("D1", "D2", "TV"):
            precision.append(
                {
                    "contrast": contrast,
                    "treatment_cell": cell_ids[min(contrast_number + 1, 5)],
                    "reference_cell": cell_ids[0],
                    "response": response,
                    "scale": "relative" if response != "TV" else "absolute",
                    "pilot_mean_difference": 0.02,
                    "pilot_sd_difference": 0.08,
                    "pilot_ci90_lower": -0.02,
                    "pilot_ci90_upper": 0.06,
                    "desired_ci_half_width": 0.05,
                    "formula_required_replicates": 16,
                    "recommended_replicates": 20 + 8 * contrast_number,
                    "exceeds_predeclared_cap": False,
                }
            )
    summary = {
        "analysis_schema_version": "2.0.0",
        "audit_status": "PASS",
        "populations": 72,
        "cells": 6,
        "seed_blocks": 12,
        "host_generations": 250,
        "stationarity_responses_passing": 29,
        "stationarity_responses_total": 30,
        "cells_passing_all_stationarity_responses": 5,
        "full_equilibrium_status": "not established",
        "resources": {
            "summed_elapsed_hours": 33.2,
            "summed_output_gib": 7.1,
            "maximum_peak_rss_mib": 180.0,
        },
        "limitations": [
            "Stationarity diagnostics are an initial screen, not proof of equilibrium.",
            "Contrasting initial communities were not tested.",
        ],
    }
    analysis.mkdir(parents=True, exist_ok=True)
    (analysis / "analysis-summary.json").write_text(
        json.dumps(summary), encoding="utf-8"
    )
    _write_tsv(analysis / "environment-trajectories.tsv", trajectories)
    _write_tsv(analysis / "run-window-summaries.tsv", windows)
    _write_tsv(analysis / "stationarity-screen.tsv", stationarity)
    _write_tsv(analysis / "precision-recommendations.tsv", precision)


@unittest.skipUnless(REPORT_DEPS, "optional report dependencies are not installed")
class Phase1SecondPilotReportingTests(unittest.TestCase):
    def test_builds_self_contained_pdf_and_editable_companion(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            analysis = root / "analysis"
            _fixture(analysis)
            artifacts = generate_second_pilot_report(
                analysis=analysis,
                design=DESIGN,
                output=root / "second-pilot.pdf",
                markdown=root / "second-pilot.md",
                assets=root / "figures",
                report_date="2026-08-27",
            )
            self.assertTrue(artifacts.pdf.read_bytes().startswith(b"%PDF-"))
            self.assertGreater(artifacts.pdf.stat().st_size, 100_000)
            markdown = artifacts.markdown.read_text(encoding="utf-8")
            self.assertIn("not definitive equilibrium", markdown)
            self.assertIn("| 10,000 |", markdown)
            self.assertIn("Stationarity screen", markdown)
            self.assertIn("Short glossary", markdown)
            self.assertEqual(len(list(artifacts.assets.glob("*.png"))), 4)


if __name__ == "__main__":
    unittest.main()
