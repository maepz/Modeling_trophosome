from __future__ import annotations

import importlib.util
import sys
import tempfile
import unittest
from pathlib import Path

from trophosome.pilot_report import (
    _design_table_rows,
    _effect_summary_rows,
    _key_findings,
    _plot_migration,
    _read_tsv,
)

REPOSITORY = Path(__file__).resolve().parents[1]
WORK = REPOSITORY / "experiments" / "work" / "trophosome"
PHASE = WORK / "p01-neutral-feedback"
ANALYSIS = PHASE / "analysis" / "s01-pilot-core-derived"
DESIGN = PHASE / "design" / "phase1-first-pilot-v210-m010-core-cells.tsv"
RUNS = PHASE / "manifests" / "phase1-first-pilot-v210-m010-runs.tsv"
ANALYSIS_SCRIPT = PHASE / "analysis" / "analyse_first_pilot_v2_1.py"


def _load_analysis_module():
    sys.path.insert(0, str(ANALYSIS_SCRIPT.parent))
    try:
        spec = importlib.util.spec_from_file_location(
            "analyse_first_pilot_v2_1_test", ANALYSIS_SCRIPT
        )
        if spec is None or spec.loader is None:
            raise RuntimeError("could not load v2.1 analysis module")
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        return module
    finally:
        sys.path.pop(0)


class Phase1V21ReportingTests(unittest.TestCase):
    def test_migration_design_and_biological_interpretation(self) -> None:
        design = _read_tsv(DESIGN)
        endpoints = _read_tsv(ANALYSIS / "run-endpoints.tsv")
        table = _design_table_rows(design)
        self.assertEqual(table[0][-1], "Migration m")
        self.assertEqual({row[-1] for row in table[1:]}, {"0.1"})

        effects = _effect_summary_rows(endpoints, design)
        self.assertIn("Migration-only baseline", effects[1][0])
        self.assertIn("stochastic exchange", effects[1][-1])
        findings = _key_findings(
            {"initial_metrics": {}},
            endpoints,
            design,
        )
        self.assertTrue(any("m=0.10" in finding for finding in findings))
        self.assertTrue(
            any("migration-only baseline" in finding for finding in findings)
        )

    @unittest.skipUnless(
        importlib.util.find_spec("matplotlib") is not None,
        "optional plotting dependency is not installed",
    )
    def test_migration_diagnostic_figure_is_created(self) -> None:
        design = _read_tsv(DESIGN)
        trajectories = _read_tsv(ANALYSIS / "environment-trajectories.tsv")
        alpha_by_cell = {row["cell_id"]: float(row["alpha"]) for row in design}
        for row in trajectories:
            if int(row["generation"]) == 0:
                row.update(
                    {
                        "post_return_richness": "",
                        "post_migration_richness": "",
                        "realized_host_feedback": "",
                        "expected_host_feedback_after_migration": "",
                    }
                )
                continue
            alpha = alpha_by_cell[row["cell_id"]]
            richness = row["environment_richness"]
            row.update(
                {
                    "post_return_richness": richness,
                    "post_migration_richness": richness,
                    "realized_host_feedback": str(alpha),
                    "expected_host_feedback_after_migration": str(0.9 * alpha),
                }
            )
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory) / "migration.png"
            _plot_migration(trajectories, design, output)
            self.assertTrue(output.is_file())
            self.assertGreater(output.stat().st_size, 10_000)

    def test_completion_gate_lists_missing_runs_before_analysis(self) -> None:
        module = _load_analysis_module()
        rows = _read_tsv(RUNS)
        with tempfile.TemporaryDirectory() as directory:
            issues = module.completion_gate_issues(
                rows,
                work=WORK,
                scratch=Path(directory),
            )
        self.assertEqual(len(rows), 60)
        self.assertGreaterEqual(len(issues), 60)
        self.assertTrue(any("completion.json" in issue for issue in issues))


if __name__ == "__main__":
    unittest.main()
