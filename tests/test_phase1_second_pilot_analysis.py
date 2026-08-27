from __future__ import annotations

import importlib.util
import json
import sys
import tempfile
import unittest
from pathlib import Path

from trophosome.config import load_config

REPOSITORY = Path(__file__).resolve().parents[1]
WORK = REPOSITORY / "experiments/work/trophosome"
PHASE = WORK / "p01-neutral-feedback"
ANALYSIS_SCRIPT = PHASE / "analysis/analyse_second_pilot.py"
RUNS = PHASE / "manifests/phase1-second-pilot-v210-m010-g250-runs.tsv"


def _load_module():
    sys.path.insert(0, str(ANALYSIS_SCRIPT.parent))
    try:
        spec = importlib.util.spec_from_file_location(
            "analyse_second_pilot_test", ANALYSIS_SCRIPT
        )
        if spec is None or spec.loader is None:
            raise RuntimeError("could not load the second-pilot analysis")
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        return module
    finally:
        sys.path.pop(0)


def _trajectory(
    seed: int,
    *,
    trend_d1: bool = False,
    d1_seed_offset: float = 0.0,
) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for generation in range(251):
        rows.append(
            {
                "seed_block_id": f"sb{seed:04d}",
                "generation": generation,
                "D0": 100.0,
                "D1": (
                    30.0 + d1_seed_offset + generation
                    if trend_d1
                    else 30.0 + d1_seed_offset
                ),
                "D2": 20.0,
                "evenness": 0.75,
                "TV": 0.1,
            }
        )
    return rows


class Phase1SecondPilotAnalysisTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.module = _load_module()

    def test_constant_late_trajectories_pass_the_stationarity_screen(self) -> None:
        _windows, screens = self.module.equilibrium_screen(
            "p01-s02-c0022",
            [_trajectory(seed) for seed in range(1, 13)],
            window_length=20,
        )
        self.assertEqual(len(screens), 5)
        self.assertTrue(all(row["stationarity_screen_pass"] for row in screens))
        self.assertTrue(all(row["approximate_combined_ess"] >= 400 for row in screens))
        self.assertTrue(
            all("not established" in row["full_equilibrium_status"] for row in screens)
        )

    def test_a_continuing_d1_trend_fails(self) -> None:
        _windows, screens = self.module.equilibrium_screen(
            "p01-s02-c0022",
            [_trajectory(seed, trend_d1=True) for seed in range(1, 13)],
            window_length=20,
        )
        by_response = {row["response"]: row for row in screens}
        self.assertFalse(by_response["D1"]["stationarity_screen_pass"])
        self.assertTrue(by_response["D2"]["stationarity_screen_pass"])

    def test_closure_analysis_accepts_twenty_matched_seed_blocks(self) -> None:
        _windows, screens = self.module.equilibrium_screen(
            "p01-s02-c0022",
            [_trajectory(seed) for seed in range(1, 21)],
            window_length=20,
        )
        self.assertEqual(len(screens), 5)
        self.assertTrue(all(row["stationarity_screen_pass"] for row in screens))

    def test_separated_diagnostic_accepts_a_stable_common_level(self) -> None:
        rows = self.module.separated_stability_diagnostic(
            "p01-s02-c0022",
            [_trajectory(seed) for seed in range(1, 13)],
            window_length=20,
            bootstrap_resamples=500,
        )
        by_response = {row["response"]: row for row in rows}
        self.assertEqual(
            by_response["D1"]["final_classification"],
            "stable_with_negligible_heterogeneity",
        )
        self.assertEqual(by_response["D1"]["within_seed_trend_status"], "stable")
        self.assertEqual(
            by_response["D1"]["between_seed_distribution_status"], "stable"
        )

    def test_separated_diagnostic_preserves_stable_seed_heterogeneity(self) -> None:
        rows = self.module.separated_stability_diagnostic(
            "p01-s02-c0022",
            [
                _trajectory(seed, d1_seed_offset=5.0 * seed)
                for seed in range(1, 13)
            ],
            window_length=20,
            bootstrap_resamples=500,
        )
        by_response = {row["response"]: row for row in rows}
        self.assertEqual(
            by_response["D1"]["final_classification"],
            "stable_with_persistent_heterogeneity",
        )
        self.assertEqual(
            by_response["D1"]["between_seed_heterogeneity_status"], "meaningful"
        )
        self.assertGreater(
            by_response["D1"]["rank_normalized_split_rhat_secondary"], 1.05
        )

    def test_separated_diagnostic_detects_continuing_change(self) -> None:
        rows = self.module.separated_stability_diagnostic(
            "p01-s02-c0022",
            [_trajectory(seed, trend_d1=True) for seed in range(1, 13)],
            window_length=20,
            bootstrap_resamples=500,
        )
        by_response = {row["response"]: row for row in rows}
        self.assertEqual(
            by_response["D1"]["final_classification"],
            "continuing_change_detected",
        )
        self.assertEqual(
            by_response["D1"]["within_seed_trend_status"], "continuing_change"
        )

    def test_precision_recommendations_use_minimum_batches_and_cap(self) -> None:
        self.assertEqual(self.module._recommended_replicates(1), (20, False))
        self.assertEqual(self.module._recommended_replicates(21), (28, False))
        self.assertEqual(self.module._recommended_replicates(96), (100, False))
        self.assertEqual(self.module._recommended_replicates(101), (100, True))

    def test_stationarity_history_identifies_persistent_assessment(self) -> None:
        self.assertEqual(
            self.module._assessment_endpoints(53),
            [212, 250],
        )
        history = [
            {
                "assessment_end_generation": generation,
                "stationarity_screen_pass": passed,
            }
            for generation, passed in (
                (80, False),
                (100, True),
                (120, False),
                (140, True),
            )
        ]
        self.assertEqual(self.module._persistent_passing_generation(history), 140)

    def test_completion_gate_lists_every_missing_population(self) -> None:
        rows = self.module.read_tsv(RUNS)
        with tempfile.TemporaryDirectory() as directory:
            issues = self.module.completion_gate_issues(
                rows,
                work=WORK,
                scratch=Path(directory),
            )
        self.assertEqual(len(rows), 120)
        self.assertGreaterEqual(len(issues), 120)
        self.assertTrue(any("completion.json" in issue for issue in issues))

    def test_completion_gate_distinguishes_file_and_resolved_config_hashes(
        self,
    ) -> None:
        rows = self.module.read_tsv(RUNS)
        run = rows[0]
        config = load_config(WORK / run["config_path"])
        resolved = config.to_dict()
        resolved_hash = self.module._resolved_config_sha256(resolved)
        self.assertNotEqual(resolved_hash, run["config_sha256"])

        with tempfile.TemporaryDirectory() as directory:
            scratch = Path(directory)
            output = scratch / run["scratch_relative_path"]
            output.mkdir(parents=True)
            for name in self.module.REQUIRED_RAW_FILES:
                if name not in {"completion.json", "resolved_config.json"}:
                    (output / name).write_bytes(b"")
            (output / "resolved_config.json").write_text(
                json.dumps(resolved), encoding="utf-8"
            )
            final = output / "final_environment_rep000.npz"
            completion = {
                "complete": True,
                "model_spec_version": self.module.MODEL_SPEC_VERSION,
                "output_schema_version": self.module.OUTPUT_SCHEMA_VERSION,
                "software_version": self.module.SOFTWARE_VERSION,
                "config_sha256": resolved_hash,
                "final_environment_sha256": {final.name: self.module.sha256(final)},
                "output_sizes": {},
            }
            (output / "completion.json").write_text(
                json.dumps(completion), encoding="utf-8"
            )
            issues = self.module.completion_gate_issues(
                [run], work=WORK, scratch=scratch
            )
        self.assertEqual(issues, ["manifest contains 1 runs; expected 120"])


if __name__ == "__main__":
    unittest.main()
