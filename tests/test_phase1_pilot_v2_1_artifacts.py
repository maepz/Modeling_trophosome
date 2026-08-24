from __future__ import annotations

import csv
import json
import subprocess
import sys
import unittest
from pathlib import Path

from trophosome.config import load_config

REPOSITORY = Path(__file__).resolve().parents[1]
WORK = REPOSITORY / "experiments" / "work" / "trophosome"
PHASE = WORK / "p01-neutral-feedback"
VARIANT = "v210-m010"


class Phase1PilotV21ArtifactTests(unittest.TestCase):
    def test_generated_variant_verifies_separately(self) -> None:
        subprocess.run(
            [
                sys.executable,
                "scripts/prepare_phase1_first_pilot_v2_1.py",
                "--verify",
            ],
            cwd=REPOSITORY,
            check=True,
            capture_output=True,
        )
        self.assertTrue((PHASE / "manifests" / "phase1-first-pilot-runs.tsv").is_file())
        self.assertTrue(
            (PHASE / "manifests" / f"phase1-first-pilot-{VARIANT}-runs.tsv").is_file()
        )

    def test_all_60_runs_have_fixed_pool_migration_and_neutral_selection(self) -> None:
        population = json.loads(
            (
                WORK / "common" / "initial-populations" / "ip001-fisher100.json"
            ).read_text(encoding="utf-8")
        )
        expected_counts = tuple(population["scaled_counts"])
        runs_path = PHASE / "manifests" / f"phase1-first-pilot-{VARIANT}-runs.tsv"
        with runs_path.open(encoding="utf-8") as handle:
            runs = list(csv.DictReader(handle, delimiter="\t"))

        self.assertEqual(len(runs), 60)
        self.assertEqual(len({row["run_id"] for row in runs}), 60)
        self.assertEqual(len({row["scratch_relative_path"] for row in runs}), 60)
        expected_seeds = {"sb0001": 666, "sb0002": 667, "sb0003": 668}
        for row in runs:
            self.assertIn(VARIANT, row["run_id"])
            self.assertIn(f"s01-pilot-{VARIANT}", row["scratch_relative_path"])
            config = load_config(WORK / row["config_path"])
            self.assertEqual(config.seed, expected_seeds[row["seed_block_id"]])
            self.assertEqual(config.replicates, 1)
            self.assertEqual(config.environment.initial_counts, expected_counts)
            self.assertEqual(config.migration.mode, "fixed_regional_pool")
            self.assertEqual(config.migration.fraction, 0.1)
            self.assertEqual(config.migration.regional_counts, expected_counts)
            self.assertFalse(config.evolution.within_host_selection)
            self.assertFalse(config.evolution.free_living_selection)

    def test_design_and_manifest_freeze_the_confirmed_variant(self) -> None:
        design_path = PHASE / "design" / f"phase1-first-pilot-{VARIANT}-cells.tsv"
        with design_path.open(encoding="utf-8") as handle:
            design = list(csv.DictReader(handle, delimiter="\t"))
        self.assertEqual(len(design), 20)
        self.assertEqual(
            {row["cell_id"] for row in design},
            {f"p01-s01-c{number:04d}" for number in range(1, 21)},
        )
        self.assertEqual({row["m"] for row in design}, {"0.1"})
        self.assertEqual(
            {row["migration_replacement_count"] for row in design},
            {"100000000"},
        )

        manifest_path = (
            PHASE / "manifests" / f"phase1-first-pilot-{VARIANT}-manifest.json"
        )
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        self.assertEqual(manifest["model_spec_version"], "2.1.0")
        self.assertEqual(manifest["software_version"], "0.7.0")
        self.assertEqual(manifest["output_schema_version"], "2.3.0")
        self.assertEqual(manifest["regional_pool"]["fraction_m"], 0.1)
        self.assertFalse(manifest["regional_pool"]["depleted_by_immigration"])
        self.assertEqual(len(manifest["runs"]), 60)

    def test_variant_paths_do_not_overlap_historical_outputs(self) -> None:
        manifest_directory = PHASE / "manifests"
        with (manifest_directory / "phase1-first-pilot-runs.tsv").open() as handle:
            historical = list(csv.DictReader(handle, delimiter="\t"))
        with (
            manifest_directory / f"phase1-first-pilot-{VARIANT}-runs.tsv"
        ).open() as handle:
            variant = list(csv.DictReader(handle, delimiter="\t"))
        self.assertTrue(
            {row["run_id"] for row in historical}.isdisjoint(
                row["run_id"] for row in variant
            )
        )
        self.assertTrue(
            {row["scratch_relative_path"] for row in historical}.isdisjoint(
                row["scratch_relative_path"] for row in variant
            )
        )


if __name__ == "__main__":
    unittest.main()
