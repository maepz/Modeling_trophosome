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


class Phase1PilotArtifactTests(unittest.TestCase):
    def test_frozen_population_and_generated_pilot_files_verify(self) -> None:
        for command in (
            [
                sys.executable,
                "scripts/freeze_fisher_logseries_population.py",
                "--verify",
            ],
            [
                sys.executable,
                "scripts/prepare_phase1_first_pilot.py",
                "--verify",
            ],
        ):
            subprocess.run(command, cwd=REPOSITORY, check=True, capture_output=True)

    def test_registered_cells_match_valid_neutral_configs(self) -> None:
        with (WORK / "registry" / "cells.csv").open() as handle:
            registered_cells = list(csv.DictReader(handle))
        expected_ids = {f"p01-s01-c{number:04d}" for number in range(1, 13)}
        self.assertTrue(
            expected_ids.issubset({cell["cell_id"] for cell in registered_cells})
        )
        cells = [
            cell for cell in registered_cells if cell["cell_id"] in expected_ids
        ]
        population = json.loads(
            (
                WORK
                / "common"
                / "initial-populations"
                / "ip001-fisher100.json"
            ).read_text(encoding="utf-8")
        )
        expected_counts = tuple(population["scaled_counts"])
        for cell in cells:
            config = load_config(WORK / cell["config_path"])
            self.assertEqual(config.environment.initial_counts, expected_counts)
            self.assertEqual(sum(config.environment.initial_counts), 1_000_000_000)
            self.assertFalse(config.evolution.within_host_selection)
            self.assertFalse(config.evolution.free_living_selection)
            self.assertEqual(config.evolution.mutation_effect_mean, 0.0)
            self.assertEqual(config.evolution.mutation_effect_sd, 0.0)
            self.assertEqual(config.output.environment_counts_mode, "final")
            self.assertEqual(config.replicates, 1)

    def test_matched_return_cells_have_identical_return_and_feedback(self) -> None:
        config_directory = (
            WORK / "p01-neutral-feedback" / "configs" / "s01-pilot"
        )
        observed = set()
        for number in range(3, 7):
            config = load_config(
                config_directory / f"p01-s01-c{number:04d}.toml"
            )
            released_per_host = round(
                config.host.carrying_capacity * config.host.escape_fraction
            )
            total_return = config.host.population_size * released_per_host
            capacity = round(
                config.environment.capacity_ratio * config.host.carrying_capacity
            )
            observed.add((total_return, total_return / (capacity + total_return)))
        self.assertEqual(observed, {(1_000_000_000, 0.5)})

    def test_each_seed_block_is_one_independently_runnable_population(self) -> None:
        manifest_directory = WORK / "p01-neutral-feedback" / "manifests"
        with (manifest_directory / "phase1-first-pilot-runs.tsv").open() as handle:
            runs = list(csv.DictReader(handle, delimiter="\t"))
        self.assertEqual(len(runs), 36)
        expected_seeds = {"sb0001": 666, "sb0002": 667, "sb0003": 668}
        for run in runs:
            config = load_config(WORK / run["config_path"])
            self.assertEqual(config.replicates, 1)
            self.assertEqual(config.seed, expected_seeds[run["seed_block_id"]])
            self.assertEqual(run["within_run_replicate_index"], "0")


if __name__ == "__main__":
    unittest.main()
