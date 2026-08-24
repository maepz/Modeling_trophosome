from __future__ import annotations

import csv
import json
import tempfile
import unittest
from concurrent.futures import ThreadPoolExecutor
from dataclasses import replace
from pathlib import Path
from unittest.mock import patch

import numpy as np

from trophosome.config import MigrationConfig, load_config
from trophosome.count_model import free_living_selection_step
from trophosome.simulation import _bounded_host_map, run_simulation

REPOSITORY = Path(__file__).resolve().parents[1]


class SimulationTests(unittest.TestCase):
    def test_parallel_submission_window_does_not_eagerly_consume_all_hosts(
        self,
    ) -> None:
        consumed = 0

        def batches():
            nonlocal consumed
            for value in range(10):
                consumed += 1
                yield (value,)

        with patch(
            "trophosome.simulation._run_host_batch",
            side_effect=lambda batch: batch,
        ):
            with ThreadPoolExecutor(max_workers=1) as executor:
                results = _bounded_host_map(executor, batches(), max_in_flight=2)
                self.assertEqual(next(results), (0,))
                self.assertEqual(consumed, 2)
                self.assertEqual(list(results), [(value,) for value in range(1, 10)])

    def test_smoke_run_writes_compact_reproducible_outputs(self) -> None:
        config = load_config(REPOSITORY / "configs" / "smoke.toml")
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory)
            summaries = run_simulation(config, output, REPOSITORY)
            self.assertEqual(len(summaries), 2)
            self.assertTrue((output / "resolved_config.json").is_file())
            self.assertTrue((output / "provenance.json").is_file())
            self.assertTrue((output / "final_environment_rep000.npz").is_file())
            self.assertTrue((output / "environment_counts.csv").is_file())
            self.assertTrue((output / "infection_counts.csv").is_file())
            self.assertTrue((output / "release_counts.csv").is_file())
            self.assertTrue((output / "migration_counts.csv").is_file())
            self.assertTrue((output / "strain_origins.csv").is_file())
            self.assertTrue((output / "strain_lineage_events.csv").is_file())
            with (output / "host_generation_summary.csv").open() as handle:
                rows = list(csv.DictReader(handle))
            self.assertEqual(len(rows), 2)
            provenance = json.loads((output / "provenance.json").read_text())
            self.assertEqual(provenance["seed"], 666)
            self.assertEqual(provenance["software_version"], "0.7.0")
            self.assertEqual(provenance["model_family"], "wright_fisher_counts")
            self.assertEqual(provenance["model_spec_version"], "2.1.0")
            self.assertEqual(provenance["output_schema_version"], "2.3.0")
            self.assertEqual(len(provenance["config_sha256"]), 64)
            with np.load(output / "final_environment_rep000.npz") as final:
                self.assertEqual(
                    int(final["counts"].sum()), config.host.carrying_capacity
                )
                self.assertIn("within_host_fitness", final.files)
                self.assertIn("free_living_fitness", final.files)
            with (output / "environment_counts.csv").open() as handle:
                header = next(csv.reader(handle))
            self.assertIn("within_host_fitness", header)
            self.assertIn("free_living_fitness", header)
            with (output / "migration_counts.csv").open() as handle:
                self.assertEqual(len(list(csv.DictReader(handle))), 0)
            with (output / "strain_origins.csv").open() as handle:
                origins = list(csv.DictReader(handle))
            self.assertEqual(len(origins), len(config.environment.initial_counts))
            self.assertEqual({row["origin_type"] for row in origins}, {"initial_focal"})

    def test_final_environment_mode_writes_only_each_replicate_endpoint(self) -> None:
        config = load_config(REPOSITORY / "configs" / "smoke.toml")
        config = replace(
            config,
            replicates=2,
            host=replace(config.host, host_generations=3),
            output=replace(config.output, environment_counts_mode="final"),
        )
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory)
            run_simulation(config, output, REPOSITORY)
            with (output / "environment_counts.csv").open() as handle:
                rows = list(csv.DictReader(handle))
            self.assertEqual(
                {(int(row["replicate"]), int(row["generation"])) for row in rows},
                {(0, 3), (1, 3)},
            )
            for replicate in range(2):
                observed = {
                    int(row["strain_id"]): int(row["count"])
                    for row in rows
                    if int(row["replicate"]) == replicate
                }
                with np.load(
                    output / f"final_environment_rep{replicate:03d}.npz"
                ) as final:
                    expected = dict(
                        zip(
                            final["genotype_ids"].tolist(),
                            final["counts"].tolist(),
                            strict=True,
                        )
                    )
                self.assertEqual(observed, expected)

    def test_zero_return_leaves_dormant_reservoir_unchanged(self) -> None:
        config = load_config(REPOSITORY / "configs" / "smoke.toml")
        config = replace(
            config,
            host=replace(
                config.host,
                host_generations=3,
                escape_fraction=0.0,
            ),
        )
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory)
            run_simulation(config, output, REPOSITORY)
            with (output / "environment_counts.csv").open() as handle:
                rows = list(csv.DictReader(handle))
            by_generation: dict[int, list[tuple[str, str]]] = {}
            for row in rows:
                by_generation.setdefault(int(row["generation"]), []).append(
                    (row["strain_id"], row["count"])
                )
            initial = by_generation[0]
            self.assertEqual(sorted(by_generation), [0, 1, 2, 3])
            self.assertTrue(
                all(composition == initial for composition in by_generation.values())
            )

    def test_free_living_selection_acts_without_host_return(self) -> None:
        config = load_config(REPOSITORY / "configs" / "smoke.toml")
        config = replace(
            config,
            environment=replace(
                config.environment,
                initial_counts=(50, 50),
                initial_within_host_fitness=(0.0, 1.0),
                initial_free_living_fitness=(1.0, 0.0),
            ),
            host=replace(
                config.host,
                host_generations=1,
                escape_fraction=0.0,
            ),
            evolution=replace(config.evolution, free_living_selection=True),
        )
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory)
            run_simulation(config, output, REPOSITORY)
            with np.load(output / "final_environment_rep000.npz") as final:
                np.testing.assert_array_equal(final["genotype_ids"], [0])
                np.testing.assert_array_equal(final["counts"], [100])
                np.testing.assert_allclose(final["within_host_fitness"], [0.0])
                np.testing.assert_allclose(final["free_living_fitness"], [1.0])

    def test_fixed_pool_migration_introduces_regional_only_strain(self) -> None:
        config = load_config(REPOSITORY / "configs" / "smoke.toml")
        config = replace(
            config,
            environment=replace(
                config.environment,
                initial_counts=(100, 0),
                initial_within_host_fitness=(1.0, 1.0),
                initial_free_living_fitness=(1.0, 1.0),
            ),
            migration=MigrationConfig(
                mode="fixed_regional_pool",
                fraction=1.0,
                regional_counts=(0, 100),
            ),
            host=replace(config.host, host_generations=1, escape_fraction=0.0),
        )
        selected_inputs: list[np.ndarray] = []

        def observe_selection(population, target_size, rng):
            selected_inputs.append(population.genotype_ids.copy())
            return free_living_selection_step(population, target_size, rng)

        config_with_selection = replace(
            config,
            evolution=replace(config.evolution, free_living_selection=True),
        )
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory)
            with patch(
                "trophosome.simulation.free_living_selection_step",
                side_effect=observe_selection,
            ):
                summaries = run_simulation(config_with_selection, output, REPOSITORY)
            with np.load(output / "final_environment_rep000.npz") as final:
                np.testing.assert_array_equal(final["genotype_ids"], [1])
                np.testing.assert_array_equal(final["counts"], [100])
            self.assertEqual(len(selected_inputs), 1)
            np.testing.assert_array_equal(selected_inputs[0], [1])
            self.assertEqual(summaries[0].migration_replacement_count, 100)
            self.assertEqual(summaries[0].realized_migration_fraction, 1.0)
            with (output / "migration_counts.csv").open() as handle:
                migration_rows = list(csv.DictReader(handle))
            self.assertEqual(
                sum(int(row["emigrant_count"]) for row in migration_rows), 100
            )
            self.assertEqual(
                sum(int(row["immigrant_count"]) for row in migration_rows), 100
            )
            with (output / "strain_origins.csv").open() as handle:
                origins = {int(row["strain_id"]): row for row in csv.DictReader(handle)}
            self.assertEqual(origins[0]["origin_type"], "initial_focal")
            self.assertEqual(origins[1]["origin_type"], "fixed_regional_pool")

    def test_results_do_not_depend_on_worker_count(self) -> None:
        config = load_config(REPOSITORY / "configs" / "smoke.toml")
        config = replace(
            config,
            migration=MigrationConfig(
                mode="fixed_regional_pool",
                fraction=0.2,
                regional_counts=(500, 1500, 3000, 5000),
            ),
        )
        parallel = replace(
            config,
            execution=replace(config.execution, workers=2, host_batch_size=1),
        )
        with tempfile.TemporaryDirectory() as first_directory:
            with tempfile.TemporaryDirectory() as second_directory:
                first = Path(first_directory)
                second = Path(second_directory)
                run_simulation(config, first, REPOSITORY)
                run_simulation(parallel, second, REPOSITORY)
                with np.load(first / "final_environment_rep000.npz") as one:
                    with np.load(second / "final_environment_rep000.npz") as two:
                        np.testing.assert_array_equal(
                            one["genotype_ids"], two["genotype_ids"]
                        )
                        np.testing.assert_array_equal(one["counts"], two["counts"])
                        np.testing.assert_array_equal(
                            one["within_host_fitness"], two["within_host_fitness"]
                        )
                        np.testing.assert_array_equal(
                            one["free_living_fitness"], two["free_living_fitness"]
                        )
                self.assertEqual(
                    (first / "strain_lineage_events.csv").read_text(),
                    (second / "strain_lineage_events.csv").read_text(),
                )
                self.assertEqual(
                    (first / "migration_counts.csv").read_text(),
                    (second / "migration_counts.csv").read_text(),
                )


if __name__ == "__main__":
    unittest.main()
