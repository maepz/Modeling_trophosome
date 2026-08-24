from __future__ import annotations

import tempfile
import unittest
from dataclasses import replace
from pathlib import Path
from unittest.mock import patch

import numpy as np

import trophosome.simulation as simulation
from trophosome.checkpointing import CheckpointMismatchError
from trophosome.config import MigrationConfig, ModelConfig, load_config

REPOSITORY = Path(__file__).resolve().parents[1]


class CheckpointRestartTests(unittest.TestCase):
    def _config(self, *, replicates: int = 1, generations: int = 3) -> ModelConfig:
        config = load_config(REPOSITORY / "configs" / "smoke.toml")
        return replace(
            config,
            replicates=replicates,
            host=replace(config.host, host_generations=generations),
            output=replace(
                config.output,
                checkpoint_interval="0.001ms",
                checkpoint_keep=2,
            ),
        )

    def _run_until_checkpoint(
        self,
        config: ModelConfig,
        output: Path,
        checkpoint_number: int,
        *,
        resume: bool = False,
    ) -> None:
        real_write = simulation.write_recovery_checkpoint
        calls = 0

        def interrupt_after_write(*args, **kwargs):
            nonlocal calls
            path = real_write(*args, **kwargs)
            calls += 1
            if calls == checkpoint_number:
                raise RuntimeError("simulated interruption")
            return path

        with patch(
            "trophosome.simulation.write_recovery_checkpoint",
            side_effect=interrupt_after_write,
        ):
            with self.assertRaisesRegex(RuntimeError, "simulated interruption"):
                simulation.run_simulation(config, output, REPOSITORY, resume=resume)

    def _assert_scientific_outputs_equal(self, first: Path, second: Path) -> None:
        first_csv = sorted(path.name for path in first.glob("*.csv"))
        second_csv = sorted(path.name for path in second.glob("*.csv"))
        self.assertEqual(first_csv, second_csv)
        for name in first_csv:
            self.assertEqual(
                (first / name).read_bytes(),
                (second / name).read_bytes(),
                name,
            )
        first_final = sorted(path.name for path in first.glob("final_environment*.npz"))
        second_final = sorted(
            path.name for path in second.glob("final_environment*.npz")
        )
        self.assertEqual(first_final, second_final)
        for name in first_final:
            with np.load(first / name) as expected, np.load(second / name) as observed:
                self.assertEqual(expected.files, observed.files)
                for field in expected.files:
                    np.testing.assert_array_equal(expected[field], observed[field])

    def test_interrupted_resume_matches_uninterrupted_run_and_truncates_rows(
        self,
    ) -> None:
        config = self._config()
        with tempfile.TemporaryDirectory() as baseline_directory:
            with tempfile.TemporaryDirectory() as resumed_directory:
                baseline = Path(baseline_directory)
                resumed = Path(resumed_directory)
                simulation.run_simulation(config, baseline, REPOSITORY)
                self._run_until_checkpoint(config, resumed, checkpoint_number=1)
                with (resumed / "environment_counts.csv").open("a") as handle:
                    handle.write("partial,row,after,checkpoint\n")

                parallel = replace(
                    config,
                    execution=replace(
                        config.execution,
                        workers=2,
                        host_batch_size=1,
                    ),
                )
                summaries = simulation.run_simulation(
                    parallel, resumed, REPOSITORY, resume=True
                )

                self.assertEqual(len(summaries), config.host.host_generations)
                self._assert_scientific_outputs_equal(baseline, resumed)
                self.assertFalse((resumed / "checkpoints").exists())
                self.assertTrue((resumed / "completion.json").is_file())

    def test_migration_enabled_resume_matches_uninterrupted_run(self) -> None:
        config = replace(
            self._config(generations=3),
            migration=MigrationConfig(
                mode="fixed_regional_pool",
                fraction=0.2,
                regional_counts=(500, 1500, 3000, 5000),
            ),
        )
        with tempfile.TemporaryDirectory() as baseline_directory:
            with tempfile.TemporaryDirectory() as resumed_directory:
                baseline = Path(baseline_directory)
                resumed = Path(resumed_directory)
                simulation.run_simulation(config, baseline, REPOSITORY)
                self._run_until_checkpoint(config, resumed, checkpoint_number=1)
                simulation.run_simulation(config, resumed, REPOSITORY, resume=True)
                self._assert_scientific_outputs_equal(baseline, resumed)

    def test_corrupt_latest_checkpoint_falls_back_and_only_two_are_retained(
        self,
    ) -> None:
        config = self._config(generations=4)
        with tempfile.TemporaryDirectory() as baseline_directory:
            with tempfile.TemporaryDirectory() as resumed_directory:
                baseline = Path(baseline_directory)
                resumed = Path(resumed_directory)
                simulation.run_simulation(config, baseline, REPOSITORY)
                self._run_until_checkpoint(config, resumed, checkpoint_number=3)
                checkpoints = sorted((resumed / "checkpoints").glob("*.npz"))
                self.assertEqual(len(checkpoints), 2)
                self.assertIn("gen000002", checkpoints[0].name)
                self.assertIn("gen000003", checkpoints[1].name)
                checkpoints[-1].write_bytes(b"corrupted checkpoint")

                simulation.run_simulation(config, resumed, REPOSITORY, resume=True)

                self._assert_scientific_outputs_equal(baseline, resumed)
                self.assertFalse((resumed / "checkpoints").exists())

    def test_final_environment_mode_resumes_to_endpoint_only_output(self) -> None:
        config = self._config(generations=3)
        config = replace(
            config,
            output=replace(config.output, environment_counts_mode="final"),
        )
        with tempfile.TemporaryDirectory() as baseline_directory:
            with tempfile.TemporaryDirectory() as resumed_directory:
                baseline = Path(baseline_directory)
                resumed = Path(resumed_directory)
                simulation.run_simulation(config, baseline, REPOSITORY)
                self._run_until_checkpoint(config, resumed, checkpoint_number=1)
                simulation.run_simulation(config, resumed, REPOSITORY, resume=True)
                self._assert_scientific_outputs_equal(baseline, resumed)
                with (resumed / "environment_counts.csv").open() as handle:
                    generations = {
                        int(line.split(",", 2)[1]) for line in handle.readlines()[1:]
                    }
                self.assertEqual(generations, {3})

    def test_scientific_configuration_mismatch_is_rejected(self) -> None:
        config = self._config()
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory)
            self._run_until_checkpoint(config, output, checkpoint_number=1)
            changed = replace(
                config,
                evolution=replace(
                    config.evolution,
                    mutation_probability=config.evolution.mutation_probability * 2,
                ),
            )
            with self.assertRaisesRegex(
                CheckpointMismatchError, "restart_config_sha256"
            ):
                simulation.run_simulation(changed, output, REPOSITORY, resume=True)

    def test_source_mismatch_is_rejected(self) -> None:
        config = self._config()
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory)
            self._run_until_checkpoint(config, output, checkpoint_number=1)
            with patch("trophosome.simulation._source_hash", return_value="0" * 64):
                with self.assertRaisesRegex(CheckpointMismatchError, "source_sha256"):
                    simulation.run_simulation(config, output, REPOSITORY, resume=True)
            with patch("trophosome.simulation.OUTPUT_SCHEMA_VERSION", "999.0.0"):
                with self.assertRaisesRegex(
                    CheckpointMismatchError, "output_schema_version"
                ):
                    simulation.run_simulation(config, output, REPOSITORY, resume=True)

    def test_repeated_resume_cycles_preserve_zero_return_invariant(self) -> None:
        config = self._config(generations=4)
        config = replace(
            config,
            host=replace(config.host, escape_fraction=0.0),
        )
        with tempfile.TemporaryDirectory() as baseline_directory:
            with tempfile.TemporaryDirectory() as resumed_directory:
                baseline = Path(baseline_directory)
                resumed = Path(resumed_directory)
                simulation.run_simulation(config, baseline, REPOSITORY)
                self._run_until_checkpoint(config, resumed, checkpoint_number=1)
                self._run_until_checkpoint(
                    config, resumed, checkpoint_number=1, resume=True
                )
                simulation.run_simulation(config, resumed, REPOSITORY, resume=True)
                self._assert_scientific_outputs_equal(baseline, resumed)

                with (resumed / "environment_counts.csv").open() as handle:
                    lines = handle.readlines()
                initial = [line for line in lines[1:] if line.startswith("0,0,")]
                for generation in range(1, config.host.host_generations + 1):
                    observed = [
                        line.replace(f"0,{generation},", "0,0,", 1)
                        for line in lines[1:]
                        if line.startswith(f"0,{generation},")
                    ]
                    self.assertEqual(observed, initial)

    def test_replicate_boundary_checkpoint_resumes_at_next_replicate(self) -> None:
        config = self._config(replicates=2, generations=2)
        config = replace(
            config,
            output=replace(config.output, checkpoint_interval="1h"),
        )
        with tempfile.TemporaryDirectory() as baseline_directory:
            with tempfile.TemporaryDirectory() as resumed_directory:
                baseline = Path(baseline_directory)
                resumed = Path(resumed_directory)
                simulation.run_simulation(config, baseline, REPOSITORY)
                self._run_until_checkpoint(config, resumed, checkpoint_number=1)
                simulation.run_simulation(config, resumed, REPOSITORY, resume=True)
                self._assert_scientific_outputs_equal(baseline, resumed)

    def test_resume_of_verified_completed_run_is_idempotent(self) -> None:
        config = self._config(generations=2)
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory)
            expected = simulation.run_simulation(config, output, REPOSITORY)
            observed = simulation.run_simulation(
                config, output, REPOSITORY, resume=True
            )
            self.assertEqual(expected, observed)
            self.assertFalse((output / "checkpoints").exists())
            with self.assertRaises(FileExistsError):
                simulation.run_simulation(config, output, REPOSITORY)


if __name__ == "__main__":
    unittest.main()
