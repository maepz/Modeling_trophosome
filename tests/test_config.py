from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

from trophosome.config import ConfigurationError, load_config

REPOSITORY = Path(__file__).resolve().parents[1]


class ConfigTests(unittest.TestCase):
    def test_smoke_config_is_valid(self) -> None:
        config = load_config(REPOSITORY / "configs" / "smoke.toml")
        self.assertEqual(config.host.infection_bottleneck, 5)
        self.assertEqual(config.environment.sampling_mode, "reservoir")
        self.assertEqual(config.model, "wright_fisher_counts")
        self.assertEqual(config.environment.capacity_ratio, 1.0)
        self.assertEqual(
            len(config.environment.initial_free_living_fitness),
            len(config.environment.initial_counts),
        )
        self.assertFalse(config.evolution.free_living_selection)
        self.assertEqual(config.output.checkpoint_interval, "1h")
        self.assertEqual(config.output.checkpoint_keep, 2)
        self.assertEqual(config.output.environment_counts_mode, "all")

    def test_checkpoint_interval_requires_a_positive_duration(self) -> None:
        text = (REPOSITORY / "configs" / "smoke.toml").read_text()
        for invalid in ('checkpoint_interval = 1', 'checkpoint_interval = "0h"'):
            with self.subTest(invalid=invalid):
                changed = text.replace('checkpoint_interval = "1h"', invalid)
                with tempfile.TemporaryDirectory() as directory:
                    path = Path(directory) / "invalid.toml"
                    path.write_text(changed)
                    with self.assertRaisesRegex(
                        ConfigurationError, "checkpoint_interval"
                    ):
                        load_config(path)

    def test_at_least_two_checkpoints_are_required(self) -> None:
        text = (REPOSITORY / "configs" / "smoke.toml").read_text()
        text = text.replace("checkpoint_keep = 2", "checkpoint_keep = 1")
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "invalid.toml"
            path.write_text(text)
            with self.assertRaisesRegex(ConfigurationError, "checkpoint_keep"):
                load_config(path)

    def test_environment_counts_mode_is_validated(self) -> None:
        text = (REPOSITORY / "configs" / "smoke.toml").read_text()
        text = text.replace(
            'host_counts_mode = "summary"',
            'environment_counts_mode = "invalid"\nhost_counts_mode = "summary"',
        )
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "invalid.toml"
            path.write_text(text)
            with self.assertRaisesRegex(
                ConfigurationError, "environment_counts_mode"
            ):
                load_config(path)

    def test_legacy_single_fitness_config_defaults_environment_to_neutral(self) -> None:
        text = (REPOSITORY / "configs" / "smoke.toml").read_text()
        text = text.replace("initial_within_host_fitness", "initial_fitness").replace(
            "initial_free_living_fitness = [0.98, 1.01, 0.99, 1.0]\n", ""
        )
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "legacy.toml"
            path.write_text(text)
            config = load_config(path)
        self.assertEqual(
            config.environment.initial_free_living_fitness,
            (1.0, 1.0, 1.0, 1.0),
        )

    def test_free_living_fitness_length_is_validated(self) -> None:
        text = (REPOSITORY / "configs" / "smoke.toml").read_text()
        text = text.replace(
            "initial_free_living_fitness = [0.98, 1.01, 0.99, 1.0]",
            "initial_free_living_fitness = [1.0]",
        )
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "invalid.toml"
            path.write_text(text)
            with self.assertRaisesRegex(
                ConfigurationError, "initial_free_living_fitness"
            ):
                load_config(path)

    def test_escape_percentage_is_rejected(self) -> None:
        text = (REPOSITORY / "configs" / "smoke.toml").read_text()
        text = text.replace("escape_fraction = 0.01", "escape_fraction = 100")
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "invalid.toml"
            path.write_text(text)
            with self.assertRaisesRegex(ConfigurationError, "escape_fraction"):
                load_config(path)

    def test_growth_factor_one_cannot_reach_larger_capacity(self) -> None:
        text = (REPOSITORY / "configs" / "smoke.toml").read_text()
        text = text.replace("growth_factor = 1.5", "growth_factor = 1.0")
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "invalid.toml"
            path.write_text(text)
            with self.assertRaisesRegex(ConfigurationError, "cannot grow"):
                load_config(path)

    def test_biological_example_warns_about_transition_count(self) -> None:
        config = load_config(REPOSITORY / "configs" / "biological-scale.example.toml")
        self.assertTrue(
            any(
                "population transitions" in item
                for item in config.feasibility_warnings()
            )
        )

    def test_finite_sampling_checks_effective_capacity(self) -> None:
        text = (REPOSITORY / "configs" / "smoke.toml").read_text()
        text = text.replace('sampling_mode = "reservoir"', 'sampling_mode = "finite"')
        text = text.replace("capacity_ratio = 1.0", "capacity_ratio = 0.1")
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "invalid.toml"
            path.write_text(text)
            with self.assertRaisesRegex(ConfigurationError, "effective capacity"):
                load_config(path)


if __name__ == "__main__":
    unittest.main()
