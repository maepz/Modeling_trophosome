from __future__ import annotations

import csv
import json
import shutil
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from trophosome.config import load_config

REPOSITORY = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPOSITORY / "scripts"))
import prepare_phase1_stage3_wave3 as design  # noqa: E402
import run_phase1_stage3_wave3 as runner  # noqa: E402


def portable_repository(root: Path) -> Path:
    initial = Path(
        "experiments/work/trophosome/common/initial-populations/ip001-fisher100.json"
    )
    destination = root / initial
    destination.parent.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(REPOSITORY / initial, destination)
    for source in design.build_files(REPOSITORY):
        target = root / source.relative_to(REPOSITORY)
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(source, target)
    scratch = root / "scratch"
    scratch.mkdir()
    layout = root / "experiments/work/trophosome/layout.local.json"
    layout.write_text(json.dumps({"scratch": str(scratch)}), encoding="utf-8")
    return scratch


class Wave3DesignTests(unittest.TestCase):
    def test_frozen_files_verify(self) -> None:
        result = subprocess.run(
            [sys.executable, "scripts/prepare_phase1_stage3_wave3.py", "--verify"],
            cwd=REPOSITORY,
            capture_output=True,
            text=True,
            check=True,
        )
        self.assertIn("Verified 116 Wave 3 files", result.stdout)

    def test_exact_bridge_cells_and_extension_are_frozen(self) -> None:
        observed = {
            (
                cell.bridge_code,
                cell.hosts,
                cell.alpha_target,
                cell.bottleneck,
                cell.migration_fraction,
                cell.extension,
            )
            for cell in design.CELLS
        }
        expected = {
            ("A1", 1000, "0.01", 1, "0.1", False),
            ("A2", 1000, "0.01", 50, "0.1", False),
            ("A3", 1000, "0.1", 1, "0.1", False),
            ("A4", 1000, "0.1", 50, "0.1", False),
            ("A5", 1000, "0.99", 1, "0.1", False),
            ("A6", 1000, "0.99", 50, "0.1", False),
            ("B1", 100, "0.1", 10, "0.01", False),
            ("B2", 100, "0.1", 10, "0.9", False),
            ("B3", 10000, "0.1", 10, "0.01", False),
            ("B4", 10000, "0.1", 10, "0.9", False),
            ("C1", 1000, "0.1", 1, "0.01", False),
            ("C2", 1000, "0.1", 1, "0.9", False),
            ("C3", 1000, "0.1", 50, "0.01", False),
            ("C4", 1000, "0.1", 50, "0.9", False),
            ("D1", 10000, "0.01", 1, "0.1", True),
            ("D2", 10000, "0.99", 1, "0.1", True),
        }
        self.assertEqual(observed, expected)
        self.assertEqual(len(design.CELLS), 16)
        self.assertEqual(sum(not cell.extension for cell in design.CELLS), 14)
        self.assertEqual(sum(cell.extension for cell in design.CELLS), 2)
        self.assertEqual(
            design.SEED_BLOCKS,
            tuple((f"sb{number:04d}", 665 + number) for number in range(1, 7)),
        )

    def test_all_96_configs_are_neutral_passage_100_fixed_pool_runs(self) -> None:
        work = REPOSITORY / "experiments/work/trophosome"
        manifest = (
            work / "p01-neutral-feedback/manifests" / f"{design.EXPERIMENT_ID}-runs.tsv"
        )
        with manifest.open(newline="", encoding="utf-8") as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        self.assertEqual(len(rows), 96)
        self.assertEqual(len({row["run_id"] for row in rows}), 96)
        self.assertEqual(
            {row["seed_block_id"] for row in rows},
            {f"sb{number:04d}" for number in range(1, 7)},
        )
        for row in rows:
            config = load_config(work / row["config_path"])
            self.assertEqual(config.host.host_generations, 100)
            self.assertLessEqual(config.host.population_size, 10_000)
            self.assertEqual(config.environment.capacity_ratio, 1)
            self.assertEqual(config.evolution.mutation_probability, 0)
            self.assertFalse(config.evolution.within_host_selection)
            self.assertFalse(config.evolution.free_living_selection)
            self.assertEqual(config.migration.mode, "fixed_regional_pool")
            self.assertEqual(
                config.migration.regional_counts,
                config.environment.initial_counts,
            )
            self.assertEqual(config.output.environment_counts_mode, "all")

    def test_dry_run_is_non_mutating(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            scratch = portable_repository(root)
            result = subprocess.run(
                [
                    sys.executable,
                    str(REPOSITORY / "scripts/run_phase1_stage3_wave3.py"),
                    "--repository",
                    str(root),
                    "--dry-run",
                ],
                capture_output=True,
                text=True,
                check=True,
            )
            self.assertIn("96 populations; 100 passages", result.stdout)
            self.assertEqual(list(scratch.iterdir()), [])

    def test_community_only_requests_tables_through_wave3(self) -> None:
        completed = subprocess.CompletedProcess([], 0)
        with (
            patch.object(
                sys, "argv", ["run_phase1_stage3_wave3.py", "--community-only"]
            ),
            patch.object(runner.subprocess, "run", return_value=completed) as launched,
        ):
            self.assertEqual(runner.main(), 0)
        command = launched.call_args.args[0]
        self.assertIn("compile_phase1_stage3_dbrda_inputs.py", command[1])
        self.assertEqual(command[command.index("--through-wave") + 1], "3")


if __name__ == "__main__":
    unittest.main()
