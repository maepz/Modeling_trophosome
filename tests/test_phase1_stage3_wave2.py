from __future__ import annotations

import csv
import json
import shutil
import subprocess
import sys
import tempfile
import unittest
from dataclasses import replace
from pathlib import Path
from unittest.mock import patch

from trophosome.config import MigrationConfig, load_config
from trophosome.simulation import run_simulation

REPOSITORY = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPOSITORY / "scripts"))
import assess_phase1_stage3_wave2_horizon as assessment  # noqa: E402
import prepare_phase1_stage3_wave2 as design  # noqa: E402
import run_phase1_stage3_wave2 as runner  # noqa: E402
from run_phase1_first_pilot import _sha256  # noqa: E402


def portable_repository(root: Path) -> tuple[Path, Path]:
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
    layout.write_text(json.dumps({"scratch": str(scratch)}))
    return root / "experiments/work/trophosome", scratch


class Wave2DesignTests(unittest.TestCase):
    def test_frozen_files_verify(self) -> None:
        result = subprocess.run(
            [sys.executable, "scripts/prepare_phase1_stage3_wave2.py", "--verify"],
            cwd=REPOSITORY,
            capture_output=True,
            text=True,
            check=True,
        )
        self.assertIn("Verified 453 Wave 2 files", result.stdout)

    def test_exact_panels_reuse_and_matched_founder_totals(self) -> None:
        self.assertEqual(len(design.CELLS), 40)
        self.assertEqual(len(design.NEW_CELLS), 34)
        self.assertEqual(len(design.REUSED_CELLS), 6)
        h_by_b = [cell for cell in design.CELLS if cell.panel == "H-by-B"]
        alpha_by_m = [cell for cell in design.CELLS if cell.panel == "alpha-by-m"]
        self.assertEqual(len(h_by_b), 12)
        self.assertEqual(len(alpha_by_m), 28)
        self.assertEqual(
            {(cell.hosts, cell.bottleneck) for cell in h_by_b},
            {(h, b) for h in (100, 1000, 10000) for b in (1, 5, 10, 50)},
        )
        self.assertEqual(
            {(cell.alpha_target, cell.migration_fraction) for cell in alpha_by_m},
            {
                (alpha, migration)
                for alpha in ("0", "0.01", "0.1", "0.99")
                for migration in ("0", "0.001", "0.01", "0.1", "0.5", "0.9", "0.99")
            },
        )
        by_hb: dict[int, list[design.Wave2Cell]] = {}
        for cell in h_by_b:
            by_hb.setdefault(cell.hosts * cell.bottleneck, []).append(cell)
            self.assertEqual(cell.total_return, design.K)
            self.assertAlmostEqual(cell.alpha, 0.5)
        self.assertEqual(
            {value for value, cells in by_hb.items() if len(cells) == 2},
            {1000, 5000, 10000, 50000},
        )
        self.assertEqual(
            {cell.reused_source for cell in design.REUSED_CELLS},
            {
                "p01-s02-c0021",
                "p01-s02-c0022",
                "p01-s02-c0024",
                "p01-s03-c0037",
                "p01-s03-c0039",
                "p01-s03-c0041",
            },
        )
        positive_feedback = {
            cell.alpha_target: cell.total_return
            for cell in alpha_by_m
            if cell.migration_fraction == "0.1"
        }
        self.assertEqual(
            positive_feedback,
            {
                "0": 0,
                "0.01": 10_100_000,
                "0.1": 111_110_000,
                "0.99": 99_000_000_000,
            },
        )

    def test_configs_freeze_neutral_dynamics_and_maximum_horizon(self) -> None:
        work = REPOSITORY / "experiments/work/trophosome"
        manifest = (
            work
            / "p01-neutral-feedback/manifests"
            / f"{design.EXPERIMENT_ID}-runs.tsv"
        )
        with manifest.open(newline="", encoding="utf-8") as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        self.assertEqual(len(rows), 408)
        self.assertEqual(len({row["run_id"] for row in rows}), 408)
        for row in rows:
            config = load_config(work / row["config_path"])
            self.assertEqual(config.host.host_generations, 1000)
            self.assertEqual(
                config.seed, dict(design.SEED_BLOCKS)[row["seed_block_id"]]
            )
            self.assertEqual(config.evolution.mutation_probability, 0)
            self.assertFalse(config.evolution.within_host_selection)
            self.assertFalse(config.evolution.free_living_selection)
            self.assertEqual(
                config.migration.regional_counts, config.environment.initial_counts
            )
            self.assertEqual(config.output.environment_counts_mode, "all")

    def test_reused_trajectories_are_complete_through_passage_100(self) -> None:
        path = (
            REPOSITORY
            / "experiments/work/trophosome/p01-neutral-feedback/design"
            / f"{design.EXPERIMENT_ID}-reused-trajectories.tsv"
        )
        with path.open(newline="", encoding="utf-8") as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        self.assertEqual(len(rows), 6 * 12 * 101)
        self.assertEqual({int(row["generation"]) for row in rows}, set(range(101)))
        self.assertEqual(len({row["cell_id"] for row in rows}), 6)
        self.assertEqual(len({row["seed_block_id"] for row in rows}), 12)

    def test_no_return_environment_is_independent_of_host_number(self) -> None:
        base = load_config(REPOSITORY / "configs/smoke.toml")
        base = replace(
            base,
            migration=MigrationConfig(
                mode="fixed_regional_pool",
                fraction=0.1,
                regional_counts=base.environment.initial_counts,
            ),
            host=replace(
                base.host,
                population_size=2,
                host_generations=3,
                escape_fraction=0,
            ),
            evolution=replace(base.evolution, mutation_probability=0),
        )
        many = replace(base, host=replace(base.host, population_size=5))
        with tempfile.TemporaryDirectory() as first_directory:
            with tempfile.TemporaryDirectory() as second_directory:
                first = Path(first_directory)
                second = Path(second_directory)
                run_simulation(base, first, REPOSITORY)
                run_simulation(many, second, REPOSITORY)
                self.assertEqual(
                    (first / "environment_counts.csv").read_bytes(),
                    (second / "environment_counts.csv").read_bytes(),
                )
                self.assertEqual(
                    (first / "migration_counts.csv").read_bytes(),
                    (second / "migration_counts.csv").read_bytes(),
                )

    def test_dry_run_is_non_mutating(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            _work, scratch = portable_repository(root)
            result = subprocess.run(
                [
                    sys.executable,
                    str(REPOSITORY / "scripts/run_phase1_stage3_wave2.py"),
                    "--repository",
                    str(root),
                    "--dry-run",
                ],
                capture_output=True,
                text=True,
                check=True,
            )
            self.assertIn("408 populations toward passage 100", result.stdout)
            self.assertEqual(list(scratch.iterdir()), [])

    def test_wave2_audit_accepts_a_verified_planned_pause(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            work = root / "work"
            scratch = root / "scratch"
            config_path = work / "configs/smoke.toml"
            config_path.parent.mkdir(parents=True)
            scratch.mkdir()
            shutil.copyfile(REPOSITORY / "configs/smoke.toml", config_path)
            output = scratch / "run"
            config = load_config(config_path)
            run_simulation(
                config,
                output,
                REPOSITORY,
                pause_after_generation=1,
            )
            run_id = "synthetic-wave2-audit"
            (output / "execution-summary.json").write_text(
                json.dumps({"run_id": run_id, "status": "paused"})
            )
            row = {
                "run_id": run_id,
                "config_path": "configs/smoke.toml",
                "config_sha256": _sha256(config_path),
                "scratch_relative_path": "run",
            }
            self.assertEqual(
                runner.state_issues([row], work=work, scratch=scratch, horizon=1),
                [],
            )


class AdaptiveDecisionTests(unittest.TestCase):
    def _stable_values(self) -> dict[tuple[str, str, int], float]:
        return {
            (cell.cell_id, seed, generation): 0.01
            for cell in design.CELLS
            for seed, _master in design.SEED_BLOCKS
            for generation in range(101)
        }

    def test_passage100_keeps_only_prespecified_anchors_when_stable(self) -> None:
        values = self._stable_values()
        with patch.object(assessment, "load_trajectories", return_value=(values, [])):
            decision, diagnostics = assessment.build_decision(REPOSITORY, 100)
        expected = {
            cell.cell_id
            for cell in design.CELLS
            if cell.panel == "alpha-by-m"
            and cell.alpha_target in {"0", "0.1"}
            and cell.migration_fraction in {"0", "0.001", "0.01"}
        }
        self.assertEqual(set(decision["selected_cell_ids"]), expected)
        self.assertEqual(decision["continuation_horizon"], 500)
        self.assertEqual(decision["selected_populations"], 72)
        self.assertEqual(len(diagnostics), 21)

    def test_unresolved_treatment_adds_it_and_its_same_m_control(self) -> None:
        values = self._stable_values()
        treatment = next(
            cell
            for cell in design.CELLS
            if cell.panel == "alpha-by-m"
            and cell.alpha_target == "0.99"
            and cell.migration_fraction == "0.01"
        )
        control = next(
            cell
            for cell in design.CELLS
            if cell.panel == "alpha-by-m"
            and cell.alpha_target == "0"
            and cell.migration_fraction == "0.01"
        )
        for seed, _master in design.SEED_BLOCKS:
            for generation in range(76, 101):
                values[(treatment.cell_id, seed, generation)] = 0.04
        with patch.object(assessment, "load_trajectories", return_value=(values, [])):
            decision, _diagnostics = assessment.build_decision(REPOSITORY, 100)
        self.assertIn(treatment.cell_id, decision["selected_cell_ids"])
        self.assertIn(control.cell_id, decision["selected_cell_ids"])
        self.assertIn("unresolved-raw-TV", decision["reasons"][treatment.cell_id])


if __name__ == "__main__":
    unittest.main()
