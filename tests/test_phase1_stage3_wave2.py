from __future__ import annotations

import csv
import json
import math
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
import analyse_phase1_stage3_wave2 as wave2_analysis  # noqa: E402
import assess_phase1_stage3_wave2_horizon as assessment  # noqa: E402
import build_phase1_stage3_wave2_report as reporting  # noqa: E402
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

    def test_summarize_only_dispatches_the_wave2_analysis_without_scratch_setup(
        self,
    ) -> None:
        completed = subprocess.CompletedProcess([], 0)
        with (
            patch.object(
                sys, "argv", ["run_phase1_stage3_wave2.py", "--summarize-only"]
            ),
            patch.object(runner.subprocess, "run", return_value=completed) as launched,
        ):
            self.assertEqual(runner.main(), 0)
        command = launched.call_args.args[0]
        self.assertIn("analyse_phase1_stage3_wave2.py", command[1])
        self.assertIn("--repository", command)

    def test_dbrda_only_dispatches_the_matrix_compiler_without_scratch_setup(
        self,
    ) -> None:
        completed = subprocess.CompletedProcess([], 0)
        with (
            patch.object(sys, "argv", ["run_phase1_stage3_wave2.py", "--dbrda-only"]),
            patch.object(runner.subprocess, "run", return_value=completed) as launched,
        ):
            self.assertEqual(runner.main(), 0)
        command = launched.call_args.args[0]
        self.assertIn("compile_phase1_stage3_dbrda_inputs.py", command[1])
        self.assertIn("--repository", command)


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


class Wave2ReportTests(unittest.TestCase):
    def test_passage100_report_is_self_contained_and_audited(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            outputs = reporting.build(
                REPOSITORY,
                pdf=root / "report.pdf",
                markdown=root / "report.md",
                assets=root / "figures",
                completion=root / "report-completion.json",
            )

            self.assertEqual(len(outputs), 5)
            self.assertTrue((root / "report.pdf").read_bytes().startswith(b"%PDF"))
            markdown = (root / "report.md").read_text(encoding="utf-8")
            self.assertIn("adaptive-decision report", markdown)
            self.assertIn("not mean that every trajectory was flat", markdown)
            self.assertIn("cannot answer the primary H-by-B comparison", markdown)
            completion = json.loads(
                (root / "report-completion.json").read_text(encoding="utf-8")
            )
            self.assertTrue(completion["complete"])
            self.assertEqual(
                completion["scope"], "passage-100-adaptive-horizon-decision"
            )
            self.assertEqual(len(completion["inputs"]), 3)
            self.assertEqual(len(completion["outputs"]), 4)
            for figure in ("late-window-tv.png", "stability-diagnostics.png"):
                self.assertGreater((root / "figures" / figure).stat().st_size, 10_000)


class Wave2AnalysisTests(unittest.TestCase):
    def test_frozen_reuse_supplies_all_72_passage100_populations(self) -> None:
        phase = REPOSITORY / "experiments/work/trophosome/p01-neutral-feedback"
        matrix = wave2_analysis._read_tsv(
            phase
            / "design/phase1-stage3-wave2-v210-adaptive-g1000-cells.tsv"
        )
        cells = {row["cell_id"]: row for row in matrix}
        trajectories, inputs = wave2_analysis._reused_rows(phase, cells)

        self.assertEqual(len(trajectories), 6 * 12 * 101)
        self.assertEqual(len(inputs), 6 * 12)
        self.assertEqual(
            len(
                {
                    (row["cell_id"], row["seed_block_id"], row["generation"])
                    for row in trajectories
                }
            ),
            len(trajectories),
        )
        self.assertTrue(all(row["source_run_id"] for row in trajectories))

    def test_environment_and_host_prefixes_are_summarized_without_later_rows(
        self,
    ) -> None:
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory)
            environment = output / "environment_counts.csv"
            with environment.open("w", newline="", encoding="utf-8") as handle:
                writer = csv.DictWriter(
                    handle,
                    fieldnames=("replicate", "generation", "strain_id", "count"),
                )
                writer.writeheader()
                for generation in range(102):
                    writer.writerow(
                        {
                            "replicate": 0,
                            "generation": generation,
                            "strain_id": 0,
                            "count": 500_000_000 + generation,
                        }
                    )
                    writer.writerow(
                        {
                            "replicate": 0,
                            "generation": generation,
                            "strain_id": 1,
                            "count": 500_000_000 - generation,
                        }
                    )
            summary = output / "host_generation_summary.csv"
            with summary.open("w", newline="", encoding="utf-8") as handle:
                writer = csv.DictWriter(
                    handle,
                    fieldnames=(
                        "replicate",
                        "host_generation",
                        "realized_host_feedback",
                        "mean_adult_richness",
                        "mean_adult_gene_diversity",
                    ),
                )
                writer.writeheader()
                for generation in range(1, 102):
                    writer.writerow(
                        {
                            "replicate": 0,
                            "host_generation": generation,
                            "realized_host_feedback": 0.5,
                            "mean_adult_richness": 2.0,
                            "mean_adult_gene_diversity": 0.5,
                        }
                    )
            (output / "pause.json").write_text(
                json.dumps({"last_completed_generation": 101}), encoding="utf-8"
            )
            run = {
                "run_id": "synthetic-run",
                "cell_id": "synthetic-cell",
                "seed_block_id": "sb0001",
            }
            cell = {
                "cell_id": "synthetic-cell",
                "cell": "synthetic",
                "panel": "H-by-B",
                "H": "100",
                "B": "10",
                "alpha_target": "0.5",
                "alpha": "0.5",
                "m": "0.1",
            }
            rows, input_record = wave2_analysis._new_run_rows(
                run,
                cell,
                output,
                {0: 500_000_000, 1: 500_000_000},
            )

            self.assertEqual(len(rows), 101)
            self.assertAlmostEqual(rows[-1]["TV"], 100 / 1_000_000_000)
            self.assertEqual(rows[-1]["realized_host_feedback"], 0.5)
            self.assertEqual(input_record["reached_generation"], 101)
            self.assertEqual(len(input_record["environment_prefix_sha256"]), 64)

    def test_complete_design_derives_all_registered_summary_tables(self) -> None:
        phase = REPOSITORY / "experiments/work/trophosome/p01-neutral-feedback"
        matrix = wave2_analysis._read_tsv(
            phase
            / "design/phase1-stage3-wave2-v210-adaptive-g1000-cells.tsv"
        )
        trajectories = []
        for cell in matrix:
            alpha = float(cell["alpha_target"])
            migration = float(cell["m"])
            hosts = int(cell["H"])
            bottleneck = int(cell["B"])
            for seed_number, (seed, _master) in enumerate(design.SEED_BLOCKS):
                for generation in range(101):
                    progress = generation / 100
                    seed_shift = (seed_number + 1) * 1e-5
                    tv = (
                        0.001
                        + alpha
                        * (1 - migration)
                        * progress
                        / math.sqrt(hosts * bottleneck)
                        + seed_shift
                    )
                    trajectories.append(
                        {
                            "run_id": f"{cell['cell_id']}-{seed}",
                            "cell_id": cell["cell_id"],
                            "cell": cell["cell"],
                            "seed_block_id": seed,
                            "panel": cell["panel"],
                            "H": hosts,
                            "B": bottleneck,
                            "HB": hosts * bottleneck,
                            "alpha_target": alpha,
                            "alpha": float(cell["alpha"]),
                            "m": migration,
                            "generation": generation,
                            "source_role": cell["initial_source_role"],
                            "source_run_id": "synthetic-source",
                            "D0": 100,
                            "shannon": math.log(30 - tv),
                            "simpson": 1 - 1 / (20 - tv),
                            "D1": 30 - tv,
                            "D2": 20 - tv,
                            "evenness": 0.8 - tv / 10,
                            "TV": tv,
                            "turnover": tv / 100,
                            "realized_host_feedback": alpha,
                            "mean_adult_richness": 3.0,
                            "mean_adult_gene_diversity": 0.5,
                        }
                    )

        tables = wave2_analysis._derive_tables(trajectories, matrix)

        self.assertEqual(len(tables["environment-trajectories-g100"]), 480 * 101)
        self.assertEqual(len(tables["run-endpoints-g100"]), 480)
        self.assertEqual(len(tables["run-tail-summaries-g100"]), 480)
        self.assertEqual(len(tables["cell-summaries-g100"]), 40 * 7)
        self.assertEqual(len(tables["h-by-b-paired-contrasts"]), 4 * 7)
        self.assertEqual(len(tables["alpha-by-m-contrasts"]), 3 * 7 * 7)
        self.assertEqual(len(tables["alpha-by-m-interactions"]), 3 * 6 * 7)
        self.assertTrue(tables["h-by-b-model-comparison"])


if __name__ == "__main__":
    unittest.main()
