from __future__ import annotations

import csv
import json
import os
import shutil
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import numpy as np

from trophosome.config import load_config
from trophosome.simulation import _output_fields

REPOSITORY = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPOSITORY / "scripts"))
import analyse_phase1_stage3_wave1 as analysis  # noqa: E402
import build_phase1_stage3_wave1_report as reporting  # noqa: E402
import prepare_phase1_stage3_wave1 as design  # noqa: E402
import run_phase1_stage3_wave1 as runner  # noqa: E402
from run_phase1_first_pilot import _sha256  # noqa: E402
from run_phase1_second_pilot import _resolved_config_sha256  # noqa: E402


def make_repository(root: Path) -> tuple[Path, Path, list[dict[str, str]]]:
    """Portable frozen inputs, independent of the user's actual scratch location."""
    for source in design.build_files(REPOSITORY):
        destination = root / source.relative_to(REPOSITORY)
        destination.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(source, destination)
    work = root / "experiments/work/trophosome"
    initial = Path("common/initial-populations/ip001-fisher100.json")
    (work / initial).parent.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(
        REPOSITORY / "experiments/work/trophosome" / initial, work / initial
    )
    scratch = root / "scratch"
    scratch.mkdir()
    (work / "layout.local.json").write_text(json.dumps({"scratch": str(scratch)}))
    for name in (
        "scripts/analyse_phase1_stage3_wave1.py",
        "scripts/prepare_phase1_stage3_wave1.py",
        "scripts/run_phase1_stage3_wave1.py",
        "src/trophosome/stage3_report.py",
        "src/trophosome/second_pilot_report.py",
        "experiments/work/trophosome/p01-neutral-feedback/analysis/analyse_first_pilot.py",
        "experiments/work/trophosome/p01-neutral-feedback/analysis/analyse_second_pilot.py",
    ):
        destination = root / name
        destination.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(REPOSITORY / name, destination)
    return runner.load_rows(root)


def make_completed_run(row: dict, work: Path, scratch: Path) -> Path:
    """Clearly synthetic stable community; no biological simulation is run."""
    config = load_config(work / row["config_path"])
    output = scratch / row["scratch_relative_path"]
    output.mkdir(parents=True, exist_ok=True)
    for name, fields in _output_fields(config).items():
        (output / name).write_text(",".join(fields) + "\n")
    initial = config.environment.initial_counts
    with (output / "environment_counts.csv").open("a") as handle:
        for generation in range(design.HOST_GENERATIONS + 1):
            for strain, count in enumerate(initial):
                handle.write(
                    f"0,{generation},{strain},{count},{count / design.K},1.0,1.0\n"
                )
    with (output / "host_generation_summary.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=_output_fields(config)["host_generation_summary.csv"]
        )
        writer.writeheader()
        writer.writerows(
            {"replicate": 0, "host_generation": g}
            for g in range(1, design.HOST_GENERATIONS + 1)
        )
    with (output / "migration_counts.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=_output_fields(config)["migration_counts.csv"]
        )
        writer.writeheader()
        writer.writerows(
            {
                "replicate": 0,
                "generation": g,
                "emigrant_count": 100000000,
                "immigrant_count": 100000000,
            }
            for g in range(1, design.HOST_GENERATIONS + 1)
        )
    np.savez_compressed(
        output / "final_environment_rep000.npz",
        genotype_ids=np.arange(100),
        counts=np.asarray(initial),
    )
    (output / "resolved_config.json").write_text(json.dumps(config.to_dict()))
    (output / "provenance.json").write_text(
        json.dumps({"source_sha256": "synthetic-test-source"})
    )
    (output / "execution-summary.json").write_text(
        json.dumps(
            {
                "run_id": row["run_id"],
                "status": "complete",
                "elapsed_seconds": 100,
                "resumed": False,
            }
        )
    )
    (output / "completion.json").write_text(
        json.dumps(
            {
                "complete": True,
                "source_sha256": "synthetic-test-source",
                "config_sha256": _resolved_config_sha256(config.to_dict()),
                "model_spec_version": design.MODEL_SPEC_VERSION,
                "software_version": design.SOFTWARE_VERSION,
                "output_schema_version": design.OUTPUT_SCHEMA_VERSION,
                "output_sizes": {
                    name: (output / name).stat().st_size
                    for name in _output_fields(config)
                },
                "final_environment_sha256": {
                    "final_environment_rep000.npz": _sha256(
                        output / "final_environment_rep000.npz"
                    )
                },
            }
        )
    )
    return output


class Stage3DesignTests(unittest.TestCase):
    def test_frozen_artifacts_and_registries(self):
        result = subprocess.run(
            [sys.executable, "scripts/prepare_phase1_stage3_wave1.py", "--verify"],
            cwd=REPOSITORY,
            capture_output=True,
            text=True,
            check=True,
        )
        self.assertIn("Verified 317 Stage 3 files", result.stdout)

    def test_exact_agreed_matrix_and_life_history(self):
        self.assertEqual(
            [
                (c.number, c.hosts, c.alpha_target, c.mutation_probability)
                for c in design.CELLS
            ],
            [
                (27 + 8 * h_index + 2 * a_index + u_index, h, a, u)
                for h_index, h in enumerate((100, 1000, 10000))
                for a_index, a in enumerate(("0.001", "0.01", "0.1", "0.99"))
                for u_index, u in enumerate(("0", "1e-10"))
            ],
        )
        self.assertEqual(design.TOTAL_OFFSPRING, 506775749305)
        for alpha, returned in zip(
            design.FEEDBACK_LEVELS,
            (1000000, 10100000, 111110000, 99000000000),
            strict=True,
        ):
            selected = [c for c in design.CELLS if c.alpha_target == alpha]
            self.assertEqual({c.total_return for c in selected}, {returned})
            for c in selected:
                self.assertEqual(
                    round(float(c.escape_fraction) * design.K) * c.hosts, returned
                )
                self.assertLess(abs(c.alpha - float(alpha)), 1.1e-6)
        self.assertEqual(
            len([c for c in design.CELLS if c.escape_range != "primary-range"]), 10
        )
        self.assertEqual(
            [
                float(c.escape_fraction)
                for c in design.CELLS
                if c.alpha_target == "0.99" and c.mutation_probability == "0"
            ],
            [0.99, 0.099, 0.0099],
        )

    def test_all_runs_match_counts_seeds_and_retention(self):
        work, _, rows = runner.load_rows(REPOSITORY)
        self.assertEqual(len(rows), 288)
        self.assertEqual(len({r["scratch_relative_path"] for r in rows}), 288)
        self.assertNotIn(design.CONTROL_ID, {r["cell_id"] for r in rows})
        self.assertEqual(len(design.matrix_rows()), 25)
        for row in rows:
            config = load_config(work / row["config_path"])
            self.assertEqual(
                config.seed, dict(design.SEED_BLOCKS)[row["seed_block_id"]]
            )
            self.assertEqual(config.host.host_generations, 100)
            self.assertLessEqual(config.host.population_size, 10000)
            self.assertEqual(config.migration.fraction, 0.1)
            self.assertEqual(
                config.migration.regional_counts, config.environment.initial_counts
            )
            self.assertFalse(config.evolution.within_host_selection)
            self.assertFalse(config.evolution.free_living_selection)
            self.assertEqual(
                config.output.host_counts_mode,
                "full" if config.host.population_size == 100 else "panel",
            )
            self.assertEqual(config.output.environment_counts_mode, "all")
            self.assertEqual(config.execution.workers, 2)

    def test_frozen_references_do_not_require_mutable_stage2_outputs(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            make_repository(root)
            self.assertEqual(design.verify_files(root), [])
            reference = (
                root
                / "experiments/work/trophosome/p01-neutral-feedback/design"
                / f"{design.EXPERIMENT_ID}-reference-endpoints.tsv"
            )
            reference.write_text(reference.read_text() + "\n")
            with self.assertRaisesRegex(ValueError, "checksum differs"):
                design.verify_files(root)

    def test_reused_references_are_passage100_with_twelve_seeds(self):
        work = REPOSITORY / "experiments/work/trophosome"
        rows, provenance = design._reference_files(work)
        self.assertEqual(len(rows), 72)
        self.assertEqual({int(r["generation"]) for r in rows}, {100})
        control = [r for r in rows if r["source_role"] == "reused-control"]
        self.assertEqual(len(control), 12)
        self.assertEqual({r["cell_id"] for r in control}, {design.CONTROL_ID})
        self.assertIn("passage-100", provenance["note"])
        self.assertTrue(all("TV_tail_sd" in r for r in rows))

    def test_dry_run_does_not_create_scratch_runs(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            _, scratch, _ = make_repository(root)
            result = subprocess.run(
                [
                    sys.executable,
                    str(REPOSITORY / "scripts/run_phase1_stage3_wave1.py"),
                    "--repository",
                    str(root),
                    "--dry-run",
                ],
                capture_output=True,
                text=True,
                check=True,
            )
            self.assertIn("288 new populations", result.stdout)
            self.assertEqual(list(scratch.iterdir()), [])

    def test_hpc_wrapper_handles_conda_unset_variables_and_forwards_arguments(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            commands, env_bin = root / "commands", root / "environment/bin"
            commands.mkdir()
            env_bin.mkdir(parents=True)
            mamba = commands / "mamba"
            mamba.write_text(
                "#!/usr/bin/env bash\n"
                '[[ "$*" == "shell hook -s bash" ]] || exit 97\n'
                "cat <<'HOOK'\n"
                'mamba() { [[ "$*" == "activate trophosome" ]] || return 98; '
                ': "$TROPHOSOME_TEST_UNSET_BACKUP"; '
                'export PATH="$FAKE_ENV_BIN:$PATH"; }\n'
                "HOOK\n"
            )
            mamba.chmod(0o755)
            python = env_bin / "python"
            python.write_text(
                "#!/usr/bin/env bash\n"
                '[[ "$OMP_NUM_THREADS" == 1 ]] || exit 99\n'
                'printf "%s\\n" "$@" > "$FAKE_CAPTURE"\n'
            )
            python.chmod(0o755)
            capture = root / "arguments.txt"
            environment = {
                **os.environ,
                "PATH": f"{commands}:{os.environ['PATH']}",
                "FAKE_ENV_BIN": str(env_bin),
                "FAKE_CAPTURE": str(capture),
                "TROPHOSOME_STAGE3_JOBS": "4",
            }
            subprocess.run(
                [
                    "bash",
                    str(REPOSITORY / "scripts/hpc/launch_phase1_stage3_wave1.sh"),
                    "--smoke-only",
                    "--dry-run",
                ],
                env=environment,
                check=True,
                capture_output=True,
            )
            self.assertEqual(
                capture.read_text().splitlines(),
                [
                    str(REPOSITORY / "scripts/run_phase1_stage3_wave1.py"),
                    "--repository",
                    str(REPOSITORY),
                    "--jobs",
                    "4",
                    "--smoke-only",
                    "--dry-run",
                ],
            )


class Stage3AuditTests(unittest.TestCase):
    def test_full_launch_is_blocked_before_smoke_and_creates_no_runs(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            _, scratch, _ = make_repository(root)
            result = subprocess.run(
                [
                    sys.executable,
                    str(REPOSITORY / "scripts/run_phase1_stage3_wave1.py"),
                    "--repository",
                    str(root),
                ],
                capture_output=True,
                text=True,
                check=False,
            )
            self.assertEqual(result.returncode, 1)
            self.assertIn("Full wave not launched", result.stdout)
            self.assertEqual(list(scratch.iterdir()), [])

    def test_source_location_is_the_committed_checkout(self):
        runner.verify_source_location(sys.executable, REPOSITORY)
        with tempfile.TemporaryDirectory() as directory:
            with self.assertRaisesRegex(RuntimeError, "not importing this frozen"):
                runner.verify_source_location(sys.executable, Path(directory))

    def test_completion_hashes_and_corrupt_output(self):
        with tempfile.TemporaryDirectory() as directory:
            work, scratch, rows = make_repository(Path(directory))
            output = make_completed_run(rows[0], work, scratch)
            self.assertEqual(
                runner.completion_issues(rows[:1], work=work, scratch=scratch), []
            )
            (output / "environment_counts.csv").write_text("damaged")
            self.assertIn(
                "committed size differs",
                runner.completion_issues(rows[:1], work=work, scratch=scratch)[0],
            )

    def test_manifest_cannot_use_raw_hash_as_completed_resolved_hash(self):
        with tempfile.TemporaryDirectory() as directory:
            work, scratch, rows = make_repository(Path(directory))
            output = make_completed_run(rows[0], work, scratch)
            path = output / "completion.json"
            completed = json.loads(path.read_text())
            completed["config_sha256"] = rows[0]["config_sha256"]
            path.write_text(json.dumps(completed))
            self.assertIn(
                "config_sha256 differs",
                runner.completion_issues(rows[:1], work=work, scratch=scratch)[0],
            )

    def test_smoke_gate_requires_three_complete_populations(self):
        with tempfile.TemporaryDirectory() as directory:
            work, scratch, rows = make_repository(Path(directory))
            self.assertFalse(
                runner.smoke_assessment(rows, work=work, scratch=scratch)["passed"]
            )
            selected = [
                r
                for r in rows
                if r["cell"] in design.SMOKE_CELLS and r["seed_block_id"] == "sb0001"
            ]
            for row in selected:
                make_completed_run(row, work, scratch)
            self.assertTrue(
                runner.smoke_assessment(rows, work=work, scratch=scratch)["passed"]
            )
            path = (
                scratch
                / selected[0]["scratch_relative_path"]
                / "execution-summary.json"
            )
            record = json.loads(path.read_text())
            record["elapsed_seconds"] = 49 * 3600
            path.write_text(json.dumps(record))
            self.assertFalse(
                runner.smoke_assessment(rows, work=work, scratch=scratch)["passed"]
            )

    def test_incomplete_report_never_builds_or_launches(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            make_repository(root)
            with patch.object(reporting, "generate_report") as builder:
                with self.assertRaisesRegex(RuntimeError, "completion gate failed"):
                    reporting.build(root)
                builder.assert_not_called()

    def test_duplicate_environment_rows_are_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "counts.csv"
            path.write_text("replicate,generation,strain_id,count\n0,0,0,1\n0,0,0,1\n")
            with self.assertRaisesRegex(ValueError, "duplicate"):
                analysis.read_trajectory(path)


class Stage3AnalysisTests(unittest.TestCase):
    def test_tv_interaction_pairs_seeds_and_uses_difference_in_differences(self):
        rows = []
        for c in design.CELLS:
            for seed, master in design.SEED_BLOCKS:
                rows.append(
                    {
                        "cell_id": c.cell_id,
                        "seed_block_id": seed,
                        "TV": float(c.alpha_target) * (1 if c.hosts == 100 else 0.5)
                        + master / 10000,
                    }
                )
        result = analysis.tv_interactions(rows[::-1])
        self.assertEqual(len(result), 18)
        selected = [r for r in result if r["H_low"] == 100 and r["mutation"] == "0"]
        self.assertTrue(all(r["mean"] < 0 for r in selected))
        self.assertTrue(
            all(
                abs(r["mean"]) < 1e-12
                for r in result
                if r["mutation"] == "mutation-modification"
            )
        )

    def test_migration_audit_uses100_not250_passages(self):
        with tempfile.TemporaryDirectory() as directory:
            work, scratch, rows = make_repository(Path(directory))
            output = make_completed_run(rows[0], work, scratch)
            analysis.audit_migration(output / "migration_counts.csv")
            with (output / "migration_counts.csv").open("a") as handle:
                handle.write("0,101,0,100000000,100000000\n")
            with self.assertRaisesRegex(ValueError, "migration generations"):
                analysis.audit_migration(output / "migration_counts.csv")

    def test_twelve_populations_define_uncertainty(self):
        self.assertEqual(analysis.mean_ci90([4.0] * 12), (4.0, 4.0, 4.0))
        with self.assertRaises(ValueError):
            analysis.mean_ci90([4.0] * 250)
        mean, low, high = analysis.mean_ci90(list(range(1, 13)))
        self.assertAlmostEqual(mean, 6.5)
        self.assertAlmostEqual(high - mean, mean - low)

    def test_seed_pairing_and_margin_classification(self):
        matrix = design.matrix_rows()
        endpoints, references = [], []
        for stage, numbers, target in (
            ("s03", range(27, 51), endpoints),
            ("s02", range(21, 27), references),
        ):
            for number in numbers:
                for seed, _ in design.SEED_BLOCKS:
                    target.append(
                        {
                            "cell_id": f"p01-{stage}-c{number:04d}",
                            "seed_block_id": seed,
                            "D0": 100.0,
                            "D1": 30.0,
                            "D2": 20.0,
                            "evenness": 0.7,
                            "TV": 0.01,
                        }
                    )
        contrasts = analysis.paired_contrasts(endpoints, references[::-1], matrix)
        self.assertTrue(
            all(r["mean"] == 0 and r["status"] == "equivalent" for r in contrasts)
        )
        categories = analysis.classify(contrasts, matrix)
        self.assertTrue(
            all(
                r["classification"] == "negligible_for_five_measured_statistics"
                for r in categories
            )
        )
        with self.assertRaisesRegex(ValueError, "duplicate"):
            analysis.paired_contrasts(endpoints + endpoints[:1], references, matrix)

    def test_full_synthetic_analysis_and_pdf(self):
        # A persistent QA destination can be requested for page-by-page visual review.
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            work, scratch, rows = make_repository(root)
            for row in rows:
                make_completed_run(row, work, scratch)
            # The report builder itself lives in the actual repository; record its
            # content in the synthetic repository so input fingerprinting is portable.
            shutil.copyfile(
                REPOSITORY / "scripts/build_phase1_stage3_wave1_report.py",
                root / "scripts/build_phase1_stage3_wave1_report.py",
            )
            with patch.object(
                reporting,
                "__file__",
                str(root / "scripts/build_phase1_stage3_wave1_report.py"),
            ):
                artifacts = reporting.build(root)
            self.assertTrue(artifacts[0].read_bytes().startswith(b"%PDF-"))
            text = artifacts[1].read_text()
            self.assertIn("Audit: `PASS`", text)
            self.assertIn("not evidence", text)
            self.assertIn("c0044", text)
            derived = (
                work
                / "p01-neutral-feedback/analysis"
                / f"{design.STAGE_DIRECTORY}-derived"
            )
            summary = json.loads((derived / "analysis-summary.json").read_text())
            self.assertEqual(summary["populations"], 288)
            self.assertEqual(
                summary["primary_populations_including_reused_control"], 300
            )
            self.assertEqual(summary["cells"], 25)
            self.assertEqual(summary["primary_endpoint"], 100)
            self.assertEqual(summary["reference_populations"], 72)
            self.assertIn("Simpson", text)
            self.assertIn("Shannon", text)
            self.assertTrue(
                json.loads((derived / "report-completion.json").read_text())["complete"]
            )
            if destination := os.environ.get("TROPHOSOME_STAGE3_QA"):
                shutil.copytree(
                    root / "output", Path(destination) / "output", dirs_exist_ok=True
                )
                shutil.copytree(
                    root / "docs", Path(destination) / "docs", dirs_exist_ok=True
                )


if __name__ == "__main__":
    unittest.main()
