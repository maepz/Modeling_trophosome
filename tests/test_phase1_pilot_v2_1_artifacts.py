from __future__ import annotations

import csv
import json
import os
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

from trophosome.config import load_config

REPOSITORY = Path(__file__).resolve().parents[1]
WORK = REPOSITORY / "experiments" / "work" / "trophosome"
PHASE = WORK / "p01-neutral-feedback"
VARIANT = "v210-m010"


class Phase1PilotV21ArtifactTests(unittest.TestCase):
    def test_hpc_launcher_activates_environment_without_mamba_run(self) -> None:
        launcher = REPOSITORY / "scripts/hpc/launch_phase1_first_pilot_v2_1.sh"
        with tempfile.TemporaryDirectory() as temporary_directory:
            temporary = Path(temporary_directory)
            command_directory = temporary / "commands"
            environment_bin = temporary / "environment" / "bin"
            command_directory.mkdir()
            environment_bin.mkdir(parents=True)
            capture_path = temporary / "python-arguments.txt"

            fake_mamba = command_directory / "mamba"
            fake_mamba.write_text(
                "#!/usr/bin/env bash\n"
                'if [[ "$1 $2 $3 $4" != "shell hook -s bash" ]]; then\n'
                "  exit 97\n"
                "fi\n"
                "cat <<'HOOK'\n"
                "mamba() {\n"
                '  if [[ "$1" != "activate" || "$2" != "trophosome" ]]; then\n'
                "    return 98\n"
                "  fi\n"
                '  export PATH="$FAKE_MAMBA_ENV_BIN:$PATH"\n'
                "}\n"
                "HOOK\n",
                encoding="utf-8",
            )
            fake_mamba.chmod(0o755)

            fake_python = environment_bin / "python"
            fake_python.write_text(
                '#!/usr/bin/env bash\nprintf \'%s\\n\' "$@" > "$FAKE_PYTHON_CAPTURE"\n',
                encoding="utf-8",
            )
            fake_python.chmod(0o755)

            environment = os.environ.copy()
            environment.update(
                {
                    "PATH": f"{command_directory}:{environment['PATH']}",
                    "FAKE_MAMBA_ENV_BIN": str(environment_bin),
                    "FAKE_PYTHON_CAPTURE": str(capture_path),
                    "TROPHOSOME_MAMBA_ENV": "trophosome",
                    "TROPHOSOME_PILOT_JOBS": "3",
                }
            )
            subprocess.run(
                ["bash", str(launcher), "--prepare-only", "--cell", "c0001"],
                cwd=REPOSITORY,
                env=environment,
                check=True,
                capture_output=True,
                text=True,
            )

            arguments = capture_path.read_text(encoding="utf-8").splitlines()
            self.assertEqual(
                arguments,
                [
                    str(REPOSITORY / "scripts/run_phase1_first_pilot_v2_1.py"),
                    "--repository",
                    str(REPOSITORY),
                    "--jobs",
                    "3",
                    "--prepare-only",
                    "--cell",
                    "c0001",
                ],
            )

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
