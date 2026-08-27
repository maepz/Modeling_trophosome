from __future__ import annotations

import csv
import importlib.util
import json
import os
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

from trophosome.config import load_config

REPOSITORY = Path(__file__).resolve().parents[1]
WORK = REPOSITORY / "experiments/work/trophosome"
PHASE = WORK / "p01-neutral-feedback"
VARIANT = "v210-m010-g250"


def _load_shared_runner():
    path = REPOSITORY / "scripts/run_phase1_first_pilot_v2_1.py"
    sys.path.insert(0, str(path.parent))
    try:
        spec = importlib.util.spec_from_file_location("shared_pilot_runner_test", path)
        if spec is None or spec.loader is None:
            raise RuntimeError("could not load shared pilot runner")
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        return module
    finally:
        sys.path.pop(0)


def _load_second_runner():
    path = REPOSITORY / "scripts/run_phase1_second_pilot.py"
    sys.path.insert(0, str(path.parent))
    try:
        spec = importlib.util.spec_from_file_location("second_pilot_runner_test", path)
        if spec is None or spec.loader is None:
            raise RuntimeError("could not load second-pilot runner")
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        return module
    finally:
        sys.path.pop(0)


class Phase1SecondPilotArtifactTests(unittest.TestCase):
    def test_generated_second_pilot_verifies(self) -> None:
        result = subprocess.run(
            [sys.executable, "scripts/prepare_phase1_second_pilot.py", "--verify"],
            cwd=REPOSITORY,
            check=True,
            capture_output=True,
            text=True,
        )
        self.assertIn("Verified 82 second-pilot files", result.stdout)

    def test_frozen_design_has_six_sentinels_and_h_is_at_most_10000(self) -> None:
        path = PHASE / "design" / f"phase1-second-pilot-{VARIANT}-cells.tsv"
        with path.open(encoding="utf-8") as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        self.assertEqual(len(rows), 6)
        self.assertEqual({int(row["host_generations"]) for row in rows}, {250})
        self.assertEqual({row["m"] for row in rows}, {"0.1"})
        self.assertLessEqual(max(int(row["H"]) for row in rows), 10_000)

        few = next(row for row in rows if row["cell"] == "c0022")
        many = next(row for row in rows if row["cell"] == "c0024")
        self.assertEqual((int(few["H"]), int(many["H"])), (100, 10_000))
        self.assertEqual(int(few["R"]), int(many["R"]))
        self.assertEqual(many["source_first_pilot_cell"], "p01-s01-c0005")
        self.assertEqual(many["host_counts_mode"], "panel")

    def test_all_72_runs_freeze_outputs_migration_and_neutral_selection(self) -> None:
        path = PHASE / "manifests" / f"phase1-second-pilot-{VARIANT}-runs.tsv"
        with path.open(encoding="utf-8") as handle:
            runs = list(csv.DictReader(handle, delimiter="\t"))
        self.assertEqual(len(runs), 72)
        self.assertEqual(len({row["run_id"] for row in runs}), 72)
        self.assertEqual(len({row["scratch_relative_path"] for row in runs}), 72)
        self.assertEqual(
            {row["seed_block_id"] for row in runs},
            {f"sb{number:04d}" for number in range(1, 13)},
        )
        for row in runs:
            config = load_config(WORK / row["config_path"])
            self.assertEqual(config.host.host_generations, 250)
            self.assertEqual(config.output.environment_counts_mode, "all")
            self.assertEqual(config.migration.mode, "fixed_regional_pool")
            self.assertEqual(config.migration.fraction, 0.1)
            self.assertEqual(
                config.migration.regional_counts,
                config.environment.initial_counts,
            )
            self.assertFalse(config.evolution.within_host_selection)
            self.assertFalse(config.evolution.free_living_selection)

    def test_manifest_freezes_resource_and_report_gates(self) -> None:
        path = PHASE / "manifests" / f"phase1-second-pilot-{VARIANT}-manifest.json"
        payload = json.loads(path.read_text(encoding="utf-8"))
        resources = payload["resource_projection_from_first_pilot"]
        self.assertLess(resources["aggregate_output_gib"], 10)
        self.assertLess(resources["aggregate_runtime_hours"], 40)
        self.assertEqual(
            payload["automatic_reporting"]["completion_gate"],
            "all 72 populations complete and internally valid",
        )
        self.assertEqual(
            payload["automatic_reporting"]["report_only_option"], "--report-only"
        )

    def test_shared_runner_prepares_a_second_pilot_scratch_manifest(self) -> None:
        manifest = PHASE / "manifests" / f"phase1-second-pilot-{VARIANT}-runs.tsv"
        with manifest.open(encoding="utf-8") as handle:
            row = next(csv.DictReader(handle, delimiter="\t"))
        runner = _load_shared_runner()
        with tempfile.TemporaryDirectory() as directory:
            scratch = Path(directory)
            runner._prepare_scratch([row], work=WORK, scratch=scratch)
            payload = json.loads(
                (scratch / row["scratch_relative_path"] / "run.json").read_text(
                    encoding="utf-8"
                )
            )
        self.assertEqual(payload["experiment_id"], row["experiment_id"])
        self.assertEqual(payload["variant_id"], VARIANT)
        self.assertEqual(payload["pilot_tier"], "second-pilot")

    def test_report_gate_uses_the_resolved_configuration_hash(self) -> None:
        manifest = PHASE / "manifests" / f"phase1-second-pilot-{VARIANT}-runs.tsv"
        with manifest.open(encoding="utf-8") as handle:
            row = next(csv.DictReader(handle, delimiter="\t"))
        runner = _load_second_runner()
        resolved = load_config(WORK / row["config_path"]).to_dict()
        resolved_hash = runner._resolved_config_sha256(resolved)
        self.assertNotEqual(resolved_hash, row["config_sha256"])
        with tempfile.TemporaryDirectory() as directory:
            scratch = Path(directory)
            output = scratch / row["scratch_relative_path"]
            output.mkdir(parents=True)
            (output / "resolved_config.json").write_text(
                json.dumps(resolved), encoding="utf-8"
            )
            for name in ("environment_counts.csv", "host_generation_summary.csv"):
                (output / name).write_bytes(b"")
            (output / "completion.json").write_text(
                json.dumps(
                    {
                        "complete": True,
                        "model_spec_version": runner.MODEL_SPEC_VERSION,
                        "output_schema_version": runner.OUTPUT_SCHEMA_VERSION,
                        "software_version": runner.SOFTWARE_VERSION,
                        "config_sha256": resolved_hash,
                        "output_sizes": {
                            "environment_counts.csv": 0,
                            "host_generation_summary.csv": 0,
                        },
                    }
                ),
                encoding="utf-8",
            )
            issues = runner._report_readiness_issues([row], work=WORK, scratch=scratch)
        self.assertEqual(issues, ["manifest contains 1 runs; expected 72"])

    def test_hpc_launcher_activates_mamba_and_forwards_report_only(self) -> None:
        launcher = REPOSITORY / "scripts/hpc/launch_phase1_second_pilot.sh"
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            commands = root / "commands"
            environment_bin = root / "environment/bin"
            commands.mkdir()
            environment_bin.mkdir(parents=True)
            capture = root / "python-arguments.txt"
            fake_mamba = commands / "mamba"
            fake_mamba.write_text(
                "#!/usr/bin/env bash\n"
                'if [[ "$1 $2 $3 $4" != "shell hook -s bash" ]]; then exit 97; fi\n'
                "cat <<'HOOK'\n"
                "mamba() {\n"
                '  [[ "$1" == "activate" && "$2" == "trophosome" ]] || return 98\n'
                '  : "$TROPHOSOME_TEST_UNSET_CONDA_BACKUP_VARIABLE"\n'
                '  export PATH="$FAKE_MAMBA_ENV_BIN:$PATH"\n'
                "}\n"
                "HOOK\n",
                encoding="utf-8",
            )
            fake_mamba.chmod(0o755)
            fake_python = environment_bin / "python"
            fake_python.write_text(
                '#!/usr/bin/env bash\nprintf \'%s\\n\' "$@" > "$FAKE_CAPTURE"\n',
                encoding="utf-8",
            )
            fake_python.chmod(0o755)
            environment = os.environ.copy()
            environment.update(
                {
                    "PATH": f"{commands}:{environment['PATH']}",
                    "FAKE_MAMBA_ENV_BIN": str(environment_bin),
                    "FAKE_CAPTURE": str(capture),
                    "TROPHOSOME_SECOND_PILOT_JOBS": "3",
                }
            )
            subprocess.run(
                ["bash", str(launcher), "--report-only"],
                cwd=REPOSITORY,
                env=environment,
                check=True,
                capture_output=True,
                text=True,
            )
            self.assertEqual(
                capture.read_text(encoding="utf-8").splitlines(),
                [
                    str(REPOSITORY / "scripts/run_phase1_second_pilot.py"),
                    "--repository",
                    str(REPOSITORY),
                    "--jobs",
                    "3",
                    "--report-only",
                ],
            )


if __name__ == "__main__":
    unittest.main()
