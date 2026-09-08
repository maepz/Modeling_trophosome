from __future__ import annotations

import csv
import gzip
import hashlib
import json
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import numpy as np

REPOSITORY = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPOSITORY / "scripts"))
import compile_phase1_stage3_dbrda_inputs as dbrda  # noqa: E402


def _synthetic_source(
    root: Path,
    run_id: str,
    *,
    ancestral_lineage: int = 0,
    work: Path | None = None,
    scratch: Path | None = None,
) -> tuple[dbrda.Sample, Path, Path, Path]:
    work = work or root / "work"
    scratch = scratch or root / "scratch"
    config = work / f"configs/{run_id}.toml"
    output = scratch / f"outputs/{run_id}"
    config.parent.mkdir(parents=True, exist_ok=True)
    output.mkdir(parents=True, exist_ok=True)
    config.write_text(f'run_id = "{run_id}"\n', encoding="utf-8")
    environment = output / "environment_counts.csv"
    environment.write_text(
        "replicate,generation,strain_id,count\n"
        + "".join(f"0,{passage},{ancestral_lineage},10\n" for passage in range(101)),
        encoding="utf-8",
    )
    lineage = output / "strain_lineage_events.csv"
    lineage.write_text(
        "replicate,generation,host_id,within_host_generation,strain_id,"
        "parent_strain_id,mutational_depth\n",
        encoding="utf-8",
    )
    marker = output / "completion.json"
    marker.write_text(
        json.dumps(
            {
                "complete": True,
                "output_sizes": {
                    environment.name: environment.stat().st_size,
                    lineage.name: lineage.stat().st_size,
                },
            }
        ),
        encoding="utf-8",
    )
    digest = hashlib.sha256(config.read_bytes()).hexdigest()
    cell = {
        "cell_id": f"synthetic-{run_id}",
        "cell": run_id,
        "H": "100",
        "B": "10",
        "f": "0.01",
        "e": "10",
        "R": "1000",
        "alpha": "0.1",
        "alpha_target": "0.1",
        "m": "0.1",
        "u": "0",
    }
    sample = dbrda.Sample(
        analysis_set="synthetic",
        panel="synthetic",
        analysis_role="factorial-treatment",
        include_primary_dbrda=True,
        cell=cell,
        seed_block_id="sb0001",
        source_role="synthetic",
        source_cell_id=cell["cell_id"],
        source_run={
            "run_id": run_id,
            "cell_id": cell["cell_id"],
            "seed_block_id": "sb0001",
            "master_seed": "1",
            "config_path": f"configs/{run_id}.toml",
            "config_sha256": digest,
            "scratch_relative_path": f"outputs/{run_id}",
        },
    )
    return sample, work, scratch, marker


class Stage3DbrdaInputTests(unittest.TestCase):
    def test_master_design_has_three_valid_analysis_subsets(self) -> None:
        phase = REPOSITORY / "experiments/work/trophosome/p01-neutral-feedback"
        sets = dbrda._resolve_samples(phase)
        self.assertEqual(
            {name: len(rows) for name, rows in sets.items()},
            {
                "wave1_h_alpha_u": 300,
                "wave2a_h_by_b": 144,
                "wave2b_alpha_by_m": 336,
            },
        )
        self.assertEqual(
            {
                name: sum(sample.include_primary_dbrda for sample in samples)
                for name, samples in sets.items()
            },
            {
                "wave1_h_alpha_u": 288,
                "wave2a_h_by_b": 144,
                "wave2b_alpha_by_m": 336,
            },
        )
        all_samples = [sample for samples in sets.values() for sample in samples]
        self.assertEqual(len(all_samples), 780)
        self.assertEqual(
            len({sample.source_run["run_id"] for sample in all_samples}), 732
        )
        for samples in sets.values():
            source_ids = [sample.source_run["run_id"] for sample in samples]
            self.assertEqual(len(source_ids), len(set(source_ids)))

    def test_wave3_extends_the_master_without_changing_earlier_subsets(self) -> None:
        phase = REPOSITORY / "experiments/work/trophosome/p01-neutral-feedback"
        sets = dbrda._resolve_samples(phase, through_wave=3)
        self.assertEqual(
            {name: len(rows) for name, rows in sets.items()},
            {
                "wave1_h_alpha_u": 300,
                "wave2a_h_by_b": 144,
                "wave2b_alpha_by_m": 336,
                "wave3_bridge": 96,
            },
        )
        all_samples = [sample for samples in sets.values() for sample in samples]
        self.assertEqual(len(all_samples), 876)
        self.assertEqual(
            len({sample.source_run["run_id"] for sample in all_samples}), 828
        )
        wave3 = sets["wave3_bridge"]
        self.assertEqual(
            {sample.seed_block_id for sample in wave3},
            {f"sb{number:04d}" for number in range(1, 7)},
        )
        self.assertEqual(
            {sample.analysis_role for sample in wave3},
            {"bridge-core-treatment", "bridge-extension-treatment"},
        )
        self.assertTrue(all(sample.include_primary_dbrda for sample in wave3))

    def test_mutant_descendants_collapse_to_their_ancestral_root(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            work = root / "work"
            scratch = root / "scratch"
            config = work / "configs/run.toml"
            output = scratch / "outputs/run"
            config.parent.mkdir(parents=True)
            output.mkdir(parents=True)
            config.write_text("synthetic = true\n", encoding="utf-8")
            environment = output / "environment_counts.csv"
            environment.write_text(
                "replicate,generation,strain_id,count\n"
                "0,99,0,10\n"
                "0,100,0,4\n"
                "0,100,100,3\n"
                "0,100,101,3\n"
                "0,101,0,10\n",
                encoding="utf-8",
            )
            lineage = output / "strain_lineage_events.csv"
            lineage.write_text(
                "replicate,generation,host_id,within_host_generation,strain_id,"
                "parent_strain_id,mutational_depth\n"
                "0,4,0,2,100,0,1\n"
                "0,70,0,9,101,100,2\n",
                encoding="utf-8",
            )
            (output / "completion.json").write_text(
                json.dumps(
                    {
                        "complete": True,
                        "output_sizes": {
                            environment.name: environment.stat().st_size,
                            lineage.name: lineage.stat().st_size,
                        },
                    }
                ),
                encoding="utf-8",
            )
            digest = hashlib.sha256(config.read_bytes()).hexdigest()
            sample = dbrda.Sample(
                analysis_set="synthetic",
                panel="synthetic",
                analysis_role="factorial-treatment",
                include_primary_dbrda=True,
                cell={"cell": "c0001"},
                seed_block_id="sb0001",
                source_role="synthetic",
                source_cell_id="p01-s00-c0001",
                source_run={
                    "run_id": "synthetic-run",
                    "cell_id": "p01-s00-c0001",
                    "seed_block_id": "sb0001",
                    "config_path": "configs/run.toml",
                    "config_sha256": digest,
                    "scratch_relative_path": "outputs/run",
                },
            )
            frequencies, provenance = dbrda._ancestral_frequencies(
                sample, work=work, scratch=scratch, capacity=10
            )
            self.assertAlmostEqual(frequencies[0], 1.0)
            self.assertAlmostEqual(float(frequencies[1:].sum()), 0.0)
            self.assertEqual(provenance["raw_strain_richness_g100"], 3)
            self.assertEqual(provenance["ancestral_lineage_richness_g100"], 1)
            self.assertEqual(provenance["retained_mutant_strains_g100"], 2)

    def test_complete_ancestral_trajectory_collapses_transient_mutants(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            work = root / "work"
            scratch = root / "scratch"
            config = work / "configs/run.toml"
            output = scratch / "outputs/run"
            config.parent.mkdir(parents=True)
            output.mkdir(parents=True)
            config.write_text("synthetic = true\n", encoding="utf-8")
            environment = output / "environment_counts.csv"
            environmental_rows = ["replicate,generation,strain_id,count"]
            for passage in range(101):
                strain = 0 if passage < 50 else 100 if passage < 80 else 101
                environmental_rows.append(f"0,{passage},{strain},10")
            environment.write_text(
                "\n".join(environmental_rows) + "\n", encoding="utf-8"
            )
            lineage = output / "strain_lineage_events.csv"
            lineage.write_text(
                "replicate,generation,host_id,within_host_generation,strain_id,"
                "parent_strain_id,mutational_depth\n"
                "0,50,0,2,100,0,1\n"
                "0,80,0,9,101,100,2\n",
                encoding="utf-8",
            )
            (output / "completion.json").write_text(
                json.dumps(
                    {
                        "complete": True,
                        "output_sizes": {
                            environment.name: environment.stat().st_size,
                            lineage.name: lineage.stat().st_size,
                        },
                    }
                ),
                encoding="utf-8",
            )
            digest = hashlib.sha256(config.read_bytes()).hexdigest()
            sample = dbrda.Sample(
                analysis_set="synthetic",
                panel="synthetic",
                analysis_role="factorial-treatment",
                include_primary_dbrda=True,
                cell={"cell": "c0001"},
                seed_block_id="sb0001",
                source_role="synthetic",
                source_cell_id="p01-s00-c0001",
                source_run={
                    "run_id": "synthetic-run",
                    "cell_id": "p01-s00-c0001",
                    "seed_block_id": "sb0001",
                    "config_path": "configs/run.toml",
                    "config_sha256": digest,
                    "scratch_relative_path": "outputs/run",
                },
            )
            trajectory, provenance = dbrda._ancestral_trajectory(
                sample, work=work, scratch=scratch, capacity=10
            )
            self.assertEqual(trajectory.shape, (101, 100))
            np.testing.assert_allclose(trajectory[:, 0], 1.0)
            np.testing.assert_allclose(trajectory[:, 1:], 0.0)
            self.assertEqual(provenance["first_passage"], 0)
            self.assertEqual(provenance["last_passage"], 100)
            self.assertEqual(provenance["passages"], 101)
            self.assertEqual(provenance["retained_mutant_strains_g100"], 1)

    def test_tv_matrix_is_symmetric_and_uses_half_l1_distance(self) -> None:
        frequencies = np.zeros((3, 100), dtype=float)
        frequencies[0, 0] = 1.0
        frequencies[1, :2] = 0.5
        frequencies[2] = frequencies[1]
        observed = dbrda.total_variation_matrix(frequencies)
        np.testing.assert_allclose(observed, observed.T)
        np.testing.assert_allclose(np.diag(observed), 0)
        self.assertAlmostEqual(observed[0, 1], 0.5)
        self.assertAlmostEqual(observed[1, 2], 0.0)

    def test_environment_trajectory_requires_every_requested_passage(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "environment_counts.csv"
            path.write_text(
                "replicate,generation,strain_id,count\n"
                "0,0,0,10\n"
                "0,1,0,7\n"
                "0,1,1,3\n"
                "0,2,1,10\n",
                encoding="utf-8",
            )
            observed = dbrda._environment_trajectory(
                path, first_passage=0, last_passage=2
            )
            self.assertEqual(sorted(observed), [0, 1, 2])
            self.assertEqual(observed[1], {0: 7, 1: 3})
            with self.assertRaisesRegex(ValueError, "missing passage.*3"):
                dbrda._environment_trajectory(path, first_passage=0, last_passage=3)

    def test_master_triplet_keeps_identical_sample_order(self) -> None:
        cell = {
            "cell_id": "p01-s03-c0001",
            "cell": "c0001",
            "H": "100",
            "B": "10",
            "f": "0.01",
            "e": "10",
            "R": "1000",
            "alpha": "0.1",
            "alpha_target": "0.1",
            "m": "0.1",
            "u": "0",
        }
        samples = []
        for number in (1, 2):
            samples.append(
                dbrda.Sample(
                    analysis_set="synthetic",
                    panel="synthetic",
                    analysis_role="factorial-treatment",
                    include_primary_dbrda=True,
                    cell=cell,
                    seed_block_id=f"sb{number:04d}",
                    source_role="synthetic",
                    source_cell_id=cell["cell_id"],
                    source_run={
                        "run_id": f"run-{number}",
                        "master_seed": str(number),
                    },
                )
            )
        frequencies = np.zeros((2, 100), dtype=float)
        frequencies[0, 0] = 1
        frequencies[1, 1] = 1
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory)
            dbrda._write_master_triplet(output, samples, frequencies)
            filenames = (
                "x-explanatory-g100.tsv",
                "y-ancestral-frequencies-g100.tsv",
                "yprime-tv-g100.tsv",
            )
            observed = []
            for filename in filenames:
                with (output / filename).open(newline="", encoding="utf-8") as handle:
                    observed.append(
                        [
                            row["sample_id"]
                            for row in csv.DictReader(handle, delimiter="\t")
                        ]
                    )
            self.assertEqual(observed[0], observed[1])
            self.assertEqual(observed[1], observed[2])
            with (output / filenames[2]).open(newline="", encoding="utf-8") as handle:
                reader = csv.reader(handle, delimiter="\t")
                header = next(reader)
            self.assertEqual(header[1:], observed[0])

    def test_prc_table_has_one_complete_trajectory_per_analysis_sample(self) -> None:
        cell = {
            "cell_id": "p01-s03-c0001",
            "cell": "c0001",
            "H": "100",
            "B": "10",
            "f": "0.01",
            "e": "10",
            "R": "1000",
            "alpha": "0.1",
            "alpha_target": "0.1",
            "m": "0.1",
            "u": "0",
        }
        sample = dbrda.Sample(
            analysis_set="synthetic",
            panel="synthetic",
            analysis_role="factorial-treatment",
            include_primary_dbrda=True,
            cell=cell,
            seed_block_id="sb0001",
            source_role="synthetic",
            source_cell_id=cell["cell_id"],
            source_run={"run_id": "run-1", "master_seed": "1"},
        )
        trajectory = np.zeros((101, 100), dtype=float)
        trajectory[:, 0] = 1
        trajectory[100, 0] = 0.25
        trajectory[100, 1] = 0.75
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory)
            metadata = dbrda._write_prc_trajectory(
                output, [sample], {"run-1": trajectory}
            )
            self.assertEqual(metadata["rows"], 101)
            self.assertEqual(metadata["population_trajectories"], 1)
            with gzip.open(
                output / metadata["path"], "rt", newline="", encoding="utf-8"
            ) as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(len(rows), 101)
            self.assertEqual(rows[0]["passage"], "0")
            self.assertEqual(rows[-1]["passage"], "100")
            self.assertEqual(
                {row["population_sample_id"] for row in rows},
                {sample.sample_id},
            )
            self.assertEqual(rows[-1]["ancestral_001"], "0.75")

    def test_source_cache_resumes_and_invalidates_changed_inputs(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            sample, work, scratch, marker = _synthetic_source(root, "run-1")
            cache = root / "cache"
            trajectory_task = dbrda.SourceTask(
                sample=sample,
                work=work,
                scratch=scratch,
                capacity=10,
                kind="trajectory",
                cache_directory=cache,
                rebuild_cache=False,
            )
            _, first, _, first_status = dbrda._process_source_task(trajectory_task)
            _, second, _, second_status = dbrda._process_source_task(trajectory_task)
            self.assertEqual(first_status, "trajectory-written")
            self.assertEqual(second_status, "trajectory-hit")
            np.testing.assert_array_equal(first, second)

            endpoint_task = dbrda.SourceTask(
                sample=sample,
                work=work,
                scratch=scratch,
                capacity=10,
                kind="endpoint",
                cache_directory=cache,
                rebuild_cache=False,
            )
            _, endpoint, _, endpoint_status = dbrda._process_source_task(endpoint_task)
            self.assertEqual(endpoint_status, "trajectory-hit")
            np.testing.assert_array_equal(endpoint, first[100])

            marker_payload = json.loads(marker.read_text(encoding="utf-8"))
            marker_payload["cache_invalidation_test"] = True
            marker.write_text(json.dumps(marker_payload), encoding="utf-8")
            _, third, _, third_status = dbrda._process_source_task(trajectory_task)
            self.assertEqual(third_status, "trajectory-written")
            np.testing.assert_array_equal(first, third)

    def test_two_workers_preserve_source_results(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            first, work, scratch, _ = _synthetic_source(
                root, "run-1", ancestral_lineage=0
            )
            second, _, _, _ = _synthetic_source(root, "run-2", ancestral_lineage=1)
            values, provenance, cache_counts, issues = (
                dbrda._compile_source_populations(
                    [first, second],
                    work=work,
                    scratch=scratch,
                    capacity=10,
                    kind="endpoint",
                    workers=2,
                    cache_directory=None,
                    rebuild_cache=False,
                    progress_every=1,
                )
            )
            self.assertEqual(issues, [])
            self.assertEqual(set(values), {"run-1", "run-2"})
            self.assertEqual(set(provenance), {"run-1", "run-2"})
            self.assertEqual(cache_counts, {"disabled": 2})
            self.assertEqual(values["run-1"][0], 1.0)
            self.assertEqual(values["run-2"][1], 1.0)

    def test_compilation_modes_write_only_the_requested_tables(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            repository = Path(directory) / "repository"
            work = repository / "experiments/work/trophosome"
            phase = work / "p01-neutral-feedback"
            scratch = Path(directory) / "scratch"
            (work / "common/initial-populations").mkdir(parents=True)
            (work / "layout.local.json").write_text(
                json.dumps({"scratch": str(scratch)}), encoding="utf-8"
            )
            (work / "common/initial-populations/ip001-fisher100.json").write_text(
                json.dumps({"scaled_counts": [10] + [0] * 99}),
                encoding="utf-8",
            )
            sample, _, _, _ = _synthetic_source(
                Path(directory),
                "run-1",
                work=work,
                scratch=scratch,
            )
            cache = scratch / "compiler-cache"
            with patch.object(
                dbrda,
                "_resolve_samples",
                return_value={"synthetic": [sample]},
            ):
                full_output = phase / "analysis/full"
                dbrda.compile_inputs(
                    repository,
                    full_output,
                    workers=1,
                    cache_directory=cache,
                    progress_every=1,
                )
                endpoint_output = phase / "analysis/endpoint"
                dbrda.compile_inputs(
                    repository,
                    endpoint_output,
                    mode="endpoint",
                    workers=1,
                    cache_directory=cache,
                    progress_every=1,
                )
                prc_output = phase / "analysis/prc"
                dbrda.compile_inputs(
                    repository,
                    prc_output,
                    mode="prc",
                    workers=1,
                    cache_directory=cache,
                    progress_every=1,
                )

            endpoint_files = {
                "x-explanatory-g100.tsv",
                "y-ancestral-frequencies-g100.tsv",
                "yprime-tv-g100.tsv",
            }
            prc_file = "prc-ancestral-trajectories-g0-g100.tsv.gz"
            for filename in endpoint_files | {prc_file}:
                self.assertTrue((full_output / filename).is_file())
            for filename in endpoint_files:
                self.assertTrue((endpoint_output / filename).is_file())
                self.assertFalse((prc_output / filename).exists())
            self.assertFalse((endpoint_output / prc_file).exists())
            self.assertTrue((prc_output / prc_file).is_file())

            endpoint_audit = json.loads(
                (endpoint_output / "community-input-audit-g100.json").read_text(
                    encoding="utf-8"
                )
            )
            prc_audit = json.loads(
                (prc_output / "community-input-audit-g100.json").read_text(
                    encoding="utf-8"
                )
            )
            self.assertEqual(endpoint_audit["mode"], "endpoint")
            self.assertEqual(endpoint_audit["performance"]["cache_hits"], 1)
            self.assertIsNone(endpoint_audit["prc_trajectory"])
            self.assertEqual(prc_audit["mode"], "prc")
            self.assertEqual(prc_audit["performance"]["cache_hits"], 1)
            self.assertIsNone(prc_audit["master_triplet"])


if __name__ == "__main__":
    unittest.main()
