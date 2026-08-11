from __future__ import annotations

import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

REPOSITORY = Path(__file__).resolve().parents[1]
SCRIPT = REPOSITORY / "scripts" / "manage_project_layout.py"


class ProjectLayoutTests(unittest.TestCase):
    def run_script(
        self, *arguments: str, check: bool = True
    ) -> subprocess.CompletedProcess[str]:
        return subprocess.run(
            [sys.executable, str(SCRIPT), *arguments],
            text=True,
            capture_output=True,
            check=check,
        )

    def test_registration_and_lookup_support_future_parameters(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            self.run_script("init", "--root", directory)
            self.run_script("init", "--root", directory)
            project = Path(directory) / "work" / "trophosome"
            self.assertTrue((project / "registry" / "cells.csv").is_file())
            self.assertFalse((project / "p03-architecture").exists())
            self.assertTrue((project / "p01-neutral-feedback").is_dir())
            layout = json.loads((project / "layout.local.json").read_text())
            self.assertEqual(
                [row["phase_id"] for row in layout["initialized_phases"]],
                ["p01"],
            )

            registered = self.run_script(
                "register-cell",
                "--root",
                directory,
                "--cell-id",
                "p03-s03-c0018",
                "--label",
                "Lobed architecture mapping",
                "--architecture",
                "arch-lobed-v1",
                "--param",
                "H=10000",
                "--param",
                "f=1e-4",
                "--param",
                "u=1e-10",
                "--param",
                "c=1",
                "--param",
                "architecture_mode=lobed",
                "--param",
                "lobe_count=20",
                "--param",
                "migration_probability=0.01",
            )
            self.assertTrue(
                (
                    Path(directory)
                    / "scratch"
                    / "trophosome"
                    / "p03-architecture"
                    / "s03-parameter-map"
                ).is_dir()
            )
            layout = json.loads((project / "layout.local.json").read_text())
            self.assertEqual(
                [row["phase_id"] for row in layout["initialized_phases"]],
                ["p01", "p03"],
            )
            self.assertIn("h10000-f1em4-u1em10-c1-archlobed", registered.stdout)

            shown = self.run_script(
                "show-cell",
                "--root",
                directory,
                "--json",
                "c0018",
            )
            payload = json.loads(shown.stdout)
            self.assertEqual(payload["cell"]["cell_id"], "p03-s03-c0018")
            parameters = {
                row["parameter_name"]: row["value"]
                for row in payload["parameters"]
            }
            self.assertEqual(parameters["lobe_count"], "20")
            self.assertEqual(parameters["migration_probability"], "0.01")
            self.assertEqual(
                Path(payload["paths"]["scratch"]).name,
                "p03-s03-c0018",
            )

            shown_from_run_id = self.run_script(
                "show-cell",
                "--root",
                directory,
                "--json",
                "p03-s03-c0018-sb0007",
            )
            self.assertEqual(
                json.loads(shown_from_run_id.stdout)["cell"]["cell_id"],
                "p03-s03-c0018",
            )

            mnemonic = self.run_script(
                "mnemonic",
                "--root",
                directory,
                "c0018",
            )
            self.assertEqual(
                mnemonic.stdout.strip(),
                "h10000-f1em4-u1em10-c1-archlobed",
            )

            cell_id = self.run_script(
                "cell-id",
                "--root",
                directory,
                mnemonic.stdout.strip(),
            )
            self.assertEqual(cell_id.stdout.strip(), "p03-s03-c0018")

    def test_short_cell_number_must_be_unique(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            for cell_id in ("p01-s01-c0042", "p02-s01-c0042"):
                self.run_script(
                    "register-cell",
                    "--root",
                    directory,
                    "--cell-id",
                    cell_id,
                )
            ambiguous = self.run_script(
                "show-cell",
                "--root",
                directory,
                "c0042",
                check=False,
            )
            self.assertEqual(ambiguous.returncode, 2)
            self.assertIn("p01-s01-c0042", ambiguous.stderr)
            self.assertIn("p02-s01-c0042", ambiguous.stderr)

    def test_reverse_mnemonic_lookup_returns_all_matches_and_filters(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            for cell_id in ("p01-s01-c0041", "p01-s02-c0042"):
                self.run_script(
                    "register-cell",
                    "--root",
                    directory,
                    "--cell-id",
                    cell_id,
                    "--param",
                    "H=100",
                )

            all_matches = self.run_script(
                "cell-id",
                "--root",
                directory,
                "h100",
            )
            self.assertEqual(
                all_matches.stdout.splitlines(),
                ["p01-s01-c0041", "p01-s02-c0042"],
            )

            filtered = self.run_script(
                "cell-id",
                "--root",
                directory,
                "--stage",
                "s02",
                "h100",
            )
            self.assertEqual(filtered.stdout.strip(), "p01-s02-c0042")

    def test_duplicate_cell_ids_are_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            arguments = (
                "register-cell",
                "--root",
                directory,
                "--cell-id",
                "p01-s01-c0001",
                "--param",
                "H=100",
            )
            self.run_script(*arguments)
            duplicate = self.run_script(*arguments, check=False)
            self.assertEqual(duplicate.returncode, 2)
            self.assertIn("already registered", duplicate.stderr)


if __name__ == "__main__":
    unittest.main()
