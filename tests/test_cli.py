from __future__ import annotations

import contextlib
import io
import json
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from trophosome.cli import build_parser, main

REPOSITORY = Path(__file__).resolve().parents[1]


class CliTests(unittest.TestCase):
    def test_top_level_help_explains_biological_cycle_and_workflow(self) -> None:
        help_text = build_parser().format_help()
        self.assertIn("environmentally acquired microbial symbionts", help_text)
        self.assertIn("Typical workflow", help_text)
        self.assertIn("trophosome validate configs/smoke.toml", help_text)

    def test_run_help_explains_results_in_plain_language(self) -> None:
        output = io.StringIO()
        with self.assertRaises(SystemExit):
            with contextlib.redirect_stdout(output):
                main(["run", "-h"])
        help_text = output.getvalue()
        self.assertIn("within-host reproduction", help_text)
        self.assertIn("host_generation_summary.csv", help_text)
        self.assertIn("Use a different output directory", help_text)
        self.assertIn("--resume", help_text)

    def test_validate_help_says_that_it_does_not_run_the_model(self) -> None:
        output = io.StringIO()
        with self.assertRaises(SystemExit):
            with contextlib.redirect_stdout(output):
                main(["validate", "-h"])
        help_text = output.getvalue()
        self.assertIn("does not simulate hosts", help_text)
        self.assertIn("does not create a results folder", help_text)

    def test_validate_prints_resolved_configuration(self) -> None:
        output = io.StringIO()
        with contextlib.redirect_stdout(output):
            status = main(["validate", str(REPOSITORY / "configs" / "smoke.toml")])
        payload = json.loads(output.getvalue())
        self.assertEqual(status, 0)
        self.assertEqual(payload["model"], "wright_fisher_counts")

    def test_run_reports_output_location(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            destination = Path(directory) / "result"
            output = io.StringIO()
            with patch("trophosome.cli.run_simulation", return_value=(object(),)):
                with contextlib.redirect_stdout(output):
                    status = main(
                        [
                            "run",
                            str(REPOSITORY / "configs" / "smoke.toml"),
                            "--output",
                            str(destination),
                            "--repository",
                            str(REPOSITORY),
                        ]
                    )
            payload = json.loads(output.getvalue())
        self.assertEqual(status, 0)
        self.assertEqual(payload["host_generations"], 1)
        self.assertEqual(payload["output"], str(destination.resolve()))


if __name__ == "__main__":
    unittest.main()
