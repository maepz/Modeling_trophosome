from __future__ import annotations

import importlib.util
import tempfile
import unittest
from pathlib import Path

from trophosome.pilot_report import (
    SLIDE28_MAIN_PURPOSES,
    _design_table_rows,
    _feedback_groups_by_host,
    _figure2_pooling_groups,
    _read_tsv,
    generate_pilot_report,
)

REPOSITORY = Path(__file__).resolve().parents[1]
ANALYSIS = (
    REPOSITORY
    / "experiments"
    / "work"
    / "trophosome"
    / "p01-neutral-feedback"
    / "analysis"
    / "s01-pilot-core-derived"
)
DESIGN = (
    REPOSITORY
    / "experiments"
    / "work"
    / "trophosome"
    / "p01-neutral-feedback"
    / "design"
    / "phase1-first-pilot-core-cells.tsv"
)
FULL_DESIGN = DESIGN.with_name("phase1-first-pilot-cells.tsv")
REPORT_DEPS = all(
    importlib.util.find_spec(module) is not None
    for module in ("matplotlib", "reportlab")
)


class PilotReportTests(unittest.TestCase):
    def test_design_table_uses_slide_28_main_purposes(self) -> None:
        rows = _design_table_rows(_read_tsv(FULL_DESIGN))
        self.assertEqual(rows[0][1], "Main purpose")
        purposes = {row[0]: row[1] for row in rows[1:]}
        self.assertEqual(purposes, SLIDE28_MAIN_PURPOSES)

    def test_feedback_series_hold_host_abundance_fixed(self) -> None:
        groups = _feedback_groups_by_host(_read_tsv(FULL_DESIGN))
        self.assertEqual(
            [[row["cell_id"].rsplit("-", 1)[-1] for row in group] for group in groups],
            [
                ["c0011", "c0016", "c0003", "c0018"],
                ["c0019", "c0020", "c0004", "c0012"],
            ],
        )
        for group in groups:
            for observed, expected in zip(
                [float(row["alpha"]) for row in group],
                [1 / 1001, 1 / 11, 0.5, 10 / 11],
                strict=True,
            ):
                self.assertAlmostEqual(observed, expected, places=12)
        self.assertEqual(
            [{int(float(row["H"])) for row in group} for group in groups],
            [{100}, {1000}],
        )

    def test_figure2_uses_feedback_and_mutation_matched_return_series(self) -> None:
        groups = _figure2_pooling_groups(_read_tsv(FULL_DESIGN))
        self.assertEqual(
            [
                (float(group[0]["alpha"]), float(group[0]["u"]))
                for group in groups
            ],
            [(0.5, 0.0), (0.5, 1e-10), (0.0909090909091, 0.0)],
        )
        self.assertEqual(
            [[row["cell_id"].rsplit("-", 1)[-1] for row in group] for group in groups],
            [
                ["c0003", "c0004", "c0005", "c0006"],
                ["c0009", "c0013", "c0014", "c0015"],
                ["c0016", "c0020", "c0017"],
            ],
        )

    @unittest.skipUnless(REPORT_DEPS, "optional report dependencies are not installed")
    def test_generates_pdf_markdown_and_diversity_glossary(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            artifacts = generate_pilot_report(
                ANALYSIS,
                DESIGN,
                root / "pilot.pdf",
                markdown_path=root / "pilot.md",
                report_date="2026-08-12",
            )
            self.assertTrue(artifacts.pdf.is_file())
            self.assertTrue(artifacts.pdf.read_bytes().startswith(b"%PDF-"))
            markdown = artifacts.markdown.read_text(encoding="utf-8")
            self.assertIn("Shannon entropy", markdown)
            self.assertIn("Simpson diversity", markdown)
            self.assertNotIn("Gini-Simpson", markdown)
            self.assertIn("Hill D1", markdown)
            self.assertIn("Appendix: glossary", markdown)
            self.assertLess(
                markdown.index("## Biological results"),
                markdown.index("## Quality control and computational feasibility"),
            )
            self.assertIn("same seed block", markdown)
            self.assertIn("Σᵢ pᵢ²", markdown)
            self.assertTrue((artifacts.assets / "endpoint-diversity.png").is_file())
            self.assertTrue((artifacts.assets / "feedback-bracket.png").is_file())


if __name__ == "__main__":
    unittest.main()
