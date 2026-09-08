#!/usr/bin/env python3
"""Build the self-contained Wave 2 passage-100 adaptive-decision report."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from datetime import UTC, datetime
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from reportlab.lib import colors  # noqa: E402
from reportlab.lib.colors import HexColor  # noqa: E402
from reportlab.lib.pagesizes import A4  # noqa: E402
from reportlab.lib.styles import ParagraphStyle, getSampleStyleSheet  # noqa: E402
from reportlab.lib.units import mm  # noqa: E402
from reportlab.platypus import (  # noqa: E402
    Image,
    PageBreak,
    Paragraph,
    SimpleDocTemplate,
    Spacer,
    Table,
    TableStyle,
)

EXPERIMENT_ID = "phase1-stage3-wave2-v210-adaptive-g1000"
DERIVED_PREFIX = EXPERIMENT_ID.replace(
    "phase1-stage3-wave2-", "s03-parameter-map-wave2-"
)
DERIVED_DIRECTORY = f"{DERIVED_PREFIX}-derived"
ACCENT = HexColor("#176B73")
INK = HexColor("#233238")
PALE = HexColor("#EAF4F4")
GREY = HexColor("#66757A")


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _reported_path(path: Path, repository: Path) -> str:
    try:
        return str(path.relative_to(repository))
    except ValueError:
        return str(path)


def _read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def _validate_inputs(
    decision: dict, diagnostics: list[dict[str, str]], cells: list[dict[str, str]]
) -> None:
    if decision.get("status") != "frozen-outcome-dependent-decision":
        raise ValueError("the passage-100 adaptive decision is not frozen")
    if decision.get("assessed_horizon") != 100:
        raise ValueError("the report requires the passage-100 decision")
    if len(cells) != 40:
        raise ValueError(f"expected 40 Wave 2 conditions, found {len(cells)}")
    if len(diagnostics) != 21:
        raise ValueError(f"expected 21 adaptive diagnostics, found {len(diagnostics)}")
    if len(decision.get("input_fingerprints", [])) != 409:
        raise ValueError("the decision must fingerprint 408 new populations and reuse")
    if decision.get("selected_populations") != 72:
        raise ValueError("the frozen decision must select 72 populations")
    if any(row["stable"] != "True" for row in diagnostics):
        raise ValueError("the committed diagnostic table contains an unstable row")


def _configure_axes(axis) -> None:
    axis.spines[["top", "right"]].set_visible(False)
    axis.grid(axis="y", alpha=0.18, linewidth=0.7)
    axis.tick_params(labelsize=8)


def _make_figures(
    diagnostics: list[dict[str, str]], assets: Path
) -> tuple[Path, Path]:
    assets.mkdir(parents=True, exist_ok=True)
    raw = [row for row in diagnostics if row["diagnostic"] == "raw-TV"]
    colours = {
        "0": "#5A6670",
        "0.01": "#2A9D8F",
        "0.1": "#E9A03B",
        "0.99": "#C8553D",
    }
    migrations = ("0", "0.001", "0.01")
    x = range(len(migrations))

    fig, axes = plt.subplots(1, 2, figsize=(9.4, 4.0), constrained_layout=True)
    for alpha in ("0", "0.01", "0.1"):
        selected = [
            next(row for row in raw if row["alpha"] == alpha and row["m"] == m)
            for m in migrations
        ]
        axes[0].plot(
            x,
            [float(row["late_mean"]) for row in selected],
            marker="o",
            linewidth=2,
            color=colours[alpha],
            label=f"alpha = {alpha}",
        )
    high = [
        next(row for row in raw if row["alpha"] == "0.99" and row["m"] == m)
        for m in migrations
    ]
    axes[1].plot(
        x,
        [float(row["late_mean"]) for row in high],
        marker="o",
        linewidth=2,
        color=colours["0.99"],
        label="alpha = 0.99",
    )
    for axis, title in zip(
        axes,
        ("Feedback alpha <= 0.1", "Strong feedback alpha = 0.99"),
        strict=True,
    ):
        _configure_axes(axis)
        axis.set_xticks(list(x), migrations)
        axis.set_xlabel("Regional replacement fraction m")
        axis.set_ylabel("Mean TV during passages 76-100")
        axis.set_title(title, fontsize=10, weight="bold")
        axis.legend(frameon=False, fontsize=8)
    axes[0].set_ylim(bottom=0)
    axes[1].set_ylim(bottom=0)
    late_path = assets / "late-window-tv.png"
    fig.savefig(late_path, dpi=220, bbox_inches="tight")
    plt.close(fig)

    ordered = sorted(raw, key=lambda row: (float(row["alpha"]), float(row["m"])))
    labels = [f"alpha={row['alpha']}, m={row['m']}" for row in ordered]
    mean = [float(row["mean_change"]) / float(row["margin"]) for row in ordered]
    low = [float(row["ci90_lower"]) / float(row["margin"]) for row in ordered]
    high_ci = [float(row["ci90_upper"]) / float(row["margin"]) for row in ordered]
    y = list(range(len(ordered)))
    fig, axis = plt.subplots(figsize=(8.8, 5.4), constrained_layout=True)
    axis.axvspan(-1, 1, color="#DDEFE9", alpha=0.8, label="Pre-specified stable range")
    axis.axvline(0, color="#53656B", linewidth=1)
    axis.axvline(-1, color="#2A9D8F", linestyle="--", linewidth=1)
    axis.axvline(1, color="#2A9D8F", linestyle="--", linewidth=1)
    for index, row in enumerate(ordered):
        colour = colours[row["alpha"]]
        axis.plot(
            [low[index], high_ci[index]],
            [index, index],
            color=colour,
            linewidth=2,
        )
        axis.scatter(mean[index], index, color=colour, s=28, zorder=3)
    axis.set_yticks(y, labels)
    axis.set_xlabel("Change between time windows / allowed stability margin")
    axis.set_title("All raw-TV confidence intervals remain inside the decision margin")
    axis.spines[["top", "right", "left"]].set_visible(False)
    axis.grid(axis="x", alpha=0.15)
    axis.legend(frameon=False, loc="lower right", fontsize=8)
    axis.invert_yaxis()
    stability_path = assets / "stability-diagnostics.png"
    fig.savefig(stability_path, dpi=220, bbox_inches="tight")
    plt.close(fig)
    return late_path, stability_path


def _styles():
    styles = getSampleStyleSheet()
    styles.add(
        ParagraphStyle(
            "ReportTitle",
            parent=styles["Title"],
            fontName="Helvetica-Bold",
            fontSize=23,
            leading=27,
            textColor=INK,
            spaceAfter=7 * mm,
        )
    )
    styles.add(
        ParagraphStyle(
            "ReportSubtitle",
            parent=styles["Normal"],
            fontSize=11,
            leading=15,
            textColor=GREY,
            spaceAfter=7 * mm,
        )
    )
    styles.add(
        ParagraphStyle(
            "Section",
            parent=styles["Heading2"],
            fontName="Helvetica-Bold",
            fontSize=15,
            leading=18,
            textColor=ACCENT,
            spaceBefore=4 * mm,
            spaceAfter=3 * mm,
        )
    )
    styles.add(
        ParagraphStyle(
            "BodyReport",
            parent=styles["BodyText"],
            fontSize=9.4,
            leading=13.2,
            textColor=INK,
            spaceAfter=2.4 * mm,
        )
    )
    styles.add(
        ParagraphStyle(
            "Caption",
            parent=styles["BodyText"],
            fontSize=8.2,
            leading=11,
            textColor=GREY,
            spaceBefore=1.5 * mm,
            spaceAfter=3 * mm,
        )
    )
    styles.add(
        ParagraphStyle(
            "Small",
            parent=styles["BodyText"],
            fontSize=7.4,
            leading=9.4,
            textColor=INK,
        )
    )
    return styles


def _page(canvas, document) -> None:
    canvas.saveState()
    width, height = A4
    canvas.setStrokeColor(PALE)
    canvas.line(18 * mm, height - 15 * mm, width - 18 * mm, height - 15 * mm)
    canvas.setFont("Helvetica", 7.5)
    canvas.setFillColor(GREY)
    canvas.drawString(18 * mm, 10 * mm, "Trophosome - Phase 1 Stage 3 Wave 2")
    canvas.drawRightString(width - 18 * mm, 10 * mm, f"Page {document.page}")
    canvas.restoreState()


def _build_pdf(
    pdf: Path,
    decision: dict,
    diagnostics: list[dict[str, str]],
    late_figure: Path,
    stability_figure: Path,
) -> None:
    pdf.parent.mkdir(parents=True, exist_ok=True)
    styles = _styles()
    doc = SimpleDocTemplate(
        str(pdf),
        pagesize=A4,
        leftMargin=18 * mm,
        rightMargin=18 * mm,
        topMargin=21 * mm,
        bottomMargin=17 * mm,
        title="Phase 1 Stage 3 Wave 2 passage-100 adaptive-decision report",
        author="Trophosome",
    )
    body = styles["BodyReport"]
    story = [
        Paragraph("Phase 1 Stage 3 Wave 2", styles["ReportTitle"]),
        Paragraph(
            "Passage-100 adaptive time-horizon report | model 2.1.0 | "
            "12 matched seed blocks",
            styles["ReportSubtitle"],
        ),
    ]
    summary = Table(
        [
            ["Primary Wave 2 design", "40 conditions; 480 population trajectories"],
            ["New HPC populations checked", "408"],
            ["Reused passage-100 trajectories", "72"],
            ["Adaptive conditions assessed", "12 low-migration conditions"],
            ["Diagnostics meeting stability rule", "21 of 21"],
            ["Continuation to passage 500", "6 anchor conditions; 72 populations"],
        ],
        colWidths=[61 * mm, 103 * mm],
    )
    summary.setStyle(
        TableStyle(
            [
                ("BACKGROUND", (0, 0), (0, -1), PALE),
                ("TEXTCOLOR", (0, 0), (-1, -1), INK),
                ("FONTNAME", (0, 0), (0, -1), "Helvetica-Bold"),
                ("FONTNAME", (1, 0), (1, -1), "Helvetica"),
                ("FONTSIZE", (0, 0), (-1, -1), 8.8),
                ("LEADING", (0, 0), (-1, -1), 11),
                ("GRID", (0, 0), (-1, -1), 0.4, colors.white),
                ("VALIGN", (0, 0), (-1, -1), "MIDDLE"),
                ("LEFTPADDING", (0, 0), (-1, -1), 6),
                ("RIGHTPADDING", (0, 0), (-1, -1), 6),
                ("TOPPADDING", (0, 0), (-1, -1), 5),
                ("BOTTOMPADDING", (0, 0), (-1, -1), 5),
            ]
        )
    )
    story.extend(
        [
            summary,
            Spacer(1, 5 * mm),
            Paragraph("Main conclusion", styles["Section"]),
            Paragraph(
                "All low-migration trajectories satisfied the pre-specified "
                "stability rule at passage 100. This means that their estimated "
                "change between passages 51-75 and 76-100 was small relative to "
                "the accepted biological margin. It does not mean that every "
                "trajectory was perfectly flat or that equilibrium has been proven.",
                body,
            ),
            Paragraph(
                "Because no optional condition failed the rule, only the six "
                "pre-specified anchors were selected for continuation: alpha=0 "
                "and alpha=0.1 at m=0, 0.001 and 0.01. These 72 populations provide "
                "a controlled long-term check through passage 500.",
                body,
            ),
            Paragraph("Scope of this report", styles["Section"]),
            Paragraph(
                "The pulled repository contains the frozen adaptive decision, "
                "21 diagnostic summaries and cryptographic fingerprints for all "
                "408 new trajectories. It does not contain the full endpoint and "
                "trajectory tables needed to analyse the H-by-B panel, the complete "
                "alpha-by-m surface, or D1, D2 and evenness. Consequently, this is "
                "an adaptive-horizon report, not the complete scientific Wave 2 "
                "report.",
                body,
            ),
            PageBreak(),
            Paragraph("Environmental displacement at low migration", styles["Section"]),
            Image(str(late_figure), width=171 * mm, height=73 * mm),
            Paragraph(
                "Figure 1. Mean total-variation distance from the initial "
                "environmental community during passages 76-100. The two panels use "
                "different vertical scales. Greater host feedback produced greater "
                "displacement. Within each feedback level, increasing regional "
                "replacement from 0 to 0.01 reduced "
                "the retained displacement, especially when feedback was strong.",
                styles["Caption"],
            ),
            Paragraph("Biological interpretation", styles["Section"]),
            Paragraph(
                "At alpha=0, the focal population remained essentially unchanged when "
                "there was no migration; the very small non-zero values at m=0.001 and "
                "0.01 reflect stochastic sampling during regional replacement. At "
                "alpha=0.99, mean late-window TV declined from 0.428 without migration "
                "to 0.312 at m=0.01. This is consistent with immigration counteracting "
                "the compositional signal created by repeated host return.",
                body,
            ),
            Paragraph(
                "These three migration levels were selected for deciding the time "
                "horizon because their environmental memory is longest. The complete "
                "seven-level migration experiment must be analysed before estimating "
                "the shape or strength of the full migration effect.",
                body,
            ),
            PageBreak(),
            Paragraph(
                "Why the trajectories passed the stability gate", styles["Section"]
            ),
            Image(str(stability_figure), width=169 * mm, height=103 * mm),
            Paragraph(
                "Figure 2. Change in raw TV between the two assessment windows, "
                "divided "
                "by the condition-specific stability margin. Points are paired means "
                "across 12 seed blocks and lines are 90% Student-t intervals. A "
                "complete "
                "interval inside -1 to +1 satisfies the frozen stability rule.",
                styles["Caption"],
            ),
            Paragraph(
                "Several confidence intervals lie above zero, particularly under "
                "strong "
                "feedback, so small directional change was still detectable. They pass "
                "because the entire interval remains smaller than the pre-specified "
                "tolerance: max(0.002 TV, 25% of the absolute late-window mean). The "
                "correct conclusion is therefore 'stable within the chosen biological "
                "margin', not 'no continuing change'.",
                body,
            ),
            PageBreak(),
            Paragraph("Frozen continuation decision", styles["Section"]),
        ]
    )
    selected_rows = [["Condition", "Feedback alpha", "Migration m", "Reason"]]
    cells = {
        "c0063": ("0", "0"),
        "c0064": ("0", "0.001"),
        "c0065": ("0", "0.01"),
        "c0077": ("0.1", "0"),
        "c0078": ("0.1", "0.001"),
        "c0079": ("0.1", "0.01"),
    }
    for cell in decision["selected_cells"]:
        alpha, migration = cells[cell]
        selected_rows.append([cell, alpha, migration, "Pre-specified anchor"])
    selected = Table(selected_rows, colWidths=[31 * mm, 36 * mm, 32 * mm, 65 * mm])
    selected.setStyle(
        TableStyle(
            [
                ("BACKGROUND", (0, 0), (-1, 0), ACCENT),
                ("TEXTCOLOR", (0, 0), (-1, 0), colors.white),
                ("FONTNAME", (0, 0), (-1, 0), "Helvetica-Bold"),
                ("FONTNAME", (0, 1), (-1, -1), "Helvetica"),
                ("FONTSIZE", (0, 0), (-1, -1), 8.2),
                ("GRID", (0, 0), (-1, -1), 0.35, colors.white),
                ("ROWBACKGROUNDS", (0, 1), (-1, -1), [colors.white, PALE]),
                ("VALIGN", (0, 0), (-1, -1), "MIDDLE"),
                ("TOPPADDING", (0, 0), (-1, -1), 5),
                ("BOTTOMPADDING", (0, 0), (-1, -1), 5),
            ]
        )
    )
    story.extend(
        [
            selected,
            Spacer(1, 4 * mm),
            Paragraph(
                "The continuation is deliberately diagnostic. It preserves both "
                "no-return controls and alpha=0.1 treatments at each eligible "
                "migration level, even though all passed the passage-100 rule. "
                "Optional alpha=0.01 "
                "and 0.99 conditions were not extended because they were stable and no "
                "same-m control was unresolved.",
                body,
            ),
            Paragraph("Audit and limitations", styles["Section"]),
            Paragraph(
                "The frozen decision records one fingerprint for the reused trajectory "
                "file and 408 fingerprints for new passage-100 trajectory prefixes. "
                "The report reads, but does not alter, those records. Twelve "
                "independently simulated populations are the replicate unit for every "
                "diagnostic.",
                body,
            ),
            Paragraph(
                "This report cannot answer whether H and B have independent effects at "
                "fixed HB, whether alpha and m interact over their full ranges, or how "
                "diversity and evenness respond. Those questions require derived "
                "tables from the HPC scratch outputs. Passage-500 results will also "
                "be selected "
                "by this frozen decision and must be labelled exploratory.",
                body,
            ),
            Paragraph("Decision rule", styles["Section"]),
            Paragraph(
                "For each condition and seed block, mean TV in passages 51-75 was "
                "subtracted from mean TV in passages 76-100. A trajectory was "
                "classified "
                "as stable when the paired 90% interval for this change lay entirely "
                "within +/-max(0.002 TV, 25% of the absolute late-window mean). "
                "Positive-feedback cells were checked both as raw TV and after "
                "subtracting their "
                "same-m alpha=0 control.",
                body,
            ),
        ]
    )
    doc.build(story, onFirstPage=_page, onLaterPages=_page)


def _build_markdown(
    markdown: Path,
    decision: dict,
    diagnostics: list[dict[str, str]],
    late_figure: Path,
    stability_figure: Path,
) -> None:
    markdown.parent.mkdir(parents=True, exist_ok=True)
    relative_late = late_figure.relative_to(markdown.parent).as_posix()
    relative_stability = stability_figure.relative_to(markdown.parent).as_posix()
    raw = [row for row in diagnostics if row["diagnostic"] == "raw-TV"]
    lines = [
        "# Phase 1 Stage 3 Wave 2 passage-100 adaptive-decision report",
        "",
        "## Main conclusion",
        "",
        "All 21 raw-TV and host-induced-TV diagnostics met the pre-specified "
        "stability rule at passage 100. Six pre-specified anchor conditions "
        f"({', '.join(decision['selected_cells'])}) nevertheless continue to "
        "passage 500, representing 72 populations.",
        "",
        "Stable means that the paired 90% interval for continuing change remained "
        "inside the chosen biological margin. It does not mean that every trajectory "
        "was flat or that equilibrium has been demonstrated.",
        "",
        "## Scope",
        "",
        "The repository contains the frozen adaptive decision, its 21 diagnostics "
        "and fingerprints for all 408 new trajectories, but not the derived endpoint "
        "tables needed for the complete H-by-B and alpha-by-m scientific report.",
        "",
        "## Late-window environmental displacement",
        "",
        f"![Late-window total-variation distance]({relative_late})",
        "",
        "Regional replacement reduced the late-window host signal over the three "
        "low-migration values used for the adaptive decision. The complete migration "
        "surface is not available in the pulled summaries.",
        "",
        "## Stability diagnostics",
        "",
        f"![Stability diagnostics]({relative_stability})",
        "",
        "| Cell | alpha | m | Late mean TV | Window change | 90% interval | Margin |",
        "|---|---:|---:|---:|---:|---:|---:|",
    ]
    for row in raw:
        lines.append(
            f"| {row['cell']} | {row['alpha']} | {row['m']} | "
            f"{float(row['late_mean']):.6f} | {float(row['mean_change']):.6f} | "
            f"[{float(row['ci90_lower']):.6f}, {float(row['ci90_upper']):.6f}] | "
            f"{float(row['margin']):.6f} |"
        )
    lines.extend(
        [
            "",
            "## Frozen continuation",
            "",
            "Only the alpha=0 and alpha=0.1 anchors at m=0, 0.001 and 0.01 "
            "continue to passage 500. Optional alpha=0.01 and alpha=0.99 cells "
            "were not extended because neither their raw nor control-adjusted "
            "diagnostics failed.",
            "",
            "## Limitation",
            "",
            "This adaptive report cannot answer the primary H-by-B comparison, the "
            "complete alpha-by-m interaction, or diversity/evenness responses. Those "
            "analyses require the HPC scratch outputs or portable derived endpoint "
            "tables.",
            "",
        ]
    )
    markdown.write_text("\n".join(lines), encoding="utf-8")


def build(
    repository: Path,
    *,
    pdf: Path | None = None,
    markdown: Path | None = None,
    assets: Path | None = None,
    completion: Path | None = None,
) -> list[Path]:
    repository = repository.resolve()
    phase = repository / "experiments/work/trophosome/p01-neutral-feedback"
    derived = phase / "analysis" / DERIVED_DIRECTORY
    decision_path = derived / "adaptive-horizon-decision-g100.json"
    diagnostics_path = derived / "adaptive-horizon-diagnostics-g100.tsv"
    cells_path = phase / "design" / f"{EXPERIMENT_ID}-cells.tsv"
    decision = json.loads(decision_path.read_text(encoding="utf-8"))
    diagnostics = _read_tsv(diagnostics_path)
    cells = _read_tsv(cells_path)
    _validate_inputs(decision, diagnostics, cells)

    pdf = (
        pdf
        or repository / "output/pdf" / f"{EXPERIMENT_ID}-passage100-report.pdf"
    )
    markdown = (
        markdown
        or repository / "docs" / f"{EXPERIMENT_ID}-passage100-report.md"
    )
    assets = (
        assets
        or repository / "docs/figures" / f"{EXPERIMENT_ID}-passage100-report"
    )
    late_figure, stability_figure = _make_figures(diagnostics, assets)
    _build_pdf(pdf, decision, diagnostics, late_figure, stability_figure)
    _build_markdown(markdown, decision, diagnostics, late_figure, stability_figure)

    inputs = [decision_path, diagnostics_path, cells_path]
    outputs = [pdf, markdown, late_figure, stability_figure]
    completion = completion or derived / "passage100-report-completion.json"
    completion.parent.mkdir(parents=True, exist_ok=True)
    completion.write_text(
        json.dumps(
            {
                "complete": True,
                "scope": "passage-100-adaptive-horizon-decision",
                "experiment_id": EXPERIMENT_ID,
                "completed_at": datetime.now(UTC).isoformat(),
                "inputs": [
                    {
                        "path": _reported_path(path, repository),
                        "sha256": _sha256(path),
                    }
                    for path in inputs
                ],
                "outputs": [
                    {
                        "path": _reported_path(path, repository),
                        "sha256": _sha256(path),
                    }
                    for path in outputs
                ],
                "limitation": (
                    "Full H-by-B, alpha-by-m and diversity analyses require derived "
                    "tables from HPC scratch outputs."
                ),
            },
            indent=2,
        )
        + "\n",
        encoding="utf-8",
    )
    print(f"Built self-contained report: {pdf}")
    print(f"Editable companion: {markdown}")
    return [*outputs, completion]


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--repository", type=Path, default=Path(__file__).resolve().parents[1]
    )
    parser.add_argument("--pdf", type=Path, help="optional PDF output path")
    parser.add_argument("--markdown", type=Path, help="optional Markdown output path")
    parser.add_argument("--assets", type=Path, help="optional figure directory")
    parser.add_argument("--completion", type=Path, help="optional audit output path")
    args = parser.parse_args()
    build(
        args.repository,
        pdf=args.pdf.resolve() if args.pdf else None,
        markdown=args.markdown.resolve() if args.markdown else None,
        assets=args.assets.resolve() if args.assets else None,
        completion=args.completion.resolve() if args.completion else None,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
