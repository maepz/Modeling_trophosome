"""Build the Phase 1 second-pilot report from audited derived tables."""

from __future__ import annotations

import csv
import html
import json
import math
import os
import tempfile
from collections import defaultdict
from collections.abc import Callable
from dataclasses import dataclass
from datetime import date
from pathlib import Path
from typing import Any


@dataclass(frozen=True)
class SecondPilotReportArtifacts:
    """Files created by :func:`generate_second_pilot_report`."""

    pdf: Path
    markdown: Path
    assets: Path


CELL_COLOURS = {
    "p01-s02-c0021": "#7B8794",
    "p01-s02-c0022": "#226C8A",
    "p01-s02-c0023": "#2B9C6A",
    "p01-s02-c0024": "#7A5195",
    "p01-s02-c0025": "#E3A23B",
    "p01-s02-c0026": "#C94A53",
}


def _configure_matplotlib() -> None:
    """Select a headless backend and a writable, reusable font cache."""

    cache = Path(tempfile.gettempdir()) / f"trophosome-matplotlib-{os.getuid()}"
    cache.mkdir(parents=True, exist_ok=True)
    os.environ.setdefault("MPLCONFIGDIR", str(cache))
    os.environ.setdefault("XDG_CACHE_HOME", str(cache / "xdg"))
    import matplotlib

    matplotlib.use("Agg", force=True)


def _read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def _short_cell(cell_id: str) -> str:
    return cell_id.rsplit("-", 1)[-1]


def _float(value: Any, default: float = math.nan) -> float:
    if value in (None, ""):
        return default
    return float(value)


def _true(value: Any) -> bool:
    return str(value).strip().lower() == "true"


def _label(row: dict[str, str]) -> str:
    return f"{_short_cell(row['cell_id'])}: {row['label']}"


def _cell_stationarity_text(cell_id: str, stationarity: list[dict[str, str]]) -> str:
    rows = [row for row in stationarity if row["cell_id"] == cell_id]
    generations = [
        int(row["persistent_stationarity_generation"])
        for row in rows
        if row.get("persistent_stationarity_generation", "") not in ("", None)
    ]
    if len(rows) != 5 or len(generations) != 5:
        return "not by 250"
    return f"by {max(generations)}"


def _load(
    analysis: Path, design_path: Path
) -> tuple[
    dict[str, Any],
    list[dict[str, str]],
    list[dict[str, str]],
    list[dict[str, str]],
    list[dict[str, str]],
    list[dict[str, str]],
]:
    required = {
        "summary": analysis / "analysis-summary.json",
        "trajectories": analysis / "environment-trajectories.tsv",
        "windows": analysis / "run-window-summaries.tsv",
        "stationarity": analysis / "stationarity-screen.tsv",
        "precision": analysis / "precision-recommendations.tsv",
        "design": design_path,
    }
    missing = [str(path) for path in required.values() if not path.is_file()]
    if missing:
        raise FileNotFoundError(
            "missing second-pilot report inputs: " + ", ".join(missing)
        )
    summary = json.loads(required["summary"].read_text(encoding="utf-8"))
    if summary.get("audit_status") != "PASS":
        raise ValueError("analysis audit is not PASS")
    design = _read_tsv(required["design"])
    if len(design) != 6:
        raise ValueError("second-pilot report requires six sentinel cells")
    return (
        summary,
        design,
        _read_tsv(required["trajectories"]),
        _read_tsv(required["windows"]),
        _read_tsv(required["stationarity"]),
        _read_tsv(required["precision"]),
    )


def _plot_trajectories(
    path: Path,
    design: list[dict[str, str]],
    trajectories: list[dict[str, str]],
) -> None:
    _configure_matplotlib()
    import matplotlib.pyplot as plt
    import numpy as np

    grouped: dict[tuple[str, int, str], list[float]] = defaultdict(list)
    for row in trajectories:
        for response in ("D1", "D2", "TV"):
            grouped[(row["cell_id"], int(row["generation"]), response)].append(
                float(row[response])
            )

    figure, axes = plt.subplots(3, 1, figsize=(9.2, 8.4), sharex=True)
    for axis, response, title in zip(
        axes,
        ("D1", "D2", "TV"),
        (
            "Common-strain diversity (Hill D1)",
            "Dominant-strain diversity (Hill D2)",
            "Compositional departure from the initial reservoir (TV)",
        ),
        strict=True,
    ):
        for cell in design:
            cell_id = cell["cell_id"]
            generations = sorted(
                generation
                for candidate, generation, metric in grouped
                if candidate == cell_id and metric == response
            )
            values = [
                grouped[(cell_id, generation, response)] for generation in generations
            ]
            median = np.asarray([np.median(item) for item in values])
            lower = np.asarray([np.quantile(item, 0.10) for item in values])
            upper = np.asarray([np.quantile(item, 0.90) for item in values])
            colour = CELL_COLOURS[cell_id]
            axis.plot(
                generations,
                median,
                color=colour,
                linewidth=1.8,
                label=_short_cell(cell_id),
            )
            axis.fill_between(generations, lower, upper, color=colour, alpha=0.12)
        axis.set_title(title, loc="left", fontsize=10.5, fontweight="bold")
        axis.grid(axis="y", color="#DCE3E8", linewidth=0.7)
        axis.spines[["top", "right"]].set_visible(False)
        axis.set_ylabel(response)
    axes[0].legend(ncol=6, frameon=False, loc="upper center", fontsize=8)
    axes[-1].set_xlabel("Host-population passage")
    figure.suptitle(
        "Environmental diversity through 250 repeated host passages",
        x=0.07,
        ha="left",
        fontsize=14,
        fontweight="bold",
        color="#173B4F",
    )
    figure.tight_layout(rect=(0, 0, 1, 0.97))
    path.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(path, dpi=200, bbox_inches="tight", facecolor="white")
    plt.close(figure)


def _plot_stationarity(
    path: Path,
    design: list[dict[str, str]],
    stationarity: list[dict[str, str]],
) -> None:
    _configure_matplotlib()
    import matplotlib.pyplot as plt
    import numpy as np
    from matplotlib.colors import ListedColormap

    responses = ("D0", "D1", "D2", "evenness", "TV")
    index = {(row["cell_id"], row["response"]): row for row in stationarity}
    matrix = np.asarray(
        [
            [
                1
                if _true(index[(cell["cell_id"], response)]["stationarity_screen_pass"])
                else 0
                for cell in design
            ]
            for response in responses
        ]
    )
    figure, axis = plt.subplots(figsize=(8.8, 4.3))
    axis.imshow(matrix, cmap=ListedColormap(["#D95C59", "#3A9D6F"]), vmin=0, vmax=1)
    for y, response in enumerate(responses):
        for x, cell in enumerate(design):
            row = index[(cell["cell_id"], response)]
            symbol = "PASS" if matrix[y, x] else "CHECK"
            axis.text(
                x,
                y - 0.08,
                symbol,
                ha="center",
                va="center",
                color="white",
                fontsize=8,
                fontweight="bold",
            )
            axis.text(
                x,
                y + 0.24,
                f"Rhat {_float(row['rank_normalized_split_rhat']):.2f}\nESS {_float(row['approximate_combined_ess']):.0f}",
                ha="center",
                va="center",
                color="white",
                fontsize=6.3,
            )
    axis.set_xticks(range(len(design)), [_short_cell(row["cell_id"]) for row in design])
    axis.set_yticks(range(len(responses)), responses)
    axis.tick_params(length=0)
    axis.set_title(
        "Stationarity screen by environmental response",
        loc="left",
        fontsize=13,
        fontweight="bold",
        color="#173B4F",
        pad=14,
    )
    axis.set_xlabel("Sentinel cell")
    for spine in axis.spines.values():
        spine.set_visible(False)
    figure.tight_layout()
    path.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(path, dpi=200, bbox_inches="tight", facecolor="white")
    plt.close(figure)


def _plot_fluctuations(
    path: Path,
    design: list[dict[str, str]],
    windows: list[dict[str, str]],
) -> None:
    _configure_matplotlib()
    import matplotlib.pyplot as plt
    import numpy as np

    responses = ("D0", "D1", "D2", "evenness", "TV")
    grouped: dict[tuple[str, str], list[float]] = defaultdict(list)
    for row in windows:
        if row["window"] == "4" and math.isfinite(_float(row["cv"])):
            grouped[(row["cell_id"], row["response"])].append(float(row["cv"]))
    matrix = np.asarray(
        [
            [np.median(grouped[(cell["cell_id"], response)]) for cell in design]
            for response in responses
        ]
    )
    figure, axis = plt.subplots(figsize=(8.8, 4.2))
    image = axis.imshow(matrix, cmap="Blues", aspect="auto")
    ceiling = max(float(np.nanmax(matrix)), 1e-12)
    for y in range(len(responses)):
        for x in range(len(design)):
            colour = "white" if matrix[y, x] > 0.55 * ceiling else "#173B4F"
            axis.text(
                x,
                y,
                f"{100 * matrix[y, x]:.2f}%",
                ha="center",
                va="center",
                fontsize=8,
                color=colour,
            )
    axis.set_xticks(range(len(design)), [_short_cell(row["cell_id"]) for row in design])
    axis.set_yticks(range(len(responses)), responses)
    axis.tick_params(length=0)
    axis.set_xlabel("Sentinel cell")
    axis.set_title(
        "Continuing fluctuation in the final diagnostic window",
        loc="left",
        fontsize=13,
        fontweight="bold",
        color="#173B4F",
        pad=14,
    )
    colourbar = figure.colorbar(image, ax=axis, shrink=0.82)
    colourbar.set_label("Median coefficient of variation")
    for spine in axis.spines.values():
        spine.set_visible(False)
    figure.tight_layout()
    path.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(path, dpi=200, bbox_inches="tight", facecolor="white")
    plt.close(figure)


def _plot_precision(path: Path, precision: list[dict[str, str]]) -> None:
    _configure_matplotlib()
    import matplotlib.pyplot as plt

    contrasts = list(dict.fromkeys(row["contrast"] for row in precision))
    responses = ("D1", "D2", "TV")
    index = {(row["contrast"], row["response"]): row for row in precision}
    figure, axes = plt.subplots(1, 3, figsize=(10.5, 4.4), sharey=True)
    for axis, response in zip(axes, responses, strict=True):
        values = [
            int(index[(contrast, response)]["recommended_replicates"])
            for contrast in contrasts
        ]
        colours = [
            "#C94A53"
            if _true(index[(contrast, response)]["exceeds_predeclared_cap"])
            else "#2C7DA0"
            for contrast in contrasts
        ]
        axis.barh(range(len(contrasts)), values, color=colours)
        axis.axvline(20, color="#6B7785", linestyle="--", linewidth=1)
        axis.set_title(response, fontsize=10.5, fontweight="bold")
        axis.set_xlabel("Recommended matched replicates")
        axis.grid(axis="x", color="#E0E6EA", linewidth=0.7)
        axis.spines[["top", "right", "left"]].set_visible(False)
        axis.tick_params(axis="y", length=0)
        for y, value in enumerate(values):
            axis.text(value + 1, y, str(value), va="center", fontsize=8)
    axes[0].set_yticks(
        range(len(contrasts)), [value.replace("-", " ") for value in contrasts]
    )
    figure.suptitle(
        "Replicates recommended for the confirmatory experiment",
        x=0.08,
        ha="left",
        fontsize=13,
        fontweight="bold",
        color="#173B4F",
    )
    figure.tight_layout(rect=(0, 0, 1, 0.93))
    path.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(path, dpi=200, bbox_inches="tight", facecolor="white")
    plt.close(figure)


def _write_markdown(
    path: Path,
    title: str,
    report_date: str,
    summary: dict[str, Any],
    design: list[dict[str, str]],
    stationarity: list[dict[str, str]],
    precision: list[dict[str, str]],
    figures: dict[str, Path],
) -> None:
    relative_figures = {
        key: figure.relative_to(path.parent).as_posix()
        for key, figure in figures.items()
    }
    cell_pass = {
        cell["cell_id"]: sum(
            _true(row["stationarity_screen_pass"])
            for row in stationarity
            if row["cell_id"] == cell["cell_id"]
        )
        for cell in design
    }
    maximum_precision = max(int(row["recommended_replicates"]) for row in precision)
    lines = [
        f"# {title}",
        "",
        f"**Report date:** {report_date}  ",
        f"**Analysis audit:** `{summary['audit_status']}`  ",
        f"**Populations:** {summary['populations']} ({summary['cells']} cells x {summary['seed_blocks']} matched seed blocks)  ",
        f"**Host passages:** {summary['host_generations']}",
        "",
        "## What this pilot asks",
        "",
        (
            "This second pilot asks whether environmental diversity appears statistically "
            "stable late in a 250-passage run, how much it continues to fluctuate, and how "
            "many independent matched replicates a later confirmatory experiment may need."
        ),
        "",
        "## Main result",
        "",
        (
            f"{summary['stationarity_responses_passing']} of "
            f"{summary['stationarity_responses_total']} cell-response screens passed; "
            f"{summary['cells_passing_all_stationarity_responses']} of 6 cells passed all five responses."
        ),
        "",
        (
            "Passing this screen supports late-run stationarity, not definitive equilibrium. "
            "The experiment starts all runs from the same environmental composition, so it "
            "does not test convergence from different starting communities."
        ),
        "",
        "## Sentinel conditions",
        "",
        "| Cell | Biological role | Hosts | Escape fraction | Total return | Mutation rate | Screen passes | Persistent screen |",
        "|---|---|---:|---:|---:|---:|---:|---:|",
    ]
    for row in design:
        lines.append(
            f"| {_short_cell(row['cell_id'])} | {row['sentinel_role']} | {int(row['H']):,} | "
            f"{float(row['f']):.3g} | {int(row['R']):,} | {float(row['u']):.3g} | "
            f"{cell_pass[row['cell_id']]}/5 | "
            f"{_cell_stationarity_text(row['cell_id'], stationarity)} |"
        )
    lines.extend(
        [
            "",
            f"![Environmental trajectories]({relative_figures['trajectories']})",
            "",
            (
                "**Figure 1.** Lines are medians across 12 independent populations; shaded "
                "bands span the 10th to 90th percentiles. TV is total-variation distance from "
                "the initial environmental composition."
            ),
            "",
            f"![Stationarity screen]({relative_figures['stationarity']})",
            "",
            (
                "**Figure 2.** A response passes only when both overlapping late-window "
                "assessments satisfy the predeclared equivalence limits, rank-normalized split "
                "R-hat is below 1.05, and approximate combined ESS is at least 400."
            ),
            "",
            f"![Continuing fluctuations]({relative_figures['fluctuations']})",
            "",
            (
                "**Figure 3.** Median within-run coefficient of variation in the last diagnostic "
                "window. Larger values mean more continuing fluctuation around the late-run level."
            ),
            "",
            "## Precision for the confirmatory experiment",
            "",
            (
                f"The largest recommendation is {maximum_precision} matched replicates. "
                "Recommendations target a 95% confidence-interval half-width of 0.05, use a "
                "minimum of 20 replicates, increase in batches of 8, and are capped at 100."
            ),
            "",
            f"![Precision recommendations]({relative_figures['precision']})",
            "",
            "**Figure 4.** Dashed lines show the minimum of 20 matched replicates.",
            "",
            "## Quality control and computational resources",
            "",
            f"- Analysis audit: `{summary['audit_status']}`.",
            f"- Summed elapsed time: {summary['resources']['summed_elapsed_hours']:.2f} hours.",
            f"- Output volume: {summary['resources']['summed_output_gib']:.2f} GiB.",
            f"- Highest recorded process-tree memory: {summary['resources']['maximum_peak_rss_mib']:.1f} MiB.",
            "- Every run was checked for 251 environmental states, constant reservoir size, migration counts, committed completion metadata, configuration hashes, and final-state checksums.",
            "",
            "## Interpretation limits",
            "",
        ]
    )
    lines.extend(f"- {limitation}" for limitation in summary["limitations"])
    lines.extend(
        [
            "",
            "## Short glossary",
            "",
            "- **D0:** strain richness; every detected strain counts equally.",
            "- **D1:** effective number of common strains; sensitive to richness and evenness.",
            "- **D2:** effective number of dominant strains; gives more weight to abundant strains.",
            "- **Evenness:** how similarly abundant the strains are.",
            "- **TV:** total-variation distance from the initial environmental composition; 0 means identical and 1 means no overlap.",
            "- **ESS:** effective sample size after accounting for temporal autocorrelation.",
            "- **R-hat:** agreement among independent replicate trajectories after splitting each trajectory; values close to 1 are preferred.",
            "",
        ]
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines), encoding="utf-8")


def _reportlab_imports() -> dict[str, Any]:
    try:
        from reportlab.lib import colors as report_colors
        from reportlab.lib.enums import TA_CENTER as report_ta_center
        from reportlab.lib.pagesizes import A4 as report_a4
        from reportlab.lib.styles import ParagraphStyle as ReportParagraphStyle
        from reportlab.lib.styles import getSampleStyleSheet as report_styles
        from reportlab.lib.units import mm as report_mm
        from reportlab.platypus import (
            Image as ReportImage,
        )
        from reportlab.platypus import (
            PageBreak as ReportPageBreak,
        )
        from reportlab.platypus import (
            Paragraph as ReportParagraph,
        )
        from reportlab.platypus import (
            SimpleDocTemplate as ReportDocument,
        )
        from reportlab.platypus import (
            Spacer as ReportSpacer,
        )
        from reportlab.platypus import (
            Table as ReportTable,
        )
        from reportlab.platypus import (
            TableStyle as ReportTableStyle,
        )
    except ImportError as exc:  # pragma: no cover - dependency guard
        raise RuntimeError("install the 'report' extra to build the PDF") from exc
    return {
        "colors": report_colors,
        "TA_CENTER": report_ta_center,
        "A4": report_a4,
        "ParagraphStyle": ReportParagraphStyle,
        "getSampleStyleSheet": report_styles,
        "mm": report_mm,
        "Image": ReportImage,
        "PageBreak": ReportPageBreak,
        "Paragraph": ReportParagraph,
        "SimpleDocTemplate": ReportDocument,
        "Spacer": ReportSpacer,
        "Table": ReportTable,
        "TableStyle": ReportTableStyle,
    }


def _write_pdf(
    path: Path,
    title: str,
    report_date: str,
    summary: dict[str, Any],
    design: list[dict[str, str]],
    stationarity: list[dict[str, str]],
    precision: list[dict[str, str]],
    figures: dict[str, Path],
) -> None:
    imports = _reportlab_imports()
    colors = imports["colors"]
    A4 = imports["A4"]
    mm = imports["mm"]
    Paragraph = imports["Paragraph"]
    ParagraphStyle = imports["ParagraphStyle"]
    SimpleDocTemplate = imports["SimpleDocTemplate"]
    Spacer = imports["Spacer"]
    PageBreak = imports["PageBreak"]
    Image = imports["Image"]
    Table = imports["Table"]
    TableStyle = imports["TableStyle"]
    getSampleStyleSheet = imports["getSampleStyleSheet"]
    TA_CENTER = imports["TA_CENTER"]

    path.parent.mkdir(parents=True, exist_ok=True)
    document = SimpleDocTemplate(
        str(path),
        pagesize=A4,
        leftMargin=17 * mm,
        rightMargin=17 * mm,
        topMargin=16 * mm,
        bottomMargin=17 * mm,
        title=title,
        author="Trophosome model reporting workflow",
    )
    base = getSampleStyleSheet()
    styles = {
        "title": ParagraphStyle(
            "second-title",
            parent=base["Title"],
            fontName="Helvetica-Bold",
            fontSize=21,
            leading=25,
            textColor=colors.HexColor("#173B4F"),
            spaceAfter=10,
        ),
        "subtitle": ParagraphStyle(
            "second-subtitle",
            parent=base["Normal"],
            fontSize=9,
            leading=13,
            textColor=colors.HexColor("#61717C"),
            spaceAfter=7,
        ),
        "h1": ParagraphStyle(
            "second-h1",
            parent=base["Heading1"],
            fontName="Helvetica-Bold",
            fontSize=14,
            leading=17,
            textColor=colors.HexColor("#173B4F"),
            spaceBefore=12,
            spaceAfter=7,
            keepWithNext=True,
        ),
        "body": ParagraphStyle(
            "second-body",
            parent=base["BodyText"],
            fontSize=9.1,
            leading=13,
            textColor=colors.HexColor("#263640"),
            spaceAfter=7,
        ),
        "small": ParagraphStyle(
            "second-small",
            parent=base["BodyText"],
            fontSize=7.4,
            leading=9.3,
            textColor=colors.HexColor("#263640"),
        ),
        "caption": ParagraphStyle(
            "second-caption",
            parent=base["BodyText"],
            fontName="Helvetica-Oblique",
            fontSize=7.7,
            leading=10.5,
            textColor=colors.HexColor("#51616B"),
            spaceAfter=9,
        ),
        "callout": ParagraphStyle(
            "second-callout",
            parent=base["BodyText"],
            fontSize=10,
            leading=14,
            textColor=colors.HexColor("#173B4F"),
            backColor=colors.HexColor("#EAF3F6"),
            borderColor=colors.HexColor("#7CB7C9"),
            borderWidth=0.7,
            borderPadding=8,
            spaceAfter=10,
        ),
        "footer": ParagraphStyle(
            "second-footer",
            parent=base["Normal"],
            fontSize=7,
            textColor=colors.HexColor("#6A7780"),
            alignment=TA_CENTER,
        ),
    }

    def page(canvas: Any, doc: Any) -> None:
        canvas.saveState()
        canvas.setStrokeColor(colors.HexColor("#D7E0E5"))
        canvas.line(17 * mm, 13 * mm, A4[0] - 17 * mm, 13 * mm)
        canvas.setFillColor(colors.HexColor("#6A7780"))
        canvas.setFont("Helvetica", 7)
        canvas.drawString(17 * mm, 8.5 * mm, "Trophosome Phase 1 second pilot")
        canvas.drawRightString(A4[0] - 17 * mm, 8.5 * mm, f"Page {doc.page}")
        canvas.restoreState()

    cell_pass = {
        cell["cell_id"]: sum(
            _true(row["stationarity_screen_pass"])
            for row in stationarity
            if row["cell_id"] == cell["cell_id"]
        )
        for cell in design
    }
    maximum_precision = max(int(row["recommended_replicates"]) for row in precision)
    cell_stationarity = "; ".join(
        f"{_short_cell(cell['cell_id'])} {_cell_stationarity_text(cell['cell_id'], stationarity)}"
        for cell in design
    )
    story: list[Any] = [
        Paragraph(html.escape(title), styles["title"]),
        Paragraph(
            f"Equilibrium-and-precision pilot | {report_date} | "
            f"Analysis audit: <b>{summary['audit_status']}</b>",
            styles["subtitle"],
        ),
        Paragraph("What this pilot asks", styles["h1"]),
        Paragraph(
            "Does environmental diversity become statistically stable late in a "
            "250-passage run, how strongly does it continue to fluctuate, and how many "
            "matched replicates should a confirmatory experiment use?",
            styles["body"],
        ),
        Paragraph(
            f"<b>{summary['stationarity_responses_passing']} of "
            f"{summary['stationarity_responses_total']}</b> cell-response screens passed; "
            f"<b>{summary['cells_passing_all_stationarity_responses']} of 6</b> cells passed "
            "all five responses.",
            styles["callout"],
        ),
        Paragraph(
            "Important interpretation: passing supports late-run <b>stationarity</b>, "
            "not definitive equilibrium. These runs all begin from the same environmental "
            "community, so convergence from contrasting starting communities is not tested.",
            styles["body"],
        ),
        Paragraph(
            "Earliest assessment from which all five response screens remained "
            f"satisfied: {html.escape(cell_stationarity)}.",
            styles["body"],
        ),
        Paragraph("The six sentinel conditions", styles["h1"]),
    ]
    table_data: list[list[Any]] = [
        [
            Paragraph("Cell", styles["small"]),
            Paragraph("Biological role", styles["small"]),
            Paragraph("Hosts", styles["small"]),
            Paragraph("Escape", styles["small"]),
            Paragraph("Return", styles["small"]),
            Paragraph("u", styles["small"]),
            Paragraph("Pass", styles["small"]),
        ]
    ]
    for row in design:
        table_data.append(
            [
                _short_cell(row["cell_id"]),
                Paragraph(html.escape(row["sentinel_role"]), styles["small"]),
                f"{int(row['H']):,}",
                f"{float(row['f']):.3g}",
                f"{int(row['R']):.2g}",
                f"{float(row['u']):.1g}",
                f"{cell_pass[row['cell_id']]}/5",
            ]
        )
    table = Table(
        table_data,
        colWidths=[16 * mm, 49 * mm, 18 * mm, 20 * mm, 21 * mm, 18 * mm, 16 * mm],
        repeatRows=1,
    )
    table.setStyle(
        TableStyle(
            [
                ("BACKGROUND", (0, 0), (-1, 0), colors.HexColor("#173B4F")),
                ("TEXTCOLOR", (0, 0), (-1, 0), colors.white),
                ("FONTNAME", (0, 0), (-1, 0), "Helvetica-Bold"),
                ("FONTNAME", (0, 1), (-1, -1), "Helvetica"),
                ("FONTSIZE", (0, 1), (-1, -1), 7.2),
                ("VALIGN", (0, 0), (-1, -1), "MIDDLE"),
                ("GRID", (0, 0), (-1, -1), 0.3, colors.HexColor("#CBD6DC")),
                (
                    "ROWBACKGROUNDS",
                    (0, 1),
                    (-1, -1),
                    [colors.white, colors.HexColor("#F2F6F8")],
                ),
                ("LEFTPADDING", (0, 0), (-1, -1), 4),
                ("RIGHTPADDING", (0, 0), (-1, -1), 4),
                ("TOPPADDING", (0, 0), (-1, -1), 4),
                ("BOTTOMPADDING", (0, 0), (-1, -1), 4),
            ]
        )
    )
    story.extend([table, Spacer(1, 5 * mm), PageBreak()])

    def report_image(key: str, width_mm: float) -> Any:
        image = Image(str(figures[key]))
        ratio = image.imageHeight / image.imageWidth
        image.drawWidth = width_mm * mm
        image.drawHeight = width_mm * mm * ratio
        return image

    figure_sections = (
        (
            "Environmental trajectories",
            "trajectories",
            "Lines are medians across 12 independent populations; shaded bands span "
            "the 10th to 90th percentiles. TV measures departure from the initial "
            "environmental composition.",
        ),
        (
            "Late-run stationarity screen",
            "stationarity",
            "PASS requires both late-window equivalence assessments, rank-normalized "
            "split R-hat below 1.05, and approximate combined ESS of at least 400.",
        ),
        (
            "Continuing fluctuation",
            "fluctuations",
            "The coefficient of variation describes within-run fluctuation during each "
            "cell's final diagnostic window; it is not variation among replicates.",
        ),
        (
            "Precision for confirmation",
            "precision",
            f"The largest recommendation is {maximum_precision} matched replicates. The "
            "target is a confidence-interval half-width of 0.05, with a minimum of 20, "
            "batches of 8, and a cap of 100.",
        ),
    )
    for number, (heading, key, caption) in enumerate(figure_sections, start=1):
        story.extend(
            [
                Paragraph(heading, styles["h1"]),
                report_image(key, 174),
                Paragraph(f"Figure {number}. {caption}", styles["caption"]),
            ]
        )
        if number in {1, 3}:
            story.append(PageBreak())

    resources = summary["resources"]
    story.extend(
        [
            Paragraph("Quality control and resources", styles["h1"]),
            Paragraph(
                f"The audit passed for all {summary['populations']} populations. It checked "
                "251 environmental states per run, constant reservoir capacity, migration "
                "counts, committed completion metadata, configuration hashes, and final-state "
                "checksums.",
                styles["body"],
            ),
            Paragraph(
                f"Summed elapsed time: <b>{resources['summed_elapsed_hours']:.2f} h</b>; "
                f"output: <b>{resources['summed_output_gib']:.2f} GiB</b>; highest process-tree "
                f"memory: <b>{resources['maximum_peak_rss_mib']:.1f} MiB</b>.",
                styles["body"],
            ),
            Paragraph("Interpretation limits", styles["h1"]),
        ]
    )
    for limitation in summary["limitations"]:
        story.append(Paragraph(f"&#8226; {html.escape(limitation)}", styles["body"]))
    story.extend(
        [
            Paragraph("Short glossary", styles["h1"]),
            Paragraph(
                "<b>D0</b> is strain richness. <b>D1</b> is the effective number of common "
                "strains. <b>D2</b> emphasizes dominant strains. <b>Evenness</b> describes "
                "similarity of strain abundances. <b>TV</b> is total-variation distance from "
                "the initial environment. <b>ESS</b> discounts temporally redundant "
                "observations. <b>R-hat</b> tests agreement among split replicate trajectories.",
                styles["body"],
            ),
        ]
    )
    document.build(story, onFirstPage=page, onLaterPages=page)


def generate_second_pilot_report(
    *,
    analysis: Path,
    design: Path,
    output: Path,
    markdown: Path,
    assets: Path,
    title: str = "Phase 1 second pilot: stationarity and precision",
    report_date: str | None = None,
) -> SecondPilotReportArtifacts:
    """Generate figures, an editable Markdown report, and a self-contained PDF."""

    summary, cells, trajectories, windows, stationarity, precision = _load(
        analysis, design
    )
    figures: dict[str, Path] = {
        "trajectories": assets / "environmental-trajectories.png",
        "stationarity": assets / "stationarity-screen.png",
        "fluctuations": assets / "continuing-fluctuations.png",
        "precision": assets / "precision-recommendations.png",
    }
    plotters: tuple[tuple[str, Callable[..., None], tuple[Any, ...]], ...] = (
        ("trajectories", _plot_trajectories, (cells, trajectories)),
        ("stationarity", _plot_stationarity, (cells, stationarity)),
        ("fluctuations", _plot_fluctuations, (cells, windows)),
        ("precision", _plot_precision, (precision,)),
    )
    for key, plotter, arguments in plotters:
        plotter(figures[key], *arguments)
    date_text = report_date or date.today().isoformat()
    _write_markdown(
        markdown,
        title,
        date_text,
        summary,
        cells,
        stationarity,
        precision,
        figures,
    )
    _write_pdf(
        output,
        title,
        date_text,
        summary,
        cells,
        stationarity,
        precision,
        figures,
    )
    return SecondPilotReportArtifacts(output, markdown, assets)
