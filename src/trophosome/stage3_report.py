"""Biologist-readable PDF and Markdown for the audited first mapping wave."""

from __future__ import annotations

import json
from datetime import date
from pathlib import Path

from trophosome.second_pilot_report import _configure_matplotlib, _read_tsv

LABELS = {
    "D0": "Strain richness",
    "D1": "Effective diversity (Hill q=1)",
    "D2": "Dominant-strain diversity (Hill q=2)",
    "evenness": "Evenness",
    "TV": "Composition distance from start",
    "shannon": "Shannon index",
    "simpson": "Simpson index (1 - concentration)",
}


def generate_report(
    *, analysis: Path, design: Path, pdf: Path, markdown: Path, assets: Path
) -> list[Path]:
    summary = json.loads((analysis / "analysis-summary.json").read_text())
    if (
        summary.get("audit_status") != "PASS"
        or summary.get("populations") != 288
        or summary.get("primary_populations_including_reused_control") != 300
        or summary.get("primary_endpoint") != 100
    ):
        raise ValueError(
            "a passed 288-new-run / 300-primary-population audit is required"
        )
    matrix = _read_tsv(design)
    new_cells = [r for r in matrix if r["source_role"] == "new-grid-cell"]
    by_cell = {r["cell_id"]: r for r in matrix}
    statistics = _read_tsv(analysis / "cell-summaries.tsv")
    trajectories = _read_tsv(analysis / "environment-trajectories.tsv")
    contrasts = _read_tsv(analysis / "paired-contrasts.tsv")
    classifications = _read_tsv(analysis / "classifications.tsv")
    references = _read_tsv(analysis / "reused-reference-endpoints.tsv")
    reference_statistics = _read_tsv(analysis / "supplementary-reference-summaries.tsv")
    reference_cells = {r["cell_id"]: r for r in references}
    interactions = _read_tsv(analysis / "tv-interactions.tsv")
    _configure_matplotlib()
    import matplotlib.pyplot as plt
    import numpy as np
    from reportlab.lib import colors
    from reportlab.lib.pagesizes import A4
    from reportlab.lib.styles import getSampleStyleSheet
    from reportlab.platypus import (
        Image,
        PageBreak,
        Paragraph,
        SimpleDocTemplate,
        Spacer,
        Table,
        TableStyle,
    )

    for directory in (pdf.parent, markdown.parent, assets):
        directory.mkdir(parents=True, exist_ok=True)
    plt.rcParams.update(
        {
            "font.size": 9,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "svg.fonttype": "none",
        }
    )
    figure_paths = []

    def save(fig, stem):
        for extension in ("png", "svg"):
            path = assets / f"{stem}.{extension}"
            fig.savefig(path, dpi=180, facecolor="white", bbox_inches="tight")
            figure_paths.append(path)
        plt.close(fig)
        return assets / f"{stem}.png"

    fig, axes = plt.subplots(1, 2, figsize=(9.8, 4.4), constrained_layout=True)
    control = next(
        r for r in statistics if r["cell"] == "c0021" and r["response"] == "TV"
    )
    for axis, u in zip(axes, ("0", "1e-10"), strict=True):
        for h, color, marker in (
            (100, "#31688e", "o"),
            (1000, "#35a779", "s"),
            (10000, "#b16d21", "^"),
        ):
            cells = [c for c in new_cells if int(c["H"]) == h and c["u"] == u]
            rows = [
                next(
                    r
                    for r in statistics
                    if r["cell_id"] == c["cell_id"] and r["response"] == "TV"
                )
                for c in cells
            ]
            mean = np.array([float(r["mean"]) for r in rows])
            axis.errorbar(
                [float(c["alpha"]) for c in cells],
                mean,
                yerr=[
                    mean - [float(r["ci90_low"]) for r in rows],
                    [float(r["ci90_high"]) for r in rows] - mean,
                ],
                marker=marker,
                color=color,
                capsize=3,
                label=f"H={h:,}",
            )
        for r in reference_statistics:
            c = reference_cells[r["cell_id"]]
            if r["response"] == "TV" and c["u"] == u:
                x, y = float(c["alpha"]), float(r["mean"])
                axis.errorbar(
                    x,
                    y,
                    yerr=[[y - float(r["ci90_low"])]],
                    fmt="D",
                    mfc="none",
                    color="#777777",
                    capsize=2,
                )
                axis.annotate(
                    c["cell"],
                    (x, y),
                    xytext=(3, 5),
                    textcoords="offset points",
                    fontsize=7,
                )
        axis.axhline(float(control["mean"]), color="#777777", linestyle=":", lw=1)
        axis.set_xscale("log")
        axis.set_xticks([0.001, 0.01, 0.1, 0.99], ["0.001", "0.01", "0.1", "0.99"])
        axis.set_title("Mutation off" if u == "0" else "Mutation u=10^-10")
        axis.set_xlabel("Feedback alpha (targets labelled; realized positions)")
        axis.set_ylabel("TV at passage 100")
        axis.legend(fontsize=8)
        axis.grid(alpha=0.15)
    map_figure = save(fig, "host-feedback-tv")

    fig, axes = plt.subplots(4, 2, figsize=(9.8, 10.0), constrained_layout=True)
    for axis, response in zip(axes.flat, LABELS, strict=False):
        selected = [r for r in statistics if r["response"] == response]
        x = np.arange(len(selected))
        means = np.array([float(r["mean"]) for r in selected])
        lows = np.array([float(r["ci90_low"]) for r in selected])
        highs = np.array([float(r["ci90_high"]) for r in selected])
        axis.errorbar(
            x,
            means,
            yerr=[means - lows, highs - means],
            fmt="none",
            color="#64748b",
            capsize=2,
        )
        axis.scatter(
            x,
            means,
            c=[
                "#777777"
                if r["cell"] == "c0021"
                else "#21718c"
                if by_cell[r["cell_id"]]["u"] == "0"
                else "#d17a2f"
                for r in selected
            ],
            s=18,
            zorder=3,
        )
        axis.set_xticks(x, [r["cell"] for r in selected], rotation=90, fontsize=7)
        axis.set_title(LABELS[response])
        axis.grid(axis="y", alpha=0.15)
    axes.flat[-1].axis("off")
    axes.flat[-1].text(
        0.05,
        0.85,
        "Each point is a mean of 12 populations.\nBars: 90% confidence interval for that mean.\nBlue: mutation-free; orange: mutation-enabled.\nGrey: reused no-return control.\nAll outcomes are measured at passage 100.\n\nCell order is not a single-factor dose series.",
        va="top",
        linespacing=1.7,
    )
    endpoint_figure = save(fig, "passage100-diversity")
    trajectory_figures = []
    for stratum, cells in (
        ("mutation-free", [r for r in new_cells if r["u"] == "0"]),
        ("mutation-enabled", [r for r in new_cells if r["u"] != "0"]),
    ):
        columns = 4
        fig, axes = plt.subplots(
            len(cells) // columns, columns, figsize=(9.8, 7.0), constrained_layout=True
        )
        for axis, cell in zip(axes.flat, cells, strict=True):
            selected = [r for r in trajectories if r["cell_id"] == cell["cell_id"]]
            for seed in sorted({r["seed_block_id"] for r in selected}):
                rows = sorted(
                    (r for r in selected if r["seed_block_id"] == seed),
                    key=lambda r: int(r["generation"]),
                )
                axis.plot(
                    [int(r["generation"]) for r in rows],
                    [float(r["TV"]) for r in rows],
                    lw=0.65,
                    alpha=0.3,
                    color="#21718c",
                )
            generations = sorted({int(r["generation"]) for r in selected})
            means = [
                np.median(
                    [float(r["TV"]) for r in selected if int(r["generation"]) == g]
                )
                for g in generations
            ]
            axis.plot(generations, means, lw=1.6, color="#172f46")
            axis.set_title(
                f"{cell['cell']} | H={int(cell['H']):,}\nalpha target={cell['alpha_target']}",
                fontsize=8,
            )
            axis.set_xlabel("Host passage", fontsize=8)
            axis.set_ylabel("TV", fontsize=8)
            axis.axvspan(51, 100, color="#b3c9d5", alpha=0.12)
            axis.grid(alpha=0.15)
        trajectory_figures.append(save(fig, f"trajectories-{stratum}"))
    chosen = [r for r in contrasts if r["contrast"].startswith("pooling-")]
    fig, axes = plt.subplots(1, 2, figsize=(9.8, 6.0), constrained_layout=True)
    for axis, response in zip(axes, ("D1", "TV"), strict=True):
        rows = [r for r in chosen if r["response"] == response]
        scale = 100 if response == "D1" else 1
        means = np.array([float(r["mean"]) for r in rows]) * scale
        lows = np.array([float(r["ci90_low"]) for r in rows]) * scale
        highs = np.array([float(r["ci90_high"]) for r in rows]) * scale
        axis.errorbar(
            means,
            range(len(rows)),
            xerr=[means - lows, highs - means],
            fmt="o",
            color="#21718c",
            ms=4,
            capsize=3,
        )
        axis.axvspan(
            -5 if response == "D1" else -0.05,
            5 if response == "D1" else 0.05,
            color="#e8eff3",
        )
        axis.axvline(0, color="#64748b", lw=0.8)
        axis.set_yticks(
            range(len(rows)),
            [f"{r['treatment'][-5:]} - {r['reference'][-5:]}" for r in rows],
            fontsize=8,
        )
        axis.invert_yaxis()
        axis.set_title(LABELS[response])
        axis.set_xlabel(
            "Relative change (%)" if response == "D1" else "Absolute difference"
        )
    contrast_figure = save(fig, "paired-comparisons")
    fig, axes = plt.subplots(1, 3, figsize=(9.8, 4.5), constrained_layout=True)
    for axis, u in zip(axes, ("0", "1e-10", "mutation-modification"), strict=True):
        rows = [r for r in interactions if r["mutation"] == u]
        mean = np.array([float(r["mean"]) for r in rows])
        axis.errorbar(
            mean,
            range(len(rows)),
            xerr=[
                mean - [float(r["ci90_low"]) for r in rows],
                [float(r["ci90_high"]) for r in rows] - mean,
            ],
            fmt="o",
            capsize=3,
            color="#21718c",
        )
        axis.axvline(0, color="#777777", lw=1)
        axis.set_yticks(
            range(len(rows)),
            [
                f"H {r['H_low']} to {r['H_high']}\na {r['alpha_low_target']} to {r['alpha_high_target']}"
                for r in rows
            ],
            fontsize=7,
        )
        axis.invert_yaxis()
        axis.set_title(
            f"Mutation u={u}"
            if u != "mutation-modification"
            else "Modification by mutation",
            fontsize=9,
        )
        axis.set_xlabel("Paired difference-in-differences", fontsize=8)
    interaction_figure = save(fig, "paired-tv-interactions")

    overview = (
        "Stage 3 part one contains 24 new conditions and 12 matched seed blocks "
        "per condition: 288 new populations. One shared c0021 control contributes "
        "12 reused populations, making 25 primary conditions and 300 primary populations. "
        "All outcomes are measured at passage 100. Five off-grid Stage 2 conditions "
        "(60 populations) provide supplementary reference points, not grid substitutes. "
        "Selection is off in both habitats. The environment receives host returns and "
        "10% replacement from a fixed regional source at each passage. "
        "These are finite-time outcomes, not evidence that every condition reached equilibrium."
    )
    interpretation = (
        "Richness counts distinct strain labels. Hill q=1 diversity is the number of equally "
        "common strains that would give the same Shannon diversity; Hill q=2 gives more "
        "weight to common strains. Evenness describes how equally abundance is shared. "
        "TV is the fraction of abundance that must be reassigned to recover the starting "
        "composition (0 means identical, 1 means entirely different)."
    )
    if summary.get("data_origin") == "synthetic-test-fixture":
        overview = "SYNTHETIC SOFTWARE TEST: not biological results. " + overview
    caution = (
        "Intervals use populations, not passages or individual hosts, as replicates. "
        "They are exploratory 90% Student-t intervals, without adjustment for multiple "
        "comparisons. The shaded bands mark working biological margins, not null-hypothesis "
        "significance thresholds. Relative diversity contrasts are calculated within each "
        "matched seed before averaging. Stage 2 references are frozen passage-100 results "
        "from sb0001-sb0012, not additional simulations. A TV contrast compares distances "
        "from the start; it does not compare unrelated mutation IDs across simulations."
    )
    lines = [
        "# Phase 1: Stage 3 part one",
        "",
        f"Report date: {date.today().isoformat()}. Audit: `PASS`.",
        "",
        overview,
        "",
        interpretation,
        "",
        "## Frozen design",
        "",
        "| Cell | # hosts H | Escape f | Mutation u | Alpha target | Alpha realized | Source / escape range |",
        "|---|---:|---:|---:|---:|---:|---|",
    ]
    for cell in matrix:
        lines.append(
            f"| {cell['cell']} | {int(cell['H']):,} | {cell['f']} | {cell['u']} | {cell['alpha_target']} | {float(cell['alpha']):.9g} | {cell['source_role']}; {cell['escape_range']} |"
        )
    lines += [
        "",
        "Alpha is host feedback before regional migration. Total return is exactly matched across H at each target; realized alpha reflects whole-cell rounding. At alpha=0.99, f is 0.99, 0.099 or 0.0099. Extended escape-range cells are labelled.",
        "",
        "## Results and uncertainty",
        "",
        caution,
        "",
    ]
    for caption, path in zip(
        (
            "Host-feedback TV response (grey diamonds: off-grid references; dotted line: shared control)",
            "Passage-100 diversity",
            "Mutation-free trajectories",
            "Mutation-enabled trajectories",
            "Matched comparisons",
            "H-by-feedback TV interactions and their modification by mutation",
        ),
        (
            map_figure,
            endpoint_figure,
            *trajectory_figures,
            contrast_figure,
            interaction_figure,
        ),
        strict=True,
    ):
        lines += [
            f"### {caption}",
            "",
            f"![{caption}](figures/{assets.name}/{path.name})",
            "",
        ]
    lines += [
        "## Exploratory classifications",
        "",
        "| Cell | Relative to the no-return control |",
        "|---|---|",
    ]
    lines += [
        f"| {r['cell']} | {r['classification'].replace('_', ' ')} |"
        for r in classifications
    ]
    lines += [
        "",
        "Mixed or unresolved does not mean no effect. Negligible means equivalence "
        "only for the five measured statistics, not identity of entire communities.",
        "",
        "## What to do next",
        "",
        f"{summary['more_seeds_flags']} cell-response means have interval half-widths larger than their working precision margin.",
        "",
        "Inspect these cells for additional seeds. Fit response surfaces separately for mutation-free and mutation-enabled cells, "
        "using whole-cell held-out validation before selecting new parameter combinations. Do not launch adaptive additions automatically. "
        "Tail summaries for passages 51-100 describe ongoing fluctuations. The paired mean-TV change from passages 51-75 to 76-100 flags cells for time-horizon review; this is not a complete stationarity or equilibrium test. "
        "Keep raw outputs until analysis and archiving have been checked.",
        "",
    ]
    markdown.write_text("\n".join(lines), encoding="utf-8")
    styles = getSampleStyleSheet()
    styles["BodyText"].fontSize = 9
    styles["BodyText"].leading = 13
    story = []

    def paragraph(text, style="BodyText"):
        story.append(Paragraph(text, styles[style]))
        story.append(Spacer(1, 7))

    def figure(path, width=490):
        from PIL import Image as PillowImage

        with PillowImage.open(path) as picture:
            w, h = picture.size
        story.append(Image(str(path), width=width, height=width * h / w))

    paragraph("Phase 1 | Stage 3 part one", "Title")
    paragraph(
        f"Exploratory results at passage 100 | {date.today().isoformat()} | Audit: PASS"
    )
    paragraph(overview)
    paragraph("The frozen design", "Heading2")
    table = [["Cell", "# hosts H", "Escape f", "Mutation u", "Alpha target"]] + [
        [c["cell"], f"{int(c['H']):,}", c["f"], c["u"], c["alpha_target"]]
        for c in matrix
    ]
    t = Table(table, colWidths=[62, 80, 90, 100, 130], repeatRows=1)
    t.setStyle(
        TableStyle(
            [
                ("BACKGROUND", (0, 0), (-1, 0), colors.HexColor("#173b4f")),
                ("TEXTCOLOR", (0, 0), (-1, 0), colors.white),
                (
                    "ROWBACKGROUNDS",
                    (0, 1),
                    (-1, -1),
                    [colors.white, colors.HexColor("#edf3f5")],
                ),
                ("FONTSIZE", (0, 0), (-1, -1), 8),
                ("BOTTOMPADDING", (0, 0), (-1, -1), 3),
            ]
        )
    )
    story.append(t)
    paragraph(
        "Alpha is feedback before migration. B=10; K=N_E=10^9; growth factor=1.2; "
        "500 steady bacterial generations; frozen 100-strain starting and regional population."
    )
    story.append(PageBreak())
    paragraph("Host abundance and feedback", "Title")
    figure(map_figure)
    paragraph(
        "Means and 90% intervals at passage 100. Grey open diamonds are supplementary Stage 2 references at their actual feedback values; the dotted line is the shared no-return control. Connecting lines guide the eye, not fitted predictions."
    )
    paragraph(caution)
    paragraph(
        "Exact total return is matched across H at each feedback target. Whole-cell rounding explains the small target/realized difference recorded in the TSV. At alpha=0.99, H=100 and H=1,000 require 99% and 9.9% release: extreme-range tests."
    )
    story.append(PageBreak())
    paragraph("How host passage changes diversity", "Title")
    paragraph(interpretation)
    figure(endpoint_figure)
    paragraph(
        "The full classification table and all five responses are included in the editable report and derived tables."
    )
    for stratum, path in zip(
        ("Mutation-free", "Mutation-enabled"), trajectory_figures, strict=True
    ):
        story.append(PageBreak())
        paragraph(f"{stratum} trajectories", "Title")
        paragraph(
            "Thin lines show the 12 individual populations; the dark line is their median. "
            "Each panel has its own vertical scale. A flat-looking mean does not establish equilibrium."
        )
        figure(path)
        paragraph(
            "Complete diversity, composition, turnover and ancestral-lineage metrics are retained for all 101 environmental states of each new population. Reused references are derived summaries; their raw records remain in Stage 2."
        )
    story.append(PageBreak())
    paragraph("Matched biological comparisons", "Title")
    paragraph(caution)
    figure(contrast_figure)
    paragraph("Next decisions", "Heading2")
    paragraph(
        f"{summary['more_seeds_flags']} cell-response means exceed the working precision margin. "
        "Consider more seeds where stochastic uncertainty is large; consider new cells where whole-cell "
        "held-out prediction is poor. These decisions remain a separate review, not an automatic launch."
    )
    paragraph(
        "No raw results or checkpoints are deleted by reporting. Archive audited outputs before changing retention."
    )
    story.append(PageBreak())
    paragraph("Does H modify the feedback effect?", "Title")
    figure(interaction_figure)
    interaction_text = "Each interaction subtracts the feedback effect at smaller H from the feedback effect at larger H, within the same seed. Negative values in the first two panels mean feedback changes TV less at larger H. The third panel subtracts the mutation-free interaction from the mutation-enabled interaction. These exploratory 90% intervals have no assigned biological equivalence threshold."
    paragraph(interaction_text)
    glossary = "For strain frequencies p_i: Shannon = -sum(p_i ln p_i), using natural logarithms. Hill D1 = exp(Shannon). Simpson = 1 - sum(p_i squared), the probability that two independent draws have different labels. Hill D2 = 1 / sum(p_i squared). Pielou evenness = Shannon / ln(D0) for D0 > 1. TV = half the sum of absolute frequency differences across the union of labels; absent strains have frequency zero. D0 counts labels; ancestral diversity merges mutant descendants with their founding strain."
    paragraph("Appendix | Diversity glossary", "Heading2")
    paragraph(glossary)
    with markdown.open("a", encoding="utf-8") as handle:
        handle.write(
            "\n"
            + interaction_text
            + "\n\n## Appendix: diversity glossary\n\n"
            + glossary
            + "\n"
        )

    def footer(canvas, document):
        canvas.setFont("Helvetica", 8)
        canvas.setFillColor(colors.HexColor("#64748b"))
        canvas.drawString(
            45, 25, "trophosome | Stage 3 wave 1 | exploratory, not confirmatory"
        )
        canvas.drawRightString(A4[0] - 45, 25, str(document.page))

    SimpleDocTemplate(
        str(pdf),
        pagesize=A4,
        rightMargin=45,
        leftMargin=45,
        topMargin=38,
        bottomMargin=42,
    ).build(story, onFirstPage=footer, onLaterPages=footer)
    return [pdf, markdown, *figure_paths]
