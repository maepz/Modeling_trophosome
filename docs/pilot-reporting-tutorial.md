# Tutorial: generate a biological pilot report

The `trophosome report` command converts standardized pilot summary tables into
a concise, self-contained PDF. It also writes an editable Markdown companion
and the plotted PNG figures used in that companion.

The command is designed for portability:

- it is part of the installed `trophosome` program;
- it does not depend on Pandoc, LaTeX, Microsoft Word or a particular Markdown
  reader;
- the PDF contains the figures, so relative image paths are not needed to read
  it;
- it discovers no-return, matched-return and mutation comparisons from the
  parameter matrix rather than from particular cell numbers; and
- it reads result summaries but never modifies raw simulations.

## 1. Install the reporting dependencies

From the repository:

```bash
python -m pip install -e '.[report]'
```

For development and tests:

```bash
python -m pip install -e '.[dev,report]'
```

## 2. Prepare the standardized analysis folder

The report command expects one analysis folder containing:

| File | Required content |
|---|---|
| `analysis-summary.json` | Experiment identity, initial diversity, audit result, seed blocks, source revision, limitations and resource totals; an optional benchmark-hardware note is used when available |
| `run-endpoints.tsv` | One row per independently simulated population, including cell ID, seed block, endpoint diversity, composition distance, mutation-fate and resource measures |
| `environment-trajectories.tsv` | One row per retained population and host passage; reserved for trajectory reporting |

The accompanying design TSV must contain one row per cell. The report currently
uses these columns:

```text
cell_id  label  H  f  R  alpha  u
```

Additional columns are allowed. New parameter columns can be added without
breaking the report; they remain in the design table for later report extensions.

The Phase 1 pilot analysis produces this data contract with:

```bash
python experiments/work/trophosome/p01-neutral-feedback/analysis/analyse_first_pilot.py
```

## 3. Generate the report

From the repository root:

```bash
trophosome report \
  --analysis experiments/work/trophosome/p01-neutral-feedback/analysis/s01-pilot-derived \
  --design experiments/work/trophosome/p01-neutral-feedback/design/phase1-first-pilot-cells.tsv \
  --output output/pdf/phase1-first-pilot-report.pdf \
  --markdown docs/phase1-first-pilot-report.md \
  --title 'Phase 1 first-pilot report' \
  --report-date 2026-08-12
```

`--report-date` is optional. Supplying it makes repeated builds display the same
date. If `--markdown` is omitted, the Markdown file is written next to the PDF.

## 4. Inspect the outputs

The command prints a JSON object listing three outputs:

1. the self-contained PDF;
2. the editable Markdown companion; and
3. the directory of PNG figure assets referenced by the Markdown.

The PDF is the portable reading copy. The Markdown and PNG files are the
transparent, version-controllable source for later edits.

The report contains:

- purpose and scope;
- the full pilot cell matrix;
- endpoint audit and computational feasibility;
- richness, Shannon entropy, Simpson diversity, Hill `D1`, Hill `D2`,
  Pielou evenness and total-variation composition distance;
- automatically detected matched-return and mutation comparisons;
- interpretation and limitations;
- reproducibility metadata; and
- a glossary of biological parameters and diversity measures.

## 5. Rebuild after adding cells

Add the new cells to the design matrix, run their independent populations, add
their rows to the standardized endpoint and trajectory tables, and rerun the
same command. Cell numbers are identifiers, not instructions to the reporter.

The reporter detects:

- no-return cells from `R = 0`;
- a matched-return series from equal positive `R`, equal `alpha`, mutation off
  and varying `H`; and
- a mutation series from equal `H`, `R` and `alpha` with varying `u`.

This means a later architecture or selection parameter can be added to the
matrix. Reports comparing those new factors will require a dedicated comparison
panel, but the current report will continue to build as long as the core columns
remain present.

## 6. Run on the HPC

Reporting is a single short post-processing job. It does not need the scheduler
or one process per simulation. Activate the same mamba environment used for the
model, install the `report` extra once, and run `trophosome report` after the
analysis tables have been collected. Generate the PDF in `data/` or another
backed-up location rather than in temporary scratch.
