# Experiments

This directory contains the portable, version-controlled part of the staged
trophosome experiments.

- `work/trophosome/` contains cell registries, designs, configurations,
  manifests, wrappers and analysis code.
- `data/trophosome/` contains only curated releases that have passed completion
  and checksum validation.
- Raw simulation outputs, active checkpoints and scheduler logs belong in the
  external scratch path recorded in the untracked `layout.local.json`.

Initialize or relocate the machine-specific directory structure with
`scripts/manage_project_layout.py`. Do not edit the normalized registries from
parallel compute jobs.

## Phase 1 first pilot

The staged 20-cell design (12 core cells followed by eight gated extension
cells) is listed in
`work/trophosome/p01-neutral-feedback/design/phase1-first-pilot-cells.tsv`.
Each cell has three independently runnable seed blocks. The corresponding 60
run configurations and output paths are indexed by
`work/trophosome/p01-neutral-feedback/manifests/phase1-first-pilot-runs.tsv`.

One seed block represents one stochastic population and has one master seed.
The model derives all infection-, host-, mutation-, release-, regulation- and
tie-breaking streams from that seed and their logical coordinates. Reusing the
same seed-block ID across cells provides the matched stochastic comparison.

After cloning or relocating the repository, recreate the empty machine-local
run directories from `layout.local.json` with:

```bash
python scripts/prepare_phase1_pilot_scratch.py --write
python scripts/prepare_phase1_pilot_scratch.py --verify
```

These commands prepare directories and `run.json` manifests only. They do not
launch simulations.

Run the complete staged workflow sequentially on a Mac with:

```bash
.venv/bin/python scripts/run_phase1_first_pilot.py
```

By default, the launcher runs or verifies all three seed blocks of the core 12,
creates the interim PDF, assesses the expansion gates, and only then runs the
eight extension cells. It overwrites the same PDF with the final 20-cell report.
The launcher uses a cost-aware cell order, resumes a valid checkpoint when an
interrupted output directory contains one, and records elapsed time, peak
process-tree memory and output size in each scratch run directory.

Use `--dry-run` to review the staged workflow. Supplying `--seed-block` or
`--cell` switches to selected-run mode for manual Mac or HPC dispatch. Distinct
cell/seed-block pairs can be assigned to separate invocations on an HPC because
they never share an output directory.

## Model 2.1 fixed-regional-pool rerun

The migration-enabled rerun retains the same 20 cells and the same three
matched seed blocks, but uses scientific model specification 2.1.0. After each
host return, exactly 10% of the one-billion-cell focal environment is exchanged
with a fixed, non-depleting regional source. The regional source has the same
starting composition as `ip001-fisher100`. Both within-host and free-living
selection remain disabled.

This is stored as the separate `v210-m010` variant. Its run IDs, configurations,
manifest and scratch paths cannot overlap the completed historical pilot. Build
or verify its 20 cell configurations and 60 run configurations with:

```bash
python scripts/prepare_phase1_first_pilot_v2_1.py --write
python scripts/prepare_phase1_first_pilot_v2_1.py --verify
```

On the HPC, the wrapper uses the mamba environment named `trophosome` and runs
eight independent populations concurrently by default. Each population retains
the two host workers specified in its TOML:

```bash
bash scripts/hpc/launch_phase1_first_pilot_v2_1.sh --prepare-only
bash scripts/hpc/launch_phase1_first_pilot_v2_1.sh --dry-run
bash scripts/hpc/launch_phase1_first_pilot_v2_1.sh
```

The launcher checks the scientific-model, software and output-schema versions,
validates every selected TOML and checksum, creates machine-local run manifests,
and refuses to start simulations from an uncommitted source tree. Repeating the
last command skips complete populations and resumes interrupted populations
that have a valid checkpoint. Runtime, peak memory and console logs are written
inside each population's scratch directory.

Set `TROPHOSOME_PILOT_JOBS` to change the number of simultaneous populations.
For example, `TROPHOSOME_PILOT_JOBS=4` uses at most eight host workers. The
historical no-migration analysis is never used for this variant.

When all 60 populations have valid model-2.1 completion markers, the launcher
automatically runs the dedicated migration-aware audit, writes standardized
analysis tables, creates report figures and knits both:

- `output/pdf/phase1-first-pilot-v210-m010-report.pdf`; and
- `docs/phase1-first-pilot-v210-m010-report.md`.

The report treats the no-return cells as the migration-only baseline, audits
the exact emigrant and immigrant totals at every passage, and shows how regional
exchange changes the expected host-derived fraction and richness. If the report
stage fails, raw completed simulations are not modified. After correcting the
reporting environment, retry only that stage with:

```bash
bash scripts/hpc/launch_phase1_first_pilot_v2_1.sh --report-only
```

An unchanged, successfully fingerprinted report is skipped on later launcher
runs. Use `--no-report` only when simulations should finish without automatic
analysis and knitting.
