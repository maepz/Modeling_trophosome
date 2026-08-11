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

The prepared 12-cell core is listed in
`work/trophosome/p01-neutral-feedback/design/phase1-first-pilot-cells.tsv`.
Each cell has three independently runnable seed blocks. The corresponding 36
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

Run the first seed block sequentially on a Mac with:

```bash
.venv/bin/python scripts/run_phase1_first_pilot.py --seed-block sb0001
```

The launcher uses a cost-aware cell order, resumes a valid checkpoint when an
interrupted output directory contains one, and records elapsed time, peak
process-tree memory and output size in each scratch run directory. Use
`--dry-run` to review the selection without launching it, or repeat `--cell`
to select individual cells. Distinct cell/seed-block pairs can be assigned to
separate invocations on an HPC because they never share an output directory.
