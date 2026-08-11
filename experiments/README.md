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
