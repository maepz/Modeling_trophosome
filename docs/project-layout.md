# Project layout and cell registry

The standard-library helper `scripts/manage_project_layout.py` creates the
manifest-driven `work`, `scratch`, and `data` hierarchy and provides quick cell
lookup by stable cell ID.

## Initialize a project

On a system where the three storage areas are below the same home directory:

```bash
python scripts/manage_project_layout.py init \
  --root "$HOME" \
  --project trophosome
```

If the storage parents differ, override them explicitly:

```bash
python scripts/manage_project_layout.py init \
  --root "$HOME" \
  --project trophosome \
  --work-root /work/my-user \
  --scratch-root /scratch/my-user \
  --data-root /data/my-user
```

Initialization is additive and idempotent. Existing registry files are not
overwritten.

The resulting top-level structure is:

```text
ROOT/
├── work/trophosome/       # wrappers, configs, logs, manifests, analyses
├── scratch/trophosome/    # replaceable run outputs and checkpoints
└── data/trophosome/       # curated releases retained or downloaded
```

Within `work/` and `scratch/`, directories are first divided by phase and then
by experimental stage. A normal initialization creates only Phase 1 and its
five stages. Phase 2 or Phase 3 directories are created automatically when the
first cell for that phase is registered. They can also be prepared explicitly:

```bash
python scripts/manage_project_layout.py init --root "$HOME" --phase p02
```

### Repository-backed layout

For a portable Mac/HPC workflow, keep experiment `work` and curated `data`
inside the model repository and place only scratch outside it:

```bash
python scripts/manage_project_layout.py init \
  --root /path/to/local-project-home \
  --project trophosome \
  --work-root /path/to/Modeling_trophosome/experiments/work \
  --data-root /path/to/Modeling_trophosome/experiments/data \
  --scratch-root /path/to/external/scratch
```

The generated `layout.local.json` contains machine-specific absolute paths and
is ignored by Git. Re-run `init` after cloning on another computer to regenerate
that local file and any empty directories. Registries, configurations,
manifests and curated releases remain portable and version controlled.

## Register a cell

```bash
python scripts/manage_project_layout.py register-cell \
  --root "$HOME" \
  --cell-id p01-s03-c0042 \
  --label "Baseline mutation mapping" \
  --group M \
  --architecture arch-panmictic-v1 \
  --initial-population ip001-fisher100 \
  --param H=10000 \
  --param B=10 \
  --param K=1e9 \
  --param growth_factor=1.2 \
  --param steady_generations=500 \
  --param f=1e-4 \
  --param e=100000 \
  --param R=1000000000 \
  --param u=1e-10 \
  --param c=1 \
  --unit H=hosts \
  --unit K=cells \
  --unit e=cells_per_host \
  --unit R=cells \
  --role e=derived
```

The script generates a mnemonic such as:

```text
h10000-f1em4-u1em10-c1-b10-k1000000000-ts500-archpanmictic
```

The compact labels mean:

| Label | Biological meaning |
|---|---|
| `h` | host abundance, `H` |
| `f` | escape fraction per host, `f` |
| `u` | whole-genome mutation probability per bacterial generation, `u` |
| `c` | environmental capacity ratio, `c` |
| `b` | infection bottleneck, `B` founder cells per host |
| `k` | within-host carrying capacity, `K`, in cells |
| `ts` | steady within-host bacterial generations per host lifespan |
| `arch` | architecture profile, for example `panmictic` or `lobed` |
| `sel` | selection profile, used when selection is activated in Phase 2 |
| `fit` | initial strain-fitness profile |

Small numbers are made filename-safe: `1em4` means \(1 \times 10^{-4}\), and
`1em10` means \(1 \times 10^{-10}\).

The numeric cell ID remains the permanent identity. Its result directory is
named only with the full cell ID, for example `p01-s03-c0042/`. The mnemonic is
kept in `work/trophosome/registry/cells.csv` for human review and lookup; it is
not included in directory names.

The open-ended parameters are stored in long form in
`work/trophosome/registry/cell_parameters.csv`. Future parameters such as
`architecture_mode`, `lobe_count`, or `migration_probability` require no
registry schema change.

For Phase 2, retain the same biological parameters and set `--selection` and
`--fitness` to versioned profile IDs. For Phase 3, set `--architecture` to the
appropriate versioned architecture profile and add any architecture-specific
values with repeated `--param NAME=VALUE` arguments.

Register cells from one planning process before submitting a local or HPC job
array. Compute jobs should read the frozen registries and manifests, not append
to the shared CSV files concurrently.

## Review a cell

```bash
python scripts/manage_project_layout.py show-cell \
  --root "$HOME" \
  p01-s03-c0042
```

If the cell number is unique in the registry, the shorter form is sufficient:

```bash
python scripts/manage_project_layout.py show-cell \
  --root "$HOME" \
  c0042
```

If `c0042` occurs in more than one phase or stage, the command reports the
matching full IDs and asks for an unambiguous one.

The lookup also accepts a cell directory path or a complete run ID:

```bash
python scripts/manage_project_layout.py show-cell \
  --root "$HOME" \
  p01-s03-c0042-sb0007
```

The parameters are not mathematically decoded from `c0042`. Instead, the cell
number is resolved to its stable full ID, which is then used as a registry key.
This is what permits later phases to add parameters without renaming older
cells or making paths unmanageably long.

## Translate between cell IDs and mnemonics

Return only the mnemonic associated with a cell ID or unique short number:

```bash
python scripts/manage_project_layout.py mnemonic \
  --root "$HOME" \
  c0042
```

Perform the reverse lookup using the exact mnemonic:

```bash
python scripts/manage_project_layout.py cell-id \
  --root "$HOME" \
  h10000-f1em4-u1em10-c1-b10-k1000000000-ts500-archpanmictic
```

The same parameter mnemonic may occur in more than one experimental stage. In
that case, `cell-id` prints every matching full ID. Use `--phase p01` and/or
`--stage s03` when a single filtered result is required.

For programmatic use:

```bash
python scripts/manage_project_layout.py show-cell \
  --root "$HOME" \
  --json \
  p01-s03-c0042
```

## List cells

```bash
python scripts/manage_project_layout.py list-cells --root "$HOME"
python scripts/manage_project_layout.py list-cells --root "$HOME" --phase p03
python scripts/manage_project_layout.py list-cells --root "$HOME" --contains lobed
```

## Registry rules

- One cell ID represents one scientific parameter combination.
- Never reuse a cell ID after changing a scientific parameter.
- One `(cell_id, replicate_id)` pair represents one stochastic run. Reusing the
  same replicate IDs across cells creates matched seed blocks for comparisons.
- One seed block has one master seed; the simulation derives the many random
  streams needed by its hosts, infections and biological processes from that
  seed and stable logical coordinates.
- The cell registry supports quick review but is not a substitute for the
  versioned TOML configuration.
- Each completed run must retain its own `resolved_config.json` and provenance.
- Each HPC job writes to a unique cell and seed-block directory.
- Recovery checkpoints are removed after a successful run, leaving the final
  environmental state and required scientific outputs.
