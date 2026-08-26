# HPC environment probe

`probe_hpc_environment.sh` creates a read-only report describing the machine,
mamba environments, Python runtime, filesystem, background-process tools and
available memory-accounting interfaces. It does not contact the network,
install software or launch a simulation.

Run it from the shell environment normally used for jobs:

```bash
bash scripts/hpc/probe_hpc_environment.sh \
  --project-root "$HOME/data/CRF_project"
```

The default report name is
`hpc_environment_probe_YYYYMMDDTHHMMSSZ.txt`. The repository ignores these
reports because they contain machine-specific hostnames and paths. Inspect the
report, transfer it to the Mac, and provide it in the task used to prepare the
HPC mamba and machine profiles.

Use `--output PATH` to select another report location and `--quiet` to suppress
the console copy.

## Launch the model-2.1 first pilot

The fixed-regional-pool pilot uses:

- scientific model specification 2.1.0;
- trophosome software 0.7.0;
- output schema 2.3.0;
- all 20 first-pilot cells and three matched seed blocks (60 populations);
- migration fraction `m = 0.1` with the fixed `ip001-fisher100` regional pool;
- neutral within-host and free-living fitness; and
- eight simultaneous populations, with two host workers each, initially.

The generated variant is named `v210-m010`. It uses separate paths from the
completed no-migration pilot.

After transferring one clean, frozen repository revision to
`/data/qiulab/CRF_project/work/Modeling_trophosome`, enter that directory and
confirm that the `trophosome` environment contains the expected software:

```bash
eval "$(mamba shell hook -s bash)"
mamba activate trophosome
trophosome --version
git status --short
```

The version should be `0.7.0`, and the Git status should be empty. Then perform
the one-time machine-local storage setup below.

### Create the machine-local storage layout

The repository stores portable experiment definitions, but raw simulation
outputs and checkpoints belong in the server's external scratch directory. A
machine-local `layout.local.json` records this mapping and is intentionally
ignored by Git. Create it once after cloning the repository or moving it to a
new machine:

```bash
TROPHOSOME_HPC_ROOT=/home/qiulab/data/CRF_project
TROPHOSOME_REPOSITORY="$(pwd -P)"

python scripts/manage_project_layout.py init \
  --root "$TROPHOSOME_HPC_ROOT" \
  --project trophosome \
  --work-root "$TROPHOSOME_REPOSITORY/experiments/work" \
  --scratch-root "$TROPHOSOME_HPC_ROOT/scratch" \
  --data-root "$TROPHOSOME_REPOSITORY/experiments/data" \
  --phase p01
```

The command is additive and safe to repeat. It should report these locations:

```text
Work:    .../Modeling_trophosome/experiments/work/trophosome
Scratch: /home/qiulab/data/CRF_project/scratch/trophosome
Data:    .../Modeling_trophosome/experiments/data/trophosome
```

Review the generated mapping with:

```bash
python -m json.tool experiments/work/trophosome/layout.local.json
```

Do not commit `layout.local.json`: its absolute paths are specific to this
server. Initialization creates the required empty directories but does not run
the model or modify the frozen pilot configurations and manifests.

Then perform the non-simulating preflight:

```bash
bash scripts/hpc/launch_phase1_first_pilot_v2_1.sh --prepare-only
bash scripts/hpc/launch_phase1_first_pilot_v2_1.sh --dry-run
```

The preflight reads the machine-local `layout.local.json`, checks all 60 frozen
configuration checksums, confirms `m = 0.1`, confirms that the focal and
regional starting compositions are identical, and creates the isolated scratch
folders. It consumes no simulation time.

Because the probed machine did not expose a batch scheduler, run the launcher
inside `tmux` so it survives a disconnected terminal:

```bash
tmux new -s trophosome-pilot-v210
bash scripts/hpc/launch_phase1_first_pilot_v2_1.sh
```

Detach with `Ctrl-b`, then `d`. Reconnect later with:

```bash
tmux attach -t trophosome-pilot-v210
```

To use another concurrency after checking the first run's CPU, memory and
filesystem behaviour:

```bash
TROPHOSOME_PILOT_JOBS=4 bash scripts/hpc/launch_phase1_first_pilot_v2_1.sh
```

The full run can be stopped with `Ctrl-c`. Run the same command again to skip
completed populations and resume valid checkpoints. Per-population `run.out`,
`run.err`, `execution-summary.json` and model outputs are stored below the
scratch path in `p01-neutral-feedback/s01-pilot-v210-m010/`.

### Automatic analysis and report

The `trophosome` environment must contain the reporting dependencies before the
simulations finish. Verify them without writing files:

```bash
python -c 'import matplotlib, reportlab; print("reporting dependencies available")'
```

Once all 60 completion markers are present, the launcher automatically:

1. verifies the software, model and output-schema versions;
2. verifies final-state checksums and committed output sizes;
3. verifies five exchanges of exactly 100 million emigrants and immigrants in
   every population;
4. creates the migration-aware endpoint, trajectory and cell-summary tables;
5. creates biological and computational figures; and
6. knits the self-contained PDF and editable Markdown report.

The resulting files are:

```text
output/pdf/phase1-first-pilot-v210-m010-report.pdf
docs/phase1-first-pilot-v210-m010-report.md
docs/figures/phase1-pilot-v210-m010-report/
```

If even one population is missing or invalid, the report is not created and the
audit lists the affected run IDs. A report failure never deletes or changes raw
simulation outputs. Retry the report without rerunning populations using:

```bash
bash scripts/hpc/launch_phase1_first_pilot_v2_1.sh --report-only
```

The report completion record fingerprints all 60 model completion files, the
design, manifest and reporting code. An unchanged report is therefore skipped
safely if the launcher is run again.

The server filesystem was 95% occupied when it was probed. Check the pilot's
total scratch use periodically:

```bash
du -sh /home/qiulab/data/CRF_project/scratch/trophosome/p01-neutral-feedback/s01-pilot-v210-m010
```

## Launch the Phase 1 second pilot

The second pilot is the long-run stationarity-and-precision screen selected from
the completed model-2.1 first pilot. It is exploratory rather than
confirmatory. The frozen matrix contains six sentinel conditions, 12 matched
seed blocks per condition, and 250 host-population passages: 72 independently
runnable populations in total.

| Cell | Biological comparison | Hosts | Escape fraction | Total return | Mutation rate |
|---|---|---:|---:|---:|---:|
| `c0021` | Migration-only, no host return | 100 | 0 | 0 | 0 |
| `c0022` | Mutation-free baseline host feedback | 100 | `1e-2` | `1e9` | 0 |
| `c0023` | Mutation enabled at baseline feedback | 100 | `1e-2` | `1e9` | `1e-10` |
| `c0024` | Many hosts at the same total return as `c0022` | 10,000 | `1e-4` | `1e9` | 0 |
| `c0025` | Weak feedback | 100 | `1e-3` | `1e8` | 0 |
| `c0026` | Strong feedback | 100 | `1e-1` | `1e10` | 0 |

All conditions use fixed-regional-pool migration with `m = 0.1`, neutral
within-host and free-living fitness, complete environmental trajectories, and
one-hour checkpoint targets. Adult strain counts are complete for the
100-host cells and use the deterministic 100-host panel for `H = 10,000`.

The first-pilot measurements project about 35.7 summed runtime hours and 7.3
GiB for the full matrix. These are linear planning estimates. The longest cell
is projected at about 2.7 hours per population; actual HPC performance should
still be checked after the first completed seed block.

### Before launching

Pull the clean repository revision, activate the same `trophosome` mamba
environment, and enter the repository:

```bash
git pull --ff-only
eval "$(mamba shell hook -s bash)"
mamba activate trophosome
cd /home/qiulab/data/CRF_project/work/Modeling_trophosome
git status --short
trophosome --version
```

The source tree must be clean and the software version must be `0.7.0`. Reuse
the machine-local `experiments/work/trophosome/layout.local.json` created for
the first pilot. If it is absent, follow **Create the machine-local storage
layout** above before continuing.

Verify the frozen configurations and review all output destinations without
running the model:

```bash
python scripts/prepare_phase1_second_pilot.py --verify
bash scripts/hpc/launch_phase1_second_pilot.sh --prepare-only
bash scripts/hpc/launch_phase1_second_pilot.sh --dry-run
```

The preflight must report 72 populations, six sentinel cells, 250 passages, 12
seed blocks, and `m = 0.1`. It also checks every configuration checksum,
selection state, starting focal/regional composition, and output-retention
mode.

Before the full launch, a non-reporting single-population check may be run and
resumed safely:

```bash
bash scripts/hpc/launch_phase1_second_pilot.sh \
  --cell c0021 --seed-block sb0001 --no-report
```

This population is part of the frozen 72-run design; it is not a disposable
extra replicate. The later full launcher will detect and skip it.

### Full launch, interruption, and resume

Run inside `tmux` because the probed server does not expose a batch scheduler:

```bash
tmux new -s trophosome-second-pilot
bash scripts/hpc/launch_phase1_second_pilot.sh
```

Eight populations run concurrently by default, with two host workers per
population. Change this only after inspecting machine load:

```bash
TROPHOSOME_SECOND_PILOT_JOBS=4 \
  bash scripts/hpc/launch_phase1_second_pilot.sh
```

Detach from `tmux` with `Ctrl-b`, then `d`, and reconnect with:

```bash
tmux attach -t trophosome-second-pilot
```

`Ctrl-c` requests an orderly stop. Run the same launcher again to skip complete
populations and resume interrupted ones from their newest valid checkpoints.
Raw results are stored below:

```text
/home/qiulab/data/CRF_project/scratch/trophosome/p01-neutral-feedback/
  s02-equilibrium-precision-v210-m010-g250/
```

### Automatic and stand-alone report

Verify the report libraries before starting the full run:

```bash
python -c 'import matplotlib, numpy, reportlab; print("reporting dependencies available")'
```

After all 72 populations complete, the launcher automatically audits the raw
outputs, creates the analysis tables and four biological figures, and writes:

```text
output/pdf/phase1-second-pilot-v210-m010-g250-report.pdf
docs/phase1-second-pilot-v210-m010-g250-report.md
docs/figures/phase1-second-pilot-v210-m010-g250-report/
```

The PDF is self-contained. The Markdown copy and PNG figures are editable. The
report explains environmental D0, D1, D2, evenness and compositional change;
late-run stationarity diagnostics; the earliest assessment from which the
screen remains satisfied; continuing fluctuations; and recommended replicate
numbers for the confirmatory experiment.

No report is created unless all 72 committed outputs pass their version,
configuration, size, final-checksum, reservoir-capacity, trajectory-completeness
and migration audits. A report error never alters the raw simulations. Build or
rebuild the report later, without launching simulations, using either:

```bash
bash scripts/hpc/launch_phase1_second_pilot.sh --report-only
```

or the direct stand-alone command:

```bash
python scripts/build_phase1_second_pilot_report.py \
  --repository "$(pwd -P)" --force
```

Without `--force`, an unchanged report is skipped using a fingerprint of the
72 completion records, design, manifest, analysis code, and report code. To run
simulations without trying the automatic report, add `--no-report`.
