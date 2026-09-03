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

## How the launchers select Python

Every maintained Phase 1 launcher can be called either before or after the
`trophosome` environment is activated. If that environment is already active,
the launcher uses its Python directly and does not require the `mamba` command
to remain in `PATH`. If it is not active, the launcher initializes mamba,
activates `trophosome`, and then uses that environment's Python. It therefore
is not necessary to deactivate the environment for `--prepare-only` or
`--dry-run`, nor to reactivate it before the full launch.

Before a simulation starts, the launcher checks for uncommitted changes under
`src/trophosome` and in `pyproject.toml`. Frozen experiment files are checked
separately by their recorded checksums. Unrelated historical results, reports,
notebooks and legacy directories may appear in `git status`, but they no longer
trigger the model-source safety error. A genuine source-code change still does.

## Phase 1 Stage 3 Wave 2

Wave 2 tests host number by infection bottleneck and host feedback by regional
exchange. Its complete passage-100 analysis has 40 conditions and 12 matched
seed blocks: 408 new populations plus 72 exact reused populations. The adaptive
continuation is bounded at passages 500 and 1,000 and is described biologically
in [`docs/phase1-stage3-wave2.md`](../../docs/phase1-stage3-wave2.md).

After pulling the frozen revision, activate the `trophosome` environment and
confirm that the maintained model source is clean. Reuse the machine-local
`layout.local.json` described below.

```bash
eval "$(mamba shell hook -s bash)"
mamba activate trophosome
cd /home/qiulab/data/CRF_project/work/Modeling_trophosome
git status --short
trophosome --version
python scripts/prepare_phase1_stage3_wave2.py --verify
```

Review and prepare the initial batch without simulating:

```bash
bash scripts/hpc/launch_phase1_stage3_wave2.sh --prepare-only
bash scripts/hpc/launch_phase1_stage3_wave2.sh --dry-run
```

The dry run must report 408 new populations toward passage 100. Run the three
included safety populations and assess their observed time and storage:

```bash
bash scripts/hpc/launch_phase1_stage3_wave2.sh --smoke-only
bash scripts/hpc/launch_phase1_stage3_wave2.sh --check-smoke
```

If the safety gate passes, start the full initial batch inside `tmux`:

```bash
tmux new -s trophosome-stage3-wave2-g100
bash scripts/hpc/launch_phase1_stage3_wave2.sh
```

The initial launcher stops every new trajectory cleanly at passage 100. When
all 408 states pass their checksum audit, it freezes the outcome-dependent
passage-100 decision. It does not start the next boundary. Review it with:

```bash
python -m json.tool \
  experiments/work/trophosome/p01-neutral-feedback/analysis/\
s03-parameter-map-wave2-v210-adaptive-g1000-derived/\
adaptive-horizon-decision-g100.json
```

Then review the authorized continuation without running it:

```bash
bash scripts/hpc/launch_phase1_stage3_wave2.sh --horizon 500 --dry-run
```

Start only those selected populations, preferably in a new `tmux` session:

```bash
tmux new -s trophosome-stage3-wave2-g500
bash scripts/hpc/launch_phase1_stage3_wave2.sh --horizon 500
```

After the passage-500 states are audited, the launcher freezes the second
decision. Inspect it and then explicitly launch the final authorized subset:

```bash
python -m json.tool \
  experiments/work/trophosome/p01-neutral-feedback/analysis/\
s03-parameter-map-wave2-v210-adaptive-g1000-derived/\
adaptive-horizon-decision-g500.json

bash scripts/hpc/launch_phase1_stage3_wave2.sh --horizon 1000 --dry-run
tmux new -s trophosome-stage3-wave2-g1000
bash scripts/hpc/launch_phase1_stage3_wave2.sh --horizon 1000
```

To rebuild a missing decision without simulating, use `--assess-only` with
`--horizon 100` or `500`. If the decision already exists, the assessor refuses
to replace a different result. Repeating a launcher safely skips a trajectory
that already reached the requested boundary and resumes a valid checkpoint for
an interrupted one.

`Ctrl-c` requests a clean stop; repeat the same command to resume. Do not edit
the model package (`src/trophosome/` or `pyproject.toml`) between boundaries,
because every continuation verifies the source checksum stored in its
checkpoint. Set `TROPHOSOME_STAGE3_WAVE2_JOBS` to change the default eight
simultaneous populations after reviewing the safety measurements.

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
confirmatory. Its parameter matrix remains frozen at six sentinel conditions
and 250 host-population passages. The original batch used 12 matched seed blocks
per condition (72 populations). The Stage 2 closure batch adds eight matched
seed blocks (`sb0013`--`sb0020`) to every condition (48 new populations), giving
20 seed blocks and 120 populations in the combined Stage 2 analysis.

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

The first-pilot measurements project about 23.8 additional summed runtime hours
and 4.9 GiB for the 48-population closure batch, or about 59.5 hours and 12.2
GiB for all 120 Stage 2 populations. These are linear planning estimates. The
longest cell is projected at about 2.7 hours per population; actual HPC
performance should still be checked after the first completed closure seed
block. The 72 completed populations must remain in scratch because the combined
report audits and analyses all 120 populations together.

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
bash scripts/hpc/launch_phase1_second_pilot_closure.sh --prepare-only
bash scripts/hpc/launch_phase1_second_pilot_closure.sh --dry-run
```

The closure preflight must report 48 populations, six sentinel cells, 250
passages, eight selected seed blocks out of 20 frozen seed blocks, and migration
fraction `m = 0.1`. It also checks every configuration checksum, selection
state, starting focal/regional composition, and output-retention mode.

Before the full launch, a non-reporting single-population check may be run and
resumed safely:

```bash
bash scripts/hpc/launch_phase1_second_pilot_closure.sh \
  --cell c0021 --seed-block sb0013 --no-report
```

This population is part of the frozen closure batch; it is not a disposable
extra replicate. The later closure launcher will detect and skip it.

### Full launch, interruption, and resume

Run inside `tmux` because the probed server does not expose a batch scheduler:

```bash
tmux new -s trophosome-second-pilot-closure
bash scripts/hpc/launch_phase1_second_pilot_closure.sh
```

Eight populations run concurrently by default, with two host workers per
population. Change this only after inspecting machine load:

```bash
TROPHOSOME_SECOND_PILOT_JOBS=4 \
  bash scripts/hpc/launch_phase1_second_pilot_closure.sh
```

Detach from `tmux` with `Ctrl-b`, then `d`, and reconnect with:

```bash
tmux attach -t trophosome-second-pilot-closure
```

`Ctrl-c` requests an orderly stop. Run the same closure launcher again to skip
complete populations and resume interrupted ones from their newest valid
checkpoints. The closure wrapper fixes the selection to `sb0013`--`sb0020`; any
additional options such as `--jobs`, `--cell`, or `--no-report` are forwarded to
the shared Stage 2 runner. Seed blocks for later confirmatory experiments must
use a separately frozen series rather than reusing these exploratory seeds.
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

After all 120 populations complete, the launcher automatically audits the raw
outputs, creates the analysis tables and six biological figures, and writes:

```text
output/pdf/phase1-second-pilot-v210-m010-g250-report.pdf
docs/phase1-second-pilot-v210-m010-g250-report.md
docs/figures/phase1-second-pilot-v210-m010-g250-report/
```

The PDF is self-contained. The Markdown copy and PNG figures are editable. The
report explains environmental D0, D1, D2, evenness and compositional change;
late-run stationarity diagnostics; the post-hoc diagnostic separating temporal
stability from persistent seed-block heterogeneity; the earliest assessment from
which the registered screen remains satisfied; continuing fluctuations; and
recommended replicate numbers for the confirmatory experiment.

No report is created unless all 120 committed outputs pass their version,
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
120 completion records, design, manifest, analysis code, and report code. To run
simulations without trying the automatic report, add `--no-report`.

## Phase 1 Stage 3: first mapping wave

This is the first **main parameter-mapping batch**, not a rerun of either pilot.
The [frozen part-one design](../../docs/phase1-stage3-wave1.md) contains 24 new cells
(`c0027`-`c0050`) x twelve matched seed blocks (`sb0001`-`sb0012`) = 288 new populations,
all at 100 passages, `m=0.1`, and at most 10,000 hosts. The shared no-return
control reuses twelve Stage 2 c0021 passage-100 outcomes, giving 25 primary
conditions and 300 primary populations. Five other Stage 2 conditions are
supplementary only. No reused pilot is rerun. Selection remains off
in both habitats. Model/software/output versions remain 2.1.0/0.7.0/2.3.0.

### 1. Update and verify the environment

The current part-one design is frozen in commit `d0fc904`; use that revision
or a later revision containing it. Connect to the server, enter the existing
repository and update it:

Finish any active Stage 2 jobs before updating this checkout. Checkpoints are
tied to their original source fingerprint; if an older run still needs to be
resumed, keep its original revision in a separate checkout until it finishes.

```bash
cd /home/qiulab/data/CRF_project/work/Modeling_trophosome
git pull --ff-only
eval "$(mamba shell hook -s bash)"
mamba activate trophosome
python -m pip install -e '.[report]'
trophosome --version
git status --short
```

The version should be `0.7.0` and Git status should be empty before launching.
If Git reports local changes or a pull conflict, stop and inspect them rather
than overwriting them. Reuse the existing `layout.local.json`; if missing,
follow [the one-time storage setup](#create-the-machine-local-storage-layout).
There is no need to regenerate the frozen TOMLs on the HPC.

### 2. Prepare the 288 isolated run directories (no simulation)

```bash
python scripts/prepare_phase1_stage3_wave1.py --verify
bash scripts/hpc/launch_phase1_stage3_wave1.sh --prepare-only
bash scripts/hpc/launch_phase1_stage3_wave1.sh --dry-run
```

Expected: 317 frozen files verified and a preflight listing 288 new populations.
All new raw outputs go under the configured scratch root, in
`p01-neutral-feedback/s03-parameter-map-wave1-v210-m010-g100/`.
Existing Stage 1 and Stage 2 outputs are untouched.
The old unlaunched g250 matrix has been replaced. Do not reuse its TOMLs or
scratch paths; verify that the launcher prints **100 passages** before launch.

### 3. Run three safety populations first

Use the site's permitted compute session/allocation. On the previously probed
server, which did not expose a scheduler, use `tmux`:

```bash
tmux new -s trophosome-stage3-wave1
bash scripts/hpc/launch_phase1_stage3_wave1.sh --smoke-only --jobs 3
```

This runs c0034 (H=100, mutation on, 99% release), c0049 (H=10,000, mutation
off), and c0050 (H=10,000, mutation on), all at alpha=0.99, with `sb0001`.
They are three of the 288 planned new populations. Detach with
`Ctrl-b`, then `d`; reconnect with:

```bash
tmux attach -t trophosome-stage3-wave1
```

### 4. Review the measured resource check

After the three safety populations finish:

```bash
bash scripts/hpc/launch_phase1_stage3_wave1.sh --check-smoke
```

Look for `"passed": true`. The check audits all three completions and estimates
full-wave runtime/storage from their measurements with a 2x safety margin.
Mutation-enabled costs interpolate between H=100 and H=10,000; mutation-free
costs scale from H=10,000. No-return and supplementary references add no new runs.
It requires at most 48 projected hours per population, less than 350 GiB
projected output and sufficient free scratch space. Confirm your user quota
and compute allocation separately. These are projections, not hard runtime
limits or a guarantee of a 48-hour entire batch. If the screen fails, review
its complete output before proceeding; do not bypass the check.
A resumed safety population also needs review because the existing timing
record covers only the most recent attempt.

### 5. Finish the first wave

In the same persistent compute session:

```bash
bash scripts/hpc/launch_phase1_stage3_wave1.sh
```

The launcher verifies the safety gate again, audits and skips the three
completed safety runs, then runs the remaining 285 populations. By default,
eight populations run simultaneously, with two host workers each. To lower
the load, use `--jobs 4` (or set `TROPHOSOME_STAGE3_JOBS=4`). Never exceed your
CPU allocation; population-manager processes also use some CPU.

`Ctrl-c` requests a graceful stop. Repeating the same command resumes from
the latest valid checkpoints and skips audited completed populations. Two
recovery checkpoints are retained while a run is active, targeting one per
hour at completed host-passage boundaries. Per-run `run.out`, `run.err`, and
`execution-summary.json` files are in scratch, not the repository.

### 6. Automatic and stand-alone reports

When all 288 new populations and the frozen references pass the audit, the
launcher automatically creates:

```text
output/pdf/phase1-stage3-wave1-v210-m010-g100-report.pdf
docs/phase1-stage3-wave1-v210-m010-g100-report.md
docs/figures/phase1-stage3-wave1-v210-m010-g100-report/
```

The self-contained PDF and editable companion include diversity at passage
100, Shannon/Simpson and Hill indices, individual TV trajectories, paired
H-by-feedback interactions and precision/late-drift checks. The report avoids
assuming equilibrium. All comparisons reuse the matching Stage 2 passage-100
endpoints, never passage 250. Off-grid references are displayed at their actual
alpha values, not relabelled as 0.1 or 0.99.

Rebuild it at any time after completion, without launching simulations:

```bash
bash scripts/hpc/launch_phase1_stage3_wave1.sh --report-only
```

Or, from the activated environment:

```bash
python scripts/build_phase1_stage3_wave1_report.py --repository "$(pwd -P)"
```

Both paths re-audit the data on every build. Do not pipe them through `head`.
`--no-report` defers automatic reporting if necessary; a report failure never
requires rerunning valid simulations. Derived analysis tables and the report
completion record live in
`experiments/work/trophosome/p01-neutral-feedback/analysis/s03-parameter-map-wave1-v210-m010-g100-derived/`.
No raw files are automatically deleted. Review the report before archiving
results or selecting the next adaptive batch.
