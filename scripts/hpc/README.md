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
mamba run -n trophosome trophosome --version
git status --short
```

The version should be `0.7.0`, and the Git status should be empty. Then perform
the non-simulating preflight:

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
mamba run -n trophosome python -c 'import matplotlib, reportlab; print("reporting dependencies available")'
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
