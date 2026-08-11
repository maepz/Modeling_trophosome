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
