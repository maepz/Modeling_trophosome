# Reproducibility workflow

## One configuration per experiment

Copy a versioned TOML file from `configs/`, change only declared inputs, and give
the output directory a unique experiment name. Do not edit model parameters in a
notebook or source file.

```bash
trophosome validate configs/my_experiment.toml
trophosome run configs/my_experiment.toml --output results/my_experiment
```

This interface runs the main exact-count prototype. Legacy V3.1 does not use it
and is not fully reproducible because some random stages and host workers create
unseeded generators. Use V3.1 only as a documented comparator until it has a
complete seed interface.

The run records:

- the fully resolved configuration;
- root seed and child seed derivation strategy;
- Git commit;
- Python, NumPy, and platform versions;
- software, model-specification and output-schema versions;
- a SHA-256 hash of the fully resolved configuration;
- host-generation summaries;
- compact environmental, infection, adult-summary, pooled-occupancy, release,
  fixed-pool migration and mutation-lineage tables;
- an initial-strain origin table identifying focal, regional and shared source
  strains;
- explicit within-host and free-living fitness metadata for initial and mutant
  strains, plus habitat-specific mean fitness summaries;
- optional full or reproducibly sampled-panel adult count tables;
- validated recovery checkpoints while a run is active and verified final
  environmental states after completion.

## Randomness

Use NumPy `Generator` objects created from `SeedSequence`. Never call the global
`np.random.seed()` inside a worker. The maintained runner derives a stream from
the seed, replicate, host generation, host ID and stochastic stage. It also
reserves strain-ID blocks by logical host, so task ordering and worker count do
not change a seeded run.

## Outputs

Raw outputs belong under `results/` or an external research-data store and are
ignored by Git. Promote only small, documented fixtures into `tests/data/`.
Record checksums for externally archived releases.

Do not pickle long-lived research results. Pickle is Python-version-sensitive and
unsafe to load from untrusted sources. The scalable path uses CSV/JSON metadata
and typed NumPy arrays.

## Checkpoints and restart

Recovery checkpoints are written atomically at completed host-generation
boundaries. The configured wall-clock interval is a target: `"1h"` means the
first safe boundary after one hour has elapsed. A checkpoint is also written at
every replicate boundary.

Each checkpoint contains the environmental strain state, active mutational
depths, logical replicate and generation coordinates, seed-scheme and
strain-identifier metadata, accumulated generation summaries, committed byte
offsets for every CSV table, the resolved configuration, and configuration and
source hashes. It is reopened and checksum-validated before older recovery
points are removed. The newest two valid checkpoints are retained while the run
is active.

Resume an interrupted run with the same configuration and output directory:

```bash
trophosome run configs/my_experiment.toml \
  --output results/my_experiment \
  --resume
```

Before continuing, resume verifies the source and state-affecting configuration,
then truncates every output table to its last committed offset. A corrupt newest
checkpoint is ignored in favour of the preceding valid checkpoint. Worker count,
host batch size, in-flight batches, checkpoint duration and checkpoint retention
may change because they do not alter the scientific transition law or logical
random streams. Scientific parameters and the planned number of generations may
not change.

After successful completion, the software verifies all output sizes and final
environmental states, writes `completion.json`, and removes all recovery
checkpoints. The final `final_environment_repNNN.npz` files remain. An idempotent
`--resume` on an already completed run verifies it and exits without rerunning
hosts.

## Validation levels

1. unit tests for invariants and edge cases;
2. stochastic regression tests comparing distributions, not identical
   trajectories;
3. exact-versus-approximate validation on small matched scenarios;
4. convergence tests over samples, replicates, and any rescaling factor;
5. biological sensitivity analysis over uncertain life-history parameters.
