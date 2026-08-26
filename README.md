# Trophosome symbiont evolution model

`trophosome` simulates how microbial symbiont strains change as they repeatedly
move between an environmental reservoir and a population of hosts. It was
developed from the association between tubeworms and their bacterial symbionts,
but its biological scope is broader than tubeworms.

## What does one simulation represent?

During each modelled host-population generation:

1. every host acquires a small sample of symbionts from a shared environmental
   reservoir;
2. those symbionts reproduce inside the host and may mutate or experience
   selection;
3. some symbionts leave the hosts and mix with the environmental population;
4. the focal environment can exchange bacteria with a fixed regional source;
5. the updated focal environment supplies the next generation of hosts.

The model records both the abundance of labelled strains and the mutations that
created new strains. This makes it possible to ask whether repeated host passage
removes diversity, creates diversity, redistributes strains, or changes the
long-term composition of the environmental population.

## Which symbioses can be studied?

The model is not restricted to tubeworms or to symbionts housed in a trophosome.
It can represent other host--microbial associations when symbionts are acquired
horizontally from an environmental source and their life cycle can reasonably be
approximated by the four stages above.

The infection bottleneck, within-host population size, time spent in the host,
host abundance, symbiont release and environmental assumptions must be
calibrated for the organism being studied. The current implementation represents
haploid microbial symbionts. The name `trophosome` is retained to acknowledge
the biological system that motivated the software.

## Start here: run the small example

Python 3.11 or newer is required. The commands below are entered in a terminal
opened in the downloaded `Modeling_trophosome` folder.

### 1. Install the software

Creating a virtual environment gives this project its own isolated Python
installation. On macOS or Linux:

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -e .
```

On Windows, activate the environment with `.venv\Scripts\activate` instead of
the `source` command.

### 2. Check the example configuration

```bash
trophosome validate configs/smoke.toml
```

Validation checks the parameter values and reports potential feasibility
problems. It does not run the simulation or create results.

### 3. Run the example

```bash
trophosome run configs/smoke.toml --output results/smoke
```

`configs/smoke.toml` is deliberately tiny and should finish quickly. Its purpose
is to confirm that the installation works; its parameter values are not proposed
as a biological scenario.

### 4. Inspect the result

Open `results/smoke/host_generation_summary.csv` in R, Python or a spreadsheet
program. It contains one row for each replicate and host-population generation.
The complete strain counts and run information are in the same results folder.

If a longer run is interrupted after a recovery checkpoint has been written,
continue it with the same configuration and output folder:

```bash
trophosome run configs/my_experiment.toml --output results/my_experiment --resume
```

Resume verifies the scientific configuration and exact source code, removes any
incomplete rows written after the last safe host-generation boundary, and then
continues from the newest valid checkpoint. CPU worker and batching settings may
be changed when moving the same run between a Mac and an HPC node.

For reminders at any time, use:

```bash
trophosome -h
trophosome validate -h
trophosome run -h
```

## Set up an experiment

Each experiment is described by a plain-text TOML file in `configs/`. Duplicate
one of the examples, give it an informative name, edit its values, validate it,
and save the run to a new results folder. Keeping one configuration and one
output folder per experiment prevents accidental mixing of parameter sets.

The supplied configurations have different purposes:

| Configuration | Intended use |
| --- | --- |
| `configs/smoke.toml` | Very small installation check |
| `configs/toy.toml` | Small exploratory or demonstration run |
| `configs/biological-scale.example.toml` | Starting point for staged scaling tests; not biologically calibrated |
| `configs/phase2-selection.example.toml` | Example with both habitat-specific fitness traits and selection |
| `configs/migration.example.toml` | Small example of exchange with a fixed regional source; selection remains off |

### Parameters most biologists will change

| Parameter | Biological meaning |
| --- | --- |
| `replicates` | Number of independent stochastic repetitions of the same experiment |
| `initial_counts` | Starting relative abundances of environmental strains |
| `capacity_ratio` | Effective environmental reservoir size relative to one within-host carrying capacity |
| `migration.mode` | Whether regional exchange is disabled or uses a fixed regional source |
| `migration.fraction` | Fraction of the focal reservoir replaced by regional immigrants in each host generation |
| `migration.regional_counts` | Fixed regional strain composition, aligned by position with `initial_counts` |
| `population_size` | Number of individual hosts in each host generation |
| `infection_bottleneck` | Number of symbiont cells that initially infect each host |
| `carrying_capacity` | Adult symbiont population size within one host |
| `growth_factor` | Multiplicative population increase per bacterial generation during host colonisation |
| `steady_generations` | Bacterial generations spent at the within-host carrying capacity |
| `host_generations` | Number of complete host-passage cycles to simulate |
| `escape_fraction` | Fraction of the adult symbiont population released by each host |
| `mutation_probability` | Probability of a strain-changing mutation per cell per bacterial generation |
| `within_host_selection` | Whether within-host reproductive success depends on within-host fitness |
| `free_living_selection` | Whether the model applies one environmental generation with free-living fitness effects |

Two units deserve particular attention:

- `escape_fraction = 0.01` means that 1% of the adult symbionts are released.
- `mutation_probability` is a probability per whole genome, or per explicitly
  defined target region, per bacterial generation. A per-site mutation rate must
  first be converted using the length of the sequence being represented.

A **within-host bacterial generation** is one reproductive step inside a host.
When free-living selection is enabled, the model also adds one free-living
bacterial generation after migration. A **host generation** is one complete
passage from environmental infection through host release and environmental
updating.

## Understand the output files

The main result tables are ordinary CSV files:

| File | Contents |
| --- | --- |
| `host_generation_summary.csv` | Environmental and mean host diversity for every replicate and host generation |
| `environment_counts.csv` | Labelled environmental counts, frequencies and fitness; every generation in `all` mode or only each replicate endpoint in `final` mode |
| `infection_counts.csv` | Strains that founded every individual host infection |
| `host_adult_summaries.csv` | Diversity and mean fitness of the adult symbionts in every host |
| `host_adult_counts.csv` | Complete adult strain counts when panel or full retention is requested |
| `pooled_host_counts_and_occupancy.csv` | Abundance of each strain across all hosts and the number of hosts containing it |
| `release_counts.csv` | Strains released by each host before environmental mixing |
| `migration_counts.csv` | Realized emigrant and immigrant counts for each strain and host generation |
| `strain_origins.csv` | Which initial strains began in the focal reservoir, fixed regional source, or both |
| `strain_lineage_events.csv` | Parentage, origin and fitness of every newly mutated strain |
| `resolved_config.json` | Exact parameter values interpreted by the software |
| `provenance.json` | Seed, software version and computing information needed for reproduction |
| `final_environment_repNNN.npz` | Final environmental strain state for each replicate |
| `completion.json` | Verification record containing output sizes and final-state checksums |

The presence and size of `host_adult_counts.csv` are controlled by
`host_counts_mode`: `summary` omits detailed adult counts, `panel` retains a
reproducible sample of hosts, and `full` retains every host. Use `panel` for most
large experiments to keep data volumes manageable.

`expected_host_feedback_after_migration` in the generation summary is the
expected fraction of host-derived bacteria remaining after exchange. It is not
an observed ancestry fraction because individual-cell provenance is not tracked
through neutral capacity regulation and migration.

While a run is active, `checkpoints/` contains the two newest validated recovery
points. `checkpoint_interval = "1h"` requests a checkpoint at the first completed
host-generation boundary after one hour has elapsed; it cannot interrupt a host
generation partway through. At every replicate boundary a checkpoint is written
regardless of elapsed time. After successful completion, the software verifies
the outputs and final states and deletes all recovery checkpoints. Failed or
interrupted runs retain them for `--resume`.

## Generate a pilot report

A completed pilot matrix can be converted into a self-contained biological PDF
with `trophosome report`. The command also creates an editable Markdown copy and
its figures. See the [pilot reporting tutorial](docs/pilot-reporting-tutorial.md)
for the standardized input tables, a complete example and instructions for
adding future experimental cells.

The maintained Phase 1 server workflows, including checkpoint/resume behavior
and automatic or stand-alone reports, are documented in the
[HPC workflow guide](scripts/hpc/README.md). The current second pilot uses six
sentinel conditions, 12 matched populations per condition, and 250 repeated
host passages.

## Model versions and scientific documentation

The command `trophosome run` uses `wright_fisher_counts`, the maintained main
model. It follows a haploid, discrete-generation Wright--Fisher process over
labelled strain counts. It supports independent within-host and free-living
fitness and selection without representing every bacterial cell as a separate
software object.

Older model versions remain in the repository for comparison:

- V3.1 is a legacy neutral endpoint approximation based on an Ewens partition.
- V1.3/V1.4 and V2.1/V2.2 are retained as historical comparison models in
  `src/project_package/`.

Before interpreting research results, consult the [formal model
specification](docs/model-specification.md), [biologist-facing distributional
validation](docs/distributional-validation.md), [model
semantics](docs/model.md), [scalability analysis](docs/scalability.md), and [V3.1
evaluation](docs/v3_1-evaluation.md).

For software development, tests and optional plotting or coalescent tools,
install the additional dependencies and run the test suite:

```bash
python -m pip install -e '.[dev,coalescent,plot]'
python -m pytest
```

## Repository layout

```text
configs/                 Versioned simulation inputs
docs/                    Model, scalability, and development documentation
experiments/work/        Versioned designs, registries, manifests and analyses
experiments/data/        Curated, verified experiment releases
scripts/                 Reproducible benchmarks and maintenance checks
src/trophosome/          Main exact-count model and tested support code
src/project_package/     Legacy V1--V3.1 comparison implementations
tests/                   Automated behavioral and regression tests
notebooks/               Exploratory analyses; not the source of tested logic
examples/                Historical example outputs (legacy)
```

## Model status

The main prototype is exact for a haploid, discrete-generation, Wright--Fisher,
infinite-alleles model with the configured habitat-specific selection and
demography.
It stores only extant strain counts during reproduction, records mutation
parentage separately, and streams host results in bounded batches. Census size
is an integer; runtime depends mainly on bacterial transitions, extant richness,
and materialised mutation events.

The focal environment is a non-depleting infection reservoir with effective
capacity `capacity_ratio * carrying_capacity`. Host releases are mixed into it
and it is returned to capacity by neutral Hamilton regulation. When migration is
enabled, a fixed number of focal cells emigrate and the same number are replaced
by immigrants sampled from an unchanged regional source. Optional free-living
Wright--Fisher selection acts after this exchange. The realised host-feedback
and migration fractions are reported each generation.

The current research priorities are:

1. design and run the neutral Phase 1 host-feedback experiments;
2. calibrate mutation and both habitat-specific fitness traits biologically;
3. use the dual-selection Phase 2 profile to examine concordant and antagonistic
   fitness effects;
4. develop a separately named common-strain/rare-mutant hybrid only where the
   exact mutation-event workload is infeasible;
5. retain V3.1 as a documented neutral endpoint comparator.

## Reproducibility and contribution

Use one TOML file per experiment, never edit parameters inside a notebook, and
keep raw results outside Git. The full workflow is in
[docs/reproducibility.md](docs/reproducibility.md); development conventions are in
[docs/development.md](docs/development.md). Bounded profiling evidence is recorded
in [docs/benchmark-results.md](docs/benchmark-results.md).

No software license has been declared yet. Add one before public redistribution.
The local release preparation checklist is in
[docs/release-checklist.md](docs/release-checklist.md).
