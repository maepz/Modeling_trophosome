# Scalability analysis and model alternatives

## Executive finding

The maintained prototype is the exact count-based Wright--Fisher model. It keeps
the forward bacterial generations required for mutation timing, drift and future
within-host selection, while removing the legacy NetworkX and trajectory-storage
costs. Its dominant workload variables are:

- `H`: number of hosts;
- `T`: bacterial transitions per host;
- `R_t`: extant labelled strain richness at transition `t`;
- `M_t`: newly materialized mutation events, approximately `N_t * mu`;
- the requested strain-level output volume.

Adult census size `N_t` is stored as integer counts, so it does not by itself
create one object per cell. It still affects the number of materialized mutation
events. The exact infinite-alleles model therefore remains practical when
`N_t*mu` and output richness are manageable, not for arbitrary census and
mutation-rate combinations.

The exact-preserving implementation improvements are compact count/fitness
arrays, capacity-managed within-host buffers, one vectorized multinomial and
binomial transition per generation, a precomputed growth schedule, endpoint-only
host retention, separate append-only lineage events, inverse-CDF reservoir
founder sampling, one-pass population aggregation, deterministic per-host random
streams and ID blocks, a persistent worker pool, and a bounded FIFO submission
window. These change runtime and memory behavior without changing the declared
transition distribution.

V3.1 remains useful as a legacy approximate neutral endpoint comparator. It
jumps from infection to the adult endpoint. Its immediate mutation workload is

```text
M = K * (1 - (1-u)^T), approximately K*u*T at small u*T,
```

where `K` is adult carrying capacity, `u` is mutation probability per bacterial
generation, and `T` includes growth and adult generations. The current V3.1 CRP
visits approximately `M` inferred mutation-derived cells. It then calculates
founder probabilities with exact Stirling integers and builds a NetworkX graph
whose size depends on represented strain richness.

Its dominant variables are:

- `H`: number of hosts;
- `M`: mutation-derived cells processed per host;
- `S`: represented strains and ancestry nodes per host;
- the number of complete graphs retained, merged and serialized.

In the tested `K=10^6`, `u=10^-3` life history, `M` was about 107,793 while
expected Ewens richness was only about 4,690. This indicates that an efficient
Ewens sampler can remove the present `M`-by-current-richness probability work
without changing V3.1's declared neutral endpoint distribution. Some exact
samplers still require roughly `M` inexpensive draws. Output must also scale with
whatever strain-level detail is actually requested.

The V1 legacy implementation keeps a NetworkX graph for every sampled bacterial
generation of every host, repeatedly copies those graphs, returns every history
from worker processes, merges graphs, and pickles metapopulation time series. Its
memory is therefore approximately proportional to `H * sum_t(R_t)` plus Python
object and multiprocessing serialization overhead. More workers can increase peak
memory even when wall time improves slightly.

## Evidence from the repository

The complete V1.4/V2.2/V3.1 scheduler matrix showed that V3.1 was faster for the
tested low-mutation defaults:

| Version | Full matrix time | Peak resident memory |
|---|---:|---:|
| V3.1 | 9 min 50 s | 15.96 GB |
| V1.4 | 39 min 58 s | 67.31 GB |
| V2.2 | 1 h 44 min 46 s | 7.64 GB |

Those jobs used `K=10,000`, `u=10^-12`, zero adult steady generations, 20 host
generations and 100 initial strains. They establish V3.1's low-mutation speed
advantage but do not exercise its high-mutation bottleneck.

A local one-host profile used infection size 10, `K=10^6`, growth factor 1.2,
50 adult generations and `u=10^-4`. About 11,336 mutation-derived cells were
inferred; the run took about 1.9 seconds. Approximately 0.98 seconds was spent in
the cell-by-cell CRP, 0.51 seconds calculating exact Stirling founder
probabilities, and 0.14 seconds converting to NetworkX. At `u=10^-3`, more than
100,000 mutation-derived cells were inferred and the run was still calculating
exact Stirling numbers after 90 seconds.

The older V1 investigation shows a separate full-trajectory failure mode:

Recorded cluster runs show the failure mode directly:

| Scenario | Peak memory | Outcome |
|---|---:|---|
| 10,000 hosts, 2 workers | 30.5 GB | stopped after ~30 min |
| 10,000 hosts, 10 workers | 62.8 GB | stopped after ~20 min |
| 100,000 hosts, 10 workers | 639.6 GB | `MemoryError` during result transfer |
| 100,000 hosts, large biological test | 1.81 TB | killed after 2.5 h |

The bounded local profile reproduced the mechanism. With one host, carrying
capacity 10,000, 50 adult generations, and mutation probability 0.001, the legacy
code stored 92 graph snapshots. Graph copying, `remove_empty_leaves`, and graph
node/edge construction dominated the profile. Raising carrying capacity from
10,000 to 100,000 at mutation probability 0.001 had nearly the same cost as
raising mutation probability from 0.001 to 0.01 at capacity 10,000: both increase
`N*mu` tenfold and took about 6.2 seconds with about 25 MB traced peak memory for a
single host.

## Specific bottlenecks and correctness risks

### V3.1

1. The current CRP rebuilds a probability array over all existing tables for
   every mutation-derived cell. It performs much more work than is required to
   sample the Ewens partition.
2. Founder occupancy probabilities construct enormous exact Stirling integers.
   This becomes a major high-richness runtime bottleneck even though the infection
   bottleneck itself is small.
3. Zero mutation asks `msprime` to simulate zero samples and crashes.
4. Rounding each original strain down and calling the remainder mutant creates
   spurious mutation-derived cells at tiny mutation probabilities.
5. The CRP, fitness propagation and host workers do not share a controlled seed
   hierarchy, so a supplied seed does not reproduce a complete run.
6. `theta=K*u`, use of census `K` as the relevant population size, and use of an
   equilibrium Ewens partition after bottleneck/growth are scientific assumptions
   requiring calibration.
7. Fitness is assigned after adult abundance is decided and therefore does not
   produce within-host selection on abundance.
8. Complete per-host graphs are still merged and full host-generation graph time
   series are pickled. Output size can replace bacterial generations as the
   memory limit.
9. Infection and environmental depletion retain the legacy mixed reservoir/
   finite-population semantics.

### V1/V2 comparison paths

1. `run_generation_of_host_pop` creates one argument dictionary per host and
   returns every host's full `results` history through `multiprocessing.Pool`.
2. The outer `sampling_rate` parameter is ignored; workers receive
   `sampling_rate=1`, maximizing retained snapshots.
3. `run_until_fixation3` copies the graph after every generation and stores those
   copies in a dictionary.
4. NetworkX stores Python dictionaries per node and edge. It is useful for final
   visualization but expensive as a simulation state.
5. Mutation supply scales as `N*mu`. At `N=1e10` and `mu=1e-4`, the expectation is
   one million new mutant cells per host per bacterial generation. Across 100,000
   hosts, explicit infinite-alleles materialization is not computationally or
   scientifically summarizable as ordinary graph objects.
6. Legacy mutation counts are drawn from the parental size `n`, not the offspring
   size. During population decline this can request more mutants than offspring.
7. Growth stops only after exceeding the threshold, so carrying capacity is
   overshot.
8. Mutation identifiers are restarted between growth and steady phases because
   the updated `new_avail_id` is not returned, allowing identifier collisions.
9. Random seeds are not reproducible: `np.random.seed()` is called without a seed
   in workers; `np.random.seed(...)` returns `None`, which is passed to msprime;
   and the command-line `random_seed` option does not seed a supplied value.
10. V2.2 correctly sets `ploidy=1`, improving on the earlier V2.1/default
    diploidy problem, but still creates samples proportional to adult census
    size. A tree sequence compresses shared ancestry, but it cannot make billions
    of requested present-day samples cheap.
11. Infection samples each host independently from one pool and then subtracts the
    combined samples. This can create negative environmental abundance.
12. A full Conda environment export is committed as `requirements.lock.txt`; it
    contains machine-specific build paths and unrelated software, so it is not a
    portable lock file.

## Alternatives, with scientific status

| Approach | Scientific status | Preserves | Loses/approximates | Best use |
|---|---|---|---|---|
| V3.1 endpoint + efficient Ewens sampler | Same approximate neutral endpoint as V3.1 if implemented with distributional-equivalence tests | bottleneck, endpoint mutation partition, bounded strain ancestry | generation-by-generation drift and selection | legacy neutral comparison; possible future separately named endpoint model |
| Streamed genotype counts (maintained main prototype) | Exact for declared haploid WF infinite-alleles model | bottleneck, mutation timing, drift, selection, growth, adult phase, escape | cell genealogy and genomes | Phase 1 inference at manageable `N*mu`; legacy validation |
| Sample-first DTWF/coalescent | Exact genealogy for a declared neutral demographic model, or Kingman approximation where valid | bottleneck and demography for released/sequenced lineages | forward selection unless extended | molecular diversity of a bounded output sample |
| Common types + rare-mutant branching process | Approximate stochastic hybrid | selection/drift for common types; establishment/extinction of rare types | interactions among many rare clones | large `N*mu` with selected mutations |
| Wright--Fisher diffusion / Dirichlet-multinomial transition | Approximate aggregate | allele-frequency moments and bottleneck drift | individual mutation genealogy, strong sweep detail | parameter sweeps and sensitivity analysis |
| Host-composition cohorts | Exact only when hosts with identical state share identical random outcomes; otherwise approximate | host heterogeneity categories | within-cohort stochastic independence | very large host populations with few host states |
| SLiM with tree-sequence recording | Exact individual simulation at simulated scale | flexible selection, demography, migration | still requires simulated individuals; scaling changes drift | validation at intermediate scale |
| Population rescaling | Approximate | selected dimensionless quantities if carefully matched | rare-event timing, clonal interference, AFS tails | sensitivity checks, never unvalidated default |

For rescaling census size by factor `C`, matching neutral diversity often requires
increasing mutation rate by `C`; matching `N*s` may require scaling selection too.
Those transformations are not universally valid during bottlenecks, strong
selection, linkage, or clonal interference. Validate against unscaled small cases.

## Recommended architecture

Keep three representations separate:

1. simulation state: compact arrays of extant genotype IDs, counts, and fitness;
2. optional lineage log/tree sequence: mutation events or sampled ancestry only;
3. analysis output: tidy summaries plus explicitly requested samples.

Never use the visualization graph as the mutable simulation state. Construct a
small graph after a run from a sampled lineage output.

## Implementation roadmap

### Phase 0: optimize and restore the exact prototype — implemented

- Replace NetworkX simulation state with compact count and fitness arrays.
- Preserve the exact Wright--Fisher reproduction, mutation and selection draws.
- Stream endpoints and outputs in bounded host batches.
- Use deterministic child random streams and strain-ID reservations.
- Add population-size, mutation-lineage and worker-invariance tests.

### Phase 1: establish a scientific baseline

- Freeze a small parameter matrix and expected summary distributions.
- Define whole-genome versus per-site mutation semantics.
- Use the non-depleting effective-reservoir assumption with configurable
  `capacity_ratio`.
- Treat escape as a physical without-replacement sample from each adult host.
- Match haploid ploidy and the within-host demographic trajectory.

### Phase 2: benchmark and validate the maintained exact prototype

- Compare the neutral exact prototype with V1.3/V1.4 across tractable bottleneck,
  capacity, adult-duration and mutation-supply grids.
- Compare replicate distributions of richness, diversity, founder survival and
  escapee composition, not individual random trajectories.
- Benchmark runtime across worker count, host count, `K`, mutation rate and
  represented richness.

### Phase 3: evaluate a separately named neutral endpoint approximation

- Replace the current CRP with an efficient Ewens partition sampler such as a
  Feller-coupling or allele-count implementation.
- Replace exact Stirling integers with a numerically stable small-state founder
  occupancy sampler.
- Add distributional-equivalence tests before changing the default.

### Phase 4: extend bounded outputs and restart support — partly implemented

- Keep full adult spectra only for validation or a reproducibly sampled panel.
- Validated restart from atomic host-generation checkpoints is implemented,
  including output truncation, source/configuration validation, two-checkpoint
  retention, and successful-run cleanup.
- Add columnar output after the table schema is frozen.

### Phase 5: deepen exact-model validation

- Compare distributions of richness, heterozygosity, mean fitness, escapee
  composition, and fixation time against matched small legacy runs.
- Add neutral analytical checks and repeated stochastic tests.
- Benchmark across host count, richness, `K`, `mu`, and worker count.

### Phase 6: implement sample-first neutral genealogy when the output is a sample

- Simulate ancestry only for released or sequenced cells.
- Use DTWF during severe recent bottlenecks/growth and coalescent farther back.
- Set `ploidy=1`, use explicit demography, and add mutations after ancestry.

### Phase 7: add a separately named selection-aware hybrid

- Track common selected genotypes forward as counts.
- Treat new rare mutants with Poisson mutation supply and branching-process
  establishment probabilities.
- Materialize only established or sampled lineages.

### Phase 8: biological calibration

- Perform simulation-based calibration and sensitivity analysis.
- Report which results are stable across exact, coalescent, and hybrid methods.
- Use the simplest validated method for each scientific question rather than one
  universal simulator.
