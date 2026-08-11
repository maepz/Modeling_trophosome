# Evaluation of model V3.1

## Short conclusion

V3.1 is now retained as a legacy neutral endpoint comparator. It was a
substantial computational improvement over the original V1.4 and V2.2 paths
because it does not reproduce every bacterial generation inside every host.
Instead, it jumps from the infection bottleneck to an expected adult population
and constructs only the mutant strain partition and a strain-level ancestry.

That is the right general direction for reaching large bacterial census sizes.
However, the present V3.1 implementation is not yet ready for biological
inference at large mutation supply. Three correctness problems should be repaired
first: zero mutation can crash, rounding creates mutant cells even when the
expected mutation supply is almost zero, and a supplied random seed does not make
the full result reproducible. At higher mutation supply, the current Chinese
restaurant process (CRP) and the exact Stirling-number calculation become the
main runtime bottlenecks.

This evaluation did not change the scientific behaviour of the legacy V3.1
function. It added documentation, tests for a separate workload estimator, and a
roadmap for repairing and validating V3.1 before extending it.

## Where the versions now fit

| Version | Main idea | Role in the project |
|---|---|---|
| V1.3 | Forward Wright--Fisher process, with within-host selection | Legacy scientific reference |
| V1.4 | Forward Wright--Fisher process, with selection switchable | Main small-scale neutral comparator |
| V2.2 | Coalescent/tree-sequence approach with haploid ploidy and within-host growth | Alternative ancestry comparator |
| V3.1 | Direct endpoint approximation plus an Ewens/CRP strain partition and a coalescent strain tree | Legacy neutral endpoint comparator |
| `trophosome/count_model.py` | Tested Wright--Fisher process over genotype counts | **Maintained main prototype**, closest in law to V1.3/V1.4 |

The command-line program runs the tested count model. It does not silently
replace or run V3.1. That boundary is deliberate because V3.1 remains an
approximate neutral endpoint with a different transition law.

## What V3.1 does for one host

In biological language, V3.1 performs the following operations:

1. It draws the infecting bacterial cells from the free-living population.
2. It calculates how many bacterial generations are needed to grow from the
   infection bottleneck to the adult carrying capacity.
3. It adds the number of generations spent at adult carrying capacity.
4. It calculates the probability that a lineage has no mutation during that
   whole interval.
5. It assigns the corresponding expected abundance to the strains that infected
   the host.
6. It treats all remaining adult cells as belonging to mutation-derived strains.
7. It divides those cells among new strains using an Ewens/Chinese restaurant
   process.
8. It generates a coalescent tree with one sampled leaf per new strain, cuts the
   tree to a sampled number of founder lineages, attaches it to the infection
   strains, and adds fitness values along the branches.
9. It samples the bacteria that escape from the adult host.

This is not simply V1.4 written more efficiently. V1.4 draws reproduction,
genetic drift, mutation and possibly selection in each bacterial generation.
V3.1 replaces most of that trajectory with an adult endpoint distribution.

## The main V3.1 quantities

Let:

- `B` be the infection bottleneck;
- `K` be adult carrying capacity;
- `g` be the multiplicative growth factor;
- `A` be the number of generations at adult carrying capacity;
- `u` be the mutation probability per bacterial cell generation.

V3.1 calculates the growth duration as

```text
T_growth = ceiling(log(K / B) / log(g))
T_total  = T_growth + A
```

and the probability that a lineage experiences at least one mutation as

```text
p_mutation = 1 - (1 - u)^T_total.
```

The approximate number of mutation-derived adult cells is therefore

```text
M = K * p_mutation.
```

When `u*T_total` is small, this is approximately `K*u*T_total`. This value is a
much better screen for V3.1 workload than `K` alone. The legacy implementation
then runs its CRP approximately `M` times.

V3.1 sets the Ewens parameter to `theta = K*u`. This evaluation records that as
the current assumption; it does not validate it. The correct relationship among
`theta`, mutation probability and effective population size depends on the
chosen haploid population-genetic model. In particular, adult census size `K`
must not automatically be assumed to equal effective population size.

## What is exact and what is approximate

| Component | Status in V3.1 | Interpretation |
|---|---|---|
| Infection draw | Exact stochastic draw for the coded multinomial reservoir assumption | The infection bottleneck is explicitly random |
| Growth duration | Deterministic calculation | Every host uses the same declared growth rule |
| Abundance of original strains | Expected value, then rounded down | Within-host drift and selection do not change these abundances |
| Total mutation-derived abundance | Derived as the remainder up to `K` | It includes accumulated rounding remainder |
| Partition among new strains | Stochastic Ewens/CRP draw conditional on the coded `M` and `theta` | A neutral equilibrium partition is being used as an endpoint approximation |
| New-strain tree | Stochastic Kingman coalescent tree over one leaf per represented strain | It is not the full cell genealogy through the actual growth history |
| Number of infection founders represented in the new-strain tree | Stochastic occupancy calculation | Current implementation uses exact Stirling integers |
| Fitness labels on the tree | Stochastic post-processing | Fitness does not affect adult abundance in V3.1 |
| Escape sample | Stochastic sample from the adult graph | Its exact meaning depends on the legacy sampling function |

The Ewens sampling formula is a neutral result. Applying it after recent growth,
a severe infection bottleneck and a finite adult phase is a modelling
approximation that needs comparison against matched V1.4 simulations. Efficient
ways to sample this same neutral distribution exist; changing the sampler need
not change the declared distribution. See Tavaré's review of the Ewens sampling
formula, the Chinese restaurant process and Feller coupling:
<https://londmathsoc.onlinelibrary.wiley.com/doi/10.1112/blms.12537>.

## Correctness and scientific issues found

### 1. A zero-mutation run can crash

When `u=0`, V3.1 produces zero new strains and then asks `msprime` for an ancestry
with zero samples. The resulting `ValueError` means the model cannot currently
represent the biologically valid no-mutation control.

**Recommended repair:** return an adult population containing only the original
infection strains when the number of mutation-derived cells is zero. Add a test
that verifies population size and escape sampling in this case.

### 2. Integer rounding creates spurious mutation-derived cells

The expected adult abundance of every original strain is individually rounded
down. All discarded fractional parts are then classified as mutation-derived.
For example, with `K=1,000` and `u=10^-12`, the biological expectation is only
about `10^-8` mutation-derived cells for the tested life history, but the code
created a new genotype with abundance two. The effect becomes larger when more
founder strains are present because each strain contributes another discarded
fraction.

**Recommended repair:** use a probability-conserving stochastic allocation. One
option is to draw the total number of unmutated cells first and then allocate
those cells among original strains with a multinomial draw. The remaining cells
are mutation-derived by construction. This preserves `K` without converting
rounding error into mutations.

### 3. A random seed does not reproduce the full run

The host function creates a generator from `seed`, but calls the nested CRP with
`rng=None`. The fitness-propagation function is also called without a seed.
Repeated runs with the same supplied seed therefore produced different strain
partitions and fitness values.

At host-population level every worker is currently given `seed=None`, so there is
no user-controlled reproducibility at all.

**Recommended repair:** create one master `SeedSequence`, derive a distinct child
stream for each host and each stochastic stage, and pass generators or explicit
integer seeds to infection, partition, ancestry, fitness and escape sampling.
The result should be reproducible independently of worker count.

### 4. The endpoint abundance model is neutral and deterministic

Original strains retain their infection frequencies in expectation. Fitness is
added only after their adult abundances have been decided. Consequently, V3.1
does not currently model within-host selection changing strain frequencies, even
though fitness values appear in its final graph. It also omits within-host drift
among original strains.

This is acceptable only if V3.1 is explicitly labelled as a neutral endpoint
model and its output questions are chosen accordingly. It is not equivalent to a
selected V1.3 run.

### 5. Fitness assumptions differ across versions

The V3.1 tree helper uses mutation effects with mean `-0.1`, standard deviation
`0.1`, and multiplicative propagation. Some V1 code uses different effect sizes
and additive updates. A performance comparison is still possible, but a
scientific comparison of mean fitness or selection is not matched unless these
assumptions are harmonised.

### 6. The Ewens parameter and endpoint approximation need validation

Two questions should be treated separately:

1. Is the coded `theta = K*u` the intended haploid neutral parameter, and should
   `K` be replaced by an effective population size?
2. Even with a calibrated `theta`, is an equilibrium Ewens partition a good
   approximation after the infection bottleneck, growth and the declared adult
   duration?

These are scientific calibration questions, not code-cleaning questions. They
should be answered with small matched simulations and analytical checks before
large V3.1 parameter sweeps.

### 7. Environmental accounting still has legacy risks

Hosts sample infection independently from one environmental graph, but the
aggregate infected counts are later subtracted. This mixes a reservoir-like
infection draw with finite-pool depletion and can subtract the wrong counts when
the set or ordering of strains differs. The project needs one explicit choice:

- a reservoir, where infection does not deplete the represented environment; or
- a finite environmental population, where all host infections are sampled
  jointly without replacement and do deplete it.

### 8. Complete graph output can become the next memory limit

V3.1 avoids storing every within-host bacterial generation, but still constructs
a NetworkX ancestry graph for each host, merges graphs and pickles full graph time
series for every host generation in the current job script. If the number of
represented strains is large, the graph and output files will dominate memory.

No implementation can store a complete strain graph in less space than the graph
itself. Large simulations therefore need bounded outputs: summaries for all
hosts, plus complete ancestry only for a declared subset of hosts or sampled
escapees.

## Performance evidence

### Existing full-matrix cluster runs

The recorded scheduler comparison covered multiple host counts and worker counts
under the job-script defaults: `K=10,000`, `u=10^-12`, zero adult steady
generations, 20 host generations and 100 starting strains.

| Version | Full matrix elapsed time | Peak resident memory |
|---|---:|---:|
| V3.1 | 9 min 50 s | 15.96 GB |
| V1.4 | 39 min 58 s | 67.31 GB |
| V2.2 | 1 h 44 min 46 s | 7.64 GB |

This supports the conclusion that V3.1's endpoint design is faster than the two
comparison paths for that matrix. It does **not** test mutation scaling: at
`u=10^-12`, almost no mutation-derived cells are expected. One earlier V3.1
scheduler run also failed because the function returned an undefined `results`
variable; the subsequently pulled code no longer has that return-value error.

### Local, mutation-focused profile

A bounded one-host profile used infection size 10, `K=10^6`, growth factor 1.2,
50 adult generations and `u=10^-4`. This implies 114 total within-host lineage
generations and approximately 11,336 mutation-derived cells. It completed in
about 1.9 seconds and produced a graph with roughly 900 nodes. The largest costs
were:

- the cell-by-cell CRP: about 0.98 seconds;
- exact Stirling-number founder probabilities: about 0.51 seconds;
- conversion to NetworkX: about 0.14 seconds.

At `u=10^-3`, the same life history implies approximately 107,793
mutation-derived cells and an expected Ewens richness of about 4,690 strains. The
legacy function was still computing exact Stirling numbers after 90 seconds and
was stopped. This confirms that high mutation supply, rather than carrying
capacity alone, is V3.1's immediate limit.

Use the new estimator before a run:

```bash
PYTHONPATH=src python scripts/estimate_v3_1_load.py \
  --infection-size 10 \
  --carrying-capacity 1000000 \
  --growth-factor 1.2 \
  --steady-generations 50 \
  --mutation-probability 0.001
```

The estimator does not simulate biology. It reports the legacy V3.1 quantities
and warns when the present implementation is likely to be expensive.

## Recommended implementation roadmap

### Stage 1: freeze and test the current V3.1 meaning

Before optimising, write a compact V3.1 specification and create small tests for:

- zero mutation;
- exactly one infection strain;
- several infection strains;
- zero adult generations;
- fixed population-size conservation at `K`;
- escape count conservation;
- identical results from identical seeds;
- statistically independent hosts from different child seeds.

This stage should not add selection or a new population-genetic approximation.

### Stage 2: repair the three correctness problems

Implement the zero-mutation branch, probability-conserving abundance allocation,
and complete random-stream propagation. Also validate growth inputs such as
`g>1` when `K>B` and reject impossible probabilities early.

These are code-level repairs, but the new stochastic abundance allocation must
still be compared with the old expected values over many replicates.

### Stage 3: validate the neutral endpoint approximation

Run V1.4 and V3.1 on a small grid where V1.4 is tractable. Use at least 100--1,000
replicates per point and compare distributions, not one trajectory.

Suggested starting grid:

- infection bottleneck: 2, 10 and 20 cells;
- carrying capacity: 100, 1,000 and 10,000 cells;
- adult duration: 0, 10 and 100 generations;
- mutation supply selected to cover `K*u*T` below 0.1, around 1, around 10 and
  around 100;
- neutral fitness in both models.

Compare adult and escapee strain richness, Simpson diversity, frequency of the
founder strains, singleton fraction, pairwise diversity and the probability that
one infection founder dominates. Predefine an acceptable discrepancy for each
summary before using V3.1 outside the tested region.

### Stage 4: replace computationally expensive calculations without changing
the neutral target distribution

1. Replace the current cell-by-cell CRP implementation with an efficient Ewens
   partition sampler, for example one based on Feller coupling or an equivalent
   allele-count formulation.
2. Replace enormous exact Stirling integers with a numerically stable founder
   occupancy sampler. Because the infection bottleneck is small, a small-state
   transition calculation or direct occupancy simulation is preferable.
3. Keep abundance arrays as the simulation state. Construct NetworkX only for
   the subset of ancestry that will actually be inspected or plotted.

This stage can be distribution-preserving under the declared Ewens endpoint
model. It is distinct from replacing V3.1 with a new approximation.

An exact computational improvement is suitable when the number of
mutation-derived cells is large but still iterable. Some exact samplers still
perform approximately one inexpensive operation per such cell. If that number
reaches millions or billions per host, use a sample-first or aggregate output
rather than claiming that an exact full partition is independent of it.

### Stage 5: bound outputs for thousands of hosts

For every host, retain a compact row of summary statistics and the released-cell
counts needed for the environmental update. Store full within-host ancestry only
for a reproducibly sampled panel of hosts or cells. Write host results in chunks
rather than returning every graph from all workers at once. Checkpoint at each
host generation.

### Stage 6: add selection only as a separately named, validated model

If the biological question requires selection to affect abundance, develop a
V3.2 or selection-aware hybrid rather than making the current neutral shortcut
appear exact. Plausible options are:

- forward counts for common selected strains plus a branching-process treatment
  of rare new mutants;
- a bounded sample-first genealogy for neutral molecular diversity, combined
  with forward trajectories for selected lineages;
- a carefully tested aggregate transition for parameter sweeps.

Each option is approximate in a different way. Validate it against V1.3/V1.4 at
small scale and report the region in which its summary statistics agree.

## Practical decision rule

V3.1 can be used now for software exploration and neutral small-scale comparison.
It should not yet be used for final large-scale biological conclusions. Promote
it to a production neutral model only when:

1. the zero-mutation, rounding and seed tests pass;
2. the Ewens parameter is scientifically defined;
3. V3.1 agrees with the chosen V1.4 summary distributions over a declared
   parameter region;
4. the faster partition and founder samplers pass distributional equivalence
   tests; and
5. output memory remains bounded as host number increases.
