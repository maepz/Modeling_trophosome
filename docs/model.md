# Model semantics

## Version roles

`wright_fisher_counts`, implemented in `src/trophosome/count_model.py`, is the
main research prototype. It preserves the declared forward Wright--Fisher law
while replacing NetworkX population state with compact labelled count arrays.
It is conceptually closest to V1.3/V1.4 and carries independent within-host and
free-living fitness traits, with a separate selection switch for each habitat.

V3.1 is a legacy neutral endpoint comparator implemented by
`project_package.run_model.run_generation_of_host_pop_v3_1`. V1.4 is the legacy
small-scale neutral Wright--Fisher comparator, V1.3 is the selection-aware legacy
comparator, and V2.2 is the coalescent/tree-sequence comparator.

## Biological cycle

For each host generation, the model applies eight stages:

1. draw `infection_bottleneck` bacterial cells from the environmental pool;
2. expand the within-host population by a fixed multiplicative growth factor;
3. cap it at `carrying_capacity` and run `steady_generations` at that size;
4. sample `escape_fraction` of the adult cells without replacement;
5. add all released cells to the non-depleted effective reservoir;
6. return the mixed representation to its fixed effective capacity by neutral
   Hamilton apportionment;
7. optionally replace a fixed fraction of the focal environment with bacteria
   sampled from a fixed regional source; and
8. optionally apply one free-living fitness-weighted Wright--Fisher transition,
   then repeat.

The infection bottleneck, growth phase, adult phase, escape, and horizontal
re-acquisition therefore remain explicit.

## V3.1 endpoint model

V3.1 keeps the infection and escape samples stochastic but does not simulate
reproduction at every bacterial generation. For each host it:

1. calculates the growth time from infection size `B` to carrying capacity `K`;
2. adds the declared adult steady duration to obtain `T` lineage generations;
3. assigns the infection strains their expected unmutated adult abundances using
   `(1-u)^T`, where `u` is mutation probability per cell generation;
4. treats the remainder up to `K` as mutation-derived cells;
5. samples an Ewens/Chinese restaurant process partition of those cells using the
   legacy assumption `theta=K*u`;
6. creates a coalescent tree with one sample per represented new strain, attaches
   it to sampled infection founders, then samples escapees.

This is a neutral endpoint approximation, not an exact acceleration of V1.4.
Drift and selection do not alter the adult abundance of infection strains, and
fitness values are assigned after abundances have been decided. The equilibrium
Ewens partition also requires validation for a population that has experienced a
recent bottleneck, growth and a finite adult phase.

The expected number of mutation-derived cells is approximately
`K * (1 - (1-u)^T)`, or `K*u*T` when `u*T` is small. This quantity determines how
many iterations the current V3.1 CRP performs. See `docs/v3_1-evaluation.md` for
the full exact/approximate distinction and validation roadmap.

## Exact-count model

`src/trophosome/count_model.py` implements a haploid Wright--Fisher process over
genotype counts. In each bacterial generation:

- parental probability is proportional to abundance times within-host fitness
  when within-host selection is enabled, or abundance alone when it is disabled;
- the next generation is one multinomial draw of the requested size;
- each offspring undergoes at most one strain-changing event, independently, with
  the configured Bernoulli probability per genome per generation;
- every event creates a unique genotype (infinite-alleles assumption);
- every mutant independently draws a within-host fitness effect and a
  free-living fitness effect from the same configured normal distribution; each
  value is parental fitness plus its habitat-specific effect, floored at the
  configured non-negative value.

This is exact for those assumptions. Census cells are not Python objects. Runtime
and memory depend mainly on within-host transitions, extant genotype richness,
and the number of newly materialized mutants, rather than directly on census
size.

The optimized implementation changes representation and execution, not the
transition law. It uses mutable capacity-managed arrays inside each host,
precomputes the size schedule, retains only the adult endpoint by default, keeps
mutation parentage in an append-only event stream, samples small reservoir
bottlenecks without a full-richness multinomial, processes hosts in bounded
batches, and reuses a persistent worker pool. Random streams and strain-ID blocks
are assigned by logical host, so changing the worker count does not change a
seeded result.

## Effective environmental reservoir

The Phase 1 default is `sampling_mode="reservoir"`. Environmental counts are
sampling weights for a spatially and temporally bounded focal effective
reservoir, not a literal census of all free-living bacteria. Infection samples
do not deplete it. Host return and optional regional exchange operate on this
focal population. Its fixed capacity is

```text
N_E = capacity_ratio * carrying_capacity.
```

If `R_t` bacteria are released, the reported realized feedback before capacity
regulation is

```text
alpha_t = R_t / (N_E + R_t).
```

The mixed focal pool is first returned to `N_E` by neutral Hamilton
apportionment. With migration disabled, this post-return state passes directly
to the optional selection stage.

With `migration.mode="fixed_regional_pool"`, the model calculates

```text
M = round(migration.fraction * N_E).
```

Exactly `M` focal cells emigrate as a without-replacement sample, and exactly
`M` immigrants are sampled with replacement from `regional_counts`. The focal
capacity therefore remains `N_E`. The regional composition is a fixed sampling
distribution: it is not depleted by immigration and does not receive or retain
the focal emigrants. Thus the model represents physical immigration and
emigration at the focal boundary, but only the fixed regional pool influences
composition across generations. The sum of `regional_counts` sets no regional
census size; only its relative strain weights matter.

The focal and regional vectors share strain positions. A zero is permitted in
either vector when the strain occurs in the other, allowing focal-only,
regional-only and initially shared strains. Their shared strain IDs and fitness
metadata make later origins unambiguous.

With `free_living_selection=true`, the post-migration focal pool undergoes one
multinomial free-living transition with probabilities proportional to abundance
times free-living fitness. With selection disabled, it proceeds unchanged to
the next infection. In either case, the stored environmental state after
generation `t` is the `pre_infection` state of generation `t+1`. The optional
`finite` infection mode is retained for explicit sensitivity experiments and
samples without replacement, but it is not the agreed Phase 1 reservoir
assumption.

## Important parameter semantics

`mutation_probability` is a whole-genome (or explicitly defined target-region)
probability per bacterial cell generation. A per-site rate must not be inserted
without conversion. If the target length is `L` and the per-site rate is `u`, a
small-rate approximation is `L*u`; the probability of at least one mutation under
independent sites is `1 - (1-u)^L`. The current Bernoulli kernel records that as
one strain-changing event, so multiple-hit distances require a future Poisson or
finite-sites extension when `L*u` is not small.

`escape_fraction` is a fraction in `[0, 1]`, not a percentage. `0.01` means 1%.

Each strain has two relative reproductive weights:
`within_host_fitness` and `free_living_fitness`. They are inherited independently
at mutation: each equals its corresponding parental value plus an independent
draw from `Normal(mutation_effect_mean, mutation_effect_sd)`, then is floored at
`fitness_floor`. The default mean effect of `-0.01` is additive, so it means a
one-percentage-point reduction when parental fitness is one; it is not a 1%
multiplicative reduction at every parental value.

`within_host_selection` and `free_living_selection` activate the two fitness
dimensions independently. Free-living selection represents exactly one
free-living Wright--Fisher generation per host-population generation, after host
return, neutral regulation and migration, but before the next infection sample.
It creates drift and selection but no environmental mutation.

## Not yet represented

- complete genomes, finite sites, recombination, or horizontal gene transfer;
- cell-level genealogy or molecular branch lengths; the lineage table records
  strain-changing mutation events, not every cell birth;
- dynamic regional demes, regional depletion, or feedback of focal emigrants
  into the regional composition;
- trophosome spatial subdivision or lobe-specific migration;
- overlapping host cohorts, variable host lifespan, or host demography;
- mutation in the free-living phase or a configurable number of free-living
  generations per host generation;
- density-dependent growth other than an imposed carrying-capacity trajectory.

These are model extensions, not software optimizations, and require explicit
scientific decisions and validation tests.
