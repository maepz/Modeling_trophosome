# Exact-count model specification

## Identity and status

- Model family: `wright_fisher_counts`
- Scientific specification: `2.1.0`
- Reference software release: `0.7.0`
- Output schema: `2.3.0`
- Status: release candidate, distributionally validated

This document defines the stochastic process. An implementation change that
violates it is a new model specification, not a performance optimization.

## State

At environmental host generation `t`, the focal effective reservoir is
represented by labelled strain counts:

```text
E_t = {(strain_id_i, count_i, w_H,i, w_E,i)}.
```

Its fixed effective capacity is

```text
N_E = round(capacity_ratio * carrying_capacity).
```

Strain identifiers are unique. Counts are positive integers. `w_H` is a finite,
non-negative within-host reproductive weight and `w_E` is an independently
inherited, finite, non-negative free-living reproductive weight. Mutation
ancestry is an append-only event relation and is not part of the active
reproductive state.

When fixed-pool migration is enabled, a second, time-invariant composition is
declared:

```text
G = {(strain_id_i, regional_count_i, w_H,i, w_E,i)}.
```

`regional_count_i` supplies relative immigration weights. `G` is not a finite
regional census and is never depleted or updated. Focal and regional input
vectors are aligned by strain position. A declared strain may initially have
zero abundance in either habitat, but not in both.

## One host-population generation

### 1. Infection

Each of `H = population_size` hosts independently receives `B =
infection_bottleneck` cells.

In the Phase 1 `reservoir` mode, cells are sampled with replacement from `E_t`.
Sampling one host does not deplete or otherwise change the reservoir. The
optional `finite` mode samples without replacement and is not the Phase 1
default.

### 2. Within-host population schedule

Starting at `B`, the requested offspring-population size is repeatedly

```text
min(K, max(N + 1, ceil(N * growth_factor)))
```

until carrying capacity `K` is reached. The model then applies
`steady_generations` transitions at size `K`.

### 3. Wright–Fisher reproduction

For parental strain `i` with abundance `n_i`, reproductive probability is

```text
p_i = n_i / sum_j(n_j)                         selection disabled
p_i = n_i * w_H,i / sum_j(n_j * w_H,j)         selection enabled.
```

For requested offspring size `N'`, pre-mutation offspring counts follow

```text
(X_1, ..., X_S) ~ Multinomial(N', p).
```

Generations are discrete and non-overlapping. There is no cell survival outside
the sampled parental contribution implicit in the Wright–Fisher draw.

### 4. Mutation

Conditional on `X_i`, the number of mutated offspring of parent `i` is

```text
M_i ~ Binomial(X_i, u),
```

independently across parental strains. Each event creates one unique strain with
count one; multiple simultaneous events never share a new strain identity. This
is a Bernoulli infinite-alleles mutation model with at most one strain-changing
event per offspring per generation.

Each mutation event records its parent strain, within-host generation and
mutational depth. Every mutant receives two fitness values:

```text
w_H,new = max(fitness_floor, w_H,parent + epsilon_H)
w_E,new = max(fitness_floor, w_E,parent + epsilon_E)

epsilon_H, epsilon_E independently ~ Normal(effect_mean, effect_sd).
```

Both habitats use the same configured effect distribution, whose default mean
is `-0.01` and standard deviation is `0.01`, but draw independently. The effect
is additive: for a parent with fitness 1, the default expected mutant fitness is
0.99. New within-host fitness affects reproduction only in later within-host
bacterial generations. New free-living fitness first acts if the strain returns
to the environmental reservoir and free-living selection is enabled.

### 5. Host release

At the adult endpoint, each host releases

```text
e = round(K * escape_fraction)
```

cells sampled without replacement from that host. The pooled release size is

```text
R = H * e.
```

### 6. Environmental return and neutral capacity regulation

All released cells are added to the unchanged effective reservoir. Before
capacity regulation, realized host feedback is

```text
alpha_t = R / (N_E + R).
```

The mixed counts are returned to `N_E` by Hamilton proportional apportionment.
The largest fractional remainders receive residual cells. Exact ties at the
allocation boundary are resolved by a dedicated seeded random stream, ensuring
that strain identifiers cannot confer a rounding advantage. This regulation is
strain-neutral and occurs whether or not later migration or free-living
selection is enabled. Call the resulting focal counts `A_t`.

### 7. Optional fixed-pool migration

When `migration.mode = "none"`, set `B_t = A_t`.

When `migration.mode = "fixed_regional_pool"`, let

```text
M = round(migration.fraction * N_E)
m_realized = M / N_E.
```

Draw the emigrant counts jointly without replacement from `A_t`:

```text
U ~ MultivariateHypergeometric(A_t, M).
```

Independently draw immigrant counts with replacement from the fixed regional
composition:

```text
V ~ Multinomial(M, regional_count / sum(regional_count)).
B_t = A_t - U + V.
```

Consequently, `sum(B_t) = N_E` exactly. Both immigration and emigration occur at
the focal boundary. The regional composition `G` nevertheless remains fixed:
focal emigrants do not alter it, and immigrant sampling does not deplete it.

Because individual ancestry is not retained through Hamilton apportionment and
migration, the summary

```text
expected_host_feedback_after_migration = alpha_t * (1 - m_realized)
```

is an expectation, not an observed lineage fraction.

### 8. Optional free-living selection

When `free_living_selection = false` (the neutral Phase 1 profile), set
`E_(t+1) = B_t`.

When `free_living_selection = true` (the Phase 2 profile), one haploid
free-living Wright--Fisher transition acts on `B_t`. For strain `i`,

```text
q_i = B_t,i * w_E,i / sum_j(B_t,j * w_E,j)
(Y_1, ..., Y_S) ~ Multinomial(N_E, q).
```

This is exactly one free-living bacterial generation per host-population
generation. It combines environmental selection and genetic drift after neutral
regulation and migration. It does not create mutations: mutation remains
confined to the within-host phase. A different number or duration of
environmental generations would be a further scientific extension rather than
a scaling optimization.

The resulting state is `E_(t+1)`. It is simultaneously the post-environmental
state of generation `t` and the pre-infection state of generation `t+1`.

## Randomness and parallel execution

Random generators are derived from the master seed and logical coordinates:
replicate, host generation, host, and stochastic stage. Strain-ID blocks are
also reserved by logical host. Changing worker count or host batch size must not
change a seeded result.

## Exactness claim

The count representation is exact for the transitions above. It does not store
individual cells, but its Multinomial, Binomial and Hypergeometric draws are the
same distributions obtained from an explicit cell implementation.

The specification does not claim to represent complete genomes, finite sites,
recombination, mutation in the environment, more than one environmental
generation per host generation, a dynamic or depleted regional metapopulation,
overlapping host cohorts or within-trophosome spatial structure.

## Version-change policy

- Array reuse, batching or output streaming that preserves every transition
  distribution does not change the model-specification version.
- A backward-compatible output column addition increments the output-schema
  minor version; removing or redefining a field increments its major version.
- Changing mutation timing, reproduction, infection, release or environmental
  regulation increments the model-specification version.
- Adding an optional scientific process requires a named profile or a new model
  family when it cannot be represented as an existing parameter value.
