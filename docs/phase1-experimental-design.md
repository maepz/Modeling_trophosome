# Phase 1 experimental design: neutral host feedback and environmental symbiont diversity

**Status:** proposed design for review. No production simulations have been
launched.

**Model:** exact-count `wright_fisher_counts` prototype, model specification
2.0.0 and output schema 2.2.0.

## Purpose

Phase 1 asks how repeated passage through tubeworm hosts changes the strain
composition and diversity of an environmental symbiont reservoir when all
evolution is neutral. In particular, the experiments must distinguish three
possible mechanisms:

1. small, unrepresentative infection samples are amplified inside hosts and
   cause existing environmental strains to be lost or become dominant;
2. mutation inside hosts creates new strains that enter the environment; and
3. host passage changes the relative abundances of existing strains without
   greatly changing how many strains are present.

The design also tests whether the environmental population approaches a stable
statistical diversity distribution, which parameters determine that distribution,
and whether many hosts returning a small number of bacteria each differ from a
few hosts returning many bacteria each.

The primary questions are:

1. Does repeated host passage change environmental strain composition relative
   to a population with no host-derived return?
2. Does passage mainly remove existing diversity, create new diversity by
   mutation, or redistribute existing diversity?
3. Which combinations of host abundance, escape, and mutation reduce or increase
   diversity, increase richness while reducing evenness, or have negligible
   effects?
4. Does environmental diversity approach a stable statistical level, and what
   controls its level, convergence time, and continuing fluctuations?
5. Which parameters and interactions explain the most variation in diversity?
6. Do host abundance and escape have separate effects, or do they operate mainly
   through the total number of bacteria returned?
7. At a fixed total return, does pooling across many hosts differ from receiving
   more bacteria from fewer hosts?

## Phase 1 assumptions

All main Phase 1 experiments use the following assumptions:

- The environmental population is a non-depleting effective reservoir.
- Infection does not remove cells from that reservoir.
- Hosts are infected independently from the same reservoir.
- Escapees from all hosts are pooled before being returned.
- Environmental capacity is fixed.
- Capacity is regulated by neutral Hamilton proportional apportionment.
- There is no free-living bacterial generation.
- There is no free-living mutation.
- Free-living selection is disabled.
- Within-host selection is disabled for the neutral Phase 1 experiments.
- All initial strains have equal reproductive fitness.
- Selection-enabled experiments are outside the main Phase 1 analysis and must
  be clearly labelled as separate sensitivity analyses or deferred to Phase 2.

## Terms and parameters

| Symbol or term | Meaning |
|---|---|
| `H` | Number of hosts in one simulated host generation. This is `population_size` in the configuration. |
| `B` | Infection bottleneck: the number of bacterial cells that initially infect one host. |
| `K` | Within-host carrying capacity: the maximum bacterial population represented inside one adult host. |
| `g` | Within-host growth factor used to expand the infection from `B` cells toward `K`. |
| `T_s` | Number of steady bacterial generations spent at carrying capacity after growth. |
| `c` | Environmental capacity ratio. It sets the effective reservoir size relative to one host's carrying capacity. |
| `N_E` | Fixed effective environmental capacity, calculated as `N_E = round(cK)`. |
| `f` | Escape fraction: the fraction of the adult within-host population released by each host. |
| `e` | Whole number of cells released by one host, calculated as `e = round(Kf)`. |
| `R` | Total cells returned by all hosts in one host generation, `R = He`. |
| \(\alpha\) | Realized host-feedback fraction before environmental capacity regulation, \(\alpha=R/(N_E+R)\). |
| `u` | Probability that one bacterial offspring becomes a new strain during one within-host bacterial generation. This is `mutation_probability`. |
| Strain | A unique labelled type in the model. Under the infinite-alleles assumption, every mutation event creates a new strain label. |
| Richness | Number of distinct strain labels present. |
| Evenness | How evenly bacterial abundance is distributed among the strains that are present. |
| Host passage | Infection, growth inside hosts, escape, pooling, and return to the environment, repeated over host generations. |
| Experimental group | A version of the simulation in which selected mechanisms are enabled or disabled. This is sometimes called an experimental “arm.” |
| Mechanistic control | A group in which one biological process is removed so that its contribution can be identified by comparison. |

The environmental population is an **effective reservoir**. Its capacity is the
number of environmental cells represented by the model for infection and host
feedback. It is not necessarily intended to equal the literal census of every
free-living symbiont in the habitat.

## What happens during one host generation?

### 1. Hosts acquire bacteria from the environment

Each of `H` hosts receives `B` bacterial cells sampled from the environmental
reservoir. Sampling is with replacement at the reservoir level, so infection of
one host does not deplete the environment or change what is available to another
host.

The infection sample may not represent the environmental population well. Rare
strains may be absent, and common strains can be over- or under-represented by
chance. This sampling effect is especially important when `B` is small.

### 2. Bacteria multiply inside each host

The `B` infection founders reproduce until the within-host population reaches
`K`. Population size grows according to the configured growth factor and is
then held at `K` for `T_s` additional bacterial generations.

Reproduction follows an exact neutral haploid Wright--Fisher process. Every
bacterial generation is a new random sample of offspring from the preceding
generation. This creates neutral genetic drift even when all strains have equal
fitness.

In mutation-enabled groups, each offspring has probability `u` of becoming a
new, uniquely labelled strain. Mutation happens only inside hosts in Phase 1.

### 3. Each host releases bacteria

At the adult endpoint, the model samples `e` cells without replacement from
each host, where

\[
e=\operatorname{round}(Kf).
\]

The exact integer `e`, rather than only the decimal escape fraction, must be
recorded in the experimental manifest. Two nominal escape fractions that round
to the same `e` describe the same treatment.

### 4. Escapees from all hosts are pooled

The pooled release contains

\[
R=He
\]

cells. Pooling may average over the chance compositions of many independent
infections. Therefore, a release assembled from many hosts may differ from an
equally large release assembled from only a few hosts.

### 5. The pooled release is added to the unchanged reservoir

Immediately before regulation, the mixed population contains `N_E + R` cells.
The fraction of that mixture that came from hosts is

\[
\alpha=\frac{R}{N_E+R}.
\]

The value of \(\alpha\) describes the immediate strength of host feedback. It is
often more interpretable than the escape fraction alone because it includes host
abundance, release per host, and environmental capacity.

### 6. The environment is returned to capacity

Hamilton proportional apportionment reduces the mixed population back to
`N_E` cells while preserving strain proportions as closely as whole-cell
counts permit. Strains with the largest fractional remainders receive any
leftover cells. Exact ties are resolved with a dedicated seeded random stream so
that strain identifiers do not receive a rounding advantage.

This step is capacity regulation, not a free-living bacterial generation. It
does not add environmental reproduction, mutation, selection, or conventional
Wright--Fisher drift.

Ignoring small whole-cell rounding effects, environmental composition follows

\[
p_{t+1}=(1-\alpha)p_t+\alpha q_t,
\]

where `p_t` is environmental strain composition and `q_t` is pooled release
composition. Host feedback is strong when \(\alpha\) is large and weak when it is
small.

### 7. The cycle is repeated

The regulated environment infects the next host generation. Repeated passage
allows small changes introduced by host return to accumulate.

### Numerical example

Suppose:

- `K = 10^9` bacteria per adult host;
- `c = 1`, giving `N_E = 10^9`;
- `H = 1,000` hosts;
- `f = 10^-4`, giving `e = 100,000` cells per host.

The total return is

\[
R=1{,}000\times100{,}000=10^8,
\]

and the realized feedback is

\[
\alpha=\frac{10^8}{10^9+10^8}=0.0909.
\]

Thus approximately 9.1% of the reservoir-plus-release mixture is host-derived
before it is reduced back to `10^9` cells.

## Experimental groups and mechanistic controls

### Group N: no host-derived return

Hosts are infected and within-host populations are simulated, but escape is set
to zero:

\[
f=0,\qquad e=0,\qquad R=0.
\]

Nothing from the hosts returns to the environmental reservoir. Under the Phase 1
reservoir assumptions, the environmental strain composition remains exactly
unchanged after its initial capacity scaling.

This is the reference group for the question: does repeated host passage change
the environmental population relative to a reservoir receiving no host-derived
bacteria?

A small no-return group with positive mutation can be retained as a calibration
control. Mutants may be observed inside hosts, but they cannot affect the
environment because none escape.

### Group A: host passage without mutation

Hosts return bacteria, but mutation is disabled:

\[
f>0,\qquad u=0.
\]

This group can change environmental composition through:

1. random sampling of the infection founders;
2. amplification of those founders to (K);
3. neutral within-host drift;
4. sampling of escapees from adults;
5. pooling across hosts; and
6. return and Hamilton regulation.

It cannot create new strains. Consequently:

- any richness loss must involve disappearance of strains that already existed;
- any change in abundance with little richness change is redistribution of
  existing diversity; and
- any increased dominance reflects amplification and neutral drift, not
  selection.

Comparing Group A with Group N measures the combined effect of infection
sampling, within-host amplification and drift, escape, pooling, and return.

### Group M: host passage with mutation

Hosts return bacteria and within-host mutation is enabled:

\[
f>0,\qquad u>0.
\]

This is the complete neutral Phase 1 passage group. New strains can arise inside
hosts, but many will be lost before reaching the environment. The retained
lineage and count outputs allow the following stages to be distinguished:

- mutations created during within-host reproduction;
- new strains surviving to the adult endpoint;
- new strains included in the escape sample;
- new strains retained in the regulated environment; and
- new strains persisting for multiple host generations.

Comparing Group M with Group A measures the incremental contribution of mutation
conditional on host passage.

### Bottleneck-size sensitivity

The primary experiment holds `B = 10`. Selected scenarios should be repeated at,
for example,

\[
B\in\{2,10,50\}.
\]

If amplification of small, unrepresentative infections is important, distortion
and diversity loss should generally be strongest at `B = 2` and weaker at
`B = 50`. This is an exploratory mechanistic sensitivity analysis, not a primary
factor in the confirmatory experiment.

### Initial-population sensitivity

Selected scenarios should begin from contrasting strain distributions. A stable
statistical diversity claim is stronger if different initial populations
eventually approach equivalent diversity distributions. If they do not, the
result should be described as initial-condition dependent or quasi-stationary.

## How the experimental groups distinguish mechanisms

The principal comparisons are:

- **Group A minus Group N:** diversity loss or redistribution produced by
  infection sampling, amplification, drift, escape, pooling, and return when no
  new strains can arise.
- **Group M minus Group A:** the additional contribution of within-host mutation.
- **Group M minus Group N:** the total effect of complete neutral host passage.
- **Different host numbers at the same (R):** the effect of pooling over many
  independent hosts, separated from total return.
- **Different bottleneck sizes:** evidence that small infection samples are the
  source of the amplification effect.

Three mechanisms can operate together:

1. **Loss:** existing strains disappear faster than new strains establish.
2. **Creation:** mutant strains establish faster than existing strains disappear.
3. **Redistribution:** strain abundances change substantially while richness
   changes little.

The analysis therefore measures all three rather than forcing each simulation
into a single explanation.

## Fixed settings and calibration choices

The recommended initial fixed settings, based on the biological-scale example,
are:

| Parameter | Proposed Phase 1 value | Reason |
|---|---:|---|
| Within-host carrying capacity, `K` | `10^9` | Biological-scale example and exact-count representation can store census size as counts. |
| Infection bottleneck, `B` | 10 | Agreed baseline; sensitivity considered separately. |
| Growth factor, `g` | 1.2 | Existing biological-scale schedule. |
| Steady bacterial generations, `T_s` | 100 | Existing biological-scale adult duration. |
| Capacity ratio, `c` | 1 | Gives `N_E = K`; sensitivity can examine other ratios. |
| Reservoir sampling | Non-depleting | Required Phase 1 assumption. |
| Within-host selection | Disabled | Required for a neutral baseline. |
| Free-living selection | Disabled | Required Phase 1 assumption. |
| Initial fitness values | All 1 | Prevents fitness differences among initial strains. |
| Mutation-effect mean and SD | 0 and 0 | Fitness effects are irrelevant when selection is disabled and zero values keep neutral metadata clear. |

With these values, a host experiences 201 within-host bacterial transitions: 101
during growth and 100 at carrying capacity. The sum of requested offspring
population sizes across one host lifetime is approximately
`1.0678 × 10^11`. Expected mutation events generated per host lifetime are
therefore approximately:

| Mutation probability (u) | Expected events per host lifetime |
|---:|---:|
| 0 | 0 |
| `10^-12` | 0.107 |
| `10^-11` | 1.07 |
| `10^-10` | 10.7 |
| `10^-9` | 106.8 |

These values are mutation-generation workloads, not expected surviving or
returned mutants. The pilot must measure survival through the adult, escape, and
environmental stages.

## Initial environmental strain population

The four-strain distribution in the biological-scale example is useful for
software demonstration but is too restricted for the main diversity experiment.
It leaves little room to distinguish richness loss from changes in evenness.

The recommended confirmatory starting population is:

- 100 labelled strains;
- one frozen Fisher log-series abundance realization using the historical
  parameter `a = 0.995`;
- integer counts scaled once to sum exactly to `N_E`;
- the identical count vector in every experimental group and replicate; and
- within-host and free-living fitness equal to 1 for every initial strain.

Starting counts should already sum to `N_E`. This prevents random tie-breaking
during initial capacity scaling from giving different replicates slightly
different starting environments.

Selected equilibrium-sensitivity runs should use:

- 100 equally abundant strains;
- an uneven 20-strain population; and
- an uneven 500-strain population.

The exact abundance vectors and their generation seed must be frozen in the
experiment manifest before confirmatory runs.

## Primary factors and proposed levels

### Host abundance

Begin with:

\[
H\in\{100,1000,10000,100000\}.
\]

These levels span the resolved scientifically relevant domain from 100 to
100,000 hosts on a logarithmic scale.

### Escape

Use the integer number released per host, `e`, as the exact design variable.
Store the corresponding `f = e/K` in each configuration.

For general mapping, begin with:

\[
f\in\{10^{-5},10^{-4},10^{-3},10^{-2}\}.
\]

Across `H = 100` to `H = 100,000`, these levels cover very weak to strong host
feedback without assuming that a linear change in escape has a linear biological
effect.

### Mutation probability

Treat `u = 0` as its own mechanistic group. Begin positive mutation calibration
with:

\[
u\in\{10^{-12},10^{-11},10^{-10},10^{-9}\}.
\]

Final levels should be selected using the number of novel strains that actually
reach the pooled release and regulated environment. Ideally, retained levels
should bracket regimes with approximately 0.1, 1, 10, and 100 newly returned or
newly established strains per environmental update.

The biological interpretation of (u) is a whole-genome mutation probability per
bacterial cell generation, not an unconverted per-site rate.

### Derived variables used for design and interpretation

The main derived variables are:

- total return `R = He`;
- feedback strength \(\alpha=R/(N_E+R)\); and
- expected per-host lifetime mutation supply
  `U = u × sum(N_j)`, where `N_j` is the within-host offspring population size in
  bacterial generation (j).

Exploratory mapping should be designed in `log H`, `log alpha`, and
`log U`, with `u = 0` kept as a separate group. This is more interpretable than
an unstructured factorial because (H), (e), and (R) are mathematically
related.

## Matched-return comparisons

The most direct test of host number versus total return holds (R) constant
while changing how many hosts contribute.

With `K = N_E = 10^9`, use the following exact comparisons:

| Total return `R` | Feedback \(\alpha\) | Host number `H` | Cells per host `e` | Escape fraction `f` |
|---:|---:|---:|---:|---:|
| `10^9` | 0.5 | 100 | `10^7` | `10^-2` |
| `10^9` | 0.5 | 1,000 | `10^6` | `10^-3` |
| `10^9` | 0.5 | 10,000 | `10^5` | `10^-4` |
| `10^9` | 0.5 | 100,000 | `10^4` | `10^-5` |

Across all four rows, total return and immediate environmental feedback are
identical. Only the number of independent hosts and the release depth per host
differ.

If the three cases have equivalent environmental outcomes, host abundance and
escape operate mainly through (R). If they differ, pooling across independent
hosts matters beyond total return.

Changing (H) at fixed (R) also changes the total number of infection founders,
(HB). This is part of the biological pooling effect: many-host treatments
sample the reservoir through more independent infections.

## Concrete first pilot matrix

This is a computational, mutation-supply, and variance pilot. It is not a
production or confirmatory experiment.

Use:

- `K = 10^9`;
- `B = 10`;
- growth factor 1.2;
- 500 steady bacterial generations;
- `c = 1`;
- the frozen 100-strain starting spectrum;
- three independently seeded replicates per cell; and
- five host generations per replicate.

| Pilot block | Host and return settings | Mutation settings | Cells |
|---|---|---|---:|
| No return | `H = 100, e = 0` | `u = 0, 10^-10` | 2 |
| Matched `R = 10^9` | All four `H,e` pairs in the table above | `u = 0` | 4 |
| Mutation bracket | `H = 100, e = 10^7, R = 10^9` | `u = 10^-12, 10^-11, 10^-10, 10^-9` | 4 |
| Feedback boundaries | Very weak `H=100, f=10^-5`; strong `H=1,000, f=10^-2` | `u = 0` | 2 |
| **Initial core** |  |  | **12** |

Five extension cells are pre-specified but are not launched automatically. Once
the mutation bracket selects an informative `u*`, three cells complete the
mutation-enabled `R=10^9` comparison at `H=1,000`, 10,000, and 100,000. Two
additional mutation-free cells compare `H=100, f=10^-3` with
`H=10,000, f=10^-5`; both have `R=10^8` and `alpha=0.090909`. This produces a
17-cell design only if the 12-cell core passes its resource and information
gates.

The pilot should measure:

- wall time per host generation;
- disk volume per output table and host generation;
- mutation events per bacterial transition and host lifetime;
- novel strains present in adults, releases, and the regulated environment;
- richness and host occupancy of the pooled adult and release populations;
- short-term changes in environmental richness, evenness, and composition;
- replicate-to-replicate variation; and
- the largest realized mutation count in any bacterial transition.

### Decision rule for expansion

Expand beyond the pilot only if all of the following hold:

1. Environmental capacity, infection size, adult size, release size, and lineage
   parentage invariants pass.
2. The no-return environment remains exactly unchanged.
3. Every output required by the agreed schema is present and internally
   consistent.
4. The largest realized per-transition mutation count remains below 25% of the
   configured materialization limit.
5. Projected runtime and storage for the next stage occupy no more than 70% of
   the available budget.
6. The positive mutation levels bracket an informative response. They must not
   all produce zero returned novelty or all produce overwhelming novelty.
7. Matched-return cells have exactly equal integer `R` and equal \(\alpha\).

If mutation levels fail criterion 6, move the mutation grid based on realized
novel release and establishment rates before increasing replicate numbers. Do
not silently rescale (K) as a computational shortcut because changing (K)
changes neutral drift and mutation supply.

## Staged execution plan

### Stage 0: scientific and software freeze

Before any production run:

1. Resolve the scientific choices listed at the end of this document.
2. Freeze the initial strain-count vector.
3. Freeze the factor ranges and biological relevance thresholds.
4. Commit or archive the exact source tree and analysis code.
5. Validate every configuration.
6. Run the maintained unit and distributional validation suite.
7. Create a machine-readable experiment manifest listing every cell, seed block,
   derived `e`, `R`, \(\alpha\), expected mutation workload, and intended
   output mode.

### Stage 1: first pilot

Run the 12-cell, three-replicate, five-generation core above. Use it only for
feasibility, mutation-supply calibration, early effect sizes, and variance
estimation. Add the five extension cells only after applying the declared
decision rule.

### Stage 2: equilibrium and precision pilot

Choose approximately six sentinel conditions representing:

- no return;
- mutation-free passage at baseline feedback;
- mutation-enabled passage at baseline feedback;
- few-host and many-host extremes at fixed (R); and
- weak and strong feedback.

Begin with 12 independent seed blocks. Use the generation rule and convergence
diagnostics below. This stage estimates the required run length and number of
replicates for confirmatory contrasts.

### Stage 3: exploratory parameter mapping

Do not run the complete `H × e × u` factorial. Instead:

1. Place 18--24 maximin space-filling cells across `log H`, `log alpha`, and
   `log U`, retaining a separate `u = 0` stratum.
2. Run 6--8 replicates per exploratory cell initially.
3. Fit an uncertainty-aware response-surface model.
4. Add 8--12 cells only where prediction error is large, response curvature is
   high, or biological classification is uncertain.
5. Reserve several cells and seed blocks for out-of-sample validation.

This sequential design locates important boundaries without paying for a large
full factorial dominated by redundant or biologically uninformative conditions.

### Stage 4: confirmatory experiment

After pilot results have frozen factor levels, margins, run lengths, and analysis
code, run new held-out seed blocks for the pre-registered comparisons:

1. host passage versus no return;
2. mutation-enabled versus mutation-free passage;
3. few versus many hosts at fixed (R);
4. whether the pooling contrast changes with mutation; and
5. equivalence of long-run diversity from contrasting initial conditions.

### Stage 5: exploratory sensitivities

Only after the main neutral results are complete, examine selected cells under:

- alternative infection bottlenecks;
- alternative initial richness and evenness;
- alternative environmental capacity ratios; and
- uncertain growth or adult-duration values.

Selection, finite environmental sampling, and free-living bacterial generations
must remain separately labelled and must not be pooled with the neutral Phase 1
results.

## Number of host generations and equilibrium

The same number of host generations is not equally informative at every feedback
level. When \(\alpha\) is small, each host generation changes only a small fraction
of environmental composition.

Define cumulative feedback exposure after `G` host generations as

\[
x_G=G[-\log(1-\alpha)]\approx G\alpha
\]

for small \(\alpha\). This gives a common time scale for comparing weak and strong
feedback.

For the first long-run screen, use

\[
G_0=\min\left(5000,\max\left(250,
\left\lceil\frac{20}{-\log(1-\alpha)}\right\rceil\right)\right).
\]

Examples are:

- approximately 250 host generations for \(\alpha=0.0909\);
- approximately 2,011 host generations for \(\alpha=0.0099\); and
- more than 5,000 generations for \(\alpha\lesssim0.004\), which triggers the
  initial cap.

The 5,000-generation cap is a computational decision, not evidence of
equilibrium. A cell that has not passed the diagnostics by the cap must be
reported as “not converged by 5,000 generations.” Scientifically important cells
can then be extended under a predeclared second-stage rule.

Do not stop individual replicates early when they appear stable. Every replicate
within an experimental cell should have the same fixed horizon so that apparent
stability does not bias the estimated equilibrium.

### Definition of statistical equilibrium

Equilibrium here means stability of statistical summaries such as richness,
effective diversity, evenness, and their fluctuations. Exact strain labels are
not expected to become stationary when mutation continually introduces new
labels.

Evaluate three consecutive non-overlapping windows, preferably five
feedback-exposure units per window. Declare equilibrium for a response only when:

1. the 90% confidence or credible interval for change across each window is
   entirely inside the predeclared equivalence margin;
2. the means of adjacent windows are equivalent;
3. rank-normalized split \(\hat R<1.05\) across replicate trajectories;
4. combined effective sample size is at least 400 after accounting for temporal
   autocorrelation;
5. the criteria remain satisfied at the next assessment; and
6. sentinel runs from contrasting initial populations converge to equivalent
   summary distributions.

Convergence time is the first generation after which all criteria remain
satisfied. Report it in both host generations and feedback-exposure units.

When `u = 0`, the process may ultimately reach fixation at one strain. A
temporary plateau above one strain should be called quasi-stationary unless the
long-run behavior is established. Record fixation time separately and treat runs
that do not fix by the maximum horizon as right-censored.

## Stochastic replicates and precision

### Independent unit

One complete seeded simulation replicate is the independent experimental unit.
Individual hosts within a replicate share the same environmental history, and
successive generations form a time series. Neither should be counted as an
independent replicate for environmental inference.

### Seed blocking

Use the same master seed and replicate indices across matched experimental
groups. This creates common stochastic blocks and can reduce variation in
planned contrasts. Each replicate index remains independent of the other
indices.

Mutant strain identifiers must not be treated as homologous across configurations.
Comparisons of mutant lineages across experimental groups should use ancestry
categories and summary distributions, not matching numeric mutant IDs.

### Precision pilot and final replicate number

Start the long-run precision pilot with 12 independent seed blocks at the
sentinel conditions. For a paired contrast with pilot standard deviation `s_d`
and desired confidence-interval half-width (h), estimate

\[
n=\left\lceil\left(\frac{1.96s_d}{h}\right)^2\right\rceil.
\]

Calculate this requirement for the co-primary common-strain diversity,
dominant-strain diversity, and composition responses. Use the largest result.

For confirmatory comparisons:

- use at least 20 replicates;
- add replicates in batches of eight;
- set a predeclared maximum of 100 replicates per cell; and
- if more than 100 are required, reconsider the biological relevance margin or
  experimental contrast rather than treating large simulation volume as a
  substitute for a meaningful effect.

Uncertainty calculations must resample whole simulation replicates. Hosts must
not be resampled as though they were independent environmental experiments.

## Environmental response variables

### Primary diversity measures

Calculate the following from labelled environmental counts:

1. **Richness, `D_0`:** number of strains present.
2. **Effective common-strain diversity, `D_1`:**
   \[
   D_1=\exp\left(-\sum_i p_i\log p_i\right).
   \]
   This is the effective number of commonly represented strains.
3. **Effective dominant-strain diversity, `D_2`:**
   \[
   D_2=\frac{1}{\sum_i p_i^2}.
   \]
   This emphasizes the most abundant strains.
4. **Pielou evenness:** Shannon diversity divided by the maximum possible value
   for the observed richness.

Using `D_0`, `D_1`, and `D_2` together distinguishes, for example, the addition of many
rare mutants from an increase in the effective number of common strains.

The existing unbiased gene-diversity summary remains a secondary response.
Calculate `D_2` directly from raw environmental frequencies so that all Hill
diversities use the same convention.

The current metrics function returns evenness 1 when only one strain remains.
That is a convenient software edge-case value but is biologically misleading.
For analysis, treat evenness as undefined when richness is one and interpret it
only jointly with richness.

### Composition and turnover

Measure compositional departure from the no-return population using total
variation distance:

\[
TV(p_t,p_0)=\frac12\sum_i|p_{i,t}-p_{i,0}|.
\]

This has a direct interpretation as the fraction of abundance that would need to
be reassigned among strains to make the two compositions identical.

Also calculate:

- Jensen--Shannon divergence as a bounded secondary composition measure;
- turnover between consecutive host generations;
- abundance of rare strains and singleton strains;
- frequency of the most abundant strain; and
- fraction of environmental abundance derived from within-host mutations.

### Continuing fluctuations at equilibrium

Within the diagnosed equilibrium window, calculate for each replicate:

- mean and median diversity;
- standard deviation and coefficient of variation;
- 5th and 95th percentiles;
- lag-one autocorrelation;
- integrated autocorrelation time;
- median generation-to-generation composition turnover; and
- frequency and duration of excursions outside the equilibrium interval.

These responses distinguish a stable mean with continuing large fluctuations
from a tightly regulated stable population.

## Mechanistic decomposition

### Exact richness accounting

For every environmental update define:

- `L_t`: strains present at generation `t` but absent at `t + 1`;
- `G_t`: strains absent at `t` but newly present at `t + 1`.

Verify:

\[
\Delta D_{0,t}=G_t-L_t.
\]

Separate `G_t` into strains newly created by mutation and previously existing
strains reappearing after being absent from a recorded snapshot, if such
reappearance is possible in the retained state.

### Tracking mutant fate

Use lineage, adult, release, and environmental tables to count mutants at each
filter:

1. generated within a host;
2. surviving to the adult population;
3. sampled among escapees;
4. retained after environmental capacity regulation; and
5. persisting for at least 1, 5, 10, and additional predeclared host generations.

Model each transition as a survival or transmission probability. This reveals
whether mutation creates substantial diversity inside hosts but little enduring
environmental diversity.

### Root-collapsed diversity

Trace each mutant through the lineage-event table to its ultimate ancestor among
the initial environmental strains. Recalculate environmental diversity after
assigning all mutant descendants to that initial root.

Compare:

- **labelled diversity**, where every mutant is a separate strain; and
- **root-collapsed diversity**, where mutant descendants are combined with their
  original ancestry.

The difference measures diversity created by mutation. Changes in root-collapsed
composition measure redistribution among the original environmental ancestries.

### Infection and pooling statistics

Use the infection-founder and pooled-host outputs to calculate:

- number and fraction of environmental strains represented among all infections;
- expected founder coverage
  `1 - (1 - p_i)^(HB)` for strain `i`, compared with observed coverage;
- variation in strain frequency among hosts;
- number of hosts occupied by each strain;
- adult and release-pool divergence from the pre-infection environment;
- release occupancy across hosts, derived from per-host release counts; and
- effective number of hosts contributing to each dominant strain.

These measurements explain why matched total returns can differ between
few-host and many-host conditions.

## Biological outcome classification

Scientific relevance margins must be approved before confirmatory runs.
Provisional margins are:

- 5% for richness, `D_1`, and `D_2`;
- 0.02 absolute units for evenness; and
- 0.05 for total-variation distance.

Use simultaneous intervals and classify each parameter combination as follows:

- **Diversity reduced:** both `D_1` and `D_2` are at least 5% below the
  no-return reference.
- **Diversity increased:** both `D_1` and `D_2` are at least 5% above the
  reference.
- **Richness increased but evenness reduced:** richness is at least 5% higher and
  evenness is at least 0.02 lower. Report the directions of `D_1` and `D_2`
  alongside this mixed outcome.
- **Negligible effect:** 90% intervals for richness, `D_1`, and `D_2` lie
  entirely within ±5%, and the interval for evenness lies within ±0.02.
- **Mixed or uncertain:** all other combinations.

Call composition meaningfully changed only if the lower simultaneous 95%
interval for total-variation distance exceeds the approved relevance threshold.
Because the no-return environmental state is deterministic, the objective is to
identify biologically meaningful change rather than merely detect any nonzero
numerical difference.

## Statistical analysis

### Confirmatory models

For each replicate, summarize the diagnosed equilibrium window before fitting
between-treatment models. This avoids treating autocorrelated generations as
independent data.

Use a blocked multivariate model containing:

- experimental group;
- total return `R` or feedback strength \(\alpha\);
- host number at fixed return, written `H | R`;
- mutation level;
- `R × u`, `H | R × u`, and other predeclared interactions; and
- seed block.

Recommended response models are:

- Gaussian models on log-transformed `D_0`, `D_1`, and `D_2`;
- beta regression for evenness when richness is greater than one;
- negative-binomial models for numbers of created, returned, and persistent
  mutant strains;
- binomial or beta-binomial models for lineage transmission proportions;
- replicate-level bootstrap intervals for composition and turnover; and
- survival or accelerated-failure-time models for convergence and fixation,
  including right-censored runs.

Use planned contrasts rather than all possible pairwise comparisons. Apply
simultaneous confidence intervals or Holm adjustment across the small family of
confirmatory tests. Equivalence conclusions should use two one-sided tests or an
equivalent interval-based procedure.

### Testing total return versus pooling

Parameterize the matched-return analysis using `R` and `H | R`, rather than
trying to fit `H`, `e`, and `R = He` simultaneously.

- The (R) term estimates the effect of changing the total number returned.
- The `H | R` term estimates the effect of redistributing the same return
  across different numbers of hosts.
- `H | R × u` tests whether pooling changes the environmental
  contribution of mutation.

If the interval for `H | R` lies entirely within the biological equivalence
margin, an (R)-only explanation is supported. Failure to reject a zero effect
is not sufficient; an explicit equivalence test is required.

### Exploratory response-surface model

Fit a heteroskedastic Gaussian-process emulator or hierarchical generalized
additive model over `log H`, `log alpha`, and `log U`, with a separate
indicator for zero mutation.

Assess prediction using held-out parameter cells and seed blocks. Add design
points sequentially where:

- cross-validated prediction error exceeds 10%;
- the estimated response is strongly nonlinear;
- uncertainty spans two biological classifications; or
- factor interactions are poorly identified.

### Parameter importance

Use the final validated response surface to estimate first-order and total Sobol
indices over a predeclared parameter domain and distribution, preferably
log-uniform for positive factors spanning orders of magnitude.

Report:

- the variation explained by each main factor;
- total influence including interactions;
- pairwise and important three-way interaction contributions;
- bootstrap uncertainty from simulation replicates; and
- uncertainty caused by fitting the response surface.

This is preferable to order-dependent sums of squares from a large factorial
ANOVA.

## Mapping the analyses to the seven questions

### 1. Does host passage change environmental composition?

Compare Groups A and M with Group N using total-variation distance, Hill
diversities, richness, evenness, and temporal turnover. The no-return reservoir
should remain exactly unchanged.

### 2. Does passage remove, create, or redistribute diversity?

Use Group A versus N for removal and redistribution, Group M versus A for
mutation, exact gained/lost richness accounting, mutant transmission filters,
and labelled versus root-collapsed diversity.

### 3. Under which parameter combinations does each outcome occur?

Apply the predeclared outcome classification to predictions from the exploratory
response surface, then verify important boundary cells with held-out simulations.

### 4. Does diversity approach a stable statistical level?

Apply the multi-window equilibrium diagnostic, initial-condition sensitivity,
and right-censored convergence analysis. Model equilibrium means, convergence
times, autocorrelation, and fluctuation magnitudes separately.

### 5. Which parameters and interactions explain the most variation?

Use validated response-surface predictions and Sobol first-order and total
indices with uncertainty intervals.

### 6. Are host number and escape independent, or mainly a total-return effect?

Use the general response surface together with the `R, H | R` matched-return
model and equivalence tests.

### 7. Does pooling across hosts matter at constant total return?

Use the exact `R = 10^9` four-host-level block and, if the expansion rule is
met, the `R = 10^8` endpoint pair. Interpret environmental differences through
founder coverage, among-host composition variance, occupancy, and release-pool
representativeness.

## Output and retention policy

Retain the agreed output specification.

| Output | Required content and use |
|---|---|
| `environment_counts.csv` | Labelled environmental counts, frequencies, and fitness metadata. Phase 1 pilot configurations use `environment_counts_mode = "final"`, so only the endpoint of each replicate is retained. |
| `infection_counts.csv` | Infection founders by host and generation. Source for founder coverage and representativeness. |
| `host_adult_summaries.csv` | Richness, diversity, and mean fitness for every adult host. |
| `host_adult_counts.csv` | Full adult counts for sentinel/mechanistic runs or a reproducibly selected host panel for large mapping runs. |
| `pooled_host_counts_and_occupancy.csv` | Total adult abundance and number of occupied hosts for each strain. |
| `release_counts.csv` | Per-host released counts. Source for pooled-release composition and release occupancy. |
| `strain_lineage_events.csv` | Mutation parentage, within-host generation, depth, and fitness metadata. Source for mutation fate and root-collapsed analysis. |
| `host_generation_summary.csv` | Environmental richness and gene diversity, adult summaries, return size, and realized feedback. |
| `resolved_config.json` | Exact resolved inputs. |
| `provenance.json` | Seed, configuration hash, source revision, software versions, model specification, output schema, Python/NumPy versions, and platform. |
| Checkpoints and final `.npz` files | Environmental states for inspection, recovery, and final verification. |

Recommended adult-count modes are:

- full adult counts for first-pilot and mechanistic runs with `H <= 100`;
- a deterministic panel of 100 hosts when `H = 1,000` during mapping;
- full counts for selected confirmatory cells when projected storage permits; and
- adult summaries for every host in every run.

Raw results should be written outside Git. After quality control, make each run
directory immutable, generate SHA-256 checksums, and then compress or transfer
it to the research data store. The other required scientific tables remain
available, but successful Phase 1 pilot runs retain only the final labelled
environmental state rather than a complete environmental composition history.
The verified final state is stored both as `environment_counts.csv` and as the
checksummed `final_environment_repNNN.npz` recovery-independent artifact.

The runner supports `environment_counts_mode = "all"` for a specifically
approved time-series analysis and `environment_counts_mode = "final"` for the
agreed pilot retention policy. Checkpoints still contain the most recent
environmental state while a run is active, regardless of this output setting.

## Reproducibility safeguards

1. Use one versioned TOML configuration per experimental cell.
2. Never change parameters inside a notebook or analysis script.
3. Maintain a machine-readable design manifest containing configuration hashes,
   seed blocks, expected derived parameters, and output locations.
4. Freeze all starting count vectors and record their checksums.
5. Use the maintained logical seed hierarchy so worker number and batch size do
   not change results.
6. Re-run a selected cell at two worker settings as a production quality-control
   check.
7. Record exact source code, model version, output-schema version, Python, NumPy,
   and platform information.
8. Run configuration validation and automated tests before the first production
   batch.
9. Check every result for environmental capacity, release size, non-negative
   counts, unique strain IDs, valid lineage parents, and consistent pooled
   abundance.
10. Reserve confirmatory seeds before exploratory analysis and do not replace
    inconvenient outcomes.

The runner now records both the Git commit and a content hash of the maintained
source, so dirty and untracked model code cannot go unnoticed. The exact source
must still be committed or archived with a source manifest before production;
the hash detects a difference but does not by itself reconstruct the source.

Validated generation-boundary restart is implemented. Active runs retain the
two newest recovery checkpoints, and `trophosome run --resume` verifies source
and configuration hashes, truncates incomplete post-checkpoint output, and
continues from the newest valid state. Recovery checkpoints are removed after a
successful verified run, leaving the final environmental state and scientific
output tables.

## Resolved Phase 1 choices

The first-pilot manifest freezes whole-genome mutation probabilities, a fully
neutral fitness profile, the 100-strain Fisher log-series vector, host abundance
from 100 to 100,000, escape from `10^-5` to `10^-2`, `B=10`, `K=10^9`, growth
factor 1.2, 500 steady bacterial generations, and primary `c=1`. Sensitivities
will consider `c=0.5` and `c=H`. Biological relevance margins are 5% for Hill
diversities, 0.02 for evenness, and 0.05 for composition distance. The computing
limits are 48 hours and 100 GB, with expansion permitted only below 70% of each
budget. Whole-cell Hamilton apportionment at very weak return is an accepted
part of the effective-reservoir interpretation. Checkpoint restart is enabled at
a user-defined interval, defaulting to one hour.

## Confirmatory and exploratory separation

### Confirmatory

- baseline host passage versus no return;
- mutation-enabled versus mutation-free passage at fixed host and return
  settings;
- few-host versus many-host comparisons at exact fixed total return;
- predeclared equivalence and outcome-classification rules;
- selected convergence and initial-condition equivalence tests; and
- new held-out seed blocks generated after pilot decisions are frozen.

### Exploratory

- global response mapping across host number, feedback, and mutation;
- adaptive placement of additional parameter cells;
- Sobol parameter and interaction importance;
- bottleneck, capacity-ratio, initial-spectrum, growth, and adult-duration
  sensitivities;
- detailed strain-frequency, occupancy, fluctuation, and lineage-persistence
  patterns; and
- any selection-enabled sensitivity, which must remain clearly separated from
  the neutral Phase 1 evidence.

## Recommended immediate next step

Commit and freeze the prepared 12-cell core, then run the maintained validation
suite and begin the Mac timing gate. The prepared manifest remains labelled
`prepared-not-launched` until the source snapshot is committed and the first
cell is explicitly started.
