# Phase 1: Stage 3, Wave 2

Frozen on 2026-09-03. This is an exploratory equivalence-and-memory experiment,
not a confirmatory test. Preparation creates configurations, manifests and the
adaptive decision machinery; it does not launch simulations.

## Questions and experiment size

Wave 2 contains two deliberately separate experiments.

1. **Host number by infection bottleneck (`H-by-B`):** can many severe
   bottlenecks collectively provide a representative host return?
2. **Host feedback by regional exchange (`alpha-by-m`):** does immigration
   erase host-induced change, or stabilize the focal environmental population?

Both experiments use 12 matched seed blocks (`sb0001`-`sb0012`) and a common
primary endpoint at passage 100. The 40-condition analysis contains 480
populations: 408 new populations from 34 conditions and 72 passage-100
environmental trajectories reused from earlier work. Five are exact prior
conditions. The alpha=0, m=0.1 environmental control reuses c0021 because, with
no return and non-depleting infection, host number cannot affect the focal
environment; its host-level outputs are not reused. Reuse does not create
additional replication.

All primary passage-100 comparisons remain valid regardless of which
trajectories later receive exploratory time extensions.

## Experiment A: host number by bottleneck size

Hold alpha=0.5, m=0.1, u=0 and total return R=10^9 cells. Escape fraction is
adjusted with H so every condition returns the same number of cells. The only
experimental factors are H and B.

| H | Escape fraction f | B=1 | B=5 | B=10 | B=50 |
|---:|---:|---|---|---|---|
| 100 | 0.01 | c0051 | c0052 | c0053* | c0054 |
| 1,000 | 0.001 | c0055 | c0056 | c0057 | c0058 |
| 10,000 | 0.0001 | c0059 | c0060 | c0061* | c0062 |

An asterisk marks an exact passage-100 reuse: c0053 reuses c0022 and c0061
reuses c0024.

Four exact total-founder comparisons test whether the product HB is sufficient:

| HB | Matched conditions |
|---:|---|
| 1,000 | c0053 (100 x 10) and c0055 (1,000 x 1) |
| 5,000 | c0054 (100 x 50) and c0056 (1,000 x 5) |
| 10,000 | c0057 (1,000 x 10) and c0059 (10,000 x 1) |
| 50,000 | c0058 (1,000 x 50) and c0060 (10,000 x 5) |

The primary analysis compares an HB-only model with a model containing
separate log(H), log(B) and interaction terms. Pairwise equivalence is assessed
from seed-paired TV log-ratios. The exploratory contour criterion is a TV ratio
from 0.8 to 1.25; the broader absolute biological margin of 0.05 is also
reported.

## Experiment B: host feedback by regional exchange

Hold H=1,000, B=10 and u=0. Cross four feedback targets with seven regional
exchange fractions:

The positive-feedback return totals use the same 10,000-cell rounding
convention as Wave 1: 10,100,000 cells for target alpha=0.01; 111,110,000 for
alpha=0.1; and 99,000,000,000 for alpha=0.99. This makes the three positive
m=0.1 reuses exact rather than approximately matched.

| Alpha target | m=0 | 0.001 | 0.01 | 0.1 | 0.5 | 0.9 | 0.99 |
|---:|---|---|---|---|---|---|---|
| 0 | c0063 | c0064 | c0065 | c0066* | c0067 | c0068 | c0069 |
| 0.01 | c0070 | c0071 | c0072 | c0073* | c0074 | c0075 | c0076 |
| 0.1 | c0077 | c0078 | c0079 | c0080* | c0081 | c0082 | c0083 |
| 0.99 | c0084 | c0085 | c0086 | c0087* | c0088 | c0089 | c0090 |

The starred m=0.1 conditions reuse c0021, c0037, c0039 and c0041,
respectively. The c0021 reuse is limited to its environmentally equivalent
no-return trajectory; the other three are exact configurations. Each migration
rate has its own alpha=0 control. This separates host-induced TV from background
change caused by regional exchange itself.

The immediate retained host signal is approximately alpha(1-m). Deliberate
near-matches include:

- about 0.001: (alpha=0.01, m=0.9) and (0.1, 0.99);
- about 0.01: (0.01, 0.001), (0.1, 0.9) and (0.99, 0.99);
- about 0.10: (0.1, 0.001) and (0.99, 0.9).

Equal outcomes along these contours would support immediate dilution as the
dominant mechanism. Differences would show that environmental memory matters
in addition to alpha(1-m).

For this panel, define **erasure** as lower mean host-induced TV relative to the
same-m alpha=0 control. Define **stabilization** separately as lower continuing
variation through time and/or lower variation among seed blocks. The registered
responses are TV at passage 100, mean and within-population SD of TV during
passages 51-100, among-seed variation, D1, D2 and evenness. Compare an
alpha(1-m)-only relationship against separate alpha, m and interaction terms
and a finite-time environmental-memory model. Every contrast preserves seed
pairing.

B=2 and B=20 are not part of this frozen batch. They require a separately
reviewed design amendment if the H-by-B surface shows important curvature.

## Frozen adaptive time horizon

One hundred passages are sufficient for the shared finite-time comparison but
cannot establish long-term stabilization when m is very small. For example,
immigration at m=0.001 has an approximate compositional half-life of 693
passages. Conversely, m=0 is a mechanistic no-immigration control and should
not be interpreted as guaranteeing a diversity equilibrium.

Every new Wave 2 TOML therefore has a maximum horizon of 1,000 passages. The
launcher initially pauses it exactly after passage 100 and retains a verified
checkpoint. A selected extension resumes the same environmental state, lineage
state and random-number trajectory; it does not rerun passages 1-100 or start a
new population.

Only alpha-by-m conditions with m in {0, 0.001, 0.01} are eligible. H-by-B and
m>=0.1 conditions stop at the primary passage-100 endpoint.

### Passage 100 to passage 500

Six central diagnostic anchors always continue:

- alpha=0 and 0.1 at m=0;
- alpha=0 and 0.1 at m=0.001; and
- alpha=0 and 0.1 at m=0.01.

The remaining eligible positive-feedback conditions continue only if raw TV or
host-induced TV remains unresolved. Host-induced TV is treatment TV minus its
same-m alpha=0 control within each seed block. A selected treatment always
brings its corresponding control.

The passage-100 assessment compares mean TV in passages 51-75 with passages
76-100. Across the 12 paired seed changes, a 90% Student-t interval must lie
fully inside:

```text
plus or minus max(0.002 TV, 25% of the absolute late-window mean)
```

Otherwise the trajectory is labelled unresolved and extended. This combines a
small absolute floor with the pre-specified 25% contour scale, avoiding an
unstable percentage criterion near TV=0.

### Passage 500 to passage 1,000

The assessment is repeated using passages 401-450 and 451-500. The alpha=0 and
0.1 anchors at m=0 and m=0.001 continue by design. Other conditions that reached
500 continue only when the same stability rule remains unresolved. The maximum
horizon is 1,000; the workflow never extends beyond it automatically.

Each decision is written once with checksums of the assessed TV trajectory
prefixes. Re-running the assessment verifies the same decision and refuses to
silently replace it if an input changed. The assessment never launches the
next stage: the researcher reviews the frozen decision and explicitly starts
the authorized continuation.

These adaptive long-run results are diagnostic and exploratory. Selection for
extension must be acknowledged when interpreting passage-500 or passage-1,000
subsets. It does not alter or filter the complete passage-100 analysis.

## Fixed biology and retained outputs

Scientific model 2.1.0, software 0.7.0 and output schema 2.3.0 are recorded.
The focal and fixed regional populations begin with the same frozen 100-strain
composition. K=N_E=10^9, capacity ratio c=1, growth factor=1.2 and 500 steady
bacterial generations are fixed. Infection samples a non-depleting reservoir;
all escapees are pooled. Mutation and both selection switches are off.

Complete labelled environmental trajectories are retained. Adult strain counts
are complete for H=100 and use the reproducible 100-host panel for H>=1,000.
All other agreed infection, adult-summary, pooled-occupancy, release, migration,
lineage, generation-summary and provenance outputs remain enabled.

Because continuation validates the model-source checksum stored in the
checkpoint, do not change `src/trophosome/` or `pyproject.toml` between the
passage-100 and later runs. Reporting scripts and derived results may be added
without changing the simulated trajectory.

## Frozen files and operation

The machine-readable design is
[`phase1-stage3-wave2-v210-adaptive-g1000-cells.tsv`](../experiments/work/trophosome/p01-neutral-feedback/design/phase1-stage3-wave2-v210-adaptive-g1000-cells.tsv).
The manifest records the complete decision rule and the checksums of 7,272
reused metric-trajectory rows. These frozen rows support Wave 2 comparisons but
do not replace the original labelled-count archives; retain the Stage 2 and Wave
1 raw outputs for strain-level reanalysis. The 408 new run configurations have
isolated scratch paths.

Use the [HPC instructions](../scripts/hpc/README.md#phase-1-stage-3-wave-2) to
prepare, inspect and run each boundary. `--prepare-only`, `--dry-run` and
`--assess-only` do not launch a simulation.
