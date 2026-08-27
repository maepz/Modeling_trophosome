# Phase 1 second pilot: stationarity and precision

**Report date:** 2026-08-27  
**Analysis audit:** `PASS`  
**Populations:** 72 (6 cells x 12 matched seed blocks)  
**Host passages:** 250

## What this pilot asks

This second pilot asks whether environmental diversity appears statistically stable late in a 250-passage run, how much it continues to fluctuate, and how many independent matched replicates a later confirmatory experiment may need.

## Main result

5 of 30 cell-response screens passed; 0 of 6 cells passed all five responses.

Passing this screen supports late-run stationarity, not definitive equilibrium. The experiment starts all runs from the same environmental composition, so it does not test convergence from different starting communities.

## Sentinel conditions

| Cell | Biological role | Hosts | Escape fraction | Total return | Mutation rate | Screen passes | Persistent screen |
|---|---|---:|---:|---:|---:|---:|---:|
| c0021 | no return | 100 | 0 | 0 | 0 | 1/5 | not by 250 |
| c0022 | baseline host passage | 100 | 0.01 | 1,000,000,000 | 0 | 1/5 | not by 250 |
| c0023 | within-host mutation | 100 | 0.01 | 1,000,000,000 | 1e-10 | 0/5 | not by 250 |
| c0024 | host pooling | 10,000 | 0.0001 | 1,000,000,000 | 0 | 1/5 | not by 250 |
| c0025 | weak feedback | 100 | 0.001 | 100,000,000 | 0 | 1/5 | not by 250 |
| c0026 | strong feedback | 100 | 0.1 | 10,000,000,000 | 0 | 1/5 | not by 250 |

![Environmental trajectories](figures/phase1-second-pilot-v210-m010-g250-report/environmental-trajectories.png)

**Figure 1.** Lines are medians across 12 independent populations; shaded bands span the 10th to 90th percentiles. TV is total-variation distance from the initial environmental composition.

![Individual seed-block trajectories](figures/phase1-second-pilot-v210-m010-g250-report/individual-seed-trajectories.png)

**Figure 2.** Thin coloured lines are all 12 independently seeded populations for the three conditions with the strongest seed-to-seed variation; the thick black line is their median. Persistent separation among seeds is expected in a stochastic model and does not, by itself, show continuing directional change.

Independent populations are not expected to converge to an identical diversity value. A stationary stochastic population can keep fluctuating, and different populations can remain apart. Longer runs are useful only for asking whether directional change within each population disappears and whether the distribution across populations stops changing—not for forcing the lines to meet. The paired late-window tests compare each seed block with itself, so a stable difference between seeds does not make those tests fail. R-hat is the component that is sensitive to persistent differences among seed histories; here it is best read as a secondary common-distribution diagnostic, not as a demand that independent biological populations become identical.

## Post-hoc separated stability diagnostic

After separating temporal trends, stability of the seed-block distribution, and the magnitude of between-seed spread, 22 of 30 responses were classified as statistically stable, 8 remained unresolved, and 0 showed clear biologically meaningful continuing change. This diagnostic is exploratory because it was added after inspection of the registered pilot results.

![Separated stability diagnostic](figures/phase1-second-pilot-v210-m010-g250-report/separated-stability-diagnostic.png)

**Figure 3.** The classification uses the agreed biological margins. First, slopes are calculated within each seed block. Second, changes in the average and spread of the seed-block distribution are tested separately. Failure to demonstrate equivalence is labelled unresolved rather than continuing change; clear change requires a 90% interval outside the biological margin. R-hat and ESS are retained as secondary diagnostics only. T? marks an unresolved trend; D? marks unresolved stability of the seed-block distribution.

![Stationarity screen](figures/phase1-second-pilot-v210-m010-g250-report/stationarity-screen.png)

**Figure 4.** The registered combined screen passes only when both overlapping late-window assessments satisfy the predeclared equivalence limits, rank-normalized split R-hat is below 1.05, and approximate combined ESS is at least 400. R-hat asks whether replicate histories are compatible with a common long-run distribution; it does not ask whether their current values are identical.

![Continuing fluctuations](figures/phase1-second-pilot-v210-m010-g250-report/continuing-fluctuations.png)

**Figure 5.** Median within-run coefficient of variation in the last diagnostic window. Larger values mean more continuing fluctuation around the late-run level.

## Precision for the confirmatory experiment

The largest recommendation is 20 matched replicates. Recommendations target a 95% confidence-interval half-width of 0.05, use a minimum of 20 replicates, increase in batches of 8, and are capped at 100.

![Precision recommendations](figures/phase1-second-pilot-v210-m010-g250-report/precision-recommendations.png)

**Figure 6.** Dashed lines show the minimum of 20 matched replicates.

## Quality control and computational resources

- Analysis audit: `PASS`.
- Summed elapsed time: 33.77 hours.
- Output volume: 7.95 GiB.
- Highest recorded process-tree memory: 106.5 MiB.
- Every run was checked for 251 environmental states, constant reservoir size, migration counts, committed completion metadata, configuration hashes, and final-state checksums.

## Interpretation limits

- Stationarity diagnostics are an initial screen, not proof of equilibrium.
- Convergence from contrasting initial strain compositions is not tested.
- Linear first-pilot cost projections are planning estimates, not guarantees.
- Confirmatory seed blocks must remain held out until final design freeze.

## Short glossary

- **D0:** strain richness; every detected strain counts equally.
- **D1:** effective number of common strains; sensitive to richness and evenness.
- **D2:** effective number of dominant strains; gives more weight to abundant strains.
- **Evenness:** how similarly abundant the strains are.
- **TV:** total-variation distance from the initial environmental composition; 0 means identical and 1 means no overlap.
- **ESS:** effective sample size after accounting for temporal autocorrelation.
- **R-hat:** whether split replicate histories appear compatible with a common long-run distribution. It does not require replicates to have identical values; values close to 1 indicate better distributional agreement.
