# Phase 1 first-pilot report: fixed regional-pool migration

**Status:** exploratory pilot report; not a confirmatory analysis  
**Report date:** 2026-08-26  
**Experiment:** `phase1-first-pilot-v210-m010`

## Purpose and scope

This pilot tested whether the simulation is computationally feasible and whether the selected parameter levels produce informative biological responses. One complete simulated population is the independent replicate. The pilot does not establish equilibrium or provide confirmatory tests.

All cells include fixed regional-pool exchange after host return and environmental regulation (m=0.10). The no-return cells therefore measure migration alone, rather than an unchanged environment.

## Pilot design

| Cell | Main purpose | # hosts H | Escape fraction f | Escapees/host e | Total return R | Host fraction alpha | Mutation rate u | Migration m |
|---|---|---|---|---|---|---|---|---|
| c0001 | No-return baseline | 100 | 0 | 0 | 0 | 0 | 0 | 0.1 |
| c0002 | Mutation occurs in hosts but nothing returns | 100 | 0 | 0 | 0 | 0 | 1.00e-10 | 0.1 |
| c0003 | Fixed-return, feedback comparison | 100 | 0.01 | 10,000,000 | 1.00e+09 | 0.5 | 0 | 0.1 |
| c0004 | Fixed-return, pooling comparison | 1,000 | 0.001 | 1,000,000 | 1.00e+09 | 0.5 | 0 | 0.1 |
| c0005 | Fixed-return, pooling comparison | 10,000 | 1.00e-04 | 100,000 | 1.00e+09 | 0.5 | 0 | 0.1 |
| c0006 | Fixed-return, pooling comparison | 100,000 | 1.00e-05 | 10,000 | 1.00e+09 | 0.5 | 0 | 0.1 |
| c0007 | Lowest positive mutation level | 100 | 0.01 | 10,000,000 | 1.00e+09 | 0.5 | 1.00e-12 | 0.1 |
| c0008 | Low mutation | 100 | 0.01 | 10,000,000 | 1.00e+09 | 0.5 | 1.00e-11 | 0.1 |
| c0009 | Intermediate mutation, pooling comparison | 100 | 0.01 | 10,000,000 | 1.00e+09 | 0.5 | 1.00e-10 | 0.1 |
| c0010 | High mutation | 100 | 0.01 | 10,000,000 | 1.00e+09 | 0.5 | 1.00e-09 | 0.1 |
| c0011 | feedback comparison | 100 | 1.00e-05 | 10,000 | 1.00e+06 | 9.99e-04 | 0 | 0.1 |
| c0012 | Strong feedback boundary | 1,000 | 0.01 | 10,000,000 | 1.00e+10 | 0.909 | 0 | 0.1 |
| c0013 | Intermediate mutation, pooling comparison | 1,000 | 0.001 | 1,000,000 | 1.00e+09 | 0.5 | 1.00e-10 | 0.1 |
| c0014 | Intermediate mutation, pooling comparison | 10,000 | 1.00e-04 | 100,000 | 1.00e+09 | 0.5 | 1.00e-10 | 0.1 |
| c0015 | Intermediate mutation, pooling comparison | 100,000 | 1.00e-05 | 10,000 | 1.00e+09 | 0.5 | 1.00e-10 | 0.1 |
| c0016 | Weaker fixed-return, pooling comparison; feedback comparison | 100 | 0.001 | 1,000,000 | 1.00e+08 | 0.091 | 0 | 0.1 |
| c0017 | Weaker fixed-return, pooling comparison | 10,000 | 1.00e-05 | 10,000 | 1.00e+08 | 0.091 | 0 | 0.1 |
| c0018 | feedback comparison | 100 | 0.1 | 100,000,000 | 1.00e+10 | 0.909 | 0 | 0.1 |
| c0019 | More hosts feedback comparison | 1,000 | 1.00e-06 | 1,000 | 1.00e+06 | 9.99e-04 | 0 | 0.1 |
| c0020 | Weaker fixed-return, pooling comparison; More hosts, feedback comparison | 1,000 | 1.00e-04 | 100,000 | 1.00e+08 | 0.091 | 0 | 0.1 |

**Colour key in the PDF:** grey marks baseline/reference purposes; purple shades distinguish the host-number H comparisons; blue shades mark the feedback-alpha comparisons; and green marks the mutation-rate u comparison. Coloured text identifies cells that also belong to another comparison series. The wording and colour assignments reproduce slide 28 of the project update.

## Biological results

- Every cell used fixed regional-pool exchange at m=0.10. This replaces 10% of the focal environment after host return, so the expected host-derived fraction carried into the next passage is multiplied by 0.90.
- No-return controls provide the migration-only baseline: their median composition distance was 7.026e-05 after five passages. Host-return effects should be interpreted relative to this background exchange.
- At fixed total return, increasing host number from 100 to 100,000 changed the median composition distance from 0.1299 to 0.0042. Total return alone therefore did not determine the short-term response.
- At the highest mutation level (1e-09), median labelled richness was 444, but the median mutant abundance fraction was only 4.21e-07. Mutation created many rare labels rather than abundant new types.
- Across the fixed-H feedback series, increasing alpha from about 0.001 to 0.909 produced these endpoint changes (H=100: median TV 3.088e-04 to 0.2286, median D1 34.2 to 28.3; H=1,000: median TV 1.287e-04 to 0.082, median D1 34.2 to 32.7).

### Summary of biological effects

| Comparison | Diversity / composition | Richness D0 | Evenness | Main interpretation |
|---|---|---|---|---|
| Migration-only baseline (c0001, c0002) | Median composition distance 0.0001 | Median D0 100 | Median 0.767 | With f=0, changes arise from stochastic exchange with the fixed regional pool; within-host events cannot return. |
| Host return, u=0 (c0003) | Redistribution (median TV 0.130); abundance-weighted change varies by seed block. | No change | Median decrease at H=100; seed-block dependent. | A small number of independent host samples can amplify an unrepresentative infection pool. |
| Increase H, fixed R, u=0 (c0003, c0004, c0005, c0006) | Redistribution falls (median TV 0.130 to 0.004). | No change in this mutation-free series | Returns toward the initial value in two of three seed blocks. | At R=1e9, H controls averaging across independent infections; u was fixed at zero. |
| Increase u (c0003, c0007, c0008, c0009, c0010) | Composition changes; median D1 and D2 increase, but direction is seed-block dependent. | Strong increase (100 to 444) | Decreases as rare mutant labels accumulate | The u dose-response was measured at H=100; the extension adds an exploratory H by u comparison. |
| Increase alpha, H=100 (c0011, c0016, c0003, c0018) | Redistribution increases overall (median TV 3.088e-04 to 0.2286); median D1 34.2 to 28.3. | No change (median D0 100 throughout). | Median 0.767 to 0.726. | H and u are fixed; f, R and alpha increase together. The alpha=0.909 endpoint uses f=0.1, outside the primary range. |
| Increase alpha, H=1,000 (c0019, c0020, c0004, c0012) | Redistribution increases overall (median TV 1.287e-04 to 0.082); median D1 34.2 to 32.7. | No change (median D0 100 throughout). | Median 0.767 to 0.758. | H and u are fixed; f, R and alpha increase together. The alpha=0.001 endpoint uses f=1e-6, outside the primary range. |
| Lower-return pooling (c0016, c0020, c0017) | Redistribution falls (median TV 0.025 to 0.003). | No change (median D0 100 throughout). | Negligible net change (median 0.766 to 0.767). | At R=1e+08 and u=0, pooling across more hosts still reduces redistribution; 3 H levels were tested. |
| Mutation-enabled pooling (c0009, c0013, c0014, c0015) | Redistribution falls (median TV 0.152 to 0.005). | Moderate, non-monotonic increase (median D0 134-149; 141 at highest H). | Small, non-monotonic change (median 0.707-0.714). | Alongside the mutation-free R=1e9 series, this explores H by u: redistribution falls in both, while rare mutant labels remain. |

The three seed blocks are matched across cells. Comparisons therefore use within-seed changes: each coloured line in the comparison figures joins the same seed block across parameter levels. Black lines show cell medians. With only three blocks, these are descriptive paired patterns, not hypothesis tests.

![Endpoint diversity](figures/phase1-pilot-v210-m010-report/endpoint-diversity.png)

### Endpoint diversity indices

Values are cell medians across independent population replicates.

| Cell | D0 | Shannon H' | Simpson | D1 | D2 | Evenness | TV |
|---|---|---|---|---|---|---|---|
| c0001 | 100 | 3.533 | 0.95497 | 34.227 | 22.206 | 0.7672 | 7.026e-05 |
| c0002 | 100 | 3.533 | 0.95497 | 34.227 | 22.206 | 0.7672 | 7.026e-05 |
| c0003 | 100 | 3.4648 | 0.95274 | 31.971 | 21.157 | 0.7524 | 0.1299 |
| c0004 | 100 | 3.507 | 0.95378 | 33.348 | 21.634 | 0.7615 | 0.0564 |
| c0005 | 100 | 3.5331 | 0.9552 | 34.231 | 22.319 | 0.7672 | 0.018 |
| c0006 | 100 | 3.5329 | 0.95499 | 34.223 | 22.218 | 0.7672 | 0.0042 |
| c0007 | 101 | 3.4795 | 0.95676 | 32.443 | 23.127 | 0.7523 | 0.1517 |
| c0008 | 103 | 3.4795 | 0.95676 | 32.444 | 23.126 | 0.7446 | 0.1516 |
| c0009 | 134 | 3.4854 | 0.95676 | 32.636 | 23.126 | 0.7116 | 0.1517 |
| c0010 | 444 | 3.4854 | 0.95676 | 32.636 | 23.125 | 0.5731 | 0.1516 |
| c0011 | 100 | 3.533 | 0.95497 | 34.227 | 22.208 | 0.7672 | 3.088e-04 |
| c0012 | 100 | 3.4885 | 0.9536 | 32.738 | 21.552 | 0.7575 | 0.082 |
| c0013 | 147 | 3.5205 | 0.95486 | 33.8 | 22.155 | 0.7069 | 0.046 |
| c0014 | 149 | 3.5255 | 0.95455 | 33.971 | 22.002 | 0.7071 | 0.0162 |
| c0015 | 141 | 3.5332 | 0.95497 | 34.232 | 22.208 | 0.7145 | 0.0045 |
| c0016 | 100 | 3.5276 | 0.9547 | 34.042 | 22.074 | 0.766 | 0.0249 |
| c0017 | 100 | 3.5334 | 0.95504 | 34.24 | 22.241 | 0.7673 | 0.0029 |
| c0018 | 100 | 3.3425 | 0.95242 | 28.289 | 21.016 | 0.7258 | 0.2286 |
| c0019 | 100 | 3.533 | 0.95497 | 34.227 | 22.205 | 0.7672 | 1.287e-04 |
| c0020 | 100 | 3.5293 | 0.95473 | 34.101 | 22.09 | 0.7664 | 0.0094 |

### Fixed regional-pool migration

![Fixed-pool migration](figures/phase1-pilot-v210-m010-report/fixed-pool-migration.png)

Panel A shows how exchange reduces the expected host-derived fraction after return and regulation. Panel B isolates the richness change caused by replacing emigrants with immigrants during passage 5. All cells use the same m, so this pilot controls for migration but does not estimate its dose-response.

### Host pooling within fixed-return series

![Matched return](figures/phase1-pilot-v210-m010-report/matched-return.png)

Colours identify matched seed blocks. Circles show u=0 and squares show u=1e-10. Filled markers with solid lines show alpha=0.5; hollow markers with dashed lines show alpha=0.09. The thicker black underlays show cell medians.

### Within-host mutation

![Mutation bracket](figures/phase1-pilot-v210-m010-report/mutation-bracket.png)

Coloured lines connect the same seed block across u; the black line is the cell median.

### Feedback strength at fixed host abundance

The same four alpha levels are evaluated separately at H=100 and H=1,000. Within each series H and u are fixed, while f, total return R and alpha increase together. Lines connect the same seed block. The extreme f values used for c0018 and c0019 are sensitivity levels outside the primary escape-rate range.

![Feedback bracket](figures/phase1-pilot-v210-m010-report/feedback-bracket.png)

### Extension pooling comparisons

The extension completes the mutation-enabled R=10^9 series and the mutation-free R=10^8 series. Lines connect the same seed block.

![Extension comparisons](figures/phase1-pilot-v210-m010-report/extension-return-comparisons.png)

## Quality control and computational feasibility

The analysis included 60 populations across 20 cells. The endpoint audit result was **PASS**. Summed runtime was 5.74 hours, diagnostic output was 4260.1 MiB, and maximum measured process-tree memory was 492.3 MiB. Every pilot run used two worker processes, five host passages and 500 steady within-host bacterial generations.

Each passage exchanged 10% of the one-billion-cell focal environment with the fixed regional source. The completion audit verified equal emigrant and immigrant totals for every population and passage.

Benchmark machine: HPC run: Linux-4.18.0-553.155.1.el8_10.x86_64-x86_64-with-glibc2.28. Timings should be recalibrated with a small batch on the target HPC before the larger experiment is launched.

![Computational feasibility](figures/phase1-pilot-v210-m010-report/computational-feasibility.png)

| # hosts H | Median runtime | Median output | Median peak RAM |
|---|---|---|---|
| 100 | 0.07 min | 0.17 MiB | 95.5 MiB |
| 1,000 | 0.34 min | 1.10 MiB | 96.5 MiB |
| 10,000 | 3.20 min | 10.78 MiB | 98.1 MiB |
| 100,000 | 31.54 min | 111.53 MiB | 102.6 MiB |

Within the mutation-free matched-return series, the fitted host-number exponent was 0.90 for runtime and 0.95 for diagnostic output. Thus, a tenfold increase in H multiplied median runtime by about 7.9 and output by about 8.8. Peak RAM was nearly independent of H because host tasks are streamed through a bounded queue.

At H=100, raising u from zero to the highest pilot level multiplied median runtime by 2.0, output by 54.9, and peak RAM by 1.15. This cost comes from materialized mutant lineages and their records.

For this implementation, a useful work model is O(P H [G S + M]), where P is the number of host passages, G is the number of within-host transitions, S is extant strain richness and M is the number of materialized mutation events. The code works with strain counts, not individual bacterial cells. Therefore K, e and R do not create one operation per cell. The independent effect of e cannot be estimated from this pilot because e changes inversely with H in the matched-return series.

**HPC starting allocation:** request two CPUs and 1 GiB RAM for each concurrent population. Begin with 8-16 concurrent populations, inspect CPU use and file-system throughput, then increase toward the CPU limit. With 100 available CPUs and two workers per population, 50 concurrent populations is the CPU ceiling. Rebenchmark before increasing workers, u, retained host detail, or the number of passages. At unchanged output settings, runtime and host-level output should grow roughly linearly with the number of host passages.

## Interpretation and next step

The pilot supports using a longer equilibrium and precision pilot before any confirmatory experiment. Three replicates per cell are sufficient for calibration and graphical ranges, but not for stable confidence intervals or equivalence tests.

## Reproducibility

- Source commit: `84938e61a3c1af84393e0ceb59cb4ce6748cc00a`
- Recorded source-tree SHA-256 values: `e278d154be073297b2fd80ec86909fbb9ce3b30ac12f101fd130d1a963fc61ab`
- Seed blocks: `sb0001, sb0002, sb0003`
- Analysis audit: `PASS`
- Complete numerical results remain in the machine-readable analysis tables.

## Appendix: glossary and diversity measures

- **Relative abundance, pᵢ:** fraction of environmental cells assigned to strain i.
- **Labelled richness, D0:** number of distinct strain labels. Every label counts once, even if it contains only one cell.
- **Shannon entropy, H':** H' = -Σᵢ pᵢ ln(pᵢ). Higher values mean that strain identity is less predictable. It is an index, not an effective number of strains.
- **Simpson diversity:** 1 - Σᵢ pᵢ². This is the probability that two independently sampled cells have different strain labels.
- **Hill D1:** D₁ = exp(H'). The number of equally abundant strains that would give the observed Shannon entropy.
- **Hill D2:** D₂ = 1 / Σᵢ pᵢ². The number of equally abundant strains that would give the observed probability of sampling the same strain twice. It gives more weight to abundant strains than D1.
- **Pielou evenness:** J = H' / ln(D₀). It compares the observed Shannon entropy with the maximum possible value for the observed richness.
- **Composition distance (TV):** TV(A,B) = ½ Σᵢ |pᵢ(A) - pᵢ(B)|. A value of 0.20 means that 20% of abundance must be reassigned among labels to make two compositions identical.
- **Cell:** one defined combination of experimental parameter values.
- **Seed block:** one independently seeded population replicate reused as a reproducible block across cells.
- **H:** number of hosts; **f:** escape fraction per host; **e:** escaping cells per host (fK); **R:** total returned cells (He); **alpha:** fraction of the pre-regulation mixture derived from hosts; **u:** whole-genome mutation probability per bacterial generation; **m:** fraction of the focal environment replaced from the fixed regional pool after each host return.
