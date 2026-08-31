# Phase 1: Stage 3, part one

Frozen on 2026-08-31. This replaces the **unlaunched** 18-cell, six-seed,
250-passage proposal. It is exploratory, not a confirmatory experiment or an
equilibrium test. No production simulations were launched during preparation.

## Biological question and experiment size

Does increasing host abundance offset the change in environmental strain
composition caused by stronger host feedback? Does within-host mutation alter
that relationship?

- Host abundance H: **100; 1,000; 10,000**.
- Feedback alpha targets: **0; 0.001; 0.01; 0.1; 0.99**.
- Whole-genome mutation probability u: **0 or 10^-10 per bacterial generation**.
- One frozen 100-strain starting population, `ip001-fisher100`.
- **12 matched seed blocks**, `sb0001`-`sb0012` (master seeds 666-677).
- **100 host passages** per new population.

There are **24 new positive-feedback cells** (`c0027`-`c0050`), hence **288 new
populations and 28,800 newly simulated passages**. Reusing the twelve Stage 2
`c0021` no-return populations gives **25 primary conditions and 300 primary
populations**. The reused control is one condition with twelve replicates,
not twelve independent controls for every treatment.

Five other Stage 2 cells contribute **60 supplementary reference outcomes**.
They are not part of the 25-cell grid and do not increase its replication.

## Exact matrix and whole-cell returns

All hosts finish at K=10^9 bacteria; N_E=10^9. If each host releases e cells,
the total return is R=H*e and alpha=R/(N_E+R), **before regional migration**.
To compare H at exactly the same total return, round the requested R to the
nearest multiple of 10,000 (ties to even), then set e=R/H and f=e/K.
This avoids slightly different total returns between host abundances.

| Alpha target | Exact R, all H | Realized alpha |
|---:|---:|---:|
| 0.001 | 1,000,000 | 0.000999000999 |
| 0.01 | 10,100,000 | 0.009999009999 |
| 0.1 | 111,110,000 | 0.099999099999 |
| 0.99 | 99,000,000,000 | 0.99 |

Each row below is run twice, once for each mutation setting. IDs in the two
cell columns correspond to u=0 and u=10^-10, respectively.

| H | Alpha target | Release fraction f | u=0 cell | u=10^-10 cell | Escape-range interpretation |
|---:|---:|---:|---|---|---|
| 100 | 0.001 | 0.00001 | c0027 | c0028 | Primary |
| 100 | 0.01 | 0.000101 | c0029 | c0030 | Primary |
| 100 | 0.1 | 0.0011111 | c0031 | c0032 | Primary |
| 100 | 0.99 | 0.99 | c0033 | c0034 | Above primary range: 99% release |
| 1,000 | 0.001 | 0.000001 | c0035 | c0036 | Below primary range |
| 1,000 | 0.01 | 0.0000101 | c0037 | c0038 | Primary |
| 1,000 | 0.1 | 0.00011111 | c0039 | c0040 | Primary |
| 1,000 | 0.99 | 0.099 | c0041 | c0042 | Above primary range: 9.9% release |
| 10,000 | 0.001 | 0.0000001 | c0043 | c0044 | Below primary range |
| 10,000 | 0.01 | 0.00000101 | c0045 | c0046 | Below primary range |
| 10,000 | 0.1 | 0.000011111 | c0047 | c0048 | Primary |
| 10,000 | 0.99 | 0.0099 | c0049 | c0050 | Primary |

The original plausible escape range was 10^-5 to 10^-2. Five H-alpha pairs
(ten mutation-specific cells) lie outside it and are labelled extended-range
mechanistic tests. Selection remains neutral in every cell.

The complete machine-readable matrix, including the reused control, is
[`phase1-stage3-wave1-v210-m010-g100-cells.tsv`](../experiments/work/trophosome/p01-neutral-feedback/design/phase1-stage3-wave1-v210-m010-g100-cells.tsv).
Both requested and realized alpha are recorded. `show-cell c0034` in the
project-layout helper resolves the new definition. The replaced c0027-c0044
definitions were never launched and are preserved in the old Git snapshot;
no pilot cell IDs or pilot records are reassigned.

## Fixed settings

Scientific model 2.1.0, software 0.7.0, output schema 2.3.0. Both selection
switches are false; all fitness values are 1 and mutation fitness effects are
0. The effective reservoir is non-depleting during infection. All escapees
are pooled; capacity is regulated neutrally by Hamilton apportionment.

B=10 founders, K=N_E=10^9, capacity ratio c=1, growth factor 1.2, and 500
steady bacterial generations remain fixed. There is no free-living generation
or mutation. The fixed regional source replaces m=0.1 of the focal population
after host return and capacity regulation. Its composition equals the frozen
initial vector. Thus the expected retained host contribution is (1-m)*alpha.
The lifetime offspring sum is 506,775,749,305: u=10^-10 gives approximately
50.68 mutation events per host lifetime. More hosts increase total mutation
supply as well as pooling; the u=0 comparisons isolate redistribution.

## Reuse of Stage 2

Reuse exactly passage **100**, never passage 250, for sb0001-sb0012:

| Cell | H | Actual alpha | u | Role |
|---|---:|---:|---:|---|
| c0021 | 100 | 0 | 0 | Shared primary no-return control |
| c0022 | 100 | 0.5 | 0 | Supplementary |
| c0023 | 100 | 0.5 | 10^-10 | Supplementary |
| c0024 | 10,000 | 0.5 | 0 | Supplementary |
| c0025 | 100 | 0.090909090909 | 0 | Supplementary, not alpha=0.1 |
| c0026 | 100 | 0.909090909091 | 0 | Supplementary, not alpha=0.99 |

The generator verifies the passed Stage 2 analysis and trajectory checksum,
then freezes 72 endpoints and passages 51-100 mean/SD summaries with provenance.
Later report rebuilds do not require mutable Stage 2 analysis files or access
to its raw scratch. These are **derived summaries**, not new labelled-count
archives: retain the original Stage 2 raw records separately if needed.
Reused seed blocks are exploratory and must not serve as held-out Stage 4
confirmation. No shared control is duplicated into extra independent data.

## Registered analysis

The primary outcome is **TV distance from the common starting population at
passage 100**. Supporting outcomes are D0, Hill D1/D2, evenness, Shannon and
Simpson (1 - sum p_i^2). Root-collapsed diversity, mutation abundance/richness
and all new-run trajectories are retained. No cross-run mutant IDs are aligned.
A TV contrast is a difference in departure from the start, not direct distance
between two independently mutated communities.

Compare adjacent H levels at fixed alpha and u; adjacent positive alpha levels
at fixed H and u; mutation on/off at every H-alpha pair; and each new condition
against the shared control. Calculate differences within each seed before
averaging. Relative changes apply to D0/D1/D2; evenness and TV use absolute
differences. Twelve paired values, not hosts or passages, define uncertainty.

The primary mechanistic interaction is a TV **difference-in-differences**:
subtract the feedback effect at smaller H from the feedback effect at larger H,
within each seed. Negative values mean the feedback effect is smaller at larger
H. Repeat in both mutation settings, then compare the two interactions. These
12 H-by-alpha and six mutation-modification contrasts are reported with 90%
Student-t intervals (df=11); they have no separately agreed biological margin.

Other contrasts use exploratory 90% intervals and working margins of 5% for
D0/D1/D2, 0.02 for evenness, and 0.05 for TV. An interval wholly within the
margin supports practical equivalence; wholly beyond a margin supports a
meaningful increase/decrease; otherwise the result is uncertain. No multiplicity
adjustment is made and no confirmatory claims should be based on these screens.
Shannon/Simpson are descriptive; no new biological margins are invented.

Passages 51-100 provide supporting within-population means and SDs. A paired
mean-TV change (76-100 minus 51-75) is compared with +/-0.05 to flag cells for
longer-run review. This limited mean-drift check does **not** prove stationarity,
constant variability or equilibrium. Inspect individual trajectories too.
Add seeds where precision is insufficient; add intermediate cells where effects
change sharply; consider longer validation runs where time dependence remains
unresolved. Additional runs require review, not automatic adaptive launching.

## Launch safeguards, output and portability

Use the [copy-paste HPC workflow](../scripts/hpc/README.md#phase-1-stage-3-first-mapping-wave).
The script paths retain `stage3_wave1` for compatibility; the experiment variant
and scratch directory are now **v210-m010-g100**. Never mix g250 scratch with it.
`--verify` checks 317 generated files and the new registry entries. The launcher
requires a clean committed source and verifies the actual imported package.

First run **c0034, c0049, c0050 with sb0001**, included in the 288. These cover
99% per-host release with mutation at H=100, and H=10,000 with mutation off/on.
Only after their audited resource check passes can the other 285 runs begin.
The check interpolates mutation-enabled costs between the two H extremes,
scales mutation-free costs by H, and applies a 2x allowance. It requires less
than 350 GiB projected output, enough free scratch, and at most 48 projected
hours **per population**, not a guaranteed whole-batch runtime. Check user
quota and allocation separately. A resumed safety run requires manual timing
review because the recorded duration covers only its latest attempt.

Default concurrency is eight populations, each with two host workers; numeric
library threads are capped at one. Checkpoints target one hour at completed
passage boundaries, keeping two recovery points while active. Valid completed
runs are audited and skipped; interrupted runs resume under unchanged inputs.

New populations retain all 101 environmental compositions, infection founders,
pooled abundance/occupancy, releases, migration, lineage events, adult summaries
and provenance. Adult strain counts are full at H=100 and use a reproducible
100-host panel at larger H. These retention choices do not change dynamics.
The report preserves compact trajectories and the final-state checksum. It
does not automatically delete raw files; audit and archive before implementing
final-state-only environmental retention. The reused summaries remain clearly
distinguished from the new complete raw records.

The automatic or stand-alone report audits all 288 new completions plus frozen
references before writing PDF, Markdown, figures and derived TSVs. Report-only
mode never simulates. No production result is generated during preparation;
end-to-end tests use explicitly synthetic data in temporary directories.
