# Phase 1 Stage 3 community-analysis input tables

The Stage 3 compiler creates one matched set of passage-100 tables plus a
trajectory table for the whole-community analyses:

- **X**, `x-explanatory-g100.tsv`: one row per analysis sample and the
  experimental parameters attached to it;
- **Y**, `y-ancestral-frequencies-g100.tsv`: the relative abundance of each of
  the 100 ancestral lineages in each sample;
- **Y′**, `yprime-tv-g100.tsv`: the square matrix of pairwise total-variation
  distances calculated from Y;
- **PRC trajectories**, `prc-ancestral-trajectories-g0-g100.tsv.gz`: ancestral-
  lineage frequencies for every analysis population at every passage from 0
  through 100.

The `sample_id` order is identical in X, the rows of Y, and both the rows and
columns of Y′. The audit fails rather than releasing a mismatched table set.

## Why there is one master matrix

One master matrix avoids maintaining several copies of the same table. The
`analysis_set` column in X identifies the three experimental subsets:

| `analysis_set` | Biological experiment | Rows in the master set | Primary db-RDA rows |
|---|---|---:|---:|
| `wave1_h_alpha_u` | Host abundance × feedback × mutation | 300 | 288 |
| `wave2a_h_by_b` | Host abundance × infection bottleneck | 144 | 144 |
| `wave2b_alpha_by_m` | Feedback × regional migration | 336 | 336 |

The 12 Wave 1 no-return populations are retained for paired comparisons but
have `include_primary_dbrda = false`. They are excluded from the primary Wave 1
factorial db-RDA because host abundance has no biological meaning when no
host-derived cells return.

Some original simulations were deliberately reused in a later experiment.
They therefore occur under a new analysis-cell identity, while
`source_run_id` records the original simulation. `source_alias_count` reveals
these cases. A source population is never duplicated within one analysis set.

Do not fit one db-RDA to the unfiltered master matrix. The three experiments
have different designs, and the master set contains cross-experiment aliases
of reused simulations. Subset X first, then use exactly those sample IDs to
subset Y and both axes of Y′.

## Definition of Y

At passage 100, environmental counts are converted to frequencies. In
mutation-enabled Wave 1 populations, each mutant is traced through
`strain_lineage_events.csv` to its frozen ancestral lineage. Its abundance is
then added to that lineage.

This gives 100 comparable biological columns in every population. A newly
created mutant ID cannot be used as a common column across simulations because
the same numeric ID in two runs does not represent the same mutation.

Y therefore describes changes in the abundance of the original lineages,
including all their mutant descendants. Mutation richness and the fate of
individual new strains remain separate mechanistic responses.

Every Y row sums to one.

The trajectory table uses the same biological response definition at every
passage. `population_sample_id` matches the endpoint `sample_id` in X and Y;
`trajectory_sample_id` uniquely identifies one population-passage observation.
It contains 101 rows per analysis population and is gzip-compressed because the
uncompressed table is large.

## Definition of Y′

For samples \(i\) and \(j\), the compiler calculates:

\[
D_{\mathrm{TV}}(i,j)=\frac{1}{2}\sum_{k=1}^{100}
|Y_{ik}-Y_{jk}|.
\]

Y′ is symmetric, has a zero diagonal, and ranges from zero to one. Because all
Y rows are frequency vectors summing to one, this distance is numerically
identical to Bray–Curtis dissimilarity.

## Creating the tables on the HPC

The raw environmental and lineage tables are stored in the machine-local
scratch directory, so the compiler must normally run on the HPC:

```bash
bash scripts/hpc/launch_phase1_stage3_wave2.sh --dbrda-only
```

This is a read-only analysis command. It does not launch simulations, change
checkpoints, or request a completion email. It reads the scratch location from
`experiments/work/trophosome/layout.local.json` and writes portable results to:

```text
experiments/work/trophosome/p01-neutral-feedback/analysis/
s03-parameter-map-dbrda-g100-derived/
```

The output directory also contains:

- `dbrda-source-provenance-g100.tsv`, linking analysis rows to original runs;
- `dbrda-input-audit-g100.json`, recording dimensions, definitions, checksums,
  and whether compilation passed.

Only this portable derived directory needs to be copied or committed. Do not
add raw scratch outputs to Git.

## Using the trajectory table for PRC

PRC requires repeated observations through time, so it cannot be reconstructed
from the passage-100 Y matrix. The extra long-format trajectory table is the
appropriate input: it avoids duplicating X/Y/Y′ for every possible comparison
while retaining all metadata required for reproducible subsetting.

Do not fit one PRC to every Stage 3 cell. Select a biologically coherent
comparison with one explicit reference treatment, retain all passages for each
selected population, and treat each population trajectory—not each passage—as
the independent experimental unit. The tutorial shows how to permute complete
trajectories within matched seed blocks.

## Reading and subsetting the master tables in R

```r
derived <- file.path(
  "experiments", "work", "trophosome", "p01-neutral-feedback", "analysis",
  "s03-parameter-map-dbrda-g100-derived"
)

X <- read.delim(
  file.path(derived, "x-explanatory-g100.tsv"),
  check.names = FALSE,
  stringsAsFactors = FALSE
)
Y_table <- read.delim(
  file.path(derived, "y-ancestral-frequencies-g100.tsv"),
  check.names = FALSE,
  stringsAsFactors = FALSE
)
D_table <- read.delim(
  file.path(derived, "yprime-tv-g100.tsv"),
  check.names = FALSE,
  stringsAsFactors = FALSE
)

stopifnot(identical(X$sample_id, Y_table$sample_id))
stopifnot(identical(X$sample_id, D_table$sample_id))
stopifnot(identical(X$sample_id, names(D_table)[-1]))

Y <- as.matrix(Y_table[, -1, drop = FALSE])
rownames(Y) <- Y_table$sample_id
D_TV <- as.matrix(D_table[, -1, drop = FALSE])
rownames(D_TV) <- D_table$sample_id

primary_flag <- tolower(as.character(X$include_primary_dbrda)) == "true"

subset_stage3 <- function(set_name) {
  keep <- X$analysis_set == set_name & primary_flag
  ids <- X$sample_id[keep]
  list(
    X = X[keep, , drop = FALSE],
    Y = Y[ids, , drop = FALSE],
    D_TV = as.dist(D_TV[ids, ids, drop = FALSE])
  )
}

wave1 <- subset_stage3("wave1_h_alpha_u")
wave2a <- subset_stage3("wave2a_h_by_b")
wave2b <- subset_stage3("wave2b_alpha_by_m")

stopifnot(nrow(wave1$X) == 288)
stopifnot(nrow(wave2a$X) == 144)
stopifnot(nrow(wave2b$X) == 336)
```

For a Hellinger RDA, transform the matching Y subset with
`vegan::decostand(Y, method = "hellinger")`. For the primary db-RDA, use the
matching `D_TV` object and condition or restrict permutations by
`seed_block_id`.

The complete, annotated workflow is in
[`phase1-stage3-community-analysis.Rmd`](phase1-stage3-community-analysis.Rmd).
