# Phase 1 Stage 3 Wave 3: bridge design

Wave 3 fills the most important gaps between the earlier factorial panels. Its
purpose is to estimate interactions that Waves 1 and 2 could not separate
cleanly, especially feedback × bottleneck, host abundance × migration, and
bottleneck × migration.

The frozen batch contains the proposed **14 core bridge conditions plus the
two-condition extension**. All 16 conditions are active. Each is run with six
matched seed blocks (`sb0001`--`sb0006`) to passage 100, giving **96 new
simulated populations**.

## Frozen conditions

| Bridge | Cell | H (hosts) | Feedback target alpha | B (founders/host) | m (regional exchange) | Role |
|---|---|---:|---:|---:|---:|---|
| A1 | c0091 | 1,000 | 0.01 | 1 | 0.1 | alpha × B core |
| A2 | c0092 | 1,000 | 0.01 | 50 | 0.1 | alpha × B core |
| A3 | c0093 | 1,000 | 0.1 | 1 | 0.1 | alpha × B core |
| A4 | c0094 | 1,000 | 0.1 | 50 | 0.1 | alpha × B core |
| A5 | c0095 | 1,000 | 0.99 | 1 | 0.1 | alpha × B core |
| A6 | c0096 | 1,000 | 0.99 | 50 | 0.1 | alpha × B core |
| B1 | c0097 | 100 | 0.1 | 10 | 0.01 | H × m core |
| B2 | c0098 | 100 | 0.1 | 10 | 0.9 | H × m core |
| B3 | c0099 | 10,000 | 0.1 | 10 | 0.01 | H × m core |
| B4 | c0100 | 10,000 | 0.1 | 10 | 0.9 | H × m core |
| C1 | c0101 | 1,000 | 0.1 | 1 | 0.01 | B × m core |
| C2 | c0102 | 1,000 | 0.1 | 1 | 0.9 | B × m core |
| C3 | c0103 | 1,000 | 0.1 | 50 | 0.01 | B × m core |
| C4 | c0104 | 1,000 | 0.1 | 50 | 0.9 | B × m core |
| D1 | c0105 | 10,000 | 0.01 | 1 | 0.1 | matched-HB extension |
| D2 | c0106 | 10,000 | 0.99 | 1 | 0.1 | matched-HB extension |

Here, H is the number of hosts, B is the number of infecting founder cells per
host, alpha is the fraction of the focal reservoir contributed by host-derived
cells immediately after return and before regional exchange, and m is the
fraction replaced from the fixed regional pool.

## Parameters held constant

- The focal reservoir and within-host carrying capacity are both (10^9)
  cells.
- The regional pool has the same frozen 100-lineage composition as the focal
  population at passage 0.
- Within-host growth lasts 500 steady bacterial generations per host passage.
- Mutation is off (`u = 0`) and selection is off in both habitats.
- Complete environmental compositions are retained at passages 0--100.
- Total host return is set exactly from alpha, then divided equally among H
  hosts. No condition exceeds 10,000 hosts.

These constraints make Wave 3 a neutral bridge experiment. They allow the new
conditions to be combined with the matching neutral populations from Waves 1
and 2 without changing the biological meaning of the response.

## Planned combined model

The bridge-augmented endpoint analysis uses the six common seed blocks and the
following predeclared terms:

```text
H + alpha + B + m + H:alpha + H:B + alpha:m + alpha:B + H:m + B:m
```

H and B should be analysed on log10 scales and continuous predictors should be
centred and scaled before fitting. The 16 Wave 3 conditions are not intended as
a stand-alone full factorial. They augment the compatible parts of Waves 1 and
2. Exact simulations reused between experimental panels must occur only once in
the fitted data; `source_run_id` provides the deduplication key.

The design is a biologically constrained bridge augmentation. It is not a claim
that an unconstrained search found a globally optimal 16-condition design.

## Safety and reproducibility

Three planned populations provide one resource measurement for each H stratum:
c0098 (H=100), c0096 (H=1,000), and c0100 (H=10,000), all with `sb0001`.
They remain part of the 96 populations. The launcher audits them, applies a 2×
resource margin, resumes valid checkpoints, and refuses a full launch until the
safety assessment passes.

The generated cell table, run manifest, checksums, seeds, and complete
configuration files are version controlled. Raw model outputs remain in the
machine-local HPC scratch tree. See the
[HPC workflow](../scripts/hpc/README.md#phase-1-stage-3-wave-3-bridge-experiment)
for the copy-paste commands.
