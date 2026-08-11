# Bounded benchmark results (2026-08-10)

These measurements are diagnostic, not publication benchmarks. They were run on
the local macOS x86_64 environment with Python 3.11, NumPy 2.4.1, and NetworkX
3.3. Peak memory is Python allocation traced by `tracemalloc`, not full process
resident memory.

## V3.1 endpoint profile

V3.1 is the legacy endpoint comparator. A one-host profile used a 10-cell inoculum,
carrying capacity `10^6`, growth factor 1.2, 50 adult generations and mutation
probability `10^-4`.

| Quantity | Result |
|---|---:|
| Calculated growth generations | 64 |
| Total lineage generations | 114 |
| Inferred mutation-derived cells | approximately 11,336 |
| Wall time | 1.91 s |
| Final adult graph | 908 nodes, 906 edges |

The leading cumulative costs were the cell-by-cell CRP (0.98 s), exact
Stirling-number founder probabilities (0.51 s), and NetworkX conversion (0.14 s).

Increasing mutation probability to `10^-3` implies approximately 107,793
mutation-derived cells but an expected Ewens richness of only about 4,690 strains.
The legacy function was still in the exact Stirling-number calculation after 90
seconds and was interrupted. This is a diagnostic upper-bound observation, not a
completed timing result.

Screen a proposed V3.1 run without simulating it using:

```bash
PYTHONPATH=src python scripts/estimate_v3_1_load.py \
  --infection-size 10 \
  --carrying-capacity 1000000 \
  --growth-factor 1.2 \
  --steady-generations 50 \
  --mutation-probability 0.001
```

## Recorded full-matrix version comparison

The pulled scheduler outputs compare V1.4, V2.2 and V3.1 over several host and
worker counts.

| Version | Full matrix elapsed time | Peak resident memory |
|---|---:|---:|
| V3.1 | 9 min 50 s | 15.96 GB |
| V1.4 | 39 min 58 s | 67.31 GB |
| V2.2 | 1 h 44 min 46 s | 7.64 GB |

The job-script defaults were carrying capacity 10,000, mutation probability
`10^-12`, no adult steady phase, 20 host generations and 100 initial strains.
The comparison supports V3.1's low-mutation speed advantage. It does not test
biologically demanding mutation supply. One earlier V3.1 job failed because of
an undefined return value; the later pulled implementation fixes that particular
error.

## Legacy single-host profile

The legacy NetworkX path used a 10-cell inoculum, growth factor 1.2, 50 adult
generations, and saved every bacterial generation.

| Carrying capacity | Mutation probability | Time (s) | Traced peak (MB) | Stored snapshots | Total nodes across snapshots |
|---:|---:|---:|---:|---:|---:|
| 10,000 | 0.001 | 0.86 | 3.04 | 92 | 4,218 |
| 10,000 | 0.01 | 6.15 | 23.70 | 92 | 34,156 |
| 100,000 | 0.001 | 6.26 | 24.81 | 104 | 34,370 |

The last two scenarios have similar mutation supply and nearly identical cost.
The profile attributed most cumulative time to graph copying, graph updates, and
`remove_empty_leaves`.

## Exact count-kernel scaling probe

The array kernel held total census size at one billion cells for 20 neutral
generations with no new mutations.

| Represented genotypes | Time (s) | Traced peak (MB) |
|---:|---:|---:|
| 10 | 0.0047 | 0.017 |
| 100 | 0.0050 | 0.022 |
| 1,000 | 0.0117 | 0.129 |
| 10,000 | 0.0505 | 1.236 |

This demonstrates the intended scaling property: census size is an integer in the
multinomial draw, while work grows with represented richness. It does not remove
the mutation-supply limit. With one million cells, 100 starting genotypes,
mutation probability `1e-4`, and 20 generations, the kernel materialized 1,979
mutants, ended with 560 extant genotypes, and used 0.078 MB traced peak memory in
0.015 seconds. Biologically large cases must still be screened using `N*mu` and
the total number of host-generation transitions.

Reproduce the array benchmark with:

```bash
python scripts/benchmark_count_kernel.py \
  --richness 10,100,1000,10000 \
  --population-size 1000000000 \
  --generations 20 \
  --output results/count-kernel.csv
```

## Exact prototype versus V3.1 and the legacy neutral count model

A toy end-to-end compute benchmark compares the maintained exact-count prototype
with V3.1 and the graph-based V1.3-style model. V1.3 itself hard-codes selection,
so the neutral legacy line uses V1.4's documented selection switch, which
otherwise retains the V1.3 graph and trajectory implementation. V3.1 is the
legacy neutral endpoint approximation and does not reproduce the same forward
transition law as the other two models.

Each timing covers one host generation: infection, within-host growth, five
steady generations, escape, host aggregation and environmental return. File
output is excluded; worker-pool startup is included. The baseline uses 8 hosts,
infection size 10, growth factor 1.5, `K=4,000`, `u=0.001`, escape fraction 0.01
and one worker. Each panel changes only its named variable, and values below are
medians of three runs.

| Scaling variable | Exact-count prototype | V3.1 endpoint | Legacy neutral path |
|---|---:|---:|---:|
| Hosts: 1 to 16 | 0.0045 to 0.0487 s | 2.395 to 2.860 s | 2.397 to 2.715 s |
| Hosts: 100, 1,000 | 0.237, 3.916 s | 4.934, 31.457 s | 3.779, 18.066 s |
| Mutation rate: 0 to 0.005 | 0.0080 to 0.0186 s | fails at 0; 2.437 to 2.868 s thereafter | 2.810 to 2.288 s |
| `K`: 250 to 16,000 | 0.0072 to 0.0274 s | 3.627 to 4.710 s | 2.156 to 2.400 s |
| `K = 10^6` | 0.345 s | preflight guard | 42.758 s |
| `K = 10^8` | preflight guard | preflight guard | preflight guard |
| Workers: 1, 2, 4 | 0.0182, 0.391, 0.422 s | 2.232, 2.849, 4.303 s | 1.813, 2.486, 3.849 s |

The worker result is deliberately not presented as parallel speedup: this toy
workload is too small to amortize process startup and serialization. One-worker
execution is the correct choice in that regime. The host, mutation and capacity
panels show the new representation's much lower baseline overhead and the
expected increase as more hosts, mutation events and extant strains are
materialized. Both legacy curves are dominated by pool/import and graph or tree
construction costs at this toy scale, so a larger, separately budgeted benchmark
is required to locate their computational scaling slopes. V3.1's zero-mutation
point is recorded as a failed run (`ValueError: Zero samples specified`), not as
a runtime measurement. Guarded points are also not runtime measurements. At the
baseline mutation rate and eight hosts, `K=10^8` implies about 6.67 million
forward mutation events for the count paths and 35.2 million mutation-derived
endpoint cells for V3.1. Those estimates exceed the benchmark's declared
one-million-event and 100,000-cell materialisation limits. V3.1 also crosses its
limit at `K=10^6` (about 267,560 mutation-derived endpoint cells), whereas the
two forward-count paths remain below their guard and were timed.

Reproduce the comparison with:

```bash
PYTHONPATH=src python scripts/benchmark_v1_3_vs_exact.py \
  --output results/benchmarks/v1_3-v3_1-vs-exact-toy.csv \
  --repeats 3
```

The added large-grid points can be requested with:

```bash
PYTHONPATH=src python scripts/benchmark_v1_3_vs_exact.py \
  --output results/benchmarks/v1_3-v3_1-vs-exact-large-grid.csv \
  --repeats 3 \
  --dimensions individual_hosts,within_host_population_size \
  --hosts 100,1000 \
  --carrying-capacities 1000000,100000000
```

## Historical versus matched-neutral regimes

The slide-13 result and the later one-generation toy result answer different
questions. The recorded slide-13 capacity sweep used 100 hosts, 10 workers,
`u=10^-12`, 20 host generations and a 100-strain Fisher-log reservoir. It also
gave selected V1.3 50 steady generations but V3.1 zero. Those recorded values
are preserved in
`examples/investigate_scaling_up/slide13_capacity_runtime.csv`: V3.1 required
3.25--3.99 seconds while selected V1.3 rose from 5.59 to 26.97 seconds over
`K=10^4` to `10^8`.

A new matched-neutral control holds the reservoir, host count, worker count,
growth, steady phase and number of host generations constant across the neutral
V1.4 compatibility path, V3.1 and the exact-count prototype. At `u=10^-12`, all
three paths remain broadly flat with capacity. Median runtimes over `K=10^4` to
`10^8` were 9.47--13.56 seconds for V3.1, 10.23--21.46 seconds for the neutral
legacy path and 2.08--4.26 seconds for the exact-count prototype. This confirms
that the historical V3.1 advantage at negligible mutation supply is real, while
separating it from the historical steady-phase mismatch.

The matched mutation-supply profile fixes `K=100,000` and varies `u`. V3.1 and
the neutral legacy path were nearly tied through an expected endpoint mutation
load of about 500 cells. At about 49,878 mutation-derived endpoint cells, V3.1
became slightly slower (10.42 versus 9.74 seconds). At about 487,944 cells, the
V3.1 run was preflight-rejected rather than reported as a runtime; the exact and
neutral forward paths completed in 2.39 and 14.18 seconds, respectively.

Run the matched controls with:

```bash
python scripts/benchmark_scalability_regimes.py \
  --profiles matched-neutral-capacity,matched-neutral-mutation \
  --output-dir results/benchmarks/regime-comparison \
  --repeats 3
```

The following command also reruns the faithful historical profile. It is much
more expensive because it executes 20 host generations per point:

```bash
python scripts/benchmark_scalability_regimes.py \
  --profiles historical-slide13,matched-neutral-capacity,matched-neutral-mutation
```
