# Distributional validation report

- Overall result: **PASS**
- Model family: `wright_fisher_counts`
- Model specification: `2.1.0`
- Software: `0.7.0`
- Output schema: `2.3.0`
- Seed: `20260810`
- Runtime: Python 3.14.7, NumPy 2.5.2
- Generated: 2026-08-24T02:03:04.512883+00:00

| Check | Result | Draws | Main diagnostic |
|---|---:|---:|---|
| neutral Wright-Fisher drift | PASS | 30,000 | observed_mean=24.05, expected_mean=24, mean_z=2.12, observed_variance=16.84, expected_variance=16.8, variance_z=0.2878 |
| fitness-weighted Wright-Fisher reproduction | PASS | 30,000 | observed_mean=27.68, expected_mean=27.69, mean_z=0.7419, observed_variance=14.86, expected_variance=14.91, variance_z=0.44 |
| free-living fitness-weighted reproduction | PASS | 30,000 | observed_mean=27.7, expected_mean=27.69, mean_z=0.1716, observed_variance=14.92, expected_variance=14.91, variance_z=0.05519 |
| independent dual-habitat mutation fitness effects | PASS | 30,000 | within_host_effect_mean=-0.01005, free_living_effect_mean=-0.01016, within_host_effect_variance=0.0004012, free_living_effect_variance=0.0004029, effect_correlation=-0.008313, maximum_z=1.44 |
| mutation-event count | PASS | 30,000 | observed_mean=4.997, expected_mean=5, mean_z=0.2776, observed_variance=4.464, expected_variance=4.5, variance_z=0.9724 |
| mutation-parent assignment | PASS | 30,000 | observed_mean=1.992, expected_mean=2, mean_z=1.054, observed_variance=1.899, expected_variance=1.92, variance_z=1.317 |
| with-replacement population sampling | PASS | 30,000 | observed_mean=4.993, expected_mean=5, mean_z=0.6062, observed_variance=4.012, expected_variance=4, variance_z=0.3761 |
| without-replacement population sampling | PASS | 30,000 | observed_mean=5, expected_mean=5, mean_z=0.04312, observed_variance=3.051, expected_variance=3.03, variance_z=0.8392 |
| optimized reservoir-founder sampling | PASS | 30,000 | observed_mean=1.597, expected_mean=1.6, mean_z=0.4746, observed_variance=1.277, expected_variance=1.28, variance_z=0.2437 |
| environmental apportionment label neutrality | PASS | 30,000 | observed_mean=29.99, expected_mean=30, mean_z=1.086, observed_variance=5.423, expected_variance=5.375, variance_z=1.09 |
| fixed-pool focal emigration | PASS | 30,000 | observed_mean=6.008, expected_mean=6, mean_z=0.7396, observed_variance=3.404, expected_variance=3.394, variance_z=0.3721 |
| fixed-pool regional immigration | PASS | 30,000 | observed_mean=3.976, expected_mean=4, mean_z=2.333, observed_variance=3.212, expected_variance=3.2, variance_z=0.4756 |
| multi-generation neutral drift distribution | PASS | 30,000 | total_variation_distance=0.006649, maximum_bin_z=1.794 |
| cell-level reference trajectory | PASS | 30,000 | richness_mean_difference=-0.004333, richness_z=0.3067, founder_cells_mean_difference=0.01447, founder_cells_z=0.4401, largest_mutant_clone_mean_difference=-0.0168, largest_mutant_clone_z=0.7201, richness_total_variation_distance=0.011, maximum_mean_z=0.7201 |

## Scope

These checks validate the declared stochastic law: neutral and selected within-host Wright–Fisher reproduction, independent dual-habitat fitness effects, free-living selection, Bernoulli infinite-alleles mutation, mutation parentage, reservoir founder sampling, finite escape sampling, fixed-pool emigration and immigration, multi-generation drift, and mutation-timing jackpot effects.

Hamilton environmental capacity regulation is checked for exact capacity, seeded label-neutral tie resolution, and no-return invariance when free-living selection is disabled. It is not yet compared with a biological equilibrium distribution; that belongs to the subsequent Phase 1 experiment.

Acceptance limits were declared in the validation program and are stored with every result. They are deliberately conservative six-standard-error release gates, supplemented by total-variation limits for complete distributions.
