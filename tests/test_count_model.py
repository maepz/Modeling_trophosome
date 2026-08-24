from __future__ import annotations

import unittest

import numpy as np

from trophosome.config import EvolutionConfig, HostConfig
from trophosome.count_model import (
    IdAllocator,
    LineageRecorder,
    PopulationState,
    fixed_pool_migration_step,
    free_living_selection_step,
    merge_populations,
    population_size_schedule,
    proportional_rescale,
    sample_population,
    simulate_within_host,
    subtract_population,
    wright_fisher_step,
)


class CountModelTests(unittest.TestCase):
    def setUp(self) -> None:
        self.population = PopulationState.from_counts([30, 70], [1.0, 1.0], [0.8, 1.2])

    def test_step_conserves_requested_population_size(self) -> None:
        evolution = EvolutionConfig(mutation_probability=0.0)
        next_state, mutants = wright_fisher_step(
            self.population,
            1_000_000,
            evolution,
            np.random.default_rng(1),
            IdAllocator(2),
        )
        self.assertEqual(next_state.size, 1_000_000)
        self.assertEqual(mutants, 0)

    def test_mutation_is_drawn_from_offspring_not_parental_size(self) -> None:
        evolution = EvolutionConfig(
            mutation_probability=1.0, max_materialized_mutants=100
        )
        recorder = LineageRecorder.from_founders(self.population.genotype_ids)
        next_state, mutants = wright_fisher_step(
            self.population,
            10,
            evolution,
            np.random.default_rng(2),
            IdAllocator(2),
            lineage_recorder=recorder,
        )
        self.assertEqual(mutants, 10)
        self.assertEqual(next_state.size, 10)
        self.assertEqual(next_state.richness, 10)
        self.assertEqual(len(recorder.events), 10)
        self.assertTrue(all(event.mutational_depth == 1 for event in recorder.events))

    def test_optimized_step_matches_reference_random_draws(self) -> None:
        evolution = EvolutionConfig(
            mutation_probability=0.2,
            mutation_effect_mean=-0.01,
            mutation_effect_sd=0.02,
            within_host_selection=True,
            max_materialized_mutants=1_000,
        )
        reference_rng = np.random.default_rng(991)
        weights = self.population.counts * self.population.within_host_fitness
        offspring = reference_rng.multinomial(200, weights / weights.sum())
        mutant_counts = reference_rng.binomial(
            offspring, evolution.mutation_probability
        )
        parent_index = np.repeat(np.arange(2), mutant_counts)
        expected_fitness = np.maximum(
            self.population.within_host_fitness[parent_index]
            + reference_rng.normal(
                evolution.mutation_effect_mean,
                evolution.mutation_effect_sd,
                size=int(mutant_counts.sum()),
            ),
            evolution.fitness_floor,
        )
        free_living_rng = np.random.default_rng(992)
        expected_free_living_fitness = np.maximum(
            self.population.free_living_fitness[parent_index]
            + free_living_rng.normal(
                evolution.mutation_effect_mean,
                evolution.mutation_effect_sd,
                size=int(mutant_counts.sum()),
            ),
            evolution.fitness_floor,
        )
        retained = offspring - mutant_counts

        observed, mutants = wright_fisher_step(
            self.population,
            200,
            evolution,
            np.random.default_rng(991),
            IdAllocator(2),
            free_living_fitness_rng=np.random.default_rng(992),
        )
        expected_counts = np.concatenate(
            (retained[retained > 0], np.ones(mutants, dtype=np.int64))
        )
        np.testing.assert_array_equal(observed.counts, expected_counts)
        np.testing.assert_allclose(
            observed.within_host_fitness[-mutants:], expected_fitness
        )
        np.testing.assert_allclose(
            observed.free_living_fitness[-mutants:], expected_free_living_fitness
        )

    def test_two_fitness_dimensions_control_only_their_own_habitat(self) -> None:
        population = PopulationState.from_counts(
            [50, 50],
            [0.0, 1.0],
            [1.0, 0.0],
        )
        within_host, _ = wright_fisher_step(
            population,
            100,
            EvolutionConfig(
                mutation_probability=0.0,
                within_host_selection=True,
            ),
            np.random.default_rng(4),
            IdAllocator(2),
        )
        np.testing.assert_array_equal(within_host.genotype_ids, [1])

        free_living = free_living_selection_step(
            population, 100, np.random.default_rng(4)
        )
        np.testing.assert_array_equal(free_living.genotype_ids, [0])

    def test_mutant_fitness_effects_are_independent_between_habitats(self) -> None:
        population = PopulationState.from_counts([1], [1.0], [1.0])
        recorder = LineageRecorder.from_founders(population.genotype_ids)
        observed, _ = wright_fisher_step(
            population,
            20,
            EvolutionConfig(
                mutation_probability=1.0,
                mutation_effect_mean=-0.01,
                mutation_effect_sd=0.02,
                max_materialized_mutants=20,
            ),
            np.random.default_rng(10),
            IdAllocator(1),
            lineage_recorder=recorder,
            free_living_fitness_rng=np.random.default_rng(11),
        )
        self.assertFalse(
            np.array_equal(observed.within_host_fitness, observed.free_living_fitness)
        )
        self.assertTrue(
            all(
                event.within_host_fitness != event.free_living_fitness
                for event in recorder.events
            )
        )

    def test_vectorized_merge_sums_shared_strains(self) -> None:
        other = PopulationState(
            np.asarray([0, 2]),
            np.asarray([5, 10]),
            np.asarray([1.0, 0.9]),
            np.asarray([0.8, 1.1]),
        )
        merged = merge_populations(self.population, other)
        np.testing.assert_array_equal(merged.genotype_ids, [0, 1, 2])
        np.testing.assert_array_equal(merged.counts, [35, 70, 10])
        np.testing.assert_allclose(merged.free_living_fitness, [0.8, 1.2, 1.1])

    def test_declared_zero_counts_preserve_aligned_strain_ids(self) -> None:
        population = PopulationState.from_counts(
            [30, 0, 70, 0],
            [1.0, 0.9, 1.1, 0.8],
            [0.8, 0.9, 1.0, 1.1],
        )
        np.testing.assert_array_equal(population.genotype_ids, [0, 2])
        np.testing.assert_array_equal(population.counts, [30, 70])
        np.testing.assert_allclose(population.within_host_fitness, [1.0, 1.1])
        np.testing.assert_allclose(population.free_living_fitness, [0.8, 1.0])

    def test_fixed_pool_migration_exchanges_equal_exact_counts(self) -> None:
        focal = PopulationState.from_counts(
            [30, 70, 0], [1.0, 1.0, 1.0], [0.9, 1.0, 1.1]
        )
        regional = PopulationState.from_counts(
            [0, 0, 100], [1.0, 1.0, 1.0], [0.9, 1.0, 1.1]
        )
        migrated, emigrants, immigrants = fixed_pool_migration_step(
            focal,
            regional,
            100,
            np.random.default_rng(7),
            np.random.default_rng(8),
        )
        self.assertEqual(migrated.size, focal.size)
        np.testing.assert_array_equal(migrated.genotype_ids, [2])
        np.testing.assert_array_equal(migrated.counts, [100])
        assert emigrants is not None
        assert immigrants is not None
        np.testing.assert_array_equal(emigrants.counts, [30, 70])
        np.testing.assert_array_equal(immigrants.genotype_ids, [2])
        self.assertEqual(emigrants.size, immigrants.size)

    def test_zero_migration_returns_an_unchanged_copy(self) -> None:
        migrated, emigrants, immigrants = fixed_pool_migration_step(
            self.population,
            self.population,
            0,
            np.random.default_rng(7),
            np.random.default_rng(8),
        )
        self.assertIsNone(emigrants)
        self.assertIsNone(immigrants)
        np.testing.assert_array_equal(
            migrated.genotype_ids, self.population.genotype_ids
        )
        np.testing.assert_array_equal(migrated.counts, self.population.counts)

    def test_proportional_rescale_has_exact_capacity(self) -> None:
        scaled = proportional_rescale(self.population, 33)
        self.assertEqual(scaled.size, 33)
        np.testing.assert_array_equal(scaled.counts, [10, 23])

    def test_proportional_rescale_randomizes_only_tied_boundary(self) -> None:
        population = PopulationState.from_counts([31, 69], [1.0, 1.0])
        first = proportional_rescale(population, 50, np.random.default_rng(1))
        second = proportional_rescale(population, 50, np.random.default_rng(2))
        self.assertEqual(first.size, 50)
        self.assertEqual(second.size, 50)
        self.assertIn(first.counts.tolist(), ([16, 34], [15, 35]))
        self.assertIn(second.counts.tolist(), ([16, 34], [15, 35]))
        self.assertFalse(np.array_equal(first.counts, second.counts))

    def test_without_replacement_sample_never_exceeds_source_counts(self) -> None:
        sample = sample_population(
            self.population, 90, np.random.default_rng(3), replace=False
        )
        source = dict(
            zip(
                self.population.genotype_ids.tolist(),
                self.population.counts.tolist(),
                strict=True,
            )
        )
        for genotype_id, count in zip(
            sample.genotype_ids.tolist(), sample.counts.tolist(), strict=True
        ):
            self.assertLessEqual(count, source[genotype_id])
        self.assertEqual(sample.size, 90)

    def test_finite_population_can_be_fully_sampled(self) -> None:
        sample = sample_population(
            self.population,
            self.population.size,
            np.random.default_rng(5),
            replace=False,
        )
        self.assertIsNone(subtract_population(self.population, sample))

    def test_billion_cell_population_can_be_sampled_without_replacement(
        self,
    ) -> None:
        population = PopulationState.from_counts(
            [500_000_000, 300_000_000, 200_000_000],
            [1.0, 1.0, 1.0],
        )
        sample = sample_population(
            population,
            10_000_000,
            np.random.default_rng(17),
            replace=False,
        )
        self.assertEqual(sample.size, 10_000_000)
        source = dict(
            zip(
                population.genotype_ids.tolist(),
                population.counts.tolist(),
                strict=True,
            )
        )
        for genotype_id, count in zip(
            sample.genotype_ids.tolist(), sample.counts.tolist(), strict=True
        ):
            self.assertLessEqual(count, source[genotype_id])

    def test_seeded_trajectory_is_reproducible_and_hits_capacity(self) -> None:
        host = HostConfig(
            population_size=1,
            infection_bottleneck=100,
            carrying_capacity=1_000,
            growth_factor=1.5,
            steady_generations=4,
            host_generations=1,
            escape_fraction=0.1,
        )
        evolution = EvolutionConfig(mutation_probability=0.001)
        first, first_summary = simulate_within_host(
            self.population,
            host,
            evolution,
            np.random.default_rng(44),
            IdAllocator(2),
        )
        second, second_summary = simulate_within_host(
            self.population,
            host,
            evolution,
            np.random.default_rng(44),
            IdAllocator(2),
        )
        self.assertEqual(first.size, 1_000)
        np.testing.assert_array_equal(first.genotype_ids, second.genotype_ids)
        np.testing.assert_array_equal(first.counts, second.counts)
        self.assertEqual(first_summary, second_summary)

    def test_endpoint_only_uses_same_trajectory(self) -> None:
        host = HostConfig(
            population_size=1,
            infection_bottleneck=100,
            carrying_capacity=1_000,
            growth_factor=1.5,
            steady_generations=4,
            host_generations=1,
            escape_fraction=0.1,
        )
        evolution = EvolutionConfig(mutation_probability=0.001)
        schedule = population_size_schedule(host)
        full, summaries = simulate_within_host(
            self.population,
            host,
            evolution,
            np.random.default_rng(71),
            IdAllocator(2),
            size_schedule=schedule,
        )
        endpoint, no_summaries = simulate_within_host(
            self.population,
            host,
            evolution,
            np.random.default_rng(71),
            IdAllocator(2),
            size_schedule=schedule,
            record_summaries=False,
        )
        np.testing.assert_array_equal(full.genotype_ids, endpoint.genotype_ids)
        np.testing.assert_array_equal(full.counts, endpoint.counts)
        self.assertTrue(summaries)
        self.assertEqual(no_summaries, ())


if __name__ == "__main__":
    unittest.main()
