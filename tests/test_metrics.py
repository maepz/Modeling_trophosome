from __future__ import annotations

import unittest

from trophosome.count_model import PopulationState
from trophosome.metrics import population_metrics


class MetricsTests(unittest.TestCase):
    def test_single_genotype_has_defined_zero_diversity(self) -> None:
        metrics = population_metrics(PopulationState.from_counts([10], [0.8], [1.2]))
        self.assertEqual(metrics.gene_diversity, 0.0)
        self.assertEqual(metrics.shannon_diversity, 0.0)
        self.assertEqual(metrics.pielou_evenness, 1.0)
        self.assertEqual(metrics.mean_within_host_fitness, 0.8)
        self.assertEqual(metrics.mean_free_living_fitness, 1.2)

    def test_equal_two_genotypes_have_maximal_evenness(self) -> None:
        metrics = population_metrics(PopulationState.from_counts([50, 50], [1, 1]))
        self.assertAlmostEqual(metrics.pielou_evenness, 1.0)
        self.assertAlmostEqual(metrics.gene_diversity, 100 / 99 * 0.5)


if __name__ == "__main__":
    unittest.main()
