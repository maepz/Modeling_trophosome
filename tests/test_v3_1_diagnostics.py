from __future__ import annotations

import math
import unittest

from trophosome.v3_1_diagnostics import (
    estimate_v3_1_load,
    ewens_expected_richness_bounds,
    growth_generations,
    mutation_carrier_probability,
)


class V31DiagnosticsTests(unittest.TestCase):
    def test_growth_duration_matches_legacy_formula(self) -> None:
        expected = math.ceil(math.log(1_000_000 / 10) / math.log(1.2))
        self.assertEqual(growth_generations(10, 1_000_000, 1.2), expected)

    def test_stable_mutation_probability_at_tiny_rate(self) -> None:
        observed = mutation_carrier_probability(1e-12, 10)
        self.assertAlmostEqual(observed, 1e-11, delta=1e-20)

    def test_zero_mutation_is_reported_without_crashing(self) -> None:
        estimate = estimate_v3_1_load(
            infection_size=10,
            carrying_capacity=1_000,
            growth_factor=2,
            steady_generations=0,
            mutation_probability=0,
        )
        self.assertEqual(estimate.mutant_cells_for_screening, 0)
        self.assertTrue(any("zero-mutant" in warning for warning in estimate.warnings))

    def test_ewens_bounds_contain_exact_expected_richness(self) -> None:
        n = 100
        theta = 2.0
        exact = sum(theta / (theta + i) for i in range(n))
        lower, upper = ewens_expected_richness_bounds(n, theta)
        self.assertLessEqual(lower, exact)
        self.assertGreaterEqual(upper, exact)

    def test_mutation_supply_warning_uses_lineage_exposure(self) -> None:
        estimate = estimate_v3_1_load(
            infection_size=10,
            carrying_capacity=1_000_000,
            growth_factor=1.2,
            steady_generations=50,
            mutation_probability=0.001,
        )
        self.assertGreater(estimate.expected_mutant_cells, 100_000)
        self.assertTrue(any("current CRP" in warning for warning in estimate.warnings))


if __name__ == "__main__":
    unittest.main()
