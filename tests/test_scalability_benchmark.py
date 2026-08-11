import unittest

from scripts.benchmark_scalability_regimes import _profile_arguments
from scripts.benchmark_v1_3_vs_exact import (
    EXACT_LABEL,
    V3_1_LABEL,
    Case,
    _endpoint_mutation_load,
    _expected_mutation_workload,
    _legacy_environment,
)


class ScalabilityBenchmarkTests(unittest.TestCase):
    def workload(self, carrying_capacity: int, model: str) -> float:
        workload, _ = _expected_mutation_workload(
            Case(
                workers=1,
                hosts=8,
                mutation_probability=0.001,
                carrying_capacity=carrying_capacity,
            ),
            model=model,
            bottleneck=10,
            growth_factor=1.5,
            steady_generations=5,
            escape_fraction=0.01,
        )
        return workload

    def test_one_million_cell_case_only_crosses_v3_guard(self) -> None:
        self.assertLess(self.workload(1_000_000, EXACT_LABEL), 1_000_000)
        self.assertGreater(self.workload(1_000_000, V3_1_LABEL), 100_000)

    def test_one_hundred_million_cell_case_crosses_both_guards(self) -> None:
        self.assertGreater(self.workload(100_000_000, EXACT_LABEL), 1_000_000)
        self.assertGreater(self.workload(100_000_000, V3_1_LABEL), 100_000)

    def test_workload_scales_with_host_generations(self) -> None:
        case = Case(1, 8, 0.001, 4_000)
        one_generation = _endpoint_mutation_load(
            case,
            bottleneck=10,
            growth_factor=1.5,
            steady_generations=5,
            escape_fraction=0.01,
        )
        twenty_generations = _endpoint_mutation_load(
            case,
            bottleneck=10,
            growth_factor=1.5,
            steady_generations=5,
            escape_fraction=0.01,
            host_generations=20,
        )
        self.assertAlmostEqual(twenty_generations, 20 * one_generation)

    def test_fisher_log_environment_is_seeded_and_has_requested_richness(self) -> None:
        first = _legacy_environment(100, 10, initial_strains=100, seed=666)
        second = _legacy_environment(100, 10, initial_strains=100, seed=666)
        first_rows = list(first.nodes(data=True))
        second_rows = list(second.nodes(data=True))
        self.assertEqual(len(first_rows), 100)
        self.assertEqual(first_rows, second_rows)

    def test_regime_presets_separate_historical_and_matched_semantics(self) -> None:
        historical = _profile_arguments("historical-slide13")
        matched = _profile_arguments("matched-neutral-capacity")
        self.assertIn("v1_3,v3_1", historical)
        self.assertIn("20", historical)
        self.assertIn("legacy,v3_1,exact", matched)
        self.assertIn("1", matched)


if __name__ == "__main__":
    unittest.main()
