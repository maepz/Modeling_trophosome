from __future__ import annotations

import unittest

from scripts.validate_distributions import run_validation


class DistributionalValidationTests(unittest.TestCase):
    def test_release_gating_distributions(self) -> None:
        failures = [
            result.name
            for result in run_validation(repetitions=2_500, seed=20260810)
            if not result.passed
        ]
        self.assertEqual(failures, [])


if __name__ == "__main__":
    unittest.main()
