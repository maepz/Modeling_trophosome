from __future__ import annotations

import unittest

from trophosome import (
    MODEL_FAMILY,
    MODEL_SPEC_VERSION,
    OUTPUT_SCHEMA_VERSION,
    __version__,
)


class PackageMetadataTests(unittest.TestCase):
    def test_release_versions_are_declared(self) -> None:
        self.assertEqual(__version__, "0.7.0")
        self.assertEqual(MODEL_FAMILY, "wright_fisher_counts")
        self.assertEqual(MODEL_SPEC_VERSION, "2.1.0")
        self.assertEqual(OUTPUT_SCHEMA_VERSION, "2.3.0")


if __name__ == "__main__":
    unittest.main()
