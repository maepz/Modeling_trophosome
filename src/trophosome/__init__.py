"""Main exact-count Wright--Fisher model for trophosome symbionts."""

from importlib.metadata import PackageNotFoundError, version

from trophosome.config import ModelConfig, load_config
from trophosome.count_model import PopulationState, simulate_within_host

MODEL_FAMILY = "wright_fisher_counts"
MODEL_SPEC_VERSION = "2.0.0"
OUTPUT_SCHEMA_VERSION = "2.1.0"

try:
    __version__ = version("trophosome-model")
except PackageNotFoundError:  # pragma: no cover - source tree without installation
    __version__ = "0+unknown"

__all__ = [
    "MODEL_FAMILY",
    "MODEL_SPEC_VERSION",
    "OUTPUT_SCHEMA_VERSION",
    "ModelConfig",
    "PopulationState",
    "__version__",
    "load_config",
    "simulate_within_host",
]
