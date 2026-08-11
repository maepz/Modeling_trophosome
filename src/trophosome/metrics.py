"""Population summaries with defined edge-case behavior."""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np

from trophosome.count_model import PopulationState


@dataclass(frozen=True)
class PopulationMetrics:
    population_size: int
    richness: int
    mean_within_host_fitness: float
    mean_free_living_fitness: float
    gene_diversity: float
    shannon_diversity: float
    pielou_evenness: float


def population_metrics(state: PopulationState) -> PopulationMetrics:
    frequencies = state.counts / state.size
    shannon = float(-np.sum(frequencies * np.log(frequencies)))
    gene_diversity = float(1.0 - np.sum(frequencies**2))
    if state.size > 1:
        gene_diversity *= state.size / (state.size - 1)
    evenness = shannon / math.log(state.richness) if state.richness > 1 else 1.0
    return PopulationMetrics(
        population_size=state.size,
        richness=state.richness,
        mean_within_host_fitness=float(
            np.average(state.within_host_fitness, weights=state.counts)
        ),
        mean_free_living_fitness=float(
            np.average(state.free_living_fitness, weights=state.counts)
        ),
        gene_diversity=gene_diversity,
        shannon_diversity=shannon,
        pielou_evenness=evenness,
    )
