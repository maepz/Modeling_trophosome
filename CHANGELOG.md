# Changelog

All notable changes to the maintained model are documented here. Software
versions follow Semantic Versioning. The scientific model specification and
output schema have independent versions recorded in every simulation.

## [Unreleased]

No changes yet.

## [0.5.0] - 2026-08-11

### Added

- Validated, atomic recovery checkpoints at completed host-generation
  boundaries.
- Wall-clock checkpoint targets such as `"30m"`, `"1h"`, and `"2h"`, with two
  validated recovery points retained by default while a run is active.
- Exact restart through `trophosome run --resume`, including active lineage
  depths, accumulated summaries, source and configuration hashes, deterministic
  strain-ID metadata, and committed output-file offsets.
- A committed `completion.json` record with output sizes and final-environment
  checksums.

### Changed

- The output schema is now `2.1.0`; generation summaries are committed as the
  run proceeds and provenance records resume events.
- Resume safely truncates incomplete rows written after the newest valid
  checkpoint and falls back to the preceding checkpoint if the newest is
  corrupted.
- Recovery checkpoints are deleted after successful output and final-state
  verification. Only `final_environment_repNNN.npz` is retained permanently as
  population state.
- Execution-only settings such as worker count and host batch size may change on
  resume without changing the seeded scientific result.
- The scientific model specification remains `wright_fisher_counts/2.0.0`.

## [0.4.0] - 2026-08-11

### Added

- Independent `within_host_fitness` and `free_living_fitness` traits for every
  strain.
- Independently sampled habitat-specific mutation fitness effects from the same
  configured normal distribution.
- Optional exact free-living Wright--Fisher selection after host return.
- A dual-selection example configuration and distributional checks for
  free-living selection and fitness-effect independence.

### Changed

- The scientific model specification is now `wright_fisher_counts/2.0.0`.
- The output schema is now `2.0.0`; environment, lineage and summary outputs use
  explicit habitat-specific fitness names.
- The neutral Phase 1 behavior remains available with
  `free_living_selection=false`.

## [0.3.0] - 2026-08-10

### Added

- Maintained `wright_fisher_counts` package and command-line interface.
- Exact count-based haploid Wright–Fisher reproduction with optional selection,
  explicit mutation timing and infinite-alleles lineage events.
- Bounded host batching, persistent workers and reproducible logical-host random
  streams.
- Versioned configurations, streaming outputs, checkpoints and run provenance.
- Analytical and cell-level distributional validation harness.
- Biologist-facing validation report with reproducible distribution figures.
- Formal model specification `wright_fisher_counts/1.0.0` and output schema
  `1.0.0`.
- Documented V1.3/V1.4, V2.2 and V3.1 legacy comparison roles.

### Changed

- The exact-count implementation is restored as the main research prototype.
- Environmental capacity regulation randomizes exact Hamilton-apportionment
  ties using a seeded stream so strain identifiers cannot receive a systematic
  rounding advantage.
- Package version is read from installed project metadata rather than duplicated
  in the source package.

### Validated

- Neutral and fitness-weighted Wright–Fisher offspring distributions.
- Bernoulli mutation counts and mutation-parent assignment.
- With- and without-replacement founder/escape sampling.
- Optimized inverse-CDF reservoir sampling.
- Multi-generation neutral drift.
- Mutation timing and jackpot-clone distributions against a cell-level reference.
- Environmental apportionment label neutrality.

[Unreleased]: https://github.com/maepz/Modeling_trophosome/compare/v0.5.0...HEAD
[0.5.0]: https://github.com/maepz/Modeling_trophosome/compare/v0.4.0...v0.5.0
[0.4.0]: https://github.com/maepz/Modeling_trophosome/compare/v0.3.0...v0.4.0
[0.3.0]: https://github.com/maepz/Modeling_trophosome/releases/tag/v0.3.0
