# Changelog

All notable changes to the maintained model are documented here. Software
versions follow Semantic Versioning. The scientific model specification and
output schema have independent versions recorded in every simulation.

## [Unreleased]

### Added

- The frozen Stage 3 part-one matrix: 24 new cells, twelve matched seed blocks,
  288 new populations, 100 passages, and a reused twelve-population no-return
  control (25 primary conditions, 300 primary populations). Crosses H=100,
  1,000, 10,000 with alpha targets 0.001, 0.01, 0.1, 0.99 and u=0/10^-10.
  Exact total return is matched across H, with target/realized alpha recorded.
  Includes three gated HPC safety runs, checkpoint/resume launching, frozen
  passage-100 Stage 2 references, seed-paired H-by-feedback TV interactions,
  Shannon/Simpson indices, late-time drift/precision checks, and
  automatic or stand-alone PDF/Markdown reporting. No scientific transition
  rules or model/output version numbers changed.

- An optional fixed, non-depleting regional source that exchanges an exact
  configured fraction of the focal environmental population after host return.
- Aligned focal and regional root-strain vectors, including focal-only,
  regional-only and shared initial strains.
- Per-generation realized emigrant and immigrant counts, initial-strain origin
  metadata, and post-return and post-migration richness summaries.
- Distributional validation of without-replacement focal emigration and
  with-replacement immigration from the fixed regional composition.
- A frozen six-sentinel, 72-population Phase 1 second-pilot workflow with
  checkpoint/resume support, 250-passage environmental trajectories, automatic
  auditing, stationarity diagnostics, precision recommendations, and a
  self-contained biological PDF report.

### Changed

- Environmental stage order is now explicit: neutral regulation after host
  return, fixed-pool exchange, then optional free-living selection.
- Software version is `0.7.0`, model specification is `2.1.0`, and output schema
  is `2.3.0`.
- Existing configuration files remain migration-free by default.

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
