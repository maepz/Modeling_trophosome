# Release checklist for 0.5.0

## Completed locally

- [x] Model family and specification version declared.
- [x] Output schema version declared.
- [x] Package version has one authoritative source in `pyproject.toml`.
- [x] Full automated test suite passes.
- [x] Distributional validation passes at 30,000 draws per analytical check.
- [x] Independent cell-level mutation-timing reference passes.
- [x] Lint gate passes for maintained Python code.
- [x] Maintained package coverage exceeds 85%.
- [x] Reproducible validation report and machine-readable results generated.
- [x] Biologist-facing validation narrative and figures generated.
- [x] Changelog, citation metadata and model specification prepared.
- [x] Wheel and source distribution build successfully.
- [x] Clean wheel import, version metadata and CLI validation succeed.
- [x] Build checks included in CI.
- [x] Generation-boundary checkpoint restart and output truncation validated.

## Required before public release

- [ ] Choose and add a software license. No license has been selected on the
  researcher's behalf.
- [ ] Review the author and citation metadata.
- [ ] Separate unrelated pre-existing working-tree changes and generated OS
  metadata from the release commit.
- [ ] Review the complete proposed diff.
- [ ] Commit through a dedicated branch and reviewed pull request.
- [ ] Tag the accepted commit `v0.5.0` and create the GitHub release.

Commit, tag, push and GitHub-release operations were intentionally not performed
during local release preparation.
