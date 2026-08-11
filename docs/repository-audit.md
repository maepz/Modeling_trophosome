# Repository audit (2026-08-10)

The checkout was approximately 402 MB. Git reported about 335 MB of loose objects.
Tracked generated material included 30 notebook-checkpoint files, 5 `.DS_Store`
files, 17 Python bytecode files, 4 egg-info files, 12 scheduler log/accounting
files, and at least 16 generated time-series example outputs.

Changes made in this reorganization:

- added portable `pyproject.toml` dependency groups and a command-line entry point;
- added a domain-specific package name (`trophosome`) while retaining the legacy
  package for notebook compatibility;
- added validated TOML configurations and recorded provenance;
- added automated tests and continuous integration;
- added model, scaling, reproducibility, and development documentation;
- added ignores for generated outputs, caches, notebooks checkpoints, and logs.

After the repository was updated, V1.4, V2.2 and V3.1 were re-inspected. The
optimized `trophosome/count_model.py` is now the maintained main prototype;
V3.1 is a legacy neutral endpoint comparator. A V3.1 workload estimator and
tests remain available without modifying the legacy V3.1 scientific
implementation.

Deliberately not performed:

- deleting or rewriting existing notebooks and research outputs;
- rewriting Git history to purge large generated objects;
- assigning a software license;
- silently changing legacy model behavior.

Those actions can cause data loss, disrupt collaborators, or require a scientific
decision, so they should be handled in a separate reviewed cleanup.
