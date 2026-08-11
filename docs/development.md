# Development guide

## Setup

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -e '.[dev,coalescent,plot]'
```

Run checks before committing:

```bash
ruff check src tests scripts
python -m pytest --cov=trophosome --cov-report=term-missing --cov-fail-under=85
python scripts/validate_distributions.py --repetitions 30000 \
  --json docs/distributional-validation.json \
  --markdown docs/distributional-validation-results.md
python scripts/plot_distributional_validation.py --repetitions 30000
python -m build
```

The automated test suite runs a smaller release-gating distributional check.
The 30,000-draw command regenerates the scientific validation evidence used for
a release candidate.

The tests also use the standard `unittest` API and can be run without pytest:

```bash
PYTHONPATH=src python -m unittest discover -s tests -p 'test_*.py'
```

## Boundaries

- Treat `trophosome/count_model.py` and the streaming runner as the maintained
  main prototype.
- Treat V1.3/V1.4, V2.2 and V3.1 as legacy comparison implementations.
- Put reusable logic in `src/trophosome/`, not notebooks.
- Keep model assumptions in typed configuration and `docs/model.md`.
- Keep plotting separate from state transitions.
- Pass random generators explicitly.
- Do not keep per-host histories by default.
- Do not introduce a new approximation without documenting its exactness status
  and adding comparison tests against a smaller exact model.

## Main-prototype optimization boundary

Performance changes to the count model must preserve its haploid, discrete-
generation Wright--Fisher transition distribution. Appropriate changes include
array reuse, vectorized draws, bounded task submission, deterministic batching,
streaming outputs, and avoiding unrequested histories. A shortcut that changes
mutation timing, drift, selection, or clone-size distributions must receive a
new model name and exact-versus-approximate validation.

## Legacy comparisons

`src/project_package/` remains importable so existing notebooks are not broken.
Do not silently rewrite V3.1 or relabel it as exact. Add regression and
distributional tests around `run_generation_of_host_pop_v3_1` before changing
its behavior. A future maintained endpoint approximation should have a separate
name and an explicit validation region. Useful legacy work remains:

1. input validation and deterministic seed management;
2. population-size-conserving endpoint abundance allocation;
3. efficient Ewens partition sampling;
4. stable founder occupancy sampling;
5. bounded summaries and optional ancestry output;
6. host-generation/environment update.

At each step, distinguish a bug fix from a change in scientific assumptions. Add
a regression or distributional equivalence test. Delete a legacy function only
after all references have migrated.

Use `scripts/estimate_v3_1_load.py` to screen proposed V3.1 parameters before a
large run. It is a diagnostic calculator, not a simulator.

## Repository hygiene

The audit found tracked notebook checkpoints, OS metadata, bytecode, egg metadata,
scheduler logs, and generated example results. `.gitignore` prevents new copies.
Existing tracked artifacts should be removed from Git history only after verifying
that scientifically important outputs are archived elsewhere. History rewriting
is a separate, coordinated operation because collaborators must re-clone.
