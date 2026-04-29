# bifo-graph tests

Test suite for `bifo-graph`. Tests are organized by scope, with markers
controlling which run by default.

## Setup

```bash
pip install -r requirements.txt -r requirements-dev.txt
```

## Running tests

From the repository root:

```bash
make test         # Fast tests only (unit + regression). Default.
make test-slow    # Slow tests (full-pipeline regressions).
make test-kg      # Tests that require a live KG backend connection.
make test-cov     # Fast tests with coverage report.
make test-all     # Unit + regression + slow (excludes requires_kg).
```

Or invoke pytest directly:

```bash
pytest                       # Same as `make test`.
pytest -m slow               # Slow tests only.
pytest -m requires_kg        # KG-dependent tests only.
pytest -m ""                 # Everything (no marker filter).
pytest tests/unit/           # All tests in tests/unit/.
pytest -k test_config        # Tests whose name contains "test_config".
```

## Layout

| Directory                    | Purpose |
|------------------------------|---------|
| `tests/unit/`                | Pure-function tests. No I/O, no subprocess, no fixtures larger than a few KB. Should run in <1s each. |
| `tests/regression/`          | Tests that run pipeline modules against checked-in fixture inputs and compare outputs to checked-in expected values. |
| `tests/regression/fixtures/` | Small (<=100 KB) inputs that exercise specific code paths. |
| `tests/regression/baselines/`| Manifests describing known-good HPC runs. JSONL format; one record per output file. Validates that downstream code changes don't silently alter pipeline output. |
| `tests/integration/`         | Tests requiring external infrastructure (live KG, full HPC compute). Marked `slow` and/or `requires_kg`. |

## Markers

| Marker           | Meaning                                               | Run by default? |
|------------------|-------------------------------------------------------|-----------------|
| (unmarked)       | Fast, isolated, runs anywhere.                        | Yes             |
| `slow`           | Wall time >5s (typically full-pipeline runs).         | No              |
| `requires_kg`    | Needs a live knowledge-graph backend connection.      | No              |

Apply with `@pytest.mark.slow` or `@pytest.mark.requires_kg` on test functions
or classes. The `--strict-markers` setting in `pyproject.toml` makes typo'd
marker names a hard error.

## Adding tests

- A new pure-function test: `tests/unit/test_<module>.py`.
- A new pipeline-output regression test: `tests/regression/test_<stage>_regression.py`.
- A new full-cohort baseline diff: `tests/integration/test_<cohort>_<scenario>.py`,
  marked `@pytest.mark.slow`.
- New fixtures: `tests/regression/fixtures/`. Keep them small; large reference
  data lives on the HPC and is referenced via baseline manifests in
  `tests/regression/baselines/`.

When adding a regression test, also commit the expected output files alongside
the test. The validation gate at every step in development is "diff the new
output against the committed expected output."
