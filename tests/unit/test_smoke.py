"""
Smoke test confirming the pytest scaffold runs.

Replace or remove once real unit tests exist. Until then, this prevents
pytest from exiting with status 5 ("no tests collected"), which Make would
treat as a failure.
"""


def test_smoke():
    """Trivial assertion; confirms pytest, conftest, and markers wire up."""
    assert True
