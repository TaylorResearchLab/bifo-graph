# bifo-graph Makefile — common test and maintenance targets.
# All test targets require a Python env with requirements-dev.txt installed.

.PHONY: help test test-slow test-kg test-cov test-all clean

help:
	@echo "bifo-graph make targets:"
	@echo "  test        Run unit + regression tests (fast, default)."
	@echo "  test-slow   Run slow tests only (full-pipeline regressions)."
	@echo "  test-kg     Run tests that require a live KG backend connection."
	@echo "  test-cov    Run fast tests with coverage report."
	@echo "  test-all    Run unit + regression + slow (excludes requires_kg)."
	@echo "  clean       Remove pytest, coverage, and __pycache__ artifacts."

test:
	pytest

# Targets that filter by marker may match zero tests (especially early in
# development). pytest exits 5 on "no tests collected", which Make treats
# as failure. The `|| [ $$? -eq 5 ]` translates exit 5 to success so that
# an empty marker bucket isn't a CI/Make failure; any other non-zero exit
# (real test failure, collection error, etc.) still propagates.

test-slow:
	pytest -m slow || [ $$? -eq 5 ]

test-kg:
	pytest -m requires_kg || [ $$? -eq 5 ]

test-cov:
	pytest --cov=pipeline --cov-report=term-missing

test-all:
	pytest -m "not requires_kg" || [ $$? -eq 5 ]

clean:
	rm -rf .pytest_cache .coverage htmlcov
	find . -type d -name __pycache__ -prune -exec rm -rf {} +
	find . -type f -name "*.pyc" -delete
