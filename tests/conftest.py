"""
Shared pytest configuration and path constants for bifo-graph tests.

Per-subdir conftest.py files (if any are added later) layer on top of this one.
"""
from pathlib import Path

# The repository root, found by walking up from this file.
PROJECT_ROOT: Path = Path(__file__).resolve().parent.parent

# Test directory roots.
TESTS_ROOT: Path = PROJECT_ROOT / "tests"
UNIT_TESTS_DIR: Path = TESTS_ROOT / "unit"
REGRESSION_TESTS_DIR: Path = TESTS_ROOT / "regression"
INTEGRATION_TESTS_DIR: Path = TESTS_ROOT / "integration"

# Paths to checked-in fixture data and run manifests.
# FIXTURES_DIR is created in subsequent steps when first needed.
# RUN_MANIFESTS_DIR holds JSONL manifests describing known-good HPC runs;
# these are byte-level diff targets for regression tests.
FIXTURES_DIR: Path = REGRESSION_TESTS_DIR / "fixtures"
RUN_MANIFESTS_DIR: Path = REGRESSION_TESTS_DIR / "run_manifests"

# Pipeline source directory, for tests that import pipeline modules.
PIPELINE_DIR: Path = PROJECT_ROOT / "pipeline"
