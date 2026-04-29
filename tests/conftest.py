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

# Paths to checked-in fixture data and baseline manifests.
# These directories are created in subsequent steps (Step 0 creates BASELINES_DIR;
# Step 1 creates FIXTURES_DIR/configs/). pathlib.Path objects don't require
# the underlying directory to exist.
FIXTURES_DIR: Path = REGRESSION_TESTS_DIR / "fixtures"
BASELINES_DIR: Path = REGRESSION_TESTS_DIR / "baselines"

# Pipeline source directory, for tests that import pipeline modules.
PIPELINE_DIR: Path = PROJECT_ROOT / "pipeline"
