"""
Tests for pipeline/check_configs.py.

Covers:
  - Success path produces expected output and exit 0.
  - Failure paths produce error output and exit 1.
  - Missing files exit 2.
  - --quiet suppresses success output but not errors.
  - "Did you mean" hints fire correctly for typos and forward-direction
    predicates.
"""
from __future__ import annotations

import importlib.util
import io
from pathlib import Path

import pytest

# ----------------------------------------------------------------------
# Module loading
# ----------------------------------------------------------------------

_THIS_DIR = Path(__file__).resolve().parent
_PROJECT_ROOT = _THIS_DIR.parent.parent
_CHECK_CONFIGS = _PROJECT_ROOT / "pipeline" / "check_configs.py"
_FIXTURES_DIR = _PROJECT_ROOT / "tests" / "regression" / "fixtures" / "configs"

_spec = importlib.util.spec_from_file_location("check_configs", _CHECK_CONFIGS)
check_configs = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(check_configs)


def fixture_path(name: str) -> Path:
    p = _FIXTURES_DIR / name
    assert p.exists(), f"fixture not found: {p}"
    return p


def run_check(
    pipeline_config_name: str = "pipeline_config_minimal.yaml",
    bifo_mapping_name: str = "bifo_mapping_minimal.yaml",
    quiet: bool = False,
) -> tuple[int, str, str]:
    """Run check_configs.check() against the named fixtures.

    Returns (exit_code, stdout_text, stderr_text).
    """
    out = io.StringIO()
    err = io.StringIO()
    rc = check_configs.check(
        pipeline_config_path=fixture_path(pipeline_config_name),
        bifo_mapping_path=fixture_path(bifo_mapping_name),
        quiet=quiet,
        out=out,
        err=err,
    )
    return rc, out.getvalue(), err.getvalue()


# ======================================================================
# Success-path tests
# ======================================================================

class TestSuccess:
    def test_minimal_consistent_pair_exits_zero(self):
        rc, out, err = run_check()
        assert rc == check_configs.EXIT_OK

    def test_success_output_has_ok_markers(self):
        rc, out, err = run_check()
        assert "[OK]" in out
        assert "[FAIL]" not in out

    def test_success_output_mentions_msigdb(self):
        rc, out, err = run_check()
        assert "MSIGDB" in out

    def test_success_output_lists_predicates(self):
        rc, out, err = run_check()
        assert "inverse_pathway_associated_with_gene" in out

    def test_success_output_includes_classifications(self):
        rc, out, err = run_check()
        assert "weak_mechanistic_or_observational" in out

    def test_success_output_concludes_with_ready_message(self):
        rc, out, err = run_check()
        assert "ready" in out.lower()

    def test_quiet_mode_suppresses_success_output(self):
        rc, out, err = run_check(quiet=True)
        assert rc == check_configs.EXIT_OK
        assert out == ""
        assert err == ""


# ======================================================================
# Failure-path tests
# ======================================================================

class TestFailure:
    def test_invalid_pipeline_config_exits_one(self):
        rc, out, err = run_check(
            pipeline_config_name="pipeline_config_invalid_extra_field.yaml"
        )
        assert rc == check_configs.EXIT_CONFIG_ERROR

    def test_invalid_pipeline_config_error_includes_path(self):
        rc, out, err = run_check(
            pipeline_config_name="pipeline_config_invalid_extra_field.yaml"
        )
        # File-path context wrapping is in place
        assert "pipeline_config_invalid_extra_field.yaml" in err

    def test_invalid_pipeline_config_error_marked_fail(self):
        rc, out, err = run_check(
            pipeline_config_name="pipeline_config_invalid_extra_field.yaml"
        )
        assert "[FAIL]" in err

    def test_predicate_missing_from_bifo_mapping_exits_one(self):
        rc, out, err = run_check(
            bifo_mapping_name="bifo_mapping_predicate_missing.yaml"
        )
        assert rc == check_configs.EXIT_CONFIG_ERROR

    def test_predicate_missing_error_names_predicate(self):
        rc, out, err = run_check(
            bifo_mapping_name="bifo_mapping_predicate_missing.yaml"
        )
        assert "inverse_pathway_associated_with_gene" in err

    def test_predicate_missing_error_grouped_by_source(self):
        rc, out, err = run_check(
            bifo_mapping_name="bifo_mapping_predicate_missing.yaml"
        )
        # The structured error display lists "source MSIGDB:"
        assert "source MSIGDB" in err

    def test_two_errors_both_reported(self):
        rc, out, err = run_check(
            pipeline_config_name="pipeline_config_two_predicates.yaml",
            bifo_mapping_name="bifo_mapping_two_errors.yaml",
        )
        assert rc == check_configs.EXIT_CONFIG_ERROR
        # Both errors mentioned in the same output
        assert "inverse_pathway_associated_with_gene" in err
        assert "inverse_has_signature_gene" in err

    def test_two_errors_numbered_separately(self):
        """Errors are numbered [1], [2] for the user to count."""
        rc, out, err = run_check(
            pipeline_config_name="pipeline_config_two_predicates.yaml",
            bifo_mapping_name="bifo_mapping_two_errors.yaml",
        )
        assert "[1]" in err
        assert "[2]" in err

    def test_failure_quiet_still_emits_error(self):
        """--quiet suppresses success output but not failure output."""
        rc, out, err = run_check(
            pipeline_config_name="pipeline_config_invalid_extra_field.yaml",
            quiet=True,
        )
        assert rc == check_configs.EXIT_CONFIG_ERROR
        # No success/header output
        assert out == ""
        # But the failure detail is still in stderr
        assert "[FAIL]" in err


# ======================================================================
# Setup-error tests (file not found)
# ======================================================================

class TestSetupErrors:
    def test_missing_pipeline_config_exits_two(self, tmp_path):
        bogus = tmp_path / "does_not_exist.yaml"
        out = io.StringIO()
        err = io.StringIO()
        rc = check_configs.check(
            pipeline_config_path=bogus,
            bifo_mapping_path=fixture_path("bifo_mapping_minimal.yaml"),
            quiet=False,
            out=out,
            err=err,
        )
        assert rc == check_configs.EXIT_SETUP_ERROR
        assert "not found" in err.getvalue().lower()

    def test_missing_bifo_mapping_exits_two(self, tmp_path):
        bogus = tmp_path / "does_not_exist.yaml"
        out = io.StringIO()
        err = io.StringIO()
        rc = check_configs.check(
            pipeline_config_path=fixture_path("pipeline_config_minimal.yaml"),
            bifo_mapping_path=bogus,
            quiet=False,
            out=out,
            err=err,
        )
        assert rc == check_configs.EXIT_SETUP_ERROR
        assert "not found" in err.getvalue().lower()


# ======================================================================
# "Did you mean" hint tests
# ======================================================================

class TestDidYouMean:
    """The pre-flight tool offers helpful suggestions for common errors."""

    def test_typo_suggestion_uses_difflib(self, tmp_path):
        """A typo'd predicate name produces a suggestion of the closest
        match in bifo_mapping."""
        # Construct a pipeline_config with a typo'd predicate
        bad_config = tmp_path / "typo.yaml"
        bad_config.write_text("""
schema_version: 1
pathway_sources:
  MSIGDB:
    enabled: true
    propagating_predicates_g_to_pw:
      - invrese_pathway_associated_with_gene
gene_sabs:
  - HGNC
burden_control_exclusions: []
""")
        out = io.StringIO()
        err = io.StringIO()
        rc = check_configs.check(
            pipeline_config_path=bad_config,
            bifo_mapping_path=fixture_path("bifo_mapping_minimal.yaml"),
            quiet=False,
            out=out,
            err=err,
        )
        assert rc == check_configs.EXIT_CONFIG_ERROR
        # The suggestion should be the correct spelling.
        err_text = err.getvalue()
        assert "Did you mean" in err_text
        assert "inverse_pathway_associated_with_gene" in err_text

    def test_forward_direction_suggests_inverse(self, tmp_path):
        """If a forward-direction predicate is listed under
        propagating_predicates_g_to_pw, suggest the inverse_ form
        if it exists in bifo_mapping."""
        # bifo_mapping_minimal has inverse_pathway_associated_with_gene
        # defined. We reference its forward counterpart.
        # First we need a bifo_mapping that ALSO defines the forward
        # predicate (so we can test the wrong-direction error path,
        # not the missing-predicate path).
        custom_mapping = tmp_path / "mapping_with_forward.yaml"
        custom_mapping.write_text("""
version: 0.0.0-fixture
edge_resolution:
  predicate_to_flow:
    pathway_associated_with_gene:
      direction: target_to_source
      classification: nonpropagating_context
    inverse_pathway_associated_with_gene:
      direction: source_to_target
      classification: weak_mechanistic_or_observational
""")
        bad_config = tmp_path / "config_with_forward.yaml"
        bad_config.write_text("""
schema_version: 1
pathway_sources:
  MSIGDB:
    enabled: true
    propagating_predicates_g_to_pw:
      - pathway_associated_with_gene
gene_sabs:
  - HGNC
burden_control_exclusions: []
""")
        out = io.StringIO()
        err = io.StringIO()
        rc = check_configs.check(
            pipeline_config_path=bad_config,
            bifo_mapping_path=custom_mapping,
            quiet=False,
            out=out,
            err=err,
        )
        assert rc == check_configs.EXIT_CONFIG_ERROR
        # The hint should suggest the inverse_ form.
        err_text = err.getvalue()
        assert "Did you mean" in err_text
        assert "inverse_pathway_associated_with_gene" in err_text
        # And explain the convention.
        assert "inverse" in err_text.lower()


# ======================================================================
# CLI tests (via main())
# ======================================================================

class TestMain:
    def test_main_with_explicit_paths(self, monkeypatch, capsys):
        """main() accepts CLI arguments and returns the right exit code."""
        rc = check_configs.main([
            "--pipeline-config", str(fixture_path("pipeline_config_minimal.yaml")),
            "--bifo-mapping", str(fixture_path("bifo_mapping_minimal.yaml")),
        ])
        assert rc == check_configs.EXIT_OK

    def test_main_quiet_flag(self, monkeypatch, capsys):
        """--quiet flag suppresses success output via the CLI."""
        rc = check_configs.main([
            "--pipeline-config", str(fixture_path("pipeline_config_minimal.yaml")),
            "--bifo-mapping", str(fixture_path("bifo_mapping_minimal.yaml")),
            "--quiet",
        ])
        assert rc == check_configs.EXIT_OK
        captured = capsys.readouterr()
        assert captured.out == ""
