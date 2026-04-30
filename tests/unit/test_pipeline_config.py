"""
Unit tests for pipeline/pipeline_config.py.

Test categories:

  TestSchema           — Pydantic schema/loader behavior, no cross-validation.
                         Uses pipeline_config fixture files.

  TestEnabledHelpers   — enabled_sources() and enabled_pathway_sabs() methods.

  TestCrossValidation  — Cross-validation against bifo_mapping fixtures.
                         Each invalid bifo_mapping fixture should produce a
                         specific, informative error.

  TestProductionConfig — Loads the canonical config/pipeline_config.yaml
                         and config/bifo_mapping_ddkg.yaml together. Serves
                         as a guard against config drift: if the production
                         pair becomes inconsistent, this test fails.
"""
from __future__ import annotations

import importlib.util
from pathlib import Path

import pytest
from pydantic import ValidationError

# ----------------------------------------------------------------------
# Module loading
# ----------------------------------------------------------------------

_THIS_DIR = Path(__file__).resolve().parent
_PROJECT_ROOT = _THIS_DIR.parent.parent
_PIPELINE_CONFIG_MODULE = _PROJECT_ROOT / "pipeline" / "pipeline_config.py"
_FIXTURES_DIR = _PROJECT_ROOT / "tests" / "regression" / "fixtures" / "configs"
_PRODUCTION_PIPELINE_CONFIG = _PROJECT_ROOT / "config" / "pipeline_config.yaml"
_PRODUCTION_BIFO_MAPPING = _PROJECT_ROOT / "config" / "bifo_mapping_ddkg.yaml"

_spec = importlib.util.spec_from_file_location(
    "pipeline_config", _PIPELINE_CONFIG_MODULE
)
pipeline_config = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(pipeline_config)


# ----------------------------------------------------------------------
# Helpers
# ----------------------------------------------------------------------

def fixture_path(name: str) -> Path:
    """Return the absolute path to a fixture YAML."""
    p = _FIXTURES_DIR / name
    assert p.exists(), f"fixture not found: {p}"
    return p


# ======================================================================
# TestSchema — loader behavior on pipeline_config.yaml fixtures
# ======================================================================

class TestSchema:
    """Schema-level loader behavior: parsing, defaults, validation."""

    def test_minimal_config_loads(self):
        """Happy-path: minimal config loads without cross-validation."""
        config = pipeline_config.load_config(
            fixture_path("pipeline_config_minimal.yaml")
        )
        assert config.schema_version == 1
        assert "MSIGDB" in config.pathway_sources
        assert "GO" in config.pathway_sources
        assert config.gene_sabs == ["HGNC"]
        assert config.burden_control_exclusions == []

    def test_minimal_config_has_correct_pathway_sources(self):
        """Verify the loaded pathway sources match fixture content."""
        config = pipeline_config.load_config(
            fixture_path("pipeline_config_minimal.yaml")
        )
        msigdb = config.pathway_sources["MSIGDB"]
        assert msigdb.sab == "MSIGDB"
        assert msigdb.enabled is True
        assert msigdb.propagating_predicates_g_to_pw == [
            "inverse_pathway_associated_with_gene"
        ]
        assert msigdb.kg_name_prefixes == ["HALLMARK_"]

        go = config.pathway_sources["GO"]
        assert go.sab == "GO"
        assert go.enabled is False
        assert go.kg_name_prefixes == []

    def test_extra_top_level_field_rejected(self):
        """Strict schema rejects unknown top-level fields."""
        with pytest.raises(ValueError) as exc_info:
            pipeline_config.load_config(
                fixture_path("pipeline_config_invalid_extra_field.yaml")
            )
        # Wrapped in our outer "failed schema validation" error message;
        # the underlying Pydantic error names the offending field.
        assert "unknown_top_level_field" in str(exc_info.value)
        # File path included for context
        assert "pipeline_config_invalid_extra_field.yaml" in str(exc_info.value)

    def test_enabled_source_with_no_predicates_rejected(self):
        """Enabled source with empty propagating_predicates_g_to_pw fails."""
        with pytest.raises(ValueError) as exc_info:
            pipeline_config.load_config(
                fixture_path("pipeline_config_invalid_enabled_no_predicates.yaml")
            )
        msg = str(exc_info.value)
        assert "MSIGDB" in msg
        assert "predicates" in msg.lower()

    def test_sab_with_whitespace_rejected(self):
        """SAB names with whitespace are rejected."""
        with pytest.raises(ValueError) as exc_info:
            pipeline_config.load_config(
                fixture_path("pipeline_config_invalid_sab_whitespace.yaml")
            )
        assert "whitespace" in str(exc_info.value).lower()

    def test_unsupported_schema_version_rejected(self):
        """Schema versions other than 1 are rejected."""
        with pytest.raises(ValueError) as exc_info:
            pipeline_config.load_config(
                fixture_path(
                    "pipeline_config_invalid_unsupported_schema_version.yaml"
                )
            )
        msg = str(exc_info.value).lower()
        assert "schema_version" in msg
        # The error suggests an action.
        assert "set" in msg or "update" in msg

    def test_pathway_source_is_frozen(self):
        """PathwaySource instances are immutable (frozen=True)."""
        config = pipeline_config.load_config(
            fixture_path("pipeline_config_minimal.yaml")
        )
        msigdb = config.pathway_sources["MSIGDB"]
        with pytest.raises(ValidationError):
            msigdb.enabled = False

    def test_pipeline_config_is_frozen(self):
        """PipelineConfig instances are immutable (frozen=True)."""
        config = pipeline_config.load_config(
            fixture_path("pipeline_config_minimal.yaml")
        )
        with pytest.raises(ValidationError):
            config.schema_version = 2


# ======================================================================
# TestEnabledHelpers — accessor methods
# ======================================================================

class TestEnabledHelpers:
    """The enabled_sources() and enabled_pathway_sabs() helper methods."""

    def test_enabled_sources_filters_disabled(self):
        """enabled_sources() returns only sources with enabled=true."""
        config = pipeline_config.load_config(
            fixture_path("pipeline_config_minimal.yaml")
        )
        enabled = config.enabled_sources()
        assert len(enabled) == 1
        assert enabled[0].sab == "MSIGDB"

    def test_enabled_pathway_sabs_returns_sorted(self):
        """enabled_pathway_sabs() returns SABs sorted alphabetically."""
        config = pipeline_config.load_config(
            fixture_path("pipeline_config_minimal.yaml")
        )
        sabs = config.enabled_pathway_sabs()
        assert sabs == ["MSIGDB"]
        # Sortedness invariant: the returned list equals its sorted self.
        assert sabs == sorted(sabs)

    def test_enabled_sources_empty_when_all_disabled(self, tmp_path):
        """If no sources are enabled, enabled_sources() returns empty list."""
        # Build a config inline using tmp_path (one of the few cases where
        # inline YAML is justified — testing a degenerate state without
        # cluttering the fixtures dir with a single-use file).
        yaml_content = """
schema_version: 1
pathway_sources:
  MSIGDB:
    enabled: false
    propagating_predicates_g_to_pw:
      - inverse_pathway_associated_with_gene
gene_sabs:
  - HGNC
burden_control_exclusions: []
"""
        cfg_path = tmp_path / "all_disabled.yaml"
        cfg_path.write_text(yaml_content)
        config = pipeline_config.load_config(cfg_path)
        assert config.enabled_sources() == []
        assert config.enabled_pathway_sabs() == []


# ======================================================================
# TestCrossValidation — checks against bifo_mapping fixtures
# ======================================================================

class TestCrossValidation:
    """Cross-validation between pipeline_config and bifo_mapping fixtures."""

    def test_minimal_consistent_pair_loads(self):
        """Cross-validation succeeds for the minimal consistent pair."""
        config = pipeline_config.load_config(
            fixture_path("pipeline_config_minimal.yaml"),
            bifo_mapping_path=fixture_path("bifo_mapping_minimal.yaml"),
        )
        # Loading didn't raise; that's the point.
        assert config.enabled_pathway_sabs() == ["MSIGDB"]

    def test_predicate_missing_from_bifo_mapping_fails(self):
        """Cross-validation reports predicates missing from bifo_mapping."""
        with pytest.raises(ValueError) as exc_info:
            pipeline_config.load_config(
                fixture_path("pipeline_config_minimal.yaml"),
                bifo_mapping_path=fixture_path(
                    "bifo_mapping_predicate_missing.yaml"
                ),
            )
        msg = str(exc_info.value)
        assert "inverse_pathway_associated_with_gene" in msg
        assert "not" in msg  # "is not defined in bifo_mapping"

    def test_wrong_direction_in_bifo_mapping_fails(self):
        """Cross-validation flags direction mismatches."""
        with pytest.raises(ValueError) as exc_info:
            pipeline_config.load_config(
                fixture_path("pipeline_config_minimal.yaml"),
                bifo_mapping_path=fixture_path(
                    "bifo_mapping_wrong_direction.yaml"
                ),
            )
        msg = str(exc_info.value)
        assert "target_to_source" in msg
        assert "source_to_target" in msg

    def test_nonpropagating_classification_fails(self):
        """Cross-validation flags predicates classified as nonpropagating."""
        with pytest.raises(ValueError) as exc_info:
            pipeline_config.load_config(
                fixture_path("pipeline_config_minimal.yaml"),
                bifo_mapping_path=fixture_path(
                    "bifo_mapping_nonpropagating_classification.yaml"
                ),
            )
        msg = str(exc_info.value)
        assert "nonpropagating_context" in msg

    def test_unknown_classification_fails_with_helpful_error(self):
        """Cross-validation flags classifications outside the known taxonomy."""
        with pytest.raises(ValueError) as exc_info:
            pipeline_config.load_config(
                fixture_path("pipeline_config_minimal.yaml"),
                bifo_mapping_path=fixture_path(
                    "bifo_mapping_unknown_classification.yaml"
                ),
            )
        msg = str(exc_info.value)
        assert "invented_unknown_class" in msg
        # The error should list the known classifications so the user
        # can fix it without consulting the source.
        assert "mechanistic" in msg

    def test_multiple_errors_collected_in_one_raise(self):
        """All errors reported in one exception, not one-at-a-time."""
        with pytest.raises(ValueError) as exc_info:
            pipeline_config.load_config(
                fixture_path("pipeline_config_two_predicates.yaml"),
                bifo_mapping_path=fixture_path(
                    "bifo_mapping_two_errors.yaml"
                ),
            )
        msg = str(exc_info.value)
        # Both error conditions should appear in a single error message.
        assert "inverse_pathway_associated_with_gene" in msg
        assert "inverse_has_signature_gene" in msg
        assert "nonpropagating_context" in msg
        assert "target_to_source" in msg

    def test_disabled_sources_skip_cross_validation(self):
        """Disabled sources don't trigger cross-validation errors.

        Even if a disabled source's predicates would fail validation,
        the loader doesn't check them. Only enabled sources are validated.
        """
        # The minimal pipeline_config has GO disabled with the predicate
        # `gene_plays_role_in_process`. The bifo_mapping_predicate_missing
        # fixture doesn't define that predicate either. But since GO is
        # disabled, that's not an error.
        # However, MSIGDB's `inverse_pathway_associated_with_gene` IS
        # missing from bifo_mapping_predicate_missing, so the load DOES fail
        # — but on the MSIGDB predicate, not the GO one.
        with pytest.raises(ValueError) as exc_info:
            pipeline_config.load_config(
                fixture_path("pipeline_config_minimal.yaml"),
                bifo_mapping_path=fixture_path(
                    "bifo_mapping_predicate_missing.yaml"
                ),
            )
        msg = str(exc_info.value)
        # The error mentions MSIGDB's predicate but should NOT mention
        # GO or its predicates.
        assert "inverse_pathway_associated_with_gene" in msg
        assert "gene_plays_role_in_process" not in msg

    def test_skip_cross_validation_when_path_is_none(self):
        """Passing bifo_mapping_path=None skips cross-validation entirely."""
        # Use a config that would FAIL cross-validation against ANY real
        # bifo_mapping (or use minimal). The point is that without a
        # bifo_mapping path, no cross-validation happens.
        config = pipeline_config.load_config(
            fixture_path("pipeline_config_minimal.yaml"),
            bifo_mapping_path=None,
        )
        assert config.enabled_pathway_sabs() == ["MSIGDB"]


# ======================================================================
# TestStructuralErrors — bifo_mapping structural-error reporting
# ======================================================================

class TestStructuralErrors:
    """The loader should produce informative errors when bifo_mapping
    is structurally malformed, not silently cascade into misleading
    'predicate not found' errors."""

    def test_bifo_mapping_missing_edge_resolution(self):
        """bifo_mapping without an edge_resolution section: clear error."""
        with pytest.raises(ValueError) as exc_info:
            pipeline_config.load_config(
                fixture_path("pipeline_config_minimal.yaml"),
                bifo_mapping_path=fixture_path(
                    "bifo_mapping_missing_edge_resolution.yaml"
                ),
            )
        msg = str(exc_info.value)
        assert "edge_resolution" in msg
        assert "missing" in msg.lower()
        # Path is included so the user knows which file to fix.
        assert "bifo_mapping_missing_edge_resolution.yaml" in msg

    def test_bifo_mapping_missing_predicate_to_flow(self):
        """bifo_mapping missing predicate_to_flow: clear error pointing
        at the right subsection."""
        with pytest.raises(ValueError) as exc_info:
            pipeline_config.load_config(
                fixture_path("pipeline_config_minimal.yaml"),
                bifo_mapping_path=fixture_path(
                    "bifo_mapping_missing_predicate_to_flow.yaml"
                ),
            )
        msg = str(exc_info.value)
        assert "predicate_to_flow" in msg

    def test_bifo_mapping_malformed_predicate(self):
        """A predicate entry missing required fields: the error names
        the offending predicate so the user can find it."""
        with pytest.raises(ValueError) as exc_info:
            pipeline_config.load_config(
                fixture_path("pipeline_config_minimal.yaml"),
                bifo_mapping_path=fixture_path(
                    "bifo_mapping_malformed_predicate.yaml"
                ),
            )
        msg = str(exc_info.value)
        # Predicate name in error
        assert "inverse_pathway_associated_with_gene" in msg
        # File path in error
        assert "bifo_mapping_malformed_predicate.yaml" in msg

    def test_bifo_mapping_not_a_mapping(self):
        """A YAML file that's a list at top level: clear 'not a mapping'
        error."""
        with pytest.raises(ValueError) as exc_info:
            pipeline_config.load_config(
                fixture_path("pipeline_config_minimal.yaml"),
                bifo_mapping_path=fixture_path(
                    "bifo_mapping_not_a_mapping.yaml"
                ),
            )
        msg = str(exc_info.value)
        assert "mapping" in msg.lower()

    def test_pipeline_config_invalid_yaml_includes_path(self, tmp_path):
        """When pipeline_config.yaml fails YAML parse, the file path is
        in the error."""
        bad = tmp_path / "broken.yaml"
        # Actually-malformed YAML (unclosed brace).
        bad.write_text("foo: {bar: [1, 2,\n")
        with pytest.raises(ValueError) as exc_info:
            pipeline_config.load_config(bad)
        msg = str(exc_info.value)
        assert "broken.yaml" in msg
        assert "yaml" in msg.lower()

    def test_pipeline_config_path_in_schema_error(self):
        """When pipeline_config.yaml has a schema violation, the file
        path is in the wrapped error message."""
        with pytest.raises(ValueError) as exc_info:
            pipeline_config.load_config(
                fixture_path("pipeline_config_invalid_extra_field.yaml")
            )
        msg = str(exc_info.value)
        # Confirms the file-path-context wrapping is in place.
        assert "pipeline_config_invalid_extra_field.yaml" in msg


# ======================================================================
# TestProductionConfig — guards against config drift
# ======================================================================

class TestProductionConfig:
    """Tests against the real config/pipeline_config.yaml and bifo_mapping.

    These tests serve a different purpose from the loader-mechanic tests
    above: they assert the canonical configuration pair is internally
    consistent. If someone edits one file and breaks consistency with the
    other, this catches it.
    """

    def test_production_pipeline_config_loads_with_cross_validation(self):
        """The shipped pipeline_config.yaml is internally consistent with
        bifo_mapping_ddkg.yaml."""
        # If this test fails, either:
        #   (a) someone edited bifo_mapping_ddkg.yaml and changed a
        #       predicate's classification or direction; or
        #   (b) someone added a predicate to pipeline_config.yaml that
        #       doesn't exist in bifo_mapping_ddkg.yaml.
        config = pipeline_config.load_config(
            _PRODUCTION_PIPELINE_CONFIG,
            bifo_mapping_path=_PRODUCTION_BIFO_MAPPING,
        )
        assert config.schema_version == 1

    def test_production_config_has_msigdb_enabled(self):
        """MSIGDB is the canonical enabled source."""
        config = pipeline_config.load_config(_PRODUCTION_PIPELINE_CONFIG)
        assert "MSIGDB" in config.pathway_sources
        assert config.pathway_sources["MSIGDB"].enabled is True

    def test_production_config_has_go_disabled(self):
        """GO is declared but disabled by default (per Step 1 design)."""
        config = pipeline_config.load_config(_PRODUCTION_PIPELINE_CONFIG)
        assert "GO" in config.pathway_sources
        assert config.pathway_sources["GO"].enabled is False

    def test_production_config_has_hgnc_in_gene_sabs(self):
        """HGNC is the canonical gene SAB."""
        config = pipeline_config.load_config(_PRODUCTION_PIPELINE_CONFIG)
        assert "HGNC" in config.gene_sabs

    def test_production_config_burden_exclusions_empty_by_default(self):
        """The shipped config has no burden exclusions applied."""
        config = pipeline_config.load_config(_PRODUCTION_PIPELINE_CONFIG)
        assert config.burden_control_exclusions == []
