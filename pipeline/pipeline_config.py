"""
pipeline/pipeline_config.py

Loader and types for config/pipeline_config.yaml, with cross-validation
against config/bifo_mapping_ddkg.yaml.

Public API:
    load_config(
        pipeline_config_path,
        bifo_mapping_path=None,
    ) -> PipelineConfig

If bifo_mapping_path is provided, the loader cross-validates that every
predicate listed under an enabled source's propagating_predicates_g_to_pw is:
    1. Present in bifo_mapping_ddkg.yaml's edge_resolution.predicate_to_flow.
    2. Has direction == 'source_to_target' (gene -> pathway).
    3. Classified as one of the propagating classes (NOT nonpropagating_context).

Cross-validation collects all errors before raising, so a single load call
reports the full set of inconsistencies.
"""
from pathlib import Path

import yaml
from pydantic import BaseModel, ConfigDict, Field, field_validator


# Predicate classifications considered propagating under BIFO. These are
# the classes that bifo_conditioning preserves when building the conditioned
# operator. Predicates outside this set (specifically nonpropagating_context)
# are dropped during conditioning.
PROPAGATING_CLASSIFICATIONS = frozenset({
    "mechanistic",
    "weak_mechanistic_or_observational",
    "observational",
    "contextual_constraint",
})

NONPROPAGATING_CLASSIFICATIONS = frozenset({
    "nonpropagating_context",
})

# Convenience: union of the two for "is this a known classification at all?"
_KNOWN_CLASSIFICATIONS = PROPAGATING_CLASSIFICATIONS | NONPROPAGATING_CLASSIFICATIONS


class PathwaySource(BaseModel):
    """One pathway source declared in pipeline_config.yaml."""

    model_config = ConfigDict(extra="forbid", frozen=True)

    sab: str  # injected by the loader from the YAML dict key
    enabled: bool = False
    propagating_predicates_g_to_pw: list[str] = Field(default_factory=list)
    kg_name_prefixes: list[str] = Field(default_factory=list)

    @field_validator("sab")
    @classmethod
    def _sab_nonempty_no_whitespace(cls, v: str) -> str:
        if not v or v.strip() != v or any(c in v for c in (" ", "\t", "\n")):
            raise ValueError(
                f"sab must be a non-empty string with no whitespace; got {v!r}"
            )
        return v

    @field_validator("propagating_predicates_g_to_pw")
    @classmethod
    def _validate_predicate_names(cls, v: list[str]) -> list[str]:
        for pred in v:
            if not pred or pred.strip() != pred:
                raise ValueError(
                    f"empty or whitespace-padded predicate name in {v!r}"
                )
        return v


class PipelineConfig(BaseModel):
    """Top-level config object loaded from pipeline_config.yaml."""

    model_config = ConfigDict(extra="forbid", frozen=True)

    schema_version: int = 1
    pathway_sources: dict[str, PathwaySource] = Field(default_factory=dict)
    gene_sabs: list[str] = Field(default_factory=lambda: ["HGNC"])
    burden_control_exclusions: list[str] = Field(default_factory=list)

    @field_validator("schema_version")
    @classmethod
    def _supported_schema_version(cls, v: int) -> int:
        if v != 1:
            raise ValueError(
                f"unsupported schema_version {v}; this loader supports "
                f"schema_version 1 only. If your pipeline_config.yaml is "
                f"newer than this loader, update bifo-graph; if you "
                f"intended to use the current schema, set 'schema_version: 1'."
            )
        return v

    @field_validator("pathway_sources")
    @classmethod
    def _enabled_must_have_predicates(
        cls, v: dict[str, PathwaySource]
    ) -> dict[str, PathwaySource]:
        for sab, source in v.items():
            if source.enabled and not source.propagating_predicates_g_to_pw:
                raise ValueError(
                    f"enabled pathway source {sab!r} has empty "
                    f"propagating_predicates_g_to_pw; an enabled source must "
                    f"declare at least one predicate"
                )
        return v

    def enabled_sources(self) -> list[PathwaySource]:
        """Return only the pathway sources with enabled=true."""
        return [s for s in self.pathway_sources.values() if s.enabled]

    def enabled_pathway_sabs(self) -> list[str]:
        """Return the SABs of enabled pathway sources, sorted alphabetically.

        Determinism matters: callers downstream may use the order to drive
        cypher UNION ordering or column ordering. Sorted output makes the
        same config produce the same downstream artifacts.
        """
        return sorted(s.sab for s in self.enabled_sources())


# -----------------------------------------------------------------------
# Cross-validation against bifo_mapping_ddkg.yaml
# -----------------------------------------------------------------------

class _BifoMappingPredicate(BaseModel):
    """Subset of bifo_mapping's predicate entry used for cross-validation.

    bifo_mapping_ddkg.yaml entries have many fields (flow, ro_curie, note,
    direction, classification, ...); we only care about direction and
    classification for our checks. extra='allow' so unknown fields don't
    break us.
    """
    model_config = ConfigDict(extra="allow")

    direction: str
    classification: str


def _load_bifo_mapping_predicates(
    bifo_mapping_path: Path,
) -> dict[str, _BifoMappingPredicate]:
    """Read bifo_mapping_ddkg.yaml and extract the predicate_to_flow map.

    Performs structural validation: raises a clear error if bifo_mapping
    is missing the expected sections, rather than returning an empty dict
    and producing misleading "predicate not found" errors downstream.

    Raises:
        FileNotFoundError: file does not exist.
        ValueError: YAML parse failure, or bifo_mapping is structurally
            malformed (missing edge_resolution, missing predicate_to_flow,
            or per-predicate entries malformed).
    """
    try:
        with open(bifo_mapping_path) as f:
            raw = yaml.safe_load(f)
    except yaml.YAMLError as e:
        raise ValueError(
            f"bifo_mapping at {bifo_mapping_path} is not valid YAML: {e}"
        ) from e

    if raw is None:
        raise ValueError(
            f"bifo_mapping at {bifo_mapping_path} is empty"
        )
    if not isinstance(raw, dict):
        raise ValueError(
            f"bifo_mapping at {bifo_mapping_path} is not a YAML mapping at "
            f"the top level (got {type(raw).__name__})"
        )

    edge_res = raw.get("edge_resolution")
    if edge_res is None:
        raise ValueError(
            f"bifo_mapping at {bifo_mapping_path} is missing the "
            f"'edge_resolution' section. This loader expects bifo_mapping "
            f"to follow the schema with edge_resolution.predicate_to_flow "
            f"containing per-predicate metadata."
        )
    if not isinstance(edge_res, dict):
        raise ValueError(
            f"bifo_mapping at {bifo_mapping_path}: 'edge_resolution' must "
            f"be a mapping; got {type(edge_res).__name__}"
        )

    pred_map = edge_res.get("predicate_to_flow")
    if pred_map is None:
        raise ValueError(
            f"bifo_mapping at {bifo_mapping_path} has 'edge_resolution' "
            f"but is missing 'edge_resolution.predicate_to_flow'. Expected "
            f"a mapping of predicate-name -> "
            f"{{direction, classification, ...}}."
        )
    if not isinstance(pred_map, dict):
        raise ValueError(
            f"bifo_mapping at {bifo_mapping_path}: "
            f"edge_resolution.predicate_to_flow must be a mapping; "
            f"got {type(pred_map).__name__}"
        )

    result: dict[str, _BifoMappingPredicate] = {}
    for name, body in pred_map.items():
        try:
            result[name] = _BifoMappingPredicate.model_validate(body)
        except Exception as e:
            raise ValueError(
                f"bifo_mapping at {bifo_mapping_path}: predicate {name!r} "
                f"is malformed: {e}"
            ) from e

    return result


def _cross_validate(
    config: PipelineConfig,
    bifo_predicates: dict[str, _BifoMappingPredicate],
    bifo_mapping_path: Path,
) -> None:
    """Cross-validate pipeline_config's predicates against bifo_mapping.

    Collects all errors before raising so users see the full set of
    inconsistencies in one shot.
    """
    errors: list[str] = []
    for source in config.enabled_sources():
        for pred in source.propagating_predicates_g_to_pw:
            if pred not in bifo_predicates:
                errors.append(
                    f"source {source.sab!r}: predicate {pred!r} is not "
                    f"defined in bifo_mapping (looked under "
                    f"edge_resolution.predicate_to_flow)"
                )
                continue
            entry = bifo_predicates[pred]
            if entry.direction != "source_to_target":
                errors.append(
                    f"source {source.sab!r}: predicate {pred!r} has direction "
                    f"{entry.direction!r} in bifo_mapping but is listed under "
                    f"propagating_predicates_g_to_pw, which requires "
                    f"direction 'source_to_target' (gene -> pathway)"
                )
            if entry.classification in NONPROPAGATING_CLASSIFICATIONS:
                errors.append(
                    f"source {source.sab!r}: predicate {pred!r} is classified "
                    f"as {entry.classification!r} in bifo_mapping (not "
                    f"propagating); cannot be used as a membership predicate"
                )
            elif entry.classification not in PROPAGATING_CLASSIFICATIONS:
                errors.append(
                    f"source {source.sab!r}: predicate {pred!r} has unknown "
                    f"classification {entry.classification!r} in bifo_mapping "
                    f"(expected one of {sorted(_KNOWN_CLASSIFICATIONS)})"
                )

    if errors:
        raise ValueError(
            f"pipeline_config is inconsistent with {bifo_mapping_path}:\n"
            + "\n".join(f"  - {e}" for e in errors)
        )


# -----------------------------------------------------------------------
# Public loader
# -----------------------------------------------------------------------

def load_config(
    pipeline_config_path: Path | str,
    bifo_mapping_path: Path | str | None = None,
) -> PipelineConfig:
    """Load and validate pipeline_config.yaml.

    Args:
        pipeline_config_path: path to pipeline_config.yaml.
        bifo_mapping_path: optional path to bifo_mapping_ddkg.yaml. If
            provided, the loader cross-validates that every enabled
            source's predicates are propagating-classified there. If None,
            cross-validation is skipped (use only in unit tests where
            bifo_mapping isn't relevant to the test).

    Raises:
        FileNotFoundError: pipeline_config.yaml or bifo_mapping_path not found.
        ValueError: YAML parse error, schema violation, structural problem
            in bifo_mapping, or cross-validation inconsistency. All
            ValueError raises include the file path and a human-readable
            description of the problem.
    """
    pipeline_config_path = Path(pipeline_config_path)

    try:
        with open(pipeline_config_path) as f:
            raw = yaml.safe_load(f) or {}
    except yaml.YAMLError as e:
        raise ValueError(
            f"pipeline_config at {pipeline_config_path} is not valid YAML: {e}"
        ) from e

    # The YAML uses SAB as the dict key under pathway_sources; inject it
    # as a `sab` field in each source body so PathwaySource can validate it.
    if "pathway_sources" in raw and raw["pathway_sources"]:
        if not isinstance(raw["pathway_sources"], dict):
            raise ValueError(
                f"pipeline_config at {pipeline_config_path}: "
                f"'pathway_sources' must be a mapping; "
                f"got {type(raw['pathway_sources']).__name__}"
            )
        raw["pathway_sources"] = {
            sab: {**(body or {}), "sab": sab}
            for sab, body in raw["pathway_sources"].items()
        }

    try:
        config = PipelineConfig.model_validate(raw)
    except Exception as e:
        # Pydantic ValidationError doesn't include file path; wrap it.
        # Re-raise as ValueError so callers handle one exception type.
        raise ValueError(
            f"pipeline_config at {pipeline_config_path} failed schema "
            f"validation:\n{e}"
        ) from e

    if bifo_mapping_path is not None:
        bifo_mapping_path = Path(bifo_mapping_path)
        bifo_predicates = _load_bifo_mapping_predicates(bifo_mapping_path)
        _cross_validate(config, bifo_predicates, bifo_mapping_path)

    return config
