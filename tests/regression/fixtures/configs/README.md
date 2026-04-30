# Test fixtures: configs

Small YAML fixtures used by `tests/unit/test_pipeline_config.py`.

Each fixture exercises a specific schema or cross-validation case.
Keeping fixtures as checked-in files (rather than inline strings in
test code) makes them inspectable, diffable, and reproducible.

## pipeline_config fixtures

| File | Purpose |
|---|---|
| `pipeline_config_minimal.yaml` | Smallest valid config (one enabled source, GO disabled). Happy path for loader and cross-validation tests. |
| `pipeline_config_invalid_extra_field.yaml` | Has an unknown top-level field. Should fail with `ValidationError` (extra='forbid'). |
| `pipeline_config_invalid_enabled_no_predicates.yaml` | Source enabled but with empty `propagating_predicates_g_to_pw`. Should fail with `ValueError` from cross-field validator. |
| `pipeline_config_invalid_sab_whitespace.yaml` | SAB key contains whitespace. Should fail with `ValueError` from SAB field validator. |
| `pipeline_config_invalid_unsupported_schema_version.yaml` | `schema_version` is not 1. Should fail with `ValueError` from schema_version validator. |
| `pipeline_config_two_predicates.yaml` | Lists two predicates so cross-validation can hit two errors at once (paired with `bifo_mapping_two_errors.yaml`). |

## bifo_mapping fixtures

These are minimal versions of `config/bifo_mapping_ddkg.yaml`, containing
only the `edge_resolution.predicate_to_flow` subsection that the loader's
cross-validation reads. Real bifo_mapping files have many more sections;
the loader uses `extra='allow'` on its predicate model, so unknown fields
are accepted but ignored.

| File | Purpose |
|---|---|
| `bifo_mapping_minimal.yaml` | Defines the predicates referenced in `pipeline_config_minimal.yaml`, classified correctly. Cross-validation passes. |
| `bifo_mapping_predicate_missing.yaml` | Doesn't define a predicate that pipeline_config references. Cross-validation should report "predicate not found". |
| `bifo_mapping_wrong_direction.yaml` | Defines the predicate but with `direction: target_to_source`. Cross-validation should fail. |
| `bifo_mapping_nonpropagating_classification.yaml` | Defines the predicate but classified as `nonpropagating_context`. Cross-validation should fail. |
| `bifo_mapping_unknown_classification.yaml` | Defines the predicate but with a classification not in the known taxonomy. Cross-validation should fail. |
| `bifo_mapping_two_errors.yaml` | Two predicates, each with a different problem. Cross-validation should report both in one error. |

## Adding fixtures

When adding a new fixture:

1. Use a name that describes what the fixture does, not how the test uses it.
2. Add a header comment explaining the fixture's purpose and the expected
   loader/validator behavior against it.
3. Add a row to this README's tables.
4. Reference from the test by relative path; never inline the content.
