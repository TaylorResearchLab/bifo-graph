# Configuring pathway sources

This document explains how to configure pathway sources for the bifo-graph
pipeline. It is the reference for editing `config/pipeline_config.yaml` and
its companion `config/bifo_mapping_ddkg.yaml`.

## Overview: two configuration files

The pipeline reads two YAML files at startup:

| File | Role | Edit when |
|---|---|---|
| `config/pipeline_config.yaml` | Policy: which sources to use, which prefixes, which seeds to exclude | You want to change pipeline behavior |
| `config/bifo_mapping_ddkg.yaml` | KG description: predicate inventory, semantic classifications | The underlying knowledge graph build changes |

These files are validated together at startup. If they disagree (for
example, `pipeline_config.yaml` references a predicate that doesn't exist
in `bifo_mapping_ddkg.yaml`), the pipeline refuses to run and reports the
inconsistency.

To check your configuration without running the full pipeline:

```bash
python -m pipeline.check_configs
```

This pre-flight tool validates both files, runs cross-validation, and
reports any problems in human-readable form.

## When to edit `pipeline_config.yaml`

This is the file most users edit. Common reasons:

- Enable or disable a pathway source (e.g., turn on Gene Ontology).
- Add a new pathway-name prefix filter (e.g., add `MIR_` for miRNA targets).
- Apply a burden-control exclusion list to seed genes.
- Adapt the pipeline to a different gene SAB (e.g., GENCODE instead of HGNC).

### Anatomy of `pipeline_config.yaml`

```yaml
schema_version: 1                        # do not change unless upgrading

pathway_sources:                         # registry of pathway sources
  MSIGDB:                                # SAB name from the KG
    enabled: true                        # source is active
    propagating_predicates_g_to_pw:      # KG predicate names: gene -> pathway
      - inverse_pathway_associated_with_gene
      - inverse_has_signature_gene
    kg_name_prefixes:                    # filter pathways by name
      - HALLMARK_
      - REACTOME_
      - WP_

  GO:
    enabled: false                       # source is inert
    propagating_predicates_g_to_pw:
      - gene_plays_role_in_process
      - gene_product_plays_role_in_biological_process
    kg_name_prefixes: []                 # empty = accept all

gene_sabs:                               # gene Concept SAB(s)
  - HGNC

burden_control_exclusions: []            # empty = no exclusion
```

### Field-by-field

**`schema_version`**: integer, must be `1`. Reserved for future schema changes.

**`pathway_sources`**: mapping where each key is a SAB name from the KG.
Each value is a record with `enabled`, `propagating_predicates_g_to_pw`,
and `kg_name_prefixes`.

- An enabled source must declare at least one propagating predicate. The
  pre-flight tool will reject an enabled source with an empty predicate list.
- A disabled source is documented but inert. The pipeline ignores its
  predicates and prefixes.

**`propagating_predicates_g_to_pw`**: list of KG predicate names that
connect a gene Concept to a pathway Concept. The direction matters: these
predicates must be the form that flows gene -> pathway. See "Direction
conventions" below.

**`kg_name_prefixes`**: list of pathway-name prefixes used to filter the
source. For MSigDB collections in DDKG/UBKG, names like `HALLMARK_HYPOXIA`
or `REACTOME_DNA_REPAIR` carry a collection prefix. Empty list means
no filter applied (accept all pathways from this source).

**`gene_sabs`**: list of SAB(s) the pipeline treats as gene Concepts.
HGNC is the canonical human gene SAB in DDKG/UBKG. Other graph builds
might use GENCODE, ENSEMBL, or another convention.

**`burden_control_exclusions`**: list of gene symbols (HGNC symbols) to
exclude from seed lists before CUI resolution. Empty by default. Use this
to remove genes commonly excluded from rare-variant burden tests due to
broad curation, polymorphism tolerance, or known recurrence as
false-positive P/LP variant calls.

## Direction conventions

The most common confusion when adding a predicate is direction.

**The rule for `propagating_predicates_g_to_pw`**: the predicate must
flow source -> target where source is a gene Concept and target is a
pathway Concept.

In DDKG/UBKG, predicates are often represented in both directions:

| Forward (pathway -> gene) | Inverse (gene -> pathway) |
|---|---|
| `pathway_associated_with_gene` | `inverse_pathway_associated_with_gene` |
| `has_signature_gene` | `inverse_has_signature_gene` |

For BIFO-PPR, you want the **inverse** form: gene -> pathway.

This is how BIFO-PPR's pathway membership query (Q5) extracts membership
edges. If you list the forward (`pathway_associated_with_gene`) form,
the cross-validator rejects it with:

```
predicate 'pathway_associated_with_gene' has direction 'target_to_source'
in bifo_mapping but is listed under propagating_predicates_g_to_pw,
which requires direction 'source_to_target' (gene -> pathway)
    --> Did you mean 'inverse_pathway_associated_with_gene'?
```

## Classification taxonomy

`bifo_mapping_ddkg.yaml` classifies every predicate as one of:

| Classification | Propagating? | Use as membership predicate? |
|---|---|---|
| `mechanistic` | yes | yes |
| `weak_mechanistic_or_observational` | yes | yes |
| `observational` | yes | yes |
| `contextual_constraint` | yes | yes |
| `nonpropagating_context` | no | no |

A predicate listed in `propagating_predicates_g_to_pw` must be one of the
four propagating classes. The pre-flight tool will reject a
`nonpropagating_context` predicate.

## Adding a new pathway source

To add a new source (example: a custom `CUSTOM_PW` SAB):

1. **Confirm the SAB exists in your KG.** Run a Cypher query against
   the KG to count Concepts with `SAB == 'CUSTOM_PW'`.

2. **Identify the membership predicate(s).** Find the KG predicate name
   that connects gene Concepts to your `CUSTOM_PW` Concepts in the
   gene -> pathway direction.

3. **Confirm the predicate is in `bifo_mapping_ddkg.yaml`.** Search
   under `edge_resolution.predicate_to_flow`. If it's not there, add it
   (consult the lab's Q5 design notes for the right `flow` value and
   `classification`).

4. **Edit `pipeline_config.yaml`.** Add a new entry under `pathway_sources`:

   ```yaml
     CUSTOM_PW:
       enabled: true
       propagating_predicates_g_to_pw:
         - inverse_custom_membership_predicate
       kg_name_prefixes:
         - CUSTOM_              # or empty list for no filter
   ```

5. **Run the pre-flight tool** to confirm consistency:

   ```bash
   python -m pipeline.check_configs
   ```

   Fix any errors it reports.

6. **Run the pipeline.** Output artifacts will now include `CUSTOM_PW`
   pathways alongside MSIGDB.

## Adapting to a non-DDKG knowledge graph

The `pipeline_config.yaml` and `bifo_mapping_ddkg.yaml` pair is specific
to one KG build. To use a different KG:

1. Build a `bifo_mapping_<kgname>.yaml` describing the new KG's predicates
   under `edge_resolution.predicate_to_flow`. Each predicate needs:
   `direction` (source_to_target or target_to_source) and `classification`
   (one of the five values in the taxonomy above).

2. Edit `pipeline_config.yaml` to reference the new SAB names, predicate
   names, and pathway-name prefixes appropriate for your KG.

3. Run the pre-flight tool with `--bifo-mapping path/to/bifo_mapping_<kgname>.yaml`
   to confirm consistency.

The pipeline does not perform automatic translation between KGs. The
mapping file you provide is taken as ground truth for what's in the KG.

## Interpreting pre-flight error output

The pre-flight tool prints two kinds of lines:

- `[OK]` lines indicate a check passed.
- `[FAIL]` lines indicate a check failed.

For cross-validation errors, errors are grouped by SAB and numbered:

```
[FAIL] Cross-validation found problems:

  source MSIGDB:
    [1] predicate 'X' is not defined in bifo_mapping ...
        --> Did you mean 'Y'?
    [2] predicate 'Z' has direction 'target_to_source' ...
        --> Did you mean 'inverse_Z'?
```

Each error names the source, the offending predicate, and (where
applicable) a "Did you mean" hint that points at the most likely fix.

## Exit codes

The pre-flight tool returns one of three exit codes:

| Code | Meaning |
|---|---|
| 0 | All checks passed |
| 1 | Configuration errors (cross-validation, schema, structure) |
| 2 | Setup errors (file not found, unreadable) |

This makes the tool composable with shell scripts:

```bash
python -m pipeline.check_configs --quiet && bash scripts/run_full_pipeline.sh
```

## Related files

- `config/pipeline_config.yaml` - the policy file you edit
- `config/bifo_mapping_ddkg.yaml` - the KG description
- `pipeline/pipeline_config.py` - the loader (Pydantic v2)
- `pipeline/check_configs.py` - the pre-flight tool

## See also

- `REPRODUCE.md` - end-to-end pipeline reproduction
- `BENCHMARK_MANIFEST.md` - benchmark configuration history
