"""
bifo_conditioning.py  v0.3.0  (final merged)

Merge provenance:
  - CoverageReport dataclass        from bifo_conditioning_corrected.py
  - _get_first column helper        from bifo_conditioning_corrected.py
  - build_directionality_only_edges from bifo_conditioning_corrected.py
  - hard fail on missing seeds      from bifo_conditioning_corrected.py
  - --alpha / --go-obo CLI args     from bifo_conditioning.py (v0.2.0)
  - graceful empty-operator warn    from bifo_conditioning.py (v0.2.0)

Fixes applied in this merge:
  BUG 1 (CRITICAL): resolve_flow entity-pair overrides hardcoded "mechanistic",
    ignoring the classification field in the YAML. Observational overrides
    (e.g. G→DS gene_associated_with_disease_or_phenotype) were propagating.
    Fixed: rule.get("classification", "mechanistic") is now read from YAML.
  BUG 2 (METHODOLOGICAL): random sparsification sampled from directionality_edges,
    not raw_concept_edges. This confounded two variables simultaneously.
    Fixed: sparsification pool is always raw_concept_edges.
  BUG 3 (GO): label heuristic ran even when namespace was available, risking
    override of a correct namespace-based assignment.
    Fixed: namespace check is primary; label heuristic runs only when
    namespace is absent.

Expected input:
  nodes CSV/TSV: node_id, label, [sab, code, code_id, concept_id, parent_concept_id]
  edges CSV/TSV: source, target, predicate, [sab, confidence, evidence_type, value]
"""

from __future__ import annotations

import argparse

def removesuffix(s, suffix):
    return s[:-len(suffix)] if suffix and s.endswith(suffix) else s
import json
import logging
import re
from collections import defaultdict
from dataclasses import dataclass, asdict
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Set, Tuple

import numpy as np
import pandas as pd
import yaml
from scipy import sparse
from sklearn.metrics import average_precision_score, roc_auc_score

logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
log = logging.getLogger("bifo")

CONCEPT_LABEL = "concept"
CODE_LABEL = "code"


def clean_node_id(raw_id: str) -> str:
    """Strip DDKG schema artifacts from Concept CUI strings.
    In DDKG/UBKG, some Concept CUIs are stored with a trailing ' CUI' suffix
    (e.g. 'MSIGDB:M39610 CUI', 'HGNC:20343 CUI'). Strip universally at load time.
    """
    return removesuffix(raw_id.strip(), ' CUI')


# ---------------------------------------------------------------------------
# Coverage reporting
# ---------------------------------------------------------------------------

@dataclass
class CoverageReport:
    total_nodes: int
    total_concept_nodes: int
    total_code_nodes: int
    resolved_code_nodes: int
    unresolved_code_nodes: int
    resolved_concept_nodes: int
    unresolved_concept_nodes: int
    concept_fallthrough_count: int
    level1_concept_coverage: float
    total_edges: int
    total_concept_edges: int
    concept_edges_with_resolved_endpoints: int
    kept_edges: int
    dropped_edges: int
    dropped_non_flow: int
    dropped_unmapped_predicate: int
    dropped_unresolved_entity: int
    level2_edge_coverage: float


@dataclass
class ConditioningResult:
    kept_edges: pd.DataFrame
    dropped_edges: pd.DataFrame
    node_entity_map: Dict[str, str]
    edge_flow_map: Dict[int, str]
    coverage: CoverageReport
    concept_fallthrough_nodes: List[str]


# ---------------------------------------------------------------------------
# I/O helpers
# ---------------------------------------------------------------------------

def load_mapping(path: str | Path) -> dict:
    with open(path, "r", encoding="utf-8") as fh:
        return yaml.safe_load(fh)


def load_table(path: str | Path) -> pd.DataFrame:
    path = Path(path)
    sep = "\t" if path.suffix.lower() in {".tsv", ".tab"} else ","
    # encoding=utf-8-sig strips BOM that Neo4j sometimes prepends to CSV exports
    return pd.read_csv(path, sep=sep, dtype=str, encoding="utf-8-sig").fillna("")


def _read_node_list(path: str | Path) -> List[str]:
    """Read node IDs from a file, one per line.
    Strips comments (# prefix lines and inline # comments), blank lines,
    and DDKG CUI suffix artifacts via clean_node_id().
    """
    ids = []
    for line in Path(path).read_text().splitlines():
        line = line.strip()
        if not line or line.startswith('#'):
            continue
        token = line.split('#')[0].strip().split()[0] if line.split('#')[0].strip() else ''
        if token:
            ids.append(clean_node_id(token))
    return ids


def _get_first(row: pd.Series, candidates: Iterable[str]) -> str:
    """Return the first non-empty value from candidate column names."""
    for c in candidates:
        if c in row.index and row[c]:
            return str(row[c])
    return ""


# ---------------------------------------------------------------------------
# GO subontology lookup
# ---------------------------------------------------------------------------

def build_go_lookup(obo_path: Optional[str]) -> Dict[str, str]:
    """
    Parse GO OBO file → dict: GO accession → subontology.
    Values: 'biological_process', 'molecular_function', 'cellular_component', 'complex'.
    Returns {} if no OBO file; GO terms fall back to PW.
    """
    lookup: Dict[str, str] = {}
    if not obo_path or not Path(obo_path).exists():
        log.warning("GO OBO not provided or not found — GO terms will default to PW.")
        return lookup

    current_id: Optional[str] = None
    current_ns: Optional[str] = None
    current_name: Optional[str] = None
    in_term = False

    with open(obo_path, encoding="utf-8") as fh:
        for line in fh:
            line = line.strip()
            if line == "[Term]":
                in_term = True
                current_id = current_ns = current_name = None
            elif line == "" and in_term:
                if current_id and current_ns:
                    ns = current_ns
                    if (ns == "cellular_component" and current_name
                            and "complex" in current_name.lower()):
                        ns = "complex"
                    lookup[current_id] = ns
                in_term = False
            elif in_term:
                if line.startswith("id: "):
                    current_id = line[4:].strip()
                elif line.startswith("namespace: "):
                    current_ns = line[11:].strip()
                elif line.startswith("name: "):
                    current_name = line[6:].strip()

    log.info("GO OBO: %d terms indexed", len(lookup))
    return lookup


def _go_entity_from_lookup(go_id: str, go_lookup: Dict[str, str]) -> Optional[str]:
    """Primary GO resolution via OBO-parsed namespace."""
    ns = go_lookup.get(go_id)
    if ns in ("biological_process", "molecular_function"):
        return "PW"
    if ns == "complex":
        return "CM"
    if ns == "cellular_component":
        return "S"
    return None


def _go_entity_from_row_heuristic(row: pd.Series, mapping: dict) -> Optional[str]:
    """
    Fallback GO resolution using namespace/label columns on the row.
    Only called when OBO lookup returns nothing.
    Namespace check runs first; label heuristic only when namespace absent.
    """
    go_cfg = (mapping.get("node_resolution", {})
              .get("sab_overrides", {})
              .get("GO", {}))
    namespace_cols = go_cfg.get("namespace_columns", [])
    label_cols = go_cfg.get("label_columns", [])

    namespace_text = _get_first(row, namespace_cols).strip().lower()
    label_text = _get_first(row, label_cols).strip().lower()

    # Namespace check is primary
    if namespace_text:
        is_cc = (namespace_text in {"cellular_component", "go_component", "component"}
                 or "cellular component" in namespace_text)
        if is_cc:
            if go_cfg.get("complex_terms_to_CM") and "complex" in label_text:
                return "CM"
            if go_cfg.get("compartment_terms_to_S"):
                return "S"
        if (namespace_text in {"biological_process", "go_process", "process"}
                or "biological process" in namespace_text):
            if go_cfg.get("biological_process_terms_to_PW"):
                return "PW"
        if (namespace_text in {"molecular_function", "go_function", "function"}
                or "molecular function" in namespace_text):
            if go_cfg.get("molecular_function_terms_to_PW"):
                return "PW"
        # Namespace present but unrecognised — fall through to default
        return None

    # Label heuristic only when namespace is absent
    if label_text:
        if go_cfg.get("complex_terms_to_CM") and "complex" in label_text:
            return "CM"
        if go_cfg.get("compartment_terms_to_S") and any(
            t in label_text for t in [
                "membrane", "nucleus", "cytoplasm", "cytosol",
                "mitochond", "extracellular", "organelle",
            ]
        ):
            return "S"

    return None


# ---------------------------------------------------------------------------
# Entity resolution
# ---------------------------------------------------------------------------

def _regex_entity(code_id: str, regex_rules: List[dict]) -> Optional[str]:
    for rule in regex_rules:
        if re.match(rule["pattern"], code_id):
            return rule["entity"]
    return None


def resolve_code_entity(
    row: pd.Series,
    mapping: dict,
    go_lookup: Dict[str, str],
) -> Optional[str]:
    sab = _get_first(row, ["sab", "SAB"])
    code = _get_first(row, ["code", "CODE"])
    code_id = _get_first(row, ["code_id", "CodeID", "node_id"])
    if not code_id and sab and code:
        code_id = f"{sab}:{code}"

    nr = mapping["node_resolution"]
    ignored: set = set(nr.get("ignored_or_auxiliary_sabs", []))
    if sab in ignored:
        return None

    # GO subontology override — before generic SAB lookup
    if sab in ("GO", "GOCC") or code_id.startswith("GO:"):
        if go_lookup and code_id:
            entity = _go_entity_from_lookup(code_id, go_lookup)
            if entity:
                return entity
        entity = _go_entity_from_row_heuristic(row, mapping)
        if entity:
            return entity
        return "PW"  # final fallback

    # ENSEMBL prefix disambiguation
    if sab == "ENSEMBL" and code:
        if code.startswith("ENSG"):
            return "G"
        if code.startswith("ENST"):
            return "RNA"
        if code.startswith("ENSP"):
            return "P"

    sab_map: dict = nr["sab_to_entity"]
    if sab in sab_map:
        return sab_map[sab]

    if code_id:
        entity = _regex_entity(code_id, nr["codeid_regex_to_entity"])
        if entity:
            return entity

    return None


def infer_concept_entities(
    nodes: pd.DataFrame,
    mapping: dict,
    go_lookup: Dict[str, str],
) -> Tuple[Dict[str, str], str, Dict[str, List[str]]]:
    priority: list = mapping["node_resolution"]["priority_order"]
    priority_rank = {entity: i for i, entity in enumerate(priority)}

    code_nodes = nodes[nodes["label"].str.lower() == CODE_LABEL].copy()
    if code_nodes.empty:
        log.warning("No Code nodes found — falling back to direct SAB column on Concept nodes.")

        # Fallback: use 'sab' column directly on Concept nodes.
        # This handles exports where each Concept row carries its own SAB
        # (e.g. single-layer exports from Neo4j without Code node rows).
        concept_nodes = nodes[nodes["label"].str.lower() == CONCEPT_LABEL].copy()
        if "sab" not in concept_nodes.columns and "SAB" not in concept_nodes.columns:
            log.warning("No 'sab' column on Concept nodes either — resolution will be empty.")
            return {}, "", {}

        sab_col = "sab" if "sab" in concept_nodes.columns else "SAB"
        id_col = "node_id" if "node_id" in concept_nodes.columns else "id"

        nr = mapping["node_resolution"]
        sab_map: dict = nr["sab_to_entity"]
        ignored: set = set(nr.get("ignored_or_auxiliary_sabs", []))
        regex_rules = nr.get("codeid_regex_to_entity", [])

        concept_resolved: Dict[str, str] = {}
        concept_to_entities: Dict[str, List[str]] = {}

        for _, row in concept_nodes.iterrows():
            concept_id = clean_node_id(str(row.get(id_col, "")))
            sab = str(row.get(sab_col, "")).strip()
            node_id_val = str(row.get("node_id", concept_id))

            if not concept_id or sab in ignored:
                continue

            entity = None

            # GO subontology override
            if sab in ("GO", "GOCC"):
                if go_lookup and concept_id:
                    entity = _go_entity_from_lookup(concept_id, go_lookup)
                if not entity:
                    entity = _go_entity_from_row_heuristic(row, mapping)
                if not entity:
                    entity = "PW"

            # ENSEMBL prefix disambiguation
            elif sab == "ENSEMBL":
                code_val = concept_id.split(":")[-1] if ":" in concept_id else concept_id
                if code_val.startswith("ENSG"):
                    entity = "G"
                elif code_val.startswith("ENST"):
                    entity = "RNA"
                elif code_val.startswith("ENSP"):
                    entity = "P"

            # Direct SAB lookup
            if not entity and sab in sab_map:
                entity = sab_map[sab]

            # Regex fallback on concept_id
            if not entity and concept_id:
                entity = _regex_entity(concept_id, regex_rules)

            if entity:
                concept_resolved[concept_id] = entity
                concept_to_entities[concept_id] = [entity]

        log.info("Direct SAB fallback resolved %d / %d concept nodes",
                 len(concept_resolved), len(concept_nodes))
        return concept_resolved, id_col, concept_to_entities

    code_nodes["resolved_entity"] = code_nodes.apply(
        resolve_code_entity, axis=1, mapping=mapping, go_lookup=go_lookup
    )

    concept_key = ""
    for candidate in ("parent_concept_id", "concept_id", "Concept", "concept",
                      "parent", "concept_cui"):
        if candidate in code_nodes.columns:
            concept_key = candidate
            break

    if not concept_key:
        log.warning("No concept linkage column found — concept resolution will be empty.")
        return {}, "", {}

    concept_to_entities: Dict[str, List[str]] = defaultdict(list)
    for _, row in code_nodes.iterrows():
        # Clean concept_id so keys match cleaned lookups in build_node_entity_map
        concept_id = clean_node_id(str(row.get(concept_key, "")))
        ent = row.get("resolved_entity", "")
        if concept_id and ent:
            concept_to_entities[concept_id].append(ent)

    concept_resolved: Dict[str, str] = {}
    for concept_id, entities in concept_to_entities.items():
        uniq = sorted(set(entities), key=lambda x: priority_rank.get(x, 999))
        if uniq:
            concept_resolved[concept_id] = uniq[0]

    return concept_resolved, concept_key, dict(concept_to_entities)


def build_node_entity_map(
    nodes: pd.DataFrame,
    mapping: dict,
    go_lookup: Dict[str, str],
) -> Tuple[Dict[str, str], CoverageReport, List[str]]:
    concept_entities, _, _ = infer_concept_entities(nodes, mapping, go_lookup)
    node_entity_map: Dict[str, str] = {}
    concept_fallthrough: List[str] = []

    total_nodes = int(nodes.shape[0])
    concept_mask = nodes["label"].str.lower() == CONCEPT_LABEL
    code_mask = nodes["label"].str.lower() == CODE_LABEL
    total_concepts = int(concept_mask.sum())
    total_codes = int(code_mask.sum())

    resolved_concepts = unresolved_concepts = 0
    resolved_code_nodes = unresolved_code_nodes = 0

    for _, row in nodes.iterrows():
        node_id = clean_node_id(str(row["node_id"]))
        label = str(row["label"]).lower()

        if label == CONCEPT_LABEL:
            concept_id = clean_node_id(_get_first(row, ["node_id", "concept_id"]))
            entity = concept_entities.get(concept_id)
            if entity:
                node_entity_map[node_id] = entity
                resolved_concepts += 1
            else:
                unresolved_concepts += 1
                concept_fallthrough.append(node_id)

        elif label == CODE_LABEL:
            entity = resolve_code_entity(row, mapping, go_lookup)
            if entity:
                node_entity_map[node_id] = entity
                resolved_code_nodes += 1
            else:
                unresolved_code_nodes += 1

    coverage = CoverageReport(
        total_nodes=total_nodes,
        total_concept_nodes=total_concepts,
        total_code_nodes=total_codes,
        resolved_code_nodes=resolved_code_nodes,
        unresolved_code_nodes=unresolved_code_nodes,
        resolved_concept_nodes=resolved_concepts,
        unresolved_concept_nodes=unresolved_concepts,
        concept_fallthrough_count=len(concept_fallthrough),
        level1_concept_coverage=(
            resolved_concepts / total_concepts if total_concepts else float("nan")
        ),
        # Edge fields populated later in condition_edges
        total_edges=0,
        total_concept_edges=0,
        concept_edges_with_resolved_endpoints=0,
        kept_edges=0,
        dropped_edges=0,
        dropped_non_flow=0,
        dropped_unmapped_predicate=0,
        dropped_unresolved_entity=0,
        level2_edge_coverage=float("nan"),
    )

    log.info(
        "Level 1 — concepts: %d / %d resolved (%.1f%%), fallthrough: %d | "
        "codes: %d / %d resolved",
        resolved_concepts, total_concepts,
        coverage.level1_concept_coverage * 100 if total_concepts else 0,
        len(concept_fallthrough),
        resolved_code_nodes, total_codes,
    )
    return node_entity_map, coverage, concept_fallthrough


# ---------------------------------------------------------------------------
# Flow resolution
# ---------------------------------------------------------------------------

def resolve_flow(
    predicate: str,
    source_entity: Optional[str],
    target_entity: Optional[str],
    mapping: dict,
) -> Tuple[Optional[str], str]:
    """
    Returns (flow_class, classification).

    FIX v0.3.0: entity-pair overrides read classification from the YAML rule.
    Previously hardcoded 'mechanistic', causing observational overrides
    (e.g. G→DS) to propagate incorrectly.
    """
    edge_cfg = mapping["edge_resolution"]

    if predicate in edge_cfg["non_flow_edges"]:
        return None, "non_flow"

    if predicate in edge_cfg["observational_edges"]:
        return "Observational Association", "observational"

    # Entity-pair overrides take precedence — classification comes from YAML
    for rule in edge_cfg.get("entity_pair_overrides", []):
        if (
            rule["predicate"] == predicate
            and rule.get("source_entity") == source_entity
            and rule.get("target_entity") == target_entity
        ):
            return rule["flow"], rule.get("classification", "mechanistic")

    pred_map = edge_cfg["predicate_to_flow"]
    if predicate in pred_map:
        item = pred_map[predicate]
        return item["flow"], item.get("classification", "mechanistic")

    return None, "unmapped"


# ---------------------------------------------------------------------------
# Core conditioning
# ---------------------------------------------------------------------------

def condition_edges(
    nodes: pd.DataFrame,
    edges: pd.DataFrame,
    mapping: dict,
    go_lookup: Dict[str, str],
) -> ConditioningResult:
    node_entity_map, coverage, concept_fallthrough = build_node_entity_map(
        nodes, mapping, go_lookup
    )
    concept_ids = set(
        [clean_node_id(str(x)) for x in nodes.loc[nodes["label"].str.lower() == CONCEPT_LABEL, "node_id"].tolist()]
    )

    kept_rows: List[dict] = []
    dropped_rows: List[dict] = []
    edge_flow_map: Dict[int, str] = {}

    total_edges = int(edges.shape[0])
    # Clean edge endpoint IDs before filtering so artifacted CUIs
    # are not silently excluded before condition_edges() sees them
    edges = edges.copy()
    edges["source"] = edges["source"].map(lambda x: clean_node_id(str(x)))
    edges["target"] = edges["target"].map(lambda x: clean_node_id(str(x)))
    concept_edge_mask = (
        edges["source"].isin(concept_ids) & edges["target"].isin(concept_ids)
    )
    concept_edges_df = edges[concept_edge_mask].copy()
    total_concept_edges = int(concept_edges_df.shape[0])

    resolved_endpoints = dropped_non_flow = dropped_unmapped = dropped_unresolved = 0

    for idx, row in concept_edges_df.iterrows():
        src, dst, predicate = (
            clean_node_id(str(row["source"])),
            clean_node_id(str(row["target"])),
            row["predicate"],
        )
        src_entity = node_entity_map.get(src)
        dst_entity = node_entity_map.get(dst)

        if src_entity and dst_entity:
            resolved_endpoints += 1

        flow, classification = resolve_flow(predicate, src_entity, dst_entity, mapping)

        row_out = row.to_dict()
        row_out.update({
            "source": src,          # overwrite with cleaned ID
            "target": dst,          # overwrite with cleaned ID
            "source_entity": src_entity or "",
            "target_entity": dst_entity or "",
            "bifo_flow": flow or "",
            "classification": classification,
        })

        if not src_entity or not dst_entity:
            row_out["drop_reason"] = "unresolved_entity"
            dropped_rows.append(row_out)
            dropped_unresolved += 1
            continue

        if classification == "non_flow":
            row_out["drop_reason"] = "non_flow"
            dropped_rows.append(row_out)
            dropped_non_flow += 1
            continue

        if classification == "unmapped":
            row_out["drop_reason"] = "unmapped_predicate"
            dropped_rows.append(row_out)
            dropped_unmapped += 1
            continue

        row_out["propagating"] = classification in {
            "mechanistic",
            "weak_mechanistic_or_observational",
        }
        kept_rows.append(row_out)
        edge_flow_map[idx] = flow or ""

    coverage.total_edges = total_edges
    coverage.total_concept_edges = total_concept_edges
    coverage.concept_edges_with_resolved_endpoints = resolved_endpoints
    coverage.kept_edges = len(kept_rows)
    coverage.dropped_edges = len(dropped_rows)
    coverage.dropped_non_flow = dropped_non_flow
    coverage.dropped_unmapped_predicate = dropped_unmapped
    coverage.dropped_unresolved_entity = dropped_unresolved
    coverage.level2_edge_coverage = (
        len(kept_rows) / total_concept_edges if total_concept_edges else float("nan")
    )

    log.info(
        "Level 2 — kept: %d / %d concept edges (%.1f%%) | "
        "dropped: non_flow=%d, unmapped=%d, unresolved=%d",
        coverage.kept_edges, total_concept_edges,
        coverage.level2_edge_coverage * 100 if total_concept_edges else 0,
        dropped_non_flow, dropped_unmapped, dropped_unresolved,
    )

    return ConditioningResult(
        kept_edges=pd.DataFrame(kept_rows) if kept_rows else pd.DataFrame(),
        dropped_edges=pd.DataFrame(dropped_rows) if dropped_rows else pd.DataFrame(),
        node_entity_map=node_entity_map,
        edge_flow_map=edge_flow_map,
        coverage=coverage,
        concept_fallthrough_nodes=concept_fallthrough,
    )


# ---------------------------------------------------------------------------
# Directionality-only control
# ---------------------------------------------------------------------------

def build_metadata_filtered_edges(
    concept_edges: pd.DataFrame, mapping: dict
) -> pd.DataFrame:
    """
    Strip metadata/hierarchy edges (CODE, is_a, part_of, etc.) but retain all
    remaining edges as propagating, including observational and unmapped edges.

    This is the metadata_filtered control arm. It tests whether removing
    annotation noise alone — without any biological admissibility filtering —
    improves signal recovery. It is not a directionality control: edges are
    not reoriented from YAML direction fields. Edge direction is inherited
    from the raw export as-is.

    Comparison to other arms:
      raw:                all concept-to-concept edges, uniform weight
      metadata_filtered:  non-metadata edges only, uniform weight, all propagating
      conditioned:        BIFO-admissible mechanistic edges only, typed weights
      random_sparsified:  random sample from raw matching conditioned edge count
    """
    non_flow = set(mapping["edge_resolution"]["non_flow_edges"])
    work = concept_edges[~concept_edges["predicate"].isin(non_flow)].copy()
    if work.empty:
        return work
    work["propagating"] = True
    work["classification"] = "metadata_filtered"
    work["bifo_flow"] = ""
    return work


# ---------------------------------------------------------------------------
# Graph construction
# ---------------------------------------------------------------------------

def concept_subgraph(
    nodes: pd.DataFrame, edges: pd.DataFrame
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    concept_nodes = nodes[nodes["label"].str.lower() == CONCEPT_LABEL].copy()
    concept_ids = {clean_node_id(str(x)) for x in concept_nodes["node_id"].tolist()}
    if edges.empty:
        return concept_nodes, pd.DataFrame()
    mask = edges["source"].isin(concept_ids) & edges["target"].isin(concept_ids)
    return concept_nodes, edges[mask].copy()


def edge_weights(edges: pd.DataFrame) -> np.ndarray:
    tier = {"experimental": 1.0, "curated": 0.75, "inferred": 0.5, "predicted": 0.25}
    weights: List[float] = []
    for _, row in edges.iterrows():
        conf = _get_first(row, ["confidence", "value"])
        if conf:
            try:
                weights.append(float(conf))
                continue
            except (ValueError, TypeError):
                pass
        evidence = _get_first(row, ["evidence_type"]).lower()
        weights.append(tier.get(evidence, 0.5))
    return np.asarray(weights, dtype=float)


def build_sparse_operator(
    concept_nodes: pd.DataFrame,
    edges: pd.DataFrame,
    propagating_only: bool = True,
) -> Tuple[sparse.csr_matrix, Dict[str, int], Dict[int, str]]:
    work = edges.copy()
    if work.empty or "source" not in work.columns:
        node_ids = [clean_node_id(str(nid)) for nid in concept_nodes["node_id"].tolist()]
        n = len(node_ids)
        node_to_idx = {nid: i for i, nid in enumerate(node_ids)}
        idx_to_node = {i: nid for nid, i in node_to_idx.items()}
        log.warning("No edges supplied to build_sparse_operator — returning empty %dx%d operator.", n, n)
        return sparse.csr_matrix((n, n), dtype=float), node_to_idx, idx_to_node
    if propagating_only and "propagating" in work.columns:
        work = work[work["propagating"] == True].copy()

    node_ids = [clean_node_id(str(nid)) for nid in concept_nodes["node_id"].tolist()]
    n = len(node_ids)
    node_to_idx = {nid: i for i, nid in enumerate(node_ids)}
    idx_to_node = {i: nid for nid, i in node_to_idx.items()}

    work["source"] = work["source"].map(lambda x: clean_node_id(str(x)))
    work["target"] = work["target"].map(lambda x: clean_node_id(str(x)))
    mask = work["source"].isin(node_to_idx) & work["target"].isin(node_to_idx)
    work = work[mask].copy()

    if work.empty:
        log.warning("No usable edges — returning empty %dx%d operator.", n, n)
        return sparse.csr_matrix((n, n), dtype=float), node_to_idx, idx_to_node

    rows = work["source"].map(lambda x: node_to_idx.get(clean_node_id(str(x)))).to_numpy()
    cols = work["target"].map(lambda x: node_to_idx.get(clean_node_id(str(x)))).to_numpy()
    # NOTE: All analyses reported in the manuscript use uniform edge weights.
    # DDKG exports do not carry confidence or evidence_type fields; edge_weights()
    # returns 0.5 for all edges in this case, which after row normalization
    # produces identical results to unit weights.
    vals = edge_weights(work)
    return (
        sparse.csr_matrix((vals, (rows, cols)), shape=(n, n), dtype=float),
        node_to_idx,
        idx_to_node,
    )


# ---------------------------------------------------------------------------
# Propagation
# ---------------------------------------------------------------------------

def row_normalize(A: sparse.csr_matrix) -> sparse.csr_matrix:
    row_sums = np.asarray(A.sum(axis=1)).flatten()
    inv = np.zeros_like(row_sums, dtype=float)
    inv[row_sums > 0] = 1.0 / row_sums[row_sums > 0]
    return sparse.diags(inv) @ A


def personalized_pagerank(
    A: sparse.csr_matrix,
    seed_idx: List[int],
    alpha: float = 0.5,
    tol: float = 1e-10,
    max_iter: int = 500,
) -> np.ndarray:
    A_tilde = row_normalize(A)
    n = A.shape[0]
    s = np.zeros(n, dtype=float)
    if seed_idx:
        s[seed_idx] = 1.0 / len(seed_idx)
    f = s.copy()
    for iteration in range(max_iter):
        nxt = (1.0 - alpha) * (A_tilde.T @ f) + alpha * s
        if np.linalg.norm(nxt - f, ord=1) < tol:
            f = nxt
            log.debug("PPR converged at iteration %d", iteration)
            break
        f = nxt
    total = f.sum()
    if total > 0:
        f /= total
    return f


# ---------------------------------------------------------------------------
# Metrics
# ---------------------------------------------------------------------------

def entropy(signal: np.ndarray) -> float:
    x = signal[signal > 0]
    return float(-(x * np.log(x)).sum())


def topk_enrichment(signal: np.ndarray, positives: Set[int], k: int) -> float:
    topk = set(np.argsort(signal)[::-1][:k].tolist())
    return float(len(topk & positives) / max(k, 1))


def functional_recovery(signal: np.ndarray, positives: Set[int]) -> Dict[str, float]:
    y_true = np.zeros(len(signal), dtype=int)
    for idx in positives:
        if idx < len(y_true):
            y_true[idx] = 1
    if len(set(y_true)) < 2:
        return {"auroc": float("nan"), "auprc": float("nan")}
    return {
        "auroc": float(roc_auc_score(y_true, signal)),
        "auprc": float(average_precision_score(y_true, signal)),
    }


def localization(signal: np.ndarray, positives: Set[int]) -> float:
    valid = [i for i in positives if i < len(signal)]
    return float(signal[valid].sum()) if valid else 0.0


def summarize_signal(signal: np.ndarray, positives: Set[int]) -> Dict[str, float]:
    return {
        **functional_recovery(signal, positives),
        "localization": localization(signal, positives),
        "topk_50": topk_enrichment(signal, positives, k=50),
        "topk_100": topk_enrichment(signal, positives, k=100),
        "entropy": entropy(signal),
    }


def random_sparsify(edges: pd.DataFrame, target_n: int, seed: int = 1) -> pd.DataFrame:
    if len(edges) <= target_n:
        return edges.copy()
    return edges.sample(n=target_n, random_state=seed).copy()


# ---------------------------------------------------------------------------
# Analysis runner
# ---------------------------------------------------------------------------

def run_analysis(
    nodes_path: str,
    edges_path: str,
    mapping_path: str,
    seed_nodes_path: str,
    heldout_nodes_path: str,
    out_json: str,
    go_obo_path: Optional[str] = None,
    alpha: float = 0.5,
    mechanistic_only: bool = False,
) -> Dict[str, object]:
    """
    Run four-arm PPR analysis.

    Parameters
    ----------
    mechanistic_only : bool
        If True, build the conditioned propagation operator using only edges
        with classification == 'mechanistic'. Observational, perturbational,
        pathway-contribution, and contextual edges are excluded from PPR.
        The full kept_edges set is still written to disk; only the operator
        used for conditioned and random-sparsified arms is filtered.
        Raw and metadata-filtered arms are unaffected.
        Use this to isolate signal from strictly mechanistic flow.
    """
    nodes = load_table(nodes_path)
    edges = load_table(edges_path)
    mapping = load_mapping(mapping_path)
    go_lookup = build_go_lookup(go_obo_path)

    seed_nodes = _read_node_list(seed_nodes_path)
    heldout_nodes = _read_node_list(heldout_nodes_path)

    # --- Conditioning ---
    log.info("=== Conditioning ===")
    cond_result = condition_edges(nodes, edges, mapping, go_lookup)

    concept_nodes = nodes[nodes["label"].str.lower() == CONCEPT_LABEL].copy()

    # Raw concept-to-concept edges — used for raw operator AND sparsification pool
    # FIX v0.3.0: sparsification must sample from raw, not directionality edges
    # concept_nodes["node_id"] is already cleaned by build_sparse_operator;
    # edges["source"/"target"] were cleaned in condition_edges, so this is consistent
    raw_concept_node_ids = {clean_node_id(str(x)) for x in concept_nodes["node_id"].tolist()}
    raw_mask = (
        edges["source"].map(lambda x: clean_node_id(str(x))).isin(raw_concept_node_ids) &
        edges["target"].map(lambda x: clean_node_id(str(x))).isin(raw_concept_node_ids)
    )
    raw_concept_edges = edges[raw_mask].copy()

    # --- Build operators ---
    log.info("=== Building operators ===")

    # Conditioned edges for PPR: full set, or mechanistic-only subset
    propagating_edges = cond_result.kept_edges
    if mechanistic_only:
        if "classification" in propagating_edges.columns:
            propagating_edges = propagating_edges[
                propagating_edges["classification"] == "mechanistic"
            ].copy()
            n_mech = len(propagating_edges)
            n_total = len(cond_result.kept_edges)
            log.info(
                "Mechanistic-only filter: %d / %d kept edges retained "
                "(%.1f%%) — Observational, Perturbational, Pathway Contribution "
                "and contextual edges excluded from PPR operator.",
                n_mech, n_total, 100 * n_mech / n_total if n_total else 0,
            )
        else:
            log.warning(
                "--mechanistic-only requested but 'classification' column absent "
                "from kept_edges — using full conditioned edge set."
            )

    A_cond, node_to_idx, _ = build_sparse_operator(
        concept_nodes, propagating_edges, propagating_only=True
    )
    A_raw, _, _ = build_sparse_operator(
        concept_nodes,
        raw_concept_edges.assign(confidence="0.5", evidence_type="", propagating=True),
        propagating_only=False,
    )
    dir_edges = build_metadata_filtered_edges(raw_concept_edges, mapping)
    A_meta, _, _ = build_sparse_operator(
        concept_nodes,
        dir_edges.assign(confidence="0.5", evidence_type=""),
        propagating_only=False,
    )

    # Sparsification target = number of propagating edges in the conditioned arm
    # When mechanistic_only=True this matches the mechanistic edge count,
    # ensuring the random control is fairly matched.
    prop_mask = (
        propagating_edges.get(
            "propagating",
            pd.Series(True, index=propagating_edges.index)
        ) == True
    ) if not propagating_edges.empty else pd.Series(dtype=bool)
    target_edge_n = int(prop_mask.sum()) if len(prop_mask) else 0

    rand_edges = random_sparsify(
        raw_concept_edges.assign(confidence="0.5", evidence_type="", propagating=True),
        target_edge_n,
        seed=1,
    )
    A_rand, _, _ = build_sparse_operator(
        concept_nodes, rand_edges, propagating_only=False
    )

    # --- Seeds and held-out ---
    seed_idx = [node_to_idx[n] for n in seed_nodes if n in node_to_idx]
    heldout_idx = {node_to_idx[n] for n in heldout_nodes if n in node_to_idx}
    missing_seeds = [n for n in seed_nodes if n not in node_to_idx]
    missing_heldout = [n for n in heldout_nodes if n not in node_to_idx]

    if missing_seeds:
        log.warning("%d seed nodes not in graph: %s", len(missing_seeds), missing_seeds[:10])
    if missing_heldout:
        log.warning("%d heldout nodes not in graph: %s", len(missing_heldout), missing_heldout[:10])

    # Hard fail: silent partial overlap would invalidate paper results
    if not seed_idx:
        raise ValueError(
            "No seed nodes found in the concept graph. "
            "Check that seed node IDs match Concept node_ids in your export."
        )
    if not heldout_idx:
        raise ValueError(
            "No held-out nodes found in the concept graph. "
            "Check that heldout node IDs match Concept node_ids in your export."
        )

    # --- Propagation ---
    log.info("=== Propagating (alpha=%.2f) ===", alpha)
    f_raw = personalized_pagerank(A_raw, seed_idx, alpha=alpha)
    f_cond = personalized_pagerank(A_cond, seed_idx, alpha=alpha)
    f_meta = personalized_pagerank(A_meta, seed_idx, alpha=alpha)
    f_rand = personalized_pagerank(A_rand, seed_idx, alpha=alpha)

    # Flow class distribution — built from propagating_edges (the actual PPR operator),
    # not from cond_result.kept_edges. When --mechanistic-only is set, this correctly
    # reports only the mechanistic subset rather than the full kept-edge composition.
    flow_dist: Dict[str, int] = defaultdict(int)
    if not propagating_edges.empty and "bifo_flow" in propagating_edges.columns:
        prop_col = propagating_edges.get(
            "propagating",
            pd.Series(True, index=propagating_edges.index)
        )
        for flow, prop in zip(propagating_edges["bifo_flow"], prop_col):
            if flow and prop:
                flow_dist[flow] += 1

    results: Dict[str, object] = {
        "parameters": {
            "alpha": alpha,
            "n_seeds_provided": len(seed_nodes),
            "n_seeds_in_graph": len(seed_idx),
            "n_heldout_provided": len(heldout_nodes),
            "n_heldout_in_graph": len(heldout_idx),
            "missing_seeds": missing_seeds,
            "missing_heldout": missing_heldout,
        },
        "coverage": asdict(cond_result.coverage),
        "concept_fallthrough_count": len(cond_result.concept_fallthrough_nodes),
        "concept_fallthrough_preview": cond_result.concept_fallthrough_nodes[:100],
        "graph_stats": {
            "raw_concept_edges": int(len(raw_concept_edges)),
            "conditioned_propagating_edges": target_edge_n,
            "flow_class_distribution": dict(flow_dist),
        },
        "raw": summarize_signal(f_raw, heldout_idx),
        "metadata_filtered": summarize_signal(f_meta, heldout_idx),
        "conditioned": summarize_signal(f_cond, heldout_idx),
        "random_sparsification_control": summarize_signal(f_rand, heldout_idx),
    }

    Path(out_json).write_text(json.dumps(results, indent=2))
    log.info("Results written to %s", out_json)

    # -----------------------------------------------------------------------
    # Persist outputs required by score_pathways.py (Analysis B)
    # Written alongside out_json using the same stem.
    # -----------------------------------------------------------------------
    out_stem = str(Path(out_json).with_suffix(""))

    # Score vectors — one per arm, indexed by node_to_idx
    np.save(f"{out_stem}_scores_cond.npy", f_cond)
    np.save(f"{out_stem}_scores_raw.npy", f_raw)
    np.save(f"{out_stem}_scores_meta.npy", f_meta)
    np.save(f"{out_stem}_scores_rand.npy", f_rand)
    log.info("Score vectors written to %s_scores_*.npy", out_stem)

    # Node index — {concept_id: matrix_row_index}
    node_index_path = f"{out_stem}_node_index.json"
    with open(node_index_path, "w") as f_ni:
        json.dump(node_to_idx, f_ni)
    log.info("Node index written to %s", node_index_path)

    # Kept (conditioned) edges — for local_bg computation in score_pathways.py
    kept_edges_path = f"{out_stem}_kept_edges.csv"
    if not cond_result.kept_edges.empty:
        cond_result.kept_edges.to_csv(kept_edges_path, index=False)
        log.info("Kept edges written to %s (%d rows)",
                 kept_edges_path, len(cond_result.kept_edges))
    else:
        log.warning("kept_edges is empty — kept_edges.csv not written")

    return results


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main() -> None:
    parser = argparse.ArgumentParser(
        description="Condition a DDKG graph to BIFO and run a four-arm propagation benchmark."
    )
    parser.add_argument("--nodes",        required=True)
    parser.add_argument("--edges",        required=True)
    parser.add_argument("--mapping",      required=True)
    parser.add_argument("--seed-nodes",   required=True)
    parser.add_argument("--heldout-nodes",required=True)
    parser.add_argument("--out-json",     required=True,
                        help="Primary results JSON. Score vectors, node index, and "
                             "kept_edges are written alongside using the same path stem. "
                             "E.g. --out-json results/run.json also writes "
                             "results/run_scores_cond.npy, results/run_node_index.json, "
                             "results/run_kept_edges.csv.")
    parser.add_argument("--go-obo",  default=None,
                        help="GO OBO file for subontology disambiguation (optional)")
    parser.add_argument("--alpha",   type=float, default=0.5,
                        help="PageRank restart probability (default 0.5)")
    parser.add_argument("--mechanistic-only", action="store_true", default=False,
                        help="Build PPR operator from mechanistic edges only "
                             "(classification == 'mechanistic'). Observational, "
                             "weak-mechanistic perturbational/correlation, pathway-contribution, "
                             "and contextual edges are excluded from propagation. "
                             "Mechanistic Perturbational Effect predicates are retained. Use to isolate "
                             "strictly mechanistic biological flow signal.")
    args = parser.parse_args()

    run_analysis(
        nodes_path=args.nodes,
        edges_path=args.edges,
        mapping_path=args.mapping,
        seed_nodes_path=args.seed_nodes,
        heldout_nodes_path=args.heldout_nodes,
        out_json=args.out_json,
        go_obo_path=args.go_obo,
        alpha=args.alpha,
        mechanistic_only=args.mechanistic_only,
    )


if __name__ == "__main__":
    main()
