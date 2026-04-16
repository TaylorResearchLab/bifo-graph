#!/usr/bin/env python3
"""
score_pathways.py — Post-propagation pathway/process scoring for BIFO benchmarks.

Takes a PPR score vector (from bifo_conditioning.py) plus the exported
nodes/edges tables and produces ranked pathways with five scoring variants
and per-pathway diagnostics.

Architecture
------------
This module is downstream of conditioning and propagation.
It does NOT re-run PPR. It consumes the score vector produced by
bifo_conditioning.py and applies pathway-level aggregation.

Five scoring variants
---------------------
1. direct       — raw PPR score on the pathway Concept node itself
2. member_mean  — mean PPR over connected member gene Concepts
3. member_max   — max PPR over connected member gene Concepts
4. degree_norm  — direct / sqrt(pathway membership degree)
                  penalizes large generic pathways
5. local_bg     — direct minus mean score of immediate graph neighbors
                  the key comparison: does BIFO shift scores relative to
                  local background?

Per-pathway diagnostics
-----------------------
- contributing_seeds : top seed nodes by score (proxy for signal origin)
- raw_vs_cond_delta  : direct_cond minus direct_raw (conditioning effect)
- degree_flag        : True if pathway degree > 90th percentile (hub warning)
- in_chd_set         : True if pathway is in the supplied CHD reference set

Membership handler map
----------------------
Pathway membership predicates and SABs are fully configurable via
MembershipSource objects. Defaults are placeholders — will be updated
once the pathway/gene query results come back from DDKG.

Usage (CLI)
-----------
  python score_pathways.py \\
    --nodes nodes.csv \\
    --edges edges_raw.csv \\
    --scores-cond scores_conditioned.npy \\
    --scores-raw  scores_raw.npy \\
    --node-index  node_index.json \\
    --seed-nodes  seed_nodes.txt \\
    --chd-pathways chd_pathway_set.txt \\
    --out-csv     pathway_scores.csv \\
    --out-json    pathway_scores.json

Usage (API)
-----------
  from score_pathways import score_pathways, DEFAULT_MEMBERSHIP_SOURCES
  scored, metrics = score_pathways(nodes, edges, scores_cond, scores_raw,
                                   node_to_idx, seed_ids, chd_set)
"""

from __future__ import annotations

import argparse

def removesuffix(s, suffix):
    return s[:-len(suffix)] if suffix and s.endswith(suffix) else s

import json
import logging
import math
from collections import defaultdict
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Dict, FrozenSet, List, Optional, Set, Tuple

import numpy as np
import pandas as pd

log = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO,
                    format="%(asctime)s %(levelname)s %(message)s")


# ---------------------------------------------------------------------------
# Membership configuration
# ---------------------------------------------------------------------------

@dataclass
class MembershipSource:
    """
    Defines how one pathway source connects to gene Concept nodes.

    forward_predicates : edge stored as pathway → gene in DDKG
    reverse_predicates : edge stored as gene → pathway in DDKG

    source_type:
        'gene_member'  — discrete pathway membership (Reactome, MSigDB)
        'annotation'   — GO-style functional annotation (less specific)
    """
    pathway_sab: str
    gene_sabs: List[str]
    forward_predicates: List[str]
    reverse_predicates: List[str]
    source_type: str = 'gene_member'

    def all_predicates(self) -> List[str]:
        return self.forward_predicates + self.reverse_predicates


# ---------------------------------------------------------------------------
# DEFAULT MEMBERSHIP CONFIG
# Predicate names verified against live DDKG query results
# (Query_A_for_score_pathways.csv, Query_B_for_score_pathways.csv).
#
# KEY FINDINGS FROM QUERIES:
#
# MSIGDB → HGNC  (two real membership predicates):
#   pathway_associated_with_gene   486,642  forward (PW→GENE)  MSigDB C2
#   has_signature_gene               7,321  forward (PW→GENE)  MSigDB Hallmarks
#
#   EXCLUDED (look like MSIGDB edges but are NOT pathway membership):
#   targets_expression_of_gene     670,694  → LINCS perturbation data, not membership
#   has_marker_gene                104,751  → cell type marker genes
#   chr_band_contains_gene          22,905  → genomic position data
#
# GO → HGNC  (annotation-style, lower count):
#   process_involves_gene           14,315  forward (PW→GENE)
#   gene_plays_role_in_process      14,315  reverse (GENE→PW)
#   location_of                      2,369  forward GO CC → HGNC
#     (cellular_component annotation, scored separately as S not PW context)
#
# REACTOME, WIKIPATHWAYS, KEGG: NO direct Concept→HGNC membership edges found.
#   These sources connect to GO terms (REACTOME→GO via has_go_term, 120,275 edges)
#   but not directly to HGNC genes in this DDKG build.
#   Reactome gene membership must be accessed via GO intermediary or
#   via a separate Reactome-specific query.
# ---------------------------------------------------------------------------
DEFAULT_MEMBERSHIP_SOURCES: List[MembershipSource] = [
    MembershipSource(
        pathway_sab='MSIGDB',
        gene_sabs=['HGNC'],
        # Verified from query: both stored as PW→GENE (forward=true)
        # pathway_associated_with_gene = MSigDB C2 curated gene sets
        # has_signature_gene            = MSigDB Hallmark gene sets
        # NOTE: targets_expression_of_gene is EXCLUDED — it is LINCS
        # perturbation data (670K edges) stored under MSIGDB SAB but
        # semantically is Perturbational Effect, not pathway membership.
        forward_predicates=[
            'pathway_associated_with_gene',
            'has_signature_gene',
        ],
        reverse_predicates=[
            'inverse_pathway_associated_with_gene',
            'inverse_has_signature_gene',
        ],
        source_type='gene_member',
    ),
    MembershipSource(
        pathway_sab='GO',
        gene_sabs=['HGNC'],
        # Verified from query: process_involves_gene stored GO→HGNC (forward=true)
        # gene_plays_role_in_process stored HGNC→GO (forward=false = reverse here)
        # location_of (GO CC → HGNC, 2,369) excluded — CC annotation, not BP membership
        # associated_with (789) excluded — too generic
        # causally_influenced_by (220) excluded — low count, UMLS artifact
        # ro, aq, qb excluded — UMLS administrative edges
        forward_predicates=[
            'process_involves_gene',
        ],
        reverse_predicates=[
            'gene_plays_role_in_process',
        ],
        source_type='annotation',  # GO is annotation-style, not discrete membership
    ),
    # REACTOME: No direct REACTOME→HGNC edges in this build.
    # Reactome connects to GO via has_go_term (120K edges REACTOME→GO).
    # Gene membership for Reactome pathways can only be obtained by traversing
    # REACTOME→GO→HGNC two-hop paths, which is NOT handled here.
    # Stub provided for future implementation when direct edges are available.
    MembershipSource(
        pathway_sab='REACTOME',
        gene_sabs=['HGNC'],
        forward_predicates=[],   # No direct edges in current build
        reverse_predicates=[],
        source_type='gene_member',
    ),
    MembershipSource(
        pathway_sab='NCC_CUSTOM',
        gene_sabs=['HGNC'],
        forward_predicates=[
            'pathway_associated_with_gene',
        ],
        reverse_predicates=[
            'inverse_pathway_associated_with_gene',
        ],
        source_type='gene_member',
    ),
]

PATHWAY_SABS: FrozenSet[str] = frozenset([
    'REACTOME', 'MSIGDB', 'GO', 'WIKIPATHWAYS', 'KEGG', 'WP',
    'NCC_CUSTOM',   # custom NCC/cilia gene sets
])

GENE_SABS: FrozenSet[str] = frozenset([
    'HGNC', 'NCBIGENE', 'ENTREZ', 'NCBI',
])


# ---------------------------------------------------------------------------
# Data structures
# ---------------------------------------------------------------------------

@dataclass
class PathwayScore:
    """All scores and diagnostics for one pathway Concept node."""
    concept_id: str
    name: str
    sab: str
    # Five scoring variants
    direct: float = 0.0
    member_mean: float = 0.0
    member_max: float = 0.0
    degree_norm: float = 0.0
    local_bg: float = 0.0
    # Raw arm scores for comparison
    direct_raw: float = 0.0
    raw_vs_cond_delta: float = 0.0
    # Diagnostics
    member_gene_count: int = 0
    degree: int = 0
    degree_flag: bool = False
    in_chd_set: bool = False
    contributing_seeds: List[str] = field(default_factory=list)

    def to_dict(self) -> dict:
        d = asdict(self)
        d['contributing_seeds'] = '|'.join(self.contributing_seeds)
        return d


# ---------------------------------------------------------------------------
# I/O helpers
# ---------------------------------------------------------------------------

def load_table(path: str) -> pd.DataFrame:
    sep = '\t' if path.endswith('.tsv') else ','
    # encoding=utf-8-sig strips BOM that Neo4j sometimes prepends to CSV exports
    return pd.read_csv(path, sep=sep, low_memory=False, encoding='utf-8-sig')


def load_scores(path: str) -> np.ndarray:
    if path.endswith('.npy'):
        return np.load(path)
    return np.loadtxt(path, dtype=float)


def load_node_index(path: str) -> Dict[str, int]:
    with open(path) as f:
        raw = json.load(f)
    if isinstance(raw, dict):
        return {clean_node_id(str(k)): int(v) for k, v in raw.items()}
    return {clean_node_id(str(item[0])): int(item[1]) for item in raw}


def _read_node_list(path: str) -> List[str]:
    """Read a node ID list, one ID per line.
    Lines starting with '#' are skipped.
    Inline comments (anything after first whitespace or '#') are stripped.
    Only the first token on each line is returned as the ID.
    """
    ids = []
    for ln in Path(path).read_text().splitlines():
        ln = ln.strip()
        if not ln or ln.startswith('#'):
            continue
        # Strip inline comment: take only the part before '#' or first whitespace
        token = ln.split('#')[0].strip().split()[0] if ln.split('#')[0].strip() else ''
        if token:
            ids.append(token)
    return ids


# ---------------------------------------------------------------------------
# Pathway node identification
# ---------------------------------------------------------------------------

def clean_node_id(raw_id: str) -> str:
    """
    Strip DDKG schema artifacts from Concept CUI strings.

    In DDKG/UBKG, some Concept CUIs (notably MSIGDB) are stored with a
    trailing ' CUI' suffix — e.g. 'MSIGDB:M39610 CUI' instead of
    'MSIGDB:M39610'. This is a graph schema artifact, not meaningful data.
    Strip it universally at ID load time so all lookups are consistent.
    """
    return removesuffix(raw_id.strip(), ' CUI')


def identify_pathway_nodes(
    nodes: pd.DataFrame,
    membership_sources: List[MembershipSource],
) -> Dict[str, str]:
    """
    Returns {concept_id: sab} for Concept nodes belonging to a pathway SAB.

    Requires a 'sab' column in nodes.csv — no edge-based fallback.
    Node IDs are cleaned via clean_node_id() to strip DDKG schema artifacts
    (e.g. trailing ' CUI' suffix on MSIGDB Concept CUIs).
    """
    all_pw_sabs = {src.pathway_sab for src in membership_sources}
    pw_nodes: Dict[str, str] = {}

    if 'sab' not in nodes.columns:
        raise ValueError(
            "nodes.csv must have a 'sab' column for pathway node identification. "
            "Ensure your DDKG export includes Code SAB per Concept node. "
            "Edge-based SAB inference is not supported — it misclassifies endpoints."
        )

    id_col = next((c for c in ['node_id', 'id', 'CUI'] if c in nodes.columns), None)
    if id_col is None:
        raise ValueError(
            "nodes.csv must have a node identifier column ('node_id', 'id', or 'CUI')."
        )

    mask = nodes['sab'].str.upper().isin(all_pw_sabs)
    for _, row in nodes[mask].iterrows():
        pw_nodes[clean_node_id(str(row[id_col]))] = str(row['sab']).upper()

    log.info("Identified %d pathway nodes from nodes.sab column", len(pw_nodes))
    return pw_nodes


# ---------------------------------------------------------------------------
# Membership graph construction
# ---------------------------------------------------------------------------


def preflight_check(nodes: pd.DataFrame) -> Dict[str, object]:
    """
    Inspect nodes.csv for multi-SAB Concepts before running scoring.

    Returns a diagnostic dict describing:
    - n_concepts          : total unique Concept IDs
    - n_multi_sab         : Concepts with more than one SAB row
    - multi_sab_examples  : up to 10 example Concept IDs and their SABs
    - pathway_sab_counts  : how many Concepts have each pathway SAB
    - recommendation      : string advice

    Run this before score_pathways() on any new export to verify that
    source-aware membership scoring will behave as expected.
    """
    id_col = next((c for c in ['node_id', 'id', 'CUI'] if c in nodes.columns), None)
    if id_col is None or 'sab' not in nodes.columns:
        return {'error': "nodes.csv missing id or sab column"}

    id_series = nodes[id_col].astype(str).map(clean_node_id)
    sab_series = nodes['sab'].str.upper()
    n_concepts = id_series.nunique()

    # Group SABs per Concept
    from itertools import groupby
    sab_by_concept: Dict[str, Set[str]] = defaultdict(set)
    for nid, sab in zip(id_series, sab_series):
        if sab:
            sab_by_concept[nid].add(sab)

    multi_sab = {nid: sabs for nid, sabs in sab_by_concept.items() if len(sabs) > 1}
    n_multi = len(multi_sab)

    # Count pathway SAB coverage
    pathway_sab_counts = (
        sab_series[sab_series.isin(PATHWAY_SABS)]
        .value_counts()
        .to_dict()
    )

    examples = [
        {'concept_id': nid, 'sabs': sorted(sabs)}
        for nid, sabs in list(multi_sab.items())[:10]
    ]

    if n_multi == 0:
        recommendation = (
            "nodes.csv has exactly one SAB per Concept. "
            "Source-aware membership scoring will be reliable."
        )
    else:
        fraction = n_multi / n_concepts
        if fraction < 0.01:
            recommendation = (
                f"{n_multi} Concepts ({100*fraction:.1f}%) have multiple SABs. "
                "Impact on membership scoring is likely small but should be verified "
                "for pathway SABs specifically."
            )
        else:
            recommendation = (
                f"WARNING: {n_multi} Concepts ({100*fraction:.1f}%) have multiple SABs. "
                "Source-aware membership scoring may be unreliable. "
                "Consider exporting nodes.csv with one representative SAB per Concept, "
                "or implementing priority-based SAB resolution before scoring."
            )

    report = {
        'n_concepts': n_concepts,
        'n_rows': len(nodes),
        'n_multi_sab_concepts': n_multi,
        'multi_sab_fraction': round(n_multi / n_concepts, 4) if n_concepts else 0,
        'multi_sab_examples': examples,
        'pathway_sab_counts': pathway_sab_counts,
        'recommendation': recommendation,
    }

    log.info("Preflight check: %d concepts, %d multi-SAB, pathway SABs: %s",
             n_concepts, n_multi, pathway_sab_counts)
    return report


def build_membership_map(
    edges: pd.DataFrame,
    node_sab_lookup: Dict[str, str],
    membership_sources: List[MembershipSource],
) -> Dict[str, List[str]]:
    """
    Returns {pathway_concept_id: [member_gene_concept_ids]}
    built from explicit typed Concept-to-Concept edges.

    Enforces MembershipSource constraints on BOTH endpoints:
      - source (or target for reverse edges) SAB must match source.pathway_sab
      - target (or source for reverse edges) SAB must be in source.gene_sabs
      - predicate must be in the source-specific forward or reverse set

    Previously this function pooled all predicates and only checked whether
    endpoints were in broad pathway/gene ID sets — source.gene_sabs was defined
    but never used. This version enforces the full triple constraint.
    """
    pred_col = next(
        (c for c in ['predicate', 'relation', 'type', 'edge_type']
         if c in edges.columns), None
    )
    if pred_col is None:
        log.warning("No predicate column in edges — membership map will be empty")
        return {}

    membership_sets: Dict[str, Set[str]] = defaultdict(set)

    for src in membership_sources:
        if not src.forward_predicates and not src.reverse_predicates:
            log.debug("Skipping %s — no predicates configured", src.pathway_sab)
            continue

        gene_sabs_upper = {s.upper() for s in src.gene_sabs}
        forward_set = set(src.forward_predicates)
        reverse_set = set(src.reverse_predicates)

        for _, row in edges.iterrows():
            a = clean_node_id(str(row['source']))
            b = clean_node_id(str(row['target']))
            pred = str(row[pred_col])

            # Forward: edge stored as pathway → gene
            if pred in forward_set:
                a_sab = node_sab_lookup.get(a, '')
                b_sab = node_sab_lookup.get(b, '')
                if a_sab == src.pathway_sab and b_sab in gene_sabs_upper:
                    membership_sets[a].add(b)

            # Reverse: edge stored as gene → pathway
            if pred in reverse_set:
                a_sab = node_sab_lookup.get(a, '')
                b_sab = node_sab_lookup.get(b, '')
                if b_sab == src.pathway_sab and a_sab in gene_sabs_upper:
                    membership_sets[b].add(a)

    # Convert sets → sorted lists (deterministic order, no duplicates)
    membership: Dict[str, List[str]] = {
        k: sorted(v) for k, v in membership_sets.items()
    }
    n_with_members = sum(1 for v in membership.values() if v)
    log.info("Membership map: %d pathways with ≥1 member gene "
             "(source-SAB constraints enforced, duplicates removed)", n_with_members)
    return membership


# ---------------------------------------------------------------------------
# Neighbor map (for local_bg score)
# ---------------------------------------------------------------------------

def build_neighbor_map(
    conditioned_edges: pd.DataFrame,
    node_to_idx: Dict[str, int],
) -> Dict[str, List[str]]:
    """
    Undirected adjacency built from the CONDITIONED edge set only.

    local_bg compares a pathway's PPR score to its neighbors in the
    biologically admissible flow graph. Using raw edges (which include
    ontology structure, code links, and non-flow edges) would produce a
    semantically muddy background — a node's neighborhood would include
    ontology hierarchy parents, metadata edges, and non-propagating context
    edges, making the local background biologically uninterpretable.

    Caller should pass the kept_edges output from bifo_conditioning.py,
    not edges_raw.csv.
    """
    neighbor_map: Dict[str, List[str]] = defaultdict(list)
    for _, row in conditioned_edges.iterrows():
        s = clean_node_id(str(row['source']))
        t = clean_node_id(str(row['target']))
        if s in node_to_idx and t in node_to_idx:
            neighbor_map[s].append(t)
            neighbor_map[t].append(s)
    return dict(neighbor_map)


# ---------------------------------------------------------------------------
# Scoring engine
# ---------------------------------------------------------------------------

def compute_pathway_scores(
    pathway_node_map: Dict[str, str],
    membership_map: Dict[str, List[str]],
    scores_cond: np.ndarray,
    scores_raw: np.ndarray,
    node_to_idx: Dict[str, int],
    nodes: pd.DataFrame,
    conditioned_edges: pd.DataFrame,
    seed_node_ids: List[str],
    chd_pathway_set: Set[str],
    membership_sources: List[MembershipSource],
) -> List[PathwayScore]:

    # Name lookup
    name_col = next(
        (c for c in ['name', 'pref_term', 'label', 'preferred_name']
         if c in nodes.columns), None
    )
    id_col = next((c for c in ['node_id', 'id', 'CUI'] if c in nodes.columns), None)
    name_lookup: Dict[str, str] = {}
    if name_col and id_col:
        name_lookup = dict(zip(
            nodes[id_col].astype(str).map(clean_node_id),
            nodes[name_col].astype(str)
        ))

    # Membership degree per pathway (for degree_norm and degree_flag)
    pathway_degrees = {pid: len(m) for pid, m in membership_map.items()}
    all_degrees = list(pathway_degrees.values()) or [0]
    p90_degree = float(np.percentile(all_degrees, 90)) if len(all_degrees) >= 10 else float('inf')

    # Neighbor map built from CONDITIONED edges only (fix 3)
    neighbor_map = build_neighbor_map(conditioned_edges, node_to_idx)

    # Top seed nodes by conditioned PPR score.
    # NOTE: These are global top seeds, not per-pathway attribution.
    # Stored as top_global_seeds (fix 5) — labeled honestly.
    # True per-pathway attribution would require path-tracing through the
    # PPR transition matrix, which is not implemented here.
    global_top_seeds: List[str] = [
        sid for sid, _ in sorted(
            [(s, float(scores_cond[node_to_idx[s]]))
             for s in seed_node_ids if s in node_to_idx],
            key=lambda x: -x[1]
        )[:5]
    ]

    results: List[PathwayScore] = []

    for concept_id, sab in pathway_node_map.items():
        if concept_id not in node_to_idx:
            continue

        idx = node_to_idx[concept_id]
        direct = float(scores_cond[idx])
        direct_raw = float(scores_raw[idx])
        members = membership_map.get(concept_id, [])

        # member_mean, member_max
        member_scores = [
            float(scores_cond[node_to_idx[m]])
            for m in members if m in node_to_idx
        ]
        member_mean = float(np.mean(member_scores)) if member_scores else 0.0
        member_max = float(np.max(member_scores)) if member_scores else 0.0

        # degree_norm: penalise large pathways
        degree = pathway_degrees.get(concept_id, 0)
        degree_norm = direct / math.sqrt(degree) if degree > 0 else direct

        # local_bg: direct minus mean conditioned-graph neighbor score
        neighbors = neighbor_map.get(concept_id, [])
        neighbor_scores = [
            float(scores_cond[node_to_idx[n]])
            for n in neighbors if n in node_to_idx
        ]
        bg_mean = float(np.mean(neighbor_scores)) if neighbor_scores else 0.0
        local_bg = direct - bg_mean

        results.append(PathwayScore(
            concept_id=concept_id,
            name=name_lookup.get(concept_id, concept_id),
            sab=sab,
            direct=direct,
            member_mean=member_mean,
            member_max=member_max,
            degree_norm=degree_norm,
            local_bg=local_bg,
            direct_raw=direct_raw,
            raw_vs_cond_delta=direct - direct_raw,
            member_gene_count=len(members),
            degree=degree,
            degree_flag=degree > p90_degree,
            in_chd_set=concept_id in chd_pathway_set,
            contributing_seeds=global_top_seeds[:],  # see top_global_seeds note above
        ))

    log.info("Scored %d pathway nodes (%d with member genes)",
             len(results), sum(1 for r in results if r.member_gene_count > 0))
    return results


# ---------------------------------------------------------------------------
# Evaluation metrics (Analysis B)
# ---------------------------------------------------------------------------

def evaluate_pathway_recovery(
    scored: List[PathwayScore],
    score_variant: str = 'degree_norm',
    k_values: Tuple[int, ...] = (10, 20, 50),
) -> Dict[str, float]:
    """
    Evaluate pathway recovery against the CHD reference set.

    Metrics
    -------
    top_k_precision        : fraction of top-k in CHD set
    top_k_enrichment_ratio : top-k CHD rate / background CHD rate
    mean_rank_chd          : mean rank of CHD pathways (conditioned)
    mean_rank_chd_raw      : mean rank of CHD pathways (raw arm)
    rank_improvement       : mean_rank_raw - mean_rank_cond (positive = better)
    """
    if not scored:
        return {}

    ranked = sorted(scored, key=lambda x: -getattr(x, score_variant, 0.0))
    ranked_raw = sorted(scored, key=lambda x: -x.direct_raw)

    n_total = len(ranked)
    n_chd = sum(1 for p in ranked if p.in_chd_set)
    bg_rate = n_chd / n_total if n_total > 0 else 0.0

    metrics: Dict[str, float] = {
        'n_pathways_scored': float(n_total),
        'n_chd_in_reference': float(n_chd),
        'background_chd_rate': bg_rate,
    }

    for k in k_values:
        top_k = ranked[:k]
        top_k_chd = sum(1 for p in top_k if p.in_chd_set)
        precision = top_k_chd / k
        enrichment = precision / bg_rate if bg_rate > 0 else float('nan')
        metrics[f'top{k}_precision'] = precision
        metrics[f'top{k}_enrichment_ratio'] = enrichment

    chd_ranks_cond = [i + 1 for i, p in enumerate(ranked) if p.in_chd_set]
    chd_ranks_raw = [i + 1 for i, p in enumerate(ranked_raw) if p.in_chd_set]

    metrics['mean_rank_chd_cond'] = (
        float(np.mean(chd_ranks_cond)) if chd_ranks_cond else float('nan')
    )
    metrics['mean_rank_chd_raw'] = (
        float(np.mean(chd_ranks_raw)) if chd_ranks_raw else float('nan')
    )
    if chd_ranks_cond and chd_ranks_raw:
        metrics['rank_improvement_cond_vs_raw'] = (
            metrics['mean_rank_chd_raw'] - metrics['mean_rank_chd_cond']
        )

    # Hub contamination: fraction of top-k that are flagged as hubs
    for k in k_values:
        top_k = ranked[:k]
        metrics[f'top{k}_hub_fraction'] = (
            sum(1 for p in top_k if p.degree_flag) / k
        )

    return metrics


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def score_pathways(
    nodes: pd.DataFrame,
    edges_raw: pd.DataFrame,
    conditioned_edges: pd.DataFrame,
    scores_cond: np.ndarray,
    scores_raw: np.ndarray,
    node_to_idx: Dict[str, int],
    seed_node_ids: List[str],
    chd_pathway_set: Set[str],
    membership_sources: Optional[List[MembershipSource]] = None,
    score_variant: str = 'degree_norm',
    min_members: int = 1,
    max_members: Optional[int] = None,
    excluded_name_patterns: Optional[List[str]] = None,
) -> Tuple[List[PathwayScore], Dict[str, float]]:
    """
    Score all pathway nodes and return (ranked_scores, metrics).

    Parameters
    ----------
    nodes                   : nodes.csv DataFrame (must have sab column)
    edges_raw               : edges_raw.csv -- for membership traversal
    conditioned_edges       : kept_edges.csv -- for local_bg
    scores_cond             : PPR vector from conditioned arm
    scores_raw              : PPR vector from raw arm
    node_to_idx             : {concept_id: matrix_index}
    seed_node_ids           : seed Concept IDs (for diagnostics)
    chd_pathway_set         : reference Concept IDs for evaluation
    membership_sources      : override DEFAULT_MEMBERSHIP_SOURCES if provided
    score_variant           : score variant for metrics ranking
    min_members             : exclude pathways with fewer mapped member genes
                              (default 1)
    max_members             : exclude pathways with more mapped member genes
                              (default None = no upper limit). Use to remove
                              giant hub pathways that dilute signal.
    excluded_name_patterns  : uppercase substrings; pathways whose names
                              contain any of these are excluded. Default
                              removes MSigDB C3 motif and miRNA sets.
    """
    sources = membership_sources or DEFAULT_MEMBERSHIP_SOURCES

    # Default exclusion patterns -- MSigDB C3 motif and miRNA sets
    if excluded_name_patterns is None:
        excluded_name_patterns = ['_Q2', '_Q3', '_Q4', '_Q5', '_Q6', 'MIR']

    pathway_node_map = identify_pathway_nodes(nodes, sources)
    if not pathway_node_map:
        log.error("No pathway nodes identified -- check nodes.csv sab column")
        return [], {}

    node_sab_lookup = build_node_sab_lookup(nodes)
    membership_map = build_membership_map(edges_raw, node_sab_lookup, sources)

    scored = compute_pathway_scores(
        pathway_node_map, membership_map,
        scores_cond, scores_raw,
        node_to_idx, nodes, conditioned_edges,
        seed_node_ids, chd_pathway_set, sources,
    )

    # Filter 1: minimum member gene count
    n_before = len(scored)
    if min_members > 0:
        scored = [p for p in scored if p.member_gene_count >= min_members]
        log.info("min_members filter (>=%d): %d -> %d pathways",
                 min_members, n_before, len(scored))

    # Filter 2: maximum member gene count
    if max_members is not None:
        n_before2 = len(scored)
        scored = [p for p in scored if p.member_gene_count <= max_members]
        log.info("max_members filter (<=%d): %d -> %d pathways",
                 max_members, n_before2, len(scored))

    # Filter 3: excluded name patterns
    if excluded_name_patterns:
        n_before2 = len(scored)
        scored = [p for p in scored
                  if not any(pat.upper() in p.name.upper()
                             for pat in excluded_name_patterns)]
        log.info("name pattern filter: %d -> %d pathways (excluded %d sets)",
                 n_before2, len(scored), n_before2 - len(scored))

    metrics = evaluate_pathway_recovery(scored, score_variant=score_variant)
    return scored, metrics



# ---------------------------------------------------------------------------
# SAB priority for multi-SAB concept resolution.
# When a Concept has multiple SAB rows, the SAB earliest in this list wins.
# Configurable — pass a custom list to build_node_sab_lookup().
DEFAULT_SAB_PRIORITY: List[str] = [
    'MSIGDB', 'REACTOME', 'GO', 'WIKIPATHWAYS', 'KEGG', 'WP',  # pathway
    'HGNC', 'NCBIGENE', 'ENTREZ', 'NCBI',                       # gene
    'UNIPROTKB', 'PR',                                           # protein
    'MONDO', 'OMIM', 'DOID',                                     # disease
]


def build_node_sab_lookup(
    nodes: pd.DataFrame,
    sab_priority: Optional[List[str]] = None,
) -> Dict[str, str]:
    """
    Build {concept_id: sab} lookup from nodes table.

    When a Concept appears with multiple SABs (multiple rows in nodes.csv),
    resolves to the SAB that appears earliest in sab_priority rather than
    using first-occurrence (export-order dependent) or last-write-wins.

    Falls back to first-occurrence only for SABs not in sab_priority.
    """
    id_col = next((c for c in ['node_id', 'id', 'CUI'] if c in nodes.columns), None)
    if id_col is None or 'sab' not in nodes.columns:
        log.warning("Cannot build node SAB lookup -- nodes.csv missing id or sab column")
        return {}

    priority = sab_priority if sab_priority is not None else DEFAULT_SAB_PRIORITY
    priority_rank: Dict[str, int] = {s.upper(): i for i, s in enumerate(priority)}
    UNRANKED = len(priority)

    id_series = nodes[id_col].astype(str).map(clean_node_id)
    sab_series = nodes['sab'].str.upper()

    n_multi = int((id_series.value_counts() > 1).sum())
    if n_multi == 0:
        return dict(zip(id_series, sab_series))

    log.info(
        "%d Concepts have multiple SAB rows -- resolving by priority "
        "(MSIGDB>REACTOME>GO>HGNC>...). Run --preflight to inspect.", n_multi
    )

    best: Dict[str, tuple] = {}
    for nid, sab in zip(id_series, sab_series):
        if not sab:
            continue
        rank = priority_rank.get(sab, UNRANKED)
        if nid not in best or rank < best[nid][1]:
            best[nid] = (sab, rank)

    return {nid: sab for nid, (sab, _) in best.items()}


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main() -> None:
    parser = argparse.ArgumentParser(
        description="Score pathway nodes from a BIFO propagation result (Analysis B).",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  Preflight check only (--nodes is the only required arg):
    python score_pathways.py --nodes nodes.csv --preflight

  Full scoring run:
    python score_pathways.py --nodes nodes.csv --edges-raw edges_raw.csv \
      --edges-conditioned kept_edges.csv --scores-cond scores_cond.npy \
      --scores-raw scores_raw.npy --node-index node_index.json \
      --seed-nodes seed_nodes.txt --out-csv pathway_scores.csv \
      --out-json pathway_scores.json
        """
    )

    # --nodes required in both modes
    parser.add_argument("--nodes",             required=True,
                        help="nodes.csv from DDKG export (must have sab column)")
    # --preflight makes all scoring arguments optional
    parser.add_argument("--preflight",         action="store_true",
                        help="Inspect nodes.csv and exit. Only --nodes required.")
    # Scoring arguments -- optional at parse time, validated below
    parser.add_argument("--edges-raw",         default=None,
                        help="edges_raw.csv -- for membership traversal")
    parser.add_argument("--edges-conditioned", default=None,
                        help="kept_edges.csv from bifo_conditioning.py -- for local_bg")
    parser.add_argument("--scores-cond",       default=None,
                        help="PPR scores conditioned arm (.npy or txt)")
    parser.add_argument("--scores-raw",        default=None,
                        help="PPR scores raw arm (.npy or txt)")
    parser.add_argument("--node-index",        default=None,
                        help="node_id->index JSON from bifo_conditioning.py")
    parser.add_argument("--seed-nodes",        default=None,
                        help="Seed concept IDs (one per line)")
    parser.add_argument("--chd-pathways",      default=None,
                        help="CHD reference pathway concept IDs (one per line)")
    parser.add_argument("--out-csv",           default=None,
                        help="Output ranked pathway scores CSV")
    parser.add_argument("--out-json",          default=None,
                        help="Output metrics JSON")
    parser.add_argument("--score-variant",     default="degree_norm",
                        choices=["direct","member_mean","member_max","degree_norm","local_bg"],
                        help="Score variant for ranking AND metrics (default: degree_norm)")
    parser.add_argument("--min-members",       type=int, default=1,
                        help="Exclude pathways with fewer than this many mapped member "
                             "genes (default 1 -- removes zero-member motif/miRNA sets)")
    parser.add_argument("--max-members",       type=int, default=None,
                        help="Exclude pathways with more than this many mapped member "
                             "genes (default: no limit). Use to remove giant hub "
                             "pathways that dilute signal (e.g. --max-members 300)")

    args = parser.parse_args()
    nodes = load_table(args.nodes)

    # --- Preflight mode ---
    if args.preflight:
        report = preflight_check(nodes)
        print(json.dumps(report, indent=2))
        raise SystemExit(0)

    # --- Scoring mode: validate required args now ---
    required_for_scoring = {
        "--edges-raw":         args.edges_raw,
        "--edges-conditioned": args.edges_conditioned,
        "--scores-cond":       args.scores_cond,
        "--scores-raw":        args.scores_raw,
        "--node-index":        args.node_index,
        "--seed-nodes":        args.seed_nodes,
        "--out-csv":           args.out_csv,
        "--out-json":          args.out_json,
    }
    missing = [k for k, v in required_for_scoring.items() if v is None]
    if missing:
        parser.error(
            "Scoring mode requires: " + ", ".join(missing) +
            "\nUse --preflight with only --nodes to inspect your export."
        )

    edges_raw         = load_table(args.edges_raw)
    conditioned_edges = load_table(args.edges_conditioned)
    scores_cond       = load_scores(args.scores_cond)
    scores_raw_arr    = load_scores(args.scores_raw)
    node_to_idx       = load_node_index(args.node_index)
    seed_ids          = _read_node_list(args.seed_nodes)
    chd_set: Set[str] = set()
    if args.chd_pathways and Path(args.chd_pathways).exists():
        chd_set = set(_read_node_list(args.chd_pathways))
        log.info("Loaded %d CHD reference pathway IDs", len(chd_set))

    scored, metrics = score_pathways(
        nodes=nodes,
        edges_raw=edges_raw,
        conditioned_edges=conditioned_edges,
        scores_cond=scores_cond,
        scores_raw=scores_raw_arr,
        node_to_idx=node_to_idx,
        seed_node_ids=seed_ids,
        chd_pathway_set=chd_set,
        score_variant=args.score_variant,
        min_members=args.min_members,
        max_members=args.max_members,
    )

    if not scored:
        log.error("No pathways scored -- check inputs")
        raise SystemExit(1)

    ranked = sorted(scored, key=lambda x: -getattr(x, args.score_variant, 0.0))

    df = pd.DataFrame([p.to_dict() for p in ranked])
    df.to_csv(args.out_csv, index=False)
    log.info("Wrote %d pathway scores to %s", len(df), args.out_csv)

    output = {
        "parameters": {"score_variant": args.score_variant,
                        "n_seeds": len(seed_ids),
                        "n_chd_reference": len(chd_set)},
        "metrics": metrics,
        "top_50": [p.to_dict() for p in ranked[:50]],
    }
    Path(args.out_json).write_text(json.dumps(output, indent=2))
    log.info("Wrote metrics to %s", args.out_json)

    print(f"\n=== Pathway Scoring Summary (score_variant={args.score_variant}) ===")
    print(f"  Pathways scored:          {len(ranked)}")
    print(f"  CHD pathways (reference): {metrics.get('n_chd_in_reference', 0):.0f}")
    print(f"  Background CHD rate:      {metrics.get('background_chd_rate', 0):.3f}")
    print(f"  Top-10 precision:         {metrics.get('top10_precision', 0):.3f}")
    print(f"  Top-10 enrichment:        {metrics.get('top10_enrichment_ratio', float('nan')):.2f}x")
    print(f"  Top-20 precision:         {metrics.get('top20_precision', 0):.3f}")
    print(f"  Mean rank CHD (cond):     {metrics.get('mean_rank_chd_cond', float('nan')):.1f}")
    print(f"  Mean rank CHD (raw):      {metrics.get('mean_rank_chd_raw', float('nan')):.1f}")
    delta = metrics.get('rank_improvement_cond_vs_raw', float('nan'))
    print(f"  Rank improvement:         {delta:+.1f} (positive = conditioning helps)")
    print(f"\n  Top 10 pathways [{args.score_variant}]:")
    for i, p in enumerate(ranked[:10], 1):
        markers = ("\u2605" if p.in_chd_set else " ") + ("H" if p.degree_flag else " ")
        print(f"    {i:2d}. [{markers}][{p.sab:10s}] {p.name[:55]:55s} "
              f"{getattr(p, args.score_variant):.5f}")
    print(f"\n  Legend: \u2605=CHD reference  H=hub (>p90 degree)")
    print(f"  NOTE: contributing_seeds = top-5 global seeds, not per-pathway attribution")


if __name__ == "__main__":
    main()
