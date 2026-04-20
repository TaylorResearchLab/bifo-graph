#!/usr/bin/env python3
"""
score_pathways.py — Post-propagation pathway/process scoring for BIFO benchmarks.

Takes a PPR score vector (from bifo_conditioning.py) plus the exported
nodes/edges tables and produces ranked pathways with five scoring variants
and per-pathway diagnostics.

Architecture
------------
This module is downstream of conditioning and propagation.
It does NOT re-run PPR for observed scores. It consumes the score vector
produced by bifo_conditioning.py and applies pathway-level aggregation.

For empirical null scoring (--n-permutations > 0), it rebuilds the
conditioned PPR operator from kept_edges.csv and runs N permutation PPR
runs with matched random seed sets drawn from eligible pathway-connected
genes. This requires no new artifacts beyond what bifo_conditioning.py
already produces.

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

Empirical null model (degree_norm only)
----------------------------------------
When --n-permutations > 0, an empirical null is computed by running PPR
N times with matched random seed sets drawn from genes that are:
  (a) present in the conditioned propagation graph
  (b) connected to at least one scored pathway via membership edges

For each pathway, the null distribution of degree_norm scores is used to
compute a finite-sample-corrected empirical p-value:
  p = (1 + count(null >= observed)) / (1 + N)
and a BH-adjusted q-value across all tested pathways.

Pathways with q < 0.05 (or user-specified threshold) are considered
statistically supported under this empirical null. This provides a
principled usability criterion: a pathway is actionable when its
propagation score exceeds what matched random seed sets typically achieve.

Per-pathway diagnostics
-----------------------
- contributing_seeds : top seed nodes by score (proxy for signal origin)
- raw_vs_cond_delta  : direct_cond minus direct_raw (conditioning effect)
- degree_flag        : True if pathway degree > 90th percentile (hub warning)
- in_chd_set         : True if pathway is in the supplied reference set
- empirical_p        : empirical p-value vs matched null (if --n-permutations > 0)
- empirical_q        : BH-adjusted q-value (if --n-permutations > 0)
- null_mean          : mean degree_norm under null (if --n-permutations > 0)
- null_sd            : std dev of degree_norm under null (if --n-permutations > 0)
- null_z             : z-score (observed - null_mean) / null_sd (if --n-permutations > 0)

Usage (CLI)
-----------
  # Observed scores only:
  python score_pathways.py \\
    --nodes nodes.csv --edges-raw edges_raw.csv \\
    --edges-conditioned kept_edges.csv \\
    --scores-cond scores_cond.npy --scores-raw scores_raw.npy \\
    --node-index node_index.json --seed-nodes seed_nodes.txt \\
    --out-csv pathway_scores.csv --out-json pathway_scores.json

  # With empirical null (1000 permutations, 4 cores):
  python score_pathways.py \\
    --nodes nodes.csv --edges-raw edges_raw.csv \\
    --edges-conditioned kept_edges.csv \\
    --scores-cond scores_cond.npy --scores-raw scores_raw.npy \\
    --node-index node_index.json --seed-nodes seed_nodes.txt \\
    --out-csv pathway_scores.csv --out-json pathway_scores.json \\
    --n-permutations 1000 --n-cores 4

Usage (API)
-----------
  from score_pathways import score_pathways, DEFAULT_MEMBERSHIP_SOURCES
  scored, metrics = score_pathways(nodes, edges_raw, conditioned_edges,
                                   scores_cond, scores_raw, node_to_idx,
                                   seed_ids, chd_set,
                                   n_permutations=1000)
"""

from __future__ import annotations

import argparse


# ---------------------------------------------------------------------------
# Parallel scoring worker
# ---------------------------------------------------------------------------
_WORKER_GLOBALS: dict = {}

def _init_worker(globals_dict: dict) -> None:
    """Initializer for Pool workers — loads shared read-only data once."""
    _WORKER_GLOBALS.update(globals_dict)


def _score_chunk(chunk: list) -> list:
    """Score a chunk of (concept_id, sab) pairs. Returns list of PathwayScore dicts."""
    g = _WORKER_GLOBALS
    scores_cond    = g['scores_cond']
    scores_raw     = g['scores_raw']
    node_to_idx    = g['node_to_idx']
    name_lookup    = g['name_lookup']
    membership_map = g['membership_map']
    neighbor_map   = g['neighbor_map']
    pathway_degrees = g['pathway_degrees']
    p90_degree     = g['p90_degree']
    chd_pathway_set = g['chd_pathway_set']
    global_top_seeds = g['global_top_seeds']

    results = []
    for concept_id, sab in chunk:
        if concept_id not in node_to_idx:
            continue
        idx        = node_to_idx[concept_id]
        direct     = float(scores_cond[idx])
        direct_raw = float(scores_raw[idx])
        members    = membership_map.get(concept_id, [])

        member_scores = [
            float(scores_cond[node_to_idx[m]])
            for m in members if m in node_to_idx
        ]
        member_mean = float(np.mean(member_scores)) if member_scores else 0.0
        member_max  = float(np.max(member_scores))  if member_scores else 0.0

        degree      = pathway_degrees.get(concept_id, 0)
        degree_norm = direct / math.sqrt(degree) if degree > 0 else direct

        neighbors = neighbor_map.get(concept_id, [])
        neighbor_scores = [
            float(scores_cond[node_to_idx[n]])
            for n in neighbors if n in node_to_idx
        ]
        bg_mean  = float(np.mean(neighbor_scores)) if neighbor_scores else 0.0
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
            contributing_seeds=global_top_seeds[:],
        ))
    return results

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
import multiprocessing as mp
import os
import pandas as pd
from scipy import sparse as sp_sparse

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
        forward_predicates=[
            'process_involves_gene',
        ],
        reverse_predicates=[
            'gene_plays_role_in_process',
        ],
        source_type='annotation',
    ),
    MembershipSource(
        pathway_sab='REACTOME',
        gene_sabs=['HGNC'],
        forward_predicates=[],
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
    'NCC_CUSTOM',
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
    # Empirical null model outputs — membership rewiring null for degree_norm
    # (None when --n-permutations not used)
    empirical_p: Optional[float] = None
    empirical_q: Optional[float] = None
    null_mean: Optional[float] = None
    null_sd: Optional[float] = None
    null_z: Optional[float] = None
    # Member-mean null outputs — stratified gene set permutation
    # (None when --n-permutations not used)
    member_mean_p: Optional[float] = None
    member_mean_q: Optional[float] = None
    member_mean_null_mean: Optional[float] = None
    member_mean_null_sd: Optional[float] = None
    member_mean_null_z: Optional[float] = None

    def to_dict(self) -> dict:
        d = asdict(self)
        d['contributing_seeds'] = '|'.join(self.contributing_seeds)
        return d


# ---------------------------------------------------------------------------
# I/O helpers
# ---------------------------------------------------------------------------

def load_table(path: str) -> pd.DataFrame:
    sep = '\t' if path.endswith('.tsv') else ','
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
    ids = []
    for ln in Path(path).read_text().splitlines():
        ln = ln.strip()
        if not ln or ln.startswith('#'):
            continue
        token = ln.split('#')[0].strip().split()[0] if ln.split('#')[0].strip() else ''
        if token:
            ids.append(token)
    return ids


# ---------------------------------------------------------------------------
# Pathway node identification
# ---------------------------------------------------------------------------

def clean_node_id(raw_id: str) -> str:
    """Strip DDKG schema artifacts from Concept CUI strings."""
    return removesuffix(raw_id.strip(), ' CUI')


def identify_pathway_nodes(
    nodes: pd.DataFrame,
    membership_sources: List[MembershipSource],
) -> Dict[str, str]:
    """Returns {concept_id: sab} for Concept nodes belonging to a pathway SAB."""
    all_pw_sabs = {src.pathway_sab for src in membership_sources}
    pw_nodes: Dict[str, str] = {}

    if 'sab' not in nodes.columns:
        raise ValueError(
            "nodes.csv must have a 'sab' column for pathway node identification."
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
    """Inspect nodes.csv for multi-SAB Concepts before running scoring."""
    id_col = next((c for c in ['node_id', 'id', 'CUI'] if c in nodes.columns), None)
    if id_col is None or 'sab' not in nodes.columns:
        return {'error': "nodes.csv missing id or sab column"}

    id_series = nodes[id_col].astype(str).map(clean_node_id)
    sab_series = nodes['sab'].str.upper()
    n_concepts = id_series.nunique()

    sab_by_concept: Dict[str, Set[str]] = defaultdict(set)
    for nid, sab in zip(id_series, sab_series):
        if sab:
            sab_by_concept[nid].add(sab)

    multi_sab = {nid: sabs for nid, sabs in sab_by_concept.items() if len(sabs) > 1}
    n_multi = len(multi_sab)

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
                "Impact on membership scoring is likely small."
            )
        else:
            recommendation = (
                f"WARNING: {n_multi} Concepts ({100*fraction:.1f}%) have multiple SABs. "
                "Source-aware membership scoring may be unreliable."
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
    """Returns {pathway_concept_id: [member_gene_concept_ids]}."""
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
            continue

        gene_sabs_upper = {s.upper() for s in src.gene_sabs}
        forward_set = set(src.forward_predicates)
        reverse_set = set(src.reverse_predicates)

        for _, row in edges.iterrows():
            a = clean_node_id(str(row['source']))
            b = clean_node_id(str(row['target']))
            pred = str(row[pred_col])

            if pred in forward_set:
                a_sab = node_sab_lookup.get(a, '')
                b_sab = node_sab_lookup.get(b, '')
                if a_sab == src.pathway_sab and b_sab in gene_sabs_upper:
                    membership_sets[a].add(b)

            if pred in reverse_set:
                a_sab = node_sab_lookup.get(a, '')
                b_sab = node_sab_lookup.get(b, '')
                if b_sab == src.pathway_sab and a_sab in gene_sabs_upper:
                    membership_sets[b].add(a)

    membership: Dict[str, List[str]] = {
        k: sorted(v) for k, v in membership_sets.items()
    }
    n_with_members = sum(1 for v in membership.values() if v)
    log.info("Membership map: %d pathways with ≥1 member gene", n_with_members)
    return membership


# ---------------------------------------------------------------------------
# Neighbor map (for local_bg score)
# ---------------------------------------------------------------------------

def build_neighbor_map(
    conditioned_edges: pd.DataFrame,
    node_to_idx: Dict[str, int],
) -> Dict[str, List[str]]:
    """Undirected adjacency from conditioned edge set only."""
    neighbor_map: Dict[str, List[str]] = defaultdict(list)
    for _, row in conditioned_edges.iterrows():
        s = clean_node_id(str(row['source']))
        t = clean_node_id(str(row['target']))
        if s in node_to_idx and t in node_to_idx:
            neighbor_map[s].append(t)
            neighbor_map[t].append(s)
    return dict(neighbor_map)


# ---------------------------------------------------------------------------
# PPR infrastructure for empirical null (Option B: rebuild from conditioned edges)
# These mirror bifo_conditioning.py exactly to guarantee identical propagation.
# ---------------------------------------------------------------------------

def _row_normalize(A: sp_sparse.csr_matrix) -> sp_sparse.csr_matrix:
    """Row-normalize sparse matrix. Identical to bifo_conditioning.py."""
    row_sums = np.asarray(A.sum(axis=1)).flatten()
    inv = np.zeros_like(row_sums, dtype=float)
    inv[row_sums > 0] = 1.0 / row_sums[row_sums > 0]
    return sp_sparse.diags(inv) @ A


def _personalized_pagerank(
    A_tilde_T: sp_sparse.csr_matrix,
    seed_idx: List[int],
    n: int,
    alpha: float = 0.5,
    tol: float = 1e-10,
    max_iter: int = 500,
) -> np.ndarray:
    """
    PPR with pre-transposed row-normalized matrix.
    Identical convergence logic to bifo_conditioning.py.
    """
    s = np.zeros(n, dtype=float)
    if seed_idx:
        s[seed_idx] = 1.0 / len(seed_idx)
    f = s.copy()
    for _ in range(max_iter):
        nxt = (1.0 - alpha) * (A_tilde_T @ f) + alpha * s
        if np.linalg.norm(nxt - f, ord=1) < tol:
            f = nxt
            break
        f = nxt
    total = f.sum()
    if total > 0:
        f /= total
    return f


def _build_operator_from_edges(
    conditioned_edges: pd.DataFrame,
    node_to_idx: Dict[str, int],
) -> sp_sparse.csr_matrix:
    """
    Rebuild sparse PPR operator from conditioned kept_edges.
    Uses only propagating=True edges (same as bifo_conditioning.py).
    Uniform edge weights — matching all reported analyses.
    """
    n = len(node_to_idx)
    work = conditioned_edges.copy()

    if 'propagating' in work.columns:
        work = work[work['propagating'] == True].copy()

    work['source'] = work['source'].astype(str).map(
        lambda x: removesuffix(x.strip(), ' CUI'))
    work['target'] = work['target'].astype(str).map(
        lambda x: removesuffix(x.strip(), ' CUI'))

    mask = work['source'].isin(node_to_idx) & work['target'].isin(node_to_idx)
    work = work[mask]

    if work.empty:
        log.warning("_build_operator_from_edges: no usable edges — returning empty operator")
        return sp_sparse.csr_matrix((n, n), dtype=float)

    rows = work['source'].map(node_to_idx).to_numpy()
    cols = work['target'].map(node_to_idx).to_numpy()
    vals = np.ones(len(rows), dtype=float)  # uniform weights — matching paper

    A = sp_sparse.csr_matrix((vals, (rows, cols)), shape=(n, n), dtype=float)
    log.info("_build_operator_from_edges: %dx%d operator from %d edges", n, n, len(rows))
    return A


# ---------------------------------------------------------------------------
# Empirical null model
# ---------------------------------------------------------------------------

def build_null_seed_universe(
    membership_map: Dict[str, List[str]],
    node_to_idx: Dict[str, int],
) -> List[str]:
    """
    Return sorted list of eligible null seed gene IDs.

    Eligible genes are those that are:
      (a) present in the conditioned propagation graph (in node_to_idx)
      (b) connected to at least one scored pathway (in membership_map values)

    This is the default matched background for null seed sampling.
    Null seeds are drawn from genes that, like the real seeds, are connected
    to pathway nodes via membership edges and are propagatable in the graph.
    """
    universe: Set[str] = set()
    for members in membership_map.values():
        for gene_id in members:
            if gene_id in node_to_idx:
                universe.add(gene_id)
    result = sorted(universe)
    log.info("Null seed universe: %d eligible pathway-connected genes", len(result))
    return result


# Module-level shared state for null permutation workers
_NULL_SHARED = {}


def _null_worker_init(shared: dict) -> None:
    """Initialize null permutation worker."""
    global _NULL_SHARED
    _NULL_SHARED = shared


def _null_worker_run(perm_seed: int) -> np.ndarray:
    """
    Single null permutation run.

    Draws a matched random seed set (size = n_seeds, without replacement
    from the eligible pathway-connected gene universe), runs PPR, and
    computes degree_norm scores for all pathway nodes.

    The degree_norm formula is IDENTICAL to _score_chunk():
        degree_norm = direct / sqrt(degree) if degree > 0 else direct
    This guarantees that observed and null scores are directly comparable.

    Returns 1D array of degree_norm scores in pathway_ids order.
    """
    g = _NULL_SHARED
    A_tilde_T       = g['A_tilde_T']
    node_to_idx     = g['node_to_idx']
    n               = g['n']
    alpha           = g['alpha']
    universe        = g['universe']
    n_seeds         = g['n_seeds']
    pathway_ids     = g['pathway_ids']
    pathway_degrees = g['pathway_degrees']

    rng = np.random.default_rng(perm_seed)

    if len(universe) < n_seeds:
        chosen = list(universe)
    else:
        chosen = list(rng.choice(universe, size=n_seeds, replace=False))

    seed_idx = [node_to_idx[g_id] for g_id in chosen if g_id in node_to_idx]

    if not seed_idx:
        return np.zeros(len(pathway_ids), dtype=float)

    f = _personalized_pagerank(A_tilde_T, seed_idx, n, alpha=alpha)

    null_scores = np.zeros(len(pathway_ids), dtype=float)
    for i, pid in enumerate(pathway_ids):
        if pid not in node_to_idx:
            continue
        idx    = node_to_idx[pid]
        direct = float(f[idx])
        deg    = pathway_degrees.get(pid, 0)
        null_scores[i] = direct / math.sqrt(deg) if deg > 0 else direct

    return null_scores


def run_empirical_null(
    conditioned_edges: pd.DataFrame,
    node_to_idx: Dict[str, int],
    pathway_ids: List[str],
    pathway_degrees: Dict[str, int],
    membership_map: Dict[str, List[str]],
    n_seeds: int,
    n_permutations: int = 1000,
    alpha: float = 0.5,
    n_cores: int = 1,
    random_seed: int = 42,
) -> np.ndarray:
    """
    Run N permutation PPR runs with matched random seed sets.

    The conditioned PPR operator is rebuilt from kept_edges.csv (Option B):
    no serialized matrix artifacts required. This ensures that the null
    model uses exactly the same propagation graph as the observed scoring.

    Parameters
    ----------
    conditioned_edges : kept_edges.csv from bifo_conditioning.py
    node_to_idx       : {concept_id: matrix_index}
    pathway_ids       : ordered list of pathway concept IDs (defines row order)
    pathway_degrees   : {pathway_id: membership_degree}
    membership_map    : {pathway_id: [member_gene_ids]}
    n_seeds           : number of seeds per permutation (matches real run)
    n_permutations    : number of permutation runs (default 1000)
    alpha             : PPR restart probability (default 0.5)
    n_cores           : parallel workers (default 1)
    random_seed       : base random seed for reproducibility

    Returns
    -------
    null_matrix : np.ndarray of shape (n_pathways, n_permutations)
        Each column is one permutation's degree_norm scores in pathway_ids order.
    """
    log.info("Building conditioned operator for null model (%d permutations, alpha=%.2f)...",
             n_permutations, alpha)

    A = _build_operator_from_edges(conditioned_edges, node_to_idx)
    A_tilde = _row_normalize(A)
    A_tilde_T = A_tilde.T.tocsr()
    n = A.shape[0]

    universe = build_null_seed_universe(membership_map, node_to_idx)
    if len(universe) == 0:
        log.error("Null seed universe is empty — cannot run permutations")
        return np.zeros((len(pathway_ids), n_permutations), dtype=float)

    if len(universe) < n_seeds:
        log.warning(
            "Null seed universe (%d genes) < seed set size (%d). "
            "Null runs will use all universe genes.",
            len(universe), n_seeds
        )

    shared = {
        'A_tilde_T':       A_tilde_T,
        'node_to_idx':     node_to_idx,
        'n':               n,
        'alpha':           alpha,
        'universe':        universe,
        'n_seeds':         n_seeds,
        'pathway_ids':     pathway_ids,
        'pathway_degrees': pathway_degrees,
    }

    # Deterministic per-permutation seeds
    rng = np.random.default_rng(random_seed)
    perm_seeds = rng.integers(0, 2**31, size=n_permutations).tolist()

    cores = min(n_cores, n_permutations) if n_cores > 1 else 1

    if cores > 1:
        log.info("Running %d null permutations across %d cores...", n_permutations, cores)
        with mp.Pool(
            processes=cores,
            initializer=_null_worker_init,
            initargs=(shared,)
        ) as pool:
            columns = pool.map(_null_worker_run, perm_seeds)
    else:
        log.info("Running %d null permutations (single core)...", n_permutations)
        _null_worker_init(shared)
        columns = [_null_worker_run(ps) for ps in perm_seeds]

    null_matrix = np.column_stack(columns)
    log.info("Null model complete: matrix shape %s", null_matrix.shape)
    return null_matrix


def compute_empirical_pvalues(
    observed: np.ndarray,
    null_matrix: np.ndarray,
) -> Dict[str, np.ndarray]:
    """
    Compute per-pathway empirical p-values and BH-adjusted q-values.

    Empirical p-value uses finite-sample correction (Phipson & Smyth 2010):
        p = (1 + count(null >= observed)) / (1 + n_permutations)
    This avoids zero p-values and is unbiased under the permutation null.

    BH correction is applied across all tested pathways to control FDR.

    Parameters
    ----------
    observed    : 1D array of observed degree_norm scores (n_pathways,)
    null_matrix : 2D array of null scores (n_pathways, n_permutations)

    Returns
    -------
    dict with arrays of length n_pathways:
        empirical_p : finite-sample-corrected empirical p-value
        empirical_q : BH-adjusted FDR q-value
        null_mean   : mean of null distribution
        null_sd     : standard deviation of null distribution
        null_z      : z-score = (observed - null_mean) / null_sd
    """
    n_pathways, n_perm = null_matrix.shape

    exceed_counts = (null_matrix >= observed[:, np.newaxis]).sum(axis=1)
    empirical_p = (1.0 + exceed_counts) / (1.0 + n_perm)

    # BH FDR correction
    n = len(empirical_p)
    order = np.argsort(empirical_p)
    ranks = np.empty(n, dtype=int)
    ranks[order] = np.arange(1, n + 1)
    empirical_q_raw = np.minimum(1.0, empirical_p * n / ranks)
    # Enforce monotonicity
    empirical_q_mono = np.minimum.accumulate(empirical_q_raw[order][::-1])[::-1]
    empirical_q = np.empty(n)
    empirical_q[order] = empirical_q_mono

    null_mean = null_matrix.mean(axis=1)
    null_sd   = null_matrix.std(axis=1)

    with np.errstate(divide='ignore', invalid='ignore'):
        null_z = np.where(null_sd > 0, (observed - null_mean) / null_sd, 0.0)

    return {
        'empirical_p': empirical_p,
        'empirical_q': empirical_q,
        'null_mean':   null_mean,
        'null_sd':     null_sd,
        'null_z':      null_z,
    }


# ---------------------------------------------------------------------------
# Membership rewiring null model (degree-preserving bipartite edge swap)
# Optimized: bridge edges pre-indexed as int32 arrays to eliminate pickle
# overhead when passing to multiprocessing workers.
# ---------------------------------------------------------------------------

def _extract_bridge_edges_indexed(
    conditioned_edges: pd.DataFrame,
    node_to_idx: Dict[str, int],
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """
    Extract and pre-index bridge and non-bridge edges as int32 arrays.

    Pre-indexing eliminates string tuple pickle overhead when passing
    data to multiprocessing workers — the dominant bottleneck for large
    graphs (600K+ bridge edges).

    Returns
    -------
    bridge_src      : int32 array of gene matrix indices
    bridge_dst      : int32 array of pathway matrix indices
    non_bridge_rows : int32 array of non-bridge source indices
    non_bridge_cols : int32 array of non-bridge target indices
    """
    prop = conditioned_edges.copy()
    if 'propagating' in prop.columns:
        prop = prop[prop['propagating'] == True]

    prop['source'] = prop['source'].astype(str).map(
        lambda x: removesuffix(x.strip(), ' CUI'))
    prop['target'] = prop['target'].astype(str).map(
        lambda x: removesuffix(x.strip(), ' CUI'))

    is_bridge = (prop['bifo_flow'] == 'Pathway Contribution') \
        if 'bifo_flow' in prop.columns \
        else pd.Series(False, index=prop.index)

    bridge_df  = prop[is_bridge]
    non_bridge = prop[~is_bridge]

    # Pre-index bridge edges
    b_src_str = bridge_df['source'].values
    b_dst_str = bridge_df['target'].values
    b_mask = np.array([
        (s in node_to_idx and t in node_to_idx)
        for s, t in zip(b_src_str, b_dst_str)
    ])
    bridge_src = np.array(
        [node_to_idx[s] for s, t, m in zip(b_src_str, b_dst_str, b_mask) if m],
        dtype=np.int32
    )
    bridge_dst = np.array(
        [node_to_idx[t] for s, t, m in zip(b_src_str, b_dst_str, b_mask) if m],
        dtype=np.int32
    )

    # Pre-index non-bridge edges
    nb_mask = (non_bridge['source'].isin(node_to_idx) &
               non_bridge['target'].isin(node_to_idx))
    non_bridge = non_bridge[nb_mask]
    non_bridge_rows = non_bridge['source'].map(node_to_idx).to_numpy(dtype=np.int32)
    non_bridge_cols = non_bridge['target'].map(node_to_idx).to_numpy(dtype=np.int32)

    log.info(
        "Bridge edges indexed: %d Pathway Contribution | Non-bridge: %d",
        len(bridge_src), len(non_bridge_rows)
    )
    return bridge_src, bridge_dst, non_bridge_rows, non_bridge_cols


def _swap_integer_edges(
    bridge_src: np.ndarray,
    bridge_dst: np.ndarray,
    n_swaps: int,
    rng: np.random.Generator,
) -> Tuple[np.ndarray, np.ndarray, dict]:
    """
    Degree-preserving edge swap on integer-indexed bridge edges.

    Swaps pathway endpoints (dst) while preserving:
      - Each pathway's in-degree (member gene count)
      - Each gene's out-degree (pathway membership count)

    Uses int64 bit-packing for O(1) edge set membership checks:
        packed = (src << 32) | dst

    Each permutation call receives a fresh copy of the arrays.
    """
    src = bridge_src.copy()
    dst = bridge_dst.copy()
    n = len(src)

    packed = (src.astype(np.int64) << 32) | dst.astype(np.int64)
    edge_set = set(packed.tolist())

    completed = 0
    attempts = 0
    max_attempts = n_swaps * 20

    while completed < n_swaps and attempts < max_attempts:
        attempts += 1
        i = int(rng.integers(n))
        j = int(rng.integers(n))
        if i == j:
            continue

        g1, p1 = int(src[i]), int(dst[i])
        g2, p2 = int(src[j]), int(dst[j])
        if p1 == p2:
            continue

        new1 = (g1 << 32) | p2
        new2 = (g2 << 32) | p1
        if new1 in edge_set or new2 in edge_set:
            continue

        edge_set.discard((g1 << 32) | p1)
        edge_set.discard((g2 << 32) | p2)
        edge_set.add(new1)
        edge_set.add(new2)
        dst[i] = p2
        dst[j] = p1
        completed += 1

    return src, dst, {
        'requested_swaps': n_swaps,
        'completed_swaps': completed,
        'attempts': attempts,
        'acceptance_rate': completed / max(attempts, 1),
    }


def _rewiring_worker_init(shared: dict) -> None:
    """Initialize membership rewiring worker."""
    global _REWIRING_SHARED
    _REWIRING_SHARED = shared


_REWIRING_SHARED = {}


def _rewiring_worker_run(perm_seed: int) -> np.ndarray:
    """
    Single membership rewiring permutation (integer-indexed, optimized).

    Workers receive only numpy arrays — no string tuples — eliminating
    pickle serialization overhead for large graphs.
    """
    g = _REWIRING_SHARED
    bridge_src_orig     = g['bridge_src']
    bridge_dst_orig     = g['bridge_dst']
    non_bridge_rows     = g['non_bridge_rows']
    non_bridge_cols     = g['non_bridge_cols']
    n                   = g['n']
    alpha               = g['alpha']
    seed_idx            = g['seed_idx']
    n_swaps             = g['n_swaps']
    pathway_indices     = g['pathway_indices']
    pathway_degrees_arr = g['pathway_degrees_arr']

    rng = np.random.default_rng(perm_seed)

    rewired_src, rewired_dst, _ = _swap_integer_edges(
        bridge_src_orig, bridge_dst_orig, n_swaps, rng
    )

    all_rows = np.concatenate([non_bridge_rows, rewired_src])
    all_cols = np.concatenate([non_bridge_cols, rewired_dst])
    all_vals = np.ones(len(all_rows), dtype=float)

    A = sp_sparse.csr_matrix(
        (all_vals, (all_rows, all_cols)), shape=(n, n), dtype=float
    )
    row_sums = np.asarray(A.sum(axis=1)).flatten()
    with np.errstate(divide='ignore', invalid='ignore'):
        inv = np.where(row_sums > 0, 1.0 / row_sums, 0.0)
    A_tilde_T = (sp_sparse.diags(inv) @ A).T.tocsr()

    f = _personalized_pagerank(A_tilde_T, seed_idx, n, alpha=alpha)

    # Vectorized degree_norm computation
    pw_scores = f[pathway_indices]
    degs = pathway_degrees_arr.astype(float)
    with np.errstate(divide='ignore', invalid='ignore'):
        null_scores = np.where(degs > 0, pw_scores / np.sqrt(degs), pw_scores)

    return null_scores


def run_membership_rewiring_null(
    conditioned_edges: pd.DataFrame,
    node_to_idx: Dict[str, int],
    pathway_ids: List[str],
    pathway_degrees: Dict[str, int],
    seed_node_ids: List[str],
    n_permutations: int = 1000,
    alpha: float = 0.5,
    n_cores: int = 1,
    random_seed: int = 42,
    n_swaps_multiplier: int = 10,
) -> np.ndarray:
    """
    Run N membership rewiring permutations (optimized integer-indexed version).

    For each permutation:
      1. Apply degree-preserving edge swaps to Pathway Contribution bridge edges
         starting fresh from the original integer-indexed edge arrays
      2. Rebuild PPR operator (non-bridge edges held fixed)
      3. Rerun PPR with real seeds (unchanged)
      4. Compute degree_norm for all scored pathways

    This null tests: does this pathway's actual gene membership receive more
    propagated mass than expected under a randomized membership structure with
    the same degree constraints?

    Valid at all seed sizes because seeds are never randomized.
    Optimized: bridge edges pre-indexed as int32 arrays to minimize pickle
    overhead in multiprocessing (dominant bottleneck for large graphs).
    """
    log.info("Extracting and pre-indexing bridge edges for rewiring null...")
    bridge_src, bridge_dst, non_bridge_rows, non_bridge_cols = \
        _extract_bridge_edges_indexed(conditioned_edges, node_to_idx)

    if len(bridge_src) == 0:
        log.error("No Pathway Contribution bridge edges found")
        return np.zeros((len(pathway_ids), n_permutations), dtype=float)

    n_swaps = len(bridge_src) * n_swaps_multiplier
    n = len(node_to_idx)
    log.info("Rewiring null: %d bridge edges, %d swaps/perm, %d perms",
             len(bridge_src), n_swaps, n_permutations)

    pathway_indices = np.array(
        [node_to_idx.get(pid, 0) for pid in pathway_ids], dtype=np.int32
    )
    pathway_degrees_arr = np.array(
        [pathway_degrees.get(pid, 0) for pid in pathway_ids], dtype=np.int32
    )

    effective_seeds = [s for s in seed_node_ids if s in node_to_idx]
    if not effective_seeds:
        log.error("No seed nodes in graph — cannot run rewiring null")
        return np.zeros((len(pathway_ids), n_permutations), dtype=float)
    if len(effective_seeds) < len(seed_node_ids):
        log.warning("Rewiring null: %d effective seeds (%d absent from graph)",
                    len(effective_seeds), len(seed_node_ids) - len(effective_seeds))

    shared = {
        'bridge_src':           bridge_src,
        'bridge_dst':           bridge_dst,
        'non_bridge_rows':      non_bridge_rows,
        'non_bridge_cols':      non_bridge_cols,
        'n':                    n,
        'alpha':                alpha,
        'seed_idx':             [node_to_idx[s] for s in effective_seeds],
        'n_swaps':              n_swaps,
        'pathway_indices':      pathway_indices,
        'pathway_degrees_arr':  pathway_degrees_arr,
    }

    rng = np.random.default_rng(random_seed)
    perm_seeds = rng.integers(0, 2**31, size=n_permutations).tolist()
    cores = min(n_cores, n_permutations) if n_cores > 1 else 1

    if cores > 1:
        log.info("Running %d rewiring permutations across %d cores...",
                 n_permutations, cores)
        with mp.Pool(
            processes=cores,
            initializer=_rewiring_worker_init,
            initargs=(shared,)
        ) as pool:
            columns = pool.map(_rewiring_worker_run, perm_seeds)
    else:
        log.info("Running %d rewiring permutations (single core)...", n_permutations)
        _rewiring_worker_init(shared)
        columns = [_rewiring_worker_run(ps) for ps in perm_seeds]

    null_matrix = np.column_stack(columns)
    log.info("Membership rewiring null complete: matrix shape %s", null_matrix.shape)
    return null_matrix


# ---------------------------------------------------------------------------
# Member-mean null model (stratified gene set permutation — no PPR rerun)
# Valid at all seed sizes. Tests whether pathway member genes carry
# disproportionate propagated signal relative to degree+score matched sets.
# Parallelized: each pathway's null distribution computed independently.
# ---------------------------------------------------------------------------

def build_gene_strata(
    eligible_genes: List[str],
    node_to_idx: Dict[str, int],
    scores_cond: np.ndarray,
    conditioned_edges: pd.DataFrame,
    n_score_bins: int = 10,
    n_degree_bins: int = 10,
) -> Dict[str, List[str]]:
    """
    Bin eligible genes by structural features only for matched null sampling.

    Matching is based on:
      (1) conditioned graph out-degree (log-binned) — controls for hub genes
      (2) pathway membership count (log-binned) — controls for broadly annotated genes

    Observed propagated score is intentionally excluded from strata to avoid
    conditioning the null on the same quantity being tested (member_mean).
    This gives a structurally matched but score-naive null.
    """
    prop = conditioned_edges.copy()
    if 'propagating' in prop.columns:
        prop = prop[prop['propagating'] == True]
    prop['source'] = prop['source'].astype(str).map(
        lambda x: removesuffix(x.strip(), ' CUI'))
    degree_counts = prop['source'].value_counts()

    # Count pathway memberships per gene across all pathways
    bridge = prop[prop.get('bifo_flow', pd.Series(dtype=str)) == 'Pathway Contribution']         if 'bifo_flow' in prop.columns else prop.iloc[0:0]
    membership_counts = bridge['source'].value_counts()         if not bridge.empty else pd.Series(dtype=int)

    genes_in_idx = [g for g in eligible_genes if g in node_to_idx]
    if not genes_in_idx:
        log.warning("build_gene_strata: no eligible genes in node index")
        return {}

    gene_degrees = np.array([float(degree_counts.get(g, 1)) for g in genes_in_idx])
    gene_memberships = np.array([float(membership_counts.get(g, 1))
                                 for g in genes_in_idx])

    eps = 1e-30
    degree_bins = np.percentile(np.log10(gene_degrees + eps),
                                np.linspace(0, 100, n_degree_bins + 1))
    membership_bins = np.percentile(np.log10(gene_memberships + eps),
                                    np.linspace(0, 100, n_score_bins + 1))

    degree_bin_ids = np.digitize(np.log10(gene_degrees + eps), degree_bins[1:-1])
    membership_bin_ids = np.digitize(np.log10(gene_memberships + eps),
                                     membership_bins[1:-1])

    strata: Dict[str, List[str]] = {}
    for gene, db, mb in zip(genes_in_idx, degree_bin_ids, membership_bin_ids):
        key = f"d{int(db)}_m{int(mb)}"
        strata.setdefault(key, []).append(gene)

    log.info("Gene strata: %d genes in %d strata (degree x membership bins)",
             len(genes_in_idx), len(strata))
    return strata


# Module-level shared state for member_mean null workers
_MM_SHARED = {}


def _mm_worker_init(shared: dict) -> None:
    global _MM_SHARED
    _MM_SHARED = shared


def _mm_worker_run_perm(perm_seed: int) -> np.ndarray:
    """
    Compute ONE permutation across ALL pathways.
    Parallelizes over permutations (1000 workers) not pathways,
    giving better load balance for large pathway universes.

    For each pathway, draws a matched null gene set without replacement
    within each stratum, avoiding duplicate genes in the null set.
    Strata are structural only (degree + membership count bins).

    Returns 1D array of length n_pathways (null member_mean for this perm).
    """
    g = _MM_SHARED
    strata          = g['strata']
    gene_to_stratum = g['gene_to_stratum']
    node_to_idx     = g['node_to_idx']
    scores_cond     = g['scores_cond']
    pathway_members = g['pathway_members']  # List[List[str]]

    all_eligible = [gg for genes in strata.values() for gg in genes]
    rng = np.random.default_rng(perm_seed)

    null_scores = np.zeros(len(pathway_members))
    for i, members in enumerate(pathway_members):
        if not members:
            continue

        # Draw without replacement across the full set to avoid duplicates
        # Track already-selected genes to prevent reuse within this pathway
        selected: Set[str] = set()
        null_genes = []

        for gene in members:
            sk = gene_to_stratum.get(gene)
            # Candidates: same stratum, not the true gene, not already selected
            candidates = [gg for gg in (strata.get(sk, []) if sk else [])
                          if gg != gene and gg not in selected]
            if not candidates:
                # Fall back to any eligible gene not yet selected
                fallback = [gg for gg in all_eligible
                            if gg != gene and gg not in selected]
                if fallback:
                    chosen = str(rng.choice(fallback))
                else:
                    chosen = gene  # degenerate: reuse true gene
            else:
                chosen = str(rng.choice(candidates))

            null_genes.append(chosen)
            selected.add(chosen)

        vals = [scores_cond[node_to_idx[gg]]
                for gg in null_genes if gg in node_to_idx]
        null_scores[i] = float(np.mean(vals)) if vals else 0.0

    return null_scores


def run_member_mean_null(
    scored: List[PathwayScore],
    membership_map: Dict[str, List[str]],
    node_to_idx: Dict[str, int],
    scores_cond: np.ndarray,
    conditioned_edges: pd.DataFrame,
    seed_node_ids: List[str],
    n_permutations: int = 1000,
    n_score_bins: int = 10,
    n_degree_bins: int = 10,
    random_seed: int = 42,
    n_cores: int = 1,
) -> None:
    """
    Compute member_mean significance via stratified gene set permutation.

    Parallelizes over permutations (not pathways) for optimal load balance:
    each worker computes one full permutation across all pathways.
    With N permutations and C cores, each core handles N/C permutations.

    No PPR reruns required — valid at all seed sizes.
    Results written directly to PathwayScore.member_mean_* fields.
    """
    eligible: Set[str] = set()
    for members in membership_map.values():
        for gg in members:
            if gg in node_to_idx:
                eligible.add(gg)
    eligible_list = sorted(eligible)

    if not eligible_list:
        log.error("member_mean null: no eligible genes found")
        return

    strata = build_gene_strata(
        eligible_list, node_to_idx, scores_cond, conditioned_edges,
        n_score_bins=n_score_bins, n_degree_bins=n_degree_bins
    )
    gene_to_stratum: Dict[str, str] = {
        gg: key for key, genes in strata.items() for gg in genes
    }

    pathways_with_members = [p for p in scored if p.member_gene_count > 0]
    pathway_members = [
        [m for m in membership_map.get(p.concept_id, []) if m in node_to_idx]
        for p in pathways_with_members
    ]

    cores = min(n_cores, n_permutations) if n_cores > 1 else 1
    log.info("member_mean null: %d pathways, %d perms, %d cores",
             len(pathways_with_members), n_permutations, cores)

    # Compute observed member_mean
    observed_vals = np.array([
        float(np.mean([scores_cond[node_to_idx[m]] for m in members]))
        if members else 0.0
        for members in pathway_members
    ])

    shared = {
        'strata':          strata,
        'gene_to_stratum': gene_to_stratum,
        'node_to_idx':     node_to_idx,
        'scores_cond':     scores_cond,
        'pathway_members': pathway_members,
    }

    # Generate deterministic per-permutation seeds
    rng = np.random.default_rng(random_seed)
    perm_seeds = rng.integers(0, 2**31, size=n_permutations).tolist()

    if cores > 1:
        with mp.Pool(processes=cores,
                     initializer=_mm_worker_init,
                     initargs=(shared,)) as pool:
            null_cols = pool.map(_mm_worker_run_perm, perm_seeds)
    else:
        _mm_worker_init(shared)
        null_cols = [_mm_worker_run_perm(ps) for ps in perm_seeds]

    # null_cols: list of n_permutations arrays each of shape (n_pathways,)
    null_matrix = np.column_stack(null_cols)  # shape: (n_pathways, n_perms)

    pval_results = compute_empirical_pvalues(observed_vals, null_matrix)

    sig_05 = 0
    for i, pathway in enumerate(pathways_with_members):
        pathway.member_mean_p         = float(pval_results['empirical_p'][i])
        pathway.member_mean_q         = float(pval_results['empirical_q'][i])
        pathway.member_mean_null_mean = float(pval_results['null_mean'][i])
        pathway.member_mean_null_sd   = float(pval_results['null_sd'][i])
        pathway.member_mean_null_z    = float(pval_results['null_z'][i])
        if pathway.member_mean_q < 0.05:
            sig_05 += 1

    log.info("member_mean null complete: %d / %d pathways significant at q < 0.05",
             sig_05, len(pathways_with_members))


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
    n_cores: int = 1,
) -> List[PathwayScore]:

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

    pathway_degrees = {pid: len(m) for pid, m in membership_map.items()}
    all_degrees = list(pathway_degrees.values()) or [0]
    p90_degree = float(np.percentile(all_degrees, 90)) if len(all_degrees) >= 10 else float('inf')

    neighbor_map = build_neighbor_map(conditioned_edges, node_to_idx)

    global_top_seeds: List[str] = [
        sid for sid, _ in sorted(
            [(s, float(scores_cond[node_to_idx[s]]))
             for s in seed_node_ids if s in node_to_idx],
            key=lambda x: -x[1]
        )[:5]
    ]

    shared = {
        'scores_cond':      scores_cond,
        'scores_raw':       scores_raw,
        'node_to_idx':      node_to_idx,
        'name_lookup':      name_lookup,
        'membership_map':   membership_map,
        'neighbor_map':     neighbor_map,
        'pathway_degrees':  pathway_degrees,
        'p90_degree':       p90_degree,
        'chd_pathway_set':  chd_pathway_set,
        'global_top_seeds': global_top_seeds,
    }

    items = list(pathway_node_map.items())
    cores = min(n_cores, len(items)) if n_cores > 1 else 1

    if cores > 1:
        chunk_size = max(1, len(items) // cores)
        chunks = [items[i:i + chunk_size] for i in range(0, len(items), chunk_size)]
        log.info("Scoring %d pathways across %d cores (%d chunks)",
                 len(items), cores, len(chunks))
        with mp.Pool(processes=cores,
                     initializer=_init_worker,
                     initargs=(shared,)) as pool:
            chunk_results = pool.map(_score_chunk, chunks)
        results: List[PathwayScore] = [r for chunk in chunk_results for r in chunk]
    else:
        _init_worker(shared)
        results = _score_chunk(items)

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
    """Evaluate pathway recovery against the reference set."""
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

    for k in k_values:
        top_k = ranked[:k]
        metrics[f'top{k}_hub_fraction'] = (
            sum(1 for p in top_k if p.degree_flag) / k
        )

    # Empirical null summary (if available)
    n_sig_05 = sum(1 for p in scored if p.empirical_q is not None and p.empirical_q < 0.05)
    n_sig_10 = sum(1 for p in scored if p.empirical_q is not None and p.empirical_q < 0.10)
    if any(p.empirical_q is not None for p in scored):
        metrics['n_significant_q05'] = float(n_sig_05)
        metrics['n_significant_q10'] = float(n_sig_10)

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
    n_cores: int = 1,
    n_permutations: int = 0,
    null_type: str = 'membership-rewiring',
    ppr_alpha: float = 0.5,
    random_seed: int = 42,
    n_swaps_multiplier: int = 10,
) -> Tuple[List[PathwayScore], Dict[str, float]]:
    """
    Score all pathway nodes and return (ranked_scores, metrics).

    Parameters
    ----------
    nodes                   : nodes.csv DataFrame (must have sab column)
    edges_raw               : edges_raw.csv -- for membership traversal
    conditioned_edges       : kept_edges.csv -- for local_bg and null model
    scores_cond             : PPR vector from conditioned arm
    scores_raw              : PPR vector from raw arm
    node_to_idx             : {concept_id: matrix_index}
    seed_node_ids           : seed Concept IDs
    chd_pathway_set         : reference Concept IDs for evaluation
    membership_sources      : override DEFAULT_MEMBERSHIP_SOURCES if provided
    score_variant           : score variant for ranking and metrics
    min_members             : exclude pathways with fewer mapped member genes
    max_members             : exclude pathways with more mapped member genes
    excluded_name_patterns  : name substrings to exclude (MSigDB motif/miRNA sets)
    n_cores                 : parallel workers
    n_permutations          : number of null permutation runs (0 = disabled)
                              When > 0, computes empirical_p, empirical_q,
                              null_mean, null_sd, null_z for each pathway.
                              Null seeds drawn from eligible pathway-connected genes.
    ppr_alpha               : PPR restart probability for null runs (default 0.5)
    random_seed             : base random seed for null permutations (default 42)
    """
    sources = membership_sources or DEFAULT_MEMBERSHIP_SOURCES

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
        n_cores=n_cores,
    )

    # Filters
    n_before = len(scored)
    if min_members > 0:
        scored = [p for p in scored if p.member_gene_count >= min_members]
        log.info("min_members filter (>=%d): %d -> %d pathways",
                 min_members, n_before, len(scored))

    if max_members is not None:
        n_before2 = len(scored)
        scored = [p for p in scored if p.member_gene_count <= max_members]
        log.info("max_members filter (<=%d): %d -> %d pathways",
                 max_members, n_before2, len(scored))

    if excluded_name_patterns:
        n_before2 = len(scored)
        scored = [p for p in scored
                  if not any(pat.upper() in p.name.upper()
                             for pat in excluded_name_patterns)]
        log.info("name pattern filter: %d -> %d pathways",
                 n_before2, len(scored))

    # Empirical null model
    if n_permutations > 0 and scored:
        log.info("Running %s null model (%d permutations)...", null_type, n_permutations)

        pathway_ids_ordered = [p.concept_id for p in scored]
        pathway_degrees_map = {p.concept_id: p.degree for p in scored}
        observed_scores = np.array([p.degree_norm for p in scored])
        # Build null universe from RETAINED scored pathways only (post-filter)
        # Excluded pathways should not contribute genes to the null universe
        scored_pathway_ids = set(pathway_ids_ordered)
        membership_map_scored = {
            pid: members
            for pid, members in membership_map.items()
            if pid in scored_pathway_ids
        }

        # Use effective seed count: seeds that actually map into the graph
        effective_seed_ids = [s for s in seed_node_ids if s in node_to_idx]
        n_effective_seeds = len(effective_seed_ids)
        if n_effective_seeds == 0:
            log.error("No seed nodes found in graph — skipping null permutations")
        else:
            if n_effective_seeds < len(seed_node_ids):
                log.warning(
                    "Null model uses %d effective seeds (%d provided, %d absent from graph)",
                    n_effective_seeds, len(seed_node_ids),
                    len(seed_node_ids) - n_effective_seeds
                )

        if n_effective_seeds > 0:
            if null_type == 'membership-rewiring':
                null_matrix = run_membership_rewiring_null(
                    conditioned_edges=conditioned_edges,
                    node_to_idx=node_to_idx,
                    pathway_ids=pathway_ids_ordered,
                    pathway_degrees=pathway_degrees_map,
                    seed_node_ids=seed_node_ids,
                    n_permutations=n_permutations,
                    alpha=ppr_alpha,
                    n_cores=n_cores,
                    random_seed=random_seed,
                    n_swaps_multiplier=n_swaps_multiplier,
                )
            else:  # seed-permutation
                null_matrix = run_empirical_null(
                conditioned_edges=conditioned_edges,
                node_to_idx=node_to_idx,
                pathway_ids=pathway_ids_ordered,
                pathway_degrees=pathway_degrees_map,
                membership_map=membership_map_scored,  # scored pathways only
            n_seeds=n_effective_seeds,
            n_permutations=n_permutations,
            alpha=ppr_alpha,
            n_cores=n_cores,
                    random_seed=random_seed,
                )

            pval_results = compute_empirical_pvalues(observed_scores, null_matrix)

            for i, p in enumerate(scored):
                p.empirical_p = float(pval_results['empirical_p'][i])
                p.empirical_q = float(pval_results['empirical_q'][i])
                p.null_mean   = float(pval_results['null_mean'][i])
                p.null_sd     = float(pval_results['null_sd'][i])
                p.null_z      = float(pval_results['null_z'][i])

            n_sig = sum(1 for p in scored if p.empirical_q is not None and p.empirical_q < 0.05)
            log.info("Empirical null complete: %d / %d pathways significant at q < 0.05",
                     n_sig, len(scored))
        # Store null run metadata for reporting — will be added to metrics after evaluate
        _null_universe_size = len(build_null_seed_universe(membership_map_scored, node_to_idx))
        _n_effective_seeds = n_effective_seeds

        # Also run member_mean stratified null (no PPR reruns — valid at all seed sizes)
        log.info("Running member_mean stratified null (%d permutations)...",
                 n_permutations)
        run_member_mean_null(
            scored=scored,
            membership_map=membership_map,
            node_to_idx=node_to_idx,
            scores_cond=scores_cond,
            conditioned_edges=conditioned_edges,
            seed_node_ids=seed_node_ids,
            n_permutations=n_permutations,
            random_seed=random_seed,
            n_cores=n_cores,
        )

    metrics = evaluate_pathway_recovery(scored, score_variant=score_variant)
    # Add null model metadata to metrics if available
    if n_permutations > 0 and scored and any(p.empirical_q is not None for p in scored):
        metrics['null_universe_size'] = float(_null_universe_size)
        metrics['n_effective_seeds'] = float(_n_effective_seeds)
        metrics['n_permutations_run'] = float(n_permutations)
    if n_permutations > 0 and scored and any(p.member_mean_q is not None for p in scored):
        n_sig_mm = sum(1 for p in scored
                       if p.member_mean_q is not None and p.member_mean_q < 0.05)
        metrics['n_significant_member_mean_q05'] = float(n_sig_mm)
    return scored, metrics


# ---------------------------------------------------------------------------
# SAB priority for multi-SAB concept resolution.
# ---------------------------------------------------------------------------
DEFAULT_SAB_PRIORITY: List[str] = [
    'MSIGDB', 'REACTOME', 'GO', 'WIKIPATHWAYS', 'KEGG', 'WP',
    'HGNC', 'NCBIGENE', 'ENTREZ', 'NCBI',
    'UNIPROTKB', 'PR',
    'MONDO', 'OMIM', 'DOID',
]


def build_node_sab_lookup(
    nodes: pd.DataFrame,
    sab_priority: Optional[List[str]] = None,
) -> Dict[str, str]:
    """Build {concept_id: sab} lookup with priority-based resolution."""
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

    log.info("%d Concepts have multiple SAB rows -- resolving by priority.", n_multi)

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
  Preflight check only:
    python score_pathways.py --nodes nodes.csv --preflight

  Full scoring run (observed scores only):
    python score_pathways.py --nodes nodes.csv --edges-raw edges_raw.csv \\
      --edges-conditioned kept_edges.csv --scores-cond scores_cond.npy \\
      --scores-raw scores_raw.npy --node-index node_index.json \\
      --seed-nodes seed_nodes.txt --out-csv pathway_scores.csv \\
      --out-json pathway_scores.json

  With empirical null model (1000 permutations, 4 cores):
    python score_pathways.py --nodes nodes.csv --edges-raw edges_raw.csv \\
      --edges-conditioned kept_edges.csv --scores-cond scores_cond.npy \\
      --scores-raw scores_raw.npy --node-index node_index.json \\
      --seed-nodes seed_nodes.txt --out-csv pathway_scores.csv \\
      --out-json pathway_scores.json \\
      --n-permutations 1000 --n-cores 4
        """
    )

    parser.add_argument("--nodes",             required=True)
    parser.add_argument("--preflight",         action="store_true")
    parser.add_argument("--edges-raw",         default=None)
    parser.add_argument("--edges-conditioned", default=None)
    parser.add_argument("--scores-cond",       default=None)
    parser.add_argument("--scores-raw",        default=None)
    parser.add_argument("--node-index",        default=None)
    parser.add_argument("--seed-nodes",        default=None)
    parser.add_argument("--chd-pathways",      default=None)
    parser.add_argument("--out-csv",           default=None)
    parser.add_argument("--out-json",          default=None)
    parser.add_argument("--score-variant",     default="degree_norm",
                        choices=["direct","member_mean","member_max","degree_norm","local_bg"])
    parser.add_argument("--min-members",       type=int, default=1)
    parser.add_argument("--max-members",       type=int, default=None)
    parser.add_argument("--n-cores",           type=int, default=1,
                        help="Parallel workers (0 = auto-detect). Used for both "
                             "observed scoring and null permutations.")
    parser.add_argument("--n-permutations",    type=int, default=0,
                        help="Number of null permutation runs for empirical p-values "
                             "(default: 0 = disabled). Recommended: 1000.")
    parser.add_argument("--null-type",          default="membership-rewiring",
                        choices=["membership-rewiring", "seed-permutation"],
                        help="Null model type (default: membership-rewiring). "
                             "membership-rewiring: keeps seeds fixed, rewires pathway "
                             "membership edges via degree-preserving swaps. Valid at "
                             "all seed sizes. "
                             "seed-permutation: draws random seed sets from eligible "
                             "pathway-connected genes. Valid only for sparse seed sets "
                             "(seed-to-universe ratio < ~5%%).")
    parser.add_argument("--n-swaps-multiplier", type=int, default=10,
                        help="Swaps per permutation = n_bridge_edges * multiplier "
                             "(default: 10). Higher values = more mixing.")
    parser.add_argument("--ppr-alpha",         type=float, default=0.5,
                        help="PPR restart probability for null runs (default: 0.5). "
                             "Should match the alpha used in bifo_conditioning.py.")
    parser.add_argument("--perm-random-seed",  type=int, default=42,
                        help="Base random seed for null permutations (default: 42).")

    args = parser.parse_args()
    nodes = load_table(args.nodes)

    if args.preflight:
        report = preflight_check(nodes)
        print(json.dumps(report, indent=2))
        raise SystemExit(0)

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
        parser.error("Scoring mode requires: " + ", ".join(missing))

    edges_raw         = load_table(args.edges_raw)
    conditioned_edges = load_table(args.edges_conditioned)
    scores_cond       = load_scores(args.scores_cond)
    scores_raw_arr    = load_scores(args.scores_raw)
    node_to_idx       = load_node_index(args.node_index)
    seed_ids          = _read_node_list(args.seed_nodes)
    chd_set: Set[str] = set()
    if args.chd_pathways and Path(args.chd_pathways).exists():
        chd_set = set(_read_node_list(args.chd_pathways))
        log.info("Loaded %d reference pathway IDs", len(chd_set))

    n_cores = args.n_cores if args.n_cores > 0 else mp.cpu_count()

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
        n_cores=n_cores,
        n_permutations=args.n_permutations,
        null_type=args.null_type,
        ppr_alpha=args.ppr_alpha,
        random_seed=args.perm_random_seed,
        n_swaps_multiplier=args.n_swaps_multiplier,
    )

    if not scored:
        log.error("No pathways scored -- check inputs")
        raise SystemExit(1)

    ranked = sorted(scored, key=lambda x: -getattr(x, args.score_variant, 0.0))

    df = pd.DataFrame([p.to_dict() for p in ranked])
    df.to_csv(args.out_csv, index=False)
    log.info("Wrote %d pathway scores to %s", len(df), args.out_csv)

    output = {
        "parameters": {
            "score_variant":      args.score_variant,
            "n_seeds_provided":   len(seed_ids),
            "n_seeds_effective":  metrics.get('n_effective_seeds', len(seed_ids)),
            "n_chd_reference":    len(chd_set),
            "n_permutations":     args.n_permutations,
            "null_type":          args.null_type,
            "ppr_alpha":          args.ppr_alpha,
            "random_seed":        args.perm_random_seed,
        },
        "metrics": metrics,
        "top_50": [p.to_dict() for p in ranked[:50]],
    }
    Path(args.out_json).write_text(json.dumps(output, indent=2))
    log.info("Wrote metrics to %s", args.out_json)

    print(f"\n=== Pathway Scoring Summary (score_variant={args.score_variant}) ===")
    print(f"  Pathways scored:          {len(ranked)}")
    print(f"  Reference pathways:       {metrics.get('n_chd_in_reference', 0):.0f}")
    print(f"  Background rate:          {metrics.get('background_chd_rate', 0):.3f}")
    print(f"  Top-10 precision:         {metrics.get('top10_precision', 0):.3f}")
    print(f"  Top-10 enrichment:        {metrics.get('top10_enrichment_ratio', float('nan')):.2f}x")
    print(f"  Top-20 precision:         {metrics.get('top20_precision', 0):.3f}")
    print(f"  Mean rank (cond):         {metrics.get('mean_rank_chd_cond', float('nan')):.1f}")
    print(f"  Mean rank (raw):          {metrics.get('mean_rank_chd_raw', float('nan')):.1f}")
    delta = metrics.get('rank_improvement_cond_vs_raw', float('nan'))
    print(f"  Rank improvement:         {delta:+.1f}")
    if args.n_permutations > 0:
        print(f"\n  Empirical null ({args.n_permutations} permutations):")
        print(f"  Null universe size:       {metrics.get('null_universe_size', 0):.0f} genes")
        print(f"  Effective seeds used:     {metrics.get('n_effective_seeds', 0):.0f}")
        print(f"  Significant (q<0.05):     {metrics.get('n_significant_q05', 0):.0f}")
        print(f"  Significant (q<0.10):     {metrics.get('n_significant_q10', 0):.0f}")
    print(f"\n  Top 10 pathways [{args.score_variant}]:")
    for i, p in enumerate(ranked[:10], 1):
        sig = ""
        if p.empirical_q is not None:
            sig = f" q={p.empirical_q:.3f}"
        markers = ("\u2605" if p.in_chd_set else " ") + ("H" if p.degree_flag else " ")
        print(f"    {i:2d}. [{markers}][{p.sab:10s}] {p.name[:50]:50s} "
              f"{getattr(p, args.score_variant):.5f}{sig}")
    print(f"\n  Legend: \u2605=reference  H=hub (>p90 degree)")
    if args.n_permutations > 0:
        print(f"  q = BH-adjusted empirical q-value vs matched null (N={args.n_permutations})")


if __name__ == "__main__":
    main()
