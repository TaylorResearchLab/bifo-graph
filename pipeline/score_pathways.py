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
    # Empirical null model outputs (None when --n-permutations not used)
    empirical_p: Optional[float] = None
    empirical_q: Optional[float] = None
    null_mean: Optional[float] = None
    null_sd: Optional[float] = None
    null_z: Optional[float] = None

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
_NULL_SHARED: dict = {}


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
    ppr_alpha: float = 0.5,
    random_seed: int = 42,
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
        log.info("Running empirical null model (%d permutations)...", n_permutations)

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

    metrics = evaluate_pathway_recovery(scored, score_variant=score_variant)
    # Add null model metadata to metrics if available
    if n_permutations > 0 and scored and any(p.empirical_q is not None for p in scored):
        metrics['null_universe_size'] = float(_null_universe_size)
        metrics['n_effective_seeds'] = float(_n_effective_seeds)
        metrics['n_permutations_run'] = float(n_permutations)
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
                             "(default: 0 = disabled). Recommended: 1000. "
                             "Null seeds drawn from eligible pathway-connected genes.")
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
        ppr_alpha=args.ppr_alpha,
        random_seed=args.perm_random_seed,
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
