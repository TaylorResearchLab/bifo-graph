#!/usr/bin/env python3
"""
baseline_enrichment.py  —  Analysis 4 + 6 for BIFO paper

Implements three conventional enrichment baselines for comparison with
BIFO-conditioned pathway scoring, plus enhanced ranking metrics (NDCG,
recall@k, average precision) for all arms.

Baselines
---------
B0  degree_overlap    — degree-weighted seed overlap (graph membership baseline)
                        OPTIONAL: requires --kept-edges; skipped if omitted
B1  seed_fisher       — hypergeometric enrichment on seed genes only
B2  neighborhood_fisher — hypergeometric on all 1-hop gene neighbors of seeds
B3  raw_ppr_gsea      — preranked GSEA-style enrichment on raw PPR scores

BIFO result
-----------
B4  bifo_full         — BIFO-conditioned pathway scoring (degree_norm)
                        read from existing pathway_scores_full.csv

Enhanced metrics (Analysis 6)
------------------------------
Added to all arms:
  - NDCG@k           — normalized discounted cumulative gain
  - recall@k         — fraction of reference set in top-k
  - average_precision — area under precision-recall curve over ranks

Usage
-----
  python baseline_enrichment.py \\
    --edges-merged  test_output/edges_merged.csv \\
    --node-index    test_output/results_node_index.json \\
    --scores-raw    test_output/results_scores_raw.npy \\
    --scores-cond   test_output/results_scores_cond.npy \\
    --bifo-scores   test_output/pathway_scores_full.csv \\
    --chd-pathways  chd_pathway_reference.txt \\
    --seed-nodes    minimal_test_run_seed_nodes.txt \\
    --kept-edges    test_output/results_kept_edges.csv \\
    --out-csv       test_output/baseline_comparison.csv \\
    --out-json      test_output/baseline_comparison.json \\
    --min-members   8 \\
    --max-members   300
"""

import argparse

def removesuffix(s, suffix):
    return s[:-len(suffix)] if suffix and s.endswith(suffix) else s

import csv
import json
import logging
import math
import numpy as np
from collections import defaultdict
from pathlib import Path
from scipy.stats import hypergeom
from typing import Dict, List, Optional, Set, Tuple

logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
log = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Utilities
# ---------------------------------------------------------------------------

def clean_node_id(raw: str) -> str:
    return removesuffix(raw.strip(), ' CUI')


def read_node_list(path: str) -> List[str]:
    ids = []
    for line in Path(path).read_text().splitlines():
        line = line.strip()
        if not line or line.startswith('#'):
            continue
        token = line.split('#')[0].strip().split()[0] if line.split('#')[0].strip() else ''
        if token:
            ids.append(clean_node_id(token))
    return ids


def load_node_index(path: str) -> Dict[str, int]:
    raw = json.load(open(path))
    return {clean_node_id(str(k)): int(v) for k, v in raw.items()}


def load_scores(path: str) -> np.ndarray:
    return np.load(path)


# ---------------------------------------------------------------------------
# Gene universe and membership construction
# ---------------------------------------------------------------------------

MEMBERSHIP_PREDS_FWD = {
    'pathway_associated_with_gene',
    'has_signature_gene',
    'process_involves_gene',
}
MEMBERSHIP_PREDS_REV = {
    'inverse_pathway_associated_with_gene',
    'inverse_has_signature_gene',
    'gene_plays_role_in_process',
}


def build_membership_and_universe(
    edges_merged_path: str,
    seed_ids: Set[str],
    min_members: int = 8,
    max_members: int = 300,
) -> Tuple[Dict[str, Set[str]], Set[str], Set[str]]:
    """
    Returns:
        membership      : {pathway_id: set of gene_ids}
        hop1_genes      : C-prefixed 1-hop gene neighbors of seeds
        gene_universe   : all C-prefixed gene-like nodes in the graph
    """
    membership: Dict[str, Set[str]] = defaultdict(set)
    hop1_neighbors: Set[str] = set()
    with open(edges_merged_path, encoding='utf-8-sig') as f:
        reader = csv.DictReader(f)
        for row in reader:
            src = clean_node_id(row['source'])
            tgt = clean_node_id(row['target'])
            pred = row['predicate'].strip()

            # Track 1-hop neighbors of seeds
            if src in seed_ids:
                hop1_neighbors.add(tgt)
            if tgt in seed_ids:
                hop1_neighbors.add(src)

            # Build membership
            if pred in MEMBERSHIP_PREDS_FWD:
                membership[src].add(tgt)
            elif pred in MEMBERSHIP_PREDS_REV:
                membership[tgt].add(src)

    hop1_genes = {n for n in hop1_neighbors
                  if n.startswith('C') and n not in seed_ids}

    # Filter membership by size first
    membership = {
        pw: genes for pw, genes in membership.items()
        if min_members <= len(genes) <= max_members
    }

    # Gene universe = union of all pathway member genes.
    pathway_member_universe = {
        gene for genes in membership.values() for gene in genes
    }
    gene_universe = pathway_member_universe - seed_ids

    log.info("Membership: %d pathways (min=%d, max=%d members)",
             len(membership), min_members, max_members)
    log.info("1-hop gene neighbors: %d", len(hop1_genes))
    log.info("Gene universe (pathway members only): %d", len(gene_universe))
    log.info("  (excludes %d seed genes from universe denominator)", len(seed_ids))

    return membership, hop1_genes, gene_universe


# ---------------------------------------------------------------------------
# Enhanced ranking metrics (Analysis 6)
# ---------------------------------------------------------------------------

def ndcg_at_k(ranked_relevance: List[int], k: int) -> float:
    k = min(k, len(ranked_relevance))
    dcg = sum(rel / math.log2(i + 2) for i, rel in enumerate(ranked_relevance[:k]))
    n_relevant = sum(ranked_relevance)
    ideal_hits = min(n_relevant, k)
    idcg = sum(1.0 / math.log2(i + 2) for i in range(ideal_hits))
    return dcg / idcg if idcg > 0 else 0.0


def recall_at_k(ranked_relevance: List[int], k: int, n_relevant: int) -> float:
    if n_relevant == 0:
        return 0.0
    return sum(ranked_relevance[:k]) / n_relevant


def average_precision(ranked_relevance: List[int]) -> float:
    hits = 0
    precision_sum = 0.0
    for i, rel in enumerate(ranked_relevance):
        if rel:
            hits += 1
            precision_sum += hits / (i + 1)
    n_relevant = sum(ranked_relevance)
    return precision_sum / n_relevant if n_relevant > 0 else 0.0


def compute_ranking_metrics(
    ranked_items: List[str],
    reference_set: Set[str],
    k_values: Tuple[int, ...] = (10, 20, 50),
) -> Dict[str, float]:
    n_total = len(ranked_items)
    n_ref = sum(1 for item in ranked_items if item in reference_set)
    bg_rate = n_ref / n_total if n_total > 0 else 0.0

    relevance = [1 if item in reference_set else 0 for item in ranked_items]

    metrics: Dict[str, float] = {
        'n_scored': float(n_total),
        'n_reference': float(n_ref),
        'background_rate': bg_rate,
    }

    for k in k_values:
        top_k_hits = sum(relevance[:k])
        precision = top_k_hits / k
        enrichment = precision / bg_rate if bg_rate > 0 else float('nan')
        metrics[f'precision_at_{k}'] = precision
        metrics[f'enrichment_at_{k}'] = enrichment
        metrics[f'recall_at_{k}'] = recall_at_k(relevance, k, n_ref)
        metrics[f'ndcg_at_{k}'] = ndcg_at_k(relevance, k)

    ref_ranks = [i + 1 for i, r in enumerate(relevance) if r]
    metrics['mean_ref_rank'] = float(np.mean(ref_ranks)) if ref_ranks else float('nan')
    metrics['average_precision'] = average_precision(relevance)

    return metrics


# ---------------------------------------------------------------------------
# B0: Degree-weighted seed overlap  (NEW)
# ---------------------------------------------------------------------------

def compute_conditioned_degree(kept_edges_path: str) -> Dict[str, int]:
    """
    Compute out-degree of each node in the BIFO-conditioned propagation graph.
    Reads results_kept_edges.csv (all kept edges, including non-propagating).
    Out-degree is used as a graph-structure weight in B0 scoring.
    """
    degree: Dict[str, int] = defaultdict(int)
    with open(kept_edges_path, encoding='utf-8-sig') as f:
        reader = csv.DictReader(f)
        for row in reader:
            src = clean_node_id(row['source'])
            degree[src] += 1
    log.info("Conditioned degree: %d nodes with outgoing edges", len(degree))
    return dict(degree)


def degree_overlap_enrichment(
    seed_ids: Set[str],
    membership: Dict[str, Set[str]],
    conditioned_degree: Dict[str, int],
) -> List[Tuple[str, float, int, int]]:
    """
    B0: Degree-weighted seed overlap baseline.

    score(p) = Σ_{g ∈ seeds ∩ members(p)} degree_conditioned(g) / √(|members(p)|)

    Uses out-degree in the BIFO-conditioned graph as a graph-structure weight.
    Unlike B1 (Fisher), no statistical test is performed.
    Unlike B3 (GSEA), no PPR propagation is used.
    Provides a lower bound on graph-guided methods and isolates the contribution
    of propagation beyond direct seed-membership overlap.

    Returns list of (pathway_id, score, overlap_count, pathway_size),
    sorted by descending score (ties broken by descending pathway size).
    """
    results = []
    for pw_id, members in membership.items():
        overlap = members & seed_ids
        if not overlap:
            score = 0.0
        else:
            score = (
                sum(conditioned_degree.get(g, 1) for g in overlap)
                / math.sqrt(len(members))
            )
        results.append((pw_id, score, len(overlap), len(members)))

    results.sort(key=lambda x: (-x[1], -x[3]))
    return results


# ---------------------------------------------------------------------------
# B1: Seed-only hypergeometric enrichment
# ---------------------------------------------------------------------------

def seed_fisher_enrichment(
    seed_ids: Set[str],
    membership: Dict[str, Set[str]],
    gene_universe: Set[str],
) -> List[Tuple[str, float, int, int]]:
    N = len(gene_universe)
    n = len(seed_ids & gene_universe)
    results = []

    for pw_id, members in membership.items():
        K = len(members)
        N_full = len(gene_universe) + len(seed_ids)
        n_full = len(seed_ids)
        k = len(members & seed_ids)
        if k == 0:
            log_p = 0.0
        else:
            log_p = float(hypergeom.logsf(k - 1, N_full, K, n_full))
        pval = float(np.exp(log_p)) if log_p > -700 else 0.0
        results.append((pw_id, pval, k, K, log_p))

    results.sort(key=lambda x: x[4])
    return [(pw, pval, k, K) for pw, pval, k, K, _ in results]


# ---------------------------------------------------------------------------
# B2: 1-hop neighborhood hypergeometric enrichment
# ---------------------------------------------------------------------------

def neighborhood_fisher_enrichment(
    hop1_genes: Set[str],
    seed_ids: Set[str],
    membership: Dict[str, Set[str]],
    gene_universe: Set[str],
) -> List[Tuple[str, float, int, int]]:
    query = (seed_ids | hop1_genes) & gene_universe
    N = len(gene_universe)
    n = len(query)
    results = []

    for pw_id, members in membership.items():
        K = len(members)
        N_full = len(gene_universe) + len(seed_ids)
        k = len(members & query)
        if k == 0:
            log_p = 0.0
        else:
            log_p = float(hypergeom.logsf(k - 1, N_full, K, n))
        pval = float(np.exp(log_p)) if log_p > -700 else 0.0
        results.append((pw_id, pval, k, K, log_p))

    results.sort(key=lambda x: x[4])
    return [(pw, pval, k, K) for pw, pval, k, K, _ in results]


# ---------------------------------------------------------------------------
# Benjamini-Hochberg correction
# ---------------------------------------------------------------------------

def bh_correct(pvalues: List[float]) -> List[float]:
    n = len(pvalues)
    if n == 0:
        return []
    indexed = sorted(enumerate(pvalues), key=lambda x: x[1])
    adjusted = [0.0] * n
    prev = 1.0
    for rank, (orig_idx, pval) in enumerate(reversed(indexed)):
        rank_from_end = n - rank
        adj = pval * n / rank_from_end
        adj = min(adj, prev)
        prev = adj
        adjusted[orig_idx] = adj
    return adjusted


# ---------------------------------------------------------------------------
# B3: Preranked GSEA-style enrichment
# ---------------------------------------------------------------------------

def preranked_gsea_enrichment(
    gene_scores: Dict[str, float],
    membership: Dict[str, Set[str]],
    gene_universe: Set[str],
    weighted_score: bool = True,
) -> List[Tuple[str, float, float]]:
    scored_genes = [(g, s) for g, s in gene_scores.items() if g in gene_universe]
    scored_genes.sort(key=lambda x: -x[1])
    ranked_genes = [g for g, _ in scored_genes]
    ranked_scores = np.array([s for _, s in scored_genes])

    N = len(ranked_genes)
    results = []

    for pw_id, members in membership.items():
        pw_in_universe = members & gene_universe
        if len(pw_in_universe) == 0:
            continue

        hits = np.array([1 if g in pw_in_universe else 0 for g in ranked_genes])
        n_hit = hits.sum()
        if n_hit == 0:
            results.append((pw_id, 0.0, 0.0))
            continue

        if weighted_score:
            hit_weights = np.where(hits == 1, np.abs(ranked_scores), 0.0)
            hit_sum = hit_weights.sum()
            if hit_sum == 0:
                hit_weights = hits.astype(float)
                hit_sum = float(n_hit)
        else:
            hit_weights = hits.astype(float)
            hit_sum = float(n_hit)

        miss_penalty = 1.0 / (N - n_hit) if N > n_hit else 0.0

        running_sum = 0.0
        max_dev = 0.0
        max_pos = 0
        for i in range(N):
            if hits[i]:
                running_sum += hit_weights[i] / hit_sum
            else:
                running_sum -= miss_penalty
            if abs(running_sum) > abs(max_dev):
                max_dev = running_sum
                max_pos = i

        leading_edge = hits[:max_pos + 1].sum() / n_hit if n_hit > 0 else 0.0
        results.append((pw_id, float(max_dev), float(leading_edge)))

    return sorted(results, key=lambda x: -x[1])


# ---------------------------------------------------------------------------
# Main analysis
# ---------------------------------------------------------------------------

def run_baseline_comparison(
    edges_merged_path: str,
    node_index_path: str,
    scores_raw_path: str,
    scores_cond_path: str,
    bifo_scores_path: str,
    chd_pathways_path: str,
    seed_nodes_path: str,
    out_csv: str,
    out_json: str,
    min_members: int = 8,
    max_members: int = 300,
    k_values: Tuple[int, ...] = (10, 20, 50),
    allowed_pathway_ids: Optional[Set[str]] = None,
    kept_edges_path: Optional[str] = None,   # NEW: required for B0
) -> Dict:

    # Load inputs
    seed_ids = set(read_node_list(seed_nodes_path))
    chd_ref = set(read_node_list(chd_pathways_path))
    node_to_idx = load_node_index(node_index_path)
    scores_raw  = load_scores(scores_raw_path)
    scores_cond = load_scores(scores_cond_path)

    log.info("Seeds: %d  CHD reference: %d", len(seed_ids), len(chd_ref))

    # Build membership and universe
    membership, hop1_genes, gene_universe = build_membership_and_universe(
        edges_merged_path, seed_ids, min_members, max_members
    )

    # Pathway names from bifo scores
    pathway_names: Dict[str, str] = {}
    pathway_sabs:  Dict[str, str] = {}
    with open(bifo_scores_path) as f:
        for row in csv.DictReader(f):
            pathway_names[row['concept_id']] = row['name']
            pathway_sabs[row['concept_id']]  = row['sab']

    # -----------------------------------------------------------------------
    # B1: Seed-only Fisher
    # -----------------------------------------------------------------------
    log.info("Running B1: Seed-only Fisher...")
    b1_results = seed_fisher_enrichment(seed_ids, membership, gene_universe)
    b1_pvals_raw = [r[1] for r in b1_results]
    b1_pvals = bh_correct(b1_pvals_raw)
    b1_ranked = [r[0] for r in b1_results]
    b1_metrics = compute_ranking_metrics(b1_ranked, chd_ref, k_values)
    b1_metrics['method'] = 'seed_fisher'
    b1_metrics['description'] = f'Hypergeometric on {len(seed_ids)} seeds'

    # -----------------------------------------------------------------------
    # B2: 1-hop neighborhood Fisher
    # -----------------------------------------------------------------------
    log.info("Running B2: 1-hop neighborhood Fisher...")
    b2_results = neighborhood_fisher_enrichment(
        hop1_genes, seed_ids, membership, gene_universe
    )
    b2_pvals_raw = [r[1] for r in b2_results]
    b2_pvals = bh_correct(b2_pvals_raw)
    b2_ranked = [r[0] for r in b2_results]
    b2_metrics = compute_ranking_metrics(b2_ranked, chd_ref, k_values)
    b2_metrics['method'] = 'neighborhood_fisher'
    b2_metrics['description'] = (
        f'Hypergeometric on {len(hop1_genes)} 1-hop gene neighbors '
        f'(+{len(seed_ids)} seeds = {len(seed_ids | hop1_genes)} total query)'
    )

    # -----------------------------------------------------------------------
    # B3: Preranked GSEA on raw PPR scores
    # -----------------------------------------------------------------------
    log.info("Running B3: Preranked GSEA on raw PPR scores...")
    gene_scores_raw = {
        gid: float(scores_raw[node_to_idx[gid]])
        for gid in gene_universe if gid in node_to_idx
    }
    b3_results = preranked_gsea_enrichment(gene_scores_raw, membership, gene_universe)
    b3_ranked = [r[0] for r in b3_results]
    b3_metrics = compute_ranking_metrics(b3_ranked, chd_ref, k_values)
    b3_metrics['method'] = 'raw_ppr_gsea'
    b3_metrics['description'] = 'Preranked GSEA on raw PPR gene scores'

    # -----------------------------------------------------------------------
    # B3b: Preranked GSEA on conditioned PPR scores (bonus arm)
    # -----------------------------------------------------------------------
    log.info("Running B3b: Preranked GSEA on conditioned PPR scores...")
    gene_scores_cond = {
        gid: float(scores_cond[node_to_idx[gid]])
        for gid in gene_universe if gid in node_to_idx
    }
    b3b_results = preranked_gsea_enrichment(gene_scores_cond, membership, gene_universe)
    b3b_ranked = [r[0] for r in b3b_results]
    b3b_metrics = compute_ranking_metrics(b3b_ranked, chd_ref, k_values)
    b3b_metrics['method'] = 'cond_ppr_gsea'
    b3b_metrics['description'] = 'Preranked GSEA on conditioned PPR gene scores'

    # -----------------------------------------------------------------------
    # B4: BIFO full-arm (read from existing scores)
    # -----------------------------------------------------------------------
    log.info("Loading B4: BIFO full-arm scores...")
    bifo_rows = sorted(
        csv.DictReader(open(bifo_scores_path)),
        key=lambda x: float(x.get('degree_norm', 0)),
        reverse=True
    )
    b4_ranked = [r['concept_id'] for r in bifo_rows]
    bifo_universe = set(b4_ranked)
    if allowed_pathway_ids is None:
        log.info("Restricting baselines to BIFO pathway universe (%d pathways)", len(bifo_universe))
        membership = {pw: genes for pw, genes in membership.items() if pw in bifo_universe}
        log.info("Membership after BIFO-universe filter: %d pathways", len(membership))
        # Rerun baselines with matched universe
        b1_results = seed_fisher_enrichment(seed_ids, membership, gene_universe)
        b1_pvals_raw = [r[1] for r in b1_results]
        b1_pvals = bh_correct(b1_pvals_raw)
        b1_ranked = [r[0] for r in b1_results]
        b1_metrics = compute_ranking_metrics(b1_ranked, chd_ref, k_values)
        b1_metrics['method'] = 'seed_fisher'
        b1_metrics['description'] = f'Hypergeometric on {len(seed_ids)} seeds (BIFO universe)'

        b2_results = neighborhood_fisher_enrichment(hop1_genes, seed_ids, membership, gene_universe)
        b2_pvals_raw = [r[1] for r in b2_results]
        b2_pvals = bh_correct(b2_pvals_raw)
        b2_ranked = [r[0] for r in b2_results]
        b2_metrics = compute_ranking_metrics(b2_ranked, chd_ref, k_values)
        b2_metrics['method'] = 'neighborhood_fisher'
        b2_metrics['description'] = f'Hypergeometric on {len(hop1_genes)} 1-hop genes (BIFO universe)'

        gene_scores_raw = {gid: float(scores_raw[node_to_idx[gid]])
                           for gid in gene_universe if gid in node_to_idx}
        b3_results = preranked_gsea_enrichment(gene_scores_raw, membership, gene_universe)
        b3_ranked = [r[0] for r in b3_results]
        b3_metrics = compute_ranking_metrics(b3_ranked, chd_ref, k_values)
        b3_metrics['method'] = 'raw_ppr_gsea'
        b3_metrics['description'] = 'Preranked GSEA on raw PPR scores (BIFO universe)'

        gene_scores_cond = {gid: float(scores_cond[node_to_idx[gid]])
                            for gid in gene_universe if gid in node_to_idx}
        b3b_results = preranked_gsea_enrichment(gene_scores_cond, membership, gene_universe)
        b3b_ranked = [r[0] for r in b3b_results]
        b3b_metrics = compute_ranking_metrics(b3b_ranked, chd_ref, k_values)
        b3b_metrics['method'] = 'cond_ppr_gsea'
        b3b_metrics['description'] = 'Preranked GSEA on conditioned PPR scores (BIFO universe)'

    b4_metrics = compute_ranking_metrics(b4_ranked, chd_ref, k_values)
    b4_metrics['method'] = 'bifo_full'
    b4_metrics['description'] = 'BIFO-conditioned degree_norm pathway scoring'

    # -----------------------------------------------------------------------
    # B0: Degree-weighted seed overlap  (NEW — optional, requires --kept-edges)
    # -----------------------------------------------------------------------
    b0_metrics = None
    b0_results: List[Tuple[str, float, int, int]] = []
    b0_ranked: List[str] = []

    if kept_edges_path:
        log.info("Running B0: Degree-weighted seed overlap...")
        conditioned_degree = compute_conditioned_degree(kept_edges_path)
        b0_results = degree_overlap_enrichment(seed_ids, membership, conditioned_degree)
        b0_ranked = [r[0] for r in b0_results]
        b0_metrics = compute_ranking_metrics(b0_ranked, chd_ref, k_values)
        b0_metrics['method'] = 'degree_overlap'
        b0_metrics['description'] = (
            'Degree-weighted seed overlap: '
            'score(p)=sum(degree_cond(g) for g in seeds∩members(p))/sqrt(|members(p)|)'
        )
        log.info("B0 non-zero pathways: %d / %d",
                 sum(1 for r in b0_results if r[1] > 0), len(b0_results))
    else:
        log.info("B0 skipped (--kept-edges not provided)")

    # -----------------------------------------------------------------------
    # Compile comparison table
    # -----------------------------------------------------------------------
    all_methods = [b1_metrics, b2_metrics, b3_metrics, b3b_metrics, b4_metrics]
    if b0_metrics:
        all_methods.insert(0, b0_metrics)   # prepend — B0 is the simplest baseline

    print("\n" + "=" * 90)
    print("BASELINE COMPARISON (Analysis 4 + 6)")
    print("=" * 90)
    print(f"  Pathway universe: {len(membership)} pathways  "
          f"(min={min_members}, max={max_members} members)")
    print(f"  CHD reference: {len(chd_ref)} pathway IDs  "
          f"(background rate: {len([r for r in b4_ranked if r in chd_ref])/len(b4_ranked):.3f})")
    print()

    header = f"  {'Method':25s}  {'P@10':>6s}  {'P@20':>6s}  {'P@50':>6s}  "
    header += f"{'R@10':>6s}  {'NDCG@10':>8s}  {'AP':>6s}  {'MeanRank':>9s}"
    print(header)
    print("  " + "-" * 87)

    for m in all_methods:
        row = (f"  {m['method']:25s}  "
               f"{m.get('precision_at_10',0):>6.3f}  "
               f"{m.get('precision_at_20',0):>6.3f}  "
               f"{m.get('precision_at_50',0):>6.3f}  "
               f"{m.get('recall_at_10',0):>6.3f}  "
               f"{m.get('ndcg_at_10',0):>8.4f}  "
               f"{m.get('average_precision',0):>6.4f}  "
               f"{m.get('mean_ref_rank',float('nan')):>9.1f}")
        print(row)

    print()
    print("  Top 5 pathways per method:")

    # Build top-5 print list — include B0 if present
    top5_methods = []
    if b0_metrics:
        top5_methods.append(('B0 degree_overlap', b0_ranked, None))
    top5_methods += [
        ('B1 seed_fisher',         b1_ranked, b1_pvals),
        ('B2 neighborhood_fisher', b2_ranked, b2_pvals),
        ('B3 raw_ppr_gsea',        b3_ranked, None),
        ('B3b cond_ppr_gsea',      b3b_ranked, None),
        ('B4 bifo_full',           b4_ranked, None),
    ]

    for label, ranked, pvals in top5_methods:
        print(f"\n  {label}:")
        for i, pw_id in enumerate(ranked[:5], 1):
            chd = '★' if pw_id in chd_ref else ' '
            name = pathway_names.get(pw_id, pw_id)[:50]
            pv_str = f"  BH_p={pvals[i-1]:.2e}" if pvals else ""
            print(f"    {i}. [{chd}] {name}{pv_str}")

    # -----------------------------------------------------------------------
    # Write outputs
    # -----------------------------------------------------------------------
    csv_rows = []
    b1_bh_map  = {b1_results[i][0]: b1_pvals[i]  for i in range(len(b1_results))}
    b2_bh_map  = {b2_results[i][0]: b2_pvals[i]  for i in range(len(b2_results))}

    # Build the methods list for CSV output — include B0 if present
    csv_methods = []
    if b0_results:
        csv_methods.append(
            ('degree_overlap', b0_ranked,
             {r[0]: r for r in b0_results})
        )
    csv_methods += [
        ('seed_fisher',         b1_ranked, {r[0]: r for r in b1_results}),
        ('neighborhood_fisher', b2_ranked, {r[0]: r for r in b2_results}),
        ('raw_ppr_gsea',        b3_ranked, {r[0]: r for r in b3_results}),
        ('cond_ppr_gsea',       b3b_ranked, {r[0]: r for r in b3b_results}),
        ('bifo_full',           b4_ranked,
         {r['concept_id']: (r['concept_id'], r.get('degree_norm',''),
                             r.get('member_gene_count',''), '')
          for r in bifo_rows}),
    ]

    for method_label, ranked, raw_results in csv_methods:
        for rank, pw_id in enumerate(ranked, 1):
            r = raw_results.get(pw_id, (pw_id, '', '', ''))
            if method_label == 'seed_fisher':
                raw_p = r[1] if len(r) > 1 and r[1] != '' else ''
                bh_p  = b1_bh_map.get(pw_id, '')
                stat_val = f"{float(raw_p):.6e}" if raw_p != '' else ''
                bh_val   = f"{float(bh_p):.6e}" if bh_p != '' else ''
                log10_bh = f"{-math.log10(max(float(bh_p), 1e-300)):.3f}" if bh_p != '' and float(bh_p) > 0 else 'inf'
            elif method_label == 'neighborhood_fisher':
                raw_p = r[1] if len(r) > 1 and r[1] != '' else ''
                bh_p  = b2_bh_map.get(pw_id, '')
                stat_val = f"{float(raw_p):.6e}" if raw_p != '' else ''
                bh_val   = f"{float(bh_p):.6e}" if bh_p != '' else ''
                log10_bh = f"{-math.log10(max(float(bh_p), 1e-300)):.3f}" if bh_p != '' and float(bh_p) > 0 else 'inf'
            elif method_label == 'degree_overlap':
                # B0: stat = raw score, no BH p-value
                stat_val = f"{float(r[1]):.6f}" if len(r) > 1 and r[1] != '' else ''
                bh_val   = ''
                log10_bh = ''
            else:
                stat_val = f"{float(r[1]):.6f}" if len(r) > 1 and r[1] != '' else ''
                bh_val = ''
                log10_bh = ''
            csv_rows.append({
                'method': method_label,
                'rank': rank,
                'concept_id': pw_id,
                'name': pathway_names.get(pw_id, pw_id),
                'sab': pathway_sabs.get(pw_id, ''),
                'in_chd_set': pw_id in chd_ref,
                'stat': stat_val,
                'bh_adjusted_p': bh_val,
                'neg_log10_bh_p': log10_bh,
                'overlap_or_es': r[2] if len(r) > 2 else '',
                'pathway_size': r[3] if len(r) > 3 else '',
            })

    with open(out_csv, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=list(csv_rows[0].keys()))
        writer.writeheader()
        writer.writerows(csv_rows)
    log.info("Wrote %d rows to %s", len(csv_rows), out_csv)

    # JSON: summary metrics
    output = {
        'parameters': {
            'min_members': min_members,
            'max_members': max_members,
            'n_pathways_universe': len(membership),
            'n_chd_reference': len(chd_ref),
            'n_seeds': len(seed_ids),
            'n_hop1_genes': len(hop1_genes),
            'n_gene_universe': len(gene_universe),
            'b0_enabled': kept_edges_path is not None,
        },
        'methods': all_methods,
    }
    Path(out_json).write_text(json.dumps(output, indent=2))
    log.info("Wrote metrics to %s", out_json)

    return output


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="BIFO Analysis 4+6: baseline enrichment comparison"
    )
    parser.add_argument("--edges-merged",  required=True)
    parser.add_argument("--node-index",    required=True)
    parser.add_argument("--scores-raw",    required=True)
    parser.add_argument("--scores-cond",   required=True)
    parser.add_argument("--bifo-scores",   required=True,
                        help="pathway_scores_full.csv from score_pathways.py")
    parser.add_argument("--chd-pathways",  required=True)
    parser.add_argument("--seed-nodes",    required=True)
    parser.add_argument("--out-csv",       required=True)
    parser.add_argument("--out-json",      required=True)
    parser.add_argument("--min-members",   type=int, default=8)
    parser.add_argument("--max-members",   type=int, default=300)
    parser.add_argument("--kept-edges",    default=None,
                        help="results_kept_edges.csv from bifo_conditioning.py. "
                             "Required for B0 (degree-weighted seed overlap). "
                             "If omitted, B0 is skipped and all other baselines run normally.")
    args = parser.parse_args()

    run_baseline_comparison(
        edges_merged_path=args.edges_merged,
        node_index_path=args.node_index,
        scores_raw_path=args.scores_raw,
        scores_cond_path=args.scores_cond,
        bifo_scores_path=args.bifo_scores,
        chd_pathways_path=args.chd_pathways,
        seed_nodes_path=args.seed_nodes,
        out_csv=args.out_csv,
        out_json=args.out_json,
        min_members=args.min_members,
        max_members=args.max_members,
        kept_edges_path=args.kept_edges,
    )


if __name__ == "__main__":
    main()
